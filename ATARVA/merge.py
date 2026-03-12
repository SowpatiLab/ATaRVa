from math import sqrt, modf
import sys, os
import pysam
import timeit as ti
# import argparse as ap
from multiprocessing import Process
from ATARVA.tamatr import reader

def merge_parser(subparsers):
    parser = subparsers.add_parser("merge", help="merging multiple ATaRVa VCF files over specified regions", description="Merge ATaRVa VCF files")
    parser._action_groups.pop()


    required = parser.add_argument_group('Required arguments')
    required.add_argument('-r', '--regions', required=True, metavar='<FILE>', help='input regions file. the regions file should be strictly in bgzipped tabix format. \
                                                                  If the regions input file is in bed format. First sort it using bedtools. Compress it using bgzip. \
                                                                  Index the bgzipped file with tabix command from samtools package.')
    required.add_argument('-i', '--vcfs', nargs='+', required=True, metavar='<FILE>', help='text file containing paths to input vcf files to be merged. The text file should list each path on a separate line. The vcf files should be strictly in bgzipped tabix format. \
                                                                  If the vcfs input file is in vcf format. First sort it using bcftools. Compress it using bgzip. \
                                                                  Index the bgzipped file with tabix command from samtools package.')

    required.add_argument('-f', '--fasta', required=True, metavar='<FILE>', help='input reference fasta file. The file should be indexed.')
    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('--contigs', nargs='+', help='contigs to get merged [chr1 chr12 chr22 ..]. If not mentioned every contigs in the region file will be merged.')
    optional.add_argument('-o', '--outname', type=str, metavar='<STR>', default='', help='name of the output file, output is in vcf format.')
    optional.add_argument('-t',  '--threads', type=int, metavar='<INT>', default=1, help='number of threads. [default: 1]')

    if (len(sys.argv) == 2) and (sys.argv[1] == 'merge'):
        parser.print_help()
        sys.exit()

    parser.set_defaults(func=merge_run)

def shredder(threads):
    sq_threads = sqrt(threads)
    decimal_val, integer_value = modf(sq_threads)
    integer_value = int(integer_value)
    if decimal_val>0:
        thread_list = [integer_value]*(integer_value+1)
        remaining_threads = threads-sum(thread_list)
        if remaining_threads>1:
            while remaining_threads>0:
                for i in range(len(thread_list[1:])):
                    thread_list[i+1]+=1
                    remaining_threads-=1
                    if remaining_threads<=0:
                        break
        else:
            thread_list[-1] += remaining_threads
    else:
        thread_list = [integer_value] + [integer_value-1]*(integer_value)
        remaining_threads = threads-sum(thread_list)
        thread_list[-1] += remaining_threads

    return thread_list

def merge_run(args):
    start_time = ti.default_timer()
    # args = parse_args()

    for arg in vars(args):
        if arg in ['func', 'command']: continue
        print (arg, getattr(args, arg))
    print('\n')
    
    out_file = sys.stdout
    if args.outname:
        if '.vcf' == args.outname[-4:]:
            out_file = f'{args.outname}'[:-4]
        elif args.outname[-1]=='/':
            out_file = args.outname + "atarva_merged"
        else:
            out_file = f'{args.outname}'
    else:
        out_file = "atarva_merged"

    vcf_list = []
    if len(args.vcfs)==1:
        with open(args.vcfs[0], 'r') as vh:
            for line in vh:
                if line[0]=='#': continue
                line = line.strip()
                vcf_list.append(line)
    else:
       vcf_list = args.vcfs 
    
    tbx  = pysam.Tabixfile(args.regions)
    total_loci = 0
    if not args.contigs:
        contigs = sorted(tbx.contigs)
        for row in tbx.fetch():
            total_loci += 1
    else:
        contigs = sorted(args.contigs)
        for each_contig in contigs:
            for row in tbx.fetch(each_contig):
                total_loci += 1
    
    # print('total_loci = ', total_loci)

    threads = args.threads
    threads = threads - 1

    
    split_point = 5000 if total_loci > 5000 else total_loci // 5
    if split_point == 0:
        split_point = 1
        partition_point = 1
    else:
        partition_point = total_loci//split_point

    split_point_chunks = 0 # to count number of split_point chunks excluding the 'minimum chunks' eg 9920 from 1 contig and 80 from another contig to add up to 10000
    fetcher = []
    line_count = 0
    current_split = []
    for each_contig in contigs:
        init = 0
        for row in tbx.fetch(each_contig):
            line_count += 1
            if init == 0:
                Row = row.split('\t')
                chrom = Row[0]
                start_coord = (int(Row[1]), int(Row[2]))
                init=1
            if split_point_chunks < partition_point-1:
                if line_count % split_point == 0:
                    end_coord = (int(row.split('\t')[1]), int(row.split('\t')[2]))
                    current_split.append([chrom, start_coord, end_coord])
                    # fetcher.append(tuple(current_split))
                    fetcher.extend(current_split)
                    split_point_chunks += 1
                    line_count = 0
                    current_split = []
                    init = 0
        if init != 0:
            end_coord = (int(row.split('\t')[1]), int(row.split('\t')[2]))
            current_split.append([chrom, start_coord, end_coord])
    # fetcher.append(tuple(current_split))
    fetcher.extend(current_split)
    tbx.close()

    # print('Length of fetcher = ', len(fetcher))
    # print('partition_point for fetcher = ', partition_point)
    # fetcher = fetcher[:2]

    region_file = args.regions
    ref_file = args.fasta
    if threads > 1:

        thread_list = shredder(threads)
        reader_thread_pool = []
        reader_threads = thread_list[0]
        partition = len(fetcher) // reader_threads
        # print('reader_threads = ', reader_threads)
        # print('partition for reading = ', partition)
        initial = 0
        track = partition
        for each_reader_thread in range(reader_threads):
            if each_reader_thread+1 == reader_threads:
                reader_contigs = fetcher[initial : ]
            else:
                reader_contigs = fetcher[initial : track]
            # print('Thread = ', each_reader_thread, ' contig length = ', len(reader_contigs))    
            thread_x = Process(target = reader, args = (out_file, region_file, ref_file, vcf_list, reader_contigs, each_reader_thread, thread_list[each_reader_thread+1]))
            thread_x.start()
            reader_thread_pool.append(thread_x)
            
            initial = track
            track += partition
    
        # joining Threads 
        for thread_x in reader_thread_pool:
            thread_x.join()
        # emptying thread_pool
        reader_thread_pool.clear()
        #sys.exit()
        out = open(f'{out_file}.vcf', 'a')
        print('Concatenating thread outputs!', file=sys.stderr)
        for tidx in range(reader_threads)[1:]:
            thread_out = f'{out_file}_thread_{tidx}.vcf'
            with open(thread_out, 'r') as fh:
                # if tidx!=0: next(fh)
                for line in fh:
                    repeat_info = line.strip().split('\t')
                    #print(*repeat_info, file=out, sep='\t')
                    out.write("\t".join(map(str, repeat_info)) + "\n")
            os.remove(thread_out)
        out.close()
        print('Concatenation completed!! ^_^', file=sys.stderr)

    else:
        reader(out_file, region_file, ref_file, vcf_list, fetcher, -1, 0)

    time_now = ti.default_timer()
    sys.stderr.write('CPU time: {} seconds\n'.format(time_now - start_time))