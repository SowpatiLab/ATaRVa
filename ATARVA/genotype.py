import sys, os, gzip
import pysam
import timeit as ti
# import argparse as ap
from multiprocessing import Process, Pool
from tqdm import tqdm

from ATARVA.version import __version__
from ATARVA.baseline import *
# from ATARVA.vcf_writer import set_info_opt_tag
from ATARVA.sub_operation_utils import set_methviz_tag

def genotype_parser(subparsers):
    """
    Parse command line arguments.
    """
    parser = subparsers.add_parser("genotype", help="tandem repeat genotyper specially designed for long read data", description="Tandem Repeat Genotyper")
    parser._action_groups.pop()

    required = parser.add_argument_group('Required arguments')
    required.add_argument('-f',  '--fasta',   required=True, metavar='<FILE>', help='input reference fasta file')
    required.add_argument('-b', '--bam', nargs='+', required=True, metavar='<FILE>', help='samples alignment files. allowed formats: SAM, BAM, CRAM')
    required.add_argument('-r', '--regions', required=True, metavar='<FILE>', help='input regions file. the regions file should be strictly in bgzipped tabix format. \
                                                                  If the regions input file is in bed format. First sort it using bedtools. Compress it using bgzip. \
                                                                  Index the bgzipped file with tabix command from samtools package.')

    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('--format', type=str, metavar='<STR>', default='bam', help='format of input alignment file. allowed options: [cram, bam, sam]. default: [bam]')
    optional.add_argument('-q', '--map-qual', type=int, metavar='<INT>', default=5, help='minimum mapping quality of the reads to be considered. [default: 5]')
    optional.add_argument('--contigs', nargs='+', help='contigs to get genotyped [chr1 chr12 chr22 ..]. If not mentioned every contigs in the region file will be genotyped.')
    optional.add_argument('--min-reads', type=int, metavar='<INT>', default=10, help='minimum read coverage after quality cutoff at a locus to be genotyped. [default: 10]')
    optional.add_argument('--max-reads', type=int, metavar='<INT>', default=None, help='maximum number of reads to be used for genotyping a locus. [default: 100]')
    optional.add_argument('--snp-dist', type=int, metavar='<INT>', default=3000, help='maximum distance of the SNP from repeat region to be considered for phasing. [default: 3000]')
    optional.add_argument('--snp-count', type=int, metavar='<INT>', default=3, help='number of SNPs to be considered for phasing (minimum value = 1). [default: 3]')
    optional.add_argument('--snp-qual', type=int, metavar='<INT>', default=20, help='minimum basecall quality at the SNP position to be considered for phasing. [default: 20]')
    optional.add_argument('--flank', type=int, metavar='<INT>', default=None, help='length of the flanking region (in base pairs) to search for insertion with a repeat in it. [default: 10]')
    optional.add_argument('--snp-read', type=float, metavar='<FLOAT>', default=0.2, help='a positive float as the minimum fraction of snp\'s read contribution to be used for phasing. [default: 0.25]')
    optional.add_argument('--meth-prob', type=float, metavar='<FLOAT>', default=0.5, help='a minimum probability cutoff for methylation. [default: 0.5]')
    optional.add_argument('--phasing-read', type=float, metavar='<FLOAT>', default=0.4, help='a positive float as the minimum fraction of total read contribution from the phased read clusters. [default: 0.4]')
    optional.add_argument('-o',  '--vcf', type=str, metavar='<FILE>', default='', help='name of the output file, output is in vcf format. [default: sys.stdout]')
    optional.add_argument('--karyotype', nargs='+', help='karyotype of the samples [XY XX]')
    optional.add_argument('-t',  '--threads', type=int, metavar='<INT>', default=1, help='number of threads. [default: 1]')
    optional.add_argument('--haplotag', type=str, metavar='<STR>', default=None, help='use haplotagged information for phasing. eg: [HP]. [default: None]')
    optional.add_argument('--decompose', action='store_true', help="write the motif-decomposed sequence to the vcf. [default: False]")
    optional.add_argument('--methviz', action='store_true', help="write the methylation encoded sequence to the vcf for visualization purpose. [default: False]")
    optional.add_argument('--amplicon', action='store_true', help="genotype mode for targeted-sequenced samples. In this mode, the default values for `max-reads` and `flank` values are 1000 and 20 respectively. [default: False]")
    optional.add_argument('--somatic', action='store_true', help="genotype mode for capturing mosaicism in samples. In this mode, default `max-reads` and `flank` values are same as amplicon mode. [default: False]")
    optional.add_argument('--read-wise', action='store_true', help="Read-wise genotyping mode for BED file with dense regions. [default: False]")
    optional.add_argument('--loci-wise', action='store_true', help="Loci-wise genotyping mode instead of Read-wise for BED file with sparse regions. [default: False]")
    optional.add_argument('-log', '--debug_mode', action='store_true', help="write the debug messages to log file. [default: False]")
    optional.add_argument('-v', '--version', action='version', version=f'ATaRVa version {__version__}')
    

    

    if (len(sys.argv) == 2) and (sys.argv[1] == 'genotype'):
        parser.print_help()
        sys.exit()
    
    parser.set_defaults(func=genotype_run)

def f_check(path):
    """
    Check if the provided FASTA file is valid.
    
    Args:
        path (str): Path to the FASTA file.
    
    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If the file is not a valid FASTA file.
        OSError: If there is an error reading the file.
    """
    try:
        f = pysam.FastaFile(path)
        f.close()
    except (FileNotFoundError, ValueError, OSError) as e:
        print(f"Error: {path} is not a valid FASTA file. {str(e)}")
        sys.exit()
    except Exception as e:
        print("An unexpected error occurred:", str(e))
        sys.exit()

def b_check(path, aln_format):
    """
    Check if the provided BAM file is valid and sorted by coordinate.
    
    Args:
        path (str): Path to the BAM file.
        aln_format (str): Format of the alignment file.
    
    Raises:
        FileNotFoundError: If the file does not exist.
        SortOrderError: If the BAM file is not sorted by coordinate.
        ValueError: If the file is not a valid BAM file.
        OSError: If there is an error reading the file.
    """

    try:
        b = pysam.AlignmentFile(path, aln_format)
        header = b.header
        if 'HD' in header and 'SO' in header['HD']:
            sort_order = header['HD']['SO']
            if sort_order == 'coordinate':
                pass
                # print(f"Alignment file sort order: {sort_order}")
            else:
                print(f"Alignment file sort order: {sort_order}. It should be sorted by \'coordinate\'!!")
                print(f"Use: samtools sort sorted_{path.split('/')[-1]} {path.split('/')[-1]}")
                sys.exit()
        else:
            print("No sort order specified in the header.")
            print(f"Use: samtools sort sorted_{path.split('/')[-1]} {path.split('/')[-1]}")
            sys.exit()

        b.close()
    except (FileNotFoundError, ValueError, OSError) as e:
        print(f"Error: {path} is not a valid alignment file. {str(e)}")
        sys.exit()
    except Exception as e:
        print("An unexpected error occurred:", str(e))
        sys.exit()

def t_check(path):
    """
    Check if the provided regions file is valid and indexed.
    
    Args:
        path (str): Path to the regions file.
        
    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If the file is not a valid tabix file.
        OSError: If there is an error reading the file.
    """
    try:
        columns = ["Chrom", "Start", "End", "Motif", "Motif length"]
        incorrect_symbols = {' ', ',', ';', ':', '|', '!', '@', '#', '$', '%', '^', '&', '*', '(', ')', '+', '=', '[', ']', '{', '}', '<', '>', '?', '/', '\\'}
        t = pysam.TabixFile(path)
        for rows in t.fetch():
            rows = rows.strip().split('\t')
            break
        if len(rows) < 5:
            print(f"Error: {path} should have at least 5 columns (chrom, start, end, motif, motif length)!!")
            sys.exit()
        else:
            compt_list = [rows[0].isalnum(), rows[1].isdigit(), rows[2].isdigit(), rows[3].isalpha(), rows[4].isdigit()]
            if not all(compt_list):
                incorrect_cols = [columns[index] for index,i in enumerate(compt_list) if not i]
                print(f"Error: The first 5 columns of {path} should be chrom, start, end, motif, motif length respectively!!")
                print(f"Incorrect columns: {', '.join(incorrect_cols)}")
                sys.exit()
            elif (len(rows) > 5) and ( len(incorrect_symbols & set(rows[5])) > 0):
                print(f"Error: The 6th column of {path} should be a alphanumeric string representing the optional annotation!!\n eg: HTT, FMR1")
                sys.exit()

        t.close()
    except (FileNotFoundError, ValueError, OSError) as e:
        print(f"Error: {path} is not a valid tabix file. {str(e)}")
        sys.exit()
    except Exception as e:
        print("An unexpected error occurred:", str(e))
        sys.exit()

def genotype_run(args):

    start_time = ti.default_timer()
    # args = parse_args()

    for arg in vars(args):
        if arg in ['func', 'help', 'command']: continue
        print (arg, getattr(args, arg))
    print('\n')

    f_check(args.fasta)

    aln_format = ''         # format of the alignment file
    if   args.format == 'cram': aln_format = 'rc'
    elif args.format == 'sam':  aln_format = 'r'
    else:            aln_format = 'rb'

    for each_bam in args.bam:
        b_check(each_bam, aln_format)
    t_check(args.regions)


    out_file = sys.stdout
    if args.vcf:
        if '.vcf' == args.vcf[-4:]:
            out_file = f'{args.vcf}'[:-4]
        elif args.vcf[-1]=='/':
            out_file = args.vcf + "atarva_tr"
        else:
            out_file = f'{args.vcf}'
    # else: out_file = f'{".".join(args.bams.split(".")[:-1])}'
    external_name = out_file

    set_info_mp_cutoff(args.meth_prob)

    set_methviz_tag(args.methviz)

    tbx  = pysam.Tabixfile(args.regions)
    total_loci = 0
    if not args.contigs:
        contigs = sorted(tbx.contigs)
        for row in tbx.fetch():
            total_loci += 1
    else:
        contigs = sorted(args.contigs)
        for each_contig in args.contigs:
            for row in tbx.fetch(each_contig):
                total_loci += 1

    if not args.karyotype:
        karyotype_list = [False]*len(args.bam)
    else:
        karyotype_list = [i=='XY' for i in args.karyotype]

    maxR = args.max_reads
    flank_length = args.flank
    if args.amplicon or args.somatic:
        if args.max_reads is None: maxR = 1000
        if args.flank is None: flank_length = 20
    else:
        if maxR is None: maxR = 100
        if flank_length is None: flank_length = 10

    threads = args.threads
    split_point = total_loci // threads # number of loci to be handled by each thread. the last thread will handle the remaining loci if total_loci is not perfectly divisible by threads
    # split_point is 0 when the total_loci is less than the number of threads
    # this is a rare case with a bed file with very few loci; all these loci will be handled by a single thread
    if split_point == 0:
        split_point = 1
        threads = 1

    
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
            if len(fetcher) < threads-1:
                if line_count % split_point == 0:
                    end_coord = (int(row.split('\t')[1]), int(row.split('\t')[2]))
                    current_split.append([chrom, start_coord, end_coord])
                    fetcher.append(tuple(current_split))
                    line_count = 0
                    current_split = []
                    init = 0
        if init != 0:
            end_coord = (int(row.split('\t')[1]), int(row.split('\t')[2]))
            current_split.append([chrom, start_coord, end_coord])
    fetcher.append(tuple(current_split))
    tbx.close()

    mbso = 0
    if (len(args.bam)>1) and (args.vcf):
        mbso = 1
    
    for kidx, each_bam in enumerate(args.bam):
        out_file = external_name
        print(f"Processing sample {each_bam.split('/')[-1]}\n")

        count = 0
        aln_file = pysam.AlignmentFile(each_bam, aln_format, reference_filename=args.fasta)
        length = 0
        for read in aln_file.fetch():
            if (read.flag & 0X400) or (read.flag & 0X100): continue 
            count+=1
            string = read.cigarstring
            length += read.query_length
            if read.has_tag('cs'):
                print("CS tag detected. Processing using CS tag...\n")
                break
            elif (string!=None) and (('X' in string) or ('=' in string)):
                print("CIGAR(X/=) tag detected. Processing using CIGAR(X/=) tag...\n")
                break
            elif read.has_tag('MD'):
                print("MD tag detected. Processing using MD tag...")
                print("Include CS tag or CIGAR tag with 'X/=' for faster processing.\n")
                break
            if count>100:
                print(f"No tags detected in {each_bam.split('/')[-1]}. Processing without tags...")
                print("Include the CS tag, MD tag, or CIGAR tag with 'X/=' for faster processing.\n")
                break
                # sys.exit()
        aln_file.close()

        amplicon = False
        somatic = False
        if args.amplicon:
            amplicon = True
            srs = True
            print('Processing in amplicon mode...')
        elif args.somatic:
            srs = True
            somatic = True
            print('Processing in somatic mode...')
            amplicon = True
        elif (args.read_wise and args.loci_wise):
            print('Error: Choose either Read-wise or Loci-wise genotyping mode!!')
            sys.exit()
        elif args.loci_wise:
            srs = True
            print('Processing in Loci-wise genotyping mode...')
        else:
            srs = False
            print('Processing in Read-wise genotyping mode...')

        if not args.vcf:
            out_file = f'{".".join(each_bam.split("/")[-1].split(".")[:-1])}'
        elif mbso or (out_file[-1]=='/'):
            out_file = out_file + '_' + ".".join(each_bam.split("/")[-1].split('.')[:-1])

        if threads > 1:
            def update(_):
                pbar.update()
            pbar = tqdm(total=threads, desc="Processing ", ascii="_>", ncols=75, bar_format="{l_bar}{bar}{n_fmt}/{total_fmt}")
            with Pool(processes=threads) as pool:
                for tidx in range(threads):
                    contig = fetcher[tidx]
                    if srs:
                        pool.apply_async(mini_cooper, 
                                         args=(each_bam, args.regions, args.fasta, aln_format, contig, args.map_qual, out_file, args.snp_qual, args.snp_count, args.snp_dist, maxR, args.min_reads, args.snp_read, args.phasing_read, tidx, flank_length, args.debug_mode, karyotype_list[kidx], args.decompose, args.haplotag, amplicon, args.meth_prob, somatic),
                                         callback=update)
                    else:
                        pool.apply_async(cooper, 
                                         args=(each_bam, args.regions, args.fasta, aln_format, contig, args.map_qual, out_file, args.snp_qual, args.snp_count, args.snp_dist, maxR, args.min_reads, args.snp_read, args.phasing_read, tidx, flank_length, args.debug_mode, karyotype_list[kidx], args.decompose, args.haplotag, amplicon, args.meth_prob, somatic),
                                         callback=update)
                pool.close()
                pool.join()
            pbar.close()
        
            out = open(f'{out_file}.vcf', 'a')
            print('Concatenating thread outputs!', file=sys.stderr)
            idx = out_file.rfind('/')
            hid_outfile = out_file[:idx+1] + '.' + out_file[idx+1:]
            for tidx in range(threads)[1:]:
                thread_out = f'{hid_outfile}_thread_{tidx}.vcf'
                with open(thread_out, 'r') as fh:
                    # if tidx!=0: next(fh)
                    for line in fh:
                        repeat_info = line.strip().split('\t')
                        print(*repeat_info, file=out, sep='\t')
                os.remove(thread_out)
            out.close()
            print('Concatenation completed!! ^_^', file=sys.stderr)

            if args.debug_mode:
                out_log = open(f'{out_file}_debug.log', 'a')
                for tidx in range(threads)[1:]:
                    thread_log_out = f'{hid_outfile}_debug_{tidx}.log'
                    with open(thread_log_out, 'r') as fh:
                        for line in fh:
                            log_info = line.strip()
                            print(log_info, file=out_log)
                    os.remove(thread_log_out)
                out_log.close()
        else:
            if srs:
                mini_cooper(each_bam, args.regions, args.fasta, aln_format, fetcher[0], args.map_qual, out_file, args.snp_qual, args.snp_count, args.snp_dist, maxR, args.min_reads, args.snp_read, args.phasing_read, -1, flank_length, args.debug_mode, karyotype_list[kidx], args.decompose, args.haplotag, amplicon, args.meth_prob, somatic)
            else:
                cooper(each_bam, args.regions, args.fasta, aln_format, fetcher[0], args.map_qual, out_file, args.snp_qual, args.snp_count, args.snp_dist, maxR, args.min_reads, args.snp_read, args.phasing_read, -1, flank_length, args.debug_mode, karyotype_list[kidx], args.decompose, args.haplotag, amplicon, args.meth_prob, somatic)

    time_now = ti.default_timer()
    sys.stderr.write('CPU time: {} seconds\n'.format(time_now - start_time))