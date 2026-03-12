import sys, os, gzip
import pysam
import timeit as ti
# import argparse as ap
from multiprocessing import Process, Pool
from tqdm import tqdm

from ATARVA.readers import fasta_check, bam_check, tabix_check
from ATARVA.version import __version__
from ATARVA.baseline import Cooper
from ATARVA.readers import *
from ATARVA.vcf_writer import set_info_mp_cutoff
from ATARVA.sub_operation_utils import set_methviz_tag

def genotype_parser(subparsers):
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
    optional.add_argument('-q', '--map-qual', type=int, metavar='<INT>', default=5, help='minimum mapping quality of a read to be considered. [default: 5]')
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
    optional.add_argument('--read-wise', action='store_true', help="Read-wise genotyping mode for BED file with dense regions. [default: False]")
    optional.add_argument('--loci-wise', action='store_true', help="Loci-wise genotyping mode instead of Read-wise for BED file with sparse regions. [default: False]")
    optional.add_argument('-log', '--debug_mode', action='store_true', help="write the debug messages to log file. [default: False]")
    optional.add_argument('-v', '--version', action='version', version=f'ATaRVa version {__version__}')


    if (len(sys.argv) == 2) and (sys.argv[1] == 'genotype'):
        parser.print_help()
        sys.exit()
    
    parser.set_defaults(func=genotype_run)


def splitfile_threads(tbx, contigs, total_loci, threads):
    """
    splits the input region file into thread-wise chunks for parallel processing

    :param tbx: pysam Tabixfile object of the input region file
    :param contigs: list of contigs to be processed
    :param total_loci: total number of loci in the input region file for the specified contigs
    :param threads: number of threads for parallel processing

    :return: list of region chunks for each thread
    """

    split_point = total_loci // threads
    # split_point is 0 when the total_loci is less than the number of threads
    # this is a rare case with a bed file with very few loci; all these loci will be handled by a single thread
    if split_point == 0: split_point = total_loci
    
    fetcher = []
    line_count = 0
    current_split = []
    for each_contig in contigs:
        init = 0
        for row in tbx.fetch(each_contig):
            line_count += 1
            if init == 0:
                fields = row.split('\t')
                chrom  = fields[0]
                start_coord = (int(fields[1]), int(fields[2]))
                init = 1
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

    return fetcher


def worker_init(bam_file, region_ranges, args, out_file, sample_idx, thread_idx):
    """
    Initializes the Cooper class for genotyping in each thread.
    
    :param bam_file: path to the alignment file
    :param region_ranges: list of regions to be processed by the thread
    :param args: command line arguments
    :param out_file: output file name
    :param sample_idx: index of the sample being processed
    :param thread_idx: index of the thread
    """
    cooper = Cooper(bam_file, region_ranges, args, out_file, sample_idx, thread_idx)


def genotype_run(args):

    start_time = ti.default_timer()
    # args = parse_args()

    for arg in vars(args):
        if arg in ['func', 'help', 'command']: continue
        print (arg, getattr(args, arg))
    print('\n')

    fasta_check(args.fasta)
    tabix_check(args.regions)

    if   args.format == 'cram': args.format = 'rc'
    elif args.format == 'sam':  args.format = 'r'
    elif args.format == 'bam':  args.format = 'rb'

    for each_bam in args.bam:
        bam_check(each_bam, args.format)

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

    temp_tbx  = pysam.Tabixfile(args.regions)
    total_loci = 0
    if not args.contigs:
        args.contigs = sorted(temp_tbx.contigs)
        for row in temp_tbx.fetch(): total_loci += 1
    else:
        args.contigs = sorted(args.contigs)
        total_loci = sum(1 for contig in args.contigs for row in temp_tbx.fetch(contig))

    if not args.karyotype:
        args.karyotype = [False]*len(args.bam)
    else:
        args.karyotype = [i=='XY' for i in args.karyotype]

    if args.amplicon:
        if args.max_reads is None: args.max_reads = 1000
        if args.flank is None: args.flank = 20
    else:
        if args.max_reads is None: args.max_reads = 100
        if args.flank is None: args.flank = 10

    split_point = total_loci // args.threads
    if split_point == 0:
        split_point = 1
        args.threads = 1

    fetcher = splitfile_threads(temp_tbx, args.contigs, total_loci, args.threads)
    temp_tbx.close()

    mbso = False    # muliple bams single output
    if len(args.bam) > 1 and (args.vcf):
        mbso = True

    for sample_idx, bam_file in enumerate(args.bam):
        out_file = external_name
        print(f"Processing sample {bam_file.split('/')[-1]}\n")
        
        srs = check_alnfile(bam_file, args)

        if not args.vcf:
            out_file = f'{".".join(bam_file.split("/")[-1].split(".")[:-1])}'
        elif mbso or (out_file[-1]=='/'):
            out_file = out_file + '_' + ".".join(bam_file.split("/")[-1].split('.')[:-1])

        if args.threads > 1:
            def update(_):
                pbar.update()
            pbar = tqdm(total=args.threads, desc="Processing ", ascii="_>", ncols=75, bar_format="{l_bar}{bar}{n_fmt}/{total_fmt}")
            with Pool(processes=args.threads) as pool:
                for thread_idx in range(args.threads):
                    region_ranges = fetcher[thread_idx]
                    
                    pool.apply_async(worker_init, args=(bam_file, region_ranges, args, sample_idx, thread_idx), callback=update)
                pool.close()
                pool.join()
            pbar.close()
        
            out = open(f'{out_file}.vcf', 'a')
            print('Concatenating thread outputs!', file=sys.stderr)
            idx = out_file.rfind('/')
            hid_outfile = out_file[:idx+1] + '.' + out_file[idx+1:]
            for tidx in range(args.threads)[1:]:
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
                for tidx in range(args.threads)[1:]:
                    thread_log_out = f'{hid_outfile}_debug_{tidx}.log'
                    with open(thread_log_out, 'r') as fh:
                        for line in fh:
                            log_info = line.strip()
                            print(log_info, file=out_log)
                    os.remove(thread_log_out)
                out_log.close()
        else:
            worker_init(bam_file, fetcher[0], args, out_file, sample_idx, 0)
    time_now = ti.default_timer()
    sys.stderr.write('CPU time: {} seconds\n'.format(time_now - start_time))
