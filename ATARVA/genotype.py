import sys, os
import pysam
import timeit as ti
from multiprocessing import Pool
from tqdm import tqdm

from ATARVA.readers import fasta_check, bam_check, tabix_check
from ATARVA.version import __version__
from ATARVA.baseline import Cooper
from ATARVA.readers import *
from ATARVA.vcf_writer import set_info_mp_cutoff
from ATARVA.sub_operation_utils import set_methviz_tag


def genotype_parser(subparsers):
    """Argument parser for the genotype sub-command of ATaRVa and run the genotyping workflow"""

    parser = subparsers.add_parser(
        "genotype",
        help        = "tandem repeat genotyper specially designed for long read data",
        description = "Tandem Repeat Genotyper"
    )
    parser._action_groups.pop()

    # ── Required arguments ────────────────────────────────────────────────────
    req = parser.add_argument_group('Required arguments')
    req.add_argument('-f', '--fasta',   required=True, metavar='<FILE>', help='input reference fasta file')
    req.add_argument('-b', '--bam',     required=True, metavar='<FILE>', nargs='+', help='sample alignment files [SAM | BAM | CRAM]')
    req.add_argument('-r', '--regions', required=True, metavar='<FILE>', help='bgzipped + tabix-indexed regions file. '
                          'To prepare: sort with bedtools → compress with bgzip → index with tabix')

    # ── Optional arguments ────────────────────────────────────────────────────
    opt = parser.add_argument_group('Optional arguments')

    # --- Input / Output ---
    opt.add_argument('-o', '--vcf',      metavar='<FILE>', default='',    help='output VCF file [default: stdout]')
    opt.add_argument('--format',         metavar='<STR>',  default='bam', help='alignment file format [cram | bam | sam] [default: bam]')
    opt.add_argument('--contigs',        metavar='<STR>',  nargs='+', help='contigs to genotype e.g. chr1 chr12 chr22. [default: all]')
    opt.add_argument('--karyotype',      metavar='<STR>',  nargs='+', help='sample karyotypes e.g. XY XX')

    # --- Read filtering ---
    opt.add_argument('-q', '--map-qual', metavar='<INT>',   type=int,   default=5,    help='minimum mapping quality [default: 5]')
    opt.add_argument('--min-reads',      metavar='<INT>',   type=int,   default=10,   help='minimum read coverage at a locus [default: 10]')
    opt.add_argument('--max-reads',      metavar='<INT>',   type=int,   default=None, help='maximum reads used per locus [default: 100]')
    opt.add_argument('--flank',          metavar='<INT>',   type=int,   default=None, help='flank length (bp) to search for insertions [default: 10]')

    # --- SNP phasing ---
    opt.add_argument('--snp-dist',       metavar='<INT>',   type=int,   default=3000, help='max SNP distance from repeat for phasing [default: 3000]')
    opt.add_argument('--snp-count',      metavar='<INT>',   type=int,   default=3,    help='number of SNPs for phasing [default: 3]')
    opt.add_argument('--snp-qual',       metavar='<INT>',   type=int,   default=20,   help='min base quality at SNP position [default: 20]')
    opt.add_argument('--snp-read',       metavar='<FLOAT>', type=float, default=0.2,  help='min SNP read fraction for phasing [default: 0.2]')
    opt.add_argument('--phasing-read',   metavar='<FLOAT>', type=float, default=0.4,  help='min fraction of reads from phased clusters [default: 0.4]')
    opt.add_argument('--haplotag',       metavar='<STR>',               default=None, help='haplotag for phasing e.g. HP [default: None]')

    # --- Methylation ---
    opt.add_argument('--meth-prob',      metavar='<FLOAT>', type=float, default=0.5, help='min methylation probability cutoff [default: 0.5]')
    opt.add_argument('--methviz',        action='store_true', help='write methylation-encoded sequence to VCF [default: False]')

    # --- Modes ---
    opt.add_argument('--amplicon',       action='store_true', help='targeted sequencing mode [max-reads=1000, flank=20] [default: False]')
    opt.add_argument('--read-wise',      action='store_true', help='read-wise genotyping for dense BED regions [default: False]')
    opt.add_argument('--loci-wise',      action='store_true', help='loci-wise genotyping for sparse BED regions [default: False]')
    opt.add_argument('--decompose',      action='store_true', help='write motif-decomposed sequence to VCF [default: False]')

    # --- Misc ---
    opt.add_argument('-t', '--threads',      metavar='<INT>',     type=int,   default=1, help='number of threads [default: 1]')
    opt.add_argument('-log', '--debug_mode', action='store_true', help='write debug messages to log file [default: False]')


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
    """
    Main function to run the genotyping workflow of ATaRVa.

    :param args: command line arguments
    """

    start_time = ti.default_timer()

    for arg in vars(args):
        if arg in ['func', 'help', 'command']: continue
        print (arg, getattr(args, arg))
    print('\n')

    # --- File checks ---
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
    external_name = out_file

    # --- parameters to be set in VCF ---
    set_info_mp_cutoff(args.meth_prob)
    set_methviz_tag(args.methviz)

    # --- process tbx ---
    temp_tbx  = pysam.Tabixfile(args.regions)
    total_loci = 0
    if not args.contigs:
        args.contigs = sorted(temp_tbx.contigs)
        for row in temp_tbx.fetch(): total_loci += 1
    else:
        args.contigs = sorted(args.contigs)
        total_loci = sum(1 for contig in args.contigs for row in temp_tbx.fetch(contig))

    split_point = total_loci // args.threads
    if split_point == 0:
        split_point = 1
        args.threads = 1

    fetcher = splitfile_threads(temp_tbx, args.contigs, total_loci, args.threads)
    temp_tbx.close()

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
