import sys
import os
import timeit as ti
from pathlib import Path
from multiprocessing import Pool
from tqdm import tqdm
import pysam

from ATARVA.readers       import fasta_check, bam_check, tabix_check
from ATARVA.version       import __version__
from ATARVA.baseline      import Cooper
from ATARVA.readers       import *
from ATARVA.vcf_writer    import set_info_mp_cutoff
from ATARVA.sub_operation_utils import set_methviz_tag

# --- Constants ---
FORMAT_MAP = {
    'cram': 'rc',
    'sam' : 'r',
    'bam' : 'rb'
}

DEFAULT_MAX_READS = 100
DEFAULT_FLANK     = 10

def genotype_parser(subparsers):
    """Argument parser for the genotype sub-command of ATaRVa"""

    parser = subparsers.add_parser(
        "genotype",
        help        = "tandem repeat genotyper for long read data",
        description = "Tandem Repeat Genotyper"
    )
    parser._action_groups.pop()

    req = parser.add_argument_group('Required arguments')
    req.add_argument('-f', '--fasta',   required=True, metavar='<FILE>', help='input reference fasta file')
    req.add_argument('-b', '--bam',     required=True, metavar='<FILE>', nargs='+',
                     help='sample alignment files [SAM | BAM | CRAM]')
    req.add_argument('-r', '--regions', required=True, metavar='<FILE>',
                     help='bgzipped + tabix-indexed regions file. '
                          'Prepare: sort with bedtools → bgzip → tabix index')

    opt = parser.add_argument_group('Optional arguments')

    # Input / Output
    opt.add_argument('-o', '--vcf',      metavar='<FILE>', default='', help='output VCF file [default: stdout]')
    opt.add_argument('--aln-format',         metavar='<STR>',  default='bam', help='alignment format [cram | bam | sam] [default: bam]')
    opt.add_argument('--rna',         action='store_true', help='if the input alignment data is RNA-seq [default: False]')
    opt.add_argument('--instability', action='store_true', help='if instability metrics should be calculated and written to VCF [default: False]')
    opt.add_argument('--contigs',        metavar='<STR>',  nargs='+', help='contigs to genotype e.g. chr1 chr12 [default: all]')
    opt.add_argument('--karyotype',      metavar='<STR>',  nargs='+', help='sample karyotypes e.g. XY XX')

    # Read filtering
    opt.add_argument('-q', '--map-qual', metavar='<INT>',   type=int,   default=5, help='minimum mapping quality [default: 5]')
    opt.add_argument('--min-reads',      metavar='<INT>',   type=int,   default=10, help='minimum read coverage at a locus [default: 10]')
    opt.add_argument('--max-reads',      metavar='<INT>',   type=int,   default=None, help='maximum reads per locus [default: 100]')
    opt.add_argument('--flank',          metavar='<INT>',   type=int,   default=None, help='flank length (bp) to search for insertions [default: 10]')

    # SNP phasing
    opt.add_argument('--snp-dist',       metavar='<INT>',   type=int,   default=3000, help='max SNP distance from repeat [default: 3000]')
    opt.add_argument('--snp-count',      metavar='<INT>',   type=int,   default=3, help='number of SNPs for phasing [default: 3]')
    opt.add_argument('--snp-qual',       metavar='<INT>',   type=int,   default=20, help='min base quality at SNP position [default: 20]')
    opt.add_argument('--snp-read',       metavar='<FLOAT>', type=float, default=0.2, help='min SNP read fraction [default: 0.2]')
    opt.add_argument('--phasing-read',   metavar='<FLOAT>', type=float, default=0.4, help='min phased read cluster fraction [default: 0.4]')
    opt.add_argument('--haplotag',       metavar='<STR>',               default=None, help='haplotag for phasing e.g. HP [default: None]')

    # Methylation
    opt.add_argument('--meth-prob',      metavar='<FLOAT>',  type=float, default=0.5, help='min methylation probability [default: 0.5]')
    opt.add_argument('--methviz',        action='store_true', help='write methylation-encoded sequence to VCF [default: False]')

    # Modes
    opt.add_argument('--read-wise',      action='store_true', help='read-wise genotyping for dense BED regions')
    opt.add_argument('--locus-wise',     action='store_true', help='locus-wise genotyping for sparse BED regions')
    opt.add_argument('--decompose',      action='store_true', help='write motif-decomposed sequence to VCF')

    # Misc
    opt.add_argument('-t',   '--threads',      metavar='<INT>', type=int, default=1, help='number of threads [default: 1]')
    opt.add_argument('-log', '--debug_mode', action='store_true', help='write debug messages to log file')
    opt.add_argument('-v',   '--version',      action='version', version=f'ATaRVa version {__version__}')

    if len(sys.argv) == 2 and sys.argv[1] == 'genotype':
        parser.print_help()
        sys.exit()

    parser.set_defaults(func=genotype_run)


def splitfile_threads(tbx, contigs, total_loci, threads):
    """
    split the input region file into per-thread chunks for parallel processing.

    :param tbx:        pysam Tabixfile object
    :param contigs:    list of contigs to process
    :param total_loci: total number of loci across all contigs
    :param threads:    number of parallel threads
    :return:           list of region chunk tuples, one per thread
    """

    split_size  = total_loci // threads
    if split_size == 0:
        split_size = total_loci
    fetcher     = []
    line_count  = 0
    cur_chunk   = []
    chrom       = None
    start_coord = None

    for contig in contigs:
        first_in_contig = True

        for row in tbx.fetch(contig):
            fields = row.split('\t')
            line_count += 1

            if first_in_contig:
                chrom       = fields[0]
                start_coord = (int(fields[1]), int(fields[2]))
                first_in_contig = False

            if len(fetcher) < threads - 1 and line_count % split_size == 0:
                end_coord = (int(fields[1]), int(fields[2]))
                cur_chunk.append([chrom, start_coord, end_coord])
                fetcher.append(tuple(cur_chunk))
                cur_chunk       = []
                line_count      = 0
                first_in_contig = True

        if not first_in_contig:
            end_coord = (int(fields[1]), int(fields[2]))
            cur_chunk.append([chrom, start_coord, end_coord])

    fetcher.append(tuple(cur_chunk)) # pad with empty chunks if needed
    return fetcher


def worker_init(bam_file, region_ranges, args, out_file, sample_idx, thread_idx):
    """
    Initialise and run a Cooper genotyping worker for a thread.

    :param bam_file:      path to the alignment file
    :param region_ranges: list of regions for this thread
    :param args:          parsed command-line arguments
    :param out_file:      output file path
    :param sample_idx:    sample index
    :param thread_idx:    thread index
    """
    Cooper(bam_file, region_ranges, args, out_file, sample_idx, thread_idx)


def resolve_output_path(vcf_arg):
    """
    Resolve and normalise the output file path from --vcf argument.
    :param vcf_arg: raw --vcf argument value
    :return: resolved base output path (no extension)
    """
    if not vcf_arg:
        return ''
    path = vcf_arg.rstrip('/')
    if path.endswith('.vcf'):
        path = path[:-4]
    return path


def sample_stem(bam_path: str) -> str:
    """Extract sample name stem from BAM file path."""
    return Path(bam_path).stem


def concat_thread_outputs(base_path: str, n_threads: int, debug: bool = False, instability: bool = False) -> None:
    """
    Concatenate per-thread VCF (and optionally log) outputs into a single file.

    :param base_path: base output path (no extension)
    :param n_threads: total number of threads
    :param debug:     whether to also concatenate debug logs
    :param instability: whether to also concatenate instability metrics
    """
    parent   = Path(base_path).parent
    stem     = Path(base_path).name
    hidden   = parent / f'.{stem}'              # hidden prefix for thread files

    # --- VCF concatenation ---
    print('\nConcatenating thread outputs...', file=sys.stderr)
    with open(f'{base_path}.vcf', 'a') as out_vcf:
        for tidx in range(1, n_threads):
            thread_vcf = f'{hidden}_thread_{tidx}.vcf'
            with open(thread_vcf) as fh:
                for line in fh:
                    print(*line.strip().split('\t'), file=out_vcf, sep='\t')
            os.remove(thread_vcf)
            os.remove(f'{hidden}_debug_{tidx}.log')
    print('Concatenation complete.', file=sys.stderr)

    # --- log concatenation ---
    if debug:
        with open(f'{base_path}_debug.log', 'a') as out_log:
            for tidx in range(1, n_threads):
                thread_log = f'{hidden}_debug_{tidx}.log'
                with open(thread_log) as fh:
                    for line in fh:
                        print(line.strip(), file=out_log)
                os.remove(thread_log)

    if instability:
        with open(f'{base_path}_instability.tsv', 'a') as out_ins:
            for tidx in range(1, n_threads):
                thread_ins = f'{hidden}_instability_{tidx}.tsv'
                with open(thread_ins) as fh:
                    for line in fh:
                        print(line.strip(), file=out_ins)
                os.remove(thread_ins)


def genotype_run(args) -> None:
    """
    Main entry point for the ATaRVa genotyping workflow.

    :param args: parsed command-line arguments
    """
    start_time = ti.default_timer()

    # --- print args ---
    for arg in vars(args):
        if arg not in ('func', 'help', 'command'):
            print(arg, getattr(args, arg), file=sys.stderr)

    # --- file validation ---
    fasta_check(args.fasta)
    tabix_check(args.regions)

    args.aln_format = FORMAT_MAP.get(args.aln_format, 'rb')

    for bam in args.bam:
        bam_check(bam, args.aln_format)

    # --- VCF output setup ---
    set_info_mp_cutoff(args.meth_prob)
    set_methviz_tag(args.methviz)

    base_output = resolve_output_path(args.vcf)

    # --- tabix region processing ---
    tbx = pysam.Tabixfile(args.regions)
    if not args.contigs:
        args.contigs = sorted(tbx.contigs)
        total_loci   = sum(1 for _ in tbx.fetch())
    else:
        args.contigs = sorted(args.contigs)
        total_loci   = sum(1 for c in args.contigs for _ in tbx.fetch(c))

    if total_loci < args.threads:
        args.threads = 1

    fetcher = splitfile_threads(tbx, args.contigs, total_loci, args.threads)
    tbx.close()

    # --- karyotype setup ---
    args.karyotype = (
        [k == 'XY' for k in args.karyotype]
        if args.karyotype
        else [False] * len(args.bam)
    )

    # --- mode defaults ---
    if args.max_reads is None: args.max_reads = DEFAULT_MAX_READS
    if args.flank     is None: args.flank     = DEFAULT_FLANK

    multi_bam_single_out = len(args.bam) > 1 and args.vcf

    # --- per-sample processing ---
    for sample_idx, bam_file in enumerate(args.bam):
        sample_name = sample_stem(bam_file)
        print(f'Processing sample: {sample_name}\n', file=sys.stderr)

        check_alnfile(bam_file, args)

        # resolve per-sample output path
        if not args.vcf:
            out_file = sample_name
        elif multi_bam_single_out or args.vcf.endswith('/'):
            out_file = f'{base_output}_{sample_name}'
        else:
            out_file = base_output

        # --- multithreaded ---
        if args.threads > 1:
            all_success = True
            errors      = []

            with tqdm(total=args.threads, desc='Processing',
                    ascii='_>', ncols=75,
                    bar_format='{l_bar}{bar}{n_fmt}/{total_fmt}') as pbar:

                def on_success(_): pbar.update()

                def on_error(e):
                    nonlocal all_success
                    all_success = False
                    errors.append(str(e))
                    pbar.update()
                    print(f'Thread failed: {e}', file=sys.stderr)

                with Pool(processes=args.threads) as pool:
                    results = []
                    for thread_idx in range(args.threads):
                        r = pool.apply_async(
                            worker_init,
                            args     = (bam_file, fetcher[thread_idx],
                                        args, out_file, sample_idx, thread_idx),
                            callback = on_success,
                            error_callback = on_error    # catch thread failures
                        )
                        results.append(r)
                    pool.close()
                    pool.join()

                # ── verify all threads completed successfully ─────────────────
                for i, r in enumerate(results):
                    try:
                        r.get()                          # raises exception if thread failed
                    except Exception as e:
                        all_success = False
                        errors.append(f'Thread {i}: {e}')

            # ── only concat if all threads succeeded ─────────────────────────
            if all_success:
                concat_thread_outputs(out_file, args.threads, debug=args.debug_mode, instability=args.instability)
            else:
                print(f'Skipping concat — {len(errors)} thread(s) failed:',
                    file=sys.stderr)
                for err in errors:
                    print(f'  ✗ {err}', file=sys.stderr)
                sys.exit(1)

        # --- single threaded ---
        else:
            worker_init(bam_file, fetcher[0], args, out_file, sample_idx, 0)

    elapsed = ti.default_timer() - start_time
    sys.stderr.write(f'\nCPU time: {elapsed:.2f} seconds\n')
