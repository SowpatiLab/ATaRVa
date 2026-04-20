import numpy as np
import pysam
import sys, os, logging

from pathlib import Path
from tqdm import tqdm
from sortedcontainers import SortedList
from collections import deque

from ATARVA.structures        import ReadLocusInfo, LocusInfo, ReadInfo, LocusVariation, ExtendedRead
from ATARVA.vcf_writer        import vcf_writer, write_homozygous_call, write_heterozygous_call, write_fail_call
from ATARVA.operation_utils   import record_homopolymers, clean_eqsign_readseq
from ATARVA.cstag_utils       import parse_cstag
from ATARVA.cigar_utils       import parse_cigar
from ATARVA.sub_operation_utils import mm_tag_extract, calculate_methylation
from ATARVA.locus_utils       import process_locus
from ATARVA.consensus         import consensus_seq_poa
from ATARVA.genotype_utils    import analyse_genotype


SKIP_MESSAGES = {
    0: 'Locus failed - insufficient reads support',
    1: 'Locus failed - insufficient reads for haplogroup',
    6: 'Locus skipped - the locus is not a repetitive sequence.'
}


def _hidden_path(out_file):
    """Return a hidden-file-prefixed path for thread output files."""
    path = Path(out_file)
    return str(path.parent / f'.{path.name}')


class PysamWarningCapture:
    """
    Context manager to capture pysam/htslib C-level stderr warnings
    and redirect them to a log file.
    """

    def __init__(self, logfile: str):
        self.logfile    = logfile
        self.log_fd     = None
        self.old_stderr = None

    def __enter__(self):
        self.log_fd     = open(self.logfile, 'a')
        self.old_stderr = os.dup(2)                    # save original stderr fd
        os.dup2(self.log_fd.fileno(), 2)               # redirect fd 2 → log file
        return self

    def __exit__(self, *args):
        sys.stderr.flush()
        os.dup2(self.old_stderr, 2)                    # restore original stderr
        os.close(self.old_stderr)
        self.log_fd.close()


class Cooper:
    """
    Independently handles genotyping a set of regions from a single BAM file.
    """

    def __init__(self, bam_file, region_ranges, args, out_file, sample_idx, thread_idx):
        """
        initialise Cooper and run genotyping.

        :param bam_file:      path to BAM/CRAM/SAM file
        :param region_ranges: list of (chrom, (start1,end1), (start2,end2)) tuples
        :param args:          parsed command-line arguments
        :param out_file:      output file base path
        :param sample_idx:    index of sample in BAM list
        :param thread_idx:    thread index (-1 = single-thread mode)
        """

        self.terminal_stderr = os.fdopen(os.dup(2), 'w')

        # --- output paths ---
        is_primary = thread_idx in (-1, 0)
        if is_primary:
            self.outfile = f'{out_file}.vcf'
            self.logfile = f'{out_file}_debug.log'
        else:
            hidden = _hidden_path(out_file)
            self.outfile = f'{hidden}_thread_{thread_idx}.vcf'
            self.logfile = f'{hidden}_debug_{thread_idx}.log'

        with PysamWarningCapture(self.logfile):
            self.tbx = pysam.Tabixfile(args.regions)
            self.bam = pysam.AlignmentFile(bam_file, args.aln_format)
            self.ref = pysam.FastaFile(args.fasta)

        self.args       = args
        self.sample_idx = sample_idx
        self.thread_idx = thread_idx
        self.karyotype  = args.karyotype[sample_idx]
        self.chrom      = None
        self.haploid    = False
        self.somatic    = False
        self.logger     = None

        self.outhandle = open(self.outfile, 'w')

        if is_primary:
            vcf_writer(self.outhandle, self.bam,
                       Path(bam_file).stem)

        # --- logging ---
        if args.debug_mode:
            Path(self.logfile).touch()
            Path(self.logfile).write_text('') # clear log file if it already exists
            logging.basicConfig(
                filename = self.logfile,
                level    = logging.DEBUG,
                format   = '%(levelname)s - %(message)s'
            )
            self.logger = logging.getLogger('ATaRVa_logger')

        self.disable_progress = False

        # --- count loci per region ---
        self.cooper_nloci  = 0
        self.range_nloci   = []

        for region_range in region_ranges:
            chrom, first_coords, last_coords = region_range
            count = 0
            for row in self.tbx.fetch(chrom, first_coords[0], last_coords[1]):
                fields     = row.split('\t')
                row_start  = int(fields[1])
                if count == 0 and row_start != first_coords[0]:
                    continue
                count += 1
                if int(fields[1]) == last_coords[0]:
                    break
            self.range_nloci.append(count)
            self.cooper_nloci += count

        # --- progress bar ---
        self.progress_bar = tqdm(
            total      = self.cooper_nloci,
            disable    = self.disable_progress,
            desc       = 'Processing ',
            file       = self.terminal_stderr, 
            ascii      = '_>',
            ncols      = 75,
            bar_format = '{l_bar}{bar}{n_fmt}/{total_fmt}'
        )

        # --- run genotyping ---
        for cidx, region_range in enumerate(region_ranges):
            chrom = region_range[0]
            bam_stem = Path(bam_file).stem
            print(f'Processing region {region_range} — sample {bam_stem} — thread {thread_idx}\n')
            self._reinitialise()
            self.cooper_readmode(region_range, cidx)

        self.bam.close()
        self.ref.close()
        self.tbx.close()
        self.outhandle.close()

    # --- state management ---

    def _reinitialise(self):
        """Reset per-region tracking state."""
        self.cooper_snp_data         = {}
        self.cooper_read_data        = {}
        self.cooper_loci_data        = {}
        self.cooper_loci_info        = {}
        self.cooper_loci_ends        = deque()
        self.cooper_loci_keys        = deque()
        self.cooper_read_ends        = deque()
        self.cooper_read_indices     = deque()
        self.prev_reads              = set()
        self.cooper_sorted_snps      = SortedList()
        self.cooper_insert_positions = set()

    # --- read processing ---

    def cooper_readmode(self, region_range: tuple, cidx: int):
        """
        Genotype a range of loci by streaming through reads.

        :param region_range: (chrom, (start1,end1), (start2,end2))
        :param cidx:         index of this region in region_ranges
        """
        chrom, first_coords, last_coords = region_range
        self.chrom   = chrom
        self.haploid = chrom in {'chrX', 'chrY', 'X', 'Y'} and self.karyotype

        region_start = first_coords[0]
        region_end   = last_coords[1]

        genotyped_count = 0
        read_index      = 0

        if not self.disable_progress:
            tqdm.write(f'> {chrom} {region_start} {region_end} '
                       f'Total loci = {self.range_nloci[cidx]}')
        with PysamWarningCapture(self.logfile):
            for raw_read in self.bam.fetch(chrom, region_start, region_end):

                if raw_read.mapping_quality < self.args.map_qual:
                    continue

                read = ExtendedRead.from_read(raw_read)

                # --- flush completed loci ---
                while self.cooper_loci_ends and read.ref_start > self.cooper_loci_ends[0]:
                    genotyped_count  += self.locus_processor()
                    self.progress_bar.update(1)

                # --- evict expired reads ---
                while self.cooper_read_ends:
                    if read.ref_start <= self.cooper_read_ends[0]:
                        break
                    if (self.cooper_loci_ends and
                            self.cooper_read_ends[0] > self.cooper_loci_ends[0]):
                        break

                    read_end = self.cooper_read_ends.popleft()
                    rindex   = self.cooper_read_indices.popleft()

                    if rindex in self.cooper_read_data:
                        for pos in self.cooper_read_data[rindex].snps:
                            if pos in self.cooper_snp_data:
                                self.cooper_snp_data[pos].cov -= 1
                                if self.cooper_snp_data[pos].cov == 0:
                                    del self.cooper_snp_data[pos]
                                    self.cooper_sorted_snps.remove(pos)
                        del self.cooper_read_data[rindex]
                        self.prev_reads.discard(rindex)

                    # evict expired insert positions
                    self.cooper_insert_positions = {
                        p for p in self.cooper_insert_positions if p > read_end
                    }

                # --- region end reached ---
                if read.ref_start > region_end:
                    while self.cooper_loci_ends:
                        genotyped_count  += self.locus_processor()
                        self.progress_bar.update(1)
                    break
                
                soft_start = read.ref_start
                soft_end   = read.ref_end
                if read.cigartuples[0][0] == 4:
                    soft_start = max(0, read.ref_start - read.cigartuples[0][1])
                if read.cigartuples[-1][0] == 4:
                    soft_end = read.ref_end + read.cigartuples[-1][1]
                if abs(soft_start - read.ref_start) > self.args.flank or abs(soft_end - read.ref_end) > self.args.flank:
                    normal_loci, softclip_loci = [[],[]]
                    for row in self.tbx.fetch(chrom, read.ref_start, read.ref_end):
                        f = row.split('\t')
                        if read.ref_start <= int(f[1]) and int(f[2]) <= read.ref_end:
                            normal_loci.append((int(f[1]), int(f[2])))
                    for row in self.tbx.fetch(chrom, soft_start, soft_end):
                        f = row.split('\t')
                        if soft_start <= int(f[1]) and int(f[2]) <= soft_end:
                            softclip_loci.append((int(f[1]), int(f[2])))

                    if len(normal_loci) < len(softclip_loci):
                        print('More loci in softclip range than normal range — use softclip loci')

                # --- assign loci to read ---
                for row in self.tbx.fetch(chrom, read.ref_start, read.ref_end):
                    fields      = row.split('\t')
                    locus_start = int(fields[1])
                    locus_end   = int(fields[2])
                    locus_len   = locus_end - locus_start

                    locus_name  = fields[5] if len(fields) > 5 else None

                    read.loci.append((locus_start, locus_end))

                    # region boundary checks
                    if locus_start < first_coords[0]:
                        continue
                    if locus_start >= last_coords[0]:
                        break
                    if not (first_coords[0] <= locus_start and locus_end <= last_coords[1]):
                        continue

                    # read must fully span the locus
                    if not (read.ref_start <= locus_start and locus_end <= read.ref_end):
                        continue

                    left_flank  = min(self.args.flank, locus_start - read.ref_start)
                    right_flank = min(self.args.flank, read.ref_end  - locus_end)

                    read.left_flanks.append(left_flank)
                    read.right_flanks.append(right_flank)
                    read.loci_coords.append((locus_start - left_flank,
                                            locus_end   + right_flank))

                    locus_key = f'{chrom}:{locus_start}-{locus_end}'
                    read.loci_keys.append(locus_key)
                    read.loci_data[locus_key] = ReadLocusInfo(halen=0, alen=0, rlen=locus_len, seq=[])

                    if locus_key not in self.cooper_loci_data:
                        self.cooper_loci_data[locus_key] = LocusVariation()
                        self.cooper_loci_info[locus_key] = LocusInfo(chrom=chrom, start=locus_start,
                                                                    end=locus_end, motif=fields[3],
                                                                    name=locus_name)
                        self.cooper_loci_ends.append(locus_end)
                        self.cooper_loci_keys.append(locus_key)

                if not read.loci_coords:
                    continue

                # --- EQX normalisation ---
                read_sequence = read.sequence
                for op, length in read.cigartuples:
                    if op in (0, 7):
                        break
                    if op in (1, 4, 8):
                        pass

                if '=' in read_sequence:
                    read_sequence = clean_eqsign_readseq(read.chrom, read.ref_start, read.cigartuples,
                                                        read_sequence, self.ref)

                # --- register read ---
                read_index   += 1
                read.index    = read_index

                self.cooper_read_ends.append(read.ref_end)
                self.cooper_read_indices.append(read.index)
                self.cooper_read_data[read.index] = ReadInfo(
                    start       = read.ref_start,
                    end         = read.ref_end,
                    snps        = set(),
                    dels        = [],
                    methylation = [],
                    mean_qual   = read.mean_qual,
                    left_flank  = read.left_flanks,
                    right_flank = read.right_flanks
                )

                # --- haplotag extraction ---
                read.haplotag = [False, None]
                if self.args.haplotag and read.has_tag(self.args.haplotag):
                    read.haplotag = [True, read.get_tag(self.args.haplotag)]

                # --- methylation extraction ---
                mod_bases = (list(read.modified_bases.items()) if read.modified_bases else [])

                if read.has_tag('cs'):
                    parse_cstag(self, read)
                else:
                    parse_cigar(self, read)

                for mods in mod_bases:
                    if mods[0][0] == 'C' and mods[0][2] == 'm':
                        mm_tag_extract(read, mods[1])
                        self.cooper_read_data[read_index].methylation = read.methylation_calls
                        break

                # --- populate loci data ---
                for locus_key, locus_read_info in read.loci_data.items():
                    if not locus_read_info.seq: continue
                    ldata = self.cooper_loci_data[locus_key]
                    ldata.reads.append(read.index)
                    ldata.read_alens[read.index]    = [locus_read_info.halen, locus_read_info.alen]
                    ldata.read_aseqs[read.index]    = locus_read_info.seq
                    ldata.read_haplotags.append(read.haplotag[1])

        # --- flush remaining loci ---
        while self.cooper_loci_ends:
            genotyped_count  += self.locus_processor()
            self.progress_bar.update(1)


    def locus_processor(self):
        """
        Genotype the next queued locus using collected read data.

        :return: 1 if locus was successfully genotyped, 0 otherwise
        """
        self.cooper_loci_ends.popleft()
        locus_key   = self.cooper_loci_keys.popleft()
        locus_data  = self.cooper_loci_data[locus_key]
        locus       = self.cooper_loci_info[locus_key]

        if locus_key not in self.cooper_loci_data:
            print(f'No reads for locus {locus_key} — skipping.')
            return 0

        # --- fetch neighbouring loci ---
        locus_data.neighbors = [
            (int(r.split('\t')[1]), int(r.split('\t')[2]))
            for r in self.tbx.fetch(self.chrom, locus.start - self.args.flank, locus.end + self.args.flank)
        ]

        # --- ensure sorted SNPs are populated ---
        if not self.cooper_sorted_snps:
            for pos in self.cooper_snp_data:
                self.cooper_sorted_snps.add(pos)

        # --- process locus ---
        process_locus(self, locus_key)

        locus_data    = self.cooper_loci_data[locus_key]
        read_seqs     = locus_data.read_aseqs
        # --- category 1 — homozygous ---
        if locus_data.hap_category == 1:
            read_seqs = [read_seqs[rid][0] for rid in locus_data.reads if read_seqs[rid][0] != '']

            if len(read_seqs) == 0:
                # homozygous deletion genotype
                ALT = '<DEL>'
                locus_data.gt_alens  = (0, 0)
                locus_data.gt_arange = '0-0,0-0'
                locus_data.gt_aseqs  = (ALT, ALT)

            elif read_seqs[0] == self.ref.fetch(locus.chrom, locus.start, locus.end):
                # homozygous reference genotype
                ALT = '.'
                locus_data.gt_alens  = (locus.length, locus.length)
                locus_data.gt_arange = f'{locus.length}-{locus.length},{locus.length}-{locus.length}'
                locus_data.gt_aseqs  = (ALT, ALT)
                ref_allele = self.ref.fetch(locus.chrom, locus.start, locus.end)
                meth_data = calculate_methylation(locus_data.reads, locus_data.read_methylation, ref_allele)
                locus_data.hap_meth_data = (meth_data, meth_data)

            else:
                ALT = read_seqs[0]
                locus_data.gt_alens      = (len(ALT), len(ALT))
                locus_data.gt_arange     = f'{len(ALT)}-{len(ALT)}'
                locus_data.gt_aseqs      = (ALT, ALT)
                meth_data = calculate_methylation(locus_data.reads, locus_data.read_methylation, ALT)
                locus_data.hap_meth_data = (meth_data, meth_data)

            write_homozygous_call(self, locus_key)
            locus_data.is_genotyped = 1

        # --- category 2 — ambiguous ---
        elif locus_data.hap_category == 2:
            analyse_genotype(self, locus_key)
            if not locus_data.is_genotyped:
                msg = SKIP_MESSAGES.get(
                    locus_data.skip_code,
                    "Locus skipped — insufficient significant SNPs.")
                # tqdm.write(msg)

        # --- category 3 — phased / heterozygous ---
        elif locus_data.hap_category == 3:
            genotypes    = []
            allele_count = {}
            ALT_seqs     = []
            phased_reads = []
            alen_lists   = []

            for hap_reads in locus_data.hap_read_sets:
                phased_reads.append(len(hap_reads))
                seqs      = [read_seqs[rid][0] for rid in hap_reads
                             if read_seqs[rid][0]]
                alen_list = [len(read_seqs[rid][0]) for rid in hap_reads]
                alen_lists.append(alen_list)

                if seqs:
                    ALT           = consensus_seq_poa(seqs)
                    allele_length = len(ALT)
                else:
                    ALT           = '<DEL>'
                    allele_length = 0

                ALT_seqs.append(ALT)
                genotypes.append(allele_length)
                key = allele_length if allele_length not in allele_count \
                      else str(allele_length)
                allele_count[key] = len(hap_reads)

            (l1, u1) = np.percentile(alen_lists[0], [2.5, 97.5])
            (l2, u2) = np.percentile(alen_lists[1], [2.5, 97.5])
            allele_range = f'{l1}-{u1},{l2}-{u2}'

            write_heterozygous_call(self, locus_key)
            locus_data.is_genotyped = 1

        # --- no category — write fail if needed ---
        else:
            if locus_data.skip_code == 0:
                write_fail_call(self, locus_key)

        # --- cleanup ---
        del self.cooper_loci_data[locus_key]
        del self.cooper_loci_info[locus_key]

        return locus_data.is_genotyped
