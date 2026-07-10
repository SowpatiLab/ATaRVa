import numpy as np
import pysam
import sys, os, logging

from pathlib import Path
from tqdm import tqdm
from sortedcontainers import SortedList
from collections import deque

from ATARVA.structures        import ReadLocusInfo, LocusInfo, ReadInfo, LocusVariation, ExtendedRead
from ATARVA.vcf_writer        import vcf_writer, write_homozygous_call, write_heterozygous_call, write_fail_call
from ATARVA.operation_utils   import clean_eqsign_readseq
from ATARVA.cstag_utils       import parse_cstag
from ATARVA.cigar_utils       import parse_cigar
from ATARVA.sub_operation_utils import mm_tag_extract, calculate_methylation, clamp_zero
from ATARVA.locus_utils       import process_locus
from ATARVA.consensus         import consensus_seq_poa
from ATARVA.genotype_utils    import analyse_genotype
from ATARVA.process_softclips import process_flank_stretches, check_flank
from ATARVA.flank_utils       import check_flank_order, process_flanks


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
            if args.instability:
                self.insfile = f'{out_file}_instability.jsonl'
        else:
            hidden = _hidden_path(out_file)
            self.outfile = f'{hidden}_thread_{thread_idx}.vcf'
            self.logfile = f'{hidden}_debug_{thread_idx}.log'
            if args.instability:
                self.insfile = f'{hidden}_instability_{thread_idx}.jsonl'

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
        if args.instability: self.ins_handle = open(self.insfile, 'w')

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
        if self.args.threads > 1: self.disable_progress = True

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
            total        = self.cooper_nloci,
            disable      = self.disable_progress,
            position     = 1 + 2*thread_idx + 1,
            desc         = 'Processing ' if args.threads == 1 else f'Processing thread {thread_idx + 1}',
            file         = self.terminal_stderr, 
            ascii        = '_>',
            ncols        = 80,
            bar_format   = '{l_bar}{bar}{n_fmt}/{total_fmt}',
            leave        = True,
            mininterval  = 1      # refresh rate that helps reduce overhead and doesn't cause stray progress bars in multi-thread mode
        )

        # --- run genotyping ---
        for cidx, region_range in enumerate(region_ranges):
            chrom = region_range[0]
            bam_stem = Path(bam_file).stem
            self._reinitialise()
            self.cooper_readmode(region_range, cidx)

        self.bam.close()
        self.ref.close()
        self.tbx.close()
        self.outhandle.close()
        if self.args.instability: self.ins_handle.close()

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
        NONREP_FLANK    = 30        # Minimum non-repetitive flank considered for locus processing from softclip region

        softclip_mode   = True
        DROP_DISTANCE   = 0 if softclip_mode else 100000    # distance beyond which reads and loci are dropped from memory

        with PysamWarningCapture(self.logfile):
            for raw_read in self.bam.fetch(chrom, region_start, region_end):

                if raw_read.mapping_quality < self.args.map_qual:
                    continue

                read = ExtendedRead.from_read(raw_read)

                start_softclip = end_softclip = 0
                if read.cigartuples[0][0]  == 4 and softclip_mode: start_softclip = read.cigartuples[0][1]
                if read.cigartuples[-1][0] == 4 and softclip_mode: end_softclip = read.cigartuples[-1][1]
                fetch_start = max(0, read.ref_start - start_softclip)
                fetch_end   = min(read.ref_end + end_softclip, last_coords[1])

                # --- flush completed loci ---
                while self.cooper_loci_ends and fetch_start - DROP_DISTANCE > self.cooper_loci_ends[0]:
                    genotyped_count  += self.locus_processor()
                    self.progress_bar.update(1)

                loci_rindices = set()
                for lk in self.cooper_loci_data: loci_rindices |= set(self.cooper_loci_data[lk].reads)
                # --- evict expired reads ---
                while self.cooper_read_ends and fetch_start - DROP_DISTANCE > self.cooper_read_ends[0]:
                    if (self.cooper_loci_ends and self.cooper_read_ends[0] > self.cooper_loci_ends[0]):
                        break
                    if self.cooper_read_indices[0] in loci_rindices: break

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
                if fetch_start > region_end:
                    while self.cooper_loci_ends:
                        genotyped_count  += self.locus_processor()
                        self.progress_bar.update(1)
                    break

                # process the SA tag to store relevant supplementary alignments
                # if read.has_tag('SA'): read.process_satag()

                # --- assign loci to read ---
                softclip_loci  = {'keys': [], 'loci': [], 'coords': [], 'flags': []}
                for row in self.tbx.fetch(chrom, fetch_start, fetch_end):
                    fields      = row.split('\t')
                    locus_start = int(fields[1])
                    locus_end   = int(fields[2])
                    locus_len   = locus_end - locus_start

                    locus_name  = fields[5] if len(fields) > 5 else None

                    read.loci.append((locus_start, locus_end))

                    # region boundary checks
                    if locus_start < first_coords[0]:
                        continue
                    if locus_start >= last_coords[1]:
                        break

                    if not (first_coords[0] <= locus_start and locus_end <= last_coords[1]):
                        continue
                    
                    # check if the locus is outside the read's reference boundaries check if it's in the softclipped region
                    softclip_result = None
                    if softclip_mode and (locus_start < read.ref_start or locus_end > read.ref_end):
                        softclip_result = check_flank(self, read, locus_start, locus_end, start_softclip, end_softclip)
                    # result structure {'upstream':   (ref_flank_start, ref_flank_end, query_flank_start, query_flank_end),
                    #                   'downstream': (ref_flank_start, ref_flank_end, query_flank_start, query_flank_end)}
                    if softclip_result is not None:
                        flank_order = check_flank_order(softclip_result)
                        if flank_order:
                            softclip_loci['keys'].append(f'{chrom}:{locus_start}-{locus_end}')
                            softclip_loci['loci'].append((chrom, locus_start, locus_end))
                            softclip_loci['coords'].append(softclip_result)
                            softclip_loci['flags'].append(None)
                        else:
                            softclip_loci['keys'].append(f'{chrom}:{locus_start}-{locus_end}')
                            softclip_loci['loci'].append((chrom, locus_start, locus_end))
                            softclip_loci['coords'].append(softclip_result)
                            softclip_loci['flags'].append('FLANK_ORDER_INVALID')
                    if not (read.ref_start <= locus_start and locus_end <= read.ref_end) and not softclip_result:
                        continue

                    left_flank  = min(self.args.flank, clamp_zero(locus_start - read.ref_start))
                    right_flank = min(self.args.flank, clamp_zero(read.ref_end  - locus_end))

                    read.left_flanks.append(left_flank)
                    read.right_flanks.append(right_flank)
                    read.loci_coords.append((locus_start - left_flank, locus_end + right_flank))

                    locus_key = f'{chrom}:{locus_start}-{locus_end}'
                    read.loci_keys.append(locus_key)
                    read.loci_data[locus_key] = ReadLocusInfo(halen=0, alen=0, rlen=locus_len, seq=[])

                    if locus_key not in self.cooper_loci_data:
                        self.cooper_loci_data[locus_key] = LocusVariation()
                        self.cooper_loci_info[locus_key] = LocusInfo(chrom=chrom, start=locus_start, end=locus_end,
                                                                     motif=fields[3], name=locus_name)
                        self.cooper_loci_ends.append(locus_end)
                        self.cooper_loci_keys.append(locus_key)

                if not read.loci_coords:
                    continue

                if softclip_loci['coords']:
                    merged_coords = []
                    if sum([flag is None for flag in softclip_loci['flags']]) >= 1:
                        merged_coords = process_flanks(softclip_loci, read.ref_start, read.ref_end)
                    if len(merged_coords) > 0:
                        process_flank_stretches(self, read, merged_coords)

                read_required = True
                for key in read.loci_keys:
                    if self.cooper_loci_data[key].depth >= self.args.max_reads and self.cooper_loci_data[key].min_read_qual >= read.mean_qual:
                        read_required = False
                    else:
                        read_required = True
                        break
                # if all loci have enough supporting reads with highest quality; this read is not processed
                if not read_required: continue

                if '=' in read.sequence:
                    read.sequence = clean_eqsign_readseq(read.chrom, read.ref_start, read.cigartuples,
                                                         read.sequence, self.ref)

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
                mod_bases = ()
                if not self.args.rna:
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
                    if ldata.depth >= self.args.max_reads:
                        if read.mean_qual > ldata.min_read_qual:
                            ldata.reads.append(read.index)
                            if self.args.instability: ldata.read_names.append(read.query_name)
                            ldata.depth += 1
                            ldata.read_alens[read.index]     = [locus_read_info.halen, locus_read_info.alen]
                            ldata.read_aseqs[read.index]     = locus_read_info.seq
                            ldata.read_haplotags[read.index] = read.haplotag[1]
                            # remove lowest quality read
                            if ldata.min_qual_read in ldata.reads:
                                ldata.reads.remove(ldata.min_qual_read)
                                del ldata.read_alens[ldata.min_qual_read]
                                del ldata.read_aseqs[ldata.min_qual_read]
                                del ldata.read_haplotags[ldata.min_qual_read]
                                ldata.depth -= 1
                            # update minimum quality read
                            ldata.min_read_qual = float('inf')
                            for r in ldata.reads:
                                if ldata.min_read_qual > self.cooper_read_data[r].mean_qual:
                                    ldata.min_read_qual = self.cooper_read_data[r].mean_qual
                                    ldata.min_qual_read = r

                    else:
                        ldata.reads.append(read.index)
                        if self.args.instability: ldata.read_names.append(read.query_name)
                        ldata.depth += 1
                        ldata.read_alens[read.index]     = [locus_read_info.halen, locus_read_info.alen]
                        ldata.read_aseqs[read.index]     = locus_read_info.seq
                        ldata.read_haplotags[read.index] = read.haplotag[1]
                        if ldata.min_read_qual > read.mean_qual:
                            ldata.min_read_qual = read.mean_qual
                            ldata.min_qual_read = read.index

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
            write_fail_call(self, locus_key)
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
            read_seqs = [read_seqs[rid][0] for rid in locus_data.reads]
            seq_counter = {}
            for seq in read_seqs:
                try: seq_counter[seq] += 1
                except KeyError: seq_counter[seq] = 1
            max_allele = max(seq_counter, key=seq_counter.get) if seq_counter else None

            if len(read_seqs) == 0 or max_allele == '':
                # homozygous deletion genotype
                ALT = '<DEL>'
                locus_data.gt_alens  = (0, 0)
                locus_data.gt_arange = '0-0,0-0'
                locus_data.gt_aseqs  = (ALT, ALT)

            elif max_allele == self.ref.fetch(locus.chrom, locus.start, locus.end):
                # homozygous reference genotype
                ALT = '.'
                locus_data.gt_alens  = (locus.length, locus.length)
                locus_data.gt_arange = f'{locus.length}-{locus.length},{locus.length}-{locus.length}'
                locus_data.gt_aseqs  = (max_allele, max_allele)
                ref_allele = self.ref.fetch(locus.chrom, locus.start, locus.end)
                meth_data = calculate_methylation(locus_data.reads, locus_data.read_methylation, ref_allele)
                locus_data.hap_meth_data = (meth_data, meth_data)

            else:
                ALT = max_allele
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

        if self.args.instability and locus_data.is_genotyped:
            instability_info = []
            hap_als = [locus_data.gt_alens[0], locus_data.gt_alens[1]]
            for i, rid in enumerate(locus_data.reads):
                read_name = locus_data.read_names[i]
                aseq = locus_data.read_aseqs[rid][0]
                alen = len(aseq)
                hap = 0
                if not self.haploid: hap = 0 if rid in locus_data.hap_read_sets[0] else 1

                methyl_Cs = None; methyl_probab = None
                if i in locus_data.read_methylation and  locus_data.read_methylation[i] is not None:
                    methyl_Cs = len(locus_data.read_methylation[i][1])
                    methyl_probab = locus_data.read_methylation[i][0]

                instability_info.append((read_name, hap, alen, aseq, methyl_Cs, methyl_probab))
            instability_info.sort(key=lambda x: x[1]) # sort by read name for consistent output
            for read_name, hap, alen, aseq, methyl_Cs, methyl_probab in instability_info:
                self.ins_handle.write(f"{locus.chrom}\t{locus.start}\t{locus.end}\t{locus.motif}\t{read_name}\t{hap}\t{hap_als[hap]}\t{alen}\t{aseq}\t{methyl_Cs}\t{methyl_probab}\n")

        # --- cleanup ---
        del self.cooper_loci_data[locus_key]
        del self.cooper_loci_info[locus_key]

        return locus_data.is_genotyped
