import numpy as np
import pysam
import logging

from pathlib import Path
from tqdm import tqdm
from sortedcontainers import SortedList
from collections import deque

from ATARVA.structures        import ReadLocusInfo, LocusInfo, ReadInfo, LocusVariation, ExtendedRead
from ATARVA.vcf_writer        import vcf_writer, write_homozygous_call, vcf_heterozygous_writer, write_fail_call
from ATARVA.operation_utils   import record_homopolymers, clean_eqsign_readseq
from ATARVA.cstag_utils       import parse_cstag
from ATARVA.cigar_utils       import parse_cigar
from ATARVA.sub_operation_utils import mm_tag_extract, calculate_methylation
from ATARVA.locus_utils       import process_locus
from ATARVA.consensus         import consensus_seq_poa
from ATARVA.genotype_utils    import analyse_genotype


SKIP_MESSAGES = {
    0: 'Locus skipped — insignificant SNPs at read-split level.',
    1: 'Locus skipped — low read contribution from significant SNPs.',
    2: 'Locus skipped — insufficient reads in phased clusters.',
    6: 'Locus skipped — wide allele distribution with single-read support.',
}


class Cooper:
    """
    Independently handles genotyping a set of regions from a single BAM file.
    """

    def __init__(self, bam_file: str, region_ranges: list,
                 args, out_file: str,
                 sample_idx: int, thread_idx: int):
        """
        initialise Cooper and run genotyping.

        :param bam_file:      path to BAM/CRAM/SAM file
        :param region_ranges: list of (chrom, (start1,end1), (start2,end2)) tuples
        :param args:          parsed command-line arguments
        :param out_file:      output file base path
        :param sample_idx:    index of sample in BAM list
        :param thread_idx:    thread index (-1 = single-thread mode)
        """

        self.tbx = pysam.Tabixfile(args.regions)
        self.bam = pysam.AlignmentFile(bam_file, args.format)
        self.ref = pysam.FastaFile(args.fasta)

        self.args       = args
        self.sample_idx = sample_idx
        self.thread_idx = thread_idx
        self.karyotype  = args.karyotype[sample_idx]
        self.chrom      = None
        self.haploid    = False
        self.somatic    = False
        self.logger     = None

        # --- output paths ---
        is_primary = thread_idx in (-1, 0)
        if is_primary:
            self.outfile = f'{out_file}.vcf'
            self.logfile = f'{out_file}_debug.log'
        else:
            hidden = _hidden_path(out_file)
            self.outfile = f'{hidden}_thread_{thread_idx}.vcf'
            self.logfile = f'{hidden}_debug_{thread_idx}.log'

        self.outhandle = open(self.outfile, 'w')

        if is_primary:
            vcf_writer(self.outhandle, self.bam,
                       Path(bam_file).stem)

        # --- logging ---
        if args.debug_mode:
            Path(self.logfile).touch()
            logging.basicConfig(
                filename = self.logfile,
                level    = logging.DEBUG,
                format   = '%(levelname)s - %(message)s'
            )
            self.logger = logging.getLogger('ATaRVa_logger')

        if args.amplicon:
            args.haplotag = None

        self.disable_progress = thread_idx != -1

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

        for raw_read in self.bam.fetch(chrom, region_start, region_end):

            if raw_read.mapping_quality < self.args.map_qual:
                continue

            read = ExtendedRead.from_read(raw_read)

            # --- flush completed loci ---
            while self.cooper_loci_ends and read.ref_start > self.cooper_loci_ends[0]:
                genotyped_count  += self.locus_processor()
                self.progress_bar.update(1)

            # --- evict stale reads ---
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
                read.loci_data[locus_key] = ReadLocusInfo(
                    halen=0, alen=0, rlen=locus_len, seq=[])

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
                qual        = read.mean_qual,
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
                    self.cooper_read_data[read_index].methylation = \
                        read.methylation_calls
                    break

            # --- populate loci data ---
            for locus_key, locus_read_info in read.loci_data.items():
                if self.args.amplicon and not locus_read_info.seq:
                    continue
                ldata = self.cooper_loci_data[locus_key]
                ldata.reads.append(read.index)
                ldata.read_alens[read.index]    = [locus_read_info.halen,
                                                   locus_read_info.alen]
                ldata.read_seqs[read.index]     = locus_read_info.seq
                ldata.read_haplotags.append(read.haplotag[1])

        # --- flush remaining loci ---
        while self.cooper_loci_ends:
            genotyped_count  += self.locus_processor()
            self.progress_bar.update(1)

        self._reinitialise()
        self.bam.close()
        self.ref.close()
        self.tbx.close()
        self.outhandle.close()


    def locus_processor(self):
        """
        Genotype the next queued locus using collected read data.

        :return: 1 if locus was successfully genotyped, 0 otherwise
        """
        self.cooper_loci_ends.popleft()
        locus_key = self.cooper_loci_keys.popleft()
        locus     = self.cooper_loci_info[locus_key]

        if locus_key not in self.cooper_loci_data:
            print(f'No reads for locus {locus_key} — skipping.')
            return 0

        # --- fetch neighbouring loci ---
        neighbours = [
            (int(r.split('\t')[1]), int(r.split('\t')[2]))
            for r in self.tbx.fetch(self.chrom, locus.start - self.args.flank, locus.end + self.args.flank)
        ]

        # --- ensure sorted SNPs are populated ---
        if not self.cooper_sorted_snps:
            for pos in self.cooper_snp_data:
                self.cooper_sorted_snps.add(pos)

        # --- process locus ---
        (category, homozygous_allele, read_indices,
         hallele_counter, skip_point,
         haplotypes, allele_lengths) = process_locus(self, locus_key, neighbours)

        read_seqs     = self.cooper_loci_data[locus_key].read_seqs
        genotyped     = 0

        print(category)
        # --- category 1 — homozygous ---
        if category == 1:
            seqs = [read_seqs[rid][0] for rid in read_indices if read_seqs[rid][0]]

            if homozygous_allele != locus.length and seqs:
                ALT              = consensus_seq_poa(seqs)
                homozygous_allele = len(ALT)
            else:
                ALT = '.' if homozygous_allele == locus.length else '<DEL>'

            lower, upper = (round(x) for x in np.percentile(np.array(allele_lengths), [2.5, 97.5]))

            ref_seq   = self.ref.fetch(locus.chrom, locus.start, locus.end)
            meth_src  = ALT if ALT != '.' else ref_seq
            meth_info = calculate_methylation(read_indices, self.cooper_loci_data[locus_key].read_methylation, meth_src)

            allele_range = f'{lower}-{upper}' if self.haploid else f'{lower}-{upper},{lower}-{upper}'

            write_homozygous_call(self, locus_key, homozygous_allele, hallele_counter, allele_range, meth_info,
                                  '.', ALT, len(read_indices), len(read_indices))
            genotyped = 1

        # --- category 2 — ambiguous ---
        elif category == 2:
            state, skip_point = analyse_genotype(self, locus_key, hallele_counter, read_indices)
            if state:
                genotyped = 1
            else:
                msg = SKIP_MESSAGES.get(
                    skip_point,
                    "Locus skipped — insufficient significant SNPs.")
                tqdm.write(msg)

        # --- category 3 — phased / heterozygous ---
        elif category == 3:
            genotypes    = []
            allele_count = {}
            ALT_seqs     = []
            phased_reads = []
            alen_lists   = []
            meth_info    = []

            for hap_reads in haplotypes:
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

            vcf_heterozygous_writer(
                self.chrom, genotypes,
                locus.start, locus.end,
                allele_count, len(read_indices),
                self.cooper_loci_info, self.ref,
                self.outhandle, '.', phased_reads,
                0, ALT_seqs, self.args.debug_mode,
                'HP', '.', hallele_counter,
                allele_range, [None], meth_info
            )
            genotyped = 1

        # --- no category — write fail if needed ---
        else:
            if skip_point == 0:
                write_fail_call(self, locus_key, skip_point)

        # --- cleanup ---
        del self.cooper_loci_data[locus_key]
        del self.cooper_loci_info[locus_key]

        return genotyped


def _hidden_path(out_file: str) -> str:
    """Return a hidden-file-prefixed path for thread output files."""
    path = Path(out_file)
    return str(path.parent / f'.{path.name}')
