import numpy as np
import pysam
import logging

from weakref import ref
from tqdm import tqdm
from sortedcontainers import SortedList
from collections import deque

from ATARVA.structures import ReadLocusInfo, LocusInfo, ReadInfo, LocusVariation, ExtendedRead
from ATARVA.vcf_writer import vcf_writer, vcf_homozygous_writer, vcf_heterozygous_writer, vcf_fail_writer
from ATARVA.operation_utils import record_homopolymers, clean_eqsign_readseq
from ATARVA.cstag_utils import parse_cstag
from ATARVA.cigar_utils import parse_cigar
from ATARVA.sub_operation_utils import mm_tag_extract
from ATARVA.locus_utils import process_locus
from ATARVA.consensus import consensus_seq_poa
from ATARVA.genotype_utils import analyse_genotype
# , methylation_calc


class Cooper:
    """
    Object to independently handle genotyping a set of regions in a single bam
    """

    tbx = None
    bam = None
    ref = None
    args = None
    logger = None
    thread_idx = None
    sample_idx = None
    karyotype = None
    haploid = None
    chrom = None

    outfile = None; outhandle = None
    logfile = None; loghandle = None

    cooper_nloci = 0
    range_nloci = []
    disable_progress = False

    cooper_snp_data = dict()       # tracking the encountered SNPs
    cooper_read_data = dict()      # tracking the variations on each read
    cooper_loci_data = dict()      # tracking the variation for each locus
    cooper_loci_info = dict()   # saving the information of each loci

    # tracking the loci
    cooper_loci_ends = deque()
    cooper_loci_keys = deque()
    
    cooper_read_ends    = deque()
    cooper_read_indices = deque()
    prev_reads = set()

    cooper_sorted_snps = SortedList()
    cooper_insert_positions = set() # tracking the repeat insertion positions globally to avoiding same insertion into multiple loci


    def reinitialize_globals(self):
        self.cooper_snp_data = dict()       # tracking the encountered SNPs
        self.cooper_read_data = dict()      # tracking the variations on each read
        self.cooper_loci_data = dict()      # tracking the variation for each locus
        self.cooper_loci_info = dict()      # saving the information of each loci

        # tracking the loci
        self.cooper_loci_ends = deque()
        self.cooper_loci_keys = deque()

        self.cooper_read_ends    = deque()
        self.cooper_read_indices = deque()
        self.prev_reads = set()

        self.cooper_sorted_snps = SortedList()
        self.cooper_sorted_ins_rpos_set = set()

        self.chrom = None


    def __init__(self, bam_file, region_ranges, args, out_file, sample_idx, thread_idx):
        """
        initialize cooper object
        :param bam_file: path to the bam file
        :param region_ranges: list of tuples of the form (chrom, (start1, end1), (start2, end2))
        :param args: arguments from the command line
        :param out_file: path to the output file
        :param sample_idx: index of the sample in the bam list
        :param thread_idx: index of the thread for parallel processing
        """

        self.tbx = pysam.Tabixfile(args.regions)
        self.bam = pysam.AlignmentFile(bam_file, args.format)
        self.ref = pysam.FastaFile(args.fasta)

        self.args = args
        self.sample_idx = sample_idx
        self.thread_idx = thread_idx
        self.karyotype = args.karyotype[sample_idx]

        if thread_idx == -1 or thread_idx == 0:
            self.outfile = f'{out_file}.vcf'
            self.logfile = f'{out_file}_debug.log'
        else:
            idx = args.outfile.rfind('/')
            hid_outfile = args.outfile[:idx+1] + '.' + args.outfile[idx+1:]
            self.outfile = f'{hid_outfile}_thread_{thread_idx}.vcf'
            self.logfile = f'{hid_outfile}_debug_{thread_idx}.log'

        self.outhandle = open(self.outfile, 'w')

        if thread_idx == -1 or thread_idx == 0:
            vcf_writer(self.outhandle, self.bam, bam_file.split("/")[-1].split('.')[0])
        
        if args.debug_mode:
            with open(self.logfile , 'w'): pass
            logging.basicConfig(filename=self.logfile, level=logging.DEBUG,
                                format='%(levelname)s - %(message)s')
            self.logger = logging.getLogger("ATaRVa_logger")

        if args.amplicon: args.haplotag = None

        self.disable_progress = thread_idx != -1

        self.range_nloci = []
        for region_range in region_ranges:
            chrom, first_coords, last_coords = region_range

            range_nloci = 0
            for row in self.tbx.fetch(chrom, first_coords[0], last_coords[1]):
                row = row.split('\t')
                row_start, row_end = int(row[1]), int(row[2])
                if range_nloci == 0 and row_start != first_coords[0]: continue
                range_nloci += 1
                if last_coords[0] == row_start:
                    self.range_nloci.append(range_nloci)
                    self.cooper_nloci += 1
                    break

        self.progress_bar = tqdm(total= self.cooper_nloci, disable = self.disable_progress,
                                 desc="Processing ", ascii="_>", ncols=75,
                                 bar_format="{l_bar}{bar}{n_fmt}/{total_fmt}")
        for cidx, region_range in enumerate(region_ranges):
            print(f"Processing region {region_range} in sample {bam_file.split('/')[-1]} and thread {thread_idx}\n\n")
            self.reinitialize_globals()
            self.cooper_readmode(region_range, cidx)


    def cooper_readmode(self, region_range, cidx):
        """
        Genotyping range of loci (limited to a chromosome) by looping through reads. 

        :param region_range: tuple of the form (chrom, (start1, end1), (start2, end2))
        :param cidx: index of the region range in the list of region ranges
        """

        chrom, first_coords, last_coords = region_range
        self.chrom = chrom
        self.haploid = (chrom in {'chrX', 'chrY', 'X', 'Y'}) and self.karyotype
        region_range_start = first_coords[0]
        region_range_end = last_coords[1]

        genotyped_loci_count = 0

        if not self.disable_progress:
            tqdm.write(f"> {chrom} {region_range_start} {region_range_end} Total loci =  {self.range_nloci[cidx]}")

        read_index = 0      # custom index for reads

        # try realigning the read with the sofclips included

        for read in self.bam.fetch(chrom, region_range_start, region_range_end):

            # skip read with low mapping quality
            if read.mapping_quality < self.args.map_qual:
                continue

            read = ExtendedRead.from_read(read)

            while self.cooper_loci_ends and read.ref_start > self.cooper_loci_ends[0]:
                genotype_status = self.locus_processor()
                genotyped_loci_count += genotype_status[0]
                self.progress_bar.update(1)

            while self.cooper_read_ends and read.ref_start > self.cooper_read_ends[0]:
                # if the read is beyond the end of the first read that was tracked
                if self.cooper_loci_ends and self.cooper_read_ends[0] > self.cooper_loci_ends[0]:
                    # if the initial read useful for the first locus being tracked then it is retained
                    break

                else:
                    # remove the read information if the current read is beyond the first read and the locus
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

                    del_ins_pos_idx = 0
                    list_rpos = sorted(sorted_global_ins_rpos_set)
                    for i in list_rpos:
                        del_ins_pos_idx+=1
                        if i > read_end: break
                    del list_rpos[:del_ins_pos_idx]
                    sorted_global_ins_rpos_set = set(list_rpos)                    


            # if the read is beyond the last locus in the bed file the loop stops
            if read.ref_start > region_range_end:
                while self.cooper_loci_ends:
                    genotype_status = self.locus_processor()
                    genotyped_loci_count += genotype_status
                    self.progress_bar.update(1)
                # process the loci left in global_loci_variation
                break

            for row in self.tbx.fetch(chrom, read.ref_start, read.ref_end):

                # adjust read start and end based on soft and hard clippings
                # soft and hard clippings do not consume the reference bases

                row = row.split('\t')
                locus_start = int(row[1])
                locus_end = int(row[2])
                locus_len = locus_end - locus_start

                # this just stores the loci encountered in the read
                read.loci.append((locus_start, locus_end))

                if (locus_start >= first_coords[0]) and (locus_end <= last_coords[1]):
                    if locus_start == first_coords[0]:
                        if locus_end == first_coords[1]: pass
                        else: continue
                    pass
                elif locus_start < first_coords[0]:
                    continue
                elif locus_start >= last_coords[0]: break

                passed_loci = False # if the loci passed in normal or amplicon mode, then write it in global variables

                # if only the read completely covers the repeat
                if read.ref_start <= locus_start and locus_end <= read.ref_end:
                    passed_loci = True

                    left_flank  = min(self.args.flank, locus_start - read.ref_start)
                    right_flank = min(self.args.flank, read.ref_end - locus_end)

                    read.left_flanks.append(left_flank)
                    read.right_flanks.append(right_flank)

                    # NOTE: storing the locus coordinates with flanks included
                    read.loci_coords.append((locus_start - left_flank, locus_end + right_flank))

                if passed_loci:
                    locus_key = f'{chrom}:{locus_start}-{locus_end}'
                    read.loci_keys.append(locus_key)
                    read.loci_data[locus_key] = ReadLocusInfo(halen=0, alen=0, rlen = locus_len, seq=[])

                    if locus_key not in self.cooper_loci_data:
                        self.cooper_loci_data[locus_key] = LocusVariation()
                        self.cooper_loci_info[locus_key] = LocusInfo(chrom=chrom, start=locus_start, end=locus_end, motif=row[3])

                        # adding the locus key when it is first encountered
                        self.cooper_loci_ends.append(locus_end)
                        self.cooper_loci_keys.append(locus_key)


            # if no repeats are covered by the read
            if len(read.loci_coords) == 0: continue

            read_index += 1
            read.index = read_index

            tmp_qpos = 0
            for op, length in read.cigartuples:
                if (op == 0) or (op == 7):
                    tmp_seq = read.sequence[tmp_qpos :tmp_qpos + length] 
                    break
                elif (op == 1) or (op == 4) or (op == 8): tmp_qpos += length 
                
            if '=' in tmp_seq:
                read_sequence = clean_eqsign_readseq(read.chrom, read.ref_start, read.cigartuples, read_sequence, ref)

            self.cooper_read_ends.append(read.ref_end)
            self.cooper_read_indices.append(read.index)
            self.cooper_read_data[read.index] = ReadInfo(start = read.ref_start,
                                                         end = read.ref_end, snps=set(),
                                                         dels = [], methyl = [],
                                                         qual = read.mean_qual,
                                                         left_flank = read.left_flanks,
                                                         right_flank = read.right_flanks)

            if self.args.haplotag:
                read.haplotag[0] = read.has_tag(self.args.haplotag)
            if read.haplotag[0]: read.haplotag[1] = read.get_tag(self.args.haplotag)

            if read.has_tag('cs'):
                read_modified_bases = list(read.modified_bases.items()) if read.modified_bases is not None else []
                if len(read_modified_bases)>0:
                    for mods in read_modified_bases:
                        if (mods[0][0]=='C') and (mods[0][2]=='m'):
                            mm_tag_extract(read, mods[1], 0, not(mods[0][1])) # last arg is bool value for strand state; forward = True, reverse = False
                            self.cooper_read_data[read_index].methyl = read.methyl_range
                            break
            else :
                parse_cigar(self, read)
                read_modified_bases = list(read.modified_bases.items()) if read.modified_bases is not None else []
                if len(read_modified_bases) > 0:
                    for mods in read_modified_bases:
                        if (mods[0][0]=='C') and (mods[0][2]=='m'):
                            mm_tag_extract(read, mods[1], 0, not(mods[0][1])) # last arg is bool value for strand state; forward = True, reverse = False
                            self.cooper_read_data[read_index].methyl = read.methyl_range
                            break

            for locus_key in read.loci_data:

                if self.args.amplicon:
                    if not read.loci_data[locus_key].seq:
                        continue

                self.cooper_loci_data[locus_key].reads.append(read.index)
                self.cooper_loci_data[locus_key].read_alens[read.index] = [read.loci_data[locus_key].halen, read.loci_data[locus_key].alen]
                self.cooper_loci_data[locus_key].read_seqs[read.index] = read.loci_data[locus_key].seq
                self.cooper_loci_data[locus_key].read_haplotags.append(read.haplotag[1])

        while self.cooper_loci_ends:
            genotype_status = self.locus_processor()
            genotyped_loci_count += genotype_status
            self.progress_bar.update(1)
                
        # if not dwrite: tqdm.write(f'\nTotal genotyped loci = {genotyped_loci_count} out of {tot_loci_list[cidx]} in {Chrom} {first_coords[0]}-{last_coords[1]}\n')
    
        self.reinitialize_globals() # reinitialize the global variables for the next region
        self.bam.close()
        self.ref.close()
        self.tbx.close()
        self.out.close()


    def locus_processor(self):
        """
        Genotypes loci based on the read data collected
        """

        genotyped_loci = 0
        popped    = self.cooper_loci_ends.popleft()
        
        locus_key = self.cooper_loci_keys.popleft()
        locus = self.cooper_loci_info[locus_key]

        neighbors = []
        for row in self.tbx.fetch(self.chrom, locus.start - self.args.flank, locus.end + self.args.flank):
            row = row.split('\t')
            neighbors.append( ( int(row[1]), int(row[2]) ) )

        if locus_key in self.cooper_loci_data:

            if len(self.cooper_sorted_snps) == 0:
                for snp_pos in self.cooper_snp_data.keys():
                    self.cooper_sorted_snps.add(snp_pos)

            category, homozygous_allele, reads_of_homozygous, hallele_counter, skip_point, haplotypes, homozygous_alens = process_locus(self, locus_key, neighbors)
            read_seqs = self.cooper_loci_data[locus_key].read_seqs

            if category == 1:
                if homozygous_allele != locus.length:
                    seqs = [seq for seq in [read_seqs[read_id][0] for read_id in reads_of_homozygous] if seq!='']
                    if len(seqs) > 0:
                        ALT = consensus_seq_poa(seqs)
                        homozygous_allele = len(ALT)
                    else:
                        ALT = '<DEL>'
                        homozygous_allele = 0
                else: ALT = '.'

                lower, upper = [round(x) for x in np.percentile(np.array(homozygous_alens), [2.5, 97.5])]

                # if ALT != '.':
                #     meth_info = methylation_calc(reads_of_homozygous, self.cooper_loci_data, locus_key, ALT)
                # else:
                #     meth_info = methylation_calc(reads_of_homozygous, self.cooper_loci_data, locus_key, ref.fetch(locus.chrom, locus.start, locus.end))

                if self.haploid:
                    allele_range = f'{lower}-{upper}'
                else:
                    allele_range = f'{lower}-{upper},{lower}-{upper}'

                vcf_homozygous_writer(self, locus_key, homozygous_allele, len(reads_of_homozygous), len(reads_of_homozygous),
                                      ALT, '.', hallele_counter, allele_range, None, meth_info)
                genotyped_loci += 1
            elif category == 2:
                state, skip_point = analyse_genotype(Chrom, locus_key, global_loci_info, global_loci_variations, global_read_variations, global_snp_positions,
                                                     hallele_counter, ref, out, sorted_global_snp_list, snpQ, snpC, snpD, snpR, phasingR, reads_of_homozygous,
                                                     male, log_bool, decomp, amplicon, somatic)
                if state: genotyped_loci += 1
                else:
                    skip_messages = {
                        0: 'Locus skipped due to insignificant snps at the level of read split.',
                        1: 'Locus skipped due to less read contribution of Significant snps.',
                        2: 'Locus skipped due to less read contribution in the phased clusters.',
                        6: f'Locus {locus_key} skipped due to wide distribution of alleles with one read supporting to it.',
                    }
                    tqdm.write(skip_messages.get(skip_point, 'Locus skipped due to less number of significant snps based on user\'s parameter.'))
            elif category == 3:
                genotypes = []
                allele_count = {}
                ALT_seqs = []
                phased_read = []
                alen_list = []
                meth_info = []
                for hap_reads in haplotypes:
                    phased_read.append(len(hap_reads))
                    seqs = [seq for seq in [read_seqs[read_id][0] for read_id in hap_reads] if seq!='']
                    alen_list.append([len(read_seqs[read_id][0]) for read_id in hap_reads])
                    if len(seqs)>0:
                        ALT = consensus_seq_poa(seqs)
                        allele_length = len(ALT)
                    else:
                        ALT = '<DEL>'
                        allele_length = 0

                    ALT_seqs.append(ALT)
                    genotypes.append(allele_length)

                    if allele_length not in allele_count:
                        allele_count[allele_length] = len(hap_reads)
                    else:
                        allele_count[str(allele_length)] = len(hap_reads)

                    # meth_info.append(methylation_calc(hap_reads, global_loci_variations, locus_key))

                lower1,upper1 = confidence_interval(alen_list[0])
                lower2,upper2 = confidence_interval(alen_list[1])
                allele_range = f'{lower1}-{upper1},{lower2}-{upper2}'

                vcf_heterozygous_writer(Chrom, genotypes, lstart, lend, allele_count, len(reads_of_homozygous), global_loci_info,
                                        ref, out, '.', phased_read, 0, ALT_seqs, log_bool, 'HP', decomp, hallele_counter,
                                        allele_range, [None], meth_info)
                genotyped_loci += 1
            else:
                if skip_point == 0:
                    vcf_fail_writer(Chrom, locus_key, global_loci_info, ref, out, len(prev_reads), skip_point)

                    
            del global_loci_variations[locus_key]
            
        return genotyped_loci
