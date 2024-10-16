from locus_utils import process_locus
from cstag_utils import parse_cstag
from cigar_utils import parse_cigar_tag
from operation_utils import update_homopolymer_coords
from genotype_utils import analyse_genotype
from vcf_writer import *

import pysam
import sys

def locus_processor(global_loci_keys, global_loci_ends, global_loci_variations, global_read_variations, global_snp_positions, prev_reads, sorted_global_snp_list, maxR, minR, ref, Chrom, global_loci_info, out, level_split, snpQ, snpC, snpD, snpR, phasingR):

    genotyped_loci = 0
    popped    = global_loci_ends.pop(0)
    locus_key = global_loci_keys.pop(0)
    allele_lengths = tuple()

    if locus_key in global_loci_variations:
        

        if sorted_global_snp_list == []:
            sorted_global_snp_list = sorted(list(global_snp_positions.keys()))
        

        prev_reads, homozygous, ambiguous, homozygous_allele, reads_of_homozygous, hallele_counter, skip_point, max_limit = process_locus(locus_key, global_loci_variations, global_read_variations, global_snp_positions, prev_reads, sorted_global_snp_list, maxR, minR)

        if homozygous:
            vcf_homozygous_writer(ref, Chrom, locus_key, global_loci_info, homozygous_allele, global_loci_variations, reads_of_homozygous, out)
            genotyped_loci += 1
        elif ambiguous:
            state, skip_point = analyse_genotype(Chrom, locus_key, global_loci_info, global_loci_variations, global_read_variations, global_snp_positions, hallele_counter, ref, out, sorted_global_snp_list, level_split, snpQ, snpC, snpD, snpR, phasingR, maxR, max_limit)
            if state: genotyped_loci += 1
            elif skip_point == 0: print('Locus skipped due to insignificant snps at the level of read split.')
            elif skip_point == 1: print('Locus skipped due to less read contribution of Significant snps.')
            elif skip_point == 2: print('Locus skipped due to less read contribution in the phased clusters.')
            else: print('Locus skipped due to less number of significant snps based on user\'s parameter.')
        else:
            if skip_point == 0:
                # print('Locus skipped due to minimum number of supporting reads')
                vcf_fail_writer(Chrom, locus_key, global_loci_info, ref, out, len(prev_reads), skip_point)
            # else:
            #     print('Locus skipped due to maximum number of supporting reads')
                
        del global_loci_variations[locus_key]
        
    return genotyped_loci




def cooper(bam_file, tbx_file, ref_file, aln_format, contigs, mapq_threshold, outfile, seq_tech, level_split, snpQ, snpC, snpD, maxR, minR, snpR, phasingR, tidx):
    # this function iterates through each contig and processes the genotypes for each locus
    tbx  = pysam.Tabixfile(tbx_file)
    bam  = pysam.AlignmentFile(bam_file, aln_format)
    ref  = pysam.FastaFile(ref_file)

    if tidx!=-1:
        out = open(f'{outfile}_thread_{tidx}.vcf', 'w')
    else: out = open(f'{outfile}.vcf', 'w')
    vcf_writer(out)

    for contig in contigs:

        Chrom, Start, End = contig

        # print(f"\nProcessing contig {contig}..\n", file=sys.stderr)
        total_loci = 0
        end_coord = 0
        genotyped_loci_count = 0
        for row in tbx.fetch(Chrom, Start[0], End[1]):
            row = row.split('\t')
            if (total_loci == 0) and (int(row[2]) != Start[1]): continue
            total_loci += 1
            end_coord = int(row[2])
            if End[0] == int(row[1]): break
        print(f"{Chrom} {Start} {End} Total loci = ", total_loci, file=sys.stderr)

        homozygous = False
        ambiguous = False

        global_snp_positions = dict()       # tracking the encountered SNPs
        global_read_variations = {}         # tracking the variations on each read
        global_loci_variations = {}         # tracking the variation for each locus
        global_loci_info = {}               # saving the information of each loci

        # tracking the loci
        global_loci_ends = []; global_loci_keys    = []        
        global_read_ends = []; global_read_indices = []

        prev_reads = set()
        sorted_global_snp_list = []

        read_index = 0
        for read in bam.fetch(Chrom, Start[0], End[1]):
        

            # skip read with low mapping quality
            if read.mapping_quality < mapq_threshold:
                continue

            read_chrom = read.reference_name
            read_start = read.reference_start
            read_end   = read.reference_end

            while len(global_loci_ends) > 0 and read_start > global_loci_ends[0]:

                genotyped_loci_count += locus_processor(global_loci_keys, global_loci_ends, global_loci_variations, global_read_variations, global_snp_positions, prev_reads, sorted_global_snp_list, maxR, minR, ref, Chrom, global_loci_info, out, level_split, snpQ, snpC, snpD, snpR, phasingR)
                

            while len(global_read_ends) > 0 and read_start > global_read_ends[0]:
                # if the read is beyond the end of the first read that was tracked
                if len(global_loci_ends) > 0 and global_read_ends[0] > global_loci_ends[0]:
                    # if the initial read useful for the first locus being tracked then it is retained
                    break
                else:

                    # remove the read information if the current read is beyond the first read and the locus
                    popped = global_read_ends.pop(0)
                    rindex = global_read_indices.pop(0)
                    if rindex in global_read_variations:
                        for pos in global_read_variations[rindex]['snps']:
                            if pos in global_snp_positions:
                                global_snp_positions[pos]['cov'] -= 1
                        del_snps = []
                        for pos in global_snp_positions: #!! modify this to terminate iteration after some condition
                            if global_snp_positions[pos]['cov'] == 0:
                                del_snps.append(pos)
                                sorted_global_snp_list.remove(pos)
                        for snp in del_snps:
                            del global_snp_positions[snp]
                        del global_read_variations[rindex]

                        if rindex in prev_reads: prev_reads.remove(rindex)


            # if the read is beyond the last locus in the bed file the loop stops
            if read_start > end_coord:
                while len(global_loci_ends) > 0:
                    genotyped_loci_count += locus_processor(global_loci_keys, global_loci_ends, global_loci_variations, global_read_variations, global_snp_positions, prev_reads, sorted_global_snp_list, maxR, minRref, Chrom, global_loci_info, out, level_split, snpQ, snpC, snpD, snpR, phasingR)
                # process the loci left in global_loci_variation
                break

            # information locally saved for each read and all the loci it covers
            read_loci_variations = {}
            # set of homopolymer positions within the reference part that is covered by the read
            homopoly_positions = {}

            # repeat loci covered by the read
            loci_coords = []; loci_keys = []
            start_check = 0
            for row in tbx.fetch(read_chrom, read_start, read_end):
                
                # adjust read start and end based on soft and hard clippings
                # soft and hard clippings do not consume the reference bases

                row = row.split('\t')
                locus_start = int(row[1]);  locus_end = int(row[2]); locus_len = locus_end-locus_start

                

                if locus_start>Start[0]:
                    pass
                elif locus_start==Start[0] and locus_end==Start[1]:
                    start_check = 1
                elif start_check!=1:
                    continue
            
                if locus_start>End[0]: break
                
                # a location which is tracked with absurdly large number of SNPs
                if locus_start >= 143242921 and locus_end <= 143276906: continue

                # if only the read completely covers the repeat
                if ( locus_start >= read_start ) & ( locus_end <= read_end ):
                    loci_coords.append((locus_start, locus_end))
                    locus_key = f'{read_chrom}:{locus_start}-{locus_end}'
                    loci_keys.append(locus_key)
                    read_loci_variations[locus_key] = {'halen': locus_len, 'alen': locus_len, 'rlen': locus_len}

                    if locus_key not in global_loci_variations:
                        global_loci_variations[locus_key] = {'rlen': locus_len, 'reads': [], 'read_allele': {}}
                        global_loci_info[locus_key] = row

                        # adding the locus key when it is first encountered
                        global_loci_ends.append(locus_end)
                        global_loci_keys.append(locus_key)

            # if no repeats are covered by the read
            if len(loci_coords) == 0: continue
        
            read_index += 1
            read_quality = read.query_qualities
            cigar_tuples = read.cigartuples
            read_sequence = read.query_sequence

            cigar_one = cigar_tuples[0]

            global_read_ends.append(read_end)
            global_read_indices.append(read_index)
            global_read_variations[read_index] = {'s': read_start, 'e': read_end, 'snps': set(), 'dels': set()}

            for each_coords in loci_coords:
                update_homopolymer_coords(ref.fetch(Chrom, each_coords[0]-100, each_coords[1]+100), each_coords[0]-100, homopoly_positions)

            if read.has_tag('cs'):
                cs_tag = read.get_tag('cs')
                parse_cstag(read_index, cs_tag, read_start, loci_keys, loci_coords, read_loci_variations, homopoly_positions, global_read_variations, global_snp_positions, read_quality, cigar_one, sorted_global_snp_list)
            else :
                parse_cigar_tag(read_index, cigar_tuples, read_start, loci_keys, loci_coords, read_loci_variations,
                                homopoly_positions, global_read_variations, global_snp_positions, read_sequence, read, ref, read_quality, sorted_global_snp_list)

            for locus_key in read_loci_variations:
                global_loci_variations[locus_key]['reads'].append(read_index)
                global_loci_variations[locus_key]['read_allele'][read_index] = (read_loci_variations[locus_key]['halen'], read_loci_variations[locus_key]['alen'])

        while len(global_loci_ends) > 0:
            genotyped_loci_count += locus_processor(global_loci_keys, global_loci_ends, global_loci_variations, global_read_variations, global_snp_positions, prev_reads, sorted_global_snp_list, maxR, minR, ref, Chrom, global_loci_info, out, level_split, snpQ, snpC, snpD, snpR, phasingR)
                
        print(f'\nTotal genotyped loci = {genotyped_loci_count}\n', file=sys.stderr)

    bam.close()
    ref.close()
    tbx.close()