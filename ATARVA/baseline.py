from ATARVA.locus_utils import process_locus
from ATARVA.cstag_utils import parse_cstag
from ATARVA.cigar_utils import parse_cigar_tag
from ATARVA.operation_utils import update_homopolymer_coords
from ATARVA.genotype_utils import analyse_genotype
from ATARVA.vcf_writer import *

import pysam
import sys
import logging

def locus_processor(global_loci_keys, global_loci_ends, global_loci_variations, global_read_variations, global_snp_positions, prev_reads, sorted_global_snp_list, maxR, minR, ref, Chrom, global_loci_info, out, snpQ, snpC, snpD, snpR, phasingR, tbx, flank, sorted_global_ins_rpos_set, log_bool, logger, male):

    genotyped_loci = 0
    popped    = global_loci_ends.pop(0)
    locus_key = global_loci_keys.pop(0)
    lstart = int(locus_key[locus_key.index(':')+1 : locus_key.index('-')])
    lend = int(locus_key[locus_key.index('-')+1:])
    near_by_loci = []
    for row in tbx.fetch(Chrom, lstart-flank, lend+flank):
        row = row.split('\t')
        near_by_loci.append( ( int(row[1]), int(row[2]) ) )

    if locus_key in global_loci_variations:
        

        if sorted_global_snp_list == []:
            sorted_global_snp_list = sorted(list(global_snp_positions.keys()))
        

        prev_reads, homozygous, ambiguous, homozygous_allele, reads_of_homozygous, hallele_counter, skip_point, max_limit = process_locus(locus_key, global_loci_variations, global_read_variations, global_snp_positions, prev_reads, sorted_global_snp_list, maxR, minR, global_loci_info, near_by_loci, sorted_global_ins_rpos_set, Chrom, lstart, lend, ref, log_bool, logger)

        if homozygous:
            vcf_homozygous_writer(ref, Chrom, locus_key, global_loci_info, homozygous_allele, global_loci_variations, len(reads_of_homozygous), out, reads_of_homozygous)
            genotyped_loci += 1
        elif ambiguous:
            state, skip_point = analyse_genotype(Chrom, locus_key, global_loci_info, global_loci_variations, global_read_variations, global_snp_positions, hallele_counter, ref, out, sorted_global_snp_list, snpQ, snpC, snpD, snpR, phasingR, maxR, max_limit, male)
            if state: genotyped_loci += 1
            elif skip_point == 0: print('Locus skipped due to insignificant snps at the level of read split.')
            elif skip_point == 1: print('Locus skipped due to less read contribution of Significant snps.')
            elif skip_point == 2: print('Locus skipped due to less read contribution in the phased clusters.')
            elif skip_point == 6: print(f'Locus {locus_key} skipped due to wide distribution of alleles with one read supporting to it.') # add vcf_fail_writer with new FILTER tag
            else: print('Locus  skipped due to less number of significant snps based on user\'s parameter.')
        else:
            if skip_point == 0:
                # print('Locus skipped due to minimum number of supporting reads')
                vcf_fail_writer(Chrom, locus_key, global_loci_info, ref, out, len(prev_reads), skip_point)
            # else:
            #     print('Locus skipped due to maximum number of supporting reads')
                
        del global_loci_variations[locus_key]
        
    return genotyped_loci

# def flank_adjustment(loci_coords, exact_loci_coords, right_flank_list, left_flank_list): # incorrect code
#     # Adjusting flanks based on the nearby repeat region
#     for idx,locus in enumerate(loci_coords):
#         current_left = left_flank_list[idx]
#         current_right = right_flank_list[idx]
#         new_left_flanks = []
#         new_right_flanks = []
#         for exact_locus in exact_loci_coords:
#             if exact_locus[0]>locus[1]: break
#             elif locus[1]>exact_locus[0]:  # after skipping locus (or) after came across the current locus, adjust the right flank dist
#                 new_right_flanks.append(locus[1]-exact_locus[0] if (locus[1]-exact_locus[0]) < current_right else current_right)
#             elif locus[0]<exact_locus[1]: # adjust the left flank dist
#                 new_left_flanks.append(exact_locus[1] - locus[0] if (exact_locus[1] - locus[0]) < current_left else current_left) # add adjusted left flank dist for each overlapping exact repeat regions with cureent repeat region's flank on left. Zero is not to adjust any dist 
#         if new_left_flanks!=[]:
#             left_flank_list[idx] -= max(new_left_flanks) # take the max value of adjusted value, so that it wont affect other exact repeats
#         if new_right_flanks!=[]:
#             right_flank_list[idx] -= max(new_right_flanks) # take the max value of adjusted value, so that it wont affect other exact repeats

#     return [right_flank_list, left_flank_list]

def flank_adjustment(loci_coords, exact_loci_coords, right_flank_list, left_flank_list, tracker_idx):
    # Adjusting flanks based on the nearby repeat region
    for idx,locus in enumerate(loci_coords):
        current_left = left_flank_list[idx]
        current_right = right_flank_list[idx]
        new_left_flanks = []
        new_right_flanks = []
        tr_idx = tracker_idx[idx]
        for e_idx,exact_locus in enumerate(exact_loci_coords):
            if exact_locus[0]>locus[1]:
                break
            if e_idx == tr_idx:
                continue
            if (e_idx>tr_idx) and (locus[1]>exact_locus[0]):  # after skipping locus (or) after came across the current locus, adjust the right flank dist
                new_right_flanks.append(locus[1]-exact_locus[0] if (locus[1]-exact_locus[0]) < current_right else current_right)
            if (e_idx<tr_idx) and (locus[0]<exact_locus[1]): # adjust the left flank dist
                new_left_flanks.append(exact_locus[1] - locus[0] if (exact_locus[1] - locus[0]) < current_left else current_left) # add adjusted left flank dist for each overlapping exact repeat regions with cureent repeat region's flank on left. Zero is not to adjust any dist 
        if new_left_flanks!=[]:
            left_flank_list[idx] -= max(new_left_flanks) # take the max value of adjusted value, so that it wont affect other exact repeats
        if new_right_flanks!=[]:
            right_flank_list[idx] -= max(new_right_flanks) # take the max value of adjusted value, so that it wont affect other exact repeats

    return [right_flank_list, left_flank_list]


def cooper(bam_file, tbx_file, ref_file, aln_format, contigs, mapq_threshold, outfile, seq_tech, snpQ, snpC, snpD, maxR, minR, snpR, phasingR, tidx, flank, log_bool, karyotype):
    # this function iterates through each contig and processes the genotypes for each locus
    tbx  = pysam.Tabixfile(tbx_file)
    bam  = pysam.AlignmentFile(bam_file, aln_format)
    ref  = pysam.FastaFile(ref_file)

    logger = 0
    if tidx!=-1:
        if tidx==0:
            out = open(f'{outfile}.vcf', 'w')
            vcf_writer(out, bam, bam_file.split("/")[-1].split('.')[0])
            log_name = f'{outfile}_debug.log'
        else:
            out = open(f'{outfile}_thread_{tidx}.vcf', 'w')
            log_name = f'{outfile}_debug_{tidx}.log'

    else:
        out = open(f'{outfile}.vcf', 'w')
        vcf_writer(out, bam, bam_file.split("/")[-1].split('.')[0])
        log_name = f'{outfile}_debug.log'

    if log_bool:
        with open(log_name, 'w'):
            pass
        logging.basicConfig(
            filename=log_name,
            level=logging.DEBUG,
            format='%(levelname)s - %(message)s'
        )
        logger = logging.getLogger("MyLogger")
    
    for contig in contigs:

        Chrom, Start, End = contig
        male = False
        if (Chrom in {'chrX', 'chrY', 'X', 'Y'}) and karyotype: male = True

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
        print(f"> {Chrom} {Start} {End} Total loci = ", total_loci, file=sys.stderr)

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
        sorted_global_ins_rpos_set = set()

        read_index = 0
        for read in bam.fetch(Chrom, Start[0], End[1]):
        

            # skip read with low mapping quality
            if read.mapping_quality < mapq_threshold:
                continue

            read_chrom = read.reference_name
            read_start = read.reference_start
            read_end   = read.reference_end

            while len(global_loci_ends) > 0 and read_start > global_loci_ends[0]:

                genotyped_loci_count += locus_processor(global_loci_keys, global_loci_ends, global_loci_variations, global_read_variations, global_snp_positions, prev_reads, sorted_global_snp_list, maxR, minR, ref, Chrom, global_loci_info, out, snpQ, snpC, snpD, snpR, phasingR, tbx, flank, sorted_global_ins_rpos_set, log_bool, logger, male)
                

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
                        del del_snps

                        if rindex in prev_reads: prev_reads.remove(rindex)
                    del_ins_pos_idx = 0
                    list_rpos = sorted(sorted_global_ins_rpos_set)
                    for i in list_rpos:
                        del_ins_pos_idx+=1
                        if i > popped: break
                    del list_rpos[:del_ins_pos_idx]
                    sorted_global_ins_rpos_set = set(list_rpos)
                    


            # if the read is beyond the last locus in the bed file the loop stops
            if read_start > end_coord:
                while len(global_loci_ends) > 0:
                    genotyped_loci_count += locus_processor(global_loci_keys, global_loci_ends, global_loci_variations, global_read_variations, global_snp_positions, prev_reads, sorted_global_snp_list, maxR, minR, ref, Chrom, global_loci_info, out, snpQ, snpC, snpD, snpR, phasingR, tbx, flank, sorted_global_ins_rpos_set, log_bool, logger, male)
                # process the loci left in global_loci_variation
                break

            # information locally saved for each read and all the loci it covers
            read_loci_variations = {}
            # set of homopolymer positions within the reference part that is covered by the read
            homopoly_positions = {}

            # repeat loci covered by the read
            loci_coords = []; loci_keys = []
            left_flank_list = []; right_flank_list = []

            for row in tbx.fetch(read_chrom, read_start, read_end):
                
                # adjust read start and end based on soft and hard clippings
                # soft and hard clippings do not consume the reference bases

                row = row.split('\t')
                locus_start = int(row[1]);  locus_end = int(row[2]); locus_len = locus_end-locus_start


                if (locus_start>=Start[0]) and (locus_end<=End[1]):
                    if locus_start==Start[0]:
                        if locus_end==Start[1]: pass
                        else: continue
                    pass
                elif locus_start<Start[0]:
                    continue
                elif locus_start>=End[0]: break
                

                # if only the read completely covers the repeat
                if ( locus_start >= read_start ) & ( locus_end <= read_end ):
                    left_flank = flank; right_flank = flank
                    if (locus_start - flank) < read_start:
                        left_flank = locus_start - read_start
                    if (locus_end + flank) > read_end:
                        right_flank  = read_end - locus_end
                    left_flank_list.append(left_flank)
                    right_flank_list.append(right_flank)

                    loci_coords.append((locus_start - left_flank, locus_end + right_flank))

                    locus_key = f'{read_chrom}:{locus_start}-{locus_end}'
                    loci_keys.append(locus_key)
                    read_loci_variations[locus_key] = {'halen': locus_len, 'alen': locus_len, 'rlen': locus_len, 'seq': []}

                    if locus_key not in global_loci_variations:
                        global_loci_variations[locus_key] = {'rlen': locus_len, 'reads': [], 'read_allele': {}, 'read_sequence': {}}
                        global_loci_info[locus_key] = row

                        # adding the locus key when it is first encountered
                        global_loci_ends.append(locus_end)
                        global_loci_keys.append(locus_key)

            # if no repeats are covered by the read
            if len(loci_coords) == 0: continue

            # right_flank_list, left_flank_list = flank_adjustment(loci_coords, exact_loci_coords, right_flank_list, left_flank_list, tracker_idx)
            # for idx,locus in enumerate(passed_loci_coords):
            #     loci_coords[idx] = (locus[0]-left_flank_list[idx], locus[1]+right_flank_list[idx]) 

            read_index += 1
            read_quality = read.query_qualities
            cigar_tuples = read.cigartuples
            read_sequence = read.query_sequence

            cigar_one = cigar_tuples[0]

            global_read_ends.append(read_end)
            global_read_indices.append(read_index)
            global_read_variations[read_index] = {'s': read_start, 'e': read_end, 'snps': set(), 'dels': set()}

            # for each_coords in loci_coords:
            #     update_homopolymer_coords(ref.fetch(Chrom, each_coords[0]-100, each_coords[1]+100), each_coords[0]-100, homopoly_positions)

            if read.has_tag('cs'):
                cs_tag = read.get_tag('cs')
                parse_cstag(read_index, cs_tag, read_start, loci_keys, loci_coords, read_loci_variations, homopoly_positions, global_read_variations, global_snp_positions, read_sequence, read_quality, cigar_one, sorted_global_snp_list, left_flank_list, right_flank_list, male)
            else :
                parse_cigar_tag(read_index, cigar_tuples, read_start, loci_keys, loci_coords, read_loci_variations,
                                homopoly_positions, global_read_variations, global_snp_positions, read_sequence, read, ref, read_quality, sorted_global_snp_list, left_flank_list, right_flank_list, male)

            for locus_key in read_loci_variations:
                global_loci_variations[locus_key]['reads'].append(read_index)
                global_loci_variations[locus_key]['read_allele'][read_index] = [read_loci_variations[locus_key]['halen'], read_loci_variations[locus_key]['alen']]
                global_loci_variations[locus_key]['read_sequence'][read_index] = read_loci_variations[locus_key]['seq']

        while len(global_loci_ends) > 0:
            genotyped_loci_count += locus_processor(global_loci_keys, global_loci_ends, global_loci_variations, global_read_variations, global_snp_positions, prev_reads, sorted_global_snp_list, maxR, minR, ref, Chrom, global_loci_info, out, snpQ, snpC, snpD, snpR, phasingR, tbx, flank, sorted_global_ins_rpos_set, log_bool, logger, male)
                
        print(f'\nTotal genotyped loci = {genotyped_loci_count} out of {total_loci} in {Chrom} {Start[0]}-{End[1]}\n', file=sys.stderr)

    bam.close()
    ref.close()
    tbx.close()
    out.close()


def mini_cooper(bam_file, tbx_file, ref_file, aln_format, contigs, mapq_threshold, outfile, seq_tech, snpQ, snpC, snpD, maxR, minR, snpR, phasingR, tidx, flank, log_bool, karyotype):
    # this function iterates through each contig and processes the genotypes for each locus
    tbx  = pysam.Tabixfile(tbx_file)
    tbx2  = pysam.Tabixfile(tbx_file)
    bam  = pysam.AlignmentFile(bam_file, aln_format)
    ref  = pysam.FastaFile(ref_file)


    logger = 0
    if tidx!=-1:
        if tidx==0:
            out = open(f'{outfile}.vcf', 'w')
            vcf_writer(out, bam, bam_file.split("/")[-1].split('.')[0])
            log_name = f'{outfile}_debug.log'
        else:
            out = open(f'{outfile}_thread_{tidx}.vcf', 'w')
            log_name = f'{outfile}_debug_{tidx}.log'

    else:
        out = open(f'{outfile}.vcf', 'w')
        vcf_writer(out, bam, bam_file.split("/")[-1].split('.')[0])
        log_name = f'{outfile}_debug.log'

    if log_bool:
        logging.basicConfig(
            filename=log_name,
            level=logging.DEBUG,
            format='%(levelname)s - %(message)s'
        )
        logger = logging.getLogger("MyLogger")

    for contig in contigs:

        Chrom, Start, End = contig
        male = False
        if (Chrom in {'chrX', 'chrY', 'X', 'Y'}) and karyotype: male = True

        # print(f"\nProcessing contig {contig}..\n", file=sys.stderr)
        total_loci = 0
        genotyped_loci_count = 0
        for row in tbx.fetch(Chrom, Start[0], End[1]):
            row = row.split('\t')
            if (total_loci == 0) and (int(row[2]) != Start[1]): continue
            total_loci += 1
            if End[0] == int(row[1]): break
        print(f"> {Chrom} {Start} {End} Total loci =  {total_loci}", file=sys.stderr)

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
        sorted_global_ins_rpos_set = set()

        read_index = 0
        
        for row in tbx.fetch(Chrom, Start[0], End[1]):

            row = row.split('\t')
            locus_start = int(row[1]);  locus_end = int(row[2]); locus_len = locus_end-locus_start


            if (locus_start>=Start[0]) and (locus_end<=End[1]):
                if locus_start==Start[0]:
                    if locus_end==Start[1]: pass
                    else: continue
                pass
            elif locus_start<Start[0]:
                continue
            elif locus_start>=End[0]: break
            
           
            for read in bam.fetch(Chrom, int(row[1]), int(row[2])):
                if read.mapping_quality < mapq_threshold:
                    continue

                read_chrom = read.reference_name
                read_start = read.reference_start
                read_end   = read.reference_end
                
                # if only the read completely covers the repeat
                if not (( locus_start >= read_start ) & ( locus_end <= read_end )):
                    continue

                
                # information locally saved for each read and all the loci it covers
                read_loci_variations = {}
                # set of homopolymer positions within the reference part that is covered by the read
                homopoly_positions = {}
    
                # repeat loci covered by the read
                loci_coords = []; loci_keys = []
 
                left_flank_list = []; right_flank_list = []

                left_flank = flank; right_flank = flank
                if (locus_start - flank) < read_start:
                    left_flank = locus_start - read_start
                if (locus_end + flank) > read_end:
                    right_flank  = read_end - locus_end
                left_flank_list.append(left_flank)
                right_flank_list.append(right_flank)

                loci_coords.append((locus_start - left_flank, locus_end + right_flank))

                locus_key = f'{read_chrom}:{locus_start}-{locus_end}'
                loci_keys.append(locus_key)
                read_loci_variations[locus_key] = {'halen': locus_len, 'alen': locus_len, 'rlen': locus_len, 'seq': []}

                if locus_key not in global_loci_variations:
                    global_loci_variations[locus_key] = {'rlen': locus_len, 'reads': [], 'read_allele': {}, 'read_sequence': {}}
                    global_loci_info[locus_key] = row
                    global_loci_ends.append(locus_end)
                    global_loci_keys.append(locus_key)

        
                read_index += 1
                read_quality = read.query_qualities
                cigar_tuples = read.cigartuples
                read_sequence = read.query_sequence

                cigar_one = cigar_tuples[0]

                global_read_ends.append(read_end)
                global_read_indices.append(read_index)
                global_read_variations[read_index] = {'s': read_start, 'e': read_end, 'snps': set(), 'dels': set()}
    
                # for each_coords in loci_coords:
                #     update_homopolymer_coords(ref.fetch(Chrom, each_coords[0]-100, each_coords[1]+100), each_coords[0]-100, homopoly_positions)
    
                if read.has_tag('cs'):
                    cs_tag = read.get_tag('cs')
                    parse_cstag(read_index, cs_tag, read_start, loci_keys, loci_coords, read_loci_variations, homopoly_positions, global_read_variations, global_snp_positions, read_sequence, read_quality, cigar_one, sorted_global_snp_list, left_flank_list, right_flank_list, male)
                else :
                    parse_cigar_tag(read_index, cigar_tuples, read_start, loci_keys, loci_coords, read_loci_variations,
                                homopoly_positions, global_read_variations, global_snp_positions, read_sequence, read, ref, read_quality, sorted_global_snp_list, left_flank_list, right_flank_list, male)
    
                for locus_key in read_loci_variations:
                    global_loci_variations[locus_key]['reads'].append(read_index)
                    global_loci_variations[locus_key]['read_allele'][read_index] = [read_loci_variations[locus_key]['halen'], read_loci_variations[locus_key]['alen']]
                    global_loci_variations[locus_key]['read_sequence'][read_index] = read_loci_variations[locus_key]['seq']
            
            while len(global_loci_ends) > 0:

                genotyped_loci_count += locus_processor(global_loci_keys, global_loci_ends, global_loci_variations, global_read_variations, global_snp_positions, prev_reads, sorted_global_snp_list, maxR, minR, ref, Chrom, global_loci_info, out, snpQ, snpC, snpD, snpR, phasingR, tbx2, flank, sorted_global_ins_rpos_set, log_bool, logger, male)

            while len(global_read_ends) > 0:

                popped = global_read_ends.pop(0)
                rindex = global_read_indices.pop(0)
                if rindex in global_read_variations:
                    for pos in global_read_variations[rindex]['snps']:
                        if pos in global_snp_positions:
                            global_snp_positions[pos]['cov'] -= 1
                    del_snps = []
                    for pos in global_snp_positions: 
                        if global_snp_positions[pos]['cov'] == 0:
                            del_snps.append(pos)
                            sorted_global_snp_list.remove(pos)
                    for snp in del_snps:
                        del global_snp_positions[snp]
                    del global_read_variations[rindex]
                    del del_snps

                    if rindex in prev_reads: prev_reads.remove(rindex)
                del_ins_pos_idx = 0
                list_rpos = sorted(sorted_global_ins_rpos_set)
                for i in list_rpos:
                    del_ins_pos_idx+=1
                    if i > popped: break
                del list_rpos[:del_ins_pos_idx]
                sorted_global_ins_rpos_set = set(list_rpos)

        print(f'\nTotal genotyped loci = {genotyped_loci_count} out of {total_loci} in {Chrom} {Start[0]}-{End[1]}\n', file=sys.stderr)

    bam.close()
    ref.close()
    tbx.close()
    tbx2.close()
    out.close()