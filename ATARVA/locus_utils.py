from ATARVA.realignment_utils import *
import sys, bisect
# import statistics

def count_alleles(locus_key, read_indices, global_loci_variations, allele_counter, hallele_counter):
    """
    Counts the read distribution for each allele length
    """
    for rindex in read_indices:
        halen, alen = global_loci_variations[locus_key]['read_allele'][rindex]

        try: allele_counter[alen] += 1
        except KeyError: allele_counter[alen] = 1

        try: hallele_counter[halen] += 1
        except KeyError: hallele_counter[halen] = 1


def record_snps(read_indices, old_reads, new_reads, global_read_variations, global_snp_positions, sorted_global_snp_list, locus_start, locus_end, snp_dist, prev_locus_end):

    read_indices = set(read_indices)
    snp_start = locus_start - snp_dist
    snp_end = locus_end + snp_dist
    prev_locus_end += snp_dist

    for rindex in read_indices:
        # if rindex not in new_reads: continue

        read_variation = global_read_variations[rindex]
        rstart = read_variation['s']
        rend   = read_variation['e']
        snps   = read_variation['snps']
        dels   = read_variation['dels']

            
        for pos in sorted_global_snp_list:
            if pos < prev_locus_end: continue
            if not (snp_start <= pos <= snp_end): continue
            if pos < rstart: continue
            if pos > rend: break
            
            if (pos not in snps) and (bisect.bisect(dels, pos) % 2 == 0):
                if 'r' in global_snp_positions[pos]: global_snp_positions[pos]['r'].add(rindex)
                else: global_snp_positions[pos]['r'] = {rindex}
                global_snp_positions[pos]['cov'] += 1


def inrepeat_ins(near_by_loci, ins_rpos, sorted_global_ins_rpos_set):
    for locus in near_by_loci:
        if locus[0] <= ins_rpos <= locus[1]:
            sorted_global_ins_rpos_set.add(ins_rpos)
            return 1
    return 0


def process_locus(locus_key, global_loci_variations, global_read_variations, global_snp_positions, prev_reads, sorted_global_snp_list, maxR, minR, global_loci_info, near_by_loci, sorted_global_ins_rpos_set, Chrom, locus_start, locus_end, ref, log_bool, logger, snp_dist, prev_locus_end, hp_code):


    ref_seq = ref.fetch(Chrom, locus_start, locus_end)
    ref_len = len(ref_seq)
    locus_tuple = (locus_start, locus_end)
    near_by_loci.remove(locus_tuple)
    
    read_tag = global_loci_variations[locus_key]['read_tag']
    category, haplotypes = [None, None]
    homozygous_allele = 0
    read_indices = global_loci_variations[locus_key]['reads']   # the read indices which cover the locus
    reads_of_homozygous = read_indices.copy()
    total_reads = len(read_indices)                             # total number of reads
    max_limit=0

    period = int(float(global_loci_info[locus_key][4]))
    new_ins_rpos_current_loci = set()

    # remove if the locus has poor coverage
    if total_reads < minR:
        # coverage of the locus is low
        prev_reads = set(read_indices)
        return [prev_reads, category, homozygous_allele, reads_of_homozygous, {}, 0, max_limit, haplotypes]
    elif total_reads > maxR:
        # coverage of the locus is high
        read_indices = read_indices[:maxR]
        read_tag = read_tag[:maxR]
        max_limit=1
    
    current_reads = set(read_indices)
    old_reads = prev_reads - current_reads
    new_reads = current_reads - prev_reads

    locus_read_allele = global_loci_variations[locus_key]['read_allele'] # extracting allele info from global_loci_variation
    locus_read_seq = global_loci_variations[locus_key]['read_sequence']
    ILR=0;PI=0;CI=0;
    for each_read in read_indices:

        query,rep_range,ins_left,ins_right, left_rpos, right_rpos = locus_read_seq[each_read] # fetching repeat seq with flanks, correct start end position and insertion coordinates

        new_start,new_end = rep_range # new coordinates same as correct corrdinates
        if new_end-new_start != locus_read_allele[each_read][0]:
            print(f'Calculated allele length is not same at locus {locus_key} where alen is {locus_read_allele[each_read][0]} and ranges is {rep_range} and end-start = {new_end-new_start}')
            sys.exit()
        
        
        sorted_left = sorted(ins_left,key = lambda x: x[0]) # sorting the coordinates so if the 1st insertion is itself a repeat, fetch seq from that position; no need for checking the successive ins (only for all left ins)
        sorted_right = sorted(ins_right,key = lambda x: x[0], reverse=True) # for right ins, there are no breaks
        sorted_left_rpos = sorted(left_rpos)
        sorted_right_rpos = sorted(right_rpos, reverse=True)


        for lid,each_tuple in enumerate(sorted_left):# checking the insertion on left, whether its a repeats or not
            ins_len = each_tuple[1]-each_tuple[0]
            
            if ins_len < period:
                if ins_len>=10: pass
                else: continue
            ins_rpos = sorted_left_rpos[lid]
            if ins_rpos in sorted_global_ins_rpos_set: continue
            elif inrepeat_ins(near_by_loci, ins_rpos, sorted_global_ins_rpos_set): continue
            else:
                test_query = query[each_tuple[0]: each_tuple[1]]
                align, pos = stripSW(Inputs(ref_seq, test_query))
                que_len = len(test_query)
                align_len = len(align)

                if align_len<=round(0.2*min([que_len,ref_len])):
                    continue
                elif (align_len >= round(0.75*ref_len)) and (align.count('|') >= round(0.75*align_len)): # when insertion is larger then the ref seq
                    ILR+=1
                    new_start = each_tuple[0] #+ pos[0]
                    for ins in sorted_left_rpos[lid:]:
                        new_ins_rpos_current_loci.add(ins)
                    break
                elif (align.count('|') >= round(0.75*align_len)) and (pos[1]>=round(0.7*que_len)) and (align_len>=round(0.45*que_len)):
                    if align_len<=0.5*que_len: PI+=1
                    else: CI+=1
                    # print('either CI or PI')
                    new_start = each_tuple[0] + pos[0]
                    for ins in sorted_left_rpos[lid:]:
                        new_ins_rpos_current_loci.add(ins)
                    break

        for rid,each_tuple in enumerate(sorted_right):
            ins_len = each_tuple[1]-each_tuple[0]
            if ins_len < period:
                if ins_len>=10: pass
                else: continue
            ins_rpos = sorted_right_rpos[rid]
            if ins_rpos in sorted_global_ins_rpos_set: continue
            elif inrepeat_ins(near_by_loci, ins_rpos, sorted_global_ins_rpos_set): continue
            else:
                test_query = query[each_tuple[0]: each_tuple[1]]
                align, pos = stripSW(Inputs(ref_seq, test_query))
                que_len = len(test_query)
                align_len = len(align)
                if align_len<=round(0.2*min([que_len,ref_len])):
                    continue
                elif (align_len >= round(0.75*ref_len)) and (align.count('|') >= round(0.75*align_len)): # when insertion is larger then the ref seq
                    ILR+=1
                    new_end = each_tuple[1] #each_tuple[0] + pos[1]
                    for ins in sorted_right_rpos[rid:]:
                        new_ins_rpos_current_loci.add(ins)
                    break
                elif (align.count('|') >= round(0.75*align_len)) and (pos[0]<=round(0.3*que_len)) and (align_len>=round(0.45*que_len)):
                    if align_len<=0.5*que_len: PI+=1
                    else: CI+=1
                    new_end = each_tuple[0] + pos[1]
                    for ins in sorted_right_rpos[rid:]:
                        new_ins_rpos_current_loci.add(ins)
                    break
                    


        locus_read_seq[each_read][0] = query[new_start:new_end] # over-writing the query seq with modified seq with/without ins
        locus_read_allele[each_read][0] = new_end-new_start # over-writing the allele length after modification

    if log_bool: logger.debug(f"{locus_key};Larger_ins={ILR};Partial_ins={PI};Complete_ins={CI}")
    sorted_global_ins_rpos_set |= new_ins_rpos_current_loci
    # recording the counts of each allele length across all reads
    allele_counter = {};  hallele_counter = {}
    count_alleles(locus_key, read_indices, global_loci_variations, allele_counter, hallele_counter)

    hap_status = False
    if hp_code:
        haplotypes = ([read_indices[i] for i in [idx for idx,i in enumerate(read_tag) if i == 1]], [read_indices[i] for i in [idx for idx,i in enumerate(read_tag) if i == 2]])
        hap_status = all([len(hap)>0 for hap in haplotypes])
    
    if hap_status & ((read_tag.count(None)/total_reads) <= 0.15): # processing haplotagged reads to write into vcf_heterozygous
        category = 3 # phased
    
    elif len(hallele_counter) == 1:
        category = 1 # homozygous
        homozygous_allele = list(hallele_counter.keys())[0]
    
    else:
        filtered_alleles = list(filter(lambda x: hallele_counter[x] > 1, hallele_counter.keys()))
        if len(filtered_alleles) == 1 and hallele_counter[filtered_alleles[0]]/total_reads >= 0.75:
            category = 1 # homozygous
            homozygous_allele = filtered_alleles[0]
            reads_of_homozygous = [rindex for rindex in global_loci_variations[locus_key]['read_allele'] if homozygous_allele == global_loci_variations[locus_key]['read_allele'][rindex][0]]
        else:
            category = 2 # ambiguous
            
    record_snps(read_indices, old_reads, new_reads, global_read_variations, global_snp_positions, sorted_global_snp_list, locus_start, locus_end, snp_dist, prev_locus_end)
    
    prev_reads = current_reads.copy()
    return [prev_reads, category, homozygous_allele, reads_of_homozygous, hallele_counter, 10, max_limit, haplotypes]