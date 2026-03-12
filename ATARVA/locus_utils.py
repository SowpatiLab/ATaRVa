import bisect
from ATARVA.realignment import *

def break_locuskey(locus_key):
    """
    Break the locus key into chrom, start and end coordinates

    :param locus_key: string in the format chrom:start-end

    :return: chrom, start and end coordinates
    """
    chrom = locus_key[:locus_key.index(':')]
    start = int(locus_key[(locus_key.index(':') + 1) : locus_key.index('-')])
    end   = int(locus_key[(locus_key.index('-') + 1) :])
    return chrom, start, end


def count_alleles(locus_key, read_indices, global_loci_variations, allele_counter, hallele_counter,alen_list):
    """
    Frequency of each allele length across all reads for a locus

    :param locus_key: string in the format chrom:start-end
    :param read_indices: list of read indices covering the locus
    :param global_loci_variations: dict containing locus-wise variations
    :param allele_counter: dict to store allele length counts
    :param hallele_counter: dict to store haplotype-wise allele length counts
    :param alen_list: list to store allele lengths for all reads
    
    :return: None (updates the allele_counter, hallele_counter and alen_list in place)
    """

    for read_index in read_indices:
        halen, alen = global_loci_variations[locus_key]['read_allele'][read_index]
        alen_list.append(halen)

        try: allele_counter[alen] += 1
        except KeyError: allele_counter[alen] = 1

        try: hallele_counter[halen] += 1
        except KeyError: hallele_counter[halen] = 1


def record_snps(read_indices, old_reads, new_reads, global_read_variations, global_snp_positions, sorted_global_snp_list, locus_start, locus_end, snp_dist, prev_locus_end):

    read_indices = set(read_indices)
    snp_start = locus_start - snp_dist
    snp_end = locus_end + snp_dist
    prev_locus_end += snp_dist

    for read_index in read_indices:
        if read_index not in new_reads: continue 

        read_variation = global_read_variations[read_index]
        rstart = read_variation['s']
        rend   = read_variation['e']
        snps   = read_variation['snps']
        dels   = read_variation['dels']

        for pos in sorted_global_snp_list:
            if not (snp_start <= pos <= snp_end): continue
            if pos < rstart: continue
            if pos > rend: break
            
            if (pos not in snps) and (bisect.bisect(dels, pos) % 2 == 0):
                if 'r' in global_snp_positions[pos]: global_snp_positions[pos]['r'].add(read_index)
                else: global_snp_positions[pos]['r'] = { read_index }
                global_snp_positions[pos]['cov'] += 1


def inrepeat_ins(near_by_loci, ins_rpos, sorted_global_ins_rpos_set):
    for locus in near_by_loci:
        if locus[0] <= ins_rpos <= locus[1]:
            sorted_global_ins_rpos_set.add(ins_rpos)
            return 1
    return 0


def subset_reads(global_read_variations, max_reads, read_indices, read_hgroup):
    """
    Subsets reads for loci with high coverage based on read quality
    
    :param global_read_variations: dict containing read-wise variations
    :param max_reads: maximum number of reads to subset
    :param read_indices: list of read indices covering the locus
    :param read_hgroup: list of haplotype tags for the reads

    :return: subsetted read indices and corresponding haplotype tags
    """
    read_qualities = {} # dict to store read quality
    for read_index in read_indices:
        read_qualities[read_index] = global_read_variations[read_index]['q']

    # sort reads based on quality and subset the top max_reads
    sorted_reads = [read_index for read_index,_ in sorted(read_qualities.items(), key = lambda x: x[1], reverse = True)]
    
    tmp_read_indices = set(sorted_reads[:max_reads])
    good_qual_read_idx = [idx for idx, i in enumerate(read_indices) if i in tmp_read_indices] # getting the index of the good_qual read-ids
    read_hgroup = [read_hgroup[i] for i in good_qual_read_idx] # extracting hp-tags of good-qual reads
    read_indices = sorted(tmp_read_indices)

    del tmp_read_indices, good_qual_read_idx

    return read_indices, read_hgroup


def process_locus(cooper, locus_key, chrom, locus_start, locus_end, near_by_loci):


    ref_allele = cooper.ref.fetch(chrom, locus_start, locus_end)
    ref_allele_length = len(ref_allele)
    locus_tuple = (locus_start, locus_end)
    near_by_loci.remove(locus_tuple)
    
    read_tag = cooper.cooper_loci_data[locus_key].read_hapgp
    category, haplotypes = [None, None]
    homozygous_allele = 0
    read_indices = cooper.cooper_loci_data[locus_key].reads   # the read indices which cover the locus
    reads_of_homozygous = read_indices.copy()
    total_reads = len(read_indices)                             # total number of reads
    max_limit = 0

    period = int(float(cooper.cooper_loci_info[locus_key][4]))
    new_ins_rpos_current_loci = set()

    # remove if the locus has poor coverage
    if total_reads < cooper.args.min_reads:
        # coverage of the locus is low
        prev_reads = set(read_indices)
        return [prev_reads, category, homozygous_allele, reads_of_homozygous, {}, 0, haplotypes, []]
    elif total_reads > cooper.args.max_reads:
        read_indices, read_tag = subset_reads(cooper, cooper.args.max_reads, read_indices, read_tag)
    
    current_reads = set(read_indices)
    old_reads = prev_reads - current_reads
    new_reads = current_reads - prev_reads

    locus_read_allele = cooper.cooper_loci_data[locus_key].read_alens # extracting allele info from global_loci_variation
    locus_read_seq = cooper.cooper_loci_data[locus_key].read_seqs

    ILR=0;PI=0;CI=0

    for read_index in read_indices:

        query, rep_range, ins_left, ins_right, left_rpos, right_rpos, read_repeat_start, read_repeat_end = locus_read_seq[read_index] # fetching repeat seq with flanks, correct start end position and insertion coordinates

        new_start, new_end = rep_range

        sorted_left  = sorted(ins_left, key = lambda x: x[0]) # sorting the coordinates so if the 1st insertion is itself a repeat, fetch seq from that position; no need for checking the successive ins (only for all left ins)
        sorted_right = sorted(ins_right, key = lambda x: x[0], reverse=True) # for right ins, there are no breaks
        sorted_left_rpos = sorted(left_rpos)
        sorted_right_rpos = sorted(right_rpos, reverse=True)

        for lid, each_tuple in enumerate(sorted_left): # checking the insertion on left, whether its a repeats or not
            ins_len = each_tuple[1]-each_tuple[0]

            if ins_len < period:
                if ins_len>=10: pass
                else: continue
            ins_rpos = sorted_left_rpos[lid]
            if ins_rpos in sorted_global_ins_rpos_set: continue
            elif inrepeat_ins(near_by_loci, ins_rpos, sorted_global_ins_rpos_set): continue
            else:
                test_query = query[each_tuple[0]: each_tuple[1]]
                align, pos = stripSW(Inputs(ref_allele, test_query), True)
                que_len = len(test_query)
                align_len = len(align)

                if align_len<=round(0.2*min([que_len, ref_allele_length])):
                    continue
                elif (align_len >= round(0.75 * ref_allele_length)) and (align.count('|') >= round(0.75 * align_len)): # when insertion is larger then the ref seq
                    ILR += 1
                    new_start = each_tuple[0] #+ pos[0]
                    for ins in sorted_left_rpos[lid:]:
                        new_ins_rpos_current_loci.add(ins)
                    break
                elif (align.count('|') >= round(0.75*align_len)) and (pos[1]>=round(0.7*que_len)) and (align_len>=round(0.45*que_len)):
                    if align_len<=0.5*que_len: PI+=1
                    else: CI+=1
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
                align, pos = stripSW(Inputs(ref_allele, test_query), True)
                que_len = len(test_query)
                align_len = len(align)
                if align_len<=round(0.2*min([que_len, ref_allele_length])):
                    continue
                elif (align_len >= round(0.75 * ref_allele_length)) and (align.count('|') >= round(0.75 * align_len)): # when insertion is larger then the ref seq
                    ILR += 1
                    new_end = each_tuple[1] #each_tuple[0] + pos[1]
                    for ins in sorted_right_rpos[rid:]:
                        new_ins_rpos_current_loci.add(ins)
                    break
                elif (align.count('|') >= round(0.75 * align_len)) and (pos[0] <= round(0.3 * que_len)) and (align_len >= round(0.45 * que_len)):
                    if align_len <= 0.5 * que_len: PI += 1
                    else: CI += 1
                    new_end = each_tuple[0] + pos[1]
                    for ins in sorted_right_rpos[rid:]:
                        new_ins_rpos_current_loci.add(ins)
                    break

        locus_read_seq[read_index][0] = query[new_start:new_end] # over-writing the query seq with modified seq with/without ins
        locus_read_allele[read_index][0] = new_end-new_start # over-writing the allele length after modification

        # Extracting the methylation probability for each read at the locus
        subseq_len = read_repeat_end - read_repeat_start
        adj_start = read_repeat_start + new_start #( new_start - rep_range[0] ) # adjusting the start position with respect to original read coordinates by adding the diff(after local-alignment)
        adj_end = read_repeat_end - (subseq_len - new_end) # adjusting the start position with respect to original read coordinates by adding the diff(after local-alignment)

        read_mod_bases  = cooper.cooper_read_data[read_index].methyl   # modified_bases data
        locus_read_meth = cooper.cooper_loci_data[locus_key].read_meth # read_wise methylation data at the locus

        # calculating average methylation probability at the locus for each read
        meth_count = 0
        meth_qual  = 0
        meth_pos   = []
        meth_encode = []
        upper_bound = cooper.args.meth_prob
        lower_bound = 1 - cooper.args.meth_prob
        for each_pos in read_mod_bases:
            if adj_start <= each_pos[0] <= adj_end:
                repeat_base_pos = each_pos[0] - adj_start # position of the CG wrt the repeat start
                current_prob = each_pos[1]/255
                if lower_bound < current_prob < upper_bound: # skip the MM if it is within the lower and upper bound; only store the extremies
                    meth_pos.append(repeat_base_pos)
                    meth_encode.append(-1) # to indicate skipped positions
                    continue

                meth_count += 1
                if current_prob >= cooper.args.meth_prob:
                    meth_pos.append(repeat_base_pos)
                    meth_encode.append(1) # encoding the probability as binary
                    meth_qual += 1
                else:
                    meth_pos.append(repeat_base_pos)
                    meth_encode.append(0) # encoding the probability as binary
            
        if meth_count > 0:
            avg_qual = meth_qual/meth_count
            locus_read_meth[read_index] = (round(avg_qual, 2), meth_encode, meth_pos) # storing meth level, position meth prob encoding and meth occurrence positions
        else:
            locus_read_meth[read_index] = None                

    if cooper.args.log_bool: cooper.logger.debug(f"{locus_key};Larger_ins={ILR};Partial_ins={PI};Complete_ins={CI}")
    sorted_global_ins_rpos_set |= new_ins_rpos_current_loci
    # recording the counts of each allele length across all reads
    allele_counter = {};  hallele_counter = {}; alen_list = []
    count_alleles(locus_key, read_indices, cooper.cooper_loci_data, allele_counter, hallele_counter, alen_list)

    if not cooper.args.amplicon:
        record_snps(read_indices, old_reads, new_reads, cooper, locus_start, locus_end, cooper.prev_locus_end)

        hap_status = False
        if cooper.args.haplotype:
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
                # reads_of_homozygous = [rindex for rindex in global_loci_variations[locus_key]['read_allele'] if homozygous_allele == global_loci_variations[locus_key]['read_allele'][rindex][0]]
            else:
                category = 2 # ambiguous
    else:
        category = 2 # ambiguous
    
    prev_reads = current_reads.copy()
    return [prev_reads, category, homozygous_allele, read_indices, hallele_counter, 10, haplotypes, alen_list]