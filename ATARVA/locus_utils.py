import bisect
from ATARVA.realignment import *

def break_locuskey(locus_key):
    """
    break the locus key into chrom, start and end coordinates

    :param locus_key: string in the format chrom:start-end
    :return: chrom, start and end coordinates
    """
    chrom = locus_key[:locus_key.index(':')]
    start = int(locus_key[(locus_key.index(':') + 1) : locus_key.index('-')])
    end   = int(locus_key[(locus_key.index('-') + 1) :])
    return chrom, start, end


def count_alleles(locus_key, read_indices, global_loci_variations, allele_counter, hallele_counter, alen_list):
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
        halen, alen = global_loci_variations[locus_key].read_alens[read_index]
        alen_list.append(halen)

        try: allele_counter[alen] += 1
        except KeyError: allele_counter[alen] = 1

        try: hallele_counter[halen] += 1
        except KeyError: hallele_counter[halen] = 1


def record_snps(read_indices, old_reads, new_reads, global_read_variations, global_snp_positions, sorted_global_snp_list, locus_start, locus_end, snp_dist):

    read_indices = set(read_indices)
    snp_start = locus_start - snp_dist
    snp_end = locus_end + snp_dist

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


def inrepeat_ins(near_by_loci, ins_rpos, sorted_ins_rpos):
    """
    if flank repeat falls in another repeat

    :param near_by_loci: list of nearby loci coordinates
    :param ins_rpos: reference position of the insertion
    :param sorted_ins_rpos: set of reference positions of insertions across all loci
    """
    for locus in near_by_loci:
        if locus[0] <= ins_rpos <= locus[1]:
            sorted_ins_rpos.add(ins_rpos)
            return 1
    return 0


def subset_reads(cooper, read_indices: list, read_haplotags: list) -> tuple:
    """
    subsets reads for loci with high coverage based on read quality.

    :param cooper:         cooper object containing read data and args
    :param read_indices:   list of read indices covering the locus
    :param read_haplotags: list of haplotype tags for the reads
    :return:               subsetted read indices and corresponding haplotype tags
    """
    # sort read indices by quality in one step — no intermediate dict
    sorted_reads = sorted(
        read_indices,
        key     = lambda i: cooper.cooper_read_data[i].mean_qual,
        reverse = True
    )

    # subset top max_reads
    top_reads = set(sorted_reads[:cooper.args.max_reads])

    # filter haplotags and indices in single pass
    read_indices, read_haplotags = zip(
        *[(i, ht) for i, ht in zip(read_indices, read_haplotags) if i in top_reads]
    ) if top_reads else ([], [])

    return sorted(read_indices), list(read_haplotags)


def process_locus(cooper, locus_key, locus_neighbors):
    """
    process the reads for a locus

    :param cooper:        cooper object containing read data and args
    :param locus_key:     string in the format chrom:start-end
    :param chrom:         chromosome name
    :param locus_start:   start coordinate of the locus
    :param locus_end:     end coordinate of the locus
    :param neighboring_loci:  list of neighboring loci coordinates
    :return: prev_reads, category, homozygous_allele, reads_of_homozygous, hallele_counter, skip_point, haplotypes, homozygous_alens
    """

    locus = cooper.cooper_loci_info[locus_key]

    ref_allele = cooper.ref.fetch(cooper.chrom, locus.start, locus.end)
    ref_length = locus.length
    locus_neighbors.remove((locus.start, locus.end))
    
    read_haplotags = cooper.cooper_loci_data[locus_key].read_haplotags
    category, haplotypes = [None, None]
    is_homozygous = False
    read_indices = cooper.cooper_loci_data[locus_key].reads   # the read indices which cover the locus
    total_reads = len(read_indices)                           # total number of reads

    pending_insertions = set()

    # remove if the locus has poor coverage
    if total_reads < cooper.args.min_reads:
        cooper.prev_reads = set(read_indices)
        return [category, is_homozygous, {}, 0, haplotypes, []]
    elif total_reads > cooper.args.max_reads:
        read_indices, read_haplotags = subset_reads(cooper, cooper.args.max_reads, read_indices, read_haplotags)
    
    current_reads = set(read_indices)
    old_reads = cooper.prev_reads - current_reads
    new_reads = current_reads - cooper.prev_reads

    locus_read_allele = cooper.cooper_loci_data[locus_key].read_alens # extracting allele info from global_loci_variation
    locus_read_seq = cooper.cooper_loci_data[locus_key].read_seqs

    ILR=0; PI=0; CI=0

    for read_index in read_indices:

        # relative_range is locus coordinates with respect to flank start
        query, relative_range, left_flank_insertions, right_flank_insertions, flank_query_start, flank_query_end = locus_read_seq[read_index]
    
        relative_start, relative_end = relative_range[0], relative_range[1]
        
        # sort based on reference position
        left_flank_insertions.sort(key=lambda x: x[0])
        right_flank_insertions.sort(key=lambda x: x[0], reverse=True)

        adj_locus_query_start, adj_locus_query_end = relative_range

        for lid, flank_insertion in enumerate(left_flank_insertions): # checking the insertion on left, whether its a repeats or not
            ins_rpos, ins_query_start, ins_query_end = flank_insertion
            insert_length = ins_query_end - ins_query_start

            if insert_length < locus.motif_length:
                if insert_length >= 10: pass # process if the the flank insert is larger than 10bp
                else: continue

            # insertion already recorded from a neighboring repeat
            if ins_rpos in cooper.cooper_insert_positions: continue

            # if the insert falls in the neighoring repeats; continue
            elif inrepeat_ins(locus_neighbors, ins_rpos, cooper.cooper_insert_positions): continue

            else:
                insert = query[ins_query_start: ins_query_end]
                alignment, align_coords = stripSW(Inputs(ref_allele, insert), True)
                alignment_length = len(alignment)
                matches = alignment.count('|')

                # not even 20% of the short sequence is covered in the alignment; skip the insert
                if alignment_length <= round(0.2*min([insert_length, ref_length])): continue

                # atleast 75% of ref is aligned and has 75% matches; adjust the coordinate of the locus
                elif (alignment_length >= round(0.75 * ref_length)) and \
                     (matches >= round(0.75 * alignment_length)):
                    ILR += 1
                    adj_locus_query_start = ins_query_start
                    for ins in left_flank_insertions[lid:]:
                        pending_insertions.add(ins)
                    break

                # atleast 75% of the aligned sequence has matches and the alignment starts within the first 30% of the insert; adjust the coordinate of the locus
                elif (matches >= round(0.75 * alignment_length)) and \
                     (align_coords[1] >= round(0.7*insert_length)) and \
                     (alignment_length >= round(0.45*insert_length)):
                    if alignment_length <= 0.5*insert_length: PI += 1
                    else: CI += 1
                    adj_locus_query_start = ins_query_start + align_coords[0]
                    for ins in left_flank_insertions[lid:]:
                        pending_insertions.add(ins[0])
                    break

        for rid, flank_insertion in enumerate(right_flank_insertions):
            ins_rpos, ins_query_start, ins_query_end = flank_insertion
            ins_len = ins_query_end - ins_query_start

            if ins_len < locus.motif_length:
                if ins_len >= 10: pass
                else: continue

            if ins_rpos in cooper.cooper_insert_positions: continue

            elif inrepeat_ins(locus_neighbors, ins_rpos, cooper.cooper_insert_positions): continue
            else:
                insert = query[ins_query_start: ins_query_end]
                alignment, align_coords = stripSW(Inputs(ref_allele, insert), True)
                alignment_length = len(alignment)
                matches = alignment.count('|')

                # not even 20% of the short sequence is covered in the alignment; skip the insert
                if alignment_length <= round(0.2 * min([ins_len, ref_length])): continue

                # atleast 75% of ref is aligned and has 75% matches; adjust the coordinate of the locus
                elif (alignment_length >= round(0.75 * ref_length)) and \
                     (matches >= round(0.75 * alignment_length)): # when insertion is larger then the ref seq
                    ILR += 1
                    adj_locus_query_end = ins_query_end
                    for ins in right_flank_insertions[rid:]:
                        pending_insertions.add(ins[0])
                    break
                elif (matches >= round(0.75 * alignment_length)) and \
                     (align_coords[0] <= round(0.3 * ins_len)) and \
                     (alignment_length >= round(0.45 * ins_len)):
                    if alignment_length <= 0.5 * ins_len: PI += 1
                    else: CI += 1
                    adj_locus_query_end = ins_query_end + align_coords[1]
                    for ins in right_flank_insertions[rid:]:
                        pending_insertions.add(ins[0])
                    break

        locus_read_seq[read_index][0] = query[adj_locus_query_start:adj_locus_query_end] # over-writing the query seq with modified seq with/without ins
        locus_read_allele[read_index][0] = adj_locus_query_end - adj_locus_query_start # over-writing the allele length after modification

        # Extracting the methylation probability for each read at the locus
        subseq_len = flank_query_end - flank_query_start
        adj_flank_query_start  = flank_query_start + adj_locus_query_start #( adjusted_query_start - rep_range[0] ) # adjusting the start position with respect to original read coordinates by adding the diff(after local-alignment)
        adj_flank_query_end    = flank_query_end - (subseq_len - adj_locus_query_end) # adjusting the start position with respect to original read coordinates by adding the diff(after local-alignment)

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
            if adj_flank_query_start <= each_pos[0] <= adj_flank_query_end:
                repeat_base_pos = each_pos[0] - adj_flank_query_start # position of the CG wrt the repeat start
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

    if cooper.args.debug_mode: cooper.logger.debug(f"{locus_key};Larger_ins={ILR};Partial_ins={PI};Complete_ins={CI}")
    cooper.cooper_insert_positions |= pending_insertions
    # recording the counts of each allele length across all reads
    allele_counter = {};  hallele_counter = {}; alen_list = []
    count_alleles(locus_key, read_indices, cooper.cooper_loci_data, allele_counter, hallele_counter, alen_list)

    if not cooper.args.amplicon:
        record_snps(read_indices, old_reads, new_reads, cooper, locus.start, locus.end)

        hap_status = False
        if cooper.args.haplotype:
            haplotypes = ([read_indices[i] for i in [idx for idx,i in enumerate(read_haplotags) if i == 1]], [read_indices[i] for i in [idx for idx,i in enumerate(read_haplotags) if i == 2]])
            hap_status = all([len(hap)>0 for hap in haplotypes])
        
        if hap_status & ((read_haplotags.count(None)/total_reads) <= 0.15): # processing haplotagged reads to write into vcf_heterozygous
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
    
    cooper.prev_reads = current_reads.copy()
    return [category, homozygous_allele, read_indices, hallele_counter, 10, haplotypes, alen_list]
