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


def count_alleles(cooper, locus_key):
    """
    frequency of each allele length across all reads for a locus

    :param cooper: Cooper object
    :param locus_key: string in the format chrom:start-end
    :return: None (updates the cooper.cooper_loci_data[locus_key].alen_frequency and cooper.cooper_loci_data[locus_key].halen_frequency in place)
    """

    locus_data = cooper.cooper_loci_data[locus_key]
    read_indices = locus_data.reads
    for read_index in read_indices:
        halen, alen = locus_data.read_alens[read_index]
        locus_data.allele_lengths.append(halen)

        try: locus_data.alen_frequency[alen] += 1
        except KeyError: locus_data.alen_frequency[alen] = 1

        try: locus_data.halen_frequency[halen] += 1
        except KeyError: locus_data.halen_frequency[halen] = 1


def record_ref_snps(cooper, new_reads, locus_start, locus_end):
    """
    update ref allele and coverage for all reads with no SNPs

    :param cooper: ATaRVa object
    :param read_indices: list of read indices covering the locus
    :param new_reads: set of read indices which are new to the locus
    :param locus_start: start coordinate of the locus
    :param locus_end: end coordinate of the locus
    :param snp_dist: distance from the locus to consider for recording SNPs
    :return: None (updates the cooper.cooper_snp_data in place)
    """

    snp_dist = cooper.args.snp_dist
    snp_start = locus_start - snp_dist
    snp_end   = locus_end + snp_dist

    for read_index in new_reads:

        read_data = cooper.cooper_read_data[read_index]

        for pos in cooper.cooper_sorted_snps:
            # if snp not in range skip
            if not (snp_start <= pos <= snp_end): continue
            
            if pos < read_data.start: continue
            if pos > read_data.end: break
            
            if (pos not in read_data.snps) and (bisect.bisect(read_data.dels, pos) % 2 == 0):
                cooper.cooper_snp_data[pos].ref.add(read_index)
                cooper.cooper_snp_data[pos].cov += 1


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


def subset_reads(cooper, locus_data):
    """
    subsets reads for loci with high coverage based on read quality.

    :param cooper:         cooper object containing read data and args
    :param read_indices:   list of read indices covering the locus
    :param read_haplotags: list of haplotype tags for the reads
    :return:               subsetted read indices and corresponding haplotype tags
    """
    read_indices   = locus_data.reads
    read_haplotags = locus_data.read_haplotags

    locus_data.raw_reads     = read_indices
    locus_data.raw_haplotags = read_haplotags
    locus_data.raw_depth     = len(read_indices)

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

    locus_data.reads          = list(read_indices)
    locus_data.read_haplotags = list(read_haplotags)
    locus_data.depth          = len(read_indices)


def process_flank_insertions(flank_insertions, ref_allele, ref_length, query, locus, locus_neighbors, insert_positions, is_left):
    """
    process insertions in flanks and adjusts locus boundaries if needed

    :param flank_insertions: list of tuples with insertion reference position and query start and end positions
    :param ref_allele: reference allele sequence for the locus
    :param ref_length: reference allele length for the locus
    :param query: read sequence for the locus
    :param locus: locus object containing locus information
    :param locus_neighbors: list of neighboring loci coordinates
    :param insert_positions: set of reference positions of insertions across all loci
    :param is_left: boolean value indicating if the flank is left or right
    :return: adjusted query start or end position, set of pending insertions, counts of larger, partial and complete insertions
    """
    
    ILR = PI = CI = 0
    pending = set()
    adj_pos = None

    ref_75 = round(0.75 * ref_length)

    for fid, (ins_rpos, ins_qs, ins_qe) in enumerate(flank_insertions):
        ins_len = ins_qe - ins_qs
        if ins_len < locus.motif_length and ins_len < 10: continue
        if ins_rpos in insert_positions:                  continue
        if inrepeat_ins(locus_neighbors, ins_rpos, insert_positions): continue

        insert          = query[ins_qs:ins_qe]
        alignment, coords = stripSW(Inputs(ref_allele, insert), True)
        align_len       = len(alignment)
        matches         = alignment.count('|')
        min_len         = min(ins_len, ref_length)

        if align_len <= round(0.2 * min_len): continue

        if align_len >= ref_75 and matches >= round(0.75 * align_len):
            ILR  += 1
            adj_pos = ins_qs if is_left else ins_qe
            pending.update(ins[0] for ins in flank_insertions[fid:])
            break

        elif (matches   >= round(0.75 * align_len) and
              align_len >= round(0.45 * ins_len)):
            if is_left and coords[1] >= round(0.7 * ins_len):
                adj_pos = ins_qs + coords[0]
            elif not is_left and coords[0] <= round(0.3 * ins_len):
                adj_pos = ins_qe + coords[1]
            else:
                continue
            PI += 1 if align_len <= 0.5 * ins_len else 0
            CI += 1 if align_len >  0.5 * ins_len else 0
            pending.update(ins[0] for ins in flank_insertions[fid:])
            break

    return adj_pos, pending, ILR, PI, CI


def assign_hap_category(locus_data):
    """
    assign category for the locus based on haplotype status and allele distribution

    :param locus_data: LocusVariation object containing locus-specific data
    :param hap_status: boolean value indicating if haplotagging is available and informative
    :param read_haplotags: list of haplotype tags for the reads
    :param total_reads: total number of reads covering the locus
    :param hallele_counter: dict with haplotype-wise allele length counts
    :return: category (1 for homozygous, 2 for ambiguous, 3 for phased), homozygous allele length if applicable
    """

    if locus_data.hap_status and (locus_data.read_haplotags.count(None) / locus_data.depth) <= 0.15:
        locus_data.hap_category = 3 # phased
        return

    locus_data.phase_mode = None
    if len(locus_data.halen_frequency) == 1:
        locus_data.hap_category = 1 # homozygous
        locus_data.homozygous_alen = next(iter(locus_data.halen_frequency))   # homozygous
        return

    filtered = [a for a, c in locus_data.halen_frequency.items() if c > 1]
    if len(filtered) == 1 and locus_data.halen_frequency[filtered[0]] / locus_data.depth >= 0.75:
        locus_data.hap_category = 1 # homozygous
        locus_data.homozygous_alen = filtered[0]
        return                   # homozygous

    locus_data.hap_category = 2 # ambiguous


def process_locus(cooper, locus_key):
    """
    process the reads for a locus

    :param cooper:        cooper object containing read data and args
    :param locus_key:     string in the format chrom:start-end
    :return: None (updates the cooper.cooper_loci_data[locus_key] in place)
    """

    locus      = cooper.cooper_loci_info[locus_key]
    locus_data = cooper.cooper_loci_data[locus_key]

    ref_allele = cooper.ref.fetch(cooper.chrom, locus.start, locus.end)
    ref_length = locus.length
    locus_data.neighbors.remove((locus.start, locus.end))

    read_indices   = locus_data.reads
    read_haplotags = locus_data.read_haplotags
    locus_data.depth    = len(read_indices)

    pending_insertions     = set()
    ILR = PI = CI          = 0

    # --- coverage check ---
    if locus_data.depth < cooper.args.min_reads:
        cooper.prev_reads = set(read_indices)
        locus_data.skipped_tag = 0
        return

    if locus_data.depth > cooper.args.max_reads:
        subset_reads(cooper, locus_data)

    current_reads = set(read_indices)
    new_reads     = current_reads - cooper.prev_reads

    upper_bound = cooper.args.meth_prob
    lower_bound = 1 - upper_bound

    # --- per read processing ---
    for read_index in read_indices:
        read_data = cooper.cooper_read_data[read_index]
        query, relative_qrange, left_ins, right_ins, fqs, fqe = locus_data.read_seqs[read_index]
        adj_qs, adj_qe = relative_qrange

        left_ins.sort( key=lambda x: x[0])
        right_ins.sort(key=lambda x: x[0], reverse=True)

        # process flanks
        new_qs, pend_l, ilr, pi, ci    = process_flank_insertions(left_ins,  ref_allele, ref_length, query, locus,
                                                                  locus_data.neighbors, cooper.cooper_insert_positions, is_left=True)
        new_qe, pend_r, ilr2, pi2, ci2 = process_flank_insertions(right_ins, ref_allele, ref_length, query, locus,
                                                                  locus_data.neighbors, cooper.cooper_insert_positions, is_left=False)

        if new_qs is not None: adj_qs = new_qs
        if new_qe is not None: adj_qe = new_qe

        ILR += ilr + ilr2
        PI  += pi  + pi2
        CI  += ci  + ci2
        pending_insertions |= pend_l | pend_r

        locus_data.read_seqs[read_index][0]   = query[adj_qs:adj_qe]
        locus_data.read_alens[read_index][0]  = adj_qe - adj_qs

        # methylation
        subseq_len     = fqe - fqs
        adj_fqs        = fqs + adj_qs
        adj_fqe        = fqe - (subseq_len - adj_qe)
        read_mod_bases = read_data.methylation

        meth_pos = []; meth_encode = []; meth_count = meth_qual = 0
        for pos, raw_prob in read_mod_bases:
            if not (adj_fqs <= pos <= adj_fqe): continue
            prob = raw_prob / 255
            meth_pos.append(pos - adj_fqs)
            if lower_bound < prob < upper_bound:
                meth_encode.append(-1)
            else:
                meth_count += 1
                is_meth     = prob >= upper_bound
                meth_encode.append(int(is_meth))
                meth_qual  += is_meth

        locus_data.read_methylation[read_index] = (round(meth_qual / meth_count, 2), meth_pos, meth_encode) if meth_count else None

    if cooper.args.debug_mode:
        cooper.logger.debug(f"{locus_key};Larger_ins={ILR};Partial_ins={PI};Complete_ins={CI}")

    cooper.cooper_insert_positions |= pending_insertions
    count_alleles(cooper, locus_key) # updates frequency of allele lengths in cooper.cooper_loci_data[locus_key]

    if cooper.args.amplicon:
        locus.hap_category = 2  # ambigous for amplicons
    else:
        record_ref_snps(cooper, new_reads, locus.start, locus.end)

        if cooper.args.haplotag:
            hap1, hap2 = [], []
            for read_idx, tag in zip(read_indices, read_haplotags):
                if tag == 1: hap1.append(read_idx)
                if tag == 2: hap2.append(read_idx)
            locus_data.haplotypes = (hap1, hap2)
            # haplotyped successfully if there are reads in both haplogroups
            locus_data.hap_status = all(locus_data.haplotypes)
            if locus_data.hap_status:
                locus_data.phase_mode = 'haplotag'

        assign_hap_category(locus_data)

    cooper.prev_reads = current_reads
    locus_data.skipped_tag = 10
