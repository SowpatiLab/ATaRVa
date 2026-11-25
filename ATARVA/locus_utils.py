from bisect import bisect_left, bisect_right
import string

from ATARVA.realignment import stripSW, Inputs
from ATARVA.process_softclips import cigar_stats
from ATARVA.decompose         import motif_decomposition

# base64 encodings
BASE64_CHARS = string.ascii_uppercase + string.ascii_lowercase + string.digits + '+//'

def break_locuskey(locus_key):
    """
    break the locus key into chrom, start and end coordinates

    :param locus_key: string in the format chrom:start-end
    :return: chrom, start and end coordinates
    """
    chrom = locus_key[:locus_key.index(':')]
    locus_key = locus_key[locus_key.index(':') + 1:]
    start, end, motif = locus_key.split('-')
    start = int(start)
    end   = int(end)
    return chrom, start, end, motif


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

    snp_dist  = cooper.args.snp_dist
    snp_start = locus_start - snp_dist
    snp_end   = locus_end + snp_dist

    # SNPs inside the locus window: computed once, not once per read
    sorted_snps = cooper.cooper_sorted_snps
    window      = sorted_snps[bisect_left(sorted_snps, snp_start):
                              bisect_right(sorted_snps, snp_end)]
    if not window:
        return

    snp_data = cooper.cooper_snp_data
    reads    = cooper.cooper_read_data

    for read_index in new_reads:
        read = reads[read_index]

        # only SNPs overlapping this read: O(log n) instead of a linear scan
        candidates = [pos for pos in window[bisect_left(window, read.start):
                            bisect_right(window, read.end)] if pos not in read.ref_snps]
        if not candidates:
            continue

        alt_snps = read.snps
        dels = read.dels

        if dels:
            for pos in candidates:
                if pos not in alt_snps and bisect_right(dels, pos) % 2 == 0:
                    d = snp_data[pos]
                    d.ref.add(read_index)
                    d.cov += 1
                    read.ref_snps.add(pos)
        else:
            # no deletions: skip the bisect entirely
            for pos in candidates:
                if pos not in alt_snps:
                    d = snp_data[pos]
                    d.ref.add(read_index)
                    d.cov += 1
                    read.ref_snps.add(pos)


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
    # read_indices, read_haplotags = zip(
    #     *[(i, ht) for i, ht in zip(read_indices, read_haplotags) if i in top_reads]
    # ) if top_reads else ([], [])

    locus_data.reads          = list(top_reads)
    locus_data.depth          = len(read_indices)


def _check_insertion_fit(cooper, ins_neighbors, insert, locus):
    """
    check if the insertion fits in the reference context

    :param cooper: ATaRVa object
    :param ins_neighbors: list of neighboring loci coordinates for the insertion
    :param insert: inserted sequence
    :return: True if the insertion fits, False otherwise
    """

    max_match_locus = None
    max_match_score = 0
    for neighbor in ins_neighbors:
        n_chrom, n_start, n_end, n_motif = neighbor
        ref_seq = cooper.ref.fetch(n_chrom, n_start, n_end)
        target  = ""
        while len(target) < len(insert): target += ref_seq
        alignment_score, _align_score, target_begin, target_end, query_begin, query_end, sCigar = stripSW(Inputs(target, insert), False)
        matches, _ = cigar_stats(sCigar, False)

        if query_end - query_begin <= round(0.2 * len(insert)) or \
           matches / (query_end - query_begin) < 0.7:
            continue
        if alignment_score > max_match_score:
            max_match_score = alignment_score
            max_match_locus = neighbor

    if max_match_locus is None:
        return False

    if max_match_locus == (locus.chrom, locus.start, locus.end, locus.motif):
        return True
    else:
        return False


def process_flank_insertions(cooper, flank_insertions, ref_allele, ref_length, query, locus, locus_neighbors, insert_positions, is_left):
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
    included_insertions = set()

    for fid, (ins_rpos, ins_qs, ins_qe) in enumerate(flank_insertions):
        ins_len = ins_qe - ins_qs
        if ins_len < locus.motif_length and ins_len < 10: continue
        if ins_rpos in insert_positions:                  continue  # if the insertion position is already recorded for another locus, skip
        if inrepeat_ins(locus_neighbors, ins_rpos, insert_positions): continue

        insert = query[ins_qs:ins_qe]
        if insert == "": continue

        ins_neighbors = [
            (r.split('\t')[0], int(r.split('\t')[1]), int(r.split('\t')[2]), r.split('\t')[3])
            for r in cooper.tbx.fetch(locus.chrom, ins_rpos - cooper.args.flank, ins_rpos + cooper.args.flank)
        ]

        check_insfit = _check_insertion_fit(cooper, ins_neighbors, insert, locus)
        if check_insfit:
            ILR  += 1
            adj_pos = ins_qs if is_left else ins_qe
            included_insertions.add(ins_rpos)

    return adj_pos, included_insertions, ILR, PI, CI


def assign_hap_category(locus_data):
    """
    assign category for the locus based on haplotype status and allele distribution

    :param locus_data: LocusVariation object containing locus-specific data
    :return: None (updates the locus_data.hap_category in place)
    """

    if locus_data.is_phased and (list(locus_data.read_haplotags.values()).count(None) / locus_data.depth) <= 0.15:
        locus_data.hap_category = 3 # phased based on haplotag
        return

    # read_aseqs = [aseq[0] for aseq in locus_data.read_aseqs.values()]
    # locus_data.phase_mode = None
    # seq_counter = {}
    # for aseq in read_aseqs:
    #     try: seq_counter[aseq] += 1
    #     except KeyError: seq_counter[aseq] = 1
    # filtered_seqs = {seq: count for seq, count in seq_counter.items() if count > 1}
    # max_freq = max(filtered_seqs.values()) if filtered_seqs else 0
    # if len(filtered_seqs) == 1 or max_freq / locus_data.depth >= 0.75:
    # if max_freq / locus_data.depth >= 0.9:
    # if len(seq_counter) == 1:
    #     locus_data.hap_category    = 1 # homozygous
    #     return

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
    locus_data.neighbors.remove((locus.start, locus.end, locus.motif))

    read_indices      = locus_data.reads

    pending_insertions  = set()
    ILR = PI = CI       = 0

    # --- coverage check ---
    if locus_data.depth < cooper.args.min_reads:
        cooper.prev_reads    = set(read_indices)
        cooper.prev_reads.clear()
        locus_data.skip_code = 0
        return

    if not cooper.args.amplicon and locus_data.depth > cooper.args.max_reads:
        subset_reads(cooper, locus_data)

    current_reads = set(read_indices)
    new_reads     = current_reads - cooper.prev_reads

    # --- methylation thresholds ---
    upper_bound = cooper.args.meth_prob
    lower_bound = 1 - upper_bound

    read_level_bool = cooper.args.read_level

    # --- per read processing ---
    for read_index in read_indices:
        read = cooper.cooper_read_data[read_index]
        query, relative_qrange, left_ins, right_ins, left_dels, right_dels, fqs, fqe = locus_data.read_aseqs[read_index]
        adj_qs, adj_qe = relative_qrange

        left_ins.sort(key=lambda x: x[0])
        right_ins.sort(key=lambda x: x[0], reverse=True)

        # process flanks
        new_qs, pend_l, ilr, pi, ci    = process_flank_insertions(cooper, left_ins,  ref_allele, ref_length, query, locus, locus_data.neighbors,
                                                                  cooper.cooper_insert_positions[read_index], is_left=True)
        new_qe, pend_r, ilr2, pi2, ci2 = process_flank_insertions(cooper, right_ins, ref_allele, ref_length, query, locus, locus_data.neighbors,
                                                                  cooper.cooper_insert_positions[read_index], is_left=False)

        if cooper.args.strict:
            if len(pend_l) > 0 or len(pend_r) > 0:
                locus_data.flank_var = True
        else:
            if new_qs is not None: adj_qs = new_qs
            if new_qe is not None: adj_qe = new_qe

        ILR += ilr + ilr2
        PI  += pi  + pi2
        CI  += ci  + ci2
        pending_insertions |= pend_l | pend_r

        current_allele                         = query[adj_qs:adj_qe]
        alen                                   = adj_qe - adj_qs
        locus_data.read_aseqs[read_index][0]   = current_allele
        locus_data.read_alens[read_index][1]   = alen
  

        # methylation
        subseq_len     = fqe - fqs
        adj_fqs        = fqs + adj_qs
        adj_fqe        = fqe - (subseq_len - adj_qe)

        current_meth_info                           = read.process_methylation_info(adj_fqs, adj_fqe, lower_bound, upper_bound)
        locus_data.read_methylation[read_index]     = current_meth_info
        cooper.cooper_insert_positions[read_index] |= pending_insertions

        if read_level_bool:
            dseq,_ = motif_decomposition(current_allele, locus.motif_length)
            if current_meth_info is not None:
                _, _, meth_scores, _ = current_meth_info
                avg_read_meth = round(sum([(i+1) for i in meth_scores]) / len(meth_scores) / 256, 2)
                base64_scaled = ''.join(BASE64_CHARS[round((i+1)/4)] for i in meth_scores)
            else:
                avg_read_meth = '.'
                base64_scaled = '.'    
            cooper.readlevel_handle.write(f"{locus_key}\t{read_index}\t{current_allele}\t{alen}\t{dseq}\t{avg_read_meth}\t{base64_scaled}\n")

    if cooper.args.debug_mode:
        cooper.logger.debug(f"{locus_key};Larger_ins={ILR};Partial_ins={PI};Complete_ins={CI}")

    count_alleles(cooper, locus_key) # updates frequency of allele lengths in cooper.cooper_loci_data[locus_key]

    record_ref_snps(cooper, new_reads, locus.start, locus.end)

    if cooper.args.haplotag and (len(set(locus_data.read_haplotag_ps.values())) == 1 and None not in set(locus_data.read_haplotag_ps.values())):
        hap1, hap2 = [], []
        for read_idx in read_indices:
            tag = locus_data.read_haplotags[read_idx]
            if tag == 1: hap1.append(read_idx)
            if tag == 2: hap2.append(read_idx)
        locus_data.hap_read_sets        = (hap1, hap2)
        locus_data.hap_alen_sets = ([locus_data.read_alens[ridx][0] for ridx in hap1], [locus_data.read_alens[ridx][0] for ridx in hap2])
        # haplotyped successfully if there are reads in both haplogroups
        locus_data.is_phased  = all(locus_data.hap_read_sets)
        locus_data.phase_mode = 'haplotag' if locus_data.is_phased else None

    assign_hap_category(locus_data)

    cooper.prev_reads.clear()
    # cooper.prev_reads.update(current_reads)
    locus_data.skip_code = 10
