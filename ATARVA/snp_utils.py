

def haplocluster_reads(cooper, locus_key):
    """
    cluster reads into haplotypes using SNP allele information.

    :param cooper:      cooper object
    :param locus_key:   locus key
    :return:            minimum SNP position used for clustering, or -1 if clustering fails
    """

    MIN_ALLELE_FRAC       = 0.2    # min fraction for allele quality calculation
    THRESHOLD_RANGES      = [(0.3, 0.7), (0.25, 0.75), (0.2, 0.8)]
    MIN_SNP_COVERAGE_FRAC = 0.6    # min fraction of reads a SNP must cover

    locus = cooper.cooper_loci_info[locus_key]
    locus_data = cooper.cooper_loci_data[locus_key]
    locus_cov  = locus_data.depth

    relevant_snps = set()    # SNPs relevant to the locus based on proximity
    for rindex in locus_data.reads:
        relevant_snps |= (cooper.cooper_read_data[rindex].snps)
    relevant_snps = sorted(list(filter(lambda x: (x in cooper.cooper_snp_data) and (cooper.cooper_snp_data[x].cov >= 3) and
                                                 (locus.start - cooper.args.snp_dist < x < locus.end + cooper.args.snp_dist),
                            relevant_snps)))

    relevant_snp_data = {}
    alt_snp_cov = {}
    for pos in relevant_snps:
        qual_check = False
        snp_data = cooper.cooper_snp_data[pos]
        max_coverage = 0

        passed_alts = set()
        for alt in snp_data.sub:
            alt_reads = snp_data.sub[alt].intersection(set(locus_data.reads))
            alt_cov   = len(alt_reads)
            if alt_cov == 0: continue
            if alt_cov > max_coverage: max_coverage = alt_cov

            avg_alt_qual = sum([snp_data.qual[idx] for idx in alt_reads])/alt_cov

            if avg_alt_qual <= cooper.args.snp_qual:
                qual_check = True
                break
            passed_alts.add(alt)
        if max_coverage == 0 or qual_check:
            continue
        alt_snp_cov[pos] = max_coverage

        relevant_snp_data[pos] = { 'cov': 0, 'alleles': {}, 'qual': {} }
        for alt in passed_alts:
            alt_reads = snp_data.sub[alt].intersection(set(locus_data.reads))
            relevant_snp_data[pos]['cov']           += len(alt_reads)
            relevant_snp_data[pos]['alleles'][alt]   = alt_reads
            relevant_snp_data[pos]['qual'][alt]      = snp_data.qual
        ref_reads = snp_data.ref.intersection(set(locus_data.reads))
        if len(ref_reads) > 0:
            relevant_snp_data[pos]['alleles']['r'] = ref_reads
            relevant_snp_data[pos]['cov']         += len(ref_reads)

    del_positions = list(filter(lambda x: relevant_snp_data[x]['cov'] < 5, relevant_snp_data.keys()))

    for pos in del_positions:
        del relevant_snp_data[pos]
    ordered_snp_on_cov = sorted(relevant_snp_data.keys(), key = lambda item : alt_snp_cov[item], reverse = True)

    sig_snp_data    = {}
    ordered_sig_snps = []

    for pos in ordered_snp_on_cov:
        snp_data  = relevant_snp_data[pos]
        pos_cov   = snp_data['cov']
        if pos_cov < 0.6 * locus_cov: break
        ordered_sig_snps.append(pos)
        sig_snp_data[pos] = { 'cov': snp_data['cov'], 'alleles': snp_data['alleles'], 'qual': snp_data['qual'] }

    min_snp_pos      = qvalue_phasing(cooper, locus, locus_data, sig_snp_data, ordered_sig_snps)
    # for tier_idx, (lower_thresh, upper_thresh) in enumerate(THRESHOLD_RANGES):

    #     sig_snp_data    = {}
    #     ordered_sig_snps = []

    #     for pos in ordered_snp_on_cov:
    #         snp_data  = relevant_snp_data[pos]
    #         pos_cov   = snp_data['cov']
    #         if pos_cov < 0.6 * locus_cov: break

    #         # count alleles where read fraction is within threshold bounds
    #         balanced_alleles = sum(
    #             lower_thresh * pos_cov <= len(reads) <= upper_thresh * pos_cov
    #             for reads in snp_data['alleles'].values()
    #         )

    #         if balanced_alleles >= 2:
    #             ordered_sig_snps.append(pos)
    #             sig_snp_data[pos] = { 'cov': snp_data['cov'], 'alleles': snp_data['alleles'], 'qual': snp_data['qual'] }

    #     if not ordered_sig_snps:
    #         if tier_idx < 2: continue
    #         return -1

    #     min_snp_pos = merge_snpreadsets(cooper, locus_data, sig_snp_data, ordered_sig_snps)

    #     if locus_data.is_genotyped or tier_idx == 2:
    #         break

    return min_snp_pos


def merge_snpreadsets(cooper, locus_data, sig_snp_data, ordered_sig_snps):
    """
    merge SNP read sets to phase reads into two haplotype clusters.

    :param cooper:              cooper object
    :param locus_data:          locus data object
    :param sig_snp_data:        significant SNP position data
    :param ordered_sig_snps:    SNP positions ordered by significance
    :return:                    [haplotypes, success, min_snp, skip_point, snp_quals, phased_reads, snp_count]
    """
    CLUSTER_JOIN_HIGH     = 0.7    # min intersection to join a cluster
    CLUSTER_JOIN_LOW      = 0.05   # max intersection to confirm non-membership

    del_positions = []
    for pos in ordered_sig_snps:
        remove_alts = []
        for alt, reads in sig_snp_data[pos]['alleles'].items():
            if len(reads) < 0.2 * sig_snp_data[pos]['cov']:
                remove_alts.append(alt)
        for alt in remove_alts:
            del sig_snp_data[pos]['alleles'][alt]
            if alt in sig_snp_data[pos]['qual']:
                del sig_snp_data[pos]['qual'][alt]
        if len(sig_snp_data[pos]['alleles']) < 2:
            del_positions.append(pos)
    for pos in del_positions:
        del sig_snp_data[pos]
        ordered_sig_snps.remove(pos)

    snps     = ordered_sig_snps[:cooper.args.snp_count]

    # --- compute quality values for top SNPs ---
    snp_quals = [max(list(sig_snp_data[pos]['qual'].values())) for pos in snps]
    snp_quals = ','.join(str(int(q)) for q in snp_quals)

    # --- compute pairwise mismatch scores between SNPs ---
    mismatch_scores = {}                  # {pos_a: {pos_b: mismatch_score}}

    for i, pos_a in enumerate(snps):
        if mismatch_scores and i == len(snps) - 1:
            break
        mismatch_scores[pos_a] = {}
        alleles_a = list(sig_snp_data[pos_a]['alleles'].values())

        for pos_b in snps[i + 1:]:
            score = sum(
                min(len(reads_b & allele_a) for allele_a in alleles_a)
                for reads_b in sig_snp_data[pos_b]['alleles'].values()
            )
            mismatch_scores[pos_a][pos_b] = score

    # --- select best SNPs for phasing ---
    # prefer positions with at least 2 zero-mismatch neighbours
    sig_snps = []
    for pos, neighbours in mismatch_scores.items():
        if sum(1 for s in neighbours.values() if s == 0) >= 2:
            sig_snps.append(pos)
            sig_snps.extend(sorted(neighbours, key=lambda p: neighbours[p]))
            break

    # fallback — pick position with lowest sum of two best scores
    if not sig_snps:
        best_pos = min(
            mismatch_scores,
            key=lambda p: sum(sorted(mismatch_scores[p].values())[:2])
        )
        sig_snps.append(best_pos)
        sig_snps.extend(
            sorted(mismatch_scores[best_pos], key=lambda p: mismatch_scores[best_pos][p])
        )

    # --- build final ordered SNP dict ---
    final_snp_dict = {
        pos: sig_snp_data[pos]['alleles']
        for pos in sig_snps
        if pos in sig_snp_data
    }
    min_snp_pos = min(final_snp_dict)

    # --- cluster reads into two haplotypes ---
    cluster1: set = set()
    cluster2: set = set()

    for alleles in final_snp_dict.values():
        for read_set in alleles.values():
            if not cluster1:
                cluster1 |= read_set
            elif not cluster2:
                cluster2 |= read_set
            elif (len(read_set & cluster1) > CLUSTER_JOIN_HIGH * len(read_set) and
                  len(read_set & cluster2) < CLUSTER_JOIN_LOW  * len(read_set)):
                cluster1 |= read_set
            elif (len(read_set & cluster2) > CLUSTER_JOIN_HIGH * len(read_set) and
                  len(read_set & cluster1) < CLUSTER_JOIN_LOW  * len(read_set)):
                cluster2 |= read_set

    # --- validate phasing coverage ---
    total_phased = len(cluster1) + len(cluster2)
    if total_phased >= cooper.args.phasing_read * locus_data.depth:
        locus_data.hap_read_sets      = (list(cluster1), list(cluster2))
        locus_data.hap_alen_sets      = ([locus_data.read_alens[ridx][0] for ridx in cluster1], [locus_data.read_alens[ridx][0] for ridx in cluster2])
        locus_data.is_genotyped       = True
        locus_data.phase_mode         = 'snp'
        locus_data.hap_category       = 0
        locus_data.n_phasing_snps     = len(final_snp_dict)
        locus_data.phasing_snp_quals  = snp_quals

        return min_snp_pos

    locus_data.is_genotyped    = False
    locus_data.skip_code     = 1
    return min_snp_pos


def score_cluster_membership(cluster, read_idx, read_indices, score_matrix):
    """
    Function to score the membership of a read in a cluster based on its supporting reads.

    Args:
        cluster (set): Set of read indices in the cluster.
        read_idx (int): Index of the read to be scored.
        read_indices (set): Set of all read indices.
        score_matrix (list): List of scores for each read pair.

    Returns:
        int: Score indicating the membership of the read in the cluster.
    """
    scores = []
    for read in cluster:
        i = read_indices.index(read)
        j = read_indices.index(read_idx)
        index = (i * len(read_indices)) + j
        scores.append(score_matrix[index])
    return sum(scores)/len(scores) if scores else 0


def qvalue_phasing(cooper, locus, locus_data, sig_snp_data, ordered_sig_snps):
    """
    Function to cluster the reads based on SNPs and their supporting reads.

    Args:
        snp_allelereads (dict): Dictionary containing SNP positions and their alleles with supporting reads.
        coverage_sorted_snps (list): List of SNP positions sorted by coverage.
        read_indices (set): Set of read indices that support the SNPs.
        args (Namespace): Command line arguments containing parameters like snpC, snpR, phasingR, snpQ.

    Returns:
        list: A list containing the final haplotypes, minimum SNP position, skip point,
    """
    del_positions = []
    for pos in ordered_sig_snps:
        remove_alts = []
        for alt, reads in sig_snp_data[pos]['alleles'].items():
            if len(reads) < 0.2 * sig_snp_data[pos]['cov']:
                remove_alts.append(alt)
        for alt in remove_alts:
            del sig_snp_data[pos]['alleles'][alt]
            if alt in sig_snp_data[pos]['qual']:
                del sig_snp_data[pos]['qual'][alt]
        if len(sig_snp_data[pos]['alleles']) < 2:
            del_positions.append(pos)
    for pos in del_positions:
        del sig_snp_data[pos]
        ordered_sig_snps.remove(pos)

    if not ordered_sig_snps: return -1

    if len(ordered_sig_snps) == 1:
        max_alt_cov = max(len(reads) for reads in sig_snp_data[ordered_sig_snps[0]]['alleles'].values())
        if max_alt_cov > 0.7 * locus_data.depth or max_alt_cov < 0.3 * locus_data.depth:
            return -1

    min_snp_pos = min(ordered_sig_snps)
    sorted_reads = sorted(locus_data.reads)

    # --- compute quality values for top SNPs ---
    # snp_quals = [max(list([max(x) for x in sig_snp_data[pos]['qual'].values()])) for pos in ordered_sig_snps]
    snp_quals = [max(list([max(list(x.values())) for x in sig_snp_data[pos]['qual'].values()])) for pos in ordered_sig_snps]
    snp_quals = ','.join(str(int(q)) for q in snp_quals)

    score_matrix = []
    for i, read_a in enumerate(sorted_reads):
        for j, read_b in enumerate(sorted_reads):
            if read_a == read_b: score_matrix.append(0); continue
            if j < i: score_matrix.append(score_matrix[j*len(sorted_reads) + i]); continue
            
            score = 0
            for snp_pos in ordered_sig_snps:
                allele_a = None; qual_a = None
                allele_b = None; qual_b = None
                for allele in sig_snp_data[snp_pos]['alleles']:
                    if read_a in sig_snp_data[snp_pos]['alleles'][allele]:
                        allele_a = allele
                        if allele != 'r': qual_a = sig_snp_data[snp_pos]['qual'][allele][read_a]
                        else: qual_a = cooper.cooper_read_data[read_a].mean_qual
                    if read_b in sig_snp_data[snp_pos]['alleles'][allele]:
                        allele_b = allele
                        if allele != 'r': qual_b = sig_snp_data[snp_pos]['qual'][allele][read_b]
                        else: qual_b = cooper.cooper_read_data[read_b].mean_qual
                if allele_a is None or allele_b is None:
                    pass
                elif allele_a == allele_b:
                    pass
                elif qual_a < cooper.args.snp_qual or qual_b < cooper.args.snp_qual:
                    pass
                else:
                    score -= (qual_a + qual_b)
            score_matrix.append(score)

    sorted_indices = sorted(range(len(score_matrix)), key=score_matrix.__getitem__)

    num_reads = len(sorted_reads)
    cluster1 = set()
    cluster2 = set()
    phased_reads = set()
    for itr, idx in enumerate(sorted_indices):
        if score_matrix[idx] == 0: break
        read_a = sorted_reads[idx // num_reads]
        read_b = sorted_reads[idx % num_reads]

        if read_a in phased_reads and read_b in phased_reads: continue

        if itr == 0:
            cluster1.add(read_a)
            cluster2.add(read_b)
            phased_reads.add(read_a); phased_reads.add(read_b)
            continue
        if read_b not in phased_reads:
            score_b1 = score_cluster_membership(cluster1, read_b, sorted_reads, score_matrix)
            score_b2 = score_cluster_membership(cluster2, read_b, sorted_reads, score_matrix)
            if   score_b1 < score_b2:
                cluster2.add(read_b)
            elif score_b1 > score_b2:
                cluster1.add(read_b)
            phased_reads.add(read_b)
        if read_a not in phased_reads:
            score_a1 = score_cluster_membership(cluster1, read_a, sorted_reads, score_matrix)
            score_a2 = score_cluster_membership(cluster2, read_a, sorted_reads, score_matrix)
            if   score_a1 < score_a2:
                cluster2.add(read_a)
            elif score_a1 > score_a2:
                cluster1.add(read_a)
            phased_reads.add(read_a)

        if len(phased_reads) == len(sorted_reads):
            break

    for read_a in sorted_reads:
        if read_a not in phased_reads:
            score_a1 = score_cluster_membership(cluster1, read_a, sorted_reads, score_matrix)
            score_a2 = score_cluster_membership(cluster2, read_a, sorted_reads, score_matrix)
            if   score_a1 < score_a2: cluster2.add(read_a)
            elif score_a1 > score_a2: cluster1.add(read_a)
            phased_reads.add(read_a)

    # --- validate phasing coverage ---
    total_phased = len(cluster1) + len(cluster2)
    if total_phased >= cooper.args.phasing_read * locus_data.depth:
        locus_data.hap_read_sets      = (list(cluster1), list(cluster2))
        locus_data.hap_alen_sets      = ([locus_data.read_alens[ridx][0] for ridx in cluster1], [locus_data.read_alens[ridx][0] for ridx in cluster2])
        locus_data.is_genotyped       = True
        locus_data.phase_mode         = 'snp'
        locus_data.hap_category       = 0
        locus_data.n_phasing_snps     = len(sig_snp_data)
        locus_data.phasing_snp_quals  = snp_quals

        return min_snp_pos

    locus_data.is_genotyped    = False
    locus_data.skip_code     = 1
    return min_snp_pos