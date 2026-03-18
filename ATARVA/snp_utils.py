import numpy as np

def haplocluster_reads(cooper, locus_key):
    """
    cluster reads into haplotypes using SNP allele information.

    :param cooper:                  cooper object
    :param locus_data:              locus data object
    :param relevant_snp_data:       dict of locus relevant SNP data at the locus
    :param ordered_snp_on_cov:      SNP positions ordered by alt allele coverage
    :return:                        [haplotypes, min_snp, skip_point, snp_quals, phased_reads, snp_count]
    """

    MIN_ALLELE_FRAC       = 0.2    # min fraction for allele quality calculation

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
        snp_data = cooper.cooper_snp_data[pos]
        
        # filter to get reads that are relevant to the locus
        ref_reads = snp_data.ref.intersection(set(locus_data.reads))
        pos_reads = ref_reads.copy()
        for alt in snp_data.sub:
            pos_reads |= snp_data.sub[alt]
        pos_reads = pos_reads.intersection(set(locus_data.reads))

        pos_cov = len(pos_reads)
        if pos_cov < 5 or pos_cov < locus_cov * MIN_SNP_COVERAGE_FRAC: continue    # at least 5 relevant reads needed at SNP position

        passed_alleles = []
        alt_covs       = []
        alt_reads      = {}
        for alt in snp_data.sub:
            alt_reads[alt] = snp_data.sub[alt].intersection(set(locus_data.reads))
            alt_cov   = len(alt_reads[alt])
            if alt_cov == 0: continue

            avg_alt_qual = np.mean([snp_data.qual[idx] for idx in alt_reads[alt]])

            if avg_alt_qual >= cooper.args.snp_qual and alt_cov >= pos_cov * MIN_ALLELE_FRAC:
                # alternate allele passes when having 20% read support and average quality above threshold
                passed_alleles.append(alt)
                alt_covs.append(alt_cov)
        
        if len(ref_reads) >= pos_cov * MIN_ALLELE_FRAC: passed_alleles.append('r')

        if len(passed_alleles) < 2: continue # if no alleles passed quality and coverage thresholds, skip the SNP

        alt_snp_cov[pos] = max(alt_covs)

        relevant_snp_data[pos] = { 'cov': 0, 'alleles': {}, 'qual': {} }
        for alt in passed_alleles:
            relevant_snp_data[pos]['cov']           += len(alt_reads[alt])
            relevant_snp_data[pos]['alleles'][alt]   = alt_reads[alt]
            relevant_snp_data[pos]['qual'][alt]      = np.mean([snp_data.qual[idx] for idx in alt_reads[alt]])
        if len(snp_data.ref) > 0:
            relevant_snp_data[pos]['alleles']['r'] = ref_reads
            relevant_snp_data[pos]['cov']         += len(ref_reads)

    ordered_snp_on_cov = sorted(relevant_snp_data.keys(), key = lambda item : alt_snp_cov[item], reverse = True)

    THRESHOLD_RANGES      = [(0.3, 0.7), (0.25, 0.75), (0.2, 0.8)]
    FAIL_RESULT           = [(), -1, 0, '', '', 0]
    MIN_SNP_COVERAGE_FRAC = 0.6    # min fraction of reads a SNP must cover

    final_haplotypes = ()
    min_snp_pos      = -1
    skip_point       = 10
    min_snp_coverage = MIN_SNP_COVERAGE_FRAC * locus_cov

    for tier_idx, (lower_thresh, upper_thresh) in enumerate(THRESHOLD_RANGES):

        sig_snp_data    = {}
        ordered_sig_snps = []

        for pos in ordered_snp_on_cov:
            snp_data  = relevant_snp_data[pos]
            pos_cov   = snp_data['cov']

            # count alleles where read fraction is within threshold bounds
            balanced_alleles = sum(
                lower_thresh * pos_cov <= len(reads) <= upper_thresh * pos_cov
                for reads in snp_data['alleles'].values()
            )

            if balanced_alleles >= 2:
                ordered_sig_snps.append(pos)
                sig_snp_data[pos] = { 'cov': snp_data['cov'], 'alleles': snp_data['alleles'], 'qual': snp_data['qual'] }

        if not ordered_sig_snps:
            if tier_idx < 2: continue
            return FAIL_RESULT

        
        min_snp_pos = merge_snpreadsets(cooper, locus_data, sig_snp_data, ordered_sig_snps)

        if locus_data.hap_status or tier_idx == 2:
            break

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
        alleles_a = list(sig_snp_data[pos_a].values())

        for pos_b in snps[i + 1:]:
            score = sum(
                min(len(reads_b & allele_a) for allele_a in alleles_a)
                for reads_b in sig_snp_data[pos_b].values()
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
        pos: sig_snp_data[pos]
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
        locus_data.haplotypes   = (list(cluster1), list(cluster2))
        locus_data.hap_status   = True
        locus_data.phase_mode   = 'snp'
        locus_data.hap_category = 0

        return min_snp_pos
