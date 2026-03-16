# --- Constants ---

THRESHOLD_RANGES      = [(0.3, 0.7), (0.25, 0.75), (0.2, 0.8)]
MIN_SNP_COVERAGE_FRAC = 0.6    # min fraction of reads a SNP must cover
MIN_ALLELE_FRAC       = 0.2    # min fraction for allele quality calculation
CLUSTER_JOIN_HIGH     = 0.7    # min intersection to join a cluster
CLUSTER_JOIN_LOW      = 0.05   # max intersection to confirm non-membership
FAIL_RESULT           = [(), -1, 0, '', '', 0]

def haplocluster_reads(snp_allele_reads, ordered_snps_by_coverage, read_indices, qual_threshold, num_snps,
                       snp_cov_threshold, phased_fraction_threshold):
    """
    cluster reads into haplotypes using SNP allele information.

    :param snp_allele_reads:          SNP position → allele → read indices
    :param ordered_snps_by_coverage:  SNP positions ordered by coverage
    :param read_indices:              all read indices at this locus
    :param qual_threshold:            minimum base quality at SNP position
    :param num_snps:                  number of SNPs to use for phasing
    :param snp_cov_threshold:         minimum SNP coverage threshold
    :param phased_fraction_threshold: min fraction of reads that must be phased
    :return:                          [haplotypes, min_snp, skip_point, snp_quals, phased_reads, snp_count]
    """

    final_haplotypes = ()
    min_snp_pos      = -1
    skip_point       = 10
    n_reads          = len(read_indices)
    min_snp_coverage = MIN_SNP_COVERAGE_FRAC * n_reads

    for tier_idx, (lower_thresh, upper_thresh) in enumerate(THRESHOLD_RANGES):

        sig_snp_data    = {}
        ordered_sig_snps = []

        for pos in ordered_snps_by_coverage:
            snp         = snp_allele_reads[pos]
            total_allele_reads = sum(len(reads) for reads in snp['alleles'].values())

            if snp['cov'] < min_snp_coverage:
                break

            # count alleles where read fraction is within threshold bounds
            balanced_alleles = sum(
                lower_thresh * total_allele_reads <= len(reads) <= upper_thresh * total_allele_reads
                for reads in snp['alleles'].values()
            )

            if balanced_alleles >= 2:
                ordered_sig_snps.append(pos)
                sig_snp_data[pos] = { 'cov': snp['cov'], 'alleles': snp['alleles'], 'qual': snp['qual'] }

        if not sig_snp_data:
            if tier_idx < 2: continue
            return FAIL_RESULT

        (final_haplotypes,
         success,
         min_snp_pos,
         skip_point,
         snp_quals,
         phased_reads,
         snp_count) = merge_snpreadsets(sig_snp_data, ordered_sig_snps, read_indices, qual_threshold, num_snps,
                                        snp_cov_threshold, phased_fraction_threshold)

        if success or tier_idx == 2:
            break

    return [final_haplotypes, min_snp_pos, skip_point, snp_quals, phased_reads, snp_count]


def merge_snpreadsets(sig_snp_data, ordered_sig_snps, read_indices, snp_qual_threshold, max_snps, min_allele_read_frac, phased_fraction_threshold):
    """
    merge SNP read sets to phase reads into two haplotype clusters.

    :param sig_snp_data:              significant SNP position data
    :param ordered_sig_snps:          SNP positions ordered by significance
    :param read_indices:              all read indices at locus
    :param snp_qual_threshold:        minimum base quality threshold
    :param max_snps:                  max SNPs to use for phasing
    :param min_allele_read_frac:      min read fraction to keep an allele
    :param phased_fraction_threshold: min fraction of reads that must be phased
    :return:                          [haplotypes, success, min_snp, skip_point, snp_quals, phased_reads, snp_count]
    """
    skip_point   = 10
    top_snps     = ordered_sig_snps[:max_snps]
    phased_reads = ''

    # --- compute quality values for top SNPs ---
    snp_qual_list = []
    for pos in top_snps:
        snp        = sig_snp_data[pos]
        total_reads = snp['cov']
        allele_quals = {
            sum(snp['qual'][r] for r in reads) / len(reads)
            for allele, reads in snp['alleles'].items()
            if allele != 'r' and len(reads) / total_reads >= MIN_ALLELE_FRAC
        }
        if allele_quals:
            snp_qual_list.append(str(int(max(allele_quals))))

    snp_quals = ','.join(snp_qual_list)
    snp_count = len(top_snps)

    if not top_snps:
        return [(), False, -1, 5, snp_quals, phased_reads, snp_count]

    # --- filter alleles by read fraction ---
    filtered_snps = {pos: dict(sig_snp_data[pos]['alleles']) for pos in top_snps}

    for pos in list(filtered_snps):
        alleles   = filtered_snps[pos]
        tot_reads = sum(len(reads) for reads in alleles.values())

        # remove alleles below read fraction threshold
        filtered_snps[pos] = {
            allele: reads
            for allele, reads in alleles.items()
            if len(reads) / tot_reads >= min_allele_read_frac
        }

        # remove positions with fewer than 2 valid alleles
        if len(filtered_snps[pos]) < 2:
            del filtered_snps[pos]

    if not filtered_snps:
        return [(), False, -1, 1, snp_quals, phased_reads, snp_count]

    # --- compute pairwise mismatch scores between SNPs ---
    snp_positions  = list(filtered_snps)
    mismatch_scores = {}                  # {pos_a: {pos_b: mismatch_score}}

    for i, pos_a in enumerate(snp_positions):
        if mismatch_scores and i == len(snp_positions) - 1:
            break
        mismatch_scores[pos_a] = {}
        alleles_a = list(filtered_snps[pos_a].values())

        for pos_b in snp_positions[i + 1:]:
            score = sum(
                min(len(reads_b & allele_a) for allele_a in alleles_a)
                for reads_b in filtered_snps[pos_b].values()
            )
            mismatch_scores[pos_a][pos_b] = score

    # --- select best SNPs for phasing ---
    # prefer positions with at least 2 zero-mismatch neighbours
    sig_snps = []
    for pos, neighbours in mismatch_scores.items():
        if sum(1 for s in neighbours.values() if s == 0) >= 2:
            sig_snps.append(pos)
            sig_snps.extend(
                sorted(neighbours, key=lambda p: neighbours[p])
            )
            break

    # fallback — pick position with lowest sum of two best scores
    if not sig_snps:
        best_pos = min(
            mismatch_scores,
            key=lambda p: sum(
                sorted(mismatch_scores[p].values())[:2]
            )
        )
        sig_snps.append(best_pos)
        sig_snps.extend(
            sorted(mismatch_scores[best_pos],
                   key=lambda p: mismatch_scores[best_pos][p])
        )

    # --- build final ordered SNP dict ---
    final_snp_dict = {
        pos: filtered_snps[pos]
        for pos in sig_snps
        if pos in filtered_snps
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
    if total_phased >= phased_fraction_threshold * len(read_indices):
        phased_reads = [len(cluster1), len(cluster2)]
        return [(cluster1, cluster2), True, min_snp_pos,
                skip_point, snp_quals, phased_reads, snp_count]

    return [(), False, -1, 2, snp_quals, phased_reads, snp_count]