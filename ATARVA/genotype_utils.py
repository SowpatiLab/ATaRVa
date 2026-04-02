import numpy as np
import warnings
from sklearn.cluster import KMeans
from threadpoolctl import threadpool_limits

from ATARVA.snp_utils import haplocluster_reads
from ATARVA.vcf_writer import *
from ATARVA.sub_operation_utils import alt_sequence
from ATARVA.somatic_utils import *

# NOTE: after genotyping and haplogrouping, before reporting the locus build the ALT sequence
#       and calculate methylation level for the haplotypes


def homozygous_call(cooper, locus_key):
    """
    genotype a homozygous locus and build the ALT sequence for the haplotype
    :param cooper:     cooper object
    :param locus_key:  key for the locus
    :return:           [bool_state, category]
    """

    locus        = cooper.cooper_loci_info[locus_key]
    locus_data   = cooper.cooper_loci_data[locus_key]
    hap_reads    = locus_data.hap_read_sets[0]
    hap_lengths  = locus_data.hap_alen_sets[0]

    lower, upper = (round(x) for x in np.percentile(np.array(hap_lengths), [2.5, 97.5]))

    ALT, allele_length, decomposed_seq, is_repetitive = alt_sequence(locus_data.read_aseqs, hap_reads, locus.motif_length)
    if not is_repetitive:
        locus_data.skip_code = 6
        return

    locus_data.gt_aseqs        = (ALT, None)
    locus_data.gt_alens        = (allele_length, None)
    locus_data.gt_decomp_seqs  = (decomposed_seq, None)
    locus_data.hap_meth_data   = (calculate_methylation(hap_reads, locus_data.read_methylation, ALT), None)
    locus_data.gt_arange       = (f'{lower}-{upper}', None)

    write_homozygous_call(cooper, locus_key)
    return


def heterozygous_call(cooper, locus_key):
    """
    genotype a heterozygous locus and build the ALT sequences for the haplotypes
    :param cooper:          cooper object
    :param locus_key:       key for the locus
    :return:                [bool_state, category]
    """

    locus      = cooper.cooper_loci_info[locus_key]
    locus_data = cooper.cooper_loci_data[locus_key]
    hap_read_sets = locus_data.hap_read_sets
    hap_alen_sets = locus_data.hap_alen_sets

    for i in range(2):
        hap_reads = hap_read_sets[i]
        hap_lengths = hap_alen_sets[i]
        ALT, allele_length, decomp_seq, is_repetitive = alt_sequence(locus_data.read_aseqs, hap_reads, locus.motif_length)
        lower, upper = (round(x) for x in np.percentile(np.array(hap_lengths), [2.5, 97.5]))
        if i == 0:
            locus_data.gt_aseqs         = (ALT, locus_data.gt_aseqs[1])
            locus_data.gt_alens         = (allele_length, locus_data.gt_alens[1])
            locus_data.gt_decomp_seqs   = (decomp_seq, locus_data.gt_decomp_seqs[1])
            locus_data.hap_meth_data    = (calculate_methylation(hap_reads, locus_data.read_methylation, ALT), locus_data.hap_meth_data[1])
            locus_data.gt_arange        = (f'{lower}-{upper}', locus_data.gt_arange[1])
        else:
            locus_data.gt_aseqs         = (locus_data.gt_aseqs[0], ALT)
            locus_data.gt_alens         = (locus_data.gt_alens[0], allele_length)
            locus_data.gt_decomp_seqs   = (locus_data.gt_decomp_seqs[0], decomp_seq)
            locus_data.hap_meth_data    = (locus_data.hap_meth_data[0], calculate_methylation(hap_reads, locus_data.read_methylation, ALT))
            locus_data.gt_arange        = (locus_data.gt_arange[0], f'{lower}-{upper}')

    write_heterozygous_call(cooper, locus_key)

    return [True, 10]


def compute_cluster_cutoff(minor_cluster, major_cluster):
    """
    compute minimum read cutoff for the minor cluster based on
    whether it overlaps with the major cluster's allele range.

    :param minor_cluster: smaller allele length cluster
    :param major_cluster: larger allele length cluster
    :return:              minimum read count cutoff
    """
    max_major    = max(major_cluster)
    tolerance    = max(max_major * 0.1, 10)        # 10% of max or at least 10bp
    lower_bound  = min(major_cluster) - tolerance
    upper_bound  = max_major          + tolerance

    # if minor cluster overlaps with major - no cutoff needed
    overlaps = any(lower_bound <= alen <= upper_bound for alen in minor_cluster)
    if overlaps: return 0.15 * (len(major_cluster) + len(minor_cluster)) # 15% of total reads in both clusters

    # min 3% of major cluster size, at least 2 reads
    ratio_cutoff = int(max(0.03, len(minor_cluster) / len(major_cluster)) * len(major_cluster))
    return max(2, ratio_cutoff)


def length_genotyper(cooper, locus_key):
    """
    genotype a locus by clustering allele lengths using KMeans.

    :param cooper:     cooper object
    :param locus_key:  key for the locus
    :return:           [bool_state, category]
    """

    MIN_READS        = 3
    MIN_CLUSTER_FRAC = 0.15
    WINDOW_FRAC      = 0.1

    locus_data = cooper.cooper_loci_data[locus_key]
    read_alens = locus_data.read_alens

    read_indices    = sorted(locus_data.reads)
    unique_alens    = set(locus_data.allele_lengths)
    singleton_alens = {alen for alen, count in locus_data.halen_frequency.items() if count == 1}

    # --- pre-compute 10% windows for each unique allele ---
    windows = { i: (round(i * (1 - WINDOW_FRAC)), round(i * (1 + WINDOW_FRAC))) for i in unique_alens }

    # --- filter singleton alleles not near any other allele ---
    main_read_ids  = []
    filtered_alens = []

    for read_index in read_indices:
        alen = read_alens[read_index][0]
        if alen in singleton_alens:
            near_other = any(lo <= alen <= hi for i, (lo, hi) in windows.items() if i != alen)
            # NOTE: filters out singleton alleles which are not within 10% of any other allele,
            #       these outlier can potentially create spurious clusters
            if not near_other: continue
        main_read_ids.append(read_index)
        filtered_alens.append(alen)

    if len(filtered_alens) < MIN_READS:
        locus_data.skip_code = 0
        return

    # --- KMeans clustering ---
    alen_array = np.array(filtered_alens).reshape(-1, 1)
    with threadpool_limits(limits=1):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=UserWarning)
            kmeans = KMeans(n_clusters=2, init='k-means++', n_init=5, random_state=0).fit(alen_array)

    # --- split into clusters in single pass ---
    c1_idx, c2_idx = [], []
    for i, label in enumerate(kmeans.labels_):
        (c1_idx if label == 0 else c2_idx).append(i)

    c1_lengths    = [filtered_alens[i] for i in c1_idx]
    c2_lengths    = [filtered_alens[i] for i in c2_idx]
    hap_read_sets = ([main_read_ids[i] for i in c1_idx], [main_read_ids[i] for i in c2_idx])

    locus_data.gt_depth = len(filtered_alens)
    # --- compute cluster cutoff ---
    min_cluster_size = MIN_CLUSTER_FRAC * len(filtered_alens)

    if c1_idx and c2_idx:
        if len(c1_idx) < min_cluster_size <= len(c2_idx):
            min_cluster_size = compute_cluster_cutoff(c1_lengths, c2_lengths)
        elif len(c2_idx) < min_cluster_size <= len(c1_idx):
            min_cluster_size = compute_cluster_cutoff(c2_lengths, c1_lengths)

    # --- haploid genotyping ---
    if cooper.haploid:
        # we assign the larger cluster to the haploid alelles
        major_idx   = 0 if len(c1_idx) >= len(c2_idx) else 1
        major_alens = c1_lengths    if major_idx == 0 else c2_lengths
        major_reads = hap_read_sets[major_idx]
        locus_data.hap_read_sets = (major_reads, [])
        locus_data.hap_alen_sets = (major_alens, None)

        if len(major_reads) < min_cluster_size:
            locus_data.skip_code = 1
            return

        homozygous_call(cooper, locus_key)
        return

    # --- diploid genotyping ---
    c1_valid = bool(c1_idx) and len(c1_idx) >= min_cluster_size
    c2_valid = bool(c2_idx) and len(c2_idx) >= min_cluster_size

    if c1_valid and c2_valid:
        locus_data.hap_read_sets = hap_read_sets
        locus_data.hap_alen_sets = (c1_lengths, c2_lengths)
        heterozygous_call(cooper, locus_key)
        return

    if c1_valid:
        locus_data.hap_read_sets = (hap_read_sets[0], [])
        locus_data.hap_alen_sets = (c1_lengths, None)
        homozygous_call(cooper, locus_key)
        return

    if c2_valid:
        locus_data.hap_read_sets = (hap_read_sets[1], [])
        locus_data.hap_alen_sets = (c2_lengths, None)
        homozygous_call(cooper, locus_key)
        return

    return


def analyse_genotype(cooper, locus_key):
    """
    genotype the locus based on the read data

    :param cooper: Cooper object
    :param locus_key: locus identifier string
    """

    locus      = cooper.cooper_loci_info[locus_key]
    locus_data = cooper.cooper_loci_data[locus_key]

    if cooper.haploid:
        length_genotyper(cooper, locus_key)
        return

    min_snp = haplocluster_reads(cooper, locus_key)

    if not locus_data.is_genotyped: # if the loci has no significant snps
        length_genotyper(cooper, locus_key)
        return

    if min_snp != -1:
        snp_left_boundary = locus.start - cooper.args.snp_dist
        min_idx = 0
        for each_spos in cooper.cooper_sorted_snps:
            if each_spos >= snp_left_boundary:
                break
            del cooper.cooper_snp_data[each_spos]
            min_idx += 1
        del cooper.cooper_sorted_snps[:min_idx]

    for i, hap_reads in enumerate(locus_data.hap_read_sets):
        hap_lengths = locus_data.hap_alen_sets[i]
        ALT, allele_length, decomp_seq, is_repetitive = alt_sequence(locus_data.read_aseqs, hap_reads, locus.motif_length)
        lower, upper = (round(x) for x in np.percentile(np.array(hap_lengths), [2.5, 97.5]))
        if i == 0:
            locus_data.gt_aseqs        = (ALT, locus_data.gt_aseqs[1])
            locus_data.gt_alens        = (allele_length, locus_data.gt_alens[1])
            locus_data.gt_decomp_seqs  = (decomp_seq, locus_data.gt_decomp_seqs[1])
            locus_data.hap_meth_data   = (calculate_methylation(hap_reads, locus_data.read_methylation, ALT), locus_data.hap_meth_data[1])
            locus_data.gt_arange       = (f'{lower}-{upper}', locus_data.gt_arange[1])
        else:
            locus_data.gt_aseqs        = (locus_data.gt_aseqs[0], ALT)
            locus_data.gt_alens        = (locus_data.gt_alens[0], allele_length)
            locus_data.gt_decomp_seqs  = (locus_data.gt_decomp_seqs[0], decomp_seq)
            locus_data.hap_meth_data   = (locus_data.hap_meth_data[0], calculate_methylation(hap_reads, locus_data.read_methylation, ALT))
            locus_data.gt_arange       = (locus_data.gt_arange[0], f'{lower}-{upper}')
    write_heterozygous_call(cooper, locus_key)

    return
