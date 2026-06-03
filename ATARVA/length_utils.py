import warnings
from sklearn.cluster import KMeans
from threadpoolctl import threadpool_limits
from sklearn.mixture import GaussianMixture
from scipy.signal import find_peaks
from hdbscan import HDBSCAN
import numpy as np

import stringzilla as sz

from ATARVA.vcf_writer import *
from ATARVA.sub_operation_utils import alt_sequence, calculate_methylation


def _assign_genotype(cooper, locus_key, locus_data, c1_idx, c2_idx, c1_lengths,
                     c2_lengths, hap_read_sets, min_cluster_size):
    """
    Shared logic to assign genotype after clustering, handles both haploid and diploid cases.
    :param cooper:            cooper object
    :param locus_key:         key for the locus
    :param locus_data:        locus data object to update
    :param c1_idx:            indices of reads in cluster 1
    :param c2_idx:            indices of reads in cluster 2
    :param c1_lengths:        allele lengths for cluster 1
    :param c2_lengths:        allele lengths for cluster 2
    :param hap_read_sets:     tuple of (cluster 1 read indices, cluster 2 read indices)
    :param min_cluster_size:  minimum read count cutoff for a valid cluster
    """

    # ── haploid ───────────────────────────────────────────────────────
    if cooper.haploid:
        major_idx   = 0 if len(c1_idx) >= len(c2_idx) else 1
        major_alens = c1_lengths if major_idx == 0 else c2_lengths
        major_reads = hap_read_sets[major_idx]

        locus_data.hap_read_sets = (major_reads, [])
        locus_data.hap_alen_sets = (major_alens, None)

        if len(major_reads) < min_cluster_size:
            locus_data.skip_code = 1
            return

        homozygous_call(cooper, locus_key)
        return

    # ── diploid ───────────────────────────────────────────────────────
    c1_valid = bool(c1_idx) and len(c1_idx) >= min_cluster_size
    c2_valid = bool(c2_idx) and len(c2_idx) >= min_cluster_size

    if c1_valid and c2_valid:
        locus_data.hap_read_sets = hap_read_sets
        locus_data.hap_alen_sets = (c1_lengths, c2_lengths)
        heterozygous_call(cooper, locus_key)
        return

    if c1_valid:
        locus_data.hap_read_sets = (hap_read_sets[0], hap_read_sets[0])
        locus_data.hap_alen_sets = (c1_lengths, c1_lengths)
        homozygous_call(cooper, locus_key)
        return

    if c2_valid:
        locus_data.hap_read_sets = (hap_read_sets[1], hap_read_sets[1])
        locus_data.hap_alen_sets = (c2_lengths, c2_lengths)
        homozygous_call(cooper, locus_key)
        return

    locus_data.skip_code = 6


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
    # if not is_repetitive:
    #     locus_data.skip_code = 6
    #     return

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

    _assign_genotype(cooper, locus_key, locus_data,
                     c1_idx, c2_idx, c1_lengths, c2_lengths,
                     hap_read_sets, min_cluster_size)


def length_genotyper_gmm(cooper, locus_key):
    """
    Genotype using Gaussian Mixture Model.
    
    :param cooper:     cooper object
    :param locus_key:  key for the locus
    """

    MIN_READS        = 3
    MIN_CLUSTER_FRAC = 0.15
    WINDOW_FRAC      = 0.1

    locus_data = cooper.cooper_loci_data[locus_key]
    read_alens = locus_data.read_alens

    read_indices    = sorted(locus_data.reads)
    unique_alens    = set(locus_data.allele_lengths)
    singleton_alens = {alen for alen, count in locus_data.halen_frequency.items() if count == 1}
    windows         = {i: (round(i * (1 - WINDOW_FRAC)), round(i * (1 + WINDOW_FRAC)))
                       for i in unique_alens}

    main_read_ids  = []
    filtered_alens = []
    for read_index in read_indices:
        alen = read_alens[read_index][0]
        if alen in singleton_alens:
            if not any(lo <= alen <= hi for i, (lo, hi) in windows.items() if i != alen):
                continue
        main_read_ids.append(read_index)
        filtered_alens.append(alen)

    if len(filtered_alens) < MIN_READS:
        locus_data.skip_code = 0
        return

    alen_array = np.array(filtered_alens).reshape(-1, 1)

    # ── GMM fitting ───────────────────────────────────────────────────
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=UserWarning)
        gmm    = GaussianMixture(
            n_components    = 2,
            covariance_type = 'full',
            random_state    = 0,
            n_init          = 5
        ).fit(alen_array)

    labels = gmm.predict(alen_array)
    probs  = gmm.predict_proba(alen_array)   # confidence per assignment

    # ── split clusters ────────────────────────────────────────────────
    c1_idx = [i for i, l in enumerate(labels) if l == 0]
    c2_idx = [i for i, l in enumerate(labels) if l == 1]

    # filter low confidence assignments
    MIN_PROB      = 0.7
    c1_idx = [i for i in c1_idx if probs[i][0] >= MIN_PROB]
    c2_idx = [i for i in c2_idx if probs[i][1] >= MIN_PROB]

    c1_lengths    = [filtered_alens[i] for i in c1_idx]
    c2_lengths    = [filtered_alens[i] for i in c2_idx]
    hap_read_sets = (
        [main_read_ids[i] for i in c1_idx],
        [main_read_ids[i] for i in c2_idx]
    )

    locus_data.gt_depth      = len(filtered_alens)
    min_cluster_size         = MIN_CLUSTER_FRAC * len(filtered_alens)

    if c1_idx and c2_idx:
        if len(c1_idx) < min_cluster_size <= len(c2_idx):
            min_cluster_size = compute_cluster_cutoff(c1_lengths, c2_lengths)
        elif len(c2_idx) < min_cluster_size <= len(c1_idx):
            min_cluster_size = compute_cluster_cutoff(c2_lengths, c1_lengths)
    
    _assign_genotype(cooper, locus_key, locus_data,
                     c1_idx, c2_idx, c1_lengths, c2_lengths,
                     hap_read_sets, min_cluster_size)



def length_genotyper_hdbscan(cooper, locus_key):
    """
    Genotype using HDBSCAN — auto cluster count, outlier aware.
    
    :param cooper:     cooper object
    :param locus_key:  key for the locus
    """

    MIN_READS        = 3
    MIN_CLUSTER_FRAC = 0.15
    WINDOW_FRAC      = 0.1

    locus_data = cooper.cooper_loci_data[locus_key]
    locus      = cooper.cooper_loci_info[locus_key]
    read_alens = locus_data.read_alens

    read_indices    = sorted(locus_data.reads)
    unique_alens    = set(locus_data.allele_lengths)
    singleton_alens = {alen for alen, count in locus_data.halen_frequency.items() if count == 1}
    windows         = {i: (round(i * (1 - WINDOW_FRAC)), round(i * (1 + WINDOW_FRAC)))
                       for i in unique_alens}

    main_read_ids  = []
    filtered_alens = []
    for read_index in read_indices:
        alen = read_alens[read_index][0]
        if alen in singleton_alens:
            if not any(lo <= alen <= hi for i, (lo, hi) in windows.items() if i != alen):
                continue
        main_read_ids.append(read_index)
        edit_dist = sz.edit_distance(cooper.ref.fetch(cooper.chrom, locus.start, locus.end), locus_data.read_aseqs[read_index][0])
        filtered_alens.append([edit_dist, alen])

    if len(filtered_alens) < MIN_READS:
        locus_data.skip_code = 0
        return

    alen_array = np.array(filtered_alens)
    alen_array_normalized = (alen_array - alen_array.mean(axis=0)) / (alen_array.std(axis=0) + 1e-8)

    # ── HDBSCAN clustering ────────────────────────────────────────────
    clusterer = HDBSCAN(
        min_cluster_size     = max(MIN_READS, int(MIN_CLUSTER_FRAC * len(filtered_alens))),
        allow_single_cluster = True    # handles homozygous
    ).fit(alen_array)

    labels     = clusterer.labels_
    unique_labels = set(labels) - {-1}   # -1 = noise/outlier

    if not unique_labels:
        locus_data.skip_code = 0
        return

    # ── build clusters — ignore noise reads ──────────────────────────
    clusters = {
        label: [i for i, l in enumerate(labels) if l == label]
        for label in unique_labels
    }

    # ── take two largest clusters ─────────────────────────────────────
    top2      = sorted(clusters, key=lambda l: len(clusters[l]), reverse=True)[:2]
    c1_idx    = clusters[top2[0]]
    c2_idx    = clusters[top2[1]] if len(top2) > 1 else []

    c1_lengths    = [alen_array[i, 1] for i in c1_idx]
    c2_lengths    = [alen_array[i, 1] for i in c2_idx]
    hap_read_sets = (
        [main_read_ids[i] for i in c1_idx],
        [main_read_ids[i] for i in c2_idx]
    )

    locus_data.gt_depth  = len(filtered_alens)
    min_cluster_size     = MIN_CLUSTER_FRAC * len(filtered_alens)

    if c1_idx and c2_idx:
        if len(c1_idx) < min_cluster_size <= len(c2_idx):
            min_cluster_size = compute_cluster_cutoff(c1_lengths, c2_lengths)
        elif len(c2_idx) < min_cluster_size <= len(c1_idx):
            min_cluster_size = compute_cluster_cutoff(c2_lengths, c1_lengths)

    _assign_genotype(cooper, locus_key, locus_data,
                     c1_idx, c2_idx, c1_lengths, c2_lengths,
                     hap_read_sets, min_cluster_size)



def length_genotyper_histogram(cooper, locus_key):
    """
    Genotype using histogram peak detection.
    
    :param cooper:     cooper object
    :param locus_key:  key for the locus
    """

    MIN_READS        = 3
    MIN_CLUSTER_FRAC = 0.15
    WINDOW_FRAC      = 0.1
    BIN_WIDTH        = 5

    locus_data = cooper.cooper_loci_data[locus_key]
    read_alens = locus_data.read_alens

    read_indices    = sorted(locus_data.reads)
    unique_alens    = set(locus_data.allele_lengths)
    singleton_alens = {alen for alen, count in locus_data.halen_frequency.items() if count == 1}
    windows         = {i: (round(i * (1 - WINDOW_FRAC)), round(i * (1 + WINDOW_FRAC)))
                       for i in unique_alens}

    main_read_ids  = []
    filtered_alens = []
    for read_index in read_indices:
        alen = read_alens[read_index][0]
        if alen in singleton_alens:
            if not any(lo <= alen <= hi for i, (lo, hi) in windows.items() if i != alen):
                continue
        main_read_ids.append(read_index)
        filtered_alens.append(alen)

    if len(filtered_alens) < MIN_READS:
        locus_data.skip_code = 0
        return

    arr = np.array(filtered_alens)

    # ── histogram peak detection ──────────────────────────────────────
    bins         = np.arange(arr.min(), arr.max() + BIN_WIDTH, BIN_WIDTH)
    counts, edges = np.histogram(arr, bins=bins)
    min_count    = max(2, round(MIN_CLUSTER_FRAC * len(filtered_alens)))
    peaks, _     = find_peaks(counts, height=min_count, distance=2)

    if len(peaks) == 0:
        locus_data.skip_code = 0
        return

    # ── assign reads to nearest peak ─────────────────────────────────
    peak_centres = [(edges[p] + edges[p + 1]) / 2 for p in peaks]
    top2_centres = sorted(peak_centres,
                          key=lambda c: counts[np.searchsorted(edges, c) - 1],
                          reverse=True)[:2]

    def nearest_peak(alen):
        return min(range(len(top2_centres)),
                   key=lambda i: abs(alen - top2_centres[i]))

    c1_idx = [i for i, a in enumerate(filtered_alens) if nearest_peak(a) == 0]
    c2_idx = [i for i, a in enumerate(filtered_alens) if nearest_peak(a) == 1]

    c1_lengths    = [filtered_alens[i] for i in c1_idx]
    c2_lengths    = [filtered_alens[i] for i in c2_idx]
    hap_read_sets = (
        [main_read_ids[i] for i in c1_idx],
        [main_read_ids[i] for i in c2_idx]
    )

    locus_data.gt_depth  = len(filtered_alens)
    min_cluster_size     = MIN_CLUSTER_FRAC * len(filtered_alens)

    if c1_idx and c2_idx:
        if len(c1_idx) < min_cluster_size <= len(c2_idx):
            min_cluster_size = compute_cluster_cutoff(c1_lengths, c2_lengths)
        elif len(c2_idx) < min_cluster_size <= len(c1_idx):
            min_cluster_size = compute_cluster_cutoff(c2_lengths, c1_lengths)

    _assign_genotype(cooper, locus_key, locus_data,
                     c1_idx, c2_idx, c1_lengths, c2_lengths,
                     hap_read_sets, min_cluster_size)
