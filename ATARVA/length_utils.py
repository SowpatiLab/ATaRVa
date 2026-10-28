import warnings
import numpy as np
import stringzilla as sz

from sklearn.cluster   import KMeans
from sklearn.mixture   import GaussianMixture
from sklearn.neighbors import KernelDensity
from scipy.signal      import peak_widths
from threadpoolctl     import threadpool_limits
from scipy.signal      import find_peaks
from hdbscan           import HDBSCAN

from ATARVA.vcf_writer import *
from ATARVA.sub_operation_utils import alt_sequence, calculate_methylation
from ATARVA.consensus import consensus_seq_poa


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

    if len (c1_idx) + len(c2_idx) >= cooper.args.phased_frac * locus_data.depth \
        and len(c1_idx) >= cooper.args.cluster_frac * locus_data.depth \
        and len(c2_idx) >= cooper.args.cluster_frac * locus_data.depth:

        if c1_valid and c2_valid:
            locus_data.hap_read_sets = hap_read_sets
            locus_data.hap_alen_sets = (c1_lengths, c2_lengths)
            locus_data.is_genotyped = True
            heterozygous_call(cooper, locus_key)
            return

        if c1_valid:
            locus_data.hap_read_sets = (hap_read_sets[0], hap_read_sets[0])
            locus_data.hap_alen_sets = (c1_lengths, c1_lengths)
            locus_data.is_genotyped = True
            homozygous_call(cooper, locus_key)
            return

        if c2_valid:
            locus_data.hap_read_sets = (hap_read_sets[1], hap_read_sets[1])
            locus_data.hap_alen_sets = (c2_lengths, c2_lengths)
            locus_data.is_genotyped = True
            homozygous_call(cooper, locus_key)
            return

    locus_data.hap_read_sets = (locus_data.reads, [])
    locus_data.hap_alen_sets = ([locus_data.read_alens[rid][0] for rid in locus_data.read_alens], [])
    locus_data.is_genotyped = True
    homozygous_call(cooper, locus_key)

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
    ucluster     = [str(a) for a in sorted([x for x in hap_lengths if x < lower or x > upper])]

    ALT, allele_length = alt_sequence(locus_data.read_aseqs, hap_reads, cooper.args.amplicon)

    locus_data.gt_aseqs        = (ALT, None)
    locus_data.gt_alens        = (allele_length, None)
    locus_data.hap_meth_data   = (calculate_methylation(hap_reads, locus_data.read_methylation, ALT), None)
    locus_data.gt_arange       = (f'{lower}-{upper}', None)
    locus_data.gt_ucluster     = ('-'.join(ucluster), None)
    if cooper.args.decompose and ALT != '<DEL>':
        decomp_seq, nonrep_fraction = motif_decomposition(ALT, locus.motif_length)
        locus_data.gt_decomp_seqs  = (decomp_seq, None)

    ucluster = [len(locus_data.read_aseqs[rid][0]) for rid in locus_data.reads if rid not in hap_reads]
    if ucluster:
        locus_data.gt_ucluster = (locus_data.gt_ucluster[0], locus_data.gt_ucluster[1], '-'.join([str(x) for x in sorted(ucluster)]))

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
    phased_reads = set()

    for i in range(2):
        hap_reads = hap_read_sets[i]
        phased_reads |= set(hap_reads)
        hap_lengths = hap_alen_sets[i]
        ALT, allele_length = alt_sequence(locus_data.read_aseqs, hap_reads, cooper.args.amplicon)
        lower, upper = (round(x) for x in np.percentile(np.array(hap_lengths), [2.5, 97.5]))
        ucluster = [str(a) for a in sorted([x for x in hap_lengths if x < lower or x > upper])]
        if i == 0:
            locus_data.gt_aseqs         = (ALT, locus_data.gt_aseqs[1])
            locus_data.gt_alens         = (allele_length, locus_data.gt_alens[1])
            locus_data.hap_meth_data    = (calculate_methylation(hap_reads, locus_data.read_methylation, ALT), locus_data.hap_meth_data[1])
            locus_data.gt_arange        = (f'{lower}-{upper}', locus_data.gt_arange[1])
            locus_data.gt_ucluster      = ('-'.join(ucluster), locus_data.gt_ucluster[1])
            if cooper.args.decompose and ALT != '<DEL>':
                decomp_seq, nonrep_fraction = motif_decomposition(ALT, locus.motif_length)
                locus_data.gt_decomp_seqs   = (decomp_seq, locus_data.gt_decomp_seqs[1])
        else:
            locus_data.gt_aseqs         = (locus_data.gt_aseqs[0], ALT)
            locus_data.gt_alens         = (locus_data.gt_alens[0], allele_length)
            locus_data.hap_meth_data    = (locus_data.hap_meth_data[0], calculate_methylation(hap_reads, locus_data.read_methylation, ALT))
            locus_data.gt_arange        = (locus_data.gt_arange[0], f'{lower}-{upper}')
            locus_data.gt_ucluster      = (locus_data.gt_ucluster[0], '-'.join(ucluster))
            if cooper.args.decompose and ALT != 'DEL':
                decomp_seq, nonrep_fraction = motif_decomposition(ALT, locus.motif_length)
                locus_data.gt_decomp_seqs   = (locus_data.gt_decomp_seqs[0], decomp_seq)

    ucluster = [len(locus_data.read_aseqs[rid][0]) for rid in locus_data.reads if rid not in phased_reads]
    if ucluster:
        locus_data.gt_ucluster = (locus_data.gt_ucluster[0], locus_data.gt_ucluster[1], '-'.join([str(x) for x in sorted(ucluster)]))

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
    unique_alens    = set([l[0] for l in read_alens.values()])
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
    unique_alens    = set([l[0] for l in read_alens.values()])
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


def collapse_condensed_tree(clusterer, n_clusters=2):
    """
    Collapse HDBSCAN condensed-tree clusters upward until n_clusters
    remain.

    Returns
    -------
    cluster_labels : np.ndarray
        Cluster label for every original observation.
        Noise = -1.
    clusters : list
        Final condensed-tree cluster IDs.
    """

    tree = clusterer.condensed_tree_.to_numpy()

    parent = tree["parent"]
    child = tree["child"]
    child_size = tree["child_size"]

    # Cluster IDs are the nodes that have children.
    cluster_nodes = set(parent)

    # Start with the clusters selected by HDBSCAN's normal flat clustering.
    selected = list(clusterer.condensed_tree_._select_clusters())

    # If already <= requested number
    if len(selected) <= n_clusters:
        clusters = selected
    else:
        clusters = selected.copy()

        # Parent -> children mapping
        children = {}
        for p, c in zip(parent, child):
            children.setdefault(p, []).append(c)

        # Repeatedly collapse sibling clusters into their parent.
        while len(clusters) > n_clusters:

            cluster_set = set(clusters)

            # Find parents whose children contain >= 2 currently
            # represented clusters.
            candidates = []

            for p, cs in children.items():
                represented = cluster_set.intersection(cs)

                if len(represented) >= 2:
                    # Number of clusters that would disappear after
                    # replacing represented children by this parent.
                    reduction = len(represented) - 1

                    candidates.append(
                        (reduction, p, represented)
                    )

            if not candidates:
                break

            # Prefer a collapse that gets us closest to n_clusters.
            candidates.sort(
                key=lambda x: (
                    abs((len(clusters) - x[0]) - n_clusters),
                    -x[0]
                )
            )

            _, p, represented = candidates[0]

            clusters = [
                c for c in clusters
                if c not in represented
            ]

            clusters.append(p)

    # ---------------------------------------------------------
    # Convert final condensed-tree clusters to observation labels
    # ---------------------------------------------------------

    # Build child -> parent mapping
    child_to_parent = dict(zip(child, parent))

    # Points are leaves in the condensed tree.
    final_clusters = set(clusters)

    labels = np.full(clusterer.labels_.shape, -1, dtype=int)

    # Map each observation to its first ancestor belonging to
    # one of the final clusters.
    for point in range(len(labels)):

        node = point

        while node not in final_clusters:

            if node not in child_to_parent:
                break

            node = child_to_parent[node]

        if node in final_clusters:
            labels[point] = clusters.index(node)

    return labels, clusters


def length_genotyper_hdbscan(cooper, locus_key):
    """
    Genotype using HDBSCAN — auto cluster count, outlier aware.

    :param cooper:     cooper object
    :param locus_key:  key for the locus
    """

    MIN_READS        = 3
    MIN_CLUSTER_FRAC = cooper.args.cluster_frac
    WINDOW_FRAC      = 0.1

    locus_data = cooper.cooper_loci_data[locus_key]
    locus      = cooper.cooper_loci_info[locus_key]
    read_alens = locus_data.read_alens

    read_indices    = sorted(locus_data.reads)
    unique_alens    = set([l[0] for l in read_alens.values()])
    singleton_alens = {alen for alen, count in locus_data.halen_frequency.items() if count == 1}
    windows         = {i: (round(i * (1 - WINDOW_FRAC)), round(i * (1 + WINDOW_FRAC)))
                       for i in unique_alens}

    main_read_ids  = []
    filtered_alens = []
    filtered_seqs  = []
    for read_index in read_indices:
        alen = read_alens[read_index][0]
        if alen in singleton_alens:
            if not any(lo <= alen <= hi for i, (lo, hi) in windows.items() if i != alen):
                continue
        main_read_ids.append(read_index)
        filtered_alens.append(alen)
        filtered_seqs.append(locus_data.read_aseqs[read_index][0])

    if len(filtered_alens) < MIN_READS:
        locus_data.skip_code = 0
        return

    # ── compute edit distance from reference allele ──────────────────
    reference_seq  = cooper.ref.fetch(locus.chrom, locus.start, locus.end)  # use reference sequence as reference
    edit_distances = [sz.edit_distance(reference_seq, seq) for seq in filtered_seqs]

    # ── normalize features for clustering ──────────────────────────────
    dist_normalized = np.array(edit_distances)
    feature_array = np.column_stack([filtered_alens, dist_normalized])

    # ── HDBSCAN clustering with 2D features ─────────────────────────────
    clusterer = HDBSCAN(
        min_cluster_size     = max(MIN_READS, int(MIN_CLUSTER_FRAC * len(filtered_alens))),
        allow_single_cluster = True
    ).fit(feature_array)

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

    # if len(clusters) > 2:
    #     collapsed = collapse_condensed_tree(clusterer, n_clusters=2)
    #     clusters = {}
    #     for i, label in enumerate(collapsed[0]):
    #         if label not in clusters:
    #             clusters[label] = []
    #         clusters[label].append(i)

    # ── take two largest clusters ─────────────────────────────────────
    top2      = sorted(clusters, key=lambda l: len(clusters[l]), reverse=True)[:2]
    c1_idx    = clusters[top2[0]]
    c2_idx    = clusters[top2[1]] if len(top2) > 1 else []

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

    _assign_genotype(cooper, locus_key, locus_data, c1_idx, c2_idx, 
                     c1_lengths, c2_lengths, hap_read_sets, min_cluster_size)


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
    unique_alens    = set([l[0] for l in read_alens.values()])
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


def score_calc(x_grid, density, initial_peaks, valleys, top_contour_widths):
    """
    Calculate the score for each peak based on prominence, width, and skewness.

    :param x_grid:             grid of x values for density estimation
    :param density:            density values corresponding to x_grid
    :param initial_peaks:      indices of detected peaks in the density
    :param valleys:            indices of valleys between peaks in the density
    :param top_contour_widths: widths of the peaks at half prominence
    :return:                    final_score, initial_prominence, area_covered
    """

    peak_density   = density[initial_peaks]
    valley_density = density[valleys]
    peak_points    = x_grid[initial_peaks]
    initial_score      = []
    initial_prominence = []
    initial_skewness   = []
    area_covered       = []
    for idx in range(len(peak_density)):
        if idx == 0:
            f_dense   = density[0]
            base_left = x_grid[0][0]

        if idx < len(peak_density)-1:
            valley_dense = valley_density[idx]
            base_right   = x_grid[valleys[idx]][0]
        else:
            valley_dense = density[-1]
            base_right   = x_grid[-1][0]

        max_point  = max(f_dense, valley_dense)
        prominence = (peak_density[idx] - max_point)# - diffs
        f_dense    = valley_dense

        flattened_x_grid = x_grid[:, 0]
        mask   = (flattened_x_grid >= base_left) & (flattened_x_grid <= base_right)
        x_vals = flattened_x_grid[mask]
        y_vals = density[mask]

        area = np.trapz(y_vals, x_vals)
        area_covered.append(area)
        min_area = 1 if area >= 0.02 else 0 # min area covered by the peak should be atleast 2% to be considered as valid peak

        current_peak = peak_points[idx][0]
        L = (current_peak - base_left); R = (base_right-current_peak) # distance between left boundary to the peak and right boundary to the peak

        eps = 1e-8
        # normalized asymmetry/skewness
        K = abs(L - R) / (L + R + eps)
        initial_skewness.append(K)

        # final score
        score = ((prominence ** 2) / ((top_contour_widths[idx] + eps) ** 2)) * min_area

        initial_score.append(score)
        initial_prominence.append(prominence)

        base_left = base_right

    initial_skewness = np.array(initial_skewness)
    median_K         = np.median(initial_skewness)
    mad_K            = np.median( np.abs( initial_skewness - median_K ) ) + 1e-8
    zK               = (initial_skewness - median_K) / mad_K
    quality_zk       = 1 / (1 + np.exp(zK))

    final_score = ( np.array(initial_score) * quality_zk )

    return final_score, initial_prominence, area_covered


def length_genotyper_kde(cooper, locus_key):
    """
    Genotype using Kernel Density Estimation (KDE) for peak detection.

    :param cooper:     cooper object
    :param locus_key:  key for the locus
    """

    MIN_READS        = 3
    MIN_CLUSTER_FRAC = 0.15
    WINDOW_FRAC      = 0.1
    BIN_WIDTH        = 5

    locus_data = cooper.cooper_loci_data[locus_key]
    locus      = cooper.cooper_loci_info[locus_key]
    read_alens = locus_data.read_alens

    read_indices    = sorted(locus_data.reads)
    unique_alens    = set([l[0] for l in read_alens.values()])
    singleton_alens = {alen for alen, count in locus_data.halen_frequency.items() if count == 1}
    windows         = {i: (round(i * (1 - WINDOW_FRAC)), round(i * (1 + WINDOW_FRAC)))
                        for i in unique_alens}

    main_read_ids  = []
    filtered_alens = []
    filtered_seqs  = []
    for read_index in read_indices:
        alen = read_alens[read_index][0]
        if alen in singleton_alens:
            if not any(lo <= alen <= hi for i, (lo, hi) in windows.items() if i != alen):
                continue
        main_read_ids.append(read_index)
        filtered_alens.append(alen)
        filtered_seqs.append(locus_data.read_aseqs[read_index][0])

    alen_units = np.array([length//locus.motif_length for length in filtered_alens])
    #### KDE with mode peaks and valley for definitive split point for each peak based on the area under the peaks
    bandwidth = 10; tot_data_points = 1000 # for amplicon, to get better density estimation and peaks
    stdev = np.std(alen_units)
    if stdev != 0:
        bandwidth = 0.5 * stdev * (len(alen_units) ** (-1/5))

    alen_units = alen_units.reshape(-1, 1)

    # Fit kde to the data
    kde = KernelDensity(kernel='gaussian', algorithm='kd_tree', metric='minkowski', bandwidth=bandwidth).fit(alen_units)
    # Evaluate the density on a grid
    x_grid      = np.linspace(alen_units.min()-50, alen_units.max()+50, tot_data_points).reshape(-1, 1)
    log_density = kde.score_samples(x_grid)
    density     = np.exp(log_density)

    # Analysing the distribution to identify the sharp narrow peaks
    initial_peaks, _   = find_peaks(density)
    original_widths    = peak_widths(density, initial_peaks)
    top_contour_widths = peak_widths(density, initial_peaks, rel_height=0.2)
    valleys, _         = find_peaks(-density)

    score, initial_prominence, area_covered = score_calc(x_grid, density, initial_peaks, valleys, top_contour_widths[0])
    narrow_peaks_idx = list(np.argsort(score)[-2:]) # taking top two peaks with more area under the curve, as the peaks with higher area will be sharper and more prominent

    top_2_prominence   = [initial_prominence[idx] for idx in narrow_peaks_idx]
    min_height_covered = min(top_2_prominence) >= 0.15 * (max(top_2_prominence)) # min peak should have atleast 15% of the max peak prominence to be considered as a valid peak
    min_area_covered   = all([area_covered[idx]>=0.05 for idx in narrow_peaks_idx]) # both peaks should have atleast 5% of the area covered to be considered as valid peaks

    if min_height_covered or min_area_covered: # any of this should be True to consider both peaks as valid peaks, otherwise only the max peak will be considered for split
        pass
    elif (narrow_peaks_idx[0] > narrow_peaks_idx[1]) and (area_covered[narrow_peaks_idx[0]] >= 0.02): # if the min_peak is on right side of the max_peak and if it has atleast 2% of total area consider both peaks as valid; to report the longer allele
        pass
    else: # if the min_peak is left side of the max_peak, then consider only the max_peak as valid and as homozygous 
        narrow_peaks_idx = [narrow_peaks_idx[1]]

    width     = sorted(original_widths[0][narrow_peaks_idx])
    top_count = len(narrow_peaks_idx)

    # Getting new peaks with analysed data
    peaks, _  = find_peaks(density, width = width)

    # Choose split 
    peak_heights = density[peaks] # extracting only the peaks frim density
    top_peaks    = peaks[np.argsort(peak_heights)[-top_count:]] # taking top two peaks
    sorted_peaks = sorted(top_peaks)

    # flattening the data and initializing the labels for each data point as -1 (unassigned)
    alen_units = alen_units.flatten()
    labels     = np.full(len(alen_units), -1)

    left = sorted_peaks[0]
    # boundaries of left peak
    peak1_left  = valleys[(valleys < left)]
    p1_start    = peak1_left[-1] if peak1_left.size > 0 else 0
    peak1_right = valleys[(valleys > left)]
    p1_end      = peak1_right[0] if peak1_right.size > 0 else len(x_grid) - 1

    p1_left_split = x_grid[p1_start][0]
    p1_right_split = x_grid[p1_end][0]

    p1_allele_bool = (alen_units >= p1_left_split) & (alen_units <= p1_right_split)
    if p1_allele_bool.sum() > 4: # atleast 5 reads should be there in cluster to consider it as valid cluster
        labels[p1_allele_bool] = 0

    if len(sorted_peaks) > 1:
        right = sorted_peaks[1]
        # boundaries of right peak
        peak2_left  = valleys[(valleys < right)]
        p2_start    = peak2_left[-1] if peak2_left.size > 0 else 0
        peak2_right = valleys[(valleys > right)]
        p2_end      = peak2_right[0] if peak2_right.size > 0 else len(x_grid) - 1

        p2_left_split  = x_grid[p2_start][0]
        p2_right_split = x_grid[p2_end][0]

        labels[(alen_units >= p2_left_split) & (alen_units <= p2_right_split)] = 1 # no min read cutoff for the longer allele

    c1_idx = [i for i, x in enumerate(labels) if x == 0]
    c2_idx = [i for i, x in enumerate(labels) if x == 1]

    c1_lengths = [filtered_alens[i] for i in c1_idx]
    c2_lengths = [filtered_alens[i] for i in c2_idx]

    cutoff = MIN_CLUSTER_FRAC * len(filtered_alens)
    br     = False

    if c1_lengths and c2_lengths:
        def confidence_interval(data):
            data = np.array(data)
            ci   = np.percentile(data, [2.5, 97.5])
            return [round(ci[0]), round(ci[1])]

        def process_conditions(alen_x, alen_y):
            nonlocal br, cutoff
            min_val, max_val = confidence_interval(alen_y)
            min_slide = 10 if min_val * 0.1 > 10 else min_val * 0.1  # max(min_val*0.1, 10)#; maximum length slide is 10
            max_slide = 10 if max_val * 0.1 > 10 else max_val * 0.1  # max(max_val*0.1, 10)#; maximum length slide is 10
            min_bound = min_val - min_slide
            max_bound = max_val + max_slide

            for min_al in alen_x:
                if min_bound <= min_al <= max_bound:
                    br = True; break

            if not br:
                cutoff = min(5, cutoff)

        if len(c1_lengths) < cutoff and len(c2_lengths) >= cutoff:
            process_conditions(c1_lengths, c2_lengths)
                               
        elif len(c2_lengths) < cutoff and len(c1_lengths) >= cutoff:
            process_conditions(c2_lengths, c1_lengths)

    cutoff = int(cutoff)

    if cooper.haploid:
        cluster_len = [len(c1_lengths), len(c2_lengths)]
        cidx        = cluster_len.index(max(cluster_len))
        if cluster_len[cidx] >= cutoff:
            major_alens = c1_idx if cidx == 0 else c2_idx # major allele cluster
            major_reads = c1_lengths if cidx == 0 else c2_lengths # major allele lengths
            locus_data.hap_read_sets = (major_reads, [])
            locus_data.hap_alen_sets = (major_alens, None)

        else:
            locus_data.hap_read_sets = (locus_data.reads, [])
            locus_data.hap_alen_sets = ([locus_data.read_alens[rid][0] for rid in locus_data.read_alens], [])
            locus_data.is_genotyped = True
            homozygous_call(cooper, locus_key)
        
    elif len(c1_lengths) >= cutoff and len(c2_lengths) >= cutoff:
        locus_data.hap_read_sets = ( [main_read_ids[i] for i in c1_idx], [main_read_ids[i] for i in c2_idx] )
        locus_data.hap_alen_sets = (c1_lengths, c2_lengths)
        locus_data.is_genotyped  = True
        heterozygous_call(cooper, locus_key)
    
    elif len(c1_lengths) >= cutoff:
        locus_data.hap_read_sets = ([main_read_ids[i] for i in c1_idx], [main_read_ids[i] for i in c1_idx])
        locus_data.hap_alen_sets = (c1_lengths, c1_lengths)
        locus_data.is_genotyped = True
        homozygous_call(cooper, locus_key)

    elif len(c2_lengths) >= cutoff:
        locus_data.hap_read_sets = ([main_read_ids[i] for i in c2_idx], [main_read_ids[i] for i in c2_idx])
        locus_data.hap_alen_sets = (c2_lengths, c2_lengths)
        locus_data.is_genotyped = True
        homozygous_call(cooper, locus_key)

    # _assign_genotype(cooper, locus_key, locus_data,
    #                  c1_idx, c2_idx, c1_lengths, c2_lengths,
    #                  ( [main_read_ids[i] for i in c1_idx], [main_read_ids[i] for i in c2_idx] ),
    #                  min_cluster_size = 1) # no min cluster size cutoff for the longer allele
