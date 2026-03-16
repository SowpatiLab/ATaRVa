import numpy as np
import warnings
from sklearn.cluster import KMeans
from threadpoolctl import threadpool_limits

from ATARVA.snp_utils import haplocluster_reads
from ATARVA.vcf_writer import *
from ATARVA.sub_operation_utils import alt_sequence
from ATARVA.somatic_utils import *


def homozygous_call(cooper, locus_key, allele_lengths, hallele_counter, read_indices):

    locus = cooper.cooper_loci_info[locus_key]
    locus_data = cooper.cooper_loci_data[locus_key]
    lower, upper = [round(x) for x in np.percentile(allele_lengths, [2.5, 97.5])]
    allele_range = f'{lower}-{upper},{lower}-{upper}'
    ALT, allele_length, decomposed_seq, is_repetitive = alt_sequence(locus_data.read_seqs, read_indices, locus.motif_length)
    if is_repetitive:
        methylation_data = calculate_methylation(read_indices, locus_data.read_methylation, ALT)
        write_homozygous_call(cooper, locus_key, allele_length, hallele_counter, allele_range, methylation_data, decomposed_seq, ALT)
    else:
        return [False, 6]
    return [True, 10]


def heterozygous_call(cooper, haplotypes, read_seqs, new_alen, locus_key, read_indices, hallele_counter):

    alen_c1 = new_alen[0]
    alen_c2 = new_alen[1]
    phased_read = ['.','.']
    chosen_snpQ = '.'
    snp_num = '.'        

    genotypes = []
    allele_count = {}
    ALT_seqs = []
    repeativity_list = []
    decomp_seq_list = []
    meth_info = []
    for hap_reads in haplotypes:
        ALT, allele_length, decomp_seq, repeativity = alt_sequence(read_seqs, hap_reads, motif_size)
        repeativity_list.append(repeativity)
        decomp_seq_list.append(decomp_seq)
        ALT_seqs.append(ALT)
        genotypes.append(allele_length)
        if allele_length not in allele_count:
            allele_count[allele_length] = len(hap_reads)
        else:
            allele_count[str(allele_length)] = len(hap_reads)

        meth_info.append(calculate_methylation(hap_reads, global_loci_variations, locus_key, ALT))

    lower1,upper1 = np.percentile(alen_c1, [2.5, 97.5])
    lower2,upper2 = np.percentile(alen_c2, [2.5, 97.5])
    allele_range = f'{lower1}-{upper1},{lower2}-{upper2}'

    if all(repeativity_list):
        vcf_heterozygous_writer(contig, genotypes, locus_start, locus_end, allele_count, len(read_indices), global_loci_info, ref, out, chosen_snpQ, phased_read, snp_num, ALT_seqs, log_bool, 'kmeans', decomp, hallele_counter, allele_range, decomp_seq_list, meth_info)
    elif any(repeativity_list):
        if repeativity_list[0]:
            allele_range = f'{lower1}-{upper1},{lower1}-{upper1}'
            write_homozygous_call(ref, contig, locus_key, global_loci_info, genotypes[0], len(haplotypes[0]), len(read_indices), out, ALT_seqs[0], log_bool, 'kmeans', decomp, hallele_counter, False, allele_range, decomp_seq_list[0], meth_info[0])
        else:
            allele_range = f'{lower2}-{upper2},{lower2}-{upper2}'
            write_homozygous_call(ref, contig, locus_key, global_loci_info, genotypes[1], len(haplotypes[1]), len(read_indices), out, ALT_seqs[1], log_bool, 'kmeans', decomp, hallele_counter, False, allele_range, decomp_seq_list[1], meth_info[1])
    else:
        return [False, 6]
    
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

    # if minor cluster overlaps with major — no cutoff needed
    overlaps = any(lower_bound <= alen <= upper_bound for alen in minor_cluster)
    if overlaps: return 0.15 * (len(major_cluster) + len(minor_cluster)) # 15% of total reads in both clusters

    # min 3% of major cluster size, at least 2 reads
    ratio_cutoff = int(max(0.03, len(minor_cluster) / len(major_cluster)) * len(major_cluster))
    return max(2, ratio_cutoff)

    
def length_genotyper(cooper, locus_key, hallele_counter, read_indices):
    """
    genotype a locus by clustering allele lengths using KMeans.

    :param cooper:          cooper object
    :param hallele_counter: allele length counter dict
    :param locus_key:       locus identifier string
    :param read_indices:    list of read indices covering the locus
    :return:                [bool_state, category]
    """
    FAIL             = [False, 6]
    MIN_READS        = 3
    MIN_CLUSTER_FRAC = 0.15
    WINDOW_FRAC      = 0.1

    locus      = cooper.cooper_loci_info[locus_key]
    locus_data = cooper.cooper_loci_data[locus_key]
    read_alens = locus_data.read_alens
    read_seqs  = locus_data.read_seqs

    read_indices  = sorted(read_indices)
    unique_alens  = set(hallele_counter)
    singleton_alens = {alen for alen, count in hallele_counter.items() if count == 1}

    # --- pre-compute 10% windows for each unique allele ---
    windows = { i: (round(i * (1 - WINDOW_FRAC)), round(i * (1 + WINDOW_FRAC))) for i in unique_alens }

    # --- filter singleton alleles not near any other allele ---
    main_read_ids  = []
    filtered_alens = []

    for read_index in read_indices:
        alen = read_alens[read_index][0]
        if alen in singleton_alens:
            near_other = any(lo <= alen <= hi for i, (lo, hi) in windows.items() if i != alen)
            if not near_other: continue
        main_read_ids.append(read_index)
        filtered_alens.append(alen)

    if len(filtered_alens) < MIN_READS: return FAIL

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

    alen_c1    = [filtered_alens[i] for i in c1_idx]
    alen_c2    = [filtered_alens[i] for i in c2_idx]
    haplotypes = ([main_read_ids[i] for i in c1_idx], [main_read_ids[i] for i in c2_idx])

    # --- compute cluster cutoff ---
    min_cluster_size = MIN_CLUSTER_FRAC * len(filtered_alens)

    if c1_idx and c2_idx:
        if len(c1_idx) < min_cluster_size <= len(c2_idx):
            min_cluster_size = compute_cluster_cutoff(alen_c1, alen_c2)
        elif len(c2_idx) < min_cluster_size <= len(c1_idx):
            min_cluster_size = compute_cluster_cutoff(alen_c2, alen_c1)

    n_total = len(read_indices)

    # --- haploid genotyping ---
    if cooper.haploid:
        major_idx   = 0 if len(c1_idx) >= len(c2_idx) else 1
        major_alens = alen_c1    if major_idx == 0 else alen_c2
        major_reads = haplotypes[major_idx]

        if len(major_reads) < min_cluster_size: return FAIL

        lower, upper  = np.percentile(major_alens, [2.5, 97.5])
        allele_range  = f'{lower}-{upper}'
        ALT, allele_length, decomposed_seq, is_repetitive = alt_sequence(read_seqs, major_reads, locus.motif_length)
        methylation_data = calculate_methylation(major_reads, locus_data.read_methylation, ALT)

        if not is_repetitive: return FAIL
        # def homozygous_call(cooper, locus_key, allele_lengths, hallele_counter, read_indices):
        homozygous_call(cooper, locus_key, major_alens, hallele_counter, major_reads)
        return [True, 10]

    # ── diploid genotyping ────────────────────────────────────────────
    c1_valid = bool(c1_idx) and len(c1_idx) >= min_cluster_size
    c2_valid = bool(c2_idx) and len(c2_idx) >= min_cluster_size

    if c1_valid and c2_valid:
        # bool_state, category = heterozygous_call(haplotypes, read_seqs, [alen_c1, alen_c2], contig, locus_key, read_indices, hallele_counter)
        return [bool_state, category]

    if c1_valid:
        bool_state, category = homozygous_call(cooper, locus_key, alen_c1, hallele_counter, haplotypes[0])
        return [bool_state, category]

    if c2_valid:
        bool_state, category = homozygous_call(cooper, locus_key, alen_c2, hallele_counter, haplotypes[1])
        return [bool_state, category]

    return FAIL


def analyse_genotype(cooper, locus_key, hallele_counter, read_indices):
    """
    genotype the locus based on the read data

    :param cooper: Cooper object
    :param locus_key: key for the locus
    :param hallele_counter: counter for haplotype alleles
    :param read_indices: indices of reads to consider
    :return: list of genotype status and category
    """

    locus = cooper.cooper_loci_info[locus_key]
    status = False

    read_seqs = cooper.cooper_loci_data[locus_key].read_seqs

    if cooper.somatic: # for somatic variant calling
        # state, skip_point, genotype_dict = correlation_clustering(read_seqs, read_indices, motif_size, global_loci_variations, locus_key)
        # if state:
        #     vcf_multizygous_writer(contig, genotype_dict, locus_start, locus_end, len(read_indices), global_loci_info, ref, out, log_bool, decomp, hallele_counter)
        return [state, skip_point]
    
    elif cooper.haploid: # for haploid and amplicon genotyping
        state, skip_point = length_genotyper(cooper, locus_key, hallele_counter, read_indices)
        return [state, skip_point]

    snp_positions = set()
    for rindex in read_indices:
        snp_positions |= (cooper.cooper_read_data[rindex].snps)

    snp_positions = sorted(list(filter(lambda x: (x in cooper.cooper_snp_data) and (cooper.cooper_snp_data[x].cov >= 3) and
                                                    (locus.start - cooper.args.snp_dist < x < locus.end + cooper.args.snp_dist),
                            snp_positions)))

    snp_allelereads = {}
    read_indices = set(read_indices)
    non_ref_snp_cov = {}
    for pos in snp_positions:
        c_point = 0
        coverage = set()
        alt_nucs = [nuc for nuc in cooper.cooper_snp_data[pos].sub]
        for alt_nuc in alt_nucs:
            reads_of_nuc = cooper.cooper_snp_data[pos].sub[alt_nuc].intersection(read_indices)
            if len(reads_of_nuc) == 0: continue
            coverage.add(len(reads_of_nuc))

            if (sum([cooper.cooper_snp_data[pos].qual[read_idx] for read_idx in reads_of_nuc])/len(reads_of_nuc)) <= cooper.args.snp_qual:
                c_point=1
                break
        if (len(coverage)==0) or (c_point==1): continue
        else: non_ref_snp_cov[pos] = max(coverage)
            
        snp_allelereads[pos] = { 'cov': 0, 'reads': set(), 'alleles': {}, 'qual': {} }
        for nuc in cooper.cooper_snp_data[pos].sub:
            snp_allelereads[pos]['alleles'][nuc] = cooper.cooper_snp_data[pos].sub[nuc].intersection(read_indices)
            snp_allelereads[pos]['cov'] += len(snp_allelereads[pos]['alleles'][nuc])
            if nuc!='r':
                snp_allelereads[pos]['qual'].update(dict([(read_idx,cooper.cooper_snp_data[pos].qual[read_idx]) for read_idx in snp_allelereads[pos]['alleles'][nuc]]))

    del_positions = list(filter(lambda x: snp_allelereads[x]['cov'] < 5, snp_allelereads.keys()))

    for pos in del_positions: del snp_allelereads[pos]

    ordered_snp_on_cov = sorted(snp_allelereads.keys(), key = lambda item : non_ref_snp_cov[item], reverse = True)

    print(f'SNP positions considered for phasing: {snp_allelereads}')
    (haplotypes,
     min_snp, skip_point,
     chosen_snpQ,
     phased_read, snp_num) = haplocluster_reads(snp_allelereads, ordered_snp_on_cov, read_indices, cooper.args.snp_qual,
                                                cooper.args.snp_count, cooper.args.snp_read, cooper.args.phasing_read)

    if haplotypes == (): # if the loci has no significant snps
        state, skip_point = length_genotyper(cooper, locus_key, hallele_counter, read_indices)
        del read_seqs
        return [state, skip_point]
    
    if min_snp != -1:
        snp_left_boundary = locus.start - cooper.args.snp_dist
        min_idx = 0
        for each_spos in cooper.cooper_sorted_snps:
            if each_spos >= snp_left_boundary:
                break
            del cooper.cooper_snp_data[each_spos]
            min_idx += 1
        del cooper.cooper_sorted_snps[:min_idx]


    genotypes = []
    allele_count = {}
    ALT_seqs = []
    alen_list = []
    meth_info = []
    for hap_reads in haplotypes:
        ALT, allele_length,_,_ = alt_sequence(read_seqs, hap_reads, False, locus.motif_length) # false for genome, to not check the repetitiveness in the sequence
        alen_list.append([len(read_seqs[read_id][0]) for read_id in hap_reads])
        ALT_seqs.append(ALT)
        genotypes.append(allele_length)
        if allele_length not in allele_count:
            allele_count[allele_length] = len(hap_reads)
        else:
            allele_count[str(allele_length)] = len(hap_reads)

        meth_info.append(calculate_methylation(hap_reads, cooper.cooper_loci_data[locus_key].read_methylation, ALT))

    del read_seqs
    lower1, upper1 = np.percentile(alen_list[0], [2.5, 97.5])
    lower2, upper2 = np.percentile(alen_list[1], [2.5, 97.5])
    allele_range = f'{lower1}-{upper1},{lower2}-{upper2}'
    vcf_heterozygous_writer(contig, genotypes, locus_start, locus_end, allele_count, len(read_indices), global_loci_info, ref, out, chosen_snpQ, phased_read, snp_num, ALT_seqs, log_bool, 'SNP', decomp, hallele_counter, allele_range, [None], meth_info)
    state = True
    return [state, skip_point]
    