import numpy as np
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, dendrogram
from ATARVA.sub_operation_utils import *

seq_dict = {'A': [0, 0], 'C': [0, 1], 'G': [1, 0], 'T': [1, 1], 'N': [1, 1]}

def correlation_clustering(read_seqs, read_indices, motif_size, global_loci_variations, locus_key):

    seq_list = []
    for read_id in read_indices:
        seq_list.append(read_seqs[read_id][0])
    # Convert sequence to bit representation using NumPy
    # if variance is zero, paiwise distance will cause error
    seq_len = [len(i) for i in seq_list]
    bit_seqs = [np.concatenate([seq_dict.get(i, [1,1]) for i in seq]) if seq!='' else np.array([1,1]) for seq in seq_list]
    # the array shouldn't contain same numbers eg [1,1,1,1], it should not have variance of '0' i.e why  'else np.array([1,1])' this is added, if the allele is a full deletion
    # then padding will -1-1-1-1-1-1 giving 0 variance

    max_seq_len = max(len(seq) for seq in bit_seqs) + 4
    # if the lengthiest sequence is a monomer eg AAAAAAAAAAAAAAAAAAAAAAA, its variance is zero
    # thats why extra '4' (any number of padding can be added) is added to the max length, so variance wont be zero

    bit_matrix = np.array([np.pad(seq, (0, max_seq_len-len(seq)), constant_values = -1) for seq in bit_seqs], dtype=np.int8)
    dist_matrix = pdist(bit_matrix, metric='correlation') # correlation formula is used for distance calculation
    linkage_matrix = linkage(dist_matrix, method='complete') # complete method is used to get compact clusters

    # 3rd col of Linkage has the distance values and 90th percentile is calculated to cut the dendrogram to get tight unique clusters
    dist_percentile = np.percentile(linkage_matrix[:, 2], 90)

    # Generating dendrogram details without plotting
    den_detail = dendrogram(linkage_matrix, count_sort='descending', color_threshold=dist_percentile, no_plot=True)

    # Grouping reads based on cluster colors
    color = den_detail['leaves_color_list']
    leaves = den_detail['leaves']
    cluster_dict = {}
    for idx, each_color in enumerate(color):
        if each_color == 'C0': continue # skipping the outlier cluster
        if each_color not in cluster_dict:
            cluster_dict[each_color] = [read_indices[leaves[idx]]] # getting the read_id from leave_index of dendrogram
        else:
            cluster_dict[each_color].append(read_indices[leaves[idx]])

    if len(cluster_dict)==0: # no valid clusters found
        return [False, 6, {}]
    
    haplotypes = list(cluster_dict.values())

    alt_seq_lens = set()
    alt_seqs = set()

    genotype_dict = {}
    for hap_reads in haplotypes:

        if len(hap_reads)<3: continue

        ALT, allele_length, decomp_seq, repeativity = alt_sequence(read_seqs, hap_reads, True, motif_size) # true for amplicon, to check the repetitiveness in the sequence


        # skip if the cluster has less reads and is not a repeat
        if not repeativity:
            continue

        if allele_length in alt_seq_lens:
            if ALT not in alt_seqs:
                alt_seqs.add(ALT)
            else:
                continue
        else:
            alt_seq_lens.add(allele_length)
            alt_seqs.add(ALT) 

        ci = confidence_interval([len(read_seqs[read_id][0]) for read_id in hap_reads])
        meth_info = methylation_calc(hap_reads, global_loci_variations, locus_key, ALT)
        if allele_length not in genotype_dict:
            genotype_dict[allele_length] = (ALT, ci, decomp_seq, len(hap_reads), meth_info)
        else:
            genotype_dict[str(allele_length)] = (ALT, ci, decomp_seq, len(hap_reads), meth_info)

    del read_seqs, alt_seqs, alt_seq_lens, haplotypes, cluster_dict, bit_matrix, dist_matrix, linkage_matrix, den_detail
    return [True, 10, genotype_dict]
