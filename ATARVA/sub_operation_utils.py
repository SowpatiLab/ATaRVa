import hdbscan
import numpy as np
import warnings
import statistics
import regex as re
import string

from ATARVA.consensus import *
from ATARVA.decompose import motif_decomposition

methviz_tag = False
def set_methviz_tag(value):
    global methviz_tag
    methviz_tag = value


# base64 encodings
BASE64_CHAR = string.ascii_uppercase + string.ascii_lowercase + string.digits + '+//'

def mm_tag_extract(read, mod_probs, meth_cutoff, frwd_strand):
    """
    record the positions with methylation calls
    
    :param read: pysam AlignedSegment object
    :param mod_probs: list of tuples with modified base positions and their respective probabilities
    :param meth_cutoff: score cutoff for calling methylation
    :param frwd_strand: boolean value for strand state; forward = True, reverse = False
    :return: none; the methylation positions are recorded in the read.methyl_range attribute
    """

    last_index = len(read.query_sequence) - 1
    if (read.methyl_start != None) and (read.methyl_end != None):
        for pos, prob in mod_probs:
            meth_chunk_start = pos     if frwd_strand else pos - 1 # to check the meth context, start index
            meth_chunk_end   = pos + 2 if frwd_strand else pos + 1 # to check the meth context, end index
            if read.methyl_start <= pos <= read.methyl_end:
                if (pos + 1 <= last_index) and (read.query_sequence[meth_chunk_start : meth_chunk_end]=='CG'):
                    read.methyl_range.append((pos, prob))


def pos_align(cg_positions, methyl_positions):
    
    cg_gaps = np.diff(cg_positions, prepend=0)
    methyl_gaps = np.diff(methyl_positions, prepend=0)
    remainder = 0
    crossed_read_pos = 0
    final_read_idx = [] # position index to take from read_pos
    
    methyl_npositions = len(methyl_gaps)
    longer_read_pos = methyl_positions >= len(cg_positions)

    for idx, cg_gap in enumerate(cg_gaps):
        tmp_len = len(final_read_idx)
        
        # Finding the start pos to begin aligning
        if longer_read_pos and (idx == 0) and (crossed_read_pos == 0): # when the read_pos are more in nu. compared to true pos
            tmp_diff = methyl_gaps[crossed_read_pos] - cg_gap

            if abs(tmp_diff) < 4: # Initial tagging of true & read meth pos can have a tolerance of 3bp on either direction
                final_read_idx.append(idx) # idx is zero here
                remainder = tmp_diff
                crossed_read_pos += 1
            else:
                jump_dist = cons_pos[idx] - read_pos[idx] # adjusting position by sliding the values from the initial position
                tmp_read_pos = [i - jump_dist for i in read_pos]

                distance_ins = [abs(cons_pos[idx] - i) for i in tmp_read_pos] # adjusted pos
                distance_nor = [abs(cons_pos[idx] - i) for i in read_pos] # normal pos
                min_dist_list = [sum(distance_nor), sum(distance_ins)]
                
                if min_dist_list.index(min(min_dist_list)) == 0: # checkpoint for finding optimal start position
                    crossed_read_pos = distance_nor.index(min(distance_nor))
                else:
                    crossed_read_pos = distance_ins.index(min(distance_ins))

                final_read_idx.append(crossed_read_pos)
                
                crossed_read_pos += 1
            continue
            
        check_pos = remainder
        check_diff = [100000] # list of diffs
        while crossed_read_pos < read_pos_len:
            check_pos += read_diff[crossed_read_pos] # incrementing to compare the diffs
            
            if abs(check_pos - true_diff) > min(check_diff): # choose the previous pos, if the current diff started to increase
                if crossed_read_pos-1 not in final_read_idx:
                    final_read_idx.append(crossed_read_pos-1)
                    remainder = (check_pos - read_diff[crossed_read_pos]) - true_diff
                break
            elif (crossed_read_pos+1 == read_pos_len) and abs(true_diff-read_diff[crossed_read_pos]) < 10: # for adding last pos
                if crossed_read_pos not in final_read_idx:
                    final_read_idx.append(crossed_read_pos)
                crossed_read_pos += 1
            else:
                check_diff.append(abs(check_pos - true_diff))
                crossed_read_pos += 1
                
        if tmp_len == len(final_read_idx):
            final_read_idx.append(-2)
            
    return final_read_idx


def encode_methylation(score_matrix, pos_matrix, cg_positions):
    
    encoded_meth = ''

    # Extracting only those positions which are within 2bp of true CG positions
    new_pos_matrix = []
    for methyl_positions in pos_matrix:
        new_pos_matrix.append(pos_align(cg_positions, methyl_positions))

    # Creating new matrix with only those positions
    new_matrix = []
    for idx,each_meth_cat in enumerate(matrix):
        pos_list = new_pos_matrix[idx]
        new_matrix.append([each_meth_cat[pos] if pos!=-2 else -2 for pos in pos_list])

    for col in zip(*new_matrix):
        col_array = np.array(col)
        mode = statistics.mode(col_array)
        if mode == -2:
            encryted_meth += '*' # adding * for skipping the positions where there is error call in few reads
        elif mode == -1:
            encryted_meth += '-' #adding - for ambiguous calls
        else:
            col_array = col_array[col_array != -1]
            col_array = col_array[col_array != -2]
            col_mean = round(np.mean(col_array), 2) * 100
            col_mean= round(col_mean/1.5625) # scaling to 0-64
            encryted_meth += BASE64_CHAR[col_mean]
    return encryted_meth


def methylation_calc(reads, read_methyl_scores, ALT_seq):

    arr  = np.frombuffer(ALT_seq.upper().encode(), dtype=np.uint8)
    c_mask = arr[:-1] == ord('C')
    g_mask = arr[1:]  == ord('G')
    cg_positions = np.where(c_mask & g_mask)[0]

    if cg_positions.size == 0:
        return [None, None, None]

    methyl_nreads = 0
    methyl_score   = 0
    encoded_string = None
    score_matrix = []
    pos_matrix = []
    for read_index in reads:
        if read_methyl_scores[read_index] is not None:
            methyl_nreads += 1
            methyl_score   += read_methyl_scores[read_index][0]
            methyl_scores  = read_methyl_scores[read_index][1]
            # for meth visualization
            if methviz_tag:
                score_matrix.append(methyl_scores)
                pos_matrix.append(read_methyl_scores[read_index][2])
            
    if methyl_nreads > 0:
        if methviz_tag:
            encrypted_meth = encode_methylation(score_matrix, pos_matrix, cg_positions)
        return [round(methyl_score/methyl_nreads, 2), methyl_nreads, encrypted_meth]
    else:
        return [None, None, None]


def dbscan(data, hap_reads):
    data = np.array(data).reshape(-1, 1)
    min_samples = max(10, round(0.2*len(data))) # min 20% of the data or 10 reads
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=FutureWarning)
        clusterer = hdbscan.HDBSCAN(min_cluster_size=min_samples)
        cluster_labels = clusterer.fit_predict(data)
    unique_labels = set(cluster_labels)
    
    if len(unique_labels)==1: # cluster case = (0), (-1)
        return [False,None,None] # proceed with Kmeans
        
    elif (len(unique_labels)==2) and (-1 in unique_labels): # cluster case = (0,-1)
        return [False,None,None] # proceed with Kmeans
        
    elif len(unique_labels)>=2: # cluster case = (0,1), (0,1,-1), (0,1,2)
        main_label = unique_labels-{-1}

        main_clusters = {}
        
        for label in main_label:
            c_label = [i for i, x in enumerate(cluster_labels) if x == label]
            alen = [data[i][0] for i in c_label]
            if len(c_label) in main_clusters:
                main_clusters[len(c_label)+1] = [c_label, alen]
            else:
                main_clusters[len(c_label)] = [c_label, alen]
            
        top2_clus_idx = [v for _,v in sorted(main_clusters.items(), reverse=True)[:2]] # getting top 2 cluster with more support

        new_haplotypes = [[hap_reads[idx] for idx in top2_clus_idx[0][0]], [hap_reads[idx] for idx in top2_clus_idx[1][0]]] # getting respective read ids

        new_alen = [top2_clus_idx[0][1], top2_clus_idx[1][1]]

        if set(new_alen[0])==set(new_alen[1]):
            return [False,None,None]
        
        return [True, new_haplotypes, new_alen]




    


def alt_sequence(read_seqs, hap_reads, amplicon, motif_size):
    seqs = [seq for seq in [read_seqs[read_id][0] for read_id in hap_reads] if seq!='']
    if len(seqs)>0:
        ALT = consensus_seq_poa(seqs)
        allele_length = len(ALT)
    else:
        ALT = '<DEL>'
        allele_length = 0

    decomp_seq = ''
    repeativity = True
    if amplicon and allele_length and (motif_size<=10):
        decomp_seq, nonrep_percent = motif_decomposition(ALT, motif_size)
        
        if nonrep_percent > 0.30: # if more than 30% of the sequence is non-repeat, repeativity = False
            repeativity = False
    return [ALT, allele_length, decomp_seq, repeativity]






