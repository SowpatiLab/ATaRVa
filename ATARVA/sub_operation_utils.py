import hdbscan
import numpy as np
import warnings
import statistics
import regex as re
from ATARVA.consensus import *
from ATARVA.decompose import motif_decomposition

methviz_tag = False
def set_methviz_tag(value):
    global methviz_tag
    methviz_tag = value

# base64 encodings
encode64_dict = {0: 'A', 1: 'B', 2: 'C', 3: 'D', 4: 'E', 5: 'F', 6: 'G', 7: 'H', 8: 'I', 9: 'J', 10: 'K', 11: 'L', 12: 'M', 13: 'N', 14: 'O', 15: 'P', 16: 'Q',17: 'R', 18: 'S', 19: 'T', 20: 'U', 21: 'V', 22: 'W',
                23: 'X', 24: 'Y', 25: 'Z', 26: 'a', 27: 'b', 28: 'c', 29: 'd', 30: 'e', 31: 'f', 32: 'g', 33: 'h', 34: 'i', 35: 'j', 36: 'k', 37: 'l', 38: 'm', 39: 'n', 40: 'o', 41: 'p', 42: 'q', 43: 'r', 44: 's',
                45: 't', 46: 'u', 47: 'v', 48: 'w', 49: 'x', 50: 'y', 51: 'z', 52: '0', 53: '1', 54: '2', 55: '3', 56: '4', 57: '5', 58: '6', 59: '7', 60: '8', 61: '9',62: '+', 63: '/', 64: '/'}

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


def mm_tag_extract(read, mod_scores, meth_cutoff, frwd_strand):
    """
    record the positions with methylation calls
    
    :param read: pysam AlignedSegment object
    :param mod_scores: list of tuples with modified base positions and their respective scores
    :param meth_cutoff: score cutoff for calling methylation
    :param frwd_strand: boolean value for strand state; forward = True, reverse = False
    :return: none; the methylation positions are recorded in the read.methyl_range attribute
    """

    last_index = len(read.query_sequence) - 1
    if (read.methyl_start != None) and (read.methyl_end!=None):
        for mod_val in mod_scores:
            position = mod_val[0]
            meth_chunk_start = position     if frwd_strand else position - 1 # to check the meth context, start index
            meth_chunk_end   = position + 2 if frwd_strand else position + 1 # to check the meth context, end index
            if read.methyl_start <= position <= read.methyl_end:
                if (position + 1 <= last_index) and (read.query_sequence[meth_chunk_start : meth_chunk_end]=='CG'):
                    read.methyl_range.append(position)

def cg_pos(seq):
    start_positions = [m.start() for m in re.finditer("CG", seq, overlapped=False)]
    return start_positions
            
def methylation_calc(hap_reads, global_loci_variations, locus_key, ALT_seq):
    if cg_pos(ALT_seq) == []:
        return [None, None, None]
    meth_reads = 0
    hap_meth = 0
    encrypted_meth = None
    locus_read_meth = global_loci_variations[locus_key]['read_meth']
    matrix = []
    pos_matrix = []
    for read_id in hap_reads:
        if locus_read_meth[read_id] is not None:
            meth_reads += 1
            hap_meth += locus_read_meth[read_id][0]
            meth_probs = locus_read_meth[read_id][1]
            # for meth visualization
            if methviz_tag:
                matrix.append(meth_probs)
                pos_matrix.append(locus_read_meth[read_id][2])
            
    if meth_reads > 0:
        if methviz_tag:
            encrypted_meth = methylation_encoding(matrix, pos_matrix, ALT_seq)
        return [round(hap_meth/meth_reads, 2), meth_reads, encrypted_meth]
    else:
        return [None, None, None]
    
def confidence_interval(data):
    data = np.array(data)
    ci = np.percentile(data, [2.5, 97.5])
    return [round(ci[0]), round(ci[1])]

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

def pos_diffs(pos_list):
    diffs = []
    for idx,pos in enumerate(pos_list):
        if idx == 0:
            diffs.append(pos)
            prev_val = pos
        else:
            diffs.append(pos - prev_val)
            prev_val = pos
    return diffs

def pos_align(cons_pos, cons_diff, read_pos):
    read_diff = pos_diffs(read_pos)
    remainder = 0
    crossed_read_pos = 0
    final_read_idx = [] # position index to take from read_pos
    read_pos_len = len(read_diff)
    longer_read_pos = read_pos_len >= len(cons_pos)

    for idx, true_diff in enumerate(cons_diff):    
        tmp_len = len(final_read_idx)
        # Finding the start pos to begin aligning
        if longer_read_pos and (idx == 0) and (crossed_read_pos == 0): # when the read_pos are more in nu. compared to true pos
            tmp_diff = read_diff[crossed_read_pos] - true_diff
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

def methylation_encoding(matrix, pos_matrix, ALT_seq):
    encryted_meth = ''
    cons_pos = cg_pos(ALT_seq) # getting the pos of CGs in the ALT sequence
    cons_diff = pos_diffs(cons_pos)
    # Extracting only those positions which are within 2bp of true CG positions
    new_pos_matrix = []
    for read_pos in pos_matrix:
        new_pos_matrix.append(pos_align(cons_pos, cons_diff, read_pos))

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
            encryted_meth += encode64_dict[col_mean]
    return encryted_meth
