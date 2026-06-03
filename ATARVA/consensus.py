import pyabpoa as pa
import statistics as stats
from collections import Counter
from itertools import chain

def seq_organiser(seqs):
    seq_len_dict = {}
    for each_seq in seqs:
        lent = len(each_seq)
        if lent in seq_len_dict:
            seq_len_dict[lent].append(each_seq)
        else:
            seq_len_dict[lent] = [each_seq]
    sorted_mode_seqs = sorted(seq_len_dict.values(), key=lambda x : len(x), reverse=True)
    flat_list = list(chain.from_iterable(sorted_mode_seqs))
    del seq_len_dict, sorted_mode_seqs
    return flat_list

def seq_alter_median(seqs):
    odd = len(seqs)%2==1
    sorted_seqs = sorted(seqs)
    median_idx = len(sorted_seqs)//2
    upper_batch = sorted_seqs[median_idx:] # 1st half
    lower_batch = sorted_seqs[:median_idx][::-1] # reversed 2nd half
    insert_idx = []
    increment = 1
    for x in list(range(len(upper_batch))): # adjusted indices to insert the lower batch seqs alternatively
        insert_idx.append(x + increment)
        increment += 1
    low_idx = 0
    break_idx = insert_idx[-1] # Have to break the follwoing loop for even numbers of the total reads, because the lower batch will have one seq lesser than upper batch
    for idx in insert_idx:
        if odd and (break_idx == idx): break
        upper_batch.insert(idx, lower_batch[low_idx]) # Direct alternative insertion of lower batch seqs into upper batch
        low_idx += 1
    del sorted_seqs, lower_batch
    return upper_batch

def consensus_seq_poa(seqs):
    if len(seqs)<7:
        cons_algrm='MF'
    else:
        cons_algrm='HB'

    median_len = len(stats.median_high(seqs))
    if median_len > 10000:
        fseqs = [seq for seq in seqs if len(seq) == median_len]
        return ''.join(Counter(bases).most_common(1)[0][0] for bases in zip(*fseqs))

    sorted_seqs = seq_alter_median(seqs)
    abpoa = pa.msa_aligner(cons_algrm=cons_algrm)
    result = abpoa.msa(sorted_seqs, out_cons=True, out_msa=False)
    return result.cons_seq[0]