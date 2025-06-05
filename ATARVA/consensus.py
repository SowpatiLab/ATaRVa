import pyabpoa as pa
import statistics as stats
from collections import Counter

def consensus_seq_poa(seqs):
    if len(seqs)<7:
        cons_algrm='MF'
    else:
        cons_algrm='HB'

    median_len = len(stats.median_high(seqs))
    if median_len > 10000:
        fseqs = [seq for seq in seqs if len(seq) == median_len]
        return ''.join(Counter(bases).most_common(1)[0][0] for bases in zip(*fseqs))

    abpoa = pa.msa_aligner(cons_algrm=cons_algrm)
    result = abpoa.msa(seqs, out_cons=True, out_msa=False)
    return result.cons_seq[0]