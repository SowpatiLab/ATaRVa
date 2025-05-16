import pyabpoa as pa

def consensus_seq_poa(seqs):
    if len(seqs)<7:
        cons_algrm='MF'
    else:
        cons_algrm='HB'

    abpoa = pa.msa_aligner(cons_algrm=cons_algrm)
    result = abpoa.msa(seqs, out_cons=True, out_msa=False)
    return result.cons_seq[0]