import sys
import statistics

def consensus_seq_poa(seqs, alen):

    tot_seq = len(seqs)
    seq_len = len(seqs[0])
    
    consensus = ''
    for idx in range(seq_len):
        current_base = []
        for each_seq in seqs:
            current_base.append(each_seq[idx])
        consensus+=statistics.mode(current_base)

    if len(consensus)!=alen:
        print('Error in the Consensus!!!')
        print('consensus_len = ', len(consensus), 'alen = ', alen)
        print(seqs)
        sys.exit()
    else:
        return consensus