from collections import Counter
from bitarray import bitarray

def convert_to_bitset(seq):
    lbit = {'A': '0', 'C': '0', 'G': '1', 'T': '1', 'N': '1'}
    rbit = {'A': '0', 'C': '1', 'G': '0', 'T': '1', 'N': '1'}
    
    lbitseq = bitarray()
    rbitseq = bitarray()
    
    for s in seq:
        lbitseq.extend(lbit.get(s, '1'))
        rbitseq.extend(rbit.get(s, '1'))
    
    return lbitseq, rbitseq

def shift_and_match(seq, motif_length):
    best_shift = motif_length
    max_matches = 0

    shift_values = set(range(1, 7)) | {motif_length - 1, motif_length, motif_length + 1}

    for shift in sorted(shift_values):
        if shift < 1:
            continue

        lbitseq, rbitseq = convert_to_bitset(seq)

        lmatch = ~(lbitseq ^ (lbitseq >> shift))
        rmatch = ~(rbitseq ^ (rbitseq >> shift))
        match = lmatch & rmatch

        match_positions = match.count()

        if match_positions > max_matches:
            max_matches = match_positions
            best_shift = shift

    return best_shift

def kmp_search_non_overlapping(text, pattern):
    def compute_lps(pattern):
        lps = [0] * len(pattern)
        length = 0
        i = 1
        while i < len(pattern):
            if pattern[i] == pattern[length]:
                length += 1
                lps[i] = length
                i += 1
            else:
                if length != 0:
                    length = lps[length - 1]
                else:
                    lps[i] = 0
                    i += 1
        return lps

    lps = compute_lps(pattern)
    result = []
    i = 0  
    j = 0  

    while i < len(text):
        if pattern[j] == text[i]:
            i += 1
            j += 1

        if j == len(pattern):  
            result.append(i - j)
            j = 0  
        elif i < len(text) and pattern[j] != text[i]:
            if j != 0:
                j = lps[j - 1]
            else:
                i += 1

    return result

def get_most_frequent_motif(sequence, motif_size):
    if len(sequence) < motif_size:  
        return sequence, None  

    repeating_units = [sequence[i:i + motif_size] for i in range(len(sequence) - motif_size + 1)]
    motif_counts = Counter(repeating_units)
    
    if not motif_counts:  
        return sequence, None  

    most_common = motif_counts.most_common(2) 
    
    primary_motif = most_common[0][0] if most_common else sequence
    secondary_motif = most_common[1][0] if len(most_common) > 1 else None

    return primary_motif, secondary_motif

def motif_decomposition(sequence, motif_size):
    best_motif_size = min(motif_size , shift_and_match(sequence, motif_size))  
    primary_motif, secondary_motif = get_most_frequent_motif(sequence, best_motif_size)
    
    positions = kmp_search_non_overlapping(sequence, primary_motif)
    
    if not positions:
        return sequence  
    
    decomposed_parts = []
    count = 1  

    if positions[0] != 0:
        decomposed_parts.append(sequence[:positions[0]])  
    
    for i in range(1, len(positions)):
        if positions[i] == positions[i - 1] + len(primary_motif):
            count += 1
        else:
            decomposed_parts.append(f"({primary_motif}){count}")
            interspersed = sequence[positions[i - 1] + len(primary_motif):positions[i]]
            if interspersed:
                if len(interspersed) >= motif_size * 2:
                    secondary_motif_size = shift_and_match(interspersed, motif_size)
                    secondary_decomposition = motif_decomposition(interspersed, secondary_motif_size)
                    decomposed_parts.append(secondary_decomposition)
                else:
                    decomposed_parts.append(interspersed)
            count = 1  

    decomposed_parts.append(f"({primary_motif}){count}")
    last_motif_end = positions[-1] + len(primary_motif)
    leftover_sequence = sequence[last_motif_end:]
    
    if leftover_sequence:
        if len(leftover_sequence) >= motif_size * 2:
            secondary_motif_size = shift_and_match(leftover_sequence, motif_size)
            secondary_decomposition = motif_decomposition(leftover_sequence, secondary_motif_size)
            decomposed_parts.append(secondary_decomposition)
        else:
            decomposed_parts.append(leftover_sequence)
    
    return "-".join(decomposed_parts)