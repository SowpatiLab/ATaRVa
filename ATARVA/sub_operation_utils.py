import string, sys
import numpy as np
from scipy import stats

from ATARVA.consensus import *
from ATARVA.decompose import motif_decomposition

methviz_tag = False
def set_methviz_tag(value):
    global methviz_tag
    methviz_tag = value


# base64 encodings
BASE64_CHARS = string.ascii_uppercase + string.ascii_lowercase + string.digits + '+//'
NO_MATCH  = -2   # sentinel for unmatched position
AMBIGUOUS = -1   # sentinel for ambiguous methylation call
METH_ENCODING = { NO_MATCH  : '*',   # error / missing call
                 AMBIGUOUS : '-' }   # ambiguous call


def mm_tag_extract(read, mod_probs):
    """
    record the positions with methylation calls
    
    :param read: pysam AlignedSegment object
    :param mod_probs: list of tuples with modified base positions and their respective probabilities
    :param meth_cutoff: score cutoff for calling methylation
    :param frwd_strand: boolean value for strand state; forward = True, reverse = False
    :return: none; the methylation positions are recorded in the read.methylation_calls attribute
    """

    last_index = len(read.query_sequence) - 1
    if (read.methyl_start != None) and (read.methyl_end != None):
        for pos, prob in mod_probs:
            meth_chunk_start = pos if read.is_forward else pos - 1 # to check the meth context, start index
            if read.methyl_start <= pos <= read.methyl_end:
                if (pos + 1 <= last_index) and (read.query_sequence[meth_chunk_start : meth_chunk_start + 2]=='CG'):
                    read.methylation_calls.append((pos, prob))


def align_cg_positions(consensus_cg_pos, read_cg_pos):
    """
    Align CpG positions from a read to consensus CpG positions.
    Returns indices into read_cg_pos that best match each consensus position.
    -2 indicates no match found for a consensus position.

    :param consensus_cg_pos: CpG positions in the consensus sequence
    :param read_cg_pos:      CpG positions in the read sequence
    :return:                 list of read indices aligned to each consensus position
    """
    consensus_gaps  = np.diff(consensus_cg_pos, prepend=0)
    read_gaps       = np.diff(read_cg_pos,       prepend=0)

    n_read_gaps     = len(read_gaps)
    read_exceeds    = n_read_gaps >= len(consensus_cg_pos)  # read has more CGs than consensus

    aligned_read_idx  = []   # matched read index for each consensus position
    read_cursor       = 0    # current position in read_gaps
    gap_remainder     = 0    # carry-over gap difference between iterations

    TOLERANCE         = 4    # max bp tolerance for initial position matching
    LAST_POS_WINDOW   = 10   # max bp tolerance for last position matching
    NO_MATCH          = -2   # sentinel for unmatched consensus position

    for cons_idx, cons_gap in enumerate(consensus_gaps):
        n_matched_before = len(aligned_read_idx)

        # --- find optimal start position ---
        if read_exceeds and cons_idx == 0 and read_cursor == 0:
            gap_diff = read_gaps[read_cursor] - cons_gap

            if abs(gap_diff) < TOLERANCE:
                # gaps are close enough — direct match
                aligned_read_idx.append(cons_idx)
                gap_remainder = gap_diff
                read_cursor  += 1
            else:
                # slide read positions to find better alignment start
                slide_dist    = consensus_cg_pos[cons_idx] - read_cg_pos[cons_idx]
                shifted_pos   = read_cg_pos - slide_dist

                dist_original = np.abs(consensus_cg_pos[cons_idx] - read_cg_pos)
                dist_shifted  = np.abs(consensus_cg_pos[cons_idx] - shifted_pos)

                # pick whichever alignment has lower total distance
                if dist_original.sum() <= dist_shifted.sum():
                    read_cursor = int(np.argmin(dist_original))
                else:
                    read_cursor = int(np.argmin(dist_shifted))

                aligned_read_idx.append(read_cursor)
                read_cursor += 1

            continue

        # --- scan read gaps to find best match for current consensus gap ---
        cumulative_gap = gap_remainder
        min_diffs      = [100_000]

        while read_cursor < n_read_gaps:
            cumulative_gap += read_gaps[read_cursor]
            current_diff    = abs(cumulative_gap - cons_gap)

            # diff increasing — previous position was better match
            if current_diff > min(min_diffs):
                best_cursor = read_cursor - 1
                if best_cursor not in aligned_read_idx:
                    aligned_read_idx.append(best_cursor)
                    gap_remainder = (cumulative_gap - read_gaps[read_cursor]) - cons_gap
                break

            # last position check
            is_last     = read_cursor + 1 == n_read_gaps
            last_in_win = abs(cons_gap - read_gaps[read_cursor]) < LAST_POS_WINDOW
            if is_last and last_in_win:
                if read_cursor not in aligned_read_idx:
                    aligned_read_idx.append(read_cursor)
                read_cursor += 1
                break

            min_diffs.append(current_diff)
            read_cursor += 1

        # --- no match found for this consensus position ---
        if len(aligned_read_idx) == n_matched_before:
            aligned_read_idx.append(NO_MATCH)

    return aligned_read_idx


def encode_methylation(score_matrix, pos_matrix, consensus_cgs):
    """
    encode methylation scores at consensus CpG positions into a compact string for visualization

    :param score_matrix: list of lists with methylation scores for each read at their respective CpG positions
    :param pos_matrix: list of lists with read CpG positions aligned to consensus CpG
    :param consensus_cgs: list of CpG positions in the consensus sequence
    :return: encoded methylation string where each character represents the consensus CpG position and its methylation status across reads
    """

    # --- align each read's CG positions to consensus ---
    aligned_pos_matrix = [
        align_cg_positions(consensus_cgs, np.array(read_cgs))
        for read_cgs in pos_matrix
    ]

    # --- build score matrix using aligned positions ---
    aligned_score_matrix = [
        [meth_scores[pos] if pos != NO_MATCH else NO_MATCH
         for pos in aligned_pos]
        for meth_scores, aligned_pos in zip(score_matrix, aligned_pos_matrix)
    ]

    # --- encode each consensus CG position ---
    encoded_meth = []
    for col in zip(*aligned_score_matrix):
        col_array = np.array(col)
        mode      = stats.mode(col_array, keepdims=True).mode[0]

        if mode in METH_ENCODING:
            encoded_meth.append(METH_ENCODING[mode])
            continue

        # filter sentinels and compute mean
        valid     = col_array[(col_array != AMBIGUOUS) & (col_array != NO_MATCH)]
        col_mean  = round(np.mean(valid) * 100)
        col_scaled = round(col_mean / 1.5625)    # scale to 0-64
        encoded_meth.append(BASE64_CHARS[col_scaled])

    return ''.join(encoded_meth)


def calculate_methylation(read_indices, read_methylation, consensus_seq):
    """
    calculate the methylation level at CG positions for a haplogroup based on the consensus sequence of the haplogroup

    :param read_indices: list of read ids supporting the locus
    :param methyl_probabilities: list of tuples with modified base positions and their respective probabilities for each read
    :param consensus_seq: the consensus sequence
    :return: average methylation level, number of reads with methylation information, encrypted methylation string for visualization
    """

    arr  = np.frombuffer(consensus_seq.upper().encode(), dtype=np.uint8)
    c_mask = arr[:-1] == ord('C')
    g_mask = arr[1:]  == ord('G')
    cg_positions = np.where(c_mask & g_mask)[0]

    if cg_positions.size == 0:
        return [None, None, None]

    methyl_nreads       = 0
    haplotype_avgmeth   = 0
    encoded_methylation = None
    score_matrix = []
    pos_matrix = []
    for read_index in read_indices:
        if read_methylation[read_index] is not None:
            methyl_nreads += 1
            haplotype_avgmeth += read_methylation[read_index][0]
            # for meth visualization
            if methviz_tag:
                pos_matrix.append(read_methylation[read_index][1])
                score_matrix.append(read_methylation[read_index][2])
            
    if methyl_nreads > 0:
        if methviz_tag:
            encoded_methylation = encode_methylation(score_matrix, pos_matrix, cg_positions)
        return [round(haplotype_avgmeth/methyl_nreads, 2), methyl_nreads, encoded_methylation]
    else:
        return [None, None, None]


def alt_sequence(read_alleles, hap_reads, motif_size):
    """
    generates alt sequence and reports characteristics

    :param read_alleles: dict with read_id as key and allele sequence as value
    :param hap_reads: list of read_ids belonging to the haplotype cluster
    :param motif_size: size of the motif for checking repetitiveness
    :return: ALT, allele_length, decomposed_seq, repetitive: alt sequence, its length, motif decomposition and repetitiveness status
    """

    # wouldn't most of the reads say deleted and only one read has some sequence ???
    ALT = ''
    seqs = [seq for seq in [read_alleles[read_id][0] for read_id in hap_reads] if seq!='']
    if seqs:
        ALT = consensus_seq_poa(seqs)
        allele_length = len(ALT)
    else:
        ALT = '<DEL>'
        allele_length = 0

    decomposed_seq = ''
    is_repetitive = True
    if allele_length and (motif_size <= 10):
        decomposed_seq, nonrep_fraction = motif_decomposition(ALT, motif_size)
        if nonrep_fraction > 0.30: is_repetitive = False
    return [ALT, allele_length, decomposed_seq, is_repetitive]
