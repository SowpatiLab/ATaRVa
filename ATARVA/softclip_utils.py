import parasail
import edlib
import ahocorasick

from ATARVA.realignment import *
from ATARVA.validation_tests import *


def build_kmer_automaton(seq: str, k: int = 10):
    """Build once per sequence X."""
    A = ahocorasick.Automaton()
    for i in range(len(seq) - k + 1):
        A.add_word(seq[i:i + k], (i, seq[i:i + k]))
    A.make_automaton()
    return A


def has_kmer_match(automaton, Y: str) -> bool:
    """Fast existence check — O(|Y|)."""
    return any(automaton.iter(Y))


def check_flank(cooper, read, locus_start, locus_end, up_softclip, down_softclip):
    """
    Check and identify flanking sequences of a locus in the softclipped regions.

    :param cooper: Cooper object containing reference and other data
    :param read: pysam read object
    :param locus_start: start coordinate of the repeat in reference
    :param locus_end: end coordinate of the repeat in reference
    :param up_softclip: length of soft clip at the start of the read
    :param down_softclip: length of soft clip at the end of the read
    :return: dict with 'upstream' and 'downstream' coordinates in read, or None if locus not in softclips
    """

    NONREP_FLANK = 30
    score_threshold = int(2 * (0.9 * NONREP_FLANK))  # a match score of 90% as threshold

    # contains reference and read coordinates of upstream and downstream flanks of a locus; initialized to None
    result = {'upstream': None, 'downstream': None}  

    # looking for the upstream flank for loci that could be in upstream softclip
    if locus_start - NONREP_FLANK < read.ref_start and up_softclip > NONREP_FLANK:
        upstream     = cooper.ref.fetch(read.chrom, locus_start - NONREP_FLANK, locus_start)
        softclip_seq = read.query_sequence[:up_softclip]
        # build once per read/sequence X
        automaton  = build_kmer_automaton(upstream)
        kmer_match = has_kmer_match(automaton, softclip_seq)
        if not kmer_match: return None

        alignment_score, target_begin, target_end, query_begin, query_end, sCigar = stripSW(Inputs(softclip_seq, upstream), False)

        if alignment_score >= score_threshold and target_end > 0:
            result['upstream'] = (locus_start - NONREP_FLANK + query_begin,
                                  locus_start - NONREP_FLANK + query_end,
                                  target_begin, target_end, alignment_score)
            # check if the downstream flank is already covered by the read
            overlap            = max((locus_end + NONREP_FLANK) - read.ref_start, target_end + (locus_end - locus_start) + NONREP_FLANK - read.query_start)

            # If upstream flank found and locus <80% covered, look for downstream flank in same softclip
            if overlap <= 0:
                softclip_seq  = read.query_sequence[target_end:up_softclip]
                downstream = cooper.ref.fetch(read.chrom, locus_end, locus_end + NONREP_FLANK)
                alignment_score_down, target_begin_down, target_end_down, query_begin_down, query_end_down, sCigar_down = stripSW(Inputs(softclip_seq, downstream), False)

                if alignment_score_down >= score_threshold:
                    result['downstream'] = (locus_end + query_begin_down,
                                            locus_end + query_end_down,
                                            target_begin_down + (target_end),
                                            target_end_down   + (target_end), alignment_score_down)
            else:
                result['downstream'] = (read.ref_start, read.ref_start, read.query_start, read.query_start)  # Indicate that downstream flank is already covered

    # Check if downstream flank region is in right softclip
    if locus_end + NONREP_FLANK > read.ref_end and down_softclip > NONREP_FLANK:
        downstream     = cooper.ref.fetch(read.chrom, locus_end, locus_end + NONREP_FLANK)
        softclip_seq   = read.query_sequence[-down_softclip:]
        automaton  = build_kmer_automaton(downstream)
        kmer_match = has_kmer_match(automaton, softclip_seq)
        if not kmer_match: return None

        alignment_score, target_begin, target_end, query_begin, query_end, sCigar = stripSW(Inputs(softclip_seq, downstream), False)

        if alignment_score >= score_threshold and target_begin < len(softclip_seq):
            result['downstream'] = (locus_end + query_begin,
                                    locus_end + query_end,
                                    read.query_end + target_begin,
                                    read.query_end + target_end, alignment_score)
            overlap              = max(read.ref_end - (locus_start - NONREP_FLANK), -(target_begin - ((locus_end - locus_start) + NONREP_FLANK)))
            # If downstream flank found and locus <80% covered, look for upstream flank in same softclip
            if overlap <= 0:
                softclip_seq = read.query_sequence[read.query_end:read.query_end + target_begin]
                upstream = cooper.ref.fetch(read.chrom, locus_start - NONREP_FLANK, locus_start)
                alignment_score_up, target_begin_up, target_end_up, query_begin_up, query_end_up, sCigar_up = stripSW(Inputs(softclip_seq, upstream), False)

                if alignment_score_up >= score_threshold:
                    result['upstream'] = (locus_start - NONREP_FLANK + query_begin_up,
                                          locus_start - NONREP_FLANK + query_end_up,
                                          read.query_end + target_begin_up,
                                          read.query_end + target_end_up, alignment_score_up)
            else:
                result['upstream'] = (read.ref_end, read.ref_end, read.query_end, read.query_end)  # Indicate that upstream flank is already covered

    return result if result['upstream'] or result['downstream'] else None