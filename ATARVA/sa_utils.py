from ATARVA.process_softclips import align_sequences, generate_md_tag, generate_cs_tag
from ATARVA.process_softclips import _cigar_tuples, _cigar_string, _collapse_mismatches
from ATARVA.process_softclips import join_cigars, join_mdtags, join_cstags, valid_md, valid_cs

import re


def flip_cigar(cigar: str) -> str:
    """
    Flip a CIGAR string so that the roles of target and query are swapped
    (i.e. convert a target-referenced CIGAR into a query-referenced one,
    or vice versa).

    In CIGAR notation:
      - 'D' (deletion) = target has a base, query doesn't  -> becomes 'I'
      - 'I' (insertion) = query has a base, target doesn't -> becomes 'D'
      - 'M', '=', 'X', 'S', 'H', 'N', 'P' consume both / are symmetric
        and stay the same.

    Example:
        >>> flip_cigar("5M2D3M1I4M")
        '5M2I3M1D4M'
    """
    swap = {'D': 'I', 'I': 'D'}

    ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    if not ops:
        raise ValueError(f"Invalid CIGAR string: {cigar!r}")

    flipped = ''.join(f"{length}{swap.get(op, op)}" for length, op in ops)
    return flipped


def _query_length(cigartuples: list[tuple[int, int]]):
    query_length = 0
    for op, length in cigartuples:
        if op in (0, 1, 4, 5, 7, 8): # M/I/S/=/X consume query
            query_length += length
    return query_length

def _ref_length(cigartuples: list[tuple[int, int]]):
    ref_length = 0
    for op, length in cigartuples:
        if op in (0, 2, 3, 7, 8): # M/D/N/=/X consume reference
            ref_length += length
    return ref_length


def _trim_upcigar(cigar_tuples, qpos, rpos, repflank_start, repflank_end):
    """
    Trim the CIGAR from upstream to exclude repeat coords.

    :param cigar_tuples: list of tuples representing the CIGAR string
    :param qpos: query position
    :param rpos: reference position
    :param repflank_start: start position of the repeat flank
    """

    for idx, (op, length) in enumerate(cigar_tuples):
        if op == 4 or op == 5: qpos += length # soft clip
        if op == 1: qpos += length # insertion
        if op == 2: # deletion
            for _ in range(length):
                if _ == 0:
                    if rpos >= repflank_end:
                        trimmed_cigar = [(2, length)] + cigar_tuples[idx+1:]
                        return rpos, qpos, trimmed_cigar
                rpos += 1
                if rpos >= repflank_end:
                    trimmed_cigar = [(2, length - (_ + 1))] + cigar_tuples[idx+1:]
                    return rpos, qpos, trimmed_cigar

        if op == 0:
            for _ in range(length):
                if _ == 0:
                    if rpos >= repflank_end:
                        trimmed_cigar = [(0, length)] + cigar_tuples[idx+1:]
                        return rpos, qpos, trimmed_cigar
                rpos += 1
                qpos += 1
                if rpos >= repflank_end:
                    trimmed_cigar = [(0, length - (_ + 1))] + cigar_tuples[idx+1:]
                    return rpos, qpos, trimmed_cigar


def _trim_downcigar(cigar_tuples, qpos, rpos, repflank_start, repflank_end):
    """
    Trim the CIGAR from downstream to exclude repeat coords.

    :param cigar_tuples: list of tuples representing the CIGAR string
    :param qpos: query position
    :param rpos: reference position
    :param repflank_end: end position of the repeat flank
    """

    for idx, (op, length) in enumerate(cigar_tuples[::-1]):
        if op == 4 or op == 5: qpos -= length # soft clip
        if op == 1: qpos -= length # insertion
        if op == 2: # deletion
            for _ in range(length):
                if _ == 0:
                    if rpos <= repflank_start:
                        trimmed_cigar = cigar_tuples[:len(cigar_tuples) - (idx + 1)] + [(2, length)]
                        return rpos, qpos, trimmed_cigar
                rpos -= 1
                if rpos <= repflank_start:
                    trimmed_cigar = cigar_tuples[:len(cigar_tuples) - (idx + 1)] + [(2, length - (_ + 1))]
                    return rpos, qpos, trimmed_cigar

        if op == 0:
            for _ in range(length):
                if _ == 0:
                    if rpos <= repflank_start:
                        trimmed_cigar = cigar_tuples[:len(cigar_tuples) - (idx + 1)] + [(0, length)]
                        return rpos, qpos, trimmed_cigar
                rpos -= 1
                qpos -= 1
                if rpos <= repflank_start:
                    trimmed_cigar = cigar_tuples[:len(cigar_tuples) - (idx + 1)] + [(0, length - (_ + 1))]
                    return rpos, qpos, trimmed_cigar


def _conf_downpos(cigar_tuples, rpos, qpos, length_threshold=30):

    for idx, (op, length) in enumerate(cigar_tuples[::-1]):
        if op == 4 or op == 5: qpos -= length # soft clip
        if op == 1: qpos -= length # insertion
        if op == 2: rpos -= length # deletion

        if op == 0:
            if length >= length_threshold:  # match
                trimmed_cigar = cigar_tuples[:len(cigar_tuples) - idx] + [(0, length - length_threshold)]
                return rpos - length_threshold, qpos -  length_threshold, trimmed_cigar
            rpos -= length
            qpos -= length

    best_match_score = 0
    best_match_rpos  = -1
    best_match_qpos  = -1
    best_match_cigar = None
    aln_length = 0
    aln_score  = 0
    prev_op = None

    for idx, op, length in enumerate(cigar_tuples[::-1]):

        if op == 4 or op == 5: qpos -= length # soft clip

        if op == 1: # insertion
            for _ in range(length):
                qpos -= 1
                if aln_length < length_threshold:
                    aln_length += 1
                    aln_score  -= 1
                    prev_op = 1
                else:
                    if prev_op == 0: aln_score -= 1
                    else: aln_score += 1
                    aln_score -= 1
                    prev_op = 1
                    if aln_score > best_match_score:
                        best_match_score = aln_score
                        best_match_rpos = rpos
                        best_match_qpos = qpos
                        best_match_cigar = cigar_tuples[:len(cigar_tuples) - (idx + 1)] + [(1, length - (_ + 1))]

        if op == 2: # deletion
            for _ in range(length):
                rpos -= 1
                if aln_length < length_threshold:
                    aln_length += 1
                    aln_score  -= 1
                    prev_op = 2
                else:
                    if prev_op == 0: aln_score -= 1
                    else: aln_score += 1
                    aln_score -= 1
                    prev_op = 2
                    if aln_score > best_match_score:
                        best_match_score = aln_score
                        best_match_rpos = rpos
                        best_match_qpos = qpos
                        best_match_cigar = cigar_tuples[:len(cigar_tuples) - (idx + 1)] + [(2, length - (_ + 1))]

        if op == 0: # match
            for _ in range(length):
                rpos -= 1
                qpos -= 1
                if aln_length < length_threshold:
                    aln_length += 1
                    aln_score  += 1
                    prev_op = 0
                else:
                    if prev_op == 0: aln_score -= 1
                    else: aln_score -= 1
                    aln_score += 1
                    prev_op = 0
                    if aln_score > best_match_score:
                        best_match_score = aln_score
                        best_match_rpos = rpos
                        best_match_qpos = qpos
                        best_match_cigar = cigar_tuples[:len(cigar_tuples) - (idx + 1)] + [(0, length - (_ + 1))]

    return best_match_rpos, best_match_qpos, best_match_cigar


def _conf_uppos(cigar_tuples, rpos, qpos, length_threshold=30):

    for idx, (op, length) in enumerate(cigar_tuples):
        if op == 4 or op == 5: qpos += length # soft clip
        if op == 1: qpos += length # insertion
        if op == 2: rpos += length # deletion

        if op == 0:
            if length >= length_threshold:  # match
                trimmed_cigar = [(0, length - length_threshold)] + cigar_tuples[idx+1:]
                return rpos + length_threshold, qpos + length_threshold, trimmed_cigar
            rpos += length
            qpos += length

    best_match_score = 0
    best_match_rpos  = -1
    best_match_qpos  = -1
    best_match_cigar = None
    aln_length = 0
    aln_score  = 0
    prev_op = None

    for idx, op, length in enumerate(cigar_tuples):

        if op == 4 or op == 5: qpos += length # soft clip

        if op == 1: # insertion
            for _ in range(length):
                qpos += 1
                if aln_length < length_threshold:
                    aln_length += 1
                    aln_score  -= 1
                    prev_op = 1
                else:
                    if prev_op == 0: aln_score -= 1
                    else: aln_score += 1
                    aln_score -= 1
                    prev_op = 1
                    if aln_score > best_match_score:
                        best_match_score = aln_score
                        best_match_rpos = rpos
                        best_match_qpos = qpos
                        best_match_cigar = [(1, length - (_ + 1))] + cigar_tuples[idx + 1:]

        if op == 2: # deletion
            for _ in range(length):
                rpos += 1
                if aln_length < length_threshold:
                    aln_length += 1
                    aln_score  -= 1
                    prev_op = 2
                else:
                    if prev_op == 0: aln_score -= 1
                    else: aln_score += 1
                    aln_score -= 1
                    prev_op = 2
                    if aln_score > best_match_score:
                        best_match_score = aln_score
                        best_match_rpos = rpos
                        best_match_qpos = qpos
                        best_match_cigar = [(2, length - (_ + 1))] + cigar_tuples[idx + 1:]

        if op == 0: # match
            for _ in range(length):
                rpos += 1
                qpos += 1
                if aln_length < length_threshold:
                    aln_length += 1
                    aln_score  += 1
                    prev_op = 0
                else:
                    if prev_op == 0: aln_score -= 1
                    else: aln_score -= 1
                    aln_score += 1
                    prev_op = 0
                    if aln_score > best_match_score:
                        best_match_score = aln_score
                        best_match_rpos = rpos
                        best_match_qpos = qpos
                        best_match_cigar = [(0, length - (_ + 1))] + cigar_tuples[idx + 1:]

    return best_match_rpos, best_match_qpos, best_match_cigar


def _collapse_matches(cigar, match_char='M'):
    """
    collapse matches within matches in CIGAR string

    :param cigar: CIGAR string (e.g., "10M1I5M8M")
    :param match_char: character to represent matches (e.g., 'M' or '=')

    :return collapsed CIGAR string
    """

    length = ''
    prev_op  = None

    match_length = 0
    new_cigar = ''

    for c in cigar:
        if c.isdigit():
            length += c
            continue
        else:
            if c == match_char:
                match_length += int(length)
            else:
                if match_length > 0:
                    new_cigar += f'{match_length}{match_char}'
                    match_length = 0
                new_cigar += f'{length}{c}'
            length = ''
    if match_length > 0:
        new_cigar += f'{match_length}{match_char}'

    return new_cigar


def process_upstreamsa(cooper, ref, read, sa_read, sa_start, sa_end, sa_cigar, repeat_flank_start, repeat_flank_end):
    """
    Process the supplementary alignments of a read to check if any of them cover the locus of interest.
    """

    sa_startclip = 0
    if sa_read.cigartuples[0][0] == 4 or sa_read.cigartuples[0][0] == 5:
        sa_startclip = sa_read.cigartuples[0][1]
    read_startclip = 0
    if read.cigartuples[0][0] == 4 or read.cigartuples[0][0] == 5:
        read_startclip = read.cigartuples[0][1]
    if sa_startclip >= read_startclip: return False

    end_softclip = 0
    if read.cigartuples[-1][0] == 4:
        end_softclip = read.cigartuples[-1][1]

    if not _query_length(read.cigartuples) == _query_length(sa_read.cigartuples):
        return False
    if read.ref_end <= sa_read.ref_end: return False

    query_len = _query_length(read.cigartuples)
    query_seq = ""
    primary = "READ"
    if len(read.query_sequence) == query_len and len(read.query_qualities) == query_len:
        query_seq = read.query_sequence
        primary = "READ"
    elif len(sa_read.query_sequence) == query_len and len(sa_read.query_qualities) == query_len:
        query_seq = sa_read.query_sequence
        primary = "SA"

    if query_seq == "":
        return False

    flank = 0
    result_up   = _trim_upcigar(read.cigartuples, 0, read.ref_start, repeat_flank_start-flank, repeat_flank_end+flank)
    result_down = _trim_downcigar(sa_cigar, _query_length(sa_cigar), sa_end, repeat_flank_start-flank, repeat_flank_end+flank)
    if result_up is None or result_down is None:
        return False
    
    pa_rpos, pa_qpos, pa_trimmed_cigar = result_up
    sa_rpos, sa_qpos, sa_trimmed_cigar = result_down

    sa_conf_rpos, sa_conf_qpos, sa_conf_cigar = sa_rpos, sa_qpos, sa_trimmed_cigar
    pa_conf_rpos, pa_conf_qpos, pa_conf_cigar = pa_rpos, pa_qpos, pa_trimmed_cigar

    sa_target = ref.fetch(cooper.chrom, sa_conf_rpos, pa_conf_rpos)
    sa_query  = query_seq[sa_conf_qpos:pa_conf_qpos]
    if sa_conf_qpos >= pa_conf_qpos:
        return False

    sa_cigarstring, score = align_sequences(sa_query, sa_target) # faster if the the query is shorter
    pa_conf_cigarstring   = _cigar_string(pa_conf_cigar)
    joined_cigar = f'{sa_conf_qpos}S' + join_cigars(sa_cigarstring, pa_conf_cigarstring)

    if primary != "READ":
        read.query_sequence  = query_seq
        read.query_qualities = sa_read.query_qualities
        read.mod_bases       = sa_read.modified_bases
    if '=' not in read.cigarstring:
        joined_cigar = joined_cigar.replace('=', 'M').replace('X', 'M')
    if 'X' not in read.cigarstring:
        joined_cigar = joined_cigar.replace('=', 'M').replace('X', 'M')
        match_char = 'M' if '=' not in read.cigarstring else '='
        joined_cigar = _collapse_matches(joined_cigar, match_char)

    read.cigarstring = joined_cigar
    read.cigartuples = _cigar_tuples(joined_cigar)
    read.ref_start   = sa_conf_rpos
    read.query_start = sa_conf_qpos
    if read.has_tag('MD'):
        sa_mdtag    = generate_md_tag(sa_cigarstring, sa_target, sa_query)
        qseq        = query_seq[pa_conf_qpos:-end_softclip] if end_softclip > 0 else query_seq[pa_conf_qpos:]
        pa_mdtag    = generate_md_tag(_cigar_string(pa_conf_cigar), ref.fetch(cooper.chrom, pa_conf_rpos, read.ref_end), qseq)
        read.md_tag = join_mdtags(sa_mdtag, pa_mdtag)
        valid = valid_md(read.md_tag, joined_cigar)
        assert valid, f"Invalid MD tag: {read.md_tag} for CIGAR: {joined_cigar}"
    if read.has_tag('cs'):
        sa_cstag = generate_cs_tag(sa_query, sa_target, sa_cigarstring)
        qseq     = query_seq[pa_conf_qpos:-end_softclip] if end_softclip > 0 else query_seq[pa_conf_qpos:]
        pa_cstag = generate_cs_tag(qseq, ref.fetch(cooper.chrom, pa_conf_rpos, read.ref_end), pa_conf_cigarstring)
        read.cs_tag = join_cstags(sa_cstag, pa_cstag)
        valid = valid_cs(read.cs_tag, joined_cigar)
        assert valid, f"Invalid CS tag: {read.cs_tag} for CIGAR: {joined_cigar}"

    return True
