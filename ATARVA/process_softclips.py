import parasail
import edlib
import ahocorasick

from ATARVA.realignment import stripSW, Inputs
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



def _mean(numbers):
    """
    calculate mean of a list of numbers

    :param numbers: list of numbers
    :return:        mean value
    """
    return float(sum(numbers)) / max(len(numbers), 1)


def generate_cs_tag(query_seq, ref_seq, cigar):
    """
    Generate CS tag from CIGAR string with X (mismatch) operations.
    CS tag format: = match, * mismatch, + insertion, - deletion

    :param query_seq:    query/read sequence
    :param ref_seq:      reference sequence for the aligned region
    :param cigartuples:  list of (op, length) cigar tuples
    :return:             CS tag string
    """
    cs   = []
    qpos = 0    # query position
    rpos = 0    # reference position
    cigartuples = _cigar_tuples(cigar)

    for op, length in cigartuples:

        # ── match / mismatch (M) ──────────────────────────────────
        if op == 0:
            match_length = 0
            for i in range(length):
                q = query_seq[qpos + i]
                r = ref_seq[rpos + i]
                if q == r:
                    match_length += 1
                else:
                    if match_length > 0:
                        cs.append(f':{match_length}')
                    cs.append(f'*{r}{q}')   # mismatch ref→query
                    match_length = 0

            qpos += length
            rpos += length

        # ── sequence match (=) ────────────────────────────────────
        elif op == 7:
            cs.append(f':{length}')
            qpos += length
            rpos += length

        # ── sequence mismatch (X) ─────────────────────────────────
        elif op == 8:
            for i in range(length):
                cs.append(f'*{ref_seq[rpos + i]}{query_seq[qpos + i]}')
            qpos += length
            rpos += length

        # ── insertion (I) ─────────────────────────────────────────
        elif op == 1:
            cs.append(f'+{query_seq[qpos:qpos + length]}')
            qpos += length

        # ── deletion (D) ─────────────────────────────────────────
        elif op == 2:
            cs.append(f'-{ref_seq[rpos:rpos + length]}')
            rpos += length

        # ── soft clip (S) — not included in CS tag ────────────────
        elif op == 4:
            qpos += length

        # ── hard clip (H) / padding (P) — skip ───────────────────
        elif op in (5, 6):
            pass

        # ── reference skip (N) ───────────────────────────────────
        elif op == 3:
            cs.append(f'-{ref_seq[rpos:rpos + length]}')
            rpos += length

    return ''.join(cs)


def generate_md_tag(cigar, reference, query):
    """
    generate MD tag from CIGAR, reference/target sequence, and query/read sequence.

    :param cigar: CIGAR string (e.g., "10M1I5M2D8M")
    :param target: reference-aligned sequence span consumed by CIGAR ref ops
    :param query: read-aligned sequence span consumed by CIGAR query ops
    :return: MD tag string (e.g., "10A5^CT8")
    """
    qpos = 0
    rpos = 0

    prev_op = None

    md = []
    match_count = 0

    i = 0
    while i < len(cigar):
        # parse length
        mismatch = False
        num = ''
        while i < len(cigar) and cigar[i].isdigit():
            num += cigar[i]
            i += 1
        if num == '': num = 0
        else: num = int(num)
        op = cigar[i]
        i += 1

        if op == '=':
            match_count += num
            qpos += num
            rpos += num

        elif op == 'X':
            if prev_op in ('X', 'D'):
                md.append('0')
            for _ in range(num):
                if match_count > 0:
                    md.append(str(match_count))
                    match_count = 0
                if prev_op == 'X' and md[-1] != '0': md.append('0')
                md.append(reference[rpos])
                prev_op = 'X'
                qpos += 1
                rpos += 1

        elif op == 'D':
            if prev_op in ('X', 'D'):
                md.append('0')
            if match_count > 0:
                md.append(str(match_count))
                match_count = 0
            md.append("^" + reference[rpos:rpos+num])
            rpos += num

        elif op == 'I':
            # insertion: only consumes query
            qpos += num

        elif op == 'M':
            # fallback if M is present
            for _ in range(num):
                if query[qpos] == reference[rpos]:
                    mismatch = False
                    match_count += 1
                else:
                    if _ == 0 and prev_op in ('X', 'D'):
                        md.append('0')
                    mismatch = True
                    if match_count > 0:
                        md.append(str(match_count))
                        match_count = 0
                    md.append(reference[rpos])
                qpos += 1
                rpos += 1

        elif op == 'S':
            # soft clip: only consumes query, not represented in MD
            qpos += num

        prev_op = op
        if mismatch: prev_op = 'X'

    if match_count > 0:
        md.append(str(match_count))

    return ''.join(md)


def join_cstags(cs_left: str, cs_right: str) -> str:
    """
    join two continuous minimap2 cs tags (short format) and normalize.

    Supports tokens:
      :<num>   match length
      *ac      substitution
      +acgt    insertion
      -acgt    deletion
      ~gt12ag  intron/splice

    Example:
      cs_left=":10+aa:5"
      cs_right=":3-tt:7"
      -> ":10+aa:8-tt:7"
    
    :param cs_left:  left cs tag fragment (e.g., ":10+aa:5")
    :param cs_right: right cs tag fragment (e.g., ":3-tt:7")
    :return:         joined cs tag (e.g., ":10+aa:8-tt:7")
    """

    operations = {':', '-', '+', '*', '=', '~', 'B'}
    left_op    = ''
    left_opdat = ''
    for x in cs_left[::-1]:
        if x in operations:
            left_op = x; break
        left_opdat = x + left_opdat

    right_op = ''
    right_opdat = ''
    for x in cs_right:
        if x in operations:
            if not right_op:
                right_op = x
            else: break
        else:
            right_opdat += x

    if left_op == right_op:
        if left_op == '*':
            return cs_left + cs_right
        elif left_op == ':' or left_op == '~' or left_op == 'B':
            merged_opdat = str(int(left_opdat) + int(right_opdat))
        else:
            merged_opdat = left_opdat + right_opdat
        return cs_left[:-len(left_opdat)-1] + left_op + merged_opdat + cs_right[len(right_opdat)+1:]
    else:
        return cs_left + cs_right


def right_md_token(md):
    i = 0
    if md[i].isdigit():
        num = ''
        while i < len(md) and md[i].isdigit():
            num += md[i]; i += 1
        return int(num), i

    elif md[i] == '^':
        i += 1
        seq = ''
        while i < len(md) and md[i].isalpha():
            seq += md[i]; i += 1
        return '^' + seq, i
    else:
        return md[i], i + 1


def left_md_token(md):
    i = len(md) - 1
    if md[i].isdigit():
        num = ''
        while i >= 0 and md[i].isdigit():
            num = md[i] + num
            i -= 1
        return int(num), i + 1 

    else:
        if i > 0 and md[i-1].isdigit():
            return md[i], i
        seq = ''
        while i >= 0 and md[i] != '^':
            seq = md[i] + seq
            i -= 1
        return '^' + seq, i


def join_mdtags(md1, md2):
    if md1 == '': return md2
    if md2 == '': return md1

    t1, t1_idx = left_md_token(md1)
    t2, t2_idx = right_md_token(md2)

    if not t1:
        return md2
    if not t2:
        return md1

    merged = ''
    # merge boundary numbers
    if isinstance(t1, int) and isinstance(t2, int):
        merged = t1 + t2
    elif not isinstance(t1, int) and not isinstance(t2, int):
        if t1[0] == '^' and t2[0] == '^':
            merged = t1 + t2[1:]
        else:
            merged = t1 + '0' + t2
    else:
        merged = str(t1) + str(t2)

    return md1[:t1_idx] + str(merged) + md2[t2_idx:]


def join_cigars(cigar_left: str, cigar_right: str) -> str:
    """
    join two continuous CIGAR strings and normalize.
    :param cigar_left:  left CIGAR string (e.g., "10M1I5M")
    :param cigar_right: right CIGAR string (e.g., "2M3D8M")
    :return:            joined CIGAR string (e.g., "10M1I7M3D8M")
    """

    operations = {'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X', 'B'}
    left_op    = ''
    left_opdat = ''
    left_index = len(cigar_left)
    for x in cigar_left[::-1]:
        if x in operations:
            if left_op == '':
                left_op = x
            else: break
        else:
            left_opdat = x + left_opdat
        left_index -= 1

    right_op = ''
    right_opdat = ''
    right_index = 0
    for x in cigar_right:
        if x in operations:
            right_op = x
            right_index += 1; break
        else:
            right_opdat += x
        right_index += 1

    if left_op == right_op:
        merged_opdat = str(int(left_opdat) + int(right_opdat))
        return cigar_left[:left_index] + merged_opdat + left_op + cigar_right[right_index:]
    else:
        return cigar_left + cigar_right


def nw_align(query: str, target: str, gap_open: int = 6, gap_extend: int = 3) -> tuple[str, int]:
    """
    Needleman-Wunsch global alignment returning CIGAR and score.

    :param query:      query sequence
    :param target:     reference sequence
    :param gap_open:   gap open penalty
    :param gap_extend: gap extend penalty
    :return:           (cigar_string, alignment_score)
    """
    score_matrix = parasail.matrix_create("ACGT", 2, -2)  # simple match/mismatch scoring
    result = parasail.nw_trace(
        query, target,
        gap_open, gap_extend,
        score_matrix        # DNA scoring matrix
    )
    return result.cigar.decode.decode('utf-8'), result.score


def align_sequences(query: str, target: str, gap_open: int = 6, gap_extend: int = 3) -> tuple[str, int]:
    """
    Align query to target using parasail if both sequences < 10kb, otherwise use edlib.
    
    :param query:      query sequence
    :param target:     reference sequence
    :param gap_open:   gap open penalty (used only for parasail)
    :param gap_extend: gap extend penalty (used only for parasail)
    :return:           (cigar_string, alignment_score)
    """

    query_len = len(query)
    target_len = len(target)
    threshold = 10000  # 10kb

    # Use parasail if both sequences are below 10kb
    if query_len < threshold and target_len < threshold:
        return nw_align(query, target, gap_open, gap_extend)

    # Use edlib for larger sequences
    else:
        alignment = edlib.align(query, target, task='path')

        cigar = alignment['cigar']
        # Estimate score from alignment (edlib doesn't provide a direct score, so we use edit distance)
        # Negative edit distance as a rough score metric
        score = -alignment['editDistance']

        return cigar, score


def _edit_distance(cigar: str):
    """
    calculate edit distance from CIGAR string (considering M/I/D/X as edit operations)

    :param cigar: CIGAR string (e.g., "10M1I5M2D8M")
    :return:      edit distance
    """

    length = ''
    edit_distance = 0

    for c in cigar:
        if c.isdigit():
            length += c
            continue
        else:
            if c in ('I', 'D', 'X'): edit_distance += int(length)
            length = ''

    return edit_distance


def _match_lengths(cigar: str):
    """
    list of match lengths from CIGAR string (considering M/= as matches)
    NOTE: if X in CIGAR it is merged into M

    :param cigar: CIGAR string (e.g., "10M1I5M2D8M")
    :return:      list of match lengths (e.g., [10, 5, 8])
    """


    length = ''
    match_break = False
    prev_op  = None
    prev_len = 0

    match_lengths = []
    match_length = 0

    for c in cigar:
        if c.isdigit():
            length += c
            continue
        else:
            if c == 'M' or c == '=':
                match_length = int(length)
                if match_break and prev_op == 'X':
                    match_length += prev_len
                match_break = False
            if c == 'X':
                if match_break: pass
                else:
                    match_length += int(length)
            if c in ('I', 'D'):
                if match_length > 0: match_lengths.append(match_length)
                match_break = True
                match_length = 0
            prev_op = c
            prev_len = int(length)
            length = ''

    return match_lengths


def _collapse_mismatches(cigar: str, match_char: str):
    """
    collapse mismatches within matches in CIGAR string

    :param cigar: CIGAR string (e.g., "10M1I5M2D8M")
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
            if c == 'M' or c == '=' or c == 'X':
                match_length += int(length)

            elif c in ('I', 'D', 'B', 'N'):
                if match_length > 0:
                    new_cigar += f'{match_length}{match_char}'
                    match_length = 0
                new_cigar += f'{length}{c}'
                match_length = 0
            elif c == 'S':
                new_cigar += f'{length}S'
            prev_op = c
            length = ''
    if match_length > 0:
        new_cigar += f'{match_length}{match_char}'

    return new_cigar


def _cigar_tuples(cigar: str):
    length = ''
    cigar_map = {'M': 0, 'I': 1, 'D': 2, 'N': 3, 'S': 4, 'H': 5, 'P': 6, '=': 7, 'X': 8, 'B': 9}
    cigar_tuples = []

    for c in cigar:
        if c.isdigit():
            length += c
            continue
        else:
            cigar_tuples.append((cigar_map[c], int(length)))
            length = ''

    return cigar_tuples


def _query_length(cigartuples: list[tuple[int, int]]):
    query_length = 0
    for op, length in cigartuples:
        if op in (0, 1, 4, 7, 8): # M/I/S/=/X consume query
            query_length += length
    return query_length


def _ref_length(cigartuples: list[tuple[int, int]]):
    ref_length = 0
    for op, length in cigartuples:
        if op in (0, 2, 3, 7, 8): # M/D/N/=/X consume reference
            ref_length += length
    return ref_length


def strip_softclip(cigar: str, dir: str) -> tuple[str, int]:
    """
    strip leading and trailing soft-clips from CIGAR string

    :param cigar: CIGAR string (e.g., "10S90M5S")
    :param dir: direction of soft-clip to strip ("left" or "right" or "both")
    :return:      stripped CIGAR string
    """

    left_softclip_length = 0
    right_softclip_length = 0

    if dir in ("right", "both") and cigar.endswith('S'):
        i = len(cigar) - 2
        while i >= 0 and cigar[i].isdigit():
            i -= 1
        right_softclip_length = int(cigar[i+1:-1])
        cigar = cigar[:i+1]

    if dir in ("left", "both"):
        length = ''
        for i, c in enumerate(cigar):
            if c.isdigit():
                continue
            elif c == 'S':
                i += 1; break
            else:
                return cigar
        cigar = cigar[i:]

    return cigar


def md_stats(md):
    matches = 0
    mismatches = 0
    deletions = 0

    i = 0
    while i < len(md):
        if md[i].isdigit():
            num = ""
            while i < len(md) and md[i].isdigit():
                num = num + md[i]
                i += 1
            matches += int(num)
        elif md[i] == '^':
            i += 1
            while i < len(md) and md[i].isalpha():
                deletions += 1
                i += 1
        else:
            mismatches += 1
            i += 1

    return matches, mismatches, deletions


def cigar_stats(cigar):
    match_len = 0
    del_len = 0

    length = ''
    op = ''
    for c in cigar:
        if c.isdigit():
            length += c
        else:
            op = c
            if op in ('M', '=', 'X'):
                match_len += int(length)
            elif op == 'D':
                del_len += int(length)
            length = ''
            op = ''

    return match_len, del_len


def valid_md(md, cigar):
    m, x, d_md = md_stats(md)
    m_cigar, d_cigar = cigar_stats(cigar)

    return (m + x == m_cigar) and (d_md == d_cigar)


def cs_stats(cs):
    """
    Extract match and mismatch counts from CS tag.
    CS tag format: :num (match), *ref_alt (mismatch), +seq (insertion), -seq (deletion)
    
    :param cs: CS tag string
    :return: tuple (matches, mismatches)
    """
    matches = 0
    mismatches = 0

    i = 0
    while i < len(cs):
        if cs[i] == ':':
            i += 1
            num = ""
            while i < len(cs) and cs[i].isdigit():
                num += cs[i]
                i += 1
            matches += int(num)
        elif cs[i] == '*':
            mismatches += 1
            i += 3  # skip ref and alt bases
        elif cs[i] == '+':
            i += 1
            while i < len(cs) and cs[i].isalpha():
                i += 1
        elif cs[i] == '-':
            i += 1
            while i < len(cs) and cs[i].isalpha():
                i += 1
        else:
            i += 1

    return matches, mismatches


def valid_cs(cs, cigar):
    """
    Validate CS tag against CIGAR string by comparing match and mismatch counts.
    
    :param cs: CS tag string
    :param cigar: CIGAR string
    :return: True if valid, False otherwise
    """
    cs_m, cs_x = cs_stats(cs)
    m_cigar, d_cigar = cigar_stats(cigar)

    return (cs_m + cs_x == m_cigar)


def detect_flank(cooper, read, locus_start, locus_end, stream):
    """
    Identify the flanking sequences of the repeat region in the read
    Args:
        read: pysam read object
        fasta: pysam fasta object
        flank_length: length of the flanking sequence to be considered
        locus_start: start coordinate of the repeat in reference
        locus_end: end coordinate of the repeat in reference

    Returns:
        [start, end, sub_cigar]
        start and end coordinates of the flanking sequence in the read
        sub_cigar: CIGAR string of the repeat alignment to the reference in the read
    """

    flank = 50
    score_threshold = int(2*(0.9*flank))

    if stream: # if its upstream (True)

        # finding the upstream flank in the read 
        upstream = cooper.ref.fetch(read.chrom, locus_start - flank, locus_start)
        alignment_score, _align_score, target_begin, target_end, query_begin, query_end, sCigar = stripSW(Inputs(read.query_sequence, upstream), False)

        # target_begin is the start of the alignment in the read
        # query_begin  is the start of the alignment in the reference flank
        # target_end   is the end of the alignment in the read
        # query_end    is the end of the alignment in the reference flank

        # if the upstream flank aligns to the already aligned region(by aligner), there is no expansion in the soft-clip, here target_end is dummy
        if target_end > read.query_start: return [-1, -1]

        # reject read if less than 70% match with the reference flank sequence
        if alignment_score < score_threshold: return [-1, -1]

        new_read_query_start = target_begin - 1
        new_read_ref_start   = locus_start - flank + query_begin - 1

        return [new_read_query_start, new_read_ref_start]

    else: # if its downstream (False)

        # finding the downstream flank in the read
        downstream = cooper.ref.fetch(read.chrom, locus_end, locus_end + flank)
        alignment_score, _align_score, target_begin, target_end, query_begin, query_end, sCigar = stripSW(Inputs(read.query_sequence, downstream), False)

        # target_begin is the start of the alignment in the read
        # query_begin  is the start of the alignment in the reference flank
        # target_end   is the end of the alignment in the read
        # query_end    is the end of the alignment in the reference flank

        # if the downstream flank aligns to the already aligned region(by aligner), there is no expansion in the soft-clip, here target_begin is dummy
        if target_begin <= read.query_end: return [-1, -1]

        # reject read if less than 70% match with the reference flank sequence
        if alignment_score < score_threshold: return [-1, -1]

        new_read_query_end = target_end - 1
        new_read_ref_end   = locus_end + query_end - 1

        return [new_read_query_end, new_read_ref_end]


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
    low_score_threshold = int(2 * (0.8 * NONREP_FLANK))  # a match score of 90% as threshold

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

        alignment_score, _align_score, target_begin, target_end, query_begin, query_end, sCigar = stripSW(Inputs(softclip_seq, upstream), False)

        if alignment_score >= score_threshold and target_end > 0 and _align_score < low_score_threshold:
            result['upstream'] = (locus_start - NONREP_FLANK + query_begin,
                                  locus_start - NONREP_FLANK + query_end,
                                  target_begin, target_end, alignment_score)
            # check if the downstream flank is already covered by the read
            overlap            = max((locus_end + NONREP_FLANK) - read.ref_start, target_end + (locus_end - locus_start) + NONREP_FLANK - read.query_start)

            # If upstream flank found and locus <80% covered, look for downstream flank in same softclip
            if overlap <= 0:
                softclip_seq  = read.query_sequence[target_end:up_softclip]
                downstream = cooper.ref.fetch(read.chrom, locus_end, locus_end + NONREP_FLANK)
                alignment_score_down, _align_score_down, target_begin_down, target_end_down, query_begin_down, query_end_down, sCigar_down = stripSW(Inputs(softclip_seq, downstream), False)

                if alignment_score_down >= score_threshold and _align_score_down < low_score_threshold:
                    result['downstream'] = (locus_end + query_begin_down,
                                            locus_end + query_end_down,
                                            target_begin_down + (target_end),
                                            target_end_down   + (target_end), alignment_score_down)
            else:
                result['downstream'] = (read.ref_start, read.ref_start, read.query_start, read.query_start, None)  # Indicate that downstream flank is already covered

    # Check if downstream flank region is in right softclip
    if locus_end + NONREP_FLANK > read.ref_end and down_softclip > NONREP_FLANK:
        downstream     = cooper.ref.fetch(read.chrom, locus_end, locus_end + NONREP_FLANK)
        softclip_seq   = read.query_sequence[-down_softclip:]
        automaton  = build_kmer_automaton(downstream)
        kmer_match = has_kmer_match(automaton, softclip_seq)
        if not kmer_match: return None

        alignment_score, _align_score, target_begin, target_end, query_begin, query_end, sCigar = stripSW(Inputs(softclip_seq, downstream), False)

        if alignment_score >= score_threshold and target_begin < len(softclip_seq) and _align_score < low_score_threshold:
            result['downstream'] = (locus_end + query_begin,
                                    locus_end + query_end,
                                    read.query_end + target_begin,
                                    read.query_end + target_end, alignment_score)
            overlap              = max(read.ref_end - (locus_start - NONREP_FLANK), -(target_begin - ((locus_end - locus_start) + NONREP_FLANK)))
            # If downstream flank found and locus <80% covered, look for upstream flank in same softclip
            if overlap <= 0:
                softclip_seq = read.query_sequence[read.query_end:read.query_end + target_begin]
                upstream = cooper.ref.fetch(read.chrom, locus_start - NONREP_FLANK, locus_start)
                alignment_score_up, _align_score_up, target_begin_up, target_end_up, query_begin_up, query_end_up, sCigar_up = stripSW(Inputs(softclip_seq, upstream), False)

                if alignment_score_up >= score_threshold and _align_score_up < low_score_threshold:
                    result['upstream'] = (locus_start - NONREP_FLANK + query_begin_up,
                                          locus_start - NONREP_FLANK + query_end_up,
                                          read.query_end + target_begin_up,
                                          read.query_end + target_end_up, alignment_score_up)
            else:
                result['upstream'] = (read.ref_end, read.ref_end, read.query_end, read.query_end, None)  # Indicate that upstream flank is already covered

    return result if result['upstream'] or result['downstream'] else None


def align_flank_stretch(cooper, read, flank_coords, flank_type):
    """
    Align a single flank stretch to reference and generate CIGAR string.

    :param cooper: Cooper object containing reference
    :param read: pysam read object
    :param flank_coords: tuple (ref_start, ref_end, query_start, query_end)
    :param flank_type: 'upstream' or 'downstream'
    :return: dict with alignment info {cigar, score, ref_start, ref_end, query_start, query_end}
    """

    if not flank_coords:
        return None

    ref_start, ref_end, query_start, query_end = flank_coords

    # Extract sequences for this stretch
    query_seq = read.query_sequence[query_start:query_end]
    ref_seq = cooper.ref.fetch(read.chrom, ref_start, ref_end)

    # Validate sequence lengths match coordinate ranges
    if len(query_seq) != (query_end - query_start) or len(ref_seq) != (ref_end - ref_start):
        print(f"Warning: Sequence lengths do not match coordinate ranges for flank stretch. "
              f"Query length: {len(query_seq)}, Expected: {query_end - query_start}; "
              f"Ref length: {len(ref_seq)}, Expected: {ref_end - ref_start}.")
        return None

    if not query_seq or not ref_seq:
        print(f"Warning: Empty sequence for flank stretch. Query sequence: '{query_seq}', Ref sequence: '{ref_seq}'.")
        return None

    # Align reference sequence (ref_start:ref_end) to query sequence (query_start:query_end)
    # Use parasail for small sequences (< 10kb each), edlib for larger sequences
    cigar_string, score = align_sequences(query_seq, ref_seq)

    return {
        'cigar': cigar_string,
        'md_tag': generate_md_tag(cigar_string, ref_seq, query_seq),
        'cs_tag': generate_cs_tag(query_seq, ref_seq, cigar_string),
        'score': score,
        'ref_start': ref_start,
        'ref_end': ref_end,
        'query_start': query_start,
        'query_end': query_end,
        'query_length': len(query_seq),
        'ref_length': len(ref_seq),
        'flank_type': flank_type,
        'gap': False
    }


def process_flank_stretches(cooper, read, softclip_loci_coords):
    """
    Process and align individual flank stretches from softclipped regions.

    :param cooper: Cooper object containing reference
    :param read: pysam read object
    :param softclip_loci: list of tuples (locus_start, locus_end)
    :param softclip_loci_coords: list of dicts with flank coordinates
    :return: dict with processed flanks and updated CIGAR info
    """

    result = {'upstream': [], 'downstream': []}
    prev_ref_end = 0
    prev_query_end = 0

    # in cs tag, for N use '~' and for B use 'B' to distinguish from insertions and deletions, which are represented by + and - respectively


    if len(softclip_loci_coords) == 0:
        return

    prev_sofclip_dir = None
    read_ref_start   = read.ref_start
    read_ref_end     = read.ref_end
    read_query_start = read.query_start
    read_query_end   = read.query_end
    for coords in softclip_loci_coords:
        ref_start   = coords['upstream'][0]
        ref_end     = coords['downstream'][1]
        query_start = coords['upstream'][2]
        query_end   = coords['downstream'][3]

        if ref_start < read_ref_start: read_ref_start = ref_start
        if ref_end > read_ref_end: read_ref_end = ref_end
        if query_start < read_query_start: read_query_start = query_start
        if query_end > read_query_end: read_query_end = query_end

        softclip_dir = 'upstream' if ref_end < read.ref_end else 'downstream'
        coords    = (ref_start, ref_end, query_start, query_end)
        alignment = align_flank_stretch(cooper, read, coords, softclip_dir)
        if softclip_dir == 'upstream' and prev_query_end == 0 and query_start > 0:
            result['upstream'].append({'cigar': f'{query_start}S', 'md_tag': '', 'cs_tag': '', 'gap': True, 'flank_type': 'upstream'})
        if prev_sofclip_dir == softclip_dir:
            gap = False
            gap_cigar = ''
            gap_md    = ''
            gap_cs    = ''
            if prev_ref_end > 0:
                gap_length = ref_start - prev_ref_end
                if gap_length > 0:
                    gap = True
                    gap_cigar += f'{gap_length}N'
                    gap_cs    += f'~{ref_start-prev_ref_end}'
                    gap_md    += f'^{cooper.ref.fetch(read.chrom, prev_ref_end, ref_start)}'

            if prev_query_end > 0 and query_start > prev_query_end:
                gap_length = query_start - prev_query_end
                if gap_length > 0:
                    gap = True
                    gap_cigar += f'{gap_length}B'
                    gap_cs    += f'B{gap_length}'

            if gap:
                result[softclip_dir].append({'gap': True, 'cigar': gap_cigar, 'md_tag': gap_md,
                                             'cs_tag': gap_cs, 'flank_type': softclip_dir})

        result[softclip_dir].append(alignment)

        prev_ref_end = ref_end
        prev_query_end = query_end
        prev_sofclip_dir = softclip_dir

    upstream_cigar = ''; upstream_ref_end = 0; upstream_query_end = 0
    for aln in result['upstream']:
        if 'ref_end' in aln:
            upstream_ref_end = aln['ref_end']
            upstream_query_end = aln['query_end']
        if upstream_cigar == '':
            upstream_cigar = aln['cigar']
        else:
            upstream_cigar = join_cigars(upstream_cigar, aln['cigar'])

    read_cigar = read.cigarstring
    upstream_md = ''; downstream_md = ''
    upstream_cs = ''; downstream_cs = ''

    if upstream_cigar != '':

        if upstream_ref_end < read.ref_start or upstream_query_end < read.query_start:
            ref_gap   = read.ref_start - upstream_ref_end
            query_gap = read.query_start - upstream_query_end
            sub_cigar = f'{ref_gap}N' if ref_gap > 0 else ''
            sub_cigar += f'{query_gap}B' if query_gap > 0 else ''
            upstream_cigar = join_cigars(upstream_cigar, sub_cigar)
            if ref_gap > 0:
                result['upstream'].append({'gap': True, 'cigar': '', 'md_tag': f'^{cooper.ref.fetch(read.chrom, upstream_ref_end, read.ref_start)}',
                                             'cs_tag': f'~{ref_gap}', 'flank_type': softclip_dir})
            if query_gap > 0:
                result['upstream'].append({'gap': True, 'cigar': '', 'md_tag': '',
                                             'cs_tag': f'B{query_gap}', 'flank_type': softclip_dir})

        if read.has_tag('MD'):
            for aln in result['upstream']:
                if 'md_tag' in aln:
                    if upstream_md == '':
                        upstream_md = aln['md_tag']
                    else:
                        upstream_md = join_mdtags(upstream_md, aln['md_tag'])
        if read.has_tag('cs'):
            for aln in result['upstream']:
                if 'cs_tag' in aln:
                    if upstream_cs == '':
                        upstream_cs = aln['cs_tag']
                    else:
                        upstream_cs = join_cstags(upstream_cs, aln['cs_tag'])

        if 'X' in read.cigarstring:
            if '=' in read.cigarstring: pass
            else: upstream_cigar = upstream_cigar.replace('=', 'M')
        else:
            if 'M' in read.cigarstring:
                upstream_cigar = upstream_cigar.replace('=', 'M')
                upstream_cigar = _collapse_mismatches(upstream_cigar, 'M')
            elif '=' in read.cigarstring:
                upstream_cigar = _collapse_mismatches(upstream_cigar, '=')
        read_cigar = join_cigars(upstream_cigar, strip_softclip(read_cigar, 'left'))

    downstream_cigar = ''
    downstream_ref_start = 0; downstream_query_start = 0
    downstream_query_end = 0
    for i, aln  in enumerate(result['downstream']):
        if i == 0:
            if 'ref_start' in aln:
                downstream_ref_start = aln['ref_start']
                downstream_query_start = aln['query_start']
        if 'query_end' in aln:
            downstream_query_end = aln['query_end']
        if downstream_cigar == '':
            downstream_cigar = aln['cigar']
        else:
            downstream_cigar = join_cigars(downstream_cigar, aln['cigar'])
    if downstream_cigar != '':

        if downstream_ref_start > read.ref_end or downstream_query_start > read.query_end:
            ref_gap   = downstream_ref_start - read.ref_end
            query_gap = downstream_query_start - read.query_end
            if query_gap > 0:
                result['downstream'] = [{'gap': True, 'cigar': '', 'md_tag': '',
                                             'cs_tag': f'B{query_gap}', 'flank_type': softclip_dir}] + result['downstream']
            if ref_gap > 0:
                result['downstream'] = [{'gap': True, 'cigar': '', 'md_tag': f'^{cooper.ref.fetch(read.chrom, read.ref_end, downstream_ref_start )}',
                                             'cs_tag': f'~{ref_gap}', 'flank_type': softclip_dir}] + result['downstream']
            sub_cigar = f'{ref_gap}N' if ref_gap > 0 else ''
            sub_cigar += f'{query_gap}B' if query_gap > 0 else ''
            downstream_cigar = join_cigars(sub_cigar, downstream_cigar)

        if read.has_tag('MD'):
            for aln in result['downstream']:
                if 'md_tag' in aln:
                    if downstream_md == '':
                        downstream_md = aln['md_tag']
                    else:
                        downstream_md = join_mdtags(downstream_md, aln['md_tag'])
        if read.has_tag('cs'):
            for aln in result['downstream']:
                if 'cs_tag' in aln:
                    if downstream_cs == '':
                        downstream_cs = aln['cs_tag']
                    else:
                        downstream_cs = join_cstags(downstream_cs, aln['cs_tag'])

        if 'X' in read.cigarstring:
            if '=' in read.cigarstring: pass
            else: downstream_cigar = downstream_cigar.replace('=', 'M')
        else:
            if 'M' in read.cigarstring:
                downstream_cigar = downstream_cigar.replace('=', 'M')
                downstream_cigar = _collapse_mismatches(downstream_cigar, 'M')
            elif '=' in read.cigarstring:
                downstream_cigar = _collapse_mismatches(downstream_cigar, '=')
        read_cigar = join_cigars(strip_softclip(read_cigar, 'right'), downstream_cigar)
        if downstream_query_end < len(read.query_sequence):
            read_cigar += f'{len(read.query_sequence) - downstream_query_end}S'

    # if _query_length(_cigar_tuples(read.cigarstring)) != _query_length(_cigar_tuples(read_cigar)):
    #     raise ValueError("Query length after processing flank stretches does not match read sequence length.\n" +
    #                     f"Length from CIGAR: {(_query_length(_cigar_tuples(read_cigar)))}, " +
    #                     f"Read sequence length: {len(read.query_sequence)}")

    read.cigarstring = read_cigar
    read.cigartuples = _cigar_tuples(read_cigar)
    read.ref_start   = read_ref_start
    read.ref_end     = read_ref_end
    read.query_start = read_query_start
    read.query_end   = read_query_end
    if read.has_tag('MD'):
        read.md_tag = join_mdtags(upstream_md, read.md_tag)
        read.md_tag = join_mdtags(read.md_tag, downstream_md)
        # cannot validate MD tag because of introduction of N and B operations which are not represented in MD tag, so we skip validation for now
        # if not validate_md_tag(read.md_tag, read.cigarstring).is_valid:
        #     raise ValueError(f"Invalid MD tag for read: {read.query_name} with reference start position: {read.ref_start}")
    if read.has_tag('cs'):
        read.cs_tag = join_cstags(upstream_cs, read.cs_tag)
        read.cs_tag = join_cstags(read.cs_tag, downstream_cs)
        # if not validate_cs_cigar(read.cs_tag, read.cigarstring).is_valid:
        #     raise ValueError(f"Invalid CS tag for read: {read.query_name} with reference start position: {read.ref_start}")
