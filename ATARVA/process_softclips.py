import parasail
import cstag
import numpy as np

from ATARVA.realignment import *


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
                md.append(reference[rpos])
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

    operations = {':', '-', '+', '*', '=', '~'}
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
        elif left_op == ':':
            merged_opdat = str(int(left_opdat) + int(right_opdat))
        else:
            merged_opdat = left_opdat + right_opdat
        return cs_left[:-len(left_opdat)-1] + left_op + merged_opdat + cs_right[len(right_opdat)+1:]
    else:
        return cs_left + cs_right


# def join_mdtags(left_md: str, right_md: str) -> str:
#     """
#     join MD tags of continuous alignments

#     Examples:
#         ["10", "2", "A", "0", "^C", "^T", "3"] -> "12A0^CT3"
#         "10A00^C2" -> "10A0^C2"
    
#     :param left_md:  left MD tag fragment (e.g., "10A0^C2")
#     :param right_md: right MD tag fragment (e.g., "12T5")
#     :return:         joined MD tag (e.g., "10A0^C2T5")
#     """

#     left_op = ''
#     left_opdat = ''
#     left_index = len(left_md)
#     for x in left_md[::-1]:
#         if x.isdigit():
#             if left_op == '':
#                 left_op = 'M'
#                 left_opdat = x + left_opdat
#             else: break
#         elif x == '^':
#             left_op = 'D'
#             left_index -= 1
#             break
#         elif x in 'ACGTN':
#             if left_op == 'M':
#                 break
#             else:
#                 left_op = 'X'
#                 left_opdat = x + left_opdat
#         left_index -= 1

#     right_op = ''
#     right_opdat = ''
#     right_index = 0
#     for x in right_md:
#         if x.isdigit():
#             if right_op == '':
#                 right_op = 'M'
#                 right_opdat += x
#             else: break
#         elif x == '^':
#             if right_op == 'M':
#                 break
#             else: right_op = 'D'
#         elif x in 'ACGTN':
#             if right_op == '':
#                 right_op = 'X'
#                 right_opdat += x
#             elif right_op == 'D':
#                 right_opdat += x
#             elif right_op == 'M':
#                 break
#         right_index += 1

#     if left_op == right_op:
#         if left_op == 'M':
#             return left_md[:left_index] + str(int(left_opdat) + int(right_opdat)) + right_md[right_index:]
#         elif left_op == 'X':
#             return left_md + '0' + right_md
#         elif left_op == 'D':
#             return left_md[:left_index] + '^' + left_opdat + right_opdat + right_md[right_index:]
#     else:
#         if left_op == 'D' and right_op == 'X':
#             return left_md + '0' + right_md
#         else:
#             return left_md + right_md

def tokenize_md(md):
    tokens = []
    i = 0
    while i < len(md):
        if md[i].isdigit():
            num = 0
            while i < len(md) and md[i].isdigit():
                num = num * 10 + int(md[i])
                i += 1
            tokens.append(num)
        elif md[i] == '^':
            i += 1
            seq = []
            while i < len(md) and md[i].isalpha():
                seq.append(md[i])
                i += 1
            tokens.append('^' + ''.join(seq))
        else:
            tokens.append(md[i])
            i += 1
    return tokens


def join_mdtags(md1, md2):
    t1 = tokenize_md(md1)
    t2 = tokenize_md(md2)

    if not t1:
        return md2
    if not t2:
        return md1

    # merge boundary numbers
    if isinstance(t1[-1], int) and isinstance(t2[0], int):
        t1[-1] += t2[0]
        t2 = t2[1:]

    merged = t1 + t2

    # rebuild string
    out = []
    for t in merged:
        out.append(str(t))
    return ''.join(out)


def join_cigars(cigar_left: str, cigar_right: str) -> str:
    """
    join two continuous CIGAR strings and normalize.
    :param cigar_left:  left CIGAR string (e.g., "10M1I5M")
    :param cigar_right: right CIGAR string (e.g., "2M3D8M")
    :return:            joined CIGAR string (e.g., "10M1I7M3D8M")
    """

    operations = {'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X'}
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

            elif c in ('I', 'D'):
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

    flank = 30

    if stream: # if its upstream (True)

        # finding the upstream flank in the read 
        upstream = cooper.ref.fetch(read.chrom, locus_start - flank, locus_start)
        alignment_score, target_begin, target_end, query_begin, query_end, sCigar = stripSW(Inputs(read.query_sequence, upstream), False)

        # target_begin is the start of the alignment in the read
        # query_begin  is the start of the alignment in the reference flank
        # target_end   is the end of the alignment in the read
        # query_end    is the end of the alignment in the reference flank

        # if the upstream flank aligns to the already aligned region(by aligner), there is no expansion in the soft-clip, here target_end is dummy
        if target_end > read.query_start: return [-1, -1]

        # reject read if less than 70% match with the reference flank sequence
        if alignment_score < int(2*(0.8*flank)): return [-1, -1]

        new_read_query_start = target_begin - 1
        new_read_ref_start   = locus_start - flank + query_begin - 1

        return [new_read_query_start, new_read_ref_start]
        
    else: # if its downstream (False)

        # finding the downstream flank in the read
        downstream = cooper.ref.fetch(read.chrom, locus_end, locus_end + flank)
        alignment_score, target_begin, target_end, query_begin, query_end, sCigar = stripSW(Inputs(read.query_sequence, downstream), False)

        # target_begin is the start of the alignment in the read
        # query_begin  is the start of the alignment in the reference flank
        # target_end   is the end of the alignment in the read
        # query_end    is the end of the alignment in the reference flank

        # if the downstream flank aligns to the already aligned region(by aligner), there is no expansion in the soft-clip, here target_begin is dummy
        if target_begin <= read.query_end: return [-1, -1]

        # reject read if less than 70% match with the reference flank sequence
        if alignment_score < int(2*(0.8*flank)): return [-1, -1]

        new_read_query_end = target_end - 1
        new_read_ref_end   = locus_end + query_end - 1

        return [new_read_query_end, new_read_ref_end]


def process_softclips(cooper, read, locus, dir):
    """
    process the soft-clipped regions of the read to process missed repeats
    :param cooper: Cooper object containing reference and other data
    :param read: pysam AlignedSegment object
    :param softclip_loci: list of tuples (locus_start, locus_end) for the soft-clipped regions to be processed

    :return: None
             Modifies the read object in place
             - updates the read cigar
             - updates read reference start and end positions
             - updates read query     start and end positions
    """

    # start_locus = sorted(softclip_loci, key=lambda x: x[0])[0]
    # end_locus = sorted(softclip_loci, key=lambda x: x[1])[-1]
    default_flank = 30

    if dir == "left":
        if read.ref_start > locus[0] + default_flank:
            # detect the flank for the first upstream locus
            new_read_query_start, new_read_ref_start = detect_flank(cooper, read, locus[0], locus[1], True)
            if new_read_query_start == -1: return
            unaligned_query = read.query_sequence[new_read_query_start:read.query_start]
            unaligned_ref   = cooper.ref.fetch(read.chrom, new_read_ref_start, read.ref_start)
            sub_cigar, alignment_score = nw_align(unaligned_query, unaligned_ref)
            aligned_matches   = _match_lengths(read.cigarstring)
            softclip_matches  = _match_lengths(sub_cigar)
            if np.mean(softclip_matches)/np.mean(aligned_matches) >= 0.6 or np.mean(softclip_matches) >= 20:
                sub_mdtag = generate_md_tag(sub_cigar, unaligned_ref, unaligned_query)
                # update read positions and cigar
                if read.has_tag('cs'):
                    sub_cs_tag = cstag.call(sub_cigar, sub_mdtag, unaligned_query)
                    read.cs_tag = join_cstags(sub_cs_tag, read.cs_tag)
                if read.has_tag('MD') and read._read.has_tag('MD'):
                    read.md_tag = join_mdtags(sub_mdtag, read.md_tag)

                read.query_start = new_read_query_start
                read.ref_start   = new_read_ref_start
                if 'X' in read.cigarstring:
                    if '=' in read.cigarstring: pass
                    else: sub_cigar = sub_cigar.replace('=', 'M')
                else:
                    if 'M' in read.cigarstring:
                        sub_cigar = sub_cigar.replace('=', 'M')
                        sub_cigar = _collapse_mismatches(sub_cigar, 'M')
                    elif '=' in read.cigarstring:
                        sub_cigar = _collapse_mismatches(sub_cigar, '=')
                read.cigarstring = join_cigars(sub_cigar, strip_softclip(read.cigarstring, 'left'))
                if new_read_query_start > 0:
                    read.cigarstring = f'{new_read_query_start}S'  + read.cigarstring
                read.cigartuples = _cigar_tuples(read.cigarstring)
                if _query_length(read.cigartuples) != len(read.query_sequence):
                    raise ValueError("Query length after softclip processing does not match read sequence length.\n" +
                                    f"Length from CIGAR: {(_query_length(read.cigartuples))}, " +
                                    f"Read sequence length: {len(read.query_sequence)}")

                if read.has_tag('MD') and not valid_md(read.md_tag, read.cigarstring):
                    raise ValueError("MD tag is not consistent with CIGAR string after softclip processing.\n" +
                                    f"MD stats: {md_stats(read.md_tag)}\nCIGAR stats: {cigar_stats(read.cigarstring)}")

            return True

    if dir == "right":
        if read.ref_end < locus[1] - default_flank:
            # detect the flank for the first downstream locus
            new_read_query_end, new_read_ref_end = detect_flank(cooper, read, locus[0], locus[1], False)
            if new_read_query_end == -1: return
            unaligned_query = read.query_sequence[read.query_end:new_read_query_end]
            unaligned_ref   = cooper.ref.fetch(read.chrom, read.ref_end, new_read_ref_end)
            sub_cigar, alignment_score = nw_align(unaligned_query, unaligned_ref)
            aligned_matches   = _match_lengths(read.cigarstring)
            softclip_matches  = _match_lengths(sub_cigar)
            if np.mean(softclip_matches)/np.mean(aligned_matches) >= 0.6 or np.mean(softclip_matches) >= 20:
                sub_mdtag = generate_md_tag(sub_cigar, unaligned_ref, unaligned_query)
                # update read positions and cigar
                if read.has_tag('cs') and read._read.has_tag('cs'):
                    sub_cs_tag = cstag.call(sub_cigar, sub_mdtag, unaligned_query)
                    read.cs_tag = join_cstags(read.cs_tag, sub_cs_tag)
                if read.has_tag('MD') and read._read.has_tag('MD'):
                    read.md_tag = join_mdtags(read.md_tag, sub_mdtag)
                    # read.md_tag = generate_md_tag(read.cigarstring, cooper.ref.fetch(read.chrom, read.ref_start, read.ref_end), read.query_sequence)

                read.query_end = new_read_query_end
                read.ref_end   = new_read_ref_end
                if 'X' in read.cigarstring:
                    if '=' in read.cigarstring: pass
                    else: sub_cigar = sub_cigar.replace('=', 'M')
                else:
                    if 'M' in read.cigarstring:
                        sub_cigar = sub_cigar.replace('=', 'M')
                        sub_cigar = _collapse_mismatches(sub_cigar, 'M')
                    elif '=' in read.cigarstring:
                        sub_cigar = _collapse_mismatches(sub_cigar, '=')
                read.cigarstring = join_cigars(strip_softclip(read.cigarstring, 'right'), sub_cigar)
                if new_read_query_end < len(read.query_sequence):
                    read.cigarstring += f'{len(read.query_sequence) - new_read_query_end}S'
                read.cigartuples = _cigar_tuples(read.cigarstring)
                if _query_length(read.cigartuples) != len(read.query_sequence):
                    raise ValueError("Query length after softclip processing does not match read sequence length.\n" +
                                    f"Length from CIGAR: {(_query_length(read.cigartuples))}, " +
                                    f"Read sequence length: {len(read.query_sequence)}")

                if _ref_length(read.cigartuples) != (read.ref_end - read.ref_start):
                    raise ValueError("Reference length after softclip processing does not match read reference span.\n" +
                                    f"Length from CIGAR: {(_ref_length(read.cigartuples))}, " +
                                    f"Read reference span: {(read.ref_end - read.ref_start)}")

                if read.has_tag('MD') and not valid_md(read.md_tag, read.cigarstring):
                    raise ValueError("MD tag is not consistent with CIGAR string after softclip processing.\n" +
                                    f"MD stats: {md_stats(read.md_tag)}\nCIGAR stats: {cigar_stats(read.cigarstring)}")

                return True
    return False
