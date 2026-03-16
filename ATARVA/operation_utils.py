from ATARVA.process_softclips import detect_flank


def clean_eqsign_readseq(chrom, ref_pos, cigar_tuples, read_seq, ref):
    """
    Convert the read sequence to replace "=" with actual bases

    :param chrom: chromosome of the read
    :param ref_pos: reference position of the read
    :param cigar_tuples: list of cigar tuples for the read
    :param read_seq: original read sequence with '=' placeholders
    :param ref: reference fasta file object
    :return: read sequence with actual bases instead of '=' placeholders
    """

    new_read_seq = []
    qpos = 0

    CIGAR_ACTIONS = {4, 1, 2, 0, 7, 8}

    for op, length in cigar_tuples:
        if op in (4, 1):                                         # S, I — copy from read
            new_read_seq.append(read_seq[qpos:qpos + length])
            qpos += length

        elif op == 2:                                            # D — advance ref only
            ref_pos += length

        elif op in (0, 7):                                       # M, = — fetch from ref
            ref_end    = ref_pos + length
            inter_ref  = ref.fetch(chrom, ref_pos, ref_end)

            if op == 7:                                          # = exact match, use ref
                new_read_seq.append(inter_ref)
            else:                                                # M resolve '=' placeholders
                segment = read_seq[qpos:qpos + length]
                new_read_seq.append(''.join(
                    inter_ref[i] if base == '=' else base
                    for i, base in enumerate(segment)
                ))

            ref_pos  = ref_end
            qpos       += length

        elif op == 8:                                            # X — mismatch, use read
            new_read_seq.append(read_seq[qpos:qpos + length])
            ref_pos += length
            qpos    += length

    return ''.join(new_read_seq)


def record_homopolymers(ref_seq, locus_start, homopolymer_positions, threshold=3):
    """
    Update the homopolymer positions for a given reference sequence and locus start position
    For each position part of the homopolymer, the length of the homopolymer from the at position is stored

    :param ref_seq: reference sequence for the locus
    :param locus_start: start position of the locus in the reference
    :param homopolymer_positions: dictionary to store the homopolymer positions
    :param threshold: minimum length of homopolymer to be considered (default is 3)    
    """

    i = 0
    while i < len(ref_seq):
        j = i
        # extend j while same base
        while j < len(ref_seq) and ref_seq[j] == ref_seq[i]:
            j += 1

        length = j - i
        if length >= threshold:
            for offset in range(length):
                homopolymer_positions[locus_start + i + offset] = length - offset


def match_jump(read, match_refend, match_qryend, match_length, repeat_index, locus_query_range,
               flank_query_range, locus_reached, locus_boundary_crossed):
    """
    process a match and update data for affected repeats

    :param read: the read being processed
    :param match_refend: the reference end position of the match
    :param match_qryend: the query end position of the match
    :param match_length: length of the match
    :param repeat_index: the current repeat index being processed
    :param locus_query_range: list of query position ranges for each locus in the read
    :param flank_query_range: list of query position ranges for the flanks of each locus in the read
    :param locus_reached: list of bools storing if reached the locus scanning the cigar
    :param locus_boundary_crossed: bools indicating if the repeat boundaries have been crossed scanning the cigar
    :return: number of repeat indices to jump after processing the match segment
    """

    match_refstart = match_refend - match_length
    r = 0 
    for r, coord in enumerate(read.loci_coords[repeat_index:]):
        current_index = r + repeat_index
        flank_start, flank_end = coord
        locus_start = flank_start + read.left_flanks[current_index]
        locus_end = flank_end - read.right_flanks[current_index]

        # if the match segment is before the start of the repeat; repeat is unaffected
        if match_refend < flank_start: break

        # if the match segment is beyond the end of the repeat; repeat is unaffected
        if match_refstart > flank_end: continue
            
        if not locus_reached[current_index]:

            if flank_start <= match_refend:
                flank_query_range[current_index][0] = match_qryend - (match_refend - flank_start)

            if flank_end <= match_refend:
                flank_query_range[current_index][1] = match_qryend - (match_refend - flank_end)

            locus_reached[current_index] = True

            if not locus_boundary_crossed[current_index][0]:
                if locus_start <= match_refend:
                    locus_query_range[current_index][0] = match_qryend - (match_refend - locus_start)
                    locus_boundary_crossed[current_index][0] = True
            if not locus_boundary_crossed[current_index][1]:
                if locus_end <= match_refend:
                    locus_query_range[current_index][1] = match_qryend - (match_refend - locus_end)
                    # the postion is greater than the end to consider inserts at the end of the repeat
                    if match_refend > locus_end: locus_boundary_crossed[current_index][1] = True

        elif flank_end <= match_refend:
            flank_query_range[current_index][1] = match_qryend - (match_refend - flank_end)

        # for storing repeat qpos ranges
        if not locus_boundary_crossed[current_index][0]:
            if locus_start <= match_refend:
                locus_query_range[current_index][0] = match_qryend - (match_refend - locus_start)
                locus_boundary_crossed[current_index][0] = True
            if locus_end <= match_refend:
                locus_query_range[current_index][1] = match_qryend - (match_refend - locus_end)
                if match_refend > locus_end: locus_boundary_crossed[current_index][1] = True

        elif not locus_boundary_crossed[current_index][1]:
            if locus_end <= match_refend:
                locus_query_range[current_index][1] = match_qryend - (match_refend - locus_end)
                if match_refend > locus_end: locus_boundary_crossed[current_index][1] = True

    jump = 0    # jump beyond the repeat where all positions are tracked
    if read.loci_coords[repeat_index + r - 1][1] < match_refend:
        for f in read.loci_coords[repeat_index:]:
            if f[1] < match_refend: jump += 1
            else: break

    return jump


def deletion_jump(read, del_refend, qpos, deletion_len, repeat_index, locus_query_range, flank_query_range, locus_reached, locus_boundary_crossed):
    """
    process a deletion and update data for affected repeats

    :param read: the read being processed
    :param del_refend: the reference end position of the deletion
    :param qpos: the query position corresponding to the reference end position of the deletion
    :param deletion_len: length of the deletion
    :param repeat_index: the current repeat index being processed
    :param locus_query_range: list of query position ranges for each locus in the read
    :param flank_query_range: list of query position ranges for the flanks of each locus in the read
    :param locus_reached: list of bools storing if reached the locus scanning the cigar
    :param locus_boundary_crossed: bools indicating if the repeat boundaries have been crossed scanning the cigar
    :return: number of repeat indices to jump after processing the deletion segment
    """

    del_refstart = del_refend - deletion_len

    r = 0   # required to be initialised outside the loop
    for r, coord in enumerate(read.loci_coords[repeat_index:]):
        current_index = r + repeat_index

        flank_start, flank_end = coord  # the locus coordinates include flanks
        locus_start = flank_start + read.left_flanks[current_index]
        locus_end = flank_end - read.right_flanks[current_index]

        # if rpos is before the start of the repeat; repeat is unaffected
        if del_refend < flank_start: break

        # actual position in the reference where the deletion is occurring
        if del_refstart > flank_end: continue

        locus_key = read.loci_keys[current_index]
        if not locus_reached[current_index]:
            # if the locus is not tracked so far
            if locus_start <= del_refend:    
                flank_query_range[current_index][0] = qpos        
                locus_reached[current_index] = True    # set tracked as true

            if locus_end < del_refend:
                flank_query_range[current_index][1] = qpos

            # for storing repeat qpos ranges
            if not locus_boundary_crossed[current_index][0]:
                if locus_start <= del_refend:
                    locus_query_range[current_index][0] = qpos
                    locus_boundary_crossed[current_index][0] = True
            if not locus_boundary_crossed[current_index][1]:
                # should probably be locus_end < del_refend because locus_end includes flank
                if locus_end < del_refend:
                    locus_query_range[current_index][1] = qpos
                    locus_boundary_crossed[current_index][1] = True

        elif flank_end < del_refend:
            flank_query_range[current_index][1] = qpos

        # for storing repeat qpos ranges
        if not locus_boundary_crossed[current_index][0]:
            if locus_start <= del_refend:
                locus_query_range[current_index][0] = qpos
                locus_boundary_crossed[current_index][0] = True
            if locus_end < del_refend:
                locus_query_range[current_index][1] = qpos
                locus_boundary_crossed[current_index][1] = True
        elif not locus_boundary_crossed[current_index][1]:
            if locus_end <= del_refend:
                locus_query_range[current_index][1] = qpos
                if del_refend > locus_end: locus_boundary_crossed[current_index][1] = True

        # if deletion does not overlap the repeat; skip
        if del_refend < locus_start or del_refstart > locus_end:
            continue
        else:
            # updating the allele with the deletion considered
            del_len = min(locus_end, del_refend) - max(locus_start, del_refstart)
            if del_refstart not in read.homopolymer_positions:
                read.loci_data[locus_key].alen  -= del_len
                read.loci_data[locus_key].halen -= del_len
            else:
                # if the deletion is only limited to the homopolymer positions
                if del_len <= read.homopolymer_positions[del_refstart]:
                    read.loci_data[locus_key].halen -= del_len
                else:
                    read.loci_data[locus_key].alen  -= del_len
                    read.loci_data[locus_key].halen -= del_len

    jump = 0    # jump beyond the repeat where all positions are tracked
    if read.loci_coords[repeat_index + r - 1][1] < del_refend:
        for coord in read.loci_coords[repeat_index:]:
            if coord[1] < del_refend: jump += 1
            else: break

    return jump


def insertion_jump(read, ins_refpos, ins_qend, insert_len, homopolymer_insert, repeat_index, locus_query_range, flank_query_range,
                   locus_reached, locus_boundary_crossed, left_flank_insertions, right_flank_insertions):
    """
    process an insertion and update data for affected repeats
    :param read: the read being processed
    :param ins_refpos: the reference position where the insertion is occurring
    :param ins_qend: the query end position of the insertion
    :param insert_len: length of the insertion
    :param homopolymer_insert: boolean indicating if the insertion is a homopolymer    :param repeat_index: the current repeat index being processed
    :param locus_query_range: list of query position ranges for each locus in the read
    :param flank_query_range: list of query position ranges for the flanks of each locus in the read
    :param locus_reached: list of bools storing if reached the locus scanning the cigar
    :param locus_boundary_crossed: bools indicating if the repeat boundaries have been crossed scanning the cigar
    :param left_flank_insertions: list of insertions in the left flank for each locus
    :param right_flank_insertions: list of insertions in the right flank for each locus
    :return: number of repeat indices to jump after processing the insertion segment
    """

    ins_qstart = ins_qend - insert_len  # start position of the insertion in the query

    r = 0   # required to be initialised outside the loop
    for r, coord in enumerate(read.loci_coords[repeat_index:]):
        current_index = r + repeat_index

        flank_start, flank_end = coord
        locus_start = flank_start + read.left_flanks[current_index]
        locus_end = flank_end - read.right_flanks[current_index]
        # if rpos is before the start of the repeat; repeat is unaffected
        if ins_refpos < flank_start: break

        # if the insertion is happening beyond, the repeat in unaffected
        if ins_refpos > flank_end: continue

        locus_key = read.loci_keys[current_index]
        if not locus_reached[current_index]:
            if flank_start <= ins_refpos:      # can this be locus_start - 1??
                flank_query_range[current_index][0] = ins_qstart
                locus_reached[current_index] = True    # set tracked as true

            if flank_end == ins_refpos:
                flank_query_range[current_index][1] = ins_qend

            # for storing repeat qpos ranges
            if not locus_boundary_crossed[current_index][0]:
                if locus_start - 1 <= ins_refpos:
                    locus_query_range[current_index][0] = ins_qstart
                    locus_boundary_crossed[current_index][0] = True
            if not locus_boundary_crossed[current_index][1]:
                if locus_end <= ins_refpos:
                    locus_query_range[current_index][1] = ins_qend
                    if ins_refpos > flank_end: locus_boundary_crossed[current_index][1] = True


        elif flank_end == ins_refpos:
            flank_query_range[current_index][1] = ins_qend

        # for storing repeat qpos ranges
        if not locus_boundary_crossed[current_index][0]:
            if locus_start <= ins_refpos:
                locus_query_range[current_index][0] = ins_qstart
                locus_boundary_crossed[current_index][0] = True
            if locus_end <= ins_refpos:
                locus_query_range[current_index][1] = ins_qend
                if ins_refpos > locus_end: locus_boundary_crossed[current_index][1] = True
        elif not locus_boundary_crossed[current_index][1]:
            if locus_end <= ins_refpos:
                locus_query_range[current_index][1] = ins_qend
                if ins_refpos > locus_end: locus_boundary_crossed[current_index][1] = True

        # read_loci_variations[locus_key][rpos] = f'I|{insertion_length}'
        if locus_start <= ins_refpos <= locus_end: # introduced to include length only if it comes inside repeat region
            if ins_refpos not in read.homopolymer_positions:
                read.loci_data[locus_key].alen  += insert_len
                read.loci_data[locus_key].halen += insert_len
            else:
                if homopolymer_insert:
                    # only if the insertion is a homopolymer; consider it as homopolymer insertion
                    read.loci_data[locus_key].halen += insert_len
                else:
                    read.loci_data[locus_key].alen  += insert_len
                    read.loci_data[locus_key].halen += insert_len

        if flank_start <= ins_refpos <= locus_start - 1: # -1 is included so ins near the start pos is not taken into account as it is already added
            try:
                left_flank_insertions[current_index].append((ins_refpos, ins_qstart, ins_qend))
            except AttributeError:
                pass
        elif locus_end + 1 <= ins_refpos <= flank_end: # +1 is included so ins near the end pos is not taken into account as it is already added
            try:
                right_flank_insertions[current_index].append((ins_refpos, ins_qstart, ins_qend))
            except AttributeError:
                pass

    jump = 0    # jump beyond the repeat where all positions are tracked
    if read.loci_coords[repeat_index + r - 1][1] < ins_refpos:
        for f in read.loci_coords[repeat_index:]:
            if f[1] < ins_refpos: jump += 1
            else: break

    return jump
