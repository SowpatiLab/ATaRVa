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
               flank_query_range, repeat_tracked, flank_tracked):
            #    qpos, match_len, loci_flank_qpos_range, flank_track, left_flank,
            #    right_flank, amp_right_flank_list, amp_left_flank_list,
            #    out_insertion_qpos_ranges_right, out_insertion_qpos_ranges_left,
            #    right_ins_rpos, left_ins_rpos, amplicon_variables):
    """
    Return the number of repeat indices to jump when scanning through a match segment
    """

    previous_rpos = match_refend - match_length
    r = 0 
    for r, coord in enumerate(read.loci_coords[repeat_index:]):
        coord_start, coord_end = coord
        current_index = r + repeat_index

        # if the match segment is before the start of the repeat; repeat is unaffected
        if match_refend < coord_start: break

        # if the match segment is beyond the end of the repeat; repeat is unaffected
        if previous_rpos > coord_end: continue
            
        locus_key = read.loci_keys[current_index]
        if not repeat_tracked[current_index]:

            if coord_start <= match_refend:
                locus_query_range[current_index][0] = match_qryend - (match_refend - coord_start)

            if coord_end <= match_refend:                
                locus_query_range[current_index][1] = match_qryend - (match_refend - coord_end)

            repeat_tracked[current_index] = True 

            # for storing repeat qpos ranges
            # ????? Recheck this
            if not flank_tracked[current_index][0]:
                if coord_start <= match_refend:
                    flank_query_range[current_index][0] = match_qryend - (match_refend - (coord_start - read.left_flanks[current_index]))
                    flank_tracked[current_index][0] = True
            if not flank_tracked[current_index][1]:
                if coord_end + read.right_flanks[current_index] <= match_refend:
                    flank_query_range[current_index][1] = match_qryend - (match_refend - (coord_end + read.right_flanks[current_index]))
                    if match_refend > coord_end - read.right_flanks[current_index]: flank_tracked[current_index][1] = True

        elif coord_end <= match_refend:
            locus_query_range[current_index][1] = match_qryend - (match_refend - coord_end)

        # for storing repeat qpos ranges
        if not flank_tracked[current_index][0]:
            if coord_start + read.left_flanks[current_index] <= match_refend:
                flank_query_range[current_index][0] = match_qryend - (match_refend - (coord_start - read.left_flanks[current_index]))
                flank_tracked[current_index][0] = True
            if coord_end + read.right_flanks[current_index] <= match_refend:
                flank_query_range[current_index][1] = match_qryend - (match_refend - (coord_end + read.right_flanks[current_index]))
                if match_refend > coord_end - read.right_flanks[current_index]: flank_tracked[current_index][1] = True
            
        elif not flank_tracked[current_index][1]:
            if coord_end + read.right_flanks[current_index] <= match_refend:
                flank_query_range[current_index][1] = match_qryend - (match_refend - (coord_end + read.right_flanks[current_index]))
                if match_refend > coord_end - read.right_flanks[current_index]: flank_tracked[current_index][1] = True

    jump = 0    # jump beyond the repeat where all positions are tracked
    if read.loci_coords[repeat_index + r - 1][1] < match_refend:
        for f in read.loci_coords[repeat_index:]:
            if f[1] < match_refend: jump += 1
            else: break

    return jump


def deletion_jump(cooper, read, del_refend, qpos, deletion_len, repeat_index, locus_query_range, flank_query_range, repeat_tracked, flank_tracked):
    """
    Return the number of repeat indices to jump when scanning through a deletion segment.
    The function tracks specifically if the deletion is segment has complete repeats in them
    or segments of the repeat is deleted.
    """

    # if amp_left_flank_list:
        # chrom, ref, query_sequence, flank_length, qpos_start, qpos_end = amplicon_variables

    # rpos - corresponds to the position in the reference after tracking the deletion
    r = 0   # required to be initialised outside the loop
    for r, coord in enumerate(read.loci_coords[repeat_index:]):
        coord_start, coord_end = coord
        current_index = r + repeat_index
        # if rpos is before the start of the repeat; repeat is unaffected
        if del_refend < coord_start: break

        # actual position in the reference where the deletion is occurring
        del_refstart = del_refend - deletion_len
        if del_refstart > coord_end: continue

        locus_key = read.loci_keys[current_index]
        if not repeat_tracked[current_index]:
            # if the locus is not tracked
            # deletion is encountered beyond
            if coord_start <= del_refend:    
                locus_query_range[current_index][0] = qpos        
                repeat_tracked[current_index] = True    # set tracked as true
                # if amp_left_flank_list:
                #     lstart, lend = ref_repeat(locus_key)
                #     Record_left_out_ins(amp_left_flank_list, r, repeat_index, chrom, ref, query_sequence, flank_length, lstart, lend, locus_qpos_range, loci_flank_qpos_range, out_insertion_qpos_ranges_left, left_ins_rpos, flank_track, qpos_start, qpos_end)

            if coord_end < del_refend:
                locus_query_range[current_index][1] = qpos
                # if amp_left_flank_list:
                #     lstart, lend = ref_repeat(locus_key)
                #     Record_right_out_ins(amp_right_flank_list, r, repeat_index, chrom, ref, query_sequence, flank_length, lstart, lend, locus_qpos_range, loci_flank_qpos_range, out_insertion_qpos_ranges_right, right_ins_rpos, flank_track, qpos_start, qpos_end)

            # for storing repeat qpos ranges
            if not flank_tracked[current_index][0]:
                if coord_start + read.left_flanks[current_index] <= del_refend:
                    flank_query_range[current_index][0] = qpos
                    flank_tracked[current_index][0] = True
            if not flank_tracked[current_index][1]:
                if coord_end - read.right_flanks[current_index] < del_refend:
                    flank_query_range[current_index][1] = qpos
                    flank_tracked[current_index][1] = True

        elif coord_end < del_refend:
            locus_query_range[current_index][1] = qpos
            # if amp_left_flank_list:
            #     lstart, lend = ref_repeat(locus_key)
            #     Record_right_out_ins(amp_right_flank_list, r, repeat_index, chrom, ref, query_sequence, flank_length, lstart, lend, locus_qpos_range, loci_flank_qpos_range, out_insertion_qpos_ranges_right, right_ins_rpos, flank_track, qpos_start, qpos_end)

        # for storing repeat qpos ranges
        if not flank_tracked[current_index][0]:
            if coord_start + read.left_flanks[current_index] <= del_refend:
                flank_query_range[current_index][0] = qpos
                flank_tracked[current_index][0] = True 
            if coord_end-read.right_flanks[current_index] < del_refend:
                flank_query_range[current_index][1] = qpos
                flank_tracked[current_index][1] = True
        elif (not flank_tracked[current_index][1]) and (coord_end - read.right_flanks[current_index] <= del_refend):
            flank_query_range[current_index][1] = qpos
            if del_refend > coord_end - read.right_flanks[current_index]: flank_tracked[current_index][1] = True

        # updating the allele with the deletion considered
        # read_loci_variations[locus_key][rpos] = f'D|{deletion_length}'

        # del_len = min(coord[1], rpos) - max(coord[0], del_pos)
        del_len = min(coord_end - read.right_flanks[current_index], del_refend) - max(coord_start+read.left_flanks[current_index], del_refstart)
        if (del_refstart >= coord_start + read.left_flanks[current_index]) and (del_refend <= coord_end - read.right_flanks[current_index]): # introduced to include length only if it comes inside repeat region
            if del_refstart not in read.homopolymer_positions:
                read.loci_data[locus_key].alen -= del_len
                read.loci_data[locus_key].halen -= del_len
            else:
                if del_len <= read.homopolymer_positions[del_refstart]:
                    # if the deletion is only limited to the homopolymer positions
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


def insertion_jump(cooper, read, ins_refpos, qpos, insert_len, repeat_index, locus_query_range, flank_query_range,
                   repeat_tracked, flank_tracked, out_insertion_qpos_ranges_left, out_insertion_qpos_ranges_right, left_ins_rpos,
                   right_ins_rpos): #, amp_right_flank_list, amp_left_flank_list, amplicon_variables):
    """
    Return the number of repeat indices to jump when scanning through a insertion segment.
    The function tracks specifically if the deletion is segment has complete repeats in them
    or segments of the repeat is deleted.
    """

    # if amp_left_flank_list:
    #     chrom, ref, query_sequence, flank_length, qpos_start, qpos_end = amplicon_variables

    r = 0   # required to be initialised outside the loop
    for r, coord in enumerate(read.loci_coords[repeat_index:]):
        current_index = r + repeat_index

        coord_start, coord_end = coord
        # if rpos is before the start of the repeat; repeat is unaffected
        if ins_refpos < coord_start: break

        # if the insertion is happening beyond, the repeat in unaffected
        if ins_refpos > coord_end: continue

        locus_key = read.loci_keys[current_index]
        if not repeat_tracked[current_index]:
            # if the locus is not tracked
            # deletion is encountered beyond
            if coord_start <= ins_refpos:
                locus_query_range[current_index][0] = qpos - insert_len
                repeat_tracked[current_index] = True    # set tracked as true
                # if amp_left_flank_list:
                #     lstart, lend = ref_repeat(locus_key)
                #     Record_left_out_ins(amp_left_flank_list, r, repeat_index, chrom, ref, query_sequence, flank_length, lstart, lend, locus_qpos_range, loci_flank_qpos_range, out_insertion_qpos_ranges_left, left_ins_rpos, flank_track, qpos_start, qpos_end)

            if coord_end == ins_refpos:
                locus_query_range[current_index][1] = qpos
                # if amp_left_flank_list:
                #     lstart, lend = ref_repeat(locus_key)
                #     Record_right_out_ins(amp_right_flank_list, r, repeat_index, chrom, ref, query_sequence, flank_length, lstart, lend, locus_qpos_range, loci_flank_qpos_range, out_insertion_qpos_ranges_right, right_ins_rpos, flank_track, qpos_start, qpos_end)

                # here jump can be done

            # for storing repeat qpos ranges
            if not flank_tracked[current_index][0]:
                if coord_start + read.left_flanks[current_index] - 1 <= ins_refpos:
                    locus_query_range[current_index][0] = qpos - insert_len
                    flank_tracked[current_index][0] = True
            if not flank_tracked[current_index][1]:
                if coord_end - read.right_flanks[current_index] <= ins_refpos:
                    locus_query_range[current_index][1] = qpos
                    if ins_refpos > coord_end - read.right_flanks[current_index]: flank_tracked[current_index][1] = True


        elif coord_end == ins_refpos:
            locus_query_range[current_index][1] = qpos
            # if amp_left_flank_list:
            #     lstart, lend = ref_repeat(locus_key)
            #     Record_right_out_ins(amp_right_flank_list, r, repeat_index, chrom, ref, query_sequence, flank_length, lstart, lend, locus_qpos_range, loci_flank_qpos_range, out_insertion_qpos_ranges_right, right_ins_rpos, flank_track, qpos_start, qpos_end)

        # for storing repeat qpos ranges
        if not flank_tracked[current_index][0]:
            if coord_start + read.left_flanks[current_index] <= ins_refpos:
                flank_query_range[current_index][0] = qpos - insert_len
                flank_tracked[current_index][0] = True
            if coord_end - read.right_flanks[current_index] <= ins_refpos:
                flank_query_range[current_index][1] = qpos
                if ins_refpos > coord_end - read.right_flanks[current_index]: flank_tracked[current_index][1] = True
        elif (not flank_tracked[current_index][1]) and (coord_end - read.right_flanks[current_index] <= ins_refpos):
            flank_query_range[current_index][1] = qpos
            if ins_refpos > coord_end - read.right_flanks[current_index]: flank_tracked[current_index][1] = True

        # read_loci_variations[locus_key][rpos] = f'I|{insertion_length}'
        if coord_start + read.left_flanks[current_index] <= ins_refpos <= coord_end - read.right_flanks[current_index]: # introduced to include length only if it comes inside repeat region
            if ins_refpos not in read.homopolymer_positions:
                read.loci_data[locus_key].alen  += insert_len
                read.loci_data[locus_key].halen += insert_len
            else:
                if len(set(insert)) == 1:
                    # only if the insertion is a homopolymer; consider it as homopolymer insertion
                    read.loci_data[locus_key].halen += insert_len
                else:
                    read.loci_data[locus_key].alen  += insert_len
                    read.loci_data[locus_key].halen += insert_len

        if coord_start <= ins_refpos <= coord_start + read.left_flanks[current_index]-1: # -1 is included so ins near the start pos is not taken into account as it is already added
            try:
                out_insertion_qpos_ranges_left[current_index].append((qpos-insert_len, qpos))
                left_ins_rpos[current_index].append(ins_refpos)
            except AttributeError:
                pass
        elif coord_end - read.right_flanks[current_index] + 1 <= ins_refpos <= coord_end: # +1 is included so ins near the end pos is not taken into account as it is already added
            try:
                out_insertion_qpos_ranges_right[current_index].append((qpos-insert_len, qpos))
                right_ins_rpos[current_index].append(ins_refpos)
            except AttributeError:
                pass
    jump = 0    # jump beyond the repeat where all positions are tracked
    if read.loci_coords[repeat_index + r - 1][1] < ins_refpos:
        for f in read.loci_coords[repeat_index:]:
            if f[1] < ins_refpos: jump += 1
            else: break

    return jump

def ref_repeat(locus_key):
    lstart = int(locus_key[locus_key.index(':')+1 : locus_key.index('-')])
    lend = int(locus_key[locus_key.index('-')+1:])
    return lstart, lend

def Record_left_out_ins(amp_left_flank_list, r, repeat_index, chrom, ref, query_sequence, flank_length, coord_start, coord_end, locus_qpos_range, loci_flank_qpos_range, out_insertion_qpos_ranges_left, left_ins_rpos, flank_track, qpos_start, qpos_end):
    needed_len = amp_left_flank_list[r+repeat_index]
    if needed_len>0:
        updated_seq_start, status = detect_flank(chrom, query_sequence, ref, flank_length, coord_start, coord_end, qpos_start, qpos_end, True)
        if status:
            locus_qpos_range[r+repeat_index][0] = updated_seq_start # setting the read flank start
            flank_track[r+repeat_index][0] = True 
            loci_flank_qpos_range[r+repeat_index][0] = updated_seq_start # setting the repeat start in read same as its flank start
            out_insertion_qpos_ranges_left[r+repeat_index] = () # restricting the left insertion ranges as empty
            left_ins_rpos[r+repeat_index] = () # restricting the left insertion ref pos as empty
        else:
            out_insertion_qpos_ranges_left[r+repeat_index] = (None,) # tagging the read as None to ignore this read for genotyping, as it's softclip did not have sufficient ref flank

def Record_right_out_ins(amp_right_flank_list, r, repeat_index, chrom, ref, query_sequence, flank_length, coord_start, coord_end, locus_qpos_range, loci_flank_qpos_range, out_insertion_qpos_ranges_right, right_ins_rpos, flank_track, qpos_start, qpos_end):
    needed_len = amp_right_flank_list[r+repeat_index]
    if needed_len>0:
        updated_seq_end, status = detect_flank(chrom, query_sequence, ref, flank_length, coord_start, coord_end, qpos_start, qpos_end, False)
        if status:
            locus_qpos_range[r+repeat_index][1] = updated_seq_end # setting the read flank end
            flank_track[r+repeat_index][1] = True
            loci_flank_qpos_range[r+repeat_index][1] = updated_seq_end # setting the repeat end in read same as its flank end
            out_insertion_qpos_ranges_right[r+repeat_index] = () # restricting the right insertion ranges as empty
            right_ins_rpos[r+repeat_index] = () # restricting the right insertion ref pos as empty
        else:
            out_insertion_qpos_ranges_right[r+repeat_index] = (None,) # tagging the read as None to ignore this read for genotyping, as it's softclip did not have sufficient ref flank