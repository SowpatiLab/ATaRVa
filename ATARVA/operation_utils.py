def update_homopolymer_coords(ref_seq, locus_start, homopoly_positions):
    """
    Record all the homopolymer stretches of at least 3 bases within the repeat coordinates
    """
    prev_N = ref_seq[0]; start = -1
    for i, n in enumerate(ref_seq[1:]):
        if n == prev_N:
            if start == -1: start = i
        else:
            if start != -1 and (i+1)-start >= 4:
                for l,c in enumerate(range(locus_start+start, locus_start+i+1)):
                    homopoly_positions[c] = (i-start+1)-l
            start = -1
        prev_N = n

    if start != -1 and (i+1)-start >= 4:
        for l,c in enumerate(range(locus_start+start, locus_start+i+1)):
            # for each position in the homopolymer stretch we record the length of the 
            # homopolymer nucleotides on the right
            homopoly_positions[c] = (i-start+1)-l


def match_jump(rpos, repeat_index, loci_coords, tracked, locus_qpos_range, qpos, match_len, loci_flank_qpos_range, flank_track, left_flank, right_flank):
    """
    Return the number of repeat indices to jump when scanning through a match segment
    """
    previous_rpos = rpos - match_len
    r = 0 
    for r,coord in enumerate(loci_coords[repeat_index:]):
        coord_start, coord_end = coord
        
        if rpos < coord_start: break
        
        if previous_rpos > coord_end: continue
            
        if not tracked[r+repeat_index]:

            if coord_start <= rpos:
                
                locus_qpos_range[r+repeat_index][0] = qpos - (rpos - coord_start)
            if coord_end < rpos:
                
                locus_qpos_range[r+repeat_index][1] = qpos - (rpos - coord_end)

            tracked[r+repeat_index] = True 

            # for storing repeat qpos ranges
            if coord_start+left_flank[r+repeat_index] <= rpos:
                loci_flank_qpos_range[r+repeat_index][0] = qpos - (rpos - coord_start)+left_flank[r+repeat_index]
                flank_track[r+repeat_index][0] = True
            if coord_end-right_flank[r+repeat_index] <= rpos:
                loci_flank_qpos_range[r+repeat_index][1] = qpos - (rpos - coord_end)-right_flank[r+repeat_index]
                if rpos > coord_end-right_flank[r+repeat_index]: flank_track[r+repeat_index][1] = True
                
                

        elif coord_end <= rpos:
            
            locus_qpos_range[r+repeat_index][1] = qpos - (rpos -coord_end)

        # for storing repeat qpos ranges
        if not flank_track[r+repeat_index][0]:
            if coord_start+left_flank[r+repeat_index] <= rpos:
                loci_flank_qpos_range[r+repeat_index][0] = qpos - (rpos -coord_start)+left_flank[r+repeat_index]
                flank_track[r+repeat_index][0] = True 
            if coord_end-right_flank[r+repeat_index] <= rpos:
                loci_flank_qpos_range[r+repeat_index][1] = qpos - (rpos - coord_end)-right_flank[r+repeat_index]
                if rpos > coord_end-right_flank[r+repeat_index]: flank_track[r+repeat_index][1] = True
            
        elif (not flank_track[r+repeat_index][1]) and (coord_end-right_flank[r+repeat_index] <= rpos):
            loci_flank_qpos_range[r+repeat_index][1] = qpos - (rpos -coord_end)-right_flank[r+repeat_index]
            if rpos > coord_end-right_flank[r+repeat_index]: flank_track[r+repeat_index][1] = True

    
    jump = 0    # jump beyond the repeat where all positions are tracked
    if loci_coords[repeat_index + r - 1][1] < rpos:
        for f in loci_coords[repeat_index:]:
            if f[1] < rpos: jump += 1
            else: break
    
    return jump


def deletion_jump(deletion_length, rpos, repeat_index, loci_keys, tracked, loci_coords, homopoly_positions, read_loci_variations, locus_qpos_range, qpos, loci_flank_qpos_range, flank_track, left_flank, right_flank):
    """
    Return the number of repeat indices to jump when scanning through a deletion segment.
    The function tracks specifically if the deletion is segment has complete repeats in them
    or segments of the repeat is deleted.
    """

    # rpos - corresponds to the position in the reference after tracking the deletion
    r = 0   # required to be initialised outside the loop
    for r, coord in enumerate(loci_coords[repeat_index:]):
        coord_start, coord_end = coord
        # if rpos is before the start of the repeat; repeat is unaffected
        if rpos < coord_start: break

        # actual position in the reference where the deletion is occurring
        del_pos = rpos - deletion_length
        if del_pos > coord_end: continue

        locus_key = loci_keys[r+repeat_index]
        if not tracked[r+repeat_index]:
            # if the locus is not tracked
            # deletion is encountered beyond
            if coord_start <= rpos:    
                locus_qpos_range[r+repeat_index][0] = qpos        
                tracked[r+repeat_index] = True    # set tracked as true

            if coord_end < rpos:
                
                locus_qpos_range[r+repeat_index][1] = qpos

            # for storing repeat qpos ranges
            if coord_start+left_flank[r+repeat_index] <= rpos:
                loci_flank_qpos_range[r+repeat_index][0] = qpos
                flank_track[r+repeat_index][0] = True
            if coord_end-right_flank[r+repeat_index] < rpos:
                loci_flank_qpos_range[r+repeat_index][1] = qpos
                flank_track[r+repeat_index][1] = True

        elif coord_end < rpos:
            
            locus_qpos_range[r+repeat_index][1] = qpos

        # for storing repeat qpos ranges
        if not flank_track[r+repeat_index][0]:
            if coord_start+left_flank[r+repeat_index] <= rpos:
                loci_flank_qpos_range[r+repeat_index][0] = qpos
                flank_track[r+repeat_index][0] = True 
            if coord_end-right_flank[r+repeat_index] < rpos:
                loci_flank_qpos_range[r+repeat_index][1] = qpos
                flank_track[r+repeat_index][1] = True
        elif (not flank_track[r+repeat_index][1]) and (coord_end-right_flank[r+repeat_index] <= rpos):
            loci_flank_qpos_range[r+repeat_index][1] = qpos
            if rpos > coord_end-right_flank[r+repeat_index]: flank_track[r+repeat_index][1] = True

        # updating the allele with the deletion considered
        # read_loci_variations[locus_key][rpos] = f'D|{deletion_length}'
        
        # del_len = min(coord[1], rpos) - max(coord[0], del_pos)
        del_len = min(coord_end-right_flank[r+repeat_index], rpos) - max(coord_start+left_flank[r+repeat_index], del_pos)
        if (rpos >= coord_start+left_flank[r+repeat_index]) and (del_pos <= coord_end-right_flank[r+repeat_index]): # introduced to include length only if it comes inside repeat region
            if del_pos not in homopoly_positions:
                read_loci_variations[locus_key]['alen'] -= del_len
                read_loci_variations[locus_key]['halen'] -= del_len
            else:
                if del_len <= homopoly_positions[del_pos]:
                    # if the deletion is only limited to the homopolymer positions
                    read_loci_variations[locus_key]['halen'] -= del_len
                else:
                    read_loci_variations[locus_key]['alen'] -= del_len
                    read_loci_variations[locus_key]['halen'] -= del_len


    jump = 0    # jump beyond the repeat where all positions are tracked
    if loci_coords[repeat_index + r - 1][1] < rpos:
        for f in loci_coords[repeat_index:]:
            if f[1] < rpos: jump += 1
            else: break

    return jump


def insertion_jump(insertion_length, insert, rpos, repeat_index, loci_keys, tracked, loci_coords, homopoly_positions, read_loci_variations, locus_qpos_range, qpos, loci_flank_qpos_range, flank_track, left_flank, right_flank, out_insertion_qpos_ranges_left, out_insertion_qpos_ranges_right, left_ins_rpos, right_ins_rpos):
    """
    Return the number of repeat indices to jump when scanning through a insertion segment.
    The function tracks specifically if the deletion is segment has complete repeats in them
    or segments of the repeat is deleted.
    """
    r = 0   # required to be initialised outside the loop
    for r, coord in enumerate(loci_coords[repeat_index:]):
        coord_start, coord_end = coord
        # if rpos is before the start of the repeat; repeat is unaffected
        if rpos < coord_start: break

        # if the insertion is happening beyond, the repeat in unaffected
        if rpos > coord_end: continue

        locus_key = loci_keys[r+repeat_index]
        if not tracked[r+repeat_index]:
            # if the locus is not tracked
            # deletion is encountered beyond
            if coord_start <= rpos:
                locus_qpos_range[r+repeat_index][0] = qpos-insertion_length
                tracked[r+repeat_index] = True    # set tracked as true
            if coord_end == rpos:
               
                locus_qpos_range[r+repeat_index][1] = qpos
                # here jump can be done

            # for storing repeat qpos ranges
            if coord_start+left_flank[r+repeat_index]-1 <= rpos:
                loci_flank_qpos_range[r+repeat_index][0] = qpos-insertion_length
                flank_track[r+repeat_index][0] = True
            if coord_end-right_flank[r+repeat_index] <= rpos:
                loci_flank_qpos_range[r+repeat_index][1] = qpos
                if rpos > coord_end-right_flank[r+repeat_index]: flank_track[r+repeat_index][1] = True


        elif coord_end == rpos:
            
            locus_qpos_range[r+repeat_index][1] = qpos

        # for storing repeat qpos ranges
        if not flank_track[r+repeat_index][0]:
            if coord_start+left_flank[r+repeat_index] <= rpos:
                loci_flank_qpos_range[r+repeat_index][0] = qpos-insertion_length
                flank_track[r+repeat_index][0] = True 
            if coord_end-right_flank[r+repeat_index] <= rpos:
                loci_flank_qpos_range[r+repeat_index][1] = qpos
                if rpos > coord_end-right_flank[r+repeat_index]: flank_track[r+repeat_index][1] = True
        elif (not flank_track[r+repeat_index][1]) and (coord_end-right_flank[r+repeat_index] <= rpos):
            loci_flank_qpos_range[r+repeat_index][1] = qpos
            if rpos > coord_end-right_flank[r+repeat_index]: flank_track[r+repeat_index][1] = True

        # read_loci_variations[locus_key][rpos] = f'I|{insertion_length}'
        if coord_start+left_flank[r+repeat_index] <= rpos <= coord_end-right_flank[r+repeat_index]: # introduced to include length only if it comes inside repeat region
            if rpos not in homopoly_positions:
                read_loci_variations[locus_key]['alen'] += insertion_length
                read_loci_variations[locus_key]['halen'] += insertion_length
            else:
                if len(set(insert)) == 1:
                    # only if the insertion is a homopolymer; consider it as homopolymer insertion
                    read_loci_variations[locus_key]['halen'] += insertion_length
                else:
                    read_loci_variations[locus_key]['alen'] += insertion_length
                    read_loci_variations[locus_key]['halen'] += insertion_length

        if coord_start <= rpos <= coord_start+left_flank[r+repeat_index]-1: # -1 is included so ins near the start pos is not taken into account as it is already added
            out_insertion_qpos_ranges_left[r+repeat_index].append((qpos-insertion_length, qpos))
            left_ins_rpos[r+repeat_index].append(rpos)
        elif coord_end-right_flank[r+repeat_index]+1 <= rpos <= coord_end: # +1 is included so ins near the end pos is not taken into account as it is already added
            out_insertion_qpos_ranges_right[r+repeat_index].append((qpos-insertion_length, qpos))
            right_ins_rpos[r+repeat_index].append(rpos)

    jump = 0    # jump beyond the repeat where all positions are tracked
    if loci_coords[repeat_index + r - 1][1] < rpos:
        for f in loci_coords[repeat_index:]:
            if f[1] < rpos: jump += 1
            else: break

    return jump