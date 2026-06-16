import bisect
import sys
import numpy as np
from ATARVA.operation_utils import match_jump, deletion_jump, insertion_jump
from ATARVA.md_utils import parse_mdtag

def subex(ref, que):

    if '=' in que:
        array = np.frombuffer(que.encode(), dtype=np.byte)
        substitution_indices = np.where(array != ord('='))[0]
        
    else:
        array1 = np.frombuffer(ref.encode(), dtype=np.byte)
        array2 = np.frombuffer(que.encode(), dtype=np.byte)
        substitution_indices = np.where(array1 != array2)[0]
        
    return substitution_indices.tolist()

def outside_locus(loci_coords, rpos):
    """
    Check if the position is outside the locus range
    """
    for each_coord in loci_coords:
        if each_coord[0] <= rpos <= each_coord[1]:
            return False
    return True

def parse_cigar_tag(read_index, cigar_tuples, read_start, loci_keys, loci_coords, read_loci_variations,
                homopoly_positions, global_read_variations, global_snp_positions, read_sequence, read, ref, read_quality, sorted_global_snp_list, left_flank_list, right_flank_list, male, hp, init_amp_var, same_read_loci, snpQ):
    
    rpos = read_start   # NOTE: The coordinates are 1 based in SAM
    qpos = 0            # starts from 0 the sub string the read sequence in python

    # chrom = read.reference_name
    repeat_index = 0
    tracked = [False] * len(loci_coords)
    # seq_len = len(read_sequence)

    locus_qpos_range = []
    loci_flank_qpos_range = []
    out_insertion_qpos_ranges_left = []
    out_insertion_qpos_ranges_right = []
    left_ins_rpos = []
    right_ins_rpos = []
    for _ in loci_coords:
        locus_qpos_range.append([0,0])
        loci_flank_qpos_range.append([0,0])
        out_insertion_qpos_ranges_left.append([])
        out_insertion_qpos_ranges_right.append([])
        left_ins_rpos.append([])
        right_ins_rpos.append([])

    flank_track = [[False,False] for _ in loci_coords]

    X_tag = False
    insertion_point = {}

    md = 0
    if read.has_tag('MD'):
        md = 1

    if sorted_global_snp_list == None:
        sorted_global_snp_list = []

    amp_right_flank_list, amp_left_flank_list, chrom, flank, qpos_start, qpos_end = init_amp_var
    amplicon_variables = []
    if amp_left_flank_list:
        amplicon_variables = [chrom, ref, read_sequence, flank, qpos_start, qpos_end]

    for c, cigar in enumerate(cigar_tuples):
        if cigar[0] == 4:
           qpos += cigar[1] 
        elif cigar[0] == 2:     # deletion
            deletion_length = cigar[1]
            if not male:
                global_read_variations[read_index]['dels'].extend([rpos, rpos+deletion_length])
            rpos += cigar[1]
            repeat_index += deletion_jump(deletion_length, rpos, repeat_index, loci_keys, tracked, loci_coords,
                                          homopoly_positions, read_loci_variations, locus_qpos_range, qpos, loci_flank_qpos_range, flank_track, left_flank_list, right_flank_list, amp_right_flank_list, amp_left_flank_list, out_insertion_qpos_ranges_right, out_insertion_qpos_ranges_left, right_ins_rpos, left_ins_rpos, amplicon_variables)
        elif cigar[0] == 1:     # insertion
            insertion_point[rpos] = cigar[1]
            # insert = read_sequence[qpos:qpos+cigar[1]]
            insert_length = cigar[1]
            qpos += cigar[1]
            repeat_index += insertion_jump(insert_length, '', rpos, repeat_index, loci_keys,
                                           tracked, loci_coords, homopoly_positions, read_loci_variations, locus_qpos_range, qpos, loci_flank_qpos_range, flank_track, left_flank_list, right_flank_list, out_insertion_qpos_ranges_left, out_insertion_qpos_ranges_right, left_ins_rpos, right_ins_rpos, amp_right_flank_list, amp_left_flank_list, amplicon_variables)
        elif cigar[0] == 0: # match (both equals & difference)
            if (not md) & (not male) & (not hp)and (amplicon_variables==[]):
            # if (not md) & (not male) & (not hp): # considering substitutions for amplicon data as well
                ref_sequence = ref.fetch(chrom, rpos, rpos+cigar[1])
                query_sequence = read_sequence[qpos:qpos+cigar[1]]
                sub_pos = []
                if len(ref_sequence)==len(query_sequence):
                    sub_pos = subex(ref_sequence, query_sequence)
                else:
                    print('Error in fetching in sequences of ref & read for substitution')
                    sys.exit()
                for each_sub in sub_pos:
                    if not outside_locus(same_read_loci, each_sub):
                        continue
                    sub_nuc = query_sequence[each_sub]
                    Q_value = read_quality[qpos+each_sub]
                    if Q_value < snpQ:
                        continue
                    global_read_variations[read_index]['snps'].add(rpos+each_sub)
                    if rpos+each_sub not in global_snp_positions:
                        global_snp_positions[rpos+each_sub] = { 'cov': 1, sub_nuc: {read_index}, 'Qval': {read_index:Q_value} }
                        bisect.insort(sorted_global_snp_list, rpos+each_sub)
                    else:
                        global_snp_positions[rpos+each_sub]['cov'] += 1
                        global_snp_positions[rpos+each_sub]['Qval'][read_index] = Q_value
                        if sub_nuc in global_snp_positions[rpos+each_sub]: 
                            global_snp_positions[rpos+each_sub][sub_nuc].add(read_index)
                            
                        else: global_snp_positions[rpos+each_sub][sub_nuc] = {read_index}

            qpos += cigar[1]; rpos += cigar[1]; match_len = cigar[1]
            repeat_index += match_jump(rpos, repeat_index, loci_keys, loci_coords,tracked, locus_qpos_range, qpos, match_len, loci_flank_qpos_range, flank_track, left_flank_list, right_flank_list, amp_right_flank_list, amp_left_flank_list, out_insertion_qpos_ranges_right, out_insertion_qpos_ranges_left, right_ins_rpos, left_ins_rpos, amplicon_variables)

        elif cigar[0] == 7: # exact match (equals)
            X_tag = True
            qpos += cigar[1]; rpos += cigar[1]; match_len = cigar[1]
            repeat_index += match_jump(rpos, repeat_index, loci_keys, loci_coords,tracked, locus_qpos_range, qpos, match_len, loci_flank_qpos_range, flank_track, left_flank_list, right_flank_list, amp_right_flank_list, amp_left_flank_list, out_insertion_qpos_ranges_right, out_insertion_qpos_ranges_left, right_ins_rpos, left_ins_rpos, amplicon_variables)

        elif cigar[0] == 8: # substitution (difference)
            X_tag = True
            if (not male) and outside_locus(same_read_loci, rpos) and (not hp) and (amplicon_variables==[]):
            # if (not male) and outside_locus(same_read_loci, rpos) and (not hp): # considering substitutions for amplicon data as well
                sub_nuc = read_sequence[qpos]
                Q_value = read_quality[qpos]
                if Q_value >= snpQ:
                    global_read_variations[read_index]['snps'].add(rpos)
                    if rpos not in global_snp_positions:
                        global_snp_positions[rpos] = { 'cov': 1, sub_nuc: {read_index}, 'Qval': {read_index:Q_value} }
                        bisect.insort(sorted_global_snp_list, rpos)
                    else:
                        global_snp_positions[rpos]['cov'] += 1
                        global_snp_positions[rpos]['Qval'][read_index] = Q_value
                        if sub_nuc in global_snp_positions[rpos]: 
                            global_snp_positions[rpos][sub_nuc].add(read_index)
                            
                        else: global_snp_positions[rpos][sub_nuc] = {read_index}
            qpos += cigar[1]; rpos += cigar[1]; match_len = cigar[1]
            repeat_index += match_jump(rpos, repeat_index, loci_keys, loci_coords,tracked, locus_qpos_range, qpos, match_len, loci_flank_qpos_range, flank_track, left_flank_list, right_flank_list, amp_right_flank_list, amp_left_flank_list, out_insertion_qpos_ranges_right, out_insertion_qpos_ranges_left, right_ins_rpos, left_ins_rpos, amplicon_variables)

    if not X_tag :
        if read.has_tag('MD'):
            if cigar_tuples[0][0] == 4: qpos = cigar_tuples[0][1]
            else: qpos=0
            MD_tag = read.get_tag('MD')
            parse_mdtag(MD_tag, qpos, read_start, global_read_variations, global_snp_positions, read_index, read_quality, read_sequence, sorted_global_snp_list, insertion_point, same_read_loci, male, hp, amplicon_variables, snpQ)

    tot_loc = len(loci_keys) - 1 # 0-based counting
    meth_start = None; meth_end = None # repeat sequence coordinates in the read for meth-calc
    for idx,each_key in enumerate(loci_keys):

        if amp_left_flank_list:
            soft_flank_status = [] # 0 means presence, 1 means absence of sufficient soft-flank
            if type(out_insertion_qpos_ranges_left[idx]) == tuple:
                soft_flank_status.append(len(out_insertion_qpos_ranges_left[idx])) # can be empty or None value
            else:
                soft_flank_status.append(0) # if its a list (without reaching the soft-clip; sufficient flank already aligned by aligner)
            if type(out_insertion_qpos_ranges_right[idx]) == tuple:
                soft_flank_status.append(len(out_insertion_qpos_ranges_right[idx])) # can be empty or None value
            else:
                soft_flank_status.append(0) # if its a list (without reaching the soft-clip; sufficient flank already aligned by aligner)

            if any(soft_flank_status): # if any of the soft-flanks are not sufficient then dont write
                pass #continue
            else:
                pass
            
        s_pos = locus_qpos_range[idx][0]
        if idx==0:
            meth_start = s_pos
        e_pos = locus_qpos_range[idx][1]
        if idx==tot_loc:
            meth_end = e_pos

        loci_flank_qpos_range[idx][0] = loci_flank_qpos_range[idx][0] - s_pos
        loci_flank_qpos_range[idx][1] = loci_flank_qpos_range[idx][1] - s_pos

        tmp_left_set = out_insertion_qpos_ranges_left[idx]
        ins_left = [(each_tuple[0]-s_pos, each_tuple[1]-s_pos) for each_tuple in tmp_left_set if len(tmp_left_set)==2 and None not in tmp_left_set]
        tmp_right_set = out_insertion_qpos_ranges_right[idx]
        ins_right = [(each_tuple[0]-s_pos, each_tuple[1]-s_pos) for each_tuple in tmp_right_set if len(tmp_right_set)==2 and None not in tmp_right_set]

        read_loci_variations[each_key]['seq'] = [read_sequence[s_pos:e_pos], loci_flank_qpos_range[idx], ins_left, ins_right, left_ins_rpos[idx], right_ins_rpos[idx], s_pos, e_pos]

    return meth_start, meth_end

    