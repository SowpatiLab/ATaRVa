import bisect
from operation_utils import match_jump, deletion_jump, insertion_jump


def parse_cstag(read_index, cs_tag, read_start, loci_keys, loci_coords, read_loci_variations,
                homopoly_positions, global_read_variations, global_snp_positions, read_sequence, read_quality, cigar_one, sorted_global_snp_list, left_flank_list, right_flank_list, male):
    """
    Parse the CS tag for a read and record the variations observed for the read also for the loci
    """
    if sorted_global_snp_list == None:
        sorted_global_snp_list = []
    operations = {':', '-', '+', '*', '=', '~'}
    rpos = read_start   # NOTE: The coordinates are 1 based in SAM
    qpos = 0            # starts from 0 the sub string the read sequence in python

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

    repeat_index = 0
    tracked = [False] * len(loci_coords)        

    flank_track = [[False,False] for _ in loci_coords]

    if cigar_one[0] == 4:
        qpos+=cigar_one[1]

    i = 0; cs_len = len(cs_tag)
    while i<cs_len:

        if cs_tag[i] == ':':        # sequence match in short CS is followed by the length of match
            match_len = '0'; i += 1
            while i < cs_len and cs_tag[i] not in operations:
                match_len += cs_tag[i]; i += 1

            match_len = int(match_len)
            qpos += match_len; rpos += match_len
            repeat_index += match_jump(rpos, repeat_index, loci_coords, tracked, locus_qpos_range, qpos, match_len, loci_flank_qpos_range, flank_track, left_flank_list, right_flank_list)

        elif cs_tag[i] == '=':      # sequence match in long CS is followed by nucs which are matching       
            match_len = 0
            while i < cs_len and cs_tag[i] not in operations:
                match_len += 1; i += 1

            qpos += match_len; rpos += match_len
            repeat_index += match_jump(rpos, repeat_index, loci_coords, tracked, locus_qpos_range, qpos, match_len, loci_flank_qpos_range, flank_track, left_flank_list, right_flank_list)

        elif cs_tag[i] == '*':      # substitution of a base; is followed by reference and substituted base
            ref_nuc = cs_tag[i+1]; sub_nuc = cs_tag[i+2]
            i += 3

            if not male:
                Q_value = read_quality[qpos]
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
            
            qpos += 1; rpos += 1; match_len = 1
            repeat_index += match_jump(rpos, repeat_index, loci_coords, tracked, locus_qpos_range, qpos, match_len, loci_flank_qpos_range, flank_track, left_flank_list, right_flank_list)

        elif cs_tag[i] == '+':      # insertion; is followed by the inserted bases
            insert = ''; insert_length = 0; i += 1
            while i < cs_len and cs_tag[i] not in operations:
                insert += cs_tag[i]; insert_length += 1 
                i += 1
            qpos += insert_length
            repeat_index += insertion_jump(insert_length, insert, rpos, repeat_index, loci_keys,
                                           tracked, loci_coords, homopoly_positions, read_loci_variations, locus_qpos_range, qpos, loci_flank_qpos_range, flank_track, left_flank_list, right_flank_list, out_insertion_qpos_ranges_left, out_insertion_qpos_ranges_right, left_ins_rpos, right_ins_rpos)

        elif cs_tag[i] == '-':      # deletion; is followed by the deleted bases
            deletion = ''; deletion_length = 0; i += 1
            while i < cs_len and cs_tag[i] not in operations:
                deletion += cs_tag[i]; deletion_length += 1
                i += 1
            if not male:
                global_read_variations[read_index]['dels'] |= set(range(rpos, rpos+deletion_length))
            rpos += deletion_length
            repeat_index += deletion_jump(deletion_length, rpos, repeat_index, loci_keys, tracked, loci_coords,
                                          homopoly_positions, read_loci_variations, locus_qpos_range, qpos, loci_flank_qpos_range, flank_track, left_flank_list, right_flank_list)
            
    for idx,each_key in enumerate(loci_keys):
        s_pos = locus_qpos_range[idx][0]
        e_pos = locus_qpos_range[idx][1]
        
        loci_flank_qpos_range[idx][0] = loci_flank_qpos_range[idx][0] - s_pos
        loci_flank_qpos_range[idx][1] = loci_flank_qpos_range[idx][1] - s_pos
        ins_left = [(each_tuple[0]-s_pos, each_tuple[1]-s_pos) for each_tuple in out_insertion_qpos_ranges_left[idx]]
        ins_right = [(each_tuple[0]-s_pos, each_tuple[1]-s_pos) for each_tuple in out_insertion_qpos_ranges_right[idx]]
        read_loci_variations[each_key]['seq'] = [read_sequence[s_pos:e_pos], loci_flank_qpos_range[idx], ins_left, ins_right, left_ins_rpos[idx], right_ins_rpos[idx]]