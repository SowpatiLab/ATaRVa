import bisect
from ATARVA.operation_utils import match_jump, deletion_jump, insertion_jump

def match_parse(sub_base, int_base, rpos, repeat_index, loci_coords, tracked, locus_qpos_range, qpos, loci_flank_qpos_range, flank_track, left_flank_list, right_flank_list):
    rpos += int(sub_base)
    qpos += int(sub_base)
    rpos += int_base
    qpos += int_base 
    match_len = int_base if int_base else int(sub_base)
    repeat_index += match_jump(rpos, repeat_index, loci_coords, tracked, locus_qpos_range, qpos, match_len, loci_flank_qpos_range, flank_track, left_flank_list, right_flank_list)
    sub_base = '0'
    int_base = 0
    return [rpos, qpos, sub_base, int_base, repeat_index]

def sub_parse(base, sub_char, male, sorted_global_snp_list, global_snp_positions, read_index, global_read_variations, read_quality, rpos, repeat_index, loci_coords, tracked, locus_qpos_range, qpos, loci_flank_qpos_range, flank_track, left_flank_list, right_flank_list, hp):
    #ref_nuc = sub_char[0]
    sub_nuc = sub_char[1]

    outside_loci = True
    for each_coord in loci_coords:
        if each_coord[0] <= rpos <= each_coord[1]:
            outside_loci = False
            break

    if (not male) and outside_loci and (not hp):
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
    rpos += base
    qpos += base
    repeat_index += match_jump(rpos, repeat_index, loci_coords, tracked, locus_qpos_range, qpos, base, loci_flank_qpos_range, flank_track, left_flank_list, right_flank_list)
    base = 0
    sub_char = ''
    return [rpos, qpos, base, sub_char, repeat_index]

def ins_parse(base, insert, rpos, repeat_index, loci_keys,
                           tracked, loci_coords, homopoly_positions, read_loci_variations, locus_qpos_range, qpos, loci_flank_qpos_range, flank_track, left_flank_list, right_flank_list, out_insertion_qpos_ranges_left, out_insertion_qpos_ranges_right, left_ins_rpos, right_ins_rpos):
    qpos += base
    repeat_index += insertion_jump(base, insert, rpos, repeat_index, loci_keys,
                           tracked, loci_coords, homopoly_positions, read_loci_variations, locus_qpos_range, qpos, loci_flank_qpos_range, flank_track, left_flank_list, right_flank_list, out_insertion_qpos_ranges_left, out_insertion_qpos_ranges_right, left_ins_rpos, right_ins_rpos)
    base = 0
    return [qpos, base, repeat_index]

def del_parse(base, male, global_read_variations, read_index, rpos, repeat_index, loci_keys, tracked, loci_coords,
                          homopoly_positions, read_loci_variations, locus_qpos_range, qpos, loci_flank_qpos_range, flank_track, left_flank_list, right_flank_list):
    if not male:
        global_read_variations[read_index]['dels'].extend([rpos, rpos+base])
    rpos += base
    repeat_index += deletion_jump(base, rpos, repeat_index, loci_keys, tracked, loci_coords,
                          homopoly_positions, read_loci_variations, locus_qpos_range, qpos, loci_flank_qpos_range, flank_track, left_flank_list, right_flank_list)
    base = 0
    return [rpos, base, repeat_index]


def parse_cstag(read_index, cs_tag, read_start, loci_keys, loci_coords, read_loci_variations,
                homopoly_positions, global_read_variations, global_snp_positions, read_sequence, read_quality, cigar_one, sorted_global_snp_list, left_flank_list, right_flank_list, male, hp):
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

    base=0
    sub_base = '0'
    int_base = 0
    sub_char = ''
    symbols = {'=',':','*','-','+'}
    subs,dels,ins,iden = [False]*4
    for i in cs_tag:
        if i in symbols:
            if i == '*':
                subs = True
                if (sub_base != '0') or (int_base != 0):  #match
                    rpos, qpos, sub_base, int_base, repeat_index = match_parse(sub_base, int_base, rpos, repeat_index, loci_coords, tracked, locus_qpos_range, qpos, loci_flank_qpos_range, flank_track, left_flank_list, right_flank_list)
                    
                elif ins: # & (base != 0):
                    qpos, base, repeat_index = ins_parse(base, '', rpos, repeat_index, loci_keys,
                           tracked, loci_coords, homopoly_positions, read_loci_variations, locus_qpos_range, qpos, loci_flank_qpos_range, flank_track, left_flank_list, right_flank_list, out_insertion_qpos_ranges_left, out_insertion_qpos_ranges_right, left_ins_rpos, right_ins_rpos)
                    
                elif dels: # & (base != 0):
                    rpos, base, repeat_index = del_parse(base, male, global_read_variations, read_index, rpos, repeat_index, loci_keys, tracked, loci_coords,
                          homopoly_positions, read_loci_variations, locus_qpos_range, qpos, loci_flank_qpos_range, flank_track, left_flank_list, right_flank_list)
                    
                elif subs & (sub_char != ''):
                    rpos, qpos, base, sub_char, repeat_index = sub_parse(base, sub_char, male, sorted_global_snp_list, global_snp_positions, read_index, global_read_variations, read_quality, rpos, repeat_index, loci_coords, tracked, locus_qpos_range, qpos, loci_flank_qpos_range, flank_track, left_flank_list, right_flank_list, hp)
                    
                ins, dels, iden = [False]*3
                
            elif i == '-':
                dels = True
                if (sub_base != '0') or (int_base != 0):  #match
                    rpos, qpos, sub_base, int_base, repeat_index = match_parse(sub_base, int_base, rpos, repeat_index, loci_coords, tracked, locus_qpos_range, qpos, loci_flank_qpos_range, flank_track, left_flank_list, right_flank_list)
                    
                elif subs & (sub_char != ''):
                    rpos, qpos, base, sub_char, repeat_index = sub_parse(base, sub_char, male, sorted_global_snp_list, global_snp_positions, read_index, global_read_variations, read_quality, rpos, repeat_index, loci_coords, tracked, locus_qpos_range, qpos, loci_flank_qpos_range, flank_track, left_flank_list, right_flank_list, hp)
                    
                elif ins: # & (base != 0):
                    qpos, base, repeat_index = ins_parse(base, '', rpos, repeat_index, loci_keys,
                           tracked, loci_coords, homopoly_positions, read_loci_variations, locus_qpos_range, qpos, loci_flank_qpos_range, flank_track, left_flank_list, right_flank_list, out_insertion_qpos_ranges_left, out_insertion_qpos_ranges_right, left_ins_rpos, right_ins_rpos)
                    
                ins, subs, iden = [False]*3
                
            elif i == '+':
                ins = True
                if (sub_base != '0') or (int_base != 0):  #match
                    rpos, qpos, sub_base, int_base, repeat_index = match_parse(sub_base, int_base, rpos, repeat_index, loci_coords, tracked, locus_qpos_range, qpos, loci_flank_qpos_range, flank_track, left_flank_list, right_flank_list)
                    
                elif subs & (sub_char != ''):
                    rpos, qpos, base, sub_char, repeat_index = sub_parse(base, sub_char, male, sorted_global_snp_list, global_snp_positions, read_index, global_read_variations, read_quality, rpos, repeat_index, loci_coords, tracked, locus_qpos_range, qpos, loci_flank_qpos_range, flank_track, left_flank_list, right_flank_list, hp)
                    
                elif dels: # & (base != 0):
                    rpos, base, repeat_index = del_parse(base, male, global_read_variations, read_index, rpos, repeat_index, loci_keys, tracked, loci_coords,
                          homopoly_positions, read_loci_variations, locus_qpos_range, qpos, loci_flank_qpos_range, flank_track, left_flank_list, right_flank_list)
                    
                dels, subs, iden = [False]*3
                
            else: #elif (i == ':') or (i == '='):
                if i == '=':
                    iden = True
                if subs & (sub_char != ''):
                    rpos, qpos, base, sub_char, repeat_index = sub_parse(base, sub_char, male, sorted_global_snp_list, global_snp_positions, read_index, global_read_variations, read_quality, rpos, repeat_index, loci_coords, tracked, locus_qpos_range, qpos, loci_flank_qpos_range, flank_track, left_flank_list, right_flank_list, hp)
                    
                elif dels: # & (base != 0):
                    rpos, base, repeat_index = del_parse(base, male, global_read_variations, read_index, rpos, repeat_index, loci_keys, tracked, loci_coords,
                          homopoly_positions, read_loci_variations, locus_qpos_range, qpos, loci_flank_qpos_range, flank_track, left_flank_list, right_flank_list)
                    
                elif ins: # & (base != 0):
                    qpos, base, repeat_index = ins_parse(base, '', rpos, repeat_index, loci_keys,
                           tracked, loci_coords, homopoly_positions, read_loci_variations, locus_qpos_range, qpos, loci_flank_qpos_range, flank_track, left_flank_list, right_flank_list, out_insertion_qpos_ranges_left, out_insertion_qpos_ranges_right, left_ins_rpos, right_ins_rpos)
                    
                subs,dels,ins = [False]*3
                

        else:
            if i.isalpha():
                if iden:
                    int_base += 1
                elif subs:
                    sub_char += i
                    base = 1
                else: #elif ins or dels:
                    base += 1
            else: 
                sub_base+=i
                
    if (sub_base != '0') or (int_base != 0):  #match
        rpos, qpos, sub_base, int_base, repeat_index = match_parse(sub_base, int_base, rpos, repeat_index, loci_coords, tracked, locus_qpos_range, qpos, loci_flank_qpos_range, flank_track, left_flank_list, right_flank_list)
        
    elif ins: # & (base != 0):
        qpos, base, repeat_index = ins_parse(base, '', rpos, repeat_index, loci_keys,
               tracked, loci_coords, homopoly_positions, read_loci_variations, locus_qpos_range, qpos, loci_flank_qpos_range, flank_track, left_flank_list, right_flank_list, out_insertion_qpos_ranges_left, out_insertion_qpos_ranges_right, left_ins_rpos, right_ins_rpos)
        
    elif dels: # & (base != 0):
        rpos, base, repeat_index = del_parse(base, male, global_read_variations, read_index, rpos, repeat_index, loci_keys, tracked, loci_coords,
              homopoly_positions, read_loci_variations, locus_qpos_range, qpos, loci_flank_qpos_range, flank_track, left_flank_list, right_flank_list)
        
    elif subs & (sub_char != ''):
        rpos, qpos, base, sub_char, repeat_index = sub_parse(base, sub_char, male, sorted_global_snp_list, global_snp_positions, read_index, global_read_variations, read_quality, rpos, repeat_index, loci_coords, tracked, locus_qpos_range, qpos, loci_flank_qpos_range, flank_track, left_flank_list, right_flank_list, hp)
            
    for idx,each_key in enumerate(loci_keys):
        s_pos = locus_qpos_range[idx][0]
        e_pos = locus_qpos_range[idx][1]
        
        loci_flank_qpos_range[idx][0] = loci_flank_qpos_range[idx][0] - s_pos
        loci_flank_qpos_range[idx][1] = loci_flank_qpos_range[idx][1] - s_pos
        ins_left = [(each_tuple[0]-s_pos, each_tuple[1]-s_pos) for each_tuple in out_insertion_qpos_ranges_left[idx]]
        ins_right = [(each_tuple[0]-s_pos, each_tuple[1]-s_pos) for each_tuple in out_insertion_qpos_ranges_right[idx]]
        read_loci_variations[each_key]['seq'] = [read_sequence[s_pos:e_pos], loci_flank_qpos_range[idx], ins_left, ins_right, left_ins_rpos[idx], right_ins_rpos[idx]]