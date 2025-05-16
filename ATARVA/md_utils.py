import bisect

def update_global_snpPos(ref_start, pos, global_read_variations, global_snp_positions, read_index, read_sequence, read_quality, sorted_global_snp_list, insertion_point, qpos, loci_coords, male, hp):

    rpos = ref_start+pos
    for ins in insertion_point:
        if ins<rpos:
            qpos+=insertion_point[ins]
        elif ins>rpos: break

    outside_loci = True
    for each_coord in loci_coords:
        if each_coord[0] <= rpos <= each_coord[1]:
            outside_loci = False
            break

    if (not male) and outside_loci and (not hp):
        Q_value = read_quality[qpos]
        sub_char = read_sequence[qpos]
        global_read_variations[read_index]['snps'].add(rpos)
        if rpos not in global_snp_positions:
            global_snp_positions[rpos] = { 'cov': 1, sub_char: {read_index}, 'Qval': {read_index:Q_value} }
            bisect.insort(sorted_global_snp_list, rpos)
        else:
            global_snp_positions[rpos]['cov'] += 1
            global_snp_positions[rpos]['Qval'][read_index] = Q_value
            if sub_char in global_snp_positions[rpos]: 
                global_snp_positions[rpos][sub_char].add(read_index)

            else: global_snp_positions[rpos][sub_char] = {read_index}



def parse_mdtag(MD_tag, qpos, ref_start, global_read_variations, global_snp_positions, read_index, read_quality, read_sequence, sorted_global_snp_list, insertion_point, loci_coords, male, hp):
        
    if sorted_global_snp_list == None:
        sorted_global_snp_list = []

    base = 0
    sub_base='0'
    sub_char=''
    pos=0

    deletion = False
    replacing = False
    
    for i in MD_tag:
    
        if deletion:
            if i.isalpha():
                base+=1
                continue
            else: deletion = False
                
        
        if i.isnumeric():
            sub_base+=i
            if sub_char != '':
                update_global_snpPos(ref_start, pos, global_read_variations, global_snp_positions, read_index, read_sequence, read_quality, sorted_global_snp_list, insertion_point, qpos, loci_coords, male, hp)
                replacing = False
                qpos+=1
                sub_char = ''
                
        elif i.isalpha():
            replacing = True

            if sub_char == '':
                base += int(sub_base)+1
                pos = base - 1
                qpos+=int(sub_base)
            else:
                base+=1
                update_global_snpPos(ref_start, pos, global_read_variations, global_snp_positions, read_index, read_sequence, read_quality, sorted_global_snp_list, insertion_point, qpos, loci_coords, male, hp)
                pos = base - 1
                qpos+=1
                sub_char = ''
                
            sub_base = ''
            sub_char += i

    
        else: #i == '^':
            deletion = True

            if replacing:
                if sub_char != '':
                    update_global_snpPos(ref_start, pos, global_read_variations, global_snp_positions, read_index, read_sequence, read_quality, sorted_global_snp_list, insertion_point, qpos, loci_coords, male, hp)
                    replacing = False
                    qpos+=1
                    sub_char = ''
            else:
                base += int(sub_base)
                qpos+=int(sub_base)
                sub_base = ''

                
    if sub_char != '':
        update_global_snpPos(ref_start, pos, global_read_variations, global_snp_positions, read_index, read_sequence, read_quality, sorted_global_snp_list, insertion_point, qpos, loci_coords, male, hp)