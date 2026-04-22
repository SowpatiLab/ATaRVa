from sortedcontainers import SortedList

from ATARVA.structures import SNP

def update_snps(cooper, read, pos, qpos, insertion_point):

    rpos = read.ref_start + pos
    for ins in insertion_point:
        if ins < rpos:
            qpos += insertion_point[ins]
        elif ins>rpos: break

    Q_value = read.query_qualities[qpos]
    sub_char = read.query_sequence[qpos]
    cooper.cooper_read_data[read.index].snps.add(rpos)
    if rpos not in cooper.cooper_snp_data:
        cooper.cooper_snp_data[rpos] = SNP(cov = 1, sub = { sub_char: {read.index} }, qual={ read.index: Q_value })
        cooper.cooper_sorted_snps.add(rpos)
    else:
        cooper.cooper_snp_data[rpos].cov += 1
        cooper.cooper_snp_data[rpos].qual[read.index] = Q_value
        if sub_char in cooper.cooper_snp_data[rpos].sub: 
            cooper.cooper_snp_data[rpos].sub[sub_char].add(read.index)

        else: cooper.cooper_snp_data[rpos].sub[sub_char] = {read.index}


def parse_mdtag(cooper, read, qpos, insertion_point):
        
    if cooper.cooper_sorted_snps == None:
        cooper.cooper_sorted_snps = SortedList()

    base = 0
    sub_base = '0'
    sub_char = ''
    pos = 0

    deletion  = False
    replacing = False
    
    for i in read.md_tag:
    
        if deletion:
            if i.isalpha():
                base += 1
                continue
            else: deletion = False

        if i.isnumeric():
            sub_base += i
            if sub_char != '':
                update_snps(cooper, read, pos, qpos, insertion_point)
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
                update_snps(cooper, read, pos, qpos, insertion_point)
                pos = base - 1
                qpos+=1
                sub_char = ''
                
            sub_base = ''
            sub_char += i

    
        else: #i == '^':
            deletion = True

            if replacing:
                if sub_char != '':
                    update_snps(cooper, read, pos, qpos, insertion_point)
                    replacing = False
                    qpos+=1
                    sub_char = ''
            else:
                base += int(sub_base)
                qpos+=int(sub_base)
                sub_base = ''

                
    if sub_char != '':
        update_snps(cooper, read, pos, qpos, insertion_point)