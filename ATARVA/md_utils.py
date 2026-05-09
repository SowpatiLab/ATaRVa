from sortedcontainers import SortedList

from ATARVA.structures import SNP


def update_snps(cooper_snp_data, cooper_sorted_snps, cooper_loci_keys, cooper_loci_info, cooper_read_data, cooper_prev_reads,
                snp_qual_cutoff, is_haploid, is_haplotag, read, pos, qpos, insertion_point, from_md=True):
    """
    update SNP information in cooper.cooper_snp_data and also for this read in read.snps

    :param cooper:             cooper object
    :param read:               read object
    :param pos:                position of SNP on reference
    :param qpos:               position of SNP on read
    :param insertion_point:    dictionary of insertion positions and lengths for this read
    :return:                   None; updates cooper.cooper_snp_data and read.snps in place
    """

    rpos = 0
    if from_md:
        rpos = read.ref_start + pos

        for ins in insertion_point:
            if ins < rpos:
                qpos += insertion_point[ins]
            # elif ins == rpos:
            #     return  # if the SNP is next to an insertion ignore it
            elif ins > rpos: break
    else: rpos = pos

    if rpos in cooper_read_data[read.index].no_snps:
        return
    for locus in cooper_loci_keys:
        if cooper_loci_info[locus].start <= rpos <= cooper_loci_info[locus].end:
            return
        if cooper_loci_info[locus].start > rpos: break

    sub_qval = read.query_qualities[qpos]
    
    if sub_qval >= snp_qual_cutoff and (not is_haploid) and (not is_haplotag):
        sub_char = read.query_sequence[qpos]
        cooper_read_data[read.index].snps.add(rpos)

        for read_index in cooper_prev_reads.copy():
            _read = cooper_read_data[read_index]
            if _read.start <= rpos <= _read.end and rpos not in _read.snps:
                cooper_prev_reads.remove(read_index)

        if rpos not in cooper_snp_data:
            cooper_snp_data[rpos] = SNP(cov = 1, sub = { sub_char: {read.index} }, qual={ read.index: sub_qval })
            cooper_sorted_snps.add(rpos)
        else:
            cooper_snp_data[rpos].cov += 1
            cooper_snp_data[rpos].qual[read.index] = sub_qval
            if sub_char in cooper_snp_data[rpos].sub: 
                cooper_snp_data[rpos].sub[sub_char].add(read.index)

            else: cooper_snp_data[rpos].sub[sub_char] = {read.index}


def parse_mdtag(cooper_snp_data, cooper_sorted_snps, cooper_loci_keys, cooper_loci_info, cooper_read_data, cooper_prev_reads,
                snp_qual_cutoff, is_haploid, is_haplotag, read, qpos, insertion_point):
    """
    parse the MD tag of a read to identify SNP positions and update cooper.cooper_snp_data and read.snps accordingly
    :param cooper:             cooper object
    :param read:               read object
    :param qpos:               current position on read sequence (0-based)
    :param insertion_point:    dictionary of insertion positions and lengths for this read
    :return:                   None; updates cooper.cooper_snp_data and read.snps in place
    """
        
    if cooper_sorted_snps == None:
        cooper_sorted_snps = SortedList()

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
                update_snps(cooper_snp_data, cooper_sorted_snps, cooper_loci_keys, cooper_loci_info, cooper_read_data, cooper_prev_reads,
                            snp_qual_cutoff, is_haploid, is_haplotag, read, pos, qpos, insertion_point, True)
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
                update_snps(cooper_snp_data, cooper_sorted_snps, cooper_loci_keys, cooper_loci_info, cooper_read_data, cooper_prev_reads,
                            snp_qual_cutoff, is_haploid, is_haplotag, read, pos, qpos, insertion_point, True)
                pos = base - 1
                qpos+=1
                sub_char = ''
                
            sub_base = ''
            sub_char += i

    
        else: #i == '^':
            deletion = True

            if replacing:
                if sub_char != '':
                    update_snps(cooper_snp_data, cooper_sorted_snps, cooper_loci_keys, cooper_loci_info, cooper_read_data, cooper_prev_reads,
                                snp_qual_cutoff, is_haploid, is_haplotag, read, pos, qpos, insertion_point, True)
                    replacing = False
                    qpos+=1
                    sub_char = ''
            else:
                base += int(sub_base)
                qpos+=int(sub_base)
                sub_base = ''

                
    if sub_char != '':
        update_snps(cooper_snp_data, cooper_sorted_snps, cooper_loci_keys, cooper_loci_info, cooper_read_data, cooper_prev_reads,
                    snp_qual_cutoff, is_haploid, is_haplotag, read, pos, qpos, insertion_point, True)
    