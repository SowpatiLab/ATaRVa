from ATARVA.operation_utils import match_jump, deletion_jump, insertion_jump
from ATARVA.structures import SNP

import bisect
import numpy as np

def parse_cstag(cooper, read):
    """
    Parse the CS tag for a read and record the variations observed for the read also for the loci

    :param cooper: ATaRVa object
    :param read_index: index of the read in the global read list
    :param read: pysam AlignedSegment object
    :param loci_keys: list of locus keys
    :param loci_coords: list of locus coordinates
    :param left_flank_list: list of left flank sequences for the loci
    :param right_flank_list: list of right flank sequences for the loci
    :param init_amp_var: list of initial amplicon variation information for the read
    :return: meth_start, meth_end: start and end positions of methylation in the
    """

    if cooper.cooper_sorted_snps == None: cooper.cooper_sorted_snps = []
    operations = {':', '-', '+', '*', '=', '~'}

    locus_query_range = np.zeros((len(read.loci_coords), 2), dtype=int)
    flank_query_range = np.zeros((len(read.loci_coords), 2), dtype=int)
    flankINS_query_range_left  = [[] for _ in read.loci_coords]
    flankINS_query_range_right = [[] for _ in read.loci_coords]
    left_ins_rpos = [[] for _ in read.loci_coords]
    right_ins_rpos = [[] for _ in read.loci_coords]
    repeat_tracked = [False for _ in read.loci_coords]
    flank_tracked = [[False,False] for _ in read.loci_coords]

    repeat_index = 0

    qpos = 0            # starts from 0 the sub string the read sequence in python
    if read.cigartuples[0][0] == 4:     # adjust of softclip
        qpos += read.cigartuples[0][1]

    i = 0; cs_len = len(read.cs_tag)
    while i < cs_len:

        if read.cs_tag[i] == ':':        # sequence match in short CS is followed by the length of match
            match_len = '0'; i += 1
            while i < cs_len and read.cs_tag[i] not in operations:
                match_len += read.cs_tag[i]; i += 1

            match_len = int(match_len)
            qpos += match_len; rpos += match_len
            repeat_index += match_jump(read, rpos, qpos, match_len, repeat_index, locus_query_range,
                                       flank_query_range, repeat_tracked, flank_tracked)

        elif read.cs_tag[i] == '=':      # sequence match in long CS is followed by nucs which are matching       
            match_len = 0
            while i < cs_len and read.cs_tag[i] not in operations:
                match_len += 1; i += 1

            qpos += match_len; rpos += match_len
            repeat_index += match_jump(read, rpos, qpos, match_len, repeat_index, locus_query_range,
                                       flank_query_range, repeat_tracked, flank_tracked)

        elif read.cs_tag[i] == '*':      # substitution of a base; is followed by reference and substituted base
            ref_nuc, sub_nuc = read.cs_tag[i+1], read.cs_tag[i+2]
            i += 3

            if not cooper.male:
                qual = read.query_qualities[qpos]
                cooper.cooper_read_data[read.index]['snps'].add(rpos)

                if rpos not in cooper.cooper_snp_data:
                    cooper.cooper_snp_data[rpos] = SNP(cov = 1, sub = { sub_nuc: {read.index} }, qual={ read.index: qual })
                    bisect.insort(cooper.cooper_sorted_snps, rpos)
                else:
                    cooper.cooper_snp_data[rpos].cov += 1
                    cooper.cooper_snp_data[rpos].qual[read.index] = qual
                    if sub_nuc in cooper.cooper_snp_data[rpos]: 
                        cooper.cooper_snp_data[rpos].sub[sub_nuc].add(read.index)

                    else: cooper.cooper_snp_data[rpos].sub[sub_nuc] = {read.index}

            qpos += 1; rpos += 1; match_len = 1
            repeat_index += match_jump(read, rpos, qpos, match_len, repeat_index, locus_query_range,
                                       flank_query_range, repeat_tracked, flank_tracked)

        elif read.cs_tag[i] == '+':      # insertion; is followed by the inserted bases
            insert = ''; insert_length = 0; i += 1
            while i < cs_len and read.cs_tag[i] not in operations:
                insert += read.cs_tag[i]; insert_length += 1 
                i += 1
            qpos += insert_length
            repeat_index += insertion_jump(insert_length, insert, rpos, repeat_index, loci_keys,
                                           tracked, loci_coords, homopoly_positions, read_loci_variations, locus_qpos_range, qpos, loci_flank_qpos_range, flank_track, left_flank_list, right_flank_list, out_insertion_qpos_ranges_left, out_insertion_qpos_ranges_right, left_ins_rpos, right_ins_rpos)

        elif read.cs_tag[i] == '-':      # deletion; is followed by the deleted bases
            deletion = ''; deletion_length = 0; i += 1
            while i < cs_len and read.cs_tag[i] not in operations:
                deletion += read.cs_tag[i]; deletion_length += 1
                i += 1
            if not cooper.male:
                cooper.READ_VARs[read_index].dels |= set(range(rpos, rpos+deletion_length))
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