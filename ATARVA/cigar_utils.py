from ATARVA.operation_utils import match_jump, deletion_jump, insertion_jump
from ATARVA.md_utils import parse_mdtag
from ATARVA.structures import SNP

import bisect
import numpy as np

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
    Check if the position is outside all locus ranges (sorted coords)
    
    :param loci_coords: list of locus coordinates
    :param rpos: position to be checked
    :return: True if the position is outside all locus ranges, False otherwise
    """

    starts = [c[0] for c in loci_coords]
    idx    = bisect.bisect_right(starts, rpos) - 1
    return idx < 0 or rpos > loci_coords[idx][1]


def parse_cigar(cooper, read):
    """
    parse the cigar to identify loci specific variations and the SNPs

    :param cooper: ATaRVa object
    :param read: pysam AlignedSegment object
    """

    rpos = read.reference_start   # NOTE: The coordinates are 1 based in SAM
    qpos = 0            # starts from 0 the sub string the read sequence in python

    chrom = read.chrom
    repeat_index = 0

    locus_query_range = np.zeros((len(read.loci_coords), 2), dtype=int)
    flank_query_range = np.zeros((len(read.loci_coords), 2), dtype=int)
    flankINS_query_range_left  = [[] for _ in read.loci_coords]
    flankINS_query_range_right = [[] for _ in read.loci_coords]
    left_ins_rpos = [[] for _ in read.loci_coords]
    right_ins_rpos = [[] for _ in read.loci_coords]
    repeat_tracked = [False for _ in read.loci_coords]
    flank_tracked = [[False,False] for _ in read.loci_coords]

    has_subop = False
    has_MD = read.has_tag('MD')

    insert_positions = {}

    for c, cigar in enumerate(read.cigartuples):
    
        if cigar[0] == 4: qpos += cigar[1]    # softclip; adjust the qpos but not rpos; also consider the possibility of match after softclip
 
        elif cigar[0] == 2:     # deletion
            deletion_len = cigar[1]
            if not cooper.male:        ## ?? Why are we doing this only for female samples
                cooper.cooper_read_data[read.index].dels.extend([rpos, rpos + deletion_len])
            rpos += deletion_len
            repeat_index += deletion_jump(cooper, read, rpos, qpos, deletion_len, repeat_index, locus_query_range, flank_query_range, repeat_tracked, flank_tracked)

        elif cigar[0] == 1:     # insertion
            insert_positions[rpos] = cigar[1]
            insert_len = cigar[1]

            qpos += insert_len
            repeat_index += insertion_jump(cooper, read, rpos, qpos, insert_len, repeat_index, locus_query_range, flank_query_range, repeat_tracked, flank_tracked,
                                           flankINS_query_range_left, flankINS_query_range_right, left_ins_rpos, right_ins_rpos)
        
        elif cigar[0] == 0: # match (both equals & difference)
            match_len = cigar[1]
            if not has_MD and (not cooper.male) and (not cooper.args.haplotag):
                ref_sequence   = cooper.ref.fetch(chrom, rpos, rpos + match_len)
                query_sequence = read.query_sequence[qpos:qpos + match_len]

                assert len(ref_sequence) == len(query_sequence), \
                    "Error: fetching sequences based on CIGAR. Please check cigar formats"

                substitutions = subex(ref_sequence, query_sequence)
                for sub_pos in substitutions:

                    if not outside_locus(read.loci_coords, rpos + sub_pos):
                        continue
                    sub_nuc  = query_sequence[sub_pos]
                    sub_qual = read.query_qualities[qpos + sub_pos]
                    if sub_qual < cooper.args.snp_qual: continue
                    cooper.cooper_read_data[read.index].snps.add(rpos + sub_pos)
                    if rpos + sub_pos not in cooper.cooper_snp_data:
                        cooper.cooper_snp_data[rpos + sub_pos] = SNP(cov = 1, sub = { sub_nuc: {read.index} }, qual={ read.index: sub_qual })
                        cooper.cooper_sorted_snps.add(rpos + sub_pos)
                    else:
                        cooper.cooper_snp_data[rpos + sub_pos].cov += 1
                        cooper.cooper_snp_data[rpos + sub_pos].qual[read.index] = sub_qual
                        if sub_nuc in cooper.cooper_snp_data[rpos + sub_pos].sub:
                            cooper.cooper_snp_data[rpos + sub_pos].sub[sub_nuc].add(read.index)
                        else:
                            cooper.cooper_snp_data[rpos + sub_pos].sub[sub_nuc] = {read.index}

            qpos += match_len; rpos += match_len
            repeat_index += match_jump(read, rpos, qpos, match_len, repeat_index, locus_query_range, flank_query_range, repeat_tracked, flank_tracked)

        elif cigar[0] == 7: # exact match (equals)
            has_subop = True
            match_len = cigar[1]
            qpos += match_len; rpos += match_len
            repeat_index += match_jump(read, rpos, qpos, match_len, repeat_index, locus_query_range, flank_query_range, repeat_tracked, flank_tracked)

        elif cigar[0] == 8: # substitution (difference)
            has_subop = True
            match_len = cigar[1]
            if (not cooper.male) and outside_locus(read.loci_coords, rpos) and (not cooper.args.haplotag):
                sub_nuc  = read.query_sequence[qpos]
                sub_qual = read.query_qualities[qpos]
                if sub_qual >= cooper.args.snp_qual:
                    cooper.cooper_read_data[read.index].snps.add(rpos)
                    if rpos not in cooper.cooper_snp_data:
                        cooper.cooper_snp_data[rpos] = SNP(cov = 1, sub = { sub_nuc: {read.index} }, qual={ read.index: sub_qual })
                        cooper.cooper_sorted_snps.add(rpos)
                    else:
                        cooper.cooper_snp_data[rpos].cov += 1
                        cooper.cooper_snp_data[rpos].qual[read.index] = sub_qual
                        if sub_nuc in cooper.cooper_snp_data[rpos].sub: 
                            cooper.cooper_snp_data[rpos].sub[sub_nuc].add(read.index)
                        else:
                            cooper.cooper_snp_data[rpos].sub[sub_nuc] = {read.index}

            qpos += match_len; rpos += match_len
            repeat_index += match_jump(read, rpos, qpos, match_len, repeat_index, locus_query_range, flank_query_range, repeat_tracked, flank_tracked)

    if not has_subop :
        if read.has_tag('MD'):
            if read.cigartuples[0][0] == 4: qpos = read.cigartuples[0][1]
            else: qpos = 0
            parse_mdtag(cooper, read, qpos, insert_positions)

    num_read_loci = len(read.loci_coords)
    for idx, locus_key in enumerate(read.loci_keys):
        locus_start = locus_query_range[idx][0]
        if idx == 0: read.methyl_start = locus_start
        
        locus_end = locus_query_range[idx][1]
        if idx == num_read_loci - 1: read.methyl_end = locus_end

        flank_query_range[idx][0] = flank_query_range[idx][0] - locus_start
        flank_query_range[idx][1] = flank_query_range[idx][1] - locus_start
        
        ins_left  = [(coords[0] - locus_start, coords[1] - locus_start) for coords in flankINS_query_range_left[idx] ]
        ins_right = [(coords[0] - locus_start, coords[1] - locus_start) for coords in flankINS_query_range_right[idx]]
        read.loci_data[locus_key].seq = [read.query_sequence[locus_start:locus_end], flank_query_range[idx],
                                         ins_left, ins_right, left_ins_rpos[idx], right_ins_rpos[idx],
                                         locus_start, locus_end]
    