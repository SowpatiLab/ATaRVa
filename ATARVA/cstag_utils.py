from ATARVA.md_utils import update_snps
from ATARVA.operation_utils import match_jump, deletion_jump, insertion_jump
from ATARVA.md_utils import update_snps

import numpy as np

def parse_cstag(cooper, read):
    """
    Parse the CS tag for a read and record the variations observed for the read also for the loci

    :param cooper: ATaRVa object
    :param read: pysam AlignedSegment object
    """

    if cooper.cooper_sorted_snps == None: cooper.cooper_sorted_snps = []
    operations = {':', '-', '+', '*', '=', '~'}

    rpos = read.ref_start   # NOTE: The coordinates are 1 based in SAM
    qpos = 0                # starts from 0 the sub string the read sequence in python

    chrom = read.chrom
    repeat_index = 0

    locus_query_range = np.zeros((len(read.loci_coords), 2), dtype=int)
    flank_query_range = np.zeros((len(read.loci_coords), 2), dtype=int)
    left_flank_insertions  = [[] for _ in read.loci_coords] # stores insertions in left flank as (rpos, qstart, qend)
    right_flank_insertions = [[] for _ in read.loci_coords] # stores insertions in right flank as (rpos, qstart, qend)
    locus_reached = [False for _ in read.loci_coords]
    locus_boundary_crossed = [[False,False] for _ in read.loci_coords]

    insert_positions = {}
    if read.cigartuples[0][0] == 4:     # adjust of softclip
        qpos += read.cigartuples[0][1]

    i = 0; cs_len = len(read.cs_tag)
    while i < cs_len:

        if read.cs_tag[i] == ':':        # sequence match in short CS is followed by the length of match
            match_len = ''; i += 1
            while i < cs_len and read.cs_tag[i] not in operations:
                match_len += read.cs_tag[i]; i += 1

            match_len = int(match_len)
            qpos += match_len; rpos += match_len
            repeat_index += match_jump(read, rpos, qpos, match_len, repeat_index, locus_query_range, flank_query_range,
                                       locus_reached, locus_boundary_crossed)

        elif read.cs_tag[i] == '=':      # sequence match in long CS is followed by nucs which are matching       
            match_len = 0
            while i < cs_len and read.cs_tag[i] not in operations:
                match_len += 1; i += 1

            qpos += match_len; rpos += match_len
            repeat_index += match_jump(read, rpos, qpos, match_len, repeat_index, locus_query_range, flank_query_range,
                                       locus_reached, locus_boundary_crossed)


        elif read.cs_tag[i] == '*':      # substitution of a base; is followed by reference and substituted base
            ref_nuc, sub_nuc = read.cs_tag[i+1], read.cs_tag[i+2]
            i += 3

            match_len = 1
            update_snps(cooper, read, rpos, qpos, insert_positions, False)

            qpos += match_len; rpos += match_len
            repeat_index += match_jump(read, rpos, qpos, match_len, repeat_index, locus_query_range, flank_query_range,
                                       locus_reached, locus_boundary_crossed)

        elif read.cs_tag[i] == '+':      # insertion; is followed by the inserted bases
            insert = ''; insert_len = 0; i += 1
            while i < cs_len and read.cs_tag[i] not in operations:
                insert += read.cs_tag[i]; insert_len += 1 
                i += 1

            insert_positions[rpos] = insert_len
            homopolymer_insert = False
            if len(set(read.query_sequence[qpos:qpos+insert_len])) == 1: homopolymer_insert = True

            qpos += insert_len
            repeat_index += insertion_jump(read, rpos, qpos, insert_len, homopolymer_insert, repeat_index, locus_query_range, flank_query_range,
                                           locus_reached, locus_boundary_crossed, left_flank_insertions, right_flank_insertions)

        elif read.cs_tag[i] == '-':      # deletion; is followed by the deleted bases
            deletion = ''; deletion_len = 0; i += 1
            while i < cs_len and read.cs_tag[i] not in operations:
                deletion += read.cs_tag[i]; deletion_len += 1
                i += 1
            if not cooper.haploid:
                cooper.cooper_read_data[read.index].dels.extend([rpos, rpos + deletion_len])
            rpos += deletion_len
            repeat_index += deletion_jump(read, rpos, qpos, deletion_len, repeat_index, locus_query_range, flank_query_range, locus_reached, locus_boundary_crossed)

            
    num_read_loci = len(read.loci_coords)
    for idx, locus_key in enumerate(read.loci_keys):
        # changing all the query coordinates to be relative to the start of the flank start
        flank_query_start = flank_query_range[idx][0]
        if idx == 0: read.methyl_start = flank_query_start

        flank_query_end = flank_query_range[idx][1]
        if idx == num_read_loci - 1: read.methyl_end = flank_query_end

        locus_query_range[idx][0] = locus_query_range[idx][0] - flank_query_start
        locus_query_range[idx][1] = locus_query_range[idx][1] - flank_query_start

        left_flank_insertions[idx]  = [(coords[0], coords[1] - flank_query_start, coords[2] - flank_query_start) for coords in left_flank_insertions[idx] ]
        right_flank_insertions[idx] = [(coords[0], coords[1] - flank_query_start, coords[2] - flank_query_start) for coords in right_flank_insertions[idx] ]
        read.loci_data[locus_key].seq = [read.query_sequence[flank_query_start:flank_query_end],
                                         locus_query_range[idx],
                                         left_flank_insertions[idx],
                                         right_flank_insertions[idx],
                                         flank_query_start, flank_query_end]