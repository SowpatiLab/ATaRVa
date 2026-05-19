import bisect
import numpy as np

from ATARVA.operation_utils import match_jump, deletion_jump, insertion_jump
from ATARVA.md_utils import parse_mdtag
from ATARVA.structures import SNP

def subex(ref, query):
    """
    compare reference and query sequence and get mismatch positions

    :param ref: reference sequence
    :param query: query sequence
    :return: list of mismatch positions
    """

    if '=' in query:
        array = np.frombuffer(query.encode(), dtype=np.byte)
        substitution_indices = np.where(array != ord('='))[0]
        
    else:
        array1 = np.frombuffer(ref.encode(), dtype=np.byte)
        array2 = np.frombuffer(query.encode(), dtype=np.byte)
        substitution_indices = np.where(array1 != array2)[0]
        
    return substitution_indices


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

    rpos = read.ref_start   # NOTE: The coordinates are 1 based in SAM
    qpos = 0                # starts from 0 the sub string the read sequence in python

    chrom = read.chrom
    repeat_index = 0

    locus_query_range = [[0, 0] for _ in read.loci_coords]
    flank_query_range = [[0, 0] for _ in read.loci_coords]
    left_flank_insertions  = [[] for _ in read.loci_coords] # stores insertions in left flank as (rpos, qstart, qend)
    right_flank_insertions = [[] for _ in read.loci_coords] # stores insertions in right flank as (rpos, qstart, qend)
    locus_reached = [False for _ in read.loci_coords]
    locus_boundary_crossed = [[False,False] for _ in read.loci_coords]

    has_subop = False
    has_MD = read.has_tag('MD')

    insert_positions = {}

    haploid              = cooper.haploid
    snp_data             = cooper.cooper_snp_data
    cooper_prev_reads    = cooper.prev_reads
    is_haplotag          = cooper.args.haplotag
    reference            = cooper.ref
    snp_qual_cutoff      = cooper.args.snp_qual
    cooper_sorted_snps   = cooper.cooper_sorted_snps
    cooper_loci_keys     = cooper.cooper_loci_keys
    cooper_loci_info     = cooper.cooper_loci_info
    cooper_snp_data      = cooper.cooper_snp_data
    cooper_read_data     = cooper.cooper_read_data
    is_haploid           = cooper.haploid

    read_index = read.index
    qseq       = read.query_sequence
    qquals     = read.query_qualities

    no_snp_range = 3

    for cigar in read.cigartuples:
        if cigar[0] == 4: qpos += cigar[1]    # softclip; adjust the qpos but not rpos; also consider the possibility of match after softclip
 
        elif cigar[0] == 2:     # deletion
            deletion_len = cigar[1]
            if not haploid:
                cooper_read_data[read_index].dels.extend([rpos, rpos + deletion_len])
                cooper_read_data[read_index].no_snps.update(range(rpos-no_snp_range, rpos + deletion_len + 1 + no_snp_range))
            rpos += deletion_len
            repeat_index += deletion_jump(read, rpos, qpos, deletion_len, repeat_index, locus_query_range, flank_query_range, locus_reached, locus_boundary_crossed)

        elif cigar[0] == 1:     # insertion
            insert_positions[rpos] = cigar[1]
            cooper_read_data[read_index].no_snps.update(range(rpos-no_snp_range, rpos + 1 + no_snp_range))
            insert_len = cigar[1]
            homopolymer_insert = False
            if len(set(qseq[qpos:qpos+insert_len])) == 1: homopolymer_insert = True

            qpos += insert_len
            repeat_index += insertion_jump(read, rpos, qpos, insert_len, homopolymer_insert, repeat_index, locus_query_range, flank_query_range,
                                           locus_reached, locus_boundary_crossed, left_flank_insertions, right_flank_insertions)
        
        elif cigar[0] == 0: # match (includes substitutions)
            match_len = cigar[1]
            if not has_MD and (not haploid) and (not is_haplotag):
                ref_sequence   = reference.fetch(chrom, rpos, rpos + match_len)
                query_sequence = qseq[qpos:qpos + match_len]

                assert len(ref_sequence) == len(query_sequence), \
                    "Error: fetching sequences based on CIGAR. Please check cigar formats"

                substitutions = subex(ref_sequence, query_sequence)
                for sub_pos in substitutions:
                    if not outside_locus(read.loci_coords, rpos + sub_pos):
                        continue
                    sub_nuc  = query_sequence[sub_pos]
                    sub_qual = qquals[qpos + sub_pos]
                    if sub_qual < snp_qual_cutoff: continue

                    sub_pos = rpos + sub_pos
                    cooper_read_data[read_index].snps.add(sub_pos)
                    if sub_pos not in snp_data:
                        snp_data[sub_pos] = SNP(cov = 1, sub = { sub_nuc: {read_index} }, qual={ read_index: sub_qual })
                        cooper_sorted_snps.add(sub_pos)
                    else:
                        snp_data[sub_pos].cov += 1
                        snp_data[sub_pos].qual[read_index] = sub_qual
                        if sub_nuc in snp_data[sub_pos].sub:
                            snp_data[sub_pos].sub[sub_nuc].add(read_index)
                        else:
                            snp_data[sub_pos].sub[sub_nuc] = {read_index}

            qpos += match_len; rpos += match_len
            repeat_index += match_jump(read, rpos, qpos, match_len, repeat_index, locus_query_range, flank_query_range, locus_reached, locus_boundary_crossed)

        elif cigar[0] == 7: # exact match (does not include substitutions)
            has_subop = True
            match_len = cigar[1]
            qpos += match_len; rpos += match_len
            repeat_index += match_jump(read, rpos, qpos, match_len, repeat_index, locus_query_range, flank_query_range, locus_reached, locus_boundary_crossed)

        elif cigar[0] == 8: # substitution (difference)
            has_subop = True
            match_len = cigar[1]
            if (not haploid) and outside_locus(read.loci_coords, rpos) and (not is_haplotag):
                sub_nuc  = qseq[qpos]
                sub_qual = qquals[qpos]
                if sub_qual >= snp_qual_cutoff:
                    cooper_read_data[read_index].snps.add(rpos)
                    if rpos not in snp_data:
                        snp_data[rpos] = SNP(cov = 1, sub = { sub_nuc: {read_index} }, qual={ read_index: sub_qual })
                        cooper_sorted_snps.add(rpos)
                    else:
                        snp_data[rpos].cov += 1
                        snp_data[rpos].qual[read_index] = sub_qual
                        if sub_nuc in snp_data[rpos].sub: 
                            snp_data[rpos].sub[sub_nuc].add(read_index)
                        else:
                            snp_data[rpos].sub[sub_nuc] = {read_index}

            qpos += match_len; rpos += match_len
            repeat_index += match_jump(read, rpos, qpos, match_len, repeat_index, locus_query_range, flank_query_range, locus_reached, locus_boundary_crossed)

    if not has_subop :
        if read.has_tag('MD'):
            if read.cigartuples[0][0] == 4: qpos = read.cigartuples[0][1]
            else: qpos = 0
            parse_mdtag(cooper_snp_data, cooper_sorted_snps, cooper_loci_keys, cooper_loci_info, cooper_read_data, cooper_prev_reads,
                        snp_qual_cutoff, is_haploid, is_haplotag, read, qpos, insert_positions)

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
        read.loci_data[locus_key].seq = [qseq[flank_query_start:flank_query_end],
                                         locus_query_range[idx],
                                         left_flank_insertions[idx],
                                         right_flank_insertions[idx],
                                         flank_query_start, flank_query_end]
