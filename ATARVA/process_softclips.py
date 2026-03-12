from ATARVA.realignment import *

def detect_flank(cooper, read, locus_start, locus_end, stream):
    """
    Identify the flanking sequences of the repeat region in the read
    Args:
        read: pysam read object
        fasta: pysam fasta object
        flank_length: length of the flanking sequence to be considered
        locus_start: start coordinate of the repeat in reference
        locus_end: end coordinate of the repeat in reference

    Returns:
        [start, end, sub_cigar]
        start and end coordinates of the flanking sequence in the read
        sub_cigar: CIGAR string of the repeat alignment to the reference in the read
    """

    flank = cooper.args.flank

    if stream: # if its upstream (True)

        # finding the upstream flank in the read 
        upstream = cooper.ref.fetch(read.chrom, locus_start - cooper.args.flank, locus_start)
        alignment_score, target_begin, target_end, query_begin, query_end, sCigar = stripSW(Inputs(read.query_sequence, upstream), False)

        # target_begin is the start of the alignment in the read
        # query_begin  is the start of the alignment in the reference flank
        # target_end   is the end of the alignment in the read
        # query_end    is the end of the alignment in the reference flank

        # reject read if less than 70% match with the reference flank sequence
        if alignment_score < int(2*(0.7*flank)): return [0, False]

        # reference position where the read (upstream flank) is starts to align
        read_reference_start = locus_start - flank + query_begin - 1

        # if the flank aligns completely this should be equal to the repeat start
        uflank_reference_end = read_reference_start + query_end - query_begin + 1

        # remaining bases in the flank that are not present in the read; these could be either deleted or substituted
        # these bases are at the end of the upstream flank
        uflank_remaining_length = flank_length - query_end
        uflank_read_end = target_end - 1

        # adding the remaining bases of the upstream flank as substitutions
        if uflank_remaining_length > 0:
            # sub_cigar += f'{uflank_remaining_length}X'
            uflank_reference_end += uflank_remaining_length
            uflank_read_end += uflank_remaining_length

        if uflank_read_end >= read.query_start:
            # if the upstream flank aligns to the already aligned region(by aligner), there is no expansion in the soft-clip
            return [uflank_read_end, False] # False means no expansion in the soft-clip, here uflank_read_end is dummy
        else:
            # upstream flank is beyond the true aligned region, so there is expansion in the soft-clip
            return [uflank_read_end+1, True] # True means expansion in the soft-clip, +1 will give start_index in read_sequence
        
    else: # if its downstream (False)

        # finding the downstream flank in the read
        downstream = cooper.ref.fetch(read.chrom, locus_end, locus_end + flank)
        alignment_score, target_begin, target_end, query_begin, query_end, sCigar = stripSW(Inputs(read.query_sequence, downstream), False)

        # target_begin is the start of the alignment in the read
        # query_begin  is the start of the alignment in the reference flank
        # target_end   is the end of the alignment in the read
        # query_end    is the end of the alignment in the reference flank

        # reject read if less than 70% match with the reference flank sequence
        if alignment_score < int(2*(0.7*flank)): return [0, False]

        # should be qual to the repeat end if the flank aligns completely
        dflank_reference_start = locus_end + query_begin - 1

        # remaining bases in the flank that are not present in the read; these could be either deleted or substituted
        # these bases are at the start of downstream flank
        dflank_remaining_length = flank - query_end # number of bases in the flank that are not aligned
        dflank_read_start = target_begin - 1

        # adding the remaining bases of the downstream flank as substitutions
        if dflank_reference_start > locus_end:
            # dflank_cigar += f'{dflank_remaining_length}X'     # adding the remaining bases as substitutions
            dflank_reference_start -= (dflank_remaining_length)
            dflank_read_start -= (dflank_remaining_length)

        if dflank_read_start <= read.query_end:
            # if the downstream flank aligns to the already aligned region(by aligner), there is no expansion in the soft-clip
            return [dflank_read_start, False] # False means no expansion in the soft-clip, here dflank_read_start is dummy
        else:
            # downstream flank is beyond the true aligned region, so there is expansion in the soft-clip
            return [dflank_read_start, True] # True means expansion in the soft-clip