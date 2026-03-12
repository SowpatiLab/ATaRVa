

def check_flank(cooper, read, locus_start, locus_end):
    flank = cooper.args.flank
    enough_soft_flank = False
    left_flank = None; right_flank = None

    # read alignment covers the locus with enough soft clipped flanks on both sides
    if read.ref_start < locus_start - flank and read.ref_end > locus_end + flank:
        enough_soft_flank = True
        left_flank  = flank
        right_flank = flank

    # for upstream soft_clip, checking the length compared to flank length
    elif locus_start < read.ref_start:
        mod_locus_start = read.ref_start
        left_flank = min(flank, locus_start - read.ref_start)

    # for downstream soft_clip, checking the length compared to flank length
    elif locus_end > read.ref_end:
        mod_locus_end = read.ref_end
        right_flank = min(flank, read.ref_end - locus_end)

    return enough_soft_flank, left_flank, right_flank