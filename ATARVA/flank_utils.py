from itertools import product


def check_flank_order(coords, chrom, start, end):
    """
    Validate that softclip flank coordinates are in reference order and that the
    corresponding query coordinates preserve the same sequential order.

    Each entry in ``softclip_loci_coords`` is expected to be a dict with
    ``'upstream'`` and ``'downstream'`` keys mapped to 4-tuples of the form:
    ``(ref_start, ref_end, query_start, query_end)``.

    :param softclip_loci_coords: list of flank coordinate dictionaries
    :param softclip_loci: list of locus coordinate tuples
    :return: True when every locus has ordered upstream/downstream reference and
             query coordinates, False otherwise
    """

    upstream   = coords['upstream']   if coords else None
    downstream = coords['downstream'] if coords else None

    if upstream is None and downstream is None:
        return False

    # if either of them is None we assume it is in order
    elif upstream is None or downstream is None:
        return True

    upstream_ref_start, upstream_ref_end, upstream_query_start, upstream_query_end = upstream
    downstream_ref_start, downstream_ref_end, downstream_query_start, downstream_query_end = downstream

    ref_in_order = (
        upstream_ref_start <= upstream_ref_end <= downstream_ref_start <= downstream_ref_end
    )
    query_in_order = (
        upstream_query_start <= upstream_query_end <= downstream_query_start <= downstream_query_end
    )

    if not ref_in_order:
        # raise ValueError(f"Reference flank coordinates are not in order for locus {chrom}:{start}-{end}")
        return False
    if not query_in_order:
        return False

    return True


def flag_misaligned_loci(softclip_loci):
    """
    Identify loci that break continuous ascending order of query coordinates.
    
    Assumes loci are in continuous ascending order in the read. Any locus whose
    query start position is <= the previous locus is marked for removal.
    Uses a greedy approach to keep the maximum number of loci in correct order.
    
    Example: query positions [1,2,3,4,5,2,3,6,7,8,1,2,10] 
             -> removes indices of [2,3,1,2] (positions that restart or go backward)
             -> keeps [1,2,3,4,5,6,7,8,10]
    
    :param softclip_loci: list of dicts with 'upstream' and 'downstream' coordinates
    :return: set of indices that should be dropped
    """
    
    if not softclip_loci['coords'] or len(softclip_loci['coords']) < 2:
        return
    
    # Extract upstream query starts (read positions)
    indices    = []
    ups_starts = []
    prev_uqe   = -1
    for i, coords in enumerate(softclip_loci['coords']):
        if softclip_loci['flags'][i] is not None: continue

        if prev_uqe == -1 and coords['upstream'] is None:
            softclip_loci['flags'][i] = 'INVALID_UPSTREAM'
            continue

        _us = False
        if coords['upstream'] is not None:
            _us = True
            _, _, uqs, uqe = coords['upstream']  # extract query start
            indices.append(i)
            ups_starts.append(uqs)
        
        elif coords['downstream'] is not None:
            _, _, dqs, dqe = coords['downstream']
            if dqs > prev_uqe:
                indices.append(i)
                ups_starts.append(prev_uqe)

        if _us: prev_uqe = uqe
    
    if len(ups_starts) == 0:
        return

    orders     = [[ups_starts[0]]]
    order_idxs = [[indices[0]]]
    for i,uqs in enumerate(ups_starts[1:]):
        j = i + 1
        check = False
        for o, order in enumerate(orders):
            if uqs >= order[-1]:
                orders[o].append(uqs)
                order_idxs[o].append(indices[j])
                check = True
        if not check:
            orders.append([uqs])
            order_idxs.append([indices[j]])

    max_order = []
    max_order_idxs = []
    for o, order in enumerate(orders):
        if len(order) > len(max_order):
            max_order = order
            max_order_idxs = order_idxs[o]

    for idx in range(len(softclip_loci['coords'])):
        if idx not in max_order_idxs and softclip_loci['flags'][idx] is None:
            softclip_loci['flags'][idx] = 'INVALID_ORDER'


def _coordinates_overlap(coord1, coord2):
    """
    Check if two coordinate tuples (start, end) overlap.

    :param coord1: tuple (start, end)
    :param coord2: tuple (start, end)
    :return: True if overlapping, False otherwise
    """
    if not coord1 or not coord2:
        return False
    start1, end1 = coord1
    start2, end2 = coord2

    # if separated by less than a kb merge them
    # if end1 < start2 and start2 - end1 <= 0:
    #     return True
    
    if start2 <= end1 and end1 <= end2:
        return True

    return not (end1 < start2 or end2 < start1)


def _merge_coordinates(coord1, coord2, upstream=True):
    """
    Merge two overlapping coordinate tuples.

    :param coord1: tuple (start, end)
    :param coord2: tuple (start, end)
    :return: merged tuple (min_start, max_end)
    """
    if not coord1:
        return coord2
    if not coord2:
        return coord1
    start_ref1, end_ref1, start_query1, end_query1 = coord1
    start_ref2, end_ref2, start_query2, end_query2 = coord2

    if upstream:
        if start_query1 < start_query2:
            return coord1
        elif start_query1 > start_query2:
            return coord2
        elif start_query1 == start_query2:
            if start_ref1 < start_ref2:
                return coord1
            else:
                return coord2
    else:
        if end_query1 > end_query2:
            return coord1
        elif end_query1 < end_query2:
            return coord2
        elif end_query1 == end_query2:
            if end_ref1 > end_ref2:
                return coord1
            else:
                return coord2
            

def check_ref_query_order_consistency(current_coords, next_coords, merged_i, merged_j):
    """
    Check if the order of query coordinates matches the order of reference 
    coordinates for current and next loci. Compares current upstream with 
    next upstream and downstream, and current downstream with next upstream 
    and downstream (for all non-None pairs).

    :param current_coords: dict with 'upstream' and 'downstream' coordinate tuples
    :param next_coords: dict with 'upstream' and 'downstream' coordinate tuples
    :return: tuple (locus_to_drop, reason) where locus_to_drop is 'current', 'next', or None
             reason is a string describing why a locus was dropped, or 'ORDER_CONSISTENT' if valid
    """
    
    if current_coords is None or next_coords is None:
        return None, 'INVALID_COORDS'
    
    curr_upstream   = current_coords.get('upstream')
    curr_downstream = current_coords.get('downstream')
    next_upstream   = next_coords.get('upstream')
    next_downstream = next_coords.get('downstream')
    
    # List to track order mismatches
    mismatch_count = 0
    total_comparisons = 0
    
    # Compare current_upstream with next_upstream
    if curr_upstream is not None and next_upstream is not None:
        curr_urs, curr_ure, curr_uqs, curr_uqe = curr_upstream
        next_urs, next_ure, next_uqs, next_uqe = next_upstream
        
        curr_ref_order   = curr_urs < next_urs
        curr_query_order = curr_uqs < next_uqs
        total_comparisons += 1
        if curr_ref_order != curr_query_order:
            mismatch_count += 1
    
    # Compare current_upstream with next_downstream
    if curr_upstream is not None and next_downstream is not None:
        curr_urs, curr_ure, curr_uqs, curr_uqe = curr_upstream
        next_drs, next_dre, next_dqs, next_dqe = next_downstream
        
        curr_ref_order   = curr_urs < next_drs
        curr_query_order = curr_uqs < next_dqs
        total_comparisons += 1
        if curr_ref_order != curr_query_order:
            mismatch_count += 1
    
    # Compare current_downstream with next_upstream
    if curr_downstream is not None and next_upstream is not None:
        curr_drs, curr_dre, curr_dqs, curr_dqe = curr_downstream
        next_urs, next_ure, next_uqs, next_uqe = next_upstream
        
        curr_ref_order   = curr_drs < next_urs
        curr_query_order = curr_dqs < next_uqs
        total_comparisons += 1
        if curr_ref_order != curr_query_order:
            mismatch_count += 1
    
    # Compare current_downstream with next_downstream
    if curr_downstream is not None and next_downstream is not None:
        curr_drs, curr_dre, curr_dqs, curr_dqe = curr_downstream
        next_drs, next_dre, next_dqs, next_dqe = next_downstream
        
        curr_ref_order   = curr_drs < next_drs
        curr_query_order = curr_dqs < next_dqs
        total_comparisons += 1
        if curr_ref_order != curr_query_order:
            mismatch_count += 1
    
    # If no comparisons were made, return as consistent
    if total_comparisons == 0:
        return None, 'INSUFFICIENT_COORDS'
    
    # If all comparisons match, order is consistent
    if mismatch_count == 0:
        return None, 'ORDER_CONSISTENT'
    
    # If there are mismatches, calculate which locus to drop based on length discrepancy
    if curr_upstream is not None and curr_downstream is not None:
        curr_urs, curr_ure, curr_uqs, curr_uqe = curr_upstream
        curr_drs, curr_dre, curr_dqs, curr_dqe = curr_downstream
        curr_ref_length = curr_dre - curr_urs
        curr_query_length = curr_dqe - curr_uqs
        curr_discrepancy = abs(curr_query_length - curr_ref_length)
    else:
        curr_discrepancy = float('inf')
    
    if next_upstream is not None and next_downstream is not None:
        next_urs, next_ure, next_uqs, next_uqe = next_upstream
        next_drs, next_dre, next_dqs, next_dqe = next_downstream
        next_ref_length = next_dre - next_urs
        next_query_length = next_dqe - next_uqs
        next_discrepancy = abs(next_query_length - next_ref_length)
    else:
        next_discrepancy = float('inf')
    
    # Drop the one with larger discrepancy
    if merged_i == True:
        return 'next', 'ORDER_MISMATCH_NEXT_WORSE'
    elif curr_discrepancy > next_discrepancy:
        return 'current', 'ORDER_MISMATCH_CURRENT_WORSE'
    elif next_discrepancy > curr_discrepancy:
        return 'next', 'ORDER_MISMATCH_NEXT_WORSE'
    else:
        # If discrepancies are equal, drop current by default
        return 'next', 'ORDER_MISMATCH_EQUAL_DROP_NEXT'


def process_flanks(softclip_loci, read_ref_start, read_ref_end):
    """
    checks the order of alignment coords of the flanks and merges them

    :param softclip_loci: list of locus coordinates
    :param softclip_loci_coords: list of dicts with 'upstream' and 'downstream' coordinates
    :param read_ref_end: end position of the read in the reference
    :return: tuple of (merged_loci, merged_coords)
    
    NOTE:
    urs = upstream reference start;   ure = upstream reference end
    drs = downstream reference start; dre = downstream reference end
    uqs = upstream query start;       uqe = upstream query end
    dqs = downstream query start;     dqe = downstream query end
    """

    merged_coords = []

    # Identify and flag loci that are in the wrong order compared to the majority
    flag_misaligned_loci(softclip_loci)

    # if there is only one qualifying locus with flanks
    if sum([flag is None for flag in softclip_loci['flags']]) == 1:
        check = False
        for i, coords in enumerate(softclip_loci['coords']):
            if softclip_loci['flags'][i] is None:
                if coords['upstream'] is not None and coords['downstream'] is not None:
                    check = True
                    merged_coords.append(coords)
                    break
                else:
                    softclip_loci['flags'][i] = 'FLANK_MISSING'
        if not check:
            return []
        else: return merged_coords
    
    elif sum([flag is None for flag in softclip_loci['flags']]) == 0:
        return []


    merged = [False] * len(softclip_loci['coords'])
    i = 0
    while i < len(softclip_loci['coords']):
        if softclip_loci['flags'][i] is not None:
            i += 1
            continue
        
        current_coords = softclip_loci['coords'][i]
        current_locus  = softclip_loci['loci'][i]

        # a locus with no upstream flank either is merged with the locus before or is skipped
        if current_coords['upstream'] is None:
            softclip_loci['flags'][i] = 'MISSING_UPFLANK'
            i += 1
            continue

        # Check for overlaps with following loci
        j = i + 1
        dir_switch = False
        while j < len(softclip_loci['coords']) and not dir_switch:
            next_coords = softclip_loci['coords'][j]
            next_locus  = softclip_loci['loci'][j]

            if softclip_loci['flags'][i] is not None:
                j = i + 1
                break
            
            if softclip_loci['flags'][j] is not None:
                j += 1
                continue

            # softclip_dir = 'upstream' if ref_end < read.ref_end else 'downstream'
            if current_coords['upstream'][1] < read_ref_end:
                check = False
                if next_coords['downstream'] is not None and next_coords['downstream'][1] > read_ref_end:
                    dir_switch = True
                    check = True
                elif next_coords['upstream'] is not None and (next_coords['upstream'][1] > read_ref_end or (next_coords['upstream'][0] == read_ref_end and next_coords['upstream'][1] == read_ref_end)):
                    check = True
                if check: continue
            
            locus_to_drop, reason = check_ref_query_order_consistency(current_coords, next_coords, merged[i], merged[j])

            if locus_to_drop == 'current':
                softclip_loci['flags'][i] = reason
                j = i + 1
                break

            elif locus_to_drop == 'next':
                softclip_loci['flags'][j] = reason
                j += 1
                continue

            curr_urs, curr_ure, curr_uqs, curr_uqe = current_coords['upstream'] if current_coords['upstream'] else (None, None, None, None)
            next_urs, next_ure, next_uqs, next_uqe = next_coords['upstream'] if next_coords['upstream'] else (None, None, None, None)
            curr_drs, curr_dre, curr_dqs, curr_dqe = current_coords['downstream'] if current_coords['downstream'] else (None, None, None, None)
            next_drs, next_dre, next_dqs, next_dqe = next_coords['downstream'] if next_coords['downstream'] else (None, None, None, None)
            
            # check if there is a missing centre flank and merge the flanks if the coordinates are in the correct order
            missing_centre = current_coords['upstream'] is not None and (current_coords['downstream'] is None or next_coords['upstream'] is None) and next_coords['downstream'] is not None
            if missing_centre:
                current_coords['downstream'] = next_coords['downstream']
                current_locus = softclip_loci['loci'][j]
                merged[i] = True; merged[j] = True
                j += 1
                continue

            # if the downstream flank of the next coordinate is missing and the downstream flank of the current coordinate is missing, loci cannot be merged
            # hence skipped with flag of FLANK_MISSING
            if next_coords['downstream'] is None:
                j += 1
                continue

            ref_overlap   = _coordinates_overlap((curr_urs, curr_dre), (next_urs, next_dre))
            query_overlap = _coordinates_overlap((curr_uqs, curr_dqe), (next_uqs, next_dqe))

            if ref_overlap or query_overlap:
                # Merge coordinates
                current_coords['upstream']   = _merge_coordinates(current_coords['upstream'], next_coords['upstream'], True)
                current_coords['downstream'] = _merge_coordinates(current_coords['downstream'], next_coords['downstream'], False)
                current_locus = softclip_loci['loci'][j]
                merged[i] = True; merged[j] = True
                j += 1
            else:
                break

        if current_coords['upstream'] is not None and current_coords['downstream'] is not None and softclip_loci['flags'][i] is None:
            merged_coords.append(current_coords)

        i = j if j > i + 1 else i + 1

    return merged_coords



# def check_ref_query_order_consistency(current_coords: dict,
#                                        next_coords:    dict,
#                                        merged_i:       bool) -> tuple[str | None, str]:
#     """
#     Check if query coordinate order matches reference coordinate order
#     across current and next loci flanking regions.

#     :param current_coords: dict with 'upstream' and 'downstream' coord tuples
#                            each tuple: (ref_start, ref_end, query_start, query_end)
#     :param next_coords:    same structure for next locus
#     :param merged_i:       if True — always drop next
#     :return:               (locus_to_drop, reason)
#                            locus_to_drop: 'current', 'next', or None
#     """
#     if current_coords is None or next_coords is None:
#         return None, 'INVALID_COORDS'

#     # ── extract coordinates ───────────────────────────────────────────
#     curr_up   = current_coords.get('upstream')
#     curr_down = current_coords.get('downstream')
#     next_up   = next_coords.get('upstream')
#     next_down = next_coords.get('downstream')

#     # ── count order mismatches across all valid pairs ─────────────────
#     pairs = [
#         (curr_up,   next_up),
#         (curr_up,   next_down),
#         (curr_down, next_up),
#         (curr_down, next_down),
#     ]

#     total      = 0
#     mismatches = 0

#     for a, b in pairs:
#         if a is None or b is None:
#             continue
#         total += 1
#         # ref_start is index 0, query_start is index 2
#         ref_order   = a[0] < b[0]
#         query_order = a[2] < b[2]
#         if ref_order != query_order:
#             mismatches += 1

#     # ── early exits ───────────────────────────────────────────────────
#     if total == 0:
#         return None, 'INSUFFICIENT_COORDS'

#     if mismatches == 0:
#         return None, 'ORDER_CONSISTENT'

#     # ── compute length discrepancy per locus ──────────────────────────
#     def discrepancy(up, down) -> float:
#         """Absolute difference between ref and query span lengths."""
#         if up is None or down is None:
#             return float('inf')
#         ref_len   = down[1] - up[0]   # down.ref_end   - up.ref_start
#         query_len = down[3] - up[2]   # down.query_end - up.query_start
#         return abs(query_len - ref_len)

#     # ── decide which locus to drop ────────────────────────────────────
#     if merged_i:
#         return 'next', 'ORDER_MISMATCH_NEXT_WORSE'

#     curr_disc = discrepancy(curr_up, curr_down)
#     next_disc = discrepancy(next_up, next_down)

#     if curr_disc > next_disc:
#         return 'current', 'ORDER_MISMATCH_CURRENT_WORSE'
#     elif next_disc > curr_disc:
#         return 'next',    'ORDER_MISMATCH_NEXT_WORSE'
#     else:
#         return 'next',    'ORDER_MISMATCH_EQUAL_DROP_NEXT'
