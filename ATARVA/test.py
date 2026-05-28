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

loci = {'keys': ['chr1:143249365-143249376', 'chr1:143249453-143249461', 'chr1:143249654-143249663', 'chr1:143249921-143249930', 'chr1:143250000-143250008', 'chr1:143250101-143250109', 'chr1:143250389-143250398', 'chr1:143250493-143250502', 'chr1:143250569-143250578', 'chr1:143251209-143251217', 'chr1:143251333-143251341', 'chr1:143251359-143251367', 'chr1:143251434-143251442', 'chr1:143251722-143251731', 'chr1:143251748-143251757', 'chr1:143251826-143251835', 'chr1:143251904-143251913', 'chr1:143252567-143252575', 'chr1:143252593-143252603', 'chr1:143252642-143252650', 'chr1:143252743-143252751', 'chr1:143253136-143253144', 'chr1:143253525-143253534', 'chr1:143253656-143253663', 'chr1:143254196-143254204', 'chr1:143254222-143254230', 'chr1:143254323-143254331', 'chr1:143254615-143254624', 'chr1:143254664-143254673', 'chr1:143254817-143254825', 'chr1:143254843-143254851', 'chr1:143254892-143254900', 'chr1:143254918-143254926', 'chr1:143254993-143255001', 'chr1:143255272-143255281', 'chr1:143255298-143255307', 'chr1:143255454-143255463', 'chr1:143255913-143255922', 'chr1:143255962-143255971', 'chr1:143256028-143256039', 'chr1:143256115-143256123', 'chr1:143256141-143256151', 'chr1:143256216-143256224', 'chr1:143256291-143256299', 'chr1:143256579-143256588', 'chr1:143256605-143256614', 'chr1:143256761-143256770', 'chr1:143257221-143257230', 'chr1:143257336-143257347', 'chr1:143257548-143257556', 'chr1:143257574-143257582', 'chr1:143257937-143257946', 'chr1:143258408-143258417', 'chr1:143258434-143258443', 'chr1:143258512-143258521', 'chr1:143258590-143258599', 'chr1:143259318-143259326', 'chr1:143259606-143259615'], 'loci': [('chr1', 143249365, 143249376), ('chr1', 143249453, 143249461), ('chr1', 143249654, 143249663), ('chr1', 143249921, 143249930), ('chr1', 143250000, 143250008), ('chr1', 143250101, 143250109), ('chr1', 143250389, 143250398), ('chr1', 143250493, 143250502), ('chr1', 143250569, 143250578), ('chr1', 143251209, 143251217), ('chr1', 143251333, 143251341), ('chr1', 143251359, 143251367), ('chr1', 143251434, 143251442), ('chr1', 143251722, 143251731), ('chr1', 143251748, 143251757), ('chr1', 143251826, 143251835), ('chr1', 143251904, 143251913), ('chr1', 143252567, 143252575), ('chr1', 143252593, 143252603), ('chr1', 143252642, 143252650), ('chr1', 143252743, 143252751), ('chr1', 143253136, 143253144), ('chr1', 143253525, 143253534), ('chr1', 143253656, 143253663), ('chr1', 143254196, 143254204), ('chr1', 143254222, 143254230), ('chr1', 143254323, 143254331), ('chr1', 143254615, 143254624), ('chr1', 143254664, 143254673), ('chr1', 143254817, 143254825), ('chr1', 143254843, 143254851), ('chr1', 143254892, 143254900), ('chr1', 143254918, 143254926), ('chr1', 143254993, 143255001), ('chr1', 143255272, 143255281), ('chr1', 143255298, 143255307), ('chr1', 143255454, 143255463), ('chr1', 143255913, 143255922), ('chr1', 143255962, 143255971), ('chr1', 143256028, 143256039), ('chr1', 143256115, 143256123), ('chr1', 143256141, 143256151), ('chr1', 143256216, 143256224), ('chr1', 143256291, 143256299), ('chr1', 143256579, 143256588), ('chr1', 143256605, 143256614), ('chr1', 143256761, 143256770), ('chr1', 143257221, 143257230), ('chr1', 143257336, 143257347), ('chr1', 143257548, 143257556), ('chr1', 143257574, 143257582), ('chr1', 143257937, 143257946), ('chr1', 143258408, 143258417), ('chr1', 143258434, 143258443), ('chr1', 143258512, 143258521), ('chr1', 143258590, 143258599), ('chr1', 143259318, 143259326), ('chr1', 143259606, 143259615)], 'coords': [{'upstream': (143249370, 143249370, 28061, 28061), 'downstream': (143249376, 143249405, 31319, 31348)}, {'upstream': None, 'downstream': (143249461, 143249489, 28200, 28230)}, {'upstream': (143249624, 143249652, 28213, 28241), 'downstream': (143249663, 143249692, 28814, 28843)}, {'upstream': (143249891, 143249920, 28577, 28606), 'downstream': (143249930, 143249959, 29081, 29110)}, {'upstream': (143249970, 143249999, 29121, 29150), 'downstream': (143250008, 143250037, 29159, 29188)}, {'upstream': None, 'downstream': (143250109, 143250138, 28327, 28356)}, {'upstream': (143250359, 143250388, 29510, 29539), 'downstream': (143250398, 143250427, 29549, 29578)}, {'upstream': (143250463, 143250492, 29614, 29643), 'downstream': (143250502, 143250531, 30664, 30693)}, {'upstream': (143250539, 143250568, 30701, 30731), 'downstream': (143250578, 143250607, 32960, 32989)}, {'upstream': None, 'downstream': (143251217, 143251246, 31403, 31432)}, {'upstream': (143251303, 143251332, 28289, 28318), 'downstream': (143251341, 143251370, 29159, 29188)}, {'upstream': (143251329, 143251357, 28214, 28242), 'downstream': (143251367, 143251396, 28717, 28746)}, {'upstream': (143251404, 143251433, 28289, 28318), 'downstream': (143251442, 143251471, 28327, 28356)}, {'upstream': (143251692, 143251721, 29510, 29539), 'downstream': (143251731, 143251760, 29549, 29578)}, {'upstream': (143251718, 143251747, 29536, 29565), 'downstream': (143251757, 143251786, 30586, 30615)}, {'upstream': (143251796, 143251825, 29614, 29643), 'downstream': (143251837, 143251864, 30666, 30693)}, {'upstream': (143251874, 143251903, 32921, 32950), 'downstream': (143251913, 143251942, 32960, 32989)}, {'upstream': None, 'downstream': (143252575, 143252604, 39412, 39441)}, {'upstream': None, 'downstream': (143252603, 143252632, 28179, 28208)}, {'upstream': (143252612, 143252641, 28188, 28217), 'downstream': (143252650, 143252679, 31478, 31507)}, {'upstream': (143252713, 143252742, 28289, 28318), 'downstream': (143252751, 143252780, 28327, 28356)}, {'upstream': (143253106, 143253134, 28214, 28242), 'downstream': (143253144, 143253173, 28717, 28746)}, {'upstream': (143253495, 143253524, 29536, 29565), 'downstream': (143253534, 143253563, 30586, 30615)}, {'upstream': None, 'downstream': (143253663, 143253692, 32932, 32960)}, {'upstream': None, 'downstream': (143254205, 143254233, 28666, 28694)}, {'upstream': (143254192, 143254221, 29121, 29150), 'downstream': (143254230, 143254259, 29159, 29188)}, {'upstream': (143254293, 143254322, 28289, 28318), 'downstream': (143254331, 143254360, 28327, 28356)}, {'upstream': (143254585, 143254614, 31163, 31192), 'downstream': (143254624, 143254653, 36970, 36999)}, {'upstream': (143254634, 143254663, 31212, 31241), 'downstream': (143254673, 143254702, 31251, 31280)}, {'upstream': (143254787, 143254816, 28113, 28142), 'downstream': (143254825, 143254854, 31403, 31432)}, {'upstream': (143254813, 143254842, 31391, 31420), 'downstream': (143254851, 143254880, 31429, 31458)}, {'upstream': (143254862, 143254891, 28188, 28217), 'downstream': (143254900, 143254929, 31478, 31507)}, {'upstream': (143254888, 143254916, 28214, 28242), 'downstream': (143254926, 143254955, 28717, 28746)}, {'upstream': None, 'downstream': (143255001, 143255030, 28327, 28356)}, {'upstream': None, 'downstream': (143255281, 143255310, 29549, 29578)}, {'upstream': (143255268, 143255297, 29536, 29565), 'downstream': (143255307, 143255336, 30586, 30615)}, {'upstream': (143255424, 143255453, 32921, 32950), 'downstream': (143255463, 143255492, 32960, 32990)}, {'upstream': (143255883, 143255912, 31163, 31192), 'downstream': (143255922, 143255951, 36970, 36999)}, {'upstream': (143255932, 143255961, 31212, 31241), 'downstream': (143255971, 143256000, 31251, 31280)}, {'upstream': (143249370, 143249370, 28061, 28061), 'downstream': (143256039, 143256068, 28067, 28096)}, {'upstream': (143256085, 143256114, 28113, 28142), 'downstream': (143256123, 143256152, 39412, 39441)}, {'upstream': None, 'downstream': (143256151, 143256180, 28179, 28208)}, {'upstream': None, 'downstream': (143256224, 143256253, 28717, 28746)}, {'upstream': (143256261, 143256290, 28289, 28318), 'downstream': (143256299, 143256328, 28327, 28356)}, {'upstream': (143256549, 143256577, 29510, 29538), 'downstream': (143256588, 143256617, 29549, 29578)}, {'upstream': (143256575, 143256604, 29536, 29565), 'downstream': (143256614, 143256643, 30586, 30615)}, {'upstream': (143256731, 143256759, 30703, 30731), 'downstream': (143256770, 143256799, 32960, 32989)}, {'upstream': None, 'downstream': (143257230, 143257259, 36970, 36999)}, {'upstream': (143257306, 143257335, 31278, 31307), 'downstream': (143257347, 143257376, 31319, 31348)}, {'upstream': (143257518, 143257547, 28289, 28318), 'downstream': (143257556, 143257585, 29159, 29188)}, {'upstream': (143257544, 143257572, 28214, 28242), 'downstream': (143257582, 143257611, 28717, 28746)}, {'upstream': (143257907, 143257936, 28577, 28606), 'downstream': (143257946, 143257975, 29523, 29552)}, {'upstream': None, 'downstream': (143258417, 143258446, 29549, 29578)}, {'upstream': (143258404, 143258433, 28603, 28632), 'downstream': (143258443, 143258472, 29107, 29136)}, {'upstream': (143258482, 143258511, 29614, 29643), 'downstream': (143258521, 143258550, 30664, 30693)}, {'upstream': (143258560, 143258589, 32921, 32950), 'downstream': (143258599, 143258628, 32960, 32989)}, {'upstream': (143259288, 143259317, 28188, 28217), 'downstream': (143259326, 143259355, 28327, 28356)}, {'upstream': (143259576, 143259605, 28577, 28606), 'downstream': (143259615, 143259644, 29549, 29578)}], 'flags': [None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None]}



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


def process_flanks(softclip_loci, read_ref_start):
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

    print('After flagging misaligned loci')
    for i, coords in enumerate(softclip_loci['coords']):
        print(f'  {softclip_loci["loci"][i]}: {coords}\t{softclip_loci["flags"][i]} ')

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
            print(i, j, softclip_loci['flags'][i], softclip_loci['flags'][j], current_locus, next_locus)

            if softclip_loci['flags'][i] is not None:
                print('A')
                j = i + 1
                break
            
            if softclip_loci['flags'][j] is not None:
                print('B')
                j += 1
                continue

            if current_locus[1] < read_ref_start and next_locus[1] > read_ref_start:
                print('C')
                dir_switch = True
                continue
            
            print(i, j, current_coords, next_coords, sep='\t')
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

merge_coords = process_flanks(loci, 143228447)

# for i, locus in enumerate(loci['loci']):
#     print(locus, loci['coords'][i], loci['flags'][i])

print(merge_coords)
