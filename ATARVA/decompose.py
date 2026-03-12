from collections import Counter
from bitarray import bitarray
import regex as re
import warnings

divisor_dict = {2:[1], 3:[1], 4:[1,2], 5:[1], 6:[1,2,3], 7:[1], 8:[1,2,4], 9:[1,3], 10:[1,2,5]}

cyclic_motif_registry = {}

def is_valid_repeat_block(s):

    if not s.startswith('('):
        return False
    if ')' not in s:
        return False
    close = s.find(')')
    if close <= 1:
        return False
    return s[close+1:].isdigit()


def get_repeat_components(s):

    close = s.index(')')
    motif = s[1:close]
    count = int(s[close+1:])
    return motif, count

def get_cyclic_variants(motif):
    n = len(motif)
    return [motif[i:] + motif[:i] for i in range(n)]

def get_canonical_motif(motif):
    return min(get_cyclic_variants(motif))

def register_and_get_motif(motif):
    if not motif or len(motif) == 0:
        return motif
    
    canonical = get_canonical_motif(motif)
    
    if canonical not in cyclic_motif_registry:
        cyclic_motif_registry[canonical] = motif
        return motif
    else:
        return cyclic_motif_registry[canonical]


def convert_to_bitset(seq):
    lbit = {'A': '0', 'C': '0', 'G': '1', 'T': '1', 'N': '1'}
    rbit = {'A': '0', 'C': '1', 'G': '0', 'T': '1', 'N': '1'}
    
    lbitseq = bitarray()
    rbitseq = bitarray()
    
    for s in seq:
        lbitseq.extend(lbit.get(s, '1'))
        rbitseq.extend(rbit.get(s, '1'))
    
    return lbitseq, rbitseq

def shift_and_match(seq):
    shift_list = []
    # best_shift = motif_length
    max_matches = 0

    shift_values = set(range(1, 11))

    for shift in sorted(shift_values):
        if shift < 1:
            continue

        lbitseq, rbitseq = convert_to_bitset(seq)

        lmatch = ~(lbitseq ^ (lbitseq >> shift))
        rmatch = ~(rbitseq ^ (rbitseq >> shift))
        match = lmatch & rmatch
        shift_list.append(match)

    return shift_list

def kmp_search_non_overlapping(text, pattern):
    def compute_lps(pattern):
        lps = [0] * len(pattern)
        length = 0
        i = 1
        while i < len(pattern):
            if pattern[i] == pattern[length]:
                length += 1
                lps[i] = length
                i += 1
            else:
                if length != 0:
                    length = lps[length - 1]
                else:
                    lps[i] = 0
                    i += 1
        return lps

    lps = compute_lps(pattern)
    result = []
    i = 0
    j = 0  

    while i < len(text):
        if pattern[j] == text[i]:
            i += 1
            j += 1

        if j == len(pattern):  
            result.append(i - j)
            j = 0  
        elif i < len(text) and pattern[j] != text[i]:
            if j != 0:
                j = lps[j - 1]
            else:
                i += 1

    return result

def get_most_frequent_motif(sequence, motif_size, motif):
    cyc_motif = get_canonical_motif(motif) if motif else None
    if len(sequence) < motif_size:  
        return sequence, None  

    repeating_units = [sequence[i:i + motif_size] for i in range(len(sequence) - motif_size + 1)]
    motif_counts = Counter(repeating_units)
    
    
    if not motif_counts:  
        return sequence, None  

    most_common = []
    for each_motif,count in motif_counts.most_common():
        if len(most_common) >=2: break
        if get_canonical_motif(each_motif) == cyc_motif:
            if motif not in most_common: most_common.append(motif)
        else:
            most_common.append(each_motif)

    primary_motif = most_common[0] if most_common else sequence
    primary_motif = register_and_get_motif(primary_motif)
    secondary_motif = most_common[1] if len(most_common) > 1 else None

    return primary_motif, secondary_motif

def max_match(shift_list, gap_regions, motif_size):
    gap_wise_shift = []
    
    for start, end in gap_regions:
        gap_length = end - start
        
        candidate_shifts = []
        
        for idx, shift_pattern in enumerate(shift_list):
            if not shift_pattern.any():
                continue
                
            shift_value = idx + 1
            
            if shift_value * 2 > gap_length:  
                continue
                
            sub_pattern = shift_pattern[start:end]
            if len(sub_pattern) == 0:
                continue
            
            total_matches = sub_pattern.count()
            
            pattern_str = sub_pattern.to01()
            max_consecutive = max(len(run) for run in pattern_str.split('0')) if '1' in pattern_str else 0
            ideal_consecutive = max_consecutive + shift_value
            
            # Score based on multiple factors
            #score = (match_density * 0.4) + (max_consecutive / gap_length * 0.4) + (1.0 / shift_value * 0.2)
            
            candidate_shifts.append((shift_value, ideal_consecutive, total_matches))
        
        candidate_shifts.sort(key=lambda x: x[1], reverse=True)
        
        best_shift = -1
               
        for shift_value, ideal_consecutive, total_matches in candidate_shifts:
            #repeat_runs = max_consecutive / shift_value
            if ( ideal_consecutive / shift_value) >= 2.0:
                best_shift = shift_value
                break
        
        gap_wise_shift.append(best_shift)
    
    return gap_wise_shift
    
slide_threshold = {1:1, 2:2, 3:3, 4:4, 5:5, 6:6, 7:7, 8:8, 9:9, 10:10}
def window_scan(shift_list, motif_size, sequence, sequential_decomp, sequential_part, overall_boundary, current_gap):
    shift_seq = shift_list[motif_size-1]#[current_gap[0]:current_gap[1]]

    slide_size = slide_threshold.get(motif_size, 8)
    # slide_size = motif_size
    
    i=current_gap[0]
    start = i; end = current_gap[1]
    start_track = start
    initial = True
    while i<(current_gap[1]-slide_size + 1):
        words = shift_seq[i:i+slide_size]

        if (words.count()/len(words)) >= 0.9:
            if initial:
                calc_start = i-motif_size if i-motif_size >= start_track else start_track
                if calc_start>-1:
                    for b in overall_boundary:
                        if not (b[0] <= calc_start < b[1]):
                            pass
                        else:
                            start = b[1]
                            break
                    else:
                        start = calc_start
                else:
                    start = i #i-motif_size if i-motif_size > -1 else i
                initial = False
            end = i+motif_size
            i+=motif_size

            continue

        else:
            i+=1
            if not initial:
                calc_end = end+motif_size
                for b in overall_boundary:
                    if (start< calc_end <= b[0]) or (calc_end > start >= b[1]):
                        pass
                    else:
                        end = b[0]
                        break
                else:
                    end = calc_end    
                if start!=end:
                    start_track = decomposer([start, end], sequence, motif_size, sequential_decomp, sequential_part, shift_list, overall_boundary)
                initial = True
    if not initial:
        calc_end = end+motif_size
        for b in overall_boundary:
            if (start< calc_end <= b[0]) or (calc_end > start >= b[1]):
                pass
            else:
                end = b[0]
                break
        else:
            end = calc_end 
        if start!=end:
            start_track = decomposer([start, end], sequence, motif_size, sequential_decomp, sequential_part, shift_list, overall_boundary)


def shift_decomp(seq, motif_size, motif, boundary, state):
    decomposed_parts = []
    count = 1  
    
    primary_motif, _ = get_most_frequent_motif(seq, motif_size, motif)
    
    positions = kmp_search_non_overlapping(seq, primary_motif)
    
    if not positions:
        return [seq], boundary


    seq = seq[positions[0] : positions[-1]+motif_size]
    
    if state:
        b1 = boundary[0]
        boundary = [b1+positions[0], b1+positions[-1]+motif_size]

        
    for i in range(1, len(positions)):
        if positions[i] == positions[i - 1] + len(primary_motif):
            count += 1
        else:
            if count > 1:
                decomposed_parts.append(f"({primary_motif}){count}")
            else:
                decomposed_parts.append(primary_motif)
            interspersed = seq[positions[i - 1] + len(primary_motif):positions[i]]
            if interspersed:
                primary_motif, secondary_motif = get_most_frequent_motif(seq, motif_size, '')
                secondary_decomp, boundary = shift_decomp(interspersed, motif_size, '', boundary, False)
                if secondary_decomp:
                    decomposed_parts.extend(secondary_decomp)

            count = 1
            
    if count > 1:
        decomposed_parts.append(f"({primary_motif}){count}")
    else:
        decomposed_parts.append(primary_motif)
    last_motif_end = positions[-1] + len(primary_motif)
    leftover_sequence = seq[last_motif_end:]
    if leftover_sequence:
        decomposed_parts.append(leftover_sequence)
    return decomposed_parts, boundary

def decomposer(tuple_bound, sequence, best_shift, sequential_decomp, sequential_part, shift_list, overall_boundary):
    processed_seq, tuple_bound = shift_decomp(sequence[tuple_bound[0]:tuple_bound[1]], best_shift, '', tuple_bound, True)
    b1 = tuple_bound[0]; b2 = tuple_bound[1]
    sequential_decomp[b1] = processed_seq
    sequential_part[b1:b2] = 0
    for id in range(len(shift_list)):
        if shift_list[id] == 0: continue
        shift_list[id][b1:b2] = 0
    overall_boundary.append(tuple_bound)
    return b2

def gap_boundaries(overall_boundary, seq_len, chunk_state):
    overall_boundary = sorted(overall_boundary)

    if chunk_state and len(overall_boundary)>=1:
        if len(overall_boundary)>1:
            chunk_boundary = overall_boundary[0][0]
            chunk_boundary_end = overall_boundary[-1][1]
        else:
            chunk_boundary = overall_boundary[0][0]
            chunk_boundary_end = seq_len
    else:
        chunk_boundary = 0
        chunk_boundary_end = seq_len
        
    last_start = overall_boundary[-1][0]
    gap_regions = []
    for id,occupied_reg in enumerate(overall_boundary):
    
        if id == 0:
            current_end = occupied_reg[1]
            if occupied_reg[0] != chunk_boundary:
                gap_regions.append([chunk_boundary, occupied_reg[0]])
            elif len(overall_boundary) == 1:
                if occupied_reg[1] <= (seq_len-1):
                    gap_regions.append([occupied_reg[1],seq_len])
                    current_end = seq_len
            
    
        else:
            if occupied_reg[0] != current_end:
                gap_regions.append([current_end, occupied_reg[0]])
            current_end = occupied_reg[1]

            
    if current_end < chunk_boundary_end:
        gap_regions.append([current_end, chunk_boundary_end])

    return gap_regions

def motif_decomposition(sequence, motif_size):

    global cyclic_motif_registry
    cyclic_motif_registry = {}
    

    shift_list = shift_and_match(sequence)
    overall_boundary = []
    sequential_decomp = {}
    sequential_part = bitarray([1]*len(sequence))
    seq_len = len(sequence)
    gap_regions = [[0, seq_len]]
    rounds = 0
    while any(sequential_part):
        if rounds == 1:
            shift_list[0][:] = 0
            
        if (rounds==0) and (motif_size==1):
            gap_wise_shift = [1]
        else:
            gap_wise_shift = max_match(shift_list, gap_regions, motif_size)

        for index, best_shift in enumerate(gap_wise_shift):
            
            current_gap = gap_regions[index]
            
            if best_shift!=-1:
                previous_ovbound = overall_boundary.copy()
                window_scan(shift_list, best_shift, sequence, sequential_decomp, sequential_part, overall_boundary, current_gap)
                if previous_ovbound == overall_boundary:
 
                    c1 = current_gap[0]; c2 = current_gap[1]
                    covered_chunks = [i for i in overall_boundary if (i[0]>=c1 and i[1]<=c2) or ((i[1]==c1) or (i[0]==c2))]

                    if covered_chunks == []:
                        sequential_part[c1 : c2] = 0
                        sequential_decomp[c1] = [sequence[c1 : c2]]
                        overall_boundary.append([c1,c2])
                        continue

                    gap_in_chunks = gap_boundaries(covered_chunks, c2, True)
                    if not gap_in_chunks:
                        if (c1 == 0) or (c2 == seq_len):
                            gap_in_chunks = [current_gap]

                    for each_gap in gap_in_chunks:
                        sequential_part[each_gap[0] : each_gap[1]] = 0
                        for shift_id in range(len(shift_list)):
                            shift_list[shift_id][each_gap[0] : each_gap[1]] = 0
                        sequential_decomp[each_gap[0]] = [sequence[each_gap[0] : each_gap[1]]]

                        overall_boundary.append([each_gap[0] , each_gap[1]])

                
            else:
                if sequence[current_gap[0] : current_gap[1]]:
                    sequential_decomp[current_gap[0]] = [sequence[current_gap[0] : current_gap[1]] ]
                sequential_part[current_gap[0] : current_gap[1]] = 0
                overall_boundary.append([current_gap[0] , current_gap[1]])

        gap_regions = gap_boundaries(overall_boundary, seq_len, False)

        rounds += 1
        if rounds >= 20 :
            warnings.warn("Max rounds reached; decomposition may be incomplete")
        
    fseq = [i for k,v in sorted(sequential_decomp.items(), key = lambda x : x[0]) for i in v]

    if len(fseq)==1:
        if '(' in fseq[0]:
            non_rep_percent = 0
        else:
            non_rep_percent=1
    else:
        fseq, non_rep_percent = refine_decomposition(fseq, motif_size, len(sequence))

    return ["-".join(fseq), non_rep_percent]

def refine_decomposition(fseq, motif_size, seq_len):

    
    new_seq_list = []
    non_repeat = 0
    
    for i in fseq:
        if '(' in i and ')' in i:  # refining only the decomposed parts
            try:
                end_point = i.index(')')
                count = int(i[end_point+1:])
                tmp_motif = i[1:end_point]
                motif_len = len(tmp_motif)
                
                if len(set(tmp_motif)) == 1:
                    new_seq_list.append(f'({tmp_motif[0]}){motif_len*count}')
                
                elif (motif_len % 2 == 0):  
                    mid_point = int(motif_len/2)
                    if tmp_motif[0: mid_point] == tmp_motif[mid_point:]:
                        check_motif = tmp_motif[0: mid_point]
                        new_seq_list.append(f'({check_motif}){count*2}')
                    else:
                        new_seq_list.append(i)
                        
                elif motif_len <= motif_size:  
                    new_seq_list.append(i)
                        
                elif (motif_len != 1) and (motif_len > motif_size) and (motif_len % motif_size == 0):
                    for div in divisor_dict.get(motif_len, [1]):
                        check_motif = tmp_motif[:div]
                        e = 0
                        pattern = "(" + check_motif + "){e<=" + str(e) + "}"
                        matches = re.finditer(pattern, tmp_motif, overlapped=False)
                        tot_rep = sum(1 for _ in matches)
                        if tot_rep * len(check_motif) == motif_len:
                            motif_count = tot_rep * count
                            new_seq_list.append(f'({check_motif}){motif_count}')
                            break
                    else:
                        new_seq_list.append(i)
                else:
                    new_seq_list.append(i)
            except (ValueError, IndexError):
                new_seq_list.append(i)
        else:

            segment = i
            decomposed_segment = []
            
            pos = 0
            while pos < len(segment):
                found_repeat = False
                
                for test_len in range(motif_size, min(6, len(segment) - pos) + 1):
                    if test_len == 0:
                        continue
                    
                    test_motif = segment[pos:pos+test_len]
                    
                    count = 1
                    next_pos = pos + test_len
                    
                    while next_pos + test_len <= len(segment) and segment[next_pos:next_pos+test_len] == test_motif:
                        count += 1
                        next_pos += test_len
                    
                    if count >= 2:
                        decomposed_segment.append(f'({test_motif}){count}')
                        pos = next_pos
                        found_repeat = True
                        break
                
                if not found_repeat:

                    remaining_len = len(segment) - pos
                    if remaining_len <= motif_size:
                        decomposed_segment.append(segment[pos:])
                        break
                    else:
                        chunk_size = min(motif_size, remaining_len)
                        decomposed_segment.append(segment[pos:pos+chunk_size])
                        pos += chunk_size
            
            if len(decomposed_segment) > 1 or (len(decomposed_segment) == 1 and '(' in decomposed_segment[0]):
                new_seq_list.extend(decomposed_segment)
            else:
                new_seq_list.append(i)
                non_repeat += len(i)

    
    dlen = len(new_seq_list) - 1
    loc = 0
    tmp_loc = 1 # different value then 'loc' to start the loop
    refined_list = []
    
    while loc < dlen:
        if tmp_loc != loc:
            current = new_seq_list[loc]
        else:
            loc += 1
            if (loc == dlen) and current:
                refined_list.append(current)
                loc += 1
                break
        tmp_loc = loc
        next_item = new_seq_list[loc+1]

        current_state = 1 if '(' in current and ')' in current else 0
        next_state = 1 if '(' in next_item and ')' in next_item else 0
        
        if current_state == 1 and next_state == 1:
            try:
                e1 = current.index(')')
                e2 = next_item.index(')')
                tmp1 = current[1:e1]
                tmp2 = next_item[1:e2]
                
                if tmp1 == tmp2:  
                    combined_count = int(current[e1+1:]) + int(next_item[e2+1:])
                    refined_list.append(f'({tmp1}){combined_count}')
                    loc += 2
                else:
                    refined_list.append(current)
                    loc += 1
            except (ValueError, IndexError):
                refined_list.append(current)
                loc += 1
                
        elif current_state == 1 and next_state == 0:
            try:
                e1 = current.index(')')
                tmp1 = current[1:e1]
                
                if len(next_item) >= len(tmp1) and tmp1 == next_item[0:len(tmp1)]:
                    combined_count = int(current[e1+1:]) + 1
                    remaining = next_item[len(tmp1):]

                    while tmp1 == remaining[0:len(tmp1)]:
                        combined_count += 1 
                        remaining = remaining[len(tmp1):]
                        if (len(remaining)==1) and tmp1 == remaining:
                            combined_count += 1
                            remaining = ''
                            break
                    
                    refined_list.append(f'({tmp1}){combined_count}')
                    
                    if remaining:
                        current = remaining
                        loc += 1
                    else:
                        loc += 2
                        
                elif tmp1 == next_item:  
                    combined_count = int(current[e1+1:]) + 1
                    refined_list.append(f'({tmp1}){combined_count}')
                    loc += 2
                else:
                    refined_list.append(current)
                    loc += 1
                    
            except (ValueError, IndexError):
                refined_list.append(current)
                loc += 1
                
        elif current_state == 0 and next_state == 1:
            try:
                e2 = next_item.index(')')
                tmp2 = next_item[1:e2]
                
                if len(current) >= len(tmp2) and tmp2 == current[len(current)-len(tmp2):]:
                    combined_count = int(next_item[e2+1:]) + 1
                    
                    beginning = current[:len(current)-len(tmp2)]
                    if beginning:
                        refined_list.append(beginning)
                    
                    refined_list.append(f'({tmp2}){combined_count}')
                    loc += 2
                elif tmp2 == current:  
                    combined_count = int(next_item[e2+1:]) + 1
                    refined_list.append(f'({tmp2}){combined_count}')
                    loc += 2
                else:
                    refined_list.append(current)
                    loc += 1
            except (ValueError, IndexError):
                refined_list.append(current)
                loc += 1
                
        else:
            current = current + next_item

    
    if loc == dlen:
        refined_list.append(new_seq_list[loc])
    
    final_merged = []
    i = 0
    while i < len(refined_list):
        if i < len(refined_list) - 1:
            current = refined_list[i]
            next_item = refined_list[i+1]
            
            if '(' in current and ')' in current and '(' in next_item and ')' in next_item:
                try:
                    e1 = current.index(')')
                    e2 = next_item.index(')')
                    tmp1 = current[1:e1]
                    tmp2 = next_item[1:e2]
                    
                    if tmp1 == tmp2:
                        combined_count = int(current[e1+1:]) + int(next_item[e2+1:])
                        final_merged.append(f'({tmp1}){combined_count}')
                        i += 2
                        continue
                except (ValueError, IndexError):
                    pass
        
        final_merged.append(refined_list[i])
        i += 1
    
    non_repeat_len = 0
    for elem in final_merged:
        if '(' in elem and ')' in elem:
            try:
                end_point = elem.index(')')
                pass
            except ValueError:
                non_repeat_len += len(elem)
        else:
            non_repeat_len += len(elem)
    
    non_rep_percent = round((non_repeat_len / seq_len), 2) if seq_len > 0 else 0
    
    return final_merged, non_rep_percent