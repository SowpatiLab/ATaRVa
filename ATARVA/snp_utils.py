# import sys

def haplocluster_reads(snp_allelereads, ordered_snp_on_cov, read_indices, snpQ, snpC, snpR, phasingR):

    threshold_range = [(0.3, 0.7),(0.25, 0.75),(0.2, 0.8)] # threshold values to get Significant_snps
    for idx,range in enumerate(threshold_range):

        final_haplos = ()
        min_snp = -1
        skip_point = 10

        r1 = range[0] # threshold value 1
        r2 = range[1] # threshold value 2

        filtered_significant_poses = {}
        ordered_split_pos = []
        for pos in ordered_snp_on_cov:
            if snp_allelereads[pos]['cov'] < (0.6*len(read_indices)): break
            if sum(r1 * sum(len(snp_allelereads[pos]['alleles'][nucs]) for nucs in snp_allelereads[pos]['alleles']) <= len(snp_allelereads[pos]['alleles'][nucs]) <= r2 * sum(len(snp_allelereads[pos]['alleles'][nucs]) for nucs in snp_allelereads[pos]['alleles']) for nucs in snp_allelereads[pos]['alleles']) >= 2:
                ordered_split_pos.append(pos)
                filtered_significant_poses[pos] = {'cov': snp_allelereads[pos]['cov'], 'alleles': snp_allelereads[pos]['alleles'], 'Qval': snp_allelereads[pos]['Qval']}


        if len(filtered_significant_poses) == 0:
            if idx < 2: #level_split:
                continue
            else:
                skip_point = 0
                return [final_haplos, min_snp, skip_point, '', '', 0]

        final_haplos, status, min_snp, skip_point, chosen_snpQ, phased_read, snp_num = merge_snpreadsets(filtered_significant_poses, ordered_split_pos, read_indices, snpQ, snpC, snpR, phasingR) # calling the phasing function
        if status: break
        if idx == 2: #level_split:
            break

    return [final_haplos, min_snp, skip_point, chosen_snpQ, phased_read, snp_num]


def merge_snpreadsets(Significant_poses, ordered_split_pos, read_indices, snpQ, snpC, snpR, phasingR):


    skip_point = 10

    sorted_alt_snp = ordered_split_pos[:snpC]
    Alt_snp_Qval = []
    for snps in sorted_alt_snp:
        total_reads = Significant_poses[snps]['cov']
        qvalue_list = set()
        for nucs in Significant_poses[snps]['alleles']:
            if (nucs!='r') and ((len(Significant_poses[snps]['alleles'][nucs])/total_reads) >= 0.2):
                qvalue_list.add(sum([Significant_poses[snps]['Qval'][r_idx] for r_idx in Significant_poses[snps]['alleles'][nucs]]) / len(Significant_poses[snps]['alleles'][nucs]))

        if qvalue_list != set():
            Alt_snp_Qval.append(str(int(max(qvalue_list))))
        
        
    chosen_snpQ = ','.join(Alt_snp_Qval)
    phased_read = ''
    snp_num = len(sorted_alt_snp)
    
    if len(sorted_alt_snp) == 0:
        skip_point = 5
        return [(), False, -1, skip_point, chosen_snpQ, phased_read, snp_num]
    

    sorted_filtered_dict = dict((snps,Significant_poses[snps]['alleles']) for snps in sorted_alt_snp)


    #to remove the nucs with lower read contribution in sig_snps
    for pos in sorted_filtered_dict:
        del_nucs = []
        tot_reads = sum([len(vals) for vals in sorted_filtered_dict[pos].values()]) # calculate the total read for that snp
        if len(sorted_filtered_dict[pos]) == 2:
            for nucs in sorted_filtered_dict[pos]:
                if len(sorted_filtered_dict[pos][nucs])/tot_reads < snpR: del_nucs.append(nucs) # if the reads in 'nucs' have reads less than 25%, delete it
        else:  # if snp has more than 2 'nucs' in it  # SUS!!!
            for nucs in sorted_filtered_dict[pos]:
                if len(sorted_filtered_dict[pos][nucs])/tot_reads < snpR: del_nucs.append(nucs) # if the reads in 'nucs' have reads less than 25%, delete it
        for nucs in del_nucs:
            del sorted_filtered_dict[pos][nucs]

    # to remove pos with less than 2 nucs
    del_pos = []
    for pos in sorted_filtered_dict:
        if len(sorted_filtered_dict[pos]) == 2:
            pass
        elif len(sorted_filtered_dict[pos]) < 2:
            del_pos.append(pos)
    for pos in del_pos:
        del sorted_filtered_dict[pos]

    if len(sorted_filtered_dict) == 0: # after removing snp based on their nuc's read coverage, if there are no snp left then go to next threshold range
        skip_point = 1
        return [(), False, -1, skip_point, chosen_snpQ, phased_read, snp_num]


    # SNP Combination method 
    poses = list(sorted_filtered_dict.keys()) # get only SNP coordinate
    pos_cluster = {}
    for idx in range(len(poses)): # Comparison starting from 1st snp to all other snp
        if (pos_cluster!={}) and idx == len(poses)-1: break 
        pos_cluster[poses[idx]] = {}
        for pos in poses[idx + 1:]:
            current_pos_values = list(sorted_filtered_dict[poses[idx]].values())  # [  {1,2,3,4,5},  {6,7,8,9,10} ]  nuc's read set for current Target SNP
            mis_score = 0
            for reads_set in sorted_filtered_dict[pos].values():  # [  {1,2,3,4,6},  {5,7,8,9,10}  ]  nuc's read set for Query SNP
                similar_reads = set()
                for i in range(len(current_pos_values)):
                    intersection = reads_set & current_pos_values[i]  # intersection is calculated for each combination of 4 read sets
                    similar_reads.add(len(intersection))   # taking the min value from a pairwise interection values
                mis_score += min(similar_reads)
            pos_cluster[poses[idx]][pos] = mis_score # pos_cluster = { 1023 : {1036 : 0, 1045 : 1, 1123 : 3,..},  1036 : { 1045: 1, 1123 : 2,..}, ..............}

    significant_snps = []
    for each_pos in pos_cluster:
        if list(pos_cluster[each_pos].values()).count(0) >= 2: # check whether any of the pos_cluster have atleast 2 zeros in their mismatch scores
            significant_snps.append(each_pos)
            significant_snps.extend(sorted(pos_cluster[each_pos].keys(), key=lambda item: pos_cluster[each_pos][item])) # if yes take that pos_cluster and proceed for clustering
            break
    if significant_snps == []: # if there are no pos_cluster with atleast 2 zeros in it
        least_mismatches = {}
        for each_pos in pos_cluster:
            least_mismatches[each_pos] = sum(sorted(val for val in pos_cluster[each_pos].values())[:2]) # take sum of 1st 2 mismatch scores from sorted pos:mis_score and append in to 'least_mismatches' dictionary
        current_pos = sorted(least_mismatches.keys(), key=lambda item: least_mismatches[item])[0] # now take the snp position with least score and proceed for clustering
        significant_snps.append(current_pos)
        significant_snps.extend(sorted(pos_cluster[current_pos].keys(), key=lambda item: pos_cluster[current_pos][item])) # [1023, 1036, 1045, 1123]



    final_ordered_dict = {k: sorted_filtered_dict[k] for k in significant_snps if k in sorted_filtered_dict} # dict of 'final snp & its nucs : reads' for haplotyping, generated from 'significant_snps' list above
    
    min_snp = min(list(final_ordered_dict.keys()))

    cluster1, cluster2 = set(), set()
    for position_keys in final_ordered_dict: # clustering starts from the position of how its arranged in the dict, The 2 read_set will be two new haplotypes
        for nuc_read in final_ordered_dict[position_keys].values():
            if not cluster1:
                cluster1 |= nuc_read
            elif not cluster2:
                cluster2 |= nuc_read
            elif (len(nuc_read & cluster1)>0.7*len(nuc_read)) and  (len(nuc_read & cluster2)<0.05*len(nuc_read)): # then successive read_set are joined based on their intersection percentage
                cluster1 |= nuc_read
            elif (len(nuc_read & cluster2)>0.7*len(nuc_read)) and  (len(nuc_read & cluster1)<0.05*len(nuc_read)):
                cluster2 |= nuc_read

    if ((len(cluster1)+len(cluster2)) >= phasingR*len(read_indices)): # after clustering the total reads in the both reads should be greater thn 50% of total supporting reads for that specific locus
        phased_read = [len(cluster1),len(cluster2)]
        return [(cluster1, cluster2), True, min_snp, skip_point, chosen_snpQ, phased_read, snp_num]
    else:
        skip_point = 2
        return [(), False, -1, skip_point, chosen_snpQ, phased_read, snp_num]