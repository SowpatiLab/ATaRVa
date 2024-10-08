from snp_utils import haplocluster_reads
from vcf_writer import *
import numpy as np
import sys
import statistics


def length_genotyper(hallele_counter, global_loci_info, global_loci_variations, locus_key, read_indices, contig, locus_start, locus_end, ref, out):
    stringency_factor = [0.5, 1, 2]

    motif_len = float(global_loci_info[locus_key][4])
    if motif_len <= 6:
        factor = stringency_factor[0]
    elif motif_len <= 15:
        factor = stringency_factor[1]
    else:
        factor = stringency_factor[2]
    
    alleles = sorted(hallele_counter.keys())
    diffs = np.diff(alleles)
    mean_diff = np.mean(diffs)
    std_diff = np.std(diffs)
    
    lower_bound = mean_diff - std_diff*factor
    upper_bound = mean_diff + std_diff*factor

    clusters = []
    current_cluster = []
    for idx,num in enumerate(alleles):
        if idx==0:
            current_cluster.append(num)
        else:
            diff = num - alleles[idx - 1]

            if lower_bound <= diff <= upper_bound:
                current_cluster.append(num)
            else:
                # Start a new cluster
                clusters.append(current_cluster)
                current_cluster = [num]
    clusters.append(current_cluster)
    
    del_cl = []
    for idx,each_cluster in enumerate(clusters):
        cl_reads = 0
        for each_num in each_cluster:
            cl_reads += hallele_counter[each_num]
        if cl_reads < len(read_indices)*0.2:
            del_cl.append(each_cluster)
    for less_cl in del_cl:
        clusters.remove(less_cl)

    while len(clusters) > 2:
        cluster_reads = []
        for each_cluster in clusters:
            cl_reads = 0
            for each_num in each_cluster:
                cl_reads += hallele_counter[each_num]
            cluster_reads.append(cl_reads)
        del_cl = clusters[cluster_reads.index(min(cluster_reads))]
        clusters.remove(del_cl)

    def genotyper(group, hallele_counter):
        allele_freq = []
        for allele in group:
            allele_freq.append(hallele_counter[allele])
        fidx = allele_freq.index(max(allele_freq))
        final_allele = group[fidx]
        return (final_allele, hallele_counter[final_allele], sum(allele_freq))

    genotypes = []
    allele_count = {}
    phased_read = []
    chosen_snpQ = '.'
    snp_num = '.'
    if len(clusters) == 1:
        homo_allele = genotyper(clusters[0], hallele_counter)
        vcf_homozygous_writer(ref, contig, locus_key, global_loci_info, homo_allele[0], global_loci_variations, homo_allele[1], out)

    elif len(clusters) == 2:
        hetero_alleles = []
        for each_cluster in clusters:
            allele = genotyper(each_cluster, hallele_counter)
            hetero_alleles.append(allele)
            genotypes.append(allele[0])
            phased_read.append(allele[2])
            if allele[0] not in allele_count:
                allele_count[allele[0]] = allele[1]
            else:
                allele_count[str(allele[0])] = allele[1]
        vcf_heterozygous_writer(contig, genotypes, locus_start, locus_end, allele_count, len(read_indices), global_loci_info, ref, out, chosen_snpQ, phased_read, snp_num)

    return [True, 10]
        
            


    

def analyse_genotype(contig, locus_key, global_loci_info,
                     global_loci_variations, global_read_variations, global_snp_positions, hallele_counter,
                     ref, out, sorted_global_snp_list, level_split, snpQ, snpC, snpD, snpR, phasingR):

    locus_start = int(global_loci_info[locus_key][1])
    locus_end = int(global_loci_info[locus_key][2])
    state = False


    read_indices = global_loci_variations[locus_key]['reads']


    snp_positions = set()
    for rindex in read_indices:
        snp_positions |= (global_read_variations[rindex]['snps'])

    snp_positions = sorted(list(filter(lambda x: (x in global_snp_positions) and (global_snp_positions[x]['cov'] >= 3) and
                                                    (locus_start - snpD < x < locus_end + snpD),
                            snp_positions)))


    snp_allelereads = {}
    read_indices = set(read_indices)
    non_ref_snp_cov = {}
    for pos in snp_positions:
        c_point=0
        coverage = set()
        non_ref_nucs = [nucleotides for nucleotides in global_snp_positions[pos] if nucleotides not in ['cov', 'Qval', 'r']]
        for each_nuc in non_ref_nucs:
            reads_of_nuc = global_snp_positions[pos][each_nuc].intersection(read_indices)
            if len(reads_of_nuc) == 0: continue
            coverage.add(len(reads_of_nuc))

            if (sum([global_snp_positions[pos]['Qval'][read_idx] for read_idx in reads_of_nuc])/len(reads_of_nuc)) <= 13:
                c_point=1
                break
        if (len(coverage)==0) or (c_point==1): continue
        else: non_ref_snp_cov[pos] = max(coverage)
            
        snp_allelereads[pos] = { 'cov': 0, 'reads': set(), 'alleles': {}, 'Qval': {} }
        for nuc in global_snp_positions[pos]:
            if (nuc == 'cov') or (nuc == 'Qval'): continue
            snp_allelereads[pos]['alleles'][nuc] = global_snp_positions[pos][nuc].intersection(read_indices)
            snp_allelereads[pos]['cov'] += len(snp_allelereads[pos]['alleles'][nuc])
            if nuc!='r':
                snp_allelereads[pos]['Qval'].update(dict([(read_idx,global_snp_positions[pos]['Qval'][read_idx]) for read_idx in snp_allelereads[pos]['alleles'][nuc]]))

    del_positions = list(filter(lambda x: snp_allelereads[x]['cov'] < 5, snp_allelereads.keys()))
    # snp_positions = list(filter(lambda x: snp_allelereads[x]['cov'] >= 5, snp_allelereads.keys()))
    for pos in del_positions:
        del snp_allelereads[pos]


    ordered_snp_on_cov = sorted(snp_allelereads.keys(), key = lambda item : non_ref_snp_cov[item], reverse = True)


    haplotypes, min_snp, skip_point, chosen_snpQ, phased_read, snp_num = haplocluster_reads(snp_allelereads, ordered_snp_on_cov, read_indices, level_split, snpQ, snpC, snpR, phasingR) # SNP ifo and supporting reads for specific locus are given to the phasing function

    if haplotypes == (): # if the loci has no significant snps
        state, skip_point = length_genotyper(hallele_counter, global_loci_info, global_loci_variations, locus_key, read_indices, contig, locus_start, locus_end, ref, out)
        return [state, skip_point]
    
    if min_snp != -1:
        min_idx = sorted_global_snp_list.index(min_snp)
        del sorted_global_snp_list[:min_idx]
        del_snps = set()
        for pos in global_snp_positions:
            if pos < min_snp: del_snps.add(pos)
        for pos in del_snps:
            del global_snp_positions[pos]

    
    locus_read_allele = global_loci_variations[locus_key]['read_allele'] # extracting allele info from global_loci_variation
    hap_alleles = ([], [])
    for idx, haplo_tuple in enumerate(haplotypes): # Getting the mode value of alleles from clusters
        for each_read in haplo_tuple:
            hap_alleles[idx].append(locus_read_allele[each_read][0])


    genotypes = []
    for alleles_set in hap_alleles: # making the set of final alleles from two clusters
        genotypes.append(statistics.mode(alleles_set))

    ref_len = int(locus_key[locus_key.index('-')+1:]) - int(locus_key[locus_key.index(':')+1 : locus_key.index('-')])
    
        

    allele_count = {}
    for index, allele in enumerate(genotypes):
        if allele not in allele_count:
            allele_count[allele] = hap_alleles[index].count(allele)
        else:
            allele_count[str(allele)] = hap_alleles[index].count(allele)

    vcf_heterozygous_writer(contig, genotypes, locus_start, locus_end, allele_count, len(read_indices), global_loci_info, ref, out, chosen_snpQ, phased_read, snp_num)
    state = True
    return [state, skip_point]
    