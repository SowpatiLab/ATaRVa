from snp_utils import haplocluster_reads
from vcf_writer import *
import numpy as np
import sys
import statistics


def length_genotyper(hallele_counter, global_loci_info, global_loci_variations, locus_key, read_indices, contig, locus_start, locus_end, ref, out):
    
    total_reads = len(read_indices)
    filtered_alleles = list(filter(lambda x: hallele_counter[x] > 1, hallele_counter.keys()))
    top_alleles = [al for al in filtered_alleles if (hallele_counter[al]/total_reads) > 0.2]
    locus_read_allele = global_loci_variations[locus_key]['read_allele'] # extracting allele info from global_loci_variation

    if len(top_alleles) == 1:
        hap_reads = [read_id for read_id in locus_read_allele if locus_read_allele[read_id][0]==top_alleles[0]]
        vcf_homozygous_writer(ref, contig, locus_key, global_loci_info, top_alleles[0], global_loci_variations, hallele_counter[top_alleles[0]], out, hap_reads)

    elif len(top_alleles)>1:
        phased_read = ['.','.']
        chosen_snpQ = '.'
        snp_num = '.'
        hetero_alleles = sorted(top_alleles, key=lambda x: hallele_counter[x], reverse=True)[:2]

        hap_reads = ([],[])
        for i,al in enumerate(hetero_alleles):
            hap_reads[i].extend([read_id for read_id in locus_read_allele if locus_read_allele[read_id][0]==al])

        allele_count = {}
        for al in hetero_alleles:
            allele_count[al] = hallele_counter[al]
        vcf_heterozygous_writer(contig, hetero_alleles, locus_start, global_loci_variations, locus_end, allele_count, len(read_indices), global_loci_info, ref, out, chosen_snpQ, phased_read, snp_num, hap_reads)
    else:
        pass # write allele distribution with only one read supporting to it in vcf

    return [True, 10]
        
            


    

def analyse_genotype(contig, locus_key, global_loci_info,
                     global_loci_variations, global_read_variations, global_snp_positions, hallele_counter,
                     ref, out, sorted_global_snp_list, level_split, snpQ, snpC, snpD, snpR, phasingR, maxR, max_limit):

    locus_start = int(global_loci_info[locus_key][1])
    locus_end = int(global_loci_info[locus_key][2])
    state = False

    if max_limit == 0:
        read_indices = global_loci_variations[locus_key]['reads']
    else:
        read_indices = global_loci_variations[locus_key]['reads'][:maxR]


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
    
    hap_reads = ([],[])
    for i,al in enumerate(genotypes):
        hap_reads[i].extend([read_id for read_id in haplotypes[i] if locus_read_allele[read_id][0]==al])

    allele_count = {}
    for index, allele in enumerate(genotypes):
        if allele not in allele_count:
            allele_count[allele] = hap_alleles[index].count(allele)
        else:
            allele_count[str(allele)] = hap_alleles[index].count(allele)

    vcf_heterozygous_writer(contig, genotypes, locus_start, global_loci_variations, locus_end, allele_count, len(read_indices), global_loci_info, ref, out, chosen_snpQ, phased_read, snp_num, hap_reads)
    state = True
    return [state, skip_point]
    