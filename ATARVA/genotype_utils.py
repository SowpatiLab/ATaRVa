from ATARVA.snp_utils import haplocluster_reads
from ATARVA.vcf_writer import *
import numpy as np
import statistics
from sklearn.cluster import KMeans
import warnings
from threadpoolctl import threadpool_limits


def length_genotyper(hallele_counter, global_loci_info, global_loci_variations, locus_key, read_indices, contig, locus_start, locus_end, ref, out, male, log_bool, decomp):

    read_indices = sorted(read_indices)
    locus_read_allele = global_loci_variations[locus_key]['read_allele']
    
    alen_with_1read = [item[0] for item in hallele_counter.items() if item[1]==1] # allele with 1 read contribution
    alen_with_gread = set(hallele_counter.keys()) - set(alen_with_1read) # allele with more than 1 read contribution
    main_read_id = []
    alen_data = []
    
    for id in read_indices:
        if locus_read_allele[id][0] in alen_with_1read: # checking if the '1 read - allele' is nearby any of other 'good read - allele'
            num = locus_read_allele[id][0]
            for i in alen_with_gread:
                if i in range(num-10, num+10): # '1 read - allele' is considered if other allele are within 10 bp on either of the side
                    alen_data.append(num)
                    main_read_id.append(id)
                    break
        else:
            alen_data.append(locus_read_allele[id][0])
            main_read_id.append(id)

    if alen_data == []:
        return [False, 6]

    data = np.array(alen_data)
    data = data.reshape(-1, 1)
    with threadpool_limits(limits=1):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=UserWarning)
            kmeans = KMeans(n_clusters=2, init='k-means++', n_init=5, random_state=0).fit(data)
    cluster_labels = kmeans.labels_  
    c1 = [i for i, x in enumerate(cluster_labels) if x == 0]
    c2 = [i for i, x in enumerate(cluster_labels) if x == 1]

    haplotypes = ([main_read_id[idx] for idx in c1], [main_read_id[idx] for idx in c2])
    cutoff = 0.15*len(alen_data) # 15%

    if male:
        cluster_len = [len(c1), len(c2)]
        cidx = cluster_len.index(max( cluster_len ))
        if cluster_len[cidx]>=cutoff:
            mac = haplotypes[cidx]
            hap_al1 = [alen_data[i] for i in [c1,c2][cidx]]
            genotypes = statistics.mode(hap_al1)
            hap_reads = [read_id for read_id in mac if locus_read_allele[read_id][0]==genotypes]
            vcf_homozygous_writer(ref, contig, locus_key, global_loci_info, genotypes, global_loci_variations, hap_al1.count(genotypes), out, hap_reads, log_bool, 'kmeans', decomp)
    
    elif (c1!=[] and len(c1)>=cutoff) and (c2!=[] and len(c2)>=cutoff):
        phased_read = ['.','.']
        chosen_snpQ = '.'
        snp_num = '.'
        hap_al1 = [alen_data[i] for i in c1]
        hap_al2 = [alen_data[i] for i in c2]
        haplo_alleles = (hap_al1,hap_al2)
        

        genotypes = []
        for alleles_set in haplo_alleles: # making the set of final alleles from two clusters
            genotypes.append(statistics.mode(alleles_set))

        hap_reads = ([],[])
        for i,al in enumerate(genotypes):
            hap_reads[i].extend([read_id for read_id in haplotypes[i] if locus_read_allele[read_id][0]==al])

        allele_count = {}
        for al in genotypes:
            allele_count[al] = hallele_counter[al]
        vcf_heterozygous_writer(contig, genotypes, locus_start, global_loci_variations, locus_end, allele_count, len(read_indices), global_loci_info, ref, out, chosen_snpQ, phased_read, snp_num, hap_reads, log_bool, 'kmeans', decomp)

    elif c1!=[] and len(c1)>=cutoff:
        hap_al1 = [alen_data[i] for i in c1]
        genotypes = statistics.mode(hap_al1)
        hap_reads = [read_id for read_id in haplotypes[0] if locus_read_allele[read_id][0]==genotypes]
        vcf_homozygous_writer(ref, contig, locus_key, global_loci_info, genotypes, global_loci_variations, hap_al1.count(genotypes), out, hap_reads, log_bool, 'kmeans', decomp)
        

    elif c2!=[] and len(c2)>=cutoff:
        hap_al2 = [alen_data[i] for i in c2]
        genotypes = statistics.mode(hap_al2)
        hap_reads = [read_id for read_id in haplotypes[1] if locus_read_allele[read_id][0]==genotypes]
        vcf_homozygous_writer(ref, contig, locus_key, global_loci_info, genotypes, global_loci_variations, hap_al2.count(genotypes), out, hap_reads, log_bool, 'kmeans', decomp)
        
    else:
        return [False, 6] # write allele distribution with only one read supporting to it in vcf
    
    return [True, 10]


def analyse_genotype(contig, locus_key, global_loci_info,
                     global_loci_variations, global_read_variations, global_snp_positions, hallele_counter,
                     ref, out, sorted_global_snp_list, snpQ, snpC, snpD, snpR, phasingR, maxR, max_limit, male, log_bool, decomp):
            
    locus_start = int(global_loci_info[locus_key][1])
    locus_end = int(global_loci_info[locus_key][2])

    state = False

    if max_limit == 0:
        read_indices = global_loci_variations[locus_key]['reads']
    else:
        read_indices = global_loci_variations[locus_key]['reads'][:maxR]

    if male: 
        state, skip_point = length_genotyper(hallele_counter, global_loci_info, global_loci_variations, locus_key, read_indices, contig, locus_start, locus_end, ref, out, male, log_bool, decomp)
        return [state, skip_point]


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


    haplotypes, min_snp, skip_point, chosen_snpQ, phased_read, snp_num = haplocluster_reads(snp_allelereads, ordered_snp_on_cov, read_indices, snpQ, snpC, snpR, phasingR) # SNP ifo and supporting reads for specific locus are given to the phasing function

    if haplotypes == (): # if the loci has no significant snps
        state, skip_point = length_genotyper(hallele_counter, global_loci_info, global_loci_variations, locus_key, read_indices, contig, locus_start, locus_end, ref, out, male, log_bool, decomp)
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

    
    hap_reads = ([],[])
    for i,al in enumerate(genotypes):
        hap_reads[i].extend([read_id for read_id in haplotypes[i] if locus_read_allele[read_id][0]==al])

    allele_count = {}
    for index, allele in enumerate(genotypes):
        if allele not in allele_count:
            allele_count[allele] = hap_alleles[index].count(allele)
        else:
            allele_count[str(allele)] = hap_alleles[index].count(allele)

    vcf_heterozygous_writer(contig, genotypes, locus_start, global_loci_variations, locus_end, allele_count, len(read_indices), global_loci_info, ref, out, chosen_snpQ, phased_read, snp_num, hap_reads, log_bool, 'SNP', decomp)
    state = True
    return [state, skip_point]
    