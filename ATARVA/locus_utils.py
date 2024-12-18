def count_alleles(locus_key, read_indices, global_loci_variations, allele_counter, hallele_counter):
    """
    Counts the read distribution for each allele length
    """
    for rindex in read_indices:
        halen, alen = global_loci_variations[locus_key]['read_allele'][rindex]

        try: allele_counter[alen] += 1
        except KeyError: allele_counter[alen] = 1

        try: hallele_counter[halen] += 1
        except KeyError: hallele_counter[halen] = 1


def record_snps(read_indices, old_reads, new_reads, global_read_variations, global_snp_positions, sorted_global_snp_list):

    for rindex in read_indices:
        if rindex not in new_reads: continue
        
        rstart = global_read_variations[rindex]['s']
        rend   = global_read_variations[rindex]['e']
        snps   = global_read_variations[rindex]['snps']
        dels   = global_read_variations[rindex]['dels']
        
        for pos in sorted_global_snp_list:
        
            if pos < rstart: continue
            if pos > rend: break
            if pos not in snps and pos not in dels:
                if 'r' in global_snp_positions[pos]: global_snp_positions[pos]['r'].add(rindex)
                else: global_snp_positions[pos]['r'] = {rindex}
                global_snp_positions[pos]['cov'] += 1


def process_locus(locus_key, global_loci_variations, global_read_variations, global_snp_positions, prev_reads, sorted_global_snp_list, maxR, minR):
    
    homozygous = False
    ambiguous = False
    homozygous_allele = 0
    reads_of_homozygous = []
    read_indices = global_loci_variations[locus_key]['reads']   # the read indices which cover the locus
    total_reads = len(read_indices)                             # total number of reads
    max_limit=0

    # remove if the locus has poor coverage
    if total_reads < minR:
        # coverage of the locus is low
        prev_reads = set(read_indices)
        return [prev_reads, homozygous, ambiguous, homozygous_allele, reads_of_homozygous, {}, 0, max_limit]
    elif total_reads > maxR:
        # coverage of the locus is high
        read_indices = read_indices[:maxR]
        max_limit=1
        # prev_reads = set(read_indices)
        # return [prev_reads, homozygous, ambiguous, homozygous_allele, reads_of_homozygous, {}, 1]
    
    current_reads = set(read_indices)
    old_reads = prev_reads - current_reads
    new_reads = current_reads - prev_reads
    
    # recording the counts of each allele length across all reads
    allele_counter = {};  hallele_counter = {}
    count_alleles(locus_key, read_indices, global_loci_variations, allele_counter, hallele_counter)
    if len(hallele_counter) == 1:
        homozygous = True
        homozygous_allele = list(hallele_counter.keys())[0]
        reads_of_homozygous = read_indices.copy()
    
    else:
        filtered_alleles = list(filter(lambda x: hallele_counter[x] > 1, hallele_counter.keys()))
        if len(filtered_alleles) == 1 and hallele_counter[filtered_alleles[0]]/total_reads >= 0.75:
            homozygous = True
            homozygous_allele = filtered_alleles[0]
            reads_of_homozygous = [rindex for rindex in global_loci_variations[locus_key]['read_allele'] if homozygous_allele == global_loci_variations[locus_key]['read_allele'][rindex][0]]
        else:
            ambiguous = True
            
    
    record_snps(read_indices, old_reads, new_reads, global_read_variations, global_snp_positions, sorted_global_snp_list)
    
    prev_reads = current_reads.copy()
    return [prev_reads, homozygous, ambiguous, homozygous_allele, reads_of_homozygous, hallele_counter, 10, max_limit]