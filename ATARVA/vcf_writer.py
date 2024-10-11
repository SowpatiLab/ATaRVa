def vcf_writer(out):

    print(*['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE'], file=out, sep='\t')

def vcf_homozygous_writer(ref, contig, locus_key, global_loci_info, homozygous_allele, global_loci_variations, reads_of_homozygous, out):

    locus_start = int(global_loci_info[locus_key][1])
    locus_end = int(global_loci_info[locus_key][2])
    

    if type(reads_of_homozygous) == list:
        reads_of_homozygous = len(reads_of_homozygous)
    
    ref_allele_length = locus_end - locus_start
    DP = len(global_loci_variations[locus_key]['reads'])

    AC = 0; AN = 1; GT = '0/0'
    if homozygous_allele != ref_allele_length:
        AC = 1
        GT = '1/1'

    # INFO = 'AC=' + str(AC) + ';AN=' + str(AN) + ';DP=' + str(DP)+ ';END=' + str(locus_end)

    INFO = 'AC=' + str(AC) + ';AN=' + str(AN) + ';MOTIF=' + str(global_loci_info[locus_key][3]) + ';END=' + str(locus_end)

    # FORMAT = 'GT:AL:SD'
    # SAMPLE = str(GT) + ':' + str(homozygous_allele) + ',' + str(homozygous_allele) + ':' + str(reads_of_homozygous)

    FORMAT = 'GT:AL:SD:PC:DP:SN:SQ'
    SAMPLE = str(GT) + ':' + str(homozygous_allele) + ',' + str(homozygous_allele) + ':' + str(reads_of_homozygous) + ':.:' + str(DP) + ':.:.'

    print(*[contig, locus_start, '.',  ref.fetch(contig, locus_start-1, locus_end), '.', 0, 'PASS', INFO, FORMAT, SAMPLE], file=out, sep='\t')
    del global_loci_info[locus_key]


def vcf_heterozygous_writer(contig, genotypes, locus_start, locus_end, allele_count, DP, global_loci_info, ref, out, chosen_snpQ, phased_read, snp_num):

    locus_key = f'{contig}:{locus_start}-{locus_end}'
    final_allele = set(genotypes)
    heterozygous_allele = ''
    AC = 'AC'
    AN = 'AN'
    GT = 'GT'
    SD = 'SD'
    PC = 'PC'

    ref_allele_length = locus_end - locus_start

    if len(final_allele) == 1:
        AN = 1
        if ref_allele_length == tuple(final_allele)[0]:
            AC = 0
            GT = '0|0'
            heterozygous_allele+=str(ref_allele_length)+','+str(ref_allele_length)
            SD = str(allele_count[ref_allele_length])+','+str(allele_count[str(ref_allele_length)])
        else:
            AC = 1; GT = '1|1'
            heterozygous_allele+=str(tuple(final_allele)[0])+','+str(tuple(final_allele)[0])
            SD = str(allele_count[tuple(final_allele)[0]])+','+str(allele_count[str(tuple(final_allele)[0])])
        PC = str(phased_read[0])+','+str(phased_read[1])
    else:
        AN = 2
        if len(set((ref_allele_length,)) & final_allele) == 1:
            AC = 1
            GT = '0|1'
            heterozygous_allele+=str(ref_allele_length)+','+str(tuple(final_allele-{ref_allele_length})[0])
            SD = str(allele_count[ref_allele_length])+','+str(allele_count[tuple(final_allele-{ref_allele_length})[0]])
            if genotypes.index(ref_allele_length) == 0:
                PC = str(phased_read[0])+','+str(phased_read[1])
            else: PC = str(phased_read[1])+','+str(phased_read[0])
        else:
            AC = 2
            GT = '1|2'
            heterozygous_allele+=str(genotypes[0])+','+str(genotypes[1])
            SD = str(allele_count[genotypes[0]])+','+str(allele_count[genotypes[1]])
            PC = str(phased_read[0])+','+str(phased_read[1])

    # INFO = 'AC='+str(AC)+';AN='+str(AN)+';DP='+str(DP)+';END='+str(locus_end)+';PC='+PC+';SN='+str(snp_num)+';SQ='+chosen_snpQ

    INFO = 'AC='+str(AC)+';AN='+str(AN)+';MOTIF=' + str(global_loci_info[locus_key][3]) + ';END='+str(locus_end)

    # FORMAT = 'GT:AL:SD'
    # SAMPLE = str(GT)+':'+heterozygous_allele+':'+SD

    FORMAT = 'GT:AL:SD:PC:DP:SN:SQ'
    SAMPLE = str(GT)+':'+heterozygous_allele+':' + SD + ':' + PC + ':' + str(DP) + ':' + str(snp_num) + ':' + chosen_snpQ

    print(*[contig, locus_start, '.',  ref.fetch(contig, locus_start-1, locus_end), '.', 0, 'PASS', INFO, FORMAT, SAMPLE], file=out, sep='\t')
    del global_loci_info[locus_key]

def vcf_fail_writer(contig, locus_key, global_loci_info, ref, out, DP, skip_point):

    locus_start = int(global_loci_info[locus_key][1])
    locus_end = int(global_loci_info[locus_key][2])

    if skip_point == 0:
        FILTER = 'LESS_READS'    
    locus_key = f'{contig}:{locus_start}-{locus_end}'
    INFO = 'AC=0;AN=0;MOTIF=' + str(global_loci_info[locus_key][3]) + ';END=' + str(locus_end)
    FORMAT = 'GT:AL:SD:PC:DP:SN:SQ'
    SAMPLE = '.:.:.:.:.:.:.'
    print(*[contig, locus_start, '.',  ref.fetch(contig, locus_start-1, locus_end), '.', 0, FILTER, INFO, FORMAT, SAMPLE], file=out, sep='\t')
    del global_loci_info[locus_key]
