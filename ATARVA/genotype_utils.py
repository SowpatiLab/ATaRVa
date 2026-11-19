import numpy as np

from ATARVA.snp_utils           import haplocluster_reads
from ATARVA.vcf_writer          import write_heterozygous_call
from ATARVA.length_utils        import length_genotyper, length_genotyper_kde, length_genotyper_hdbscan
from ATARVA.decompose           import motif_decomposition, longest_pure_repeat
from ATARVA.sub_operation_utils import alt_sequence, calculate_methylation

# NOTE: after genotyping and haplogrouping, before reporting the locus build the ALT sequence
#       and calculate methylation level for the haplotypes


def analyse_genotype(cooper, locus_key):
    """
    genotype the locus based on the read data

    :param cooper: Cooper object
    :param locus_key: locus identifier string
    """

    locus      = cooper.cooper_loci_info[locus_key]
    locus_data = cooper.cooper_loci_data[locus_key]

    if cooper.haploid:
        length_genotyper(cooper, locus_key)
        return

    min_snp = haplocluster_reads(cooper, locus_key)

    if not locus_data.is_genotyped: # if the loci has no significant snps
        if cooper.args.amplicon:
            length_genotyper_kde(cooper, locus_key)
        if not locus_data.is_genotyped: # if the loci has significant snps but not enough to haplotype
            length_genotyper_hdbscan(cooper, locus_key)
        return

    if min_snp != -1:
        snp_left_boundary = locus.start - cooper.args.snp_dist
        min_idx = 0
        for each_spos in cooper.cooper_sorted_snps:
            if each_spos >= snp_left_boundary:
                break
            del cooper.cooper_snp_data[each_spos]
            min_idx += 1
        del cooper.cooper_sorted_snps[:min_idx]

    phased_reads = set()
    for i, hap_reads in enumerate(locus_data.hap_read_sets):
        hap_reads   = sorted(hap_reads)
        phased_reads |= set(hap_reads)
        hap_lengths = locus_data.hap_alen_sets[i]
        ALT, allele_length = alt_sequence(locus_data.read_aseqs, hap_reads, cooper.args.amplicon, cooper.poa_aligner)
        lower, upper = (round(x) for x in np.percentile(np.array(hap_lengths), [2.5, 97.5]))
        ucluster = [str(a) for a in sorted([x for x in hap_lengths if x < lower or x > upper])]
        if i == 0:
            locus_data.gt_aseqs        = (ALT, locus_data.gt_aseqs[1])
            locus_data.gt_alens        = (allele_length, locus_data.gt_alens[1])
            locus_data.hap_meth_data   = (calculate_methylation(cooper, hap_reads, locus_data.read_methylation, ALT), locus_data.hap_meth_data[1])
            locus_data.gt_arange       = (f'{lower}-{upper}', locus_data.gt_arange[1])
            locus_data.gt_ucluster     = ('-'.join(ucluster), locus_data.gt_ucluster[1])
            if cooper.args.decompose and ALT != '<DEL>':
                decomp_seq, nonrep_fraction = motif_decomposition(ALT, locus.motif_length)
                lpm, lpcn  = longest_pure_repeat(decomp_seq, locus.motif)
                locus_data.gt_lpm          = (lpm,  locus_data.gt_lpm[1])
                locus_data.gt_lpcn         = (lpcn, locus_data.gt_lpcn[1])
                locus_data.gt_decomp_seqs  = (decomp_seq, locus_data.gt_decomp_seqs[1])
        else:
            locus_data.gt_aseqs        = (locus_data.gt_aseqs[0], ALT)
            locus_data.gt_alens        = (locus_data.gt_alens[0], allele_length)
            locus_data.hap_meth_data   = (locus_data.hap_meth_data[0], calculate_methylation(cooper, hap_reads, locus_data.read_methylation, ALT))
            locus_data.gt_arange       = (locus_data.gt_arange[0], f'{lower}-{upper}')
            locus_data.gt_ucluster     = (locus_data.gt_ucluster[0], '-'.join(ucluster))
            if cooper.args.decompose and ALT != '<DEL>':
                decomp_seq, nonrep_fraction = motif_decomposition(ALT, locus.motif_length)
                lpm, lpcn  = longest_pure_repeat(decomp_seq, locus.motif)
                locus_data.gt_lpm          = (locus_data.gt_lpm[0],  lpm)
                locus_data.gt_lpcn         = (locus_data.gt_lpcn[0], lpcn)
                locus_data.gt_decomp_seqs  = (locus_data.gt_decomp_seqs[0], decomp_seq)

    ucluster = [len(locus_data.read_aseqs[rid][0]) for rid in locus_data.reads if rid not in phased_reads]
    if ucluster:
        locus_data.gt_ucluster = (locus_data.gt_ucluster[0], locus_data.gt_ucluster[1], '-'.join([str(x) for x in sorted(ucluster)]))

    write_heterozygous_call(cooper, locus_key)

    return
