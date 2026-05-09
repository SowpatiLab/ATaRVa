import numpy as np

from ATARVA.snp_utils import haplocluster_reads
from ATARVA.vcf_writer import *
from ATARVA.sub_operation_utils import alt_sequence
from ATARVA.somatic_utils import *
from ATARVA.length_utils import *

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

    for i, hap_reads in enumerate(locus_data.hap_read_sets):
        hap_reads = sorted(hap_reads)
        hap_lengths = locus_data.hap_alen_sets[i]
        ALT, allele_length, decomp_seq, is_repetitive = alt_sequence(locus_data.read_aseqs, hap_reads, locus.motif_length)
        lower, upper = (round(x) for x in np.percentile(np.array(hap_lengths), [2.5, 97.5]))
        if i == 0:
            locus_data.gt_aseqs        = (ALT, locus_data.gt_aseqs[1])
            locus_data.gt_alens        = (allele_length, locus_data.gt_alens[1])
            locus_data.gt_decomp_seqs  = (decomp_seq, locus_data.gt_decomp_seqs[1])
            locus_data.hap_meth_data   = (calculate_methylation(hap_reads, locus_data.read_methylation, ALT), locus_data.hap_meth_data[1])
            locus_data.gt_arange       = (f'{lower}-{upper}', locus_data.gt_arange[1])
        else:
            locus_data.gt_aseqs        = (locus_data.gt_aseqs[0], ALT)
            locus_data.gt_alens        = (locus_data.gt_alens[0], allele_length)
            locus_data.gt_decomp_seqs  = (locus_data.gt_decomp_seqs[0], decomp_seq)
            locus_data.hap_meth_data   = (locus_data.hap_meth_data[0], calculate_methylation(hap_reads, locus_data.read_methylation, ALT))
            locus_data.gt_arange       = (locus_data.gt_arange[0], f'{lower}-{upper}')
    write_heterozygous_call(cooper, locus_key)

    return
