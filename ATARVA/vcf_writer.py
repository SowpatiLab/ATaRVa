import sys
import pysam
import numpy as np
from ATARVA.decompose import motif_decomposition


INFO_MP_CUTOFF = 0.5
def set_info_mp_cutoff(val):
    global INFO_MP_CUTOFF
    INFO_MP_CUTOFF = val


def vcf_writer(out, bam, bam_name):
    """
    Initialize the VCF header and write to output file

    :param out: output file handle
    :param bam: pysam AlignmentFile object for the input BAM file
    :param bam_name: sample name for the BAM file, to be added in the VCF header
    """

    vcf_header = pysam.VariantHeader()

    # command
    vcf_header.add_line(f"##command=ATaRVa_0.6.0 {' '.join(sys.argv)}")

    for contig in bam.header['SQ']:
        vcf_header.contigs.add(contig['SN'], length=contig['LN'])
    #sample_name
    vcf_header.add_sample(bam_name)

    # FILTER
    vcf_header.filters.add('LESS_READS', number=None, type=None, description="Read depth below threshold")

    # INFO
    vcf_header.info.add("AC",    number='A', type="Integer",  description="Number of alternate alleles in called genotypes")
    vcf_header.info.add("AN",    number=1,   type="Integer",  description="Number of alleles in called genotypes")
    vcf_header.info.add("MOTIF", number=1,   type="String",   description="Repeat motif")
    vcf_header.info.add("START", number=1,   type="Integer",  description="Start position of the repeat region in 0-based coordinate system")
    vcf_header.info.add("END",   number=1,   type="Integer",  description="End position of the repeat region")
    vcf_header.info.add("ID",    number=1,   type="String",   description="Locus identifier tag")
    vcf_header.info.add("REFCN", number=1,   type="Integer",  description="Reference allele copy number")
    vcf_header.info.add("CT",    number=1,   type="String",   description="Cluster type")
    vcf_header.info.add("MPC",   number=1,   type="String",   description=f"{INFO_MP_CUTOFF}")
    vcf_header.info.add("AFD",   number=1,   type="String",   description="Allele length frequency for the reads supporting the locus")

    # FORMAT
    vcf_header.formats.add("GT", number=1,   type="String",   description="Genotype")
    vcf_header.formats.add("AL", number=2,   type="Integer",  description="Allele length in base pairs")
    vcf_header.formats.add("CN", number=2,   type="Integer",  description="Motif copy number for each allele")
    vcf_header.formats.add("AR", number='.', type="String",   description="Allele length range")
    vcf_header.formats.add("SD", number='.', type="Integer",  description="Number of reads supporting each haplogroup")
    vcf_header.formats.add("DP", number=1,   type="Integer",  description="Total supporting reads for the repeat locus")
    vcf_header.formats.add("SN", number='.', type="Integer",  description="Number of informative SNPs used for phasing")
    vcf_header.formats.add("SQ", number='.', type="Float",    description="Average Phred-scale base call quality of each informative SNP across the supporting reads")
    vcf_header.formats.add("MA", number='.', type="Float",    description="Average  methylation level for each allele")
    vcf_header.formats.add("MR", number='.', type="Integer",  description="Number of informative reads for methylation scoring for each allele")
    vcf_header.formats.add("DS", number='A', type="String",   description="Motif decomposed sequence for each allele sequence")
    vcf_header.formats.add("MV", number='.', type="String",   description="Base methylation score encoded for visualization for each allele")

    out.write(str(vcf_header))


def write_fail_call(cooper, locus_key):
    """
    entry for vcf failed call

    :param cooper: Cooper object containing the locus information
    :param locus_key: key for the locus in the format 'chrom:start-end'
    :param skip_point: integer indicating the reason for failure, to be added in the FILTER column of the VCF
    """

    locus = cooper.cooper_loci_info[locus_key]
    locus_data = cooper.cooper_loci_data[locus_key]
    refcn = str(locus.length // locus.motif_length)

    if locus.name is not None:
        optional_tag = f';ID={locus.name}'
    else:
        optional_tag = ';ID=.'

    if locus_data.skip_code == 0: FILTER = 'LESS_READS'

    locus_key = f'{locus.chrom}:{locus.start}-{locus.end}'

    INFO = 'AC=0;AN=0;MOTIF=' + str(locus.motif) + ';START=' + str(locus.start) + ';END=' + str(locus.end) + optional_tag + ';REFCN='+refcn
    FORMAT = 'GT:AL:CN:AR:SD:DP:SN:SQ:MA:MR:DS:MV'
    SAMPLE = '.:.:.:.:.:.:.:.:.:.:.:.'

    print(*[locus.chrom, locus.start + 1, '.',  cooper.ref.fetch(locus.chrom, locus.start, locus.end), '.', 0, FILTER, INFO, FORMAT, SAMPLE], file=cooper.outhandle, sep='\t')


def write_homozygous_call(cooper, locus_key):
    """
    write a homozygous VCF record for a given locus.

    :param cooper:           cooper object
    :param locus_key:        locus identifier string
    """
    locus = cooper.cooper_loci_info[locus_key]
    locus_data = cooper.cooper_loci_data[locus_key]

    # --- locus optional tag ---
    optional_tag = f';ID={locus.name}' if locus.name else ';ID=.'

    # --- methylation fields ---
    if locus_data.hap_meth_data is not None and locus_data.hap_meth_data[0] is not None:
        meth_avg_prob, meth_read_count, meth_vis_enc = locus_data.hap_meth_data[0]
    else:
        meth_avg_prob, meth_read_count, meth_vis_enc = None, None, None

    meth_prob_str  = str(meth_avg_prob)   if meth_avg_prob   is not None else '.'
    meth_reads_str = str(meth_read_count) if meth_read_count is not None else '.'
    meth_vis_str   = meth_vis_enc         if meth_vis_enc    is not None else '.'

    # homozygous — duplicate values for diploid FORMAT consistency

    allele = locus_data.gt_aseqs[0] if locus_data.gt_aseqs else None

    # --- ref / alt ---
    ref_seq    = cooper.ref.fetch(cooper.chrom, locus.start, locus.end)
    is_alt     = allele and ref_seq != allele
    is_seq_alt = is_alt and not allele.startswith('<')

    AC  = 2       if is_alt else 0
    GT  = '1/1'   if is_alt else '0/0'
    ALT = allele  if is_alt else '.'

    # --- copy number fields ---
    allele_length = 0 if allele == '<DEL>' else len(allele)
    ref_units  = locus.length  // locus.motif_length
    motif_copy = allele_length // locus.motif_length

    # --- decomposed sequence ---
    if cooper.args.decompose and is_seq_alt and locus.motif_length <= 10:
        decomposed_seq = locus_data.gt_decomp_seqs[0] if locus_data.gt_decomp_seqs[0] else None
    else: decomposed_seq = '.'

    # --- INFO field ---
    INFO = (
        f'AC={AC};AN=2'
        f';MOTIF={locus.motif}'
        f';START={locus.start}'
        f';END={locus.end}'
        f'{optional_tag}'
        f';REFCN={ref_units}'
    )

    if cooper.args.debug_mode:
        alen_dist   = sorted(locus_data.halen_frequency.items(), key=lambda x: x[1], reverse=True)
        INFO += f';CT=<TAG>'

    # --- SAMPLE field ---
    FORMAT = 'GT:AL:CN:AR:SD:DP:SN:SQ:MA:MR:DS:MV'

    allele_length = 0 if allele == '<DEL>' else len(allele)
    depth         = locus_data.depth
    hap_depth     = depth
    length_GT     = f'{allele_length}' if cooper.haploid else f'{allele_length},{allele_length}'
    units_GT      = f'{motif_copy}' if cooper.haploid else f'{motif_copy},{motif_copy}'
    allele_range  = f'{locus_data.gt_arange[0]}' if cooper.haploid else f'{locus_data.gt_arange[0]},{locus_data.gt_arange[1]}'
    GT            = GT if not cooper.haploid else GT[0]  # convert to haploid GT if needed
    MA = f'{meth_prob_str}'  if cooper.haploid else f'{meth_prob_str},{meth_prob_str}'
    MV = f'{meth_vis_str}'   if cooper.haploid else f'{meth_vis_str},{meth_vis_str}'
    MR = f'{meth_reads_str}' if cooper.haploid else f'{meth_reads_str},{meth_reads_str}'
    SAMPLE = (
            f'{GT}'
            f':{length_GT}'
            f':{units_GT}'
            f':{allele_range}'
            f':{hap_depth}'
            f':{depth}'
            f':.:.'
            f':{MA}'
            f':{MR}'
            f':{decomposed_seq}'
            f':{MV}'
        )

    # --- write VCF record ---
    print(cooper.chrom, locus.start + 1, '.', ref_seq, ALT, 0, 'PASS', INFO, FORMAT, SAMPLE, file=cooper.outhandle, sep='\t')


def write_heterozygous_call(cooper, locus_key):
    """
    write a heterozygous VCF record for a given locus.

    :param cooper:           cooper object
    :param locus_key:        locus identifier string
    """

    locus = cooper.cooper_loci_info[locus_key]
    locus_data = cooper.cooper_loci_data[locus_key]

    if locus.name is not None:
        optional_tag = f';ID={locus.name}'
    else: optional_tag = ';ID=.'

    length_GT = ''
    units_GT = ''
    AC  = 'AC'
    AN  = 2
    GT  = 'GT'
    SD  = 'SD'
    ALT = '.'

    ref_alen   = locus.end - locus.start
    ref_allele = cooper.ref.fetch(cooper.chrom, locus.start, locus.end)
    ref_units = str(ref_alen // locus.motif_length)

    meth_prob   = []
    meth_reads  = []
    meth_viztag = []
    for hap_methyl in locus_data.hap_meth_data:
        meth_prob.append(str(hap_methyl[0]) if hap_methyl[0] is not None else '.') #methylation probability
        meth_reads.append(str(hap_methyl[1]) if hap_methyl[1] is not None else '.') #number of methylated reads
        meth_viztag.append(hap_methyl[2] if hap_methyl[2] is not None else '.') #methylation visual encoding

    if locus_data.gt_aseqs[0] == locus_data.gt_aseqs[1]: # if the two alleles are the same, make it a homozygous call
        allele = locus_data.gt_aseqs[0]
        if ref_allele == allele:
            AC = 0; GT = '0|0'
        else:
            AC = 2; GT = '1|1'
            ALT = allele
        allele_length = 0 if allele == '<DEL>' else len(allele)
        length_GT += f'{allele_length},{allele_length}'
        units_GT += f'{allele_length//locus.motif_length},{allele_length//locus.motif_length}'
        SD = f'{len(locus_data.hap_read_sets[0])},{len(locus_data.hap_read_sets[1])}'
        allele_range = f'{locus_data.gt_arange[0]},{locus_data.gt_arange[1]}'

    else:
        if ref_allele in locus_data.gt_aseqs:
            ref_index = 0; allele_index = 1
            if locus_data.gt_aseqs[0] != ref_allele:
                allele_index = 0
                ref_index = 1
            allele = locus_data.gt_aseqs[allele_index]
            AC = 1
            GT = '0|1'
            allele_length = 0 if allele == '<DEL>' else len(allele)
            length_GT += f'{ref_alen},{allele_length}'
            units_GT += f'{ref_units},{allele_length//locus.motif_length}'

            SD = f'{len(locus_data.hap_read_sets[ref_index])},{len(locus_data.hap_read_sets[allele_index])}'
            ALT = allele
            if ref_index == 1:
                meth_prob = meth_prob[::-1] # reverse the meth_prob to keep the order consistent with GT
                meth_reads = meth_reads[::-1]
                meth_viztag = meth_viztag[::-1]
            allele_range = f'{locus_data.gt_arange[ref_index]},{locus_data.gt_arange[allele_index]}'
        else:
            AC = '1,1'
            GT = '1|2'
            ALT1 = locus_data.gt_aseqs[0]
            ALT2 = locus_data.gt_aseqs[1]
            length_GT += f'{len(ALT1)},{len(ALT2)}'
            units_GT += f'{len(ALT1)//locus.motif_length},{len(ALT2)//locus.motif_length}'

            SD = f'{len(locus_data.hap_read_sets[0])},{len(locus_data.hap_read_sets[1])}'

            ALT = f'{ALT1},{ALT2}'
            allele_range = f'{locus_data.gt_arange[0]},{locus_data.gt_arange[1]}'
    MA = ','.join(meth_prob)
    MR = ','.join(meth_reads)
    MV = ','.join(meth_viztag)

    INFO = f'AC={AC};AN={AN};MOTIF={locus.motif};START={locus.start};END={locus.end}{optional_tag};REFCN={ref_units}'

    decomposed_seqs = '.,.'
    if cooper.args.decompose:
        decomposed_seqs = f'{locus_data.gt_decomp_seqs[0]},{locus_data.gt_decomp_seqs[1]}'

    num_snps = 0
    snp_quals = ''
    if locus_data.phase_mode == 'snp' and locus_data.phasing_snp_quals:
        num_snps  = locus_data.n_phasing_snps
        snp_quals = locus_data.phasing_snp_quals
    FORMAT = 'GT:AL:CN:AR:SD:DP:SN:SQ:MA:MR:DS:MV'
    SAMPLE = (
            f'{GT}'
            f':{length_GT}'
            f':{units_GT}'
            f':{allele_range}'
            f':{SD}'
            f':{str(locus_data.depth)}'
            f':{num_snps}'
            f':{snp_quals}'
            f':{MA}'
            f':{MR}'
            f':{decomposed_seqs}'
            f':{MV}'
        )

    print(*[cooper.chrom, locus.start + 1, '.',  ref_allele, ALT, 0, 'PASS', INFO, FORMAT, SAMPLE], file=cooper.outhandle, sep='\t')


def vcf_multizygous_writer(contig, genotype_dict, locus_start, locus_end, DP, global_loci_info, ref, out, log_bool, decomp, hallele_counter):

    locus_key = f'{contig}:{locus_start}-{locus_end}'

    tag = "correlation_clustering"

    if len(global_loci_info[locus_key]) > 5:
        optional_tag = f';ID={global_loci_info[locus_key][5]}'
    else:
        optional_tag = ';ID=.'

    motif_size = int(float(global_loci_info[locus_key][4]))

    GT_dict = {}
    gt_idx = 0
    ref_allele_length = locus_end - locus_start
    refcn = str(ref_allele_length // int(float(global_loci_info[locus_key][4])))
    ref_seq = ref.fetch(contig, locus_start, locus_end)
    for each_genotype in genotype_dict:
        current_gt = genotype_dict[each_genotype]
        if int(each_genotype) == ref_allele_length:
            if ref_seq == current_gt[0]:
                GT_dict[0] = (current_gt[0], str(each_genotype), current_gt[3], f'{current_gt[1][0]}-{current_gt[1][1]}', current_gt[4][0], current_gt[4][1], current_gt[2], current_gt[4][2])
            else:
                gt_idx += 1
                GT_dict[gt_idx] = (current_gt[0], str(each_genotype), current_gt[3], f'{current_gt[1][0]}-{current_gt[1][1]}', current_gt[4][0], current_gt[4][1], current_gt[2], current_gt[4][2])
        else:
            gt_idx += 1
            GT_dict[gt_idx] = (current_gt[0], str(each_genotype), current_gt[3], f'{current_gt[1][0]}-{current_gt[1][1]}', current_gt[4][0], current_gt[4][1], current_gt[2], current_gt[4][2])
    del genotype_dict

    GT = []
    if gt_idx> 0:
        AN = (gt_idx + 1) if (0 in GT_dict) else gt_idx
    else:
        AN = 2
    # AN = (gt_idx + 1) if gt_idx>0 else 2
    AC = ','.join(['1']*gt_idx) if gt_idx>0 else 0
    ALT = []
    AL = []
    AR = []
    SD = []
    MA = []
    MR = []
    deseq = []
    MV = []
    for gt_key in sorted(GT_dict.keys()):
        GT.append(str(gt_key))
        if gt_key != 0:
            ALT.append(GT_dict[gt_key][0])
            deseq.append(GT_dict[gt_key][6] if GT_dict[gt_key][6] else '.')
        AL.append(GT_dict[gt_key][1])
        SD.append(str(GT_dict[gt_key][2]))
        AR.append(GT_dict[gt_key][3])
        MA.append(str(GT_dict[gt_key][4]) if GT_dict[gt_key][4] is not None else '.')
        MR.append(str(GT_dict[gt_key][5]) if GT_dict[gt_key][5] is not None else '.')
        MV.append(str(GT_dict[gt_key][7]) if GT_dict[gt_key][5] is not None else '.')
    del GT_dict

    GT = '/'.join(GT)
    ALT = ','.join(ALT) if ALT else '.'
    CN = ','.join([str(i // motif_size) for i in AL])
    AL = ','.join(AL)
    AR = ','.join(AR)
    SD = ','.join(SD)
    MA = ','.join(MA)
    MR = ','.join(MR)
    MV = ','.join(MV)

    if log_bool:
        eac = sorted(hallele_counter.items(), key = lambda x: x[1], reverse=True)
        INFO = 'AC='+str(AC)+';AN='+str(AN)+';MOTIF=' + str(global_loci_info[locus_key][3]) + ';START=' + str(locus_start) + ';END='+str(locus_end) + optional_tag + ';REFCN='+refcn + ';CT=' + tag + ';EAC=' + str(eac)
    else:
        INFO = 'AC='+str(AC)+';AN='+str(AN)+';MOTIF=' + str(global_loci_info[locus_key][3]) + ';START=' + str(locus_start) + ';END='+str(locus_end) + optional_tag + ';REFCN='+refcn

    if decomp:
        deseq = ','.join(deseq) if deseq else '.'
    else:
        deseq = '.'

    FORMAT = 'GT:AL:CN:AR:SD:DP:SN:SQ:MA:MR:DS:MV'
    SAMPLE = GT + ':' + AL + ':' + CN + ':' + AR + ':' + SD + ':' + str(DP) + ':.:.:' + MA + ':' + MR + ':' + deseq + ':' + MV

    print(*[contig, locus_start+1, '.',  ref_seq, ALT, 0, 'PASS', INFO, FORMAT, SAMPLE], file=out, sep='\t')

    del GT, ALT, AL, CN, AR, SD, MA, MR, MV, deseq, global_loci_info[locus_key]
