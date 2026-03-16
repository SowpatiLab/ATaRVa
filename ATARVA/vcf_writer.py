import sys
import pysam
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
    vcf_header.info.add("AC", number='A', type="Integer", description="Number of alternate alleles in called genotypes")
    vcf_header.info.add("AN", number=1, type="Integer", description="Number of alleles in called genotypes")
    vcf_header.info.add("MOTIF", number=1, type="String", description="Repeat motif")
    vcf_header.info.add("START", number=1, type="Integer", description="Start position of the repeat region in 0-based coordinate system")
    vcf_header.info.add("END", number=1, type="Integer", description="End position of the repeat region")
    vcf_header.info.add("ID", number=1, type="String", description="Locus identifier tag")
    vcf_header.info.add("REFCN", number=1, type="Integer", description="Reference allele copy number")
    vcf_header.info.add("CT", number=1, type="String", description="Cluster type")
    vcf_header.info.add("EAC", number=1, type="String", description="Each Allele Count")
    vcf_header.info.add("MPC", number=1, type="String", description=f"{INFO_MP_CUTOFF}")

    # FORMAT
    vcf_header.formats.add("GT", number=1, type="String", description="Genotype")
    vcf_header.formats.add("AL", number=2, type="Integer", description="Allele length in base pairs")
    vcf_header.formats.add("CN", number=2, type="Integer", description="Motif copy number for each allele")
    vcf_header.formats.add("AR", number='.', type="String", description="Allele length range")
    vcf_header.formats.add("SD", number='.', type="Integer", description="Number of reads supporting for the alleles")
    vcf_header.formats.add("PC", number=2, type="Integer", description="Number of reads in the phased cluster for each allele")
    vcf_header.formats.add("DP", number=1, type="Integer", description="Number of the supporting reads for the repeat locus")
    vcf_header.formats.add("SN", number='.', type="Integer", description="Number of SNPs used for phasing")
    vcf_header.formats.add("SQ", number='.', type="Float", description="Phred-scale qualities of the SNPs used for phasing")
    vcf_header.formats.add("MA", number='.', type="Float", description="Mean methylation level for each allele")
    vcf_header.formats.add("MR", number='.', type="Integer", description="Number of reads providing methylation info for each allele")
    vcf_header.formats.add("DS", number='A', type="String", description="Motif decomposed sequence")
    vcf_header.formats.add("MV", number='.', type="String", description="Visual methylation encodings for the alleles")

    out.write(str(vcf_header))


def write_fail_call(cooper, locus_key, skip_point):
    """
    entry for vcf failed call

    :param cooper: Cooper object containing the locus information
    :param locus_key: key for the locus in the format 'chrom:start-end'
    :param skip_point: integer indicating the reason for failure, to be added in the FILTER column of the VCF
    """

    locus = cooper.cooper_loci_info[locus_key]
    refcn = str(locus.length // locus.motif_length)

    if locus.name is not None:
        optional_tag = f';ID={locus.name}'
    else:
        optional_tag = ';ID=.'

    if skip_point == 0: FILTER = 'LESS_READS'
          
    locus_key = f'{locus.chrom}:{locus.start}-{locus.end}'

    INFO = 'AC=0;AN=0;MOTIF=' + str(locus.motif) + ';START=' + str(locus.start) + ';END=' + str(locus.end) + optional_tag + ';REFCN='+refcn
    FORMAT = 'GT:AL:CN:AR:SD:DP:SN:SQ:MA:MR:DS:MV'
    SAMPLE = '.:.:.:.:.:.:.:.:.:.:.:.'

    print(*[locus.chrom, locus.start + 1, '.',  cooper.ref.fetch(locus.chrom, locus.start, locus.end), '.', 0, FILTER, INFO, FORMAT, SAMPLE], file=cooper.outhandle, sep='\t')


def write_homozygous_call(cooper, locus_key, allele_length, hallele_counter, allele_range,
                          methylation_data, decomposed_seq, ALT_seq, reads_len, depth):
    """
    write a homozygous VCF record for a given locus.

    :param cooper:           cooper object
    :param locus_key:        locus identifier string
    :param allele_length:    length of the called allele
    :param hallele_counter:  allele length counter dict
    :param allele_range:     allele length range string
    :param methylation_data: tuple of (prob, read_count, vis_encoding)
    :param decomposed_seq:   motif-decomposed sequence
    :param ALT_seq:          alternative allele sequence
    :param reads_len:        number of reads supporting the call
    :param depth:            total read depth at locus
    """
    locus = cooper.cooper_loci_info[locus_key]

    # --- locus optional tag ---
    optional_tag = f';ID={locus.name}' if locus.name else ';ID=.'

    # --- methylation fields ---
    meth_avg_prob, meth_read_count, meth_vis_enc = methylation_data

    meth_prob_str  = str(meth_avg_prob)   if meth_avg_prob   is not None else '.'
    meth_reads_str = str(meth_read_count) if meth_read_count is not None else '.'
    meth_vis_str   = meth_vis_enc         if meth_vis_enc    is not None else '.'

    # homozygous — duplicate values for diploid FORMAT consistency
    MA = f'{meth_prob_str},{meth_prob_str}'
    MV = f'{meth_vis_str},{meth_vis_str}'

    # --- ref / alt ---
    ref_seq    = cooper.ref.fetch(cooper.chrom, locus.start, locus.end)
    is_alt     = ALT_seq and ref_seq != ALT_seq
    is_seq_alt = is_alt and not ALT_seq.startswith('<')

    AC  = 2     if is_alt else 0
    GT  = '1/1' if is_alt else '0/0'
    ALT = ALT_seq if is_alt else '.'

    # --- copy number fields ---
    ref_cn     = locus.length  // locus.motif_length
    motif_copy = allele_length // locus.motif_length

    # --- decomposed sequence ---
    if cooper.args.decompose and is_seq_alt and locus.motif_length <= 10:
        deseq = decomposed_seq or motif_decomposition(ALT, locus.motif_length)[0]
    else: deseq = '.'

    # --- INFO field ---
    INFO = (
        f'AC={AC};AN=2'
        f';MOTIF={locus.motif}'
        f';START={locus.start}'
        f';END={locus.end}'
        f'{optional_tag}'
        f';REFCN={ref_cn}'
    )

    if cooper.args.debug_mode:
        eac   = sorted(hallele_counter.items(), key=lambda x: x[1], reverse=True)
        INFO += f';CT=<TAG>;EAC={eac}'

    # --- SAMPLE field ---
    FORMAT = 'GT:AL:CN:AR:SD:DP:SN:SQ:MA:MR:DS:MV'

    if not cooper.haploid:
        SAMPLE = (
            f'{GT}'
            f':{allele_length},{allele_length}'
            f':{motif_copy},{motif_copy}'
            f':{allele_range}'
            f':{reads_len}'
            f':{depth}'
            f':.:.'
            f':{MA}'
            f':{meth_reads_str}'
            f':{deseq}'
            f':{MV}'
        )
    else:
        SAMPLE = (
            f'{GT[0]}'
            f':{allele_length}'
            f':{motif_copy}'
            f':{allele_range}'
            f':{reads_len}'
            f':{depth}'
            f':.:.'
            f':{meth_prob_str}'
            f':{meth_reads_str}'
            f':{deseq}'
            f':{MV}'
        )

    # --- write VCF record ---
    print(cooper.chrom, locus.start + 1, '.', ref_seq, ALT, 0, 'PASS', INFO, FORMAT, SAMPLE, file=cooper.outhandle, sep='\t')


def vcf_heterozygous_writer(contig, genotypes, locus_start, locus_end, allele_count, DP, global_loci_info, ref, out,
                            chosen_snpQ, phased_read, snp_num, ALT_reads, log_bool, tag, decomp, hallele_counter, allele_range,
                            decomp_seq, meth_info):

    locus_key = f'{contig}:{locus_start}-{locus_end}'

    if len(global_loci_info[locus_key]) > 5:
        optional_tag = f';ID={global_loci_info[locus_key][5]}'
    else:
        optional_tag = ';ID=.'

    motif_size = int(float(global_loci_info[locus_key][4]))

    final_allele = set(genotypes)
    heterozygous_allele = ''
    AC = 'AC'
    AN = 2
    GT = 'GT'
    SD = 'SD'
    PC = 'PC'
    ALT = '.'
    alt_seqs = []

    ref_allele_length = locus_end - locus_start
    refcn = str(ref_allele_length // int(float(global_loci_info[locus_key][4])))

    meth_prob = []
    meth_reads = []
    meth_vistag = []
    for each_meth in meth_info:
        meth_prob.append(str(each_meth[0]) if each_meth[0] is not None else '.') #methylation probability
        meth_reads.append(str(each_meth[1]) if each_meth[1] is not None else '.') #number of methylated reads
        meth_vistag.append(each_meth[2] if each_meth[2] is not None else '.') #methylation visual encoding

    if len(final_allele) == 1:

        if ref_allele_length == tuple(final_allele)[0]:
            AC = 0
            GT = '0|0'
            heterozygous_allele+=str(ref_allele_length)+','+str(ref_allele_length)
            SD = str(allele_count[ref_allele_length])+','+str(allele_count[str(ref_allele_length)])
            alt_seqs.append('')
        else:
            AC = 2; GT = '1|1'
            heterozygous_allele+=str(tuple(final_allele)[0])+','+str(tuple(final_allele)[0])
            SD = str(allele_count[tuple(final_allele)[0]])+','+str(allele_count[str(tuple(final_allele)[0])])

            ALT = ALT_reads[0]
            if ALT[0]!='<': alt_seqs.append(ALT)
            else: alt_seqs.append('')
        PC = str(phased_read[0])+','+str(phased_read[1])
        MA = ','.join(meth_prob)
        MR = ','.join(meth_reads)
        MV = ','.join(meth_vistag)
    else:

        if len(set((ref_allele_length,)) & final_allele) == 1:
            AC = 1
            GT = '0|1'
            heterozygous_allele+=str(ref_allele_length)+','+str(tuple(final_allele-{ref_allele_length})[0])
            SD = str(allele_count[ref_allele_length])+','+str(allele_count[tuple(final_allele-{ref_allele_length})[0]])
            if genotypes.index(ref_allele_length) == 0:
                PC = str(phased_read[0])+','+str(phased_read[1])

                alt_seqs.append(None) # dummy added for ref, to keep the length of alt_seqs as 2
                ALT = ALT_reads[1]
                if ALT[0]!='<': alt_seqs.append(ALT)
                else: alt_seqs.append('')
                MA = ','.join(meth_prob)
                MR = ','.join(meth_reads)
                MV = ','.join(meth_vistag)
            else:
                PC = str(phased_read[1])+','+str(phased_read[0])

                ALT = ALT_reads[0]
                if ALT[0]!='<': alt_seqs.append(ALT)
                else: alt_seqs.append('')
                alt_seqs.append(None) # dummy added for ref, to keep the length of alt_seqs as 2
                allele_range = ','.join(allele_range.split(',')[::-1]) # reverse the allele range to keep the order consistent with GT
                MA = ','.join(meth_prob[::-1]) # reverse the meth_prob to keep the order consistent with GT
                MR = ','.join(meth_reads[::-1])
                MV = ','.join(meth_vistag[::-1])
        else:
            AC = '1,1'
            GT = '1|2'
            heterozygous_allele+=str(genotypes[0])+','+str(genotypes[1])
            SD = str(allele_count[genotypes[0]])+','+str(allele_count[genotypes[1]])
            PC = str(phased_read[0])+','+str(phased_read[1])

            ALT1 = ALT_reads[0]
            if ALT1[0]!='<': alt_seqs.append(ALT1)
            else: alt_seqs.append('')
                
            ALT2 = ALT_reads[1]
            if ALT2[0]!='<': alt_seqs.append(ALT2)
            else: alt_seqs.append('')

            ALT = ALT1 + ',' + ALT2
            MA = ','.join(meth_prob)
            MR = ','.join(meth_reads)
            MV = ','.join(meth_vistag)


    if PC == '.,.': PC = '.' # due to length genotyper
    if log_bool:
        eac = sorted(hallele_counter.items(), key = lambda x: x[1], reverse=True)
        INFO = 'AC='+str(AC)+';AN='+str(AN)+';MOTIF=' + str(global_loci_info[locus_key][3]) + ';START=' + str(locus_start) + ';END='+str(locus_end) + optional_tag + ';REFCN='+ refcn + ';CT=' + tag + ';EAC=' + str(eac)
    else:
        INFO = 'AC='+str(AC)+';AN='+str(AN)+';MOTIF=' + str(global_loci_info[locus_key][3]) + ';START=' + str(locus_start) + ';END='+str(locus_end) + optional_tag + ';REFCN='+ refcn

    deseq = '.,.'
    if decomp:
        motif_size = int(float(global_loci_info[locus_key][4]))
        
        if motif_size>10:
            deseq = ','.join(['.']*len(alt_seqs))
        else:
            ds = []
            for index,iseq in enumerate(alt_seqs):
                if iseq:
                    if all(decomp_seq):
                        ds.append(decomp_seq[index])
                    else:
                        i_deseq,_ = motif_decomposition(iseq, motif_size)
                        ds.append(i_deseq)
                elif iseq=='':
                    ds.append('.')
            deseq = ','.join(ds)

    motif_copy = ','.join([str(int(i) // motif_size) for i in heterozygous_allele.split(',')])
     
    FORMAT = 'GT:AL:CN:AR:SD:DP:SN:SQ:MA:MR:DS:MV'
    SAMPLE = str(GT)+':'+heterozygous_allele+':' + motif_copy + ':' + allele_range + ':' + SD + ':' + str(DP) + ':' + str(snp_num) + ':' + chosen_snpQ + ':' + MA + ':' + MR + ':' + deseq + ':' + MV

    del ALT_reads
    del alt_seqs

    print(*[contig, locus_start+1, '.',  ref.fetch(contig, locus_start, locus_end), ALT, 0, 'PASS', INFO, FORMAT, SAMPLE], file=out, sep='\t')


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