import sys
import pysam
import time as ti

from ATARVA.version     import __version__
from ATARVA.locus_utils import break_locuskey

SA_DROP_DISTANCE = 1000000    # distance beyond which supplementary alignments are dropped from memory
FORMAT  = 'GT:AL:CN:LPM:LPS:AR:SD:DP:SN:SQ:MA:MR:DS:MV:OL:FV:PS'
def fail_format(depth=0):
    return f'.:.:.:.:.:.:.:{depth}:.:.:.:.:.:.:.:.:.'

INFO_MP_CUTOFF = 0.5
def set_info_mp_cutoff(val):
    global INFO_MP_CUTOFF
    INFO_MP_CUTOFF = val


def vcf_writer(args, out, bam, bam_name):
    """
    Initialize the VCF header and write to output file

    :param out: output file handle
    :param bam: pysam AlignmentFile object for the input BAM file
    :param bam_name: sample name for the BAM file, to be added in the VCF header
    """

    vcf_header = pysam.VariantHeader()
    vcf_header.add_line('##fileformat=VCFv4.5')
    vcf_header.add_line('##fileDate={}'.format(ti.strftime("%Y%m%d")))
    vcf_header.add_line('##source=ATaRVa-v{__version__}')
    if args.fasta:
        vcf_header.add_line('##reference={}'.format(args.fasta))
    else:
        vcf_header.add_line('##reference={}'.format(bam.header.get('SQ', [{}])[0].get('SN', 'unknown')))

    # command
    vcf_header.add_line(f"##command=ATaRVa_{__version__} {' '.join(sys.argv)}")

    for contig in bam.header['SQ']:
        vcf_header.contigs.add(contig['SN'], length=contig['LN'])

    #sample_name
    vcf_header.add_sample(bam_name)

    # FILTER
    vcf_header.filters.add('LESS_READS', number=None, type=None, description="Read depth below threshold")
    vcf_header.filters.add('FAIL',       number=None, type=None, description="Reads not encountered at the locus")

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
    vcf_header.formats.add("GT",   number=1,   type="String",   description="Genotype")
    vcf_header.formats.add("AL",   number=2,   type="Integer",  description="Allele length in base pairs")
    vcf_header.formats.add("CN",   number=2,   type="Integer",  description="Motif copy number for each allele")
    vcf_header.formats.add("LPM",  number=2,   type="String",   description="Longest pure motif repeat")
    vcf_header.formats.add("LPS",  number=2,   type="String",   description="Copy number of longest pure motif repeat")
    vcf_header.formats.add("AR",   number='.', type="String",   description="Allele length range")
    vcf_header.formats.add("SD",   number='.', type="Integer",  description="Number of reads supporting each haplogroup")
    vcf_header.formats.add("DP",   number=1,   type="Integer",  description="Total supporting reads for the repeat locus")
    vcf_header.formats.add("SN",   number='.', type="Integer",  description="Number of informative SNPs used for phasing")
    vcf_header.formats.add("SQ",   number='.', type="Float",    description="Average Phred-scale base call quality of each informative SNP across the supporting reads")
    vcf_header.formats.add("MA",   number='.', type="Float",    description="Average  methylation level for each allele")
    vcf_header.formats.add("MR",   number='.', type="Integer",  description="Number of informative reads for methylation scoring for each allele")
    vcf_header.formats.add("DS",   number='A', type="String",   description="Motif decomposed sequence for each allele sequence")
    vcf_header.formats.add("MV",   number='.', type="String",   description="Base methylation score encoded for visualization for each allele")
    vcf_header.formats.add("OL",   number='1', type="String",   description="The outlier lengths observed in the phased reads and unphased reads for the locus")
    vcf_header.formats.add("FV",   number='1', type="String",   description="Flag indicating presence of variants in the flank of the repeat locus")
    vcf_header.formats.add("PS",   number='1', type="String",   description="Phase Set assigned in the phasing process if available in alignment file")

    out.write(str(vcf_header))


def ordered_write(cooper, locus_key, record):
    """
    Write the VCF record for a locus in order based on the region iterator

    :param cooper: Cooper object containing the locus information
    :param locus_key: key for the locus in the format 'chrom:start-end'
    :param record: list of VCF record fields to be written
    """
    
    if locus_key == cooper.next_locus_key:
        print(*record, file=cooper.outhandle, sep='\t')
        cooper.next_locus_key = next(cooper.region_itr, None)
        while cooper.next_locus_key in cooper.genotyped_loci:
            print(*cooper.genotyped_loci[cooper.next_locus_key], file=cooper.outhandle, sep='\t')
            del cooper.genotyped_loci[cooper.next_locus_key]
            cooper.next_locus_key = next(cooper.region_itr, None)
    else:
        chrom,   start,  end,  motif = break_locuskey(locus_key)
        _chrom, _start, _end, _motif = break_locuskey(cooper.next_locus_key)
        while start - _end > SA_DROP_DISTANCE:
            if cooper.next_locus_key in cooper.genotyped_loci:
                print(*cooper.genotyped_loci[cooper.next_locus_key], file=cooper.outhandle, sep='\t')
                del cooper.genotyped_loci[cooper.next_locus_key]
            else:
                if cooper.args.gvcf:
                    print(cooper.chrom, _start + 1, '.', cooper.ref.fetch(_chrom, _start, _end), '.', 0, 'FAIL',
                        f'AC=0;AN=0;MOTIF={_motif};START={_start};END={_end};ID=.;REFCN={(_end-_start)//len(_motif)}',
                        FORMAT, fail_format(), file=cooper.outhandle, sep='\t')
            cooper.next_locus_key = next(cooper.region_itr, None)
            _chrom, _start, _end, _motif = break_locuskey(cooper.next_locus_key)
        cooper.genotyped_loci[locus_key] = record


def write_fail_call(cooper, locus_key):
    """
    entry for vcf failed call

    :param cooper: Cooper object containing the locus information
    :param locus_key: key for the locus in the format 'chrom:start-end'
    :param skip_point: integer indicating the reason for failure, to be added in the FILTER column of the VCF
    """

    FILTER = ''
    locus = cooper.cooper_loci_info[locus_key]
    depth = 0
    if locus_key not in cooper.cooper_loci_data:
        FILTER = 'LESS_READS'
    else:
        locus_data = cooper.cooper_loci_data[locus_key]
        refcn = str(locus.length // locus.motif_length)

        if locus_data.skip_code == 0:   FILTER = 'LESS_READS'
        depth = locus_data.depth

    # --- locus optional tag ---
    optional_tag = f';ID={locus.name}' if locus.name else ';ID=.'

    INFO   = 'AC=0;AN=0;MOTIF=' + str(locus.motif) + ';START=' + str(locus.start) + ';END=' + str(locus.end) + optional_tag + ';REFCN='+refcn
    SAMPLE = fail_format(depth)

    record = [locus.chrom, locus.start + 1, '.',  cooper.ref.fetch(locus.chrom, locus.start, locus.end), '.', 0, FILTER, INFO, FORMAT, SAMPLE]
    ordered_write(cooper, locus_key, record)


def write_homozygous_call(cooper, locus_key):
    """
    write a homozygous VCF record for a given locus.

    :param cooper:           cooper object
    :param locus_key:        locus identifier string
    """

    locus      = cooper.cooper_loci_info[locus_key]
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
    decomposed_seq = '.'
    lpm    = '.'
    lpm_cn = '.'
    if cooper.args.decompose and locus.motif_length <= 10:
        decomposed_seq = locus_data.gt_decomp_seqs[0] if locus_data.gt_decomp_seqs[0] else '.'
        lpm    = locus_data.gt_lpm[0]  if locus_data.gt_lpm[0]  else '.'
        lpm_cn = locus_data.gt_lpcn[0] if locus_data.gt_lpcn[0] else '.'

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
    FV = '1' if locus_data.flank_var else '0'
    OL = ''
    PS = '.'
    for uclust in locus_data.gt_ucluster:
        if uclust is None: OL += '.,'
        elif len(uclust) == 0: OL += '.,'
        else: OL += uclust + ','
    OL = OL.rstrip(',')
    if locus_data.hap_category == 3 and locus_data.is_phased:
        PS = list(set(locus_data.read_haplotag_ps.values())) if locus_data.read_haplotag_ps is not None else '.'
        if PS != '.': PS = PS[0]
        SAMPLE += f':{PS}'
    SAMPLE = (
            f'{GT}'
            f':{length_GT}'
            f':{units_GT}'
            f':{lpm}'
            f':{lpm_cn}'
            f':{allele_range}'
            f':{hap_depth}'
            f':{depth}'
            f':.:.'
            f':{MA}'
            f':{MR}'
            f':{decomposed_seq}'
            f':{MV}'
            f':{OL}'
            f':{FV}'
            f':{PS}'
        )

    # --- write VCF record ---
    record = [cooper.chrom, locus.start + 1, '.', ref_seq, ALT, 0, 'PASS', INFO, FORMAT, SAMPLE]
    ordered_write(cooper, locus_key, record)


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
    for hap_meth in locus_data.hap_meth_data:
        meth_prob.append(str(hap_meth[0])  if hap_meth[0] is not None else '.') #methylation probability
        meth_reads.append(str(hap_meth[1]) if hap_meth[1] is not None else '.') #number of methylated reads
        meth_viztag.append(hap_meth[2]     if hap_meth[2] is not None else '.') #methylation visual encoding

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
            if locus_data.gt_aseqs[0] != ref_allele:
                allele_index = 0; ref_index = 1
                GT = '1|0'
                allele = locus_data.gt_aseqs[allele_index]
                allele_length = 0 if allele == '<DEL>' else len(allele)
                length_GT += f'{allele_length},{ref_alen}'
                units_GT += f'{allele_length//locus.motif_length},{ref_units}'
            else:
                ref_index = 0; allele_index = 1
                GT = '0|1'
                allele = locus_data.gt_aseqs[allele_index]
                allele_length = 0 if allele == '<DEL>' else len(allele)
                length_GT += f'{ref_alen},{allele_length}'
                units_GT += f'{ref_units},{allele_length//locus.motif_length}'
            AC = 1

            ALT = allele
            SD = f'{len(locus_data.hap_read_sets[0])},{len(locus_data.hap_read_sets[1])}'
            allele_range = f'{locus_data.gt_arange[0]},{locus_data.gt_arange[1]}'
        else:
            AC = '1,1'
            GT = '1|2'
            ALT1 = locus_data.gt_aseqs[0]
            ALT2 = locus_data.gt_aseqs[1]
            ALT_LEN1 = len(ALT1) if ALT1 != "<DEL>" else 0
            ALT_LEN2 = len(ALT2) if ALT2 != "<DEL>" else 0
            length_GT += f'{ALT_LEN1},{ALT_LEN2}'
            units_GT += f'{ALT_LEN1//locus.motif_length},{ALT_LEN2//locus.motif_length}'

            SD = f'{len(locus_data.hap_read_sets[0])},{len(locus_data.hap_read_sets[1])}'

            ALT = f'{ALT1},{ALT2}'
            allele_range = f'{locus_data.gt_arange[0]},{locus_data.gt_arange[1]}'
    MA = ','.join(meth_prob)
    MR = ','.join(meth_reads)
    MV = ','.join(meth_viztag)

    INFO = f'AC={AC};AN={AN};MOTIF={locus.motif};START={locus.start};END={locus.end}{optional_tag};REFCN={ref_units}'

    decomposed_seqs = '.,.'
    lpm = '.,.'
    lpm_cn = '.,.'
    if cooper.args.decompose:
        decomposed_seqs = f'{locus_data.gt_decomp_seqs[0]},{locus_data.gt_decomp_seqs[1]}'
        lpm    = f'{locus_data.gt_lpm[0]},{locus_data.gt_lpm[1]}'
        lpm_cn = f'{locus_data.gt_lpcn[0]},{locus_data.gt_lpcn[1]}'

    num_snps = 0
    snp_quals = ''
    if locus_data.phase_mode == 'snp' and locus_data.phasing_snp_quals:
        num_snps  = locus_data.n_phasing_snps
        snp_quals = locus_data.phasing_snp_quals
    if snp_quals == '': snp_quals = '.'
    FV = '1' if locus_data.flank_var else '0'
    OL = ''
    for uclust in locus_data.gt_ucluster:
        if uclust is None: OL += '.,'
        elif len(uclust) == 0: OL += '.,'
        else: OL += uclust + ','
    OL = OL.rstrip(',')
    PS = '.'
    if locus_data.hap_category == 3 and locus_data.is_phased:
        PS = list(set(locus_data.read_haplotag_ps.values())) if locus_data.read_haplotag_ps is not None else '.'
        if PS != '.':
            PS = PS[0]
        SAMPLE += f':{PS}'
    SAMPLE = (
            f'{GT}'
            f':{length_GT}'
            f':{units_GT}'
            f':{lpm}'
            f':{lpm_cn}'
            f':{allele_range}'
            f':{SD}'
            f':{str(locus_data.depth)}'
            f':{num_snps}'
            f':{snp_quals}'
            f':{MA}'
            f':{MR}'
            f':{decomposed_seqs}'
            f':{MV}'
            f':{OL}'
            f':{FV}'
            f':{PS}'
        )

    record = [cooper.chrom, locus.start + 1, '.',  ref_allele, ALT, 0, 'PASS', INFO, FORMAT, SAMPLE]
    ordered_write(cooper, locus_key, record)
