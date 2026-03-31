from dataclasses import dataclass, field
import pysam
import numpy as np


@dataclass(slots=True)
class ReadLocusInfo:
    halen: int  = 0
    alen:  int  = 0
    rlen:  int  = 0
    seq:   list = field(default_factory=list)

@dataclass(slots=True)
class LocusInfo:
    chrom:        str = ''
    start:        int = 0
    end:          int = 0
    motif:        str = ''
    length:       int = 0
    motif_length: int = 0
    name:         str = None

    def __post_init__(self):
        self.length       = self.end - self.start
        self.motif_length = len(self.motif)

@dataclass(slots=True)
class ReadInfo:
    start:       int   = 0
    end:         int   = 0
    snps:        set   = field(default_factory=set)
    dels:        set   = field(default_factory=list)
    methylation: list  = field(default_factory=list)
    mean_qual:   int   = 0
    left_flank:  list  = field(default_factory=list)
    right_flank: list  = field(default_factory=list)


@dataclass(slots=True)
class LocusVariation:
    reads:             list  = field(default_factory=list)
    read_haplotags:    list  = field(default_factory=list)
    depth:             int   = 0
    read_alens:        dict  = field(default_factory=dict)
    read_seqs:         dict  = field(default_factory=dict)
    read_methylation:  dict  = field(default_factory=dict)  # read index -> (total methylation probability, methylation positions, methylation encodings)
    alen_frequency:    dict  = field(default_factory=dict)
    halen_frequency:   dict  = field(default_factory=dict)
    allele_lengths:    list  = field(default_factory=list)

    neighbors:        set   = field(default_factory=set)   # positions of neighbouring loci with shared reads

    # updated when locus has high coverage and reads are subset
    raw_depth:         int   = 0
    raw_reads:         list  = field(default_factory=list)
    raw_haplotags:     list  = field(default_factory=list)

    # data for genotyping
    haplotypes:            tuple = ([], [])       # (hap1 read indices, hap2 read indices)
    haplotype_lengths:     tuple = (None, None)   # (lengths in hap1 reads, lengths in hap2 reads)
    haplotype_alleles:     tuple = (None, None)   # (allele1, allele2) allele sequences for the genotype call
    haplotype_alens:       tuple = (None, None)   # (alen1, alen2) allele lengths for the genotype call
    haplotype_arange:      tuple = (None, None)   # allele length range for the genotype call, in the format 'lower1-upper1,lower2-upper2'
    decomposed_alleles:    tuple = (None, None)   # (decomposed allele1, decomposed allele2) motif decomposition of the allele sequences
    haplotype_methyldata:  tuple = (None, None)   # (hap1 meth data, hap2 meth data) methylation data for the genotype call, in the format (average methylation level, number of reads with methylation info, encoded methylation string for visualization)
    hap_status:            bool  = False          # whether haplotagging was successful
    num_snps:              int   = 0              # number of SNPs used for phasing, if phased by SNPs
    snp_quals:             str   = ''             # comma-separated string of SNP quality scores, if phased by SNPs
    hap_category:          int   = None           # 0 or 1
    hap_depth:             int   = 0              # reads used for haplogrouping after filtering out
    phase_mode:            str   = None           # 0: phased by SNPs, 1: phased by length
    homozygous_alen:       int   = None           # if hap_category is 1, the allele length of the homozygous locus
    fail_code:             int   = 10             # whether genotyping failed for the locus
    genotyped:             int   = 0              # whether the locus was genotyped successfully


@dataclass(slots=True)
class SNP:
    cov:  int   = 0
    sub:  dict  = field(default_factory=dict)
    qual: dict  = field(default_factory=dict)
    ref:  set   = field(default_factory=set)


class ExtendedRead(pysam.AlignedSegment):

    @classmethod
    def from_read(cls, read: pysam.AlignedSegment):
        """Factory method to create ExtendedRead from a pysam read"""
        ext        = cls()
        ext._read  = read

        # copy core pysam attributes
        ext.query_name            = read.query_name
        ext.query_sequence        = read.query_sequence
        ext.flag                  = read.flag
        ext.reference_id          = read.reference_id
        ext.reference_start       = read.reference_start
        ext.mapping_quality       = read.mapping_quality
        ext.cigar                 = read.cigar
        ext.next_reference_id     = read.next_reference_id
        ext.next_reference_start  = read.next_reference_start
        ext.template_length       = read.template_length
        ext.query_qualities       = read.query_qualities
        ext.tags                  = read.tags

        # custom attributes
        ext.index                  = None
        ext.chrom                  = read.reference_name
        ext.ref_start              = read.reference_start
        ext.ref_end                = read.reference_end
        ext.query_start            = read.query_alignment_start
        ext.query_end              = read.query_alignment_end
        ext.sequence               = read.query_sequence
        ext.mean_qual              = int(np.mean(read.query_qualities)) if read.query_qualities else 0
        ext.loci                   = []
        ext.loci_coords            = []
        ext.loci_keys              = []
        ext.loci_data              = dict()
        ext.homopolymer_positions  = dict()
        ext.left_flanks            = []
        ext.right_flanks           = []
        ext.haplotag               = [False, None]

        ext.methyl_start           = 0
        ext.methyl_end             = 0
        ext.methylation_calls      = []     # list of tuples (position, probability) for methylation calls in the read

        return ext
