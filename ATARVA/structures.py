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

    def process_methylation_info(self, adj_fqs, adj_fqe, lower=0.2, upper=0.8):
        """
        Process methylation data for storing in locus data for genotyping and visualization
        :param adj_fqs: adjusted flank start position of the read sequence supporting the locus after accounting for softclips and indels
        :param adj_fqe: adjusted flank end position of the read sequence supporting the locus after accounting for softclips and indels
        :param lower: lower threshold for methylation probability to consider as methylated
        :param upper: upper threshold for methylation probability to consider as unmethylated
        :return: tuple of (average methylation level, list of methylation positions, list of methylation encodings for visualization)
        
        Note: The methylation encoding is
              -   -1 for uncertain methylation status (probability between lower and upper), 
              -    1 for methylated (probability >= upper)
              -    0 for unmethylated (probability <= lower)
        """
        read_mod_bases = self.methylation

        meth_pos = []; meth_scores = []; meth_calls = [];
        meth_count = meth_qual = 0
        for pos, raw_prob in read_mod_bases:
            if not (adj_fqs <= pos <= adj_fqe): continue
            prob = raw_prob / 255
            meth_pos.append(pos - adj_fqs)
            meth_scores.append(raw_prob)
            if lower < prob < upper:
                meth_calls.append(-1)
            else:
                meth_count += 1
                is_meth     = prob >= upper
                meth_scores.append(raw_prob)
                meth_calls.append(int(is_meth))
                meth_qual  += is_meth

        if meth_count: return (round(meth_qual / meth_count, 2), meth_pos, meth_scores, meth_calls)
        else: return None


@dataclass(slots=True)
class LocusVariation:
    reads:             list  = field(default_factory=list)  # list of all informative read indices for the locus
    read_haplotags:    list  = field(default_factory=list)  # the haplotag assigned to each read in reads 
    depth:             int   = 0                            # depth of the locus, updated when reads are subset for high coverage loci  
    read_alens:        dict  = field(default_factory=dict)  # read index -> allele length dict for the reads supporting the locus
    read_aseqs:        dict  = field(default_factory=dict)  # read index -> allele sequence dict for the reads supporting the locus, in the format (allele sequence, allele length) 
    read_methylation:  dict  = field(default_factory=dict)  # read index -> (total methylation probability, methylation positions, methylation scores, methylation calls(0 or 1))
    halen_frequency:   dict  = field(default_factory=dict)  # allele length -> count dict for the reads supporting the locus
    alen_frequency:    dict  = field(default_factory=dict)  # allele length ignoring indels in homopolymers
    allele_lengths:    list  = field(default_factory=list)  # list of allele lengths for the reads supporting the locus

    neighbors:        set   = field(default_factory=set)    # positions of neighbouring loci with shared reads

    # updated when locus has high coverage and reads are subset
    raw_depth:         int   = 0                            # original depth of the locus before subsetting
    raw_reads:         list  = field(default_factory=list)  # original list of read indices supporting the locus before subsetting 
    raw_haplotags:     list  = field(default_factory=list)  # original list of haplotags for the reads supporting the locus before subsetting

    # data for genotyping
    hap_read_sets:     tuple = ([], [])       # (hap1 read indices, hap2 read indices)
    hap_alen_sets:     tuple = (None, None)   # (lengths in hap1 reads, lengths in hap2 reads)
    hap_meth_data:     tuple = (None, None)   # (hap1 meth data, hap2 meth data) methylation data for the genotype call, in the format (average methylation level, number of reads with methylation info, encoded methylation string for visualization)
    is_genotyped:      int   = 0              # whether the locus was genotyped successfully
    hap_category:      int   = None           # 1: homozygous, 2: ambiguous, 3: haplotag
    gt_depth:          int   = 0              # reads used for haplogrouping after filtering out
    gt_aseqs:          tuple = (None, None)   # (allele1, allele2) allele sequences for the genotype call
    gt_alens:          tuple = (None, None)   # (alen1, alen2) allele lengths for the genotype call
    gt_arange:         tuple = (None, None)   # allele length range for the genotype call, in the format 'lower1-upper1,lower2-upper2'
    gt_decomp_seqs:    tuple = (None, None)   # (decomposed allele1, decomposed allele2) motif decomposition of the allele sequences
    is_phased:         bool  = False          # whether haplotagging was successful
    skip_code:         int   = 10             # whether genotyping failed for the locus
    n_phasing_snps:    int   = 0              # number of SNPs used for phasing, if phased by SNPs
    phasing_snp_quals: str   = ''             # comma-separated string of SNP quality scores, if phased by SNPs
    phase_mode:        str   = None           # 0: phased by SNPs, 1: phased by length    


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
        ext.query_name                = read.query_name
        ext.query_sequence            = read.query_sequence
        ext.flag                      = read.flag
        ext.reference_id              = read.reference_id
        ext.reference_start           = read.reference_start
        ext.ref_start                 = read.reference_start
        ext.ref_end                   = read.reference_end
        ext.query_start               = read.query_alignment_start
        ext.query_end                 = read.query_alignment_end
        ext.mapping_quality           = read.mapping_quality
        ext.cigar                     = read.cigar
        ext.next_reference_id         = read.next_reference_id
        ext.next_reference_start      = read.next_reference_start
        ext.template_length           = read.template_length
        ext.query_qualities           = read.query_qualities
        ext.tags                      = read.tags
        ext.has_tag                   = read.has_tag
        if read.has_tag('cs'):
            ext.cs_tag            = read.cs_tag
        if read.has_tag('MD'):
            ext.md_tag            = read.get_tag('MD')

        # custom attributes
        ext.index                  = None
        ext.chrom                  = read.reference_name
        ext.ref_start              = read.reference_start
        ext.ref_end                = read.reference_end
        ext.query_start            = read.query_alignment_start
        ext.query_end              = read.query_alignment_end
        ext.sequence               = read.query_sequence
        ext.mean_qual              = int(sum(read.query_qualities)/len(read.query_qualities)) if read.query_qualities else 0
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
