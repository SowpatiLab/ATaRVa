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
class ReadInfo:
    start:       int   = 0
    end:         int   = 0
    snps:        set   = field(default_factory=set)
    dels:        set   = field(default_factory=list)
    methyl:      list  = field(default_factory=list)
    qual:        int   = 0
    left_flank:  list  = field(default_factory=list)
    right_flank: list  = field(default_factory=list)


@dataclass(slots=True)
class LocusVariation:
    reads:        list  = field(default_factory=list)
    read_alens:   dict  = field(default_factory=dict)
    read_seqs:    dict  = field(default_factory=dict)
    read_hapgp:   list  = field(default_factory=list)
    read_meth:    dict  = field(default_factory=dict)


@dataclass(slots=True)
class SNP:
    cov:  int   = 0
    sub:  dict  = field(default_factory=dict)
    qual: dict  = field(default_factory=dict)


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
        ext.hp_tag                 = [False, None]

        ext.methyl_start           = 0
        ext.methyl_end             = 0
        ext.methyl_range           = []

        return ext
