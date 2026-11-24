# Genotype Usage

## Table of contents

* [Usage](#usage)
* [Required Inputs: Overview and Preprocessing](#required-inputs-overview-and-preprocessing)
* [Advanced options](#advanced-options)
* [Examples](#examples)
* [Output format](#output-format)

## Usage
The help message and available subcommands can be accessed using

```bash
$ atarva -h
#  or
$ atarva --help
```
which gives the following output

```
ATaRVa - Analysis of Tandem Repeat Variants
Sowpati Lab & Dashnow Lab

Usage:
    atarva [OPTIONS] <COMMAND>

Commands:
  genotype  Tandem Repeat Genotyper
  merge     Merge ATaRVa VCF files

Options:
  -h, --help     Print help
  -v, --version  Print version
```

## Required Inputs: Overview and Preprocessing
`atarva genotype` accepts read alignments and a TR region catalogue of interest and output a variant file with genotype information, including the consensus sequence, allele length, and decomposed motifs.

Overview of the ATaRVa workflow:
1. ATaRVa processes the input BAM/CRAM file read-wise, assuming that most long reads span multiple TR loci.
2. After identifying the flank alignments and adjustment of any flanking indels in each read, ATaRVa calculates allele length for a locus across all spanning reads.
3. Next, ATaRVa phases reads into haplotypes. ATaRVa attempts to group reads using nearby informative SNVs. If SNP-based genotyping is unsuccessful, it uses HDBSCAN clustering based on allele length and edit distance from the reference (KDE clustering for amplicon mode). If the input file is already phased the input haplotag can be specified with --haplotag.
4. It derives consensus allele sequences for each haplotype using partial order alignment, decomposes each TR allele into motif-level representations, and outputs the results in VCF format.

The help message and available options can be accessed using

```bash
$ atarva genotype -h
#  or
$ atarva genotype --help
```
which gives the following output

```
usage: atarva genotype [-h] -f <FILE> -b <FILE> [<FILE> ...] -r <FILE> [-o <FILE>] [--aln-format <STR>]
                       [--use-softclip] [--contigs CONTIGS [CONTIGS ...]] [--karyotype <STR> [<STR> ...]]
                       [--map-qual <INT>] [--min-reads <INT>] [--max-reads <INT>] [--flank <INT>] [--snp-dist <INT>]
                       [--snp-count <INT>] [--snp-qual <INT>] [--phased-frac <FLOAT>] [--cluster-frac <FLOAT>]
                       [--haplotag <STR>] [--meth-prob <FLOAT>] [--methviz] [--locus-wise] [--decompose]
                       [--amplicon] [--strict] [--gvcf] [-t <INT>] [-log] [-v]

Required arguments:
  -f <FILE>, --fasta <FILE>
                        input reference fasta file
  -b <FILE> [<FILE> ...], --bam <FILE> [<FILE> ...]
                        sample alignment files [SAM | BAM | CRAM]
  -r <FILE>, --regions <FILE>
                        bgzipped + tabix-indexed regions file. Prepare: sort with bedtools → bgzip → tabix index
```
The details of required input files are given below:

### Reference genome
#### `-f or --fasta`
**Expects**: *FILE*<br>
**Default**: *None*<br>
The `-f` or `--fasta` option is used to specify the input FASTA file. The corresponding index file (`.fai`) should be in the same directory. ATaRVa uses [pysam](https://github.com/pysam-developers/pysam)'s `FastaFile` parser to read the input FASTA file.

### Alignment file
#### `-b or --bam`
**Expects**: *FILE*<br>
**Default**: *None*<br>
The `-b` or `--bam` option is used to specify one or more input alignment files in the same format. ATaRVa accepts any of the three alignment formats: SAM, BAM, or CRAM. The alignment file should be sorted by coordinates. The format should be specified using the `--aln-format` option. The corresponding index file (`.bai` or `.csi`) should be located in the same directory. An alignment file can be sorted and indexed using the following commands:

```bash
# to sort the alignment file
$ samtools sort -o sorted_output.bam input.bam

# to generate .bai index file
$ samtools index -b sorted_output.bam
```

An alignment file containing at least one of the following tags is preferred for faster processing: `MD` tag, `CS` tag, or a `CIGAR` string with `=/X` operations.

- The CS tag is generated using the --cs option when aligning reads with the [minimap2](https://github.com/lh3/minimap2) aligner. (`--cs=short` is prefered over `--cs=long`)
- The MD tag can be generated using the --MD option in minimap2.

If the alignment files were generated without any of these tags, you can generate the `MD` tag by running the following command to 

```bash
# input: reference genome fasta file & alignment file
# output: an alignment file with MD tag in it

# for generating MD tag
$ samtools calmd -b aln.bam ref.fa > aln_md.bam
```
### Region file
#### `-r or --regions`
**Expects**: *FILE*<br>
**Default**: *None*<br>
The `-r` or `--regions` option is used to specify the input TR regions file. ATaRVa requires a sorted, bgzipped BED file of TR repeat regions, along with its corresponding tabix-indexed file. The BED file should contain the following columns:

1. Chromosome name where TR is located
2. Start position of the TR
3. End position of the TR
4. Repeat motif
5. Motif length
6. Locus ID (optional)

Below is an example of a repeat region BED file. **NOTE: The BED file should either have no header or a header that starts with `#` symbol. The .gz and .tbi files should be in same directory**

| #CHROM | START | END | MOTIF | MOTIF_LEN |
|--------|-------|-----|-------|-----------|
| chr1   | 10000 | 10467 | TAACCC | 6    |
| chr1   | 10481 | 10497 | GCCC | 4      |
| chr2   | 10005 | 10173 | CCCACACACCACA | 13 |
| chr2   | 10174 | 10604 | ACCCTA | 6    |
| chr17  | 60483 | 60491 | AGA    | 3    |

To sort, bgzip, and index the BED file, use the following commands:

#### Sort
```bash
# input: Unsorted bed file
# output: Sorted bed file

# Sorting the BED file using sort
$ sort -k1,1 -k2,2n input.bed > sorted_output.bed
# or using bedtools
$ bedtools sort -i input.bed > sorted_output.bed
```
#### Bgzip
```bash
# input: Sorted bed file
# output: bgzipped bed file

# To keep the original file unchanged and generate separate gz file
$ bgzip -c sorted_output.bed > sorted_output.bed.gz
# or to compress the original file; converts sorted_output.bed to sorted_output.bed.gz
$ bgzip sorted_output.bed
```
#### Index
```bash
# input: bgzipped bed file
# output: tabix indexed file (.tbi)

# install samtools to use tabix
$ tabix -p bed sorted_output.bed.gz
```

## Advanced options
```
Optional arguments:
  -o <FILE>, --vcf <FILE>
                        output VCF file [default: stdout]
  --aln-format <STR>    alignment format [cram | bam | sam] [default: bam]
  --use-softclip        use softclip information [default: False]
  --contigs CONTIGS [CONTIGS ...]
                        contigs to genotype e.g. chr1 chr12 [default: all]
  --karyotype <STR> [<STR> ...]
                        sample karyotypes e.g. XY XX
  --map-qual <INT>      minimum mapping quality [default: 5]
  --min-reads <INT>     minimum read coverage at a locus [default: 10]
  --max-reads <INT>     maximum reads per locus [default: 100; 5000 for amplicon mode]
  --flank <INT>         flank length (bp) to search for insertions [default: 10]
  --snp-dist <INT>      max SNP distance from repeat [default: 3000]
  --snp-count <INT>     number of SNPs for phasing [default: 3]
  --snp-qual <INT>      min base quality at SNP position [default: 20]
  --phased-frac <FLOAT>
                        min fraction of total spanning reads that are phased [default: 0.6]
  --cluster-frac <FLOAT>
                        min fraction of total spanning reads that are phased to one cliuster[default: 0.2, amplicon mode: 0.1]
  --haplotag <STR>      haplotag for phasing e.g. HP [default: None]
  --meth-prob <FLOAT>   min methylation probability [default: 0.5]
  --methviz             write methylation-encoded sequence to VCF [default: False]
  --locus-wise          read-wise genotyping for dense BED regions
  --decompose           write motif-decomposed sequence to VCF
  --amplicon            genotype mode for target-amplified high depth sequencing of pathogenic regions. This mode first uses KDE modeling for genotyping [default: False]
  --strict              strict mode: stick to the coordinates in the reference catalog. If flank is greater than 0 and there are indels in the flanking regions, the locus is flagged for format FV as 1. [Default: False]
  --gvcf                gvcf mode: write gVCF output. Include all loci from the reference catalog in the VCF. [Default: False]
  -t <INT>, --threads <INT>
                        number of threads [default: 1]
  -log, --debug_mode    write debug messages to log file
  -v, --version         show program's version number and exit
```

<div align=center>
  <img src="../lib/workflow.png" alt="Command-driven workflow" width="700"/>
  <p><i>Workflow of Commands </i></p>
</div>

The details of each option are given below:

### `--aln-format`
**Expects**: *STRING*<br>
**Default**: *bam*<br>
This option sets the format of the alignment file. The default format is BAM. Specify the input format as `sam` for SAM files, `cram` for CRAM files, or `bam` for BAM files.

### `--use-softclip`
**Expects**: *BOOLEAN*<br>
**Default**: *False*<br>
This option enables the use of softclip information and supplementary alignments in the reads for genotyping. When set to `True`, ATaRVa will consider the softclipped portions and supplementary alignments of the reads during genotyping.

### `--contigs`
**Expects**: *STRING*<br>
**Default**: *None*<br>
Specify the chromosome(s) for genotyping; repeat loci on all other chromosomes will be skipped. If no chromosomes are mentioned, repeats on all chromosomes in the BED file will be genotyped. eg: `--contigs chr1 chr12 chr22` will genotype only the repeat loci in these mentioned chromosomes in the BED file.

### `--karyotype`
**Expects**: *STRING*<br>
**Default**: *XX*<br>
Karyotype of the samples eg. XX or XY.

### `--map-qual`
**Expects**: *INTEGER*<br>
**Default**: *5*<br>
Minimum mapping quality for the reads to be considered. All reads with a mapping quality below the specified value will be excluded during genotyping.

### `--min-reads`
**Expects**: *INTEGER*<br>
**Default**: *10*<br>
Minimum number of the supporting reads required to genotype a locus. If the number of reads is less than this value, the locus will be skipped.

### `--max-reads`
**Expects**: *INTEGER*<br>
**Default**: *100*<br>
Maximum number of supporting reads allowed for a locus to be genotyped. If the number of reads exceeds this limit, only this specified number of reads will be used for genotyping the locus.

### `--flank`
**Expects**: *INTEGER*<br>
**Default**: *10*<br>
The number of base pairs in the flanking regions to be used for realignment.

### `--snp-dist`
**Expects**: *INTEGER*<br>
**Default**: *3000*<br>
Maximum base pair (bp) distance from the flanks of the repeat locus to fetch SNPs from each read considered for phasing.

### `--snp-count`
**Expects**: *INTEGER*<br>
**Default**: *3*<br>
Maximum number of SNPs to be used for read clustering and phasing.

### `--snp-qual`
**Expects**: *INTEGER*<br>
**Default**: *20*<br>
Minimum Q value of the SNPs to be used for phasing.

### `--phased-frac`
**Expects**: *FLOAT*<br>
**Default**: *0.6*<br>
Minimum fraction of reads required in both clusters relative to the total supporting reads for the repeat locus after phasing.

### `--cluster-frac`
**Expects**: *FLOAT*<br>
**Default**: *0.2*<br>
Minimum fraction reads in a single cluster relative to the total supporting reads for the repeat locus after phasing. If the fraction of reads in a single cluster is below this value, the locus will be skipped.

### `--haplotag`
**Expects**: *STRING*<br>
**Default**: *None*<br>
Specify the haplotype tag to utilize phased information for genotyping. eg `HP`
<div align=center>
  <img src="../lib/haplotag.png" alt="HP clustering" width="500"/>
  <p><i>Haplotag based clustering of reads</i></p>
</div>

### `--meth-prob`
**Expects**: *FLOAT*<br>
**Default**: *0.5*<br>
Minimum probability value of methylation call to be considered for calculation, in the supporting reads of the repeat locus.

### `--methviz`
**Default**: *False*<br>
Calculates site-level methylation levels for CG bases(5mC) within the repeat region and writes them as a base64-encoded MV tag in the VCF SAMPLE column.

### `--locus-wise`

Genotyping mode for BED files with sparse regions across chromosomes or for BED files containing a small number of loci(<500). In this mode, loci are genotyped independently rather than using a read-wise approach.

### `--decompose`
**Default**: *False*<br>
Decomposes called allele sequences into repeat blocks and intervening sequence.
Each reported repeat block represents consecutive copies of a motif;
interruptions can occur between blocks, rather than being counted as part of
a perfect repeat. This helps describe the structure of an allele, not just
its overall length.
<div align=center>
  <img src="../lib/MD.png" alt="Motif decomposition" width="500"/>
  <p><i>Overview of TR motif decomposition </i></p>
</div>
<br>
See [Output format](#output-format) for details.

### `--amplicon`
**Default**: *False*<br>
Genotyping mode for targeted sequencing data. In this mode, the default values for `max-reads` and `flank` values are 1000 and 20 respectively.

### `--strict`
**Default**: *False*<br>
When set to `True`, ATaRVa will not considered any flanking indels for genotyping of the locus.

### `--gvcf`
**Default**: *False*<br>
When set to `True`, ATaRVa output a gVCF file with all repeat loci as in the input catalog, including those that are not genotyped.

### `-t or --threads`
**Expects**: *INTEGER*<br>
**Default**: *1*<br>
Number of threads to use for the process.

### `-log or --debug_mode`
**Expects**: *BOOLEAN*<br>
**Default**: *False*<br>
When set to `True`, ATaRVa will write debug messages to a log file.

### `-v or --version`
Prints the version info of ATaRVa.

## Examples
The following examples assume the input reference genome is in `FASTA` format and is named ref.fa, the alignment file is in `BAM` format and is named input.bam, and the TR regions file is in `BED` format and is named regions.bed.gz.

### Basic usage
To run ATaRVa with default parameters, use the following command:
```bash
$ atarva genotype -f ref.fa --bam input.bam -r regions.bed.gz
```
### With karyotype
To run ATaRVa with sex chromosome karyotype, use the following command:
```bash
$ atarva genotype -f ref.fa --bam input.bam -r regions.bed.gz --karyotype XY
```
With multiple bams:
```bash
$ atarva genotype -f ref.fa --bam input1.bam input2.bam -r regions.bed.gz --karyotype XY XX
```
### With haplotag
To run ATaRVa on haplotagged alignment file, use the following command:
```bash
$ atarva genotype -f ref.fa --bam input.bam -r regions.bed.gz --haplotag HP
```
### With amplicon
To run ATaRVa on targeted sequencing file, use the following command:
```bash
$ atarva genotype -f ref.fa --bam input.bam -r regions.bed.gz --amplicon
```
### Stringent parameter usage
To run ATaRVa with stringent parameters, use the following command:
```bash
$ atarva genotype -q 30 --snp-count 5 --snp-qual 25 --min-reads 20 -t 32 -fi ref.fa --bam input.bam -r regions.bed.gz
# The above command with --snp-count 5 will use a maximum of five heterozygous SNPs to provide accurate genotypes, but only when phasing is based on SNPs and not on length.
```
### Genotyping TRs from specific chromosome/s
To genotype TRs from specific chromosomes only, run ATaRVa with the following command:
```bash
$ atarva genotype --contigs chr9 chr15 chr17 chrX -t 32 -f ref.fa --bam input.bam -r regions.bed.gz
```
### For input alignment file other than `bam`
```bash
# input cram file
$ atarva genotype --aln-format cram -f ref.fa --bam input.cram -r regions.bed.gz

# input sam file
$ atarva genotype --aln-format sam -f ref.fa --bam input.sam -r regions.bed.gz
```
### Usage in docker
To run ATaRVa in docker container, use the following command:
```bash
$ docker run -i -t --rm -v /path_of_necessary_files/:/folder_name atarva:latest genotype -f /folder_name/ref.fa --bam /folder_name/input.bam -r /folder_name/regions.bed.gz
``` 

In all the above examples, the output of ATaRVa is saved to input.vcf unless -o is specified.

## Output format

### `-o or --vcf`
**Expects**: *STRING (to be used as filename)*<br>
**Default**: *Input Alignment Filename + .vcf*<br>
Without `-o`, each input alignment produces a VCF named after the alignment's
filename stem, such as `input.vcf` for `input.bam`. With multiple input
alignments and an explicit output name, the sample stem is appended to that
name for each output.

ATaRVa writes one sample per VCF and declares VCF version 4.5. The
[VCF specification](https://github.com/samtools/hts-specs) describes the
general format. ATaRVa-specific fields are described below.

### Record columns

| Field | Description |
|-------|-------------|
| CHROM | Chromosome containing the repeat region. |
| POS | One-based start position of the repeat region: catalog `START` + 1. |
| ID | Always `.`. A catalog locus ID is instead recorded in INFO `ID`. |
| REF | Reference sequence fetched for the catalog interval. |
| ALT | Full sequence of each called alternate (non-reference) allele. `.` indicates a reference or missing call. `<DEL>` indicates a deletion of the entire repeat locus for that allele. |
| QUAL | Set to `0`, not a calculated genotype quality. |
| FILTER | `PASS` for a successful call; `LESS_READS` for insufficient read depth; or `FAIL` for an unencountered locus. |
| INFO | Locus-level annotations. |
| FORMAT | Colon-separated sample field names. |
| SAMPLE | Colon-separated values corresponding to the FORMAT fields. |

Failed records contain missing genotype values and a read depth when available.
Use `--gvcf` to include ungenotyped catalog loci. Failed records can also
appear without `--gvcf`.

### INFO fields

| INFO | Description |
|------|-------------|
| AC | How many times each ALT allele appears in the sample’s genotype, in the same order as the ALT alleles. |
| AN | Total number of alleles called in the sample’s genotype. |
| MOTIF | Repeat motif from the catalog. |
| START | Zero-based catalog start position. |
| END | Catalog end position, the exclusive end of the BED interval. |
| ID | Optional locus ID from the catalog, or `.` if not provided. |
| REFCN | Reference interval length divided by catalog motif length, rounded down to an integer. |

### FORMAT fields

`.` indicates a missing value. Multiple values are separated by commas.

| FORMAT | Description |
|--------|-------------|
| GT | Called genotype. |
| AL | Length in bases of each called allele. `<DEL>` is assigned length `0`. |
| CN | Allele length divided by the catalog motif length, rounded down to the nearest integer. |
| LPM | Motif of the longest pure repeat segment found in the decomposition. See [motif decomposition](#motif-decomposition-longest-pure-motif-lpm-and-longest-pure-segment-lps). |
| LPS | Longest Pure Segment: number of consecutive, perfect copies of the LPM motif. See [motif decomposition](#motif-decomposition-longest-pure-motif-lpm-and-longest-pure-segment-lps). |
| AR | Lower and upper allele-length bounds for each read cluster, formatted as `lower-upper`. Bounds are calculated from the 2.5th and 97.5th percentiles and rounded. |
| SD | Number of reads supporting each haplogroup or cluster. |
| DP | Number of supporting spanning reads at the locus. |
| SN | Number of informative SNPs used for SNP-based phasing. |
| SQ | SNP quality values used by the phasing implementation. |
| MA | Mean of the per-read methylation averages for reads with methylation information in the haplogroup, rounded to two decimal places. |
| MR | Number of haplogroup reads with methylation information used for MA. |
| DS | Motif-decomposed representation of each allele sequence called in the sample. This can include the reference allele if called. See [motif decomposition](#motif-decomposition-longest-pure-motif-lpm-and-longest-pure-segment-lps). |
| MV | Encoded per-CpG methylation information for visualization, when `--methviz` is enabled and information is available. |
| OL | Hyphen-separated read allele lengths outside the cluster's AR bounds, followed where available by lengths of unphased reads. Groups are comma-separated. |
| FV | `1` when a flanking variant was detected, otherwise `0`. |
| PS | Phase set from input read tags for a phased call. |



### Motif decomposition, Longest Pure Motif (LPM) and Longest Pure Segment (LPS)

With `--decompose`, ATaRVa decomposes called allele sequences into repeat
segments and intervening sequence. A repeat segment is represented as
`(motif)copy_count`, where the count is of consecutive, perfect copies.
The allele as a whole can contain interruptions or different repeat motifs:
these are represented as separate segments, separated by hyphens, rather
than counted within one uninterrupted run. For example, the the sequence ACACACTACAC would be
decomposed as `(AC)3-T-(AC)2` with the interruption `T` between two
blocks of the same motif. Its longest pure segment is `(AC)3`, so `LPM=AC`
and `LPS=3`, not 5. LPS does not combine repeat copies across the
interruption.

The implementation chooses the decomposed repeat segment with the greatest
length in bases (motif length multiplied by copy count), not necessarily the
largest copy count. If equally long segments tie, it prefers one whose motif
is a cyclic rotation of the catalog motif; otherwise it takes the first
longest segment. A selected motif that is a cyclic rotation of the catalog
motif is reported as the catalog motif. If there are no decomposed repeat
segments, the implementation reports the catalog motif and counts its longest
exact contiguous run in the undecomposed sequence, which can be zero.
