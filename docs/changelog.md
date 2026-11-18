## Changelog

### v0.9.5
* Fixed missing entries bug in merged VCF file. Added length-based merging option in VCF merging.
* Added the LPM and LPS tags in format field of the output VCF file to match with v0.7.1. LPM - Longest Pure Motif, LPS - Longest Pure Copy Number.
* Removed locus mode. Splitting the catalog to exclude large regions without repeats consequently skipping reads in those regions.
* Added a `--gvcf` option to output a gVCF file with all repeat loci as in the input catalog.
* Filtering out outliers using HDBSCAN in non-amplicon mode before generating the consensus sequence.

### v0.9.4
* Fixed read contribution thresholds for each haplogroup to match with ATaRVA v0.7.1.
* Introduce `ref_snps` in read_data structure and optimized `record_ref_snps` function to reduce the time taken.
* Clearing out supplementary alignments data as the reads are deleted.
* Optimized Q value phasing using numpy array operations to reduce the time taken.
* Threshold for sequence lengths considered for abPOA consensus sequence generation. Limit of 10kb on median length and 5kb difference between the longest and median sequences.

### v0.9.3
* New VCF merging strategy implemented
* Removed the `--rna` and `--instability` options.

### v0.9.2
* Changed `--phasing-read` to `--phased-frac` to specify the minimum fraction of phased reads required for haplotyping. The default value is set to 0.6.
* Changed `--snp-read` to `--cluster-frac` to specify the minimum fraction of reads required for clustering. The default value is set to 0.2.
* Removed processing of softclipped and supplementary alignments as defaul. Users can enable this feature using the `--use-softclip` option.
* Included the processing of supplementary alignments along with softclipped portions of the reads.
* Reintroduced `--locus-wise` and `--amplicon` options for genotyping mode selection.
* Debugged the merge functionality.

### v0.9.1
* Includes a `--strict` mode to ignore the variants in the flank and flag a locus with flanking variations using FV.
* Inlcuded a option to skip using softclipped portions of the reads for genotyping using `--skip-softclip` option.
* Include phase set number for haplotagged BAMs

### v0.7.1+ext0.01
* Reports PS (phase set) in the VCF for haplotagged BAM inputs.
* Fixed - Report SQ value as '.' for loci with no informative SNPs.

### v0.7.1+ext
* First release. Extended version of ATaRVa with additional features as below.
* Default processing of the softclipped portions of the read to search for the presence of repeat loci from the catalogs by searching for the flanking sequences.
* Implements a slightly different version of read-haplogrouping based on the informative SNPs. The reads are heirarchically grouped based on the minimising the Qvalue penalty of the informative SNPs.
* Added a new option `--rna` to enable genotyping of tandem repeats from RNA-seq data. This option is designed parse the CIGAR string of the typicaly RNA-seq alignment which include introns and other splicing events encoded as N.
* Added a new option `--instability` which outputs an additional file with allele information from each read at a locus
* Changed `--format` to `--aln-format` to specify the input format for alignment file.

### v0.7.1
* Fixed NumPy version dependency issues across Python versions
* Updated the KDE clustering score calculation formula by incorporating a skewness-based parameter

### v0.7.0
* Replaced K-means clustering in `--amplicon` mode and non-SNP regions with KDE-based clustering and edit-distance–based HDBSCAN.
* Forced all loci into `haplotyping` mode to ensure clustering is performed, even when most reads support a single allele.
* Added `LPM` tag in the VCF_SAMPLE column to report longest pure repeat motif and its copy number.
* Updated ALT allele assignment to prioritize sequence comparison against the reference sequence instead of relying solely on allele length.
* Improved consensus sequence generation by modifying sequence ordering:
  * WGS mode now orders sequences based on the mode of allele lengths.
  * `--amplicon` mode orders sequences alternately from both sides of the median allele length.

### v0.6.0
* Fixed an incorrect code modification that caused `--amplicon` mode to produce incorrect results in previous version(V0.5.0)
* Fixed bugs in `loci-wise` mode related to storing SNP info
* Introduced subcommands to separate operating modes [`genotype` & `merge`]
* Added `MV` tag in the VCF_SAMPLE column to report base-wise methylation level for visualization purpose
* Changed the name of `MM` tag into `MA`
* Improvised the mean methylation level calculation in `MA` tag
* Increased the default `--snp-qual` from 13 to 20
* Added `CN` and `REFCN` tags to the VCF to report motif copy number in the sample and INFO fields, respectively.

### v0.5.0
* Changed the VCF-START column into 1-based coordinate system
* Included `START` tag in VCF-INFO column with 0-based coordinate system
* Added `MR` tags in the VCF_SAMPLE column to report the supporting read count for mean methylation level
* Added a confirmation step to check `MM` extraction from the reverse strand
* Added a `loci-wise` flag to perform region-wise genotyping (instead of the default read-wise mode) for BED files with sparse regions
* Improved Motif-decomposition script to maintain consistent representation of a motif (cyclic variation check)

### v0.4.0
* Added `MM` tag in VCF_SAMPLE column for mean methylation level
* Modified `AR` tag in VCF-SAMPLE column with central 95% allele range
* Implemented DBSCAN clustering in `amplicon` mode to check for multiple clusters
* Fixed bugs in decomposition function [#8](https://github.com/SowpatiLab/ATaRVa/issues/8)

### v0.3.1
* Added checkpoint in amplicon mode for non-repeatedness in ALT sequence
* Refined Motif-decomposition sequence for motif breaks
* Added `AR` tag in VCF-SAMPLE column for allele range

### v0.3.0
* Added `--amplicon` mode for targeted sequencing data
* Added function to convert eqx read sequence
* Improved Outlier cleaning in K-Means clustering
* Implemented De-novo motif identification in motif-decomposition
* Added optional tag `ID` in INFO field if BED input has additional column

### v0.2.0
* Added `--haplotag` argument to enable the use of haplotag information for genotyping.
* Fixed bugs in SNP-based clustering.
* Replaced the use of the mode function with a consensus-based approach for final allele derivation.
* Removed `PC` tag from the FORMAT field of the output VCF.

### v0.1.2
* Modified input arguments.

### v0.1.1
* Added a Mac OS compatible <code>.so</code> file.

### v0.1
* First release.
