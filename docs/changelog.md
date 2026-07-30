## Changelog
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