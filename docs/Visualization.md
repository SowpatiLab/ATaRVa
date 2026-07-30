## Visualization
For visualizing motif structure and methylation levels of TR alleles, use [**VisuaMiTRa**](https://github.com/SowpatiLab/visuamitra.git), an auxiliary visualization tool that runs in a web browser and accepts ATaRVa VCF files as input. 

To include motif decomposition and methylation information in the VCF, run ATaRVa with the `--decompose` and `--methviz` flags. Motif decomposition can still be generated from the default VCF even if `--decompose` is not specified; however, the `--methviz` flag is required for visualizing methylation levels.

VisuaMiTRa allows users to compare motif structure and methylation levels within the same allele in single-sample mode and also supports multi-sample analysis.

<div align="center">
  <img src="/lib/vis_overview.png" alt="Visualization of TRs" width="800">
  <p><i>Overview of visuamitra</i></p>
</div>