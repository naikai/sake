# Quality Control
Yu-Jui Ho and Toby Aicher, M Hammell Lab  
`r Sys.Date()`  


Following succesful data upload, the gene count expression table must be normalized before proceeding with quality control and analysis. We provide several options for data normalization, variance stabilizing transformations, and quality control (QC). QC analysis is generally aimed at identifying very noisy samples due to technical issues in library construction and/or sample preparation. Two simple data metrics can help identify samples unusually low transcript counts or coverage rates, as described below.

## Normalization 

Normalization allows you to compare read counts between samples and detect differentially expressed genes by accounting for sequencing depth (library size). Three methods to normalize scRNA-seq data are provided.

- **Reads per million mapped (RPM) normalization**: normalize each sample gene count based on the total reads derived from all annotated genes in the library. 

- **DESeq-like normalization**: noramlize the gene count using a method introduced in the [DESeq](http://bioconductor.org/packages/release/bioc/html/DESeq.html) package.

- **Upper quartile normalization**: normalize each sample gene count based on the total counts attributed to the upper quartile of expressed genes in the library.

<img src="Figures/Ting/Normalization.png" width="700px" height="450px" />

## Transformation 

Transformation mitigates the effects of extreme values in scRNA-seq data (from the large dynamic range measured) and adjusts for mean-variance dependency (the observation that genes with higher expression often have larger aboslute variance across samples). Two methods to transform scRNA-Seq data are provided.

- **Variance stablizing transformation (VST)**: transform the gene count using a method intoduced in the DEseq package.

- **Log transformation**: transform the gene count by *log2(count+1)*

<img src="Figures/Ting/Transformation.png" width="700px" height="450px" />

## Filtering 
After data normalization or transformation, there are two simple data metrics to help identify samples that deviate from the majority of samples with respect to total gene counts sampled and total genes covered by at least one read. These problematic samples can be selected and removed before downstream analysis.

### Read Counts Distribution

The read counts distribution summerizes the total number of read counts annotated to known genes for each sample. The read count distribution should be as uniform as possible, and the normalization method you choose will affect this result. Distinctively low read counts may indicate technical issues such as RNA degradation, low amplification rate, or low sequencing efficiency.

### Gene Coverage

Gene coverage summarizes the total number of genes with at least one read in each sample. This number should be relatively stable across libraries. Low gene coverage in a sample can indicate poor quality of a single-cell library. However, the number of expressed genes may also be altered based on biological differences between cell types or experimental conditions, or intrinsic heterogeneity among cell populations. Selection of outlier samples for pre-analysis removal is optional, but can improve downstream analysis.

*** 

An example QC plot displaying read distribution (x-axis) and gene coverage (y-axis) will be used for identifying potential problematic samples. 

<img src="Figures/Ting/QC_plot.png" width="700px" height="450px" />

Samples with relatively low total transcript counts and gene coverage rates usually represent degraded or poorly amplified libraries. These can be identified visually and removed from the sample set before proceeding with downstream analyses. User can use a selection box to highlight samples in the left bottom corner. Samples within the selection box will be shown on the top right table. User can then click on each library identified as problematic in the table at right. Clicking the *"Submit"*  button will discard all highlighted samples from further analysis in the SAKE pipeline. 

In the example case shown above, 3 samples have been selected for removal: `CTC_plt_MP4_13, CTC_plt_MP4_14`, and `CTC_plt_MP4_17` based on very low gene coverage/library size rates (less than 1500 genes sampled and fewer than ~2.2 million reads sequenced). 

*Note: Users have to be cautious when removing samples. Outliers may represent a special cell type with very different transcriptomic patterns relative to all other cells. *


***

Continue on the next section [Gene List Filtering](Genelist_Filtering.Rmd)

