# Quality Control
Yu-Jui Ho and Toby Aicher  
`r Sys.Date()`  


You can normalize or transform your raw data after succesfully uploading it. We provide two simple data metrics to help identify samples that deviate from the majority of samples, which may indicate potential low-quality samples that can be removed before downstream analysis.

## Normalization 

Normalization allows you to compare read counts between samples and detect differentially expressed genes by accounting for technical variability between samples and adjusting for their sequencing depth (library size). Two methods to normalize scRNA-seq data are provided.

- **Reads per millions (RPM) normalization**: normalize the gene count from each sample according to their total read count. 

- **DESeq-like normalization**: noramlize the gene count using a method introduced in the [DESeq](http://bioconductor.org/packages/release/bioc/html/DESeq.html) package.

- **Upper quartile normalization**: noramlize the gene count based on upper quartile value from each library

<img src="Figures/Ting/Normalization.png" width="700px" height="450px" />

## Transformation 

Transformation allows you to take into account the presence of extreme values in scRNA-seq data and to adjust for mean variance dependency (the observation that genes with higher expression often have larger variances across samples). Two methods to transform scRNA-Seq data are provided.

- **Variance stablizing transformation (VST)**: transform the gene count using a method intoduced in the DEseq package.

- **Log transformation**: transform the gene count by log2(count+1)

<img src="Figures/Ting/Transformation.png" width="700px" height="450px" />

## Filtering 
After data normalization or transformation, there are two simple data metrics to help identify samples that deviate from the majority of samples, which can be selected and removed before downstream analysis.

### Read distribution

Read distribution summerizes the total number of read counts for each sample. The read count distribution should be as uniform as possible, and the normalization method you choose will affect the result. Distinctively low read counts may indicate RNA degradation or low sequencing efficiency.

### Gene coverage

Gene coverage summarizes the total number of genes with at least one read in each sample. This number should be relatively stable across libraries. Low gene coverage in a sample can indicate poor quality of a single-cell library. However, the number of expressed genes may be altered based on the difference between cell types, experimental conditions or sequencing protocals, intrinsic heterogeneity among cell population, among other reasons. User should be aware of these factors and decide how to treat outlier samples for downstream analysis.

*** 

An example QC plot displaying read distribution (x-axis) and gene coverage (y-axis) will be used for identifying potential problematic samples. 

<img src="Figures/Ting/QC_plot.png" width="700px" height="450px" />

Samples with relatively low total transcript counts and gene coverage rates usually represent degraded or poorly amplified libraries. These can be identified visually and removed from the sample set before proceeding with downstream analyses.
User will be able to use a selection box to highlight samples in the left bottom corner. Samples within the selection box will be shown on the top right table. User can then inspect each of them and click on the ones they want to remove from downstream analysis. Selected samples will be displayed on the bottom right, and user can hit `Submit` button when they are ready to discard these samples. 

In this case, we will just demostrate how you can select and remove 3 samples `CTC_plt_MP4_13, CTC_plt_MP4_14`, and `CTC_plt_MP4_17` just for demonstration purpose. We actually did not remove any samples for downstream analysis in order to compare our clustering results to the published results. 

*Note: Users have to be cautious when removing samples. Sometimes they may represent a special cell type which has very different transcriptomic pattern than all the other cells. *


***

Continue on the next section [Filtering](Filtering.Rmd)

