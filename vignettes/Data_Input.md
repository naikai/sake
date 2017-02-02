# Data Input
Yu-Jui Ho and Toby Aicher, M Hammell Lab  
`r Sys.Date()`  

There are three different ways you can upload files into `sake` pacakge. It can either be 

- New Gene Count Data File 
- Pre-loaded Demo Data (for testing)
- Saved run result (for advanced users)

## New Gene Count Data file

The gene count data file is one of the most common formats for RNA-seq assays. Each row represents the expression value for a gene (in raw counts or normalized RPM, TPM, etc.), and each column represents a single sample. The package requires the first row to be the header, containing unique names for each sample. The first column is required to be the names/ID for each gene/transcript. 

While the gene count file is expected to be tab-demlited, you can specify other characters used to separate fields.

- Tab (\\t) - Default setting for `.txt` or `.out` file
- Comma (,) - Usually used by `.csv` file
- Semicolon (;) - Less common

Example data sets should look like this

Gene    |MEF-1    |MEF-10   |MEF-11   |MEF-12   |MEF-2
--------|---------|---------|---------|---------|-----
Gm15772 |1493.562 |1714.470 |1178.217 |1858.733 |1199.904
Dnajc3  |75.209   |67.320   |291.554  |49.924   |166.867
Mdn1    |29.288   |7.819    |82.620   |1.262    |0.214
Mfap1b  |4.796    |1.335    |4.308    |0.000    |0.748
Zglp1   |1.939    |78.381   |3.385    |0.541    |3.205
Gm12359 |1.225    |13.159   |1.846    |0.000    |0.320
Gm16039 |0.408    |2.861    |0.154    |0.360    |56.406
Gm11149 |0.204    |0.000    |0.000    |0.000    |0.000

## Pre-loaded Demo Data

There are several pre-loaded gene expression datasets from published single-cell studies available for learning how to use the SAKE package. These include  one study exploring neuronal differentiation over a time-course^[Treutlein *et al*, Dissecting direct reprogramming from fibroblast to neuron using single-cell RNA-seq, Nature, 2016] as well as a second study evaluating circulating tumor cells in a pancreatic cancer mouse model^[Ting *et al*, Single-Cell RNA Sequencing Identifies Extracellular Matrix Gene Expression by Pancreatic Circulating Tumor Cells, Cell Reports, 2014]. These datasets allow the user to reproduce the analysis results presented in in the SAKE paper.


An example screenshot shows selection of the pre-loaded data set downloaded from GEO and published in Ting *et al* (2014). 

<img src="Figures/Ting/preload.png" width="800px" height="450px" />


A successfully loaded data will look like this  

<img src="Figures/Ting/loaded_file.png" width="800px" height="450px" />

***



## Saved run result

Users familiar with SAKE have the option to run most of the computationally intensive portions of the SAKE clustering  algorithm on their own clustered compute servers, then upload these results to our web host for interactive analysis of the results. This is especially useful when the sample sizes of the single-cell study become too large for this web host to analyze in real time (greater than ~200 cells). The previous run results can be saved in [.rda](https://stat.ethz.ch/R-manual/R-devel/library/base/html/load.html) format and then loaded back to `sake`. Saved run results will include data analysis of NMF, t-SNE, and DESeq2 (if specified). Uploading these results to the `sake` server allows for interactive figure generation. 

Example saved data can be downloaded [here](../data/Ting.NMF5.rda)

*** 

Continue on the next section [Quality Control](Quality_Control.Rmd)

*** 

