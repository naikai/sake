<!-- README.md is generated from README.Rmd. Please edit that file -->
sake
====

[![Build Status](https://travis-ci.com/naikai/sake.svg?token=qigAqQi4xmKjKDqnm97n&branch=master)](https://travis-ci.com/naikai/sake) [![codecov](https://codecov.io/gh/naikai/sake/branch/master/graph/badge.svg?token=WEipAvcFMf)](https://codecov.io/gh/naikai/sake)

**S**ingle-cell RNA-Seq **A**nalysis and **K**lustering **E**valuation
----------------------------------------------------------------------

The aim of `sake` is to provide a user-friendly tool for easy analysis of NGS Single-Cell transcriptomic data

<img src="vignettes/Figures/SAKE_workflow.png" width="1024px" height="647px" />

**Flowchart of SAKE package and example analysis results**: **a)** Analysis workflow for analyzing single-cell RNA-Seq data. **b)** Quality Controls to compare total sequenced reads to total gene transcripts detected. **c)** Sample correlation heat map plot **d)** A heat map of sample assignment from NMF run, with dark red indicating high confidence in cluster assignments **e)** t-SNE plot to compare NMF assigned groups with t-SNE projections. **f)** A table of NMF identified features (genes defining each cluster) and a box plot of gene expression distributions across NMF assigned groups. **g)** Summary table for GO term enrichment analysis for each NMF assigned group.

Installation
------------

First we will install some prerequisite libraries before installing `sake`

For **Centos** (tested on 6.9)

``` bash
sudo yum install openssl-devel libcurl-devel libpng-devel libxml2-devel libxslt

# Require `gcc` >= 4.6 
sudo yum install centos-release-scl
sudo yum install devtoolset-3-toolchain
scl enable devtoolset-3 bash
```

For **Ubuntu** (tested on 16.10)

``` bash
sudo apt-get install libcurl4-openssl-dev libpng-dev libxslt-dev libssl-dev libxml2-dev xsltproc
```

For **Mac OS** (tested on Sierra 10.12.1)

``` bash
# follow instructions to install brew on MAC
http://brew.sh
# install required packages 
brew install curl openssl libpng libxslt libxml2 
# update gcc
brew install gcc48
```

#### Download and [Install R](http://cran.wustl.edu)

-   Download and [Install Rstudio](https://www.rstudio.com/products/rstudio/download/) (Suggested but not required)

#### To install the latest development version from GitHub (around 30 minutes)

``` r
source("http://bioconductor.org/biocLite.R")
biocLite(c("annotate", "AnnotationHub", "biomaRt", "DESeq2", "gage", "gageData", "GO.db", "pathview", "plotly", "DT"))

install.packages("devtools")
devtools::install_github("renozao/pkgmaker", ref="develop")
devtools::install_github("naikai/sake") # You may want to consider "devtools::install_github("naikai/sake", CC=gcc-7) to flag for use with the GCC compiler" 
```

#### Quick update on how to install SAKE. 2023-02-08 (inputs from @savytskanatalia)

# I reinstalled older versions of the packages: 
# plotly 4.9.2,  htmltools 0.5.0, htmlwidgets 1.5.2, 
# promises 1.1.1, crosstalk 1.1.0.1, and DT 0.16 
remove.packages("plotly")
remove.packages("htmltools")
remove.packages("htmlwidgets")
remove.packages("promises")
remove.packages("crosstalk")
remove.packages("DT")
install.packages("https://cran.r-project.org/src/contrib/Archive/crosstalk/crosstalk_1.1.0.1.tar.gz", repo=NULL, type="source")
install.packages("https://cran.r-project.org/src/contrib/Archive/promises/promises_1.1.1.tar.gz", repo=NULL, type="source")

install.packages("https://cran.r-project.org/src/contrib/Archive/plotly/plotly_4.9.2.tar.gz", repo=NULL, type="source")
install.packages("https://cran.r-project.org/src/contrib/Archive/htmltools/htmltools_0.5.0.tar.gz", repo=NULL, type="source")
install.packages("https://cran.r-project.org/src/contrib/Archive/htmlwidgets/htmlwidgets_1.5.2.tar.gz", repo=NULL, type="source")
install.packages("https://cran.r-project.org/src/contrib/Archive/DT/DT_0.16.tar.gz", repo=NULL, type="source")
# gage_2.36.0

# go even below the versions
# shiny 1.4.0.2
# shinydashboard 0.7.1
# shinythemes 1.1.2

remove.packages("shiny")
remove.packages("shinydashboard")
remove.packages("shinythemes")

install.packages("https://cran.r-project.org/src/contrib/Archive/shiny/shiny_1.4.0.2.tar.gz", repo=NULL, type="source")
install.packages("https://cran.r-project.org/src/contrib/Archive/shinydashboard/shinydashboard_0.7.1.tar.gz", repo=NULL, type="source")
install.packages("https://cran.r-project.org/src/contrib/Archive/shinythemes/shinythemes_1.1.2.tar.gz", repo=NULL, type="source")




```
> sessionInfo()
R version 4.0.4 (2021-02-15)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19042)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
[4] LC_NUMERIC=C                           LC_TIME=English_United States.1252    

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods   base    

other attached packages:
[1] sake_0.4.0          dplyr_1.0.5         NMF_0.30.1          Biobase_2.50.0      BiocGenerics_0.36.0 cluster_2.1.0      
[7] rngtools_1.5        pkgmaker_0.32.2.900 registry_0.5-1    

loaded via a namespace (and not attached):
 [1] viridis_0.5.1       pkgload_1.2.0       viridisLite_0.3.0   foreach_1.5.1       shiny_1.5.0         assertthat_0.2.1  
 [7] BiocManager_1.30.12 remotes_2.2.0       sessioninfo_1.1.1   pillar_1.5.1        glue_1.4.2          digest_0.6.27      
[13] RColorBrewer_1.1-2  promises_1.1.1      colorspace_2.0-0    htmltools_0.5.1.1   httpuv_1.5.5        plyr_1.8.6        
[19] pkgconfig_2.0.3     devtools_2.3.2      purrr_0.3.4         xtable_1.8-4        scales_1.1.1        processx_3.5.0    
[25] later_1.1.0.1       tibble_3.1.0        generics_0.1.0      ggplot2_3.3.3       usethis_2.0.1       ellipsis_0.3.1    
[31] cachem_1.0.4        withr_2.4.1         cli_2.3.1           magrittr_2.0.1      crayon_1.4.1        mime_0.10          
[37] memoise_2.0.0       ps_1.6.0            fs_1.5.0            fansi_0.4.2         doParallel_1.0.16   pkgbuild_1.2.0    
[43] tools_4.0.4         data.table_1.14.0   prettyunits_1.1.1   lifecycle_1.0.0     matrixStats_0.58.0  gridBase_0.4-7    
[49] stringr_1.4.0       munsell_0.5.0       callr_3.6.0         compiler_4.0.4      rlang_0.4.10        grid_4.0.4        
[55] iterators_1.0.13    rstudioapi_0.13     testthat_3.0.2      gtable_0.3.0        codetools_0.2-18    DBI_1.1.1          
[61] curl_4.3            reshape2_1.4.4      R6_2.5.0            gridExtra_2.3       fastmap_1.1.0       utf8_1.2.1        
[67] rprojroot_2.0.2     dendextend_1.14.0   desc_1.3.0          stringi_1.5.3       Rcpp_1.0.6          vctrs_0.3.6        
[73] tidyselect_1.1.0  
```


Usage
-----

``` r
library(sake)
shiny::runApp(system.file("sake", package="sake"))
```

Getting Started
---------------

Please follow the links to briefly walk you through the functions of `sake` package.

-   [Data input](vignettes/Data_Input.Rmd)
-   [Quality control](vignettes/Quality_Control.Rmd)
-   [Filtering](vignettes/Filtering.Rmd)
-   [Run NMF](vignettes/NMF.Rmd)
-   [Visualization](vignettes/Visualization.Rmd)
-   [DE and Enrichment](vignettes/DE_Enrich.Rmd)

-   Example1 - Ting et al [HTML](vignettes/Ting.Rmd) [PDF](vignettes/Ting_pdf.pdf)
-   Example2 - Treutlein et al [HTML](vignettes/Treutlein.Rmd) [PDF](vignettes/Treutlein_pdf.pdf)

Copying & Distribution
----------------------

Please note that this project is released with a [Contributor Code of Conduct](CONDUCT.md). By participating in this project you agree to abide by its terms.

