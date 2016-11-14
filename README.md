<!-- README.md is generated from README.Rmd. Please edit that file -->
sake
====

[![Build Status](https://travis-ci.com/naikai/sake.svg?token=qigAqQi4xmKjKDqnm97n&branch=master)](https://travis-ci.com/naikai/sake) [![codecov](https://codecov.io/gh/naikai/sake/branch/master/graph/badge.svg?token=WEipAvcFMf)](https://codecov.io/gh/naikai/sake)

### **S**ingle-cell RNA-Seq **A**nalysis and **K**lustering **E**valuation

The aim of `sake` is to provide a user-friendly tool for easy analysis of NGS Single-Cell transcriptomic data

### Installation Guide

First we will install some prerequisite libraries before installing `sake`

For **Centos** (tested on 6.9)

``` bash
sudo yum install openssl-devel libcurl-devel libpng-devel libxml2-devel libxslt

# Require `gcc` >= 4.6 
sudo yum install centos-release-scl
sudo yum install devtoolset-3-toolchain
scl enable devtoolset-3 bash
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

Download and [Install R](http://cran.wustl.edu)
Download and [Install RStudio](https://www.rstudio.com/products/rstudio/download/) (Suggested but not Required)

#### To install the latest development version from GitHub (takes around 30-40 minutes)

``` r
source("http://bioconductor.org/biocLite.R")
biocLite(c("annotate", "AnnotationHub", "biomaRt", "DESeq2", "gage", "gageData", "GO.db", "pathview"))

install.packages("devtools")
devtools::install_github("renozao/pkgmaker", ref="develop")
devtools::install_github("naikai/sake", ref="package-installation")
```

#### Or you can install package through [packrat](https://rstudio.github.io/packrat/) (takes around 7-10 minutes)

### Usage

``` r
library(sake)
shiny::runApp(system.file("sake", package="sake"))
```

### Getting Started

These are the instructions that will walk you through the functions of `sake` package.

-   [Data input](vignettes/Data_Input.Rmd)
-   [Data Metrics](vignettes/Data_Metrics.Rmd)
-   [Filtering](vignettes/Filtering.Rmd)
-   [Run NMF](vignettes/NMF.Rmd)
-   [Visualization](vignettes/Visualization.Rmd)

### Copying & Distribution

Please note that this project is released with a [Contributor Code of Conduct](CONDUCT.md). By participating in this project you agree to abide by its terms.

### Web version

Feel free to try out the web verison of the tool at url [sake](http://sake.mhammell.tools)
