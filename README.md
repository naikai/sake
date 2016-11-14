# sake 

[![Build Status](https://travis-ci.com/naikai/sake.svg?token=qigAqQi4xmKjKDqnm97n&branch=master)](https://travis-ci.com/naikai/sake)
[![codecov](https://codecov.io/gh/naikai/sake/branch/master/graph/badge.svg?token=WEipAvcFMf)](https://codecov.io/gh/naikai/sake)

### **S**ingle-cell RNA-Seq **A**nalysis and **K**lustering **E**valuation
The aim of `sake` is to provide a user-friendly tool for easy analysis of NGS Single-Cell transcriptomic data

### Installation Guide

#### First we will install these libraries before installing `sake` 
For Centos 6.9
```
sudo yum install openssl-devel libcurl-devel libpng-devel libxml2-devel libxslt
```

##### Require `gcc` >= 4.6 
```
sudo yum install centos-release-scl
sudo yum install devtoolset-3-toolchain
scl enable devtoolset-3 bash
```


To install the latest development version from GitHub:
```R
source("http://bioconductor.org/biocLite.R")
biocLite(c("annotate", "AnnotationHub", "biomaRt", "DESeq2", "gage", "gageData", "GO.db", "pathview"))

install.packages("devtools")
devtools::install_github("renozao/pkgmaker", ref="develop")
devtools::install_github("naikai/sake", ref="package-installation")
```

### Usage 
```R
library(sake)
shiny::runApp(system.file("sake", package="sake"))
```

### Getting Started
These are the instructions that will walk you through the functions of `sake` package. 

- [Data input](vignettes/Data_Input.Rmd)
- [Test input](DESCRIPTION)
- [plot](Rplot.png)

### Copying & Distribution
Please note that this project is released with a [Contributor Code of Conduct](CONDUCT.md). By participating in this project you agree to abide by its terms.

### Web version 
Feel free to try out the web verison of the tool at url [sake](http://sake.mhammell.tools)
