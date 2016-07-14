# icash 

[![Build Status](https://travis-ci.com/naikai/icash.svg?token=qigAqQi4xmKjKDqnm97n&branch=master)](https://travis-ci.com/naikai/icash)
[![codecov](https://codecov.io/gh/naikai/icash/branch/master/graph/badge.svg?token=WEipAvcFMf)](https://codecov.io/gh/naikai/icash)

### **S**ingle-cell RNA-Seq **A**nalysis and **K**lustering **E**valuation
The aim of `sake` is to provide a user-friendly tool for easy analysis of NGS Single-Cell transcriptomic data

### Version 0.2.1.4

### Installation Guide
To install the latest development version from GitHub:
```R
source("http://bioconductor.org/biocLite.R")
biocLite(c("annotate", "AnnotationHub", "biomaRt", "DESeq2", "gage", "gageData", "GO.db", "pathview"))

install.packages("devtools")
devtools::install_github("naikai/sake", auth_token = "438a6ba6acec1e0a0b3550986f83f42d88b941f9")
```

### Usage 
```R
library(sake)
shiny::runApp(system.file("sake", package="sake"))
```

### Copying & Distribution
Please note that this project is released with a [Contributor Code of Conduct](CONDUCT.md). By participating in this project you agree to abide by its terms.

### Analysis pipeline 
