# icash 

[![Build Status](https://travis-ci.com/naikai/icash.svg?token=qigAqQi4xmKjKDqnm97n&branch=master)](https://travis-ci.com/naikai/icash)
[![codecov](https://codecov.io/gh/naikai/icash/branch/master/graph/badge.svg?token=WEipAvcFMf)](https://codecov.io/gh/naikai/icash)

### **I**nteractive **C**lustering Tool for **A**ssessing **S**ingle-Cell **H**eterogeneity
The aim of `icash` is to provide an user-friendly tools for easy analyzing NGS Single-Cell transcriptomic data

### Version 0.2.1.4

### Installation Guide
To install the latest development version from GitHub:
```R
source("http://bioconductor.org/biocLite.R")
biocLite(c("annotate", "AnnotationHub", "biomaRt", "DESeq2", "gage", "gageData", "GO.db", "pathview"))

install.packages("devtools")
devtools::install_github("naikai/icash", auth_token = "438a6ba6acec1e0a0b3550986f83f42d88b941f9")
```

### Usage 
```R
library(icash)
shiny::runApp(system.file("icash", package="icash"))
```

### Copying & Distribution
Please note that this project is released with a [Contributor Code of Conduct](CONDUCT.md). By participating in this project you agree to abide by its terms.

### Analysis pipeline 

Modify the file
