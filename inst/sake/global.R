library(sake)
library(networkD3)

options(shiny.maxRequestSize=1000*1024^2)
options(shiny.jquery.version = 1)
raw_fd <- system.file('extdata', 'rawdata', package="sake")
filenames<-list.files(path=raw_fd)
# filenames<-list.files(path="./extdata", pattern="\\.txt$")
geneset_fd <- system.file("extdata", "geneset", package="sake")
genefile_fd <- system.file("extdata", "genefile", package="sake")
pre_geneset <- load_geneset(geneset_fd, ".gmt")
pre_genelist <- load_genelist(genefile_fd, ".txt")
# `%then%` <- shiny:::`%OR%`
`%then%` <- function(a, b) {
  if (is.null(a)) b else a
}
