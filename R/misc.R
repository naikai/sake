#' load predefined gene list files for heatmap
#'
#' Load predefined gene list files for heatmap (clustering) result
#' @param folder Where the files are
#' @param pattern file extensions for the ones you want to load
#' @keywords lodad
#' @export
#' @examples
#' load_genelist("data", ".txt")

load_genelist <- function(folder, pattern=".txt"){
  genelist <- list()
  filename <- dir(file.path(folder), pattern)
  id <- sapply(strsplit(filename, "\\."), function(x) x[2])
  for(i in 1:length(filename)){
    genelist[[id[i]]] <- read.table(file.path(folder, filename[i]), header=T, sep="\t", stringsAsFactors=FALSE)
  }
  return(genelist)
}


#' load predefined gene set files for GSEA
#'
#' Load predefined gene list files for gene set enrichment analysis (GSEA)
#' @param folder Where the files are
#' @param pattern file extensions for the ones you want to load
#' @keywords lodad
#' @export
#' @examples
#' load_geneset("data", ".txt")
load_geneset <- function(folder, pattern=".gmt"){
  geneset <- list()
  filename <- dir(file.path(folder), pattern)
  id <- sapply(strsplit(filename, "\\."), function(x) paste(x[1], x[2], sep="."))
  for(i in 1:length(filename)){
    geneset[[id[i]]] <- file.path(folder, filename[i])
  }
  return(geneset)
}



#' Run DEseq2
#'
#' This function allows you to run DESeq2 if provided required data
#' @keywords DESeq2
#' @export
#' @examples
#' run.DESeq2()
run.DESeq2 <- function(sub.data, sub.expinfo, design, workers=16)
{
  register(BiocParallel::MulticoreParam(workers))
  colData <- sub.expinfo
  contrast.var <- as.character(design)[2]
  ddsfeatureCounts <- DESeqDataSetFromMatrix(countData=sub.data, colData=colData, design=design)
  ### Set betaPrior=FALSE to go with MLE LFC to get simple LFC = (avg in group2/ avg in group1)
  dds <- DESeq(ddsfeatureCounts, parallel=T, betaPrior=FALSE)
  return(dds)
}


#' Wrapper function to run the icahs Shiny App
#' @export
run_icash <- function()
{
  shiny::runApp(system.file('icash', package='icash'))
}

