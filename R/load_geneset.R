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