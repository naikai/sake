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