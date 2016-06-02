#' Select top n genes 
#'
#' @param data numerical data 
#' @param n number to select
#' @param whole whether to return the whole sets
#' @keywords top n percent
#' @export
#' @examples
#' top_genes(data)

top_genes <- function(data, n=20, whole=F, name="Exp1"){
  data <- as.data.table(data)
  data.summary <- data[, length(num_classes), by="string_gene_symbols"][order(-V1), ]
  setnames(data.summary, c("string_gene_symbols", name))

  if (whole){
    return(data.summary)
  }else{
    return(head(data.summary, n=n))
  }
}