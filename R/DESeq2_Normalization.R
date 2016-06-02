#' This function allows you to VST your data 
#' @param countdata integer data frame
#' @keywords vst
#' @export
#' @examples
#' vst(countdata)

DESeq_Normalization <- function(countdata){
   	require(DESeq2)
   	condition <- factor(rep("Sample", ncol(countdata)))
   	countdata <- newCountDataSet(countdata,condition )
   	countdata <- estimateSizeFactors( countdata )
   	normalized.countdata <- counts(countdata, normalized=TRUE)
   	return(exprs(normalized.countdata))
}
