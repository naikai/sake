#' This function allows you to VST your data 
#' @param countdata integer data frame
#' @keywords vst
#' @export
#' @examples
#' vst(countdata)

vst <- function(countdata){
   library(DESeq)
      condition <- factor(rep("Tumour", ncol(countdata)))
      countdata <- newCountDataSet(countdata,condition )
      countdata <- estimateSizeFactors( countdata )
      cdsBlind <- DESeq::estimateDispersions( countdata, method="blind")
      vstdata <- varianceStabilizingTransformation( cdsBlind )
      return(exprs(vstdata))
}
