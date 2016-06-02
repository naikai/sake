#' A Cat Function
#'
#' This function allows you to run DESeq2 if provided required data 
#' @keywords DESeq DESeq2
#' @export
#' @examples
#' run.DESeq2()

run.DESeq2 <- function(sub.data, sub.expinfo, design, workers=16)
{
   register(MulticoreParam(workers))
      colData <- sub.expinfo
      contrast.var <- as.character(design)[2] 
      ddsfeatureCounts <- DESeqDataSetFromMatrix(countData=sub.data, colData=colData, design=design)
### Set betaPrior=FALSE to go with MLE LFC to get simple LFC = (avg in group2/ avg in group1)
      dds <- DESeq(ddsfeatureCounts, parallel=T, betaPrior=FALSE)
      return(dds)
}


