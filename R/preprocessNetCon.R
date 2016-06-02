#' Preprocess the connecitivty matrix based on mean connectivity and binary switch 
#' 1. First try to remove columns (rows) that have median expression 75% quantile
#' 2. Convert the matrix into binary based on quantile connectivity'
#'
#' @param data Input expression data 
#' @param tao  threshold for hard transformation 
#' @param beta parameter for soft power transformation 
#' @param num_features Number of top high connective genes to be used 
#' @keywords co-expression, network, connectivity 
#' @export
#' @examples
#' preprocessNetCon(sub_jaccard_dist, rmv.filter=0.5, binary.filter=0.75, plot=F)
# preprocessNetCon  
preprocessNetCon <- function(sub_jaccard_dist, rmv.filter=0.5, binary.filter=0.75, plot=F){
   # filter low mean connectivity genes
   # keep.idx <- colMedians(sub_jaccard_dist) > median(sub_jaccard_dist)
   keep.idx <- colSums(sub_jaccard_dist) >= median(colSums(sub_jaccard_dist))
   core_sub_jaccard_dist <- sub_jaccard_dist[keep.idx, keep.idx]
   if(plot){
      myHeatmap.3(core_sub_jaccard_dist, type = "heatmap", Colv=NULL, Rowv=NULL, scale="none", cexRow=0.5, cexCol=0.5)
   }

   # binary switch filter
   quantile(core_sub_jaccard_dist)
   top75.idx <- core_sub_jaccard_dist > quantile(core_sub_jaccard_dist)[4]
   bin_core_sub_jaccard_dist <- core_sub_jaccard_dist
   bin_core_sub_jaccard_dist[top75.idx] = 1
   if(plot){
      myHeatmap.3(bin_core_sub_jaccard_dist, type = "heatmap", Colv=NULL, Rowv=NULL, scale="none", cexRow=0.5, cexCol=0.5)
   }

   return(core_sub_jaccard_dist)
}
