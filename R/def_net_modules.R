#' Identify the the edge of the modules (blocks) within the net connectivity matrix
#' 1. Run for loops to detect the block
#' 2. Detect edge genes in each block and their gene idx
#' 3. Go back to the original 100x100 matrix and use the edge genes (idx) to define block
#'
#' @param data Input expression data
#' @param tao  threshold for hard transformation
#' @param beta parameter for soft power transformation
#' @param num_features Number of top high connective genes to be used
#' @keywords co-expression, network, connectivity
#' @export
#' @examples
#' def_net_modules(data, tao=0.5, beta=10, num_feature=0.01)
# def_net_modules
def_net_modules <- function(core_sub_jaccard_dist, sub_jaccard_dist, min.connectivity=0.2, min.size=5){
   old <- 1
   size <- 1
   blk_cnt <- 0
   # res <- list()
   res <- rep(NULL, 6)
   genenames <- colnames(core_sub_jaccard_dist)

   for (current in 2:dim(core_sub_jaccard_dist)[1]){
      # if (core_sub_jaccard_dist[current, current-1]==1){
      if (core_sub_jaccard_dist[current, current-1]>=min.connectivity){
         size <- size + 1
      }else{
         if(size>=min.size){
            blk_cnt <- blk_cnt + 1
            # Add another filter based on average connectivity of the block
            avg.connectivity <- sum(core_sub_jaccard_dist[old:current-1, old:current-1]) / ((current-1-old+1)*(current-1-old+1-1))
            # res[[blk_cnt]] <- c(old, current-1, size, genenames[old], genenames[current-1], avg.connectivity)
            res <- rbind(res, c(old, current-1, size, genenames[old], genenames[current-1], avg.connectivity))
         }
         size <- 1
         old <- current
      }
   }
   # last check #
   if (size >= min.size){
      blk_cnt <- blk_cnt + 1
      # res[[blk_cnt]] <- c(old, current, size, genenames[old], genenames[current])
      avg.connectivity <- sum(core_sub_jaccard_dist[old:current, old:current]) / ((current-old+1)*(current-old+1-1))
      res <- rbind(res, c(old, current, size, genenames[old], genenames[current], avg.connectivity))
   }

   colnames(res) <- c("Start", "End", "Size", "Start_Gene", "End_Gene", "Avg_Connectivity")
   return(res)
}
