#' Identify the the edge of the modules (blocks) within the net connectivity matrix
#' Run for loops to detect the block
#' Detect edge genes in each block and their gene idx
#' Go back to the original 100x100 matrix and use the edge genes (idx) to define block
#'
#' Extract groups (sample clustering) from NMF run results 
#' X = W x H
#' W is the feature matrix, H is the sample matrix
#' After NMF run, use this function to select import features and assign groups for the two matrix

#' @param core_sub_jaccard_dist 
#' @keywords NMF feature
#' @export
#' @examples
#' summarize_module(core_sub_jaccard_dist, sub_jaccard_dist, min.connectivity=0.2, min.size=5)
# summarize_module
summarize_module2 <- function(core_sub_jaccard_dist, old, current, blk_cnt, genenames, original_names, res_summary, res_genelist, size, min.size=5, start_idx){
   if(size>=min.size){
      blk_cnt <- blk_cnt + 1
      # Add another filter based on average connectivity of the block
      avg.connectivity <- sum(core_sub_jaccard_dist[old:(current-1), old:(current-1)]) / ((current-1-old+1)*(current-1-old+1-1))
      # extract outer gene idx from the original 100x100 matrix #
      original_idx <- match(c(genenames[old], genenames[current-1]), original_names) + start_idx - 1
      res_summary <<- rbind(res_summary, c(old, current-1, size, genenames[old], genenames[current-1], avg.connectivity, original_idx[1], original_idx[2]))
      res_genelist[[blk_cnt]] <<- genenames[old:(current-1)]
   }
}
