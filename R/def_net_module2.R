#' Identify the the edge of the modules (blocks) within the net connectivity matrix
#'
#' Run for loops to detect the block
#' Detect edge genes in each block and their gene idx
#' Go back to the original 100x100 matrix and use the edge genes (idx) to define block

#' @param core_sub_jaccard_dist preprocessed sub_jaccard_dist matrix (after applying filter, binarized the data)
#' @param sub_jaccard_dist zoomed-in small block size of jaccard_dist
#' @param min.connectivity threshold for minimum connectivity (set to 75\% quartile in each block)   
#' @param min.size threshold for minimum block size 
#' @param start_idx specifiy what is the index for the starting gene in the block 
#' @param allow.gap allow how man gaps between each sub-modules 
#' @keywords co-expression, network, connectivity
#' @export
#' @examples
#' def_net_modules2(core_sub_jaccard_dist, sub_jaccard_dist, min.connectivity=0.2, min.size=5)
def_net_module2 <- function(core_sub_jaccard_dist, sub_jaccard_dist, min.connectivity=0.2, min.size=5, start_idx=1, allow.gap=1){

   summarize_module <- function(){
      if(size>=min.size){
         blk_cnt <<- blk_cnt + 1
         # Add another filter based on average connectivity of the block
         avg.connectivity <- sum(core_sub_jaccard_dist[old:(current-1), old:(current-1)]) / ((current-1-old+1)*(current-1-old+1-1))
         # extract outer gene idx from the original 100x100 matrix #
         original_idx <- match(c(genenames[old], genenames[current-1]), original_names) + start_idx - 1
         res_summary <<- rbind(res_summary, c(old, current-1, size, genenames[old], genenames[current-1], avg.connectivity, original_idx[1], original_idx[2]))
         res_genelist[[blk_cnt]] <<- genenames[old:(current-1)]
      }
   }

   old <- 1
   size <- 1
   blk_cnt <- 0
   current_gap <- 0 # parameter to track how many gaps in between modules 
   res <- list()
   res_genelist <- list()
   res_summary <- rep(NULL, 8)
   genenames <- colnames(core_sub_jaccard_dist)
   original_names <- colnames(sub_jaccard_dist)

   for (current in 2:dim(core_sub_jaccard_dist)[1]){
   	# Allow for 1 (default) gap between each sub-modules 
      if (core_sub_jaccard_dist[current, current-1]>=min.connectivity){
         size <- size + 1 + current_gap
         current_gap <- 0
      }else{
         # current_gap = current_gap + 1 
         # if(current_gap <= allow.gap){
         #    # which means allowing this gap, and there are other connecting nodes within current modules 
         #    if(any(core_sub_jaccard_dist[current, old:current] >= min.connectivity)){
         #       next;
         #    }else{
         #       # summarize_module(core_sub_jaccard_dist, old, current, blk_cnt, genenames, original_names, res_summary, res_genelist, size, min.size=min.size, start_idx=start_idx)
         #       summarize_module()
         #    }
         # }else{
         #    # summarize_module(core_sub_jaccard_dist, old, current, blk_cnt, genenames, original_names, res_summary, res_genelist, size, min.size=min.size, start_idx=start_idx)
         #    summarize_module()
         # }
         summarize_module()
         size <- 1 
         old <- current 
         current_gap <- 0 
      }
   }
   # last check #
   summarize_module()
   # summarize_module(core_sub_jaccard_dist, old, current, blk_cnt, genenames, original_names, res_summary, res_genelist, size, min.size=min.size, start_idx=start_idx)

   if(!is.null(res_summary)){
      colnames(res_summary) <- c("Start", "End", "Size", "Start_Gene", "End_Gene", "Avg_Connectivity", "Original_Start", "Original_End")
   }

   res[['summary']] <- res_summary
   res[['genelist']] <- res_genelist
   return(res)
}
