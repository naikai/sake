#' Loop through each of the modules and connect the adjacent ones which is only allow.gap away from each other
#'
#' @param final_summary summary table
#' @param final_genelist summary gene list
#' @param allow.gap allow how man gaps between each sub-modules
#' @keywords co-expression, network, connectivity
#' @export
#' @examples
#' connect_gap_modules(final_summary, final_genelist)
# connect_gap_modules

plot_modules <- function(core_sub_jaccard_dist, final_genelist, allow.gap=1){

      if(sum(which.idx)>0){
         final_summary <- rbind(final_summary, res_summary[which.idx, ])
         final_genelist <- c(final_genelist, unlist(res_genelist[which.idx]))

         new_core_sub_jaccard_dist <- matrix(0, ncol=dim(core_sub_jaccard_dist)[1], nrow=dim(core_sub_jaccard_dist)[1], dimnames = dimnames(core_sub_jaccard_dist))
         for (j in which.idx){
            idx <- res_summary[j,1]:res_summary[j,2]
            new_core_sub_jaccard_dist[idx,idx] = core_sub_jaccard_dist[idx, idx]
            ### Can use bit vector calculation for faster speed? ###

            # save genelist for each modules
            filename <- paste0(paste(res_summary[j,c(4,5,7,8)], collapse="-"), ".txt")
            write.table(data.frame(Gene=unlist(res_genelist[j])), filename, quote=F, row.names=F)
         }
         bin_new_core_sub_jaccard_dist <- replace(new_core_sub_jaccard_dist, new_core_sub_jaccard_dist>0, 1)
         myHeatmap.3(bin_new_core_sub_jaccard_dist, type = "heatmap", Colv=NULL, Rowv=NULL, scale="none")
         myHeatmap.3(new_core_sub_jaccard_dist, type = "heatmap", Colv=NULL, Rowv=NULL, scale="none")
      }else{
          plot(1, main="No res fits the criteria", col="red")
      }
}