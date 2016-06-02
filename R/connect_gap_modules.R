#' Loop through each of the modules and connect the adjacent ones which is only allow.gap away from each other
#'
#' @param final_summary summary table
#' @param final_genelist summary gene list
#' @param allow.gap allow how man gaps between each sub-modules
#' @keywords co-expression, network, connectivity
#' @export
#' @examples
#' connect_gap_modules(final_summary, final_genelist, allow.gap=1)
# connect_gap_modules

connect_gap_modules <- function(final_summary, final_genelist, allow.gap=1){
   # create a column 'Delete' to specify which rows to be removed after we check through all the rows
   final_summary$Delete <- 0
   old <- 1

   for(current in 2:nrow(final_summary)){
      if( (final_summary[current, "Original_Start"] - final_summary[old, "Original_End"]) <= allow.gap ){
        # print("Ready to merged")
        # print(paste("Current Starts", final_summary[current, "Original_Start"], "Old ends", final_summary[old, "Original_End"]))

        # update these metrics for final_summary
        final_summary[current, "Size"] <- final_summary[current, "End"] -final_summary[old, "Start"] + 1
        new.avg.connectivity <- (final_summary[old, "Avg_Connectivity"] * final_summary[old, "Size"] + final_summary[current, "Avg_Connectivity"] * final_summary[current, "Size"]) / final_summary[current, "Size"]
        final_summary[current, "Avg_Connectivity"] <- new.avg.connectivity
        final_summary[current, "Start"] <- final_summary[old, "Start"]
        final_summary[current, "Start_Gene"] <- final_summary[old, "Start_Gene"]
        final_summary[current, "Original_Start"] <- final_summary[old, "Original_Start"]
        final_summary[old, "Delete"] <- 1

        ### update final_genelist
        final_genelist[[current]] <- unlist(c(final_genelist[[old]], final_genelist[[current]]))
      }
      old <- current
   }
   ### Old keep the ones that is not be Deleted, and then remove 'Delete' column
   idx <- final_summary$Delete == 1
   final_summary <- subset(final_summary, Delete != 1, select = -c(Delete))
   final_genelist[idx] <- NULL

   res <- list()
   res[['summary']] <- final_summary
   res[['genelist']] <- final_genelist
   return(res)
}