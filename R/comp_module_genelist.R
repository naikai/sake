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

comp_module_genelist <- function(final_comp_genelist, folderName0){
   pdf(paste0("Overlap", folderName0, ".pdf"), height=10, width=10)

   num <- length(final_comp_genelist) 
   for(i in 2:num){
      comb <- combn(num, i)
      apply(comb, 2, function(x) {
         sub_final_comp_genelist <- final_comp_genelist[x]
         venn(sub_final_comp_genelist)
      })
   }
   dev.off()
}