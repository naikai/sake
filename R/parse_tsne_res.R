#' Extract and parse tsne results 
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @export
#' @examples
#' cat_function()
# plot tSNE result 
parse_tsne_res <- function(tsne_out){
   data <- as.data.frame(tsne_out$Y)
   num.samples <- dim(data)[2]
   colnames(data) <- c("x", "y")
   data$Names <- rownames(data)

   return(data)
}


