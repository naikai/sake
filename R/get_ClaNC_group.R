#' A Cat Function
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @export
#' @examples
#' cat_function()
get_ClaNC_group <- function(Build_out){
     ### Extract features for each group ### 
     idx <- apply(Build_out$cntrds, 1, function(x) which.max(abs(x-Mode(x))))
          summary <- sapply(1:max(idx), function(x) names(idx[idx==x]))
            colnames(summary) <- paste0("Group", 1:ncol(summary))
              return(summary)
}


