#' Adjust label cex based on its number 
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @export
#' @examples
#' assign_label_cex(100)

assign_label_cex <- function(n.samples){
   if (n.samples < 40){
      label.cex <- 1.2
   }else if(n.samples < 70){
      label.cex <- 1.0
   }else if(n.samples < 100){
      label.cex <- 0.7
   }else if (n.samples < 200){
      label.cex <- 0.6
   }else if (n.samples < 300){
      label.cex <- 0.5
   }else{
      label.cex <- 0.4
   }

   return(label.cex)
}


