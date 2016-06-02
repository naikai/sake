#' Select top number from a vector
#'
#' Select top # from a vector
#' @param data numerical data
#' @param n top number
#' @param bottom whether to select from bottom?
#' @keywords select 
#' @export
#' @examples
#' select.top.n(data, 100, bottom=F)

select.top.n <- function(data, n, bottom=F){
   if(bottom){
      data <- head(sort(data, decreasing = F), n=n)
   }else{
      data <- head(sort(data, decreasing = T), n=n)
   }
   return(data)
}

