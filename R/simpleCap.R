#' Capitalize the first letter 
#'
#' This function allows you to express your love of cats.
#' @param genelist Gene list you want to convert 
#' @keywords standardize
#' @export
#' @examples
#' simpleCap("meat")

simpleCap <- function(x) {
   s <- strsplit(x, " ")[[1]]
   paste(toupper(substring(s, 1,1)), substring(s, 2), sep="", collapse=" ")
}
