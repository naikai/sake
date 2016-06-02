#' Select top n percentage of data 
#'
#' Select top n pct from a vector
#' @param data numerical data 
#' @param n percentage to select 
#' @keywords top n percent
#' @export
#' @examples
#' select.n.pct(c(1,23,1,412,51,231,516,1,23,13,3,5,1), 10)
#' [1] 412 516

select.n.pct <- function(data, n){
     data[data > quantile(data, prob=1-n/100)]
}


