#' generate silhouette plot 
#'
#' This function allows you to express your love of cats.
#' @param res NMF run result
#' @keywords silhouette
#' @export
#' @examples
#' nmf_silhouette_plot(res)

nmf_silhouette_plot <- function(res, type="consensus", silorder=F){
	require(NMF)
	si <- silhouette(res, what=type)
	plot(si)
}
