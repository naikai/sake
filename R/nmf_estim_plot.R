#' generate individual plot for estimating K
#'
#' This function allows you to express your love of cats.
#' @param res NMF run result
#' @keywords silhouette
#' @export
#' @examples
#' nmf_silhouette_plot(res)

nmf_estim_plot <- function(estim.r){
	require(NMF)

	nmf_rank <- estim.r$measures[, 1]
	for(i in 2:ncol(estim.r$measures)){
		ylabel <- colnames(estim.r$measures)[i]
  		# plot(nmf_rank, estim.r$measures[, i], type="o", xlab="Rank", ylab=ylabel)
  		a <- ggplot(estim.r$measures, aes_string(x="rank", y=ylabel)) + theme_bw() + geom_point(size=4) + geom_line() +
  				theme(axis.text=element_text(size=15), axis.title=element_text(size=15)) + 
  				ggtitle(ylabel) + theme(plot.title = element_text(lineheight=.8, size=15, face="bold"))
  		print(a)
  	}
}

