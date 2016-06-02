#' Transform the countdata by upper quantile value from each sample
#'
#' Transform the countdata using the quantile value from each sample
#' @param data countdata
#' @keywords UQ, quantile, normalization
#' @export
#' @examples
#' RPM()

# upper quartile normalization, add option to remove all zeros in the sample first 
UQ <- function(data, remove.zero=T){
	require(magrittr)
	if(remove.zero){
		upper_quartile <- apply(data, 2, function(x) quantile(x[x!=0])[4])
	}else{
		upper_quartile <- apply(data, 2, function(x) quantile(x)[4])
	}
	data.upperQ <- data %>% t %>% divide_by(upper_quartile) %>% t %>% multiply_by(1000)
	return(data.upperQ)
}


