#' Manually specify coloring for provided groups 
#'
#' Assign colors to data, if more unique data than provided colors, will impute the missing ones and fill them in
#' @param data data
#' @param group.color colors for each factorial group
#' @keywords color manual
#' @export
#' @examples
#' create.manual.color(c(1,2,3,1,2,3,2,3), c("red", "blue", "green"))

create.manual.color <- function(data, group.color)
{
	num_uniq_data <- length(unique(data))
	if(num_uniq_data > length(group.color)){
		require(RColorBrewer)
		group.color <- colorRampPalette(group.color)(num_uniq_data)
	}
	data <- as.numeric(factor(data))
	return(group.color[data])
}


