#' A Cat Function
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @export
#' @examples
#' cat_function()

create.brewer.color <- function(data, num=8, name="Set1")
{
	if(name=="naikai"){
		my_pallete <- c(
                        rgb(236,67,35, maxColorValue=255),
                        rgb(56,146,208, maxColorValue=255),
                        rgb(230,235,88, maxColorValue=255),
                        rgb(116,187,88, maxColorValue=255),
                        rgb(196,95,46, maxColorValue=255),
                        rgb(203,77,202, maxColorValue=255),
                        rgb(118,220,65, maxColorValue=255),
                        rgb(115,113,206, maxColorValue=255))
		return(create.manual.color(data, my_pallete))
	}else if(name=="naikai2"){
		my_pallete <- c(
                        rgb(10,10,10, maxColorValue=255),
                        rgb(230,235,88, maxColorValue=255),
                        rgb(118,220,65, maxColorValue=255),
                        rgb(203,77,202, maxColorValue=255),
                        rgb(196,95,46, maxColorValue=255),
                        rgb(116,187,88, maxColorValue=255),
                        rgb(236,67,35, maxColorValue=255),
                        rgb(56,146,208, maxColorValue=255),
                        rgb(115,113,206, maxColorValue=255))
		return(create.manual.color(data, my_pallete))
	}else{
		groupCodes <- as.factor(data)
		colorCodes <- colorRampPalette(brewer.pal(num, name))(num)
		color.idx <- match(groupCodes, levels(groupCodes))
		label.color <- colorCodes[color.idx]
		return(label.color)
	}
}


