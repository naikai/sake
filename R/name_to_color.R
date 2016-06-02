#' Convert name to color
#'
#' This function convert the sample_name in the data into different colors 
#' @param name usually sample name, separate with "_" 
#' @keywords color
#' @export
#' @examples
#' name_to_color()
# convert sample name to color
name_to_color <- function(name, split_pattern="\\_", num_color=1,
						  ColScheme=c("naikai", "naikai2", "Set1", "Set2", "Set3", "Dark2", "Pastel1", "Pastel2", "Paired", "Accent")
						 )
{
	require(magrittr)

	col.color <- list()
	ColSideColors <- NULL
	col.color.name <- list()
	ColSideColors.name <- NULL

	column.names <- strsplit(name, split=split_pattern)
	max_num_breaks <- max(sapply(column.names, length))

	if(num_color != length(ColScheme)){
		stop(paste0("Error: num_color:", num_color, "is different than num of ColScheme:", length(ColScheme), "\n"))
	}
	if(max_num_breaks > num_color){
		warning("There are more possible breaks in the name than the num_color specified. We will use only the num_color")
	}

	for(i in 1:num_color){
		num_names <- sapply(column.names, function(x) x[i]) %>% unique %>% length
		col.color[[i]] <- create.brewer.color(sapply(column.names, function(x) x[i]), num_names, ColScheme[i])
		ColSideColors <- cbind(ColSideColors, col.color[[i]])
		col.color.name[[i]] <- sapply(column.names, function(x) x[i])
		ColSideColors.name <- cbind(ColSideColors.name, col.color.name[[i]])
	}

	res <- list()
	if(num_color==1){
		if(length(unique(as.character(ColSideColors[,1])))==1){
			res[["color"]] <- NULL
			res[["name"]] <- NULL
		}else{
			res[["color"]] <- as.matrix(ColSideColors[, 1])
			res[["name"]] <- as.matrix(ColSideColors.name[, 1])
		}
	}else if(num_color>1){
		num_color <- ifelse(num_color > max_num_breaks, max_num_breaks, num_color)
		res[["color"]] <- as.matrix(ColSideColors[, 1:num_color])
		res[["name"]] <- as.matrix(ColSideColors.name[, 1:num_color])
	}else{
		stop("Error: num_color must be greater than 0")
	}

	return(res)
}


