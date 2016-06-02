#' Extract data by MAD value 
#'
#' This function allows you to express your love of cats.
#' @param data Original gene expression data 
#' @param topN How many genes from TopMAD list
#' @keywords cats
#' @export
#' @examples
#' cat_function()
extract_data_by_mad <- function (data, topN, by="row", type="data"){
	by.idx <- ifelse(by=="row", 1, 2)
	data.mad.genes <- apply(data, by.idx, mad) %>% 
					  select.top.n(., topN, bottom=F) %>% 
					  names
	if(type=="data"){
		if(by=="row"){
			data <- data[rownames(data) %in% data.mad.genes, ]
		}else if(by=="col"){
			data <- data[ ,colnames(data) %in% data.mad.genes]
		}else{
			stop("Wrong parameters, by can be 'row' or 'col'")
		}
		return(data)
	}else if(type=="genes"){
		return(data.mad.genes)
	}else{
		stop("Wrong parameters, type can be 'data' or 'genes'")
	}
}


