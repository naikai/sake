#' use biomaRT to convert human gene to mouse gene 
#'
#' This function allows you to express your love of cats.
#' @param genelist Gene list you want to convert 
#' @keywords standardize
#' @export
#' @examples
#' gene_convert_human_to_mouse("TP53")

gene_convert_human_to_mouse <- function(genelist){
	require(biomaRt)

	human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
	mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
	res <- getLDS(attributes = c("hgnc_symbol"),
					filters = "hgnc_symbol",
					values = genelist,
					mart = human,
					attributesL = c("mgi_symbol"),
					martL = mouse
					)
	unmatched <- genelist[! genelist %in% res[,1] ]
	unmatched <- cbind(unmatched, sapply(tolower(unmatched), simpleCap))
	colnames(unmatched) <- colnames(res)
	res <- rbind(res, unmatched)

	# resort the result to match original order 
	idx <- match(genelist, res[,1])
	return(res[idx, ])
}