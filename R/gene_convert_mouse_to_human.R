#' use biomaRT to convert mouse gene to human gene 
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords standardize
#' @export
#' @examples
#' standardize(x)

gene_convert_mouse_to_human <- function(genelist){
	require(biomaRt)

	human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
	mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
	res <- getLDS(attributes = c("mgi_symbol"),
					filters = "mgi_symbol",
					values = genelist,
					mart = mouse,
					attributesL = c("hgnc_symbol"),
					martL = human
					)
	unmatched <- genelist[! genelist %in% res[,1] ]
	unmatched <- cbind(unmatched, toupper(unmatched))
	colnames(unmatched) <- colnames(res)
	res <- rbind(res, unmatched)

	# resort the result to match original order 
	idx <- match(genelist, res[,1])
	return(res[idx, ])
}