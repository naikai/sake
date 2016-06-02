#' Check rownames from data and auto detect which platform it is from
#'
#' Basically we want to convert those affy metrix array data and just
#' leave normal gene_symbols intact
#' @param data Gene list you want to convert
#' @keywords standardize
#' @export
#' @examples
#' detect_genecnt_platform("TP53")

detect_genecnt_ID <- function(data, sum.method="mean"){
	require(magrittr)
	require(pathview)

	# Map molecular data onto Entrez Gene IDs
	id.map.sym2eg <- id2eg(ids = rownames(data), category = "SYMBOL", org="Hs")
	data.entrez <- mol.sum(mol.data = data, id.map = id.map.sym2eg, sum.method = sum.method)

	# deseq2.fc <- gene.entrez[, 1]
	# exp.fc=deseq2.fc
	# head(exp.fc)
	# out.suffix="deseq2"
	return(data.entrez)
}
