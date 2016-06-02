#' use hgu133a.db or hgu133plus2.db to convert from affymetrix probe_id to human gene_symbol
#'
#' Conversion between affy probe_ID to gene_symbol
#' @param genelist Gene list you want to convert
#' @keywords standardize
#' @export
#' @examples
#' gene_convert_human_to_mouse("TP53")

affy_probe_to_gene_symbol <- function(probe_id, GPL, sel.columns="SYMBOL"){
	require("annotate")
	PROBES<- as.character(probe_id)

	if(GPL==96 || GPL=="GPL96"){
		require(hgu133a.db)
		OUT <- select(hgu133a.db, PROBES, sel.columns)
	}else if(GPL==97 || GPL=="GPL97"){
		require(hgu133b.db)
		OUT <- select(hgu133b.db, PROBES, sel.columns)
	}else if(GPL==570 || GPL=="GPL570"){
		require(hgu133plus2.db)
		OUT <- select(hgu133plus2.db, PROBES, sel.columns)
	}else if(GPL==571 || GPL=="GPL571"){
		require(hgu133a2.db)
		OUT <- select(hgu133a2.db, PROBES, sel.columns)
	}

	# for now, just select the first gene_symbols if there are multiples
	# we can fix it later by comparing the expression values and select the best one? Ray. 2015-10-30
	# OUT2 %>% group_by(PROBEID) %>% summarise(paste(SYMBOL, collapse=","))

	return(OUT)
}
