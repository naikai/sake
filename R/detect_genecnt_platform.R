#' Check rownames from data and auto detect which platform it is from
#'
#' Basically we want to convert those affy metrix array data and just
#' leave normal gene_symbols intact
#' @param data Gene list you want to convert
#' @keywords standardize
#' @export
#' @examples
#' detect_genecnt_platform("TP53")

detect_genecnt_platform <- function(data, method="mean"){
	require("annotate")
	require(magrittr)

	PROBES <- as.character(rownames(data))
	num.genes <- nrow(data)

	# for now just use simple nrow as check point
	if(num.genes == 22283){
  	  require(hgu133a.db) # 22283 probes
	  GPL="GPL96"
	  print ("Guess it's hgu133a")
	}else if(num.genes == 22645){
  	  require(hgu133b.db) # 22645 probes
	  GPL="GPL97"
	  print ("Guess it's hgu133b")
	}else if(num.genes == 54675 || num.genes == 44137){
  	  require(hgu133plus2.db) # 54675 probes
  	  # hard fix here for now, need to add mapped rownames to all these probe.annotation and then decide which one 
	  GPL="GPL570"
	  print ("Guess it's hgu133plus2")
	}else if(num.genes == 22277){
  	  require(hgu133a2.db) 
	  GPL="GPL571"
	  print ("Guess it's hgu133a2")
	}else{
	  print ("Guess its normal RNA-Seq, do nothing")
	  return(data)
	}

	if(method == "mean"){
	  select.fun <- function(x) rowMeans(x)
	}else if(method == "sd"){
	  select.fun <- function(x) rowSds(x)
	}else if(method == "max"){
	  select.fun <- function(x) rowMax(x)
	}else if(method == "min"){
	  select.fun <- function(x) rowMin(x)
	}else if(method == "iqr"){
	  select.fun <- function(x) apply(x, 1, IQR)
	}
	# There are two levels of mutliple mapping
	# one is that each probe_ids may have multiple mapped gene_symbols
	# For now, we just select the first gene_symbols if there are multiples
	gene.symbols <- affy_probe_to_gene_symbol(PROBES, GPL) %>%
	                filter(!duplicated(PROBEID)) %>%
	                dplyr::select(SYMBOL) %>%
	                unlist

	# Another place is that there may be multiple probes designed for each gene
	# Now we select the probes with the max(mean expression across samples)
	# Mayb fix it later by adding more different filter selections? Ray. 2015-10-30
	data.gene <- data %>% as.data.table %>%
	                      '['(, RowExp := select.fun(.SD)) %>%
	                      '['(, SYMBOL := gene.symbols) %>% group_by(., SYMBOL) %>%
	                      filter(., which.max(RowExp))

	# Remove rows without gene_symbols (NAs), then convert back to data.frame
	final.data <- data.gene %>%
	                  filter(!is.na(SYMBOL)) %>%
	                  setDF %>%
	                  set_rownames(.$SYMBOL) %>%
	                  '['(, !(colnames(.) %in% c("RowExp", "SYMBOL")))

	return(final.data)
}
