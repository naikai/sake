#' use hgu133a.db or hgu133plus2.db to convert from affymetrix probe_id to human gene_symbol
#'
#' Conversion between affy probe_ID to gene_symbol
#' @param genelist Gene list you want to convert
#' @keywords standardize
#' @export
#' @examples
#' affy_probe_to_gene_symbol("TP53")
affy_probe_to_gene_symbol <- function(probe_id, GPL, sel.columns="SYMBOL"){
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
  return(OUT)
}


#' Automatically convert genecnt into Entrez ID
#'
#' Standize gene id into Entrez ID
#' @param data Gene id list you want to convert
#' @keywords conversion
#' @export
#' @examples
#' detect_genecnt_id("TP53")
detect_genecnt_id <- function(data, sum.method="mean"){
  # Map molecular data onto Entrez Gene IDs
  id.map.sym2eg <- id2eg(ids = rownames(data), category = "SYMBOL", org="Hs")
  data.entrez <- mol.sum(mol.data = data, id.map = id.map.sym2eg, sum.method = sum.method)
  # deseq2.fc <- gene.entrez[, 1]
  # exp.fc=deseq2.fc
  # head(exp.fc)
  # out.suffix="deseq2"
  return(data.entrez)
}

#' Check rownames from data and auto detect which platform it is from
#'
#' Basically we want to convert those affy metrix array data and just
#' leave normal gene_symbols intact
#' @param data Gene list you want to convert
#' @keywords standardize
#' @import data.table
#' @export
#' @examples
#' detect_genecnt_platform("TP53")
detect_genecnt_platform <- function(data, method="mean"){
  PROBES <- as.character(rownames(data))
  num.genes <- nrow(data)
  # for now just use simple nrow as check point
  if(num.genes == 22283){
    library(hgu133a.db) # 22283 probes
    GPL="GPL96"
    print ("Guess it's hgu133a")
  }else if(num.genes == 22645){
    library(hgu133b.db) # 22645 probes
    GPL="GPL97"
    print ("Guess it's hgu133b")
  }else if(num.genes == 54675 || num.genes == 44137){
    library(hgu133plus2.db) # 54675 probes
    # hard fix here for now, need to add mapped rownames to all these probe.annotation and then decide which one
    GPL="GPL570"
    print ("Guess it's hgu133plus2")
  }else if(num.genes == 22277){
    library(hgu133a2.db)
    GPL="GPL571"
    print ("Guess it's hgu133a2")
  }else{
    print ("Guess its normal RNA-Seq, do nothing")
    return(data)
  }

  if(method == "mean"){
    select.fun <- function(x) rowMeans(x)
  }else if(method == "sd"){
    select.fun <- function(x) apply(x, 1, sd)
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
                  dplyr::filter(!duplicated(PROBEID)) %>%
                  dplyr::select(SYMBOL) %>%
                  unlist

  # Another place is that there may be multiple probes designed for each gene
  # Now we select the probes with the max(mean expression across samples)
  # Mayb fix it later by adding more different filter selections? Ray. 2015-10-30
  data.gene <- data %>% as.data.table() %>%
                        '['(, RowExp := select.fun(.SD)) %>%
                        '['(, SYMBOL := gene.symbols) %>%
                        dplyr::group_by(SYMBOL) %>%
                        dplyr::filter(RowExp == max(RowExp))

  # Remove rows without gene_symbols (NAs), then convert back to data.frame
  final.data <- data.gene %>%
                dplyr::filter(!is.na(SYMBOL)) %>%
                as.data.frame() %>% 
                setDF %>%
                magrittr::set_rownames(.$SYMBOL) %>%
                '['(, !(colnames(.) %in% c("RowExp", "SYMBOL")))

  return(final.data)
}


#' Extract ensembl_id using gene symbol
#'
#' Can be extened by providing the "columan name" that we want to match in snpinfo #
#' @param genelist input genelist
#' @param ensembl where to extract emsembl information
#' @keywords ensembl, annotate
#' @export
#' @examples
#' ext_gene_function()
ext_gene_function <- function(genelist, ensembl, attributes=c('ensembl_gene_id', 'hgnc_symbol', 'description', 'gene_biotype')){
    genelist <- as.data.table(genelist)
    setnames(genelist, "hgnc_symbol")
    my_ensembl_gene_id <- getBM(attributes=attributes,
                                filters = 'hgnc_symbol',
                                values = genelist,
                                mart = ensembl)
    my_ensembl_gene_id <- as.data.table(my_ensembl_gene_id[grepl("^ENSG", my_ensembl_gene_id$ensembl_gene_id), ])
    merged <- merge(genelist, my_ensembl_gene_id, by='hgnc_symbol', all=T, allow.cartesian=TRUE)

    merged$GeneCard <- paste0("<a href='http://www.genecards.org/cgi-bin/carddisp.pl?gene=",
                             merged$hgnc_symbol, "'>", merged$hgnc_symbol, "</a>")
    merged$Ensembl <- ifelse(is.na(merged$ensembl_gene_id), NA,
                      paste("<a href='http://www.ensembl.org/Homo_sapiens/Gene/Summary?g=",
                      merged$ensembl_gene_id, "'>", merged$ensembl_gene_id, "</a>", sep=""))
    # Drop unwanted columns
    merged$ensembl_gene_id <- NULL
    # Reorder the column order
    # merged <- as.data.frame(merged)[c(3,1,2,4)]
    return(merged)
}


#' Use biomaRT to convert human gene to mouse gene
#'
#' For the human gene symbols that we can not find corresponding mouse symbols,
#' we will simply convert the first letter into Capital and return it
#' @param genelist Gene list you want to convert
#' @keywords standardize
#' @export
#' @examples
#' gene_convert_human_to_mouse("TP53")
gene_convert_human_to_mouse <- function(genelist){
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

#' use biomaRT to convert mouse gene to human gene
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords standardize
#' @export
#' @examples
#' gene_convert_human_to_mouse("Trp53")
gene_convert_mouse_to_human <- function(genelist){
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
