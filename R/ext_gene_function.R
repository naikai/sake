#' Extract ensembl_id using gene symbol 
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @export
#' @examples
#' cat_function()
# can be extened by providing the "columan name" that we want to match in snpinfo #
ext_gene_function <- function(genelist, ensembl, attributes=c('ensembl_gene_id', 'hgnc_symbol', 'description', 'gene_biotype')){
    # xli <- head(genelist, n=30)
    # res <- queryMany(genelist, scopes='symbol', fields=c('entrezgene', 'go'), species='human')
    # ensembl <- ensembl()
    # filters = listFilters(ensembl)
    # attributes <- listAttributes(ensembl)
    genelist <- as.data.table(genelist)
    setnames(genelist, "hgnc_symbol")
    my_ensembl_gene_id <- getBM(attributes=attributes,
                                filters = 'hgnc_symbol',
                                values = genelist,
                                mart = ensembl)
    # colnames(my_ensembl_gene_id) <- c("ensembl_gene_id", "string_gene_symbols", "description", "gene_biotype")
    my_ensembl_gene_id <- as.data.table(my_ensembl_gene_id[grepl("^ENSG", my_ensembl_gene_id$ensembl_gene_id), ])
    merged <- merge(genelist, my_ensembl_gene_id, by='hgnc_symbol', all=T, allow.cartesian=TRUE)

    merged$GeneCard <- paste0("<a href='http://www.genecards.org/cgi-bin/carddisp.pl?gene=",
                             merged$hgnc_symbol, "'>", merged$hgnc_symbol, "</a>")
    merged$Ensembl <- ifelse(is.na(merged$ensembl_gene_id), NA,
                      paste("<a href='http://www.ensembl.org/Homo_sapiens/Gene/Summary?g=",
                      merged$ensembl_gene_id, "'>", merged$ensembl_gene_id, "</a>", sep=""))
    # Drop unwanted columns
    merged$ensembl_gene_id <- NULL
    # merged$hgnc_symbol <- NULL
    # Reorder the column order
    # merged <- as.data.frame(merged)[c(3,1,2,4)]
    return(merged)
}
