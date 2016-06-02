#' create centroid based on the grouping info 
#'
#' This function allows you to express your love of cats.
#' @param data Input data
#' @param groups Groups for each samples 
#' @keywords centroid reshape
#' @export
#' @examples
#' create_centroid()
create_centroid <- function(data, groups, scale="none"){
   require(reshape2)
   require(reshape)
      cdata <- scale.data(data, scale=scale)
      cdata$gene <- rownames(cdata)
      ddata <- melt(cdata, id="gene")
      ddata$group <- rep(groups, each=nrow(data))
      centroid <- cast(ddata, group ~ gene, mean)
      return(t(centroid))
}

