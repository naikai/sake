#' A Cat Function
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @export
#' @examples
#' cat_function()
plot_gene_expression <- function (data, gene, groups, type="boxplot", data.return=F, ylim=c(0,20), title=NULL, las=1, cex.axis=1){
   data <- as.numeric(data[rownames(data) == gene, ])
      data <- data.frame(groups=factor(groups), Expression=data)
      if (is.null(title)){
         title <- paste("Gene:", gene)
      }

   if (type=="boxplot"){
      p <- ggplot(data, aes(x=groups, y=Expression)) + geom_boxplot(aes(fill = factor(groups)))
         p <- p + ggtitle(title) + theme_bw() + ylim(ylim) + xlab("")
         print(p)
   }else if (type=="base"){
      par(cex.axis=cex.axis)
         boxplot(Expression~groups, data=data, main=title, xlab="Group", ylab="Expression", las=las)
   }
   if(data.return)
      return(p)
}


