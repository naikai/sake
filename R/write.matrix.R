#' Own version of writing matrix output 
#' @param countdata integer data frame
#' @keywords write.matrix 
#' @export
#' @examples
#' my.write.matrix(data)

my.write.matrix <- function(x, file = "", sep = "\t", col.names=T,
   append=F,
   row.names=F,
   justify = c( "none", "left", "right"),
   pval=NULL,
   newline="\n",
   names=NULL ) 
{
   justify = match.arg( justify )
   x <- as.matrix(x)
   p <- ncol(x)
   cnames <- colnames(x)
   rnames <- rownames(x)

   if ( !is.null(pval) ) {
      x[,pval] <- format.pval( as.numeric(x[,pval]) )
   }
   if ( col.names && !is.null(cnames) )
   {
      x <- format(rbind( cnames, x ), justify=justify)
   }
   if ( row.names )
   {
      p <- p+1
      if ( col.names && !is.null(cnames) ) {
         rnames <- if (is.null(names)) c("",rnames) else c(names, rnames)
      }
      x <- cbind( format(rnames,justify=justify), x )
   }
   cat( t(x), file=file, sep=c(rep(sep,p - 1), newline), append=append )
}
