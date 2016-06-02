#' Scale data base on specify metrics
#'
#' This function allows you to scale the data based on your metrics of interest
#' @param x data 
#' @param scale By 'row' or 'column' 
#' @param na.rm Remove NAs. Default is TRUE
#' @param method What kind of metrics? Default is 'median'. (Can be 'mean', 'median', 'mode', 'max', 'min', etc)
#' @keywords scale 
#' @export
#' @examples
#' scale_data()
scale_data <- function(x, scale="row", na.rm=TRUE, method="median"){
   require(matrixStats)
   
   if(method=="mean"){
      rowCal <- rowMeans
      colCal <- colMeans
   }else if (method=="median"){
      rowCal <- rowMedians
      colCal <- colMedians
   }

   if (scale=="column" | scale=="col" | scale=="both") {
         x <- sweep(x, 2L, colCal(as.matrix(x), na.rm = na.rm), check.margin = FALSE)
         sx <- apply(x, 2L, sd, na.rm = na.rm)
         x <- sweep(x, 2L, sx, "/", check.margin = FALSE)
   }

   if (scale=="row" | scale=="both") {
         x <- sweep(x, 1L, rowCal(as.matrix(x), na.rm = na.rm), check.margin = FALSE)
         sx <- apply(x, 1L, sd, na.rm = na.rm)
         x <- sweep(x, 1L, sx, "/", check.margin = FALSE)
   }

   if(scale=="none"){
      x <- x 
   }

   return(x)
}
