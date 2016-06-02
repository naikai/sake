#' Transform the countdata into Reads per million reads 
#'
#' Transform the countdata into Reads per million reads 
#' @param data countdata
#' @param mapped_reads mapped reads in each library
#' @keywords RPM
#' @export
#' @examples
#' RPM()

RPM <- function(data, mapped_reads){
     data <- t(t(data)/mapped_reads*1000000) 
}


