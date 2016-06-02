#' Remove genes that with too many constant 0 expression across samples 
#'
#' This function allows you to remove those data entries with constant variance
#' You can specify the percentage of 0s as the threshold 
#' @param data Input data set
#' @param by Default by "row", can change to "column"
#' @keywords rmv_constant_var
#' @export
#' @examples
#' rmv_constant_var(data, by="row")
#' Remove data entries with constant variance 

rmv_constant_0 <- function(data, by="row", pct=0.75, minimum=0){
   alt_by=""
   if(by=="row"){
      num_all_var <- dim(data)[1]
      data <- data[apply(data, 1, function(x) sum(x<=minimum)<=length(x)*pct), ]
      num_rmv_var <- dim(data)[1]
      alt_by="col"
   }else if(by=="col"){
      num_all_var <- dim(data)[2]
      data <- data[, apply(data, 2, function(x) sum(x<=minimum)<=length(x)*pct)]
      num_rmv_var <- dim(data)[2]
      alt_by="row"
   }
   print(sprintf ("Original data: %d %s, Removed %d because these %s have values below or equal to %s in more than %d percent of all %s", num_all_var, by, num_all_var-num_rmv_var, by, minimum, pct*100, alt_by))

   return(data)
}
