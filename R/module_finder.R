#' Identify the modules in Data within co-expression network
#'
#' Specify the top percentile of correlation to be used 
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @export
#' @examples
#' cat_function()
# plot tSNE result 
module_finder <- function(data, p.value=0.05, tao=0.5, beta=3, num_features=100){
   require(HiClimR)
   require(magrittr)

   num_all_var <- dim(data)[1]
   # data <- rmv_constant_var(data)
   data <- rmv_constant_0(data)
   num_rmv_var <- dim(data)[1]
   print(sprintf ("Original %d transcripts, removed %d because they have constant variation", num_all_var, num_all_var-num_rmv_var))

   # Use fastCor function to compute the correlation matrix
   xcor <- t(data) %>% fastCor %>% abs 

   # net connectivity can be calculate using two transformation 
   # 1. hard threshold, define a tao (signum function)
   # 2. soft threahold, define a beta (power function)
   if(thresh=="hard"){
     net_dist <- xcor > tao 
   }else if(thresh=="soft"){
     net_dist <- xcor^beta 
   }

   sum_connectivity <- colSums(net_dist)
   # selected_features <- order(sum_connectivity, decreasing=T) %>% 
   #                      head(., n=num_features) %>% 
   #                      colnames(net_dist)[.]
   



   return(rawdata)
}
