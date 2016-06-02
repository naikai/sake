#' Identify the informative genes in Data with top high connectivity in the co-expression network
#'
#' Code is adapted and rewrote based on Wange et al, BMC Bioinformatics, 2014 
#' 'Improving the sensitivity of sample clustering by leveraging gene co-expression networks in variable selection'
#'
#' @param data Input expression data 
#' @param tao  threshold for hard transformation 
#' @param beta parameter for soft power transformation 
#' @param num_features Number of top high connective genes to be used 
#' @keywords co-expression, network, connectivity 
#' @export
#' @examples
#' NetConnectivity(data, tao=0.5, beta=10, num_feature=0.01)
# NetConnectivity  
NetConnectivity <- function(data, thresh='soft', tao=0.7, beta=1, method="pearson", diag.zero=T){
   # require(HiClimR)
   require(magrittr)

   # Use fastCor function to compute the correlation matrix
   # if(method=="Pearson"){
   #    system.time(xcor0 <- t(data) %>% fastCor %>% abs)
   # }else if (method=="Spearman"){
   #    system.time(xcor0 <- t(data) %>% cor(., method="spearman") %>% abs)
   # }
   # Use cor in WGCNA 
   system.time(xcor0 <- t(data) %>% WGCNA::cor(., method=method) %>% abs)

   if(diag.zero){
      system.time(diag(xcor0) <- 0) #for crossprod later
   }

   # net connectivity can be calculate using two transformation
   # 1. hard threshold, define a tao (signum function)
   # 2. soft threahold, define a beta (power function)
   if(thresh=="hard"){
      net_dist <- (xcor0 > tao) * 1
   }else if(thresh=="soft"){
      net_dist <- xcor0^beta
   }

   return(net_dist)
}
