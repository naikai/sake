#' A Cat Function
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @export
#' @examples
#' cat_function()
predict_from_centroid <- function(centroids, data, scale="none", method="spearman"){
# match centroids with data. Added Aug 18, 2015
   idx <- match(rownames(centroids), rownames(data))
      centroids <- centroids[!is.na(idx), ] 

      groups <- rep(NULL, ncol(data))
      if (dim(data)[2] > 1)
         data <- scale.data(data, scale=scale)

# corr.func <- function(x, data, centroids) {
#   which.max(cor(cbind(data[,x], centroids), method="spearman", use="pairwise.complete.obs")[-1, 1])
# }
            for (samples in 1:ncol(data)){
               groups.idx <- which.max(cor(cbind(data[,samples], centroids), method=method, use="pairwise.complete.obs")[-1,1])
                  groups[samples] <- colnames(centroids)[groups.idx]
            }
   return(groups)
}

