#' A Cat Function
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @export
#' @examples
#' cat_function()
# DESeq normalization factor 
norm_factors <- function(mat) {
   nz <- apply(mat, 1, function(row) !any(round(row) == 0))
      mat_nz <- mat[nz,]
      p <- ncol(mat)
      geo_means <- exp(apply(mat_nz, 1, function(row) (1/p) * sum(log(row)) ))
      s <- sweep(mat_nz, 1, geo_means, `/`)

      apply(s, 2, median)
}

