#' load predefined gene list files for heatmap
#'
#' Load predefined gene list files for heatmap (clustering) result
#' @param folder Where the files are
#' @param pattern file extensions for the ones you want to load
#' @keywords lodad 
#' @export
#' @examples
#' cor_mtest("data", ".txt")

cor_mtest <- function(mat, ...) {
    mat <- as.matrix(mat)
    n <- ncol(mat)
    p.mat<- matrix(NA, n, n)
    diag(p.mat) <- 0
    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            tmp <- cor.test(mat[, i], mat[, j], ...)
            p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
        }
    }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}