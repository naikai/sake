#' Classify samples based on the feed in centroids
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords centroids classify 
#' @export
#' @examples
#' cat_function()
centroid_classify <- function(x, cv_data, groups, method="spearman") {
   require(irr)
   credit_train <- data.frame(cv_data[, -x])
      credit_test  <- data.frame(cv_data[, x])
      groups_train <- groups[-x]
      groups_test  <- groups[x]
#   prop.table(table(groups_train))
#   prop.table(table(groups_test))
      credit_centroid <- create_centroid(credit_train, groups_train, scale="row")
      credit_pred  <- predict_from_centroid(credit_centroid, credit_test, scale="row", method=method)
      credit_actual <- groups_test
      prop <- sum(diag(table(credit_actual, credit_pred)))/length(x)
      kappa <- kappa2(data.frame(credit_actual, credit_pred))$value
      return(data.frame(prop=prop, kappa=kappa))
}


