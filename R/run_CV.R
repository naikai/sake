#' A Cat Function
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @export
#' @examples
#' cat_function()
run_CV <- function(cv_data, groups, nrun=10, k.fold=10, method="spearman"){
   require(caret)
   kappa_result <- NULL
      prop_result <- NULL
      for (i in 1:nrun){
         folds <- createFolds(cv_data, k=k.fold)
# lapply(folds, function(x) prop.table(table(groups[x])))
            cv_results <- sapply(folds, centroid_classify, cv_data=cv_data, groups=groups, method=method)
            kappa_result[i] <- mean(unlist(cv_results[1, ]))
            prop_result[i] <- mean(unlist(cv_results[2, ]))
      }
   return(data.frame(kappa=kappa_result, accuracy=prop_result))
}
