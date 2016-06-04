#' Classify samples based on the feed in centroids
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords centroids classify 
#' @export
#' @examples
#' cat_function()
centroid_classify <- function(x, cv_data, groups, method="spearman") {
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


#' Create centroid based on the grouping info 
#'
#' This function allows you to express your love of cats.
#' @param data Input data
#' @param groups Groups for each samples 
#' @keywords centroid reshape
#' @export
#' @examples
#' create_centroid(mtcars)
create_centroid <- function(data, groups, scale="none"){
   cdata <- scale_data(data, scale=scale)
   cdata$gene <- rownames(cdata)
   ddata <- reshape2::melt(cdata, id="gene")
   ddata$group <- rep(groups, each=nrow(data))
   centroid <- reshape2::cast(ddata, group ~ gene, mean)
   return(t(centroid))
}


#' Run Cross-Validation 
#'
#' @param cv_data 
#' @keywords cross-validation
#' @export
#' @examples
#' run_CV(cv_data)
run_CV <- function(cv_data, groups, nrun=10, k.fold=10, method="spearman"){
   kappa_result <- NULL
   prop_result <- NULL
   for (i in 1:nrun){
      folds <- caret::createFolds(cv_data, k=k.fold)
      # lapply(folds, function(x) prop.table(table(groups[x])))
      cv_results <- sapply(folds, centroid_classify, cv_data=cv_data, groups=groups, method=method)
      kappa_result[i] <- mean(unlist(cv_results[1, ]))
      prop_result[i] <- mean(unlist(cv_results[2, ]))
   }
   return(data.frame(kappa=kappa_result, accuracy=prop_result))
}


#' Predict group assignment from centroids
#'
#' @param centroids centroids for each group
#' @param data data to be classified
#' @keywords centroid
#' @export
#' @examples
#' predict_from_centroid()
predict_from_centroid <- function(centroids, data, scale="none", method="spearman"){
   idx <- match(rownames(centroids), rownames(data))
   centroids <- centroids[!is.na(idx), ] 

   groups <- rep(NULL, ncol(data))
   if (dim(data)[2] > 1)
      data <- scale.data(data, scale=scale)

   for (samples in 1:ncol(data)){
      groups.idx <- which.max(cor(cbind(data[,samples], centroids), method=method, use="pairwise.complete.obs")[-1,1])
      groups[samples] <- colnames(centroids)[groups.idx]
   }
   return(groups)
}

