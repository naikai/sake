#' Extract features from NMF run results 
#'
#' X = W x H
#' W is the feature matrix, H is the sample matrix
#' After NMF run, use this function to select import features and assign groups for the two matrix

#' @param res NMF run results
#' @keywords NMF feature
#' @export
#' @examples
#' nmf_extract_feature(data)

# Extract the feature Genes in each group
nmf_extract_feature <- function(res){
   library(NMF)

   feature.score <- featureScore(res)
   write.csv(feature.score, res.featureScore)

   # extracted features for each group
   extract.feature <- extractFeatures(res)

   data.extract.feature <- NULL
   for (i in 1:length(extract.feature)){
      idx <- extract.feature[[i]]
      Genes <- names(feature.score)[idx]
      featureScore <- feature.score[idx]
      Group <- i
      data.extract.feature <- rbind(data.extract.feature, cbind(Genes, featureScore, Group))
   }

   colnames(data.extract.feature) <- c("Genes", "featureScore", "Group")
   write.csv(data.extract.feature, res.extract.featureScore, row.names=F)
}
