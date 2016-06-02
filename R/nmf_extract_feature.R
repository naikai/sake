#' Extract features from NMF run results 
#'
#' X = W x H
#' W is the feature matrix, H is the sample matrix
#' After NMF run, use this function to select import features and assign groups for the two matrix

#' @param res NMF run results
#' @param rawdata Original expression matrix, if we want to rank the genes by 'MAD' value across samples 
#' @param manual.num How many genes to select from each of the group
#' @param method 'total', 'default', or 'rank' 
#' @param math Default is 'mad': median absolute deviation
#' @keywords NMF feature
#' @export
#' @examples
#' nmf_extract_feature(nmf_res)

# Extract the feature Genes in each group
nmf_extract_feature <- function(res, rawdata=NULL, manual.num=0, method="default", math="mad", FScutoff=0.9){
   require(NMF)

   feature.score <- featureScore(res)
   predict.feature <- predict(res, what="features", prob=T)

   data.feature <- data.frame(Gene=names(feature.score),
                              featureScore=feature.score,
                              Group=predict.feature$predict,
                              prob=predict.feature$prob, 
                              stringsAsFactors=FALSE)
   if(method=="total"){
      print("return all featureScores")
   }else{
      if(method=="default"){
         # extracted features for each group
         if (manual.num==0){
            extract.feature <- extractFeatures(res) 
         }else if (manual.num>0 && manual.num<=length(featureNames(res))){
            extract.feature <- extractFeatures(res, manual.num) 
         }else{
               stop("wrong number of (manual num) features ")
         }

         data.feature <- extract.feature %>% 
                     lapply(., function(x) data.feature[x, ]) %>% 
                     rbindlist %>% 
                     as.data.frame.matrix
      }else if(method=="rank"){
         if(is.null(rawdata)){
            stop("error: need to provide original expression data if method is 'rank' ")
         }
         data.feature <- cbind(data.feature, math=apply(log2(rawdata+1), 1, math)) %>%
                              filter(featureScore>=FScutoff) %>%
                              arrange(Group, dplyr::desc(math), dplyr::desc(prob)) 
         if(manual.num>0){
            data.feature <- group_by(data.feature, Group) %>%
                              top_n(manual.num, math)
                              # filter(min_rank(desc(math))<=manual.num) 
         }
      }
   }

   return(data.feature)
}
