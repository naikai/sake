#' Extract features from NMF run results 
#'
#' X = W x H
#' W is the feature matrix, H is the sample matrix
#' After NMF run, use this function to select import features and assign groups for the two matrix

#' @param res NMF run results
#' @keywords NMF feature
#' @export
#' @examples
#' nmf_summary(data)

# Extract the feature Genes in each group
nmf_summary <- function(res, class=NULL, save.data=F, filename="nmf.summary.txt"){
   summary.class <- NULL
   if(!is.null(class)){
      summary.class <- as.data.frame(NMF::summary(res, class=ann$Subtype))
   }else{
      summary.class <- as.data.frame(NMF::summary(res))
   }

   if(save.data){
      write.csv(summary.class, filename)
   }

   return(summary.class)
}
