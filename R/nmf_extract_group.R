#' Extract groups (sample clustering) from NMF run results 
#'
#' X = W x H
#' W is the feature matrix, H is the sample matrix
#' After NMF run, use this function to select import features and assign groups for the two matrix

#' @param res NMF run results
#' @keywords NMF feature
#' @export
#' @examples
#' nmf_extract_group(nmf_res)

# Extract the feature Genes in each group
nmf_extract_group <- function(res, type="consensus", matchConseOrder=F){
    data <- NULL

    if(type=="consensus"){
      predict.consensus <- predict(res, what="consensus")
      silhouette.consensus <- silhouette(res, what="consensus")
      # It turns out the factor levels is the NMF_assigned_groups from consensus matrix
      # that matches the original sampleNames(res) order
      # The attributes(a.predict.consensus)$iOrd is the idx order for it to match the
      # order of the samples in consensusmap(res). It is just for displaying
      # Therefore, the merged data frame sampleNames(res) + a.predict.consensus is the final
      # consensus results.
      data <- data.frame(Sample_ID=sampleNames(res),
                         nmf_subtypes = predict.consensus,
                         sil_width = signif(silhouette.consensus[, "sil_width"], 3))
      # If we want to display as we see in consensusmap, we just need to reoder everything.
      # Now re-order data to match consensusmap sample order 
      if(matchConseOrder){
        sample.order <- attributes(predict.consensus)$iOrd
        data <- data[sample.order, ]
      }
    }else if(type=="samples"){
      predict.samples <- predict(res, what="samples", prob=T)
      silhouette.samples <- silhouette(res, what="samples")
      data <- data.frame(Sample_ID=names(predict.samples$predict),
                         nmf_subtypes = predict.samples$predict,
                         sil_width = signif(silhouette.samples[, "sil_width"], 3),
                         prob = signif(predict.samples$prob, 3))
    }else{
      stop(paste("Wrong type:", type, "Possible options are: 'consensus', 'samples' "))
    }
    return(data)
}
