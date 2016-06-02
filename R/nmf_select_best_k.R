#' Select the best k from multiple consensus NMF run 
#'
#' X = W x H
#' W is the feature matrix, H is the sample matrix
#' After NMF run, use Cophenatic index to select the best k 
#' highest value after k=2 

#' @param res NMF run results
#' @keywords NMF feature
#' @export
#' @examples
#' nmf_extract_group(nmf_res)

# Extract the feature Genes in each group
nmf_select_best_k <- function(res, type="consensus"){
    data <- NULL

    if(type=="consensus"){
      predict.consensus <- predict(res, what="consensus")
      data <- data.frame(Sample_ID=sampleNames(res),
                         nmf_subtypes = predict.consensus)
      # If we want to display as we see in consensusmap, we just need to reoder everything.
      # Now re-order data to match consensusmap sample order 
      if(matchConseOrder){
        sample.order <- attributes(predict.consensus)$iOrd
        data <- data[sample.order, ]
      }
    }else if(type=="samples"){
      predict.samples <- predict(res, what="samples", prob=T)
      data <- data.frame(Sample_ID=names(predict.samples$predict),
                         nmf_subtypes = predict.samples$predict,
                         prob = predict.samples$prob)
    }else{
      stop(paste("Wrong type:", type))
    }
    return(data)
}
