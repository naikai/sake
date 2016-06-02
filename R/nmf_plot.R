#' Extract features from NMF run results 
#'
#' X = W x H
#' W is the feature matrix, H is the sample matrix
#' After NMF run, use this function to select import features and assign groups for the two matrix

#' @param res NMF run results
#' @param mode What kind of plot? (result, feature, consensus)
#' @keywords NMF feature
#' @export
#' @examples
#' nmf_extract_feature(data)

# Generate summary plot 
nmf_plot <- function(res, type="consensus", subsetRow=TRUE, save.image=F, hclustfun="average", silorder=F, add_original_name=T){
   if(save.image)
      pdf(res.pdf, width=18, height=15)

   if(type=="result"){
      print(plot(res))
   }else{
      si <- silhouette(res, what=type)

      if(type=="features"){
         if(silorder){
            basismap(res, Rowv = si, subsetRow=subsetRow)
         }else{
            basismap(res, subsetRow = subsetRow)
         }
      }else if(type=="samples"){
         if(silorder){
            coefmap(res, Colv = si)
         }else{
            coefmap(res)
         }
      }else if(type=="consensus"){
         if(add_original_name){
            colnames(res@consensus) <- sampleNames(res)
            rownames(res@consensus) <- sampleNames(res)
         }
         consensusmap(res, hclustfun=hclustfun)
      }
   }

   if(save.image)
      dev.off()
}
