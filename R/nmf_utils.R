#' Function to run NMF with different options
#'
#' This function allows you to perform matrix factorization using Non-negative matrix factorization (NMF) method
#' User need to provide filtered/ordered gene expression data and client info (expinfo) files
#' Workflow ###
#' 1 estimate how many clusters (k) ###
#' 1.1 estim.r <- nmf(data, 2:10, nrun=100)
#' 1.1 random shuffle the row (genes) for each column (sample) - estim.random.r
#' 1.2 generate plot(estim.r), consensusmap(estim.r), plot(estim.r, estim.random.r)
#' 2 After looking through the graph and decide k # This can be automated?
#' 2.1 res <- nmf(data, k, nrun=100?)
#' 2.2 compare the original.data and fitted.data (heatmap.2) # ignore this.
#' 2.3 summary (res), summary(res, target=data), summary(res, class=assigned.subtypes)
#' 2.4 plot basismap(res) and coefmap(res)
#' 2.5 extract metagene-specific features
#' @param data Input data sets
#' @param prefix Prefix for file output. Default is 'NMF'
#' @param cluster Estimated rank (clusters) in the data sets. Default is 3
#' @param run How many runs to perform? Default is 100
#' @param algorithm Which algorithms for NMF? (bruent, lee, nsNMF, KL, Frobenius, offset, ls-nmf, pe-nmf, siNMF). Default is brunet
#' @param mode Which modules to run? (Estim or Real). Default is real
#' @keywords NMF cluster
#' @export
#' @examples
#' myNMF(data, mode="estim", cluster=3, nrun=20)
myNMF <- function(data, prefix="NMF", cluster=3, top=1500, nrun=100, norm=F, algorithm="brunet", mode="real", seed=123211){
   if(mode=="estim"){
      r <- 2:cluster
      estim.r <- nmf(data, r, algorithm, .opt="vp24", nrun=nrun, seed=seed, maxIter=5000)
      res <- estim.r
   }else if(mode=="compare"){
      res.multi.method <- nmf(data, cluster, list("brunet", "lee", "ns"), nrun=nrun, .opt="vtp24", seed=seed)
      # compare(res.multi.method)
      # print(compare(res.multi.method))
      res <- res.multi.method
   }else if(mode=="real"){
      # option 't' will toggle error track function
      res <- nmf(data, cluster, algorithm, .opt="vtp24", nrun=nrun, seed=seed, maxIter=5000)
   }

   return(res)
}


#' Get NMF run summary stats
#'
#' X = W x H
#' W is the feature matrix, H is the sample matrix
#' After NMF run, use this function to select import features and assign groups for the two matrix
#' @param res NMF run results
#' @keywords NMF feature
#' @export
#' @examples
#' nmf_summary(data)
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


#' Generate individual plot for estimating K
#'
#' This function allows you to express your love of cats.
#' @param res NMF run result
#' @keywords silhouette
#' @export
#' @examples
#' nmf_estim_plot(estim.r)
nmf_estim_plot <- function(estim.r){
   nmf_rank <- estim.r$measures[, 1]
   for(i in 2:ncol(estim.r$measures)){
      ylabel <- colnames(estim.r$measures)[i]
      # plot(nmf_rank, estim.r$measures[, i], type="o", xlab="Rank", ylab=ylabel)
      a <- ggplot(estim.r$measures, aes_string(x="rank", y=ylabel)) + theme_bw() + geom_point(size=4) + geom_line() +
            theme(axis.text=element_text(size=15), axis.title=element_text(size=15)) +
            ggtitle(ylabel) + theme(plot.title = element_text(lineheight=.8, size=15, face="bold"))
      print(a)
   }
}


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
nmf_extract_feature <- function(res, rawdata=NULL, manual.num=0, method="default", math="mad", FScutoff=0.9){
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


#' Extract groups (sample clustering) from NMF run results
#'
#' X = W x H
#' W is the feature matrix, H is the sample matrix
#' After NMF run, use this function to select import features and assign groups for the two matrix
#' @param res NMF run results
#' @keywords NMF groups
#' @export
#' @examples
#' nmf_extract_group(nmf_res)
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


#' Summary plot for NMF result
#'
#' X = W x H
#' W is the feature matrix, H is the sample matrix
#' After NMF run, use this function to select import features and assign groups for the two matrix
#' @param res NMF run results
#' @param mode What kind of plot? (result, feature, consensus)
#' @keywords NMF feature
#' @export
#' @examples
#' nmf_plot(nmf_res)
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


#' Silhouette plot for NMF run
#'
#' This function allows you to express your love of cats.
#' @param res NMF run result
#' @keywords silhouette
#' @export
#' @examples
#' nmf_silhouette_plot(res)
nmf_silhouette_plot <- function(res, type="consensus", silorder=F){
   si <- silhouette(res, what=type)
   plot(si)
}


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
#' nmf_select_best_k(nmf_res)
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
