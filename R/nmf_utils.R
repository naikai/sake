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
myNMF <- function(data, prefix="NMF", cluster=3, nrun=100, norm=F, ncores=8, algorithm="brunet", mode="real", seed=123211){

   if(mode=="estim"){
      r <- 2:cluster
      estim.r <- nmf(data, r, algorithm, .opt=paste0("vp", ncores), nrun=nrun, seed=seed, maxIter=5000)
      res <- estim.r
   }else if(mode=="compare"){
      res.multi.method <- nmf(data, cluster, list("brunet", "lee", "ns"), nrun=nrun, .opt=paste0("vtp", ncores), seed=seed)
      # compare(res.multi.method)
      res <- res.multi.method
   }else if(mode=="real"){
      # option 't' will toggle error track function
      res <- nmf(data, cluster, algorithm, .opt=paste0("vtp", ncores), nrun=nrun, seed=seed, maxIter=5000)
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
         # data.feature <- cbind(data.feature, math=apply(log2(rawdata+1), 1, math)) %>%
         data.feature <- cbind(data.feature, math=apply(rawdata, 1, math)) %>%
                              filter(featureScore>=FScutoff) %>%
                              arrange(Group, dplyr::desc(math), dplyr::desc(prob))
         if(manual.num>0){
            data.feature <- group_by(data.feature, Group) %>%
                              top_n(manual.num)
                              # top_n(manual.num, math)
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


#' calculate weighted silhouette_width as threshold for running iter_NMF
#'
#' @param nmf_groups grouping infor after NMF run
weighted.silhouette <- function(nmf_groups){
  prop <- table(nmf_groups$nmf_subtypes) %>% prop.table
  group_by(inmf_groups, nmf_subtypes) %>% summarise(m = median(sil_width))
}


#' Run iterative NMF for K=2 until it cannot no longer separate samples
#'
#' @param rawdata gene count table
#' @param nrun number of runs for each NMF
#' @param min.sample number of required minimum samples in each subclusters
#' @param mad number of top MAD genes
#' @param silhouette min required threshold for silhouette index after NMF run
#' @keywords NMF iterative
#' @export
#' @examples
#' iter_nmf(rawdata, nrun=10, mad=5000)
iter_nmf <- function(rawdata, nrun=10, round=1, min.sample = 10, mad = 5000,
                     silhouette = 0.95, w.silhouette = 0.7, w.min.sample = 10, ncores=8) {
  repeat {
    cat(paste('\n\n### Running for iterNMF in round', round, "\n"))
    print("Dim of rawdata")
    print(dim(rawdata))

    data <- rmv_constant_0(rawdata, pct = 0.95) %>%
      extract_data_by_math(., topN = mad, math = "mad")

    nmf_res <- myNMF(data, cluster = 2, nrun=nrun, ncores = ncores)
    nmf_groups <- nmf_extract_group(nmf_res, type="samples")

    cur.silhouette <- median(nmf_groups$sil_width)
    print(paste('median silhouette is', signif(cur.silhouette, 3)))

    run_iter <- FALSE
    if (cur.silhouette < silhouette) {
      print("check here")
      if((weighted.silhouette(nmf_groups) > w.silhouette) & (ncol(data) > w.min.sample)){
        print(paste("weighted.silhouette is", weighted.silhouette, ". Proceed"))
        run_iter <- TRUE
      }else{
        print(paste("silhouette.consensus", cur.silhouette, "is less than:", silhouette, "don't proceed"))
        break
      }
    }else{
      # if both of the sub_clusters have greater than 6 samples
      if( sum(table(nmf_groups$nmf_subtypes) > 6) == 2){
        print("Fit min samples criteria in subclusters")
        # prop <- table(nmf_groups$nmf_subtypes) %>% prop.table
        run_iter <- TRUE
      }else{
        print("Sample distribution in subclusters is so different, probably due to outliers. Stop clustering here")
      }
    }

    if(run_iter){
      res <- list()
      for(grp in c(1,2)){
        idx <- which(nmf_groups$nmf_subtypes == grp)

        # less than 10(default) samples, don't proceed
        if(length(idx) <= min.sample){
          print(paste("subdata has", length(idx), "samples, which is less or equal than", min.sample, "samples, don't proceed"))
          # it turns out we need to have another NMF run in order to get NMF_res for the summary later
          # for cluster == 2, samples have to be greater than 3.
          if(length(idx) == 1){
            subdata <- rawdata[, c(idx, idx, idx)] %>% rmv_constant_0(.)
          }else if(length(idx) == 2){
            subdata <- rawdata[, c(idx, idx)] %>% rmv_constant_0(.)
          }else{
            subdata <- rawdata[, idx] %>% rmv_constant_0(.)
          }
          # sometimes we will get samples with all zeros in the data
          res[[grp]] <- myNMF(subdata[1:100, ] + 0.0001, cluster=2, nrun=2, ncores = ncores)
        }else{
          # repeat iter_NMF
          subdata <- rawdata[, idx] %>% rmv_constant_0(.)
          res[[grp]] <- iter_nmf(subdata, round = paste0(round, "-", grp),
                                 nrun=nrun, min.sample=min.sample,
                                 mad=mad, silhouette = silhouette,
                                 ncores = ncores)
        }
      }
      res = list(res = res)
      return(res)
    }else{
      break
    }
  }
  return(nmf_res)
}

extract_centroid_from_nmf_res <- function(res, data, method = "mean"){
  nmf_groups <- nmf_extract_group(res)
  data[, colnames(data) %in% nmf_groups$Sample_ID] %>%
    getcentroid(., method = method)
}

#' Get centroid from data frame
#'
#' @param df gene count table
#' @param method mean or median for calculating centroid
#' @keywords NMF centroid iterative
#' @export
getcentroid <- function(df, method = "mean") {
  # only 1 group in this nmf_group
  if(is.null(ncol(df))){
    return(df)
  }else{
    if(method == "mean"){
      return(rowMeans(df))
    }else if(method == "median"){
      return(rowMedians(df))
    }else{
      stop("Please check you have specified the correct method")
    }
  }
}

#' Calculate silhouette index between two data frame
#'
#' @param clust1 data in cluster1
#' @param clust2 data in cluster2
#' @param method method for calculating distance metrics
#' @keywords NMF iterative silhouette
#' @export
clust_silhouette <- function(clust1, clust2, method="euclidean"){
  sil.data <- list()
  sil.data$data <- cbind(clust1, clust2)
  sil.data$clustering <- c(rep(1, ncol(clust1)), rep(2, ncol(clust2)))
  sil.dist <- dist(t(sil.data$data), method = method)
  silhouette(sil.data, sil.dist)
}

#' Extract, cleanup, and merge iter_NMF results
#'
#' We will only get grouping info for now
#' @param iter_nmf vairous NMF res in a list format
#' @param rawdata gene count table for calculating centroids based on the clustering results
#' @param cor.method method for correlation calculation, default: spearman correlation
#' @param cor.thresh correlation threshold for merging subclusters
#' @param ctroid.topN number of top Genes to extract for centroid comparison
#' @param ctroid.math method for getting centroid, default: mean
#' @keywords NMF iterative
#' @export
#' @examples
#' clean_iterNMF(inmf_res, rawdata)
clean_iterNMF <- function(inmf_res, rawdata,
                          cor.method = "spearman", cor.thresh = 0.9,
                          ctriod.topN = 5000, ctroid.math = "mean",
                          return.all = FALSE) {

  bb <- unlist(inmf_res)
  names(bb) <- names(bb) %>% make.unique()

  # the number of feature selected for centroid will affect how we merge subclusters
  # using spearman.cor, purpose to rmv_constant_0 and select top rowMean 5K genes
  clust_data <- rawdata %>%
    rmv_constant_0(., pct=0.95) %>%
    extract_data_by_math(., topN = ctriod.topN, math = ctroid.math)

  clust_centroid <- bb %>%
    sapply(., function(x) extract_centroid_from_nmf_res(res = x, data = clust_data))

  cor_cluster <- clust_centroid %>%
    cor(., method=cor.method) %>%
    reshape2::melt() %>%
    filter(value != 1 & duplicated(value)) %>%
    dplyr::arrange(Var1, Var2)

  if(cor.thresh > 0){
    cor_cluster <- cor_cluster %>% filter(value > cor.thresh)
  }

  ### Check number of samples in each sub-clusters
  cc <- names(bb) %>% as_tibble() %>%
    dplyr::select(names = value ) %>%
    dplyr::mutate(res = bb[names]) %>%
    dplyr::mutate(nmf_groups = purrr::map(res, ~ nmf_extract_group(.)))

  if(nrow(cor_cluster) > 0){
    ### Now loop through cor_cluster and merge sub-clusters
    res <- list()
    for(i in 1:nrow(cor_cluster)){
      var1 <- cor_cluster[i, "Var1"] %>% as.character()
      var2 <- cor_cluster[i, "Var2"] %>% as.character()
      res[[var1]] <- c(res[[var1]], var2)
      res[[var2]] <- c(res[[var2]], var1)
    }
    # now loop through 5 times and assume this will stablize?
    res2 <- res
    for(i in 1:5){
      for(clust1 in names(res2)){
        sub_clust <- res[[clust1]]
        for(clust2 in sub_clust){
          res2[[clust1]] <- c(res2[[clust1]], res2[[clust2]])
        }
        res2[[clust1]] <- unique(res2[[clust1]])
      }
    }
    #
    dd <- cc %>% mutate(remove = 0)
    for(clust in names(res2)){
      idx1 <- which(dd$names == clust)
      if(dd[idx1, "remove"] == 0){

        for(subclust in res2[[clust]]){
          if(clust != subclust){
            idx2 <- which(dd$names == subclust)
            dd$res[idx1] <- list(c(dd$res[idx1], dd$res[idx2]))
            dd$nmf_groups[idx1] <- rbind(dd$nmf_groups[idx1][[1]], dd$nmf_groups[idx2][[1]]) %>% list
            dd[idx2, "remove"] <- 1
          }
        }

      }
    }
  }else{
    # no any subclusters are close to each other
    dd <- cc %>% mutate(remove=0)
  }

  ### Remove merged subclusters
  ee <- filter(dd, remove == 0) %>%
    mutate(grp = 1:nrow(.)) %>%
    mutate(ori_sample_id = purrr::map(nmf_groups, ~ .$Sample_ID),
           new_sample_id = purrr::map2(ori_sample_id, grp, ~ paste(paste0("iNMF", .y), .x, sep="_")))

  if(return.all){
    final_res <- list(res=ee,
                      cor_cluster = cor_cluster,
                      groups = cc,
                      clust_centroid = clust_centroid)
  }else{
    final_res <- list(res=ee)
  }
  return(final_res)
}



# check correlation between subclusters
#' @param iter_nmf_res vairous NMF res in a list format
#' @param cor_cluster correlation between subcluster
#' @param clust_centroid centroids in each subcluster
#' @keywords NMF iterative
#' @export
#' @examples
plot_cor_subiNMF <- function(inmf_res, cor_cluster, clust_centroid){

  for(i in 1:nrow(cor_cluster)){
    z1 <- inmf_res %>% filter(names == cor_cluster$Var1[i]) %>% .$nmf_groups %>% .[[1]]
    z2 <- inmf_res %>% filter(names == cor_cluster$Var2[i]) %>% .$nmf_groups %>% .[[1]]

    z1_name <-
      gsub(".*_", "", z1$Sample_ID) %>%
      table %>%
      data.frame %>%
      mutate(Freq = paste(., Freq, sep="-")) %>%
      dplyr::select(Freq) %>%
      unlist() %>%
      paste(., collapse=";")

    z2_name <-
      gsub(".*_", "", z2$Sample_ID) %>%
      table %>%
      data.frame %>%
      mutate(Freq = paste(., Freq, sep="-")) %>%
      dplyr::select(Freq) %>%
      unlist() %>%
      paste(., collapse=";")

    c.cor <- cor(clust_centroid[, cor_cluster$Var1[i]], clust_centroid[, cor_cluster$Var2[i]], method = "spearman")
    smoothScatter(clust_centroid[, cor_cluster$Var1[i]], clust_centroid[, cor_cluster$Var2[i]],
                  # colramp=viridis,
                  colramp = colorRampPalette(c("white", blues9)),
                  main = signif(c.cor,3), xlab = z1_name, ylab = z2_name)
  }
}

