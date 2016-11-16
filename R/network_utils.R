#' Compute module gene list
#'
#' @param final_comp_genelist gene list
#' @param folderName0 file folder
#' @keywords co-expression, network, connectivity
#' @export
#' @examples
#' comp_module_genelist(final_comp_genelist, folderName0="folder")
comp_module_genelist <- function(final_comp_genelist, folderName0){
   pdf(paste0("Overlap", folderName0, ".pdf"), height=10, width=10)

   num <- length(final_comp_genelist)
   for(i in 2:num){
      comb <- combn(num, i)
      apply(comb, 2, function(x) {
         sub_final_comp_genelist <- final_comp_genelist[x]
         venn(sub_final_comp_genelist)
      })
   }
   dev.off()
}


#' Loop through each of the modules and connect the adjacent ones which is only allow.gap away from each other
#'
#' @param final_summary summary table
#' @param final_genelist summary gene list
#' @param allow.gap allow how man gaps between each sub-modules
#' @keywords co-expression, network, connectivity
#' @export
#' @examples
#' connect_gap_modules(final_summary, final_genelist, allow.gap=1)
connect_gap_modules <- function(final_summary, final_genelist, allow.gap=1){
   # create a column 'Delete' to specify which rows to be removed after we check through all the rows
   final_summary$Delete <- 0
   old <- 1

   for(current in 2:nrow(final_summary)){
      if( (final_summary[current, "Original_Start"] - final_summary[old, "Original_End"]) <= allow.gap ){
        # print("Ready to merged")
        # print(paste("Current Starts", final_summary[current, "Original_Start"], "Old ends", final_summary[old, "Original_End"]))

        # update these metrics for final_summary
        final_summary[current, "Size"] <- final_summary[current, "End"] -final_summary[old, "Start"] + 1
        new.avg.connectivity <- (final_summary[old, "Avg_Connectivity"] * final_summary[old, "Size"] + final_summary[current, "Avg_Connectivity"] * final_summary[current, "Size"]) / final_summary[current, "Size"]
        final_summary[current, "Avg_Connectivity"] <- new.avg.connectivity
        final_summary[current, "Start"] <- final_summary[old, "Start"]
        final_summary[current, "Start_Gene"] <- final_summary[old, "Start_Gene"]
        final_summary[current, "Original_Start"] <- final_summary[old, "Original_Start"]
        final_summary[old, "Delete"] <- 1

        ### update final_genelist
        final_genelist[[current]] <- unlist(c(final_genelist[[old]], final_genelist[[current]]))
      }
      old <- current
   }
   ### Old keep the ones that is not be Deleted, and then remove 'Delete' column
   idx <- final_summary$Delete == 1
   final_summary <- subset(final_summary, Delete != 1, select = -c(Delete))
   final_genelist[idx] <- NULL

   res <- list()
   res[['summary']] <- final_summary
   res[['genelist']] <- final_genelist
   return(res)
}


#' Identify the the edge of the modules (blocks) within the net connectivity matrix
#'
#' \preformatted{
#' 1. Run for loops to detect the block
#' 2. Detect edge genes in each block and their gene idx
#' 3. Go back to the original 100x100 matrix and use the edge genes (idx) to define block
#' }
#' @param data Input expression data
#' @param tao  threshold for hard transformation
#' @param beta parameter for soft power transformation
#' @param num_features Number of top high connective genes to be used
#' @keywords co-expression, network, connectivity
#' @export
#' @examples
#' def_net_modules(data, tao=0.5, beta=10, num_feature=0.01)
def_net_modules <- function(core_sub_jaccard_dist, sub_jaccard_dist, min.connectivity=0.2, min.size=5){
   old <- 1
   size <- 1
   blk_cnt <- 0
   # res <- list()
   res <- rep(NULL, 6)
   genenames <- colnames(core_sub_jaccard_dist)

   for (current in 2:dim(core_sub_jaccard_dist)[1]){
      # if (core_sub_jaccard_dist[current, current-1]==1){
      if (core_sub_jaccard_dist[current, current-1]>=min.connectivity){
         size <- size + 1
      }else{
         if(size>=min.size){
            blk_cnt <- blk_cnt + 1
            # Add another filter based on average connectivity of the block
            avg.connectivity <- sum(core_sub_jaccard_dist[old:current-1, old:current-1]) / ((current-1-old+1)*(current-1-old+1-1))
            # res[[blk_cnt]] <- c(old, current-1, size, genenames[old], genenames[current-1], avg.connectivity)
            res <- rbind(res, c(old, current-1, size, genenames[old], genenames[current-1], avg.connectivity))
         }
         size <- 1
         old <- current
      }
   }
   # last check #
   if (size >= min.size){
      blk_cnt <- blk_cnt + 1
      # res[[blk_cnt]] <- c(old, current, size, genenames[old], genenames[current])
      avg.connectivity <- sum(core_sub_jaccard_dist[old:current, old:current]) / ((current-old+1)*(current-old+1-1))
      res <- rbind(res, c(old, current, size, genenames[old], genenames[current], avg.connectivity))
   }

   colnames(res) <- c("Start", "End", "Size", "Start_Gene", "End_Gene", "Avg_Connectivity")
   return(res)
}


#' Identify the the edge of the modules (blocks) within the net connectivity matrix
#'
#' Run for loops to detect the block
#' Detect edge genes in each block and their gene idx
#' Go back to the original 100x100 matrix and use the edge genes (idx) to define block
#' @param core_sub_jaccard_dist preprocessed sub_jaccard_dist matrix (after applying filter, binarized the data)
#' @param sub_jaccard_dist zoomed-in small block size of jaccard_dist
#' @param min.connectivity threshold for minimum connectivity (set to 75\% quartile in each block)
#' @param min.size threshold for minimum block size
#' @param start_idx specifiy what is the index for the starting gene in the block
#' @param allow.gap allow how man gaps between each sub-modules
#' @keywords co-expression, network, connectivity
#' @export
#' @examples
#' def_net_modules2(core_sub_jaccard_dist, sub_jaccard_dist, min.connectivity=0.2, min.size=5)
def_net_module2 <- function(core_sub_jaccard_dist, sub_jaccard_dist, min.connectivity=0.2, min.size=5, start_idx=1, allow.gap=1){

   summarize_module <- function(){
      if(size>=min.size){
         blk_cnt <<- blk_cnt + 1
         # Add another filter based on average connectivity of the block
         avg.connectivity <- sum(core_sub_jaccard_dist[old:(current-1), old:(current-1)]) / ((current-1-old+1)*(current-1-old+1-1))
         # extract outer gene idx from the original 100x100 matrix #
         original_idx <- match(c(genenames[old], genenames[current-1]), original_names) + start_idx - 1
         res_summary <<- rbind(res_summary, c(old, current-1, size, genenames[old], genenames[current-1], avg.connectivity, original_idx[1], original_idx[2]))
         res_genelist[[blk_cnt]] <<- genenames[old:(current-1)]
      }
   }

   old <- 1
   size <- 1
   blk_cnt <- 0
   current_gap <- 0 # parameter to track how many gaps in between modules
   res <- list()
   res_genelist <- list()
   res_summary <- rep(NULL, 8)
   genenames <- colnames(core_sub_jaccard_dist)
   original_names <- colnames(sub_jaccard_dist)

   for (current in 2:dim(core_sub_jaccard_dist)[1]){
      # Allow for 1 (default) gap between each sub-modules
      if (core_sub_jaccard_dist[current, current-1]>=min.connectivity){
         size <- size + 1 + current_gap
         current_gap <- 0
      }else{
         summarize_module()
         size <- 1
         old <- current
         current_gap <- 0
      }
   }
   # last check #
   summarize_module()
   # summarize_module(core_sub_jaccard_dist, old, current, blk_cnt, genenames, original_names, res_summary, res_genelist, size, min.size=min.size, start_idx=start_idx)

   if(!is.null(res_summary)){
      colnames(res_summary) <- c("Start", "End", "Size", "Start_Gene", "End_Gene", "Avg_Connectivity", "Original_Start", "Original_End")
   }

   res[['summary']] <- res_summary
   res[['genelist']] <- res_genelist
   return(res)
}



#' Identify the modules in Data within co-expression network
#'
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords module, network, co-expression
#' @export
#' @examples
#' module_finder(mtcars)
module_finder <- function(data, p.value=0.05, tao=0.5, beta=3, num_features=100){
   data <- rmv_constant_0(data)
   # Use fastCor function to compute the correlation matrix
   xcor <- t(data) %>% fastCor %>% abs
   # net connectivity can be calculate using two transformation
   # 1. hard threshold, define a tao (signum function)
   # 2. soft threahold, define a beta (power function)
   if(thresh=="hard"){
     net_dist <- xcor > tao
   }else if(thresh=="soft"){
     net_dist <- xcor^beta
   }

   sum_connectivity <- colSums(net_dist)
   return(sum_connectivity)
}


#' Identify the informative genes in Data with top high connectivity in the co-expression network
#'
#' Code is adapted and rewrote based on Wange et al, BMC Bioinformatics, 2014
#' 'Improving the sensitivity of sample clustering by leveraging gene co-expression networks in variable selection'
#'
#' @param data Input expression data
#' @param tao  threshold for hard transformation
#' @param beta parameter for soft power transformation
#' @keywords co-expression, network, connectivity
#' @export
#' @examples
#' NetConnectivity(mtcars, tao=0.5)
NetConnectivity <- function(data, thresh='soft', tao=0.7, beta=1, method="pearson", diag.zero=T){
   # Use cor in WGCNA
   xcor0 <- t(data) %>% WGCNA::cor(., method=method) %>% abs
   if(diag.zero){
      system.time(diag(xcor0) <- 0) #for crossprod later
   }
   # net connectivity can be calculate using two transformation
   # 1. hard threshold, define a tao (signum function)
   # 2. soft threahold, define a beta (power function)
   if(thresh=="hard"){
      net_dist <- (xcor0 > tao) * 1
   }else if(thresh=="soft"){
      net_dist <- xcor0^beta
   }
   return(net_dist)
}


#' Loop through each of the modules and connect the adjacent ones which is only allow.gap away from each other
#'
#' @param final_summary summary table
#' @param final_genelist summary gene list
#' @param allow.gap allow how man gaps between each sub-modules
#' @keywords co-expression, network, connectivity
#' @export
#' @examples
#' plot_modules(core_sub_jaccard_dist, final_genelist)
plot_modules <- function(core_sub_jaccard_dist, final_genelist, allow.gap=1){
      if(sum(which.idx)>0){
         final_summary <- rbind(final_summary, res_summary[which.idx, ])
         final_genelist <- c(final_genelist, unlist(res_genelist[which.idx]))

         new_core_sub_jaccard_dist <- matrix(0, ncol=dim(core_sub_jaccard_dist)[1], nrow=dim(core_sub_jaccard_dist)[1], dimnames = dimnames(core_sub_jaccard_dist))
         for (j in which.idx){
            idx <- res_summary[j,1]:res_summary[j,2]
            new_core_sub_jaccard_dist[idx,idx] = core_sub_jaccard_dist[idx, idx]
            ### Can use bit vector calculation for faster speed? ###

            # save genelist for each modules
            filename <- paste0(paste(res_summary[j,c(4,5,7,8)], collapse="-"), ".txt")
            write.table(data.frame(Gene=unlist(res_genelist[j])), filename, quote=F, row.names=F)
         }
         bin_new_core_sub_jaccard_dist <- replace(new_core_sub_jaccard_dist, new_core_sub_jaccard_dist>0, 1)
         myHeatmap.3(bin_new_core_sub_jaccard_dist, type = "heatmap", Colv=NULL, Rowv=NULL, scale="none")
         myHeatmap.3(new_core_sub_jaccard_dist, type = "heatmap", Colv=NULL, Rowv=NULL, scale="none")
      }else{
          plot(1, main="No res fits the criteria", col="red")
      }
}


#' Preprocess the connecitivty matrix based on mean connectivity and binary switch
#'
#' \preformatted{
#' 1. First try to remove columns (rows) that have median expression 75% quantile
#' 2. Convert the matrix into binary based on quantile connectivity'
#' }
#'
#' @param data Input expression data
#' @param tao  threshold for hard transformation
#' @param beta parameter for soft power transformation
#' @param num_features Number of top high connective genes to be used
#' @keywords co-expression, network, connectivity
#' @export
#' @examples
#' preprocessNetCon(sub_jaccard_dist, rmv.filter=0.5, binary.filter=0.75, plot=F)
preprocessNetCon <- function(sub_jaccard_dist, rmv.filter=0.5, binary.filter=0.75, plot=F){
   # filter low mean connectivity genes
   keep.idx <- colSums(sub_jaccard_dist) >= median(colSums(sub_jaccard_dist))
   core_sub_jaccard_dist <- sub_jaccard_dist[keep.idx, keep.idx]
   if(plot){
      myHeatmap.3(core_sub_jaccard_dist, type = "heatmap", Colv=NULL, Rowv=NULL, scale="none", cexRow=0.5, cexCol=0.5)
   }

   # binary switch filter
   quantile(core_sub_jaccard_dist)
   top75.idx <- core_sub_jaccard_dist > quantile(core_sub_jaccard_dist)[4]
   bin_core_sub_jaccard_dist <- core_sub_jaccard_dist
   bin_core_sub_jaccard_dist[top75.idx] = 1
   if(plot){
      myHeatmap.3(bin_core_sub_jaccard_dist, type = "heatmap", Colv=NULL, Rowv=NULL, scale="none", cexRow=0.5, cexCol=0.5)
   }

   return(core_sub_jaccard_dist)
}


#' Preprocess the connecitivty matrix based on mean connectivity and binary switch
#' 1. First try to remove columns (rows) that have median expression 75% quantile
#' 2. Convert the matrix into binary based on quantile connectivity'
#'
#' @param data Input expression data
#' @param tao  threshold for hard transformation
#' @param beta parameter for soft power transformation
#' @param num_features Number of top high connective genes to be used
#' @keywords co-expression, network, connectivity
#' @export
#' @examples
#' run_CoExpression(expdata)
run_CoExpression <- function(expdata, thresh="soft", tao=0.7, beta=10, method="Pearson", folderName0="Module", avg_con_thresh=0.00015, save.data=F){
   # Start the clock!
   ptm <- proc.time()

   # save all the results into this folder
   if(!dir.exists(folderName0)){
      dir.create(folderName0)
   }
   setwd(folderName0)

   # Get NetConnectivit
   net_dist <- NetConnectivity(expdata, thresh=thresh, tao=tao, beta=beta, method=method)

   # if(method=="Pearson"){
   #    system.time(xcor0 <- t(expdata) %>% fastCor %>% abs)
   # }else if (method=="Spearman"){
   #    system.time(xcor0 <- t(expdata) %>% cor(., method="spearman") %>% abs)
   # }
   # system.time(diag(xcor0) <- 0) #for crossprod later

   ## Subset gene sets for performance testing,
   ## will use whole gene networks once the pipeline is done
   # n <- 10000
   # xcor <- xcor0[1:n, 1:n]
   # n <- 5000
   # xcor <- xcor0[1:n, 1:n]
   # xcor <- xcor0
   # clean up
   # rm(xcor0)

   # # net connectivity can be calculate using two transformation
   # # 1. hard threshold, define a tao (signum function)
   # # 2. soft threahold, define a beta (power function)
   # if(thresh=="hard"){
   #    net_dist <- (xcor > tao) * 1
   # }else if(thresh=="soft"){
   #    net_dist <- xcor^beta
   #    # stick with hard threshod for now.
   #    # avg_con_thresh <- tao^beta/(2-tao^beta)
   # }

   ## Intersected connectivity is the crossproduct of net_dist matrix
   invisible(gc())
   system.time(inter_conn <- crossprod(net_dist) )

   ### Final part of Jaccard index to include sum of connectivity ###
   invisible(gc())
   sum_connectivity <- colSums(net_dist)
   net_colnames <- colnames(net_dist)
   net_rowsize <- dim(net_dist)[1]
   net_colsize <- dim(net_dist)[2]
   rm(net_dist)

   # It should work by matrix division
   # First generate an index table to expanding the sum_connectivity into a matrix
   # Use [i,j] as index and inser [value] into this matrix
   # [i,j,value]
   system.time({
      print ("### Calculating Jaccard index")
      tmp <- cbind(expand.grid(sum_connectivity, sum_connectivity)) %>%
            data.table %>% '['(, values := rowSums(.SD))
      idx_table<- expand.grid(1:net_rowsize, 1:net_colsize) %>%
            cbind(., tmp) %>% as.matrix
      sum_dist <- matrix(0, nrow=net_rowsize, ncol=net_colsize)
      sum_dist[idx_table[, 1:2]] <- idx_table[, 5]
      jaccard_dist <- inter_conn / (sum_dist - inter_conn)
      jaccard_dist[is.na(jaccard_dist)] <- 0
      diag(jaccard_dist) <- 0
   })

   colnames(jaccard_dist) <- net_colnames
   rownames(jaccard_dist) <- net_colnames
   # clean up
   rm(sum_dist)
   rm(idx_table)
   rm(inter_conn)

   ### Optimize time to run hierarchical clustering on the final Jaccard similarity index matrix ###
   invisible(gc())
   pdf(paste0(folderName0,'genes_hierarchical_on_jaccard_index.pdf'), height=10, width=14)
   y <- HiClimR(jaccard_dist)
   dev.off()

   sample_hclust_order <- y$order
   sample_hclust_names <- y$labels  ## because it will remove input with zero variance
   sample_hclust_num <- length(sample_hclust_order)
   rm(y)
   ### Can store y in symmetric matrices. Ray ###

   ### Store final gene list for comparison. ###
   blocks <- c(100,150,200,300,500)
   final_comp_genelist <- vector("list", length(blocks))

   system.time({
      print("### Identifying modules for each block")
      for (block_idx in 1:length(blocks)){
         block_size <- blocks[block_idx]
         folderName <- paste0(folderName0, "-", "top", block_size)

         # save all the results into this folder
         if(!dir.exists(folderName)){
            dir.create(folderName)
         }
         setwd(folderName)
         dir.create("GeneList")
         dir.create("Modules_List")

         ### extract subsets of them and run module detection ###
         invisible(gc())
         # initialize a large number to be filled in, later only select the non-NULL ones
         N <- 20000
         fill_idx <- 1
         summary_colnames <- c("Start", "End", "Size", "Start_Gene", "End_Gene", "Avg_Connectivity", "Original_Start", "Original_End")
         final_summary <- data.frame(matrix(nrow=N, ncol=8)) %>% set_colnames(summary_colnames)
         final_genelist <- vector("list", N)

         pdf(paste0('top', block_size, 'genes_after_hierarchical_rawOrder.pdf'), height=14, width=14)
         for (i in 1:ceiling(sample_hclust_num/block_size)){
            start <- (i-1)*block_size+1
            end <- ifelse(i*block_size < sample_hclust_num, i*block_size, sample_hclust_num)
            if( (sample_hclust_num-end+1) < 50){
               # warning(paste0("Skipped for ", start, "-", end, " because it has less than 10 genes"))
               warning(paste0("Merged the next few (less than 50) genes into this block ", start, "-", end, " .Total ", sample_hclust_num, " genes"))
               next
            }

            plot(1, main=paste(start, ":", end))
            # sub.genes <- colnames(jaccard_dist)[sample_hclust_order][start:end]
            sub.genes <- sample_hclust_names[sample_hclust_order][start:end]
            sub_jaccard_dist <- jaccard_dist[sub.genes, sub.genes]
            myHeatmap.3(sub_jaccard_dist, type = "heatmap", Colv=NULL, Rowv=NULL, scale="none")
            setwd("GeneList/")
            write.table(data.frame(Gene=sub.genes), paste0("Gene",start,"-",end,".txt"), quote=F, row.names=F)

            ### Modules are defined here ###
            ### Clean up the 100x100 matrix and define module edges and size ###
            core_sub_jaccard_dist <- preprocessNetCon(sub_jaccard_dist, rmv.filter=0.5, binary.filter=0.75, plot=T)
            res <- def_net_module2(core_sub_jaccard_dist, sub_jaccard_dist, min.connectivity = quantile(core_sub_jaccard_dist)[4], min.size=5, start_idx=start)
            res_summary <- res[['summary']]
            res_genelist <- res[['genelist']]
            setwd("../Modules_List/")

            if(!is.null(res_summary)){
               which.idx <- which(as.numeric(res_summary[, "Avg_Connectivity"])>=avg_con_thresh)
               new_core_sub_jaccard_dist <- matrix(0, ncol=dim(core_sub_jaccard_dist)[1], nrow=dim(core_sub_jaccard_dist)[1], dimnames = dimnames(core_sub_jaccard_dist))

               if(sum(which.idx)>0){
                  for(j in which.idx){
                     final_summary[fill_idx, ] <- res_summary[j, ]
                     final_genelist[[fill_idx]] <- res_genelist[j]
                     fill_idx <- fill_idx + 1

                     # The following is for visual (physical) inspection of the module identification
                     idx <- res_summary[j,1]:res_summary[j,2]
                     new_core_sub_jaccard_dist[idx,idx] = core_sub_jaccard_dist[idx, idx]
                     # save gene names for each modules
                     filename <- paste0(paste(res_summary[j,c(4,5,7,8)], collapse="-"), ".txt")
                     write.table(data.frame(Gene=unlist(res_genelist[j])), filename, quote=F, row.names=F)
                  }
                   bin_new_core_sub_jaccard_dist <- replace(new_core_sub_jaccard_dist, new_core_sub_jaccard_dist>0, 1)
                   myHeatmap.3(bin_new_core_sub_jaccard_dist, type = "heatmap", Colv=NULL, Rowv=NULL, scale="none")
                   myHeatmap.3(new_core_sub_jaccard_dist, type = "heatmap", Colv=NULL, Rowv=NULL, scale="none")
               }else{
                  plot(1, main="No res fits the criteria", col="red")
               }
            }
            setwd("../")
         }
         dev.off()

         final_summary <- final_summary[1:(fill_idx-1), ]
         final_genelist <- final_genelist[1:(fill_idx-1)]
         # transform these column into numeric data type
         final_summary[, c(1:3,6:8)] <- sapply(final_summary[, c(1:3,6:8)], as.numeric)

         # save(final_summary, file="Rdata_final_summary")
         # save(final_genelist, file="Rdata_final_genelist")

         # connect those modules that have small gap in between
         final_res <- connect_gap_modules(final_summary, final_genelist, allow.gap=1)

         ### Now rerun through the final_summary table to map the start-end genes in module to original jaccard_matrix
         ### Then define module this way?

         ### final result sorted by module avg.connectivity
         final_summary[order(as.numeric(final_summary[,"Avg_Connectivity"]), decreasing = T),] %>%
               subset(as.numeric(.[,"Size"])>=10) %>%
               head(., n=20)
         sum(as.numeric(final_summary[,"Size"]))
         final_summary <- final_summary[order(as.numeric(final_summary[, "Avg_Connectivity"]), decreasing = T), ]
         write.table(final_summary, paste0("Modules_Summary_Block_top", block_size, "genes.txt"), row.names = F, quote=F, sep="\t")
         write.table(data.frame(Gene=unlist(final_genelist)), paste0("Modules_GeneList_Block_top", block_size, "genes.txt"), row.names = F, quote=F, sep="\t")

         ### For drawing Venn diagram to see how much overlap between diff block size ###
         final_comp_genelist[[block_idx]] <- unlist(final_genelist)

         # back to coexpression folder
         setwd("../")
      }
   })

   ### Plot overlap Venn diagram ###
   require(gplots)
   names(final_comp_genelist) <- paste0("Top", blocks, "\n", sapply(final_comp_genelist, length), " genes")
   comp_module_genelist(final_comp_genelist, folderName0)

   # stop the clock
   proc.time() - ptm

   if(save.data){
      save.image(paste0(".RData", folderName0))
   }
   setwd("../")
   return(0)
}

#' Identify the the edge of the modules (blocks) within the net connectivity matrix
#'
#' Run for loops to detect the block
#' Detect edge genes in each block and their gene idx
#' Go back to the original 100x100 matrix and use the edge genes (idx) to define block
#'
#' @param core_sub_jaccard_dist jaccard distance for sub core genes
#' @keywords network, co-expression
#' @export
#' @examples
#' summarize_module(core_sub_jaccard_dist, sub_jaccard_dist, min.connectivity=0.2, min.size=5)
summarize_module2 <- function(core_sub_jaccard_dist, old, current, blk_cnt, genenames, original_names, res_summary, res_genelist, size, min.size=5, start_idx){
   if(size>=min.size){
      blk_cnt <- blk_cnt + 1
      # Add another filter based on average connectivity of the block
      avg.connectivity <- sum(core_sub_jaccard_dist[old:(current-1), old:(current-1)]) / ((current-1-old+1)*(current-1-old+1-1))
      # extract outer gene idx from the original 100x100 matrix #
      original_idx <- match(c(genenames[old], genenames[current-1]), original_names) + start_idx - 1
      res_summary <<- rbind(res_summary, c(old, current-1, size, genenames[old], genenames[current-1], avg.connectivity, original_idx[1], original_idx[2]))
      res_genelist[[blk_cnt]] <<- genenames[old:(current-1)]
   }
}



#' Extract hub genes in each module by using WGCNA
#'
#' Weighted gene correlation network analysis
#'
#' @param datExpr gene count data
#' @param MM.cutoff module membership cutoff
#' @param math.cutoff.quantile which quantile to use as min cutoff
#' @param math default is 'mad': can be 'mean', 'median', 'iqr'
#' @keywords wgcna, co-expression, network
#' @export
#' @examples
#' wgcna_ext_hubgenes(datExpr)
wgcna_ext_hubgenes <- function(datExpr, MM.cutoff=0.7, math.cutoff.quantile=3, math="mad", top.math=1500, plot=F, save.data=F, prefix="MM_vs_MAD"){
   MM.cutoff <- MM.cutoff
   math.Expr <- log(apply(datExpr, 2, math))
   Exp.cutoff <- signif(quantile(math.Expr)[math.cutoff.quantile], digits = 3)

   mad.genes <- extract_data_by_mad(datExpr, topN=top.math, by="col", type="genes")

   # OVData.top1500mad.genes <- rownames(top.mad.data)

   total.num <- 0
   total.genes <- data.frame(Gene=NULL, Cor=NULL, Module=NULL)

   for(module in modNames){
     column = match(module, modNames);
     moduleGenes = moduleColors==module;
     current.module <- data.frame(Gene=colnames(datExpr)[moduleGenes],
                                   Exp=math.Expr[moduleGenes],
                                   MM=abs(geneModuleMembership[moduleGenes, column]),
                                   Module=module)

     mad.module <- subset(current.module, Gene %in% mad.genes)
     Cor <- signif(cor(current.module$Exp, current.module$MM), 2)
     Corp <- signif(corPvalueStudent(Cor, sum(is.finite(current.module$Exp) & is.finite(current.module$MM))), 2)

     if(plot){
         a <- ggplot(current.module, aes(x=MM, y=Exp, label=Gene)) +
               geom_point(alpha=0.5, colour=module) +
               geom_text(size=1.5, hjust = 0, nudge_x = 0.005) + theme_bw() +
               geom_vline(xintercept = MM.cutoff, colour="tomato") +
               geom_hline(yintercept = Exp.cutoff, colour="tomato") +
               xlab(paste("Module Membership in", module, "module")) +
               ylab("Log (expression)") +
               ggtitle(paste("Log expression vs. Module membership\n", module, "module:", sum(moduleGenes), "genes\n", "cor=", Cor, "p=", Corp)) +
               theme(plot.title = element_text(size=16, face="bold")) +
               theme(axis.text=element_text(size=16), axis.title=element_text(size=16, face="bold")) +
               theme(strip.text.x = element_text(size=16, face="bold")) +
               theme(axis.title.y=element_text(margin=margin(0,20,0,10))) +
               theme(axis.title.x=element_text(margin=margin(20,0,10,0)))

         # add red dots only if we have genes that are within TopMAD
         if(nrow(mad.module)>0){
           a <- a + geom_point(alpha=0.7, aes(x=MM, y=Exp), colour="red", data=mad.module) +
                  geom_text(size=1.5, hjust = 0, nudge_x = 0.005, data=mad.module, colour="red")
         }
         (gg <- ggplotly(a))
         # print(gg)
         htmltools::tagList(gg)
      }

      idx <- (abs(current.module$MM) >= MM.cutoff) & (current.module$Exp >= Exp.cutoff)
      total.num <- total.num + sum(idx)
      total.genes <- rbind(total.genes, current.module[idx,])
   }

   print(paste("total", total.num, "genes have cor above", MM.cutoff, "and Log(Expression) above", Exp.cutoff))
   if(save.data){
      # filename <- paste0(cancer, "_TopMAD", topN, "_Module_Genes_MMcutoff", MM.cutoff, "_Expcutoff", Exp.cutoff, ".txt")
      filename <- paste0(prefix, "_Module_Genes_MMcutoff", MM.cutoff, "_Expcutoff", Exp.cutoff, ".txt")
      write.table(total.genes, filename, sep="\t", quote=F, row.names=F)
   }

   return(total.genes)
}
