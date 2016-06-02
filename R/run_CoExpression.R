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
#' preprocessNetCon(sub_jaccard_dist, rmv.filter=0.5, binary.filter=0.75, plot=F)
# preprocessNetCon  
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
