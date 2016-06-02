#' t-SNE run
#'
#' This function allows you to run t-SNE using multi-cores 
#' @param data expression data
#' @param iter how many iterations 
#' @param perplexity perplexity
#' @param dims how many dimensions
#' @param cores how many cores to run t-SNE
#' @keywords t-SNE 
#' @export
#' @examples
#' run_tsne(mtcars, iter=200, perplexity = 2, cores=4)
# run tSNE 
run_tsne <- function(data, iter=10, perplexity=30, dims=2, cores=1){
   require(Rtsne)
   tsne_out <- NULL
   ### single-core execution 
   if(cores==1){
      cost <- 1000
      for(i in 1:iter){
         temp <- Rtsne(as.matrix(t(data)), perplexity=perplexity, dims=dims)
         min.cost <- temp$itercosts[length(temp$itercosts)]
         if( min.cost < cost){
            cost <- min.cost
            tsne_out <- temp 
         }
      }
   }else if(round(cores) > 1){
      require(snowfall)

      sfInit(parallel = TRUE, cpus=cores, type="SOCK")
      sfExport('data', 'perplexity', 'dims')
      sfLibrary(Rtsne)

      processInput <- function(j){
        res <- Rtsne(as.matrix(t(data)), perplexity=perplexity, dims=dims)
        return(res)
      }
      tmp_out <- vector("list", length = iter)
      system.time(tmp_out <- sfLapply(seq(1:iter), processInput))
      tsne_out <- lapply(tmp_out, function(x) min(x$itercost)) %>% which.min %>% tmp_out[[.]]
      sfStop()
   }else{
      stop("Error: wrong number of cores, it must be a positive integer.")
   }

   rownames(tsne_out$Y) <- colnames(data)
   return(tsne_out)
}

