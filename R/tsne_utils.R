#' Use multi-core to run t-SNE
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


#' Extract and parse tsne results 
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @export
#' @examples
#' tsne_out <- run_tsne(mtcar, perplexity=3)
#' parse_tsne_res(tsne_out)
parse_tsne_res <- function(tsne_out){
   data <- as.data.frame(tsne_out$Y)
   num.samples <- dim(data)[2]
   colnames(data) <- c("x", "y")
   data$Names <- rownames(data)

   return(data)
}


#' Plot t-SNE results
#'
#' This function allows you to plot t-SNE results
#' @param tsne_out t-SNE results from run_tSNE
#' @param color color scheme for each dot
#' @param alpha alpha value for each dot
#' @param title figure title
#' @param brewer brewer color
#' @keywords tsne
#' @export
#' @examples
#' tsne_out <- run_tsne(mtcar, perplexity=3)
#' tsne_out <- parse_tsne_res(tsne_out)
#' plot_tsne(tsne_out)
plot_tsne <- function(tsne_out, color=NULL, alpha=1, title="tSNE", brewer="naikai", add.centroids=F, add.line=F, add.errorbar=F, add.label=F, label.size=3, conf=F, add.legend=F, save.plot=T, real.plot=T, point.size=3, centroid.size=6){
	data <- parse_tsne_res(tsne_out)
  if(is.null(color)){
    color <- rep("red", dim(tsne_out$Y)[1])
  }
	data$color <- color 
	filename <- title 
	min.cost <- signif(tsne_out$itercosts[length(tsne_out$itercosts)], digits=2)
	title <- paste(title, "\nmin.cost=", min.cost)
	# add manual color scheme 
	colors <- create.brewer.color(data$color, length(unique(color)), brewer)
	# add plot centroid function 
      if (add.centroids){
         if(add.line){
            gg <- merge(data, aggregate(cbind(mean.x=x,mean.y=y)~color, data, mean), by="color")
               a <- ggplot(gg, aes(x=x, y=y, colour=color)) + geom_point(size=3) + ggtitle(title) + theme_bw() + 
               geom_point(aes(x=mean.x, y=mean.y), size=centroid.size) +
               geom_segment(aes(x=mean.x, y=mean.y, xend=x, yend=y))
         }else if (add.errorbar){
            centroids <- aggregate(cbind(x,y)~color, data, mean)
               f <- function(z) sd(z)/sqrt(length(z))  # function to calculate std.err
               if (conf)
                  f <- function(z) qt(0.025,df=length(z)-1, lower.tail=F)* sd(z)/sqrt(length(z)) 
                     se <- aggregate(cbind(se.x=x, se.y=y)~color, data, f)
                     centroids <- merge(centroids,se, by="color")  # add std.err column to centroids
                     print(centroids)

                     a <- ggplot(data, aes(x=x, y=y, colour=color)) + geom_point(size=point.size) + ggtitle(title) + theme_bw() + 
                     geom_point(data=centroids, size=6) + 
                     geom_errorbar(data=centroids, aes(ymin=y-se.y, ymax=y+se.y), width=0.1)+
                     geom_errorbarh(data=centroids, aes(xmin=x-se.x,xmax=x+se.x), height=0.1)
         }else{
            centroids <- aggregate(cbind(x,y)~color, data, mean)
               a <- ggplot(data, aes(x=x, y=y, colour=color)) + geom_point(size=point.size) + ggtitle(title) + theme_bw() + 
               geom_point(data=centroids, size=centroid.size)
         }
      }else{
         a <- ggplot(data, aes(x=x, y=y, colour=color)) + geom_point(aes(text=Names), colour=colors, size=point.size, alpha=alpha) + ggtitle(title) + theme_bw()
      }

      if (add.label){
         a <- a + geom_text(data=data, aes(label=Names), hjust=0.5, vjust=2, size=label.size, colour=colors)
      }
      if (!add.legend){
         a <- a + theme(legend.position="none")
      }
      # rename x and y label 
      a <- a + xlab("Component1") + ylab("Component2")
      if (save.plot)
      pdf(paste0(filename, ".pdf"), height=13, width=13)
      if(real.plot)
         print(a)
      if(save.plot)
      dev.off()
   return(a)
}
