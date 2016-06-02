#' function to run NMF with different options 
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
#' myNMF(data, expinfo, method="estim", k=3, nrun=100, )

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
