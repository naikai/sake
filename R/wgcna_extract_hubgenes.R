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

# Extract the feature Genes in each group
wgcna_extract_hubgenes <- function(datExpr, MM.cutoff=0.7, math.cutoff.quantile=3, math="mad", top.math=1500, plot=F, save.data=F, prefix="MM_vs_MAD"){
   require(WGCNA)
   require(ggplot2)

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
