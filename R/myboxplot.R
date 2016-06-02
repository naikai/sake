#' A Cat Function
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @export
#' @examples
#' cat_function()
myboxplot <- function(data, prefix="Boxplot", method="vioplot", title="", ylim=c(0.6,0.85), save.plot=F, axis.text.size=14, axis.title.size=16, legend.position="right", vio.add.method="dot"){
# prefix <- paste0("NMF.ClaNC.Top",n,"MAD.Sel", select.features)
data.m <- melt(data)
groups <- colnames(data)

if(save.plot)
pdf(paste0(prefix, ".pdf"))

ray <- ggplot(data.m, aes(x=factor(variable), y=value)) + 
theme_bw() + ylim(ylim) + xlab("") + ylab("") + 
guides(fill=guide_legend(title=NULL)) + 
ggtitle(title) + 
theme(legend.position=legend.position) + 
theme(axis.text=element_text(size=axis.text.size), axis.title=element_text(size=axis.title.size, face="bold"))

if(method=="boxplot"){
   ray <- ray + geom_boxplot(aes(fill=variable)) + geom_jitter()
}else if (method=="vioplot"){
   ray <- ray + geom_violin(aes(fill=variable)) 
   if (vio.add.method=="none"){
      ray <- ray 
   }else if(vio.add.method=="dot"){
         ray <- ray + geom_boxplot(width=.1, fill="black") + stat_summary(fun.y=median, geom="point", col="white")
         }else if(vio.add.method=="boxplot"){
            ray <- ray + geom_boxplot(width=.1) 
         }else if(vio.add.method=="jitter"){
            ray <- ray + geom_violin(aes(fill=variable)) + geom_jitter() 
         }
   }
   print(ray)

   if(save.plot)
      dev.off()
}


