#' A Cat Function
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @export
#' @examples
#' myHeatmap.3()

myd3Heatmap <- function (data, color=rev(brewer.pal(8, "RdYlBu")), type="cor", save.image=F, 
                                width=18, height=18, title="", file.prefix="heatmap", 
                                ColSideColors=NULL, ColSideColors.name=NULL, ColSideColorsSize=1, 
                                RowSideColors=NULL, RowSideColors.name=NULL, RowSideColorsSize=1, 
                                dendrogram="both", dendro.ord="auto", 
                                reorderfun= function(d, w) reorder(d, w), 
                                Colv=NULL, Rowv=NULL, na.rm=T,
                                cor.method="pearson", dist.m="euclidean", hclust.m="ward.D2", 
                                scale="row", scale.method="mean", scale.first=F,
                                notecex=1.4, notecol="black", 
                                show_grid=FALSE, anim_duration=0.1,
                                col.legend=F, row.legend=F, labRow=T, labCol=T, 
                                cexRow=1, cexCol=1, lhei=c(1.5,7), lwid=c(1,7)){
  require(vcd)
  require(RColorBrewer)
  require(d3heatmap)
  require(dendextend)
  if(save.image)
    pdf(paste0(file.prefix, ".pdf"), width=width, height=height)

  # for sample labeling 
  # label.cex <- assign_label_cex(dim(data)[2])

  if(labCol)
  labCol = colnames(data)
  if (ncol(data) > 150)
  labCol = rep("", ncol(data))

  ### scale data if type is heatmap 
  if(type=="cor"){
    data <- cor(data, method=cor.method)
    labRow = labCol
  }else if(type=="heatmap"){
    if(scale.first){
      data <- scale_data(data, scale=scale, method=scale.method)
    }

    if(dendro.ord=="auto"){
      if(dendrogram=="both" | dendrogram=="row"){
        Rowv <- data %>% dist(method=dist.m) %>% hclust(method=hclust.m) %>% as.dendrogram 
        RowMean <- rowMeans(data, na.rm = na.rm)
        Rowv <- reorderfun(Rowv, RowMean)
      }
      if(dendrogram=="both" | dendrogram=="col" | dendrogram=="column"){
        Colv <- data %>% t %>% dist(method=dist.m) %>% hclust(method=hclust.m) %>% as.dendrogram 
        ColMean <- colMeans(data, na.rm = na.rm)
        Colv <- reorderfun(Colv, ColMean)
      }
    }else if(dendro.ord=="manual"){
      if(is.dendrogram(Rowv)){
        Rowv <- Rowv
      }else{
        if(!is.null(Rowv) && !is.na(Rowv) && !identical(Rowv, FALSE)){
          warning("Invalid value for Rowv, ignoreing and assign NULL to it")
          Rowv <- NULL
        }
      }
      if(is.dendrogram(Colv)){
        Colv <- Colv
      }else{
        if(!is.null(Colv) && !is.na(Colv) && !identical(Colv, FALSE)){
          warning("Invalid value for Colv, ignoreing and assign NULL to it")
          Colv <- NULL
        }
      }
    }else{
      stop("unknown dendro, must be either 'auto' or 'manual'")
    }

    if(labRow)
      labRow = rownames(data)
    if (nrow(data) > 150)
      labRow = rep("", nrow(data))
  }else{
    stop("Unknown type! it must be either 'cor' or 'heatmap'")
  }

  # remove NA data 
  if(na.rm){
    data <- na.omit(data)
  }

  # decide whether to show ColSideColors (RowSideColors)
  if(!is.null(RowSideColors) & !is.null(ColSideColors)){
    # print("both")
    a <- d3heatmap(data, main=title, trace='none', key=TRUE, col=color, 
              dendrogram=dendrogram, na.rm=na.rm,
              Colv=Colv, labCol=labCol, ColSideColors=ColSideColors, ColSideColorsSize=ColSideColorsSize, 
              Rowv=Rowv, labRow=labRow, RowSideColors=RowSideColors, RowSideColorsSize=RowSideColorsSize,
              scale=scale, lwid=lwid, lhei=lhei, cexCol=cexCol, cexRow=cexRow, 
              margins=c(height, width), show_grid=show_grid,
              anim_duration=anim_duration )
  }else if(!is.null(RowSideColors)){
    # print("RowSideColors")
    a <- d3heatmap(data, main=title, trace='none', key=TRUE, col=color, 
              dendrogram=dendrogram, na.rm=na.rm,  
              Colv=Colv, labCol=labCol, 
              Rowv=Rowv, labRow=labRow, RowSideColors=RowSideColors, RowSideColorsSize=RowSideColorsSize,
              scale=scale, lwid=lwid, lhei=lhei, cexCol=cexCol, cexRow=cexRow, 
              margins=c(height, width), show_grid=show_grid,
              anim_duration=anim_duration )
  }else if(!is.null(ColSideColors)){
    # print("ColSideColors")
    a <- d3heatmap(data, main=title, trace='none', key=TRUE, col=color, 
              dendrogram=dendrogram, na.rm=na.rm,
              Colv=Colv, labCol=labCol, ColSideColors=ColSideColors, ColSideColorsSize=ColSideColorsSize, 
              Rowv=Rowv, labRow=labRow, 
              scale=scale, lwid=lwid, lhei=lhei, cexCol=cexCol, cexRow=cexRow, 
              margins=c(height, width), show_grid=show_grid,   
              anim_duration=anim_duration )
  }else{
    # print("No SideColors")
    a <- d3heatmap(data, main=title, trace='none', key=TRUE, col=color, 
              dendrogram=dendrogram, na.rm=na.rm,
              Colv=Colv, labCol=labCol, 
              Rowv=Rowv, labRow=labRow, 
              scale=scale, lwid=lwid, lhei=lhei, cexCol=cexCol, cexRow=cexRow, 
              margins=c(height, width), show_grid=show, 
              anim_duration=anim_duration )
  }

  if(save.image)
    dev.off()
    
  return(a)
}