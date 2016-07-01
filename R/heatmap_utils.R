#' A Cat Function
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @export
myHeatmap.3 <- function (data, color=rev(brewer.pal(8, "RdYlBu")), type="cor", save.image=F,
                                width=18, height=18, title="", file.prefix="heatmap",
                                ColSideColors=NULL, ColSideColors.name=NULL, ColSideColorsSize=1,
                                RowSideColors=NULL, RowSideColors.name=NULL, RowSideColorsSize=1,
                                reorderfun= function(d, w) reorder(d, w),
                                Colv=NULL, Rowv=NULL, na.rm=T,
                                cor.method="pearson", dist.m="euclidean", hclust.m="ward.D2",
                                scale="row", scale.method="mean", scale.first=F,
                                dendrogram="both", dendro.ord="auto",
                                notecex=1.4, notecol="black", col.legend=F, row.legend=F,
                                labRow=T, labCol=T, cexRow=1.1, cexCol=1.1, lhei=c(1.3,7), lwid=c(1,7))
{
  if(save.image)
    pdf(paste0(file.prefix, ".pdf"), width=width, height=height)

  # for sample labeling
  label.cex <- assign_label_cex(dim(data)[2])

  # ColSideColors
  if(is.null(ColSideColors)){
    ColSideColors = rep("white", ncol(data))
  }
  ColSideColors <- as.matrix(ColSideColors)

  if(type=="cor"){
    data <- cor(data, method=cor.method)
    # cellnote <- format(data, digits=1, nsmall=2)
    cellnote <- format(round(data, digits=2), nsmall=2)
    if (dim(data)[2]>150)
      notecex=0.01

    Colv <- data %>% t %>% dist(method=dist.m) %>% hclust(method=hclust.m) %>% as.dendrogram
    Rowv <- Colv
    # Rowv <- data %>% dist(method=dist.m) %>% hclust(method=hclust.m) %>% as.dendrogram
    heatmap.3(data, main=title, trace="none", scale="none", key=TRUE, col=color,
            cellnote=cellnote, notecex=notecex, notecol=notecol,
            ColSideColors=ColSideColors, ColSideColorsSize=ColSideColorsSize, Colv=Colv, Rowv=Rowv,
            lwid=lwid, lhei=lhei,
            labCol=rownames(data), labRow=rownames(data), cexRow=label.cex, cexCol=label.cex,
            margins=c(height, width) )

  }else if(type=="heatmap"){
    ### This needs to be careful, by default myHeatmap.3 use manual feed in Colv and Rowv, so need to add dendro="auto" when call this functions
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
      if(dendextend::is.dendrogram(Rowv)){
        Rowv <- Rowv
      }else{
        if(!is.null(Rowv) && !is.na(Rowv) && !identical(Rowv, FALSE)){
          warning("Invalid value for Rowv, ignoreing and assign NULL to it")
          Rowv <- NULL
        }
      }
      if(dendextend::is.dendrogram(Colv)){
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

    # cexRow = assign_label_cex(nrow(data))
    if(labRow)
      labRow = rownames(data)
    if (nrow(data) > 150)
      labRow=F
    if(labCol)
      labCol = colnames(data)
    if (ncol(data) > 150)
      labCol=F

    if(!is.null(RowSideColors)){
      heatmap.3(data, main=title, trace='none', key=TRUE, col=color,
                dendro=dendrogram,
                Colv=Colv, labCol=labCol, ColSideColors=ColSideColors, ColSideColorsSize=ColSideColorsSize,
                Rowv=Rowv, labRow=labRow, RowSideColors=RowSideColors, RowSideColorsSize=RowSideColorsSize,
                scale=scale, lwid=lwid, lhei=lhei, cexCol=cexCol, cexRow=cexRow,
                margins=c(height, width) )
    }else{
      heatmap.3(data, main=title, trace='none', key=TRUE, col=color,
                dendro=dendrogram,
                Colv=Colv, labCol=labCol, ColSideColors=ColSideColors, ColSideColorsSize=ColSideColorsSize,
                Rowv=Rowv, labRow=labRow,
                scale=scale, lwid=lwid, lhei=lhei, cexCol=cexCol, cexRow=cexRow,
                margins=c(height, width) )
    }

    if (col.legend){
      if(!is.null(ColSideColors.name)){
        legend.name <- NULL
        legend.color <- NULL
        for(i in 1:ncol(ColSideColors.name)){
          legend.name <- c(legend.name, "", unique(unlist(ColSideColors.name[, i])))
          legend.color <- c(legend.color, "white", unique(unlist(ColSideColors[, i])))
        }
        legend("topright", legend=legend.name, pch=19, ncol=1, cex=0.8, col=legend.color)
      }
    }
    # Rotate the legend
    if (row.legend){
      if(!is.null(RowSideColors.name)){
        g <- vcd::grid_legend("bottomleft", labels=unique(RowSideColors.name),
              draw=F, col=unique(as.character(RowSideColors)), pch=19, gp=grid::gpar(cex=0.8),
              hgap=grid::unit(0.6, "lines"), vgap=grid::unit(0.3, "lines") )
        grid::grid.draw(grid::grobTree(g, vp=grid::viewport(x=.1, y=.08, angle=90)))
      }
    }
  }

  if(save.image)
    dev.off()
}



#' Custom version of d3heatmap
#'
#' Add more option to control d3heatmap for shiny app
#' @param data data.frame for plot
#' @keywords heatmap, d3heatmap
#' @export
#' @examples
#' myd3Heatmap(mtcars)
myd3Heatmap <- function (data, color=rev(brewer.pal(8, "RdYlBu")), type="cor",
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
      if(dendextend::is.dendrogram(Rowv)){
        Rowv <- Rowv
      }else{
        if(!is.null(Rowv) && !is.na(Rowv) && !identical(Rowv, FALSE)){
          warning("Invalid value for Rowv, ignoreing and assign NULL to it")
          Rowv <- NULL
        }
      }
      if(dendextend::is.dendrogram(Colv)){
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
    a <- d3heatmap(data, main=title, trace='none', key=TRUE, col=color,
              dendrogram=dendrogram, na.rm=na.rm,
              Colv=Colv, labCol=labCol, ColSideColors=ColSideColors, ColSideColorsSize=ColSideColorsSize,
              Rowv=Rowv, labRow=labRow, RowSideColors=RowSideColors, RowSideColorsSize=RowSideColorsSize,
              scale=scale, lwid=lwid, lhei=lhei, cexCol=cexCol, cexRow=cexRow,
              margins=c(height, width), show_grid=show_grid,
              anim_duration=anim_duration )
  }else if(!is.null(RowSideColors)){
    a <- d3heatmap(data, main=title, trace='none', key=TRUE, col=color,
              dendrogram=dendrogram, na.rm=na.rm,
              Colv=Colv, labCol=labCol,
              Rowv=Rowv, labRow=labRow, RowSideColors=RowSideColors, RowSideColorsSize=RowSideColorsSize,
              scale=scale, lwid=lwid, lhei=lhei, cexCol=cexCol, cexRow=cexRow,
              margins=c(height, width), show_grid=show_grid,
              anim_duration=anim_duration )
  }else if(!is.null(ColSideColors)){
    a <- d3heatmap(data, main=title, trace='none', key=TRUE, col=color,
              dendrogram=dendrogram, na.rm=na.rm,
              Colv=Colv, labCol=labCol, ColSideColors=ColSideColors, ColSideColorsSize=ColSideColorsSize,
              Rowv=Rowv, labRow=labRow,
              scale=scale, lwid=lwid, lhei=lhei, cexCol=cexCol, cexRow=cexRow,
              margins=c(height, width), show_grid=show_grid,
              anim_duration=anim_duration )
  }else{
    a <- d3heatmap(data, main=title, trace='none', key=TRUE, col=color,
              dendrogram=dendrogram, na.rm=na.rm,
              Colv=Colv, labCol=labCol,
              Rowv=Rowv, labRow=labRow,
              scale=scale, lwid=lwid, lhei=lhei, cexCol=cexCol, cexRow=cexRow,
              margins=c(height, width), show_grid=show,
              anim_duration=anim_duration )
  }
  return(a)
}
