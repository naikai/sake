#' Adjust label cex based on its number
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @export
#' @examples
#' assign_label_cex(100)
assign_label_cex <- function(n.samples){
   if (n.samples < 40){
      label.cex <- 1.2
   }else if(n.samples < 70){
      label.cex <- 1.0
   }else if(n.samples < 100){
      label.cex <- 0.7
   }else if (n.samples < 200){
      label.cex <- 0.6
   }else if (n.samples < 300){
      label.cex <- 0.5
   }else{
      label.cex <- 0.4
   }

   return(label.cex)
}


#' load predefined gene list files for heatmap
#'
#' Load predefined gene list files for heatmap (clustering) result
#' @param folder Where the files are
#' @param pattern file extensions for the ones you want to load
#' @keywords lodad
#' @export
#' @examples
#' cor_mtest("data", ".txt")
cor_mtest <- function(mat, ...) {
    mat <- as.matrix(mat)
    n <- ncol(mat)
    p.mat<- matrix(NA, n, n)
    diag(p.mat) <- 0
    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            tmp <- cor.test(mat[, i], mat[, j], ...)
            p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
        }
    }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}


#' Create brewer color
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @export
#' @examples
#' create.brewer.color()
create.brewer.color <- function(data, num=8, name="Set1")
{
	if(name=="naikai"){
		my_pallete <- c(
                        rgb(236,67,35, maxColorValue=255),
                        rgb(56,146,208, maxColorValue=255),
                        rgb(230,235,88, maxColorValue=255),
                        rgb(116,187,88, maxColorValue=255),
                        rgb(196,95,46, maxColorValue=255),
                        rgb(203,77,202, maxColorValue=255),
                        rgb(118,220,65, maxColorValue=255),
                        rgb(115,113,206, maxColorValue=255))
		return(create.manual.color(data, my_pallete))
	}else if(name=="naikai2"){
		my_pallete <- c(
                        rgb(10,10,10, maxColorValue=255),
                        rgb(230,235,88, maxColorValue=255),
                        rgb(118,220,65, maxColorValue=255),
                        rgb(203,77,202, maxColorValue=255),
                        rgb(196,95,46, maxColorValue=255),
                        rgb(116,187,88, maxColorValue=255),
                        rgb(236,67,35, maxColorValue=255),
                        rgb(56,146,208, maxColorValue=255),
                        rgb(115,113,206, maxColorValue=255))
		return(create.manual.color(data, my_pallete))
	}else{
		groupCodes <- as.factor(data)
		uniq_num <- length(levels(groupCodes))
		colorCodes <- colorRampPalette(brewer.pal(num, name))(uniq_num)
		color.idx <- match(groupCodes, levels(groupCodes))
		label.color <- colorCodes[color.idx]
		return(label.color)
	}
}



#' Manually specify coloring for provided groups
#'
#' Assign colors to data, if more unique data than provided colors, will impute the missing ones and fill them in
#' @param data data
#' @param group.color colors for each factorial group
#' @keywords color manual
#' @export
#' @examples
#' create.manual.color(c(1,2,3,1,2,3,2,3), c("red", "blue", "green"))
create.manual.color <- function(data, group.color)
{
	num_uniq_data <- length(unique(data))
	if(num_uniq_data > length(group.color)){
		require(RColorBrewer)
		group.color <- colorRampPalette(group.color)(num_uniq_data)
	}
	data <- as.numeric(factor(data))
	return(group.color[data])
}


#' Custom version of boxplot
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @export
#' @examples
#' myboxplot()
myboxplot <- function(data, prefix="Boxplot", method="vioplot", title="",
                      ylim=c(0.6,0.85), save.plot=F,
                      axis.text.size=14, axis.title.size=16,
                      legend.position="right", vio.add.method="dot")
{
  data.m <- reshape2::melt(data)
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


#' Convert name to color
#'
#' This function convert the sample_name in the data into different colors
#' @param name usually sample name, separate with "_"
#' @keywords color
#' @export
#' @examples
#' name_to_color(c("group1", "group1", "group2", "group2"))
name_to_color <- function(name, split_pattern="\\_", num_color=1,
						  ColScheme=c("naikai", "naikai2", "Set1", "Set2", "Set3", "Dark2", "Pastel1", "Pastel2", "Paired", "Accent")
						 )
{
	col.color <- list()
	ColSideColors <- NULL
	col.color.name <- list()
	ColSideColors.name <- NULL

	column.names <- strsplit(name, split=split_pattern)
	max_num_breaks <- max(sapply(column.names, length))

	if(num_color != length(ColScheme)){
		warning(paste0("Error: num_color:", num_color, "is different than num of ColScheme:", length(ColScheme), "\n"))
	  ColScheme <- ColScheme[1:num_color]
	}
	if(max_num_breaks > num_color){
		warning("There are more possible breaks in the name than the num_color specified. We will use only the num_color")
	}

	for(i in 1:num_color){
		num_names <- sapply(column.names, function(x) x[i]) %>% unique %>% length
		col.color[[i]] <- create.brewer.color(sapply(column.names, function(x) x[i]), num_names, ColScheme[i])
		ColSideColors <- cbind(ColSideColors, col.color[[i]])
		col.color.name[[i]] <- sapply(column.names, function(x) x[i])
		ColSideColors.name <- cbind(ColSideColors.name, col.color.name[[i]])
	}

	res <- list()
	if(num_color==1){
		if(length(unique(as.character(ColSideColors[,1])))==1){
			res[["color"]] <- NULL
			res[["name"]] <- NULL
		}else{
			res[["color"]] <- as.matrix(ColSideColors[, 1])
			res[["name"]] <- as.matrix(ColSideColors.name[, 1])
		}
	}else if(num_color>1){
		num_color <- ifelse(num_color > max_num_breaks, max_num_breaks, num_color)
		res[["color"]] <- as.matrix(ColSideColors[, 1:num_color])
		res[["name"]] <- as.matrix(ColSideColors.name[, 1:num_color])
	}else{
		stop("Error: num_color must be greater than 0")
	}

	return(res)
}




#' Plot expression for each gene
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @export
#' @examples
#' cat_function()
plot_gene_expression <- function (data, gene, groups, type="boxplot", data.return=f, ylim=c(0,20), title=null, las=1, cex.axis=1){
	data <- as.numeric(data[rownames(data) == gene, ])
	data <- data.frame(groups=factor(groups), Expression=data)
	if (is.null(title)){
		title <- paste("Gene:", gene)
	}

	if (type=="boxplot"){
		p <- ggplot(data, aes(x=groups, y=Expression)) + geom_boxplot(aes(fill = factor(groups)))
		p <- p + ggtitle(title) + theme_bw() + ylim(ylim) + xlab("")
		print(p)
	}else if (type=="base"){
		par(cex.axis=cex.axis)
		boxplot(Expression~groups, data=data, main=title, xlab="Group", ylab="Expression", las=las)
	}

	if(data.return)
	return(p)
}

#' simple PCA plot
#'
#' This is a simple version of PCA plot
pcaplot <- function(data, takelog=FALSE, center=TRUE, scale=FALSE, pca_x=1, pca_y=2, group=NULL, opacity=0.6, plotsize=7, labsize=4, legend=TRUE){
  # remove genes with zero read
  idx <- which(rowSums(data) == 0)
  if (length(idx>0)){
    data <- data[-idx, ]
  }

  if(takelog){ 
    if(any(data < 0)){
      warning("ERROR: data contain negative values. Can not perform log transformation")
    }else{
      data <- log2(data + 0.1) 
    }
  }

  if(!is.null(group)){
    if(length(group) == ncol(data)){
      # sort the data according to order of the group
      idx <- order(group)
      group <- group[idx]
      data <- data[, idx]
      color = create.brewer.color(group, length(unique(group)), "naikai")
    }else{
      warning("ERROR: length of group is not the same as number of columns in the input data")
      group = "D"
      legend = FALSE
    }
  }else{
    group = "D"
    color = "black"
    legend = FALSE
  }

  #pc <- prcomp(t(data), center=center, scale.=scale)
  pc <- prcomp(t(data), center=center, scale.= FALSE)
  vars <- pc$sdev^2
  vars <- signif(vars/sum(vars) * 100, 2)
  idx <- c(pca_x, pca_y)

  projection <- as.data.frame(pc$x) %>% 
    dplyr::select(idx) %>% 
    set_colnames(c("pca_x", "pca_y")) %>% 
    mutate(group = group, 
           color = color,
           Sample = rownames(pc$x)) %>%
    arrange(group)

  title = paste(paste0("PC", idx), paste0("(", vars[idx], "% Var)"))

  t <- list( size=labsize, color=toRGB("grey50") )
  p <- plot_ly(projection, x = ~pca_x, y = ~pca_y, type="scatter", mode="markers",
   color = ~group, colors = ~unique(color), 
   text = ~paste(Sample, '</br>Group:', group), hoverinfo="text", textposition="top middle",
   marker = list(size=plotsize, opacity=opacity))

  p %>% layout(showlegend = legend,
   xaxis = list(title = title[pca_x], zeroline=FALSE),
   yaxis = list(title = title[pca_y], zeroline=FALSE))
}


#' Reorder dendrogram based on mean expression for RNA-Seq
#' @param d dendrogram
#' @param w dendrogram
#' @keywords reorder
#' @export
#' @examples
#' reorderfun()
reorderfun <- function(d, w) reorder(d, w)


#' Screen shot of plotly object
#'
#' @param x a plotly object
#' @param file output filename
#' @export
screenshot_plotly = function(x, file = 'webshot.png', ...) {
  file2 = normalizePath(file, mustWork = FALSE)
  d = tempfile()
  dir.create(d)
  owd = setwd(d)
  on.exit({
    setwd(owd)
    unlink(d, recursive = TRUE)
  }, add = TRUE)
  htmlwidgets::saveWidget(plotly::as.widget(x), 'index.html', FALSE)
  file.copy(webshot::webshot('index.html', ...), file2)
  file
}
