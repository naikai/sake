#' A Cat Function
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @export
#' @examples
#' cat_function()
# Compare two grouping result
compare_groups <- function(group1, group2, file.prefix="Compare.group", title="", plot=T, save.image=F, add.legend=T, label.size=5, title.size=36){
	if (class(group1)=="integer" | class(group1)=="numeric")
	group1 <- paste0("Group", group1)
	if (class(group2)=="integer" | class(group2)=="numeric")
	group2 <- paste0("Group", group2)

	### Calculate the percentge first before feed into ggplot2 ### 
	### This is based on group1. 
	summary <- as.matrix(table(group1, group2))
	sums <- rowSums(summary)
	summary.pct <- apply(t(summary), 1, function(x) x/sums) * 100
	summary.pct <- melt(summary.pct)

	### Stacked Bar using ggplot2
	ray <- ggplot(data=summary.pct, aes(x=group1, y=value, fill=group2)) + 
				geom_bar(stat="identity") + labs(y = "Percentage (%)") + 
				scale_fill_brewer(palette="Set1") + 
				xlab('') + labs(fill="")
	### This part is specific for manually colorign the groups to match pam50 coloring 
	# xlab('') + scale_fill_manual(values=c("#E41A1C", "#fccde5", "#2166ac", "#a6cee3", "#33A02C")) # NMFK5-vs.PAM50
	# xlab('') + scale_fill_manual(values=c("#33A02C", "#fccde5", "#a6cee3", "#2166ac", "#E41A1C")) # PAM50-vs-NMFK5
	
	### Add the percentage within each bar ### 
	### Get ggplot data
	naikai <- ggplot_build(ray)$data[[1]]
	### Create values for text positions. 
	naikai$position = (naikai$ymax + naikai$ymin)/2
	### round up numbers and convert to character. 
	foo <- round(summary.pct[order(summary.pct$group1), "value"], digits=2)
	### Create a column for text
	naikai$label <- paste0(foo, "%")
	### Plot again
	ray <- ray + annotate(x = naikai$x, y = naikai$position, label = naikai$label, geom = "text", size=label.size) 
	ray <- ray + ggtitle(title) + theme(plot.title=element_text(face="bold", size=title.size))

	# legend 
	if (!add.legend){
		ray <- ray + theme(legend.position="none")
	}
	
	# option to plot it or not
	if(plot)
		print(ray)

	# option to save the plot or not 
	if (save.image){
		pdf(paste0(file.prefix, ".pdf"), height=10, width=12)
		print(ray)
		dev.off()
	}

	# or just simply return the plot 
	return (ray)
}


