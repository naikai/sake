#' Sample correlation plot
#' @export
corModuleUI <- function(id) {
  ns <- NS(id)
  tagList(
    plotOutput(ns('sampleCorPlot'))
  )
}

#' @export
corModule <- function(input, output, session, data,
                      tl_cex=0.2, number_cex=0.01, type="upper",
                      diag=FALSE, height=500){
  output$sampleCorPlot <- renderPlot({
    n <- 2
    withProgress(message = 'Calculating correlation', value = 0, {
      incProgress(1/n, detail = "Takes around 10 seconds")
      M <- cor(data())
      p.mat <- cor_mtest(data())
    })
    col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA", "#79AEDD", "#FFFFFF", "#2E9988", "#2B4444"))

    withProgress(message = 'Plotting..', value = 0, {
      incProgress(1/n, detail = "Takes around 10-20 seconds")
      corrplot(M, method="color", col=rev(col(200)),
               type=type, order="hclust",
               addCoef.col = "black", # Add coefficient of correlation
               tl.col="black", tl.srt=45, tl.cex=tl_cex(), #Text label color and rotation
               number.cex = number_cex(),
               p.mat = p.mat, sig.level = 0.01, insig = "blank",
               diag=diag
      )
    })
  }, height=height)
}

