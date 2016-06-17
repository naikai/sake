#' @export
featureUI <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(width=6,
           box(title="Sample Correlation", width=NULL, solidHeader=TRUE, status="info", height = "600px",
               fluidRow(
                 column(width=4, numericInput(ns("cor_sam_lab_cex"), label = "text label size", min=0.5, max=1, value=0.7, step = 0.1)),
                 column(width=4, numericInput(ns("cor_num_lab_cex"), label = "cor label size", min=0, max=1, value=0.4, step = 0.1)),
                 column(width=2, br(),
                        actionButton(ns("runSamCor"), " Plot!  ", icon("play-circle"),
                                     style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                 )
               ),
               fluidRow(
                 column(width=12, corModuleUI(ns("sample")))
               )
           )
    ),
    column(width=6,
           box(title="Gene Network", width=NULL, solidHeader=TRUE, status="info", height = "600px",
               fluidRow(
                 column(width=2, br(),
                        actionButton(ns("runGeneCor"), " Plot!  ", icon("play-circle"),
                                     style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                 )
               ),
               fluidRow(
                 column(width=12, corModuleUI(ns("gene")))
               )
           )
    )
  )
}

#' @export
feature <- function(input, output, session, data){
  ns <- session$ns
  callModule(corModule, ns("sample"), data = reactive({ t(data()) }),
                          tl_cex = reactive(input$cor_sam_lab_cex-0.3),
                          number_cex = reactive(0.001))
}

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
                      tl_cex = 0.2, number_cex = 0.01, type = "upper",
                      diag = FALSE, height = 500){
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

