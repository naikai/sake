#' @export
featureUI <- function(id) {
  ns <- NS(id)
  fluidRow(
    column(width=6,
           box(title="Sample Correlation", width=NULL, solidHeader=TRUE, status="info", height = "620px",
               fluidRow(
                 column(width=3, numericInput(ns("cor_sam_lab_cex"), label = "text label size", min=0.5, max=1, value=0.7, step = 0.1)),
                 column(width=3, numericInput(ns("cor_num_lab_cex"), label = "cor label size", min=0, max=1, value=0.4, step = 0.1)),
                 column(width=3, selectInput(ns("cor_type"), label = "Plot type", choices = c("upper", "lower", "full"), selected = "full")),
                 column(width=2, br(),
                        actionButton(ns("runSamCor"), " Plot!  ", icon("play-circle"),
                                     style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                 )
               ),
               fluidRow(
                 bsModal(ns("modalExample"), "Warning! Your sample size if above 200, it will take longer than usual. Are you sure you want to continue?",
                         trigger = "", size = "small",
                         actionButton(ns("cor_forcego_yes"), 'Yes'),
                         actionButton(ns("cor_forcego_no"), 'No')
                         )
               ),
               fluidRow(
                 # column(width=12, corModuleUI(ns("sample")))
                 column(width=12,  plotOutput(ns('sampleCorPlot')))
               ),
               fluidRow(
                 column(width=12, bsAlert("cor_alert"))
               )
           )
    ),
    column(width=6,
           box(title="Gene Network", width=NULL, solidHeader=TRUE, status="info", height = "620px",
               fluidRow(
                 column(width=2, br(),
                        actionButton(ns("runGeneCor"), " Plot!  ", icon("play-circle"),
                                     style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                 )
               ),
               fluidRow(
                 column(width=12, corModuleUI(ns("gene")))
               ),
               fluidRow(
                 bsModal(ns("modalExample2"), "Warning! Your sample size if above 200, it will take longer than usual. Are you sure you want to continue?", ns("runGeneCor"), size = "small",
                         actionButton(ns("corgene_forcego"), 'Run'))
               )
           )
    )
  )
}

#' @export
feature <- function(input, output, session, data){
  go <- reactiveValues(run=FALSE)

  observeEvent(input$runSamCor, {
    # if(ncol(data()) > 200){
    if(nrow(data()) > 200){
      print('greater than 200')
      toggleModal(session, "modalExample", toggle = "open")
      # createAlert(session, "cor_alert", "cor_alert1", title = "Warning", style = "warning",
      #             content = "Data contains more than 200 samples, It might take longer than usual to run. <br>
      #             Please make sure you are aware of this and still want to run",
      #             append = TRUE)
    }else{
      go$run <- TRUE
    }
  })

  observeEvent(input$cor_forcego_yes, {
    print('force go yes')
    go$run <- TRUE
    toggleModal(session, "modalExample", toggle = "close")
  })
  observeEvent(input$cor_forcego_no, {
    print('force go no')
    toggleModal(session, "modalExample", toggle = "close")
  })

  observe({
    if(go$run == FALSE) return()

    isolate({
      print('run')
      print(paste('gorun:', go$run))

      # tl_cex <- reactive(input$cor_sam_lab_cex)
      # number_cex <- reactive(input$cor_num_lab_cex)
      # type <- reactive(input$cor_type)
      tl_cex <- input$cor_sam_lab_cex
      number_cex <- input$cor_num_lab_cex
      type <- input$cor_type

      # callModule(corModule, "sample", reactive({ data() }),
      #            tl_cex = tl_cex,
      #            number_cex = number_cex,
      #            type = type)

      diag <- FALSE
      output$sampleCorPlot <- renderPlot({
        if(type == "full"){
          diag = TRUE
        }
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
                   tl.col="black", tl.srt=45, tl.cex=tl_cex, #Text label color and rotation
                   number.cex = number_cex,
                   p.mat = p.mat, sig.level = 0.01, insig = "blank",
                   diag=diag
          )
        })
      }, height=500)

      go$run <- FALSE
    })
  })

  # observeEvent(input$runGeneCor, {
  #   callModule(corModule, "gene", reactive({ t(data()) }),
  #              tl_cex = reactive(input$cor_sam_lab_cex-0.3),
  #              number_cex = reactive(0.001),
  #              type = reactive(input$cor_type))
  # })
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
                      tl_cex, number_cex, type,
                      diag = FALSE, height = 480){
  if(type() == "full"){
    diag = TRUE
  }
  print('running cormodule .. ')

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
               type=type(), order="hclust",
               addCoef.col = "black", # Add coefficient of correlation
               tl.col="black", tl.srt=45, tl.cex=tl_cex(), #Text label color and rotation
               number.cex = number_cex(),
               p.mat = p.mat, sig.level = 0.01, insig = "blank",
               diag=diag
      )
    })
  }, height=height)
}

