#' @export
featureUI <- function(id, title) {
  ns <- NS(id)
  box(title=title, width=NULL, solidHeader=TRUE, status="info", height = "700px",
      fluidRow(
        uiOutput(ns("ui_coropt1")),
        uiOutput(ns("ui_coropt2"))
      ),
      fluidRow(
        bsModal(ns("modalExample"), "Warning!\nYour sample size is above 200, it will take longer than usual.\nAre you sure you want to continue?",
                trigger = "", size = "small",
                actionButton(ns("cor_forcego_yes"), 'Yes'),
                actionButton(ns("cor_forcego_no"), 'No')
        )
      ),
      fluidRow(
        uiOutput(ns("ui_cordownload"))
      ),
      fluidRow(
        # column(width=12, corModuleUI(ns("sample")))
        column(width=12,  plotOutput(ns('sampleCorPlot')))
      )
  )
  # column(width=6,
  #        box(title="Gene Network", width=NULL, solidHeader=TRUE, status="info", height = "700px",
  #            fluidRow(
  #              column(width=2, br(),
  #                     actionButton(ns("runGeneCor"), " Plot!  ", icon("play-circle"),
  #                                  style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
  #              )
  #            ),
  #            fluidRow(
  #              column(width=12, corModuleUI(ns("gene")))
  #            ),
  #            fluidRow(
  #              bsModal(ns("modalExample2"), "Warning! Your sample size is above 200, it will take longer than usual. Are you sure you want to continue?", ns("runGeneCor"), size = "small",
  #                      actionButton(ns("corgene_forcego"), 'Run'))
  #            )
  #        )
  # )
  # )
}

#' @export
feature <- function(input, output, session, data){
  ns <- session$ns
  go <- reactiveValues(run=FALSE)

  observeEvent(data(), {
    output$ui_coropt1 <- renderUI({
      if(ncol(data()) > 200) return(NULL)
      tagList(
        column(width=3, numericInput(ns("cor_sam_lab_cex"), label = "Sample label size", min=0.5, max=1, value=0.4, step = 0.1)),
        column(width=3, numericInput(ns("cor_num_lab_cex"), label = "Corr label size", min=0, max=1, value=0.1, step = 0.1))
      )
    })
    output$ui_coropt2 <- renderUI({
      tagList(
        column(width=3, selectInput(ns("cor_type"), label = "Plot type", choices = c("upper", "lower", "full"), selected = "full")),
        column(width=2, br(),
               actionButton(ns("runSamCor"), " Plot!  ", icon("play-circle"),
                            style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
        )
      )
    })
  })

  observeEvent(input$runSamCor, {
    if(ncol(data()) > 200){
      toggleModal(session, "modalExample", toggle = "open")
    }else{
      go$run <- TRUE
    }
  })

  observeEvent(input$cor_forcego_yes, {
    go$run <- TRUE
    toggleModal(session, "modalExample", toggle = "close")
  })
  observeEvent(input$cor_forcego_no, {
    toggleModal(session, "modalExample", toggle = "close")
  })

  observe({
    if(go$run == FALSE) return()

    isolate({
      # these reactive generate plot before click run
      # tl_cex <- reactive(input$cor_sam_lab_cex)
      # number_cex <- reactive(input$cor_num_lab_cex)
      # type <- reactive(input$cor_type)
      if(ncol(data()) > 200){
        tl_cex <- 0.005
        number_cex <- 0.001
      }else{
        tl_cex <- input$cor_sam_lab_cex
        number_cex <- input$cor_num_lab_cex
      }

      type <- input$cor_type
      diag <- ifelse(type=="full", TRUE, FALSE)
      col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA", "#79AEDD", "#FFFFFF", "#2E9988", "#2B4444"))

      withProgress(message = 'Calculating correlation', value = 0, {
        incProgress(1/2, detail = "Takes around 10 seconds")
        M <- cor(data())
        p.mat <- cor_mtest(data())
      })

      output$sampleCorPlot <- renderPlot({
        withProgress(message = 'Plotting..', value = 0, {
          incProgress(1/2, detail = "Takes around 10-20 seconds")
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

      output$ui_cordownload <- renderUI({
        tagList(
          column(width=4, br(),
                 downloadButton(ns("dl_corrplot"), "Download", class="butt")
          )
        )
      })

      output$dl_corrplot <- downloadHandler(
        filename <- function() {
          paste("corrplot.pdf")
        },
        content = function(file) {
          pdf(file, width=input$estim_pdf_w, height=input$estim_pdf_h)
          corrplot(M, method="color", col=rev(col(200)),
                   type=type, order="hclust",
                   addCoef.col = "black", # Add coefficient of correlation
                   tl.col="black", tl.srt=45, tl.cex=tl_cex, #Text label color and rotation
                   number.cex = number_cex,
                   p.mat = p.mat, sig.level = 0.01, insig = "blank",
                   diag=diag
          )
          dev.off()
        }
      )

      # separate to module later when we figure out how to isolate these reactive expression
      # callModule(corModule, "sample", reactive({ data() }),
      #            tl_cex = tl_cex,
      #            number_cex = number_cex,
      #            type = type)
      go$run <- FALSE
    })
  })
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
                      tl_cex = 0.4, number_cex = 0.5, type = "full",
                      diag = FALSE, height = 480){
  if(type() == "full"){
    diag = TRUE
  }
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

