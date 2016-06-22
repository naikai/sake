#' @export
networkUI <- function(id, title) {
  ns <- NS(id)
  box(title=title, width=NULL, solidHeader=TRUE, status="info", height = "700px",
      fluidRow(
        uiOutput(ns("ui_coropt1")),
        uiOutput(ns("ui_coropt2")),
        uiOutput(ns("ui_coropt3")),
        uiOutput(ns("ui_cordownload"))
      ),
      fluidRow(
        bsModal(ns("modalExample"), "Warning!\nYour gene size is above 1000, it will take longer than usual.\nAre you sure you want to continue?",
                trigger = "", size = "small",
                actionButton(ns("cor_forcego_yes"), 'Yes'),
                actionButton(ns("cor_forcego_no"), 'No')
        )
      ),
      fluidRow(
        column(width=12, simpleNetworkOutput(ns("simplenetwork")))
      )
  )
}

#' @export
network <- function(input, output, session, data){
  ns <- session$ns
  go <- reactiveValues(run=FALSE)

  observeEvent(data(), {
    output$ui_coropt1 <- renderUI({
      tagList(
        column(width=3, selectInput(ns("network_type"),  label = "Network type",
                                    choices = c("Simple", "Force"), selected = "Simple")),
        column(width=3, numericInput(ns("opacity"), label = "Opacity", min=0.05, max=1, value=0.8, step = 0.05))
      )
    })
    output$ui_coropt3 <- renderUI({
      tagList(
        column(width=2, br(),
               actionButton(ns("runSamCor"), " Plot!  ", icon("play-circle"),
                            style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
        )
      )
    })
  })

  observeEvent(input$network_type, {
    output$ui_coropt2 <- renderUI({
      if(input$network_type == 'Simple'){
        column(width=3, sliderInput(ns("tao"), label = "Hard threshold", min=0.7, max=1, value=0.8, step = 0.01, ticks=F))
      }else if(input$network_type == 'Force'){
        column(width=3, sliderInput(ns("power"), label = "Power", min=2, max=10, value=6, step = 1, ticks=F))
      }
    })
  })

  observeEvent(input$runSamCor, {
    if(nrow(data()) > 1000){
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
        withProgress(message = 'Calculating correlation and network info', value = 0, {
          incProgress(1/2, detail = "Takes around 10 ~ 20 seconds")

          M <- t(data()) %>% WGCNA::cor(., nThreads = 4) %>% abs
          diag(M) <- 0
          net_dist <- (M > input$tao) * 1
          networkData <- net_dist %>% reshape2::melt(.) %>% dplyr::filter(value==1)
          networkData <- networkData[1:(nrow(networkData)/2), ]

          if(input$network_type == "simple"){
            networkData <- networkData %>% dplyr::select(-value)
            print(paste('dim networkData', dim(networkData)))
          }else if(input$network_type == "force"){
            Links = networkData
            Nodes = 2
          }
        })

      output$simplenetwork <- renderSimpleNetwork({
        withProgress(message = 'Plotting..', value = NULL, {
          validate(
            need(nrow(networkData) > 10, "Number of network connections is not enough")
          )
          simpleNetwork(networkData, opacity = input$opacity, zoom = TRUE)
        })
      })

      output$ui_cordownload <- renderUI({
        tagList(
          column(width=4, br(),
                 downloadButton(ns("dl_corrplot"), "Download", class="butt")
          )
        )
      })

      output$dl_corrplot <- downloadHandler(
        filename <- function() {
          paste("network.html")
        },
        content = function(file) {
          simpleNetwork(networkData, zoom = TRUE, opacity = input$opacity) %>%
            saveNetwork(file = "network.html")
        }
      )
      go$run <- FALSE
    })
  })
}


