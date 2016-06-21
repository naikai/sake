#' @export
networkUI <- function(id, title) {
  ns <- NS(id)
  box(title=title, width=NULL, solidHeader=TRUE, status="info", height = "700px",
      fluidRow(
        uiOutput(ns("ui_coropt1")),
        uiOutput(ns("ui_coropt2"))
      ),
      fluidRow(
        bsModal(ns("modalExample"), "Warning!\nYour gene size is above 1000, it will take longer than usual.\nAre you sure you want to continue?",
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
        column(width=12, simpleNetworkOutput(ns("simplenetwork")) )
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
        column(width=3, numericInput(ns("tao"), label = "Hard threshold", min=0.5, max=1, value=0.8, step = 0.05)),
        column(width=3, numericInput(ns("opacity"), label = "Opacity", min=0.05, max=1, value=0.8, step = 0.05))
      )
    })
    output$ui_coropt2 <- renderUI({
      tagList(
        column(width=2, br(),
               actionButton(ns("runSamCor"), " Plot!  ", icon("play-circle"),
                            style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
        )
      )
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
      if(nrow(data()) > 1000){
        tl_cex <- 0.005
        number_cex <- 0.001
      }else{
        tl_cex <- input$cor_sam_lab_cex
        number_cex <- input$cor_num_lab_cex
      }

      type <- input$cor_type
      # diag <- ifelse(type=="full", TRUE, FALSE)
      col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA", "#79AEDD", "#FFFFFF", "#2E9988", "#2B4444"))

      withProgress(message = 'Calculating correlation', value = 0, {
        incProgress(1/2, detail = "Takes around 20 ~ 30 seconds")
        # M <- cor(data())
        method <- "pearson"

        M <- t(data()) %>% WGCNA::cor(., method=method, nThreads = 4) %>% abs
        diag(M) <- 0
        net_dist <- (M > input$tao) * 1
        networkData <- net_dist %>% melt %>% dplyr::filter(value==1) %>% dplyr::select(-value)
        networkData <- networkData[1:(nrow(networkData)/2), ]
        print(dim(networkData))
      })

      lett <- function(num){
        return(sample(LETTERS, num, replace = TRUE))
      }

      output$simplenetwork <- renderSimpleNetwork({
        withProgress(message = 'Plotting..', value = NULL, {
          validate(
            need(nrow(networkData) > 10, "Number of network connections is not enough")
          )
          # networkData <- cbind(
          #   expand.grid(lett(30), 1:30) %>% tidyr::unite(src, Var1, Var2, remove = TRUE, sep = ""),
          #   expand.grid(lett(30), 1:30) %>% tidyr::unite(target, Var1, Var2, remove = TRUE, sep = "")
          # )
          # print(dim(networkData))
          simpleNetwork(networkData, opacity = input$opacity, zoom = TRUE)
          # simpleNetwork(networkData, zoom = FALSE, opacity = input$opacity)
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
            saveNetwork(file = 'network.html')
        }
      )
      go$run <- FALSE
    })
  })
}


