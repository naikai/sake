#' @export
networkUI <- function(id, title) {
  ns <- NS(id)
  box(title=title, width=NULL, solidHeader=TRUE, status="info", height = "700px",
      fluidRow(
        uiOutput(ns("ui_coropt1")),
        uiOutput(ns("ui_coropt2")),
        uiOutput(ns("ui_coropt3"))
      ),
      fluidRow(
        uiOutput(ns("ui_cordownload"))
      ),
      fluidRow(
        bsModal(ns("modal1"), "Warning!\nYour gene size is above 1000, it will take longer than usual.\nAre you sure you want to continue?",
                trigger = "", size = "small",
                actionButton(ns("cor_forcego_yes"), 'Yes'),
                actionButton(ns("cor_forcego_no"), 'No')
        ),
        bsModal(ns("modal2"), "Stop!\nThis module will not work for gene size over 2000, please select a smaller number.",
                trigger = "", size = "small"
        )
      ),
      fluidRow(
        uiOutput(ns("ui_network"))
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
        column(width=2, br(), actionButton(ns("runSamCor"), " Plot!  ", icon("play-circle"), class = 'act'))
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
    if(nrow(data()) > 2000){
      toggleModal(session, "modal2", toggle = "open")
    }else if(nrow(data()) > 1000){
      toggleModal(session, "modal1", toggle = "open")
    }else{
      go$run <- TRUE
    }
  })

  observeEvent(input$cor_forcego_yes, {
    go$run <- TRUE
    toggleModal(session, "modal1", toggle = "close")
  })
  observeEvent(input$cor_forcego_no, {
    toggleModal(session, "modal1", toggle = "close")
  })

  observe({
    if(go$run == FALSE) return()

    isolate({
          M <- t(data()) %>% WGCNA::cor(., nThreads = 4) %>% abs
          diag(M) <- 0
          networkData <- NULL
          Links <- NULL

          if(input$network_type == "Simple"){
            withProgress(message = 'Calculating correlation and network info', value = 0, {
              incProgress(1/2, detail = "Takes around 10 ~ 20 seconds")
              net_dist <- (M > input$tao) * 1
              networkData <- net_dist %>%
                              reshape2::melt(.) %>%
                              dplyr::filter(value==1) %>%
                              dplyr::select(-value)
              print('dim networkData')
              print(dim(networkData))
            })

            output$simplenetwork <- renderSimpleNetwork({
              if(input$stop_network) return()
              validate(
                need(!is.null(networkData), "Network Data not available") %then%
                  need(nrow(networkData) > 10, "Number of network connections is not enough") %then%
                  need(nrow(networkData) < 20000, "Number of network connections is too high (greater than 20000)")
              )
              withProgress(message = 'Plotting..', value = NULL, {
                simpleNetwork(networkData, opacity = isolate(input$opacity), zoom = TRUE)
              })
            })

            output$ui_network <- renderUI({
              tagList(
                column(width=12, simpleNetworkOutput(ns("simplenetwork")))
              )
            })
          }else if(input$network_type == "Force"){
            withProgress(message = 'Calculating correlation and network info', value = 0, {
              incProgress(1/4, detail = "Transformation")
              cutoff <- 0.3
              net_dist <- M ^ input$power
              idx <- rowSums(net_dist>cutoff) != 0
              validate(
                need(sum(idx) > 0, "Did not find enough gene-gene connection from parameters you specify")
              )
              net_dist <- net_dist[idx, idx]

              incProgress(2/4, detail = "Defining Nodes")
              Nodes = data.frame(name = rownames(net_dist),
                                 group = 1,
                                 size = rowSums(net_dist),
                                 row.names = NULL)

              incProgress(3/4, detail = "Defining links between nodes")
              Links <- net_dist %>%
                        reshape2::melt(.) %>%
                        dplyr::filter(value>=cutoff) %>%
                        set_colnames(c("source", "target", "value")) %>%
                        dplyr::mutate(source = match(.$source, Nodes$name) - 1) %>%
                        dplyr::mutate(target = match(.$target, Nodes$name) - 1)
              print('dim Links')
              print(dim(Links))
            })

            output$forcenetwork <- renderForceNetwork({
              if(input$stop_network) return()
              validate(
                need(!is.null(Links), "Network Data not available") %then%
                need(nrow(Links) > 10, "Number of network connections is not enough") %then%
                  need(nrow(Links) < 20000, "Number of network connections is too high (greater than 20000)")
              )
              withProgress(message = 'Plotting..', value = NULL, {
                forceNetwork(Links = Links, Nodes = Nodes,
                             Source = "source", Target = "target",
                             Value = "value", NodeID = "name",
                             Nodesize = "size", Group = "group",
                             opacity = isolate(input$opacity), opacityNoHover = 0.3,
                             zoom = TRUE)
              })
            })
            output$ui_network <- renderUI({
              tagList(
                column(width=12, forceNetworkOutput(ns("forcenetwork")))
              )
            })
          }

      output$ui_cordownload <- renderUI({
        tagList(
          column(width=3, br(),
                 downloadButton(ns("dl_network"), "Download", class="dwnld")
          ),
          column(width=2, br(),
                 actionButton(ns("stop_network"), "Stop", class="stop")
          )
        )
      })

      output$dl_network <- downloadHandler(
        filename <- function() {
          paste("network.html")
        },
        content = function(file) {
          if(input$network_type == "Simple"){
            simpleNetwork(networkData, zoom = TRUE, opacity = input$opacity) %>%
              saveNetwork(file = file)
          }else if (input$network_type == "Force"){
            forceNetwork(Links = Links, Nodes = Nodes,
                         Source = "source", Target = "target",
                         Value = "value", NodeID = "name",
                         Nodesize = "size", Group = "group",
                         opacity = isolate(input$opacity), opacityNoHover = 0.3,
                         zoom = TRUE) %>%
              saveNetwork(file = file)
          }
        }
      )
      go$run <- FALSE
    })
  })
}


