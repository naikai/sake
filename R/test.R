library(shiny)
library(plotly)

adds <- mtcars[1:5, ]
adds[, c(1, 6)] <- NA
rownames(adds) <- paste(rownames(adds), "(NA)")
some_null_cars <- rbind(mtcars, adds)

ui <- fluidPage(
  radioButtons("plotType", "Plot Type:", choices = c("ggplotly", "plotly")),
  plotlyOutput("plot"),
  verbatimTextOutput("brush")
)

server <- function(input, output, session) {

  output$plot <- renderPlotly({
    # use the key aesthetic/argument to help uniquely identify selected observations
    key <- row.names(some_null_cars)
    if (identical(input$plotType, "ggplotly")) {
      p <- ggplot(some_null_cars, aes(x = mpg, y = wt, colour = factor(vs), key = key)) +
        geom_point() + theme_bw()
      ggplotly(p) %>% layout(dragmode = "select")
    } else {
      plot_ly(some_null_cars, x = mpg, y = wt, key = key, mode = "markers", type = "scattergl", source = "scatter")
      # %>%
      #   layout(dragmode = "select")
    }
  })

  output$brush <- renderPrint({
    d <- event_data("plotly_selected")
    if (is.null(d)) "Click and drag events (i.e., select/lasso) appear here (double-click to clear)" else d
  })
}

shinyApp(ui, server)
