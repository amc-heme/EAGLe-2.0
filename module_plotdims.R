library(shiny)
library(shinyjs)
sliderUI <- function(id, min, max, value, label) {
  ns <- NS(id)
  
  shinyjs::useShinyjs()
  hidden(
  sliderInput(ns("slider"), label = label,
              min = min, max = max, value = value
    )
  )
}
switchUI <- function(id, label, value, right){
  ns <- NS(id)
  materialSwitch(ns("hidedims"), label = label, value = value, right = right)
}

sliderServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    sliderscale <- reactive({
      input$slider
    })
    return(sliderscale)
  })
}

switchServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    observe({
      shinyjs::toggle("slider", condition = input$hidedims)
    })
  })
}

sliderApp <- function() {
  ui <- fluidPage(sliderUI("slider"))
  server <- function(input, output, session) {
    sliderServer("slider")
  }
  shinyApp(ui, server)
}
sliderApp()
