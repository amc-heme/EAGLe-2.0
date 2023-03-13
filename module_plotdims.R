library(shiny)

sliderUI <- function(id, min, max, value) {
  ns <- NS(id)
  
  sliderInput(ns("slider"), "Slider",
              min = min, max = max, value = value
  )
}

sliderServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    sliderscale <- reactive({
      input$slider
    })
    
    return(sliderscale)
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
