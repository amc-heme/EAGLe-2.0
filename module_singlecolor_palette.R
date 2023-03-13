library(shiny) 
library(colourpicker)
library(esquisse)

colorUI <- function(id, choices) {
  ns <- NS(id)
  
  colourInput(
    ns("colorpalette"),
    label = "Choose a color",
    value = NULL,
    showColour = ("both"),
    palette = ("square"),
    allowedCols = NULL,
    allowTransparent = FALSE,
    returnName = FALSE,
    closeOnClick = FALSE
  )
}

colorServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    singlecolor <- reactive({
      input$colorpalette
    })
    
    return(singlecolor)
  })
}


colorApp <- function() {
  ui <- fluidPage(colorUI("color"))
  server <- function(input, output, session) {
    colorServer("color")
  }
  shinyApp(ui, server)
}
colorApp()
