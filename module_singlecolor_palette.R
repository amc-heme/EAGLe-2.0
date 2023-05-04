library(shiny) 
library(colourpicker)
library(esquisse)

colorUI <- function(id, label, value) {
  ns <- NS(id)
  
  colourInput(
    ns("colorpalette"),
    label = label,
    value = value,
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
      print("value of color palette within color module")
      print(input$colorpallete)
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
