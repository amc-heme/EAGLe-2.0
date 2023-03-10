library(shiny) 
library(colourpicker)
library(circlize)
library(ggsci)
library(scales)
library(esquisse)

paletteUI <- function(id, choices) {
  ns <- NS(id)

    palettePicker(ns("palette"), 
       "Choose a color palette",
      choices = list(
        "Viridis" = list(
          "viridis" = viridis_pal(option = "viridis")(5),
          "magma" = viridis_pal(option = "magma")(5),
          "mako" = viridis_pal(option = "mako")(5),
          "plasma" = viridis_pal(option = "plasma")(5),
          "cividis" = viridis_pal(option = "cividis")(5))
      )
    )
}

paletteServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    color <- reactive(scale_color_viridis_d(option = (input$palette)))
    
  })
}


paletteApp <- function() {
  ui <- fluidPage(
    paletteUI("color1")
  )
  server <- function(input, output, session) {
    paletteServer("color1")
  }
  shinyApp(ui, server)  
}
paletteApp()
