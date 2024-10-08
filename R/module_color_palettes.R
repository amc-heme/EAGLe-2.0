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
                    "viridis" = viridis_pal(option = "viridis")(12),
                    "magma" = viridis_pal(option = "magma")(12),
                    "mako" = viridis_pal(option = "mako")(12),
                    "plasma" = viridis_pal(option = "plasma")(12),
                    "cividis" = viridis_pal(option = "cividis")(12)
                  )
                ))
}

paletteServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    color <- reactive({
      input$palette
    })
    
    return(color)
  })
}


paletteApp <- function() {
  ui <- fluidPage(paletteUI("palette"))
  server <- function(input, output, session) {
    paletteServer("palette")
  }
  shinyApp(ui, server)
}
paletteApp()
