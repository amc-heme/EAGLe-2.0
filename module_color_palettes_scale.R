library(shiny) 
library(colourpicker)
library(circlize)
library(ggsci)
library(scales)
library(esquisse)

colorUI <- function(id) {
  ns <- NS(id)
  tagList(
    palettePicker(ns("Palette"), 
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
    
  )
}

colorServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    eventReactive(get(input$PaletteChoicesQC), {
      if(input$PaletteChoicesQC == "viridis") {
        scale_fill_viridis_d(option = "viridis")
      } else if(input$PaletteChoicesQC == "cividis") {
        scale_fill_viridis_d(option = "cividis")
      } else if(input$PaletteChoicesQC == "magma") {
        scale_fill_viridis_d(option = "magma")
      } else if(input$PaletteChoicesQC == "plasma") {
        scale_fill_viridis_d(option = "plasma")
      }else if(input$PaletteChoicesQC == "inferno") {
        scale_fill_viridis_d(option = "inferno")
      }
    }) 
  })
}

colorApp <- function() {
  ui <- fluidPage(
    colorUI("color2")
  )
  server <- function(input, output, session) {
    colorServer("color2")
  }
  shinyApp(ui, server)  
}
colorApp()
