library(shiny) 
library(colourpicker)
library(circlize)
library(colourpicker)
library(ggsci)
library(scales)
library(esquisse)
paletteUI <- function(id) {
  tagList(
    palettePicker(
      inputId = (NS(id,"PaletteChoicesQC")),
      label = "Choose a color palette",
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

paletteServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    colorpalettechoices <-
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
    #PCA plot color palette choices for color
    colorchoicesQC <-
      eventReactive(input$PaletteChoicesQC, {
        if(input$PaletteChoicesQC == "viridis") {
          scale_color_viridis_d(option = "viridis")
        } else if(input$PaletteChoicesQC == "cividis") {
          scale_color_viridis_d(option = "cividis")
        } else if(input$PaletteChoicesQC == "magma") {
          scale_color_viridis_d(option = "magma")
        } else if(input$PaletteChoicesQC == "plasma") {
          scale_color_viridis_d(option = "plasma")
        }else if(input$PaletteChoicesQC == "inferno") {
          scale_color_viridis_d(option = "inferno")
        }
      })
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
