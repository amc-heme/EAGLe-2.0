
library(shiny)

pcaUI <- function(id){

  tagList(
    ns <- shiny::NS(id),
    
    materialSwitch(
        inputId = ns(
          "PCAplots"),
      label =
        "PCA",
      value =
        FALSE,
      right =
        TRUE
    ),

    materialSwitch(
      inputId = ns(
        "PCAscreeplots"),
      label =
        "Scree",
      value =
        FALSE,
      right =
        TRUE
    ),

    materialSwitch(
      inputId = ns(
        "multiqc"),
      label =
        "MultiQC",
      value =
        FALSE,
      right =
        TRUE
    )
    #palette choices for PCA plots
    # palettePicker(
    #   inputId = ns("PaletteChoicesQC"),
    #   label = "Choose a color palette",
    #   choices = list(
    #     "Viridis" = list(
    #       "viridis" = viridis_pal(option = "viridis")(5),
    #       "magma" = viridis_pal(option = "magma")(5),
    #       "mako" = viridis_pal(option = "mako")(5),
    #       "plasma" = viridis_pal(option = "plasma")(5),
    #       "cividis" = viridis_pal(option = "cividis")(5))
    #   )
    # )
  )
}

pcaServer<- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      # colorpalettechoices <-
      #   eventReactive(input$PaletteChoicesQC, {
      #     if(input$PaletteChoicesQC == "viridis") {
      #       scale_fill_viridis_d(option = "viridis")
      #     } else if(input$PaletteChoicesQC == "cividis") {
      #       scale_fill_viridis_d(option = "cividis")
      #     } else if(input$PaletteChoicesQC == "magma") {
      #       scale_fill_viridis_d(option = "magma")
      #     } else if(input$PaletteChoicesQC == "plasma") {
      #       scale_fill_viridis_d(option = "plasma")
      #     }else if(input$PaletteChoicesQC == "inferno") {
      #       scale_fill_viridis_d(option = "inferno")
      #     }
      #   }) 
      # #PCA plot color palette choices for color
      # colorchoicesQC <-
      #   eventReactive(input$PaletteChoicesQC, {
      #     if(input$PaletteChoicesQC == "viridis") {
      #       scale_color_viridis_d(option = "viridis")
      #     } else if(input$PaletteChoicesQC == "cividis") {
      #       scale_color_viridis_d(option = "cividis")
      #     } else if(input$PaletteChoicesQC == "magma") {
      #       scale_color_viridis_d(option = "magma")
      #     } else if(input$PaletteChoicesQC == "plasma") {
      #       scale_color_viridis_d(option = "plasma")
      #     }else if(input$PaletteChoicesQC == "inferno") {
      #       scale_color_viridis_d(option = "inferno")
      #     }
      #   })
    }
  )
}

pcaApp <- function() {
  ui <- fluidPage(
    pcaUI("PCA1")
  )
  server<- function(input, output, session) {
    pcaServer("PCA1")
  }
  shinyApp(ui = ui, server = server)
}

pcaApp()
