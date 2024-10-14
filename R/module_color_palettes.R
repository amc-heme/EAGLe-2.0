library(RColorBrewer)

paletteUI <- function(id, choices) {
  ns <- NS(id)
  
  palettePicker(ns("palette"),
                "Choose a color palette",
                choices = list(
                  "RColorBrewer" = list(
                    "Set1" = brewer.pal(9, "Set1"),
                    "Set2" = brewer.pal(8, "Set2"),
                    "Dark2" = brewer.pal(8, "Dark2"),
                    "Paired" = brewer.pal(12, "Paired"),
                    "Pastel1" = brewer.pal(9, "Pastel1")
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
