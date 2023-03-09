library(shiny)

DEUI <- function(id) {
  tagList(
    
  )
}

DEServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    
  })
}

DEApp <- function() {
  ui <- fluidPage(
    DEUI("hist1")
  )
  server <- function(input, output, session) {
    DEServer("hist1")
  }
  shinyApp(ui, server)  
}
DEApp()
