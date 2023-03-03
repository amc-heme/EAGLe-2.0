##practice space for creating modules for each analysis within the EAGLe app

library(shiny)
radio_button_UI<- function(id) {
  ns <- shiny::NS(id)
  tagList(
    radioButtons(inputId = ns("XaxisVar_CDgene"), h4("X axis variable"),
                 choices = list("Value" = "xvalue",
                                "Gene" = "xgene", "Class" = "xclass"),selected = "xgene")
  )
}

radio_button_Server <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      observeEvent(input$XaxisVar_CDgene, {
        if(input$XaxisVar_CDgene == "xvalue") {
          mychoices <- c("Gene" = "ygene")
        } else if(input$XaxisVar_CDgene=="xgene") {
          mychoices <- c("Value" = "yvalue")
        }else if(input$XaxisVar_CDgene == "xclass") {
          mychoices <- c("Value" = "yvalue")
        }
      })
    }
  )
}

radio_buttonApp<- function() {
  ui <- fluidPage(
    radio_button_UI("radio1")
  )
  server<- function(input, output, session) {
    radio_button_Server("radio1")
  }
  shinyApp(ui = ui, server = server)
}
radio_buttonApp()
# 
# library(shiny)
# histogramUI <- function(id) {
#   tagList(
#     selectInput(NS(id, "var"), "Variable", choices = names(mtcars)),
#     numericInput(NS(id, "bins"), "bins", value = 10, min = 1),
#     plotOutput(NS(id, "hist"))
#   )
# }
# histogramServer <- function(id) {
#   moduleServer(id, function(input, output, session) {
#     data <- reactive(mtcars[[input$var]])
#     output$hist <- renderPlot({
#       hist(data(), breaks = input$bins, main = input$var)
#     }, res = 96)
#   })
# }
# histogramApp <- function() {
#   ui <- fluidPage(
#     histogramUI("hist1")
#   )
#   server <- function(input, output, session) {
#     histogramServer("hist1")
#   }
#   shinyApp(ui, server)  
# }
# histogramApp()
