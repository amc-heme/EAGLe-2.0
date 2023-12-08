# library(shiny)
# library(shinyjs)
# sliderUI <- function(id, min, max, value, label) {
#   ns <- NS(id)
#   tagList(
#     materialSwitch(ns("hidedims"), label = "Custom Plot Dimensions", value = FALSE, right = TRUE),
#     shinyjs::hidden(
#       sliderInput(ns("slider"), label = label,
#                   min = min, max = max, value = value
#       )
#     )
#   )
# }
# 
# sliderServer <- function(id) {
#   moduleServer(id, function(input, output, session) {
#     sliderscale <- reactive({
#       input$slider
#     })
#     observe({
#       shinyjs::toggle("slider", condition = input$hidedims)
#     })
#     return(sliderscale)
#   })
# }
# 
# sliderApp <- function() {
#   ui <- fluidPage(sliderUI("slider"))
#   server <- function(input, output, session) {
#     sliderServer("slider")
#   }
#   shinyApp(ui, server)
# }
# sliderApp()
