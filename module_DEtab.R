#DE tab module

DEtab_UI <- function(id)  {
  ns <- NS(id)
  fluidPage(
    theme =
      shinytheme("flatly"),
    titlePanel(
      "Differential Expression Tables and Plots"
    ),#end title
  tagList(
        DE_UI(ns("DEt1")),
        DE_plots_UI(ns("volplot"))
    )
  )
}

DEtab_Server <- function(id) {
  moduleServer(id, function(input, output, session) {
  
    DE_Server("DEt1")
    DE_plots_Server("volplot")
      
  })
}

DE_tab <- function() {
  ui <- fluidPage(
    DEtab_UI("DEta1")
  )
  server <- function(input, output, session) {
    DEtab_Server("DEta1")
  }
  shinyApp(ui, server)
}
DE_tab()
