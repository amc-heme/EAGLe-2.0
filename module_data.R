
data_UI <- function(id) {
  ns <- NS(id)
  fluidPage(
    theme =
      shinytheme("flatly"),
    titlePanel("Choose a Dataset"
               ),
    sidebarLayout(
      sidebarPanel(
          selectInput(
            ns("datainput"), label = NULL,
             choices = c("Pei et al, 2020" = "Cancer_Discovery",
                         "Ye et al, 2016" = "Ye_16",
                         "Ye et al, 2020" = "Ye_20",
                         "Pollyea Ven/Aza, 2018" = "Venaza",
                         "Lagadiou, 2013" = "Lagadinou",
                         "TCGA-LAML" = "TCGA",
                         "BEAT-AML" = "BEAT",
                         "Lee et al, 2018" = "Lee",
                         "Human Protein Atlas" = "HPA"), selected = "Cancer_Discovery"
          )
      ),
      mainPanel(
        
      )
    )
  )
}

data_Server <- function(id) {
  moduleServer(id, function(input, output, session){
    
user_choice <- reactive({input$datainput})
observe({print(user_choice())})
return(user_choice)
  })
}



