config.CD <- "config_CD.R"


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
             choices = c("Pei et al, 2020" = "CancerDiscovery",
                         "Ye et al, 2016" = "Ye16",
                         "Ye et al 2020" = "Ye20",
                         "Pollyea Ven/Aza 2018" = "venaza",
                         "Lagadiou 2013" = "stemcell",
                         "TCGA-LAML" = "tcga",
                         "BEAT-AML" = "beat",
                         "Target pediatric AML" = "target",
                         "Human Protein Atlas" = "hpa",
                         "Lee et al 2018" = "lee")
          )
      ),
      mainPanel(
        
      )
    )
  )
}

data_Server <- function(id) {
  moduleServer(id, function(input, output, session){
    data_file_values <- list(
      "CancerDiscovery" = config.CD
    )
    
    config <- eventReactive(input$datainput, {
      data_file_values[[input$dataunput]]
    })
    return(config)
  })
  }
