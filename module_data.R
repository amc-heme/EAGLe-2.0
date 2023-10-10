
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
                         "Lee et al, 2018" = "Lee"), selected = NULL
          )
      ),
      mainPanel(
        
      )
    )
  )
}

data_Server <- function(id, dataset) {
  moduleServer(id, function(input, output, session){
    #all of this needs to be changed. Read in dataset info from dataset yaml
data_list <- dataset

  
  dds_object <- eventReactive(input$datainput, {
    dataset_dds <- data_list[[input$datainput]]
    dataset_dds
  })
    return(dds_object)
   
  })
}



