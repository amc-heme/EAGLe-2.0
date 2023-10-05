#each config is a glabalData module that returns the appropriate data files when chosemn

#need to use global config for paths to dds files
# config.CD <- source("~/Documents/GitHub/EAGLe-2.0/config_CD.R")
# config.Ye16 <- source("~/Documents/GitHub/EAGLe-2.0/config_Ye16.R")

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

data_Server <- function(id) {
  moduleServer(id, function(input, output, session){
    #all of this needs to be changed. Read in dataset info from dataset yaml
# reactive({
#   if(datainput == "CancerDiscovery") {
#     return(config.CD)
#   } else if(datainput == "Ye16") {
#     return(config.Ye16)
#   } else {
#     return(NULL)
#   }
# })
    
    data_file_values <- list(
      "Cancer_Discovery" = ,
      "Ye_16" = 
    )
    config_choice <- eventReactive(input$datainput, {
      data_file_values[[input$datainput]]

    })
   #  print("what is config_choice")
   # print(config_choice)
    #return(config_choice)

  })
}
#config_files <- reactiveValues()

