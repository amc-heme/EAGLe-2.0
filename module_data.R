config.CD <- source("~/Documents/GitHub/EAGLe-2.0/config_CD.R")
config.Ye16 <- source("~/Documents/GitHub/EAGLe-2.0/config_Ye16.R")

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
                         "Ye et al, 2020" = "Ye20",
                         "Pollyea Ven/Aza, 2018" = "venaza",
                         "Lagadiou, 2013" = "stemcell",
                         "TCGA-LAML" = "tcga",
                         "BEAT-AML" = "beat",
                         "Target pediatric AML" = "target",
                         "Human Protein Atlas" = "hpa",
                         "Lee et al, 2018" = "lee"), selected = NULL
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
      "CancerDiscovery" = config.CD,
      "Ye16" = config.Ye16
    )
    
    config <- eventReactive(input$datainput, {
      data_file_values[[input$datainput]]
    })
    #reactive_config <- reactiveValues()
    
    ## load in data files
    # observe({
    #   config_choice <- config_x()
    #   #load objects for the chosen config file
    #   t2g <- read_rds(config_choice$t2g_file)
    #   metadata <- read_rds(config_choice$metadata_file)
    #   batch <- config_choice$batch
    #   #samples <- config$samples
    #   num_PCs <- nrow(metadata)
    #   dds <- read_rds(config_choice$dds_file)
    #   dds.res <- read_rds(config_choice$dds_res_file)
    #   vsd <- read_rds(config_choice$vsd_file)
    #   vsd.pca <- read_rds(config_choice$vsd.pca_file)
    #   vst <- read_rds(config_choice$vst_file)
    #   vst.goi <- read_rds(config_choice$vst.goi_file)
    #   qc <-load_multiqc(config_choice$qc, sections="raw")
    #   var_1 <- config_choice$var_1
    #   var_2 <- config_choice$var_2
    # })
    return(config)
  })
  }
