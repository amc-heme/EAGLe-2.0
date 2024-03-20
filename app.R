#load libraries
library(shinythemes)
library(thematic)
library(shiny)
library(kableExtra)
library(viridisLite)
library(magrittr)
library(ggplot2)
library(tidyr)
library(viridis)
library(tximport)
library(cowplot)
library(TidyMultiqc)
library(dplyr)
library(tidyverse)
library(data.table)
library(ggpubr)
library(RColorBrewer)
library(janitor)
library(reactlog)
library(DT)
library(ggrepel)
library(DESeq2)
library(ggiraph)
library(shinyWidgets)
library(shinyjs)
library(fgsea)
library(plotly)
library(BiocParallel)
library(ComplexHeatmap)
library(InteractiveComplexHeatmap)
library(circlize)
library(colourpicker)
library(colorRamp2)
library(ggsci)
library(scales)
library(esquisse)
library(ggprism)
library(shinycssloaders)
library(patchwork)
library(spatstat.utils)
library(yaml)
library(readr)
library(sva)
library(formulaic)
library(shinydashboard)


options(
  shiny.fullstacktrace = TRUE
)
## Data ####
# Config information for each dataset
dataset_config <- 
  read_yaml("./data.yaml")

raw_data_p <- getwd()
print(raw_data_p)
 
# Read RDS files from yaml
CD <- read_rds(paste0(raw_data_p, dataset_config[["Cancer_Discovery"]]$data_path))
Ye16 <- read_rds(paste0(raw_data_p, dataset_config[["Ye_16"]]$data_path))
Ye20 <- read_rds(paste0(raw_data_p, dataset_config[["Ye_20"]]$data_path))
Venaza <- read_rds(paste0(raw_data_p, dataset_config[["Venaza"]]$data_path))
Lagadinou <- read_rds(paste0(raw_data_p, dataset_config[["Lagadinou"]]$data_path))
Lee <- read_rds(paste0(raw_data_p, dataset_config[["Lee"]]$data_path))
BEAT <- read_rds(paste0(raw_data_p, dataset_config[["BEAT"]]$data_path))
TCGA <- read_rds(paste0(raw_data_p, dataset_config[["TCGA"]]$data_path))
HPA <- read_rds(paste0(raw_data_p, dataset_config[["HPA"]]$data_path))
dataset <- list(CD, Ye16, Ye20, Venaza, Lagadinou, Lee, BEAT, TCGA, HPA)

# Read QC RDS files from yaml
CD.qc <- read_rds(paste0(raw_data_p, dataset_config[["Cancer_Discovery"]]$qc_path))
Ye16.qc <- read_rds(paste0(raw_data_p, dataset_config[["Ye_16"]]$qc_path))
Ye20.qc <- read_rds(paste0(raw_data_p, dataset_config[["Ye_20"]]$qc_path))
Venaza.qc <- read_rds(paste0(raw_data_p, dataset_config[["Venaza"]]$qc_path))
Lagadinou.qc <- read_rds(paste0(raw_data_p, dataset_config[["Lagadinou"]]$qc_path))
Lee.qc <- read_rds(paste0(raw_data_p, dataset_config[["Lee"]]$qc_path))
BEAT.qc <- read_rds(paste0(raw_data_p, dataset_config[["BEAT"]]$qc_path))
TCGA.qc <- read_rds(paste0(raw_data_p, dataset_config[["TCGA"]]$qc_path))
HPA.qc <- read_rds(paste0(raw_data_p, dataset_config[["HPA"]]$qc_path))
dataset.qc <- list(CD.qc, Ye16.qc, Ye20.qc, Venaza.qc, Lagadinou.qc, Lee.qc, BEAT.qc, TCGA.qc, HPA.qc)

# Add dataset names to list generated
names(dataset) <- 
  names(dataset_config)

# Add dataset names to qc list generated
names(dataset.qc) <- 
  names(dataset_config)

# UI ####
ui <-
  navbarPage("EAGLe 2.0",
             
             #Dataset tab ####
             tabsetPanel(
               id = "page",
               type = "hidden",
               
               tabPanelBody("landing_page",
                            data_UI("data1")),
               
               tabPanelBody("content",
                            actionButton(
                              "change_data",
                              icon = icon("refresh"),
                              class = "btn-xs",
                              label = "Change Dataset",
                              style = "position: absolute; right: 40px"
                            ),      
                              #QC Menu ####
                              tabsetPanel(
                                tabPanel("QC",
                                      QC_UI("QC1")),
                              #DESeq Menu ####
                              tabPanel("Differential Expression",
                                      DE_UI("DEtab1")),
                              
                              #GSEA menu ####
                              tabPanel("GSEA",
                                      GSEA_UI("GSEA1")),
                              
                              #Gene expression analysis ####
                              tabPanel("Gene Expression",
                                      goi_UI("GOI1"))
                            ))
             ))

#Server ####
server <- 
  function(input, output, session) {
    # check to make sure server is initiating
    print("Initializing renderPlots")
    
    options(shiny.reactlog = TRUE)
  
  ##Data tab ####
    dataset_choice <- data_Server("data1")
    
    
    observeEvent(dataset_choice$close_tab(), {
      updateTabsetPanel(session, "page", "content")
    })
    
    observeEvent(input$change_data, {
      updateTabsetPanel(session, "page", "landing_page")
    })
    #reactive statement to return dataset species
    data_species <- reactive({
      dataset_config[[dataset_choice$user_dataset()]]$species
    })
    
    observe({
      print(data_species())
    })
    
  ## dds object
    dataset_dds <- dds.file_Server("dds1", dataset, dataset_choice)
    
  ## vst table
    vst <- vst_Server("vst1", dataset_dds, dataset_choice)
  ## qc object
    
    qc_table <- qc.file_Server("qct1", dataset.qc, dataset_choice)
  ## QC tab ####
  
    QC_Server("QC1",dataset_dds, dataset_choice, qc_table)
    
  ## GOI tab####
    goi_Server("GOI1", dataset_choice, dataset_dds, vst)
 
  # ##DESEq #####
    DE_res <- DE_Server("DEtab1", data_species, dataset_dds, dataset_choice) 

  # ##GSEA output ####
    GSEA_Server("GSEA1", dataset_choice, DE_res)
  #  
  # ##GOI pathway output ####
  #  pathway_Server("pathway1", dds, t2g)
  #   
  } #end server

# Run the application 
shinyApp(ui = ui, server = server)
