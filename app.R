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
library(yaml)
options(
  shiny.fullstacktrace = TRUE
)
## Data ####
# Config information for each dataset
dataset_config <- 
  read_yaml("./data.yaml")

raw_data_p <- getwd()

 
# Read RDS files from yaml
CD <- read_rds(paste0(raw_data_p, dataset_config[["Cancer_Discovery"]]$data_path))
Ye16 <- read_rds(paste0(raw_data_p, dataset_config[["Ye_16"]]$data_path))
Ye20 <- read_rds(paste0(raw_data_p, dataset_config[["Ye_20"]]$data_path))
Venaza <- read_rds(paste0(raw_data_p, dataset_config[["Venaza"]]$data_path))
Lagadinou <- read_rds(paste0(raw_data_p, dataset_config[["Lagadinou"]]$data_path))
Lee <- read_rds(paste0(raw_data_p, dataset_config[["Lee"]]$data_path))
BEAT <- read_rds(paste0(raw_data_p, dataset_config[["BEAT"]]$data_path))
dataset <- list(CD, Ye16, Ye20, Venaza, Lagadinou, Lee, BEAT)

# Add dataset names to list generated
names(dataset) <- 
  names(dataset_config)


# UI ####
ui <-
  navbarPage(
    "EAGLe: Cancer Discovery",
    #Dataset tab ####
    tabPanel(
      "Dataset",
      data_UI("data1")
    ),
    #QC Menu ####
    tabPanel(
              "QC",
              QC_UI("QC1")
    ),
    # 
    #Gene expression analysis ####
    tabPanel( 
              "Gene Expression",
              goi_UI("GOI1")
    ),
    
    #DESeq Menu ####
      tabPanel("Differential Expression",
             DE_UI("DEtab1")
             ),
    
    #GSEA menu ####
    tabPanel("GSEA",
             GSEA_UI("GSEA1")
    ),

    #GOI plots ####
    tabPanel("GSEA Pathway/Gene Visualization",
     pathway_UI("pathway1")
  )
  )


#Server ####
server <- 
  function(input, output, session) {
    # check to make sure server is initiating
    print("Initializing renderPlots")
    
    options(shiny.reactlog = TRUE)
  
  
  ##Data tab ####
    dataset_choice <- data_Server("data1")
    
    #reactive statement to return dataset species
    data_species <- reactive({
      dataset_config[[dataset_choice()]]$species
    })
    observe({
      print(data_species())
    })
    
  ## dds object
    dataset_dds <- dds.file_Server("dds1", dataset, dataset_choice)

  ## QC tab ####
  
    #QC_Server("QC1", dataset_choice, dataset_dds)
    
  ## GOI tab####  
    # goi_Server("GOI1", dataset_dds)
 
  # ##DESEq #####
     DE_Server("DEtab1",data_species, dataset_dds) 
 
  # ##GSEA output ####
  #   GSEA_Server("GSEA1", dds, t2g) 
  #  
  #   
  # ##GOI pathway output ####
  #  pathway_Server("pathway1", dds, t2g)
  #   
  } #end server

# Run the application 
shinyApp(ui = ui, server = server)
