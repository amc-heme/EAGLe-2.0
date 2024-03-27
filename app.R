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
library(waiter)

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
                                id = "tabs",
                              
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
             ),
             mainPanel(
                   div(
                     tags$li(
                       class = "dataset_info_button",
                       style = 'display:inline-block;',
                       style = 'float:top;',
                       dropMenu(
                         dropdownButton(
                           circle = TRUE, status = "info", icon = icon("info"),  size = 'sm', width = "50px",
                           tooltip = tooltipOptions(title = "Click for more information on datasets")),
                         h6("Summary of datasets:"),
                         h6("Jordan M0-M5: "),
                         h6("Pei et al. have previously shown that ROS-low enriched LSCs from primitive and monocytic AML 
             differ significantly in their transcriptome and metabolism. 
             ROS-low LSCs from monocytic AML in particular, are less dependent on BCL2 therefore are more likely to be resistant to Venetoclax-based therapy.
             This section of analysis shows gene expression in primitive vs. monocytic ROS-low LSCs."),
                         h6("Ye et al. 2016: "),
                         h6("LSC gene expression from blood, bone marrow, spleen, gonadal adipose tissue, and normal bone marrow in mice."),
                         h6("Ye et al. 2020: "),
                         h6("Transcriptomes of LSCs from bone marrow and liver in mice were compared to better understand 
             the biology of liver LSCs. A bcCML model (BCR-ABL + Nup98-Hoxa9) was used. 
             Liver and bone marrow samples were combined from 3 mice in 3 different cohorts before sequencing."),
                         h6("Pollyea et al., Nature, 2018: "),
                         h6("Time zero pheresis was taken from 3 patients, then collected again at 6hr and 24hr post ven/aza treatment."),
                         h6( "Lagadinou et al., Cell Stem Cell, 2013: "),
                         h6("BCL-2 inhibition targets oxidative phosphorylation and selectively eradicates quiescent human leukemia stem cells. 
             ROS high and ROS low LSCs were either treated with 5ul of PTL or not treated before sequencing."),
                         h6("Lee et al., Nature, 2018: "),
                         h6("Genome wide expression from 12 AML patient samples with either prior complete remission, or no prior complete remission."),
                         h6("TCGA: "),
                         h6("The Cancer Genome Atlas (TCGA) project published gene expression profiles of 151 primary AML patients 
             along with their mutational profiles and clinical characteristics. This section of analysis shows gene expression 
             in the TCGA-AML dataset parsed by various mutational/clinical variables including the French-American-British 
             subtypes, karyotype, RAS mutation status, and NPM1 mutation status."),
                         h6("BEAT-AML: "),
                         h6("The BEAT-AML project published gene expression profiles of ~400 primary AML patients 
          along with their mutational profiles and clinical characteristics. 
          This section of analysis shows gene expression in the BEAT-AML dataseit parsed by various mutational/clinical variables 
             including the French-American-British subtypes, Venetoclax response, and disease stage (de novo vs. relapse).")
                         
                       )))
                   
             ))

#Server ####
server <- 
  function(input, output, session) {
    # check to make sure server is initiating
    print("Initializing renderPlots")
    
    options(shiny.reactlog = TRUE)
  
    #reactive container for reset button(change dataset action button)
    reset_trigger <- reactive({
      input$change_data
      })
    
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
  
    QC_Server("QC1", dataset_dds, dataset_choice, qc_table, reset_trigger)
    
  ## GOI tab####
    goi_Server("GOI1", dataset_choice, dataset_dds, vst)
 
  # ##DESEq #####
    DE_res <- DE_Server("DEtab1", data_species, dataset_dds, dataset_choice, reset_trigger) 

  # ##GSEA output ####
    GSEA_Server("GSEA1", dataset_choice, DE_res, reset_trigger)
    
    observeEvent(input$change_data, {
      updateTabsetPanel(session, "tabs", selected = "QC")
    })

  } #end server

# Run the application 
shinyApp(ui = ui, server = server)
