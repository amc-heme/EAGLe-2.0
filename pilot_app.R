#load libraries
library(shinythemes)
library(thematic)
library(shiny)
library(magrittr)
library(ggplot2)
library(tidyr)
library(viridis)
library(cowplot)
library(TidyMultiqc)
library(dplyr)
library(tidyverse)
library(data.table)
library(ggpubr)
library(RColorBrewer)
library(janitor)
library(reactlog)
#options(shiny.reactlog = TRUE)
#reactlogShow(time = TRUE)

#load data
meta_lut_ven <- readRDS("/Users/stephanie/Documents/GitHub/EAGLe-2.0/data/meta_lut_ven.Rds")
qcdt<-load_multiqc("/Users/stephanie/Documents/GitHub/EAGLe-2.0/data/multiqc_data.json", sections="raw") 
# exp.jordan.m0m5 <- read.table("/Users/stephanie/Documents/App-1/ShinyLessons/data/vstlimma.csv",sep = ",",header = TRUE)
# label.jordan.m0m5 <- read.csv("/Users/stephanie/Documents/App-1/ShinyLessons/data/Jordan_M0M5_ROSlow_label.csv", header = TRUE) %>% as_tibble()
# eval.jordan.m0m5 <- as.logical(gene %in% exp.jordan.m0m5$ext_gene)

ui <-
  navbarPage(
  "EAGLe",
  navbarMenu(
    "Cancer Discovery",
    tabPanel(
      "QC",
      fluidPage(
        theme =
          shinytheme(
            "flatly"
            ),
        titlePanel(
          "MultiQC Analysis"
        ), 
        sidebarLayout(
          sidebarPanel(
            selectInput(
              "QCvar",
              label=
                "Choose a QC Variable to Display",
              choices =
                c("raw.salmon.percent_mapped", "raw.salmon.num_mapped", "raw.star.uniquely_mapped_percent", "raw.star.uniquely_mapped"),
              selected =
                "raw.salmon.percent_mapped"
            ) #end selectInput
          ), #end sidebarPanel
          mainPanel(
            tabsetPanel(
              type =
                "tabs",
              tabPanel(
                "QC Plot",
                plotOutput(
                  "QCplot"
                  )
              ), #end tabPanel
              tabPanel(
                "Table1",
                tableOutput(
                  "table1"
                  )
              ), #end tabPanel
              tabPanel(
                "Summary1",
                textOutput(
                  "summary1"
                  )
                ) #end tabPanel
              ) #end tabsetPanel
            ) #end mainPanel
          ) #end sidebarLayout
        ) #end fluidPage
      ), #end QC tabPanel
    tabPanel(
      "Gene Centric Analysis"
      ), #end Genecentric tabPanel
    tabPanel(
      "DESeq Analysis"
      ), #end DE tabPanel
    tabPanel(
      "GSEA"
      ) #end GSEA tabPanel
  ), #end CD navbarmenu
  navbarMenu("BEAT-AML",
    tabPanel(
      "Gene Centric Plots",
      fluidPage(
        theme =
          shinytheme(
            "flatly"
            ),
        titlePanel(
          "BEAT-AML"
        ), #end title
        sidebarLayout(
          sidebarPanel(
            selectInput(
              "BlastVar",
              label =
                "Choose a Variable to Display",
              choices =
                c("SampleID", "PercentBlastsInBM", "PercentBlastsInPB"),
              selected =
                "PercentBlastsInBM"
            ) #end selectInput
          ), #end sidebarPanel
          mainPanel(
            tabsetPanel(
              type =
                "tabs",
              tabPanel(
                "BlastPercentagePlot",
                plotOutput(
                  "BlastPlot"
                  )
              ), #end tabPanel
              tabPanel(
                "Table2",
                tableOutput(
                  "table2"
                  ) 
              ),#end tabPanel
              tabPanel(
                "Summary2",
                textOutput(
                  "summary2"
                  )
                ) #end tabPanel
              ) #end tabsetPanel
            ) #end mainPanel
          ) #end sidebarLayout
        ) #end fluidPage
      ), #end Gene centric tabPanel
    tabPanel(
      "DESeq and GSEA Anlaysis"
      ) #end tabPanel
    ) #end BEAT navbarMenu
  ) #end navbarPage

server <- 
  function(input, output, session) {
  
  print("Initializing renderPlots")
  options(shiny.reactlog = TRUE)
  
  output$QCplot <- renderPlot({
    QCdata <- switch(input$QCvar,
                     "raw.salmon.percent_mapped"= qcdt$raw.salmon.percent_mapped,
                     "raw.salmon.num_mapped" = qcdt$raw.salmon.num_mapped,
                     "raw.star.uniquely_mapped_percent" = qcdt$raw.star.uniquely_mapped_percent,
                     "raw.star.uniquely_mapped" = qcdt$raw.star.uniquely_mapped,
                     "metadata.sample_id"= qcdt$metadata.sample_id)
    
    
   QCplot <- ggplot(qcdt, aes(metadata.sample_id, QCdata)) +
        geom_point(size = 3) +
        theme_cowplot (12) +
        theme(axis.text.x =
                element_text(angle= 60, hjust = 1))
      print(QCplot)
      QCplot
    }) #end renderPlot
  
  
   output$BlastPlot <-renderPlot({
     BlastData <- switch (input$BlastVar,
                 "SampleID" = meta_lut_ven$SampleID,
                 "PercentBlastsInBM" = meta_lut_ven$PercentBlastsInBM,
                 "PercentBlastsInPB" = meta_lut_ven$PercentBlastsInPB)
     
      ggplot(meta_lut_ven, aes(AUC, BlastData)) +
        geom_point() +
        theme_light() +
        geom_hline(yintercept=60, linetype="dashed", color = "red") +
        geom_vline(xintercept=255, linetype="dashed", color = "blue") +
        geom_vline(xintercept=60, linetype="dashed", color = "blue") +
        geom_vline(xintercept=221, linetype="dashed", color = "green") +
        geom_vline(xintercept=94, linetype="dashed", color = "green") +
        ggtitle(label = "PercentBlastsInBM and PercentBlastsInPB vs Ven AUC; 60% blast cutoff = 77")
    }) #end renderPlot
   print("complete")
} #end server


# Run the application 
shinyApp(ui = ui, server = server)


