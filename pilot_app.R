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
vst.goi <- readRDS("/Users/stephanie/Documents/GitHub/EAGLe-2.0/data/vst.goi.rds")
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
      "Gene Centric Analysis",
      fluidPage(
        theme=
          shinytheme("flatly"),
        titlePanel(
            "Gene Centric Analysis"
            ),
          sidebarLayout(
            sidebarPanel(
              selectizeInput(
                "VSTCDgenechoice",
                label=
                  "Choose a gene for analysis",
                choices =
                  NULL,
                selected = NULL,
                options = list(maxItems = 5)
                ),
              radioButtons("XaxisVar_CDgene", h3("X axis variable"),
                           choices = list("Value" = 1, "Class" = 2,
                                          "Gene" = 3),selected = 3),
              radioButtons("YaxisVar_CDgene", h3("Y axis variable"),
                           choices = list("Value" = 4, "Class" = 5)
              #                             "Gene" = 6),selected = 4),
              # radioButtons("FillVar_CDgene", h3("Fill variable"),
              #              choices = list("Value" = 7, "Class" = 8,
              #                             "Gene" = 9), selected = 8)
              ),
            mainPanel(
              tabsetPanel(
                type =
                  "tabs",
                tabPanel(
                  "Gene Centric Plot",
                  plotOutput(
                    "VSTCDplot"
                    )
                  ),
                tabPanel(
                  "VST Summary",
                  textOutput(
                    "VSTCDsummary"
                    )
                  )
                )
              )
            )
        )
      ), #end Genecentric tabPanel
    tabPanel(
      "DESeq Analysis"
      ), #end DE tabPanel
    tabPanel(
      "GSEA"
      )#end GSEA tabPanel
    ),#end gene centric tabPanel #end CD navbarmenu
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
    QCdata <- switch(
      input$QCvar,
      "raw.salmon.percent_mapped" = qcdt$raw.salmon.percent_mapped,
      "raw.salmon.num_mapped" = qcdt$raw.salmon.num_mapped,
      "raw.star.uniquely_mapped_percent" = qcdt$raw.star.uniquely_mapped_percent,
      "raw.star.uniquely_mapped" = qcdt$raw.star.uniquely_mapped,
      "metadata.sample_id" = qcdt$metadata.sample_id
    )
    
    
    QCplot <- ggplot(qcdt, aes(metadata.sample_id, QCdata)) +
      geom_point(size = 3) +
      theme_cowplot (12) +
      theme(axis.text.x =
              element_text(angle = 60, hjust = 1))
    print(QCplot)
    QCplot
  }) #end renderPlot
#gene centric reactive plot output for selectizeInput 
  updateSelectizeInput(session,"VSTCDgenechoice", choices = vst.goi$ext_gene, server = TRUE)
#x axis output for gene centric plot in CD
  xvar_CDgene <- reactiveValues(data = vst.goi)
  observeEvent(input$XaxisVar_CDgene, {
    if (input$XaxisVar_CDgene == "1") {
      xvar_CDgene$data <- vst.goi$value
    }
    if (input$XaxisVar_CDgene == "2") {
      xvar_CDgene$data <- vst.goi$class
    }
    if (input$XaxisVar_CDgene == "3") {
      xvar_CDgene$data <- input$VSTCDgenechoice
    }
  })
#y axis output for gene centric plot in CD
  yvar_CDgene <- reactiveValues(data = vst.goi)
  observeEvent(input$YaxisVar_CDgene, {
    if (input$YaxisVar_CDgene == "4") {
      yvar_CDgene$data <- vst.goi$value
    }
    if (input$YaxisVar_CDgene == "5") {
      yvar_CDgene$data <- vst.goi$class
    }
    if (input$YaxisVar_CDgene == "6") {
      yvar_CDgene$data <- input$VSTCDgenechoice
    }
  })

  output$VSTCDplot <-
    renderPlot({
      #build a color palette
      colors <-
        colorRampPalette(c("red", "blue", "purple", "green", "yellow"))(12)
      
      ggplot(vst.goi,
             aes(
               x = xvar_CDgene$data,
               y = yvar_CDgene$data,
               fill = class
             )) +
        geom_boxplot(outlier.shape = NA) +
        scale_fill_manual(values = colors) +
        scale_color_manual(values = colors) +
        geom_point(alpha = 0.5,
                   position = position_jitterdodge(jitter.width = 0.2),
                   aes(color = class)) +
        theme_light() +
        ylab("") +
        xlab("") +
        ggtitle("Gene Expression:Sensitive vs Resistant")
    }) #end render plot
  
  output$BlastPlot <- renderPlot({
    BlastData <- switch (
      input$BlastVar,
      "SampleID" = meta_lut_ven$SampleID,
      "PercentBlastsInBM" = meta_lut_ven$PercentBlastsInBM,
      "PercentBlastsInPB" = meta_lut_ven$PercentBlastsInPB
    )
    
    ggplot(meta_lut_ven, aes(AUC, BlastData)) +
      geom_point() +
      theme_light() +
      geom_hline(yintercept = 60,
                 linetype = "dashed",
                 color = "red") +
      geom_vline(xintercept = 255,
                 linetype = "dashed",
                 color = "blue") +
      geom_vline(xintercept = 60,
                 linetype = "dashed",
                 color = "blue") +
      geom_vline(xintercept = 221,
                 linetype = "dashed",
                 color = "green") +
      geom_vline(xintercept = 94,
                 linetype = "dashed",
                 color = "green") +
      ggtitle(label = "PercentBlastsInBM and PercentBlastsInPB vs Ven AUC; 60% blast cutoff = 77")
  }) #end renderPlot
  } #end server
# Run the application 
 shinyApp(ui = ui, server = server)


