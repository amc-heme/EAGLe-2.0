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
dds.res <- readRDS("/Users/stephanie/Documents/GitHub/EAGLe-2.0/data/DEtable.rds")
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
              radioButtons("XaxisVar_CDgene", h4("X axis variable"),
                           choices = list("Value" = "xvalue", "Class" = "xclass",
                                          "Gene" = "xgene"),selected = "xgene"),
              radioButtons("YaxisVar_CDgene", h4("Y axis variable"),
                           choices = list("Value" = "yvalue", "Class" = "yclass",
                                          "Gene" = "ygene"),selected = "yvalue"),
              radioButtons("FillVar_CDgene", h4("Fill variable"),
                           choices = list("Class" = "fillclass", "Gene" = "fillgene"), selected = "fillclass")
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
    tabPanel("DESeq2 Analysis",
             fluidPage(
               theme =
                 shinytheme("flatly"),
               titlePanel("DESeq2 Table and Plot"),
               #end title
               sidebarLayout(
                 sidebarPanel( 
                   selectizeInput( #pick genes to filter DE table by
                     "DECDgenechoice",
                     label=
                       "Choose a gene or group of genes from DE table",
                     choices =
                       NULL,
                     selected = NULL,
                     options = list(maxItems = 5)
                   ),
                   
                   hr(),
                   
                   sliderInput( # OR pick a padj value range to search table by
                     "CDpadj_slider",
                     label = h4(
                       "Select padj value range"
                       ),
                     min = 0,
                     max = 8,
                     value = c(0, 1)
                     ),
                   
                   hr(),
                   
                   sliderInput( #filter by expression
                     "CDlog2foldchangeslider",
                     label = h4(
                       "Select log2 fold change range"
                     ),
                     min = -4,
                     max = 5,
                     value = c(0, 5)
                   )
                   ),
                 mainPanel(
                   tabsetPanel(
                     type =
                       "tabs",
                     tabPanel(
                       "DE Table",
                       tableOutput(
                         "DETable"
                       )
                     ),
                     tabPanel(
                       "DE Plot",
                       plotOutput(
                         "DE Volcano Plot"
                         )
                       )
                     )
                   ) #end mainPanel
                 ) #end sidebarlayout
               ) #end fluidPage
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
  
  
##QC-Cancer Discovery output
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
  
##Gene Centric-Cancer Discovery output
 updateSelectizeInput(session,"VSTCDgenechoice", choices = vst.goi$ext_gene, server = TRUE)
 
 datavst<-
   reactive({
     vst.goi %>% 
       dplyr::filter(ext_gene %in% input$VSTCDgenechoice)
   })
  #x axis output
 xvar_CDgene <-
   eventReactive(input$XaxisVar_CDgene, {
     if (input$XaxisVar_CDgene == "xvalue") {
       "value"
     } else if (input$XaxisVar_CDgene == "xclass") {
       "class"
     } else if (input$XaxisVar_CDgene == "xgene") {
       "ext_gene"
     }
   })
  #y axis output f
 yvar_CDgene <-
   eventReactive(input$YaxisVar_CDgene, {
     if (input$YaxisVar_CDgene == "yvalue") {
       "value"
     } else if (input$YaxisVar_CDgene == "yclass") {
       "class"
     } else if (input$YaxisVar_CDgene == "ygene") {
       "ext_gene"
     }
   })
  #fill output
 fillvar_CDgene <-
   eventReactive(input$FillVar_CDgene, {
     if (input$FillVar_CDgene == "fillclass") {
       "class"
     } else if (input$FillVar_CDgene == "fillgene") {
       "ext_gene"
     }
   })
 #plot output
  output$VSTCDplot <-
    renderPlot({
 #build a color palette
      colors <-
        colorRampPalette(c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2"))(5)
      ggplot(datavst(),
             aes(
               x = .data[[xvar_CDgene()]],
               y =  .data[[yvar_CDgene()]],
               fill = .data[[fillvar_CDgene()]]
             )) +
        geom_boxplot(outlier.shape = NA) +
        scale_fill_manual(values = colors) +
        scale_color_manual(values = colors) +
        geom_point(alpha = 0.5,
                   position = position_jitterdodge(jitter.width = 0.2),
                   aes(color = ext_gene)) + #this needs to be reactive too
        theme_light() +
        ylab("") +
        xlab("") +
        ggtitle("Gene Expression:Sensitive vs Resistant")
    }) #end render plot

##Blast Percent plots- BEAT-AML output
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

  
## DESeq2- Cancer Discovery outputs
  updateSelectizeInput(session,"DECDgenechoice", choices = dds.res$Gene, server = TRUE)
  
  dataDECD<-
    reactive({
      dds.res %>% 
        dplyr::filter(Gene %in% input$DECDgenechoice)
    })
  DEpadj <-
    eventReactive({
      dds.res %>% 
        dplyr::filter(padj %in% input$CDpadj_slider)
    })
  DElog2 <-
    eventReactive({
      dds.res %>% 
        dplyr::filter(log2FoldChange %in% input$CDlog2foldchangeslider)
    })
  output$DETable <- renderTable({
    
  })
  } #end server
# Run the application 
 shinyApp(ui = ui, server = server)


