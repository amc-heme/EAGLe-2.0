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
#load data
meta_lut_ven <- readRDS("/Users/stephanie/Documents/App-1/ShinyLessons/data/meta_lut_ven.Rds")
qcdt<-load_multiqc("/Users/stephanie/Library/Mobile Documents/com~apple~CloudDocs/Desktop/InformaticsProjectfiles/multiqc_data.json", sections="raw") 
# exp.jordan.m0m5 <- read.csv("/Users/stephanie/Documents/App-1/ShinyLessons/data/vstlimma.csv", header = FALSE)
# label.jordan.m0m5 <- read.csv("/Users/stephanie/Documents/App-1/ShinyLessons/data/Jordan_M0M5_ROSlow_label.csv", header = TRUE) %>% as_tibble()
# eval.jordan.m0m5 <- as.logical(gene %in% exp.jordan.m0m5$V2)
# gene <- params$gene
# # isolate data for CD gene centric plot
# gene.jordan.m0m5 <- exp.jordan.m0m5 %>%
#   dplyr::filter(V2== "ext_gene"| V2== gene )
# # transpose data for CD gene centric plot
# gene.jordan.m0m5 <- t(as.data.frame(exp.jordan.m0m5))
# #rownames(gene.jordan.m0m5) <- NULL
# #gene.jordan.m0m5 <- as_data_frame(gene.jordan.m0m5) %>%
#   dplyr::rename(Specimen = 1, gene = (2:4)) %>%
#   dplyr::filter(Specimen != "ensembl_gene_id")%>%
#   dplyr::filter(Specimen != "ext_gene")
# gene.jordan.m0m5$gene1 <- as.numeric(gene.jordan.m0m5$gene1)
# 
# AML<- c("AML011114A", "AML052214", "AML011114B", "AML011514", "AML030514", "AML041514", "AML011612", "AML020604", "AML041206", "AML062717", "AML100510", "AML100814")
# gene.jordan.m0m5$AML<- AML
# gene.jordan.m0m5<-gene.jordan.m0m5[c("AML", "Specimen", "gene1", "gene2","gene3")] 
# 
# # Add labels to table for CD gene centric plot
# df.jordan.m0m5 <-inner_join(gene.jordan.m0m5, label.jordan.m0m5, by = "AML")
# df.jordan.m0m5 <-df.jordan.m0m5[c("FAB", "Specimen.y", "Specimen.x", "AML", "gene1", "gene2","gene3")]

ui <- navbarPage(
  "EAGLe",
  navbarMenu(
    "Cancer Discovery",
    tabPanel(
      "QC",
      fluidPage(
        theme =
          shinytheme("flatly"),
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
            )
          ),
          mainPanel(
            tabsetPanel(
              type = "tabs",
              tabPanel(
                "Plot",
                plotOutput("QCplot")
              ),
              tabPanel(
                "Table",
                tableOutput("table")
              ),
              tabPanel(
                "Summary",
                textOutput("summary")
              )
            )
          )
        )
      )
    ),
    tabPanel("Gene Centric Analysis"),
    # fluidPage(
    #   theme= shinytheme("flatly"),
    #   titlePanel("Gene Centric Analysis"),
    #   sidebarLayout(
    #     sidebarPanel(
    #       selectizeInput(
    #         "gene_choice",
    #         label= "Choose a gene",
    #         choices = exp.jordan.m0m5$V2,
    #         selected = NULL, options = list(maxItems = 5)),
    #     hr(),
    #     selectInput(
    #       "gene_meta",
    #       label= "Choose a metadata variable",
    #       choices = c("FAB", "Specimen", "SRR", "AML"),
    #       selected = "FAB"),
    #     ),
    #     mainPanel(
    #       tabsetPanel(
    #         type = "tabs",
    #         tabPanel("Plot", plotOutput("plotgraph3")),
    #         tabPanel("Table", tableOutput("table")),
    #         tabPanel("Summary",textOutput("summary")))
    #       )
    #     )
    #   )
    # ),
    tabPanel("DESeq Analysis"),
    tabPanel("GSEA")
  ),
  
  navbarMenu(
    "BEAT-AML",
    tabPanel(
      "Gene Centric Analysis",
      fluidPage(
        theme =
          shinytheme("flatly"),
        titlePanel(
          "BEAT-AML"
        ),
        sidebarLayout(
          sidebarPanel(
            selectInput(
              "Blastvar",
              label =
                "Choose a Variable to Display",
              choices =
                c("SampleID", "PercentBlastsInBM", "PercentBlastsInPB"),
              selected = "PercentBlastsInBM"
            )
          ),
          mainPanel(
            tabsetPanel(
              type =
                "tabs",
              tabPanel(
                "Plot",
                plotOutput("Blastplot")
              ),
              tabPanel(
                "Table",
                tableOutput("table")
              ),
              tabPanel(
                "Summary",
                textOutput(
                  "summary")
              )
            )
          )
        )
      )
    ),
    tabPanel("DESeq and GSEA Anlaysis")
  ) #end BEAT drop down
) #end navbarPage

server <- function(input, output, session) {
  output$QCplot <- 
    renderPlot({
      print("computing QCplot")
      QCdata <- switch(input$QCvar,
                       "raw.salmon.percent_mapped"= qc$raw.salmon.percent_mapped,
                       "raw.salmon.num_mapped" = qc$raw.salmon.num_mapped,
                       "raw.star.uniquely_mapped_percent" = qc$raw.star.uniquely_mapped_percent,
                       "raw.star.uniquely_mapped" = qc$raw.star.uniquely_mapped,
                       "metadata.sample_id"= qc$metadata.sample_id)
      
      ggplot(qcdt, aes(metadata.sample_id, QCdata)) +
        geom_point(size = 3) +
        theme_cowplot (12) +
        theme(axis.text.x =
                element_text(angle= 60, hjust = 1))
    })
  output$Blastplot <- 
    renderPlot({
      print("computing blastplot")
      Blastdata <- switch(input$Blastvar,
                          "SampleID" = meta_lut_ven$SampleID,
                          "PercentBlastsInBM" = meta_lut_ven$PercentBlastsInBM,
                          "PercentBlastsInPB" = meta_lut_ven$PercentBlastsInPB)
      print("data:")
      print(Blastdata)
      
      ggplot(meta_lut_ven, aes(AUC, Blastdata)) +
        geom_point() +
        theme_light() +
        geom_hline(yintercept=60, linetype="dashed", color = "red") +
        geom_vline(xintercept=255, linetype="dashed", color = "blue") +
        geom_vline(xintercept=60, linetype="dashed", color = "blue") +
        geom_vline(xintercept=221, linetype="dashed", color = "green") +
        geom_vline(xintercept=94, linetype="dashed", color = "green") +
        ggtitle(label = "PercentBlastsInBM and PercentBlastsInPB vs Ven AUC; 60% blast cutoff = 77")
    })
  #updateSelectizeInput(session,"gene_choice", choices = exp.jordan.m0m5$V2, server = TRUE) 
  
  # output$plotgraph3 <- renderPlot({
  # 
  #   
  # data3 <- switch(input$gene_choice,
  #                "Gene1" = df.jordan.m0m5$gene1,
  #                "Gene2" = df.jordan.m0m5$gene2,
  #                "Gene3" = df.jordan.m0m5$gene3)
  # meta_dat <- switch(input$gene_meta,
  #                     "FAB"= df.jordan.m0m5$FAB,
  #                     "Specimen"= df.jordan.m0m5$Specimen.y,
  #                     "SRR"= df.jordan.m0m5$Specimen.x,
  #                     "AML" = df.jordan.m0m5$AML)
  # 
  # ggplot(df.jordan.m0m5, aes(x= meta_dat, y= data3, fill = meta_dat)) +
  #   geom_boxplot(alpha = 1.5) +
  #   scale_fill_viridis(discrete= TRUE) +
  #   geom_point(x= meta_dat, y= data3) +
  #   theme_cowplot() +
  #   geom_hline(yintercept=10, linetype="dashed", color = "red") +
  #   ylab("Norm Counts") +
  #   xlab(input$gene_meta) +
  #   ggtitle("")
  # })
}
# Run the application 
shinyApp(ui,server)


