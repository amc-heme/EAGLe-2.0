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
library(DT)
library(ggrepel)
library(DESeq2)
library(shinyWidgets)
library(shinyjs)

#options(shiny.reactlog = TRUE)
#reactlogShow(time = TRUE)

#Data ####
meta_lut_ven <- readRDS("/Users/stephanie/Documents/GitHub/EAGLe-2.0/data/meta_lut_ven.Rds")
qcdt<-load_multiqc("/Users/stephanie/Documents/GitHub/EAGLe-2.0/data/multiqc_data.json", sections="raw") 
vst.goi <- readRDS("/Users/stephanie/Documents/GitHub/EAGLe-2.0/data/vst.goi.rds")
dds.res <- readRDS("/Users/stephanie/Documents/GitHub/EAGLe-2.0/data/DEtable.rds")
metadata <- read.table(file = "/Users/stephanie/Documents/GitHub/EAGLe-2.0/data/SampleSheetJordanLab.txt")
nonvsd.pca <- readRDS("/Users/stephanie/Documents/GitHub/EAGLe-2.0/data/nonvsd.pca.rds")
vsd.pca <- readRDS("/Users/stephanie/Documents/GitHub/EAGLe-2.0/data/vsd.pca.rds")
bcvsd.pca <- readRDS("/Users/stephanie/Documents/GitHub/EAGLe-2.0/data/bcvsd.pca.rds")
nonvsd.variance <- readRDS("/Users/stephanie/Documents/GitHub/EAGLe-2.0/data/nonvsd.variance.rds")
vsd.variance <- readRDS("/Users/stephanie/Documents/GitHub/EAGLe-2.0/data/vsd.variance.rds")
bcvsd.variance <- readRDS("/Users/stephanie/Documents/GitHub/EAGLe-2.0/data/bcvsd.variance.rds")
nonvsd2.pca <- readRDS("/Users/stephanie/Documents/GitHub/EAGLe-2.0/data/nonvsd2.pca.rds")
vsd2.pca <- readRDS("/Users/stephanie/Documents/GitHub/EAGLe-2.0/data/vsd2.pca.rds")
bcvsd2.pca <- readRDS("/Users/stephanie/Documents/GitHub/EAGLe-2.0/data/bcvsd2.pca.rds")

#load pathways
pathways.hallmark <- gmtPathways("/Users/stephanie/Documents/GitHub/EAGLe-2.0/data/gmt_pathway_files copy/h.all.v7.4.symbols.gmt")
pathways.GOall <- gmtPathways("/Users/stephanie/Library/Mobile Documents/com~apple~CloudDocs/Desktop/InformaticsProjectfiles/gmt_pathway_files/c5.go.v2022.1.Hs.symbols.gmt")
pathways.GOmolec <- gmtPathways("/Users/stephanie/Documents/GitHub/EAGLe-2.0/data/gmt_pathway_files copy/c5.go.mf.v7.4.symbols.gmt")
pathways.GOcellcomp <- gmtPathways("/Users/stephanie/Documents/GitHub/EAGLe-2.0/data/gmt_pathway_files copy/c5.go.cc.v2022.1.Hs.symbols.gmt") 
pathways.GObio <- gmtPathways("/Users/stephanie/Documents/GitHub/EAGLe-2.0/data/gmt_pathway_files copy/c5.go.bp.v7.4.symbols.gmt")
pathways.TFtargets <-gmtPathways("/Users/stephanie/Documents/GitHub/EAGLe-2.0/data/gmt_pathway_files copy/c3.tft.v2022.1.Hs.symbols.gmt")
pathways.allReg <- gmtPathways("/Users/stephanie/Documents/GitHub/EAGLe-2.0/data/gmt_pathway_files copy/c3.all.v2022.1.Hs.symbols.gmt")
pathways.Wiki <- gmtPathways("/Users/stephanie/Documents/GitHub/EAGLe-2.0/data/gmt_pathway_files copy/c2.cp.wikipathways.v2022.1.Hs.symbols.gmt")
pathways.Reactome <-gmtPathways("/Users/stephanie/Documents/GitHub/EAGLe-2.0/data/gmt_pathway_files copy/c2.cp.reactome.v2022.1.Hs.symbols.gmt")
pathways.KEGG <- gmtPathways("/Users/stephanie/Documents/GitHub/EAGLe-2.0/data/gmt_pathway_files copy/c2.cp.kegg.v2022.1.Hs.symbols.gmt")
pathways.Positional <-gmtPathways("/Users/stephanie/Documents/GitHub/EAGLe-2.0/data/gmt_pathway_files copy/c1.all.v2022.1.Hs.symbols.gmt")
pathways.Biocarta <-gmtPathways("/Users/stephanie/Documents/GitHub/EAGLe-2.0/data/gmt_pathway_files copy/c2.cp.biocarta.v2022.1.Hs.symbols.gmt")
pathways.lsc <- gmtPathways("/Users/stephanie/Documents/GitHub/EAGLe-2.0/data/gmt_pathway_files copy/lsc_sigs.gmt")
pathways.aeg <- gmtPathways("/Users/stephanie/Documents/GitHub/EAGLe-2.0/data/gmt_pathway_files copy/aeg_genesets_20220602.gmt")

names(pathways.aeg)[10] <- "PM_Primitive_Blast"
names(pathways.aeg)[9] <- "PM_Monocytic_Blast"
pathways.aegGOBP <- c(pathways.aeg, pathways.GObio)
# UI ####
ui <-
  navbarPage(
  "EAGLe: Cancer Discovery",
  navbarMenu( #QC Menu ####
    "QC",
    tabPanel(  # PCA plots ####
      "PCA Plots",
      fluidPage(
        theme =
          shinytheme(
            "flatly"
            ),
        titlePanel(
          "QC Analysis: PCA Plots"
        ), 
        sidebarLayout(
          sidebarPanel(
            
            radioButtons(
              "PCAvar",
              h4(
                "Choose PCA plot"
                ),
              choices =
                list("counts PCA", "VST PCA", "VST + batch corrected PCA"),
              selected =
                "counts PCA"
            ),
            
            downloadButton("downloadPlotPCA", label = "Download Plot"),

          ), #end sidebarPanel
          mainPanel(
            plotOutput(
              "PCAplot"
              )
              )
            ) 
            )
          ),
    tabPanel(  # PCA Scree Plots ####
               "PCA Scree Plots",
               fluidPage(
                 theme =
                   shinytheme(
                     "flatly"
                   ),
                 titlePanel(
                   "QC Analysis: PCA Scree Plots"
                 ), 
                 sidebarLayout(
                   sidebarPanel(
                     radioButtons(
                       "PCAvarscree",
                       h4(
                         "Choose PCA for Scree Plot"
                       ),
                       choices =
                         list("counts PCA", "VST PCA", "VST + batch corrected PCA"),
                       selected =
                         "counts PCA"
                     ),
                     
                     downloadButton(
                       "downloadPlotscree",
                       label =
                         "Download Plot"
                     )
                   ), #end sidebarPanel
                   mainPanel(
                         plotOutput(
                           "PCAvarplot"
                         )
                       )
                   )
                 )
               ), 
    tabPanel( #MultiQC Plots ####
      "MultiQC",
      fluidPage(
        theme = 
          shinytheme("flatly"),
        titlePanel(
          "QC Analysis: MultiQC Plots"
        ),
        sidebarLayout(
          sidebarPanel(
            selectInput(
            "QCvar",
            label=
              "Choose MultiQC test",
            choices =
              c("% mapped reads", "# mapped reads", "% uniquely mapped reads", "# uniquely mapped reads"),
            selected =
              "% mapped reads"
          ) #end selectInput
        ), #end sidebar panel
        mainPanel(
            plotOutput(
              "QCplot"
              )
            )
          )
        )
        )
      ),
    tabPanel( #Gene centric analysis ####
      "Gene Centric Analysis",
      fluidPage(
        theme=
          shinytheme("flatly"),
        titlePanel(
            "Gene Centric Analysis"
            ),
          sidebarLayout(
            sidebarPanel(
              useShinyjs(),
              selectizeInput(
                "VSTCDgenechoice",
                label=
                  "Choose a gene for analysis",
                choices =
                  NULL,
                selected = NULL,
                options = list(maxItems = NULL)
                ),
              radioButtons("XaxisVar_CDgene", h4("X axis variable"),
                           choices = list("Value" = "xvalue", "Class" = "xclass",
                                          "Gene" = "xgene"),selected = "xgene"),
              radioButtons("YaxisVar_CDgene", h4("Y axis variable"),
                           choices = list("Value" = "yvalue", "Class" = "yclass",
                                          "Gene" = "ygene"),selected = "yvalue"),
              radioButtons("FillVar_CDgene", h4("Fill variable"),
                           choices = list("Class" = "fillclass", "Gene" = "fillgene"), selected = "fillclass"),
              
            hr(),
            
            h3(
              "Aesthetics:"
              ),
            hr(),
            materialSwitch(
              inputId =
                "genefacetbutton",
              label =
                "Facet Grid",
              value =
                FALSE,
              right =
                TRUE
            ),
           hr(),
           #Palettes from colorBrewer
           selectInput("PaletteChoices", "Choose a color palette", choices =
                         c("Dark2", "Paired", "Set1"), selected = "Dark2"),
           hr(),
           sliderInput("geneheightslider", "Adjust plot height",
                                   min = 200, max = 1200, value = 600
           ),
           hr(),
           sliderInput("genewidthslider", "Adjust plot width",
                       min = 200, max = 1200, value = 800
           ),
           hr(),
          
           downloadButton("downloadGenePlot", label = "Download Plot"),
           
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
    navbarMenu("DESeq Analysis",# DESeq Menu ####
               tabPanel("DESeq Analysis: Table", #DESeq table ####
                        fluidPage(
                          theme =
                            shinytheme("flatly"),
                          titlePanel(
                            "DESeq Table"
                            ),#end title
                          sidebarLayout(
                            sidebarPanel( 
                              radioButtons("padjbutton", h4("padj Value"), 
                               choices = list("<= 0.01" = "sigvar1", "<= 0.05" = "sigvar5", "All" = "allvar"), selected = "allvar"),
                              
                              hr(),
                              radioButtons("DiffExpButton", h4("Differential Expression"),
                               choices = list("Up" = "DEup", "Down" = "DEdown", "No" = "DEno", "All" = "DEall"), selected = "DEall"),
                              hr(),
                              downloadButton("downloadDEtable", label = "Download Table"),
                              ),
                            mainPanel(
                              DTOutput(
                                "DETable"
                                )
                              )
                          )
                        )
               ),
               tabPanel("DESeq Volcano Plot", #DESeq volcano plot ####
                        fluidPage(
                          theme =
                            shinytheme("flatly"),
                          titlePanel(
                            "DESeq Analysis: Volcano Plot"
                          ),#end title
                          sidebarLayout(
                            sidebarPanel( 
                              radioButtons("sigvaluesbutton", h4("padj Value"), 
                                           choices = list("<= 0.01" = "sigvar0.01", "<= 0.05" = "sigvar0.05", "All" = "allvar2"), selected = "allvar2"),
                              
                              hr(),
                              selectInput("PaletteChoicesDE", "Choose a color palette", choices =
                                            c("Dark2", "Paired", "Set1"), selected = "Dark2"),
                              sliderInput("volheightslider", "Adjust plot height",
                                          min = 200, max = 1000, value = 400
                              ),
                              sliderInput("volwidthslider", "Adjust plot width",
                                          min = 200, max = 1000, value = 600
                              ),
                              # radioButtons("DiffExpButton", h4("Differential Expression"),
                              #              choices = list("Up" = "DEup", "Down" = "DEdown", "No" = "DEno", "All" = "DEall"), selected = "DEall")
                              downloadButton(
                                          "downloadDEVolcano",
                                          label =
                                            "Download Plot"
                                        )
                            ),
                            mainPanel(
                              plotOutput(
                                "DEVolcanoPlot"
                              )
                            )
                          )
                        )
               ),
               tabPanel("DESeq MA Plot", #DESeq MA Plot####
                        fluidPage(
                          theme =
                            shinytheme("flatly"),
                          titlePanel(
                            "DESeq Analysis: MA Plot"
                          ),#end title
                          sidebarLayout(
                            sidebarPanel( 
                              # radioButtons("padjbutton", h4("padj Value"), 
                              #              choices = list("<= 0.01" = "sigvar1", "<= 0.05" = "sigvar5", "All" = "allvar"), selected = "allvar"),
                              # 
                              # hr(),
                              # radioButtons("DiffExpButton", h4("Differential Expression"),
                              #              choices = list("Up" = "DEup", "Down" = "DEdown", "No" = "DEno", "All" = "DEall"), selected = "DEall")
                              downloadButton(
                                "downloadDEMA",
                                label =
                                  "Download Plot"
                              )
                            ),
                            mainPanel(
                              plotOutput(
                                "DEMAPlot"
                              )
                            )
                          )
                        )
               )
              
               ), 
#GSEA menu ####
 navbarMenu( 
   "GSEA"
 ),
#WGCNA menu####
 navbarMenu(
   "WGCNA"
 )
)
#Server ####
server <- 
  function(input, output, session) {
  
  print("Initializing renderPlots")
  options(shiny.reactlog = TRUE)
  
  
##QC-MultiQC plots####
  output$QCplot <- renderPlot({
    QCdata <- switch(
      input$QCvar,
      "% mapped reads" = qcdt$raw.salmon.percent_mapped,
      "# mapped reads" = qcdt$raw.salmon.num_mapped,
      "% uniquely mapped reads" = qcdt$raw.star.uniquely_mapped_percent,
      "# uniquely mapped reads" = qcdt$raw.star.uniquely_mapped,
      "Sample_ID" = qcdt$metadata.sample_id
    )
 
 Sample_ID <-qcdt$metadata.sample_id 
  
  ggplot(
    qcdt,
    aes(
      x = Sample_ID,
      y = QCdata
      )) +
    geom_point(alpha = 0.5) +
    theme_cowplot (12) +
    theme(axis.text.x =
            element_text(angle = 60, hjust = 1))
  }) #end renderPlot
  
  # PCA plots ####
 
  PCAdata <-
    eventReactive(input$PCAvar, {
      if (input$PCAvar == "counts PCA") {
        nonvsd.pca
      } else if (input$PCAvar == "VST PCA") {
        vsd.pca
      } else if (input$PCAvar == "VST + batch corrected PCA") {
        bcvsd.pca
      # } else if (input$PCAvar == "PCA variance") {
      #   PCA_variance
       }
    })
  
  # PCA Scree data ####
  #non VSD PCA variance
  PC_var <-data.frame(PC =paste0("PC", 1:12),variance = (((nonvsd2.pca$sdev) ^ 2 / sum((nonvsd2.pca$sdev) ^ 2)) * 100))
  lorder <-as.vector(outer(c("PC"), 1:12, paste, sep = ""))
  PC_var$PC <-factor(PC_var$PC,levels = lorder)
  
  # VSD PCA variance
  PC_var_VST <- data.frame(PC =paste0("PC", 1:12),variance =(((vsd2.pca$sdev) ^ 2 / sum((vsd2.pca$sdev) ^ 2)) * 100))
  lorder_VST <- as.vector(outer(c("PC"), 1:12, paste, sep = ""))
  PC_var_VST$PC <-factor(PC_var_VST$PC,levels = lorder_VST)
  
  #batch corrected PCA variance
  PC_var_bc <-data.frame(PC =paste0("PC", 1:12),variance =(((bcvsd2.pca$sdev) ^ 2 / sum((bcvsd2.pca$sdev) ^ 2)) * 100))
  lorder_bc <-as.vector(outer(c("PC"), 1:12, paste, sep = ""))
  PC_var_bc$PC <-factor(PC_var_bc$PC,levels = lorder_bc)
  
  PC_var_data <-
    eventReactive(input$PCAvarscree, {
      if (input$PCAvarscree == "counts PCA") {
        PC_var
      } else if (input$PCAvarscree == "VST PCA") {
        PC_var_VST
      } else if (input$PCAvarscree == "VST + batch corrected PCA") {
        PC_var_bc
      }
    })
  # reactive function for plot title
  PCA_title <- 
    reactive({
        if (input$PCAvar == "counts PCA") {
          print("counts PCA")
        } else if (input$PCAvar == "VST PCA") {
          print("VST PCA")
        } else if (input$PCAvar == "VST + batch corrected PCA") {
          print("VST + batch corrected PCA")
        }
    })
  
  #define objects for defining x label to include % variance of PC1
  pc1var <- paste("PC1", (round(nonvsd.variance[3, 1] * 100, 1)), "% variance")
  pc1varvsd <- paste("PC1", (round(vsd.variance[3, 1] * 100, 1)), "% variance")
  pc1varbcvsd <- paste("PC1", (round(bcvsd.variance[3, 1] * 100, 1)), "% variance")
  
  # add reactive expression for x label for PC1 
  variance_PC1 <-
    eventReactive(input$PCAvar, {
      if (input$PCAvar == "counts PCA") {
        pc1var
      } else if (input$PCAvar == "VST PCA") {
        pc1varvsd
      } else if (input$PCAvar == "VST + batch corrected PCA") {
        pc1varbcvsd
      }
    })
  #new objects with calculation only to use in calculation of PC2 % variance
  pc1 <- round(nonvsd.variance[3, 1] * 100, 1)
  pc12 <- round(vsd.variance[3, 1] * 100, 1)
  pc13 <- round(bcvsd.variance[3, 1] * 100, 1)
  
  #create objects for defining Y label of PC2
  pc2var <-
    paste("PC2", (round(nonvsd.variance[3, 2] * 100 - pc1, 1)), "% variance")
  pc2varvsd <-
    paste("PC2", (round(vsd.variance[3, 2] * 100 - pc12, 1)), "% variance")
  pc2varbcvsd <-
    paste("PC2", (round(bcvsd.variance[3, 2] * 100 - pc13, 1)), "% variance")

    #reactive expression for adding y labels for pc2
  variance_PC2 <-
    reactive({
      if (input$PCAvar == "counts PCA") {
        pc2var
      } else if (input$PCAvar == "VST PCA") {
        pc2varvsd
      } else if (input$PCAvar == "VST + batch corrected PCA") {
        pc2varbcvsd
      }
    })
  
  
  output$PCAplot <- renderPlot ({
    colors <-
      colorRampPalette(
        c(
          "#1B9E77",
          "#D95F02",
          "#7570B3",
          "#E7298A",
          "#66A61E",
          "#E6AB02",
          "#A6761D",
          "#666666"
        ))(2)
   ggplot(PCAdata(), aes(x = PC1, y = PC2, fill = batch, shape = condition)) + 
      geom_point(size = 5) + 
      scale_shape_manual(values = c(21, 24), name = '') +
      scale_fill_manual(values = colors) +
      theme_cowplot() + 
      theme(plot.background = element_rect(fill = "#FFFFFF", colour = "#FFFFFF")) +
      theme(panel.background = element_rect(fill = "#FFFFFF", colour = "#FFFFFF")) +
      xlab(variance_PC1()) + 
      ylab(variance_PC2()) +
      ggtitle(PCA_title()) +
      guides(fill=guide_legend(override.aes = list(color=colors))) +
      geom_text_repel(aes(label=sample_name),hjust=0, vjust=0)
 
  })
  
  #PCA plots download ####
  output$downloadPlotPCA <- downloadHandler(
    filename = function() { paste(input$PCAvar, '.png', sep='') },
    content = function(file) {
      ggsave(file, device = "png", width = 8, height = 6, units = "in",dpi = 72)
    }
  )
  # PCA scree plot ####
  output$PCAvarplot <- renderPlot ({
    ggplot(PC_var_data(),
           aes(x = PC,
               y = variance,
               group = 2)) +
      geom_point(size = 2) +
      geom_line() +
      theme_cowplot() +
      labs(x = "PC",
           y = "% Variance") +
      labs(title =
             "")
  })
  #PCA Scree download ####
  output$downloadPlotscree <- downloadHandler(
    filename = function() { paste(input$PCAvarscree, '.png', sep='') },
    content = function(file) {
      ggsave(file, device = "png", width = 8, height = 6, units = "in",dpi = 72)
    }
  )
##Gene Centric output ####
 updateSelectizeInput(session,"VSTCDgenechoice", choices = vst.goi$ext_gene, server = TRUE)
 
 datavst<-
   reactive({
     vst.goi %>% 
       dplyr::filter(ext_gene %in% input$VSTCDgenechoice)
   })
 

#make sure duplicate selections are not allowed with radio buttons
 observeEvent(input$XaxisVar_CDgene, {
   if(input$XaxisVar_CDgene == "xvalue") {
     mychoices <- c("Class" = "yclass", "Gene" = "ygene")
      }else if(input$XaxisVar_CDgene=="xclass") {
     mychoices <- c("Value"="yvalue", "Gene"="ygene")
     } else if(input$XaxisVar_CDgene=="xgene") {
     mychoices <- c("Value" = "yvalue", "Class" = "yclass")
     }
   updateRadioButtons(session, "YaxisVar_CDgene", choices = mychoices)
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
 Gene_facet <- 
   eventReactive(input$genefacetbutton, {
     if(input$genefacetbutton == TRUE) {
       facet_grid(ext_gene ~ class, scales = 'free') 
     } else(NULL)
   })
colorpalettechoices <- 
  eventReactive(input$PalletteChoices, {
    if(input$PaletteChoices == "Dark") {
      scale_color_brewer(palette = "Dark2")
    } else if(input$PaletteChoices == "PurpleGreen") {
      scale_color_brewer(palette = "PRGn")
    } else if(input$PaletteChoices == "RedBlue") {
      scale_color_brewer(palette = "RdBu")
    } else if(input$PaletteChoices == "YellowGreenBlue") {
      scale_color_brewer(palette = "YlGnBu")
    }
  })
 #plot output
  output$VSTCDplot <-
    renderPlot(
      width = function() input$genewidthslider,
      height = function() input$geneheightslider,
      {
 #build a color palette
      # colors <-
      #   colorRampPalette(c("dodgerblue4",
      #                      "darkolivegreen4",
      #                      "darkorchid3",
      #                      "goldenrod1"))(10)
      ggplot(datavst(),
             aes(
               x = .data[[xvar_CDgene()]],
               y =  .data[[yvar_CDgene()]],
               fill = .data[[fillvar_CDgene()]]
             )) +
        geom_boxplot(outlier.shape = NA) +
        scale_fill_brewer(palette = input$PaletteChoices) +
        scale_color_brewer(palette = input$PaletteChoices) +
        geom_point(alpha = 0.5,
                   position = position_jitterdodge(jitter.width = 0.2),
                   aes(color = ext_gene)) + #this needs to be reactive too
        theme_light() +
        Gene_facet() +
        ylab("") +
        xlab("") +
        ggtitle("Gene Expression:Sensitive vs Resistant")
    }) #end render plot

  output$downloadGenePlot <- downloadHandler(
    filename = function() { paste('GeneCentricPlot','.png', sep='') },
    content = function(file) {
      ggsave(file, device = "png", width = 8,
             height = 8, dpi = 72)
    }
  )

  #DESEq #####
  
## DESeq2- Cancer Discovery outputs

#   DElog2 <-
#     reactive({
#       dds.res[dds.res$log2FoldChange >= input$CDlog2foldchangeslider & dds.res$log2FoldChange <= input$CDlog2foldchangeslider, ]
#     })

  #function for sidebar input to create filtered DE table and associated volcano plot
  CD_DE_DT <- 
    reactive({
      if (input$padjbutton == "sigvar1" &
          input$DiffExpButton == "DEup") {
        dds.res %>%
          dplyr::filter(DiffExp == "up" & padj <= 0.01)
      } else if (input$padjbutton == "sigvar5" &
                 input$DiffExpButton == "DEup") {
        dds.res %>%
          dplyr::filter(DiffExp == "up" & padj <= 0.05)
      } else if (input$padjbutton == "sigvar1" &
                 input$DiffExpButton == "DEdown") {
        dds.res %>%
          dplyr::filter(DiffExp == "down" & padj <= 0.01)
      } else if (input$padjbutton == "sigvar5" &
                 input$DiffExpButton == "DEdown") {
        dds.res %>%
          dplyr::filter(DiffExp == "down" & padj <= 0.05)
      } else if (input$padjbutton == "sigvar1" &
                 input$DiffExpButton == "DEno") {
        dds.res %>%
          dplyr::filter(DiffExp == "no" & padj <= 0.01)
      } else if (input$padjbutton == "sigvar5" &
                 input$DiffExpButton == "DEno") {
        dds.res %>%
          dplyr::filter(DiffExp == "no" & padj <= 0.05)
      } else if (input$padjbutton == "allvar" &
                 input$DiffExpButton == "DEall") {
        dds.res
      } else if (input$padjbutton == "allvar" &
                 input$DiffExpButton == "DEup") {
        dds.res %>%
          dplyr::filter(DiffExp == "up" & padj >= 0)
      } else if (input$padjbutton == "allvar" &
                 input$DiffExpButton == "DEdown") {
        dds.res %>%
          dplyr::filter(DiffExp == "down" & padj >= 0)
      } else if (input$padjbutton == "allvar" &
                 input$DiffExpButton == "DEno") {
        dds.res %>%
          dplyr::filter(DiffExp == "no" & padj >= 0)
      } else if (input$padjbutton == "sigvar1" &
                 input$DiffExpButton == "DEall") {
        dds.res %>%
          dplyr::filter(DiffExp == c("up", "down", "no") & padj <= 0.01)
      } else if (input$padjbutton == "sigvar5" &
                 input$DiffExpButton == "DEall") {
        dds.res %>%
          dplyr::filter(DiffExp == c("up", "down", "no") & padj <= 0.05)
      }
    })
 
  vol_sig_values <- 
    reactive({
      if(input$sigvaluesbutton == "sigvar0.05") {
        dds.res %>% 
          dplyr::filter(padj <= 0.05)
      } else if(input$sigvaluesbutton == "sigvar0.01") {
        dds.res %>% 
          dplyr::filter(padj <= 0.01)
      }else if(input$sigvaluesbutton == "allvar2") {
        dds.res %>% 
          dplyr::filter(padj > 0)
      }
      })
  
   output$DETable <-
     renderDataTable({
       CD_DE_DT()
       
     })
        
   output$downloadDEtable <- downloadHandler(
     filename = function() { paste("DESeqTable", '.csv', sep='') },
     content = function(file) {
       write.csv(CD_DE_DT(),file)
     }
   )
   output$DEVolcanoPlot <-
     renderPlot( 
       width = function() input$volwidthslider,
                 height = function() input$volheightslider,
                 {
       # colors <-
       #   colorRampPalette(
       #     c(
       #       "#1B9E77",
       #       "#D95F02",
       #       "#7570B3",
       #       "#E7298A",
       #       "#66A61E",
       #       "#E6AB02",
       #       "#A6761D",
       #       "#666666"
       #     )
       #   )(3)
       # 
       ggplot(vol_sig_values(), aes(
         x = log2FoldChange,
         y = -log10(padj),
         col = DiffExp
       )) +
         geom_point(alpha = 0.5) +
         theme_light() +
         scale_colour_brewer(palette = input$PaletteChoicesDE) +
         ggtitle("DE Volcano Plot") +
         geom_text_repel(
           max.overlaps = 15,
           aes(label = ifelse(
             padj < 5e-20 &
               abs(log2FoldChange) >= 0.5,
             as.character(Gene),
             ""
           )),
           hjust = 0,
           vjust = 0
         ) +
         coord_cartesian(xlim = c(-10, 7))
     })
   
   output$downloadDEVolcano <- downloadHandler(
     filename = function() { paste(input$sigvaluesbutton, '.png', sep='') },
     content = function(file) {
       ggsave(file, device = "png", width = 8, height = 6, units = "in",dpi = 72)
     }
   )
   output$DEMAPlot <- renderPlot ({
     
     ggmaplot(
       dds.res,
       fdr = 0.05,
       fc = (2 ^ 1),
       size = 1.5,
       alpha = 0.7,
       palette =  
         colorRampPalette(
         c(
           "#1B9E77",
           "#D95F02",
           "#7570B3",
           "#E7298A",
           "#66A61E",
           "#E6AB02",
           "#A6761D",
           "#666666"
         )
       )(3),
       legend = NULL,
       top = TRUE,
       title = "DE MA Plot",
       ggtheme = ggplot2::theme_light()
     )
     
   })
   output$downloadDEMA <- downloadHandler(
     filename = function() { paste('DESeqMAplot', '.png', sep='') },
     content = function(file) {
       ggsave(file, device = "png", width = 8, height = 6, units = "in",dpi = 72)
     }
   )
  } #end server
# Run the application 
 shinyApp(ui = ui, server = server)


