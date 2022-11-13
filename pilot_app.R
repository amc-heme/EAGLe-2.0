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
library(fgsea)
library(WGCNA)
library(plotly)
library(BiocParallel)
library(ComplexHeatmap)
#options(shiny.reactlog = TRUE)
#reactlogShow(time = TRUE)

#Data ####
meta_lut_ven <- readRDS("/Users/stephanie/Documents/GitHub/EAGLe-2.0/data/meta_lut_ven.Rds")
qcdt<-load_multiqc("/Users/stephanie/Documents/GitHub/EAGLe-2.0/data/multiqc_data.json", sections="raw") 
vst.goi <- readRDS("/Users/stephanie/Documents/GitHub/EAGLe-2.0/data/vst.goi.rds")
#DESeq data table
dds.res <- readRDS("/Users/stephanie/Documents/GitHub/EAGLe-2.0/data/DEtable.rds")
#sample metadata table
metadata <- read.table(file = "/Users/stephanie/Documents/GitHub/EAGLe-2.0/data/SampleSheetJordanLab.txt")
# DE table with singscore
dds.resscore <- readRDS("/Users/stephanie/Documents/GitHub/EAGLe-2.0/data/dds.resscore.rds")
#tables for PCA
vsd.pca <- readRDS("/Users/stephanie/Documents/GitHub/EAGLe-2.0/data/vsd.pca.rds")
bcvsd.pca <- readRDS("/Users/stephanie/Documents/GitHub/EAGLe-2.0/data/bcvsd.pca.rds")
#tables for variance
vsd.variance <- readRDS("/Users/stephanie/Documents/GitHub/EAGLe-2.0/data/vsd.variance.rds")
bcvsd.variance <- readRDS("/Users/stephanie/Documents/GitHub/EAGLe-2.0/data/bcvsd.variance.rds")

vsd2.pca <- readRDS("/Users/stephanie/Documents/GitHub/EAGLe-2.0/data/vsd2.pca.rds")
bcvsd2.pca <- readRDS("/Users/stephanie/Documents/GitHub/EAGLe-2.0/data/bcvsd2.pca.rds")
#GSEA data table
fgseaResTidy <- readRDS("/Users/stephanie/Documents/GitHub/EAGLe-2.0/data/fgseaResTidy.rds")
fgseaResTidyAll <- readRDS("/Users/stephanie/Documents/GitHub/EAGLe-2.0/data/fgseaResTidyAll.rds")
fgseaResTidyMolec <- readRDS("/Users/stephanie/Documents/GitHub/EAGLe-2.0/data/fgseaResTidyMolec.rds")
fgseaResTidyCC <- readRDS("/Users/stephanie/Documents/GitHub/EAGLe-2.0/data/fgseaResTidyCC.rds")
fgseaResTidyBio <- readRDS("/Users/stephanie/Documents/GitHub/EAGLe-2.0/data/fgseaResTidyBio.rds")
fgseaResTidyTF <- readRDS("/Users/stephanie/Documents/GitHub/EAGLe-2.0/data/fgseaResTidyTF.rds")
fgseaResTidyReg <- readRDS("/Users/stephanie/Documents/GitHub/EAGLe-2.0/data/fgseaResTidyReg.rds")
fgseaResTidyWiki <- readRDS("/Users/stephanie/Documents/GitHub/EAGLe-2.0/data/fgseaResTidyWiki.rds")
fgseaResTidyReactome <- readRDS("/Users/stephanie/Documents/GitHub/EAGLe-2.0/data/fgseaResTidyReactome.rds")
ranks <- readRDS("/Users/stephanie/Documents/GitHub/EAGLe-2.0/data/ranks.rds")
#load pathways
pathways.hallmark <- gmtPathways("/Users/stephanie/Documents/GitHub/EAGLe-2.0/data/gmt_pathway_files copy/h.all.v7.4.symbols.gmt")
pathways.GOall <- gmtPathways("/Users/stephanie/Documents/GitHub/EAGLe-2.0/data/gmt_pathway_files copy/c5.go.v2022.1.Hs.symbols.gmt")
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
                list("VST PCA", "VST + batch corrected PCA"),
              selected =
                "VST PCA"
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
                         list("VST PCA", "VST + batch corrected PCA"),
                       selected =
                         "VST PCA"
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
                              h4("log2FoldChange = Prim/Mono"),
                              radioButtons("padjbutton", h4("padj Value"), 
                               choices = list("<= 0.01" = "sigvar1", "<= 0.05" = "sigvar5", "All" = "allvar"), selected = "allvar"),
                              
                              hr(),
                              radioButtons("DiffExpButton", h4("Differential Expression"),
                               choices = list("Up" = "DEup", "Down" = "DEdown", "No" = "DEno", "All" = "DEall"), selected = "DEall"),
                              
                              hr(),
                              materialSwitch(
                                inputId =
                                  "singscorebutton",
                                label =
                                  "Singscore Regression",
                                value =
                                  FALSE,
                                right =
                                  TRUE
                              ),
                              
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
                              h4("log2FoldChange = Prim/Mono"),
                              radioButtons("sigvaluesbutton", h4("padj Value"), 
                                           choices = list("<= 0.01" = "sigvar0.01", "<= 0.05" = "sigvar0.05", "All" = "allvar2"), selected = "allvar2"),
                              
                              hr(),
                              selectInput("PaletteChoicesDE", "Choose a color palette", choices =
                                            c("Dark2", "Paired", "Set1"), selected = "Dark2"),
                              # sliderInput("volheightslider", "Adjust plot height",
                              #             min = 200, max = 1000, value = 400
                              # ),
                              # sliderInput("volwidthslider", "Adjust plot width",
                              #             min = 200, max = 1000, value = 600
                              # ),
                              # radioButtons("DiffExpButton", h4("Differential Expression"),
                              #              choices = list("Up" = "DEup", "Down" = "DEdown", "No" = "DEno", "All" = "DEall"), selected = "DEall")
                              downloadButton(
                                          "downloadDEVolcano",
                                          label =
                                            "Download Plot"
                                        )
                            ),
                            mainPanel(
                              plotlyOutput(
                                "DEVolcanoPlot"
                              ),
                              # textOutput(
                              #   "gene_name"
                              # )
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
                              h4("log2FoldChange = Prim/Mono"),
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
                              plotlyOutput(
                                "DEMAPlot"
                              )
                            )
                          )
                        )
               )
              
               ), 
#GSEA menu ####
navbarMenu("GSEA", 
  tabPanel("GSEA", ####GSEAtables
            fluidPage(
              theme =
                shinytheme("flatly"),
              titlePanel(
                "GSEA"
              ),#end title
              sidebarLayout(
                sidebarPanel( 
                 selectInput("filechoice", label = "Choose gmt file to load pathways",
                             choices = c(Hallmark = "hallmark", GOall = "goall",GOmolecular = "GOmolec", 
                                         GOcellcomp = "GOcellcomp", GObio = "GObio", TFtargets = "TFtargets",
                                         allRegular = "allReg", Wiki = "wiki", Reactome = "reactome", KEGG = "KEGG",
                                         Positional = "positional", Biocarta = "biocarta", lsc = "lsc", aeg = "aeg")),
                             
                             
                  
                  hr(),
                  selectInput("gseachoice", "Choose a visualization tool",
                              choices =
                                c("fgsea Table" = "fgseaTable", "Waterfall Plot" = "rankedplot", "Moustache Plot" = "moustache",
                                 "Enrichment Plot" = "eplot", "Volcano Plot" = "volcanoplot", "Heatmap" = "heatmap")
                              ),
                 
                 conditionalPanel(
                   condition = "input.gseachoice == 'fgseaTable'",
                   downloadButton(
                     "downloadfgsea",
                     label =
                       "Download Table"
                   )
                 ),
                 
                 conditionalPanel(
                   condition = "input.gseachoice == 'rankedplot'",
                   downloadButton(
                     "downloadranks",
                     label =
                       "Download Plot"
                   ),
                  hr(),
                  sliderInput("howmanypathways", "Choose How Many Pathways to Rank",
                              min = 5, max = 60, value = 25
                  ),
                  
                   sliderInput("rankedheightslider", "Adjust plot height",
                               min = 200, max = 1000, value = 400
                   ),
                   hr(),
                   
                   sliderInput("rankedwidthslider", "Adjust plot width",
                               min = 200, max = 1000, value = 600
                   )
                 ),
                 conditionalPanel(
                   condition = "input.gseachoice == 'moustache'",
                   selectInput("PaletteChoicesMoustache", "Choose a color palette", choices =
                                 c("Dark2", "Paired", "Set1"), selected = "Dark2"),
                   downloadButton(
                     "downloadmoustache",
                     label =
                       "Download Plot"
                   )
                   ),
              
                 conditionalPanel(
                   condition = "input.gseachoice == 'eplot'",

                   
                   radioButtons("topupordownbutton", h4("Top Ranked Up or Down Pathway"), 
                                choices = list("Top Ranked Up Pathway" = "topup", "Top Ranked Down Pathway" = "topdown"), selected = "topup"),
                   downloadButton(
                     "downloadeplot",
                     label =
                       "Download Plot"
                   )
                 ),
                 conditionalPanel(
                   condition = "input.gseachoice == 'volcanoplot'",
                   downloadButton(
                     "downloadvolcano",
                     label =
                       "Download Plot"
                   )
                 ),
                 conditionalPanel(
                   condition = "input.gseachoice == 'heatmap'",
                   downloadButton(
                     "downloadheatmap",
                     label =
                       "Download Plot"
                   )
                 ),
                ),
                              

                mainPanel(
                  conditionalPanel(
                    condition = "input.gseachoice == 'fgseaTable'",
                  DTOutput(
                    "fgseaTable"
                  )
                  ),
                  conditionalPanel(
                    condition = "input.gseachoice == 'rankedplot'",
                    plotOutput(
                    "GSEAranked"
                  )
                  ),
                  conditionalPanel(
                    condition = "input.gseachoice == 'moustache'",
                  plotOutput(
                    "GSEAMoustache"
                  )
                  ),
                  conditionalPanel(
                    condition = "input.gseachoice == 'eplot'",
                    plotOutput(
                      "GSEAenrichment"
                    )
                  ),
                  conditionalPanel(
                    condition = "input.gseachoice == 'volcanoplot'",
                  plotOutput(
                    "GSEAvolcano"
                  )
                  )
                  )
                )
              )
            
  ),
  tabPanel("Gene Centric Pathway Analysis",
           fluidPage(
             theme =
               shinytheme("flatly"),
             titlePanel(
               "GSEA"
             ),#end title
             sidebarLayout(
               sidebarPanel( 
                 selectizeInput(
                   "Pathwaygenechoice",
                   label=
                     "Choose a gene or group of genes",
                   choices =
                     NULL,
                   selected = NULL,
                   options = list(maxItems = NULL)
                 ),
                 hr(),
                 
                 selectInput("genefilechoice", "Choose gmt file to load pathways containing the gene or genes of interest",
                 choices = c(Hallmark = "hallmark", GOall = "GOall", GOmolecular = "GOmolec", 
                             GOcellcomp = "GOcellcomp", GObio = "GObio", TFtargets = "TFtargets",
                             allRegular = "allReg", Wiki = "wiki", Reactome = "reactome", KEGG = "KEGG",
                             Positional = "positional", Biocarta = "biocarta", lsc = "lsc", aeg = "aeg"))
                 ),
               mainPanel(
                 plotOutput(
                   "PathwaysGenePlot"
                 )
               )
           )
             )
   )
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
       if (input$PCAvar == "VST PCA") {
        vsd.pca
      } else if (input$PCAvar == "VST + batch corrected PCA") {
        bcvsd.pca
       }
    })
  
  # PCA Scree data ####
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
      if (input$PCAvarscree == "VST PCA") {
        PC_var_VST
      } else if (input$PCAvarscree == "VST + batch corrected PCA") {
        PC_var_bc
      }
    })
  # reactive function for plot title
  PCA_title <- 
    reactive({
       if (input$PCAvar == "VST PCA") {
          print("VST PCA")
        } else if (input$PCAvar == "VST + batch corrected PCA") {
          print("VST + batch corrected PCA")
        }
    })
  
  #define objects for defining x label to include % variance of PC1
  pc1varvsd <- paste("PC1", (round(vsd.variance[3, 1] * 100, 1)), "% variance")
  pc1varbcvsd <- paste("PC1", (round(bcvsd.variance[3, 1] * 100, 1)), "% variance")
  
  # add reactive expression for x label for PC1 
  variance_PC1 <-
    eventReactive(input$PCAvar, {
     if (input$PCAvar == "VST PCA") {
        pc1varvsd
      } else if (input$PCAvar == "VST + batch corrected PCA") {
        pc1varbcvsd
      }
    })
  #new objects with calculation only to use in calculation of PC2 % variance

  pc12 <- round(vsd.variance[3, 1] * 100, 1)
  pc13 <- round(bcvsd.variance[3, 1] * 100, 1)
  
  #create objects for defining Y label of PC2

  pc2varvsd <-
    paste("PC2", (round(vsd.variance[3, 2] * 100 - pc12, 1)), "% variance")
  pc2varbcvsd <-
    paste("PC2", (round(bcvsd.variance[3, 2] * 100 - pc13, 1)), "% variance")

    #reactive expression for adding y labels for pc2
  variance_PC2 <-
    reactive({
     if (input$PCAvar == "VST PCA") {
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
  
  CD_DE_DT_sing <- 
    reactive({
      if (input$padjbutton == "sigvar1" &
          input$DiffExpButton == "DEup" & input$singscorebutton == TRUE) {
        dds.resscore %>%
          dplyr::filter(DiffExp == "up" & padj <= 0.01)
      } else if (input$padjbutton == "sigvar5" &
                 input$DiffExpButton == "DEup" & input$singscorebutton == TRUE) {
        dds.resscore %>%
          dplyr::filter(DiffExp == "up" & padj <= 0.05)
      } else if (input$padjbutton == "sigvar1" &
                 input$DiffExpButton == "DEdown" & input$singscorebutton == TRUE) {
        dds.resscore %>%
          dplyr::filter(DiffExp == "down" & padj <= 0.01)
      } else if (input$padjbutton == "sigvar5" &
                 input$DiffExpButton == "DEdown" & input$singscorebutton == TRUE) {
        dds.resscore %>%
          dplyr::filter(DiffExp == "down" & padj <= 0.05)
      } else if (input$padjbutton == "sigvar1" &
                 input$DiffExpButton == "DEno" & input$singscorebutton == TRUE) {
        dds.resscore %>%
          dplyr::filter(DiffExp == "no" & padj <= 0.01)
      } else if (input$padjbutton == "sigvar5" &
                 input$DiffExpButton == "DEno" & input$singscorebutton == TRUE) {
        dds.resscore %>%
          dplyr::filter(DiffExp == "no" & padj <= 0.05)
      } else if (input$padjbutton == "allvar" &
                 input$DiffExpButton == "DEall"& input$singscorebutton == TRUE ) {
        dds.resscore
      } else if (input$padjbutton == "allvar" &
                 input$DiffExpButton == "DEup" & input$singscorebutton == TRUE) {
        dds.resscore %>%
          dplyr::filter(DiffExp == "up" & padj >= 0)
      } else if (input$padjbutton == "allvar" &
                 input$DiffExpButton == "DEdown" & input$singscorebutton == TRUE) {
        dds.resscore %>%
          dplyr::filter(DiffExp == "down" & padj >= 0)
      } else if (input$padjbutton == "allvar" &
                 input$DiffExpButton == "DEno" & input$singscorebutton == TRUE) {
        dds.resscore %>%
          dplyr::filter(DiffExp == "no" & padj >= 0)
      } else if (input$padjbutton == "sigvar1" &
                 input$DiffExpButton == "DEall" & input$singscorebutton == TRUE) {
        dds.resscore %>%
          dplyr::filter(DiffExp == c("up", "down", "no") & padj <= 0.01)
      } else if (input$padjbutton == "sigvar5" &
                 input$DiffExpButton == "DEall" & input$singscorebutton == TRUE) {
        dds.resscore %>%
          dplyr::filter(DiffExp == c("up", "down", "no") & padj <= 0.05)
      }
    })
  #object for volcano plot data using DE and singscore tables
  vol_sig_values <- 
    reactive({
      if(input$sigvaluesbutton == "sigvar0.05" ) {
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
  #output DE table with adjustment for singscore
   output$DETable <-
     renderDataTable({
       if(input$singscorebutton == TRUE) {
         CD_DE_DT_sing()
       } else if(input$singscorebutton == FALSE) {
         CD_DE_DT()
       }
       
     })
 
   # download DE table
   output$downloadDEtable <- downloadHandler(
     filename = function() { paste("DESeqTable", '.csv', sep='') },
     content = function(file) {
       write.csv(CD_DE_DT(),file)
     }
   )
   
   output$DEVolcanoPlot <-
     renderPlotly( 
       # width = function() input$volwidthslider,
       #           height = function() input$volheightslider,
                 {
      colors <- c(magma(15)[9], "grey", viridis(15)[10])
      p <- ggplot(vol_sig_values(), aes(
         x = log2FoldChange,
         y = -log10(padj),
         col = DiffExp,
         text = Gene
       )) +
         geom_point(size = 1, alpha = 0.5) +
         theme_light() +
         scale_colour_manual(values = colors) +
         ggtitle("DE Volcano Plot") +
         # geom_text_repel(
         #   max.overlaps = 15,
         #   aes(label = ifelse(
         #     padj < 5e-20 &
         #       abs(log2FoldChange) >= 0.5,
         #     as.character(Gene),
         #     ""
         #   )),
         #   hjust = 0,
         #   vjust = 0
         # ) +
         coord_cartesian(xlim = c(-10, 7))
      ggplotly(p)
     })
   # output$gene_name <- renderText({
   #   if(input$point_hover) {
   #     print(dds.res$Gene)
   #   }
   # })
   output$downloadDEVolcano <- downloadHandler(
     filename = function() { paste(input$sigvaluesbutton, '.png', sep='') },
     content = function(file) {
       ggsave(file, device = "png", width = 8, height = 6, units = "in",dpi = 72)
     }
   )
   output$DEMAPlot <- renderPlotly ({
     
     ma <- ggmaplot(
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
     ggplotly(ma)
   })
   output$downloadDEMA <- downloadHandler(
     filename = function() { paste('DESeqMAplot', '.png', sep='') },
     content = function(file) {
       ggsave(file, device = "png", width = 8, height = 6, units = "in",dpi = 72)
     }
   )
   
####GSEA output ####
#run GSEA for chosen pathway input
  

   #make an object to hold the values of the selectInput for gsea pathway choices
   gsea_file_values <- list("hallmark" = pathways.hallmark,
                            "goall" = pathways.GOall,
                            "GOmolec" = pathways.GOmolec, 
                            "GOcellcomp" = pathways.GOcellcomp,
                            "GObio" = pathways.GObio,
                            "TFtargets" = pathways.TFtargets,
                            "allReg" = pathways.allReg,
                            "wiki" = pathways.Wiki,
                            "reactome" = pathways.Reactome,
                            "KEGG" = pathways.KEGG,
                            "positional" = pathways.Positional,
                            "biocarta" = pathways.Positional,
                            "lsc" = pathways.lsc,
                            "aeg" = pathways.aeg)
   #reactive expression to run fgsea and load results table for each chosen pathway
   gseafile <-
     reactive({
       pathwaygsea <- gsea_file_values[[input$filechoice]]
          fgseaRes <- fgsea::fgsea(pathways = pathwaygsea, stats = ranks, nproc = 1)
          fgseaResTidy <- fgseaRes %>%
            as_tibble() %>%
            arrange(desc(NES))
          fgseaResTidy
       })
 
  #filter Res table for chosen pathway to show in a waterfall plot
   gseafile_waterfall <-
     reactive({
       pathwaygsea <- gsea_file_values[[input$filechoice]]
         fgseaRes <-
           fgsea::fgsea(pathways = pathwaygsea,
                        stats = ranks,
                        nproc = 1)
         fgseaResTidy <- fgseaRes %>%
           as_tibble() %>%
           arrange(desc(NES))
         top15 <- fgseaResTidy %>% 
                   top_n(n = input$howmanypathways, wt = NES)
                 bottom15 <- fgseaResTidy %>%
                   top_n(n = -(input$howmanypathways), wt = NES)
                 fgseaResTidy <- rbind.data.frame(top15, bottom15)
                 fgseaResTidy
       })
  
  #table output for Res tables
   output$fgseaTable <- renderDataTable({
     if (input$gseachoice == "fgseaTable") {
       gseafile()
   }
   })
   # 
   #### GSEA pathway ranks waterfall plot ####

   output$GSEAranked <- renderPlot(
     width = function() input$rankedwidthslider,
     height = function() input$rankedheightslider,
     {
     if(input$gseachoice == "rankedplot") {
       ggplot(gseafile_waterfall(), aes(reorder(pathway, NES), NES)) +
         geom_col(aes(fill=padj<0.05)) +
         coord_flip() +
         labs(x="Pathway", y="Normalized Enrichment Score",
              title="") +
         theme_minimal()
     }
   })

   #### GSEA moustache plot ####
 
   toplotMoustache <-
     reactive({
       pathwaygsea <- gsea_file_values[[input$filechoice]]
       fgseaRes <-
         fgsea::fgsea(pathways = pathwaygsea,
                      stats = ranks,
                      nproc = 1)
       fgseaResTidy <- fgseaRes %>%
         as_tibble() %>%
       arrange(desc(NES))
         toplotMoustache <-
           cbind.data.frame(fgseaResTidy$pathway,
                            fgseaResTidy$NES,
                            fgseaResTidy$padj,
                            fgseaResTidy$pval)
         colnames(toplotMoustache) <- c("pathway", "NES", "padj", "pval")
         toplotMoustache <- toplotMoustache %>%
           mutate(., sig = ifelse(padj <= 0.05, 'yes', 'no'))
     })

    output$GSEAMoustache <- renderPlot({
      if(input$gseachoice == "moustache") {

      m <- ggplot(toplotMoustache(), aes(x = NES, y = padj, color = sig)) + #reactive for each pathway
        geom_point() +
        theme_minimal() +
        xlab('NES') +
        scale_colour_brewer(palette = input$PaletteChoicesMoustache) +
        ylab('adjusted p-value') +
        ggtitle("") + #reactive
        geom_text_repel(aes(label=ifelse(padj<0.05,as.character(pathway),"")),hjust=0,vjust=0)
        coord_cartesian(xlim = c(-3, 3), ylim = c(-0.1, 1))
        print(m)
      }
    })
  #GSEA Enrichment Plots ####
    output$GSEAenrichment <- renderPlot ({
      pathwaygsea <- gsea_file_values[[input$filechoice]]
      if(input$topupordownbutton == "topup") {
      top.UP.path <- as.character(fgseaResTidy[1,1])
      plotEnrichment(pathwaygsea[[top.UP.path]],
                     ranks) + labs(title=top.UP.path)
      } else if(input$topupordownbutton == "topdown") {
      top.DOWN.path <- as.character(fgseaResTidy[nrow(fgseaResTidy), 1])
      plotEnrichment(pathwaygsea[[top.DOWN.path]],
                     ranks) + labs(title=top.DOWN.path)
      }
    })
   
    
 # GSEA Volcano plot ####
    #need to figure out how to access the list of lists to get gene names reactively
    dds.res.pathways <- reactive({
      dds.res %>%
        mutate(., react_path = ifelse(Gene %in% c(input$fgseaTable_rows_selected), 'yes', 'no'))
      dds.res.pathways$react_path <- factor(dds.res.pathways$react_path, levels = c('no','yes'))
    })


    output$GSEAvolcano <- renderPlot ({
      if(input$gseachoice == "volcanoplot") {
      
        
        colors <- c("grey", viridis(15)[10])
        ggplot(
          data = (dds.res.pathways() %>%
                    arrange(., (react_path))),
          aes(
            x = log2FoldChange,
            y = -log10(padj),
            col = react_path
          )
        ) +
          theme_light() +
          geom_point() +
          scale_colour_manual(values = colors) +
          geom_text_repel(
            max.overlaps = 1500,
            aes(label = ifelse(
              Gene %in% (input$fgseaTable_rows_selected) & log2FoldChange > 2.5,
              as.character(Gene),
              ""
            )),
            hjust = 0,
            vjust = 0
          ) +
          theme(
            plot.title = element_text(color = "black", size = 14, face = "bold"),
            axis.title.x = element_text(color = "black", size = 14, face =
                                          "bold"),
            axis.title.y = element_text(color = "black", size = 14, face =
                                          "bold"),
            axis.text.x = element_text(size = 14),
            axis.text.y = element_text(size = 14)
          ) +
          ggtitle("") +
          xlab("log2FoldChange")
        print("plot loaded")
      }
    })
  
      

   # #gene list for gene centric pathway analysis
   updateSelectizeInput(session,"Pathwaygenechoice", choices = dds.res$Gene, server = TRUE)
  } #end server
# Run the application 
 shinyApp(ui = ui, server = server)


