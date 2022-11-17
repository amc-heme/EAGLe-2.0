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
library(InteractiveComplexHeatmap)
library(circlize)
library(colourpicker)
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
vstlimma <- readRDS("/Users/stephanie/Documents/GitHub/EAGLe-2.0/data/vstlimma.rds")

names(pathways.aeg)[10] <- "PM_Primitive_Blast"
names(pathways.aeg)[9] <- "PM_Monocytic_Blast"
pathways.aegGOBP <- c(pathways.aeg, pathways.GObio)
# UI ####
ui <-
  navbarPage(
    "EAGLe: Cancer Discovery",
    tabPanel( #QC Menu ####
              "QC",
              fluidPage(
                theme =
                  shinytheme(
                    "flatly"
                  ),
                titlePanel(
                  "QC Analysis Plots"
                ), 
                sidebarLayout(
                  sidebarPanel(
                    materialSwitch(
                      inputId =
                        "PCAplots",
                      label =
                        "PCA",
                      value =
                        FALSE,
                      right =
                        TRUE
                    ),
                    materialSwitch(
                      inputId =
                        "PCAscreeplots",
                      label =
                        "Scree",
                      value =
                        FALSE,
                      right =
                        TRUE
                    ),
                    
                    condition = "input.PCAplots == true",
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
                    hr(),
                    
                    colourInput(
                      "PCAcolor1",
                      label = "Choose 1st color",
                      value = "#009292",
                      showColour = ("both"),
                      palette = ("square"),
                      allowedCols = NULL,
                      allowTransparent = FALSE,
                      returnName = FALSE,
                      closeOnClick = FALSE
                    ),
                    hr(),
                    
                    colourInput(
                      "PCAcolor2",
                      label = "Choose 2nd color",
                      value = "#FFFF6D",
                      showColour = ("both"),
                      palette = ("square"),
                      allowedCols = NULL,
                      allowTransparent = FALSE,
                      returnName = FALSE,
                      closeOnClick = FALSE
                    ),
                    hr(),
                    
                    downloadButton("downloadPlotPCA", label = "Download Plot"),
                    
                    conditionalPanel(
                      condition = "input.PCAscreeplots == true",
                      
                      downloadButton("downloadPlotscree",
                                     label =
                                       "Download Plot")
                    ) #end conditional
                  ),
                  mainPanel(
                    conditionalPanel(condition = "input.PCAplots == true",
                                     plotOutput("PCAplot")),
                    conditionalPanel(condition = "input.PCAscreeplots == true",
                                     plotOutput("PCAvarplot"))
                  )
                )
              )
    ), 
    
    
    # tabPanel( #MultiQC Plots ####
    #   "MultiQC",
    #   fluidPage(
    #     theme = 
    #       shinytheme("flatly"),
    #     titlePanel(
    #       "QC Analysis: MultiQC Plots"
    #     ),
    #     sidebarLayout(
    #       sidebarPanel(
    #         selectInput(
    #         "QCvar",
    #         label=
    #           "Choose MultiQC test",
    #         choices =
    #           c("% mapped reads", "# mapped reads", "% uniquely mapped reads", "# uniquely mapped reads"),
    #         selected =
    #           "% mapped reads"
    #       ) #end selectInput
    #     ), #end sidebar panel
    #     mainPanel(
    #         plotOutput(
    #           "QCplot"
    #           )
    #         )
    #       )
    #     )
    #     )
    #   ),
    tabPanel( #Gene centric analysis ####
              "Gene Expression",
              fluidPage(
                theme=
                  shinytheme("flatly"),
                titlePanel(
                  "Gene Expression Plot"
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
                                 choices = list("Value" = "xvalue",
                                                "Gene" = "xgene"),selected = "xgene"),
                    radioButtons("YaxisVar_CDgene", h4("Y axis variable"),
                                 choices = list("Value" = "yvalue",
                                                "Gene" = "ygene"),selected = "yvalue"),
                    
        
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
                    # selectInput("PaletteChoices", "Choose a color palette", choices =
                    #               c("Dark2", "Paired", "Set1"), selected = "Dark2"),
                    colourInput(
                      "genecolor1",
                      label = "Choose 1st color",
                      value = "#009292",
                      showColour = ("both"),
                      palette = ("square"),
                      allowedCols = NULL,
                      allowTransparent = FALSE,
                      returnName = FALSE,
                      closeOnClick = FALSE
                    ),
                    colourInput(
                      "genecolor2",
                      label = "Choose 2nd color",
                      value = "#ffff6d",
                      showColour = ("both"),
                      palette = ("square"),
                      allowedCols = NULL,
                      allowTransparent = FALSE,
                      returnName = FALSE,
                      closeOnClick = FALSE
                    ),
                    colourInput(
                      "genecolor3",
                      label = "Choose 3rd color",
                      value = "#490092",
                      showColour = ("both"),
                      palette = ("square"),
                      allowedCols = NULL,
                      allowTransparent = FALSE,
                      returnName = FALSE,
                      closeOnClick = FALSE
                    ),
                    
                    
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
                    plotOutput(
                      "VSTCDplot"
                    )
                  )
                )
              )
    ),
    
    tabPanel("DESeq Analysis",# DESeq Menu ####
             fluidPage(
               theme =
                 shinytheme("flatly"),
               titlePanel(
                 "DESeq Table and Plots"
               ),#end title
               sidebarLayout(
                 sidebarPanel( 
                   materialSwitch(
                     inputId =
                       "DESeqtable",
                     label =
                       "DE Table",
                     value =
                       FALSE,
                     right =
                       TRUE
                   ),
                   materialSwitch(
                     inputId =
                       "singscorebutton",
                     label =
                       "DE Table w/ Monocytic Contribution",
                     value =
                       FALSE,
                     right =
                       TRUE
                   ),
                   materialSwitch(
                     inputId =
                       "DESeqvolcano",
                     label =
                       "Volcano Plot",
                     value =
                       FALSE,
                     right =
                       TRUE
                   ),
                   materialSwitch(
                     inputId =
                       "DESeqMA",
                     label =
                       "MA Plot",
                     value =
                       FALSE,
                     right =
                       TRUE
                   ),
                   materialSwitch(
                     inputId =
                       "DESeqHeat",
                     label =
                       "Heatmap",
                     value =
                       FALSE,
                     right =
                       TRUE
                   ),
                   hr(),
                   
                   radioButtons("padjbutton", label = "Filter DE tables by padj", 
                                choices = list("<= 0.01" = "sigvar1", "<= 0.05" = "sigvar5", "All" = "allvar"), selected = "allvar"), 
                   hr(),
                   
                   h4("Aesthetics"),
                   
                   colourInput(
                     "volcanocolor1",
                     label = "Choose 1st color",
                     value = "#009292",
                     showColour = ("both"),
                     palette = ("square"),
                     allowedCols = NULL,
                     allowTransparent = FALSE,
                     returnName = FALSE,
                     closeOnClick = FALSE
                   ),
                   colourInput(
                     "volcanocolor2",
                     label = "Choose 2nd color",
                     value = "grey",
                     showColour = ("both"),
                     palette = ("square"),
                     allowedCols = NULL,
                     allowTransparent = FALSE,
                     returnName = FALSE,
                     closeOnClick = FALSE
                   ),
                   colourInput(
                     "volcanocolor3",
                     label = "Choose 3rd color",
                     value = "#490092",
                     showColour = ("both"),
                     palette = ("square"),
                     allowedCols = NULL,
                     allowTransparent = FALSE,
                     returnName = FALSE,
                     closeOnClick = FALSE
                   ),
                   hr(),
                   h4("Table and Plot Downloads"),
                   downloadButton("downloadDEtable", label = "Download DE Table"),
                   hr(),
                   
                   downloadButton(
                     "downloadDEVolcano",
                     label =
                       "Download Volcano Plot"
                   ),
                   hr(),
                   downloadButton(
                     "downloadDEMA",
                     label =
                       "Download MA Plot"
                   )
                   
                 ),
                 mainPanel(
                   conditionalPanel(
                     condition = "input.DESeqtable == true",
                     DTOutput(
                       "DETable"
                     )
                   ),
                   conditionalPanel(
                     condition = "input.DESeqvolcano == true",
                     plotlyOutput(
                       "DEVolcanoPlot"
                     )
                   ),
                   conditionalPanel(
                     condition = "input.DESeqMA == true",
                     plotlyOutput(
                       "DEMAPlot"
                     )
                   ),
                   conditionalPanel(
                     condition = "input.DESeqHeat == true",
                     InteractiveComplexHeatmapOutput(heatmap_id = 
                                                       "ht"
                     )
                   )
                 )
               )
             )
    ),
    #GSEA menu ####
    tabPanel("GSEA",  ####GSEAtables
             fluidPage(
               theme =
                 shinytheme("flatly"),
               titlePanel(
                 "GSEA"
               ),#end title
               sidebarLayout(
                 sidebarPanel( 
                   materialSwitch(
                     inputId =
                       "fgseaTable",
                     label =
                       "Table",
                     value =
                       FALSE,
                     right =
                       TRUE
                   ),
                   materialSwitch(
                     inputId =
                       "rankedplot",
                     label =
                       "Waterfall",
                     value =
                       FALSE,
                     right =
                       TRUE
                   ),
                   materialSwitch(
                     inputId =
                       "moustache",
                     label =
                       "Moustache",
                     value =
                       FALSE,
                     right =
                       TRUE
                   ),
                   materialSwitch(
                     inputId =
                       "eplot",
                     label =
                       "Enrichment",
                     value =
                       FALSE,
                     right =
                       TRUE
                   ),
                   materialSwitch(
                     inputId =
                       "volcanoplot",
                     label =
                       "Volcano",
                     value =
                       FALSE,
                     right =
                       TRUE
                   ),
                   materialSwitch(
                     inputId =
                       "heatmap",
                     label =
                       "Heatmap",
                     value =
                       FALSE,
                     right =
                       TRUE
                   ),
                   
                   
                   selectInput("filechoice", label = "Choose gmt file to load pathway sets",
                               choices = c(Hallmark = "hallmark", GOall = "goall",GOmolecular = "GOmolec", 
                                           GOcellcomp = "GOcellcomp", GObio = "GObio", TFtargets = "TFtargets",
                                           allRegular = "allReg", Wiki = "wiki", Reactome = "reactome", KEGG = "KEGG",
                                           Positional = "positional", Biocarta = "biocarta", lsc = "lsc", aeg = "aeg")),
                   
                   
                   
                   hr(),
                   
                   selectizeInput(
                     "pathwaylist",
                     label=
                       "Choose a pathway",
                     choices =
                       NULL,
                     selected = NULL,
                     options = list(maxItems = NULL)
                   ),
                   
                   conditionalPanel(
                     condition = "input.fgseaTable == true",
                     downloadButton(
                       "downloadfgsea",
                       label =
                         "Download GSEA Table"
                     )
                   ),
                   
                   conditionalPanel(
                     condition = "input.rankedplot == true",
                  
                     h4("Waterfall Plot Specific Options"),
                     
                     sliderInput("howmanypathways", "Choose How Many Pathways to Rank",
                                 min = 5, max = 60, value = 15
                     ),
                     
                     sliderInput("rankedheightslider", "Adjust plot height",
                                 min = 200, max = 1000, value = 400
                     ),
                     hr(),
                     
                     sliderInput("rankedwidthslider", "Adjust plot width",
                                 min = 200, max = 1000, value = 600
                     ),
                     colourInput(
                       "waterfallcolor1",
                       label = "Choose 1st color",
                       value = "#009292",
                       showColour = ("both"),
                       palette = ("square"),
                       allowedCols = NULL,
                       allowTransparent = FALSE,
                       returnName = FALSE,
                       closeOnClick = FALSE
                     ),
                     colourInput(
                       "waterfallcolor2",
                       label = "Choose 2nd color",
                       value = "#490092",
                       showColour = ("both"),
                       palette = ("square"),
                       allowedCols = NULL,
                       allowTransparent = FALSE,
                       returnName = FALSE,
                       closeOnClick = FALSE
                     ),
                     downloadButton(
                       "downloadranks",
                       label =
                         "Download Waterfall Plot"
                     ),
                   ),
                   conditionalPanel(
                     condition = "input.moustache == true",
                     h4("Moustache Plot Specific Options"),
                     
                     colourInput(
                       "choice1color",
                       label = "Choose 1st color",
                       value = "#009292",
                       showColour = ("both"),
                       palette = ("square"),
                       allowedCols = NULL,
                       allowTransparent = FALSE,
                       returnName = FALSE,
                       closeOnClick = FALSE
                     ),
                     
                     colourInput(
                       "choice2color",
                       label = "Choose 2nd color",
                       value = "#000000",
                       showColour = ("both"),
                       palette = ("square"),
                       allowedCols = NULL,
                       allowTransparent = FALSE,
                       returnName = FALSE,
                       closeOnClick = FALSE
                     ),
    
                     downloadButton(
                       "downloadmoustache",
                       label =
                         "Download Moustache Plot"
                     )
                   ),
                   
                   conditionalPanel(
                     condition = "input.eplot == true",
                     
                     h4("Enrichment Plot Specific Options"),
                     radioButtons("topupordownbutton", h4("Top Ranked Up or Down Pathway"), 
                                  choices = list("Top Ranked Up Pathway" = "topup", "Top Ranked Down Pathway" = "topdown"), selected = "topup"),
                     downloadButton(
                       "downloadeplot",
                       label =
                         "Download Enrichment Plot"
                     )
                   ),
                   conditionalPanel(
                     condition = "input.volcanoplot == true",
                     h4("Volcano Plot Specific Options"),
                     colourInput(
                       "gseavolcolor1",
                       label = "Choose 1st color",
                       value = "grey",
                       showColour = ("both"),
                       palette = ("square"),
                       allowedCols = NULL,
                       allowTransparent = FALSE,
                       returnName = FALSE,
                       closeOnClick = FALSE
                     ),
                     colourInput(
                       "gseavolcolor2",
                       label = "Choose 2nd color",
                       value = "#490092",
                       showColour = ("both"),
                       palette = ("square"),
                       allowedCols = NULL,
                       allowTransparent = FALSE,
                       returnName = FALSE,
                       closeOnClick = FALSE
                     ),
                     colourInput(
                       "gseavolcolor3",
                       label = "Choose 3rd color",
                       value = "#009292",
                       showColour = ("both"),
                       palette = ("square"),
                       allowedCols = NULL,
                       allowTransparent = FALSE,
                       returnName = FALSE,
                       closeOnClick = FALSE
                     ),
                     hr(),
                     downloadButton(
                       "downloadvolcano",
                       label =
                         "Download Volcano Plot"
                     )
                   ),
                   conditionalPanel(
                     condition = "input.heatmap == true",
                     
                     downloadButton(
                       "downloadheatmap",
                       label =
                         "Download Heatmap"
                     )
                   ),
                 ),
                 
                 
                 mainPanel(
                   textOutput(
                     "genelist"
                   ),
                   conditionalPanel(
                     condition = "input.fgseaTable == true",
                     DTOutput(
                       "fgseaTable"
                     )
                   ),
                   conditionalPanel(
                     condition = "input.rankedplot == true",
                     plotOutput(
                       "GSEAranked"
                     )
                   ),
                   conditionalPanel(
                     condition = "input.moustache == true",
                     plotOutput(
                       "GSEAMoustache"
                     )
                   ),
                   conditionalPanel(
                     condition = "input.eplot == true",
                     plotOutput(
                       "GSEAenrichment"
                     )
                   ),
                   conditionalPanel(
                     condition = "input.volcanoplot == true",
                     plotOutput(
                       "GSEAvolcano"
                     )
                   ),
                   conditionalPanel(
                     condition = "input.heatmap == true",
                     InteractiveComplexHeatmapOutput(heatmap_id = "htgsea")
                   )
                 )
               )
             )
             
    ),
    tabPanel("GSEA Pathway/Gene Visualization",
             fluidPage(
               theme =
                 shinytheme("flatly"),
               titlePanel(
                 "GSEA:Pathway/Gene Visualization"
               ),#end title
               h4("Positive NES is upregulated in Primitive cells and negative NES is upregulated in Monocytic cells"),
               sidebarLayout( 
                 sidebarPanel( 
                   selectizeInput(
                     "Pathwaygenechoice",
                     label=
                       "Choose a gene of interest",
                     choices =
                       NULL,
                     selected = NULL,
                     options = list(maxItems = 1)
                   ),
                   hr(),
                   
                   selectInput("genefilechoice", "Choose gmt file to load pathways containing the gene of interest",
                               choices = c(Hallmark = "hallmark", GOall = "GOall", GOmolecular = "GOmolec", 
                                           GOcellcomp = "GOcellcomp", GObio = "GObio", TFtargets = "TFtargets",
                                           allRegular = "allReg", Wiki = "wiki", Reactome = "reactome", KEGG = "KEGG",
                                           Positional = "positional", Biocarta = "biocarta", lsc = "lsc", aeg = "aeg")),
                   
                   hr(),
                   sliderInput("goiheightslider", "Adjust plot height",
                               min = 200, max = 1000, value = 400
                   ),
                   hr(),
                   
                   sliderInput("goiwidthslider", "Adjust plot width",
                               min = 200, max = 1000, value = 600
                   )
                 ),
                 mainPanel(
                   plotOutput(
                     "PathwaysGenePlot"
                   )
                 )
               )
             )
             
    ),
    
    #WGCNA menu####
    tabPanel(
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
      eventReactive(input$PCAvar, {
        if (input$PCAvar == "VST PCA") {
          PC_var_VST
        } else if (input$PCAvar == "VST + batch corrected PCA") {
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
        c(input$PCAcolor1, input$PCAcolor2)
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
    PCA_var_title <- 
      reactive({
        if (input$PCAvar == "VST PCA") {
          print("VST PC variance")
        } else if (input$PCAvar == "VST + batch corrected PCA") {
          print("VST + batch corrected PC variance")
        }
      })
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
               PCA_var_title())
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
        mychoices <- c("Gene" = "ygene")
      } else if(input$XaxisVar_CDgene=="xgene") {
        mychoices <- c("Value" = "yvalue")
      }
      updateRadioButtons(session, "YaxisVar_CDgene", choices = mychoices)
    })
    
    
    #x axis output
    xvar_CDgene <-
      eventReactive(input$XaxisVar_CDgene, {
        if (input$XaxisVar_CDgene == "xvalue") {
          "value"
        } else if (input$XaxisVar_CDgene == "xgene") {
          "ext_gene"
        }
      })
    #y axis output f
    yvar_CDgene <-
      eventReactive(input$YaxisVar_CDgene, {
        if (input$YaxisVar_CDgene == "yvalue") {
          "value"
        } else if (input$YaxisVar_CDgene == "ygene") {
          "ext_gene"
        }
      })
    #fill output
    # fillvar_CDgene <-
    #   eventReactive(input$FillVar_CDgene, {
    #     if (input$FillVar_CDgene == "fillclass") {
    #       "class"
    #     } else if (input$FillVar_CDgene == "fillgene") {
    #       "ext_gene"
    #     }
    #   })
    Gene_facet <- 
      eventReactive(input$genefacetbutton, {
        if(input$genefacetbutton == TRUE) {
          facet_grid(ext_gene ~ class, scales = 'free') 
        } else(NULL)
      })
    # colorpalettechoices <- 
    #   eventReactive(input$PalletteChoices, {
    #     if(input$PaletteChoices == "Dark") {
    #       scale_color_brewer(palette = "Dark2")
    #     } else if(input$PaletteChoices == "PurpleGreen") {
    #       scale_color_brewer(palette = "PRGn")
    #     } else if(input$PaletteChoices == "RedBlue") {
    #       scale_color_brewer(palette = "RdBu")
    #     } else if(input$PaletteChoices == "YellowGreenBlue") {
    #       scale_color_brewer(palette = "YlGnBu")
    #     }
    #   })
    #plot output
    output$VSTCDplot <-
      renderPlot(
        width = function() input$genewidthslider,
        height = function() input$geneheightslider,
        {
          #build a color palette
          colors <-
            colorRampPalette(c(input$genecolor1, input$genecolor2, input$genecolor3))(10)
          ggplot(datavst(),
                 aes(
                   x = .data[[xvar_CDgene()]],
                   y =  .data[[yvar_CDgene()]],
                   fill = class
                 )) +
            geom_boxplot(outlier.shape = NA) +
            scale_fill_manual(values = colors) +
            scale_color_manual(values = colors) +
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
        if (input$padjbutton == "sigvar1") {
          dds.res %>%
            dplyr::filter(padj <= 0.01)
        } else if (input$padjbutton == "sigvar5") {
          dds.res %>%
            dplyr::filter(padj <= 0.05)
        } else if (input$padjbutton == "sigvar1") {
          dds.res %>%
            dplyr::filter(padj <= 0.01)
        } else if (input$padjbutton == "sigvar5") {
          dds.res %>%
            dplyr::filter(padj <= 0.05)
        } else if (input$padjbutton == "sigvar1") {
          dds.res %>%
            dplyr::filter(padj <= 0.01)
        } else if (input$padjbutton == "sigvar5") {
          dds.res %>%
            dplyr::filter(padj <= 0.05)
        } else if (input$padjbutton == "allvar") {
          dds.res
        } else if (input$padjbutton == "sigvar1") {
          dds.res %>%
            dplyr::filter(padj <= 0.01)
        } else if (input$padjbutton == "sigvar5") {
          dds.res %>%
            dplyr::filter(padj <= 0.05)
        } 
      })
    
    CD_DE_DT_sing <- 
      reactive({
        if (input$padjbutton == "sigvar1" & input$singscorebutton == TRUE) {
          dds.resscore %>%
            dplyr::filter(padj <= 0.01)
        } else if (input$padjbutton == "sigvar5" & input$singscorebutton == TRUE) {
          dds.resscore %>%
            dplyr::filter(padj <= 0.05)
        } else if (input$padjbutton == "sigvar1" & input$singscorebutton == TRUE) {
          dds.resscore %>%
            dplyr::filter(padj <= 0.01)
        } else if (input$padjbutton == "sigvar5" & input$singscorebutton == TRUE) {
          dds.resscore %>%
            dplyr::filter(padj <= 0.05)
        } else if (input$padjbutton == "sigvar1" & input$singscorebutton == TRUE) {
          dds.resscore %>%
            dplyr::filter(padj <= 0.01)
        } else if (input$padjbutton == "sigvar5" & input$singscorebutton == TRUE) {
          dds.resscore %>%
            dplyr::filter(padj <= 0.05)
        } else if (input$padjbutton == "allvar" & input$singscorebutton == TRUE ) {
          dds.resscore
        } else if (input$padjbutton == "sigvar1" & input$singscorebutton == TRUE) {
          dds.resscore %>%
            dplyr::filter(padj <= 0.01)
        } else if (input$padjbutton == "sigvar5" & input$singscorebutton == TRUE) {
          dds.resscore %>%
            dplyr::filter(padj <= 0.05)
        }
      })
    #object for volcano plot data using DE and singscore tables
    # vol_sig_values <- 
    #   reactive({
    #     if(input$sigvaluesbutton == "sigvar0.05" ) {
    #       dds.res %>% 
    #         dplyr::filter(padj <= 0.05)
    #     } else if(input$sigvaluesbutton == "sigvar0.01") {
    #       dds.res %>% 
    #         dplyr::filter(padj <= 0.01)
    #     }else if(input$sigvaluesbutton == "allvar2") {
    #       dds.res %>% 
    #         dplyr::filter(padj > 0)
    #     } 
    #     })
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
    #DE Volcano Plot ####
    
    output$DEVolcanoPlot <-
      renderPlotly({
        colors <- c(input$volcanocolor1, input$volcanocolor2, input$volcanocolor3)
        p <- ggplot(dds.res, aes(
          x = `log2FoldChange(Prim/Mono)`,
          y = -log10(padj),
          col = DiffExp,
          text = Gene
        )) +
          geom_point(size = 1, alpha = 0.5) +
          theme_light() +
          scale_colour_manual(values = colors) +
          ggtitle("DE Volcano Plot") +
          coord_cartesian(xlim = c(-10, 7))
        ggplotly(p)
      })
    
    output$downloadDEVolcano <- downloadHandler(
      filename = function() { paste(input$sigvaluesbutton, '.png', sep='') },
      content = function(file) {
        ggsave(file, device = "png", width = 8, height = 6, units = "in",dpi = 72)
      }
    )
    
    #DE MA Plot ####
    output$DEMAPlot <- renderPlotly ({
      
      ma <- ggmaplot(
        dds.res,
        fdr = 0.05,
        fc = 2 ^ 1,
        size = 1.5,
        alpha = 0.7,
        palette =  
          c(input$volcanocolor1, input$volcanocolor2, input$volcanocolor3),
        legend = NULL,
        top = TRUE,
        title = "DE MA Plot",
        ggtheme = ggplot2::theme_light())
      ggplotly(ma)
    })
    
    output$downloadDEMA <- downloadHandler(
      filename = function() { paste('DESeqMAplot', '.png', sep='') },
      content = function(file) {
        ggsave(file, device = "png", width = 8, height = 6, units = "in",dpi = 72)
      }
    )
    #DE Heatmap ####
    dds.mat <- dds.res %>%
      dplyr::filter(padj < 0.05 & abs(`log2FoldChange(Prim/Mono)`) >= 2)
    
    vst.mat <- vstlimma %>%
      dplyr::filter(., ensembl_gene_id %in% dds.mat$ensembl_gene_id) %>%
      column_to_rownames(., var = "ensembl_gene_id") %>%
      dplyr::select(.,-ext_gene) %>%
      as.matrix()
    rownames(vst.mat) = dds.mat$Gene
    vst.mat <- t(scale(t(vst.mat)))
    
    vst.mat <- head(vst.mat, n = 100)
    f1 = colorRamp2(seq(min(vst.mat), max(vst.mat), length = 3), c("blue", "#EEEEEE", "red"))
    ht = draw(ComplexHeatmap::Heatmap(
      vst.mat,
      name = "z scaled expression",
      col = f1,
      row_names_gp = gpar(fontsize = 4),
      row_km = 2,
      top_annotation = HeatmapAnnotation(class = anno_block(gp = gpar(fill = c("darkorange1", "blueviolet")),
                                                            labels = c("prim", "mono"), 
                                                            labels_gp = gpar(col = "white", fontsize = 10))),
      column_km = 2, 
      column_title = NULL,
      row_title = NULL
    ))
    makeInteractiveComplexHeatmap(input, output, session, ht, "ht")
    
    
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
    
    gene_gsea_file_values <- list("hallmark" = pathways.hallmark,
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
          dplyr::select(., -pval,-log2err, -ES) %>% 
          arrange(desc(NES))
        fgseaResTidy
      })
 
    
    observe({
      pathwaygsea <- gsea_file_values[[input$filechoice]]
      updateSelectizeInput(session,"pathwaylist", choices = names(pathwaygsea), server = TRUE)})
    
  
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
          dplyr::select(., -pval,-log2err, -ES) %>% 
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
      if (input$fgseaTable == TRUE) {
        gseafile()
      }
    })
    # 
    #### GSEA pathway ranks waterfall plot ####
    
    output$GSEAranked <- renderPlot(
      width = function() input$rankedwidthslider,
      height = function() input$rankedheightslider,
      {
        if(input$rankedplot == TRUE) {
          colors <- c(input$waterfallcolor1, input$waterfallcolor2)
          ggplot(gseafile_waterfall(), aes(reorder(pathway, NES), NES)) +
            geom_col(aes(fill= padj < 0.05)) +
            scale_fill_manual(values = colors) +
            coord_flip() +
            labs(x="Pathway", y="Normalized Enrichment Score (Positive NES = Upregulated in Primitive
Negative NES = Upregulated in Monocytic)",
                 title="Pathway NES from GSEA") +
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
        fgseaResTidy <-
          fgseaRes %>%
          as_tibble() %>%
          dplyr::select(., -pval,-log2err, -ES) %>% 
          arrange(desc(NES))
        toplotMoustache <-
          cbind.data.frame(fgseaResTidy$pathway,
                           fgseaResTidy$NES,
                           fgseaResTidy$padj)
        colnames(toplotMoustache) <- c("pathway", "NES", "padj")
        toplotMoustache <- toplotMoustache %>%
          mutate(., sig = ifelse(padj <= 0.05, 'yes', 'no'))
      })
    
    output$GSEAMoustache <- renderPlot(
      # width = function()
      #   input$mwidthslider,
      # height = function()
      #   input$mheightslider,
      {
        if (input$moustache == TRUE) {
          colors <- c(input$choice1color, input$choice2color)
          m <-
            ggplot(toplotMoustache(), aes(x = NES, y = padj, color = sig)) +
            geom_point() +
            theme_minimal() +
            xlab('NES') +
            scale_colour_manual(values = colors) +
            ylab('adjusted p-value') +
            ggtitle("Pathways from GSEA") + #reactive
            # geom_text_repel(
            #   max.overlaps = 10,
            #   aes(label = ifelse(
            #     padj < 0.05, as.character(fgseaRes$pathway), ""
            #   )),
            #   hjust = 0,
            #   vjust = 0
            # )
            coord_cartesian(xlim = c(-3, 3), ylim = c(-0.1, 1))
          print(m)
        }
      }
    )
    
    #GSEA Enrichment Plots ####
    output$GSEAenrichment <- renderPlot ({
      pathwaygsea <- gsea_file_values[[input$filechoice]]
      fgseaRes <-
        fgsea::fgsea(pathways = pathwaygsea,
                     stats = ranks,
                     nproc = 1)
      fgseaResTidy <- fgseaRes %>%
        as_tibble() %>%
        dplyr::select(., -pval,-log2err, -ES) %>% 
        arrange(desc(NES))
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
      pathwaygsea <- gsea_file_values[[input$filechoice]]
      p <-
        unlist((pathwaygsea[names(pathwaygsea) %in% input$pathwaylist]))
      dds.res.pathways <- dds.res %>%
        mutate(., react_path = ifelse(Gene %in% p, 'yes', 'no'))
    })
    
    
    output$GSEAvolcano <- renderPlot ({
      colors <- c(input$gseavolcolor1, input$gseavolcolor2, input$gseavolcolor3)
      if (input$volcanoplot == TRUE) {
        ggplot(
          data = (dds.res.pathways() %>% arrange(., (react_path))),
          aes(
            x = `log2FoldChange(Prim/Mono)`,
            y = -log10(padj),
            col = react_path
          )
        ) +
          theme_light() +
          geom_point() +
          scale_colour_manual(values = colors) +
          # geom_text_repel(
          #   max.overlaps = 1500,
          #   aes(
          #     label = ifelse(
          #       Gene %in% p & `log2FoldChange(Prim/Mono)` > 1.5,
          #       as.character(Gene),
          #       ""
          #     )
          #   ),
          #   hjust = 0,
          #   vjust = 0
        # ) +
        ggtitle("") +
          xlab("log2foldchange")
      }
    })
   
    #GSEA heatmap ####
   
    observeEvent(input$pathwaylist, {
      pathwaygsea <- gsea_file_values[[input$filechoice]]
      
      p <- unlist((pathwaygsea[names(pathwaygsea) %in% input$pathwaylist]))
      
      dds.sig <- dds.res %>%
        dplyr::filter(padj < 0.05 & abs(`log2FoldChange(Prim/Mono)`) >= 0.5)
      
      vst.myc <- vstlimma %>% 
        mutate(., Hallmark_myc = ifelse(ext_gene %in% p, 'yes', 'no')) %>% 
        dplyr::filter(Hallmark_myc == "yes")
      
      vstgsea.mat <- vst.myc %>%
        dplyr::filter(., ensembl_gene_id %in% dds.sig$ensembl_gene_id) %>%
        column_to_rownames(., var = "ext_gene") %>%
        dplyr::select(.,-ensembl_gene_id, -Hallmark_myc) %>%
        as.matrix()
      
      vstgsea.mat <- t(scale(t(vstgsea.mat)))
      
      htgsea = draw(ComplexHeatmap::Heatmap(
        vstgsea.mat,
        name = "gseaht_title()",
        row_names_gp = gpar(fontsize = 6),
        column_title = NULL,
        row_title = NULL))
      
      makeInteractiveComplexHeatmap(input, output, session, htgsea, "htgsea")
    })
    #Gene Centeric pathways analysis plots ####
    genecentricgseaplot <- reactive({
      genepathwaygsea <- (gene_gsea_file_values[[input$genefilechoice]])
      fgseaRes <-
        fgsea::fgsea(pathways = genepathwaygsea,
                     stats = ranks,
                     nproc = 1)
      fgseaResTidy <- fgseaRes %>%
        as_tibble() %>%
        arrange(desc(NES))
      goi_paths <- genepathwaygsea %>% keep(grepl(input$Pathwaygenechoice, genepathwaygsea))
      goi_paths <- list(grep(input$Pathwaygenechoice, genepathwaygsea))
      goi_paths <- fgseaResTidy %>%
        dplyr::filter(grepl(input$Pathwaygenechoice, leadingEdge)) 
      GOI <- input$Pathwaygenechoice
      goi_paths$GOI <- "Yes"
      nongoi_paths <- fgseaResTidy %>%
        dplyr::filter(!grepl(input$Pathwaygenechoice, leadingEdge))  
      nongoi_paths$GOI <- "No"
      allgoi_paths <- rbind.data.frame(goi_paths, nongoi_paths)
    })
    
    
    output$PathwaysGenePlot <- renderPlot(
      width = function() input$goiwidthslider,
      height = function() input$goiheightslider,
      {
        
        ggplot(genecentricgseaplot(), aes(
          x = NES,
          y = NES,
          color = (padj < 0.05)
        )) +
          geom_boxplot()  +
          facet_wrap( ~ GOI, scales = "free") +
          theme_light() +
          geom_hline(yintercept = 0, linetype = "dashed")
      })
    
    # #gene list for gene centric pathway analysis
    updateSelectizeInput(session,"Pathwaygenechoice", choices = dds.res$Gene, server = TRUE)
  } #end server
# Run the application 
shinyApp(ui = ui, server = server)

