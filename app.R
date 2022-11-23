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
library(ggsci)
library(scales)
library(esquisse)
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
                    # palettePicker(
                    #   inputId = "PaletteChoices", 
                    #   label = "Choose a color palette", 
                    #   choices = list(
                    #     "Viridis" = list(
                    #       "viridis" = viridis_pal(option = "viridis")(5),
                    #       "magma" = viridis_pal(option = "magma")(5),
                    #       "mako" = viridis_pal(option = "mako")(5),
                    #       "plasma" = viridis_pal(option = "plasma")(5),
                    #       "cividis" = viridis_pal(option = "cividis")(5))
                    #   )
                    # ),
                    # colorPicker(
                    #   inputId = "col1QC",
                    #   label = "PCA plot: Choose 1st color:",
                    #   choices = c(scales::viridis_pal(option = "viridis")(5)[1:5], scales::brewer_pal(palette = "Dark2")(5)[1:5]),
                    #   textColor = "white"
                    # ),
                    # colorPicker(
                    #   inputId = "col2QC",
                    #   label = "PCA plot: Choose 2nd color:",
                    #   choices = c(scales::viridis_pal(option = "viridis")(5)[2:5],scales::brewer_pal(palette = "Dark2")(5)[1:5]),
                    #   textColor = "white"
                    # ),
                    colourInput(
                      "col1QC",
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
                      "col2QC",
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
                    
                    downloadButton("downloadPlotPCA", label = "Download PCA Plot"),
                    
                    hr(),
                    
                    conditionalPanel(
                      condition = "input.PCAscreeplots == true",
                      
                      downloadButton("downloadPlotscree",
                                     label =
                                       "Download Scree Plot")
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
    tabPanel( #Gene expression analysis ####
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
                                                "Gene" = "xgene", "Class" = "xclass"),selected = "xgene"),
                    radioButtons("YaxisVar_CDgene", h4("Y axis variable"),
                                 choices = list("Value" = "yvalue",
                                                "Gene" = "ygene"),selected = "yvalue"),
                    radioButtons("FillVar_CDgene", h4("color by:"),
                                 choices = list("Gene" = "fillgene",
                                                "Class" = "fillclass"),selected = "fillclass"),
                    radioButtons("PrimMonobutton", h4("Show only prim or mono gene expression"),
                                choices = list("Show Comparison" = "comparison", "Prim" = "prim", "Mono" = "mono"), selected = "comparison"),
                    hr(),
                    
                    materialSwitch("genefacetbutton", label = "Facet", value = FALSE, right = TRUE),
                    
        
                    hr(),
                    palettePicker(
                      inputId = "PaletteChoicesGene", 
                      label = "Choose a color palette", 
                      choices = list(
                        "Viridis" = list(
                          "viridis" = viridis_pal(option = "viridis")(5),
                          "magma" = viridis_pal(option = "magma")(5),
                          "mako" = viridis_pal(option = "mako")(5),
                          "plasma" = viridis_pal(option = "plasma")(5),
                          "cividis" = viridis_pal(option = "cividis")(5))
                      )
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
                       "DE Table without Monocytic Contribution",
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
                   conditionalPanel(
                     condition = "input.DESeqtable == true",
                    h4("DE Table Specific Options"),
                    
                   radioButtons("padjbutton", label = "Filter DE tables by padj", 
                                choices = list("<= 0.01" = "sigvar1", "<= 0.05" = "sigvar5", "All" = "allvar"), selected = "allvar"),
                   hr(),
                   downloadButton("downloadDEtable", label = "Download DE Table")
                   ),
                   
                   hr(),
 
                   conditionalPanel(
                     condition = "input.DESeqvolcano == true",
                     h4("Volcano Plot Specific Options"),
                     
                     colorPicker(
                       inputId = "col1",
                       label = "Volcano plot: Choose 1st color:",
                       choices = c(scales::viridis_pal(option = "viridis")(5), scales::brewer_pal(palette = "Dark2")(5)),
                       textColor = "white"
                     ),

                     colorPicker(
                       inputId = "col2",
                       label = "Volcano plot: Choose 2nd color:",
                       choices = c(scales::viridis_pal(option = "viridis")(5)[2:5],scales::brewer_pal(palette = "Dark2")(5)),
                       textColor = "white"
                     ),
                     hr(),
                     downloadButton(
                       "downloadDEVolcano",
                       label =
                         "Download Volcano Plot"
                     )
                   ),
                   hr(),
                   
                   conditionalPanel(
                     condition = "input.DESeqMA ==true", 
                     h4("MA Plot Specific Options"),
                     
                     colorPicker(
                       inputId = "MAcol1",
                       label = "MA plot: Choose 1st color:",
                       choices = c(scales::viridis_pal(option = "viridis")(5), scales::brewer_pal(palette = "Dark2")(5)),
                       textColor = "white"
                     ),
         
                     colorPicker(
                       inputId = "MAcol2",
                       label = "MA plot: Choose 2nd color:",
                       choices = c(scales::viridis_pal(option = "viridis")(5)[2:5],scales::brewer_pal(palette = "Dark2")(5)),
                       textColor = "white"
                     ),
                     hr(),
                     downloadButton(
                       "downloadDEMA",
                       label =
                         "Download MA Plot"
                     )
                   ),
                   hr(),
                   
                   conditionalPanel(
                     condition = "input.DESeqHeat == true",
                     h4("Heatmap Specific Options"),
                     
                     colourInput(
                       "heatcolor1",
                       label = "Choose 1st color",
                       value = "red",
                       showColour = ("both"),
                       palette = ("square"),
                       allowedCols = NULL,
                       allowTransparent = FALSE,
                       returnName = FALSE,
                       closeOnClick = FALSE
                     ),
                     colourInput(
                       "heatcolor2",
                       label = "Choose 2nd color",
                       value = "blue",
                       showColour = ("both"),
                       palette = ("square"),
                       allowedCols = NULL,
                       allowTransparent = FALSE,
                       returnName = FALSE,
                       closeOnClick = FALSE
                     )

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
                   h4("Choose gmt file to load pathway sets"),
                   
                   selectInput("filechoice", label = NULL,
                               choices = c(Hallmark = "hallmark", GOall = "goall",GOmolecular = "GOmolec", 
                                           GOcellcomp = "GOcellcomp", GObio = "GObio", TFtargets = "TFtargets",
                                           allRegular = "allReg", Wiki = "wiki", Reactome = "reactome", KEGG = "KEGG",
                                           Positional = "positional", Biocarta = "biocarta", lsc = "lsc", aeg = "aeg")),
                   h4("Choose a plot type"),
                   
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
                   
                   hr(),
                   # palettePicker(
                   #   inputId = "PaletteChoicesGSEA", 
                   #   label = "Choose a color palette", 
                   #   choices = list(
                   #     "Viridis" = list(
                   #       "viridis" = viridis_pal(option = "viridis")(5),
                   #       "magma" = viridis_pal(option = "magma")(5),
                   #       "mako" = viridis_pal(option = "mako")(5),
                   #       "plasma" = viridis_pal(option = "plasma")(5),
                   #       "cividis" = viridis_pal(option = "cividis")(5))
                   #   )
                   # ),
                   
                   conditionalPanel(
                     condition = "input.fgseaTable == true",
                     downloadButton(
                       "downloadfgsea",
                       label =
                         "Download GSEA Table"
                     )
                   ),
                   hr(),
                   #GSEA Waterfall####
                   conditionalPanel(
                     condition = "input.rankedplot == true",
                  
                     h4("Waterfall Plot Specific Options"),
                     hr(),
                     colourInput(
                       "colWF",
                       label = "Choose color",
                       value = "#490092",
                       showColour = ("both"),
                       palette = ("square"),
                       allowedCols = NULL,
                       allowTransparent = FALSE,
                       returnName = FALSE,
                       closeOnClick = FALSE
                     ),
                     hr(), 
                     
                     sliderInput("howmanypathways", "Choose How Many Pathways to Rank",
                                 min = 5, max = 60, value = 15
                     ),
                     hr(),
                     
                     sliderInput("rankedheightslider", "Adjust plot height",
                                 min = 200, max = 1000, value = 400
                     ),
                     hr(),
                     
                     sliderInput("rankedwidthslider", "Adjust plot width",
                                 min = 200, max = 1000, value = 600
                     ),
    
                     downloadButton(
                       "downloadranks",
                       label =
                         "Download Waterfall Plot"
                     ),
                   ),
                   hr(),
                   #GSEA Moustache ####
                   conditionalPanel(
                     condition = "input.moustache == true",
                     h4("Moustache Plot Specific Options"),
                     #h5("Choose pathway(s) to label on plot"),
                     hr(),
                     # colorPicker(
                     #   inputId = "colMoustache",
                     #   label = "Moustache plot: Choose color:",
                     #   choices = c(scales::viridis_pal(option = "viridis")(5)[1:5], scales::brewer_pal(palette = "Dark2")(5)[1:5]),
                     #   textColor = "white"
                     # ),
                     colourInput(
                       "colMoustache",
                       label = "Choose color",
                       value = "#490092",
                       showColour = ("both"),
                       palette = ("square"),
                       allowedCols = NULL,
                       allowTransparent = FALSE,
                       returnName = FALSE,
                       closeOnClick = FALSE
                     ),
                    hr(), 
    
                     downloadButton(
                       "downloadmoustache",
                       label =
                         "Download Moustache Plot"
                     )
                   ),
                   
                   hr(),
                   #GSEA Enrichment plot ####
                   conditionalPanel(
                     condition = "input.eplot == true",
                     
                     h4("Enrichment Plot Specific Options"),
      
                     radioButtons("topupordownbutton", h5("Enrichment Plot Choices"), 
                                  choices = list("Top Ranked Up Pathway" = "topup", "Top Ranked Down Pathway" = "topdown", "Pathway of Choice:" = "eplotpath"), selected = "topup"),
                     h5("Choose a specific pathway"),
                     selectizeInput(
                       "pathwaylisteplot",
                       label=
                         NULL,
                       choices =
                         NULL,
                       selected = NULL ,
                       options = list(maxItems = 1)
                     ),
                     hr(),
                     downloadButton(
                       "downloadeplot",
                       label =
                         "Download Enrichment Plot"
                     )
                   ),
                   hr(),
                   #GSEA Volcano ####
                   conditionalPanel(
                     condition = "input.volcanoplot == true",
                     
                     h4("Volcano Plot Specific Options"),
                     selectizeInput(
                       "pathwaylist",
                       label=
                         "Choose a specific pathway(s) to view on volcano plot",
                       choices =
                         NULL,
                       selected = NULL,
                       options = list(maxItems = 1)
                     ),
                
                     colourInput(
                       "colVolcano",
                       label = "Choose color",
                       value = "#490092",
                       showColour = ("both"),
                       palette = ("square"),
                       allowedCols = NULL,
                       allowTransparent = FALSE,
                       returnName = FALSE,
                       closeOnClick = FALSE
                     ),
                     
                     downloadButton(
                       "downloadvolcano",
                       label =
                         "Download Volcano Plot"
                     )
                   ),
                   hr(),
                   #GSEA Heatmap ####
                   conditionalPanel(
                     condition = "input.heatmap == true",
                     h4("Heatmap Specific Options"),

                     selectizeInput(
                       "pathwaylistht",
                       label=
                         "Choose a specific pathway to view genes on heatmap",
                       choices =
                         NULL,
                       selected = NULL ,
                       options = list(maxItems = 1)
                     ),
                     colourInput(
                       "GSEAheatcolor1",
                       label = "Choose 1st color",
                       value = "red",
                       showColour = ("both"),
                       palette = ("square"),
                       allowedCols = NULL,
                       allowTransparent = FALSE,
                       returnName = FALSE,
                       closeOnClick = FALSE
                     ),

                     colourInput(
                       "GSEAheatcolor2",
                       label = "Choose 2nd color",
                       value = "blue",
                       showColour = ("both"),
                       palette = ("square"),
                       allowedCols = NULL,
                       allowTransparent = FALSE,
                       returnName = FALSE,
                       closeOnClick = FALSE
                     ),
                     
                     hr(),
                     
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
                   
                   selectInput("genefilechoice", "Choose gmt file to plot pathways containing the gene of interest",
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
    
    #color palette choices ####
    colorpalettechoices <-
      eventReactive(input$PaletteChoices, {
        if(input$PaletteChoices == "viridis") {
           scale_fill_viridis_d(option = "viridis")
        } else if(input$PaletteChoices == "cividis") {
          scale_fill_viridis_d(option = "cividis")
        } else if(input$PaletteChoices == "magma") {
          scale_fill_viridis_d(option = "magma")
        } else if(input$PaletteChoices == "plasma") {
          scale_fill_viridis_d(option = "plasma")
        }else if(input$PaletteChoices == "inferno") {
          scale_fill_viridis_d(option = "inferno")
        }
      })

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
      colors <- c(input$col1QC, input$col2QC)
      ggplot(PCAdata(), aes(x = PC1, y = PC2, fill = batch, shape = condition)) + 
        geom_point(size = 5) + 
        scale_shape_manual(values = c(21, 24), name = '') +
        scale_fill_manual(values = colors) +
        theme_cowplot(font_size = 18) + 
        theme(axis.title = element_text(face = "bold"), title = element_text(face = "bold")) +
        theme(plot.background = element_rect(fill = "#FFFFFF", colour = "#FFFFFF")) +
        theme(panel.background = element_rect(fill = "#FFFFFF", colour = "#FFFFFF")) +
        xlab(variance_PC1()) + 
        ylab(variance_PC2()) +
        ggtitle(PCA_title()) +
        guides(fill=guide_legend(override.aes = list(color = colors))) +
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
        theme_cowplot(font_size = 18) +
        #theme_light(base_size = 18) +
        theme(axis.title = element_text(face = "bold"), title = element_text(face = "bold")) +
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
    
    datavst <-
      reactive({
        if(input$PrimMonobutton == "comparison") {
          vst.goi %>% 
            dplyr::filter(ext_gene %in% input$VSTCDgenechoice)
        } else if(input$PrimMonobutton == "prim") {
          vst.goi %>%
            dplyr::filter(ext_gene %in% input$VSTCDgenechoice) %>%
            dplyr::filter(class == "sensitive")
        } else if(input$PrimMonobutton == "mono") {
          vst.goi %>%
            dplyr::filter(ext_gene %in% input$VSTCDgenechoice) %>%
            dplyr::filter(class == "resistant")
        }
      })
    
    
    #make sure duplicate selections are not allowed with radio buttons
    observeEvent(input$XaxisVar_CDgene, {
      if(input$XaxisVar_CDgene == "xvalue") {
        mychoices <- c("Gene" = "ygene")
      } else if(input$XaxisVar_CDgene=="xgene") {
        mychoices <- c("Value" = "yvalue")
      }else if(input$XaxisVar_CDgene == "xclass") {
        mychoices <- c("Gene" = "ygene", "Value" = "yvalue")
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
        } else if(input$XaxisVar_CDgene == "xclass") {
          "class"
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
          facet_grid(cols = vars(class))
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
    
    colorpalettechoicesfgene <-
      eventReactive(input$PaletteChoicesGene, {
        if(input$PaletteChoicesGene == "viridis") {
          scale_fill_viridis_d(option = "viridis")
        } else if(input$PaletteChoicesGene == "cividis") {
          scale_fill_viridis_d(option = "cividis")
        } else if(input$PaletteChoicesGene == "magma") {
          scale_fill_viridis_d(option = "magma")
        } else if(input$PaletteChoicesGene == "plasma") {
          scale_fill_viridis_d(option = "plasma")
        }else if(input$PaletteChoicesGene == "inferno") {
          scale_fill_viridis_d(option = "inferno")
        }
      })
    colorpalettechoicesgene <-
      eventReactive(input$PaletteChoicesGene, {
        if(input$PaletteChoicesGene == "viridis") {
          scale_color_viridis_d(option = "viridis")
        } else if(input$PaletteChoicesGene == "cividis") {
          scale_color_viridis_d(option = "cividis")
        } else if(input$PaletteChoicesGene == "magma") {
          scale_color_viridis_d(option = "magma")
        } else if(input$PaletteChoicesGene == "plasma") {
          scale_color_viridis_d(option = "plasma")
        }else if(input$PaletteChoicesGene == "inferno") {
          scale_fill_viridis_d(option = "inferno")
        }
      })
    #plot output
    output$VSTCDplot <-
      renderPlot(
        width = function() input$genewidthslider,
        height = function() input$geneheightslider,
        res = 120,
        {

          ggplot(datavst(),
                 aes(
                   x = .data[[xvar_CDgene()]],
                   y =  .data[[yvar_CDgene()]],
                   fill = .data[[fillvar_CDgene()]]
                 )) +
            geom_boxplot(outlier.shape = NA) +
            Gene_facet() +
            stat_compare_means(aes(label=..p.adj..)) +
            # scale_colour_manual(values = colors) +
            # scale_fill_manual(values = colors) +
            colorpalettechoicesfgene() +
            colorpalettechoicesgene() +
            geom_point(alpha = 0.5,
                       position = position_jitterdodge(jitter.width = 0.2),
                       aes(color = class)) + #this needs to be reactive too
            theme_light() +
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
    
    colorpalettechoicesDE <-
      eventReactive(input$PaletteChoicesDE, {
        if(input$PaletteChoicesDE == "viridis") {
          scale_color_viridis_d(option = "viridis")
        } else if(input$PaletteChoicesDE == "cividis") {
          scale_color_viridis_d(option = "cividis")
        } else if(input$PaletteChoicesDE == "magma") {
          scale_color_viridis_d(option = "magma")
        } else if(input$PaletteChoicesDE == "plasma") {
          scale_color_viridis_d(option = "plasma")
        } else if(input$PaletteChoicesDE == "inferno") {
          scale_color_viridis_d(option = "inferno")
        }
      })

    #DE Volcano Plot ####
    
    output$DEVolcanoPlot <-
      renderPlotly({
        colors <- c(input$col1, "grey", input$col2)
        p <- ggplot(dds.res, aes(
          x = `log2FoldChange(Prim/Mono)`,
          y = -log10(padj),
          col = DiffExp,
          text = Gene
        )) +
          geom_point(size = 1, alpha = 0.5) +
          theme_light() +
          scale_color_manual(values = colors) +
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
        colors <- c(input$col1, "grey", input$col2)
      ma <-
        ggplot(dds.res,
               aes(
                 x = log2(baseMean),
                 y = `log2FoldChange(Prim/Mono)`,
                 col = DiffExp
               )) +
        geom_point(alpha = 0.8, size = 0.5) +
        geom_hline(aes(yintercept = 0)) +
        scale_color_manual(values = colors) +
        theme_light() +
        ylim(c(
          min(dds.res$`log2FoldChange(Prim/Mono)`),
          max(dds.res$`log2FoldChange(Prim/Mono)`)
        )) +
        ggtitle("DE MA Plot") +
        xlab("log2 Mean Expression") +
        ylab("Log2 Fold Change")
      
      ggplotly(ma)
    })
    
    output$downloadDEMA <- downloadHandler(
      filename = function() { paste('DESeqMAplot', '.png', sep='') },
      content = function(file) {
        ggsave(file, device = "png", width = 8, height = 6, units = "in",dpi = 72)
      }
    )
    
    heatpalettechoicesDE <-
      eventReactive(input$PaletteChoicesDE, {
        if(input$PaletteChoicesDE == "viridis") {
          viridis_pal(option = "viridis")(10)
        } else if(input$PaletteChoicesDE == "cividis") {
          viridis_pal(option = "cividis")(10)
        } else if(input$PaletteChoicesDE == "magma") {
          viridis_pal(option = "magma")(10)
        } else if(input$PaletteChoicesDE == "plasma") {
          viridis_pal(option = "plasma")(10)
        } else if(input$PaletteChoicesDE == "inferno") {
          viridis_pal(option = "inferno")(10)
        }
      })
    #DE Heatmap ####
    observe({
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
    colors = colorRamp2(c(-2, 0, 2), c(input$heatcolor1, "white", input$heatcolor2))
    ht = draw(ComplexHeatmap::Heatmap(
      vst.mat,
      name = "z scaled expression",
      col = colors,
      row_names_gp = gpar(fontsize = 4),
      row_km = 2,
      top_annotation = HeatmapAnnotation(class = anno_block(gp = gpar(fill = c("white", "white")),
                                                            labels = c("prim", "mono"), 
                                                            labels_gp = gpar(col = "black", fontsize = 10))),
      column_km = 2, 
      column_title = NULL,
      row_title = NULL
    ))
    makeInteractiveComplexHeatmap(input, output, session, ht, "ht")
    })
    
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
        fgseaRes <- fgsea::fgsea(pathways = pathwaygsea, stats = ranks, nproc = 5)
        fgseaResTidy <- fgseaRes %>%
          as_tibble() %>%
          dplyr::select(., -pval,-log2err, -ES) %>% 
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
     colorpalettechoicesGSEAfill <-
    eventReactive(input$PaletteChoicesGSEA, {
      if(input$PaletteChoicesGSEA == "viridis") {
        scale_fill_viridis_d(option = "viridis")
      } else if(input$PaletteChoicesGSEA == "cividis") {
        scale_fill_viridis_d(option = "cividis")
      } else if(input$PaletteChoicesGSEA == "magma") {
        scale_fill_viridis_d(option = "magma")
      } else if(input$PaletteChoicesGSEA == "plasma") {
        scale_fill_viridis_d(option = "plasma")
      }else if(input$PaletteChoicesGSEA == "inferno") {
        scale_fill_viridis_d(option = "inferno")
      }
    })
    colorpalettechoicesGSEA <-
      eventReactive(input$PaletteChoicesGSEA, {
        if(input$PaletteChoicesGSEA == "viridis") {
          scale_color_viridis_d(option = "viridis")
        } else if(input$PaletteChoicesGSEA == "cividis") {
          scale_color_viridis_d(option = "cividis")
        } else if(input$PaletteChoicesGSEA == "magma") {
          scale_color_viridis_d(option = "magma")
        } else if(input$PaletteChoicesGSEA == "plasma") {
          scale_color_viridis_d(option = "plasma")
        }else if(input$PaletteChoicesGSEA == "inferno") {
          scale_fill_viridis_d(option = "inferno")
        }
      })
    
    #### GSEA pathway ranks waterfall plot ####
    
    output$GSEAranked <- renderPlot(
      width = function()
        input$rankedwidthslider,
      height = function()
        input$rankedheightslider,
      {
        if (input$rankedplot == TRUE) {
          colors <- c("grey", input$colWF)
          ggplot(gseafile_waterfall(), aes(reorder(pathway, NES), NES)) +
            geom_col(aes(fill = padj < 0.05)) +
            scale_fill_manual(values = colors) +
            coord_flip() +
            labs(
              x = "Pathway",
              y = "Normalized Enrichment Score (Positive NES = Upregulated in Primitive
Negative NES = Upregulated in Monocytic)",
              title = "Pathway NES from GSEA"
            ) +
            theme_minimal(base_size = 16) +
             theme(axis.title = element_text(face = "bold"), title = element_text(face = "bold")) 
        }
      }
    )
    
    #### GSEA moustache plot ####
    gseamoustache_title <-
      eventReactive(input$filechoice, {
        print(input$filechoice)
      })
    
    
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
    # updateSelectizeInput(session,"pathwaylistmoustache", choices = fgseaResTidy$pathway, server = TRUE)
    output$GSEAMoustache <- renderPlot(
      # width = function()
      #   input$mwidthslider,
      # height = function()
      #   input$mheightslider,
      {
        if (input$moustache == TRUE) {
          colors <- c('grey', input$colMoustache)
          m <-
            ggplot(toplotMoustache(), aes(x = NES, y = padj, color = sig)) +
            geom_point() +
            theme_minimal(base_size = 18) +
            theme(axis.title = element_text(face = "bold"), title = element_text(face = "bold")) +
            xlab('NES') +
            scale_color_manual(values = colors) +
            ylab('adjusted p-value') +
            ggtitle("Pathways from GSEA") + 
            #geom_text() +
            geom_text_repel(colour = "black", aes(label= ifelse(padj <0.05, as.character(pathway), ""), hjust=0,vjust=0)) +
            coord_cartesian(xlim = c(-3, 3), ylim = c(-0.1, 1)) 
          print(m)
        }
      }
    )
    
    #GSEA Enrichment Plots ####
    observe({
      pathwaygsea <- gsea_file_values[[input$filechoice]]
      updateSelectizeInput(session,"pathwaylisteplot", choices = names(pathwaygsea), server = TRUE)})
    
    gseaeplot_title <-
      eventReactive(input$pathwaylisteplot, {
        print(input$pathwaylisteplot)
      })
    
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
      } else if(input$topupordownbutton == "eplotpath") {
        plotEnrichment(pathwaygsea[[input$pathwaylisteplot]],
                       ranks) + labs(title= input$pathwaylisteplot)
      }
    })
    # output$GSEAenrichmentpath <- renderPlot ({
    #   if(input$eplot == TRUE) {
    #   pathwaygsea <- gsea_file_values[[input$filechoice]]
    #   p <-
    #     unlist((pathwaygsea[names(pathwaygsea) %in% input$pathwaylisteplot]))
    #   plotEnrichment(pathwaygsea[[input$pathwaylisteplot]],
    #                  ranks) + labs(title= input$pathwaylisteplot)
    #   }
    # })
    
    # GSEA Volcano plot ####
    observe({
      pathwaygsea <- gsea_file_values[[input$filechoice]]
      updateSelectizeInput(session,"pathwaylist", choices = names(pathwaygsea), server = TRUE)})
    
    dds.res.pathways <- reactive({
      pathwaygsea <- gsea_file_values[[input$filechoice]]
      p <-
        unlist((pathwaygsea[names(pathwaygsea) %in% input$pathwaylist]))

      dds.res.pathways <- dds.res %>%
        mutate(., genes_in_pathway = ifelse(Gene %in% p, 'yes', 'no'))
 print(dds.res.pathways)
    })
    
    gseavol_title <-
      eventReactive(input$pathwaylist, {
        paste(input$pathwaylist)
      })
   
    output$GSEAvolcano <- renderPlot ({
      colors <- 
        c("grey", input$colVolcano)
      if (input$volcanoplot == TRUE) {
        ggplot(
          data = (dds.res.pathways() %>% arrange(., (genes_in_pathway))),
          aes(
            x = `log2FoldChange(Prim/Mono)`,
            y = -log10(padj),
            col = genes_in_pathway[]
          )
        ) +
          theme_light(base_size = 14) +
          theme(axis.title = element_text(face = "bold"), title = element_text(face = "bold")) +
          geom_point() +
          scale_color_manual(values = colors) +
        geom_text_repel(
          max.overlaps = 1500,
          colour = "black",
          aes(
            label = ifelse(
              genes_in_pathway == 'yes' & `log2FoldChange(Prim/Mono)` > 1.5,
              as.character(Gene),
              ""
            )
          ),
          hjust = 0,
          vjust = 0
        ) +
        ggtitle(gseavol_title()) +
          xlab("log2foldchange")
      }
    })
   
    #GSEA heatmap ####
    observe({
      pathwaygsea <- gsea_file_values[[input$filechoice]]
      updateSelectizeInput(session,"pathwaylistht", choices = names(pathwaygsea), server = TRUE)})
    
    gseaht_title <-
      eventReactive(input$pathwaylistht, {
       print(input$pathwaylistht)
      })
    
    observeEvent(input$pathwaylistht, {
      if(input$heatmap == TRUE) {
      pathwaygsea <- gsea_file_values[[input$filechoice]]
      
      p <- unlist((pathwaygsea[names(pathwaygsea) %in% input$pathwaylistht]))
      
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
      colors = colorRamp2(c(-2, 0, 2), c(input$GSEAheatcolor1, "white", input$GSEAheatcolor2))
      htgsea = draw(ComplexHeatmap::Heatmap(
        vstgsea.mat,
        name = paste(gseaht_title(), fontsize = 6),
        col = colors,
          #"paste(input$pathwaylist, sep = ",")",
        row_names_gp = gpar(fontsize = 6),
        column_km = 2,
        top_annotation = HeatmapAnnotation(class = anno_block(gp = gpar(fill = c("white", "white")),
                                                              labels = c("prim", "mono"), 
                                                              labels_gp = gpar(col = "black", fontsize = 10))),
        column_title = NULL,
        row_title = NULL))
      
      makeInteractiveComplexHeatmap(input, output, session, htgsea, "htgsea")
      }
    })
    gsea_gene_title <- 
      eventReactive(input$Pathwaygenechoice, {
        paste(input$Pathwaygenechoice)
      })
    #Gene Centeric pathways analysis plots ####
    genecentricgseaplot <- reactive({
      #load chosen pathway file based on reactive input 
      genepathwaygsea <- (gene_gsea_file_values[[input$genefilechoice]])
      #load fgsea table data for chosen pathway
      fgseaRes <-
        fgsea::fgsea(pathways = genepathwaygsea,
                     stats = ranks,
                     nproc = 1)
      #create tidy table
      fgseaResTidy <- fgseaRes %>%
        as_tibble() %>%
        arrange(desc(NES))
  

      #create object for storing pathways that contain the chosen gene 
      goi_paths <- genepathwaygsea %>% keep(grepl(input$Pathwaygenechoice, genepathwaygsea))
      goi_paths <- list(grep(input$Pathwaygenechoice, genepathwaygsea))
      #filter gsea table for pathways in which the GOI is in the leading edge
      goi_paths <- fgseaResTidy %>%
        dplyr::filter(grepl(input$Pathwaygenechoice, leadingEdge)) %>%
        mutate(., class = ifelse(NES <0, 'Mono', 'Prim'))

      #create object for gene reactive input
      GOI <- input$Pathwaygenechoice
      #make a column that says yes if goi in that pathway
      goi_paths$GOI <- "Pathways that do include the gene"
      #filter gsea table to find pathways that do not include the GOI in the leading edge
      nongoi_paths <- fgseaResTidy %>%
        dplyr::filter(!grepl(input$Pathwaygenechoice, leadingEdge)) %>%
        mutate(., class = ifelse(NES <0, 'Mono', 'Prim'))

      #put no for pathways that do not contain the goi
      nongoi_paths$GOI <- "Pathways that do NOT inlcude the gene"
      #bind the two filtered data frames into one for plotting
      allgoi_paths <- rbind.data.frame(goi_paths, nongoi_paths)
    })
    
    output$PathwaysGenePlot <- renderPlot(
      width = function() input$goiwidthslider,
      height = function() input$goiheightslider,
      {
        
        ggplot(genecentricgseaplot(), aes(
          x = class,
          y = NES,
          color = (padj < 0.05)
        )) +
          geom_boxplot()  +
          facet_wrap( ~ GOI, scales = "free") +
          theme_light(base_size = 18) +
          theme(axis.title = element_text(face = "bold"), title = element_text(face = "bold")) +
          geom_hline(yintercept = 0, linetype = "dashed")
      })
    
    # #gene list for gene centric pathway analysis
    updateSelectizeInput(session,"Pathwaygenechoice", choices = dds.res$Gene, server = TRUE)
  } #end server
# Run the application 
shinyApp(ui = ui, server = server)

