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
library(plotly)
library(BiocParallel)
library(ComplexHeatmap)
library(InteractiveComplexHeatmap)
library(circlize)
library(colourpicker)
library(ggsci)
library(scales)
library(esquisse)
library(ggprism)
library(shinycssloaders)
library(patchwork)
#Data ####
#load in data and metadata
meta_lut_ven <- readRDS("data/meta_lut_ven.Rds")
qcdt<-load_multiqc("data/multiqc_data.json", sections="raw") 
vst.goi <- readRDS("data/vst.goi.rds")
#DESeq data table
dds.res <- readRDS("data/DEtable.rds")
#sample metadata table
metadata <- read.table(file = "data/SampleSheetJordanLab.txt")
# DE table with singscore
dds.resscore <- readRDS("data/dds.resscore.rds")
#tables for PCA
vsd.pca <- readRDS("data/vsd.pca.rds")
bcvsd.pca <- readRDS("data/bcvsd.pca.rds")
#tables for variance
vsd.variance <- readRDS("data/vsd.variance.rds")
bcvsd.variance <- readRDS("data/bcvsd.variance.rds")

vsd2.pca <- readRDS("data/vsd2.pca.rds")
bcvsd2.pca <- readRDS("data/bcvsd2.pca.rds")
#GSEA data table
ranks <- readRDS("data/ranks.rds")
#load molecular pathways for GSEA
pathways.hallmark <- gmtPathways("data/gmt_pathway_files copy/h.all.v7.4.symbols.gmt")
pathways.GOall <- gmtPathways("data/gmt_pathway_files copy/c5.go.v2022.1.Hs.symbols.gmt")
pathways.GOmolec <- gmtPathways("data/gmt_pathway_files copy/c5.go.mf.v7.4.symbols.gmt")
pathways.GOcellcomp <- gmtPathways("data/gmt_pathway_files copy/c5.go.cc.v2022.1.Hs.symbols.gmt") 
pathways.GObio <- gmtPathways("data/gmt_pathway_files copy/c5.go.bp.v7.4.symbols.gmt")
pathways.TFtargets <-gmtPathways("data/gmt_pathway_files copy/c3.tft.v2022.1.Hs.symbols.gmt")
pathways.allReg <- gmtPathways("data/gmt_pathway_files copy/c3.all.v2022.1.Hs.symbols.gmt")
pathways.Wiki <- gmtPathways("data/gmt_pathway_files copy/c2.cp.wikipathways.v2022.1.Hs.symbols.gmt")
pathways.Reactome <-gmtPathways("data/gmt_pathway_files copy/c2.cp.reactome.v2022.1.Hs.symbols.gmt")
pathways.KEGG <- gmtPathways("data/gmt_pathway_files copy/c2.cp.kegg.v2022.1.Hs.symbols.gmt")
pathways.Positional <-gmtPathways("data/gmt_pathway_files copy/c1.all.v2022.1.Hs.symbols.gmt")
pathways.Biocarta <-gmtPathways("data/gmt_pathway_files copy/c2.cp.biocarta.v2022.1.Hs.symbols.gmt")
pathways.lsc <- gmtPathways("data/gmt_pathway_files copy/lsc_sigs.gmt")
pathways.aeg <- gmtPathways("data/gmt_pathway_files copy/aeg_genesets_20220602.gmt")
vstlimma <- readRDS("data/vstlimma.rds")

names(pathways.aeg)[10] <- "PM_Primitive_Blast"
names(pathways.aeg)[9] <- "PM_Monocytic_Blast"
#pathways.aegGOBP <- c(pathways.aeg, pathways.GObio)
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
                  sidebarPanel( #create toggle switches for each plot
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
                    materialSwitch(
                      inputId =
                        "multiqc",
                      label =
                        "MultiQC",
                      value =
                        FALSE,
                      right =
                        TRUE
                    ),  
                    #palette choices for PCA plots
                     paletteUI("palette"),
                   
                    conditionalPanel(
                    condition = "input.PCAplots == true",
                    radioButtons( #choose type of PCA plot
                      "PCAvar",
                      h4(
                        "Choose PCA plot"
                      ),
                      choices =
                        list("VST PCA", "VST + batch corrected PCA"),
                      selected =
                        "VST PCA"
                    ), 
                    downloadButton("downloadPlotPCA", label = "Download PCA Plot")
                    ),
                    
                    hr(),
                    
                    conditionalPanel(
                      condition = "input.PCAscreeplots == true",
                      
                      downloadButton("downloadPlotscree",
                                     label =
                                       "Download Scree Plot")
                    ),
                    hr(),
                    conditionalPanel(
                      condition = "input.multiqc == true",
                      selectInput( #choose type of multiqc test to visualize
                                "QCvar",
                                label=
                                  "Choose MultiQC test",
                                choices =
                                  c("% mapped reads", "# mapped reads", "% uniquely mapped reads", "# uniquely mapped reads"),
                                selected =
                                  "% mapped reads"
                              ) #end selectInput
                    )
                  ),
                  mainPanel(
                    conditionalPanel(
                      condition = "input.PCAplots == true",
                                     plotOutput("PCAplot")),
                    conditionalPanel(
                      condition = "input.PCAscreeplots == true",
                                     plotOutput("PCAvarplot")),
                    conditionalPanel(
                      condition = "input.multiqc == true",
                                     plotOutput("QCplot"))
                  )
                )
              )
    ), 
    
    tabPanel( #Gene expression analysis ####
              "Gene Expression",
              fluidPage(
                theme=
                  shinytheme("flatly"),
                titlePanel(
                  "Gene Expression Plot"
                ),
                h6("*p values are indicated for the comparison of gene expression between prim and mono samples"),
                sidebarLayout(
                  sidebarPanel(
                    useShinyjs(), #this is needed for javascript functions
                    selectizeInput( #gene choice dropdown menu
                      "VSTCDgenechoice",
                      label=
                        "Choose a gene for analysis",
                      choices =
                        NULL,
                      selected = NULL,
                      options = list(maxItems = NULL)
                    ), #options for axis variables, fill variable, and plot filters
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
                    #add a facet toggle switch
                    materialSwitch("genefacetbutton", label = "Facet", value = FALSE, right = TRUE),
                    hr(),
                   
                    #add palette choices for boxplot colors
                    paletteUI("palette2"),
                   
                    hr(), #js functions to hide plot dimensions until selected
                    materialSwitch("hidedims", "Custom plot dimensions", value = FALSE, right = TRUE),
            
                    #plot dimension input
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
                  
                  mainPanel( #add loading spinner
                    shinycssloaders::withSpinner(
                      plotOutput(
                      "VSTCDplot"
                    )
                    )
                  )
                )
              )
    ),
    
    tabPanel("Differential Expression",# DESeq Menu ####
             fluidPage(
               theme =
                 shinytheme("flatly"),
               titlePanel(
                 "Differential Expression Tables and Plots"
               ),#end title
               sidebarLayout(
                 sidebarPanel( #toggle buttons for choosing each plot type
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
                    #option to filter table by padj
                   radioButtons("padjbutton", label = "Filter DE tables by padj", 
                                choices = list("<= 0.01" = "sigvar1", "<= 0.05" = "sigvar5", "All" = "allvar"), selected = "allvar"),
                   hr(),
                   downloadButton("downloadDEtable", label = "Download DE Table")
                   ),
                   
                   hr(),
 
                   conditionalPanel(
                     condition = "input.DESeqvolcano == true",
                     h4("Volcano Plot Specific Options"),
                     #color palette choice for volcano plot
                     colorUI("color"),
                     colorUI("color2"),
                     
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
                     #color palette choice for MA plot
                     colorUI("color3"),
                     colorUI("color4"),
                    
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
                     #color palette choices for heatmap
                     colorUI("color5"),
                     colorUI("color6"),
                     # colourInput(
                     #   "heatcolor1",
                     #   label = "Choose 1st color",
                     #   value = "red",
                     #   showColour = ("both"),
                     #   palette = ("square"),
                     #   allowedCols = NULL,
                     #   allowTransparent = FALSE,
                     #   returnName = FALSE,
                     #   closeOnClick = FALSE
                     # ),
                     # colourInput(
                     #   "heatcolor2",
                     #   label = "Choose 2nd color",
                     #   value = "blue",
                     #   showColour = ("both"),
                     #   palette = ("square"),
                     #   allowedCols = NULL,
                     #   allowTransparent = FALSE,
                     #   returnName = FALSE,
                     #   closeOnClick = FALSE
                     # )

                       )
                     ),
                 
                 mainPanel(
                   conditionalPanel( 
                     condition = "input.DESeqtable == true",
                     shinycssloaders::withSpinner( #add loading spinner
                       DTOutput(
                       "DETable"
                     )
                     )
                   ),
                   conditionalPanel(
                     condition = "input.DESeqvolcano == true",
                     shinycssloaders::withSpinner( #add loading spinner
                       plotlyOutput(
                       "DEVolcanoPlot"
                     )
                     )
                   ),
                   conditionalPanel(
                     condition = "input.DESeqMA == true",
                     shinycssloaders::withSpinner( #add loading spinner
                       plotlyOutput(
                       "DEMAPlot"
                     )
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
                   #dropdown menu for molecular pathways 
                   selectInput("filechoice", label = NULL,
                               choices = c(Hallmark = "hallmark", GOall = "goall",GOmolecular = "GOmolec", 
                                           GOcellcomp = "GOcellcomp", GObio = "GObio", TFtargets = "TFtargets",
                                           allRegular = "allReg", Wiki = "wiki", Reactome = "reactome", KEGG = "KEGG",
                                           Positional = "positional", Biocarta = "biocarta", lsc = "lsc", aeg = "aeg")),
                   h4("Choose a plot type"),
                   #toggle buttons to load table and plots
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
                     #color palette choices for waterfall plot
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
                     #slider scale to choose how many pathways to load
                     sliderInput("howmanypathways", "Choose How Many Pathways to Rank",
                                 min = 5, max = 50, value = 15
                     ),
                     hr(),
                     #js function to hide plot dimensions until selected
                   materialSwitch("hidedimsWF", "Custom plot dimensions", value = FALSE, right = TRUE),
                   
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
                     
                     hr(),
                   #color palette choices for muostache plot
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
                     #options for type of enrichment plot to load
                     radioButtons("topupordownbutton", h5("Enrichment Plot Choices"), 
                                  choices = list("Top Ranked Up Pathway" = "topup", "Top Ranked Down Pathway" = "topdown", "Pathway of Choice:" = "eplotpath"), selected = "topup"),
                     h5("Choose a specific pathway"),
                     #dropdown menu for specific pathway choices that is reactive to the pathway set dropdown menu
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
                     #dropdown list of specific pathways for volcano plot that is reactive to the pathway set dropdown menu
                     selectizeInput(
                       "pathwaylist",
                       label=
                         "Choose a specific pathway(s) to view on volcano plot",
                       choices =
                         NULL,
                       selected = NULL,
                       options = list(maxItems = 1)
                     ),
                #color palette choices for volcano plot
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
                   #dropdown menu of specific pathways for heatmap, reactive to pathway sets dropdown
                     selectizeInput(
                       "pathwaylistht",
                       label=
                         "Choose a specific pathway to view genes on heatmap",
                       choices =
                         NULL,
                       selected = NULL ,
                       options = list(maxItems = 1)
                     ),
                   #color palette choices for heatmap
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
                     )
                     
                   ),
                 ),
                 
                 
                 mainPanel(
          
                   conditionalPanel(
                     condition = "input.fgseaTable == true",
                        DTOutput( #add loading spinners
                       "fgseaTable"
                        )
                   ),
                   conditionalPanel(
                     condition = "input.rankedplot == true",
                     shinycssloaders::withSpinner( #add loading spinners
                       plotOutput(
                       "GSEAranked"
                     )
                     )
                   ),
                   conditionalPanel(
                     condition = "input.moustache == true",
                     shinycssloaders::withSpinner( #add loading spinners
                       plotOutput(
                       "GSEAMoustache"
                     )
                     )
                   ),
                   conditionalPanel(
                     condition = "input.eplot == true",
                     shinycssloaders::withSpinner( #add loading spinners
                       plotOutput(
                       "GSEAenrichment"
                     )
                     )
                   ),
                   conditionalPanel(
                     condition = "input.volcanoplot == true",
                     shinycssloaders::withSpinner( #add loading spinners
                       plotOutput(
                       "GSEAvolcano"
                     )
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
    tabPanel("GSEA Pathway/Gene Visualization", #GOI plots ####
             fluidPage(
               theme =
                 shinytheme("flatly"),
               titlePanel(
                 "GSEA: Interrogation of pathways containing a gene of interest"
               ),#end title
               h4("Positive NES is upregulated in Primitive cells and negative NES is upregulated in Monocytic cells"),
               sidebarLayout( 
                 sidebarPanel( 
                   useShinyjs(),
                   #gene list dropdown menu of all genes from DE table
                   selectizeInput(
                     "Pathwaygenechoice",
                     label=
                       "Choose a gene of interest",
                     choices =
                       NULL,
                     selected = NULL,
                     options = list(maxItems = 1) #only one gene can be selected at once
                   ),
                   hr(),
                   #pathway set dropdown list
                   selectInput("genefilechoice", "Choose gmt file to plot pathways containing the gene of interest",
                               choices = c(Hallmark = "hallmark", GOall = "GOall", GOmolecular = "GOmolec", 
                                           GOcellcomp = "GOcellcomp", GObio = "GObio", TFtargets = "TFtargets",
                                           allRegular = "allReg", Wiki = "wiki", Reactome = "reactome", KEGG = "KEGG",
                                           Positional = "positional", Biocarta = "biocarta", lsc = "lsc", aeg = "aeg")),
                   
                   hr(),
                   #color palette options for GOI
                   palettePicker(
                     inputId = "PaletteChoicesGP",
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
                   #js function to hide plot dimension options until selected
                   materialSwitch("hidedimsGP", "Custom plot dimensions", value = FALSE, right = TRUE),
                   
                   sliderInput("goiheightslider", "Adjust plot height",
                               min = 200, max = 1000, value = 400
                   ),
                   hr(),
                   
                   sliderInput("goiwidthslider", "Adjust plot width",
                               min = 200, max = 1000, value = 600
                   ),
                   downloadButton(
                     "downloadGOI",
                     label =
                       "Download Plot"
                   )
                 ),
                 mainPanel(
                   shinycssloaders::withSpinner( #add loading spinners
                     plotOutput(
                     "PathwaysGenePlot"
                   )
                   )
                 )
               )
             )
             
    ),
    
 
    
  )
#Server ####
server <- 
  function(input, output, session) {
    # check to make sure server is initiating
    print("Initializing renderPlots")
    
    options(shiny.reactlog = TRUE)
  
   
  ###QC-MultiQC plots####
    #reactive function for multiqc plot title
    QC_title <- 
      reactive({
        if (input$QCvar == "% mapped reads") {
          print("% mapped reads per sample")
        } else if (input$QCvar == "# mapped reads") {
          print("# mapped reads per sample")
        } else if (input$QCvar == "% uniquely mapped reads") {
          print("% uniquely mapped reads per sample")
        } else if (input$QCvar == "# uniquely mapped reads") {
          print("# uniquely mapped reads per sample")
        }
      })
    #multiQC plot
    output$QCplot <- renderPlot({
      #create object for reactive data input based on user choice of multiqc test option
      QCdata <- switch(
        input$QCvar,
        "% mapped reads" = qcdt$raw.salmon.percent_mapped,
        "# mapped reads" = qcdt$raw.salmon.num_mapped,
        "% uniquely mapped reads" = qcdt$raw.star.uniquely_mapped_percent,
        "# uniquely mapped reads" = qcdt$raw.star.uniquely_mapped,
        "Sample_ID" = qcdt$metadata.sample_id
      )
     #need to make an object to tell the ggplot what the sample names are
      Sample_ID <-qcdt$metadata.sample_id

      ggplot(
        qcdt,
        aes(
          x = Sample_ID,
          y = QCdata
        )) +
        geom_point() +
        theme_cowplot (font_size = 18) +
        ggtitle(QC_title()) +
        theme(axis.title = element_text(face = "bold"), title = element_text(face = "bold"), axis.text.x =
                element_text(angle = 60, hjust = 1)) 
    }) #end renderPlot

    # PCA plots ####
#function to tell ggplot which data set to use depending on user input
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
    #write functions and store in object to calculate the % variance for each PC
    PC_var_VST <- data.frame(PC =paste0("PC", 1:12),variance =(((vsd2.pca$sdev) ^ 2 / sum((vsd2.pca$sdev) ^ 2)) * 100))
    lorder_VST <- as.vector(outer(c("PC"), 1:12, paste, sep = ""))
    PC_var_VST$PC <-factor(PC_var_VST$PC,levels = lorder_VST)
    
    #batch corrected PCA variance
    PC_var_bc <-data.frame(PC =paste0("PC", 1:12),variance =(((bcvsd2.pca$sdev) ^ 2 / sum((bcvsd2.pca$sdev) ^ 2)) * 100))
    lorder_bc <-as.vector(outer(c("PC"), 1:12, paste, sep = ""))
    PC_var_bc$PC <-factor(PC_var_bc$PC,levels = lorder_bc)
    #function to tell ggplot which data set to use for the scree plots
    PC_var_data <-
      eventReactive(input$PCAvar, {
        if (input$PCAvar == "VST PCA") {
          PC_var_VST
        } else if (input$PCAvar == "VST + batch corrected PCA") {
          PC_var_bc
        }
      })
    # reactive function for PCA plot title
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
    
    colorpaletteQC <- 
      paletteServer("palette")
  
    #PCA plot output
    output$PCAplot <- renderPlot ({
      ggplot(PCAdata(), aes(x = PC1, y = PC2, shape = condition, color = batch, fill = batch)) + 
        geom_point(size = 5) + 
        scale_shape_manual(values = c(21, 24), name = '') +
        scale_fill_viridis_d(option = colorpaletteQC()) + #scale_fill_manual reactive function
        scale_color_viridis_d(option = colorpaletteQC()) + #scale_color manual reactive function
        theme_cowplot(font_size = 18) + 
        theme(axis.title = element_text(face = "bold"), title = element_text(face = "bold")) +
        theme(plot.background = element_rect(fill = "#FFFFFF", colour = "#FFFFFF")) +
        theme(panel.background = element_rect(fill = "#FFFFFF", colour = "#FFFFFF")) +
        xlab(variance_PC1()) +  #reactive xlab
        ylab(variance_PC2()) + #reactive y lab
        ggtitle(PCA_title()) + #reactive title
        geom_text_repel(colour = "black", aes(label=sample_name),hjust=0, vjust=0)
      
    })
    
    #PCA plots download ####
    output$downloadPlotPCA <- downloadHandler(
      filename = function() { paste(input$PCAvar, '.png', sep='') },
      content = function(file) {
        ggsave(file, device = "png", width = 8, height = 6, units = "in",dpi = 72)
      }
    )
    #reactive function for scree plots title
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
        theme(axis.title = element_text(face = "bold"), title = element_text(face = "bold")) +
        labs(x = "PC",
             y = "% Variance") +
        labs(title =
               PCA_var_title()) #reactive title
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
    #reactive function for for filtering vst data table based on user input 
    datavst <-
      reactive({
        if(input$PrimMonobutton == "comparison") {
          vst.goi %>% 
            dplyr::filter(ext_gene %in% input$VSTCDgenechoice)
        } else if(input$PrimMonobutton == "prim") {
          vst.goi %>%
            dplyr::filter(ext_gene %in% input$VSTCDgenechoice) %>%
            dplyr::filter(class == "prim")
        } else if(input$PrimMonobutton == "mono") {
          vst.goi %>%
            dplyr::filter(ext_gene %in% input$VSTCDgenechoice) %>%
            dplyr::filter(class == "mono")
        }
      })
    
    
    #make sure duplicate selections are not allowed with radio buttons
    observeEvent(input$XaxisVar_CDgene, {
      if(input$XaxisVar_CDgene == "xvalue") {
        mychoices <- c("Gene" = "ygene")
      } else if(input$XaxisVar_CDgene=="xgene") {
        mychoices <- c("Value" = "yvalue")
      }else if(input$XaxisVar_CDgene == "xclass") {
        mychoices <- c("Value" = "yvalue")
      }
      updateRadioButtons(session, "YaxisVar_CDgene", choices = mychoices)
    })
    
    
    #x axis reactive output based on radio buttons
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
    #y axis reactive output based on radio buttons
    yvar_CDgene <-
      eventReactive(input$YaxisVar_CDgene, {
        if (input$YaxisVar_CDgene == "yvalue") {
          "value"
        } else if (input$YaxisVar_CDgene == "ygene") {
          "ext_gene"
        }
      })

    
    #fill reactive output based on radio buttons
    fillvar_CDgene <-
      eventReactive(input$FillVar_CDgene, {
        if (input$FillVar_CDgene == "fillclass") {
          "class"
        } else if (input$FillVar_CDgene == "fillgene") {
          "ext_gene"
        }
      })
    # facet toggle switch function to turn faceting on or off
    Gene_facet <-
      eventReactive(input$genefacetbutton, {
        if(input$genefacetbutton == TRUE) {
          facet_grid(cols = vars(class))
        } else(NULL)
      })
  
    #call in palette module for plot
    colorpaletteGene <- 
      paletteServer("palette2")
    
    # function for adding padj values to plot, position needs to change when x and y variables change for readability
    sig_label_position <- reactive({
      value <- vst.goi$value
      if(input$XaxisVar_CDgene == "xvalue") {
        geom_text(aes(x = max(value), label = paste("p=",format(padj, digit = 1, scientific = T))),check_overlap = T) 
      } else if(input$XaxisVar_CDgene == "xgene") {
        geom_text(aes(y = max(value), label = paste("p=",format(padj, digit = 1, scientific = T))),check_overlap = T) 
      } else if(input$XaxisVar_CDgene == "xclass") {
        geom_text(aes(y = max(value), label = paste("p=",format(padj, digit = 1, scientific = T))),check_overlap = T) 
      }
    })
    #reactive wrapper for showing the plot dimensions options or hiding them based on toggle selection
    observe({
      toggle(id = "geneheightslider", condition = input$hidedims)
      toggle(id ="genewidthslider", condition = input$hidedims)
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
            Gene_facet() + #reactive faceting
            scale_fill_viridis_d(option = colorpaletteGene()) + #reactive  scale_fill_manual from module
            scale_color_viridis_d(option = colorpaletteGene()) + #reactive scale_color_manual from module
            geom_point(alpha = 0.5,
                       position = position_jitterdodge(jitter.width = 0.2),
                       aes(color = class)) + 
            theme_light() +
            sig_label_position() + # function for adjusted pvalues position and format on plot
            ylab("") +
            xlab("") +
            ggtitle("Gene Expression: Prim vs Mono")
        }) #end render plot
    
    output$downloadGenePlot <- downloadHandler(
      filename = function() { paste('GeneCentricPlot','.png', sep='') },
      content = function(file) {
        ggsave(file, device = "png", width = 8,
               height = 8, dpi = 72)
      }
    )
    
    #DESEq #####
    
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
    #function for filtering DE object with monocytic contribution regressed out based on padj value chosen by user 
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

    #output DE table with adjustment for singscore reactive to toggle swich 
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
    colorDE <- 
      colorServer("color")
    
    color2DE <-
      colorServer("color2")
    
    output$DEVolcanoPlot <-
      renderPlotly({
        colors <- c(colorDE(), "grey", color2DE()) #object for colors on volcano based on user input
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
    #call for module_singlcolor_palette
    color3DE <- 
      colorServer("color3")
    
    color4DE <-
      colorServer("color4")
    
    output$DEMAPlot <- renderPlotly ({
        colors <- c(color3DE(), "grey", color4DE()) #object for color choices dependent on user input
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
    
  
    #DE Heatmap ####
    color5DE <- 
      colorServer("color5")
    
    color6DE <-
      colorServer("color6")
    #interactive heatmap needs to be wrapped in a reactive function to work
    observe({
      #filter DE object for only significantly differentially experessed genes
     dds.mat <- dds.res %>%
      dplyr::filter(padj < 0.05 & abs(`log2FoldChange(Prim/Mono)`) >= 2)
    #filter vst counts matrix by sig expressed genes
    vst.mat <- vstlimma %>%
      dplyr::filter(., ensembl_gene_id %in% dds.mat$ensembl_gene_id) %>%
      column_to_rownames(., var = "ensembl_gene_id") %>%
      dplyr::select(.,-ext_gene) %>%
      as.matrix()
    rownames(vst.mat) = dds.mat$Gene
    vst.mat <- t(scale(t(vst.mat)))
    #only show the first 100 genes for visualization in this example(can change)
    vst.mat <- head(vst.mat, n = 100)
    #create a colorRamp function based on user input in color palette choices
    colors = colorRamp2(c(-2, 0, 2), c(color5DE(), "white", color6DE()))
    #create heatmap object
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
    #object for pathway choice files to use with the goi pathway input
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
        fgseaRes <- fgsea::fgsea(pathways = pathwaygsea, stats = ranks, nproc = 10)
        fgseaResTidy <- fgseaRes %>%
          as_tibble() %>%
          dplyr::select(., -pval,-log2err, -ES) %>% 
          arrange(desc(NES))
        fgseaResTidy
      })
 
    output$downloadfgsea <- downloadHandler(
      filename = function() { paste("GSEATable", '.csv', sep='') },
      content = function(file) {
        write.csv(gseafile(),file)
      }
    )
    #filter Res table for chosen pathway to show in a waterfall plot
    gseafile_waterfall <-
      reactive({
        pathwaygsea <- gsea_file_values[[input$filechoice]]
        fgseaRes <-
          fgsea::fgsea(pathways = pathwaygsea,
                       stats = ranks,
                       nproc = 10)
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
   
    #### GSEA pathway ranks waterfall plot ####
    
    output$GSEAranked <- renderPlot(
      width = function()
        input$rankedwidthslider,
      height = function()
        input$rankedheightslider,
      {
         if (input$rankedplot == TRUE) {
           #color object reactive to user choice from palette
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
    
    #download button output- waterfall plot
    output$downloadranks <- downloadHandler(
      filename = function() { paste("Waterfall Plot", '.png', sep='') },
      content = function(file) {
        ggsave(file, device = "png", width = 8, height = 6, units = "in",dpi = 72)
      }
    )
    #reactive wrapper for js functions to hide or show plot dimensions based on toggle switch
    observe({
      toggle(id = "rankedheightslider", condition = input$hidedimsWF)
      toggle(id ="rankedwidthslider", condition = input$hidedimsWF)
    })
    #### GSEA moustache plot ####
    #reactive title for moustache plot
    gseamoustache_title <-
      eventReactive(input$filechoice, {
        print(input$filechoice)
      })
   
    #function for making moustache reactive to pathway choice
    #fgsea analysis run on pathway chosen
    #resulting table is put into tidy format and uncessary rows are removed for plotting
    #new column added to state significance 
    toplotMoustache <-
      reactive({
        pathwaygsea <- gsea_file_values[[input$filechoice]]
        fgseaRes <-
          fgsea::fgsea(pathways = pathwaygsea,
                       stats = ranks,
                       nproc = 10)
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
      {
        if (input$moustache == TRUE) {
          #color object reactive to user input from plalette choice
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
            #only label pathways that sig
            geom_text_repel(colour = "black", aes(label= ifelse(padj <0.05, as.character(pathway), ""), hjust=0,vjust=0)) +
            coord_cartesian(xlim = c(-3, 3), ylim = c(-0.1, 1)) 
          print(m)
        }
      }
    )
    #download button output= moustache plot 
    output$downloadmoustache <- downloadHandler(
      filename = function() { paste("Moustache Plot", '.png', sep='') },
      content = function(file) {
        ggsave(file, device = "png", width = 8, height = 6, units = "in",dpi = 72)
      }
    )
    #GSEA Enrichment Plots ####
    #reactive expression for specific pathway choice
    observe({
      pathwaygsea <- gsea_file_values[[input$filechoice]]
      updateSelectizeInput(session,"pathwaylisteplot", choices = names(pathwaygsea), server = TRUE)})
    #reactive expression for plot title based on enrichment plot type
    gseaeplot_title <-
      eventReactive(input$pathwaylisteplot, {
        print(input$pathwaylisteplot)
      })
    #fgsea run for selected pathway, tidy results table, filter for top up or down pathway to plot, or plot pathway of choice
    output$GSEAenrichment <- renderPlot ({
      pathwaygsea <- gsea_file_values[[input$filechoice]]
      fgseaRes <-
        fgsea::fgsea(pathways = pathwaygsea,
                     stats = ranks,
                     nproc = 10)
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
   #download button output- enrichment plot
    output$downloadeplot <- downloadHandler(
      filename = function() { paste("Enrichment Plot", '.png', sep='') },
      content = function(file) {
        ggsave(file, device = "png", width = 8, height = 6, units = "in",dpi = 72)
      }
    )
    # GSEA Volcano plot ####
    #reactive expression for selection of specific pathway by user input
    observe({
      pathwaygsea <- gsea_file_values[[input$filechoice]]
      updateSelectizeInput(session,"pathwaylist", choices = names(pathwaygsea), server = TRUE)})
    #function to select only genes from the DE object that are found in the chocen pathway
    dds.res.pathways <- reactive({
      pathwaygsea <- gsea_file_values[[input$filechoice]]
      p <-
        unlist((pathwaygsea[names(pathwaygsea) %in% input$pathwaylist]))

      dds.res.pathways <- dds.res %>%
        mutate(., genes_in_pathway = ifelse(Gene %in% p, 'yes', 'no'))
 print(dds.res.pathways)
    })
    #reactive title for volcano based on specific pathway choice
    gseavol_title <-
      eventReactive(input$pathwaylist, {
        paste(input$pathwaylist)
      })
   
    output$GSEAvolcano <- renderPlot ({
      #color object reactive to user input from palette chpice
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
          aes( #only label is gene is in pathway and sig expression 
            label = ifelse(
              genes_in_pathway == 'yes' & `log2FoldChange(Prim/Mono)` > 1.5,
              as.character(Gene),
              ""
            )
          ),
          hjust = 0,
          vjust = 0
        ) +
        ggtitle(gseavol_title()) + #reactive title
          xlab("log2foldchange")
      }
    })
   #download button output- gsea volcano plot
    output$downloadvolcano <- downloadHandler(
      filename = function() { paste("Volcano Plot", '.png', sep='') },
      content = function(file) {
        ggsave(file, device = "png", width = 8, height = 6, units = "in",dpi = 72)
      }
    )
    #GSEA heatmap ####
    #reactive expression for selected pathway choice for heatmap
    observe({
      pathwaygsea <- gsea_file_values[[input$filechoice]]
      updateSelectizeInput(session,"pathwaylistht", choices = names(pathwaygsea), server = TRUE)})
    #reactive expression of title for heatmap
    gseaht_title <-
      eventReactive(input$pathwaylistht, {
       print(input$pathwaylistht)
      })
    #interactive heatmap must be wrapped in a reactive expression
    observeEvent(input$pathwaylistht, {
      if(input$heatmap == TRUE) {
      pathwaygsea <- gsea_file_values[[input$filechoice]]
      
      p <- unlist((pathwaygsea[names(pathwaygsea) %in% input$pathwaylistht]))
      #filter DE object for significance 
      dds.sig <- dds.res %>%
        dplyr::filter(padj < 0.05 & abs(`log2FoldChange(Prim/Mono)`) >= 0.5)
      #filter vst counts matrix for genes in pathway
      vst.myc <- vstlimma %>% 
        mutate(., pathwayheat = ifelse(ext_gene %in% p, 'yes', 'no')) %>% 
        dplyr::filter(pathwayheat == "yes")
      #create matrix for heatmap using genes in significant DE that are in pathway of choice
      vstgsea.mat <- vst.myc %>%
        dplyr::filter(., ensembl_gene_id %in% dds.sig$ensembl_gene_id) %>%
        column_to_rownames(., var = "ext_gene") %>%
        dplyr::select(.,-ensembl_gene_id, -pathwayheat) %>%
        as.matrix()
      #transform and scale and transform back
      vstgsea.mat <- t(scale(t(vstgsea.mat)))
      #color function buliding a colorRamp palette based on user input from palette choices
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
  
    #Gene Centeric pathways analysis plots ####
    #reactive title
    gsea_gene_title <- 
      eventReactive(input$Pathwaygenechoice, {
        paste(input$Pathwaygenechoice)
      })
    #reactive color palette
    colorchoicesGP <-
      eventReactive(input$PaletteChoicesGP, {
        if(input$PaletteChoicesGP == "viridis") {
          scale_color_viridis_d(option = "viridis")
        } else if(input$PaletteChoicesGP == "cividis") {
          scale_color_viridis_d(option = "cividis")
        } else if(input$PaletteChoicesGP == "magma") {
          scale_color_viridis_d(option = "magma")
        } else if(input$PaletteChoicesGP == "plasma") {
          scale_color_viridis_d(option = "plasma")
        }else if(input$PaletteChoicesGP == "inferno") {
          scale_color_viridis_d(option = "inferno")
        }
      })
    #reactive wrapper for js function to hide or show plot dimension options
    observe({
      toggle(id = "goiheightslider", condition = input$hidedimsGP)
      toggle(id ="goiwidthslider", condition = input$hidedimsGP)
    })
    
    genecentricgseaplot <- reactive({
      #load chosen pathway file based on reactive input 
      genepathwaygsea <- (gene_gsea_file_values[[input$genefilechoice]])
      #load fgsea table data for chosen pathway
      fgseaRes <-
        fgsea::fgsea(pathways = genepathwaygsea,
                     stats = ranks,
                     nproc = 10)
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
      goi_paths$GOI <- "pathways with GOI"
      #filter gsea table to find pathways that do not include the GOI in the leading edge
      nongoi_paths <- fgseaResTidy %>%
        dplyr::filter(!grepl(input$Pathwaygenechoice, leadingEdge)) %>%
        mutate(., class = ifelse(NES <0, 'Mono', 'Prim'))

      #put no for pathways that do not contain the goi
      nongoi_paths$GOI <- "pathways NOT with GOI"
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
          colorchoicesGP() +
          facet_wrap( ~ GOI, scales = "free") +
          theme_light(base_size = 18) +
          theme(axis.title = element_text(face = "bold"), title = element_text(face = "bold")) +
          geom_hline(yintercept = 0, linetype = "dashed") +
          plot_annotation(
            title = "Pathways with and without Gene of Interest",
            theme =
              theme(
                plot.title =
                  element_text(
                    face = "bold",
                    hjust = 0.5,
                    size = 16
                      )
                  )
              )
      })
    output$downloadGOI <- downloadHandler(
      filename = function() { paste("Gene of Interest Plot", '.png', sep='') },
      content = function(file) {
        ggsave(file, device = "png", width = 8, height = 6, units = "in",dpi = 72)
      }
    )
    # #gene list for gene centric pathway analysis
    updateSelectizeInput(session,"Pathwaygenechoice", choices = dds.res$Gene, server = TRUE)
  } #end server
# Run the application 
shinyApp(ui = ui, server = server)
