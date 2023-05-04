#load libraries
library(shinythemes)
library(thematic)
library(shiny)
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
#Data ####
#load in data and metadata
# load analysis functions 
#meta_lut_ven <- readRDS("data/meta_lut_ven.Rds")
vst.goi <- readRDS("data/vst.goi.rds")
#bcvsd.pca <- readRDS("data/bcvsd.pca.rds")
#read in data from config file
source("~/Documents/GitHub/EAGLe-2.0/config.R")
base_dir <- config$base_dir
t2g_hs <- read.table(file = config$t2g_hs_file, sep = "\t", header = T)
metadata <- read.table(file = config$metadata_file, header = TRUE, sep = "\t")
sample_id <- config$sample_id
samples <- data.frame(SRR = metadata$SRR, batch = metadata$batch, condition = metadata$sample_type, sample_name = metadata$DESeq_Sample_Name)
samples <- samples[order(samples$SRR),]
design <- config$design_formula
salm_dirs <- sapply(sample_id, function(id) file.path(base_dir, id, 'quant.sf'))
tx2gene <- t2g_hs[,c(1,2)]
colnames(tx2gene) <- c('TXNAME', 'GENEID')
txi <- tximport(salm_dirs, type = 'salmon', tx2gene = tx2gene, ignoreTxVersion = TRUE)
ddsTxi <- DESeqDataSetFromTximport(txi, colData = samples, design = design)
ddsTxi.filt <- ddsTxi[rowMins(counts(ddsTxi)) > 5, ]
dds <- DESeq(ddsTxi.filt)
dds.res <- data.frame(results(dds)) %>%
  rownames_to_column(., var = 'ensembl_gene_id') %>%
  dplyr::select(., ensembl_gene_id, baseMean, log2FoldChange, padj) %>%
  left_join(unique(dplyr::select(t2g_hs, c(ensembl_gene_id, ext_gene))), ., by = 'ensembl_gene_id') %>%
  dplyr::rename(., Gene = ext_gene) %>%
  mutate(., DiffExp = ifelse(padj < 0.05 & log2FoldChange >= 0.5, 'up',
                             ifelse(padj < 0.05 & log2FoldChange <= -0.5, 'down', 'no'))) %>%
  na.omit(.)
vsd <- vst(ddsTxi, blind = F)
qc<-load_multiqc("data/multiqc_data.json", sections="raw") 

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
              QC_UI("QC1")
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
                   #materialSwitch("hidedims", "Custom plot dimensions", value = FALSE, right = TRUE),
                    
                    
                    #plot dimension input
                    sliderUI("plotheightslider", 200, 1200, 600, "Adjust Plot Height"),
                  
                    #  hr(),
                    # 
                    sliderUI("plotwidthslider", 200, 1200, 800, "Adjust Plot Width"),
                
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
             DE_UI("DEtab1")
              
                   # materialSwitch(
                   #   inputId =
                   #     "singscorebutton",
                   #   label =
                   #     "DE Table without Monocytic Contribution",
                   #   value =
                   #     FALSE,
                   #   right =
                   #     TRUE
                   # ),
                   
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
                     colorUI("color7", "Choose color for plot", "#FF0000"),
                     
                     hr(), 
                     #slider scale to choose how many pathways to load
                     sliderInput("howmanypathways", "Choose How Many Pathways to Rank",
                                 min = 5, max = 50, value = 15
                     ),
                     # sliderUI("howmanypathways", 5, 50, 15, "Choose How Many Pathways to Rank"), 
                     
                     hr(),
                     #js function to hide plot dimensions until selected
                   materialSwitch("hidedimsWF", "Custom plot dimensions", value = FALSE, right = TRUE),
                   
                     # sliderInput("rankedheightslider", "Adjust plot height",
                     #             min = 200, max = 1000, value = 400
                     # ),
                   #call in mosule UI for height slider
                   sliderUI("rankedheightslider", 200, 1000, 400, "Adjust Plot Height"),
                   
                     hr(),
                     
                     # sliderInput("rankedwidthslider", "Adjust plot width",
                     #             min = 200, max = 1000, value = 600
                     # ),
                   #call in module UI for width slider
                   sliderUI("rankedwidthslider", 200, 1000, 600, "Adjust Plot Width"),
    
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
                    colorUI("color8", "Choose color for plot", "#FF0000"),
                   
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
                     colorUI("color9", "Choose color for plot", "#FF0000"),
                     
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
        
                   #color palette choices for heatmap
             colorUI("color10", "Choose 1st color", "#FF0000"),
             
             colorUI("color11", "Choose 2nd color", "#0000FF"),
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
                  paletteUI("paletteGOI"),
                   
                   hr(),
                   #js function to hide plot dimension options until selected
                   materialSwitch("hidedimsGP", "Custom plot dimensions", value = FALSE, right = TRUE),
                   
                   sliderUI("goiheightslider", 200, 1000, 400, "Adjust plot height"
                   ),
                   hr(),
                   
                   sliderUI("goiwidthslider", 200, 1000, 600, "Adjust plot width"
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
  ## QC tab ####
   QC_Server("QC1",colorpaletteQC)
    
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
    
    geneheight <-
    sliderServer("plotheightslider")
    
    genewidth <-
    sliderServer("plotwidthslider")
    
    
    #reactive wrapper for showing the plot dimensions options or hiding them based on toggle selection
    # observe({
    #   toggle(id = geneheight(), condition = input$hidedims)
    #   toggle(id = genewidth(), condition = input$hidedims)
    # })
    #plot output
    output$VSTCDplot <-
      renderPlot(
        width = function() genewidth(), #input$genewidthslider,
        height = function() geneheight(), #input$geneheightslider,
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
    DE_Server("DEtab1", dds, vsd)
 
    #function for filtering DE object with monocytic contribution regressed out based on padj value chosen by user 
  

    #output DE table with adjustment for singscore reactive to toggle swich 
    # output$DETable <-
    #   renderDataTable({
    #     if(input$singscorebutton == TRUE) {
    #       CD_DE_DT_sing()
    #     } else if(input$singscorebutton == FALSE) {
    #       CD_DE_DT()
    #     }
    #     
    #   })
    
    # download DE table
    # output$downloadDEtable <- downloadHandler(
    #   filename = function() { paste("DESeqTable", '.csv', sep='') },
    #   content = function(file) {
    #     write.csv(CD_DE_DT(),file)
    #   }
    # )
    # 

    # output$downloadDEVolcano <- downloadHandler(
    #   filename = function() { paste(input$sigvaluesbutton, '.png', sep='') },
    #   content = function(file) {
    #     ggsave(file, device = "png", width = 8, height = 6, units = "in",dpi = 72)
    #   }
    # )
    
    
    # output$downloadDEMA <- downloadHandler(
    #   filename = function() { paste('DESeqMAplot', '.png', sep='') },
    #   content = function(file) {
    #     ggsave(file, device = "png", width = 8, height = 6, units = "in",dpi = 72)
    #   }
    # )

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
    #call in singlecolor_palette module for plot
    colorWF <- 
      colorServer("color7")
    rankedheight <-
      sliderServer("rankedheightslider")
    rankedwidth <-
      sliderServer("rankedwidthslider")
    
    output$GSEAranked <- renderPlot(
      width = function()
        rankedwidth(),
      height = function()
        rankedheight(),
      {
         if (input$rankedplot == TRUE) {
           #color object reactive to user choice from palette
            colors <- c("grey", colorWF())
        
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
    #call in singlecolor_module for plot
    colorM <- 
      colorServer("color8")
    output$GSEAMoustache <- renderPlot(
      {
        if (input$moustache == TRUE) {
          #color object reactive to user input from plalette choice
          colors <- c('grey', colorM())
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
   #call in singlecolor_module for plot
    colorVol <- 
      colorServer("color9")
    output$GSEAvolcano <- renderPlot ({
      #color object reactive to user input from palette chpice
      colors <- 
        c("grey", colorVol())
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
    #call in singlecolor_module for plot palette
    colorGHeat <- 
      colorServer("color10")
    colorGHeat2 <-
      colorServer("color11")
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
      colors = colorRamp2(c(-2, 0, 2), c(colorGHeat(), "white", colorGHeat2()))
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
    #reactive color palette from module
    colorGOI <-
      paletteServer("paletteGOI")
    goiheights <- 
      sliderServer("goiheightslider")
    goiwidths <- 
      sliderServer("goiwidthslider")
    #reactive wrapper for js function to hide or show plot dimension options
    # observeEvent(goiheights(), {
    #   toggle(id = "goiheightslider", condition = input$hidedimsGP)
    #   toggle(id ="goiwidthslider", condition = input$hidedimsGP)
    # })
    # 
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
      width = function() goiwidths(),
      height = function() goiheights(),
      {
        ggplot(genecentricgseaplot(), aes(
          x = class,
          y = NES,
          color = (padj < 0.05)
        )) +
          geom_boxplot()  +
          scale_color_viridis_d(option = colorGOI()) +
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
