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
ens2gene_hs <- 
t2g_hs <- read.table(file = config$t2g_hs_file, sep = "\t", header = T)
ens2gene_HS <- t2g_hs[,c(2,3)]
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
              goi_UI("GOI1")
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
    tabPanel("GSEA",  
             GSEA_UI("GSEA1")
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
    
  ## GOI tab####  
   goi_Server("GOI1", vst.goi)
    
  ##DESEq #####
   DE_Server("DEtab1", dds, vsd)
 
   
    #output DE table with adjustment for singscore reactive to toggle switch 
    # output$DETable <-
    #   renderDataTable({
    #     if(input$singscorebutton == TRUE) {
    #       CD_DE_DT_sing()
    #     } else if(input$singscorebutton == FALSE) {
    #       CD_DE_DT()
    #     }
    #     
    #   })
  
   
####GSEA output ####
    
  GSEA_Server("GSEA1", dds, ens2gene_HS)
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
