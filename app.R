#load libraries
library(shinythemes)
library(thematic)
library(shiny)
library(kableExtra)
library(viridisLite)
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
options(
  shiny.fullstacktrace = TRUE
)
#Data ####
#load in data and metadata
# load analysis functions 
vst.goi <- readRDS("data/vst.goi.rds")
#read in data from config file
source("~/Documents/GitHub/EAGLe-2.0/config.R")
base_dir <- config$base_dir
t2g_hs <- read.table(file = config$t2g_hs_file, sep = "\t", header = T)
ens2gene_HS <- t2g_hs[,c(2,3)]
metadata <- read_rds(config$metadata_file)
batch <- config$batch
sample_id <- config$sample_id
# samples <- data.frame(SRR = metadata$SRR, batch = metadata$batch, condition = metadata$sample_type, sample_name = metadata$DESeq_Sample_Name)
# samples <- samples[order(samples$SRR),]
# design <- config$design_formula
# salm_dirs <- sapply(sample_id, function(id) file.path(base_dir, id, 'quant.sf'))
# tx2gene <- t2g_hs[,c(1,2)]
# colnames(tx2gene) <- c('TXNAME', 'GENEID')
# txi <- tximport(salm_dirs, type = 'salmon', tx2gene = tx2gene, ignoreTxVersion = TRUE)
# ddsTxi <- DESeqDataSetFromTximport(txi, colData = samples, design = design)
# ddsTxi.filt <- ddsTxi[rowMins(counts(ddsTxi)) > 5, ]
# dds <- DESeq(ddsTxi.filt)
# dds.res <- data.frame(results(dds)) %>%
#   rownames_to_column(., var = 'ensembl_gene_id') %>%
#   dplyr::select(., ensembl_gene_id, baseMean, log2FoldChange, padj) %>%
#   left_join(unique(dplyr::select(t2g_hs, c(ensembl_gene_id, ext_gene))), ., by = 'ensembl_gene_id') %>%
#   dplyr::rename(., Gene = ext_gene) %>%
#   mutate(., DiffExp = ifelse(padj < 0.05 & log2FoldChange >= 0.5, 'up',
#                              ifelse(padj < 0.05 & log2FoldChange <= -0.5, 'down', 'no'))) %>%
#   na.omit(.)
dds <- read_rds(config$dds_file)
dds.res <- read_rds(config$dds_res_file)
vsd <- read_rds(config$vsd_file)
vsd.pca <- read_rds(config$vsd.pca_file)
vst <- read_rds(config$vst_file)
# vsd <- vst(ddsTxi, blind = F)
qc<-load_multiqc(config$qc, sections="raw") 
var_1 <- config$var_1
var_2 <- config$var_2

#load molecular pathways for GSEA
# pathways.hallmark <- gmtPathways("data/gmt_pathway_files copy/h.all.v7.4.symbols.gmt")
# pathways.GOall <- gmtPathways("data/gmt_pathway_files copy/c5.go.v2022.1.Hs.symbols.gmt")
# pathways.GOmolec <- gmtPathways("data/gmt_pathway_files copy/c5.go.mf.v7.4.symbols.gmt")
# pathways.GOcellcomp <- gmtPathways("data/gmt_pathway_files copy/c5.go.cc.v2022.1.Hs.symbols.gmt") 
# pathways.GObio <- gmtPathways("data/gmt_pathway_files copy/c5.go.bp.v7.4.symbols.gmt")
# pathways.TFtargets <-gmtPathways("data/gmt_pathway_files copy/c3.tft.v2022.1.Hs.symbols.gmt")
# pathways.allReg <- gmtPathways("data/gmt_pathway_files copy/c3.all.v2022.1.Hs.symbols.gmt")
# pathways.Wiki <- gmtPathways("data/gmt_pathway_files copy/c2.cp.wikipathways.v2022.1.Hs.symbols.gmt")
# pathways.Reactome <-gmtPathways("data/gmt_pathway_files copy/c2.cp.reactome.v2022.1.Hs.symbols.gmt")
# pathways.KEGG <- gmtPathways("data/gmt_pathway_files copy/c2.cp.kegg.v2022.1.Hs.symbols.gmt")
# pathways.Positional <-gmtPathways("data/gmt_pathway_files copy/c1.all.v2022.1.Hs.symbols.gmt")
# pathways.Biocarta <-gmtPathways("data/gmt_pathway_files copy/c2.cp.biocarta.v2022.1.Hs.symbols.gmt")
# pathways.lsc <- gmtPathways("data/gmt_pathway_files copy/lsc_sigs.gmt")
# pathways.aeg <- gmtPathways("data/gmt_pathway_files copy/aeg_genesets_20220602.gmt")
#vstlimma <- readRDS("data/vstlimma.rds")

# names(pathways.aeg)[10] <- "PM_Primitive_Blast"
# names(pathways.aeg)[9] <- "PM_Monocytic_Blast"
#pathways.aegGOBP <- c(pathways.aeg, pathways.GObio)
# UI ####
ui <-
  navbarPage(
    "EAGLe: Cancer Discovery",
    #QC Menu ####
    tabPanel( 
              "QC",
              QC_UI("QC1")
    ), 
    
    #Gene expression analysis ####
    tabPanel( 
              "Gene Expression",
              goi_UI("GOI1")
    ),
    
    #DESeq Menu ####
      tabPanel("Differential Expression",
             DE_UI("DEtab1")
             ),
    
    #GSEA menu ####
    tabPanel("GSEA",  
             GSEA_UI("GSEA1")
    ),
    
    #GOI plots ####
    tabPanel("GSEA Pathway/Gene Visualization", 
     pathway_UI("pathway1")
  )
  )


#Server ####
server <- 
  function(input, output, session) {
    # check to make sure server is initiating
    print("Initializing renderPlots")
    
    options(shiny.reactlog = TRUE)
  ## QC tab ####
    QC_Server("QC1", vsd, vsd.pca, metadata, var_1, batch) #ready for new data
    
  ## GOI tab####  
    goi_Server("GOI1", vst.goi)
    
  ##DESEq #####
    DE_Server("DEtab1", dds.res, vst) #ready for new data
   
  ##GSEA output ####
    GSEA_Server("GSEA1", dds, ens2gene_HS, dds.res, vst)
   
    
  ##GOI pathway output ####
    pathway_Server("pathway1", dds, ens2gene_HS, dds.res)
    
  } #end server

# Run the application 
shinyApp(ui = ui, server = server)
