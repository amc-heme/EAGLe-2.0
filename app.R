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
#read in data from config file
 

#sample_id <- config$sample_id
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


# UI ####
ui <-
  navbarPage(
    "EAGLe: Cancer Discovery",
    #Dataset tab ####
    tabPanel(
      "Dataset",
      data_UI("data1")
    ),
    #QC Menu ####
    tabPanel(
              "QC",
              QC_UI("QC1")
    ),
    # 
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
  ##Data tab ####
    config <- data_Server("data1")
    
    observeEvent(config,{
      t2g <- config$t2g
      metadata <- config$metadata
      batch <- config$batch
      #samples <- config$samples
      num_PCs <- config$num_PCs
      dds <- config$dds
      dds.res <- config$dds_res
      vsd <- config$vsd
      vsd.pca <- config$vsd.pca
      vst <- config$vst
      vst.goi <- config$vst.goi
      qc <- config$qc
      var_1 <- config$var_1
      var_2 <- config$var_2
    })
    # 
  ## QC tab ####
    QC_Server("QC1", vsd, vsd.pca, metadata, batch, var_1, var_2) #ready for new data
    
  ## GOI tab####  
    goi_Server("GOI1", vst.goi)
    
  ##DESEq #####
    DE_Server("DEtab1", dds.res, vst) #ready for new data
   
  ##GSEA output ####
    GSEA_Server("GSEA1", dds, t2g, dds.res, vst) #ready for new data
   
    
  ##GOI pathway output ####
   pathway_Server("pathway1", dds, t2g, dds.res)
    
  } #end server

# Run the application 
shinyApp(ui = ui, server = server)
