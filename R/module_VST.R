##VST module
#read in data from config file
# source("~/Documents/GitHub/EAGLe-2.0/config.R")
# base_dir <- config$base_dir
# samples <- config$samples
# sample_id <- config$sample_id
# t2g_hs <- read.table(file = config$t2g_hs_file, sep = "\t", header = T)
# tx2gene <- t2g_hs[,c(1,2)]
# colnames(tx2gene) <- c('TXNAME', 'GENEID')
# salm_dirs <- sapply(sample_id, function(id) file.path(base_dir, id, 'quant.sf'))

VSD_UI <- function(id) {
  ns <- NS(id)
  tagList(
    DTOutput(ns("vsd"))
  )
  
}

VSD_Server <- function(id) {
  moduleServer(id, function(input, output, session) {
      # 
      # txi <- tximport(salm_dirs, type = 'salmon', tx2gene = tx2gene, ignoreTxVersion = TRUE)
      # ddsTxi <- DESeqDataSetFromTximport(txi, colData = samples, design = ~batch + condition)
      # vsd <- vst(ddsTxi, blind = F)
      vst <- data.frame(assay(vsd)) %>%
        rownames_to_column(., var = "ensembl_gene_id") %>%
        left_join(unique(dplyr::select(t2g_hs, c(ensembl_gene_id, ext_gene))), ., by = 'ensembl_gene_id') %>%
        na.omit(.)


    output$vsd <- renderDT(
      vst
    )
      
  })
}

VSD_App <- function() {
  ui <- fluidPage(
    VSD_UI("VSD1")
  )
  server <- function(input, output, session) {
    VSD_Server("VSD1")
  }
  shinyApp(ui, server)
}
VSD_App()
