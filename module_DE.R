##Differential Expression ##
library(tximport)
source("~/Documents/GitHub/EAGLe-2.0/config.R")
base_dir <- config$base_dir
t2g_hs <- read.table(file = config$t2g_hs_file, sep = "\t", header = T)
metadata <- read.table(file = config$metadata_file, header = TRUE, sep = "\t")
sample_id <- config$sample_id
samples <- config$samples

DE_UI <- function(id) {
        tagList(
        materialSwitch(
          inputId =
            (NS(id,"DESeqtable")),
          label =
            "DE Table",
          value =
            FALSE,
          right =
            TRUE
        ),
        DTOutput(NS(id,"results"))
        )
}

DE_Server <- function(id) {
  moduleServer(id, function(input, output, session) {

  run_DE <- function() {
    salm_dirs <- sapply(sample_id, function(id) file.path(base_dir, id, 'quant.sf'))
    tx2gene <- t2g_hs[,c(1,2)]
    colnames(tx2gene) <- c('TXNAME', 'GENEID')
    txi <- tximport(salm_dirs, type = 'salmon', tx2gene = tx2gene, ignoreTxVersion = TRUE)
    ddsTxi <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ batch + condition)
    ddsTxi.filt <- ddsTxi[rowMins(counts(ddsTxi)) > 5, ]
    dds <- DESeq(ddsTxi.filt)
  
    dds.res<- data.frame(results(dds)) %>%
      rownames_to_column(., var = 'ensembl_gene_id') %>%
      dplyr::select(., ensembl_gene_id, baseMean, log2FoldChange, padj) %>%
      left_join(unique(dplyr::select(t2g_hs, c(ensembl_gene_id, ext_gene))), ., by = 'ensembl_gene_id') %>%
      dplyr::rename(., Gene = ext_gene) %>%
      mutate(., DiffExp = ifelse(padj < 0.05 & log2FoldChange >= 0.5, 'up', 
                                 ifelse(padj < 0.05 & log2FoldChange <= -0.5, 'down', 'no'))) %>%
      na.omit(.)
    dds.res
  }

  output$results <- renderDataTable({
    if(input$DESeqtable == TRUE) {
    run_DE()
    }
  })
  DEres_reactive <- reactive({
    run_DE()
  })
  list(dds.res = DEres_reactive)
  })
}

DE_App <- function() {
  ui <- fluidPage(
    DE_UI("DE1")
  )
  server <- function(input, output, session) {
    DE_Server("DE1")
  }
  shinyApp(ui, server)
}
DE_App()
