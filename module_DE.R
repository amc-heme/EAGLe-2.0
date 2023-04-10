##Differential Expression ##
library(tximport)
base_dir <- "~/Documents/CD files/salmon"
t2g_hs <- read.table(file = "~/Documents/CD files/HS_transcript_to_gene.txt", sep = "\t", header = T)
metadata <- read.table(file = "~/Documents/CD files/SampleSheet.txt", header = TRUE, sep = "\t")

DE_UI <- function(id) {
        tagList(
        DTOutput(NS(id,"results"))
        )
}

DE_Server <- function(id) {
  moduleServer(id, function(input, output, session) {

  run_DE <- function(base_dir, t2g_hs, metadata) {
    sample_id <- c("SRR9265370", "SRR9265373", "SRR9265371", "SRR9265372", "SRR9265363", "SRR9265364", "SRR9265366", "SRR9265367", "SRR9265369", "SRR9265365", "SRR9265368", "SRR9265374")
    salm_dirs <- sapply(sample_id, function(id) file.path(base_dir, id, 'quant.sf'))
    samples <- data.frame(SRR = metadata$SRR, batch = metadata$batch, condition = metadata$sample_type, sample_name = metadata$DESeq_Sample_Name)
    samples <- samples[order(samples$SRR),]
    tx2gene <- t2g_hs[,c(1,2)]
    colnames(tx2gene) <- c('TXNAME', 'GENEID')
    txi <- tximport(salm_dirs, type = 'salmon', tx2gene = tx2gene, ignoreTxVersion = TRUE)
    ddsTxi <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ batch + condition)
    ddsTxi.filt <- ddsTxi[rowMins(counts(ddsTxi)) > 5, ]
    dds <- DESeq(ddsTxi)
  
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
    run_DE(base_dir, tx2geneFile, metadata)
  })
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
