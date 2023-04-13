##Differential Expression tab##
library(tximport)
source("~/Documents/GitHub/EAGLe-2.0/config.R")
base_dir <- config$base_dir
t2g_hs <- read.table(file = config$t2g_hs_file, sep = "\t", header = T)
metadata <- read.table(file = config$metadata_file, header = TRUE, sep = "\t")
sample_id <- config$sample_id
samples <- config$samples

DE_UI <- function(id) {
  ns <- NS(id)
  fluidPage(
    theme =
      shinytheme("flatly"),
    titlePanel(
      "Differential Expression Tables and Plots"
    ),#end title
      sidebarLayout(
        sidebarPanel(
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
        hr(),
        
        materialSwitch(
          inputId =
            (ns("DESeqvolcano")),
          label =
            "Volcano Plot",
          value =
            FALSE,
          right =
            TRUE
        ),
        hr(),
        materialSwitch(
          inputId =
            (ns("DESeqMA")),
          label =
            "MA Plot",
          value =
            FALSE,
          right =
            TRUE
        )
        )
        ),
        mainPanel(
        DTOutput(ns("results")),
        girafeOutput(ns("volplot")),
        girafeOutput(ns("MAplot"))
        )
      )
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

  }
  
  
  dds.res<- data.frame(results(run_DE())) %>%
    rownames_to_column(., var = 'ensembl_gene_id') %>%
    dplyr::select(., ensembl_gene_id, baseMean, log2FoldChange, padj) %>%
    left_join(unique(dplyr::select(t2g_hs, c(ensembl_gene_id, ext_gene))), ., by = 'ensembl_gene_id') %>%
    dplyr::rename(., Gene = ext_gene) %>%
    mutate(., DiffExp = ifelse(padj < 0.05 & log2FoldChange >= 0.5, 'up', 
                               ifelse(padj < 0.05 & log2FoldChange <= -0.5, 'down', 'no'))) %>%
    na.omit(.)
  
  output$results <- renderDataTable({
    if(input$DESeqtable == TRUE) {
    dds.res
    }
  })
  
  output$volplot <- 
    renderGirafe({
      colors <- c(magma(15)[9], "grey", viridis(15)[10] )#object for colors on volcano based on user input called from palette module
      if(input$DESeqvolcano == TRUE) { #only create plot if the  volcano switch is toggled
        p<- ggplot(dds.res, aes( #call in the DE results from the DE module
          x = `log2FoldChange`,
          y = -log10(padj),
          col = DiffExp,
          tooltip = Gene
        )) +
          geom_point_interactive(size = 1, alpha = 0.5) +
          theme_light() +
          scale_color_manual(values = colors) +
          ggtitle("DE Volcano Plot") +
          coord_cartesian(xlim = c(-10, 7))
        
        girafe(code = print(p))
      }
    })
  
  output$MAplot <- 
    renderGirafe ({
      colors <- c(magma(15)[9], "grey", viridis(15)[10] )#object for colors on volcano based on user input called from palette module
      if(input$DESeqMA == TRUE) { #only call plot if the MA plot switch is toggled
        ma <- ggplot(dds.res, #call in the DE results from the DE module
                     aes(
                       x = log2(baseMean),
                       y = `log2FoldChange`,
                       col = DiffExp,
                       tooltip = Gene
                     )) +
          geom_point_interactive(alpha = 0.8, size = 0.5) +
          geom_hline(aes(yintercept = 0)) +
          scale_color_manual(values = colors) +
          theme_light() +
          ylim(c(
            min(dds.res$`log2FoldChange`),
            max(dds.res$`log2FoldChange`)
          )) +
          ggtitle("DE MA Plot") +
          xlab("log2 Mean Expression") +
          ylab("Log2 Fold Change")
        
        girafe(code = print(ma))
      }
    })
  
  # DEres_reactive <- reactive({
  #  run_DE()
  # })
  list(dds.res)
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
