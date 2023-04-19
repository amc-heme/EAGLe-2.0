## QC tab module
library(viridis)
library(ggplot2)
library(DESeq2)
source("~/Documents/GitHub/EAGLe-2.0/config.R")
base_dir <- config$base_dir
samples <- config$samples
sample_id <- config$sample_id
t2g_hs <- read.table(file = config$t2g_hs_file, sep = "\t", header = T)
tx2gene <- t2g_hs[,c(1,2)]
colnames(tx2gene) <- c('TXNAME', 'GENEID')
salm_dirs <- sapply(sample_id, function(id) file.path(base_dir, id, 'quant.sf'))
txi <- tximport(salm_dirs, type = 'salmon', tx2gene = tx2gene, ignoreTxVersion = TRUE)
ddsTxi <- DESeqDataSetFromTximport(txi, colData = samples, design = ~batch + condition)
vsd <- vst(ddsTxi, blind = F)

QC_UI <- function(id) {
  ns <- NS(id)
  fluidPage(
    theme =
      shinytheme("flatly"),
    titlePanel(
      "QC Analysis Plots"
    ),#end title
    sidebarLayout(
      sidebarPanel(
        materialSwitch(
          inputId =
            ns("PCAplots"),
          label =
            "PCA",
          value =
            FALSE,
          right =
            TRUE
        ),
        materialSwitch(
          inputId =
            ns("PCAscreeplots"),
          label =
            "Scree",
          value =
            FALSE,
          right =
            TRUE
        ),
        materialSwitch(
          inputId =
            ns("multiqc"),
          label =
            "MultiQC",
          value =
            FALSE,
          right =
            TRUE
        ),

        conditionalPanel(
          ns = ns,
          condition = "input['PCAplots'] == TRUE", 
          radioButtons( #choose type of PCA plot
            ns("PCAvar"),
            h4(
              "Choose PCA plot"
            ),
            choices =
              list("VST PCA", "VST + batch corrected PCA"),
            selected =
              "VST PCA"
          )
        ),
        
        conditionalPanel(
          ns = ns,
          condition = "input['multiqc'] == TRUE",
          selectInput( #choose type of multiqc test to visualize
            ns("QCvar"),
            label=
              "Choose MultiQC test",
            choices =
              c("% mapped reads", "# mapped reads", "% uniquely mapped reads", "# uniquely mapped reads"),
            selected =
              "% mapped reads"
          )
        ),
        paletteUI("palette")
      ),
        mainPanel(

            plotOutput(ns("PCAplot"))
         
          #   plotOutput(ns("PCAvarplot")),
        
          #   plotOutput(ns("QCplot"))
        )
      )
    )
}

QC_Server <- function(id, colorpaletteQC) {
  moduleServer(id, function(input, output, session) {
    
   #run pca on vsd
    vsd.pca <- data.frame(prcomp(t(assay(vsd)))$x) %>% 
      as_tibble(rownames = "SRR") %>% 
      left_join(., as_tibble(colData(vsd)))
    #data frame for variance
    vsd.pca.var <- data.frame(summary(prcomp(t(assay(vsd))))$importance) 
    #determine % variance of pc1 and pc2
    pc1var = round(vsd.pca.var[3,1] * 100, 1)
    pc2var = round(vsd.pca.var[3,2] * 100 - pc1var, 1)
    
    #batch corrected PCA
    assay(vsd) <- limma::removeBatchEffect(assay(vsd),
                                           batch = samples$batch, 
                                           design = model.matrix(~condition, data = samples))
    bcvsd.pca <- data.frame(prcomp(t(assay(vsd)))$x) %>% 
      as_tibble(rownames = "SRR") %>% 
      left_join(., as_tibble(colData(vsd)))
    #call in color palette server for use in plot
    colorpaletteQC <- 
      paletteServer("palette")
    
    #reactive funtion to choose which PCA plot is loaded
    PCAdata <-
      eventReactive(input$PCAvar, {
        if (input$PCAvar == "VST PCA") {
          vsd.pca
        } else if (input$PCAvar == "VST + batch corrected PCA") {
          bcvsd.pca
        }
      })
    #reactive expression to change the title of the PCA plot based on which PCA is loaded
    PCA_title <- 
      reactive({
        if (input$PCAvar == "VST PCA") {
          print("VST PCA")
        } else if (input$PCAvar == "VST + batch corrected PCA") {
          print("VST + batch corrected PCA")
        }
      })
    output$PCAplot <- renderPlot ({
     if(input$PCAplots == TRUE) {
     pca <- ggplot(PCAdata(), aes(x = PC1, y = PC2, shape = condition, color = batch, fill = batch)) + 
        geom_point(size = 5) + 
        scale_shape_manual(values = c(21, 24), name = '') +
        scale_fill_viridis_d(option = "viridis") + #scale_fill_manual reactive function
        scale_color_viridis_d(option = "viridis") + #scale_color manual reactive function
        theme_cowplot(font_size = 18) + 
        theme(axis.title = element_text(face = "bold"), title = element_text(face = "bold")) +
        theme(plot.background = element_rect(fill = "#FFFFFF", colour = "#FFFFFF")) +
        theme(panel.background = element_rect(fill = "#FFFFFF", colour = "#FFFFFF")) +
        xlab(paste('PC1 =', pc1var, '% variance')) +
        ylab(paste('PC2 =', pc2var, '% variance')) +
        ggtitle(PCA_title()) + #reactive title
        geom_text_repel(colour = "black", aes(label=sample_name),hjust=0, vjust=0)
      print(pca)
     }
    })
  })
}

QC_App <- function() {
  ui <- fluidPage(
    QC_UI("QC1")
  )
  server <- function(input, output, session) {
    QC_Server("QC1")
  }
  shinyApp(ui, server)
}
QC_App()
  