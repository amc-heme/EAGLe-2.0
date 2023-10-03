## QC tab module
library(viridis)
library(ggplot2)
library(DESeq2)
library(TidyMultiqc)
library(cowplot)

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
        # materialSwitch(
        #   inputId =
        #     ns("PCAscreeplots"),
        #   label =
        #     "Scree",
        #   value =
        #     FALSE,
        #   right =
        #     TRUE
        # ),
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
        
        # conditionalPanel(
        #   ns = ns,
        #   condition = "input['PCAplots'] == true", 
        #   radioButtons( #choose type of PCA plot
        #     ns("PCAvar"),
        #     h4(
        #       "Choose PCA plot"
        #     ),
        #     choices =
        #       list("VST PCA", "VST + batch corrected PCA"),
        #     selected =
        #       "VST PCA"
        #   )
        # ),
        
        conditionalPanel(
          ns = ns,
          condition = "input['multiqc'] == true",
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
        paletteUI(ns("palette"))
      ),
      mainPanel(
        conditionalPanel(
          ns = ns,
          condition = "input['PCAplots'] == true",
          plotOutput(ns("PCAplot")) 
        ),
        # conditionalPanel(
        #   ns = ns,
        #   condition = "input['PCAscreeplots'] == true",
        #   plotOutput(ns("PCAvarplot"))
        # ),
        conditionalPanel(
          ns = ns,
          condition = "input['multiqc'] == true",
          plotOutput(ns("QCplot"))
        )
      )
    )
  )
}

QC_Server <- function(id, GlobalData) {
  moduleServer(id, function(input, output, session) {
    
    #run pca on vsd
    # vsd.pca <- reactive({
    #   data.frame(prcomp(t(assay(vsd())))$x) %>%
    #   as_tibble(rownames = "SRR") %>%
    #   left_join(., as_tibble(colData(vsd)))
    # })
    #data frame for variance
    vsd.pca.var <- reactive({
      data.frame(summary(prcomp(t(assay(vsd()))))$importance) 
    })
    #determine % variance of pc1 and pc2
    
    pc1var = reactive({round(vsd.pca.var()[3,1] * 100, 1)})
    pc2var = reactive({round(vsd.pca.var()[3,2] * 100 - pc1var, 1)})
    
    scree.pca <- reactive({
      prcomp(t(assay(vsd())))
    })
    #batch corrected PCA
    # assay(vsd) <- limma::removeBatchEffect(assay(vsd),
    #                                        batch = batch, 
    #                                        design = model.matrix(~condition, data = metadata))
    # 
    # bcvsd.pca <-
    #   data.frame(prcomp(t(assay(vsd)))$x) %>%
    #   as_tibble(rownames = "SRR") %>%
    #   left_join(., as_tibble(colData(vsd))) %>%
    #   dplyr::select(SRR, batch, condition, sample_name, everything())
    
    #data frame for variance
    #bcpca.var <- data.frame(summary(prcomp(t(assay(vsd))))$importance)
    
    #determine % variance of pc1 and pc2
    # bc1var = round(bcpca.var[3,1] * 100, 1)
    # bc2var = round(bcpca.var[3,2] * 100 - bc1var, 1)
    
    #pca for batch corrected scree plot
    #scree.bcpca <- prcomp(t(assay(vsd)))
    
    # create functions for calling the batch corrected or vsd % variance for x and y labels
    # xlab_PC1 <- reactive ({
    #   if(input$PCAvar == 'VST PCA') {
    #     pc1var
    #   } else if(input$PCAvar == 'VST + batch corrected PCA') {
    #     bc1var
    #   }
    # })
    
    # ylab_PC2 <- reactive ({
    #   if(input$PCAvar == 'VST PCA') {
    #     pc2var
    #   } else if(input$PCAvar == 'VST + batch corrected PCA') {
    #     bc2var
    #   }
    # })
    #call in color palette server for use in plot
    colorpaletteQC <- 
      paletteServer("palette")
    
    #reactive funtion to choose which PCA plot is loaded
    # PCAdata <-
    #   eventReactive(input$PCAvar, {
    #     if (input$PCAvar == "VST PCA") {
    #       vsd.pca
    #     } else if (input$PCAvar == "VST + batch corrected PCA") {
    #       bcvsd.pca
    #     }
    #   })
    #reactive expression to change the title of the PCA plot based on which PCA is loaded
    # PCA_title <- 
    #   reactive({
    #     if (input$PCAvar == "VST PCA") {
    #       print("VST PCA")
    #     } else if (input$PCAvar == "VST + batch corrected PCA") {
    #       print("VST + batch corrected PCA")
    #     }
    #   })
    output$PCAplot <- renderPlot ({
      if(input$PCAplots == TRUE) {
        pca <- ggplot(vsd.pca(), aes(x = PC1, y = PC2, shape = var_1(), color = batch(), fill = batch())) + 
          geom_point(size = 5) + 
          scale_shape_manual(values = c(21, 24), name = '') +
          scale_fill_viridis_d(option = colorpaletteQC()) + #scale_fill_manual reactive function
          scale_color_viridis_d(option = colorpaletteQC()) + #scale_color manual reactive function
          theme_cowplot(font_size = 18) + 
          theme(axis.title = element_text(face = "bold"), title = element_text(face = "bold")) +
          theme(plot.background = element_rect(fill = "#FFFFFF", colour = "#FFFFFF")) +
          theme(panel.background = element_rect(fill = "#FFFFFF", colour = "#FFFFFF")) +
          xlab(paste('PC1 =', pc1var(), '% variance')) + #reactive x lab for % variance
          ylab(paste('PC2 =', pc2var(), '% variance')) + #reactive y lab for % variance
          ggtitle("PCA") + 
          geom_text_repel(colour = "black", aes(label=sample_name),hjust=0, vjust=0)
        print(pca)
      }
    })
    
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
    
    #create object for reactive data input based on user choice of multiqc test option
    QCdata <- reactive ({
      if(input$QCvar == "% mapped reads") {
        qc()$raw.salmon.percent_mapped
      } else if(input$QCvar == "# mapped reads") {
        qc()$raw.salmon.num_mapped
      } else if(input$QCvar == "% uniquely mapped reads") {
        qc()$raw.star.uniquely_mapped_percent
      } else if(input$QCvar == "# uniquly mapped reads") {
        qc()$raw.star.uniquely_mapped
      }
    })
    
    Sample_ID <-reactive({
      qc()$metadata.sample_id
    })
    
    output$QCplot <- renderPlot ({
      if(input$multiqc == TRUE) {
        ggplot(
          qc(),
          aes(
            x = Sample_ID,
            y = QCdata()
          )) +
          geom_point() +
          theme_cowplot (font_size = 18) +
          ggtitle(QC_title()) +
          theme(axis.title = element_text(face = "bold"), title = element_text(face = "bold"), axis.text.x =
                  element_text(angle = 60, hjust = 1)) 
      }
    })
    
    # PCA Scree data ####
    # VSD PCA variance
    #write functions and store in object to calculate the % variance for each PC
    # PC_var_VST <- data.frame(PC =paste0("PC", 1:num_PCs),variance =(((scree.pca$sdev) ^ 2 / sum((scree.pca$sdev) ^ 2)) * 100))
    # lorder_VST <- as.vector(outer(c("PC"), 1:num_PCs, paste, sep = ""))
    # PC_var_VST$PC <-factor(PC_var_VST$PC,levels = lorder_VST)
    
    #batch corrected PCA variance
    # PC_var_bc <-data.frame(PC =paste0("PC", 1:12),variance =(((scree.bcpca$sdev) ^ 2 / sum((scree.bcpca$sdev) ^ 2)) * 100))
    # lorder_bc <-as.vector(outer(c("PC"), 1:12, paste, sep = ""))
    # PC_var_bc$PC <-factor(PC_var_bc$PC,levels = lorder_bc)
    #function to tell ggplot which data set to use for the scree plots
    # PC_var_data <-
    #   eventReactive(input$PCAvar, {
    #     if (input$PCAvar == "VST PCA") {
    #       PC_var_VST
    #     } else if (input$PCAvar == "VST + batch corrected PCA") {
    #       PC_var_bc
    #     }
    #   })
    #reactive function for scree plots title
    # PCA_var_title <- 
    #   reactive({
    #     if (input$PCAvar == "VST PCA") {
    #       print("VST PC variance")
    #     } else if (input$PCAvar == "VST + batch corrected PCA") {
    #       print("VST + batch corrected PC variance")
    #     }
    #   })
    
    # PCA scree plot ####
    # output$PCAvarplot <- renderPlot ({
    #   if(input$PCAscreeplots == TRUE) {
    #     
    #     ggplot(PC_var_VST,
    #            aes(x = PC,
    #                y = variance,
    #                group = 2)) +
    #       geom_point(size = 2) +
    #       geom_line() +
    #       theme_cowplot(font_size = 18) +
    #       theme(axis.title = element_text(face = "bold"), title = element_text(face = "bold")) +
    #       labs(x = "PC",
    #            y = "% Variance") +
    #       labs(title =
    #              "PCA variability plot")
    #   }
    # })
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
