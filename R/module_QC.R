## QC tab module
library(viridis)
library(ggplot2)
library(DESeq2)
library(TidyMultiqc)
library(cowplot)
datasets.pca <- 
  read_yaml("./data.yaml")

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
        
        conditionalPanel(
          ns = ns,
          condition = "input['multiqc'] == true",
          selectInput( #choose type of multiqc test to visualize
            ns("QCvar"),
            label=
              "Choose MultiQC table and plots",
            choices =
              c("Data Table", "% mapped reads"),
            selected =
              "Data Table"
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
          DTOutput(ns("QCdt"))
        )
      )
    )
  )
}

QC_Server <- function(id, dataset_dds, dataset_choice, qc_table) {
  moduleServer(id, function(input, output, session) {
##PCA plot functions ####
    # function for creating pca from vsd
    
    batch_vsd <- function(dds, dataset) {
      
      dds.file <- dds 
 
      meta <- colData(dds.file) #get colData from dds
      print("colnames colData:")
      print(head(meta)) #make sure it's correct
      
      batch_variable <- datasets.pca[[dataset]]$batch_var #get batch variable from yaml
      batch_variable <- gsub("\"", "", batch_variable) #remove quotes
      print("batch_var:")
      print(batch_variable)
      
      
      pca.id <- datasets.pca[[dataset]]$ID #define ID column using yaml
  
      has_batch <- grepl("yes", datasets.pca[[dataset]]$batch)
      print(has_batch)
      
      if(has_batch) {
        # variance stabilize the counts table
        vsd <- 
          vst(dds.file, blind = F)
        
        batch1 <- meta[, batch_variable] #define batch variable from colData
        print("batch1:")
        print(batch1)
        
        condition_variable <- datasets.pca[[dataset]]$PCA_var

        # model.form <- create.formula(outcome.name = NULL,
        #                              input.names = condition_variable,
        #                              dat = meta)
        #condition_variable <- as.formula(paste("~", condition_variable))
        condition_dat <- colData(vsd)[[condition_variable]]
        design_matrix <- model.matrix(~condition_dat)
        print("model:")
        print(head(design_matrix))
        
        assay(vsd) <- limma::removeBatchEffect(assay(vsd),
                                               batch = batch1,
                                               design = design_matrix)
        
        # assay(vsd) <- ComBat(assay(vsd), batch = batch1, mod=model.matrix(model.form,
        #                                                                   data = colData(vsd))) #remove batch effect
        
        vpca <- 
          data.frame(prcomp(t(assay(vsd)))$x) %>%
          as_tibble(rownames = "ID") %>%
          left_join(., as_tibble(colData(vsd)), by = c("ID" = pca.id))
        
      } else {
        
        # variance stabilize the counts table
        vsd <- 
          #print(head(vst(dds(), blind = F)))
          vst(dds.file, blind = F)
        
        vpca <- 
          data.frame(prcomp(t(assay(vsd)))$x) %>%
          as_tibble(rownames = "ID") %>%
          left_join(., as_tibble(colData(vsd)), by = c("ID" = pca.id))
      }
      print(head(vpca))
      return(vpca)
    }
 
    # vsd.pca for plotting, only run function when PCAplots button is selected
    vsd.pca <- eventReactive(input$PCAplots, { 
      batch_vsd(dataset_dds(), dataset_choice())
    })
    
    
    vsd_fun <- function(dds) {
      dds.file <- dds
      vsd <- 
        vst(dds.file, blind = F)
      return(vsd)
    }
   
    #determine % variance of pc1 and pc2
    pc_variance <- function(dds) {
      vsd <- vsd_fun(dds)
      vsd.pca.var <- prcomp(t(assay(vsd)))
      pc_variance <- data.frame(pc1 = round(vsd.pca.var$sdev[1] ^ 2 / sum(vsd.pca.var$sdev ^ 2)*100,1),
                                pc2 = round(vsd.pca.var$sdev[2] ^ 2 / sum(vsd.pca.var$sdev ^ 2)*100,1))
     return(pc_variance)
    }
    
    pc_variance_df <- reactive({
      pc_variance(dataset_dds())
    })
    
    pc1 <- reactive({
      pc_variance_df()$pc1
      })
    pc2 <- reactive({
      pc_variance_df()$pc2
      })
    
    
    #variable to use for color distinctions on plot
    color_var <- function(dataset, dds) {
      color_v <- datasets.pca[[dataset]]$PCA_var
      meta <- colData(dds)
      color_var <- meta[, color_v]
      return(color_var)
    }
    
    #variable to use for shape distinctions on plot
    shape_var <- function(dataset, dds) {
      shape_v <- datasets.pca[[dataset]]$PCA_shape
      meta <- colData(dds)
      shape_var <- meta[, shape_v]
      return(shape_var)
    } #need to change this to only include shapes if value is not NULL
     
    #run color function after dataset is chosen by user
    pca_color <- reactive({
      pca_color <- color_var(dataset_choice(), dataset_dds())
      print("pca_color")
      print(pca_color)
    })
    
    #run shape function after dataset is chosen by user
    pca_shape <- reactive({
      pca_shape <- shape_var(dataset_choice(), dataset_dds())
      print("pca_shape")
      print(pca_shape)
    })
    
    #call in color palette server for use in plot
    colorpaletteQC <- 
      paletteServer("palette")


##PCA plot ####
    output$PCAplot <- renderPlot ({
      if(input$PCAplots == TRUE) {
        pca <- ggplot(vsd.pca(), aes(x = PC1, y = PC2, color = pca_color(), shape = pca_shape(), fill = pca_color())) +
          geom_point(size = 5) + 
          scale_shape() +
          scale_fill_viridis_d(option = colorpaletteQC()) + #scale_fill_manual reactive function
          scale_color_viridis_d(option = colorpaletteQC()) + #scale_color manual reactive function
          theme_cowplot(font_size = 18) + 
          theme(axis.title = element_text(face = "bold"), title = element_text(face = "bold")) +
          theme(plot.background = element_rect(fill = "#FFFFFF", colour = "#FFFFFF")) +
          theme(panel.background = element_rect(fill = "#FFFFFF", colour = "#FFFFFF")) +
           xlab(paste('PC1 =', pc1(), '% variance')) + #reactive x lab for % variance
           ylab(paste('PC2 =', pc2(), '% variance')) + #reactive y lab for % variance
          ggtitle("PCA on VSD") +
          geom_text_repel(colour = "black", aes(label= ID,hjust=0, vjust=0)) +
          labs(color = "Condition", shape = "Batch") +
          guides(fill = FALSE)
        print(pca)
      }
    })
    
    
 #MultiQC plot function ####   
    
    
    #reactive function for multiqc plot title
    # QC_title <- 
    #   reactive({
    #     if (input$QCvar == "% mapped reads") {
    #       print("% mapped reads per sample")
    #     } else if (input$QCvar == "# mapped reads") {
    #       print("# mapped reads per sample")
    #     } else if (input$QCvar == "% uniquely mapped reads") {
    #       print("% uniquely mapped reads per sample")
    #     } else if (input$QCvar == "# uniquely mapped reads") {
    #       print("# uniquely mapped reads per sample")
    #     }
    #   })
    # 
    # #create object for reactive data input based on user choice of multiqc test option
    # QCdata <- reactive ({
    #   if(input$QCvar == "% mapped reads") {
    #     qc()$raw.salmon.percent_mapped
    #   } else if(input$QCvar == "# mapped reads") {
    #     qc()$raw.salmon.num_mapped
    #   } else if(input$QCvar == "% uniquely mapped reads") {
    #     qc()$raw.star.uniquely_mapped_percent
    #   } else if(input$QCvar == "# uniquly mapped reads") {
    #     qc()$raw.star.uniquely_mapped
    #   }
    # })
    # 
   
    
    output$QCdt <- renderDataTable({
      if (input$multiqc == TRUE) {
        # venaza dataset does not have the same column names-need to add a condition for that dataset
        # BEAT and TCGA do not have multiqc data
        multiqc_res <- qc_table()
        
        multiqc_res <- multiqc_res %>% 
          mutate_at(vars(-Sample), ~round(., digits = 3))
        
        multiqc_res <- DT::datatable(multiqc_res,
                                     options = list(scrollX = TRUE))
       
        multiqc_res
      }
    })
    

  
    # output$QCplot <- renderPlot ({
    #   if(input$multiqc == TRUE) {
    #     ggplot(
    #       qc_table(),
    #       aes(
    #         x = Sample,
    #         y = QCdata()
    #       )) +
    #       geom_point() +
    #       theme_cowplot (font_size = 18) +
    #       ggtitle(QC_title()) +
    #       theme(axis.title = element_text(face = "bold"), title = element_text(face = "bold"), axis.text.x =
    #               element_text(angle = 60, hjust = 1)) 
    #   }
    # })
    # ggplot(multiqc_general_stats,
    #        aes(x = Sample, y = Salmon_percent_mapped, color = Source)) +
    #   geom_point() +
    #   theme_cowplot() +
    #   #scale_color_manual(color = colors) +
    #   theme(
    #     axis.title = element_text(face = "bold"),
    #     title = element_text(face = "bold"),
    #     axis.text.x =
    #       element_text(angle = 60, hjust = 1)
    #   ) +
    #   ggtitle("Salmon: % Mapped Reads ") +
    #   xlab("Sample") +
    #   ylab("% mapped") 
    
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

