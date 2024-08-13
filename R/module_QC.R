## QC tab module

datasets.pca <- 
  read_yaml("./data.yaml")

QC_UI <- function(id) {
  ns <- NS(id)
  fluidPage(
    theme =
      shinytheme("flatly"),
    titlePanel(
      "QC Analysis Plots"
    ),
    sidebarLayout(
      sidebarPanel(
        materialSwitch(
          inputId =
            ns("PCAplots"),
          label =
            "PCA plot",
          value =
            FALSE,
          right =
            TRUE
        ),
        materialSwitch(
          inputId =
            ns("PCAscreeplots"),
          label =
            "Scree plot",
          value =
            FALSE,
          right =
            TRUE
        ),
        materialSwitch(
          inputId =
            ns("multiqc"),
          label =
            "MultiQC table",
          value =
            FALSE,
          right =
            TRUE
        ),
        paletteUI(ns("palette")),
        #dropdown menu containing download buttons                                  
        dropdownButton(
          inputId = ns("download_menu"),
          label = "Download",
          icon = icon("sliders"),
          status = "primary",
          circle = FALSE,
          downloadButton(
            ns("downloadPCA"),
            label =
              "PCA"), 
          
          downloadButton(
            ns("downloadScree"),
            label =
              "Scree"),
          
          downloadButton(
            ns("downloadMultiQC"),
            label = 
              "MultiQC")
        )     
    ),
       
      mainPanel(
        conditionalPanel(
          ns = ns,
          condition = "input['PCAplots'] == true",
          shinycssloaders::withSpinner(
            plotOutput(ns("PCAplot")) 
          )
        ),
        conditionalPanel(
          ns = ns,
          condition = "input['PCAscreeplots'] == true",
          shinycssloaders::withSpinner(
          plotOutput(ns("PCAvarplot"))
          )
        ),
        conditionalPanel(
          ns = ns,
          condition = "input['multiqc'] == true",
          
          DTOutput(ns("QCdt")),
          textOutput(ns("QCdt2"))
          
        )
      )
    )
  )
}

QC_Server <- function(id, dataset_dds, dataset_choice, qc_table, reset_trigger) {
  moduleServer(id, function(input, output, session) {
##PCA plot functions ####
    # function for creating pca from vsd
   
    batch_vsd <- function(dds, dataset) {
      
      dds.file <- dds 
 
      meta <- colData(dds.file) #get colData from dds
 
      batch_variable <- datasets.pca[[dataset]]$batch_var #get batch variable from yaml
      batch_variable <- gsub("\"", "", batch_variable) #remove quotes
    
      pca.id <- datasets.pca[[dataset]]$ID #define ID column using yaml
  
      has_batch <- grepl("yes", datasets.pca[[dataset]]$batch)
      
      if(has_batch) {
        # variance stabilize the counts table
        vsd <- 
          vst(dds.file, blind = F)
        batch1 <- meta[, batch_variable] #define batch variable from colData
        condition_variable <- datasets.pca[[dataset]]$PCA_var
        condition_dat <- colData(vsd)[[condition_variable]]
        design_matrix <- model.matrix(~condition_dat)
        print("model:")
        print(head(design_matrix))

        assay(vsd) <- limma::removeBatchEffect(assay(vsd),
                                               batch = batch1,
                                               design = design_matrix)
        vpca <- 
          data.frame(prcomp(t(assay(vsd)))$x) %>%
          as_tibble(rownames = "ID") %>%
          left_join(., as_tibble(colData(vsd)), by = c("ID" = pca.id))
        
      } else{
        # variance stabilize the counts table
        vsd <- 
          vst(dds.file, blind = F)
        
        vpca <- 
          data.frame(prcomp(t(assay(vsd)))$x) %>%
          as_tibble(rownames = "ID") %>%
          left_join(., as_tibble(colData(vsd)), by = c("ID" = pca.id))
        
      } 
      return(vpca)
    }
 
    # vsd.pca for plotting, only run function when PCAplots button is selected
    vsd.pca <- eventReactive(input$PCAplots, { 
      batch_vsd(dataset_dds(), dataset_choice$user_dataset())
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
    color_var <- function(dataset, dds, model) {
      if(dataset %in% c("Cancer_Discovery", "Ye_16", "Ye_20", "Venaza", "Lagadinou", "Lee")){
        color_v <- datasets.pca[[dataset]]$PCA_var
        meta <- colData(dds)
        color_var <- meta[, color_v]
        color_var <- factor(color_var)
      } else if(dataset %in% c("BEAT_quantile", "BEAT_FAB", "BEAT_Denovo.Relapse", "TCGA_FAB", "TCGA_NPM1", "TCGA_RAS")){
        color_v <- model
        meta <- colData(dds)
        color_var <- meta[, color_v]
        color_var <- factor(color_var)
      }
      return(color_var)
    }
    
    #variable to use for shape distinctions on plot
    shape_var <- function(dataset, dds, model) {
      if(dataset %in% c("Cancer_Discovery", "Ye_16", "Ye_20", "Venaza", "Lagadinou", "Lee")) {
        shape_v <- datasets.pca[[dataset]]$PCA_shape
        meta <- colData(dds)
        shape_var <- meta[, shape_v]
        shape_var <- factor(shape_var)
        
      } else if(dataset %in% c("BEAT_quantile", "BEAT_FAB", "BEAT_Denovo.Relapse", "TCGA_FAB", "TCGA_NPM1", "TCGA_RAS")){
        shape_v <- model
        meta <- colData(dds)
        shape_var <- meta[, shape_v]
        shape_var <- factor(shape_var)
      }
      return(shape_var)
    } 
     
    color_label <- function(dataset, model) {
      if(dataset %in% c("Cancer_Discovery", "Ye_16", "Ye_20", "Venaza", "Lagadinou", "Lee")) {
        color_l <- datasets.pca[[dataset]]$PCA_var
      } else if(dataset %in% c("BEAT_quantile", "BEAT_FAB", "BEAT_Denovo.Relapse", "TCGA_FAB", "TCGA_NPM1", "TCGA_RAS")) {
        color_l <- model
      }
      color_l
    }
    
    shape_label <- function(dataset, model) {
      if(dataset %in% c("Cancer_Discovery", "Ye_16", "Ye_20", "Venaza", "Lagadinou", "Lee")){
        shape_l <- datasets.pca[[dataset]]$PCA_shape
      } else if(dataset %in% c("BEAT_quantile", "BEAT_FAB", "BEAT_Denovo.Relapse", "TCGA_FAB", "TCGA_NPM1", "TCGA_RAS")) {
        shape_l <- model
      }
     
      shape_l
    }
    #run color function after dataset is chosen by user
    pca_color <- reactive({
      pca_color <- color_var(dataset_choice$user_dataset(), dataset_dds(), dataset_choice$user_model())
    })
    
    #run shape function after dataset is chosen by user
    pca_shape <- reactive({
      pca_shape <- shape_var(dataset_choice$user_dataset(), dataset_dds(), dataset_choice$user_model())
    })
    
    
    color_legend <- reactive({
      color_label(dataset_choice$user_dataset(), dataset_choice$user_model())
    })
    
    shape_legend <- reactive({
      shape_label(dataset_choice$user_dataset(), dataset_choice$user_model())
    })
    #call in color palette server for use in plot
    colorpaletteQC <- 
      paletteServer("palette")


##PCA plot ####
    output$PCAplot <- renderPlot({
      if(input$PCAplots == TRUE) {
        pca <- ggplot(vsd.pca(), aes(x = PC1, y = PC2, color = pca_color(), shape = pca_shape(), fill = pca_color())) +
          geom_point(size = 5) + 
          scale_shape() +
          scale_fill_viridis_d(option = colorpaletteQC()) + #scale_fill_manual reactive function
          scale_color_viridis_d(option = colorpaletteQC()) + #scale_color manual reactive function
          theme_cowplot(font_size = 14) + 
          theme(axis.title = element_text(face = "bold"), title = element_text(face = "bold")) +
          theme(plot.background = element_rect(fill = "#FFFFFF", colour = "#FFFFFF")) +
          theme(panel.background = element_rect(fill = "#FFFFFF", colour = "#FFFFFF")) +
           xlab(paste('PC1 =', pc1(), '% variance')) + #reactive x lab for % variance
           ylab(paste('PC2 =', pc2(), '% variance')) + #reactive y lab for % variance
          ggtitle("PCA on VSD") +
          labs(color = color_legend(), shape = shape_legend()) + #rename legends
          guides(fill = FALSE) #hide fill legend 
        print(pca)
      }
    })
 
    #download PCA plot ####
    output$downloadPCA <- downloadHandler(
      filename = paste("PCA Plot", '.png', sep=''),
      content = function(file) {
        pca <- ggplot(vsd.pca(), aes(x = PC1, y = PC2, color = pca_color(), shape = pca_shape(), fill = pca_color())) +
          geom_point(size = 5) + 
          scale_shape() +
          scale_fill_viridis_d(option = colorpaletteQC()) + #scale_fill_manual reactive function
          scale_color_viridis_d(option = colorpaletteQC()) + #scale_color manual reactive function
          theme_cowplot(font_size = 14) + 
          theme(axis.title = element_text(face = "bold"), title = element_text(face = "bold")) +
          theme(plot.background = element_rect(fill = "#FFFFFF", colour = "#FFFFFF")) +
          theme(panel.background = element_rect(fill = "#FFFFFF", colour = "#FFFFFF")) +
          xlab(paste('PC1 =', pc1(), '% variance')) + #reactive x lab for % variance
          ylab(paste('PC2 =', pc2(), '% variance')) + #reactive y lab for % variance
          ggtitle("PCA on VSD") +
          labs(color = color_legend(), shape = shape_legend()) + #rename legends
          guides(fill = FALSE) #hide fill legend 
        ggsave(pca, file = file, device = "png", width = 8, height = 6, units = "in",dpi = 100)
      }
    )
    
 #MultiQC table function ####   

    multiqc_res <- reactive({
      multiqc_res <- qc_table()
      
      multiqc_res <- multiqc_res %>%
        mutate_at(vars(-Sample), ~ round(., digits = 3))
      
    })
    
    output$QCdt <- renderDataTable({
      if (input$multiqc == TRUE & dataset_choice$user_dataset() %in%
          c("Cancer_Discovery",
            "Ye_16",
            "Ye_20",
            "Venaza",
            "Lagadinou",
            "Lee")) {
        # venaza dataset does not have the same column names-need to add a
        # condition for that dataset
        
        multiqc_res <- qc_table()
        
        multiqc_res <- multiqc_res %>%
          mutate_at(vars(-Sample), ~ round(., digits = 3))
        
        multiqc_res_DT <- DT::datatable(multiqc_res(),
                                        options = list(scrollX = TRUE))
        
        multiqc_res_DT
      }
    })
    
    output$QCdt2 <- renderText({
      if (input$multiqc == TRUE & dataset_choice$user_dataset() %in%
          c("BEAT", "TCGA")) {
        paste("MultiQC data is unavailable for this dataset")
      }
    })
    
       
    #download multiqc table
        output$downloadMultiQC <- downloadHandler(
          filename = paste("MultiQC_table", '.csv', sep=''),
          content = function(file) {
            shiny::req(!dataset_choice$user_dataset() %in% c("BEAT_quantile", "BEAT_FAB", "BEAT_Denovo.Relapse", "TCGA_FAB", "TCGA_NPM1", "TCGA_RAS"))
            write.csv(multiqc_res(), file = file)
          }
        )
        
        observe({
          toggleState("downloadMultiQC", !dataset_choice$user_dataset() %in% c("BEAT_quantile", "BEAT_FAB", "BEAT_Denovo.Relapse", "TCGA_FAB", "TCGA_NPM1", "TCGA_RAS"))
        })
 
    scree_variance_fun <- function(dds) {
      vsd <- vsd_fun(dds)
      scree.var <- prcomp(t(assay(vsd)))
      num_pcs <- length(scree.var$sdev)
      scree_variance <- data.frame(pc = paste0("PC", 1:num_pcs), variance = (((scree.var$sdev) ^ 2 / sum((scree.var$sdev) ^ 2)) * 100))
      print("scree var:")
      print(head(scree_variance))
      lorder <- as.vector(outer(c("PC"), 1:12, paste, sep = ""))
      scree_variance$pc <- factor(scree_variance$pc, levels = lorder)
      return(scree_variance)
    }
    
    scree_variance_df <- reactive({
      scree_df <- scree_variance_fun(dataset_dds())
      print("scree df:")
      print(head(scree_df))
      scree_df
    })

    # PCA scree plot ####
    output$PCAvarplot <- renderPlot ({
      if(input$PCAscreeplots == TRUE) {

        ggplot(scree_variance_df(),
               aes(x = pc,
                   y = variance,
                   group = 2)) +
          geom_point(size = 2) +
          geom_line() +
          theme_cowplot(font_size = 14) +
          theme(axis.title = element_text(face = "bold"), title = element_text(face = "bold")) +
          labs(x = "PC",
               y = "% Variance") +
          labs(title =
                 "PCA variability plot")
      }
    })
    
    #download scree plot ####
    output$downloadScree <- downloadHandler(
      filename = paste("Scree Plot", '.png', sep=''),
      content = function(file) {
        scree <- ggplot(scree_variance_df(),
               aes(x = pc,
                   y = variance,
                   group = 2)) +
          geom_point(size = 2) +
          geom_line() +
          theme_cowplot(font_size = 14) +
          theme(axis.title = element_text(face = "bold"), title = element_text(face = "bold")) +
          theme(plot.background = element_rect(fill = "#FFFFFF", colour = "#FFFFFF")) +
          theme(panel.background = element_rect(fill = "#FFFFFF", colour = "#FFFFFF")) +
          labs(x = "PC",
               y = "% Variance") +
          labs(title =
                 "PCA variability plot")
        ggsave(scree,file = file, device = "png", width = 8, height = 6, units = "in",dpi = 100)
      }
    )
    
    observe({
      req(reset_trigger())
      updateMaterialSwitch(session, "PCAplots", value = FALSE)
      updateMaterialSwitch(session, "PCAscreeplots", value = FALSE)
      updateMaterialSwitch(session, "multiqc", value = FALSE)
    })
  })
}

