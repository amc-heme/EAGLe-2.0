datasets <- 
  read_yaml("./data.yaml")
raw_data <- getwd()
t2g_hs <- read_rds("data/t2g_hs.rds")
t2g_mm <- read_rds("data/t2g_mm.rds")

goi_UI <- function(id) {
  ns <- NS(id)
  fluidPage(
    theme=
      shinytheme("flatly"),
    titlePanel(
      "Gene Expression"
    ),
    sidebarLayout(
      sidebarPanel(
        selectizeInput( #gene choice dropdown menu
          ns("VSTCDgenechoice"),
          label=
            "Choose a gene for analysis",
          choices =
            NULL,
          selected = NULL,
          options = list(maxItems = 1)
        ), 
        
        #add palette choices for boxplot colors
        paletteUI(ns("palette2")),
        
        #js functions to hide plot dimensions until selected
        
        materialSwitch(ns("hidedims"),
                       label = "Custom Plot Dimensions",
                       value = FALSE,
                       right = TRUE),
        
        shinyjs::hidden(sliderInput(ns("plotheightslider"),
                                    "Adjust Plot Height", 200, 1200, 600)),
        
        shinyjs::hidden(sliderInput(ns("plotwidthslider"),
                                    "Adjust Plot Width", 200, 1200, 800)),
        
        #dropdown menu containing download button
        dropdownButton(
          inputId = ns("download_menu"),
          label = "Download",
          icon = icon("sliders"),
          status = "primary",
          circle = FALSE,
          downloadButton(
            ns("downloadGenePlot"),
            label = 
              NULL
          )
        )
      ), 
      
      mainPanel( #add loading spinner
        shinycssloaders::withSpinner(
          plotOutput(
            ns("VSTCDplot")
          ) 
        )
      )
    )
  )
}

goi_Server <- function(id, dataset_choice, dataset_dds, vst) {
  moduleServer(id, function(input, output, session) {

    ##Gene Centric output ####
    
    observe({
      gene_choices <- vst()$ext_gene_ensembl
      updateSelectizeInput(
        session,
        "VSTCDgenechoice",
        choices = gene_choices,
        selected = NULL,
        server = TRUE
      )
    })
    
    
    vst.goi.create <- function(dataset, dataset_model, dds, vst, gene) {
      
      if(dataset %in% c("Ye_16", "Venaza", "Lagadinou", "BEAT_quantile", "BEAT_FAB", "BEAT_Denovo.Relapse", "TCGA_FAB", "TCGA_NPM1", "TCGA_RAS")){
        
        cond_var <- dataset_model
        meta <- colData(dds)
        cond <- meta[, cond_var]
        vst.goi <- vst %>%
          dplyr::filter(ext_gene_ensembl %in% gene) %>%
          dplyr::select(., -ensembl_gene_id) %>%
          t(.) %>%
          row_to_names(row_number = 1) %>%
          as.data.frame(.) %>%
          rownames_to_column(var = "Sample") %>% 
          dplyr::mutate(condition = cond) %>%
          dplyr::rename("ext_gene_ensembl" = gene)
        
        vst.goi$ext_gene_ensembl <- as.numeric(vst.goi$ext_gene_ensembl)
        
      } else{
        cond_var <- datasets[[dataset]]$PCA_var
        meta <- colData(dds)
        cond <- meta[, cond_var]
        vst.goi <- vst %>%
          dplyr::filter(ext_gene_ensembl %in% gene) %>%
          dplyr::select(., -ensembl_gene_id) %>%
          t(.) %>%
          row_to_names(row_number = 1) %>%
          as.data.frame(.) %>%
          rownames_to_column(var = "Sample") %>% 
          dplyr::mutate(condition = cond) %>%
          dplyr::rename("ext_gene_ensembl" = gene)
        
        vst.goi$ext_gene_ensembl <- as.numeric(vst.goi$ext_gene_ensembl)
      }
      return(vst.goi)
    }
    
  
    vst.gene <- reactive({
      req(input$VSTCDgenechoice)
      
      vst.goi <- vst.goi.create(dataset_choice$user_dataset(), dataset_choice$user_model(), dataset_dds(), vst(), input$VSTCDgenechoice)
      gene_choice <- input$VSTCDgenechoice
      vst.goi
    })
 
    #call in palette module for plot
    colorpaletteGene <- 
      paletteServer("palette2")

    observe({
            shinyjs::toggle("plotwidthslider", condition = input$hidedims)
          })
    
    observe({
      shinyjs::toggle("plotheightslider", condition = input$hidedims)
    })
    #plot output
    output$VSTCDplot <-
      renderPlot(
        width = function() input$plotwidthslider,
        height = function() input$plotheightslider,
        
        {
  
          ggplot(vst.gene(),
                 aes(
                   x = condition,
                   y = ext_gene_ensembl,
                   fill = condition
                 )) +
            geom_boxplot() +
            scale_fill_viridis_d(option = colorpaletteGene()) + #reactive  scale_fill_manual from module
            scale_color_viridis_d(option = colorpaletteGene()) + #reactive scale_color_manual from module
            geom_point(alpha = 0.5,
                       position = position_jitterdodge(jitter.width = 0.2),
                       aes(color = condition)) +
            theme_cowplot(font_size = 14) +
            theme(
              axis.title = element_text(face = "bold"),
              title = element_text(face = "bold"),
              axis.text.x = element_text(face = "bold", angle = 60, hjust = 1), 
              axis.text.y = element_text(face = "bold")) +
            theme(panel.background = element_rect(fill = "#FFFFFF", colour = "#FFFFFF")) +
            ylab(input$VSTCDgenechoice) +
            xlab("") +
            ggtitle("Gene Expression")
        }) #end render plot
    
    output$downloadGenePlot <- downloadHandler(
      filename = paste('GeneCentricPlot','.png', sep=''),
      content = function(file) {
        goi.plot <- ggplot(vst.gene(),
               aes(
                 x = condition,
                 y = ext_gene_ensembl,
                 fill = condition
               )) +
          geom_boxplot() +
          scale_fill_viridis_d(option = colorpaletteGene()) + #reactive  scale_fill_manual from module
          scale_color_viridis_d(option = colorpaletteGene()) + #reactive scale_color_manual from module
          geom_point(alpha = 0.5,
                     position = position_jitterdodge(jitter.width = 0.2),
                     aes(color = condition)) +
          theme_cowplot(font_size = 14) +
          theme(
            axis.title = element_text(face = "bold"),
            title = element_text(face = "bold"),
            axis.text.x = element_text(face = "bold", angle = 60, hjust = 1), 
            axis.text.y = element_text(face = "bold")) +
          theme(panel.background = element_rect(fill = "#FFFFFF", colour = "#FFFFFF")) +
          ylab(input$VSTCDgenechoice) +
          xlab("") +
          ggtitle("Gene Expression")
        ggsave(goi.plot, file = file, device = "png", width = 8,
               height = 8, dpi = 100, bg = "white")
      }
    )
    
  })
}

