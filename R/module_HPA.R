datasets <- 
  read_yaml("./data.yaml")
raw_data <- getwd()

t2g_hs <- read_rds("~/Documents/GitHub/EAGLe-2.0/data/t2g_hs.rds")
t2g_mm <- read_rds("~/Documents/GitHub/EAGLe-2.0/data/t2g_mm.rds")

HPA_UI <- function(id) {
  ns <- NS(id)
  fluidPage(
    theme=
      shinytheme("flatly"),
    titlePanel(
      "Normal Tissue Expression"
    ),
    # dropMenu(
    #   dropdownButton(circle = TRUE, status = 'info', icon = icon('info'), size = 'sm',
    #                  width = '50px',
    #                  tooltip = tooltipOptions(title = "Information"))),
    sidebarLayout(
      sidebarPanel(
        selectizeInput( #gene choice dropdown menu
          ns("VSTgenechoice"),
          label=
            "Choose a gene for analysis",
          choices =
            NULL,
          selected = NULL,
          options = list(maxItems = 1)
        ), 
        
        #add palette choices for boxplot colors
        paletteUI(ns("palette13")),
        
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
            ns("downloadHPAPlot"),
            label = 
              NULL
          )
        )
      ), 
      
      mainPanel( #add loading spinner
        shinycssloaders::withSpinner(
          plotOutput(
            ns("HPAplot")
          ) 
        )
      )
    )
  )
}

HPA_Server <- function(id, dds.HPA, vst.HPA) {
  moduleServer(id, function(input, output, session) {

    ##Gene Centric output ####
    
    observe({
      gene_choices <- vst.HPA()$ext_gene_ensembl
      updateSelectizeInput(
        session,
        "VSTgenechoice",
        choices = gene_choices,
        selected = NULL,
        server = TRUE
      )
    })
    
    vst.HPA.create <- function(dataset, dds, vst, gene) {
    
        cond_var <- datasets[["HPA"]]$PCA_var
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
      
      return(vst.goi)
    }
    
    vst.gene <- reactive({
      req(input$VSTgenechoice)
      
      vst.goi <- vst.HPA.create(datasets[["HPA"]], dds.HPA(), vst.HPA(), input$VSTgenechoice)
      gene_choice <- input$VSTCDgenechoice
      vst.goi
    })
    
    colorpaletteHPA <- 
      paletteServer("palette13")
    
    
    observe({
      shinyjs::toggle("plotwidthslider", condition = input$hidedims)
    })
    
    observe({
      shinyjs::toggle("plotheightslider", condition = input$hidedims)
    })
    
    output$HPAplot <-
      renderPlot(
        width = function() input$plotwidthslider,
        height = function() input$plotheightslider,
        {
    #plot output
      ggplot(vst.gene(),
                 aes(
                   x = condition,
                   y = ext_gene_ensembl,
                   #color = condition,
                   fill = condition
                 )) +
            geom_boxplot() +
            #Gene_facet() + #reactive faceting
            scale_fill_viridis_d(option = colorpaletteHPA()) + #reactive  scale_fill_manual from module
            scale_color_viridis_d(option = colorpaletteHPA()) + #reactive scale_color_manual from module
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
            #scale_y_continuous(breaks = seq(12, 20, by = 1)) +
            #sig_label_position() + # function for adjusted pvalues position and format on plot
            ylab(input$VSTgenechoice) +
            xlab("") +
            ggtitle("Normal Tissue Expression")
        }) #end render plot
    
    output$downloadHPAPlot <- downloadHandler(
      filename = paste('HPAPlot','.png', sep=''),
      content = function(file) {
        ggsave(file, device = "png", width = 10,
               height = 8, dpi = 100, bg = "white")
      }
    )
    
  })
}

