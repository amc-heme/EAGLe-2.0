##Differential Expression tab##
library(colorRamp2)
source("~/Documents/GitHub/EAGLe-2.0/config.R")


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
            (ns("DESeqtable")),
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
        ),
        materialSwitch(
          inputId =
            (ns("DESeqHeat")),
          label =
            "Heatmap",
          value =
            FALSE,
          right =
            TRUE
        ),
         conditionalPanel(
           ns = ns,
           condition = "input.DESeqtable == true",
          h4("DE Table Specific Options"),
          #option to filter table by padj
         radioButtons(ns("padjbutton"), label = "Filter DE tables by padj",
                      choices = list("<= 0.01" = "sigvar1", "<= 0.05" = "sigvar5", "All" = "allvar"), selected = "allvar"),
         hr(),
         downloadButton(ns("downloadDEtable"), label = "Download DE Table")
         ),

         hr(),
        conditionalPanel(
          ns = ns,
          condition = "input.DESeqvolcano == true",
          h4("Volcano Plot Specific Options"),
          #color palette choice for volcano plot
          colorUI(ns("color"), "Choose 1st color", "#0000FF"),
          colorUI(ns("color2"), "Choose 2nd color", "028a0f"),

          hr(),
          downloadButton(
            ns("downloadDEVolcano"),
            label =
              "Download Volcano Plot"
          )
        ),
        hr(),
        conditionalPanel(
          ns = ns,
          condition = "input.DESeqMA == true",
          h4("MA Plot Specific Options"),
          #color palette choice for MA plot
          colorUI(ns("color3"), "Choose 1st color", "#0000FF"),
          colorUI(ns("color4"), "Choose 2nd color", "028a0f"),

          hr(),
          downloadButton(
            ns("downloadDEMA"),
            label =
              "Download MA Plot"
          )
        ),
        hr(),
        conditionalPanel(
          ns = ns,
          condition = "input.DESeqHeat == true",
          h4("Heatmap Specific Options"),
          #color palette choices for heatmap
          colorUI(ns("color5"),"Choose 1st color", "#0000FF"),
          colorUI(ns("color6"), "Choose 2nd color", "#FF0000"),
              )
        )
             ),
        mainPanel(
          conditionalPanel(
            ns = ns,
            condition = "input.DESeqtable == true",
            DTOutput(ns("results"))
          ),
          conditionalPanel(
            ns = ns,
            condition = "input.DESeqvolcano == true",
            girafeOutput(ns("volplot"))
          ),
          conditionalPanel(
            ns = ns,
            condition = "input.DESeqMA == true",
            girafeOutput(ns("MAplot"))
          ),
          conditionalPanel(
            ns = ns,
            condition = "input.DESeqHeat == true",
            InteractiveComplexHeatmapOutput(heatmap_id = 
                                              ns("ht")
            )
          )
     
        )
      
  )
  )
}

DE_Server <- function(id, dds, vsd) {
  moduleServer(id, function(input, output, session) {
 #DE Table ####
 
    dds.res <- data.frame(results(dds)) %>%
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
  
  #create objects for color palettes from the palette module
  colorDE <-  
    colorServer("color")
  
  color2DE <-
    colorServer("color2")
  #Volcano Plot ####
  output$volplot <- 
    renderGirafe({
      #colors <- c(colorDE(), "grey",color2DE()) #object for colors on volcano based on user input
      colors <- c(colorDE(), "grey", color2DE())
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
  #MA Plot ####
  color3DE <-
    colorServer("color3")
  
  color4DE <-
    colorServer("color4")
  
  output$MAplot <- 
    renderGirafe ({
      print("values of colors in DE module")
      print(color3DE())
      print(color4DE())
      #colors <- c(color3DE(), "grey",color4DE()) #object for colors on volcano based on user input
      colors <- c(color3DE(), "grey", color4DE())#object for colors on volcano based on user input called from palette module
      if(input$DESeqMA == TRUE) { #only call plot if the MA plot switch is toggled
        ma <- ggplot(dds.res, 
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
  

    #Heatmap ####
  color5DE <-
    colorServer("color5")

  color6DE <-
    colorServer("color6")
  #object for batch corrected vsd matrix
  assay(vsd) <-
    limma::removeBatchEffect(assay(vsd),
                             batch = samples$batch,
                             design = model.matrix(~ condition, data = samples))
  # data frame for batch corrected vsd matrix
  vstlimma <-
    data.frame(assay(vsd)) %>%
    rownames_to_column(., var = "ensembl_gene_id") %>%
    left_join(unique(dplyr::select(t2g_hs, c(
      ensembl_gene_id, ext_gene
    ))), ., by = 'ensembl_gene_id') %>% na.omit(.)
  
  #interactive heatmap needs to be wrapped in a reactive function to work
  observe({
    #filter DE object for only significantly differentially expressed genes
    dds.mat <- dds.res %>%
      dplyr::filter(padj < 0.05 & abs(`log2FoldChange`) >= 2)
    
    #filter vst counts matrix by sig expressed genes
    vst.mat <- vstlimma %>%
      dplyr::filter(., ensembl_gene_id %in% dds.mat$ensembl_gene_id) %>%
      column_to_rownames(., var = "ensembl_gene_id") %>%
      dplyr::select(.,-ext_gene) %>%
      as.matrix()
    
    rownames(vst.mat) = dds.mat$Gene
    vst.mat <- t(scale(t(vst.mat)))
    #only show the first 100 genes for visualization in this example(can change)
    vst.mat <- head(vst.mat, n = 100)
    #create a colorRamp function based on user input in color palette choices
    colors = colorRamp2(c(-2, 0, 2), c(color5DE(), "white", color6DE()))
    #create heatmap object
    if(input$DESeqHeat == TRUE) {
      ht = draw(ComplexHeatmap::Heatmap(
        vst.mat,
        name = "z scaled expression",
        col = colors,
        row_names_gp = gpar(fontsize = 4),
        row_km = 2,
        top_annotation = HeatmapAnnotation(class = anno_block(gp = gpar(fill = c("white", "white")),
                                                              labels = c("prim", "mono"), 
                                                              labels_gp = gpar(col = "black", fontsize = 10))),
        column_km = 2, 
        column_title = NULL,
        row_title = NULL
      ))
      makeInteractiveComplexHeatmap(input, output, session, ht, "ht")
    }
  })
  
  # download DE table
  output$downloadDEtable <- downloadHandler(
    filename = function() { paste("DESeqTable", '.csv', sep='') },
    content = function(file) {
      write.csv(CD_DE_DT(),file)
    }
  )

  #download Volcano
  output$downloadDEVolcano <- downloadHandler(
    filename = function() { paste(input$sigvaluesbutton, '.png', sep='') },
    content = function(file) {
      ggsave(file, device = "png", width = 8, height = 6, units = "in",dpi = 72)
    }
  )
  
  #download MA
  output$downloadDEMA <- downloadHandler(
    filename = function() { paste('DESeqMAplot', '.png', sep='') },
    content = function(file) {
      ggsave(file, device = "png", width = 8, height = 6, units = "in",dpi = 72)
    }
  )
  })
}

DE_App <- function() {
  ui <- fluidPage(
    DE_UI("DE1")
  )
  server <- function(input, output, session) {
    DE_Server("DE1", dds, vsd)
  }
  shinyApp(ui, server)
}
DE_App()

