
PW_UI <- function(id) {
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
          shinyWidgets::checkboxGroupButtons(
            inputId = (ns("pwc")),
            label = "Choose a Pairwise Comparison",
            choices = c("1", "2", "3", "4"),
            #choices = NULL,
            selected = NULL,
            status = "primary"
          ),
          materialSwitch(
            inputId =
              (ns("PWDESeqtable")),
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
              (ns("PWvolcano")),
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
              (ns("PWMA")),
            label =
              "MA Plot",
            value =
              FALSE,
            right =
              TRUE
          ),
          materialSwitch(
            inputId =
              (ns("PWHeat")),
            label =
              "Heatmap",
            value =
              FALSE,
            right =
              TRUE
          ),
  
          hr(),
          conditionalPanel(
            ns = ns,
            condition = "input.PWvolcano == true",
            h4("Volcano Plot Specific Options"),
            #color palette choice for volcano plot
            colorUI(ns("color10"), "Choose 1st color", "#0000FF"),
            colorUI(ns("color11"), "Choose 2nd color", "028a0f"),
            
            hr(),
            downloadButton(
              ns("downloadPWVolcano"),
              label =
                "Download Volcano Plot"
            )
          ),
          hr(),
          conditionalPanel(
            ns = ns,
            condition = "input.PWMA == true",
            h4("MA Plot Specific Options"),
            #color palette choice for MA plot
            colorUI(ns("color12"), "Choose 1st color", "#0000FF"),
            colorUI(ns("color13"), "Choose 2nd color", "028a0f"),
            
            hr(),
            downloadButton(
              ns("downloadPWMA"),
              label =
                "Download MA Plot"
            )
          ),
          hr(),
          conditionalPanel(
            ns = ns,
            condition = "input.PWHeat == true",
            h4("Heatmap Specific Options"),
            #color palette choices for heatmap
            colorUI(ns("color14"),"Choose 1st color", "#0000FF"),
            colorUI(ns("color15"), "Choose 2nd color", "#FF0000"),
          )
        )
      ),
      mainPanel(
        conditionalPanel(
          ns = ns,
          condition = "input.PWDEtable == true",
          DTOutput(ns("PWresults"))
        ),
        conditionalPanel(
          ns = ns,
          condition = "input.PWvolcano == true",
          girafeOutput(ns("PWvolplot"))
        ),
        conditionalPanel(
          ns = ns,
          condition = "input.PWMA == true",
          girafeOutput(ns("PWMAplot"))
        ),
        conditionalPanel(
          ns = ns,
          condition = "input.PWHeat == true",
          InteractiveComplexHeatmapOutput(heatmap_id = ns("pwht"))
        )
      )
    )
  )
  }


PW_Server <- function(id, data_species, dataset_dds, dataset_choice) {
  moduleServer(id, function(input, output, session) {
    #all of the pairwise tests need to run on the fly
    
    
    # create a data module for those, and do the contrast here
    ## update pairwise choice for each dataset
    observe({
      if(dataset_choice() == "Cancer_Discovery") {
        updateCheckboxGroupButtons(
          session = session,
          inputId = "pwc",
          choices = c("1", "2"),
          selected = NULL
        )
      } else if(dataset_choice() == "Ye_16"){
        updateCheckboxGroupButtons(
          session = session,
          inputId = "pwc",
          choices = c("blood", "bone marrow", "gonadal adipose tissue", "normal bone marrow", "spleen"),
          selected = NULL
        )
      } else if(dataset_choice() == "Venaza"){
        updateCheckboxGroupButtons(
          session = session,
          inputId = "pwc",
          choices = c("1", "2"),
          selected = NULL
        )
      }
    })

    ## run pairwise contrast depending on pairwise choice for each dataset
    #DE Table ####
    
    ## for BEAT dataset, the ensembl id's need to be modified to work:
    # dds.res1$ensembl_gene_id <- str_sub(dds.res1$ensembl_gene_id, end=-4) 
    
    # reactive to switch between mouse or human t2g
    dds.res <- reactive({
      if(data_species() == "human") {
        
        res <- data.frame(results(dataset_dds())) %>%
          rownames_to_column(., var = 'ensembl_gene_id') %>%
          dplyr::select(., ensembl_gene_id, baseMean, log2FoldChange, padj) %>%
          left_join(unique(dplyr::select(t2g_hs, c(ensembl_gene_id, ext_gene))), ., by = 'ensembl_gene_id') %>%
          dplyr::rename(., Gene = ext_gene) %>%
          mutate(., DiffExp = ifelse(padj < 0.05 & log2FoldChange >= 0.5, 'up',
                                     ifelse(padj < 0.05 & log2FoldChange <= -0.5, 'down', 'no'))) %>%
          na.omit(.)
        
      } else if(data_species() == "mouse") {
        
        res <- data.frame(results(dataset_dds())) %>%
          rownames_to_column(., var = 'ensembl_gene_id') %>%
          dplyr::select(., ensembl_gene_id, baseMean, log2FoldChange, padj) %>%
          left_join(unique(dplyr::select(t2g_mm, c(ensembl_gene_id, ext_gene))), ., by = 'ensembl_gene_id') %>%
          dplyr::rename(., Gene = ext_gene) %>%
          mutate(., DiffExp = ifelse(padj < 0.05 & log2FoldChange >= 0.5, 'up',
                                     ifelse(padj < 0.05 & log2FoldChange <= -0.5, 'down', 'no'))) %>%
          na.omit(.)
      }
      res
    })
    
    output$PWresults <- renderDataTable({
      if (input$PWDESeqtable == TRUE) {
        
        dds.res()
      }
    })
    
    #create objects for color palettes from the palette module
    colorDE <-  
      colorServer("color10")
    
    color2DE <-
      colorServer("color11")
    #Volcano Plot ####
    output$PWvolplot <- 
      renderGirafe({
        #colors <- c(colorDE(), "grey",color2DE()) #object for colors on volcano based on user input
        colors <- c(colorDE(), "grey", color2DE())
        if(input$PWvolcano == TRUE) { #only create plot if the  volcano switch is toggled
          p<- ggplot(dds.res(), aes( #call in the DE results from the DE module
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
      colorServer("color12")
    
    color4DE <-
      colorServer("color13")
    
    output$PWMAplot <- 
      renderGirafe ({
        #colors <- c(color3DE(), "grey",color4DE()) #object for colors on volcano based on user input
        colors <- c(color3DE(), "grey", color4DE())#object for colors on volcano based on user input called from palette module
        if(input$PWMA == TRUE) { #only call plot if the MA plot switch is toggled
          ma <- ggplot(dds.res(), 
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
              min(dds.res()$`log2FoldChange`),
              max(dds.res()$`log2FoldChange`)
            )) +
            ggtitle("DE MA Plot") +
            xlab("log2 Mean Expression") +
            ylab("Log2 Fold Change")
          
          girafe(code = print(ma))
        }
      })
    
    
    #Heatmap ####
    #interactive heatmap needs to be wrapped in a reactive function to work
    observe({
      req(input$DPWHeat)
      ns <- NS(id)
      
      color5DE <-
        colorServer("color14")
      
      color6DE <-
        colorServer("color15")
      #object for batch corrected vsd matrix
      # assay(vsd) <-
      #   limma::removeBatchEffect(assay(vsd),
      #                            batch = samples$batch,
      #                            design = model.matrix(~ condition, data = samples))
      # # data frame for batch corrected vsd matrix
      # vstlimma <-
      #   data.frame(assay(vsd)) %>%
      #   rownames_to_column(., var = "ensembl_gene_id") %>%
      #   left_join(unique(dplyr::select(t2g_hs, c(
      #     ensembl_gene_id, ext_gene
      #   ))), ., by = 'ensembl_gene_id') %>%
      #   na.omit(.)
      
      #filter DE object for only significantly differentially expressed genes
      dds.mat <- dds.res %>%
        dplyr::filter(padj < 0.05 & abs(`log2FoldChange`) >= 2)
      
      #filter vst counts matrix by sig expressed genes
      vst.mat <- vst %>%
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
      # if (isTruthy(input$DESeqHeat)) {
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
      makeInteractiveComplexHeatmap(input, output, session, ht, heatmap_id = ns("ht"))
      #}
    })
    
    # download DE table
    output$downloadPWDEtable <- downloadHandler(
      filename = function() { paste("DESeqTable", '.csv', sep='') },
      content = function(file) {
        write.csv(CD_DE_DT(),file)
      }
    )
    
    #download Volcano
    output$downloadPWVolcano <- downloadHandler(
      filename = function() { paste(input$sigvaluesbutton, '.png', sep='') },
      content = function(file) {
        ggsave(file, device = "png", width = 8, height = 6, units = "in",dpi = 72)
      }
    )
    
    #download MA
    output$downloadPWMA <- downloadHandler(
      filename = function() { paste('DESeqMAplot', '.png', sep='') },
      content = function(file) {
        ggsave(file, device = "png", width = 8, height = 6, units = "in",dpi = 72)
      }
    )
  })
}






















