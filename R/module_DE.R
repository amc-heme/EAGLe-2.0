##Differential Expression tab##
library(colorRamp2)
library(InteractiveComplexHeatmap)
datasets <- 
  read_yaml("./data.yaml")
raw_data <- getwd()
t2g_hs <- read_rds("~/Documents/GitHub/EAGLe-2.0/data/t2g_hs.rds")
t2g_mm <- read_rds("~/Documents/GitHub/EAGLe-2.0/data/t2g_mm.rds")

#choose relevant comparisons for the user to choose from and have them run on the fly
DE_UI <- function(id) {
  ns <- NS(id)
  fluidPage(
    theme =
      shinytheme("flatly"),
    useWaiter(),
    #waiter_on_busy(),
    useShinyjs(),
    titlePanel(
      "Differential Expression Tables and Plots"
    ),#end title
      sidebarLayout(
        sidebarPanel(
          # shinyjs::useShinyjs(),
          # id = "side-panel-de",
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
        #color palette choice for plots
        
        h4("Color Palettes:"),
        colorUI(ns("color"), "Choose 1st color", "#184C52"),
        colorUI(ns("color2"), "Choose 2nd color", "#D53031"),
         hr(),
        conditionalPanel(
          ns = ns,
          condition = "input.DESeqtable == true",
          downloadButton(
            ns("downloadDESeq"),
            label =
              "Download DEG Table"
          )
        ),
        # conditionalPanel(
        #   ns = ns,
        #   condition = "input.DESeqHeat == true",
        #   #js function to hide plot dimensions until selected
        #   materialSwitch(ns("hidedimsHM"), "Custom plot dimensions",
        #                  value = FALSE, right = TRUE),
        #   
        #   shinyjs::hidden(
        #     sliderInput(ns("hmheightslider"),
        #                 "Adjust Plot Height", 200, 1200, 600)),
        #   
        #   shinyjs::hidden(
        #     sliderInput(ns("hmwidthslider"),
        #                 "Adjust Plot Width", 200, 1200, 800)),
        #   
        # )
      )
    ),
        mainPanel(
          conditionalPanel(
            ns = ns,
            condition = "input.DESeqtable == true",
            shinycssloaders::withSpinner(
            DTOutput(ns("results"))
            )
          ),
          conditionalPanel(
            ns = ns,
            condition = "input.DESeqvolcano == true",
            shinycssloaders::withSpinner(
            girafeOutput(ns("volplot"))
            )
          ),
          conditionalPanel(
            ns = ns,
            condition = "input.DESeqMA == true",
            shinycssloaders::withSpinner(
            girafeOutput(ns("MAplot"))
            )
          ),
          conditionalPanel(
            ns = ns,
            condition = "input.DESeqHeat == true",
            plotlyOutput(ns("ht"))
        )
      )
    )
  )
}

DE_Server <- function(id, data_species, dataset_dds, dataset_choice, reset_trigger, vst) {
  moduleServer(id, function(input, output, session) {
    waiter <- Waiter$new()
 
  runDETest_GSEA <- function(dataset, dds, model, comparison) {
    # return dds if LRT is chosen
    if(comparison == "LRT" & dataset %in% c("Cancer_Discovery","Ye_16", "Ye_20", "Venaza",
                                              "Lagadinou", "Lee")) {
      return(results(dds, tidy = TRUE))
    } else if(comparison != "LRT" & dataset %in% c("Cancer_Discovery","Ye_16", "Ye_20", "Venaza",
                                                   "Lagadinou", "Lee")){
      #extract counts and metadata from preloaded dds object
      dds_counts <- counts(dds)
      meta <- colData(dds)
      print("colnames colData:")
      print(head(meta))
      #extract individual levels from the comparison choice
      levels <- unlist(strsplit(comparison, "_vs_"))
      print("levels:")
      print(levels)
      print("eval model:")
      model_term <- as.formula(paste("~", model))
      print(model_term)
      ddsTxi_dds <- DESeqDataSetFromMatrix(dds_counts, colData = meta, design = model_term)
      dds.wald <- DESeq(ddsTxi_dds, test = "Wald") 
      contrasts <- c(model, levels)
      results_df_GSEA <- results(dds.wald, contrast = contrasts, tidy = TRUE)
      #if BEAT or TCGA
     # results_df_GSEA$row <- str_sub(results_df_GSEA$row, end=-4) 
      print("res_GSEA:")
      print(head(results_df_GSEA))
      return(results_df_GSEA)
    } else if(comparison != "LRT" & dataset %in% c("BEAT", "TCGA")) {
      #extract counts and metadata from preloaded dds object
      dds_counts <- counts(dds)
      meta <- colData(dds)
      print("colnames colData:")
      print(head(meta))
      #extract individual levels from the comparison choice
      levels <- unlist(strsplit(comparison, "_vs_"))
      print("levels:")
      print(levels)
      print("eval model:")
      model_term <- as.formula(paste("~", model))
      print(model_term)
      ddsTxi_dds <- DESeqDataSetFromMatrix(dds_counts, colData = meta, design = model_term)
      dds.wald <- DESeq(ddsTxi_dds, test = "Wald") 
      contrasts <- c(model, levels)
      results_df_GSEA <- results(dds.wald, contrast = contrasts, tidy = TRUE)
      #if BEAT or TCGA
      results_df_GSEA$row <- str_sub(results_df_GSEA$row, end=-4) 
      print("res_GSEA:")
      print(head(results_df_GSEA))
      return(results_df_GSEA)
    }
  } #this is what needs to be sent to GSEA#

  runDETest <- function(dds, model, comparison) {
    # return dds if LRT is chosen
    if(comparison == "LRT") {
      return(results(dds))
    } else{
    #extract counts and metadata from preloaded dds object
    dds_counts <- counts(dds)
    meta <- colData(dds)
    print("colnames colData:")
    print(head(meta))
    #extract individual levels from the comparison choice
    levels <- unlist(strsplit(comparison, "_vs_"))
    print("levels:")
    print(levels)
    print("eval model:")
    model_term <- as.formula(paste("~", model))
    print(model_term)
    ddsTxi_dds <- DESeqDataSetFromMatrix(dds_counts, colData = meta, design = model_term)
    dds.wald <- DESeq(ddsTxi_dds, test = "Wald") 
    contrasts <- c(model, levels)
    print("this is contrasts")
    print(contrasts)
    results_df <- results(dds.wald, contrast = contrasts)
    return(results_df)
    }
  } 
 
  
  generateRes <- function(dataset, de_results) {
    #mouse or human?
    is_hs <- grepl("t2g_hs", datasets[[dataset]]$t2g)
    
    if (is_hs & dataset %in% c("BEAT", "TCGA")) {
      res <- data.frame(de_results) %>%
        rownames_to_column(., var = 'ensembl_gene_id')
      
      res$ensembl_gene_id <- str_sub(res$ensembl_gene_id, end = -4)
      
      res <- res %>%
        dplyr::select(., ensembl_gene_id, baseMean, log2FoldChange, padj) %>%
        left_join(unique(dplyr::select(t2g_hs, c(
          ensembl_gene_id, ext_gene
        ))), ., by = 'ensembl_gene_id') %>%
        dplyr::rename(., Gene = ext_gene) %>%
        mutate(., DiffExp = ifelse(
          padj < 0.05 & log2FoldChange >= 0.5,
          'up',
          ifelse(padj < 0.05 &
                   log2FoldChange <= -0.5, 'down', 'no')
        )) %>%
        na.omit(.)
      
    } else if (is_hs & dataset %in% c("Cancer_Discovery", "Venaza",
                                      "Lagadinou", "Lee")) {
      res <- data.frame(de_results) %>%
        rownames_to_column(., var = 'ensembl_gene_id') %>%
        dplyr::select(., ensembl_gene_id, baseMean, log2FoldChange, padj) %>%
        left_join(unique(dplyr::select(t2g_hs, c(
          ensembl_gene_id, ext_gene
        ))), ., by = 'ensembl_gene_id') %>%
        dplyr::rename(., Gene = ext_gene) %>%
        mutate(., DiffExp = ifelse(
          padj < 0.05 & log2FoldChange >= 0.5,
          'up',
          ifelse(padj < 0.05 &
                   log2FoldChange <= -0.5, 'down', 'no')
        )) %>%
        na.omit(.)
      
    } else{
      res <- data.frame(de_results) %>%
        rownames_to_column(., var = 'ensembl_gene_id') %>%
        dplyr::select(., ensembl_gene_id, baseMean, log2FoldChange, padj) %>%
        left_join(unique(dplyr::select(t2g_mm, c(
          ensembl_gene_id, ext_gene
        ))), ., by = 'ensembl_gene_id') %>%
        dplyr::rename(., Gene = ext_gene) %>%
        mutate(., DiffExp = ifelse(
          padj < 0.05 & log2FoldChange >= 0.5,
          'up',
          ifelse(padj < 0.05 &
                   log2FoldChange <= -0.5, 'down', 'no')
        )) %>%
        na.omit(.)
    }
    return(res)
  }
    

  dds_result <- reactiveVal(NULL)
 
  #have DE run when runDE button is clicked
  observeEvent(dataset_choice$close_tab(), {
   
   waiter_show()

      dds_result(runDETest(dataset_dds(), dataset_choice$user_model(), dataset_choice$user_PW()))

    waiter_hide()
  })
  


  observe({
    if(!is.null(dds_result())) { #only run after runDE action button has been selected
     
      dds.res <- generateRes(dataset_choice$user_dataset(), dds_result())
    
      output$results <- renderDataTable({
        if (input$DESeqtable == TRUE) {
          DT::datatable(dds.res,
                        options = list(scrollX = TRUE))
       }
    })
      
      output$downloadDESeq <- downloadHandler(
        filename = function() { paste("DEGTable", '.csv', sep='')},
        content = function(file) {
          write.csv(dds.res,file)
        }
      )
  }
})

  #create objects for color palettes from the palette module
  colorDE <-  
    colorServer("color")
  
  color2DE <-
    colorServer("color2")
  #Volcano Plot ####
  observe({
    if(!is.null(dds_result())) {
     
      res.vol <- generateRes(dataset_choice$user_dataset(), dds_result())
  
  output$volplot <- 
    renderGirafe({
     
      #colors <- c(viridis(15)[10], "grey", magma(15)[9])
      colors <- c(colorDE(), "grey", color2DE())
      if(input$DESeqvolcano == TRUE) { #only create plot if the  volcano switch is toggled
        p<- ggplot(res.vol, aes( 
          x = `log2FoldChange`,
          y = -log10(padj),
          col = DiffExp,
          tooltip = Gene
        )) +
          geom_point_interactive(size = 1, alpha = 0.5) +
          theme_cowplot(font_size = 14) +
          theme(axis.title = element_text(face = "bold"), title = element_text(face = "bold")) +
          scale_color_manual(values = colors) +
          ggtitle("DE Volcano Plot") +
          coord_cartesian(xlim = c(-10, 7))
        
        girafe(code = print(p))
      }
    })
    }
  })

  
  observe({
    if(!is.null(dds_result())) {
     
      res.ma <- generateRes(dataset_choice$user_dataset(), dds_result())  
      
  output$MAplot <- 
    renderGirafe ({
      #colors <- c(viridis(15)[10], "grey", magma(15)[9])
      colors <- c(colorDE(), "grey", color2DE())
      #colors <- c(color3DE(), "grey", color4DE())#object for colors on volcano based on user input called from palette module
      if(input$DESeqMA == TRUE) { #only call plot if the MA plot switch is toggled
        ma <- ggplot(res.ma, 
                     aes(
                       x = log2(baseMean),
                       y = `log2FoldChange`,
                       col = DiffExp,
                       tooltip = Gene
                     )) +
          geom_point_interactive(alpha = 0.8, size = 0.5) +
          geom_hline(aes(yintercept = 0)) +
          scale_color_manual(values = colors) +
          theme_cowplot(font_size = 14) +
          theme(axis.title = element_text(face = "bold"), title = element_text(face = "bold")) +
          ylim(c(
            min(res.ma$`log2FoldChange`),
            max(res.ma$`log2FoldChange`)
          )) +
          ggtitle("DE MA Plot") +
          xlab("log2 Mean Expression") +
          ylab("Log2 Fold Change")
        
        girafe(code = print(ma))
      }
    })
    }
  })
  

    #Heatmap ####
  
  output$ht <- renderPlotly({
    
    req(input$DESeqHeat)
    ns <- NS(id)
    
    res.hm <-
      generateRes(dataset_choice$user_dataset(), dds_result())
    print("res.hm")
    print(head(res.hm))
    print(class(res.hm))
    #filter DE object for only significantly differentially expressed genes
    dds.mat <- res.hm %>%
      dplyr::filter(padj < 0.05 & abs(`log2FoldChange`) >= 2)
    
    #filter vst counts matrix by sig expressed genes
    vst.mat <- vst() %>%
      dplyr::filter(., ensembl_gene_id %in% dds.mat$ensembl_gene_id) %>%
      column_to_rownames(., var = "ensembl_gene_id") %>%
      dplyr::select(., -ext_gene_ensembl) %>%
      as.matrix()
    
    rownames(vst.mat) = dds.mat$Gene
    
    vst.mat <- t(scale(t(vst.mat)))
    #only show the first 100 genes for visualization in this example(can change)
    vst.mat <- head(vst.mat, n = 100)
    print("vst.mat:")
    print(head(vst.mat))
    #create a colorRamp function based on user input in color palette choices
    colors.hm <- colorRamp(c(colorDE(), "white", color2DE()))
    #create heatmap object
    ht <- plot_ly(
      x = colnames(vst.mat),
      y = rownames(vst.mat),
      z = vst.mat,
      colorbar = list(len=1, limits = c(-2, 2)),
      colors = colors.hm,
      type = "heatmap"
    )
    # ht <- 
    #   ComplexHeatmap::Heatmap(
    #     vst.mat,
    #     name = "z scaled expression",
    #     col = colors.hm,
    #     row_names_gp = gpar(fontsize = 4),
    #     row_km = 2,
    #     # top_annotation = HeatmapAnnotation(class = anno_block(gp = gpar(fill = c("white", "white")),
    #     #                                                       labels = c("prim", "mono"),
    #     #                                                       labels_gp = gpar(col = "black", fontsize = 10))),
    #     column_km = 2,
    #     column_title = NULL,
    #     row_title = NULL
    #     
    #   )
    # ht <- draw(ht)
    ht
    #makeInteractiveComplexHeatmap(input, output, session, ht, heatmap_id = ns("ht"))
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
  
  DDS4GSEA <- reactiveVal(NULL)
  
  observeEvent(dataset_choice$close_tab(), { #only run DE for GSEA if the action button is pushed
    waiter$show(
      #spin_fading_circles()
    )

      DDS4GSEA(runDETest_GSEA(dataset_choice$user_dataset(), dataset_dds(),
                              dataset_choice$user_model(),
                              dataset_choice$user_PW()))

   waiter$hide()
  })

     #only generate the res table if all choices have been selected and runDE is pushed
    res4GSEA <- reactive({
      if(!is.null(dds_result())) {
      generateRes(dataset_choice$user_dataset(), dds_result())
    }
  })
    observe({
      req(reset_trigger())
      updateMaterialSwitch(session, "DESeqtable", value = FALSE)
      updateMaterialSwitch(session, "DESeqvolcano", value = FALSE)
      updateMaterialSwitch(session, "DESeqMA", value = FALSE)
      updateMaterialSwitch(session, "DESeqHeat", value = FALSE)
    })
  list(
    res_tidy = reactive(DDS4GSEA()), #send tidy version of dds.res to GSEA for ranks function
    dds_res = res4GSEA #send dds.res to GSEA for use in volcano plot
  )
  })
  }

# need a way to send the results of the de test to GSEA for the volcano plot as well. 
# dds_results function requires the runDEtest function to work. The DE table sent to
# GSEA is in tidy format and will not work for creating a res table
#figure out how to stop dds_result from running in generateGSEAres before button is pushed
