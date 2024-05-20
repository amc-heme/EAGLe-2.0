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
    #dropMenu(
      dropdownButton(
        circle = TRUE,
        status = 'info',
        icon = icon('info'),
        size = 'sm',
        width = '50px',
        tooltip =
          tooltipOptions(title = "Information"),
     # ), 
      h6("1."),
      h6("-"),
      h6("-"),
      br(),
      h6("2."),
      
      br(),
      h6("3."),
      h6("-"),
      h6("-")
      ),
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
        colorUI(ns("color"), "Choose 1st color", "#273F52"),
        colorUI(ns("color2"), "Choose 2nd color", "#D53031"),
         hr(),
        
        #dropdown menu containing download buttons
        dropdownButton(
          inputId = ns("download_menu"),
          label = "Download",
          icon = icon("sliders"),
          status = "primary",
          circle = FALSE,
          downloadButton(
            ns("downloadDESeq"),
            label =
              "DEG Table"
          ),
          downloadButton(
            ns("downloadDEVol"),
            label = 
              "Volcano"
          ),
          downloadButton(
            ns("downloadDEMA"),
            label = 
              "MA"
          ),
          downloadButton(
            ns("downloadDEHM"),
            label =
              "Heatmap"
          )
        )
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
          uiOutput(ns("reactiveText")),
          
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
            plotlyOutput(ns("ht")),
            textOutput(ns("htwarn"))
        )
      )
    )
  )
}

DE_Server <- function(id, data_species, dataset_dds, dataset_choice, reset_trigger, vst, vst_hm) {
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
      print("levelsGSEA:")
      print(levels)
      print("eval model:")
      model_term <- as.formula(paste("~", model))
      print(model_term)
      ddsTxi_dds <- DESeqDataSetFromMatrix(dds_counts, colData = meta, design = model_term)
      dds.wald <- DESeq(ddsTxi_dds, test = "Wald") 
      contrasts <- c(model, levels)
      print("contrasts for GSEA:")
      print(contrasts)
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
    print("levelsDE:")
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
  #function for creating reactive text for each dataset to let the user know which variables were chosen and what the plots are depicting
  text_generator <- function(dataset, comparison) {
    if(dataset == "Cancer_Discovery") {
      var_value <- unlist(strsplit("mono_vs_prim", "_vs_"))
    } else if(dataset == "Ye_20") {
      var_value <- unlist(strsplit("bone_marrow_vs_liver", "_vs_"))
    } else if(dataset == "Lee") {
      var_value <- unlist(strsplit("prior_cr_vs_no_prio_cr", "_vs_"))
    } else {
      var_value <- unlist(strsplit(comparison, "_vs_"))
    }
    print("var_value")
    print(var_value)
    return(var_value)
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
  
  # render reactive text to explain to the user which variables are being shown for each dataset in the plots
  output$reactiveText <- renderUI({
    data_text <-
      text_generator(dataset_choice$user_dataset(), dataset_choice$user_PW())
    text <-
      paste(
        "The DEG table and plots below compare",
        data_text[1],
        "vs",
        data_text[2],
        "samples in the",
        dataset_choice$user_dataset(),
        "dataset"
      )
    bold_text <-
      paste("<div style='font-size: 18px;'>",
            text,
            "</div>")
    HTML(bold_text)
  })

  # render DEG table ####

  observe({
    if(!is.null(dds_result())) { #only run after runDE action button has been selected
     
      dds.res <- generateRes(dataset_choice$user_dataset(), dds_result())
    
      output$results <- renderDataTable({
        if (input$DESeqtable == TRUE) {
          DT::datatable(dds.res,
                        options = list(scrollX = TRUE))
       }
    })
      
    # Download DE table ####
      output$downloadDESeq <- downloadHandler(
        filename = paste("DEGTable", '.csv', sep=''),
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
# Download Volcano plot ####
  output$downloadDEVol <- downloadHandler(
    filename = paste("DE Volcano", '.png', sep=''),
    content = function(file) {
      colors <- c(colorDE(), "grey", color2DE())
      res.vol <- generateRes(dataset_choice$user_dataset(), dds_result())
      p <- ggplot(data=res.vol, aes(x=log2FoldChange, y=-log10(padj), col = DiffExp)) + 
        geom_point() +
        theme_cowplot(font_size = 14) +
        scale_colour_manual(values = colors) +
        theme(axis.title = element_text(face = "bold"), title = element_text(face = "bold")) +
        theme(plot.background = element_rect(fill = "#FFFFFF", colour = "#FFFFFF")) +
        theme(panel.background = element_rect(fill = "#FFFFFF", colour = "#FFFFFF")) +
        ggtitle("DE Volcano Plot") +
        coord_cartesian(xlim = c(-10, 7))
      ggsave(p, file = file, device = "png", width = 8, height = 6, units = "in",dpi = 100)
    }
  )
  
 ##MA Plot ####
  
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
   # Download MA plot ####
  output$downloadDEMA <- downloadHandler(
    filename = paste("DE MA", '.png', sep=''),
    content = function(file) {
      colors <- c(colorDE(), "grey", color2DE())
      res.ma <- generateRes(dataset_choice$user_dataset(), dds_result())
      m <- ggplot(res.ma, 
                   aes(
                     x = log2(baseMean),
                     y = `log2FoldChange`,
                     col = DiffExp
                   )) +
        geom_point() +
        geom_hline(aes(yintercept = 0)) +
        scale_color_manual(values = colors) +
        theme_cowplot(font_size = 14) +
        theme(axis.title = element_text(face = "bold"), title = element_text(face = "bold")) +
        theme(plot.background = element_rect(fill = "#FFFFFF", colour = "#FFFFFF")) +
        theme(panel.background = element_rect(fill = "#FFFFFF", colour = "#FFFFFF")) +
        ylim(c(
          min(res.ma$`log2FoldChange`),
          max(res.ma$`log2FoldChange`)
        )) +
        ggtitle("DE MA Plot") +
        xlab("log2 Mean Expression") +
        ylab("Log2 Fold Change")
      
      ggsave(m, file = file, device = "png", width = 8, height = 6, units = "in",dpi = 100)
    }
  )
    
    #Heatmap ####

  output$ht <- renderPlotly({
    req(input$DESeqHeat)
    ns <- NS(id)
    
    #determine which column needed for cluster annotations based on model choice 
    if(dataset_choice$user_dataset() %in% c("Ye_16", "Venaza", "Lagadinou", "BEAT", "TCGA")){
      dds.file <- dataset_dds()
      cond_var <- dataset_choice$user_model()
      meta <- colData(dds.file)
      cond <- meta[, cond_var]
    } else {
      dds.file <- dataset_dds()
      meta <- colData(dds.file)
      cond_var <- datasets[[dataset_choice$user_dataset()]]$PCA_var
      cond <- meta[, cond_var]
    }
    
    #generate dds results table 
    res.hm <-
      generateRes(dataset_choice$user_dataset(), dds_result())
    
    # print("res.hm")
    # print(head(res.hm))
    # print(class(res.hm))
    #filter DE object for only significantly differentially expressed genes and
    #take the top 50 most highly expressed genes for visualization
    dds.mat <- res.hm %>%
      dplyr::filter(padj < 0.05 & abs(`log2FoldChange`) >= 2) %>% 
      dplyr::arrange(desc(abs(`log2FoldChange`))) %>% 
      slice(1:50)
    
    # print("dds.mat:")
    # print(dds.mat)
    
    if(nrow(dds.mat) == 0) {
      return(NULL)
    }
    #filter vst counts matrix by sig expressed genes
    vst.mat <- vst() %>%
      dplyr::filter(., ensembl_gene_id %in% dds.mat$ensembl_gene_id) %>%
      distinct(ext_gene_ensembl, .keep_all = TRUE) %>% 
      column_to_rownames(., var = "ensembl_gene_id") %>%
      dplyr::select(., -ext_gene_ensembl) %>%
      as.matrix()
    
    rownames(vst.mat) = dds.mat$Gene
    
    vst.mat <- t(scale(t(vst.mat)))
    #only show the first 100 genes for visualization in this example(can change)
   
    # print("vst.mat:")
    # print(head(vst.mat))
    #create a colorRamp function based on user input in color palette choices
    colors.hm <- c(colorDE(), "#FFFFFF", color2DE())
    
    # k_number <-
    #   if(dataset_choice$user_PW() == "LRT") {
    #     datasets[[dataset_choice$user_dataset()]]$k
    #   } else{
    #     datasets[[dataset_choice$user_dataset()]]$k_PW
    #   }
    #   
    # 
    #  k_number <- as.numeric(k_number)
    # print("k_number:")
    # print(k_number)
    #create heatmap object
    ht <- heatmaply(
      vst.mat,
      #k_col = k_number,
      row_text_angle = 45,
      height = 600,
      width = 600,
      colors = colors.hm,
      dendrogram = "column",
      show_dendrogram = TRUE,
      col_side_colors = cond,#adds labels to clusters
      showticklabels = c(FALSE, FALSE)# removes column and row labels
    )
   
    ht
 
  })

  
  output$htwarn <- renderText({
    req(input$DESeqHeat)
    
    res.hm <- generateRes(dataset_choice$user_dataset(), dds_result())
    dds.mat <- res.hm %>% 
      dplyr::filter(padj < 0.05 & abs(`log2FoldChange`) >= 2) %>% 
      dplyr::arrange(desc(abs(`log2FoldChange`))) %>% 
      slice(1:50)
    
    if(nrow(dds.mat) == 0){
      return("No significant DEGs found")
    }else {
      return(NULL)
    }
    
  })

# download DE Heatmap ####
  output$downloadDEHM <- downloadHandler(
    filename = paste("DE_Heatmap", '.png', sep=''),
    content = function(file) {
      #generate dds results table 
      res.hm <-
        generateRes(dataset_choice$user_dataset(), dds_result())
      print("res-class:")
      print(class(res.hm))
      dds.mat <- res.hm %>%
        dplyr::filter(padj < 0.05 & abs(`log2FoldChange`) >= 2) %>% 
        dplyr::arrange(desc(abs(`log2FoldChange`))) %>% 
        slice(1:50)
      print("res-dds.mat:")
      print(class(dds.mat))
      #filter vst counts matrix by sig expressed genes
      vst.mat <- vst() %>%
        dplyr::filter(., ensembl_gene_id %in% dds.mat$ensembl_gene_id) %>%
        distinct(ext_gene_ensembl, .keep_all = TRUE) %>% 
        column_to_rownames(., var = "ensembl_gene_id") %>%
        dplyr::select(., -ext_gene_ensembl) %>%
        as.matrix()
      print("vst.mat:")
      print(class(vst.mat))
      rownames(vst.mat) = dds.mat$Gene
      
      vst.mat <- t(scale(t(vst.mat)))
     
      #create a colorRamp function based on user input in color palette choices
      colors.hm <- colorRamp2(c(-2, 0, 2), c(colorDE(), "white", color2DE()))
      ht = ComplexHeatmap::Heatmap(
        vst.mat,
        name = "z scaled expression",
        col = colors.hm,
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 6),
        column_title = NULL,
        row_title = "Top DEGs"
      ) 
      png(file)
      # draw heatmap object
      draw(ht)
      dev.off()
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
