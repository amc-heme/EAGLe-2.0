
PW_UI <- function(id) {
  ns <- NS(id)
  fluidPage(
    theme =
      shinytheme("flatly"),
    useShinyjs(),
    titlePanel(
      "Differential Expression Tables and Plots"
    ),#end title
    sidebarLayout(
      sidebarPanel(
        tagList(
    
        #   shinyWidgets::checkboxGroupButtons(
        #     inputId = (ns("pwc")),
        #     label = "Choose a Pairwise Comparison: LSC datasets",
        #     choices = c(1,2),
        #     selected = NULL,
        #     status = "primary"
        # 
        # ),
        #for BEAT and TCGA:
        
        selectInput(
            ns("DEmodel"),
            label = "Choose a metadata variable for DE design",
            choices = NULL
          ),
        selectInput(
          ns("pwc"),
          label = "Choose a Pairwise Comparison",
          choices = NULL
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
    observe({
      shinyjs::toggle(id = "pwc", condition = dataset_choice() %in% c("Ye_16", "Venaza", "Lagadinou", "BEAT","TCGA"))
    })
    observe({
      shinyjs::toggle(id = "DEmodel", condition = dataset_choice() %in% c("BEAT", "TCGA"))
    })
    #all of the pairwise tests need to run on the fly from dds.wald
  # ye_16: 
    #[1] "Intercept"                                   
    # [2] "Source_blood_vs_bone_marrow"                 
    # [3] "Source_gonadal_adipose_tissue_vs_bone_marrow"
    # [4] "Source_normal_bm_vs_bone_marrow"             
    # [5] "Source_spleen_vs_bone_marrow"    
  # venaza
    # [1] "Intercept"                 "condition_24hr_vs_control"
    # [3] "condition_6hr_vs_control"  "patient_FJ_vs_1020"       
    # [5] "patient_RB_vs_1020"  
  # Lagadinou
    # [1] "Intercept"                             
    # [2] "Treatment_high_PTL_5uM_vs_high_no_drug"
    # [3] "Treatment_low_no_drug_vs_high_no_drug" 
    # [4] "Treatment_low_PTL_5uM_vs_high_no_drug" 
    # [5] "Batch"    
  # BEAT: resultsNames(dds.wald1) (the large datasets need to be run completely on the fly from 
    #ddstxi in order to extract the correct metadata)
    # [1] "Intercept"         "quantile_q2_vs_q1" "quantile_q3_vs_q1"
    # [4] "quantile_q4_vs_q1"
    
    # create a data module for those, and do the contrast here
    ## update pairwise choice for each dataset
    # observe({
    #   if(dataset_choice() == "Ye_16"){
    #     updateCheckboxGroupButtons(
    #       session = session,
    #       inputId = "pwc",
    #       choices = c("blood_vs_bm", "gat_vs_bm", "normBM_vs_bm", "spleen_vs_bm") ,
    #       selected = NULL
    #     )
    #   } else if(dataset_choice() == "Venaza"){
    #     updateCheckboxGroupButtons(
    #       session = session,
    #       inputId = "pwc",
    #       choices = c("24hr_vs_control", "6hr_vs_control"),
    #       selected = NULL
    #     )
    #   } else if(dataset_choice() == "Lagadinou") {
    #     updateCheckboxGroupButtons(
    #     session = session,
    #     inputId = "pwc", 
    #     choices = c("high_PTL_5uM_vs_high_no_drug", "low_no_drug_vs_high_no_drug", "low_PTL_5uM_vs_high_no_drug")
    # )
    #   }
    # })
    observe({
      if(dataset_choice() == "BEAT"){
        updateSelectInput(
          session = session,
          inputId = "DEmodel",
          choices = c("ven response quantile", "FAB morphology", "Denovo vs relapse") ,
          selected = NULL
        )
      } else if(dataset_choice() == "TCGA"){
        updateSelectInput(
          session = session,
          inputId = "DEmodel",
          choices = c("FAB morphology", "Molecular classification", "RAS mutation", "NPM1 mutation"),
          selected = NULL
        )
      } 
    })
    observe({
      if(dataset_choice() == "Ye_16"){
        updateSelectInput(
          session = session,
          inputId = "pwc",
          choices = c("blood_vs_bm", "gat_vs_bm", "normBM_vs_bm", "spleen_vs_bm") ,
          selected = NULL
        )
      } else if(dataset_choice() == "Venaza"){
        updateSelectInput(
          session = session,
          inputId = "pwc",
          choices = c("24hr_vs_control", "6hr_vs_control"),
          selected = NULL
        )
      } else if(dataset_choice() == "Lagadinou") {
        updateSelectInput(
          session = session,
          inputId = "pwc", 
          choices = c("high_PTL_5uM_vs_high_no_drug", "low_no_drug_vs_high_no_drug", "low_PTL_5uM_vs_high_no_drug"),
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






















