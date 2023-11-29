##Differential Expression tab##
library(colorRamp2)
library(InteractiveComplexHeatmap)
datasets <- 
  read_yaml("./data.yaml")
raw_data <- getwd()
t2g_hs <- read_rds("~/Documents/GitHub/EAGLe-2.0/data_rds_files/t2g_hs.rds")
t2g_mm <- read_rds("~/Documents/GitHub/EAGLe-2.0/data_rds_files/t2g_mm.rds")

#choose relevant comparisons for the user to choose from and have them run on the fly
DE_UI <- function(id) {
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
          selectInput(
            ns("DEmodel"),
            label = "Choose a metadata variable for DE design",
            choices = NULL
          ),
          selectInput(
            ns("pwc"),
            label = "Choose to run LRT OR a pairwise comparison",
            choices = "LRT"
          ),
          actionButton(
            (ns("runDE")),
            "Run DESeq"
          ),
          hr(),
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
        actionButton(
          (ns("sendGSEA")),
          "Send results to GSEA"
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
            InteractiveComplexHeatmapOutput(heatmap_id = ns("ht"))
        )
      )
    )
  )
}

DE_Server <- function(id, data_species, dataset_dds, dataset_choice) {
  moduleServer(id, function(input, output, session) {
    #Model Choice ####

 #    observe({
 #      shinyjs::toggle(id = "DEmodel", condition = dataset_choice() %in% c("BEAT", "TCGA"))
 #    })
 #    observe({
 #      if(dataset_choice() == "BEAT"){
 #        updateSelectInput(
 #          session = session,
 #          inputId = "DEmodel",
 #          choices = c("ven response quantile", "FAB morphology", "Denovo vs relapse") ,
 #          selected = NULL
 #        )
 #      } else if(dataset_choice() == "TCGA"){
 #        updateSelectInput(
 #          session = session,
 #          inputId = "DEmodel",
 #          choices = c("FAB morphology", "Molecular classification", "RAS mutation", "NPM1 mutation"),
 #          selected = NULL
 #        )
 #      } 
 #    })
 #    #LRT or Pairwise choices ####
 #    observe({
 #      shinyjs::toggle(id = "pwc", condition = dataset_choice() %in% c("Ye_16", "Venaza", "Lagadinou", "BEAT","TCGA"))
 #    })
 #    observe({
 #      if(dataset_choice() == "Ye_16"){
 #        updateSelectInput(
 #          session = session,
 #          inputId = "pwc",
 #          choices = c("LRT", "blood_vs_bm", "gat_vs_bm", "normBM_vs_bm", "spleen_vs_bm") ,
 #          selected = NULL
 #        )
 #      } else if(dataset_choice() == "Venaza"){
 #        updateSelectInput(
 #          session = session,
 #          inputId = "pwc",
 #          choices = c("LRT", "24hr_vs_control", "6hr_vs_control"),
 #          selected = NULL
 #        )
 #      } else if(dataset_choice() == "Lagadinou") {
 #        updateSelectInput(
 #          session = session,
 #          inputId = "pwc", 
 #          choices = c("LRT","high_PTL_5uM_vs_high_no_drug", "low_no_drug_vs_high_no_drug", "low_PTL_5uM_vs_high_no_drug")
 #        )
 #      } else if(dataset_choice() == "BEAT" & input$DEmodel == "ven response quantile") {
 #        updateSelectInput(
 #          session = session,
 #          inputId = "pwc",
 #          choices = c("quantile_q2_vs_q1", "quantile_q3_vs_q1", "quantile_q4_vs_q1"),
 #          selected = NULL
 #        )
 #      } else if(dataset_choice() == "BEAT" & input$DEmodel == "FAB morphology") {
 #        updateSelectInput(
 #          session = session,
 #          inputId = "pwc",
 #          choices = c("M0_vs_M5", "M1_vs_M5", "M3_vs_M5", "M4_vs_M5", "M5b_vs_M5"),
 #          selected = NULL
 #        )
 #      } else if(dataset_choice() == "BEAT" & input$DEmodel == "Denovo vs relapse") {
 #        updateSelectInput(
 #          session = session,
 #          inputId = "pwc",
 #          choices = "LRT",
 #          selected = NULL
 #        )
 #      } else if(dataset_choice() == "TCGA" & input$DEmodel == "FAB morphology") {
 #        updateSelectInput(
 #          session = session,
 #          inputId = "pwc",
 #          choices = c("M0_vs_M5", "M1_vs_M5", "M2_vs_M5", "M3_vs_M5", "M4_vs_M5"),
 #          selected = NULL
 #        )
 #      } else if(dataset_choice() == "TCGA" & input$DEmodel == "Molecular classification") {
 #        updateSelectInput(
 #          session = session,
 #          inputId = "pwc",
 #          choices = c("Normal vs BCR-ABL1", "Normal vs Complex", "Normal vs RUNX1-RUNX1T1", "Normal vs MLL translocation",
 #                      "Normal vs Poor Risk Cytogenetic Abnormality", "Normal vs PML-RARA", "Normal vs NUP98 translocation",
 #                      "Normal vs CBFB-MYH11", "Normal vs Intermediate Risk Cytogenetic Abnormality"),
 #          selected = NULL
 #        )
 #      } else if(dataset_choice() == "TCGA" & input$DEmodel == "RAS mutation") {
 #        updateSelectInput(
 #          session = session,
 #          inputId = "pwc",
 #          choices = "LRT",
 #          selected = NULL
 #        )
 #      } else if(dataset_choice() == "TCGA" & input$DEmodel == "NPM1 mutation") {
 #        updateSelectInput(
 #          session = session,
 #          inputId = "pwc",
 #          choices = "LRT",
 #          selected = NULL
 #        )
 #      }
 #    })
 #  #Run action button ####
 #    observe({
 #      shinyjs::toggle(id = "runDE", condition = dataset_choice() %in% c("Ye_16", "Venaza", "Lagadinou", "BEAT","TCGA"))
 #    })
 #  # run DE for pairwise data
 #    #extract the counts
 #    reactive({
 #      if(dataset_choice() %in% c("Ye_16", "Venaza", "Lagadinou","BEAT", "TCGA")) {
 #      dds_counts <- counts(dataset_dds())
 #      #extract the metadata
 #      meta <- colData(dataset_dds())
 #      #creat DESeq object from dds
 #      ddsTxi_dds <- DESeqDataSetFromMatrix(dds_counts, colData = meta, design = ~ input$DEmodel)
 #      #run wald
 #      dds.wald <- DESeq(ddsTxi_dds, test="Wald")
 #      res <- results(dds.wald, contrast=c(input$DEmodel,"M5b", "M5"))
 #      }
 #    })
 #  
 #    
 # #DE Table ####
 # 
 # ## for BEAT dataset, the ensembl id's need to be modified to work:
 # # dds.res1$ensembl_gene_id <- str_sub(dds.res1$ensembl_gene_id, end=-4) 
 # 
 # #reactive to creat dds.res object for table based on the datasets,
 #    #accounting for species differences and DE model/pairwise choices
 #    dds.res <- reactive({
 #      if(dataset_choice() %in% c("Cancer_Discovery", "Lee")) {
 # 
 #        res <- data.frame(results(dataset_dds())) %>%
 #          rownames_to_column(., var = 'ensembl_gene_id') %>%
 #          dplyr::select(., ensembl_gene_id, baseMean, log2FoldChange, padj) %>%
 #          left_join(unique(dplyr::select(t2g_hs, c(ensembl_gene_id, ext_gene))), ., by = 'ensembl_gene_id') %>%
 #          dplyr::rename(., Gene = ext_gene) %>%
 #          mutate(., DiffExp = ifelse(padj < 0.05 & log2FoldChange >= 0.5, 'up',
 #                                     ifelse(padj < 0.05 & log2FoldChange <= -0.5, 'down', 'no'))) %>%
 #          na.omit(.)
 # 
 #      } else if(dataset_choice() %in% "Ye_20") {
 # 
 #        res <- data.frame(results(dataset_dds())) %>%
 #          rownames_to_column(., var = 'ensembl_gene_id') %>%
 #          dplyr::select(., ensembl_gene_id, baseMean, log2FoldChange, padj) %>%
 #          left_join(unique(dplyr::select(t2g_mm, c(ensembl_gene_id, ext_gene))), ., by = 'ensembl_gene_id') %>%
 #          dplyr::rename(., Gene = ext_gene) %>%
 #          mutate(., DiffExp = ifelse(padj < 0.05 & log2FoldChange >= 0.5, 'up',
 #                                     ifelse(padj < 0.05 & log2FoldChange <= -0.5, 'down', 'no'))) %>%
 #          na.omit(.)
 #      } else if(dataset_choice() %in% c("Venaza", "Lagadinou") & input$DEmodel == "LRT") {
 #        res <- data.frame(results(dataset_dds())) %>%
 #          rownames_to_column(., var = 'ensembl_gene_id') %>%
 #          dplyr::select(., ensembl_gene_id, baseMean, log2FoldChange, padj) %>%
 #          left_join(unique(dplyr::select(t2g_hs, c(ensembl_gene_id, ext_gene))), ., by = 'ensembl_gene_id') %>%
 #          dplyr::rename(., Gene = ext_gene) %>%
 #          mutate(., DiffExp = ifelse(padj < 0.05 & log2FoldChange >= 0.5, 'up',
 #                                     ifelse(padj < 0.05 & log2FoldChange <= -0.5, 'down', 'no'))) %>%
 #          na.omit(.)
 #      } else if(dataset_choice() %in% c("Venaza", "Lagadinou")) {
 # 
 #        dds_counts <- counts(dataset_dds())
 #        #extract the metadata
 #        meta <- colData(dataset_dds())
 #        #creat DESeq object from dds
 #        ddsTxi_dds <- DESeqDataSetFromMatrix(dds_counts, colData = meta, design = ~ input$DEmodel)
 #        #run wald
 #        dds.wald <- DESeq(ddsTxi_dds, test="Wald")
 # 
 #        res <- data.frame(results(dds.wald, contrast=c("input$DEmodel","a","b"))) %>%
 #          rownames_to_column(., var = 'ensembl_gene_id')
 # 
 #        res <- res %>%
 #          left_join(unique(dplyr::select(t2g_hs, c(ensembl_gene_id, ext_gene))), ., by = 'ensembl_gene_id') %>%
 #          dplyr::rename(., Gene = ext_gene) %>%
 #          mutate(., DiffExp = ifelse(padj < 0.05 & log2FoldChange >= 0.5, 'up',
 #                                     ifelse(padj < 0.05 & log2FoldChange <= -0.5, 'down', 'no'))) %>%
 #          na.omit(.)
 #      } else if(dataset_choice() %in% c("BEAT", "TCGA")) {
 #        dds_counts <- counts(dataset_dds())
 #        #extract the metadata
 #        meta <- colData(dataset_dds())
 #        #creat DESeq object from dds
 #        ddsTxi_dds <- DESeqDataSetFromMatrix(dds_counts, colData = meta, design = ~ input$DEmodel)
 #        #run wald
 #        dds.wald <- DESeq(ddsTxi_dds, test="Wald")
 # 
 #        res <- data.frame(results(dds.wald, contrast=c("input$DEmodel","a","b"))) %>%
 #          rownames_to_column(., var = 'ensembl_gene_id')
 # 
 #        res$ensembl_gene_id <- str_sub(res$ensembl_gene_id, end=-4)
 # 
 #        res <- res %>%
 #          left_join(unique(dplyr::select(t2g_hs, c(ensembl_gene_id, ext_gene))), ., by = 'ensembl_gene_id') %>%
 #          dplyr::rename(., Gene = ext_gene) %>%
 #          mutate(., DiffExp = ifelse(padj < 0.05 & log2FoldChange >= 0.5, 'up',
 #                                     ifelse(padj < 0.05 & log2FoldChange <= -0.5, 'down', 'no'))) %>%
 #          na.omit(.)
 #      } else if(dataset_choice() %in% "Ye_16" & input$DEmodel == "LRT") {
 #        res <- data.frame(results(dataset_dds())) %>%
 #          rownames_to_column(., var = 'ensembl_gene_id') %>%
 #          dplyr::select(., ensembl_gene_id, baseMean, log2FoldChange, padj) %>%
 #          left_join(unique(dplyr::select(t2g_mm, c(ensembl_gene_id, ext_gene))), ., by = 'ensembl_gene_id') %>%
 #          dplyr::rename(., Gene = ext_gene) %>%
 #          mutate(., DiffExp = ifelse(padj < 0.05 & log2FoldChange >= 0.5, 'up',
 #                                     ifelse(padj < 0.05 & log2FoldChange <= -0.5, 'down', 'no'))) %>%
 #          na.omit(.)
 #      } else if(dataset_choice() %in% "Ye_16") {
 # 
 #        dds_counts <- counts(dataset_dds())
 #        #extract the metadata
 #        meta <- colData(dataset_dds())
 #        #creat DESeq object from dds
 #        ddsTxi_dds <- DESeqDataSetFromMatrix(dds_counts, colData = meta, design = ~ input$DEmodel)
 #        #run wald
 #        dds.wald <- DESeq(ddsTxi_dds, test="Wald")
 # 
 #        res <- data.frame(results(dds.wald, contrast=c("input$DEmodel","a","b"))) %>%
 #          rownames_to_column(., var = 'ensembl_gene_id')
 # 
 #        res <- res %>%
 #          left_join(unique(dplyr::select(t2g_mm, c(ensembl_gene_id, ext_gene))), ., by = 'ensembl_gene_id') %>%
 #          dplyr::rename(., Gene = ext_gene) %>%
 #          mutate(., DiffExp = ifelse(padj < 0.05 & log2FoldChange >= 0.5, 'up',
 #                                     ifelse(padj < 0.05 & log2FoldChange <= -0.5, 'down', 'no'))) %>%
 #          na.omit(.)
 #      }
 #      res
 #    })


  DEModelChoices <- function(session, dataset) {
    choices <- switch(dataset,
      "Cancer_Discovery" = "LRT",
      "Ye_16" = "Source",
      "Ye_20" = "LRT",
      "Venaza" = "Treatment",
      "Lagadinou" = "Treatment",
      "BEAT" = c("ven response quantile", "FAB morphology", "Denovo vs relapse"),
      "TCGA" = c("FAB morphology", "Molecular classification", "RAS mutation", "NPM1 mutation"),
      "Lee" = "LRT",
       default = character(0)
    )
    
    updateSelectizeInput(
      session,
      "DEmodel",
      choices = choices,
      selected = NULL
      )
  }
    
  PWChoices <- function(session, dataset, model) {
    choices <- switch(
      paste(dataset, model, sep = "_"),
      "Cancer_Discovery_LRT" = "LRT",
      "Ye_16_Source" = c("LRT", "blood_vs_bone_marrow", "gonadal_adipose_tissue_vs_bone_marrow", "normal_bm_vs_bone_marrow", "spleen_vs_bone_marrow"),
      "Ye_20_LRT" = "LRT",
      "Venaza_Treatment" = c("LRT", "24hr_vs_control", "6hr_vs_control"),
      "Lagadinou_Treatment" = c("LRT","high_PTL_5uM_vs_high_no_drug", "low_no_drug_vs_high_no_drug", "low_PTL_5uM_vs_high_no_drug"),
      "BEAT_ven response quantile" = c("quantile_q2_vs_q1", "quantile_q3_vs_q1", "quantile_q4_vs_q1"),
      "BEAT_FAB morphology" = c("M0_vs_M5", "M1_vs_M5", "M3_vs_M5", "M4_vs_M5", "M5b_vs_M5"),
      "BEAT_Denovo vs relapse" = "LRT",
      "TCGA_FAB morphology" = c("M0_vs_M5", "M1_vs_M5", "M2_vs_M5", "M3_vs_M5", "M4_vs_M5"),
      "TCGA_Molecular classification" = c("Normal vs BCR-ABL1", "Normal vs Complex", "Normal vs RUNX1-RUNX1T1", "Normal vs MLL translocation",
                                                                "Normal vs Poor Risk Cytogenetic Abnormality", "Normal vs PML-RARA", "Normal vs NUP98 translocation",
                                                                "Normal vs CBFB-MYH11", "Normal vs Intermediate Risk Cytogenetic Abnormality"),
      "TCGA_RAS mutation" = "LRT",
      "TCGA_NPM1 mutation" = "LRT",
      "Lee_LRT" = "LRT",
      default = character(0)
    )
    updateSelectInput(
      session,
      "pwc",
      choices = choices,
      selected = NULL
    )
  }
  
  runDETest <- function(dds, model, comparison) {
    # return dds if LRT is chosen
    if(comparison == "LRT") {
      return(results(dds))
    } else {
    #extract counts and metadata from preloaded dds object
    dds_counts <- counts(dds)
    meta <- colData(dds)
    print("colnames colData:")
    print(colnames(meta))
    #extract individual levels from the comparison choice
    levels <- unlist(strsplit(comparison, "_vs_"))
    print("levels:")
    print(levels)
    print("eval model:")
    print(eval(expr(model)))
    ddsTxi_dds <- DESeqDataSetFromMatrix(dds_counts, colData = meta, design = ~ eval(expr(model)))
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
    
    if(is_hs) {
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
    } else {
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
    
  # observe({
  #   shinyjs::toggle(id = "DEmodel", condition = dataset_choice() %in% c("Ye_16", "Venaza", "Lagadinou", "BEAT", "TCGA"))
  # })

  
  observe({
    DEModelChoices(session, dataset_choice())
  })
  

  # observe({
  #   shinyjs::toggle(id = "pwc", condition = dataset_choice() %in% c("Ye_16", "Venaza", "Lagadinou", "BEAT", "TCGA"))
  # })
  # 
  observe({
    PWChoices(session, dataset_choice(), input$DEmodel)
  })
  
  # observe({
  #   shinyjs::toggle(id = "runDE", condition = dataset_choice() %in% c("Ye_16","Venaza", "Lagadinou", "BEAT", "TCGA"))
  # })
  # 
  dds_result <- reactiveVal(NULL)

  #have DE run when runDE button is clicked
  observeEvent(input$runDE, {
    dds_result(runDETest(dataset_dds(), input$DEmodel, input$pwc))
  })

  observe({
    if(!is.null(dds_result())) {
      dds.res <- generateRes(dataset_choice(), dds_result())

      output$results <- renderDataTable({
        if (input$DESeqtable == TRUE) {
          dds.res
        }
      })
    }
  })
# problem: DE object can be passed from function as input for another function
# Warning: Error in as.vector: no method for coercing this S4 class to a vector

  
  #create objects for color palettes from the palette module
  colorDE <-  
    colorServer("color")
  
  color2DE <-
    colorServer("color2")
  #Volcano Plot ####
  observe({
    if(!is.null(dds_result())) {
      res.vol <- generateRes(dataset_choice(), dds_result())
  
  output$volplot <- 
    renderGirafe({
      #colors <- c(colorDE(), "grey",color2DE()) #object for colors on volcano based on user input
      colors <- c(colorDE(), "grey", color2DE())
      if(input$DESeqvolcano == TRUE) { #only create plot if the  volcano switch is toggled
        p<- ggplot(res.vol, aes( 
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
    }
  })
  #MA Plot ####
  color3DE <-
    colorServer("color3")
  
  color4DE <-
    colorServer("color4")
  
  observe({
    if(!is.null(dds_result())) {
      res.ma <- generateRes(dataset_choice(), dds_result())  
      
  output$MAplot <- 
    renderGirafe ({
     
      colors <- c(color3DE(), "grey", color4DE())#object for colors on volcano based on user input called from palette module
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
          theme_light() +
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
  #interactive heatmap needs to be wrapped in a reactive function to work
  observe({
    req(input$DESeqHeat)
    ns <- NS(id)
    
    color5DE <-
      colorServer("color5")
    
    color6DE <-
      colorServer("color6")
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

