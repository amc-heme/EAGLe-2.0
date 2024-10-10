## GSEA tab module
datasets <- 
  read_yaml("./data.yaml")
raw_data <- getwd()
t2g_hs <- read_rds(paste0(raw_data, "/data/t2g_hs.rds"))
t2g_mm <- read_rds(paste0(raw_data, "/data/t2g_mm.rds"))
#GSEA data table

#load molecular pathways for GSEA
pathways.hallmark <- gmtPathways("gmt_pathway_files/h.all.v7.4.symbols.gmt")
pathways.GOall <- gmtPathways("gmt_pathway_files/c5.go.v2022.1.Hs.symbols.gmt")
pathways.GOmolec <- gmtPathways("gmt_pathway_files/c5.go.mf.v7.4.symbols.gmt")
pathways.GOcellcomp <- gmtPathways("gmt_pathway_files/c5.go.cc.v2022.1.Hs.symbols.gmt") 
pathways.GObio <- gmtPathways("gmt_pathway_files/c5.go.bp.v7.4.symbols.gmt")
pathways.TFtargets <-gmtPathways("gmt_pathway_files/c3.tft.v2022.1.Hs.symbols.gmt")
pathways.allReg <- gmtPathways("gmt_pathway_files/c3.all.v2022.1.Hs.symbols.gmt")
pathways.Wiki <- gmtPathways("gmt_pathway_files/c2.cp.wikipathways.v2022.1.Hs.symbols.gmt")
pathways.Reactome <-gmtPathways("gmt_pathway_files/c2.cp.reactome.v2022.1.Hs.symbols.gmt")
pathways.KEGG <- gmtPathways("gmt_pathway_files/c2.cp.kegg.v2022.1.Hs.symbols.gmt")
pathways.Positional <-gmtPathways("gmt_pathway_files/c1.all.v2022.1.Hs.symbols.gmt")
pathways.Biocarta <-gmtPathways("gmt_pathway_files/c2.cp.biocarta.v2022.1.Hs.symbols.gmt")
pathways.lsc <- gmtPathways("gmt_pathway_files/lsc_sigs.gmt")
pathways.aeg <- gmtPathways("gmt_pathway_files/aeg_genesets_20220602.gmt")

GSEA_UI <- function(id) {
  ns <- NS(id)
  
  fluidPage(
    theme =
      shinytheme("flatly"),
    titlePanel(
      "GSEA"
    ),#end title
    sidebarLayout(
      sidebarPanel( 
        h4("Choose gmt file to load pathway sets"),
        #dropdown menu for molecular pathways 
        selectInput(ns("filechoice"), label = NULL,
                    choices = c(Hallmark = "hallmark", GOall = "goall", GOmolecular = "GOmolec", 
                                GOcellcomp = "GOcellcomp", GObio = "GObio", TFtargets = "TFtargets",
                                allRegular = "allReg", Wiki = "wiki", Reactome = "reactome", KEGG = "KEGG",
                                Positional = "positional", Biocarta = "biocarta", lsc = "lsc", aeg = "aeg")),
        h4("Choose a plot type"),
        #toggle buttons to load table and plots
        materialSwitch(
          inputId =
            ns("fgseaTable"),
          label =
            "Table",
          value =
            FALSE,
          right =
            TRUE
        ),
        materialSwitch(
          inputId =
            ns("rankedplot"),
          label =
            "Waterfall",
          value =
            FALSE,
          right =
            TRUE
        ),
        materialSwitch(
          inputId =
            ns("moustache"),
          label =
            "Moustache",
          value =
            FALSE,
          right =
            TRUE
        ),
        materialSwitch(
          inputId =
            ns("eplot"),
          label =
            "Enrichment",
          value =
            FALSE,
          right =
            TRUE
        ),
        materialSwitch(
          inputId =
            ns("volcanoplot"),
          label =
            "Volcano",
          value =
            FALSE,
          right =
            TRUE
        ),
        materialSwitch(
          inputId =
            ns("heatmap"),
          label =
            "Heatmap",
          value =
            FALSE,
          right =
            TRUE
        ),
        #color palette choices for all plots
        colorUI(ns("color7"), "Choose a color for plots", "#D53031"),
    
        
        hr(),
        #GSEA Waterfall####
        conditionalPanel(
          ns = ns,
          condition = "input.rankedplot == true",
         
          h4("Waterfall Plot Specific Options"),
      
  
          #slider scale to choose how many pathways to load
          sliderInput(ns("howmanypathways"), "Choose How Many Pathways to Rank",
                      min = 5, max = 50, value = 15
          )
        ),

         hr(),
        #GSEA Enrichment plot ####
        conditionalPanel(
          ns = ns,
          condition = "input.eplot == true",
       
          h4("Enrichment Plot Specific Options"),
        
          #options for type of enrichment plot to load
          radioButtons(ns("topupordownbutton"), label = "Enrichment Plot Choices",
                       choices = list("Top Ranked Up Pathway" = "topup",
                                      "Top Ranked Down Pathway" = "topdown", 
                                      "Pathway of Choice:" = "eplotpath"),
                       selected = "topup"),
       
          #dropdown menu for specific pathway choices that is reactive to the pathway set dropdown menu
          selectizeInput(
            ns("pathwaylisteplot"),
            label=
              "Choose a specific pathway",
            choices =
              NULL,
            selected = NULL ,
            options = list(maxItems = 1)
          )
        ),
      hr(),
        #GSEA Volcano ####
        conditionalPanel(
          ns = ns,
          condition = "input.volcanoplot == true",
 
          h4("Volcano Plot Specific Options"),
   
          #dropdown list of specific pathways for volcano plot that is reactive to the pathway set dropdown menu
          selectizeInput(
            ns("pathwaylist"),
            label=
              "Choose a specific pathway(s) to view on volcano plot",
            choices =
              NULL,
            selected = NULL,
            options = list(maxItems = 1)
          ),
          tags$style(
            type = 'text/css',
            ".selectize-input { word-wrap : break-word;}
             .selectize-input { word-break: break-word;}
             .selectize-dropdown {word-wrap : break-word;} "
          )
        ),
    
        #GSEA Heatmap ####
        conditionalPanel(
          ns = ns,
          condition = "input.heatmap == true",
          h4("Heatmap Specific Options"),

          #color palette choices for heatmap
          colorUI(ns("colorHM1"), "Choose 1st color","#273F52"),

          colorUI(ns("colorHM2"), "Choose 2nd color", "#D53031"),
          #dropdown menu of specific pathways for heatmap, reactive to pathway sets dropdown
          selectizeInput(
            ns("pathwaylistht"),
            label=
              "Choose a specific pathway to view genes on heatmap",
            choices =
              NULL,
            selected = NULL ,
            options = list(maxItems = 1)
          )
        ),
      #dropdown menu containing download buttons ####
      dropdownButton(
        inputId = ns("download_menu"),
        label = "Download",
        icon = icon("sliders"),
        status = "primary",
        circle = FALSE,
        downloadButton(
          ns("downloadfgsea"),
          label =
            "Pathway List"
        ),
        downloadButton(
          ns("downloadranks"),
          label =
            "Waterfall"
        ),
        downloadButton(
          ns("downloadmoustache"),
          label = 
            "Moustache"
        ),
        downloadButton(
          ns("downloadeplot"),
          label =
            "Enrichment"
        ),
        downloadButton(
          ns("downloadvolcano"),
          label =
            "Volcano"
        ),
        downloadButton(
          ns("downloadheatmap"),
          label = 
            "Heatmap"
        )
      )
      ),

      mainPanel(
        uiOutput(ns("reactiveTextGSEA")),
        conditionalPanel(
          ns = ns,
          condition = "input.fgseaTable == true",
          shinycssloaders::withSpinner(
          DTOutput( #add loading spinners
            ns("fgseaTable")
          )
          )
        ),
        conditionalPanel(
          ns = ns,
          condition = "input.rankedplot == true",
          shinycssloaders::withSpinner( #add loading spinners
            plotOutput(
              ns("GSEAranked"),
              height = "100%"
            )
          )
        ),
        conditionalPanel(
          ns = ns,
          condition = "input.moustache == true",
          shinycssloaders::withSpinner( #add loading spinners
            girafeOutput(
              ns("GSEAMoustache")
            )
          )
        ),
        conditionalPanel(
          ns = ns,
          condition = "input.eplot == true",
          shinycssloaders::withSpinner( #add loading spinners
            plotOutput(
              ns("GSEAenrichment")
            )
          )
        ),
        conditionalPanel(
          ns = ns,
          condition = "input.volcanoplot == true",
          shinycssloaders::withSpinner( #add loading spinners
            girafeOutput(
              ns("GSEAvolcano")
            )
          )
        ),
        conditionalPanel(
          ns = ns,
          condition = "input.heatmap == true",
          shinycssloaders::withSpinner(
          plotlyOutput(ns("htgsea"))
          ),
          uiOutput(ns("htwarnGSEA"))
        )
      )
    )
  )
}


GSEA_Server <- function(id, dataset_choice, DE_res, reset_trigger, vst, dataset_dds) {
  moduleServer(id, function(input, output, session) {
 
    #run GSEA for chosen pathway input
   
    #make an object to hold the values of the selectInput for gsea pathway choices
    gsea_file_values <- list("hallmark" = pathways.hallmark,
                             "goall" = pathways.GOall,
                             "GOmolec" = pathways.GOmolec, 
                             "GOcellcomp" = pathways.GOcellcomp,
                             "GObio" = pathways.GObio,
                             "TFtargets" = pathways.TFtargets,
                             "allReg" = pathways.allReg,
                             "wiki" = pathways.Wiki,
                             "reactome" = pathways.Reactome,
                             "KEGG" = pathways.KEGG,
                             "positional" = pathways.Positional,
                             "biocarta" = pathways.Biocarta,
                             "lsc" = pathways.lsc,
                             "aeg" = pathways.aeg)
    #function for using correct t2g
    generateEnsGene <- function(dataset) {
      #mouse or human?
      is_hs <- grepl("t2g_hs", datasets[[dataset]]$t2g)
      
    if(is_hs) {
      ens2gene <- t2g_hs[, c(2,3)]
      } else{
        ens2gene <- t2g_mm[, c(2,3)]
      }
      return(ens2gene)
    }
    #ens2gene object based on user specified dataset choice
    ens2gene <- reactive({
      generateEnsGene(dataset_choice$user_dataset())
    })
    # function for calculating ranks
    pathway_ranks <- function(DE_results, ensgene) {
     
      #add the human name of the gene to the last column
      result <- inner_join(DE_results, ensgene, by = c("row" = "ensembl_gene_id"))
      colnames(result)[8] <- 'HS_Symbol'
      # select only the human gene symbol and the 'stat' from the results,
      result <- result %>% 
        dplyr::select(HS_Symbol, stat) %>% 
        na.omit() %>% 
        distinct() %>% 
        group_by(HS_Symbol) %>% 
        summarize(stat = mean(stat))
      
      #reconfigure the data and return ranks
      ranks <- deframe(result)
      return(sort(ranks))
    }
    
    ranks <- reactive({
      pathway_ranks(DE_res$res_tidy(), ens2gene())
    })
    
    #function for creating reactive text for each dataset to let the user know which variables were chosen and what the plots are depicting
    text_generator <- function(dataset, comparison) {
      if(dataset == "Cancer_Discovery") {
        var_value <- unlist(strsplit("primitive_vs_monocytic", "_vs_"))
      } else if(dataset == "Ye_20") {
        var_value <- unlist(strsplit("liver_vs_bone_marrow", "_vs_"))
      } else if(dataset == "Lee") {
        var_value <- unlist(strsplit("prior complete remission_vs_no_prior complete remission", "_vs_"))
      } else {
        var_value <- unlist(strsplit(comparison, "_vs_"))
      }
      return(var_value)
    }
    # render reactive text to explain to the user which variables are being shown for each dataset in the plots
    output$reactiveTextGSEA <- renderUI({
      data_text <-
        text_generator(dataset_choice$user_dataset(), dataset_choice$user_PW())
      text <-
        paste(
          "The table and plots below show gene set enrichment analysis results
          from the comparison of",
          data_text[1],
          "vs",
          data_text[2],
          "samples in the",
          dataset_choice$user_dataset(),
          "dataset.  Hover cursor over points on the plots for gene or pathway names.")
        text2 <-  paste(
          data_text[1],
          "= positive NES",
          ",",
          data_text[2],
          "= negative NES."
        )
      large_text_style <-
        paste("<div style='font-size: 18px;'>",
              text,
              "</div><br><div>",
              "<div style='font-size: 18px;'>",
              text2,
              "</div>")
      HTML(large_text_style)
    })
    
    fgseaRes <- reactive({
      pathwaygsea <- gsea_file_values[[input$filechoice]]
      fgsea::fgsea(pathways = pathwaygsea,
                   stats = ranks(),
                   nproc = 10)
    })
    #reactive expression to run fgsea and load results table for each chosen pathway
    gseafile <-
      eventReactive(input$filechoice,{
        # pathwaygsea <- gsea_file_values[[input$filechoice]]
        # fgseaRes <- fgsea::fgsea(pathways = pathwaygsea, stats = ranks(), nproc = 10)
        fgseaResTidy <- fgseaRes() %>%
          as_tibble() %>%
          dplyr::select(., -pval,-log2err, -ES) %>% 
          arrange(desc(NES))
        fgseaResTidy
      })
    
    #download pathway list ####
    output$downloadfgsea <- downloadHandler(
      filename = paste("GSEATable", '.csv', sep = ''),
      content = function(file) {
        fwrite(gseafile(),file, sep=",", sep2=c("", " ", ""))
      }
    )
    
      
    #filter Res table for chosen pathway to show in a waterfall plot
    gseafile_waterfall <-
      reactive({
        fgseaResTidy <- fgseaRes() %>%
          as_tibble() %>%
          dplyr::select(., -pval,-log2err, -ES) %>% 
          arrange(desc(NES))
        top15 <- fgseaResTidy %>% 
          top_n(n = input$howmanypathways, wt = NES)
        bottom15 <- fgseaResTidy %>%
          top_n(n = -(input$howmanypathways), wt = NES)
        fgseaResTidy <- rbind.data.frame(top15, bottom15)
        fgseaResTidy
      })
    
    #table output for Res tables
    output$fgseaTable <- renderDataTable({
      if (input$fgseaTable == TRUE) {
        
        gsea_table <- gseafile()
        DT::datatable(gsea_table,
                     options = list(scrollX = TRUE))
      }
    })
    
    #### GSEA pathway ranks waterfall plot ####
    #call in singlecolor_palette module for plot
    colorWF <- 
      colorServer("color7")
 
    output$GSEAranked <- renderPlot(
      width = 800,
      height = 600,
      {
        if (input$rankedplot == TRUE) {
          #color object reactive to user choice from palette
          colors <- c("grey", colorWF())
          
            ggplot(gseafile_waterfall(), aes(reorder(pathway, NES), NES)) +
            geom_col(aes(fill = padj < 0.05)) +
            scale_fill_manual(values = colors) +
            coord_flip() +
            labs(
              x = "Pathway",
              y = "Normalized Enrichment Score",
              title = "Pathway NES from GSEA"
            ) +
            theme_cowplot(font_size = 14) +
            theme(axis.title = element_text(face = "bold"), title = element_text(face = "bold"))
        }
      }
    )
    
    #download button output- waterfall plot ####
    output$downloadranks <- downloadHandler(
      filename = function() { paste("Waterfall Plot", '.png', sep='') },
      content = function(file) {
        
        colors <- c("grey", colorWF())
        
        w <- ggplot(gseafile_waterfall(), aes(reorder(pathway, NES), NES)) +
          geom_col(aes(fill = padj < 0.05)) +
          scale_fill_manual(values = colors) +
          coord_flip() +
          labs(
            x = "Pathway",
            y = "Normalized Enrichment Score",
            title = "Pathway NES from GSEA"
          ) +
          theme_cowplot(font_size = 14) +
          theme(axis.title = element_text(face = "bold"), title = element_text(face = "bold")) +
          theme(plot.background = element_rect(fill = "#FFFFFF", colour = "#FFFFFF")) +
          theme(panel.background = element_rect(fill = "#FFFFFF", colour = "#FFFFFF"))  
        
        ggsave(w, file = file, device = "png", width = 8, height = 6, units = "in",dpi = 100)
      }
    )
    #reactive wrapper for js functions to hide or show plot dimensions based on toggle switch
    observe({
      toggle(id = "rankedheightslider", condition = input$hidedimsWF)
      toggle(id ="rankedwidthslider", condition = input$hidedimsWF)
    })
    #### GSEA moustache plot ####
    #reactive title for moustache plot
    gseamoustache_title <-
      eventReactive(input$filechoice, {
        print(input$filechoice)
      })
    
    #function for making moustache reactive to pathway choice
    #fgsea analysis run on pathway chosen
    #resulting table is put into tidy format and uncessary rows are removed for plotting
    #new column added to state significance 
    toplotMoustache <-
      reactive({
        # pathwaygsea <- gsea_file_values[[input$filechoice]]
        # fgseaRes <-
        #   fgsea::fgsea(pathways = pathwaygsea,
        #                stats = ranks(),
        #                nproc = 10)
        fgseaResTidy <-
          fgseaRes() %>%
          as_tibble() %>%
          dplyr::select(., -pval,-log2err, -ES) %>% 
          arrange(desc(NES))
        toplotMoustache <-
          cbind.data.frame(fgseaResTidy$pathway,
                           fgseaResTidy$NES,
                           fgseaResTidy$padj)
        colnames(toplotMoustache) <- c("pathway", "NES", "padj")
        toplotMoustache <- toplotMoustache %>%
          mutate(., sig = ifelse(padj <= 0.05, 'yes', 'no'))
        
      })
    #call in singlecolor_module for plot
    colorM <- 
      colorServer("color7")
    output$GSEAMoustache <- renderGirafe(
      {
        if (input$moustache == TRUE) {
          #color object reactive to user input from plalette choice
          colors <- c('grey', colorM())
          ma <-
            ggplot(toplotMoustache(), aes(x = NES, y = padj, color = sig, tooltip = pathway)) +
            geom_point_interactive(alpha = 0.8, size = 0.5) +
            theme_cowplot(font_size = 14) +
            theme(axis.title = element_text(face = "bold"), title = element_text(face = "bold")) +
            xlab('NES') +
            scale_color_manual(values = colors) +
            ylab('adjusted p-value') +
            ggtitle("Pathways from GSEA") +
            coord_cartesian(xlim = c(-3, 3), ylim = c(-0.1, 1)) 
          girafe(code = print(ma))
        }
      }
    )
    #download button output - moustache plot ####
    output$downloadmoustache <- downloadHandler(
      filename = paste("Moustache Plot", '.png', sep=''),
      content = function(file) {
        
        colors <- c('grey', colorM())
        mous <-
          ggplot(toplotMoustache(), aes(x = NES, y = padj, color = sig)) +
          geom_point() +
          theme_cowplot(font_size = 14) +
          theme(axis.title = element_text(face = "bold"), title = element_text(face = "bold")) +
          theme(plot.background = element_rect(fill = "#FFFFFF", colour = "#FFFFFF")) +
          theme(panel.background = element_rect(fill = "#FFFFFF", colour = "#FFFFFF")) +
          xlab('NES') +
          scale_color_manual(values = colors) +
          ylab('adjusted p-value') +
          ggtitle("Pathways from GSEA") +
          geom_text_repel(colour = "black", aes(label= ifelse(padj <0.05, as.character(pathway), ""), hjust=0,vjust=0)) +
          coord_cartesian(xlim = c(-3, 3), ylim = c(-0.1, 1)) 
        
      ggsave(mous, file = file, device = "png", width = 8, height = 6, units = "in",dpi = 100)
      }
    )
    
    #GSEA Enrichment Plots ####
    #reactive expression for specific pathway choice
    observe({
      pathwaygsea <- gsea_file_values[[input$filechoice]]
      updateSelectizeInput(session,"pathwaylisteplot", choices = names(pathwaygsea), server = TRUE)})
    
    #reactive expression for plot title based on enrichment plot type
    gseaeplot_title <-
      eventReactive(input$pathwaylisteplot, {
        print(input$pathwaylisteplot)
      })
    #fgsea run for selected pathway, tidy results table, filter for top up or down pathway to plot, or plot pathway of choice
    output$GSEAenrichment <- renderPlot ({
       pathwaygsea <- gsea_file_values[[input$filechoice]]
      # fgseaRes <-
      #   fgsea::fgsea(pathways = pathwaygsea,
      #                stats = ranks(),
      #                nproc = 10)
      fgseaResTidy <- fgseaRes() %>%
        as_tibble() %>%
        dplyr::select(., -pval,-log2err, -ES) %>% 
        arrange(desc(NES))
      if(input$topupordownbutton == "topup") {
        top.UP.path <- as.character(fgseaResTidy[1,1])
        plotEnrichment(pathwaygsea[[top.UP.path]],
                       ranks()) + labs(title=top.UP.path)
      } else if(input$topupordownbutton == "topdown") {
        top.DOWN.path <- as.character(fgseaResTidy[nrow(fgseaResTidy), 1])
        plotEnrichment(pathwaygsea[[top.DOWN.path]],
                       ranks()) + labs(title=top.DOWN.path)
      } else if(input$topupordownbutton == "eplotpath") {
        plotEnrichment(pathwaygsea[[input$pathwaylisteplot]],
                       ranks()) + labs(title= input$pathwaylisteplot)
      }
    })
    #download button output- enrichment plot ####
    output$downloadeplot <- downloadHandler(
      filename = paste("Enrichment Plot", '.png', sep=''),
      content = function(file) {
        # pathwaygsea <- gsea_file_values[[input$filechoice]]
        # fgseaRes <-
        #   fgsea::fgsea(pathways = pathwaygsea,
        #                stats = ranks(),
        #                nproc = 10)
        fgseaResTidy <- fgseaRes() %>%
          as_tibble() %>%
          dplyr::select(., -pval,-log2err, -ES) %>% 
          arrange(desc(NES))
        if(input$topupordownbutton == "topup") {
          top.UP.path <- as.character(fgseaResTidy[1,1])
          e <- plotEnrichment(pathwaygsea[[top.UP.path]],
                         ranks()) + labs(title=top.UP.path)
        } else if(input$topupordownbutton == "topdown") {
          top.DOWN.path <- as.character(fgseaResTidy[nrow(fgseaResTidy), 1])
          e <- plotEnrichment(pathwaygsea[[top.DOWN.path]],
                         ranks()) + labs(title=top.DOWN.path)
        } else if(input$topupordownbutton == "eplotpath") {
          e <- plotEnrichment(pathwaygsea[[input$pathwaylisteplot]],
                         ranks()) + labs(title= input$pathwaylisteplot)
        }
        ggsave(e, file = file, device = "png", width = 8, height = 6, units = "in",dpi = 100)
      }
    )
    
    
    # GSEA Volcano plot ####
    
    #drop down list of pathways within chosen pathway set
    observe({
      pathwaygsea1 <- gsea_file_values[[input$filechoice]]
      updateSelectizeInput(session,"pathwaylist", choices = names(pathwaygsea1), server = TRUE)
    })
   
    #create an object of the genes in the pathway and the genes included in that pathway that are also DE
    pathway.genes <- reactive({
      
      req(DE_res$dds_res())
      
      dds.res.v <- #This should probably be filtered for only DE genes
        DE_res$dds_res()
      
      pathwaygsea2 <- gsea_file_values[[input$filechoice]]
      
      p1 <-
        unlist((pathwaygsea2[names(pathwaygsea2) %in% input$pathwaylist]))

      dds.res.pathways <- dds.res.v %>%
        mutate(., genes_in_pathway = ifelse(Gene %in% p1, 'yes', 'no'))
      
      dds.res.pathways 
    })

   
 
    #call in singlecolor_module for plot
    v_col <- colorServer("color7")
    
    output$GSEAvolcano <- renderGirafe ({
      
      req(pathway.genes())
      
      #reactive title for volcano based on specific pathway choice
      gseavol_title <-
          paste(input$pathwaylist)
    
      
      colors <-
        c("grey", v_col())
     
      #if (input$volcanoplot == TRUE) {
       v <- ggplot(
          data = (pathway.genes() %>% arrange(., (genes_in_pathway))),
          aes(
            x = log2FoldChange,
            y = -log10(padj),
            col = genes_in_pathway[],
            tooltip = Gene
          )
        ) +
         theme_cowplot(font_size = 14) +
          theme(axis.title = element_text(face = "bold"), title = element_text(face = "bold")) +
          geom_point_interactive(size = 1, alpha = 0.5) +
          scale_color_manual(values = colors) +
          ggtitle(gseavol_title) + #reactive title
          xlab("log2foldchange")
        girafe(code = print(v))
      })
  
    #download button output- volcano plot ####
    output$downloadvolcano <- downloadHandler(
      filename = paste("Volcano Plot", '.png', sep=''),
      content = function(file) {
        
        gseavol_title <-
          paste(input$pathwaylist)
  
        colors <-
          c("grey", v_col())
        
        #if (input$volcanoplot == TRUE) {
        v <- ggplot(
          data = (pathway.genes() %>% arrange(., (genes_in_pathway))),
          aes(
            x = log2FoldChange,
            y = -log10(padj),
            col = genes_in_pathway[]
          )
        ) +
          theme_cowplot(font_size = 14) +
          theme(axis.title = element_text(face = "bold"), title = element_text(face = "bold")) +
          theme(plot.background = element_rect(fill = "#FFFFFF", colour = "#FFFFFF")) +
          theme(panel.background = element_rect(fill = "#FFFFFF", colour = "#FFFFFF")) +
          geom_point() +
          scale_color_manual(values = colors) +
        geom_text_repel(
          max.overlaps = 1500,
          colour = "black",
          aes( #only label is gene is in pathway and sig expression
            label = ifelse(
              genes_in_pathway == 'yes' & log2FoldChange > 1.5,
              as.character(Gene),
              ""
            )
          ),
          hjust = 0,
          vjust = 0
        ) +
        ggtitle(gseavol_title) + #reactive title
          xlab("log2foldchange")
        ggsave(v, file = file, device = "png", width = 8, height = 6, units = "in",dpi = 100)
      })

    
    
    #GSEA heatmap ####
    
    #check for NULL output
    null_check <- reactive({
      req(input$heatmap)
      req(DE_res$dds_res())
      req(input$pathwaylistht)
      
      res.hmg <- DE_res$dds_res() 
      
      #filter DE object for significance
      dds.sig <- res.hmg %>%
        dplyr::filter(padj < 0.05 & abs(`log2FoldChange`) >= 0.5)
      
      if(nrow(dds.sig) == 0) {
        return(TRUE)
      }
      
      pathwaygsea4 <- gsea_file_values[[input$filechoice]]
      
      
      p2 <- unlist((pathwaygsea4[names(pathwaygsea4) %in% input$pathwaylistht]))
      
      #filter vst counts matrix for genes in pathway
      vst.myc <- vst() %>%
        mutate(., pathwayheat = ifelse(ext_gene_ensembl %in% p2, 'yes', 'no')) %>%
        dplyr::filter(pathwayheat == "yes") %>% 
        dplyr::rename("Gene" = "ext_gene_ensembl")
      
      #create matrix for heatmap using genes in significant DE that are in pathway of choice
      vstgsea.mat <- vst.myc %>%
        dplyr::filter(., ensembl_gene_id %in% dds.sig$ensembl_gene_id) %>%
        distinct(Gene, .keep_all = TRUE) %>% 
        column_to_rownames(var = "Gene") %>%
        dplyr::select(.,-ensembl_gene_id, -pathwayheat) %>%
        as.matrix()
      
      if(nrow(vstgsea.mat) == 0) {
        return(TRUE)
      }
      
      return(FALSE)
    })
    #reactive expression for selected pathway choice for heatmap
    observe({
      req(input$filechoice)
      req(input$heatmap)
      
       pathwaygsea3 <- gsea_file_values[[input$filechoice]]
      # fgseaRes3 <-
      #   fgsea::fgsea(pathways = pathwaygsea3,
      #                stats = ranks(),
      #                nproc = 10)
      fgseaResTidy3 <- fgseaRes() %>%
        as_tibble() %>%
        dplyr::select(., -pval,-log2err, -ES) %>% 
        arrange(desc(NES))
      top.up <- as.character(fgseaResTidy3[1,1])
      updateSelectizeInput(session,
                           "pathwaylistht",
                           choices = names(pathwaygsea3),
                           selected = top.up,
                           server = TRUE)})

    #call in single color_module for plot palette
    colorGHeat <-
      colorServer("colorHM1")
    colorGHeat2 <-
      colorServer("colorHM2")
    #interactive heatmap must be wrapped in a reactive expression
    output$htgsea <- renderPlotly({
      if(null_check()) {
        return(NULL)
      }
      
      req(input$heatmap)
      req(DE_res$dds_res())
      req(input$pathwaylistht)
      
      #determine which column needed for cluster annotations based on model choice 
        dds.file <- dataset_dds()
        meta <- colData(dds.file)
        cond_var <- datasets[[dataset_choice$user_dataset()]]$PCA_var
        cond <- meta[, cond_var]
      
      
      res.hmg <- DE_res$dds_res() 

      #filter DE object for significance
      dds.sig <- res.hmg %>%
        dplyr::filter(padj < 0.05 & abs(`log2FoldChange`) >= 0.5)

      pathwaygsea4 <- gsea_file_values[[input$filechoice]]
     
      
      p2 <- unlist((pathwaygsea4[names(pathwaygsea4) %in% input$pathwaylistht]))

      #filter vst counts matrix for genes in pathway
      vst.myc <- vst() %>%
        mutate(., pathwayheat = ifelse(ext_gene_ensembl %in% p2, 'yes', 'no')) %>%
        dplyr::filter(pathwayheat == "yes") %>% 
        dplyr::rename("Gene" = "ext_gene_ensembl")
      
      #create matrix for heatmap using genes in significant DE that are in pathway of choice
      vstgsea.mat <- vst.myc %>%
        dplyr::filter(., ensembl_gene_id %in% dds.sig$ensembl_gene_id) %>%
        distinct(Gene, .keep_all = TRUE) %>% 
        column_to_rownames(var = "Gene") %>%
        dplyr::select(.,-ensembl_gene_id, -pathwayheat) %>%
        as.matrix()

      #transform and scale and transform back
      vstgsea.mat <- t(scale(t(vstgsea.mat)))
      #color function buliding a colorRamp palette based on user input from palette choices
      colors.hmg <- c(colorGHeat(), "white", colorGHeat2())
 
      ht <- heatmaply(
        vstgsea.mat,
        #k_col = k_number,
        row_text_angle = 45,
        height = 600,
        width = 600,
        # colorbar = list(len=1, limits = c(-2, 2)),
        colors = colors.hmg,
        dendrogram = "column",
        show_dendrogram = TRUE,
        col_side_colors = cond,
        showticklabels = c(FALSE, FALSE)
      )

      ht
    })
 
    output$htwarnGSEA <- renderUI({
      if(null_check()){
        div(
          style = "text-align: center; font-size: 18px; color: red;",
       "No significant DEGs found in pathway"
        )
      } else {
        return(NULL)
      }
    })

    # download DE Heatmap ####
    output$downloadheatmap <- downloadHandler(
      filename = paste("GSEA_Heatmap", '.png', sep=''),
      content = function(file) {
        #generate dds results table 
        res.hmg <- DE_res$dds_res() 
 
        dds.sig <- res.hmg %>%
          dplyr::filter(padj < 0.05 & abs(`log2FoldChange`) >= 0.5)

        if(nrow(dds.sig) == 0) {
          return(NULL)
        }
        
        pathwaygsea4 <- gsea_file_values[[input$filechoice]]
        
        
        p2 <- unlist((pathwaygsea4[names(pathwaygsea4) %in% input$pathwaylistht]))
  
        
        #filter vst counts matrix for genes in pathway
        vst.myc <- vst() %>%
          mutate(., pathwayheat = ifelse(ext_gene_ensembl %in% p2, 'yes', 'no')) %>%
          dplyr::filter(pathwayheat == "yes") %>% 
          dplyr::rename("Gene" = "ext_gene_ensembl")
    
        vstgsea.mat <- vst.myc %>%
          dplyr::filter(., ensembl_gene_id %in% dds.sig$ensembl_gene_id) %>%
          distinct(Gene, .keep_all = TRUE) %>% 
          column_to_rownames(var = "Gene") %>%
          dplyr::select(.,-ensembl_gene_id, -pathwayheat) %>%
          as.matrix()
     
        vstgsea.mat <- t(scale(t(vstgsea.mat)))
        #create a colorRamp function based on user input in color palette choices
        colors.hmg <- colorRamp2(c(-2, 0, 2), c(colorGHeat(), "white", colorGHeat2()))
        ht = ComplexHeatmap::Heatmap(
          vstgsea.mat,
          name = "z scaled expression",
          col = colors.hmg,
          row_names_gp = gpar(fontsize = 6),
          column_names_gp = gpar(fontsize = 6),
          column_title = NULL,
          row_title = "DEGs in Pathway"
        ) 
        png(file)
        # draw heatmap object
        draw(ht)
        dev.off()
      }
    )
    
    observe({
      req(reset_trigger())
      
      updateMaterialSwitch(session, "fgseaTable", value = FALSE)
      updateMaterialSwitch(session, "rankedplot", value = FALSE)
      updateMaterialSwitch(session, "moustache", value = FALSE)
      updateMaterialSwitch(session, "eplot", value = FALSE)
      updateMaterialSwitch(session, "volcanoplot", value = FALSE)
      updateMaterialSwitch(session, "heatmap", value = FALSE)
      
    })
    
  })
}


  