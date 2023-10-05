## GSEA tab module

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
        
        conditionalPanel(
          ns = ns,
          condition = "input.fgseaTable == true",
          downloadButton(
            ns("downloadfgsea"),
            label =
              "Download GSEA Table"
          )
        ),
        hr(),
        #GSEA Waterfall####
        conditionalPanel(
          ns = ns,
          condition = "input.rankedplot == true",
          
          h4("Waterfall Plot Specific Options"),
          hr(),
          #color palette choices for waterfall plot
          colorUI(ns("color7"), "Choose color for plot", "#FF0000"),
          
          hr(), 
          #slider scale to choose how many pathways to load
          sliderInput(ns("howmanypathways"), "Choose How Many Pathways to Rank",
                      min = 5, max = 50, value = 15
          ),
          
          hr(),
          #js function to hide plot dimensions until selected
         # materialSwitch(ns("hidedimsWF"), "Custom plot dimensions", value = FALSE, right = TRUE),
          
          sliderUI(ns("rankedheightslider"), 200, 1000, 400, "Adjust Plot Height"),
          
          hr(),
          
          #call in module UI for width slider
          sliderUI(ns("rankedwidthslider"), 200, 1000, 600, "Adjust Plot Width"),
          
          downloadButton(
            ns("downloadranks"),
            label =
              "Download Waterfall"
          )
        ),
        hr(),
        
        #GSEA Moustache ####
        conditionalPanel(
          ns = ns,
          condition = "input.moustache == true",
          h4("Moustache Plot Specific Options"),
          
          hr(),
          #color palette choices for muostache plot
          colorUI(ns("color8"), "Choose color for plot", "#FF0000"),
          
          hr(), 
          # 
          # downloadButton(
          #   ns("downloadmoustache"),
          #   label =
          #     "Download Moustache Plot"
          # )
        ),
        
        hr(),
        #GSEA Enrichment plot ####
        conditionalPanel(
          ns = ns,
          condition = "input.eplot == true",
          
          h4("Enrichment Plot Specific Options"),
          #options for type of enrichment plot to load
          radioButtons(ns("topupordownbutton"), h5("Enrichment Plot Choices"), 
                       choices = list("Top Ranked Up Pathway" = "topup", "Top Ranked Down Pathway" = "topdown", "Pathway of Choice:" = "eplotpath"), selected = "topup"),
          h5("Choose a specific pathway"),
          #dropdown menu for specific pathway choices that is reactive to the pathway set dropdown menu
          selectizeInput(
            ns("pathwaylisteplot"),
            label=
              NULL,
            choices =
              NULL,
            selected = NULL ,
            options = list(maxItems = 1)
          ),
          
          hr(),
          downloadButton(
            ns("downloadeplot"),
            label =
              "Download Enrich. Plot"
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
          #color palette choices for volcano plot
          colorUI(ns("color9"), "Choose color for plot", "#FF0000"),
          
          # downloadButton(
          #   ns("downloadvolcano"),
          #   label =
          #     "Download Volcano Plot"
          # )
        ),
        hr(),
        #GSEA Heatmap ####
        conditionalPanel(
          ns = ns, 
          condition = "input.heatmap == true",
          h4("Heatmap Specific Options"),
          
          #color palette choices for heatmap
          colorUI(ns("color10"), "Choose 1st color", "#FF0000"),
          
          colorUI(ns("color11"), "Choose 2nd color", "#0000FF"),
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
        )
      ),
      
      
      mainPanel(
        
        conditionalPanel(
          ns = ns,
          condition = "input.fgseaTable == true",
          DTOutput( #add loading spinners
            ns("fgseaTable")
          )
        ),
        conditionalPanel(
          ns = ns,
          condition = "input.rankedplot == true",
          shinycssloaders::withSpinner( #add loading spinners
            plotOutput(
              ns("GSEAranked")
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
          InteractiveComplexHeatmapOutput(heatmap_id = ns("htgsea"))
        )
      )
    )
  )
}


GSEA_Server <- function(id, dds, t2g) {
  moduleServer(id, function(input, output, session) {
    #run GSEA for chosen pathway input
    ens2gene <- reactive({
      t2g[,c(2,3)]
    })
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

    # Extract the dds results in a tidy format
    res <- reactive({
      results(dds, tidy = TRUE)
    })
    
    # Add the human name of the gene to the last column, because that's what all of the pathways are annotated using
    res <- inner_join(res(), ens2gene(), by = c("row" = "ensembl_gene_id"))
    colnames(res())[8] <- 'HS_Symbol'
    
    # Select only the human gene symbol and the 'stat' from the results, remove NAs, and average the test stat for duplicate gene symbols
    res2 <- res %>%
      dplyr::select(HS_Symbol, stat) %>%
      na.omit() %>%
      distinct() %>%
      group_by(HS_Symbol) %>%
      summarize(stat = mean(stat))
    
    # Reconfigure the data
    ranks <- deframe(res2)
    ranks <- sort(ranks)
    #reactive expression to run fgsea and load results table for each chosen pathway
    gseafile <-
      eventReactive(input$filechoice,{
        pathwaygsea <- gsea_file_values[[input$filechoice]]
        fgseaRes <- fgsea::fgsea(pathways = pathwaygsea, stats = ranks, nproc = 10)
        fgseaResTidy <- fgseaRes %>%
          as_tibble() %>%
          dplyr::select(., -pval,-log2err, -ES) %>% 
          arrange(desc(NES))
        fgseaResTidy
      })
    
    output$downloadfgsea <- downloadHandler(
      filename = function() { paste("GSEATable", '.csv', sep='')},
      content = function(file) {
        write.csv(gseafile(),file)
      }
    )
  
    #filter Res table for chosen pathway to show in a waterfall plot
    gseafile_waterfall <-
      reactive({
        pathwaygsea <- gsea_file_values[[input$filechoice]]
        fgseaRes <-
          fgsea::fgsea(pathways = pathwaygsea,
                       stats = ranks,
                       nproc = 10)
        fgseaResTidy <- fgseaRes %>%
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
        gseafile()
      }
    })
    
    #### GSEA pathway ranks waterfall plot ####
    #call in singlecolor_palette module for plot
    colorWF <- 
      colorServer("color7")
    rankedheight <-
      sliderServer("rankedheightslider")
    rankedwidth <-
      sliderServer("rankedwidthslider")
    
    output$GSEAranked <- renderPlot(
      width = function()
        rankedwidth(),
      height = function()
        rankedheight(),
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
            theme_minimal(base_size = 16) +
            theme(axis.title = element_text(face = "bold"), title = element_text(face = "bold"))
        }
      }
    )
    
    #download button output- waterfall plot
    output$downloadranks <- downloadHandler(
      filename = function() { paste("Waterfall Plot", '.png', sep='') },
      content = function(file) {
        ggsave(file, device = "png", width = 8, height = 6, units = "in",dpi = 72)
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
        pathwaygsea <- gsea_file_values[[input$filechoice]]
        fgseaRes <-
          fgsea::fgsea(pathways = pathwaygsea,
                       stats = ranks,
                       nproc = 10)
        fgseaResTidy <-
          fgseaRes %>%
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
      colorServer("color8")
    output$GSEAMoustache <- renderGirafe(
      {
        if (input$moustache == TRUE) {
          #color object reactive to user input from plalette choice
          colors <- c('grey', colorM())
          ma <-
            ggplot(toplotMoustache(), aes(x = NES, y = padj, color = sig, tooltip = pathway)) +
            geom_point_interactive(alpha = 0.8, size = 0.5) +
            theme_minimal(base_size = 18) +
            theme(axis.title = element_text(face = "bold"), title = element_text(face = "bold")) +
            xlab('NES') +
            scale_color_manual(values = colors) +
            ylab('adjusted p-value') +
            ggtitle("Pathways from GSEA") +
            #only label pathways that sig
            #geom_text_repel(colour = "black", aes(label= ifelse(padj <0.05, as.character(pathway), ""), hjust=0,vjust=0)) +
            coord_cartesian(xlim = c(-3, 3), ylim = c(-0.1, 1)) 
          girafe(code = print(ma))
        }
      }
    )
    #download button output= moustache plot 
    # output$downloadmoustache <- downloadHandler(
    #   filename = function() { paste("Moustache Plot", '.png', sep='') },
    #   content = function(file) {
    #     ggsave(file, device = "png", width = 8, height = 6, units = "in",dpi = 72)
    #   }
    # )
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
      fgseaRes <-
        fgsea::fgsea(pathways = pathwaygsea,
                     stats = ranks,
                     nproc = 10)
      fgseaResTidy <- fgseaRes %>%
        as_tibble() %>%
        dplyr::select(., -pval,-log2err, -ES) %>% 
        arrange(desc(NES))
      if(input$topupordownbutton == "topup") {
        top.UP.path <- as.character(fgseaResTidy[1,1])
        plotEnrichment(pathwaygsea[[top.UP.path]],
                       ranks) + labs(title=top.UP.path)
      } else if(input$topupordownbutton == "topdown") {
        top.DOWN.path <- as.character(fgseaResTidy[nrow(fgseaResTidy), 1])
        plotEnrichment(pathwaygsea[[top.DOWN.path]],
                       ranks) + labs(title=top.DOWN.path)
      } else if(input$topupordownbutton == "eplotpath") {
        plotEnrichment(pathwaygsea[[input$pathwaylisteplot]],
                       ranks) + labs(title= input$pathwaylisteplot)
      }
    })
    #download button output- enrichment plot
    output$downloadeplot <- downloadHandler(
      filename = function() { paste("Enrichment Plot", '.png', sep='') },
      content = function(file) {
        ggsave(file, device = "png", width = 8, height = 6, units = "in",dpi = 72)
      }
    )
    # GSEA Volcano plot ####
    #reactive expression for selection of specific pathway by user input
    observe({
      pathwaygsea <- gsea_file_values[[input$filechoice]]
      updateSelectizeInput(session,"pathwaylist", choices = names(pathwaygsea), server = TRUE)})
    #function to select only genes from the DE object that are found in the chocen pathway
    dds.res.pathways <- reactive({
      pathwaygsea <- gsea_file_values[[input$filechoice]]
      p <-
        unlist((pathwaygsea[names(pathwaygsea) %in% input$pathwaylist]))
      
      dds.res.pathways <- dds.res %>%
        mutate(., genes_in_pathway = ifelse(Gene %in% p, 'yes', 'no'))
      #print(dds.res.pathways)
    })
    #reactive title for volcano based on specific pathway choice
    gseavol_title <-
      eventReactive(input$pathwaylist, {
        paste(input$pathwaylist)
      })
    #call in singlecolor_module for plot
    colorVol <- 
      colorServer("color9")
    output$GSEAvolcano <- renderGirafe ({
      #color object reactive to user input from palette chpice
      colors <- 
        c("grey", colorVol())
      if (input$volcanoplot == TRUE) {
       v <- ggplot(
          data = (dds.res.pathways() %>% arrange(., (genes_in_pathway))),
          aes(
            x = log2FoldChange,
            y = -log10(padj),
            col = genes_in_pathway[],
            tooltip = Gene
          )
        ) +
          theme_light(base_size = 14) +
          theme(axis.title = element_text(face = "bold"), title = element_text(face = "bold")) +
          geom_point_interactive(size = 1, alpha = 0.5) +
          scale_color_manual(values = colors) +
          # geom_text_repel(
          #   max.overlaps = 1500,
          #   colour = "black",
          #   aes( #only label is gene is in pathway and sig expression 
          #     label = ifelse(
          #       genes_in_pathway == 'yes' & log2FoldChange > 1.5,
          #       as.character(Gene),
          #       ""
          #     )
          #   ),
          #   hjust = 0,
          #   vjust = 0
          # ) +
          ggtitle(gseavol_title()) + #reactive title
          xlab("log2foldchange")
        girafe(code = print(v))
      }
    })
    #download button output- gsea volcano plot
    # output$downloadvolcano <- downloadHandler(
    #   filename = function() { paste("Volcano Plot", '.png', sep='') },
    #   content = function(file) {
    #     ggsave(file, device = "png", width = 8, height = 6, units = "in",dpi = 72)
    #   }
    # )
    #GSEA heatmap ####
    #reactive expression for selected pathway choice for heatmap
    observe({
      pathwaygsea <- gsea_file_values[[input$filechoice]]
      updateSelectizeInput(session,"pathwaylistht", choices = names(pathwaygsea), server = TRUE)})
    #reactive expression of title for heatmap
    gseaht_title <-
      eventReactive(input$pathwaylistht, {
        print(input$pathwaylistht)
      })
    #call in single color_module for plot palette
    colorGHeat <- 
      colorServer("color10")
    colorGHeat2 <-
      colorServer("color11")
    #interactive heatmap must be wrapped in a reactive expression
    observeEvent(input$pathwaylistht, {
      if(input$heatmap == TRUE) {
        pathwaygsea <- gsea_file_values[[input$filechoice]]
        
        p <- unlist((pathwaygsea[names(pathwaygsea) %in% input$pathwaylistht]))
        #filter DE object for significance 
        dds.sig <- dds.res %>%
          dplyr::filter(padj < 0.05 & abs(`log2FoldChange`) >= 0.5)
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
        #filter vst counts matrix for genes in pathway
        vst.myc <- vst %>% 
          mutate(., pathwayheat = ifelse(ext_gene %in% p, 'yes', 'no')) %>% 
          dplyr::filter(pathwayheat == "yes")
        #create matrix for heatmap using genes in significant DE that are in pathway of choice
        vstgsea.mat <- vst.myc %>%
          dplyr::filter(., ensembl_gene_id %in% dds.sig$ensembl_gene_id) %>%
          column_to_rownames(., var = "ext_gene") %>%
          dplyr::select(.,-ensembl_gene_id, -pathwayheat) %>%
          as.matrix()
        #transform and scale and transform back
        vstgsea.mat <- t(scale(t(vstgsea.mat)))
        #color function buliding a colorRamp palette based on user input from palette choices
        colors = colorRamp2(c(-2, 0, 2), c(colorGHeat(), "white", colorGHeat2()))
        htgsea = draw(ComplexHeatmap::Heatmap(
          vstgsea.mat,
          name = paste(gseaht_title(), fontsize = 6),
          col = colors,
          #"paste(input$pathwaylist, sep = ",")",
          row_names_gp = gpar(fontsize = 6),
          column_km = 2,
          top_annotation = HeatmapAnnotation(class = anno_block(gp = gpar(fill = c("white", "white")),
                                                                labels = c("prim", "mono"), 
                                                                labels_gp = gpar(col = "black", fontsize = 10))),
          column_title = NULL,
          row_title = NULL))
        
        makeInteractiveComplexHeatmap(input, output, session, htgsea, "htgsea")
      }
    })
    
  })
}

# GSEA_App <- function() {
#   ui <- fluidPage(
#     GSEA_UI("GSEA1")
#   )
#   server <- function(input, output, session) {
#     GSEA_Server("GSEA1", dds, ens2gene_HS)
#   }
#   shinyApp(ui, server)
# }
# GSEA_App()

  