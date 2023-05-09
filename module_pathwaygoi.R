## pathway GOI module

pathway_UI <- function(id) {
  ns <- NS(id)
  fluidPage(
    theme =
      shinytheme("flatly"),
    titlePanel(
      "GSEA: Interrogation of pathways containing a gene of interest"
    ),#end title
    h4("Positive NES is upregulated in Primitive cells and negative NES is upregulated in Monocytic cells"),
    sidebarLayout( 
      sidebarPanel( 
        useShinyjs(),
        #gene list dropdown menu of all genes from DE table
        selectizeInput(
          ns("Pathwaygenechoice"),
          label=
            "Choose a gene of interest",
          choices =
            NULL,
          selected = NULL,
          options = list(maxItems = 1) #only one gene can be selected at once
        ),
        hr(),
        #pathway set dropdown list
        selectInput(ns("genefilechoice"), "Choose gmt file to plot pathways containing the gene of interest",
                    choices = c(Hallmark = "hallmark", GOall = "GOall", GOmolecular = "GOmolec", 
                                GOcellcomp = "GOcellcomp", GObio = "GObio", TFtargets = "TFtargets",
                                allRegular = "allReg", Wiki = "wiki", Reactome = "reactome", KEGG = "KEGG",
                                Positional = "positional", Biocarta = "biocarta", lsc = "lsc", aeg = "aeg")),
        
        hr(),
        #color palette options for GOI
        paletteUI(ns("paletteGOI")),
        
        hr(),
    
        
        sliderUI(ns("goiheightslider"), 200, 1000, 400, "Adjust plot height"
        ),
        hr(),
        
        sliderUI(ns("goiwidthslider"), 200, 1000, 600, "Adjust plot width"
        ),
        downloadButton(
          ns("downloadGOI"),
          label =
            "Download Plot"
        )
      ),
      mainPanel(
        shinycssloaders::withSpinner( #add loading spinners
          plotOutput(
            ns("PathwaysGenePlot")
          )
        )
      )
    )
  )

}

pathway_Server <- function(id, dds, ens2gene_HS, dds.res) {

    moduleServer(id, function(input,output,session) {
    
    #object for pathway choice files to use with the goi pathway input
    gene_gsea_file_values <- list("hallmark" = pathways.hallmark,
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
                                  "biocarta" = pathways.Positional,
                                  "lsc" = pathways.lsc,
                                  "aeg" = pathways.aeg)
    
    #reactive title
    gsea_gene_title <- 
      eventReactive(input$Pathwaygenechoice, {
        paste(input$Pathwaygenechoice)
      })
    
    #reactive color palette from module
    colorGOI <-
      paletteServer("paletteGOI")
    
    #height and width sliders
    goiheights <- 
      sliderServer("goiheightslider")
    goiwidths <- 
      sliderServer("goiwidthslider")
   
    # Extract the dds results in a tidy format
    res <- results(dds, tidy = TRUE)
    
    # Add the human name of the gene to the last column, because that's what all of the pathways are annotated using
    res <- inner_join(res, ens2gene_HS, by = c("row" = "ensembl_gene_id"))
    colnames(res)[8] <- 'HS_Symbol'
    
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
    
    genecentricgseaplot <- reactive({
      #load chosen pathway file based on reactive input 
      genepathwaygsea <- (gene_gsea_file_values[[input$genefilechoice]])
      #load fgsea table data for chosen pathway
      fgseaRes <-
        fgsea::fgsea(pathways = genepathwaygsea,
                     stats = ranks,
                     nproc = 10)
      #create tidy table
      fgseaResTidy <- fgseaRes %>%
        as_tibble() %>%
        arrange(desc(NES))
      
      #create object for storing pathways that contain the chosen gene 
      goi_paths <- genepathwaygsea %>% keep(grepl(input$Pathwaygenechoice, genepathwaygsea))
      goi_paths <- list(grep(input$Pathwaygenechoice, genepathwaygsea))
      #filter gsea table for pathways in which the GOI is in the leading edge
      goi_paths <- fgseaResTidy %>%
        dplyr::filter(grepl(input$Pathwaygenechoice, leadingEdge)) %>%
        mutate(., class = ifelse(NES <0, 'Mono', 'Prim'))
      
      #create object for gene reactive input
      GOI <- input$Pathwaygenechoice
      #make a column that says yes if goi in that pathway
      goi_paths$GOI <- "pathways with GOI"
      #filter gsea table to find pathways that do not include the GOI in the leading edge
      nongoi_paths <- fgseaResTidy %>%
        dplyr::filter(!grepl(input$Pathwaygenechoice, leadingEdge)) %>%
        mutate(., class = ifelse(NES <0, 'Mono', 'Prim'))
      
      #put no for pathways that do not contain the goi
      nongoi_paths$GOI <- "pathways NOT with GOI"
      #bind the two filtered data frames into one for plotting
      allgoi_paths <- rbind.data.frame(goi_paths, nongoi_paths)
    })
    
    output$PathwaysGenePlot <- renderPlot(
      width = function() goiwidths(),
      height = function() goiheights(),
      {
        ggplot(genecentricgseaplot(), aes(
          x = class,
          y = NES,
          color = (padj < 0.05)
        )) +
          geom_boxplot()  +
          scale_color_viridis_d(option = colorGOI()) +
          facet_wrap( ~ GOI, scales = "free") +
          theme_light(base_size = 18) +
          theme(axis.title = element_text(face = "bold"), title = element_text(face = "bold")) +
          geom_hline(yintercept = 0, linetype = "dashed") +
          plot_annotation(
            title = "Pathways with and without Gene of Interest",
            theme =
              theme(
                plot.title =
                  element_text(
                    face = "bold",
                    hjust = 0.5,
                    size = 16
                  )
              )
          )
      })
    
    output$downloadGOI <- downloadHandler(
      filename = function() { paste("Gene of Interest Plot", '.png', sep='') },
      content = function(file) {
        ggsave(file, device = "png", width = 8, height = 6, units = "in",dpi = 72)
      }
    )
    # #gene list for gene centric pathway analysis
    updateSelectizeInput(session,"Pathwaygenechoice", choices = dds.res$Gene, server = TRUE)
    
  })
}



}