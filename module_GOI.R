# Gene centric module
vst.goi <- readRDS("data/vst.goi.rds")
goi_UI <- function(id) {
  ns <- NS(id)
  fluidPage(
    theme=
      shinytheme("flatly"),
    titlePanel(
      "Gene Expression Plot"
    ),
    h6("*p values are indicated for the comparison of gene expression between prim and mono samples"),
    sidebarLayout(
      sidebarPanel(
        useShinyjs(), #this is needed for javascript functions
        selectizeInput( #gene choice dropdown menu
          ns("VSTCDgenechoice"),
          label=
            "Choose a gene for analysis",
          choices =
            NULL,
          selected = NULL,
          options = list(maxItems = NULL)
        ), #options for axis variables, fill variable, and plot filters
        radioButtons(ns("XaxisVar_CDgene"), h4("X axis variable"),
                     choices = list("Value" = "xvalue",
                                    "Gene" = "xgene", "Class" = "xclass"),selected = "xgene"),
        radioButtons(ns("YaxisVar_CDgene"), h4("Y axis variable"),
                     choices = list("Value" = "yvalue",
                                    "Gene" = "ygene"),selected = "yvalue"),
        radioButtons(ns("FillVar_CDgene"), h4("color by:"),
                     choices = list("Gene" = "fillgene",
                                    "Class" = "fillclass"),selected = "fillclass"),
        # radioButtons(ns("PrimMonobutton"), h4("Show only prim or mono gene expression"),
        #              choices = list("Show Comparison" = "comparison", "Prim" = "prim", "Mono" = "mono"), selected = "comparison"),
        hr(),
        #add a facet toggle switch
        materialSwitch(ns("genefacetbutton"), label = "Facet", value = FALSE, right = TRUE),
        hr(),
        
        #add palette choices for boxplot colors
        paletteUI(ns("palette2")),
        
        hr(), #js functions to hide plot dimensions until selected
        #materialSwitch("hidedims", "Custom plot dimensions", value = FALSE, right = TRUE),
        
        
        #plot dimension input
        sliderUI(ns("plotheightslider"), 200, 1200, 600, "Adjust Plot Height"),
        
        #  hr(),
        # 
        sliderUI(ns("plotwidthslider"), 200, 1200, 800, "Adjust Plot Width"),
        
        hr(),
        
        downloadButton(ns("downloadGenePlot"), label = "Download Plot"),
        
      ),
      
      mainPanel( #add loading spinner
        shinycssloaders::withSpinner(
          plotOutput(
            ns("VSTCDplot")
          )
        )
      )
    )
  )
}

goi_Server <- function(id, dds, t2g) {
  moduleServer(id, function(input, output, session) {
    ##Gene Centric output ####
    updateSelectizeInput(session,"VSTCDgenechoice", choices = vst.goi$ext_gene, server = TRUE)
    #create dds results table for use in the table generated for the plot
    dds.res <- data.frame(results(dds)) %>%
      rownames_to_column(., var = 'ensembl_gene_id') %>%
      dplyr::select(., ensembl_gene_id, baseMean, log2FoldChange, padj) %>%
      left_join(unique(dplyr::select(t2g, c(ensembl_gene_id, ext_gene))), ., by = 'ensembl_gene_id') %>%
      dplyr::rename(., Gene = ext_gene) %>%
      mutate(., DiffExp = ifelse(padj < 0.05 & log2FoldChange >= 0.5, 'up',
                                 ifelse(padj < 0.05 & log2FoldChange <= -0.5, 'down', 'no'))) %>%
      na.omit(.)
    #extract counts from dds file to use in vst
    dds.counts <- reactive({
      counts(dds)
      })
    vsd <- reactive({
      vst(dds.counts(), blind = F) 
    })
    vst <- data.frame(assay(vsd())) %>% 
      rownames_to_column(., var = "ensembl_gene_id") %>% 
      left_join(unique(dplyr::select(t2g_mm, c(ensembl_gene_id, ext_gene))), ., by = 'ensembl_gene_id') %>% 
      na.omit(.)
    #create a data table filtered for only mono sample type
    res <- (label.jordan.m0m5) %>%
      mutate(., class = condition) %>%
      dplyr::select(., SRR, class) %>%
      filter(., class == "mono")
    
    #join DE res and VST counts matrix to create data table with class, padj, and gene expression values for each sample
    vst.goi <- as_tibble(vst) %>%
      melt(.) %>%
      mutate(., class = ifelse(variable %in% res$SRR, 'mono', 'prim')) %>% 
      dplyr::filter(ensembl_gene_id %in% dds.res$ensembl_gene_id) %>%
      left_join(unique(dplyr::select(dds.res, c(
        ensembl_gene_id, padj
      ))), ., by = 'ensembl_gene_id') 
    
    #factor class and variable(sample id)
    vst.goi$class <-
      factor(vst.goi$class, levels = c('prim', 'mono'))
    
    vst.goi$variable <- factor(vst.goi$variable)
    
    #reactive function for for filtering vst data table based on user input 
    datavst <-
      reactive({
        vst.goi %>% 
             dplyr::filter(ext_gene %in% input$VSTCDgenechoice)
        # if(input$PrimMonobutton == "comparison") {
        #   vst.goi %>% 
        #     dplyr::filter(ext_gene %in% input$VSTCDgenechoice)
        # } else if(input$PrimMonobutton == "prim") {
        #   vst.goi %>%
        #     dplyr::filter(ext_gene %in% input$VSTCDgenechoice) %>%
        #     dplyr::filter(class == "prim")
        # } else if(input$PrimMonobutton == "mono") {
        #   vst.goi %>%
        #     dplyr::filter(ext_gene %in% input$VSTCDgenechoice) %>%
        #     dplyr::filter(class == "mono")
        # }
      })
    
    
    #make sure duplicate selections are not allowed with radio buttons
    observeEvent(input$XaxisVar_CDgene, {
      if(input$XaxisVar_CDgene == "xvalue") {
        mychoices <- c("Gene" = "ygene")
      } else if(input$XaxisVar_CDgene=="xgene") {
        mychoices <- c("Value" = "yvalue")
      }else if(input$XaxisVar_CDgene == "xclass") {
        mychoices <- c("Value" = "yvalue")
      }
      updateRadioButtons(session, "YaxisVar_CDgene", choices = mychoices)
    })
    
    
    #x axis reactive output based on radio buttons
    xvar_CDgene <-
      eventReactive(input$XaxisVar_CDgene, {
        if (input$XaxisVar_CDgene == "xvalue") {
          "value"
        } else if (input$XaxisVar_CDgene == "xgene") {
          "ext_gene"
        } else if(input$XaxisVar_CDgene == "xclass") {
          "class"
        }
      })
    #y axis reactive output based on radio buttons
    yvar_CDgene <-
      eventReactive(input$YaxisVar_CDgene, {
        if (input$YaxisVar_CDgene == "yvalue") {
          "value"
        } else if (input$YaxisVar_CDgene == "ygene") {
          "ext_gene"
        }
      })
    
    
    #fill reactive output based on radio buttons
    fillvar_CDgene <-
      eventReactive(input$FillVar_CDgene, {
        if (input$FillVar_CDgene == "fillclass") {
          "class"
        } else if (input$FillVar_CDgene == "fillgene") {
          "ext_gene"
        }
      })
    # facet toggle switch function to turn faceting on or off
    Gene_facet <-
      eventReactive(input$genefacetbutton, {
        if(input$genefacetbutton == TRUE) {
          facet_grid(cols = vars(class))
        } else(NULL)
      })
    
    #call in palette module for plot
    colorpaletteGene <- 
      paletteServer("palette2")
    
    # function for adding padj values to plot, position needs to change when x and y variables change for readability
    sig_label_position <- reactive({
      value <- vst.goi$value
      if(input$XaxisVar_CDgene == "xvalue") {
        geom_text(aes(x = max(value), label = paste("p=",format(padj, digit = 1, scientific = T))),check_overlap = T) 
      } else if(input$XaxisVar_CDgene == "xgene") {
        geom_text(aes(y = max(value), label = paste("p=",format(padj, digit = 1, scientific = T))),check_overlap = T) 
      } else if(input$XaxisVar_CDgene == "xclass") {
        geom_text(aes(y = max(value), label = paste("p=",format(padj, digit = 1, scientific = T))),check_overlap = T) 
      }
    })
    
    geneheight <-
      sliderServer("plotheightslider")
    
    genewidth <-
      sliderServer("plotwidthslider")
    

    #plot output
    output$VSTCDplot <-
      renderPlot(
        width = function() genewidth(), #input$genewidthslider,
        height = function() geneheight(), #input$geneheightslider,
        res = 120,
        {
          ggplot(datavst(),
                 aes(
                   x = .data[[xvar_CDgene()]],
                   y =  .data[[yvar_CDgene()]],
                   fill = .data[[fillvar_CDgene()]]
                 )) +
            geom_boxplot(outlier.shape = NA) +
            Gene_facet() + #reactive faceting
            scale_fill_viridis_d(option = colorpaletteGene()) + #reactive  scale_fill_manual from module
            scale_color_viridis_d(option = colorpaletteGene()) + #reactive scale_color_manual from module
            geom_point(alpha = 0.5,
                       position = position_jitterdodge(jitter.width = 0.2),
                       aes(color = class)) + 
            theme_light() +
            sig_label_position() + # function for adjusted pvalues position and format on plot
            ylab("") +
            xlab("") +
            ggtitle("Gene Expression")
        }) #end render plot
    
    output$downloadGenePlot <- downloadHandler(
      filename = function() { paste('GeneCentricPlot','.png', sep='') },
      content = function(file) {
        ggsave(file, device = "png", width = 8,
               height = 8, dpi = 72)
      }
    )
    
    
  })
}

