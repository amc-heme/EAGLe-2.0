# DE volcano and MA plots
library(ggiraph)
DE_plots_UI <- function(id) {
  ns <- NS(id)
  tagList(
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
    # colorUI("col1", "Choose 1st color", "#0000FF"),
    # colorUI("col2", "Choose 2nd color", "028a0f"),
    girafeOutput(ns("volplot")),
    girafeOutput(ns("MAplot"))
  )
}

DE_plots_Server <- function(id, DEres) {
  moduleServer(id, function (input, output, session) {
    DEres <- reactive({
      DE_Server("DE1")
    })
    # color1 <- reactive({
    #   colorServer("col1")
    # })
    # color2 <- reactive({
    #   colorServer("col2") 
    # })
   
    output$volplot <- 
    
      renderGirafe({
        dds.res <- DEres()$dds.res
        colors <- c(magma(15)[9], "grey", viridis(15)[10] )#object for colors on volcano based on user input called from palette module
        if(input$DESeqvolcano == TRUE) { #only create plot if the  volcano switch is toggled
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
    
    output$MAplot <- 
      renderGirafe ({
        dds.res <- DEres()$dds.res
        colors <- c(magma(15)[9], "grey", viridis(15)[10] )#object for colors on volcano based on user input called from palette module
        if(input$DESeqMA == TRUE) { #only call plot if the MA plot switch is toggled
          ma <- ggplot(dds.res(), #call in the DE results from the DE module
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
  })
} 
  
DE_plots_App <- function() {
  ui <- fluidPage(
    DE_plots_UI("DE1")
  )
  server <- function(input, output, session) {
    DE_plots_Server("DE1")
  }
  shinyApp(ui, server)
}
DE_plots_App()  
  
  