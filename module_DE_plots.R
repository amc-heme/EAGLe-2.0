# DE volcano and MA plots

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

    plotlyOutput(ns("MAplot")),
    plotlyOutput(ns("volplot"))
  )
}

DE_plots_Server <- function(id, run_DE) {
  moduleServer(id, function (input, output, session) {
    DEres <- reactive({
      run_DE()
    })
    output$volplot <-
      renderPlotly({
        colors <- c("red", "grey", "green")#object for colors on volcano based on user input called from palette module
        if(input$DESeqvolcano == TRUE) { #only create plot if the  volcano switch is toggled
        p<- ggplot(DEres(), aes( #call in the DE results from the DE module
          x = `log2FoldChange(Prim/Mono)`,
          y = -log10(padj),
          col = DiffExp,
          text = Gene
        )) +
          geom_point(size = 1, alpha = 0.5) +
          theme_light() +
          scale_color_manual(values = colors) +
          ggtitle("DE Volcano Plot") +
          coord_cartesian(xlim = c(-10, 7))
        ggplotly(p)
        }
      })
    
    output$MAplot <- 
      renderPlotly ({
        colors <- c("red", "grey", "green") #object for color choices dependent on user input called from palette module
        if(input$DESeqMA == TRUE) { #only call plot if the MA plot switch is toggled
          ma<-ggplot(DEres(), #call in the DE results from the DE module
                 aes(
                   x = log2(baseMean),
                   y = `log2FoldChange(Prim/Mono)`,
                   col = DiffExp
                 )) +
          geom_point(alpha = 0.8, size = 0.5) +
          geom_hline(aes(yintercept = 0)) +
          scale_color_manual(values = colors) +
          theme_light() +
          ylim(c(
            min(DEres()$`log2FoldChange(Prim/Mono)`),
            max(dDEres()$`log2FoldChange(Prim/Mono)`)
          )) +
          ggtitle("DE MA Plot") +
          xlab("log2 Mean Expression") +
          ylab("Log2 Fold Change")
        
        ggplotly(ma)
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
  
  