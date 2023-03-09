##practice space for creating modules for each analysis within the EAGLe app

library(shiny)
# load files
vst.goi <- readRDS("data/vst.goi.rds")
geneexpressionUI<- function(id) {
  tagList(
    selectizeInput( #gene choice dropdown menu
      NS(id,"VSTCDgenechoice"),
      label=
        "Choose a gene for analysis",
      choices =
        NULL,
      selected = NULL,
      options = list(maxItems = NULL)
    ), #options for axis variables, fill variable, and plot filters
    radioButtons(NS(id, "XaxisVar_CDgene"), h4("X axis variable"),
                 choices = list("Value" = "xvalue",
                                "Gene" = "xgene", "Class" = "xclass"),selected = "xgene"),
    radioButtons(NS(id,"YaxisVar_CDgene"), h4("Y axis variable"),
                 choices = list("Value" = "yvalue",
                                "Gene" = "ygene"),selected = "yvalue"),
    radioButtons(NS(id,"FillVar_CDgene"), h4("color by:"),
                 choices = list("Gene" = "fillgene",
                                "Class" = "fillclass"),selected = "fillclass"),
    radioButtons(NS(id,"PrimMonobutton"), h4("Show only prim or mono gene expression"),
                 choices = list("Show Comparison" = "comparison", "Prim" = "prim", "Mono" = "mono"), selected = "comparison"),
    hr(),
    #add a facet toggle switch
    materialSwitch(NS(id,"genefacetbutton"), label = "Facet", value = FALSE, right = TRUE),
    hr(),
    #add palette choices for boxplot colors
    palettePicker(
      NS(id,"PaletteChoicesGene"),
      label = "Choose a color palette",
      choices = list(
        "Viridis" = list(
          "viridis" = viridis_pal(option = "viridis")(5),
          "magma" = viridis_pal(option = "magma")(5),
          "mako" = viridis_pal(option = "mako")(5),
          "plasma" = viridis_pal(option = "plasma")(5),
          "cividis" = viridis_pal(option = "cividis")(5))
      )
    ),

    hr(), #js functions to hide plot dimensions until selected
    materialSwitch(NS(id,"hidedims"), "Custom plot dimensions", value = FALSE, right = TRUE),

    #plot dimension input
    sliderInput(NS(id,"geneheightslider"), "Adjust plot height",
                min = 200, max = 1200, value = 600
    ),
    hr(),
    sliderInput(NS(id,"genewidthslider"), "Adjust plot width",
                min = 200, max = 1200, value = 800
    ),
    hr(),

    downloadButton(NS(id, "downloadGenePlot"), label = "Download Plot"),
      plotOutput(NS(id,"VSTCDplot"))

  )
}

geneexpressionServer<- function(id) {
  moduleServer(id,function(input, output, session) {
      updateSelectizeInput(session,"VSTCDgenechoice", choices = vst.goi$ext_gene, server = TRUE)
      #reactive function for for filtering vst data table based on user input
      datavst <-
        reactive({
          if(input$PrimMonobutton == "comparison") {
            vst.goi %>%
              dplyr::filter(ext_gene %in% input$VSTCDgenechoice)
          } else if(input$PrimMonobutton == "prim") {
            vst.goi %>%
              dplyr::filter(ext_gene %in% input$VSTCDgenechoice) %>%
              dplyr::filter(class == "prim")
          } else if(input$PrimMonobutton == "mono") {
            vst.goi %>%
              dplyr::filter(ext_gene %in% input$VSTCDgenechoice) %>%
              dplyr::filter(class == "mono")
          }
        })

      observeEvent(input$XaxisVar_CDgene, {
        if(input$XaxisVar_CDgene == "xvalue") {
          mychoices <- c("Gene" = "ygene")
        } else if(input$XaxisVar_CDgene=="xgene") {
          mychoices <- c("Value" = "yvalue")
        }else if(input$XaxisVar_CDgene == "xclass") {
          mychoices <- c("Value" = "yvalue")
        }
        #updateRadioButtons(session, "YaxisVar_CDgene", choices = mychoices)
      })
    # y axis reactive output based on radio buttons
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
     #facet toggle switch function to turn faceting on or off
      Gene_facet <-
        eventReactive(input$genefacetbutton, {
          if(input$genefacetbutton == TRUE) {
            facet_grid(cols = vars(class))
          } else(NULL)
        })

     # reactive function to tell ggplot which color palette to use for fill based on user input
      colorpalettechoicesfgene <-
        eventReactive(input$PaletteChoicesGene, {
          if(input$PaletteChoicesGene == "viridis") {
            scale_fill_viridis_d(option = "viridis")
          } else if(input$PaletteChoicesGene == "cividis") {
            scale_fill_viridis_d(option = "cividis")
          } else if(input$PaletteChoicesGene == "magma") {
            scale_fill_viridis_d(option = "magma")
          } else if(input$PaletteChoicesGene == "plasma") {
            scale_fill_viridis_d(option = "plasma")
          }else if(input$PaletteChoicesGene == "inferno") {
            scale_fill_viridis_d(option = "inferno")
          }
        })
      #reactive function for scale_Color_manual based on palette choice
      colorpalettechoicesgene <-
        eventReactive(input$PaletteChoicesGene, {
          if(input$PaletteChoicesGene == "viridis") {
            scale_color_viridis_d(option = "viridis")
          } else if(input$PaletteChoicesGene == "cividis") {
            scale_color_viridis_d(option = "cividis")
          } else if(input$PaletteChoicesGene == "magma") {
            scale_color_viridis_d(option = "magma")
          } else if(input$PaletteChoicesGene == "plasma") {
            scale_color_viridis_d(option = "plasma")
          }else if(input$PaletteChoicesGene == "inferno") {
            scale_fill_viridis_d(option = "inferno")
          }
        })
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
      #reactive wrapper for showing the plot dimensions options or hiding them based on toggle selection
      observe({
        toggle(id = "geneheightslider", condition = input$hidedims)
        toggle(id ="genewidthslider", condition = input$hidedims)
      })
      #plot output
      output$VSTCDplot <-
        renderPlot(
          width = function() input$genewidthslider,
          height = function() input$geneheightslider,
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
              colorpalettechoicesfgene() + #reactive  scale_fill_manual
              colorpalettechoicesgene() + #reactive scale_color_manual
              geom_point(alpha = 0.5,
                         position = position_jitterdodge(jitter.width = 0.2),
                         aes(color = class)) +
              theme_light() +
              sig_label_position() + # function for adjusted pvalues position and format on plot
              ylab("") +
              xlab("") +
              ggtitle("Gene Expression: Prim vs Mono")
          }) #end render plot

      output$downloadGenePlot <- downloadHandler(
        filename = function() { paste('GeneCentricPlot','.png', sep='') },
        content = function(file) {
          ggsave(file, device = "png", width = 8,
                 height = 8, dpi = 72)
        }
     )
    }
  )
}

geneexpressionApp <- function() {
  ui <- fluidPage(
    geneexpressionUI("Gene1")
  )
  
  server<- function(input, output, session) {
    geneexpressionServer("Gene1")
  }
  
  shinyApp(ui, server)
}

geneexpressionApp()

##example from masteringshiny chapter 19 for reference

# library(shiny)
# histogramUI <- function(id) {
#   tagList(
#     selectInput(NS(id, "var"), "Variable", choices = names(mtcars)),
#     numericInput(NS(id, "bins"), "bins", value = 10, min = 1),
#     plotOutput(NS(id, "hist"))
#   )
# }
# histogramServer <- function(id) {
#   moduleServer(id, function(input, output, session) {
#     data <- reactive(mtcars[[input$var]])
#     output$hist <- renderPlot({
#       hist(data(), breaks = input$bins, main = input$var)
#     }, res = 96)
#   })
# }
# histogramApp <- function() {
#   ui <- fluidPage(
#     histogramUI("hist1")
#   )
#   server <- function(input, output, session) {
#     histogramServer("hist1")
#   }
#   shinyApp(ui, server)
# }
# 
# histogramApp()
