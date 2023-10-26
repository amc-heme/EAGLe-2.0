
PW_UI <- function(id) {
  ns <- NS(id)
  fluidPage(
    theme =
      shinytheme("flatly"),
    titlePanel(
      "Differential Expression Tables and Plots"
    ),#end title
    sidebarLayout(
      sidebarPanel(
        tagList(
          checkboxGroupButtons(
            inputId = ns("PWChoices"),
            label = "Choose a Pairwise Comparison",
            choices = c("choice1", "choice2", "choice3"),
            status = "primary"
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
            colorUI(ns("color"), "Choose 1st color", "#0000FF"),
            colorUI(ns("color2"), "Choose 2nd color", "028a0f"),
            
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
            colorUI(ns("color3"), "Choose 1st color", "#0000FF"),
            colorUI(ns("color4"), "Choose 2nd color", "028a0f"),
            
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
            colorUI(ns("color5"),"Choose 1st color", "#0000FF"),
            colorUI(ns("color6"), "Choose 2nd color", "#FF0000"),
          )
        )
      ),
      mainPanel(
        conditionalPanel(
          ns = ns,
          condition = "input.PWDEtable == true",
          DTOutput(ns("results"))
        ),
        conditionalPanel(
          ns = ns,
          condition = "input.PWvolcano == true",
          girafeOutput(ns("pwvolplot"))
        ),
        conditionalPanel(
          ns = ns,
          condition = "input.PWMA == true",
          girafeOutput(ns("pwMAplot"))
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






















