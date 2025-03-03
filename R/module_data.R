datasets <- 
  read_yaml("./data.yaml")
data_UI <- function(id) {
  ns <- NS(id)
  fluidPage(
    theme =
      shinytheme("flatly"),
    useShinyjs(),
    useWaiter(),
 
    titlePanel("Choose a dataset, model term, and comparison for DE and visualizations"
               ),
   
    sidebarLayout(
      sidebarPanel(
        
        tagList(
          
          selectInput(
            ns("datainput"), label = NULL,
             choices = c("Pei et al, 2020" = "Cancer_Discovery",
                         "Ye et al, 2016" = "Ye_16",
                         "Ye et al, 2020" = "Ye_20",
                         "Pollyea Ven/Aza, 2018" = "Venaza",
                         "Lagadiou, 2013" = "Lagadinou",
                         "BEAT AML: Ven Response Quantile" = "BEAT_quantile",
                         "BEAT AML: FAB Blast Morphology" = "BEAT_FAB",
                         "BEAT AML: Denovo vs Relapse" = "BEAT_Denovo.Relapse",
                         "TCGA: FAB" = "TCGA_FAB",
                         "TCGA: RAS mutation" = "TCGA_RAS",
                         "TCGA: NPM1 mutation" = "TCGA_NPM1",
                         
                         
                         
                         "Lee et al, 2018" = "Lee"), selected = "Cancer_Discovery"
          ),
          selectizeInput(
            ns("DEmodel"),
            label = "Metadata value to be tested:",
            choices = NULL,
            selected = NULL,
            options = list(maxItems = 1)
          ),
          selectizeInput(
            ns("pwc"),
            label = "Choose to run LRT OR a pairwise comparison for DE testing",
            choices = NULL,
            selected = NULL,
            options = list(maxItems = 1)
          ),
          hr(),
          actionButton(
            (ns("runDE1")),
            "Continue"
          ),
          tags$head( #css to move the pop up window to the top and middle of the screen
            tags$style(
              HTML(".shiny-notification {
             position:fixed;
             top: calc(50%);
             left: calc(40%);
             }
             "
              )
            )
          )
        )
      ),
      mainPanel(
        uiOutput(ns("Data_text")),
        textOutput(ns("LRT_alert")),
        bsAlert("alert")
      )
    )
  )
}

data_Server <- function(id) {
  moduleServer(id, function(input, output, session){
    shinyjs::disable("runDE1")
    output$Data_text <- renderUI({
      if (input$datainput == "Cancer_Discovery") {
        intro <- paste(
          "Pei et al. have previously shown that ROS-low enriched LSCs
             from primitive and monocytic AML differ significantly in their
             transcriptome and metabolism.
             ROS-low LSCs from monocytic AML in particular, are less dependent
             on BCL2 therefore are more likely to be resistant to
             Venetoclax-based therapy. The analysis will be comparing ROS-low enriched LSC's
             from monocytic cells vs primitive cells"
        )
        info <-
          paste(
            "**GSEA will not be run with LRT on datasets with more than 2 conditions.
        If you are interested in pathway analysis, please choose a pairwise comparison"
          )
       
       intro_style <-
         paste("<div style='font-size: 18px;'>",
               intro,
               "</div><br><br><div>",
               "<div style='font-size: 18px;'>",
               info,
               "</div>")
        HTML(intro_style)
      } else if (input$datainput == "Ye_16") {
       intro <-  paste(
          "Comparison of transcriptomes of LSC's from blood, bone marrow, spleen,
           gonadal adipose tissue, and normal bone marrow in mice."
        )
       intro_style <- paste("<div style='font-size: 18px;'>", intro, "</div>")
       HTML(intro_style)
      } else if (input$datainput == "Ye_20") {
        intro <- paste(
          "Transcriptomes of LSCs from bone marrow and liver in mice
               were compared to better understand the biology of liver LSCs.
               A bcCML model (BCR-ABL + Nup98-Hoxa9) was used. Liver and bone
               marrow samples were combined from 3 mice in 3 different cohorts
               before sequencing."
        )
        info <-
          paste(
            "**GSEA will not be run with LRT on datasets with more than 2 conditions.
        If you are interested in pathway analysis, please choose a pairwise comparison"
          )
        
        intro_style <-
          paste("<div style='font-size: 18px;'>",
                intro,
                "</div><br><br><div>",
                "<div style='font-size: 18px;'>",
                info,
                "</div>")
        HTML(intro_style)
      } else if (input$datainput == "Venaza") {
        intro <- paste(
          "Time zero pheresis was taken from 3 patients, then collected
                 again at 6hr and 24hr post ven/aza treatment. This analysis will compare
                 24hr or 6h vs time zero (control)"
        )
        info <-
          paste(
            "**GSEA will not be run with LRT on datasets with more than 2 conditions.
        If you are interested in pathway analysis, please choose a pairwise comparison"
          )
        
        intro_style <-
          paste("<div style='font-size: 18px;'>",
                intro,
                "</div><br><br><div>",
                "<div style='font-size: 18px;'>",
                info,
                "</div>")
        HTML(intro_style)
      } else if (input$datainput == "Lagadinou") {
       intro <- paste(
          "BCL-2 inhibition targets oxidative phosphorylation and
               selectively eradicates quiescent human leukemia stem cells.
               ROS high and ROS low LSCs were either treated with 5ul of PTL
               or not treated before sequencing."
        )
       info <-
         paste(
           "**GSEA will not be run with LRT on datasets with more than 2 conditions.
        If you are interested in pathway analysis, please choose a pairwise comparison"
         )
       
       intro_style <-
         paste("<div style='font-size: 18px;'>",
               intro,
               "</div><br><br><div>",
               "<div style='font-size: 18px;'>",
               info,
               "</div>")
       HTML(intro_style)
      } else if (input$datainput %in% c("TCGA_FAB", "TCGA_NPM1", "TCGA_RAS")) {
        intro <- paste(
          "The Cancer Genome Atlas (TCGA) project published gene expression
              profiles of 151 primary AML patients along with their mutational
              profiles and clinical characteristics. The TCGA-AML dataset parsed
              by various mutational/clinical variables including the 
              French-American-British subtypes, karyotype, RAS mutation status,
              and NPM1 mutation status."
        )
        info <-
          paste(
            "**GSEA will not be run with LRT on datasets with more than 2 conditions.
        If you are interested in pathway analysis, please choose a pairwise comparison"
          )
        
        intro_style <-
          paste("<div style='font-size: 18px;'>",
                intro,
                "</div><br><br><div>",
                "<div style='font-size: 18px;'>",
                info,
                "</div>")
        HTML(intro_style)
      } else if (input$datainput %in% c("BEAT_quantile", "BEAT_FAB", "BEAT_Denovo.Relapse")) {
        intro <- paste(
          "The BEAT-AML project published gene expression profiles of ~400
              primary AML patients along with their mutational profiles and clinical
              characteristics. The BEAT-AML dataset was parsed by various 
              mutational/clinical variables including the French-American-British
              subtypes, Venetoclax response quantile, and disease stage (de novo vs. relapse).
          <br/> Venetoclax response quantiles were split with q1 being the top 25% of samples with the highest 
          AUC values and q4 being the bottom 25% of samples with the lowest AUC values."
        )
        info <-
          paste(
            "**GSEA will not be run with LRT on datasets with more than 2 conditions.
        If you are interested in pathway analysis, please choose a pairwise comparison"
          )
        
        intro_style <-
          paste("<div style='font-size: 18px;'>",
                intro,
                "</div><br><br><div>",
                "<div style='font-size: 18px;'>",
                info,
                "</div>")
        HTML(intro_style)
      } else{
        intro <- paste(
          "Genome wide expression from 12 AML patient samples with either
              prior complete remission, or no prior complete remission."
        )
        info <-
          paste(
            "**GSEA will not be run with LRT on datasets with more than 2 conditions.
        If you are interested in pathway analysis, please choose a pairwise comparison"
          )
        
        intro_style <-
          paste("<div style='font-size: 18px;'>",
                intro,
                "</div><br><br><div>",
                "<div style='font-size: 18px;'>",
                info,
                "</div>")
        HTML(intro_style)
      }
    })
    #disable the "continue" button whenever the user changes dataset choice
    observeEvent(input$datainput,{
      shinyjs::disable("runDE1")
    })
    #enable the "continue" button after the dataset information loads and a model and pairwise choice have been made to prevent the app crashing
    observe({
      req(input$DEmodel)
      req(input$pwc)
      shinyjs::enable("runDE1")
    })
    #remove the gsea message when the "continue" button is selected
    observeEvent(input$runDE1, {
      removeNotification("gsea_message")
    })
    # model choice ####
    DEModelChoices <- function(dataset) {
      m_choices <- switch(
        dataset,
        "Cancer_Discovery" = "condition",
        "Ye_16" = "Source",
        "Ye_20" = "Source",
        "Venaza" = "condition",
        "Lagadinou" = "Treatment",
        "BEAT_quantile" = "quantile",
        "BEAT_FAB" = "FAB_BlastMorphology",
        "BEAT_Denovo.Relapse" = "Denovo.Relapse",
        "TCGA_FAB" = "FAB",
        "TCGA_RAS" = "RAS_mut",
        "TCGA_NPM1" = "NPM1_mut",
        "Lee" = "prior_cr",
        default = character(0)
      )
    }


    observe({
      model_choices <- DEModelChoices(input$datainput)

      updateSelectizeInput(session,
                           "DEmodel",
                           choices = model_choices,
                           selected = NULL)
    })

    
    
    PWChoices <- function(dataset, model) { 
      if(dataset %in% c("Cancer_Discovery", "Ye_16", "Ye_20", "Venaza", "Lagadinou", "Lee")) {
      choices <- switch(
        paste(dataset, model, sep = "_"),
        "Cancer_Discovery_condition" = "LRT",
        "Ye_16_Source" = c("LRT", "blood_vs_bone_marrow", "gonadal_adipose_tissue_vs_bone_marrow", "normal_bm_vs_bone_marrow", "spleen_vs_bone_marrow"),
        "Ye_20_Source" = "LRT",
        "Venaza_condition" = c("LRT", "24hr_vs_control", "6hr_vs_control"),
        "Lagadinou_Treatment" = c("LRT","high_PTL_5uM_vs_high_no_drug", "low_no_drug_vs_high_no_drug", "low_PTL_5uM_vs_high_no_drug"),
        "Lee_prior_cr" = "LRT",
        default = character(0)
      )
      } else{
      choices <- switch(
        paste(dataset, sep = "_"),
        "BEAT_quantile" = c("q2_vs_q1", "q3_vs_q1", "q4_vs_q1"),
        "BEAT_FAB" = c("M0_vs_M5", "M1_vs_M5", "M3_vs_M5", "M4_vs_M5", "M5b_vs_M5"),
        "BEAT_Denovo.Relapse" = "Relapse_vs_Denovo",
        "TCGA_FAB" = c("M0_vs_M5", "M1_vs_M5", "M2_vs_M5", "M3_vs_M5", "M4_vs_M5", "M6_vs_M5", "M7_vs_M5"),
        "TCGA_RAS" = "wt_vs_mut",
        "TCGA_NPM1" = "wt_vs_mut",
        default = character(0)
      )
      }
    }
    
    observe({
      pwc_choices<- PWChoices(input$datainput, input$DEmodel)
      
      updateSelectizeInput(
        session,
        "pwc",
        choices = pwc_choices,
        selected = NULL
      )
    })
    

    
    user_choice <- reactive({
      input$datainput
    })
    
    user_model_choice <- reactive({
      input$DEmodel
    })
    
    user_PW_choice <- reactive({
      input$pwc
    })
  

    close_page <- reactive({
    
      input$runDE1
  
    })

    list(user_dataset = user_choice,
         user_model = user_model_choice,
         user_PW = user_PW_choice,
         close_tab = close_page
         )
  })
}



