
data_UI <- function(id) {
  ns <- NS(id)
  fluidPage(
    theme =
      shinytheme("flatly"),
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
                         "TCGA-LAML" = "TCGA",
                         "BEAT-AML" = "BEAT",
                         "Lee et al, 2018" = "Lee",
                         "Human Protein Atlas" = "HPA"), selected = "Cancer_Discovery"
          ),
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
        )
      ),
      mainPanel(
        
      )
    )
  )
}

data_Server <- function(id, dataset_choice) {
  moduleServer(id, function(input, output, session){
    
    # model choice ####
    DEModelChoices <- function(session, dataset) {
      choices <- switch(
        dataset,
        "Cancer_Discovery" = "LRT",
        "Ye_16" = "Source",
        "Ye_20" = "LRT",
        "Venaza" = "condition",
        "Lagadinou" = "Treatment",
        "BEAT" = c("quantile", "FAB_BlastMorphology", "Denovo.Relapse"),
        "TCGA" = c("FAB", "RAS_mut", "NPM1_mut"),
        "Lee" = "LRT",
        default = character(0)
      )
      
      updateSelectizeInput(session,
                           "DEmodel",
                           choices = choices,
                           selected = NULL)
    }
    
    PWChoices <- function(session, dataset, model) {
      choices <- switch(
        paste(dataset, model, sep = "_"),
        "Cancer_Discovery_LRT" = "LRT",
        "Ye_16_Source" = c("LRT", "blood_vs_bone_marrow", "gonadal_adipose_tissue_vs_bone_marrow", "normal_bm_vs_bone_marrow", "spleen_vs_bone_marrow"),
        "Ye_20_LRT" = "LRT",
        "Venaza_condition" = c("LRT", "24hr_vs_control", "6hr_vs_control"),
        "Lagadinou_Treatment" = c("LRT","high_PTL_5uM_vs_high_no_drug", "low_no_drug_vs_high_no_drug", "low_PTL_5uM_vs_high_no_drug"),
        "BEAT_quantile" = c("q2_vs_q1", "q3_vs_q1", "q4_vs_q1"),
        "BEAT_FAB_BlastMorphology" = c("M0_vs_M5", "M1_vs_M5", "M3_vs_M5", "M4_vs_M5", "M5b_vs_M5"),
        "BEAT_Denovo.Relapse" = "Relapse_vs_Denovo",
        "TCGA_FAB" = c("M0_vs_M5", "M1_vs_M5", "M2_vs_M5", "M3_vs_M5", "M4_vs_M5", "M6_vs_M5", "M7_vs_M5"),
        "TCGA_RAS_mut" = "wt_vs_mut",
        "TCGA_NPM1_mut" = "wt_vs_mut",
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
    
    user_model_choice <- observe({
      DEModelChoices(session, input$datainput)
    })

    user_PW_choice <- observe({
      PWChoices(session, input$datainput, input$DEmodel)
    })
    
    user_choice <- reactive({
      input$datainput
    })
    
    observe({
      print(user_choice())
    })
    
    list(user_dataset = user_choice,
         user_model = user_model_choice,
         user_PW = user_PW_choice
         )
  })
}



