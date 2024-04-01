
data_UI <- function(id) {
  ns <- NS(id)
  fluidPage(
    theme =
      shinytheme("flatly"),
    useShinyjs(),
    useWaiter(),
 
    titlePanel("Choose a dataset, model term, and comparison for DE and visualizations"
               ),
    # div(
    #   tags$li(
    #     class = "dataset_info_button",
    #     style = 'display:inline-block;',
    #     style = 'float:left;',
    #     dropMenu(
    #       dropdownButton(
    #         circle = TRUE, status = "info", icon = icon("info"),  size = 'sm', width = "50px",
    #         tooltip = tooltipOptions(title = "Click for more information on datasets")),
    #       h6("Summary of datasets:"),
    #       h6("Jordan M0-M5: "),
    #       h6("Pei et al. have previously shown that ROS-low enriched LSCs from primitive and monocytic AML 
    #          differ significantly in their transcriptome and metabolism. 
    #          ROS-low LSCs from monocytic AML in particular, are less dependent on BCL2 therefore are more likely to be resistant to Venetoclax-based therapy.
    #          This section of analysis shows gene expression in primitive vs. monocytic ROS-low LSCs."),
    #       h6("Ye et al. 2016: "),
    #       h6("LSC gene expression from blood, bone marrow, spleen, gonadal adipose tissue, and normal bone marrow in mice."),
    #       h6("Ye et al. 2020: "),
    #       h6("Transcriptomes of LSCs from bone marrow and liver in mice were compared to better understand 
    #          the biology of liver LSCs. A bcCML model (BCR-ABL + Nup98-Hoxa9) was used. 
    #          Liver and bone marrow samples were combined from 3 mice in 3 different cohorts before sequencing."),
    #       h6("Pollyea et al., Nature, 2018: "),
    #       h6("Time zero pheresis was taken from 3 patients, then collected again at 6hr and 24hr post ven/aza treatment."),
    #       h6( "Lagadinou et al., Cell Stem Cell, 2013: "),
    #       h6("BCL-2 inhibition targets oxidative phosphorylation and selectively eradicates quiescent human leukemia stem cells. 
    #          ROS high and ROS low LSCs were either treated with 5ul of PTL or not treated before sequencing."),
    #       h6("Lee et al., Nature, 2018: "),
    #       h6("Genome wide expression from 12 AML patient samples with either prior complete remission, or no prior complete remission."),
    #       h6("TCGA: "),
    #       h6("The Cancer Genome Atlas (TCGA) project published gene expression profiles of 151 primary AML patients 
    #          along with their mutational profiles and clinical characteristics. This section of analysis shows gene expression 
    #          in the TCGA-AML dataset parsed by various mutational/clinical variables including the French-American-British 
    #          subtypes, karyotype, RAS mutation status, and NPM1 mutation status."),
    #       h6("BEAT-AML: "),
    #       h6("The BEAT-AML project published gene expression profiles of ~400 primary AML patients 
    #       along with their mutational profiles and clinical characteristics. 
    #       This section of analysis shows gene expression in the BEAT-AML dataseit parsed by various mutational/clinical variables 
    #          including the French-American-British subtypes, Venetoclax response, and disease stage (de novo vs. relapse).")
    #       
    #     ))),
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
                         "Lee et al, 2018" = "Lee"), selected = "Cancer_Discovery"
          ),
          selectizeInput(
            ns("DEmodel"),
            label = "Choose a metadata variable from dataset",
            choices = NULL,
            selected = NULL ,
            options = list(maxItems = 1)
          ),
          selectizeInput(
            ns("pwc"),
            label = "Choose to run LRT OR a pairwise comparison for DE testing",
            choices = NULL,
            selected = NULL ,
            options = list(maxItems = 1)
          ),
          hr(),
         
          actionButton(
            (ns("runDE1")),
            "Continue"
          ),
        )
      ),
      mainPanel(
        textOutput(ns("Data_text"))
      )
    )
  )
}

data_Server <- function(id) {
  moduleServer(id, function(input, output, session){
    
    output$Data_text <- renderText({
      if (input$datainput == "Cancer_Discovery") {
        paste(
          "Pei et al. have previously shown that ROS-low enriched LSCs
             from primitive and monocytic AML differ significantly in their
             transcriptome and metabolism.
             ROS-low LSCs from monocytic AML in particular, are less dependent
             on BCL2 therefore are more likely to be resistant to
             Venetoclax-based therapy."
        )
      } else if (input$datainput == "Ye_16") {
        paste(
          "Comparison of ranscriptomes of LSC's from blood, bone marrow, spleen,
           gonadal adipose tissue, and normal bone marrow in mice."
        )
      } else if (input$datainput == "Ye_20") {
        paste(
          "Transcriptomes of LSCs from bone marrow and liver in mice
               were compared to better understand the biology of liver LSCs.
               A bcCML model (BCR-ABL + Nup98-Hoxa9) was used. Liver and bone
               marrow samples were combined from 3 mice in 3 different cohorts
               before sequencing."
        )
      } else if (input$datainput == "Venaza") {
        paste(
          "Time zero pheresis was taken from 3 patients, then collected
                 again at 6hr and 24hr post ven/aza treatment."
        )
      } else if (input$datainput == "Lagadinou") {
        paste(
          "BCL-2 inhibition targets oxidative phosphorylation and
               selectively eradicates quiescent human leukemia stem cells.
               ROS high and ROS low LSCs were either treated with 5ul of PTL
               or not treated before sequencing."
        )
      } else if (input$datainput == "TCGA") {
        paste(
          "The Cancer Genome Atlas (TCGA) project published gene expression
              profiles of 151 primary AML patients along with their mutational
              profiles and clinical characteristics. The TCGA-AML dataset parsed
              by various mutational/clinical variables including the 
              French-American-British subtypes, karyotype, RAS mutation status,
              and NPM1 mutation status."
        )
      } else if (input$datainput == "BEAT") {
        paste(
          "The BEAT-AML project published gene expression profiles of ~400
              primary AML patients along with their mutational profiles and clinical
              characteristics. The BEAT-AML dataset was parsed by various 
              mutational/clinical variables including the French-American-British
              subtypes, Venetoclax response, and disease stage (de novo vs. relapse)."
        )
      } else{
        paste(
          "Genome wide expression from 12 AML patient samples with either
              prior complete remission, or no prior complete remission."
        )
      }
    })
    # model choice ####
    DEModelChoices <- function(dataset) {
      m_choices <- switch(
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
    }
   

    observe({
      model_choices <- DEModelChoices(input$datainput)
      
      updateSelectizeInput(session,
                           "DEmodel",
                           choices = model_choices,
                           selected = NULL)
    })
    
    
    
    PWChoices <- function(dataset, model) { #this needs to not populate if model = "LRT"
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



