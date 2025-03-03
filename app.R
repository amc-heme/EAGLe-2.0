

options(
  shiny.fullstacktrace = TRUE
)
## Data ####
# Config information for each dataset
dataset_config <- 
  read_yaml("./data.yaml")

raw_data_p <- getwd()

#load the datasets
load_datasets <- function(raw_data_p, dataset_config) {
  dataset_keys <- names(dataset_config)
  
  dataset_path <- lapply(dataset_keys, function(key) {
    read_rds(paste0(raw_data_p, dataset_config[[key]]$data_path))
  })
  return(dataset_path)
}

dataset <- load_datasets(raw_data_p, dataset_config)

#load the multiqc files
load_qc <- function(raw_data_p, dataset_config) {
  dataset_keys <- names(dataset_config)
  
  dataqc_path <- lapply(dataset_keys, function(key) {
    read_rds(paste0(raw_data_p, dataset_config[[key]]$qc_path))
  })
  return(dataqc_path)
}

dataset.qc <- load_qc(raw_data_p, dataset_config)

 # Add dataset names to list generated
names(dataset) <-
  names(dataset_config)

# Add dataset names to qc list generated
names(dataset.qc) <-
  names(dataset_config)

# UI ####
ui <-
  navbarPage("EAGLe2",
              tabsetPanel(
               id = "page",
               type = "hidden",
             # landing page UI
               tabPanelBody("landing_page",
                            waiter::useWaiter(),
                            data_UI("data1")),
            # app content 
               tabPanelBody("content",
                            # button to return to landing page from main app
                            actionButton(
                              "change_data",
                              icon = icon("refresh"),
                              class = "btn-xs",
                              label = "Change Dataset",
                              style = "position: absolute; right: 40px"
                            ),      
                              #QC Menu ####
                              tabsetPanel(
                                id = "tabs",
                              
                                tabPanel("QC",
                                      QC_UI("QC1")),
                              #DESeq Menu ####
                              tabPanel("Differential Expression",
                                      DE_UI("DEtab1")),
                              
                              #GSEA menu ####
                              tabPanel("GSEA",
                                      GSEA_UI("GSEA1")),
                              
                              #Gene expression analysis ####
                              tabPanel("Gene Expression",
                                      goi_UI("GOI1")),
                              tabPanel("Normal Tissue",
                                       HPA_UI("HPA1"))
                            ))
             ))
 
#Server ####
server <- 
  function(input, output, session) {

    options(shiny.reactlog = TRUE)
    
    #reactive container for reset button(change dataset action button)
    reset_trigger <- reactive({
      input$change_data
      })
    
  ##Data tab ####
    dataset_choice <- data_Server("data1")
    
    #hide GSEA tab if LRT is chosen in place of a pairwise comparison
    observe(
      if(dataset_choice$user_PW() == "LRT" & dataset_choice$user_dataset() %in% 
         c("Ye_16", "Venaza","Lagadinou", "TCGA", "BEAT")) {
        hideTab(inputId = "tabs", target = "GSEA")
      } else {
        showTab(inputId = "tabs", target = "GSEA")
      }
    )
    
   # close landing page and open app 
    observeEvent(dataset_choice$close_tab(), {
      updateTabsetPanel(session, "page", "content")
    })
    
    #close main app and return to landing page
    observeEvent(input$change_data, {
      updateTabsetPanel(session, "page", "landing_page")
    })

    #reactive statement to return dataset species
    data_species <- reactive({
      dataset_config[[dataset_choice$user_dataset()]]$species
    })
 
  ## dds object
  #dataset = links path to chosen dataset DESeq object
  #dataset_choice = user selected dataset from data server
    dataset_dds <- dds.file_Server("dds1", dataset, dataset_choice)
  #   
  #   #dds file for HPA needs to be made each time the app loads to show the normal tissue plot
     dds.HPA <- HPAdds.file_Server("HPAdds1", dataset)
  # ## vst table
  #   #dataset_dds = DEseq object returned by dds.file server
  #   #dataset_choice = user selected dataset from data server
     vst <- vst_Server("vst1", dataset_dds, dataset_choice)

     vst.HPA <- HPAvst_Server("HPAvst1", dataset_choice)

  # ## qc objects 
  #   #dataset.qc = opens path to stored qc file for chosen dataset
  #   #dataset_choice = user selected dataset from data server
     qc_table <- qc.file_Server("qct1", dataset.qc, dataset_choice)
  
  # ## QC tab ####
  #   #dataset_dds = DEseq object returned by dds.file server
  #   #dataset_choice = user selected dataset from data server
  #   #qc_table = stored multiqc table
  #   # reset trigger = clears all previous selections and returns to landing page
    QC_Server("QC1", dataset_dds, dataset_choice, qc_table, reset_trigger)

  # ## GOI tab####
  #   #dataset_dds = DEseq object returned by dds.file server
  #   #dataset_choice = user selected dataset from data server
  #   #vst = vst table 
     goi_Server("GOI1", dataset_choice, dataset_dds, vst)

  # # ##DESEq #####
  #   # data_species = needed to select appropiate t2g table
  #   #dataset_dds = DEseq object returned by dds.file server
  #   #dataset_choice = user selected dataset from data server
  #   #reset trigger = clears all previous selections and returns to landing page
  #   #vst = vst table 
     DE_res <- DE_Server("DEtab1", data_species, dataset_dds, dataset_choice, reset_trigger, vst, vst_hm)

  # # ##GSEA output ####
  #   #dataset_choice = user selected dataset from data server
  #   #DE_res = DE server returns the DE results table in tidy format 
  #   #reset trigger = clears all previous selections and returns to landing page
  #   #vst = vst table 
     GSEA_Server("GSEA1", dataset_choice, DE_res, reset_trigger, vst, dataset_dds)

  # # ##Human Protein Atlas- normal tissue tab ####
  #   # dds.HPA = HPA dds file 
  #   # vst.HPA = vst table created for HPA dataset
     HPA_Server("HPA1", vst.HPA, dds.HPA)

    #returns user to QC tab after switching datasets
    observeEvent(input$change_data, {
      updateTabsetPanel(session, "tabs", selected = "QC")
    })

  } #end server

# Run the application 
shinyApp(ui = ui, server = server)
