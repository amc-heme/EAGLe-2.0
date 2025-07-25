
dataset_config2 <- 
  read_yaml("./datasets.yaml")
# Extract vector of human-readable labels from the config file
dataset_labels <-
  sapply(
    dataset_config2,
    function(dataset){
      dataset$label
    }
  ) |> 
  unname()

# Construct vector of dataset choices using config file and above information
# (human-readable labels as labels, machine-readable labels as values)
dataset_choices <- names(dataset_config2)
names(dataset_choices) <- dataset_labels


geneModuleUI <- function(id) {
  ns <- NS(id)
  tagList(
    # Simple CSS style sheet for components
    includeCSS("./style.css"),
    useShinyjs(),
    
    mainPanel(
      h1("Shiny EAGLe v7.0."),
      h4(strong("Please read before use:")),
      dropMenu(
        dropdownButton(circle = TRUE, status = 'info', icon = icon('info'), size = 'sm',
                       right = T, width = '50px',
                       tooltip = tooltipOptions(title = "Information")),
        
        h6("1. To search for a gene symbol:"),
        h6("- Click on the gene symbol menu, hit the delete (mac) or backspace (windows) key, and begin typing."),
        h6("- Genes are searchable by both gene symbol and ensembl gene id. They are shown with their corresponding species and ortholog."),
        h6("- A star next to a dataset indicates that the gene was tested for differential expression while datasets without a star indicate that expression of the gene was too low to pass the filter for testing."),
        br(),
        h6("2. If a gene symbol does not appear in the autocomplete list, it does not appear in any of the datasets."),
        
        br(),
        h6("3. A note about datasets:"),
        h6("- If a gene does not appear in a dataset, or that dataset is deselected, that dataset's section will also not appear in the report."),
        h6("- All genes are labelled as being either a mouse or human gene. The report generated will only
           include datasets involving the chosen species.")
        
      ),
      br(),
      selectizeInput(
        inputId = ns("gene"), 
        label = "Enter a Gene Symbol:", 
        choices = NULL, 
        selected = NULL, 
        multiple = FALSE
      ),
      tags$head( #css to move the pop up window to the top and middle of the screen
        tags$style(
          HTML(".shiny-notification {
             position:fixed;
             top: calc(10%);
             left: calc(40%);
             }
             "
          )
        )
      ),
      hidden(
        div(
          id = ns("available_datasets_ui"),
          # Shows the datsets that contain a gene, once a gene is entered
          div(
            tags$h4(
              "Available datasets for selected gene",
              class = "center container-header",
              style = "display:inline-block;"
            ),
            tags$li(
              class = "dataset_info_button",
              style = 'display:inline-block;',
              style = 'float:right;',
              dropMenu(
                dropdownButton(
                  circle = TRUE, status = "info", icon = icon("info"),  size = 'sm', width = "50px",
                  tooltip = tooltipOptions(title = "Click for more information on datasets")),
                h6("Summary of datasets:"),
                h6("Jordan M0-M5: "),
                h6("Pei et al. have previously shown that ROS-low enriched LSCs from primitive and monocytic AML 
             differ significantly in their transcriptome and metabolism. 
             ROS-low LSCs from monocytic AML in particular, are less dependent on BCL2 therefore are more likely to be resistant to Venetoclax-based therapy.
             This section of analysis shows gene expression in primitive vs. monocytic ROS-low LSCs."),
                h6("Ye et al. 2016: "),
                h6("LSC gene expression from blood, bone marrow, spleen, gonadal adipose tissue, and normal bone marrow in mice."),
                h6("Ye et al. 2020: "),
                h6("Transcriptomes of LSCs from bone marrow and liver in mice were compared to better understand 
             the biology of liver LSCs. A bcCML model (BCR-ABL + Nup98-Hoxa9) was used. 
             Liver and bone marrow samples were combined from 3 mice in 3 different cohorts before sequencing."),
                h6("Pollyea et al., Nature, 2018: "),
                h6("Time zero pheresis was taken from 3 patients, then collected again at 6hr and 24hr post ven/aza treatment."),
                h6( "Lagadinou et al., Cell Stem Cell, 2013: "),
                h6("BCL-2 inhibition targets oxidative phosphorylation and selectively eradicates quiescent human leukemia stem cells. 
             ROS high and ROS low LSCs were either treated with 5ul of PTL or not treated before sequencing."),
                h6("Lee et al., Nature, 2018: "),
                h6("Genome wide expression from 12 AML patient samples with either prior complete remission, or no prior complete remission."),
                h6("TCGA: "),
                h6("The Cancer Genome Atlas (TCGA) project published gene expression profiles of 151 primary AML patients 
             along with their mutational profiles and clinical characteristics. This section of analysis shows gene expression 
             in the TCGA-AML dataset parsed by various mutational/clinical variables including the French-American-British 
             subtypes, karyotype, RAS mutation status, and NPM1 mutation status."),
                h6("BEAT-AML: "),
                h6("The BEAT-AML project published gene expression profiles of ~400 primary AML patients 
          along with their mutational profiles and clinical characteristics. 
          This section of analysis shows gene expression in the BEAT-AML dataset parsed by various mutational/clinical variables 
             including the French-American-British subtypes, Venetoclax response, and disease stage (de novo vs. relapse).")
                
              ))),
          
          uiOutput(
            outputId = ns("available_datasets")
          ),
          # Menu to select datasets to include in the report
          tags$b("Choose datasets to include in report:"),
          tags$p("(all available datasets selected by default)"),
          shinyWidgets::pickerInput(
            inputId = ns("select_datasets"),
            label = NULL,
            choices = dataset_choices,
            selected = character(0),
            multiple = TRUE
          )
        )
      ),
      downloadButton(ns("report"), "Generate PDF Report"),
   
    )
  )
}

geneModuleServer <- function(id, gene_input, genes_by_dataset, 
                             DEgenes_by_dataset, raw_data_p) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    #named vector of identifiers for drop down list (valid_inputs) and associated values (ensembl_gene_id)
    gene_choices <- setNames(gene_input$ensembl_gene_id, paste(
      gene_input$Gene,
      "|",
      gene_input$ensembl_gene_id,
      "|",
      gene_input$Species,
      "|",
      gene_input$Ortholog
      ))
    
    observe({
      updateSelectizeInput(
        session = session, 
        inputId = 'gene', 
        choices = c("", gene_choices),
        selected = NULL, 
        server = TRUE
      )
    })
    
    
    # 2. Preview/Select Datasets in Report -----------------------------------
    ## 2.1. Show datasets that have the selected gene ####
    gene_present <-
      reactive(
        label = "Show available datasets for gene",
        {
          # Avoids errors in the event the gene is not defined
          req(input$gene)
          
          # Construct boolean vector based on whether or not the selected
          # gene is in the gene list for each dataset
          gene_present <-
            sapply(
              names(genes_by_dataset),
              # Function takes a second parameter, gene. input$gene is passed
              # to this using the third sapply parameter
              function(dataset_key, gene){
                gene %in% genes_by_dataset[[dataset_key]]
              },
              input$gene
            )
          
          names(gene_present) <- 
            names(genes_by_dataset)
          
          gene_present
        })
    
    ## 2.2. show datasets that have gene in DE res table ####
    
    gene_DE <-
      reactive(
        label = "Show datasets where gene appears in the DE results",
        {
          # Avoids errors in the event the gene is not defined
          req(input$gene)
          
          # Construct boolean vector based on whether or not the selected
          # gene is in the gene list for each dataset
          gene_DE <-
            sapply(
              names(DEgenes_by_dataset),
              # Function takes a second parameter, gene. input$gene is passed
              # to this using the third sapply parameter
              function(dataset_key, gene){
                gene %in% DEgenes_by_dataset[[dataset_key]]
              },
              input$gene
            )
          
          names(gene_DE) <- 
            names(DEgenes_by_dataset)
          
          gene_DE
        })
    
    ## 2.3. UI to display datasets with the gene ####
    output$available_datasets <-
      renderUI({
        req(c(gene_present(), gene_DE()))
        
        # For each dataset, display label and whether gene is present
        # on the same line
        tagList(
          
          div(
            class = "container-body dataset-container",
            lapply(
              names(dataset_config2),
              function(data_key, dataset_config2){
                div(
                  class = "dataset-row",
                  # add unicode mouse icon for mouse datasets
                  if(dataset_config2[[data_key]]$species == "mouse") {
                    HTML(paste0(
                      '<span>&#x1F42D;</span>',
                      ' ',
                      tags$b(
                        dataset_config2[[data_key]]$label, 
                        style = "display: inline-block"
                      )
                    ))
                  } else {
                    # Label for dataset
                    tags$b(
                      dataset_config2[[data_key]]$label, 
                      style = "display: inline-block"
                    )
                  },
                  # Indicator: check if gene is present, "slash" if not
                  if (gene_present()[data_key] == TRUE){
                    tags$p(
                      icon(
                        name = "check-circle",
                        # "fas" class gives circle with solid background
                        # https://fontawesome.com/v5/icons/check-circle?f=classic&s=solid
                        class = "yes fas"
                      ),
                      style = "display: inline-block"
                    )
                  } else{
                    tags$p(
                      icon(
                        name = "minus-circle",
                        class = "no"
                      ),
                      style = "display: inline-block"
                    )
                  },
                  if(gene_DE()[data_key] == TRUE) {
                    print("star")
                    tags$p(
                      icon(
                        name = "star",
                        # "fas" class gives star with solid background
                        # https://fontawesome.com/v5/icons/star?f=classic&s=solid
                        class = "fas fa-star"
                      ),
                      style = "display: inline-block"
                    )
                  } else{
                    NULL
                  }
                )
              },
              dataset_config2
            )
          )
        )
      })
    
    ## 2.4. Show/hide available dataset display for gene ####
    # Shows/hides datasets using an animation
    observe({
      target_id <- "available_datasets_ui"
      
      if (isTruthy(input$gene)){
        print("show element")
        showElement(
          id = target_id,
          anim = TRUE
        )
      } else {
        print("Hide element")
        hideElement(
          id = target_id,
          anim = TRUE
        )
      }
    })
    #create a new function that looks for genes that are present, that are also in the dds.res table.
    #if they are not in both, need to include another symbol next to datasets (yellow) to alert
    #the user that the gene was not tested and no padj values will be shown
    
    ## 2.5. Select datasets to include ####
    observe({
      req(gene_present())
      
      # Select all datasets available for a gene, and disable datasets that 
      # are not available
      updatePickerInput(
        session = session,
        inputId = "select_datasets",
        choices = dataset_choices,
        selected = dataset_choices[gene_present() == TRUE],
        choicesOpt = 
          list(
            # `disabled` takes a boolean vector of length equal to the number
            # of choices. When an entry is TRUE, the choice at the corresponding
            # index is disabled.
            # gene_present() is already a boolean vector
            disabled = !gene_present(),
            # CSS style to apply to disabled choices (gray out choice)
            style = ifelse(
              !gene_present(),
              yes = "color: #88888888;",
              no = ""
            )
          ),
        options =
          list(
            "selected-text-format" = "count > 5",
            "actions-box" = TRUE,
            # # Placeholder
            "none-selected-text" = "No Datasets Selected"
          )
      )
    })
    
    #showNotification popup when the user selects a gene that has multiple ensembl gene id's
    observeEvent(input$gene,{
      #object for the gene symbol associated with the ensembl id chosen in the dropdown menu
      find_gene <- gene_input$Gene[gene_input$ensembl_gene_id == input$gene]
      
      total_count <- 0
      # loop for counting the total number of times an ensembl id is found in all of the datasets
      for(dataset_key in names(genes_by_dataset)) {
        count_for_ensembl <- sum(input$gene %in% genes_by_dataset[[dataset_key]])
        total_count <- total_count + count_for_ensembl
      }
      # how many times the gene symbol occurs 
      count_for_symbol <- sum(gene_input$Gene %in% find_gene)
      #need to see if the number of times the gene symbol occurs is > then the number of datasets in which it is found
      if(total_count < count_for_symbol) {
        
        showNotification("The selected gene has multiple loci, therefore there are more than one ensembl gene id for this gene symbol. Please make sure that the ensembl
                          id selected matches the locus of interest.", type = "message", duration = 15)
      }
      # print(total_count)
      # print(count_for_symbol)
    })
    # 3. Download Report -----------------------------------------------------
    output$report <- 
      downloadHandler(
        filename = 
          function() {
            paste0(input$gene, "_EAGLe_Analysis_v7.0.pdf")
          },
        content = 
          function(file) {
            #creates a loading spinner while pdf renders
            showModal(modalDialog("Loading", footer=NULL))
            on.exit(removeModal())
            
            # Copy the report file to a temporary directory before processing 
            # it, in case we don't have write permissions to the current 
            # working dir (which can happen when deployed).
            tempReport <- 
              file.path(tempdir(), "EAGLe_Analysis_v7.0.Rmd")
            
            file.copy(
              "EAGLe_Analysis_v7.0.Rmd", 
              tempReport, 
              overwrite = TRUE
            )
            #object used to display the gene symbol for the corresponding ensembl id selected in the drop down menu
            gene_choices_symbol <- gene_input$Gene[gene_input$ensembl_gene_id == input$gene][1]  #only show the first instance since the gene will have duplicates for each dataset
            # Set up parameters to pass to Rmd document
            params <- 
              list(
                gene = input$gene,
                gene_symbol = gene_choices_symbol,
                raw_data = raw_data_p,
                gene_present = gene_present(),
                gene_DE = gene_DE(),
                selected_datasets = input$select_datasets,
                datasets = dataset_config2
              )
            
            
            # Set up pandoc
            Sys.setenv(
              RSTUDIO_PANDOC = 
                "/Applications/RStudio.app/Contents/MacOS/quarto/bin/tools"
            )
            
            # Knit the document, passing in the `params` list, and eval it in a
            # child of the global environment (this isolates the code in the 
            # document from the code in this app).
            rmarkdown::render(
              tempReport,
              output_file = file,
              params = params,
              envir = new.env(parent = globalenv())
            )
          }
      )
  })
}