qc.file_Server <- function(id, dataset.qc, dataset_choice) {
  moduleServer(id, function(input, output, session){
    
    qc_list <- dataset.qc
    
    qc_object <- reactive({
      qc_list[[dataset_choice$user_dataset()]]
    })
    
    observe({print(head(qc_object()))})
    return(qc_object)
  })
}