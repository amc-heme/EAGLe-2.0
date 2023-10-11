dds.file_Server <- function(id, dataset, dataset_choice) {
  moduleServer(id, function(input, output, session){
    
    data_list <- dataset
    
      dds_object <- reactive({
        dataset_dds <- data_list[[dataset_choice()]]
        dataset_dds
      })
      
      return(dds_object)
  })
}
