dds.file_Server <- function(id, dataset, dataset_choice) {
  moduleServer(id, function(input, output, session){
    
    data_list <- dataset
    
      dds_object <- reactive({
        data_list[[dataset_choice$user_dataset()]]
      })
    
      return(dds_object)
  })
}
