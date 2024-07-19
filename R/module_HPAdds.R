HPAdds.file_Server <- function(id, dataset) {
  moduleServer(id, function(input, output, session){

    data_list <- dataset
    
    HPAdds_object <- reactive({
      data_list[["HPA"]]
    })
    
    return(HPAdds_object)
  })
}
