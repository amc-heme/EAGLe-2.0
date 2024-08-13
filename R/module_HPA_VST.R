# vst module

raw_data <- getwd()

HPAvst_Server <- function(id, dataset_choice) {
  moduleServer(id, function(input, output, session) {
 
    vst.HPA <- eventReactive(dataset_choice$close_tab(), {
      read_rds(paste0(raw_data, "/data/HPA_VST.rds"))
    })
    return(vst.HPA)
  })
}