# vst module

raw_data <- getwd()

HPAvst_Server <- function(id, dataset_choice, dds.HPA) {
  moduleServer(id, function(input, output, session) {
 
    vst.HPA <- eventReactive(dataset_choice$close_tab(), {
      precomputed_vst_table(dds.HPA())
    })
    return(vst.HPA)
  })
}