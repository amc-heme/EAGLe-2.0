# vst module
#dataset yaml
datasets <- 
  read_yaml("./data.yaml")
raw_data <- getwd()
#read in t2g files for human and mouse
t2g_hs <- read_rds("~/Documents/GitHub/EAGLe-2.0/data/t2g_hs.rds")
t2g_mm <- read_rds("~/Documents/GitHub/EAGLe-2.0/data/t2g_mm.rds")

HPAvst_Server <- function(id, dds.HPA) {
  moduleServer(id, function(input, output, session) {
    # function for creating vst file
    vst.create <- function(dds, dataset){
      
      dds.file <- dds
      
      vsd <-
        vst(dds.file, blind = F)
      
        HPAvst.table <- 
          data.frame(assay(vsd), check.names = FALSE) %>% 
          janitor::clean_names()  %>%
          rownames_to_column(., var = "ensembl_gene_id") %>%
          left_join(unique(dplyr::select(t2g_hs, c(
            ensembl_gene_id, ext_gene
          ))), ., by = 'ensembl_gene_id') %>%
          dplyr::mutate(ext_gene_ensembl = case_when(ext_gene == "" ~ ensembl_gene_id,
                                                     TRUE ~ ext_gene)) %>%  #if any blank ext_gene name, add ensembl gene id instead
          dplyr::select(-ext_gene) %>%
          dplyr::select(ext_gene_ensembl, everything()) %>%
          na.omit(.)
     
        HPAvst.table
    }
    vst.HPA <- reactive({
      vst.create(dds.HPA(), datasets[["HPA"]])
    })
    return(vst.HPA)
  })
}