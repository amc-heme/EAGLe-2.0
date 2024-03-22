# vst module
datasets <- 
  read_yaml("./data.yaml")
raw_data <- getwd()
t2g_hs <- read_rds("~/Documents/GitHub/EAGLe-2.0/data/t2g_hs.rds")
t2g_mm <- read_rds("~/Documents/GitHub/EAGLe-2.0/data/t2g_mm.rds")

vst_Server <- function(id, dataset_dds, dataset_choice) {
  moduleServer(id, function(input, output, session) {
    
    vst.create <- function(dds, dataset){
      
      dds.file <- dds
      
      vsd <- 
        vst(dds.file, blind = F)
      
      #mouse or human?
      is_hs <- grepl("t2g_hs", datasets[[dataset]]$t2g)
      
      if(is_hs & dataset %in% c("Cancer_Discovery", "Venaza", "Lagadinou", "Lee")){
        
        vst.table <- data.frame(assay(vsd), check.names = FALSE) %>% 
          rownames_to_column(., var = "ensembl_gene_id") %>% 
          left_join(unique(dplyr::select(t2g_hs, c(ensembl_gene_id, ext_gene))), ., by = 'ensembl_gene_id') %>%
          na.omit(.)
      } else if(is_hs & dataset %in% c("BEAT", "TCGA")) {
        
          vst.table <- data.frame(assay(vsd), check.names = FALSE) %>% 
            rownames_to_column(., var = "ensembl_gene_id") %>% 
            mutate(ensembl_gene_id = sub("\\..*", "", ensembl_gene_id)) %>%
            left_join(unique(dplyr::select(t2g_hs, c(ensembl_gene_id, ext_gene))), ., by = 'ensembl_gene_id') %>%
            na.omit(.) 
      } else {
        
        vst.table <- data.frame(assay(vsd), check.names = FALSE) %>% 
          rownames_to_column(., var = "ensembl_gene_id") %>% 
          left_join(unique(dplyr::select(t2g_mm, c(ensembl_gene_id, ext_gene))), ., by = 'ensembl_gene_id') %>%
          na.omit(.)
      }
      print("vst.table:")
      print(head(vst.table))
      vst.table
    }
    vst.global <- reactive({
      vst.create(dataset_dds(), dataset_choice$user_dataset())
    })
    return(vst.global)
  })
}