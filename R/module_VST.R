# vst module
#dataset yaml
datasets <- 
  read_yaml("./data.yaml")
raw_data <- getwd()
#read in t2g files for human and mouse
t2g_hs <- read_rds(paste0(raw_data, "/data/t2g_hs.rds"))
t2g_mm <- read_rds(paste0(raw_data, "/data/t2g_mm.rds"))

vst_Server <- function(id, dataset_dds, dataset_choice) {
  moduleServer(id, function(input, output, session) {
    
    # function for creating vst file
    vst.create <- function(dds, dataset){
      
      dds.file <- dds
      
      vsd <-
        vst(dds.file, blind = F)
      
      #mouse or human?
      is_hs <- grepl("t2g_hs", datasets[[dataset]]$t2g)
      
      if (is_hs &
          dataset %in% c("Cancer_Discovery", "Venaza", "Lagadinou", "Lee")) {
        
        vst.table <- data.frame(assay(vsd), check.names = FALSE) %>% 
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
      } else if (is_hs & dataset %in% c("BEAT_quantile", "BEAT_FAB", "BEAT_Denovo.Relapse", "TCGA_FAB", "TCGA_NPM1", "TCGA_RAS")) {
        vst.table <- data.frame(assay(vsd), check.names = FALSE) %>%
          janitor::clean_names()  %>%
          rownames_to_column(., var = "ensembl_gene_id") %>% 
          mutate(ensembl_gene_id = sub("\\..*", "", ensembl_gene_id)) %>% #remove anything after the decimal to match with t2g ensembl ids
          left_join(unique(dplyr::select(t2g_hs, c(
            ensembl_gene_id, ext_gene
          ))), ., by = 'ensembl_gene_id') %>%
          dplyr::mutate(ext_gene_ensembl = case_when(ext_gene == "" ~ ensembl_gene_id,
                                                     TRUE ~ ext_gene)) %>%
          dplyr::select(-ext_gene) %>%
          dplyr::select(ext_gene_ensembl, everything()) %>%
          na.omit(.)
      } else {
        vst.table <- data.frame(assay(vsd), check.names = FALSE) %>%
          janitor::clean_names()  %>%
          rownames_to_column(., var = "ensembl_gene_id") %>%
          left_join(unique(dplyr::select(t2g_mm, c(
            ensembl_gene_id, ext_gene
          ))), .,
          by = 'ensembl_gene_id') %>%
          dplyr::mutate(ext_gene_ensembl = case_when(ext_gene == "" ~ ensembl_gene_id,
                                                     TRUE ~ ext_gene)) %>%
          dplyr::select(-ext_gene) %>%
          dplyr::select(ext_gene_ensembl, everything()) %>%
          na.omit(.)
      }
      vst.table
    }
    vst.global <- reactive({
      vst.create(dataset_dds(), dataset_choice$user_dataset())
    })
    return(vst.global)
  })
}