# VST table with condition as column names for use in heatmaps

#dataset yaml
# datasets <- 
#   read_yaml("./data.yaml")
# raw_data <- getwd()
# #read in t2g files for human and mouse
# # t2g_hs <- read_rds("~/Documents/GitHub/EAGLe-2.0/data/t2g_hs.rds")
# # t2g_mm <- read_rds("~/Documents/GitHub/EAGLe-2.0/data/t2g_mm.rds")
# 
# vsthm_Server <- function(id, data_species, dataset_dds, dataset_choice) {
#   moduleServer(id, function(input, output, session) {
#     
#     t2g_create <- function(dataset_species) {
#       if(dataset_species == "human"){
#         t2gene <- read_rds("~/Documents/GitHub/EAGLe-2.0/data/t2g_hs.rds")
#       }else {
#         t2gene <- read_rds("~/Documents/GitHub/EAGLe-2.0/data/t2g_mm.rds")
#       }
#       return(t2gene)
#     }
#   
#     
#     # function for creating vst file
#     vst.create.hm <- function(dds, dataset, dataset_model){
#       
#       t2g <- t2g_create(data_species())
#       print("t2g:")
#       print(head(t2g))
#       
#       dds.file <- dds
#       
#       vsd <-
#         vst(dds.file, blind = F)
#       
#       #mouse or human?
#       is_hs <- grepl("t2g_hs", datasets[[dataset]]$t2g)
#       
#       if (is_hs &
#           dataset %in% c("Cancer_Discovery", "Ye_20", "Lee")) {
#         
#         meta <- colData(dds.file)
#         cond_var <- datasets[[dataset]]$PCA_var
#         cond <- meta[, cond_var]
#         print("cond")
#         print(cond)
#         vst.table <- data.frame(assay(vsd), check.names = FALSE) %>% 
#           janitor::clean_names()  %>%
#           rownames_to_column(., var = "ensembl_gene_id") %>%
#           left_join(unique(dplyr::select(t2g, c(
#             ensembl_gene_id, ext_gene
#           ))), ., by = 'ensembl_gene_id') %>%
#           dplyr::mutate(ext_gene_ensembl = case_when(ext_gene == "" ~ ensembl_gene_id,
#                                                      TRUE ~ ext_gene)) %>%  #if any blank ext_gene name, add ensembl gene id instead
#           dplyr::select(-ext_gene) %>%
#           dplyr::select(ext_gene_ensembl, everything()) %>%
#           na.omit(.)
#         
#         
#         
#         colnames(vst.table)[3:ncol(vst.table)] <- cond
#         
#       } else if (is_hs & dataset %in% c("BEAT", "TCGA")) {
#         vst.table <- data.frame(assay(vsd), check.names = FALSE) %>%
#           janitor::clean_names()  %>%
#           rownames_to_column(., var = "ensembl_gene_id") %>% 
#           mutate(ensembl_gene_id = sub("\\..*", "", ensembl_gene_id)) %>% #remove anything after the decimal to match with t2g ensembl ids
#           left_join(unique(dplyr::select(t2g, c(
#             ensembl_gene_id, ext_gene
#           ))), ., by = 'ensembl_gene_id') %>%
#           dplyr::mutate(ext_gene_ensembl = case_when(ext_gene == "" ~ ensembl_gene_id,
#                                                      TRUE ~ ext_gene)) %>%
#           dplyr::select(-ext_gene) %>%
#           dplyr::select(ext_gene_ensembl, everything()) %>%
#           na.omit(.)
#         
#         meta <- colData(dds.file)
#         cond_var <- dataset_model
#         cond <- meta[, cond_var]
#         
#         colnames(vst.table)[3:ncol(vst.table)] <- cond
#         
#       } else{
#         
#         vst.table <- data.frame(assay(vsd), check.names = FALSE) %>% 
#           janitor::clean_names()  %>%
#           rownames_to_column(., var = "ensembl_gene_id") %>%
#           left_join(unique(dplyr::select(t2g, c(
#             ensembl_gene_id, ext_gene
#           ))), ., by = 'ensembl_gene_id') %>%
#           dplyr::mutate(ext_gene_ensembl = case_when(ext_gene == "" ~ ensembl_gene_id,
#                                                      TRUE ~ ext_gene)) %>%  #if any blank ext_gene name, add ensembl gene id instead
#           dplyr::select(-ext_gene) %>%
#           dplyr::select(ext_gene_ensembl, everything()) %>%
#           na.omit(.)
#         
#         meta <- colData(dds.file)
#         cond_var <- dataset_model
#         cond <- meta[, cond_var]
#         
#         colnames(vst.table)[3:ncol(vst.table)] <- cond
#       } 
#       print("vst.table.hm:")
#       print(head(vst.table))
#       vst.table
#     }
#     vst.hm <- reactive({
#       vst.create.hm(dataset_dds(), dataset_choice$user_dataset(), dataset_choice$user_model())
#     })
#     return(vst.hm)
#   })
# }