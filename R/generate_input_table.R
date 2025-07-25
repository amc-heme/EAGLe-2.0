# Load dataset paths from datasets.yaml
datasets <- 
  read_yaml("./datasets.yaml")
raw_data_p <- getwd()
mouse <- read.table(file = paste0(raw_data_p,"/new_raw_data/Mm_TranscriptToGene_Orthologs.txt"), sep = "\t", header = T)
human <- read.table(file = paste0(raw_data_p,"/new_raw_data/Hs_TranscriptToGene_Orthologs.txt"), sep = "\t", header = T)
#load rds files from datasets.yaml
#extract variables needed for gene input drop down menu
dataset_identifiers <-
  lapply(
    # Iterate through keys (names) in the datasets list
    names(datasets),
    function(data_key){
      # Read RDS file using the path defined for the current dataset
      dataset <- read_rds(paste0(raw_data_p, datasets[[data_key]]$vst_path))
      
      #filter each dataset for only genes that contain human/mouse orthologs based on the ortholog gene LUT
      if(datasets[[data_key]]$species == "mouse") {
        dataset <- dataset %>% 
          dplyr::select("ensembl_gene_id") %>% 
          mutate(Species = "mouse") %>% 
          left_join(mouse %>% 
                      dplyr::filter(ensembl_gene_id %in% dataset$ensembl_gene_id) %>% 
                      dplyr::select(ensembl_gene_id, Gene, Ortholog), by = "ensembl_gene_id") 
        #print(length(unique(dataset$Gene)) == nrow(dataset))
      } else if(datasets[[data_key]]$species == "human") {
        dataset <- dataset %>% 
          #dplyr::filter(ensembl_gene_id %in% human$ensembl_gene_id) %>% 
          dplyr::select("ensembl_gene_id") %>%
          mutate(Species = "human") %>% 
          left_join(human %>% 
                      dplyr::filter(ensembl_gene_id %in% dataset$ensembl_gene_id) %>% 
                      dplyr::select(ensembl_gene_id, Gene, Ortholog), by = "ensembl_gene_id")
      }
      return(dataset)
    }
  )
# Add dataset names to list generated
names(dataset_identifiers) <- 
  names(datasets)
#combine mouse and human gene lists

dataset_identifiers <- bind_rows(dataset_identifiers, .id = "Dataset")

  
#save table
write.table(dataset_identifiers, paste0(raw_data_p, "/new_raw_data/input_lookup_table.txt"), quote = F, sep = "\t", row.names = F)
