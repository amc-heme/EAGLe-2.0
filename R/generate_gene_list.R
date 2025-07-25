# Load dataset paths from datasets.yaml
datasets <- 
  read_yaml("./datasets.yaml")
raw_data_p <- getwd()
#mouse_human <- read.table(file = paste0(raw_data_p,"/new_raw_data/input_lookup_table.txt"), sep = "\t", header = T)
# Generate list of genes in each dataset
dataset_genes <-
  lapply(
    # Iterate through keys (names) in the datasets list
    names(datasets),
    function(data_key){
      # Read RDS file using the path defined for the current dataset
      dataset <- read_rds(paste0(raw_data_p, datasets[[data_key]]$vst_path))
   
      # Extract/return gene names
      dataset$ensembl_gene_id
    }
  )

# Add dataset names to list generated
names(dataset_genes) <- 
  names(datasets)

# Save list to file
# Open connection, write comments
con <- file("./genes_by_dataset.yaml", open = "w")
cat(
  "# Available genes by dataset",
  "# DO NOT MODIFY BY HAND!",
  "# Use generate_gene_list.R to modify.",
  sep = "\n",
  file = con
  )
# Write data, close connection
write_yaml(
  dataset_genes,
  con
)
close(con)
