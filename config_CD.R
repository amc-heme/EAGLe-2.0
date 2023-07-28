# DE config 
config <- list(
  base_dir = "~/Documents/CD files/salmon",
  t2g_hs_file = "gmt_pathway_files/HS_transcript_to_gene.txt",
  metadata_file = "data_rds_files/CD_ROSlow_metadata.rds",
  samples = data.frame(SRR = metadata$SRR, batch = metadata$batch, condition = metadata$sample_type, sample_name = metadata$DESeq_Sample_Name),
  dds_file = "data_rds_files/CD_ROSlow_dds.rds",
  dds_res_file = "data_rds_files/CD_ROSlow_dds.res.rds",
  vsd_file = "data_rds_files/CD_ROSlow_vsd.rds",
  vst_file = "data_rds_files/CD_vst_batchcorrected.rds",
  vsd.pca_file = "data_rds_files/CD_ROSlow_vsd.pca.rds",
  vst.goi_file = "data_rds_files/CD_ROSlow_vst.goi.rds",
  batch = metadata$batch,
  qc = "data/multiqc_data.json",
  var_1 = samples$condition,
  var_2 = samples$batch,
  sample_id = c("SRR9265370", "SRR9265373", "SRR9265371", "SRR9265372", "SRR9265363", "SRR9265364", "SRR9265366", "SRR9265367", "SRR9265369", "SRR9265365", "SRR9265368", "SRR9265374"),
  design_formula = ~batch + condition
)
