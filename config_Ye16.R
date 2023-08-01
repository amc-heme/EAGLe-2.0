config_ye16 <- list(
  t2g_file = "data_rds_files/t2g_mm.rds",
  metadata_file = "data_rds_files/Ye_2016_meta.rds",
  #samples = metadata,
  dds_file = "data_rds_files/CD_ROSlow_dds.rds",
  dds_res_file = "data_rds_files/Ye_16_dds.res.rds",
  vsd_file = "data_rds_files/Ye_16_vsd.rds",
  vst_file = "data_rds_files/Ye_16_vst.rds",
  vsd.pca_file = "data_rds_files/Ye_16_vsd.pca.rds",
  vst.goi_file = "data_rds_files/Ye_16_vst.goi.rds",
  num_PC = nrow(metadata),
  qc = "data/multiqc_data.json",
  batch = NULL,
  var_1 = metadata$Source,
  var_2 = NULL)


