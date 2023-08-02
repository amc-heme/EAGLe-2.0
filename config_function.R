create_config <- function(datainput) {
  reactiveValues(
    t2g_file = {
    if (datainput == "CancerDiscovery") "data_rds_files/t2g_hs.rds"
      else if (datainput == "Ye16") "data_rds_files/t2g_mm.rds"
    },
    metadata_file = {
      if (datainput == "CancerDiscovery") "data_rds_files/CD_ROSlow_metadata.rds"
     else if (datainput == "Ye16") "data_rds_files/Ye_2016_meta.rds"
    },
    dds_file = {
      if (datainput == "CancerDiscovery") "data_rds_files/CD_ROSlow_dds.rds"
      else if (datainput == "Ye16")  "data_rds_files/Ye_16_dds.rds"
    },
    dds_res_file = {
      if (datainput == "CancerDiscovery") "data_rds_files/CD_ROSlow_dds.res.rds"
      else if (datainput == "Ye16") "data_rds_files/Ye_16_dds.res.rds"
    },
      vsd_file = {
        if (datainput == "CancerDiscovery") "data_rds_files/CD_ROSlow_vsd.rds"
        else if (datainput == "Ye16") "data_rds_files/Ye_16_vsd.rds"
      },
      vst_file = {
        if (datainput == "CancerDiscovery") "data_rds_files/CD_vst_batchcorrected.rds"
          else if (datainput == "Ye16") "data_rds_files/Ye_16_vst.rds"
      },
      vsd.pca_file = {
        if (datainput == "CancerDiscovery") "data_rds_files/CD_ROSlow_vsd.pca.rds"
        else if (datainput == "Ye16") "data_rds_files/Ye_16_vsd.pca.rds"
      },
      vst.goi_file = {
        if (datainput == "CancerDiscovery") "data_rds_files/CD_ROSlow_vst.goi.rds"
        else if (datainput == "Ye16") "data_rds_files/Ye_16_vst.goi.rds"
      })
}
     
      