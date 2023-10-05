
#   num_PC = nrow(metadata),
#   var_1 = samples$condition,
#   var_2 = samples$batch)

globalDataServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    # list(
    #   vsd = reactive({readRDS("data_rds_files/CD_ROSlow_vsd.rds")}),
    #   vsd.pca = reactive({readRDS("data_rds_files/CD_ROSlow_vsd.pca.rds")})
    # )
    stash = reactiveValues()
    stash$t2g = readRDS("data_rds_files/t2g_hs.rds")
    #stash$metadata = readRDS("data_rds_files/CD_ROSlow_metadata.rds")
    stash$dds = readRDS("data_rds_files/CD_ROSlow_dds.rds")
    #stash$dds_res = readRDS("data_rds_files/CD_ROSlow_dds.res.rds")
    stash$vsd = readRDS("data_rds_files/CD_ROSlow_vsd.rds")
    #stash$vst = readRDS("data_rds_files/CD_vst_batchcorrected.rds")
    #stash$vsd.pca = readRDS("data_rds_files/CD_ROSlow_vsd.pca.rds")
    #stash$vst.goi = readRDS("data_rds_files/CD_ROSlow_vst.goi.rds")
    stash$qc = "data/multiqc_data.json"
    return(stash)
    # return(list(get_t2g = reactive(stash$t2g), get_metadata = reactive(stash$metadata), get_dds = reactive(stash$dds),
    #             get_dds.res = reactive(stash$dds_res), get_vsd = reactive(stash$vsd), get_vst = reactive(stash$vst),
    #             get_vsd.pca = reactive(stash$vsd.pca), get_vst.goi = reactive(stash$vst.goi), get_qc = reactive(stash$qc)))
  })
}
