
# globalDataServer <- function(id) {
#   moduleServer(id, function(input, output, session) {
# 
#     stash = reactiveValues()
#     stash$t2g = readRDS("data_rds_files/t2g_hs.rds") 
#     #stash$metadata = readRDS("data_rds_files/CD_ROSlow_metadata.rds")
#     stash$dds = readRDS("data_rds_files/CD_ROSlow_dds.rds")
#     #stash$dds_res = readRDS("data_rds_files/CD_ROSlow_dds.res.rds")
#     stash$vsd = readRDS("data_rds_files/CD_ROSlow_vsd.rds")
#     #stash$vst = readRDS("data_rds_files/CD_vst_batchcorrected.rds")
#     #stash$vsd.pca = readRDS("data_rds_files/CD_ROSlow_vsd.pca.rds")
#     #stash$vst.goi = readRDS("data_rds_files/CD_ROSlow_vst.goi.rds")
#     stash$qc = "data/multiqc_data.json"
#     
#     # return(list(get_t2g = reactive(stash$t2g), get_metadata = reactive(stash$metadata), get_dds = reactive(stash$dds),
#     #             get_dds.res = reactive(stash$dds_res), get_vsd = reactive(stash$vsd), get_vst = reactive(stash$vst),
#     #             get_vsd.pca = reactive(stash$vsd.pca), get_vst.goi = reactive(stash$vst.goi), get_qc = reactive(stash$qc)))
#   })
# }
Datasets:
  Cancer_Discovery:
    path: "~/Documents/GitHub/EAGLe-2.0/data_rds_files/CD_ROSlow_dds.rds"
    species: human
  Ye_16:
    path: "~/Documents/GitHub/EAGLe-2.0/data_rds_files/Ye_16_dds.rds"
    species: mouse
  Ye_20:
    path: "~/Documents/GitHub/EAGLe-2.0/data_rds_files/Ye_20_dds.rds"
    species: mouse
  Venaza:
    path: "~/Documents/GitHub/EAGLe-2.0/data_rds_files/venaza_dds.rds"
    species: human
  Lagadinou: 
    path: "~/Documents/GitHub/EAGLe-2.0/data_rds_files/EL_AML_dds.rds"
    species: human
  Lee:
    path: "~/Documents/GitHub/EAGLe-2.0/data_rds_files/Lee_dds.rds"
    species: human
  BEAT:
    path: "~/Documents/GitHub/EAGLe-2.0/data_rds_files/BEAT_AML_dds.rds"
    species: human
Transcript_to_gene:
  mouse: "~/Documents/GitHub/EAGLe-2.0/data_rds_files/t2g_mm.rds"
  human: "~/Documents/GitHub/EAGLe-2.0/data_rds_files/t2g_hs.rds"
  
