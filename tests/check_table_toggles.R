# Run from the project root: Rscript tests/check_table_toggles.R
suppressPackageStartupMessages(source("R/00_libraries.R"))
source("R/01_dataset_config.R")
for (file in c("module_precomputed", "module_singlecolor_palette", "module_color_palettes",
               "module_DE", "module_GSEA", "module_QC", "module_IRF2BP2")) source(paste0("R/", file, ".R"))
events <- list()
record_suspended <- function(value, output) {
  result <- tryCatch(force(value), shiny.silent.error = function(e) e)
  stopifnot(inherits(result, "shiny.silent.error"))
  events[[length(events) + 1]] <<- list(output = output, suspended = TRUE)
}
record_table <- function(value, output) {
  result <- jsonlite::fromJSON(value, simplifyVector = FALSE)
  stopifnot(!is.null(result$x), !is.null(result$x$container))
  events[[length(events) + 1]] <<- list(output = output, suspended = FALSE, payload = result)
}
config <- load_analysis_config()
objects <- lapply(config[c("T0", "KO_response", "BEAT_quantile")], function(x) load_precomputed_analysis(getwd(), x))
key <- reactiveVal("T0")
run <- reactiveVal(1)
object <- reactive(objects[[key()]])
choice <- list(user_dataset = reactive(key()), user_model = reactive(config[[key()]]$PCA_var),
               user_PW = reactive("Sensitive_vs_Resistant"), close_tab = reactive(run()))
testServer(DE_Server, args = list(data_species = reactive("human"), dataset_dds = object, dataset_choice = choice,
  reset_trigger = reactive(0), vst = reactive(precomputed_vst_table(object())), vst_hm = NULL), {
  session$setInputs(DESeqtable = FALSE)
  record_suspended(output$results, "DE table initially hidden")
  session$setInputs(DESeqtable = TRUE)
  record_table(output$results, "DE table visible")
  session$setInputs(DESeqtable = FALSE)
  record_suspended(output$results, "DE table switched off")
  key("KO_response")
  run(2)
  session$flushReact()
  record_suspended(output$results, "DE table hidden after dataset change")
  session$setInputs(DESeqtable = TRUE)
  record_table(output$results, "DE table visible after dataset change")
})
key("T0")
pathways.hallmark <- list(example = "UTS2")
testServer(GSEA_Server, args = list(dataset_choice = choice,
  DE_res = list(res_tidy = reactive(precomputed_results(object(), tidy = TRUE)),
                dds_res = reactive(precomputed_display_results(precomputed_results(object())))),
  reset_trigger = reactive(0), vst = reactive(precomputed_vst_table(object())), dataset_dds = object), {
  session$setInputs(filechoice = "hallmark", fgseaTable = FALSE, heatmap = FALSE)
  record_suspended(output$fgseaResults, "GSEA table initially hidden")
  session$setInputs(fgseaTable = TRUE)
  record_table(output$fgseaResults, "GSEA table visible")
  session$setInputs(fgseaTable = FALSE)
  record_suspended(output$fgseaResults, "GSEA table switched off")
})
qc <- reactiveVal(data.frame(Sample = "test", aligned = 0.95))
testServer(QC_Server, args = list(dataset_dds = object, dataset_choice = choice,
  qc_table = qc, reset_trigger = reactive(0)), {
  session$setInputs(multiqc = FALSE)
  record_suspended(output$QCdt, "QC table initially hidden")
  session$setInputs(multiqc = TRUE)
  record_table(output$QCdt, "QC table visible")
  qc(NULL)
  session$flushReact()
  record_suspended(output$QCdt, "QC data unavailable")
  stopifnot(grepl("unavailable", output$QCdt2$html))
  qc(data.frame(Sample = "test2", aligned = 0.9))
  session$flushReact()
  record_table(output$QCdt, "QC table available again")
  session$setInputs(multiqc = FALSE)
  record_suspended(output$QCdt, "QC table switched off")
})
key("BEAT_quantile")
testServer(IRF2BP2_Server, args = list(dataset_choice = choice, dataset_dds = object,
  dataset_config = config, raw_data = getwd()), {
  session$setInputs(group_by = "quantile")
  record_suspended(output$score_table, "Signature table before calculation")
  session$setInputs(calculate = 1)
  record_table(output$score_table, "Signature table after calculation")
  key("T0")
  session$flushReact()
  record_suspended(output$score_table, "Signature table after dataset reset")
})
cat("PASS: all four DT outputs avoid null payloads when hidden, disabled, unavailable or reset\n")
# Check the actual installed client code, using the actual serialized server outputs above.
events_path <- tempfile(fileext = ".json")
jsonlite::write_json(events, events_path, auto_unbox = TRUE, null = "null")
node <- Sys.which("node")
stopifnot(nzchar(node))
status <- system2(node, c(shQuote("tests/check_dt_client.js"),
                         shQuote(system.file("htmlwidgets/datatables.js", package = "DT")), shQuote(events_path)))
unlink(events_path)
stopifnot(status == 0)
