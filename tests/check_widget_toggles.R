# Run from the project root: Rscript tests/check_widget_toggles.R
suppressPackageStartupMessages(source("R/00_libraries.R"))
source("R/01_dataset_config.R")
for (file in c("module_precomputed", "module_singlecolor_palette", "module_DE", "module_GSEA")) {
  source(paste0("R/", file, ".R"))
}
# A disabled ggiraph output must abort silently instead of serializing x:null.
expect_suspended <- function(value) {
  result <- tryCatch(force(value), shiny.silent.error = function(e) e)
  stopifnot(inherits(result, "shiny.silent.error"))
}
expect_widget <- function(value) {
  result <- jsonlite::fromJSON(value, simplifyVector = FALSE)
  stopifnot(!is.null(result$x), !is.null(result$x$html))
}
object <- load_precomputed_analysis(getwd(), load_analysis_config()$T0)
results <- precomputed_display_results(precomputed_results(object))
choice <- list(user_dataset = reactive("T0"), user_model = reactive("KO_Sensitivity"),
               user_PW = reactive("Sensitive_vs_Resistant"), close_tab = reactive(1))
testServer(DE_Server, args = list(data_species = reactive("human"), dataset_dds = reactive(object),
  dataset_choice = choice, reset_trigger = reactive(0), vst = reactive(precomputed_vst_table(object)), vst_hm = NULL), {
  session$setInputs(DESeqvolcano = FALSE, DESeqMA = FALSE,
                    `color-colorpalette` = "#273F52", `color2-colorpalette` = "#D53031")
  expect_suspended(output$volplot)
  expect_suspended(output$MAplot)
  session$setInputs(DESeqvolcano = TRUE, DESeqMA = TRUE)
  expect_widget(output$volplot)
  expect_widget(output$MAplot)
  session$setInputs(DESeqvolcano = FALSE, DESeqMA = FALSE)
  expect_suspended(output$volplot)
  expect_suspended(output$MAplot)
})
# A small real pathway keeps the enrichment calculation fast.
pathways.hallmark <- list(example = "UTS2")
testServer(GSEA_Server, args = list(dataset_choice = choice,
  DE_res = list(res_tidy = reactive(precomputed_results(object, tidy = TRUE)), dds_res = reactive(results)),
  reset_trigger = reactive(0), vst = reactive(precomputed_vst_table(object)), dataset_dds = reactive(object)), {
  session$setInputs(filechoice = "hallmark", moustache = FALSE, volcanoplot = FALSE,
                    pathwaylist = "example", heatmap = FALSE, `color7-colorpalette` = "#D53031")
  expect_suspended(output$GSEAMoustache)
  expect_suspended(output$GSEAvolcano)
  session$setInputs(moustache = TRUE, volcanoplot = TRUE)
  expect_widget(output$GSEAMoustache)
  expect_widget(output$GSEAvolcano)
  session$setInputs(moustache = FALSE, volcanoplot = FALSE)
  expect_suspended(output$GSEAMoustache)
  expect_suspended(output$GSEAvolcano)
})
cat("PASS: all ggiraph plot switches suspend safely and emit valid widgets when enabled\n")
