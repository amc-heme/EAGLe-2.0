# Run from the project root: Rscript tests/check_gsea_heatmaps.R
suppressPackageStartupMessages(source("R/00_libraries.R"))
source("R/01_dataset_config.R")
for (file in c("module_precomputed", "module_expression_heatmap", "module_singlecolor_palette",
               "module_DE", "module_GSEA")) source(paste0("R/", file, ".R"))
local({
  # Reproduce Positron's restriction even on a host that normally allows forks.
  trace("mcfork", where = asNamespace("parallel"), tracer = quote(stop("Fork forbidden by regression test")), print = FALSE)
  on.exit(untrace("mcfork", where = asNamespace("parallel")), add = TRUE)
  object <- load_precomputed_analysis(getwd(), load_analysis_config()$T0)
  results <- precomputed_display_results(precomputed_results(object))
  expression <- reactiveVal(precomputed_vst_table(object))
  symbols <- isolate(expression())$ext_gene_ensembl[
    match(results$ensembl_gene_id[results$padj < 0.05 & abs(results$log2FoldChange) >= 0.5],
          isolate(expression())$ensembl_gene_id)]
  symbols <- unique(symbols[!is.na(symbols)])
  stopifnot("UTS2" %in% symbols, length(symbols) > 10)
  # Real T0 expression and results; small gene sets exercise actual enrichment quickly.
  original_pathways <- pathways.hallmark
  on.exit(assign("pathways.hallmark", original_pathways, envir = .GlobalEnv), add = TRUE)
  assign("pathways.hallmark", list(single = "UTS2", multiple = head(symbols, 10), empty = "NOT_A_GENE"),
         envir = .GlobalEnv)
  choice <- list(user_dataset = reactive("T0"), user_model = reactive("KO_Sensitivity"),
                 user_PW = reactive("Sensitive_vs_Resistant"), close_tab = reactive(1))
  object_r <- reactive(object)
  testServer(GSEA_Server, args = list(dataset_choice = choice,
    DE_res = list(res_tidy = reactive(precomputed_results(object, tidy = TRUE)), dds_res = reactive(results)),
    reset_trigger = reactive(0), vst = expression, dataset_dds = object_r), {
    session$setInputs(filechoice = "hallmark", fgseaTable = TRUE, heatmap = FALSE,
                      `colorHM1-colorpalette` = "#273F52", `colorHM2-colorpalette` = "#D53031")
    # The real fgsea computation must succeed with a fork-forbidding host/backend.
    stopifnot(nrow(fgseaRes()) == 2, length(output$fgseaResults) > 0)
    session$setInputs(heatmap = TRUE, pathwaylistht = "single")
    stopifnot(nrow(gsea_heatmap_data()) == 1, !null_check(), length(output$htgsea) > 0)
    session$setInputs(pathwaylistht = "multiple")
    stopifnot(nrow(gsea_heatmap_data()) > 1, length(output$htgsea) > 0,
              identical(rownames(gsea_heatmap_layout()$ordered_annotation), colnames(gsea_heatmap_layout()$ordered)))
    session$setInputs(pathwaylistht = "empty")
    stopifnot(null_check(), grepl("No significant DEGs", output$htwarnGSEA$html))
    invisible(output$htgsea)
    constant <- isolate(expression())
    constant[constant$ext_gene_ensembl == "UTS2", -(1:2)] <- 2
    expression(constant)
    session$setInputs(pathwaylistht = "single")
    stopifnot(null_check(), grepl("variable expression", output$htwarnGSEA$html))
    invisible(output$htgsea)
  })
  cat("PASS: real enrichment without forks; GSEA single/multiple/empty/constant heatmaps\n")
  # Exercise the DE heatmap with exactly one qualifying gene, then a constant row.
  object_results <- S4Vectors::metadata(object)$results
  object_results$padj[object_results$Gene != "UTS2"] <- 1
  S4Vectors::metadata(object)$results <- object_results
  object_r <- reactive(object)
  expression(precomputed_vst_table(object))
  testServer(DE_Server, args = list(data_species = reactive("human"), dataset_dds = object_r,
    dataset_choice = choice, reset_trigger = reactive(0), vst = expression, vst_hm = NULL), {
    session$setInputs(DESeqHeat = TRUE, `color-colorpalette` = "#273F52", `color2-colorpalette` = "#D53031")
    stopifnot(nrow(de_heatmap_data()) == 1, length(output$ht) > 0)
    constant <- isolate(expression())
    constant[constant$ext_gene_ensembl == "UTS2", -(1:2)] <- 2
    expression(constant)
    stopifnot(nrow(de_heatmap_data()) == 0, grepl("variable expression", output$htwarn$html))
    invisible(output$ht)
  })
  cat("PASS: DE single-gene and no-usable-expression heatmaps\n")
})
