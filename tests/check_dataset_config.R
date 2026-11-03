# Run from the project root: Rscript tests/check_dataset_config.R
source("R/01_dataset_config.R")
source("R/module_precomputed.R")
all_config <- yaml::read_yaml("datasets.yaml")
config <- load_analysis_config()
stopifnot(identical(sort(list.files("R", "[.]R$"))[[1]], "00_libraries.R"),
          identical(config$Cancer_Discovery$path, all_config[["Jordan_M0-M5"]]$path),
          identical(config$BEAT_quantile$path, all_config$BEAT_AML$path),
          !file.exists("data.yaml"))
for (key in names(all_config)) {
  for (field in c("path", "vst_path", "dds_path", "metadata_path", "qc_path", "signature_score_path")) {
    path <- all_config[[key]][[field]]
    if (!is.null(path)) stopifnot(file.exists(paste0(getwd(), path)))
  }
}
for (key in names(config)) {
  entry <- config[[key]]
  if (!isTRUE(entry$precomputed)) {
    stopifnot(!is.null(entry$qc_path))
    next
  }
  stopifnot(is.null(entry$qc_path))
  object <- load_precomputed_analysis(getwd(), entry)
  values <- SummarizedExperiment::assay(object)
  meta <- as.data.frame(SummarizedExperiment::colData(object))
  stopifnot(ncol(values) == nrow(meta), identical(colnames(values), meta[[entry$ID]]),
            !anyNA(meta[[entry$PCA_var]]))
  original <- readRDS(paste0(getwd(), entry$path))
  comparisons <- if (is.null(entry$comparisons)) NA_character_ else names(entry$comparisons)
  for (comparison in comparisons) {
    comparison <- if (is.na(comparison)) NULL else comparison
    results <- precomputed_results(object, comparison = comparison)
    ranked <- precomputed_results(object, tidy = TRUE, comparison = comparison)
    suffix <- if (is.null(comparison)) NULL else entry$comparisons[[comparison]]
    padj_col <- if (is.null(suffix)) "padj" else paste0("padj_", suffix)
    lfc_col <- if (is.null(suffix)) "log2FoldChange" else paste0("log2FC_", suffix)
    stopifnot(identical(results$padj, original[[padj_col]]),
              identical(results$log2FoldChange, original[[lfc_col]]))
    if (identical(entry$gsea_rank, "log2FoldChange")) {
      stopifnot(identical(ranked$stat, original[[lfc_col]]),
                all(is.na(precomputed_display_results(results)$baseMean)))
    }
    if (identical(entry$gsea_rank, "signed_sqrt_lrt")) {
      stopifnot(identical(ranked$stat, sign(original$log2FoldChange) * sqrt(pmax(original$stat, 0))))
    }
  }
  # Gene annotations follow VST IDs even when the result table has another order.
  vtable <- precomputed_vst_table(object)
  stopifnot(identical(vtable$ensembl_gene_id, rownames(values)),
            isTRUE(all.equal(as.matrix(vtable[, -(1:2)]), values, check.attributes = FALSE)))
  cat("PASS:", key, "metadata, VST, saved contrasts and ranking\n")
}
stopifnot(!any(c("M6_vs_M5", "M7_vs_M5") %in% names(config$TCGA_FAB$comparisons)))
cat("PASS: shared configuration paths and QC availability\n")
