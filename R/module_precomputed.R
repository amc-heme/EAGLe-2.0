# Adapters for analyses supplied as results, metadata and a VST matrix.
load_precomputed_analysis <- function(raw_data, config) {
  values <- readRDS(paste0(raw_data, config$vst_path))
  meta <- as.data.frame(readRDS(paste0(raw_data, config$metadata_path)))
  res <- readRDS(paste0(raw_data, config$path))
  if (is.matrix(values)) {
    symbols <- res$Gene[match(rownames(values), res$ensembl_gene_id)]
  } else {
    symbols <- values$ext_gene
    ids <- values$ensembl_gene_id
    values <- as.matrix(values[, setdiff(names(values), c("ensembl_gene_id", "ext_gene")), drop = FALSE])
    rownames(values) <- ids
  }
  # Match the existing heatmap convention: retain the first row for repeated gene IDs.
  unique_gene <- !duplicated(rownames(values))
  symbols <- symbols[unique_gene]
  values <- values[unique_gene, , drop = FALSE]
  sample_ids <- if (is.null(config$metadata_id)) meta[[config$ID]] else {
    if (config$metadata_id == ".rownames") rownames(meta) else meta[[config$metadata_id]]
  }
  # Older tables prefix digit-leading UUIDs with X, while metadata retains the UUID.
  normalize_id <- function(x) sub("^X(?=[0-9])", "", x, perl = TRUE)
  sample_match <- match(normalize_id(colnames(values)), normalize_id(sample_ids))
  keep <- !is.na(sample_match)
  stopifnot(!anyDuplicated(sample_ids), sum(keep) == nrow(meta),
            !anyDuplicated(sample_match[keep]), !anyDuplicated(rownames(values)),
            any(res$ensembl_gene_id %in% rownames(values)))
  values <- values[, keep, drop = FALSE]
  meta <- meta[sample_match[keep], , drop = FALSE]
  meta[[config$ID]] <- colnames(values)
  rownames(meta) <- colnames(values)
  meta[[config$PCA_var]] <- factor(meta[[config$PCA_var]])
  SummarizedExperiment::SummarizedExperiment(
    assays = list(vst = values), colData = S4Vectors::DataFrame(meta),
    rowData = S4Vectors::DataFrame(symbol = symbols),
    metadata = list(precomputed = TRUE, results = res, config = config)
  )
}

is_precomputed_analysis <- function(x) {
  isTRUE(S4Vectors::metadata(x)$precomputed)
}

analysis_vsd <- function(x) {
  if (is_precomputed_analysis(x)) return(x)
  DESeq2::vst(x, blind = FALSE)
}

precomputed_vst_table <- function(x) {
  values <- SummarizedExperiment::assay(x)
  symbols <- SummarizedExperiment::rowData(x)$symbol
  symbols[is.na(symbols) | symbols == ""] <- rownames(values)[is.na(symbols) | symbols == ""]
  data.frame(ext_gene_ensembl = symbols, ensembl_gene_id = rownames(values),
             values, check.names = FALSE, row.names = NULL)
}

precomputed_results <- function(x, tidy = FALSE, signed_lrt = FALSE, comparison = NULL) {
  config <- S4Vectors::metadata(x)$config
  res <- S4Vectors::metadata(x)$results
  if (!is.null(config$comparisons)) {
    suffix <- config$comparisons[[comparison]]
    if (is.null(suffix)) stop("No saved results for comparison: ", comparison)
    res$padj <- res[[paste0("padj_", suffix)]]
    res$log2FoldChange <- res[[paste0("log2FC_", suffix)]]
    res$DiffExp <- ifelse(res$padj < 0.05 & res$log2FoldChange >= 0.5, "up",
                         ifelse(res$padj < 0.05 & res$log2FoldChange <= -0.5, "down", "no"))
  }
  if (tidy) {
    # Some saved pairwise summaries have fold changes and padj, but no test statistic.
    if (identical(config$gsea_rank, "log2FoldChange")) res$stat <- res$log2FoldChange
    if (signed_lrt || identical(config$gsea_rank, "signed_sqrt_lrt")) {
      res$stat <- sign(res$log2FoldChange) * sqrt(pmax(res$stat, 0))
    }
    for (column in c("baseMean", "lfcSE", "pvalue")) {
      if (!column %in% names(res)) res[[column]] <- NA_real_
    }
    res <- res[, c("ensembl_gene_id", "baseMean", "log2FoldChange", "lfcSE",
                   "stat", "pvalue", "padj")]
    names(res)[1] <- "row"
  }
  res
}

precomputed_display_results <- function(res) {
  if (!"baseMean" %in% names(res)) res$baseMean <- NA_real_
  if (!"DiffExp" %in% names(res)) {
    res$DiffExp <- ifelse(res$padj < 0.05 & res$log2FoldChange >= 0.5, "up",
                         ifelse(res$padj < 0.05 & res$log2FoldChange <= -0.5, "down", "no"))
  }
  res <- res[, c("ensembl_gene_id", "Gene", "baseMean", "log2FoldChange", "padj", "DiffExp")]
  res[!is.na(res$padj) & !is.na(res$log2FoldChange), , drop = FALSE]
}
