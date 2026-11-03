# Shared preparation for DE and pathway heatmaps, including downloaded images.
prepare_expression_heatmap <- function(vst, results, min_lfc, genes = NULL, max_genes = Inf) {
  selected <- results[!is.na(results$padj) & is.finite(results$log2FoldChange) &
                        results$padj < 0.05 & abs(results$log2FoldChange) >= min_lfc, , drop = FALSE]
  selected <- selected[order(abs(selected$log2FoldChange), decreasing = TRUE), , drop = FALSE]
  selected <- selected[!duplicated(selected$ensembl_gene_id), , drop = FALSE]
  rows <- vst[match(selected$ensembl_gene_id, vst$ensembl_gene_id, nomatch = 0), , drop = FALSE]
  symbols <- as.character(rows$ext_gene_ensembl)
  missing <- is.na(symbols) | symbols == ""
  symbols[missing] <- rows$ensembl_gene_id[missing]
  keep <- !duplicated(symbols)
  if (!is.null(genes)) keep <- keep & symbols %in% genes
  rows <- rows[keep, , drop = FALSE]
  symbols <- symbols[keep]
  values <- as.matrix(rows[, setdiff(names(rows), c("ensembl_gene_id", "ext_gene_ensembl")), drop = FALSE])
  rownames(values) <- symbols
  if (nrow(values) == 0 || ncol(values) < 2) return(values[integer(0), , drop = FALSE])
  # Constant rows produce NaNs when z-scored and cannot be clustered.
  finite <- rowSums(!is.finite(values)) == 0
  values <- values[finite, , drop = FALSE]
  deviations <- apply(values, 1, stats::sd)
  values <- values[is.finite(deviations) & deviations > 0, , drop = FALSE]
  if (nrow(values) == 0) return(values)
  values <- values[seq_len(min(nrow(values), max_genes)), , drop = FALSE]
  t(scale(t(values)))
}

expression_heatmap_layout <- function(values, dds, condition) {
  meta <- as.data.frame(SummarizedExperiment::colData(dds))
  sample_ids <- rownames(meta)
  indices <- match(colnames(values), sample_ids)
  # VST tables made from DESeq2 objects use janitor-cleaned sample names.
  if (anyNA(indices)) indices <- match(colnames(values), janitor::make_clean_names(sample_ids))
  stopifnot(!anyNA(indices), !anyDuplicated(indices))
  annotation <- data.frame(Condition = factor(meta[[condition]][indices]))
  rownames(annotation) <- colnames(values)
  names(annotation) <- condition
  row_tree <- if (nrow(values) > 1) as.dendrogram(stats::hclust(stats::dist(values))) else FALSE
  column_tree <- if (ncol(values) > 1) as.dendrogram(stats::hclust(stats::dist(t(values)))) else FALSE
  rows <- if (isFALSE(row_tree)) seq_len(nrow(values)) else order.dendrogram(row_tree)
  columns <- if (isFALSE(column_tree)) seq_len(ncol(values)) else order.dendrogram(column_tree)
  list(values = values, annotation = annotation, row_tree = row_tree, column_tree = column_tree,
       ordered = values[rows, columns, drop = FALSE],
       ordered_annotation = annotation[columns, , drop = FALSE])
}

# Use the same symmetric scale for the browser and saved images.
expression_heatmap_colors <- function(values, colors) {
  limit <- max(2, abs(values), na.rm = TRUE)
  map <- circlize::colorRamp2(c(-limit, 0, limit), c(colors[1], "#FFFFFF", colors[3]))
  palette <- substr(map(seq(-limit, limit, length.out = 257)), 1, 7)
  list(limit = limit, map = map, palette = palette,
       plotly = lapply(seq_along(palette), function(i) list((i - 1) / (length(palette) - 1), palette[i])))
}

interactive_expression_heatmap <- function(layout, colors) {
  values <- layout$ordered
  scale <- expression_heatmap_colors(values, colors)
  if (nrow(values) == 1) {
    # heatmaply's side-color layout does not support a single gene reliably.
    labels <- matrix(paste("Gene:", rownames(values), "<br>Sample:", colnames(values),
                           "<br>Condition:", layout$ordered_annotation[[1]]), nrow = 1)
    plot <- plotly::plot_ly(x = colnames(values), y = rownames(values), z = unname(values),
                           text = labels, hoverinfo = "text+z", type = "heatmap", colorscale = scale$plotly,
                           zmin = -scale$limit, zmax = scale$limit, zauto = FALSE,
                           height = 600, width = 600)
    return(plotly::layout(plot, xaxis = list(showticklabels = FALSE), yaxis = list(showticklabels = FALSE)))
  }
  heatmaply::heatmaply(values, plot_method = "plotly", row_text_angle = 45, height = 600, width = 600,
                      colors = scale$palette, limits = c(-scale$limit, scale$limit),
                      dendrogram = "none", Rowv = FALSE, Colv = FALSE,
                      col_side_colors = layout$ordered_annotation,
                      showticklabels = c(FALSE, FALSE))
}

save_expression_heatmap <- function(file, type, layout, colors, title) {
  scale <- expression_heatmap_colors(layout$values, colors)
  condition <- names(layout$annotation)[1]
  groups <- levels(layout$annotation[[1]])
  annotation_colors <- setNames(list(setNames(grDevices::hcl.colors(length(groups), "Dark 3"), groups)), condition)
  annotation <- ComplexHeatmap::HeatmapAnnotation(
    df = layout$annotation, col = annotation_colors, show_annotation_name = FALSE,
    annotation_legend_param = list(title_gp = grid::gpar(fontsize = 12), labels_gp = grid::gpar(fontsize = 12))
  )
  plot <- ComplexHeatmap::Heatmap(
    layout$values, name = "z scaled expression", col = scale$map,
    cluster_rows = layout$row_tree, cluster_columns = layout$column_tree,
    row_dend_reorder = FALSE, column_dend_reorder = FALSE,
    row_names_gp = grid::gpar(fontsize = 12), column_names_gp = grid::gpar(fontsize = 12),
    column_title = NULL, row_title = title, top_annotation = annotation,
    heatmap_legend_param = list(title_gp = grid::gpar(fontsize = 12), labels_gp = grid::gpar(fontsize = 12))
  )
  if (type == "png") grDevices::png(file, width = 780, height = 780, units = "px", pointsize = 12)
  else svglite::svglite(file, width = 8, height = 7, pointsize = 12)
  on.exit(grDevices::dev.off(), add = TRUE)
  ComplexHeatmap::draw(plot)
  invisible(file)
}
