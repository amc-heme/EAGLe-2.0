# Run from the project root: Rscript tests/check_heatmaps.R
source("R/module_expression_heatmap.R")
# The pure preparation helper does not depend on Shiny or fitted DESeq2 objects.
vst <- data.frame(ensembl_gene_id = c("g1", "g2", "g3", "g4", "g5"),
                  ext_gene_ensembl = c("A", "B", "constant", "invalid", "A"),
                  s3 = c(3, 8, 2, NA, 99), s1 = c(1, 4, 2, 2, 99), s2 = c(2, 6, 2, 3, 99))
res <- data.frame(ensembl_gene_id = c("g4", "missing", "g2", "g3", "g1", "g5"),
                  padj = 0.01, log2FoldChange = c(10, 9, 8, 7, 6, 5))
values <- prepare_expression_heatmap(vst, res, min_lfc = 2)
stopifnot(identical(rownames(values), c("B", "A")), all(is.finite(values)),
          isTRUE(all.equal(unname(values["A", ]), c(1, -1, 0))),
          isTRUE(all.equal(unname(values["B", ]), c(1, -1, 0))))
# Empty selections, all-constant data, and one-sample data cannot be z-scored.
stopifnot(nrow(prepare_expression_heatmap(vst, res, 2, genes = "absent")) == 0,
          nrow(prepare_expression_heatmap(vst, res, 2, genes = "constant")) == 0,
          nrow(prepare_expression_heatmap(vst[, c(1, 2, 3)], res, 2)) == 0)
one <- prepare_expression_heatmap(vst, res, 2, genes = "A")
stopifnot(identical(dim(one), c(1L, 3L)), all(is.finite(one)))
# Metadata is matched by sample ID, even when its row order is different.
meta <- data.frame(condition = c("first", "second", "third"), row.names = c("s1", "s2", "s3"))
object <- SummarizedExperiment::SummarizedExperiment(
  assays = list(vst = matrix(0, 1, 3, dimnames = list("dummy", rownames(meta)))),
  colData = S4Vectors::DataFrame(meta))
layout <- expression_heatmap_layout(values, object, "condition")
stopifnot(identical(as.character(layout$annotation$condition), c("third", "first", "second")),
          identical(rownames(layout$ordered_annotation), colnames(layout$ordered)))
# A single gene must render and export without trying to cluster its one row.
colors <- c("#2166AC", "white", "#D53031")
for (matrix in list(one, values)) {
  layout <- expression_heatmap_layout(matrix, object, "condition")
  widget <- interactive_expression_heatmap(layout, colors)
  stopifnot(inherits(widget, "plotly"), length(plotly::plotly_build(widget)$x$data) > 0)
  for (type in c("png", "svg")) {
    file <- tempfile(fileext = paste0(".", type))
    device <- grDevices::dev.cur()
    save_expression_heatmap(file, type, layout, colors, "Test")
    stopifnot(file.info(file)$size > 0, identical(grDevices::dev.cur(), device))
    unlink(file)
  }
}
cat("PASS: gene/sample alignment, empty/constant/non-finite data, single-gene rendering, and PNG/SVG exports\n")
# An asymmetric data range must not move white away from zero.
asymmetric <- matrix(c(-1, 0, 3, -0.2, 0.1, 0.2), nrow = 2, byrow = TRUE,
                     dimnames = list(c("A", "B"), c("s3", "s1", "s2")))
scale <- expression_heatmap_colors(asymmetric, colors)
stopifnot(scale$limit == 3, substr(scale$map(0), 1, 7) == "#FFFFFF", scale$palette[129] == "#FFFFFF")
negative <- grDevices::col2rgb(substr(scale$map(-0.5), 1, 7))
positive <- grDevices::col2rgb(substr(scale$map(0.5), 1, 7))
stopifnot(negative[3, 1] > negative[1, 1], positive[1, 1] > positive[3, 1])
for (matrix in list(asymmetric, asymmetric[1, , drop = FALSE])) {
  layout <- expression_heatmap_layout(matrix, object, "condition")
  traces <- plotly::plotly_build(interactive_expression_heatmap(layout, colors))$x$data
  traces <- Filter(function(x) identical(x$type, "heatmap") && identical(dim(x$z), dim(matrix)), traces)
  stopifnot(length(traces) == 1)
  trace <- traces[[1]]
  stops <- if (is.matrix(trace$colorscale)) trace$colorscale else do.call(rbind, trace$colorscale)
  midpoint <- stops[abs(as.numeric(stops[, 1]) - 0.5) < 1e-12, 2]
  stopifnot(trace$zmin == -3, trace$zmax == 3, length(midpoint) == 1,
            tolower(gsub(" ", "", midpoint)) %in% c("#ffffff", "rgba(255,255,255,1)"))
}
cat("PASS: single/multiple-gene browser scales and export maps use white at zero, blue below, red above\n")
