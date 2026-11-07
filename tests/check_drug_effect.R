# Run from the project root: Rscript tests/check_drug_effect.R
suppressPackageStartupMessages(source("R/00_libraries.R"))
source("R/01_dataset_config.R")
for (file in list.files("R", "^module.*[.]R$", full.names = TRUE)) source(file)
config <- load_analysis_config()
entry <- config$drug_effect
stopifnot(entry$type == "bulk", entry$precomputed, is.null(entry$qc_path))
object <- load_precomputed_analysis(getwd(), entry)
original <- readRDS(paste0(getwd(), entry$path))
values <- readRDS(paste0(getwd(), entry$vst_path))
stopifnot(identical(assay(object), values), ncol(object) == 44,
          identical(colnames(object), colData(object)$Sample),
          all(table(colData(object)$Treatment) == 22),
          identical(precomputed_results(object), original),
          identical(precomputed_results(object, tidy = TRUE)$stat, original$stat))
for (file in c("genes_by_dataset.yaml", "DEres_by_dataset.yaml")) {
  stopifnot(setequal(read_yaml(file)$drug_effect, original$ensembl_gene_id))
}
search <- read.table("new_raw_data/input_lookup_table.txt", sep = "\t", header = TRUE)
stopifnot(all(rownames(object) %in% search$ensembl_gene_id))
choice <- list(user_dataset = reactive("drug_effect"), user_model = reactive("Treatment"),
               user_PW = reactive("KO_vs_T0"), close_tab = reactive(1))
object_r <- reactive(object)
vst_r <- reactive(precomputed_vst_table(object))
testServer(data_Server, {
  session$setInputs(datainput = "drug_effect", DEmodel = "Treatment", pwc = "KO_vs_T0")
  stopifnot(DEModelChoices("drug_effect") == "Treatment", PWChoices("drug_effect", "Treatment") == "KO_vs_T0",
            grepl("44 samples", output$Data_text$html))
})
testServer(vst_Server, args = list(dataset_dds = object_r, dataset_choice = choice), {
  stopifnot(nrow(vst.global()) == nrow(values), ncol(vst.global()) == 46)
})
testServer(QC_Server, args = list(dataset_dds = object_r, dataset_choice = choice,
  qc_table = reactive(NULL), reset_trigger = reactive(0)), {
  session$setInputs(PCAplots = TRUE, PCAscreeplots = TRUE, multiqc = TRUE)
  stopifnot(nrow(vsd.pca()) == 44, length(pca_color()) == 44, length(pca_shape()) == 44,
            color_legend() == "Treatment", shape_legend() == "Batch",
            nrow(scree_variance_df()) == 44, grepl("unavailable", output$QCdt2$html))
})
testServer(DE_Server, args = list(data_species = reactive("human"), dataset_dds = object_r,
  dataset_choice = choice, reset_trigger = reactive(0), vst = vst_r, vst_hm = NULL), {
  session$setInputs(DESeqtable = TRUE, DESeqvolcano = TRUE, DESeqMA = TRUE, DESeqHeat = TRUE,
                    `color-colorpalette` = "#2166AC", `color2-colorpalette` = "#D53031")
  stopifnot(identical(dds_result(), original), identical(DDS4GSEA()$stat, original$stat),
            length(output$results) > 0, length(output$volplot) > 0, length(output$MAplot) > 0,
            nrow(de_heatmap_data()) > 0, all(is.finite(de_heatmap_data())), length(output$ht) > 0)
})
testServer(goi_Server, args = list(dataset_choice = choice, dataset_dds = object_r, vst = vst_r), {
  session$setInputs(VSTCDgenechoice = "SLC41A2", plotwidthslider = 800, plotheightslider = 600)
  stopifnot(nrow(vst.gene()) == 44, nlevels(vst.gene()$condition) == 2)
})
testServer(GSEA_Server, args = list(dataset_choice = choice,
  DE_res = list(res_tidy = reactive(precomputed_results(object, tidy = TRUE)),
                dds_res = reactive(precomputed_display_results(original))),
  reset_trigger = reactive(0), vst = vst_r, dataset_dds = object_r), {
  session$setInputs(filechoice = "hallmark", fgseaTable = TRUE, heatmap = FALSE,
                    `colorHM1-colorpalette` = "#2166AC", `colorHM2-colorpalette` = "#D53031")
  stopifnot(all(is.finite(ranks())), any(ranks() < 0), any(ranks() > 0),
            nrow(fgseaRes()) > 0, length(output$fgseaResults) > 0)
  sig <- original$Gene[!is.na(original$padj) & original$padj < 0.05 & abs(original$log2FoldChange) >= 0.5]
  pathway <- names(pathways.hallmark)[which(vapply(pathways.hallmark, function(x) sum(x %in% sig) > 1, logical(1)))[1]]
  session$setInputs(heatmap = TRUE, pathwaylistht = pathway)
  stopifnot(nrow(gsea_heatmap_data()) > 0, length(output$htgsea) > 0)
})
testServer(IRF2BP2_Server, args = list(dataset_choice = choice, dataset_dds = object_r,
  dataset_config = config, raw_data = getwd()), {
  session$setInputs(group_by = "Treatment", calculate = 1)
  stopifnot(nrow(score_data()) == 44, all(is.finite(score_data()$IRF2BP2_score)),
            length(unique(score_data()$Group)) == 2, length(output$score_plot) > 0,
            file.info(output$download_scores)$size > 0)
})
cat("PASS: drug-effect data, selection, QC, DE plots, gene expression, real GSEA/heatmaps, and IRF2BP2 scoring\n")
# Knit the full report with the new dataset selected, including its summary and expression plot.
report <- new.env(parent = globalenv())
report$params <- list(gene = original$ensembl_gene_id[1], gene_symbol = original$Gene[1], raw_data = getwd(),
                      datasets = read_yaml("datasets.yaml"), gene_present = c(drug_effect = TRUE),
                      gene_DE = c(drug_effect = TRUE), selected_datasets = "drug_effect")
knitr::opts_knit$set(root.dir = getwd())
figures <- file.path(tempdir(), "drug-effect-figures/")
knitr::opts_chunk$set(error = FALSE, fig.path = figures)
output <- file.path(tempdir(), "drug-effect-report.md")
knitr::knit("EAGLe_Analysis_v7.0.Rmd", output = output, envir = report, quiet = TRUE)
text <- readLines(output)
stopifnot(any(grepl("2.7. Drug-induced changes: KO versus T0 across all samples", text)),
          !any(grepl("not found in any of the Bulk", text)),
          length(list.files(figures, "png$")) == 1)
cat("PASS: bulk report summary, drug-effect section and expression figure\n")
