# The gene set and GSVA settings follow the existing BEAT-AML/TCGA scoring scripts.
read_irf2bp2_signature <- function(path) {
  signature <- read.delim(path, fill = TRUE, quote = "", comment.char = "", na.strings = c("", "NA"))
  symbols <- ifelse(!is.na(signature$gene_symbol) & nzchar(signature$gene_symbol),
                    signature$gene_symbol, signature$old_symbol)
  symbols <- toupper(trimws(symbols))
  unique(symbols[!is.na(symbols) & nzchar(symbols)])
}

irf2bp2_expression <- function(config, dds, raw_data) {
  expression <- readRDS(paste0(raw_data, config$vst_path))
  if (is.matrix(expression)) {
    values <- expression
    if (all(grepl("^ENSG", rownames(values)))) {
      annotation <- S4Vectors::metadata(dds)$results
      symbols <- annotation$Gene[match(rownames(values), annotation$ensembl_gene_id)]
    } else symbols <- rownames(values)
  } else {
    symbols <- expression$ext_gene
    values <- as.matrix(expression[, setdiff(names(expression), c("ensembl_gene_id", "ext_gene")), drop = FALSE])
  }
  samples <- match(make.names(colnames(dds)), make.names(colnames(values)))
  stopifnot(!anyNA(samples), !anyDuplicated(samples))
  values <- values[, samples, drop = FALSE]
  colnames(values) <- colnames(dds)
  symbols <- toupper(trimws(symbols))
  keep <- !is.na(symbols) & nzchar(symbols) & symbols != "IRF2BP2_SIGNATURE" & rowSums(!is.finite(values)) == 0
  values <- values[keep, , drop = FALSE]
  symbols <- symbols[keep]
  # Average VST values for loci mapping to the same symbol, as in the original scripts.
  groups <- unique(symbols)
  values <- rowsum(values, symbols, reorder = FALSE) / tabulate(match(symbols, groups), nbins = length(groups))
  deviations <- apply(values, 1, stats::sd)
  values <- values[is.finite(deviations) & deviations > 0, , drop = FALSE]
  values[order(rownames(values)), , drop = FALSE]
}

calculate_irf2bp2 <- function(config, dds, raw_data = getwd()) {
  stopifnot(identical(config$species, "human"), ncol(dds) >= 2)
  signature <- read_irf2bp2_signature(file.path(raw_data, "data/IRF2BP2_signature.tsv"))
  if (!is.null(config$signature_score_path)) {
    saved <- readRDS(paste0(raw_data, config$signature_score_path))
    samples <- match(make.names(colnames(dds)), make.names(names(saved$scores)))
    stopifnot(!anyNA(samples), !anyDuplicated(samples))
    scores <- saved$scores[samples]
    names(scores) <- colnames(dds)
    matched <- saved$matched_genes
    source <- "Saved BEAT-AML scores (original 193-sample scoring cohort)."
  } else {
    if (!requireNamespace("GSVA", quietly = TRUE)) stop("Install GSVA to calculate the IRF2BP2 signature score.")
    values <- irf2bp2_expression(config, dds, raw_data)
    matched <- intersect(signature, rownames(values))
    if (!length(matched)) stop("No variable IRF2BP2 signature genes were found in this dataset.")
    param <- GSVA::gsvaParam(exprData = values, geneSets = list(IRF2BP2 = signature),
                            kcdf = "Gaussian", maxDiff = TRUE)
    result <- GSVA::gsva(param, verbose = FALSE, BPPARAM = BiocParallel::SerialParam())
    scores <- result["IRF2BP2", ]
    source <- "Calculated from the selected dataset's VST expression."
  }
  stopifnot(length(scores) == ncol(dds), all(is.finite(scores)))
  list(scores = data.frame(Sample = colnames(dds), IRF2BP2_score = as.numeric(scores)),
       matched_genes = matched, signature_genes = signature, source = source)
}

irf2bp2_group_choices <- function(meta, config) {
  candidates <- unique(c(config$PCA_var, "KO_Sensitivity", "Treatment", "Treatment_Sensitivity", "Batch", "htb",
                         "condition", "Source", "prior_cr", "quantile", "FAB", "FAB_BlastMorphology",
                         "Denovo.Relapse", "RAS_mut", "NPM1_mut", "MClass", "Tissue"))
  candidates <- intersect(candidates, names(meta))
  candidates[vapply(meta[candidates], function(x) length(unique(x[!is.na(x)])) > 1, logical(1))]
}

IRF2BP2_UI <- function(id) {
  ns <- NS(id)
  fluidPage(
    titlePanel("IRF2BP2 Signature Score"),
    p("The 107-gene signature is scored with GSVA (Gaussian kernel, maxDiff = TRUE), following the existing scoring analyses."),
    p("Higher scores indicate greater enrichment of the signature genes. Scores depend on the scoring cohort; compare groups within a dataset."),
    sidebarLayout(
      sidebarPanel(
        actionButton(ns("calculate"), "Load / Calculate Score"),
        selectInput(ns("group_by"), "Group samples by:", choices = NULL),
        br(), downloadButton(ns("download_scores"), "Download Scores"),
        br(), br(), downloadButton(ns("download_plot"), "Download Plot"),
        br(), br(), downloadButton(ns("download_signature"), "Download Signature Genes")
      ),
      mainPanel(textOutput(ns("score_summary")),
                shinycssloaders::withSpinner(plotOutput(ns("score_plot"), height = "500px")),
                DT::DTOutput(ns("score_table")))
    )
  )
}

IRF2BP2_Server <- function(id, dataset_choice, dataset_dds, dataset_config, raw_data) {
  moduleServer(id, function(input, output, session) {
    scored <- reactiveVal(NULL)
    observeEvent(dataset_choice$user_dataset(), {
      scored(NULL)
      key <- dataset_choice$user_dataset()
      meta <- as.data.frame(SummarizedExperiment::colData(dataset_dds()))
      choices <- irf2bp2_group_choices(meta, dataset_config[[key]])
      updateSelectInput(session, "group_by", choices = choices, selected = choices[1])
    })
    observeEvent(input$calculate, {
      key <- dataset_choice$user_dataset()
      if (dataset_config[[key]]$species != "human") {
        showNotification("This human signature score is available for human datasets.", type = "message")
        return()
      }
      result <- tryCatch(withProgress(message = "Calculating IRF2BP2 signature score", value = 0.5, {
        calculate_irf2bp2(dataset_config[[key]], dataset_dds(), raw_data)
      }), error = function(e) {
        showNotification(conditionMessage(e), type = "error", duration = 15)
        NULL
      })
      if (!is.null(result)) scored(list(dataset = key, result = result))
    })
    current_score <- reactive({
      req(scored(), identical(scored()$dataset, dataset_choice$user_dataset()))
      scored()$result
    })
    score_data <- reactive({
      result <- current_score()$scores
      req(input$group_by)
      meta <- as.data.frame(SummarizedExperiment::colData(dataset_dds()))
      req(input$group_by %in% names(meta))
      result$Group <- as.character(meta[[input$group_by]][match(result$Sample, rownames(meta))])
      result
    })
    output$score_summary <- renderText({
      key <- dataset_choice$user_dataset()
      if (dataset_config[[key]]$species != "human") return("This signature is defined for human datasets.")
      if (is.null(scored()) || !identical(scored()$dataset, key)) return("Click Load / Calculate Score to view this dataset.")
      result <- current_score()
      paste(dataset_config[[key]]$label, "|", nrow(result$scores), "samples |",
            length(result$matched_genes), "of", length(result$signature_genes), "signature genes used.", result$source)
    })
    score_plot <- reactive({
      ggplot(score_data(), aes(x = Group, y = IRF2BP2_score, fill = Group)) +
        geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.15, alpha = 0.6) +
        theme_light() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
        labs(title = "IRF2BP2 Signature Score", subtitle = dataset_config[[dataset_choice$user_dataset()]]$label,
             x = input$group_by, y = "GSVA score")
    })
    output$score_plot <- renderPlot({ print(score_plot()) })
    output$score_table <- DT::renderDT({ DT::datatable(score_data(), rownames = FALSE, options = list(scrollX = TRUE)) })
    output$download_scores <- downloadHandler(
      filename = function() paste0(dataset_choice$user_dataset(), "_IRF2BP2_scores.csv"),
      content = function(file) write.csv(score_data(), file, row.names = FALSE))
    output$download_plot <- downloadHandler(
      filename = function() paste0(dataset_choice$user_dataset(), "_IRF2BP2_score.png"),
      content = function(file) ggsave(file, plot = score_plot(), width = 9, height = 6, dpi = 200))
    output$download_signature <- downloadHandler(
      filename = function() "IRF2BP2_signature.tsv",
      content = function(file) file.copy(file.path(raw_data, "data/IRF2BP2_signature.tsv"), file))
  })
}
