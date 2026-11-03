# IRF2BP2 signature

The 107-gene list is copied from `BEAT-AML/docs/irf2bbp2_genes.tsv` in the user's existing analysis repository. Gene symbols use `gene_symbol`, falling back to `old_symbol`, and are uppercased and deduplicated.

Scoring follows `BEAT-AML/analysis/IRF2BP2_sig_scoring.Rmd` and `TCGA/analysis/IRF2BP2_sig_score.Rmd`: GSVA with `kcdf = "Gaussian"` and `maxDiff = TRUE` on VST expression, averaging rows mapping to the same symbol. The original scripts name a variable `ssgsea_scores`, but call GSVA rather than ssGSEA. All mapped variable genes form the background; the signature alone is not used as the expression background. Calculations use a serial backend.

`new_raw_data/BEAT_IRF2BP2_scores.rds` preserves the exact `IRF2BP2_SIGNATURE` row from `BEAT_vst2.rds` in saved Git stash `811c2da5524292468ec57b9767881158a171b95a` (`IRF2BP2_score_addition`). These scores use the original 193-sample cohort and 100 signature genes. BEAT FAB and disease-stage views select matching sample scores from that saved cohort. The stash is unchanged.

Other human datasets are scored within their selected cohorts. Coverage is reported in the app. Scores are cohort-dependent and are not calibrated for direct numerical comparisons between datasets. Mouse ortholog scoring is not implemented. Score plots do not use the gene-level differential-expression p-values.
