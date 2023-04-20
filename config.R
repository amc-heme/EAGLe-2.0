# DE config 
config <- list(
  base_dir = "~/Documents/CD files/salmon",
  t2g_hs_file = "~/Documents/CD files/HS_transcript_to_gene.txt",
  metadata_file = "~/Documents/CD files/SampleSheet.txt",
  sample_id = c("SRR9265370", "SRR9265373", "SRR9265371", "SRR9265372", "SRR9265363", "SRR9265364", "SRR9265366", "SRR9265367", "SRR9265369", "SRR9265365", "SRR9265368", "SRR9265374"),
  design_formula = ~batch + condition
)
