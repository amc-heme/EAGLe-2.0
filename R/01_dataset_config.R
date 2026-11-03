# Reports use the canonical dataset keys. Interactive modules retain their UI IDs.
# Both views read the same paths and plotting settings from datasets.yaml.
load_analysis_config <- function(path = "datasets.yaml") {
  config <- yaml::read_yaml(path)
  config <- Filter(function(x) !is.null(x$interactive_id), config)
  names(config) <- vapply(config, function(x) x$interactive_id, character(1))
  config
}
