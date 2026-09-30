source(file.path("config", "art_config.R"))
source(file.path("R", "00_utils.R"))
source_art_files(".")

run_art_cv(
  config = art_config,
  run_label = "primary_cv",
  models = art_config$primary_models,
  genres = names(art_config$genres),
  folds_to_run = seq_len(art_config$n_folds),
  geometry_method = "residual",
  anisotropic = TRUE
)

