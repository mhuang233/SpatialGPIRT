source(file.path("config", "art_config.R"))
source(file.path("R", "00_utils.R"))
source_art_files(".")

run_honest_descriptive_analysis(
  config = art_config,
  models = art_config$primary_models,
  geometry_method = "residual",
  anisotropic = TRUE
)

