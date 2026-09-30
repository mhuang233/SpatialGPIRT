source(file.path("config", "art_config.R"))
source(file.path("R", "00_utils.R"))
source_art_files(".")

# These runs are intentionally separate and resumable. They reuse the same
# deterministic respondent splits and anchor/target allocations.

run_art_cv(
  config = art_config,
  run_label = "sensitivity_isotropic",
  models = "sgp_irt",
  geometry_method = "residual",
  anisotropic = FALSE
)

run_art_cv(
  config = art_config,
  run_label = "sensitivity_longer_range_prior",
  models = "sgp_irt",
  geometry_method = "residual",
  anisotropic = TRUE,
  a_A = 2,
  b_A = 1
)

run_art_cv(
  config = art_config,
  run_label = "sensitivity_shorter_range_prior",
  models = "sgp_irt",
  geometry_method = "residual",
  anisotropic = TRUE,
  a_A = 2,
  b_A = 0.25
)

run_art_cv(
  config = art_config,
  run_label = "sensitivity_raw_response_mds",
  models = "sgp_irt",
  geometry_method = "raw",
  anisotropic = TRUE
)

run_art_cv(
  config = art_config,
  run_label = "sensitivity_plus_type",
  models = c("sgp_irt", art_config$supplemental_models),
  geometry_method = "residual",
  anisotropic = TRUE
)
