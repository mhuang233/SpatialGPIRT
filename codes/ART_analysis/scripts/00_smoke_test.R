source(file.path("config", "art_config.R"))
source(file.path("R", "00_utils.R"))
source_art_files(".")

art_config$results_dir <- file.path("results", "smoke_test")
art_config$n_folds <- 2L
art_config$chains <- 2L
art_config$parallel_chains <- 2L
art_config$iter_warmup <- 250L
art_config$iter_sampling <- 250L
art_config$refresh <- 50L
art_config$max_scoring_draws <- 100L
art_config$mds_bootstrap_replicates <- 20L

####### Test run #########
run_art_cv(
  config = art_config,
  run_label = "one_fold",
  models = art_config$primary_models,
  genres = "SciFiFantasy",
  folds_to_run = 1L,
  geometry_method = "residual",
  anisotropic = TRUE
)

