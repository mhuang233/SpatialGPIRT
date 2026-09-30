art_config <- list(
  data_file = file.path("data", "psy.csv"),
  results_dir = "results",
  stan_dir = "stan",

  # Reproduces the N = 321 cohort described in the current manuscript.
  # Alternatives: "filter_flag" (N = 306 here), "all_complete" (N = 337).
  cohort_rule = "manuscript_321",

  genres = c(
    SciFiFantasy = "SciFiFantasy",
    SocialPoliticalCommentary = "SocialPoliticalCommentary",
    Suspense = "Suspense"
  ),
  vocabulary_column = "Vocab_Sum",

  # Honest validation design.
  n_folds = 5L,
  geometry_fraction_of_outer_training = 0.50,
  anchor_fraction = 0.50,
  split_seed = 20260730L,

  # MDS based on residual item associations from the geometry respondents.
  mds_dimensions = 2L,
  residual_correlation_shrinkage = 0.05,
  car_starting_neighbors = 3L,
  mds_bootstrap_replicates = 200L,

  # Shared priors.
  sigma_gamma_prior = 1.0,
  a_A = 2.0,
  b_A = 0.5,

  # CmdStanR defaults for reportable fits.
  chains = 4L,
  parallel_chains = 4L,
  threads_per_chain = 1L,
  iter_warmup = 1000L,
  iter_sampling = 1000L,
  adapt_delta = 0.99,
  max_treedepth = 15L,
  grainsize = 1L,
  refresh = 100L,
  max_scoring_draws = 500L,
  theta_grid = seq(-4.5, 4.5, length.out = 81L),

  primary_models = c("ind_2pl", "car_2pl", "sgp_irt"),
  supplemental_models = c("car_1pl", "car_3pl"),

  save_fit_objects = TRUE,
  skip_completed = TRUE
)

