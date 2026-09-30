# Run only the published staircase 1PLUS, 2PLUS, and 3PLUS methods on the
# completed SGP-IRT primary five-fold CV design.
#
# IMPORTANT:
#   * This script does not call run_art_cv() and cannot fit any SGP-IRT,
#     independent-IRT, or matched-CAR variant.
#   * It requires the completed primary-CV geometry, split, and anchor/target
#     files and stops if any are missing.
#   * It uses the published staircase ir_spat() source with spat_model = "car"
#     and itemtype equal to exactly "1PLUS", "2PLUS", or "3PLUS". If the old
#     package cannot be installed on current R, the identical upstream
#     R/functions.R file is loaded directly from staircase-source/staircase-main.
#   * The original ART adaptation is preserved: raw Vocab_Sum defines the
#     package-required True_Species groups <40, 40--50, and >50.
#
# Run from the top-level art_revision directory:
#   source("scripts/06_run_exact_plus_primary_cv.R")

source(file.path("config", "art_config.R"))
source(file.path("R", "00_utils.R"))
source_art_files(".")

required_runtime_packages <- c("dismo", "dplyr", "plyr", "rgeos", "rstan")
missing_runtime_packages <- required_runtime_packages[
  !vapply(required_runtime_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_runtime_packages) > 0L) {
  stop(
    "Missing runtime dependencies for the exact staircase source: ",
    paste(missing_runtime_packages, collapse = ", "),
    call. = FALSE
  )
}

load_exact_ir_spat <- function() {
  if (requireNamespace("staircase", quietly = TRUE)) {
    return(list(
      fun = getExportedValue("staircase", "ir_spat"),
      version = as.character(utils::packageVersion("staircase")),
      implementation = "installed staircase::ir_spat"
    ))
  }

  source_candidates <- c(
    file.path("staircase-source", "staircase-main", "R", "functions.R"),
    file.path("staircase-main", "R", "functions.R")
  )
  source_path <- source_candidates[file.exists(source_candidates)][1L]
  assert_true(
    length(source_path) == 1L && !is.na(source_path),
    paste0(
      "The staircase package is unavailable and the exact source file was ",
      "not found. Expected staircase-source/staircase-main/R/functions.R."
    )
  )

  # Reproduce the bindings declared in staircase's published NAMESPACE. This
  # bypasses package lazy loading only; functions.R itself is not modified.
  source_environment <- new.env(parent = globalenv())
  imported_bindings <- list(
    "%>%" = getExportedValue("dplyr", "%>%"),
    case_when = getExportedValue("dplyr", "case_when"),
    distinct = getExportedValue("dplyr", "distinct"),
    left_join = getExportedValue("dplyr", "left_join"),
    mutate = getExportedValue("dplyr", "mutate"),
    "." = getExportedValue("plyr", "."),
    voronoi = getExportedValue("dismo", "voronoi"),
    gIntersection = getExportedValue("rgeos", "gIntersection"),
    stan = getExportedValue("rstan", "stan"),
    dist = getExportedValue("stats", "dist"),
    model.matrix = getExportedValue("stats", "model.matrix"),
    model.response = getExportedValue("stats", "model.response")
  )
  list2env(imported_bindings, envir = source_environment)
  sys.source(source_path, envir = source_environment, keep.source = TRUE)
  assert_true(
    exists("ir_spat", envir = source_environment, inherits = FALSE),
    "The published staircase functions.R file did not define ir_spat()."
  )

  description_path <- file.path(
    dirname(dirname(source_path)),
    "DESCRIPTION"
  )
  version <- if (file.exists(description_path)) {
    unname(read.dcf(description_path, fields = "Version")[1L, 1L])
  } else {
    "0.0.1-source"
  }
  list(
    fun = get("ir_spat", envir = source_environment, inherits = FALSE),
    version = version,
    implementation = paste0("direct published source: ", source_path)
  )
}

staircase_runtime <- load_exact_ir_spat()
EXACT_IR_SPAT <- staircase_runtime$fun
STAIRCASE_VERSION <- staircase_runtime$version
STAIRCASE_IMPLEMENTATION <- staircase_runtime$implementation
message("Exact PLUS implementation: ", STAIRCASE_IMPLEMENTATION)

# On Windows, the legacy RStan interface used by staircase can lose its PSOCK
# worker connection when chains are launched in parallel, producing
# serialize(..., node$con): error writing to connection.  Run the same four
# chains sequentially; this changes wall-clock scheduling, not the posterior.
options(mc.cores = 1L)
rstan::rstan_options(auto_write = TRUE)

EXACT_PLUS_MODELS <- c(
  exact_1plus = "1PLUS",
  exact_2plus = "2PLUS",
  exact_3plus = "3PLUS"
)

assert_true(
  identical(unname(EXACT_PLUS_MODELS), c("1PLUS", "2PLUS", "3PLUS")),
  "The exact PLUS allow-list was altered."
)

primary_dir <- file.path(art_config$results_dir, "primary_cv")
split_manifest_path <- file.path(primary_dir, "split_manifest.csv")
assert_true(
  file.exists(split_manifest_path),
  paste0("Missing completed primary-CV split manifest: ", split_manifest_path)
)

prepared <- prepare_art_data(
  path = art_config$data_file,
  genres = art_config$genres,
  vocabulary_column = art_config$vocabulary_column,
  cohort_rule = art_config$cohort_rule,
  audit_dir = NULL
)

splits <- make_three_way_splits(
  prepared = prepared,
  n_folds = art_config$n_folds,
  geometry_fraction = art_config$geometry_fraction_of_outer_training,
  seed = art_config$split_seed
)

expected_split_manifest <- do.call(rbind, lapply(splits, function(split) {
  rbind(
    data.frame(fold = split$fold, analysis_row = split$geometry,
               role = "geometry"),
    data.frame(fold = split$fold, analysis_row = split$calibration,
               role = "calibration"),
    data.frame(fold = split$fold, analysis_row = split$test,
               role = "test")
  )
}))
expected_split_manifest$subject_source_row <-
  prepared$data$.source_row[expected_split_manifest$analysis_row]

saved_split_manifest <- utils::read.csv(
  split_manifest_path,
  stringsAsFactors = FALSE
)
split_columns <- c("fold", "analysis_row", "role", "subject_source_row")
assert_true(
  all(split_columns %in% names(saved_split_manifest)),
  "The saved primary-CV split manifest has an unexpected schema."
)
sort_manifest <- function(x) {
  x <- x[, split_columns, drop = FALSE]
  x[order(x$fold, x$analysis_row, x$role), , drop = FALSE]
}
rownames(expected_split_manifest) <- NULL
saved_split_manifest <- sort_manifest(saved_split_manifest)
expected_split_manifest <- sort_manifest(expected_split_manifest)
rownames(saved_split_manifest) <- NULL
rownames(expected_split_manifest) <- NULL
assert_true(
  identical(saved_split_manifest, expected_split_manifest),
  paste0(
    "The current data/configuration do not recreate the completed SGP-IRT ",
    "five-fold split. No PLUS model was fitted."
  )
)

primary_genres <- names(art_config$genres)
primary_folds <- seq_len(art_config$n_folds)

required_primary_files <- unlist(lapply(primary_folds, function(fold) {
  unlist(lapply(primary_genres, function(genre) {
    base <- file.path(primary_dir, sprintf("fold_%02d", fold), genre)
    c(
      file.path(base, "geometry", "geometry.rds"),
      file.path(base, "anchor_target_manifest.csv"),
      file.path(base, "sgp_irt", "prediction_metrics.csv")
    )
  }), use.names = FALSE)
}), use.names = FALSE)

missing_primary_files <- required_primary_files[!file.exists(required_primary_files)]
assert_true(
  length(missing_primary_files) == 0L,
  paste0(
    "The completed SGP-IRT primary-CV artifacts are incomplete. Missing: ",
    paste(missing_primary_files, collapse = ", "),
    ". No PLUS model was fitted."
  )
)

vocabulary_group <- function(x) {
  factor(
    ifelse(x < 40, 1L, ifelse(x <= 50, 2L, 3L)),
    levels = 1:3
  )
}

make_exact_plus_data <- function(Y_calibration, raw_vocabulary, coords) {
  Y_calibration <- as.matrix(Y_calibration)
  coords <- as.matrix(coords)
  I <- nrow(Y_calibration)
  J <- ncol(Y_calibration)
  assert_true(
    identical(dim(coords), c(J, 2L)),
    "Saved MDS coordinates do not match the genre item count."
  )
  groups <- vocabulary_group(raw_vocabulary)
  assert_true(
    !anyNA(groups) && length(unique(as.integer(groups))) == 3L,
    "All three legacy vocabulary groups must occur in each calibration set."
  )
  data.frame(
    subject_id = rep(seq_len(I), each = J),
    item_id = rep(seq_len(J), times = I),
    response = as.integer(as.vector(t(Y_calibration))),
    lon = rep(coords[, 1L], times = I),
    lat = rep(coords[, 2L], times = I),
    True_Species = factor(rep(groups, each = J), levels = 1:3)
  )
}

anchor_target_from_manifest <- function(path, Y_test) {
  manifest <- utils::read.csv(path, stringsAsFactors = FALSE)
  required <- c("subject_source_row", "item_index", "role")
  assert_true(
    all(required %in% names(manifest)),
    paste0("Unexpected anchor/target manifest schema: ", path)
  )
  source_rows <- as.integer(rownames(Y_test))
  anchors <- vector("list", length(source_rows))
  targets <- vector("list", length(source_rows))
  for (i in seq_along(source_rows)) {
    piece <- manifest[manifest$subject_source_row == source_rows[i], , drop = FALSE]
    assert_true(
      nrow(piece) == ncol(Y_test),
      "The saved anchor/target manifest does not match the test response matrix."
    )
    anchors[[i]] <- sort(as.integer(piece$item_index[piece$role == "anchor"]))
    targets[[i]] <- sort(as.integer(piece$item_index[piece$role == "target"]))
    assert_true(
      length(anchors[[i]]) > 0L && length(targets[[i]]) > 0L &&
        length(intersect(anchors[[i]], targets[[i]])) == 0L,
      "Invalid saved anchor/target allocation."
    )
  }
  list(anchors = anchors, targets = targets)
}

extract_exact_plus_draws <- function(fit, itemtype, J, max_draws, seed) {
  extracted <- rstan::extract(
    fit,
    pars = c("difficulty", "diff_species", "alpha", "eta", "sigma_abil"),
    permuted = TRUE,
    inc_warmup = FALSE
  )
  difficulty <- as.matrix(extracted$difficulty)
  diff_species <- as.matrix(extracted$diff_species)
  alpha <- as.matrix(extracted$alpha)
  eta <- as.matrix(extracted$eta)
  sigma_abil <- as.numeric(extracted$sigma_abil)
  assert_true(ncol(difficulty) == J, "PLUS difficulty draws have the wrong item count.")
  assert_true(ncol(diff_species) == 3L, "PLUS species-difficulty draws must have three groups.")
  assert_true(length(sigma_abil) == nrow(difficulty), "Malformed PLUS ability-scale draws.")
  rows <- thin_draw_rows(nrow(difficulty), max_draws, seed)
  if (itemtype == "1PLUS") {
    alpha[,] <- 1
  }
  if (itemtype != "3PLUS") {
    eta[,] <- 0
  }
  list(
    difficulty = difficulty[rows, , drop = FALSE],
    diff_species = diff_species[rows, , drop = FALSE],
    alpha = alpha[rows, , drop = FALSE],
    eta = eta[rows, , drop = FALSE],
    sigma_abil = sigma_abil[rows]
  )
}

integrate_exact_plus_person <- function(
    y,
    group,
    anchor_items,
    target_items,
    draws,
    standard_normal_grid) {
  S <- nrow(draws$difficulty)
  Q <- length(standard_normal_grid)
  T <- length(target_items)
  prior_log_weight <- stats::dnorm(standard_normal_grid, log = TRUE)
  probability_by_draw <- matrix(NA_real_, S, T)

  for (s in seq_len(S)) {
    theta_grid <- draws$sigma_abil[s] * standard_normal_grid
    alpha_anchor <- as.numeric(draws$alpha[s, anchor_items])
    beta_anchor <-
      as.numeric(draws$difficulty[s, anchor_items]) +
      as.numeric(draws$diff_species[s, group])
    guessing_s <- as.numeric(draws$eta[s, group])

    eta_anchor <- sweep(
      outer(theta_grid, alpha_anchor),
      2L,
      alpha_anchor * beta_anchor,
      "-"
    )
    p_anchor <- guessing_s + (1 - guessing_s) * stats::plogis(eta_anchor)
    p_anchor <- clip_probability(p_anchor)
    y_anchor <- as.numeric(y[anchor_items])
    log_likelihood <- rowSums(
      sweep(log(p_anchor), 2L, y_anchor, "*") +
        sweep(log1p(-p_anchor), 2L, 1 - y_anchor, "*")
    )
    log_weight <- prior_log_weight + log_likelihood
    weight <- exp(log_weight - log_sum_exp(log_weight))

    alpha_target <- as.numeric(draws$alpha[s, target_items])
    beta_target <-
      as.numeric(draws$difficulty[s, target_items]) +
      as.numeric(draws$diff_species[s, group])
    eta_target <- sweep(
      outer(theta_grid, alpha_target),
      2L,
      alpha_target * beta_target,
      "-"
    )
    p_target <- guessing_s + (1 - guessing_s) * stats::plogis(eta_target)
    probability_by_draw[s, ] <- colSums(p_target * weight)
  }
  probability_by_draw
}

score_exact_plus_test <- function(
    fit,
    itemtype,
    Y_test,
    raw_vocabulary_test,
    anchor_target,
    theta_grid,
    max_draws,
    seed) {
  Y_test <- as.matrix(Y_test)
  groups <- as.integer(vocabulary_group(raw_vocabulary_test))
  draws <- extract_exact_plus_draws(
    fit = fit,
    itemtype = itemtype,
    J = ncol(Y_test),
    max_draws = max_draws,
    seed = seed
  )
  output <- vector("list", nrow(Y_test))
  for (i in seq_len(nrow(Y_test))) {
    target_items <- anchor_target$targets[[i]]
    probability_by_draw <- integrate_exact_plus_person(
      y = Y_test[i, ],
      group = groups[i],
      anchor_items = anchor_target$anchors[[i]],
      target_items = target_items,
      draws = draws,
      standard_normal_grid = theta_grid
    )
    intervals <- t(apply(probability_by_draw, 2L, shortest_interval))
    output[[i]] <- data.frame(
      subject_source_row = as.integer(rownames(Y_test)[i]),
      test_subject_index = i,
      item_index = target_items,
      item = colnames(Y_test)[target_items],
      y = as.integer(Y_test[i, target_items]),
      p_mean = colMeans(probability_by_draw),
      p_median = apply(probability_by_draw, 2L, stats::median),
      p_lower = intervals[, "lower"],
      p_upper = intervals[, "upper"],
      n_anchor = length(anchor_target$anchors[[i]])
    )
  }
  do.call(rbind, output)
}

exact_plus_diagnostics <- function(fit, fold, genre, model, wall_seconds) {
  fit_summary <- summary(fit)$summary
  sampler <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  divergences <- sum(vapply(
    sampler,
    function(x) sum(x[, "divergent__"]),
    numeric(1)
  ))
  treedepth_hits <- sum(vapply(
    sampler,
    # staircase::ir_spat() does not expose rstan's control argument, so the
    # exact package implementation uses rstan's default max_treedepth of 10.
    function(x) sum(x[, "treedepth__"] >= 10L),
    numeric(1)
  ))
  data.frame(
    fold = fold,
    genre = genre,
    model = model,
    itemtype = unname(EXACT_PLUS_MODELS[model]),
    max_rhat = max(fit_summary[, "Rhat"], na.rm = TRUE),
    min_n_eff = min(fit_summary[, "n_eff"], na.rm = TRUE),
    divergences = divergences,
    treedepth_hits = treedepth_hits,
    wall_seconds = wall_seconds,
    chains = art_config$chains,
    iter_warmup = art_config$iter_warmup,
    iter_sampling = art_config$iter_sampling,
    max_treedepth_setting = 10L,
    adapt_delta_setting = 0.8
  )
}

write_exact_plus_summary <- function() {
  metric_paths <- list.files(
    primary_dir,
    pattern = "^prediction_metrics\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )
  metrics <- do.call(rbind, lapply(
    metric_paths,
    utils::read.csv,
    stringsAsFactors = FALSE
  ))
  keep_models <- c("sgp_irt", names(EXACT_PLUS_MODELS))
  metrics <- metrics[metrics$model %in% keep_models, , drop = FALSE]
  expected_counts <- table(metrics$genre, metrics$model)
  assert_true(
    all(primary_genres %in% rownames(expected_counts)) &&
      all(keep_models %in% colnames(expected_counts)) &&
      all(expected_counts[primary_genres, keep_models, drop = FALSE] ==
            art_config$n_folds),
    "The SGP-IRT/exact-PLUS summary does not contain five folds per genre and model."
  )

  metric_columns <- c(
    "mean_log_score", "brier", "auc", "accuracy",
    "calibration_intercept", "calibration_slope", "ece_10",
    "wall_seconds"
  )
  model_summary <- mean_se_summary(
    metrics,
    value_columns = metric_columns,
    group_columns = c("genre", "model")
  )
  reference <- metrics[metrics$model == "sgp_irt", , drop = FALSE]
  comparison <- metrics[metrics$model %in% names(EXACT_PLUS_MODELS), , drop = FALSE]
  paired <- merge(
    comparison,
    reference,
    by = c("fold", "genre"),
    suffixes = c("_model", "_reference")
  )
  differences <- do.call(rbind, lapply(metric_columns, function(metric) {
    data.frame(
      genre = paired$genre,
      fold = paired$fold,
      model = paired$model_model,
      reference_model = "sgp_irt",
      metric = metric,
      difference_model_minus_reference =
        paired[[paste0(metric, "_model")]] -
        paired[[paste0(metric, "_reference")]]
    )
  }))
  difference_summary <- mean_se_summary(
    differences,
    value_columns = "difference_model_minus_reference",
    group_columns = c("genre", "model", "reference_model", "metric")
  )

  output_dir <- ensure_dir(file.path(
    art_config$results_dir,
    "exact_plus_primary_cv_summary"
  ))
  write_csv_atomic(metrics, file.path(output_dir, "all_fold_metrics.csv"))
  write_csv_atomic(model_summary, file.path(output_dir, "model_summary.csv"))
  write_csv_atomic(differences, file.path(output_dir, "paired_differences.csv"))
  write_csv_atomic(
    difference_summary,
    file.path(output_dir, "paired_difference_summary.csv")
  )
  invisible(output_dir)
}

run_manifest <- data.frame(
  implementation = STAIRCASE_IMPLEMENTATION,
  package_version = STAIRCASE_VERSION,
  itemtypes = paste(unname(EXACT_PLUS_MODELS), collapse = ";"),
  spatial_model = "car",
  cv_run = "primary_cv",
  folds = paste(primary_folds, collapse = ";"),
  genres = paste(primary_genres, collapse = ";"),
  true_species_rule = "Vocab_Sum <40; 40<=Vocab_Sum<=50; Vocab_Sum >50",
  chains = art_config$chains,
  iter_warmup = art_config$iter_warmup,
  iter_sampling = art_config$iter_sampling,
  max_scoring_draws = art_config$max_scoring_draws,
  stringsAsFactors = FALSE
)
write_csv_atomic(
  run_manifest,
  file.path(primary_dir, "exact_plus_run_manifest.csv")
)

for (fold in primary_folds) {
  split <- splits[[fold]]
  for (genre in primary_genres) {
    message("Exact PLUS primary CV: fold ", fold, ", genre ", genre)
    fold_genre_dir <- file.path(
      primary_dir,
      sprintf("fold_%02d", fold),
      genre
    )
    geometry <- readRDS(file.path(
      fold_genre_dir,
      "geometry",
      "geometry.rds"
    ))
    coords <- geometry$mds$coords
    Y <- genre_response_matrix(prepared, genre)
    Y_calibration <- Y[split$calibration, , drop = FALSE]
    Y_test <- Y[split$test, , drop = FALSE]
    raw_vocabulary <- prepared$data[[prepared$vocabulary_column]]
    plus_data <- make_exact_plus_data(
      Y_calibration = Y_calibration,
      raw_vocabulary = raw_vocabulary[split$calibration],
      coords = coords
    )
    anchor_target <- anchor_target_from_manifest(
      file.path(fold_genre_dir, "anchor_target_manifest.csv"),
      Y_test
    )

    for (model_index in seq_along(EXACT_PLUS_MODELS)) {
      model_name <- names(EXACT_PLUS_MODELS)[model_index]
      itemtype <- unname(EXACT_PLUS_MODELS[model_index])
      model_dir <- ensure_dir(file.path(fold_genre_dir, model_name))
      metric_path <- file.path(model_dir, "prediction_metrics.csv")
      if (isTRUE(art_config$skip_completed) && file.exists(metric_path)) {
        message("  Reusing completed ", itemtype)
        next
      }
      seed <- art_config$split_seed +
        100000L * fold +
        1000L * match(genre, primary_genres) +
        100L * model_index
      set.seed(seed)
      started <- Sys.time()
      fit <- EXACT_IR_SPAT(
        data = plus_data,
        spat_model = "car",
        itemtype = itemtype,
        abil = "subject_id",
        diff = "item_id",
        y = "response",
        coords = c("lon", "lat"),
        iter = art_config$iter_warmup + art_config$iter_sampling,
        warmup = art_config$iter_warmup,
        chains = art_config$chains,
        loglik = TRUE,
        refresh = art_config$refresh,
        seed = seed
      )
      wall_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
      if (isTRUE(art_config$save_fit_objects)) {
        saveRDS(fit, file.path(model_dir, "fit.rds"))
      }
      predictions <- score_exact_plus_test(
        fit = fit,
        itemtype = itemtype,
        Y_test = Y_test,
        raw_vocabulary_test = raw_vocabulary[split$test],
        anchor_target = anchor_target,
        theta_grid = art_config$theta_grid,
        max_draws = art_config$max_scoring_draws,
        seed = seed + 1L
      )
      metrics <- prediction_metrics(predictions)
      metrics$wall_seconds <- wall_seconds
      metrics$fold <- fold
      metrics$genre <- genre
      metrics$model <- model_name
      metrics$n_geometry <- length(split$geometry)
      metrics$n_calibration <- length(split$calibration)
      metrics$n_test <- length(split$test)
      metrics$geometry_method <- "residual"
      metrics$anisotropic <- NA
      metrics$a_A <- NA
      metrics$b_A <- NA

      item_metrics <- item_level_prediction_metrics(predictions)
      item_metrics$fold <- fold
      item_metrics$genre <- genre
      item_metrics$model <- model_name
      calibration <- calibration_curve_data(predictions)
      if (nrow(calibration) > 0L) {
        calibration$fold <- fold
        calibration$genre <- genre
        calibration$model <- model_name
      }
      diagnostics <- exact_plus_diagnostics(
        fit = fit,
        fold = fold,
        genre = genre,
        model = model_name,
        wall_seconds = wall_seconds
      )

      write_csv_atomic(
        predictions,
        file.path(model_dir, "target_predictions.csv")
      )
      write_csv_atomic(metrics, metric_path)
      write_csv_atomic(
        item_metrics,
        file.path(model_dir, "item_prediction_metrics.csv")
      )
      write_csv_atomic(
        calibration,
        file.path(model_dir, "calibration_curve.csv")
      )
      write_csv_atomic(
        diagnostics,
        file.path(model_dir, "exact_plus_mcmc_diagnostics.csv")
      )
    }
  }
}

summary_dir <- write_exact_plus_summary()
message(
  "Exact 1PLUS/2PLUS/3PLUS primary CV complete. Results: ",
  summary_dir
)
