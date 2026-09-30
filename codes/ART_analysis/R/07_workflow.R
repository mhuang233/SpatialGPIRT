model_data_for_fit <- function(
    model_name,
    base_data,
    geometry,
    config,
    anisotropic = TRUE,
    a_A = config$a_A,
    b_A = config$b_A) {
  if (model_name == "ind_2pl") {
    return(base_data)
  }
  if (model_name == "sgp_irt") {
    return(add_sgp_data(
      base_data = base_data,
      coords = geometry$coords,
      anisotropic = anisotropic,
      a_A = a_A,
      b_A = b_A
    ))
  }
  if (model_name %in% c("car_2pl", "car_1pl", "car_3pl")) {
    return(add_car_data(base_data, geometry$graph))
  }
  stop("Unknown ART model: ", model_name, call. = FALSE)
}

load_or_construct_geometry <- function(
    Y_geometry,
    vocabulary_geometry,
    compiled_ind_model,
    output_dir,
    seed,
    config,
    method) {
  saved_path <- file.path(output_dir, "geometry.rds")
  if (isTRUE(config$skip_completed) && file.exists(saved_path)) {
    saved <- readRDS(saved_path)
    return(list(
      coords = saved$mds$coords,
      graph = saved$graph,
      residual = saved$residual,
      probability = saved$posterior_probability,
      mds = saved$mds,
      geometry_fit = NULL,
      method = saved$method
    ))
  }
  construct_item_geometry(
    Y_geometry = Y_geometry,
    vocabulary_geometry = vocabulary_geometry,
    compiled_ind_model = compiled_ind_model,
    output_dir = output_dir,
    seed = seed,
    config = config,
    method = method
  )
}

annotate_and_write_diagnostics <- function(
    diagnostic_result,
    fold,
    genre,
    model,
    output_dir,
    config = NULL) {
  convergence <- diagnostic_result$convergence
  convergence$fold <- fold
  convergence$genre <- genre
  convergence$model <- model
  if (!is.null(config)) {
    convergence$chains <- config$chains
    convergence$iter_warmup <- config$iter_warmup
    convergence$iter_sampling <- config$iter_sampling
    convergence$adapt_delta <- config$adapt_delta
    convergence$max_treedepth_setting <- config$max_treedepth
  }
  posterior <- diagnostic_result$posterior
  posterior$fold <- fold
  posterior$genre <- genre
  posterior$model <- model
  write_csv_atomic(convergence, file.path(output_dir, "mcmc_diagnostics.csv"))
  write_csv_atomic(posterior, file.path(output_dir, "posterior_summaries.csv"))
  list(convergence = convergence, posterior = posterior)
}

run_art_fold_genre <- function(
    prepared,
    genre,
    split,
    compiled_models,
    config,
    output_dir,
    models = config$primary_models,
    geometry_method = "residual",
    anisotropic = TRUE,
    a_A = config$a_A,
    b_A = config$b_A) {
  ensure_dir(output_dir)
  Y <- genre_response_matrix(prepared, genre)
  vocabulary <- prepared$data[[prepared$vocabulary_column]]
  Y_geometry <- Y[split$geometry, , drop = FALSE]
  Y_calibration <- Y[split$calibration, , drop = FALSE]
  Y_test <- Y[split$test, , drop = FALSE]

  geometry_dir <- file.path(output_dir, "geometry")
  geometry_seed <- config$split_seed +
    100000L * split$fold +
    1000L * match(genre, names(prepared$genres))
  geometry <- load_or_construct_geometry(
    Y_geometry = Y_geometry,
    vocabulary_geometry = vocabulary[split$geometry],
    compiled_ind_model = compiled_models$ind_2pl,
    output_dir = geometry_dir,
    seed = geometry_seed,
    config = config,
    method = geometry_method
  )
  if (!is.null(geometry$geometry_fit)) {
    geometry_diagnostics <- save_fit_diagnostics(
      fit = geometry$geometry_fit,
      model_name = "ind_2pl",
      output_dir = file.path(geometry_dir, "nonspatial_geometry_fit"),
      max_treedepth = config$max_treedepth,
      item_names = colnames(Y),
      subject_source_rows = prepared$data$.source_row[split$geometry]
    )
    annotate_and_write_diagnostics(
      diagnostic_result = geometry_diagnostics,
      fold = split$fold,
      genre = genre,
      model = "geometry_ind_2pl",
      output_dir = file.path(geometry_dir, "nonspatial_geometry_fit"),
      config = config
    )
  }

  scaling <- standardize_from_training(
    vocabulary[split$calibration],
    vocabulary[split$test]
  )
  write_csv_atomic(
    data.frame(
      fold = split$fold,
      genre = genre,
      center = scaling$center,
      scale = scaling$scale
    ),
    file.path(output_dir, "vocabulary_scaling.csv")
  )
  base_data <- long_irt_data(
    Y = Y_calibration,
    vocabulary_z = scaling$training,
    grainsize = config$grainsize,
    sigma_gamma_prior = config$sigma_gamma_prior
  )
  anchor_target <- make_anchor_target_split(
    Y_test,
    anchor_fraction = config$anchor_fraction,
    seed = geometry_seed + 17L
  )
  anchor_manifest <- do.call(rbind, lapply(seq_len(nrow(Y_test)), function(i) {
    data.frame(
      fold = split$fold,
      genre = genre,
      subject_source_row = as.integer(rownames(Y_test)[i]),
      item_index = seq_len(ncol(Y_test)),
      item = colnames(Y_test),
      role = ifelse(
        seq_len(ncol(Y_test)) %in% anchor_target$anchors[[i]],
        "anchor",
        "target"
      )
    )
  }))
  write_csv_atomic(anchor_manifest, file.path(output_dir, "anchor_target_manifest.csv"))

  results <- vector("list", length(models))
  names(results) <- models
  for (model_index in seq_along(models)) {
    model_name <- models[model_index]
    assert_true(model_name %in% names(compiled_models), paste0(
      "Model was not compiled: ", model_name
    ))
    model_dir <- file.path(output_dir, model_name)
    metric_path <- file.path(model_dir, "prediction_metrics.csv")
    if (isTRUE(config$skip_completed) && file.exists(metric_path)) {
      results[[model_name]] <- list(
        skipped = TRUE,
        metrics = utils::read.csv(metric_path, stringsAsFactors = FALSE)
      )
      next
    }
    stan_data <- model_data_for_fit(
      model_name = model_name,
      base_data = base_data,
      geometry = geometry,
      config = config,
      anisotropic = anisotropic,
      a_A = a_A,
      b_A = b_A
    )
    model_seed <- geometry_seed + 100L * model_index
    fit <- sample_art_model(
      model = compiled_models[[model_name]],
      model_name = model_name,
      stan_data = stan_data,
      output_dir = model_dir,
      seed = model_seed,
      config = config
    )
    diagnostic_result <- save_fit_diagnostics(
      fit = fit,
      model_name = model_name,
      output_dir = model_dir,
      max_treedepth = config$max_treedepth,
      item_names = colnames(Y),
      subject_source_rows = prepared$data$.source_row[split$calibration]
    )
    diagnostic_result <- annotate_and_write_diagnostics(
      diagnostic_result = diagnostic_result,
      fold = split$fold,
      genre = genre,
      model = model_name,
      output_dir = model_dir,
      config = config
    )
    predictions <- score_test_respondents(
      fit = fit,
      Y_test = Y_test,
      vocabulary_test_z = scaling$new,
      subject_source_rows = as.integer(rownames(Y_test)),
      anchor_target = anchor_target,
      theta_grid = config$theta_grid,
      max_draws = config$max_scoring_draws,
      seed = model_seed + 1L
    )
    metrics <- prediction_metrics(predictions)
    metrics$wall_seconds <- diagnostic_result$convergence$wall_seconds
    metrics$fold <- split$fold
    metrics$genre <- genre
    metrics$model <- model_name
    metrics$n_geometry <- length(split$geometry)
    metrics$n_calibration <- length(split$calibration)
    metrics$n_test <- length(split$test)
    metrics$geometry_method <- geometry_method
    metrics$anisotropic <- anisotropic
    metrics$a_A <- a_A
    metrics$b_A <- b_A

    item_metrics <- item_level_prediction_metrics(predictions)
    item_metrics$fold <- split$fold
    item_metrics$genre <- genre
    item_metrics$model <- model_name
    calibration <- calibration_curve_data(predictions)
    if (nrow(calibration) > 0L) {
      calibration$fold <- split$fold
      calibration$genre <- genre
      calibration$model <- model_name
    } else {
      calibration <- data.frame(
        bin = integer(),
        n = integer(),
        mean_predicted = numeric(),
        observed = numeric(),
        fold = integer(),
        genre = character(),
        model = character()
      )
    }

    write_csv_atomic(predictions, file.path(model_dir, "target_predictions.csv"))
    write_csv_atomic(metrics, metric_path)
    write_csv_atomic(item_metrics, file.path(model_dir, "item_prediction_metrics.csv"))
    write_csv_atomic(calibration, file.path(model_dir, "calibration_curve.csv"))
    results[[model_name]] <- list(
      skipped = FALSE,
      fit = fit,
      metrics = metrics,
      diagnostics = diagnostic_result,
      predictions = predictions
    )
  }
  invisible(list(
    geometry = geometry,
    vocabulary_scaling = scaling,
    anchor_target = anchor_target,
    models = results
  ))
}

write_split_manifest <- function(prepared, splits, path) {
  manifest <- do.call(rbind, lapply(splits, function(split) {
    rbind(
      data.frame(fold = split$fold, analysis_row = split$geometry, role = "geometry"),
      data.frame(
        fold = split$fold,
        analysis_row = split$calibration,
        role = "calibration"
      ),
      data.frame(fold = split$fold, analysis_row = split$test, role = "test")
    )
  }))
  manifest$subject_source_row <- prepared$data$.source_row[manifest$analysis_row]
  write_csv_atomic(manifest, path)
}

run_art_cv <- function(
    config,
    run_label = "primary_cv",
    models = config$primary_models,
    genres = names(config$genres),
    folds_to_run = seq_len(config$n_folds),
    geometry_method = "residual",
    anisotropic = TRUE,
    a_A = config$a_A,
    b_A = config$b_A,
    force_recompile = FALSE) {
  audit_dir <- file.path(config$results_dir, "data_audit")
  prepared <- prepare_art_data(
    path = config$data_file,
    genres = config$genres,
    vocabulary_column = config$vocabulary_column,
    cohort_rule = config$cohort_rule,
    audit_dir = audit_dir
  )
  splits <- make_three_way_splits(
    prepared = prepared,
    n_folds = config$n_folds,
    geometry_fraction = config$geometry_fraction_of_outer_training,
    seed = config$split_seed
  )
  run_dir <- file.path(config$results_dir, run_label)
  ensure_dir(run_dir)
  write_software_manifest(run_dir, config)
  write_csv_atomic(
    data.frame(
      run_label = run_label,
      models = paste(models, collapse = ";"),
      genres = paste(genres, collapse = ";"),
      folds = paste(folds_to_run, collapse = ";"),
      geometry_method = geometry_method,
      anisotropic = anisotropic,
      a_A = a_A,
      b_A = b_A
    ),
    file.path(run_dir, "design_manifest.csv")
  )
  write_split_manifest(prepared, splits, file.path(run_dir, "split_manifest.csv"))

  compile_names <- unique(c("ind_2pl", models))
  compiled <- compile_art_models(
    stan_dir = config$stan_dir,
    models = compile_names,
    force_recompile = force_recompile
  )
  for (fold in folds_to_run) {
    for (genre in genres) {
      message("ART CV: fold ", fold, ", genre ", genre)
      run_art_fold_genre(
        prepared = prepared,
        genre = genre,
        split = splits[[fold]],
        compiled_models = compiled,
        config = config,
        output_dir = file.path(
          run_dir,
          sprintf("fold_%02d", fold),
          genre
        ),
        models = models,
        geometry_method = geometry_method,
        anisotropic = anisotropic,
        a_A = a_A,
        b_A = b_A
      )
    }
  }
  summarize_primary_cv(
    cv_dir = run_dir,
    output_dir = file.path(config$results_dir, paste0(run_label, "_summary"))
  )
}

make_honest_two_way_split <- function(
    prepared,
    geometry_fraction = 0.5,
    seed = 20260730L) {
  vocabulary <- prepared$data[[prepared$vocabulary_column]]
  # The honest split is stratified only by the external vocabulary covariate,
  # never by ART responses.
  strata <- make_balance_strata(vocabulary)
  set.seed(seed)
  geometry <- integer(0)
  for (level in levels(strata)) {
    candidates <- which(strata == level)
    n_geometry <- max(
      1L,
      min(length(candidates) - 1L, round(geometry_fraction * length(candidates)))
    )
    geometry <- c(geometry, sample(candidates, n_geometry))
  }
  calibration <- setdiff(seq_len(nrow(prepared$data)), geometry)
  list(
    fold = 0L,
    geometry = sort(geometry),
    calibration = sort(calibration),
    test = integer(0L)
  )
}

fit_honest_genre_models <- function(
    prepared,
    genre,
    split,
    compiled,
    config,
    output_dir,
    models,
    geometry_method = "residual",
    anisotropic = TRUE,
    a_A = config$a_A,
    b_A = config$b_A) {
  ensure_dir(output_dir)
  Y <- genre_response_matrix(prepared, genre)
  vocabulary <- prepared$data[[prepared$vocabulary_column]]
  seed <- config$split_seed + 1000L * match(genre, names(prepared$genres))
  geometry <- load_or_construct_geometry(
    Y_geometry = Y[split$geometry, , drop = FALSE],
    vocabulary_geometry = vocabulary[split$geometry],
    compiled_ind_model = compiled$ind_2pl,
    output_dir = file.path(output_dir, "geometry"),
    seed = seed,
    config = config,
    method = geometry_method
  )
  if (!is.null(geometry$geometry_fit)) {
    geometry_diagnostic <- save_fit_diagnostics(
      fit = geometry$geometry_fit,
      model_name = "ind_2pl",
      output_dir = file.path(output_dir, "geometry", "nonspatial_geometry_fit"),
      max_treedepth = config$max_treedepth,
      item_names = colnames(Y),
      subject_source_rows = prepared$data$.source_row[split$geometry]
    )
    annotate_and_write_diagnostics(
      geometry_diagnostic,
      fold = 0L,
      genre = genre,
      model = "geometry_ind_2pl",
      output_dir = file.path(output_dir, "geometry", "nonspatial_geometry_fit"),
      config = config
    )
  }

  scaling <- standardize_from_training(vocabulary[split$calibration])
  base_data <- long_irt_data(
    Y = Y[split$calibration, , drop = FALSE],
    vocabulary_z = scaling$training,
    grainsize = config$grainsize,
    sigma_gamma_prior = config$sigma_gamma_prior
  )
  write_csv_atomic(
    data.frame(center = scaling$center, scale = scaling$scale),
    file.path(output_dir, "vocabulary_scaling.csv")
  )

  fits <- list()
  for (model_index in seq_along(models)) {
    model_name <- models[model_index]
    model_dir <- file.path(output_dir, model_name)
    stan_data <- model_data_for_fit(
      model_name,
      base_data,
      geometry,
      config,
      anisotropic,
      a_A,
      b_A
    )
    fit <- sample_art_model(
      model = compiled[[model_name]],
      model_name = model_name,
      stan_data = stan_data,
      output_dir = model_dir,
      seed = seed + 100L * model_index,
      config = config
    )
    diagnostic <- save_fit_diagnostics(
      fit = fit,
      model_name = model_name,
      output_dir = model_dir,
      max_treedepth = config$max_treedepth,
      item_names = colnames(Y),
      subject_source_rows = prepared$data$.source_row[split$calibration]
    )
    annotate_and_write_diagnostics(
      diagnostic,
      fold = 0L,
      genre = genre,
      model = model_name,
      output_dir = model_dir,
      config = config
    )
    ggplot2::ggsave(
      file.path(model_dir, "selected_trace_plots.pdf"),
      selected_trace_plot(fit, model_name, genre),
      width = 10,
      height = 6.5
    )
    ppc <- posterior_predictive_item_rates(
      fit = fit,
      Y = Y[split$calibration, , drop = FALSE],
      vocabulary_z = scaling$training,
      max_draws = config$max_scoring_draws,
      seed = seed + 100L * model_index + 1L
    )
    write_csv_atomic(
      ppc,
      file.path(model_dir, "posterior_predictive_item_rates.csv")
    )
    ggplot2::ggsave(
      file.path(model_dir, "posterior_predictive_item_rates.pdf"),
      plot_posterior_predictive_rates(ppc, genre, model_name),
      width = 5.5,
      height = 5
    )
    fits[[model_name]] <- fit
  }

  stability <- bootstrap_mds_stability(
    residual = geometry$residual,
    reference_coords = geometry$coords,
    shrinkage = config$residual_correlation_shrinkage,
    replicates = config$mds_bootstrap_replicates,
    seed = seed + 77L
  )
  write_csv_atomic(stability$item, file.path(output_dir, "mds_item_stability.csv"))
  write_csv_atomic(stability$summary, file.path(output_dir, "mds_stability_summary.csv"))

  sgp_fit <- fits$sgp_irt
  if (!is.null(sgp_fit)) {
    map_data <- item_map_data(sgp_fit, geometry$coords, colnames(Y))
    write_csv_atomic(map_data, file.path(output_dir, "sgp_irt", "item_map_data.csv"))
    map_plot <- plot_item_maps(map_data, genre)
    ggplot2::ggsave(
      file.path(output_dir, "sgp_irt", "difficulty_and_uncertainty_maps.pdf"),
      map_plot,
      width = 10,
      height = 4.8
    )

    omni <- rbind(
      posterior_semivariogram(sgp_fit, geometry$coords, "beta", FALSE, seed = seed),
      posterior_semivariogram(sgp_fit, geometry$coords, "f", FALSE, seed = seed + 1L)
    )
    directional <- rbind(
      posterior_semivariogram(sgp_fit, geometry$coords, "beta", TRUE, seed = seed),
      posterior_semivariogram(sgp_fit, geometry$coords, "f", TRUE, seed = seed + 1L)
    )
    write_csv_atomic(omni, file.path(output_dir, "sgp_irt", "semivariogram_omni.csv"))
    write_csv_atomic(
      directional,
      file.path(output_dir, "sgp_irt", "semivariogram_directional.csv")
    )
    ggplot2::ggsave(
      file.path(output_dir, "sgp_irt", "semivariogram_omni.pdf"),
      plot_semivariogram(omni, paste0(genre, ": omnidirectional")),
      width = 7.5,
      height = 4.8
    )
    ggplot2::ggsave(
      file.path(output_dir, "sgp_irt", "semivariogram_directional.pdf"),
      plot_semivariogram(directional, paste0(genre, ": directional")),
      width = 8.5,
      height = 5.2
    )

    curves <- posterior_icc_tic(sgp_fit, seed = seed)
    write_csv_atomic(curves, file.path(output_dir, "sgp_irt", "icc_tic_data.csv"))
    ggplot2::ggsave(
      file.path(output_dir, "sgp_irt", "icc_tic.pdf"),
      plot_icc_tic(curves, genre),
      width = 9,
      height = 4.2
    )
  }
  ggplot2::ggsave(
    file.path(output_dir, "mds_bootstrap_stability.pdf"),
    plot_mds_stability(geometry$coords, stability$item, genre),
    width = 5.5,
    height = 5
  )
  invisible(list(
    geometry = geometry,
    stability = stability,
    fits = fits,
    vocabulary_scaling = scaling
  ))
}

run_honest_descriptive_analysis <- function(
    config,
    models = config$primary_models,
    geometry_method = "residual",
    anisotropic = TRUE,
    force_recompile = FALSE) {
  prepared <- prepare_art_data(
    path = config$data_file,
    genres = config$genres,
    vocabulary_column = config$vocabulary_column,
    cohort_rule = config$cohort_rule,
    audit_dir = file.path(config$results_dir, "data_audit")
  )
  split <- make_honest_two_way_split(
    prepared,
    geometry_fraction = 0.5,
    seed = config$split_seed
  )
  output_dir <- file.path(config$results_dir, "honest_descriptive")
  ensure_dir(output_dir)
  write_software_manifest(output_dir, config)
  manifest <- rbind(
    data.frame(analysis_row = split$geometry, role = "geometry"),
    data.frame(analysis_row = split$calibration, role = "calibration")
  )
  manifest$subject_source_row <- prepared$data$.source_row[manifest$analysis_row]
  write_csv_atomic(manifest, file.path(output_dir, "split_manifest.csv"))

  compiled <- compile_art_models(
    stan_dir = config$stan_dir,
    models = unique(c("ind_2pl", models)),
    force_recompile = force_recompile
  )
  results <- list()
  for (genre in names(config$genres)) {
    message("Honest descriptive ART analysis: ", genre)
    results[[genre]] <- fit_honest_genre_models(
      prepared = prepared,
      genre = genre,
      split = split,
      compiled = compiled,
      config = config,
      output_dir = file.path(output_dir, genre),
      models = models,
      geometry_method = geometry_method,
      anisotropic = anisotropic
    )
  }

  map_plots <- lapply(names(results), function(genre) {
    result <- results[[genre]]
    if (is.null(result$fits$sgp_irt)) {
      return(NULL)
    }
    Y <- genre_response_matrix(prepared, genre)
    map_data <- item_map_data(
      result$fits$sgp_irt,
      result$geometry$coords,
      colnames(Y)
    )
    plot_item_maps(map_data, genre)
  })
  map_plots <- Filter(Negate(is.null), map_plots)
  if (length(map_plots) > 0L) {
    ggplot2::ggsave(
      file.path(output_dir, "ART_MDS_item_difficulty_maps_faceted.pdf"),
      patchwork::wrap_plots(map_plots, ncol = 1L),
      width = 10,
      height = 4.8 * length(map_plots)
    )
  }

  variogram_plots <- lapply(names(results), function(genre) {
    result <- results[[genre]]
    if (is.null(result$fits$sgp_irt)) {
      return(NULL)
    }
    seed <- config$split_seed + 1000L * match(genre, names(prepared$genres))
    omni <- rbind(
      posterior_semivariogram(
        result$fits$sgp_irt,
        result$geometry$coords,
        "beta",
        FALSE,
        seed = seed
      ),
      posterior_semivariogram(
        result$fits$sgp_irt,
        result$geometry$coords,
        "f",
        FALSE,
        seed = seed + 1L
      )
    )
    plot_semivariogram(omni, genre)
  })
  variogram_plots <- Filter(Negate(is.null), variogram_plots)
  if (length(variogram_plots) > 0L) {
    ggplot2::ggsave(
      file.path(output_dir, "ART_MDS_semivariograms_faceted.pdf"),
      patchwork::wrap_plots(variogram_plots, ncol = 1L),
      width = 8,
      height = 4.8 * length(variogram_plots)
    )
  }

  posterior <- read_result_csvs(output_dir, "posterior_summaries.csv")
  diagnostics <- read_result_csvs(output_dir, "mcmc_diagnostics.csv")
  write_csv_atomic(posterior, file.path(output_dir, "all_posterior_summaries.csv"))
  write_csv_atomic(diagnostics, file.path(output_dir, "all_mcmc_diagnostics.csv"))

  covariance_parameters <- posterior$parameter %in% c(
    "sigma_f", "sigma_u", "ell", "nugget_ratio", "anisotropy_ratio", "rho"
  )
  vocabulary_effect <- posterior$parameter == "gamma"
  write_csv_atomic(
    posterior[covariance_parameters, , drop = FALSE],
    file.path(output_dir, "spatial_parameter_summary.csv")
  )
  write_csv_atomic(
    posterior[vocabulary_effect, , drop = FALSE],
    file.path(output_dir, "vocabulary_effect_summary.csv")
  )
  invisible(results)
}
