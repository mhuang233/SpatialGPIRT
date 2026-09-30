pair_geometry <- function(coords) {
  coords <- as.matrix(coords)
  pairs <- which(upper.tri(matrix(0, nrow(coords), nrow(coords))), arr.ind = TRUE)
  delta1 <- coords[pairs[, 2L], 1L] - coords[pairs[, 1L], 1L]
  delta2 <- coords[pairs[, 2L], 2L] - coords[pairs[, 1L], 2L]
  angle <- atan2(delta2, delta1) %% pi
  direction <- ifelse(
    angle < pi / 8 | angle >= 7 * pi / 8,
    "MDS1 direction",
    ifelse(
      angle >= 3 * pi / 8 & angle < 5 * pi / 8,
      "MDS2 direction",
      "Diagonal directions"
    )
  )
  data.frame(
    item1 = pairs[, 1L],
    item2 = pairs[, 2L],
    distance = sqrt(delta1^2 + delta2^2),
    angle = angle,
    direction = direction
  )
}

variogram_bin_definition <- function(pair_data, n_bins = NULL) {
  n_pairs <- nrow(pair_data)
  n_bins <- n_bins %||% max(3L, min(7L, floor(sqrt(n_pairs))))
  breaks <- unique(stats::quantile(
    pair_data$distance,
    probs = seq(0, 1, length.out = n_bins + 1L),
    names = FALSE
  ))
  if (length(breaks) < 3L) {
    breaks <- seq(min(pair_data$distance), max(pair_data$distance), length.out = 4L)
  }
  breaks[1L] <- -Inf
  breaks[length(breaks)] <- Inf
  cut(
    pair_data$distance,
    breaks = breaks,
    include.lowest = TRUE,
    labels = FALSE
  )
}

semivariogram_from_values <- function(
    coords,
    values,
    directional = FALSE,
    n_bins = NULL) {
  pairs <- pair_geometry(coords)
  pairs$bin <- variogram_bin_definition(pairs, n_bins)
  pairs$semivariance <- 0.5 * (
    values[pairs$item1] - values[pairs$item2]
  )^2
  if (!directional) {
    pairs$direction <- "Omnidirectional"
  }
  groups <- interaction(pairs$bin, pairs$direction, drop = TRUE)
  pieces <- split(pairs, groups)
  do.call(rbind, lapply(pieces, function(piece) {
    data.frame(
      bin = piece$bin[1L],
      direction = piece$direction[1L],
      distance = mean(piece$distance),
      semivariance = mean(piece$semivariance),
      n_pairs = nrow(piece)
    )
  }))
}

posterior_semivariogram <- function(
    fit,
    coords,
    variable = c("beta", "f"),
    directional = FALSE,
    max_draws = 1000L,
    seed = 1L) {
  variable <- match.arg(variable)
  values <- indexed_draws(fit, variable, nrow(coords))
  rows <- thin_draw_rows(nrow(values), max_draws, seed)
  values <- values[rows, , drop = FALSE]
  reference <- semivariogram_from_values(
    coords,
    colMeans(values),
    directional = directional,
    n_bins = if (directional) 3L else NULL
  )
  reference$plugin_posterior_mean_semivariance <- reference$semivariance
  key <- paste(reference$bin, reference$direction)
  semivariance_draws <- matrix(
    NA_real_,
    nrow = nrow(values),
    ncol = nrow(reference)
  )
  for (s in seq_len(nrow(values))) {
    current <- semivariogram_from_values(
      coords,
      values[s, ],
      directional = directional,
      n_bins = if (directional) 3L else NULL
    )
    current_key <- paste(current$bin, current$direction)
    semivariance_draws[s, ] <- current$semivariance[match(key, current_key)]
  }
  intervals <- t(apply(semivariance_draws, 2L, shortest_interval))
  reference$posterior_mean <- colMeans(semivariance_draws, na.rm = TRUE)
  reference$lower <- intervals[, "lower"]
  reference$upper <- intervals[, "upper"]
  reference$component <- variable
  reference
}

plot_semivariogram <- function(data, title = NULL) {
  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = distance,
      y = posterior_mean,
      color = direction,
      fill = direction
    )
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper),
      alpha = 0.16,
      linewidth = 0,
      show.legend = FALSE
    ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_line(
      ggplot2::aes(y = plugin_posterior_mean_semivariance),
      linewidth = 0.6,
      linetype = 2
    ) +
    ggplot2::geom_point(ggplot2::aes(size = n_pairs)) +
    ggplot2::facet_wrap(~component, scales = "free_y") +
    ggplot2::labs(
      title = title,
      x = "MDS distance",
      y = "Semivariance",
      color = NULL,
      size = "Pairs",
      caption = "Solid: posterior expected semivariogram; dashed: variogram of posterior mean."
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")
}

item_map_data <- function(fit, coords, item_names) {
  beta <- indexed_draws(fit, "beta", length(item_names))
  f <- indexed_draws(fit, "f", length(item_names))
  beta_interval <- t(apply(beta, 2L, shortest_interval))
  f_interval <- t(apply(f, 2L, shortest_interval))
  data.frame(
    item = item_names,
    MDS1 = coords[, 1L],
    MDS2 = coords[, 2L],
    beta_mean = colMeans(beta),
    beta_lower = beta_interval[, "lower"],
    beta_upper = beta_interval[, "upper"],
    beta_width = beta_interval[, "upper"] - beta_interval[, "lower"],
    f_mean = colMeans(f),
    f_lower = f_interval[, "lower"],
    f_upper = f_interval[, "upper"],
    f_width = f_interval[, "upper"] - f_interval[, "lower"]
  )
}

plot_item_maps <- function(map_data, genre) {
  mean_plot <- ggplot2::ggplot(
    map_data,
    ggplot2::aes(MDS1, MDS2, color = beta_mean)
  ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_text(
      ggplot2::aes(label = sub("^[^0-9]+", "", item)),
      nudge_y = 0.035,
      size = 2.5,
      check_overlap = TRUE
    ) +
    ggplot2::scale_color_viridis_c(option = "D") +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = paste0(genre, ": posterior mean difficulty"),
      x = "MDS dimension 1",
      y = "MDS dimension 2",
      color = expression(E(beta[j] ~ "|" ~ Y))
    ) +
    ggplot2::theme_minimal(base_size = 11)

  uncertainty_plot <- ggplot2::ggplot(
    map_data,
    ggplot2::aes(MDS1, MDS2, color = beta_width)
  ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_text(
      ggplot2::aes(label = sub("^[^0-9]+", "", item)),
      nudge_y = 0.035,
      size = 2.5,
      check_overlap = TRUE
    ) +
    ggplot2::scale_color_viridis_c(option = "C") +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = paste0(genre, ": 95% interval width"),
      x = "MDS dimension 1",
      y = "MDS dimension 2",
      color = "Width"
    ) +
    ggplot2::theme_minimal(base_size = 11)
  mean_plot + uncertainty_plot
}

posterior_icc_tic <- function(
    fit,
    theta_grid = seq(-3, 3, length.out = 201L),
    vocabulary_z = 0,
    max_draws = 1000L,
    seed = 1L) {
  alpha <- indexed_draws(fit, "alpha")
  beta <- indexed_draws(fit, "beta", ncol(alpha))
  gamma <- indexed_draws(fit, "gamma", 1L)[, 1L]
  guessing <- indexed_draws(fit, "guessing", ncol(alpha))
  rows <- thin_draw_rows(nrow(alpha), max_draws, seed)
  alpha <- alpha[rows, , drop = FALSE]
  beta <- beta[rows, , drop = FALSE]
  gamma <- gamma[rows]
  guessing <- guessing[rows, , drop = FALSE]

  icc_draw <- matrix(0, nrow = length(rows), ncol = length(theta_grid))
  tic_draw <- matrix(0, nrow = length(rows), ncol = length(theta_grid))
  for (s in seq_along(rows)) {
    # A one-row subset of posterior::draws_matrix can retain matrix-like
    # attributes.  Coerce each draw explicitly so outer() always returns the
    # intended length(theta_grid)-by-J matrix rather than a higher-order array.
    alpha_s <- as.numeric(alpha[s, , drop = TRUE])
    beta_s <- as.numeric(beta[s, , drop = TRUE])
    guessing_s <- as.numeric(guessing[s, , drop = TRUE])
    assert_true(
      length(alpha_s) == ncol(alpha) &&
        length(beta_s) == ncol(alpha) &&
        length(guessing_s) == ncol(alpha),
      "ICC/TIC posterior draw dimensions do not match the number of items."
    )
    eta <- outer(as.numeric(theta_grid), alpha_s) -
      matrix(beta_s, nrow = length(theta_grid), ncol = ncol(alpha), byrow = TRUE) +
      as.numeric(vocabulary_z) * as.numeric(gamma[s])
    base_probability <- stats::plogis(eta)
    probability <- sweep(base_probability, 2L, 1 - guessing_s, "*")
    probability <- sweep(probability, 2L, guessing_s, "+")
    icc_draw[s, ] <- rowMeans(probability)

    # For a 3PL curve this is the exact Bernoulli information in theta.
    derivative <- sweep(
      base_probability * (1 - base_probability),
      2L,
      (1 - guessing_s) * alpha_s,
      "*"
    )
    tic_draw[s, ] <- rowSums(
      derivative^2 / (probability * (1 - probability) + 1e-12)
    )
  }
  icc_interval <- t(apply(icc_draw, 2L, shortest_interval))
  tic_interval <- t(apply(tic_draw, 2L, shortest_interval))
  rbind(
    data.frame(
      theta = theta_grid,
      curve = "Average endorsement probability",
      mean = colMeans(icc_draw),
      lower = icc_interval[, "lower"],
      upper = icc_interval[, "upper"]
    ),
    data.frame(
      theta = theta_grid,
      curve = "Test information",
      mean = colMeans(tic_draw),
      lower = tic_interval[, "lower"],
      upper = tic_interval[, "upper"]
    )
  )
}

plot_icc_tic <- function(curves, genre) {
  ggplot2::ggplot(
    curves,
    ggplot2::aes(theta, mean)
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper),
      fill = "#2C7FB8",
      alpha = 0.18
    ) +
    ggplot2::geom_line(color = "#2C7FB8", linewidth = 0.9) +
    ggplot2::facet_wrap(~curve, scales = "free_y", nrow = 1L) +
    ggplot2::labs(
      title = genre,
      x = expression(theta),
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11)
}

plot_mds_stability <- function(coords, item_stability, genre) {
  plot_data <- merge(
    data.frame(item = rownames(coords), coords, row.names = NULL),
    item_stability,
    by = "item",
    sort = FALSE
  )
  ggplot2::ggplot(plot_data, ggplot2::aes(MDS1, MDS2)) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = MDS1 - 1.96 * mds1_sd,
        xend = MDS1 + 1.96 * mds1_sd,
        yend = MDS2
      ),
      color = "grey60"
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        y = MDS2 - 1.96 * mds2_sd,
        yend = MDS2 + 1.96 * mds2_sd,
        xend = MDS1
      ),
      color = "grey60"
    ) +
    ggplot2::geom_point(size = 2.5, color = "#2C7FB8") +
    ggplot2::geom_text(
      ggplot2::aes(label = sub("^[^0-9]+", "", item)),
      nudge_y = 0.035,
      size = 2.5
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = paste0(genre, ": bootstrap MDS stability"),
      x = "MDS dimension 1",
      y = "MDS dimension 2"
    ) +
    ggplot2::theme_minimal(base_size = 11)
}

selected_trace_plot <- function(fit, model_name, genre) {
  variables <- switch(
    model_name,
    sgp_irt = c(
      "sigma_a", "sigma_f", "sigma_u", "ell",
      "nugget_ratio", "anisotropy_ratio", "gamma"
    ),
    car_2pl = c(
      "sigma_a", "sigma_f", "sigma_u", "rho", "nugget_ratio", "gamma"
    ),
    ind_2pl = c("sigma_a", "sigma_beta", "gamma"),
    car_1pl = c("sigma_f", "sigma_u", "rho", "nugget_ratio", "gamma"),
    car_3pl = c(
      "sigma_a", "sigma_f", "sigma_u", "rho", "nugget_ratio", "gamma"
    )
  )
  draws <- posterior::as_draws_df(fit$draws(variables = variables))
  value_columns <- setdiff(
    names(draws),
    c(".chain", ".iteration", ".draw")
  )
  long <- tidyr::pivot_longer(
    as.data.frame(draws),
    cols = dplyr::all_of(value_columns),
    names_to = "parameter",
    values_to = "value"
  )
  ggplot2::ggplot(
    long,
    ggplot2::aes(.iteration, value, group = .chain, color = factor(.chain))
  ) +
    ggplot2::geom_line(linewidth = 0.25, alpha = 0.7) +
    ggplot2::facet_wrap(~parameter, scales = "free_y") +
    ggplot2::labs(
      title = paste(genre, model_name, "selected trace plots"),
      x = "Post-warmup iteration",
      y = NULL,
      color = "Chain"
    ) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(legend.position = "bottom")
}

posterior_predictive_item_rates <- function(
    fit,
    Y,
    vocabulary_z,
    max_draws = 500L,
    seed = 1L) {
  Y <- as.matrix(Y)
  vocabulary_z <- as.numeric(vocabulary_z)
  I <- nrow(Y)
  J <- ncol(Y)
  theta <- indexed_draws(fit, "theta", I)
  alpha <- indexed_draws(fit, "alpha", J)
  beta <- indexed_draws(fit, "beta", J)
  gamma <- as.numeric(indexed_draws(fit, "gamma", 1L)[, 1L])
  guessing <- indexed_draws(fit, "guessing", J)
  assert_true(ncol(theta) == I, "Theta draws do not match the number of respondents.")
  assert_true(ncol(alpha) == J, "Alpha draws do not match the number of items.")
  assert_true(ncol(beta) == J, "Beta draws do not match the number of items.")
  assert_true(ncol(guessing) == J, "Guessing draws do not match the number of items.")
  assert_true(
    nrow(alpha) == nrow(theta) &&
      nrow(beta) == nrow(theta) &&
      nrow(guessing) == nrow(theta) &&
      length(gamma) == nrow(theta),
    "Posterior variables do not contain the same number of draws."
  )
  rows <- thin_draw_rows(nrow(theta), max_draws, seed)
  replicated_rate <- matrix(NA_real_, nrow = length(rows), ncol = J)
  expected_rate <- matrix(NA_real_, nrow = length(rows), ncol = J)
  set.seed(seed)

  for (s_index in seq_along(rows)) {
    s <- rows[s_index]
    theta_s <- as.numeric(theta[s, , drop = TRUE])
    alpha_s <- as.numeric(alpha[s, , drop = TRUE])
    beta_s <- as.numeric(beta[s, , drop = TRUE])
    guessing_s <- as.numeric(guessing[s, , drop = TRUE])
    eta <- tcrossprod(theta_s, alpha_s)
    eta <- eta -
      matrix(beta_s, nrow = I, ncol = J, byrow = TRUE) +
      matrix(
        vocabulary_z * gamma[s],
        nrow = I,
        ncol = J,
        byrow = FALSE
      )
    assert_true(
      identical(dim(eta), c(I, J)),
      "The reconstructed predictor is not an I-by-J matrix."
    )
    probability <- stats::plogis(eta) *
      matrix(1 - guessing_s, nrow = I, ncol = J, byrow = TRUE) +
      matrix(guessing_s, nrow = I, ncol = J, byrow = TRUE)
    expected_rate[s_index, ] <- colMeans(probability)
    replicated_rate[s_index, ] <- colMeans(matrix(
      stats::rbinom(I * J, 1L, as.vector(probability)),
      nrow = I,
      ncol = J
    ))
  }
  intervals <- t(apply(replicated_rate, 2L, shortest_interval))
  observed_rate <- colMeans(Y)
  lower_probability <- colMeans(
    sweep(replicated_rate, 2L, observed_rate, "<=")
  )
  upper_probability <- colMeans(
    sweep(replicated_rate, 2L, observed_rate, ">=")
  )
  data.frame(
    item_index = seq_len(J),
    item = colnames(Y),
    observed_rate = observed_rate,
    posterior_expected_rate = colMeans(expected_rate),
    replicated_lower = intervals[, "lower"],
    replicated_upper = intervals[, "upper"],
    posterior_predictive_p_two_sided = pmin(
      1,
      2 * pmin(lower_probability, upper_probability)
    )
  )
}

plot_posterior_predictive_rates <- function(data, genre, model_name) {
  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = posterior_expected_rate,
      y = observed_rate,
      label = sub("^[^0-9]+", "", item)
    )
  ) +
    ggplot2::geom_abline(
      intercept = 0,
      slope = 1,
      linetype = 2,
      color = "grey50"
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = replicated_lower,
        xend = replicated_upper,
        yend = observed_rate
      ),
      color = "#2C7FB8",
      alpha = 0.65
    ) +
    ggplot2::geom_point(color = "#2C7FB8", size = 2) +
    ggplot2::geom_text(nudge_y = 0.02, size = 2.5, check_overlap = TRUE) +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(
      title = paste(genre, model_name, "item-rate posterior predictive check"),
      x = "Posterior expected / replicated-rate interval",
      y = "Observed endorsement rate"
    ) +
    ggplot2::theme_minimal(base_size = 10)
}

read_result_csvs <- function(root, filename) {
  paths <- list.files(
    root,
    pattern = paste0("^", filename, "$"),
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(paths) == 0L) {
    return(data.frame())
  }
  do.call(rbind, lapply(paths, utils::read.csv, stringsAsFactors = FALSE))
}

mean_se_summary <- function(data, value_columns, group_columns) {
  key <- do.call(
    interaction,
    c(unname(data[, group_columns, drop = FALSE]), list(drop = TRUE))
  )
  pieces <- split(data, key)
  do.call(rbind, lapply(pieces, function(piece) {
    base <- piece[1L, group_columns, drop = FALSE]
    for (column in value_columns) {
      value <- piece[[column]]
      base[[paste0(column, "_mean")]] <- mean(value, na.rm = TRUE)
      base[[paste0(column, "_se")]] <- stats::sd(value, na.rm = TRUE) /
        sqrt(sum(is.finite(value)))
    }
    base$n_folds <- nrow(piece)
    base
  }))
}

summarize_primary_cv <- function(
    cv_dir,
    output_dir,
    reference_model = "sgp_irt") {
  ensure_dir(output_dir)
  metrics <- read_result_csvs(cv_dir, "prediction_metrics.csv")
  diagnostics <- read_result_csvs(cv_dir, "mcmc_diagnostics.csv")
  posterior <- read_result_csvs(cv_dir, "posterior_summaries.csv")
  assert_true(nrow(metrics) > 0L, "No CV prediction metrics were found.")

  metric_columns <- c(
    "mean_log_score", "brier", "auc", "accuracy",
    "calibration_intercept", "calibration_slope", "ece_10",
    "wall_seconds"
  )
  summary <- mean_se_summary(
    metrics,
    value_columns = metric_columns,
    group_columns = c("genre", "model")
  )

  reference <- metrics[metrics$model == reference_model, , drop = FALSE]
  comparisons <- metrics[metrics$model != reference_model, , drop = FALSE]
  if (nrow(reference) > 0L && nrow(comparisons) > 0L) {
    paired <- merge(
      comparisons,
      reference,
      by = c("fold", "genre"),
      suffixes = c("_model", "_reference")
    )
    differences <- do.call(rbind, lapply(metric_columns, function(metric) {
      value <- paired[[paste0(metric, "_model")]] -
        paired[[paste0(metric, "_reference")]]
      data.frame(
        genre = paired$genre,
        fold = paired$fold,
        model = paired$model_model,
        reference_model = reference_model,
        metric = metric,
        difference_model_minus_reference = value
      )
    }))
    difference_summary <- mean_se_summary(
      differences,
      value_columns = "difference_model_minus_reference",
      group_columns = c("genre", "model", "reference_model", "metric")
    )
  } else {
    paired <- data.frame()
    differences <- data.frame(
      genre = character(),
      fold = integer(),
      model = character(),
      reference_model = character(),
      metric = character(),
      difference_model_minus_reference = numeric()
    )
    difference_summary <- data.frame()
  }

  write_csv_atomic(metrics, file.path(output_dir, "all_fold_metrics.csv"))
  write_csv_atomic(summary, file.path(output_dir, "model_summary.csv"))
  write_csv_atomic(differences, file.path(output_dir, "paired_differences.csv"))
  if (nrow(difference_summary) > 0L) {
    write_csv_atomic(
      difference_summary,
      file.path(output_dir, "paired_difference_summary.csv")
    )
  }
  if (nrow(diagnostics) > 0L) {
    write_csv_atomic(
      diagnostics,
      file.path(output_dir, "all_mcmc_diagnostics.csv")
    )
  }
  spatial_signal <- data.frame()
  gain_signal <- data.frame()
  if (nrow(posterior) > 0L) {
    write_csv_atomic(
      posterior,
      file.path(output_dir, "all_posterior_summaries.csv")
    )
    selected <- posterior[
      posterior$model == reference_model &
        posterior$parameter %in% c(
          "sigma_f", "sigma_u", "nugget_ratio",
          "ell", "anisotropy_ratio"
        ),
      ,
      drop = FALSE
    ]
    if (nrow(selected) > 0L) {
      # Older posterior_summaries.csv files may contain NA indices because the
      # former parser left the closing bracket in names such as ell[1].  Repair
      # those indices from the preserved variable name before reshaping.
      ell_rows <- selected$parameter == "ell"
      missing_ell_index <- ell_rows & is.na(selected$index)
      if (any(missing_ell_index) && "variable" %in% names(selected)) {
        recovered <- sub("^ell\\[([0-9]+)\\]$", "\\1", selected$variable[missing_ell_index])
        recovered[!grepl("^[0-9]+$", recovered)] <- NA_character_
        selected$index[missing_ell_index] <- suppressWarnings(as.integer(recovered))
      }
      selected$parameter_label <- ifelse(
        ell_rows & !is.na(selected$index),
        paste0("ell_", selected$index),
        selected$parameter
      )
      # There should be one value per fold, genre, and parameter label.  This
      # aggregation is defensive against stale duplicated summary rows.
      selected <- stats::aggregate(
        mean ~ fold + genre + parameter_label,
        data = selected,
        FUN = mean
      )
      selected$spatial_value <- selected$mean
      spatial_signal <- stats::reshape(
        selected[, c("fold", "genre", "parameter_label", "spatial_value")],
        idvar = c("fold", "genre"),
        timevar = "parameter_label",
        direction = "wide"
      )
      names(spatial_signal) <- sub("^spatial_value\\.", "", names(spatial_signal))
      if ("nugget_ratio" %in% names(spatial_signal)) {
        spatial_signal$spatial_fraction <- 1 - spatial_signal$nugget_ratio
      }
      log_score_gain <- differences[
        differences$metric == "mean_log_score",
        ,
        drop = FALSE
      ]
      # Positive values mean SGP-IRT has the higher held-out log score.
      log_score_gain$sgp_log_score_gain <-
        -log_score_gain$difference_model_minus_reference
      if (nrow(log_score_gain) > 0L) {
        gain_signal <- merge(
          log_score_gain[, c(
            "fold", "genre", "model", "reference_model", "sgp_log_score_gain"
          )],
          spatial_signal,
          by = c("fold", "genre"),
          all.x = TRUE
        )
      }
      write_csv_atomic(
        spatial_signal,
        file.path(output_dir, "spatial_signal_by_fold.csv")
      )
      if (nrow(gain_signal) > 0L) {
        write_csv_atomic(
          gain_signal,
          file.path(output_dir, "predictive_gain_vs_spatial_signal.csv")
        )
      }
    }
  }
  invisible(list(
    metrics = metrics,
    summary = summary,
    paired = differences,
    paired_summary = difference_summary,
    diagnostics = diagnostics,
    posterior = posterior,
    spatial_signal = spatial_signal,
    gain_signal = gain_signal
  ))
}
