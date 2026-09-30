# The primary curves are conditional on vocabulary_z = 0, corresponding to a
# respondent at the calibration-sample mean vocabulary score. For posterior
# draw s and item j,
#
#   p_j^(s)(theta, x) = logit^{-1}{alpha_j^(s) theta
#                                  - beta_j^(s) + gamma^(s) x},
#
#   average endorsement = J_g^{-1} sum_j p_j^(s)(theta, x),
#   test information     = sum_j [alpha_j^(s)]^2 p_j^(s)(theta, x)
#                                              {1-p_j^(s)(theta, x)}.

source(file.path("config", "art_config.R"))
source(file.path("R", "00_utils.R"))
source_art_files(".")

required_packages <- c("cmdstanr", "posterior", "ggplot2", "patchwork")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Install the required R packages first: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

# ---------------------------------------------------------------------------
# User-adjustable post-processing settings
# ---------------------------------------------------------------------------

THETA_GRID <- seq(-3, 3, length.out = 301L)
VOCABULARY_Z <- 0
MAX_POSTERIOR_DRAWS <- 2000L
POSTERIOR_PROBABILITY <- 0.95
INFORMATION_THRESHOLD <- 2
CURVE_SEED <- art_config$split_seed + 905L

genre_order <- c(
  "SciFiFantasy",
  "SocialPoliticalCommentary",
  "Suspense"
)

genre_labels <- c(
  SciFiFantasy = "Sci-Fi/Fantasy",
  SocialPoliticalCommentary = "Social/Political Commentary",
  Suspense = "Suspense"
)

genre_colors <- c(
  SciFiFantasy = "#1B9E77",
  SocialPoliticalCommentary = "#D95F02",
  Suspense = "#7570B3"
)

fit_root <- file.path(art_config$results_dir, "honest_descriptive")
output_dir <- ensure_dir(file.path(fit_root, "genre_operating_characteristics"))

fit_paths <- setNames(
  file.path(fit_root, genre_order, "sgp_irt", "fit.rds"),
  genre_order
)

missing_fits <- fit_paths[!file.exists(fit_paths)]
assert_true(
  length(missing_fits) == 0L,
  paste0(
    "Missing saved SGP-IRT fit(s): ",
    paste(missing_fits, collapse = ", "),
    ". Run scripts/02_run_honest_descriptive.R first."
  )
)

interval_matrix <- function(draw_matrix, probability) {
  t(apply(
    draw_matrix,
    2L,
    shortest_interval,
    probability = probability
  ))
}

summarize_scalar_draws <- function(x, prefix, probability) {
  summary <- summarize_vector(x, probability = probability)
  result <- unlist(summary[1L, c("mean", "sd", "median", "lower", "upper")])
  names(result) <- paste0(prefix, "_", names(result))
  result
}

compute_sgp_genre_operating_curves <- function(
    fit,
    genre,
    theta_grid,
    vocabulary_z,
    max_draws,
    probability,
    information_threshold,
    seed) {
  alpha_all <- indexed_draws(fit, "alpha")
  J <- ncol(alpha_all)
  beta_all <- indexed_draws(fit, "beta", J)
  gamma_all <- as.numeric(indexed_draws(fit, "gamma", 1L)[, 1L])

  assert_true(J >= 1L, paste0(genre, ": no item parameters were found."))
  assert_true(
    nrow(alpha_all) == nrow(beta_all) &&
      nrow(alpha_all) == length(gamma_all),
    paste0(genre, ": alpha, beta, and gamma draw counts do not agree.")
  )
  assert_true(
    length(vocabulary_z) == 1L && is.finite(vocabulary_z),
    "VOCABULARY_Z must be one finite standardized vocabulary value."
  )

  rows <- thin_draw_rows(nrow(alpha_all), max_draws, seed)
  alpha <- alpha_all[rows, , drop = FALSE]
  beta <- beta_all[rows, , drop = FALSE]
  gamma <- gamma_all[rows]

  Q <- length(theta_grid)
  S <- length(rows)
  endorsement_draws <- matrix(NA_real_, nrow = S, ncol = Q)
  information_draws <- matrix(NA_real_, nrow = S, ncol = Q)

  for (s in seq_len(S)) {
    # Explicit numeric coercion avoids the one-row draws_matrix dimension bug
    # encountered in earlier versions of the ART post-processing code.
    alpha_s <- as.numeric(alpha[s, , drop = TRUE])
    beta_s <- as.numeric(beta[s, , drop = TRUE])
    gamma_s <- as.numeric(gamma[s])

    assert_true(
      length(alpha_s) == J && length(beta_s) == J,
      paste0(genre, ": posterior item-vector dimensions are inconsistent.")
    )

    eta <- outer(as.numeric(theta_grid), alpha_s) -
      matrix(beta_s, nrow = Q, ncol = J, byrow = TRUE) +
      vocabulary_z * gamma_s
    p <- stats::plogis(eta)

    endorsement_draws[s, ] <- rowMeans(p)
    information_draws[s, ] <- rowSums(
      sweep(p * (1 - p), 2L, alpha_s^2, "*")
    )
  }

  endorsement_intervals <- interval_matrix(
    endorsement_draws,
    probability = probability
  )
  information_intervals <- interval_matrix(
    information_draws,
    probability = probability
  )

  curve_data <- rbind(
    data.frame(
      genre = genre,
      genre_label = unname(genre_labels[genre]),
      theta = theta_grid,
      curve = "Average endorsement probability",
      posterior_mean = colMeans(endorsement_draws),
      posterior_median = apply(endorsement_draws, 2L, stats::median),
      lower = endorsement_intervals[, "lower"],
      upper = endorsement_intervals[, "upper"],
      vocabulary_z = vocabulary_z,
      n_items = J,
      n_posterior_draws = S
    ),
    data.frame(
      genre = genre,
      genre_label = unname(genre_labels[genre]),
      theta = theta_grid,
      curve = "Test information",
      posterior_mean = colMeans(information_draws),
      posterior_median = apply(information_draws, 2L, stats::median),
      lower = information_intervals[, "lower"],
      upper = information_intervals[, "upper"],
      vocabulary_z = vocabulary_z,
      n_items = J,
      n_posterior_draws = S
    )
  )

  # Draw-wise summaries propagate posterior uncertainty in alpha, beta, and
  # gamma. The peak of the posterior mean curve is recorded separately because
  # it is the quantity most naturally described in the manuscript.
  peak_index_by_draw <- max.col(information_draws, ties.method = "first")
  peak_information_by_draw <- information_draws[
    cbind(seq_len(S), peak_index_by_draw)
  ]
  peak_theta_by_draw <- theta_grid[peak_index_by_draw]
  peak_se_by_draw <- 1 / sqrt(pmax(peak_information_by_draw, 1e-12))

  theta_zero_index <- which.min(abs(theta_grid))
  theta_50_index_by_draw <- apply(
    endorsement_draws,
    1L,
    function(x) which.min(abs(x - 0.5))
  )
  theta_50_by_draw <- theta_grid[theta_50_index_by_draw]

  mean_information_curve <- colMeans(information_draws)
  mean_endorsement_curve <- colMeans(endorsement_draws)
  mean_peak_index <- which.max(mean_information_curve)
  threshold_index <- which(mean_information_curve >= information_threshold)

  summary_values <- c(
    genre = genre,
    genre_label = unname(genre_labels[genre]),
    n_items = J,
    vocabulary_z = vocabulary_z,
    n_posterior_draws = S,
    theta_grid_min = min(theta_grid),
    theta_grid_max = max(theta_grid),
    posterior_mean_curve_peak_theta = theta_grid[mean_peak_index],
    posterior_mean_curve_peak_information = mean_information_curve[mean_peak_index],
    posterior_mean_curve_peak_se = 1 / sqrt(mean_information_curve[mean_peak_index]),
    posterior_mean_endorsement_at_theta_0 = mean_endorsement_curve[theta_zero_index],
    posterior_mean_information_at_theta_0 = mean_information_curve[theta_zero_index],
    posterior_mean_info_threshold = information_threshold,
    posterior_mean_info_threshold_lower_theta = if (length(threshold_index)) {
      min(theta_grid[threshold_index])
    } else {
      NA_real_
    },
    posterior_mean_info_threshold_upper_theta = if (length(threshold_index)) {
      max(theta_grid[threshold_index])
    } else {
      NA_real_
    },
    summarize_scalar_draws(
      peak_information_by_draw,
      "drawwise_peak_information",
      probability
    ),
    summarize_scalar_draws(
      peak_theta_by_draw,
      "drawwise_peak_theta",
      probability
    ),
    summarize_scalar_draws(
      peak_se_by_draw,
      "drawwise_peak_se",
      probability
    ),
    summarize_scalar_draws(
      endorsement_draws[, theta_zero_index],
      "endorsement_at_theta_0",
      probability
    ),
    summarize_scalar_draws(
      information_draws[, theta_zero_index],
      "information_at_theta_0",
      probability
    ),
    summarize_scalar_draws(
      theta_50_by_draw,
      "theta_at_average_endorsement_0_5",
      probability
    )
  )

  summary_data <- as.data.frame(
    as.list(summary_values),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  numeric_columns <- setdiff(names(summary_data), c("genre", "genre_label"))
  summary_data[numeric_columns] <- lapply(
    summary_data[numeric_columns],
    function(x) as.numeric(as.character(x))
  )

  list(curves = curve_data, summary = summary_data)
}

genre_results <- lapply(seq_along(genre_order), function(k) {
  genre <- genre_order[k]
  message("Reading saved SGP-IRT posterior: ", genre)
  fit <- readRDS(fit_paths[[genre]])
  compute_sgp_genre_operating_curves(
    fit = fit,
    genre = genre,
    theta_grid = THETA_GRID,
    vocabulary_z = VOCABULARY_Z,
    max_draws = MAX_POSTERIOR_DRAWS,
    probability = POSTERIOR_PROBABILITY,
    information_threshold = INFORMATION_THRESHOLD,
    seed = CURVE_SEED + k
  )
})

curve_data <- do.call(rbind, lapply(genre_results, `[[`, "curves"))
summary_data <- do.call(rbind, lapply(genre_results, `[[`, "summary"))

curve_data$genre <- factor(curve_data$genre, levels = genre_order)
curve_data$genre_label <- factor(
  curve_data$genre_label,
  levels = unname(genre_labels[genre_order])
)

write_csv_atomic(
  curve_data,
  file.path(output_dir, "ART_SGP_IRT_genre_ICC_TIC_curves.csv")
)
write_csv_atomic(
  summary_data,
  file.path(output_dir, "ART_SGP_IRT_genre_ICC_TIC_summary.csv")
)

endorsement_data <- curve_data[
  curve_data$curve == "Average endorsement probability",
  ,
  drop = FALSE
]
information_data <- curve_data[
  curve_data$curve == "Test information",
  ,
  drop = FALSE
]

common_theme <- ggplot2::theme_minimal(base_size = 13) +
  ggplot2::theme(
    legend.position = "bottom",
    legend.title = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(7, 10, 7, 10),
    axis.title = ggplot2::element_text(size = 13),
    axis.text = ggplot2::element_text(size = 11)
  )

p_endorsement <- ggplot2::ggplot(
  endorsement_data,
  ggplot2::aes(
    x = theta,
    y = posterior_mean,
    color = genre,
    fill = genre
  )
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = lower, ymax = upper),
    alpha = 0.12,
    color = NA
  ) +
  ggplot2::geom_line(linewidth = 1.0) +
  ggplot2::scale_color_manual(
    values = genre_colors,
    breaks = genre_order,
    labels = unname(genre_labels[genre_order])
  ) +
  ggplot2::scale_fill_manual(
    values = genre_colors,
    breaks = genre_order,
    labels = unname(genre_labels[genre_order])
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25),
    expand = ggplot2::expansion(mult = c(0.01, 0.02))
  ) +
  ggplot2::labs(
    x = expression(theta),
    y = "Average endorsement probability"
  ) +
  common_theme

p_information <- ggplot2::ggplot(
  information_data,
  ggplot2::aes(
    x = theta,
    y = posterior_mean,
    color = genre,
    fill = genre
  )
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = lower, ymax = upper),
    alpha = 0.12,
    color = NA
  ) +
  ggplot2::geom_line(linewidth = 1.0) +
  ggplot2::scale_color_manual(
    values = genre_colors,
    breaks = genre_order,
    labels = unname(genre_labels[genre_order])
  ) +
  ggplot2::scale_fill_manual(
    values = genre_colors,
    breaks = genre_order,
    labels = unname(genre_labels[genre_order])
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0.02, 0.06))
  ) +
  ggplot2::labs(
    x = expression(theta),
    y = "Test information"
  ) +
  common_theme

combined_plot <- (
  p_endorsement + p_information +
    patchwork::plot_layout(guides = "collect", widths = c(1, 1))
) &
  ggplot2::theme(legend.position = "bottom")

pdf_path <- file.path(output_dir, "ART_SGP_IRT_genre_ICC_TIC.pdf")
png_path <- file.path(output_dir, "ART_SGP_IRT_genre_ICC_TIC.png")

ggplot2::ggsave(
  pdf_path,
  combined_plot,
  width = 12,
  height = 5.8,
  units = "in",
  device = grDevices::cairo_pdf
)
ggplot2::ggsave(
  png_path,
  combined_plot,
  width = 12,
  height = 5.8,
  units = "in",
  dpi = 400,
  bg = "white"
)

message("SGP-IRT genre operating-characteristic outputs written to: ", output_dir)
message("  ", pdf_path)
message("  ", png_path)
message("  ", file.path(output_dir, "ART_SGP_IRT_genre_ICC_TIC_curves.csv"))
message("  ", file.path(output_dir, "ART_SGP_IRT_genre_ICC_TIC_summary.csv"))
message("No Stan model was compiled and no MCMC chain was run.")

