# Add matched 1PLUS, 2PLUS, and 3PLUS models to the corrected ART analysis.
#
# Run this script from the top-level art_revision directory, after placing
# psy.csv in data/psy.csv. Existing completed primary-CV and descriptive fits
# are reused because art_config$skip_completed is TRUE by default.
#
# Reporting labels:
#   car_1pl = matched 1PLUS  (fixed item discrimination, no guessing)
#   car_2pl = matched 2PLUS  (item discrimination, no guessing)
#   car_3pl = matched 3PLUS  (item discrimination and item guessing)
#
# All three PLUS models use the same honest respondent splits, residual-MDS
# coordinates, continuous standardized vocabulary covariate, proper CAR plus
# nugget difficulty prior, identification, common priors, MCMC settings, and
# held-out scoring code as the other ART competitors.

source(file.path("config", "art_config.R"))
source(file.path("R", "00_utils.R"))
source_art_files(".")

# Set either flag to FALSE if that part of the analysis is not needed.
RUN_HONEST_CV <- TRUE
RUN_HONEST_DESCRIPTIVE <- TRUE

comparison_models <- c(
  "ind_2pl",
  "car_1pl",
  "car_2pl",
  "car_3pl",
  "sgp_irt"
)

model_labels <- c(
  ind_2pl = "Independent 2PL",
  car_1pl = "1PLUS",
  car_2pl = "2PLUS",
  car_3pl = "3PLUS",
  sgp_irt = "SGP-IRT"
)

add_model_labels <- function(data) {
  if (nrow(data) == 0L) {
    return(data)
  }
  if ("model" %in% names(data)) {
    data$model_label <- unname(model_labels[as.character(data$model)])
    data$model_order <- match(as.character(data$model), comparison_models)
  }
  if ("reference_model" %in% names(data)) {
    data$reference_model_label <- unname(
      model_labels[as.character(data$reference_model)]
    )
  }
  data
}

read_if_present <- function(path) {
  if (!file.exists(path)) {
    return(data.frame())
  }
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

write_plus_cv_outputs <- function(config) {
  summary_dir <- file.path(config$results_dir, "primary_cv_summary")
  assert_true(
    dir.exists(summary_dir),
    "The primary CV summary directory was not created."
  )

  summary <- add_model_labels(read_if_present(
    file.path(summary_dir, "model_summary.csv")
  ))
  fold_metrics <- add_model_labels(read_if_present(
    file.path(summary_dir, "all_fold_metrics.csv")
  ))
  paired <- add_model_labels(read_if_present(
    file.path(summary_dir, "paired_difference_summary.csv")
  ))
  diagnostics <- add_model_labels(read_if_present(
    file.path(summary_dir, "all_mcmc_diagnostics.csv")
  ))
  posterior <- add_model_labels(read_if_present(
    file.path(summary_dir, "all_posterior_summaries.csv")
  ))

  keep_summary <- summary$model %in% comparison_models
  summary <- summary[keep_summary, , drop = FALSE]
  summary <- summary[order(summary$genre, summary$model_order), , drop = FALSE]
  write_csv_atomic(
    summary,
    file.path(summary_dir, "ART_plus_comparison_model_summary.csv")
  )

  keep_fold <- fold_metrics$model %in% comparison_models
  fold_metrics <- fold_metrics[keep_fold, , drop = FALSE]
  fold_metrics <- fold_metrics[
    order(fold_metrics$genre, fold_metrics$fold, fold_metrics$model_order),
    ,
    drop = FALSE
  ]
  write_csv_atomic(
    fold_metrics,
    file.path(summary_dir, "ART_plus_comparison_fold_metrics.csv")
  )

  if (nrow(paired) > 0L) {
    paired <- paired[paired$model %in% comparison_models, , drop = FALSE]
    paired <- paired[
      order(paired$genre, paired$metric, paired$model_order),
      ,
      drop = FALSE
    ]
    write_csv_atomic(
      paired,
      file.path(summary_dir, "ART_plus_comparison_paired_differences.csv")
    )
  }

  if (nrow(diagnostics) > 0L) {
    diagnostics <- diagnostics[
      diagnostics$model %in% comparison_models,
      ,
      drop = FALSE
    ]
    diagnostics$diagnostic_pass <- with(
      diagnostics,
      is.finite(max_rhat) & max_rhat < 1.01 &
        is.finite(min_ess_bulk) & min_ess_bulk >= 400 &
        is.finite(min_ess_tail) & min_ess_tail >= 400 &
        divergences == 0 &
        treedepth_hits == 0 &
        is.finite(min_ebfmi) & min_ebfmi > 0.30
    )
    diagnostics <- diagnostics[
      order(diagnostics$genre, diagnostics$fold, diagnostics$model_order),
      ,
      drop = FALSE
    ]
    write_csv_atomic(
      diagnostics,
      file.path(summary_dir, "ART_plus_comparison_diagnostics.csv")
    )
  }

  if (nrow(posterior) > 0L) {
    guessing <- posterior[
      posterior$model == "car_3pl" & posterior$parameter == "guessing",
      ,
      drop = FALSE
    ]
    guessing <- guessing[order(guessing$genre, guessing$fold, guessing$index), ]
    write_csv_atomic(
      guessing,
      file.path(summary_dir, "ART_3PLUS_guessing_summaries.csv")
    )
  }

  invisible(list(
    model_summary = summary,
    fold_metrics = fold_metrics,
    paired_differences = paired,
    diagnostics = diagnostics,
    posterior = posterior
  ))
}

write_plus_descriptive_outputs <- function(config) {
  output_dir <- file.path(config$results_dir, "honest_descriptive")
  diagnostics <- add_model_labels(read_if_present(
    file.path(output_dir, "all_mcmc_diagnostics.csv")
  ))
  posterior <- add_model_labels(read_if_present(
    file.path(output_dir, "all_posterior_summaries.csv")
  ))

  if (nrow(diagnostics) > 0L) {
    diagnostics <- diagnostics[
      diagnostics$model %in% comparison_models,
      ,
      drop = FALSE
    ]
    diagnostics$diagnostic_pass <- with(
      diagnostics,
      is.finite(max_rhat) & max_rhat < 1.01 &
        is.finite(min_ess_bulk) & min_ess_bulk >= 400 &
        is.finite(min_ess_tail) & min_ess_tail >= 400 &
        divergences == 0 &
        treedepth_hits == 0 &
        is.finite(min_ebfmi) & min_ebfmi > 0.30
    )
    diagnostics <- diagnostics[
      order(diagnostics$genre, diagnostics$model_order),
      ,
      drop = FALSE
    ]
    write_csv_atomic(
      diagnostics,
      file.path(output_dir, "ART_plus_descriptive_diagnostics.csv")
    )
  }

  if (nrow(posterior) > 0L) {
    posterior <- posterior[
      posterior$model %in% comparison_models,
      ,
      drop = FALSE
    ]
    posterior <- posterior[
      order(posterior$genre, posterior$model_order, posterior$parameter,
            posterior$index),
      ,
      drop = FALSE
    ]
    write_csv_atomic(
      posterior,
      file.path(output_dir, "ART_plus_descriptive_posterior_summaries.csv")
    )
    write_csv_atomic(
      posterior[
        posterior$model == "car_3pl" & posterior$parameter == "guessing",
        ,
        drop = FALSE
      ],
      file.path(output_dir, "ART_3PLUS_descriptive_guessing_summaries.csv")
    )
  }

  invisible(list(diagnostics = diagnostics, posterior = posterior))
}

if (isTRUE(RUN_HONEST_CV)) {
  # This uses the same run label as the completed primary analysis. With
  # skip_completed = TRUE, existing independent-2PL, 2PLUS, and SGP-IRT
  # fits are loaded/skipped; only missing 1PLUS and 3PLUS fits are sampled.
  run_art_cv(
    config = art_config,
    run_label = "primary_cv",
    models = comparison_models,
    genres = names(art_config$genres),
    folds_to_run = seq_len(art_config$n_folds),
    geometry_method = "residual",
    anisotropic = TRUE
  )
  write_plus_cv_outputs(art_config)
}

if (isTRUE(RUN_HONEST_DESCRIPTIVE)) {
  # Geometry and calibration respondents remain disjoint. Existing model fits
  # in results/honest_descriptive are reused when present.
  run_honest_descriptive_analysis(
    config = art_config,
    models = comparison_models,
    geometry_method = "residual",
    anisotropic = TRUE
  )
  write_plus_descriptive_outputs(art_config)
}

message(
  "PLUS real-data comparison complete. See results/primary_cv_summary/",
  "ART_plus_comparison_model_summary.csv and the associated diagnostics."
)
