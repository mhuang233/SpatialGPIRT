variable_draw_summaries <- function(
    fit,
    variables,
    probability = 0.95) {
  draws <- draw_matrix(fit, variables)
  summaries <- lapply(seq_len(ncol(draws)), function(k) {
    value <- draws[, k]
    interval <- shortest_interval(value, probability)
    variable <- colnames(draws)[k]
    base <- sub("\\[.*$", "", variable)
    index_text <- sub("^.*\\[|\\]$", "", variable)
    index <- if (grepl("\\[", variable)) {
      suppressWarnings(as.integer(index_text))
    } else {
      NA_integer_
    }
    data.frame(
      variable = variable,
      parameter = base,
      index = index,
      mean = mean(value),
      sd = stats::sd(value),
      median = stats::median(value),
      lower = unname(interval["lower"]),
      upper = unname(interval["upper"])
    )
  })
  do.call(rbind, summaries)
}

fit_convergence_summary <- function(
    fit,
    monitored_variables,
    max_treedepth,
    wall_seconds = NA_real_) {
  parameter_summary <- fit$summary(variables = monitored_variables)
  sampler <- posterior::as_draws_array(fit$sampler_diagnostics())
  sampler_names <- dimnames(sampler)[[3L]]
  divergences <- if ("divergent__" %in% sampler_names) {
    sum(sampler[, , "divergent__"])
  } else {
    NA_real_
  }
  treedepth_hits <- if ("treedepth__" %in% sampler_names) {
    sum(sampler[, , "treedepth__"] >= max_treedepth)
  } else {
    NA_real_
  }
  bfmi <- if ("energy__" %in% sampler_names) {
    vapply(seq_len(dim(sampler)[2L]), function(chain) {
      energy <- sampler[, chain, "energy__"]
      mean(diff(energy)^2) / stats::var(energy)
    }, numeric(1))
  } else {
    NA_real_
  }

  data.frame(
    max_rhat = max(parameter_summary$rhat, na.rm = TRUE),
    rhat_99_percentile = stats::quantile(
      parameter_summary$rhat,
      0.99,
      na.rm = TRUE,
      names = FALSE
    ),
    min_ess_bulk = min(parameter_summary$ess_bulk, na.rm = TRUE),
    min_ess_tail = min(parameter_summary$ess_tail, na.rm = TRUE),
    divergences = divergences,
    treedepth_hits = treedepth_hits,
    min_ebfmi = min(bfmi, na.rm = TRUE),
    median_ebfmi = stats::median(bfmi, na.rm = TRUE),
    wall_seconds = wall_seconds
  )
}

default_monitored_variables <- function(model_name) {
  switch(
    model_name,
    sgp_irt = c(
      "theta", "alpha", "beta", "f", "u", "gamma", "sigma_a",
      "sigma_f", "sigma_u", "ell", "nugget_ratio", "anisotropy_ratio"
    ),
    car_2pl = c(
      "theta", "alpha", "beta", "f", "u", "gamma", "sigma_a",
      "sigma_f", "sigma_u", "rho", "nugget_ratio"
    ),
    ind_2pl = c(
      "theta", "alpha", "beta", "gamma", "sigma_a", "sigma_beta"
    ),
    car_1pl = c(
      "theta", "beta", "f", "u", "gamma", "sigma_f", "sigma_u",
      "rho", "nugget_ratio"
    ),
    car_3pl = c(
      "theta", "alpha", "beta", "f", "u", "guessing", "gamma",
      "sigma_a", "sigma_f", "sigma_u", "rho", "nugget_ratio"
    ),
    stop("No monitored-variable list for ", model_name, call. = FALSE)
  )
}

save_fit_diagnostics <- function(
    fit,
    model_name,
    output_dir,
    max_treedepth,
    item_names = NULL,
    subject_source_rows = NULL) {
  ensure_dir(output_dir)
  wall_path <- file.path(output_dir, "wall_seconds.txt")
  wall_seconds <- if (file.exists(wall_path)) {
    suppressWarnings(as.numeric(readLines(wall_path, n = 1L)))
  } else {
    NA_real_
  }
  variables <- default_monitored_variables(model_name)
  convergence <- fit_convergence_summary(
    fit = fit,
    monitored_variables = variables,
    max_treedepth = max_treedepth,
    wall_seconds = wall_seconds
  )
  posterior <- variable_draw_summaries(fit, variables)

  if (!is.null(item_names)) {
    item_parameters <- posterior$parameter %in% c(
      "alpha", "beta", "f", "u", "guessing"
    )
    valid_index <- item_parameters & !is.na(posterior$index)
    posterior$item <- NA_character_
    posterior$item[valid_index] <- item_names[posterior$index[valid_index]]
  }
  if (!is.null(subject_source_rows)) {
    subject_parameter <- posterior$parameter == "theta" & !is.na(posterior$index)
    posterior$subject_source_row <- NA_integer_
    posterior$subject_source_row[subject_parameter] <-
      subject_source_rows[posterior$index[subject_parameter]]
  }
  write_csv_atomic(convergence, file.path(output_dir, "mcmc_diagnostics.csv"))
  write_csv_atomic(posterior, file.path(output_dir, "posterior_summaries.csv"))
  list(convergence = convergence, posterior = posterior)
}

diagnostic_pass <- function(diagnostics) {
  with(
    diagnostics,
    max_rhat < 1.01 &&
      divergences == 0 &&
      treedepth_hits == 0 &&
      min_ebfmi > 0.30
  )
}
