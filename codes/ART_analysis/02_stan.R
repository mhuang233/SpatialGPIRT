art_model_files <- function(stan_dir = "stan") {
  c(
    sgp_irt = file.path(stan_dir, "sgp_irt_2pl.stan"),
    car_2pl = file.path(stan_dir, "car_irt_2pl.stan"),
    ind_2pl = file.path(stan_dir, "ind_irt_2pl.stan"),
    car_1pl = file.path(stan_dir, "car_irt_1pl.stan"),
    car_3pl = file.path(stan_dir, "car_irt_3pl.stan")
  )
}

check_art_dependencies <- function() {
  packages <- c(
    "cmdstanr", "posterior", "ggplot2", "patchwork",
    "dplyr", "tidyr", "readr"
  )
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    stop(
      "Install the required R packages first: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  version <- tryCatch(cmdstanr::cmdstan_version(), error = function(e) NULL)
  assert_true(!is.null(version), "CmdStan is not installed or configured for CmdStanR.")
  invisible(version)
}

compile_art_models <- function(
    stan_dir = "stan",
    models = names(art_model_files(stan_dir)),
    force_recompile = FALSE,
    quiet = FALSE) {
  check_art_dependencies()
  files <- art_model_files(stan_dir)
  assert_true(all(models %in% names(files)), "Unknown model requested for compilation.")
  selected <- files[models]
  missing <- selected[!file.exists(selected)]
  assert_true(length(missing) == 0L, paste0(
    "Missing Stan files: ", paste(missing, collapse = ", ")
  ))

  compiled <- lapply(selected, function(path) {
    cmdstanr::cmdstan_model(
      stan_file = path,
      cpp_options = list(stan_threads = TRUE),
      force_recompile = force_recompile,
      quiet = quiet
    )
  })
  names(compiled) <- names(selected)
  compiled
}

long_irt_data <- function(Y, vocabulary_z, grainsize, sigma_gamma_prior) {
  Y <- as.matrix(Y)
  I <- nrow(Y)
  J <- ncol(Y)
  assert_true(length(vocabulary_z) == I, "Vocabulary length does not match respondents.")
  assert_true(all(Y %in% c(0L, 1L)), "Response matrix is not binary.")
  list(
    N = I * J,
    I = I,
    J = J,
    y = as.integer(as.vector(t(Y))),
    subj_id = rep(seq_len(I), each = J),
    item_id = rep(seq_len(J), times = I),
    p = 1L,
    X = matrix(rep(vocabulary_z, each = J), ncol = 1L),
    grainsize = as.integer(grainsize),
    sigma_gamma_prior = sigma_gamma_prior
  )
}

squared_distance_matrices <- function(coords) {
  coords <- as.matrix(coords)
  assert_true(ncol(coords) == 2L, "The current ART Stan model requires two MDS coordinates.")
  lapply(seq_len(2L), function(d) {
    delta <- outer(coords[, d], coords[, d], "-")
    delta^2
  })
}

add_sgp_data <- function(base_data, coords, anisotropic, a_A, b_A) {
  distance <- squared_distance_matrices(coords)
  c(
    base_data,
    list(
      sqdist1 = distance[[1L]],
      sqdist2 = distance[[2L]],
      n_A = if (isTRUE(anisotropic)) 2L else 1L,
      a_A = a_A,
      b_A = b_A
    )
  )
}

add_car_data <- function(base_data, graph) {
  c(
    base_data,
    graph[c(
      "degree",
      "N_edges",
      "node1",
      "node2",
      "normalized_adjacency_eigenvalues",
      "log_degree_sum"
    )]
  )
}

model_initial_values <- function(model_name, I, J, n_A = 2L, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  common <- list(
    z_theta = stats::rnorm(I, 0, 0.25),
    gamma = rep(0, 1L)
  )
  switch(
    model_name,
    ind_2pl = c(common, list(
      z_log_alpha = rep(0, J),
      sigma_a = 0.25,
      z_beta = stats::rnorm(J, 0, 0.05),
      sigma_beta = 0.75
    )),
    sgp_irt = c(common, list(
      z_log_alpha = rep(0, J),
      sigma_a = 0.25,
      z_f = stats::rnorm(J, 0, 0.05),
      z_u = stats::rnorm(J, 0, 0.05),
      sigma_f = 0.50,
      sigma_u = 0.25,
      A_pow = rep(4, n_A)
    )),
    car_2pl = c(common, list(
      z_log_alpha = rep(0, J),
      sigma_a = 0.25,
      f = rep(0, J),
      z_u = stats::rnorm(J, 0, 0.05),
      sigma_f = 0.50,
      sigma_u = 0.25,
      rho = 0.50
    )),
    car_1pl = c(common, list(
      f = rep(0, J),
      z_u = stats::rnorm(J, 0, 0.05),
      sigma_f = 0.50,
      sigma_u = 0.25,
      rho = 0.50
    )),
    car_3pl = c(common, list(
      z_log_alpha = rep(0, J),
      sigma_a = 0.25,
      f = rep(0, J),
      z_u = stats::rnorm(J, 0, 0.05),
      sigma_f = 0.50,
      sigma_u = 0.25,
      rho = 0.50,
      guessing = rep(0.10, J)
    )),
    stop("No initializer for model: ", model_name, call. = FALSE)
  )
}

sample_art_model <- function(
    model,
    model_name,
    stan_data,
    output_dir,
    seed,
    config,
    init = NULL) {
  ensure_dir(output_dir)
  fit_path <- file.path(output_dir, "fit.rds")
  if (isTRUE(config$skip_completed) && file.exists(fit_path)) {
    return(readRDS(fit_path))
  }

  init <- init %||% function(chain_id = 1L) {
    model_initial_values(
      model_name = model_name,
      I = stan_data$I,
      J = stan_data$J,
      n_A = stan_data$n_A %||% 2L,
      seed = seed + chain_id
    )
  }

  started <- Sys.time()
  csv_dir <- ensure_dir(file.path(output_dir, "cmdstan_csv"))
  fit <- model$sample(
    data = stan_data,
    seed = seed,
    chains = config$chains,
    parallel_chains = min(config$parallel_chains, config$chains),
    threads_per_chain = config$threads_per_chain,
    iter_warmup = config$iter_warmup,
    iter_sampling = config$iter_sampling,
    adapt_delta = config$adapt_delta,
    max_treedepth = config$max_treedepth,
    refresh = config$refresh,
    output_dir = csv_dir,
    init = init
  )
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  writeLines(format(elapsed, digits = 15), file.path(output_dir, "wall_seconds.txt"))
  if (isTRUE(config$save_fit_objects)) {
    fit$save_object(file = fit_path)
  }
  fit
}

draw_matrix <- function(fit, variables = NULL) {
  posterior::as_draws_matrix(fit$draws(variables = variables))
}

indexed_draws <- function(fit, variable, length_expected = NULL) {
  draws <- draw_matrix(fit, variable)
  pattern <- paste0("^", variable, "\\[([0-9]+)\\]$")
  columns <- grep(pattern, colnames(draws), value = TRUE)
  if (length(columns) == 0L && variable %in% colnames(draws)) {
    result <- matrix(draws[, variable], ncol = 1L)
    colnames(result) <- variable
    return(result)
  }
  index <- as.integer(sub(pattern, "\\1", columns))
  columns <- columns[order(index)]
  result <- as.matrix(draws[, columns, drop = FALSE])
  if (!is.null(length_expected)) {
    assert_true(
      ncol(result) == length_expected,
      paste0("Unexpected length for ", variable, ".")
    )
  }
  result
}

thin_draw_rows <- function(n_draws, max_draws, seed) {
  if (n_draws <= max_draws) {
    return(seq_len(n_draws))
  }
  set.seed(seed)
  sort(sample.int(n_draws, max_draws))
}
