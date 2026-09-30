# Core data-generating, fitting, diagnostic, and evaluation functions for the
# reviewer-facing SGP-IRT simulation study.

required_simulation_packages <- function() {
  pkgs <- c("cmdstanr", "posterior", "MASS", "deldir")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    stop(
      "Install the following packages before running the simulations: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

scenario_grid <- function() {
  data.frame(
    scenario = c(
      "gp_aniso_low_nugget",
      "gp_aniso_high_nugget",
      "car_truth",
      "matern_rough",
      "piecewise",
      "iid_no_spatial"
    ),
    sigma_f = c(1.00, 0.65, 0.35, 1.00, 0.60, 0.00),
    sigma_u = c(0.15, 0.90, 0.15, 0.15, 0.10, 1.00),
    ell1 = c(0.15, 0.15, NA, 0.15, 0.18, NA),
    ell2 = c(0.45, 0.45, NA, 0.35, 0.35, NA),
    nu = c(NA, NA, NA, 0.50, NA, NA),
    rho = c(NA, NA, 0.95, NA, NA, NA),
    jump_size = c(NA, NA, NA, NA, 1.00, NA),
    stringsAsFactors = FALSE
  )
}

get_scenario <- function(scenario) {
  grid <- scenario_grid()
  ans <- grid[grid$scenario == scenario, , drop = FALSE]
  if (nrow(ans) != 1L) {
    stop("Unknown scenario: ", scenario, call. = FALSE)
  }
  ans
}

standardize_ability <- function(z) {
  centered <- z - mean(z)
  scale <- sqrt(mean(centered^2))
  if (!is.finite(scale) || scale <= 0) {
    stop("Cannot standardize a degenerate ability draw.", call. = FALSE)
  }
  centered / scale
}

generate_quasi_uniform_coords <- function(J) {
  stopifnot(J >= 4L)
  side <- ceiling(sqrt(J))
  candidates <- expand.grid(
    x = (seq_len(side) - 0.5) / side,
    y = (seq_len(side) - 0.5) / side
  )
  keep <- sample(seq_len(nrow(candidates)), J, replace = FALSE)
  coords <- as.matrix(candidates[keep, , drop = FALSE])
  jitter_size <- 0.30 / side
  coords <- coords + matrix(
    runif(2L * J, -jitter_size, jitter_size),
    nrow = J,
    ncol = 2L
  )
  coords[, 1L] <- pmin(1, pmax(0, coords[, 1L]))
  coords[, 2L] <- pmin(1, pmax(0, coords[, 2L]))
  colnames(coords) <- c("s1", "s2")
  coords
}

make_sqdist_mats <- function(coords) {
  stopifnot(is.matrix(coords), ncol(coords) == 2L)
  list(
    sqdist1 = outer(coords[, 1L], coords[, 1L], "-")^2,
    sqdist2 = outer(coords[, 2L], coords[, 2L], "-")^2
  )
}

make_delaunay_graph <- function(coords) {
  stopifnot(is.matrix(coords), ncol(coords) == 2L)
  vor <- deldir::deldir(
    x = coords[, 1L],
    y = coords[, 2L],
    rw = c(range(coords[, 1L]), range(coords[, 2L]))
  )
  edges <- unique(vor$dirsgs[, c("ind1", "ind2"), drop = FALSE])
  J <- nrow(coords)
  W <- matrix(0, nrow = J, ncol = J)
  for (r in seq_len(nrow(edges))) {
    j1 <- edges$ind1[r]
    j2 <- edges$ind2[r]
    W[j1, j2] <- 1
    W[j2, j1] <- 1
  }
  diag(W) <- 0
  if (any(rowSums(W) == 0)) {
    stop("Delaunay graph contains an isolated item.", call. = FALSE)
  }
  W
}

rbf_correlation <- function(coords, ell) {
  d <- make_sqdist_mats(coords)
  R <- exp(
    -0.5 * (
      d$sqdist1 / ell[1L]^2 +
        d$sqdist2 / ell[2L]^2
    )
  )
  diag(R) <- diag(R) + 1e-8
  R
}

matern_correlation <- function(coords, ell, nu) {
  stopifnot(nu > 0)
  d <- make_sqdist_mats(coords)
  r <- sqrt(
    d$sqdist1 / ell[1L]^2 +
      d$sqdist2 / ell[2L]^2
  )
  if (abs(nu - 0.5) < 1e-12) {
    R <- exp(-r)
  } else {
    scaled <- sqrt(2 * nu) * r
    R <- matrix(1, nrow(coords), nrow(coords))
    positive <- r > 0
    R[positive] <-
      2^(1 - nu) / gamma(nu) *
      scaled[positive]^nu *
      besselK(scaled[positive], nu = nu)
  }
  diag(R) <- diag(R) + 1e-8
  R
}

draw_from_covariance <- function(Sigma) {
  U <- chol(Sigma)
  drop(t(U) %*% rnorm(nrow(Sigma)))
}

draw_car_component <- function(W, rho, sigma_f) {
  degree <- rowSums(W)
  Q <- diag(degree) - rho * W
  diag(Q) <- diag(Q) + 1e-8
  U <- chol(Q)
  sigma_f * drop(backsolve(U, rnorm(nrow(W))))
}

stratified_cell_holdout <- function(subj_id, J, test_fraction = 0.20) {
  respondents <- sort(unique(subj_id))
  n_test_each <- max(1L, min(J - 1L, floor(test_fraction * J)))
  test_idx <- unlist(
    lapply(
      respondents,
      function(i) {
        eligible <- which(subj_id == i)
        sample(eligible, n_test_each, replace = FALSE)
      }
    ),
    use.names = FALSE
  )
  sort(test_idx)
}

simulate_irt_dgp <- function(
    scenario,
    I = 200L,
    J = 50L,
    p = 1L,
    sigma_a = 0.35,
    gamma = 0.50,
    test_fraction = 0.20,
    seed = 1L) {
  stopifnot(I >= 10L, J >= 4L, p >= 1L)
  cfg <- get_scenario(scenario)
  set.seed(seed)

  coords <- generate_quasi_uniform_coords(J)
  W <- make_delaunay_graph(coords)
  ell <- c(cfg$ell1, cfg$ell2)

  if (scenario %in% c("gp_aniso_low_nugget", "gp_aniso_high_nugget")) {
    R <- rbf_correlation(coords, ell)
    f_true <- cfg$sigma_f * draw_from_covariance(R)
  } else if (scenario == "car_truth") {
    f_true <- draw_car_component(W, cfg$rho, cfg$sigma_f)
  } else if (scenario == "matern_rough") {
    R <- matern_correlation(coords, ell, cfg$nu)
    f_true <- cfg$sigma_f * draw_from_covariance(R)
  } else if (scenario == "piecewise") {
    R <- rbf_correlation(coords, ell)
    smooth <- cfg$sigma_f * draw_from_covariance(R)
    region <- as.integer(coords[, 1L] > 0.5)
    jump <- cfg$jump_size * (region - mean(region))
    f_true <- smooth + jump
  } else if (scenario == "iid_no_spatial") {
    f_true <- rep(0, J)
  } else {
    stop("Scenario is not implemented: ", scenario, call. = FALSE)
  }

  u_true <- rnorm(J, mean = 0, sd = cfg$sigma_u)
  beta_true <- f_true + u_true

  theta_true <- standardize_ability(rnorm(I))
  log_alpha_true <- sigma_a * rnorm(J)
  alpha_true <- exp(log_alpha_true)

  if (length(gamma) == 1L) {
    gamma <- rep(gamma, p)
  }
  if (length(gamma) != p) {
    stop("gamma must have length one or p.", call. = FALSE)
  }

  N <- I * J
  subj_id <- rep(seq_len(I), each = J)
  item_id <- rep(seq_len(J), times = I)
  X <- matrix(runif(N * p, min = -1, max = 1), nrow = N, ncol = p)
  eta <- alpha_true[item_id] * theta_true[subj_id] -
    beta_true[item_id] +
    drop(X %*% gamma)
  prob <- plogis(eta)
  y <- rbinom(N, size = 1L, prob = prob)

  test_idx <- stratified_cell_holdout(
    subj_id = subj_id,
    J = J,
    test_fraction = test_fraction
  )
  train_idx <- setdiff(seq_len(N), test_idx)
  sqdist <- make_sqdist_mats(coords)

  nugget_ratio <- if (cfg$sigma_f^2 + cfg$sigma_u^2 > 0) {
    cfg$sigma_u^2 / (cfg$sigma_f^2 + cfg$sigma_u^2)
  } else {
    NA_real_
  }

  truth <- list(
    scenario = scenario,
    I = I,
    J = J,
    p = p,
    theta = theta_true,
    log_alpha = log_alpha_true,
    alpha = alpha_true,
    f = f_true,
    u = u_true,
    beta = beta_true,
    gamma = gamma,
    sigma_a = sigma_a,
    sigma_f = cfg$sigma_f,
    sigma_u = cfg$sigma_u,
    sigma_beta = if (scenario == "iid_no_spatial") cfg$sigma_u else NA_real_,
    ell = ell,
    rho = cfg$rho,
    nu = cfg$nu,
    jump_size = cfg$jump_size,
    nugget_ratio = nugget_ratio,
    anisotropy_ratio = if (all(is.finite(ell))) ell[1L] / ell[2L] else NA_real_,
    coords = coords,
    W = W,
    X = X,
    y = y,
    prob = prob,
    subj_id = subj_id,
    item_id = item_id,
    train_idx = train_idx,
    test_idx = test_idx,
    seed = seed
  )

  list(truth = truth, sqdist = sqdist)
}

make_common_stan_data <- function(
    sim,
    grainsize = NULL,
    a_A = 2,
    b_A = 0.1,
    sigma_gamma_prior = 1) {
  tr <- sim$truth
  train <- tr$train_idx
  test <- tr$test_idx
  if (is.null(grainsize)) {
    grainsize <- max(1L, floor(length(train) / 100L))
  }

  list(
    N_train = length(train),
    N_test = length(test),
    I = tr$I,
    J = tr$J,
    p = tr$p,
    y_train = as.integer(tr$y[train]),
    subj_train = as.integer(tr$subj_id[train]),
    item_train = as.integer(tr$item_id[train]),
    X_train = tr$X[train, , drop = FALSE],
    y_test = as.integer(tr$y[test]),
    subj_test = as.integer(tr$subj_id[test]),
    item_test = as.integer(tr$item_id[test]),
    X_test = tr$X[test, , drop = FALSE],
    grainsize = as.integer(grainsize),
    a_A = a_A,
    b_A = b_A,
    sigma_gamma_prior = sigma_gamma_prior
  )
}

make_stan_data <- function(sim, model, ...) {
  model <- match.arg(model, c("sgp", "car", "ind"))
  common <- make_common_stan_data(sim, ...)
  if (model == "sgp") {
    common$sqdist1 <- sim$sqdist$sqdist1
    common$sqdist2 <- sim$sqdist$sqdist2
  } else if (model == "car") {
    common$a_A <- NULL
    common$b_A <- NULL
    W <- sim$truth$W
    degree <- as.numeric(rowSums(W))
    edges <- which(upper.tri(W) & W != 0, arr.ind = TRUE)
    normalized_W <- W / sqrt(outer(degree, degree))
    eigenvalues <- eigen(
      normalized_W,
      symmetric = TRUE,
      only.values = TRUE
    )$values
    common$degree <- degree
    common$N_edges <- nrow(edges)
    common$node1 <- as.integer(edges[, 1L])
    common$node2 <- as.integer(edges[, 2L])
    common$normalized_adjacency_eigenvalues <-
      pmin(1, pmax(-1, eigenvalues))
    common$log_degree_sum <- sum(log(degree))
  } else {
    common$a_A <- NULL
    common$b_A <- NULL
  }
  common
}

compile_matched_models <- function(
    stan_dir = "stan",
    models = c("sgp", "car", "ind"),
    force_recompile = FALSE) {
  required_simulation_packages()
  models <- match.arg(models, c("sgp", "car", "ind"), several.ok = TRUE)
  filenames <- c(
    sgp = "sgp_irt_2pl.stan",
    car = "car_irt_2pl.stan",
    ind = "ind_irt_2pl.stan"
  )
  out <- vector("list", length(models))
  names(out) <- models
  for (model in models) {
    stan_file <- file.path(stan_dir, filenames[[model]])
    if (!file.exists(stan_file)) {
      stop("Stan file not found: ", stan_file, call. = FALSE)
    }
    out[[model]] <- cmdstanr::cmdstan_model(
      stan_file = stan_file,
      cpp_options = list(stan_threads = TRUE),
      force_recompile = force_recompile
    )
  }
  out
}

fit_matched_model <- function(
    model_object,
    stan_data,
    seed,
    chains = 4L,
    parallel_chains = chains,
    threads_per_chain = 1L,
    iter_warmup = 1000L,
    iter_sampling = 1000L,
    adapt_delta = 0.99,
    max_treedepth = 15L,
    refresh = 0L) {
  timing <- system.time({
    fit <- model_object$sample(
      data = stan_data,
      seed = seed,
      chains = chains,
      parallel_chains = parallel_chains,
      threads_per_chain = threads_per_chain,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth,
      refresh = refresh,
      save_warmup = FALSE
    )
  })
  list(fit = fit, elapsed_seconds = unname(timing[["elapsed"]]))
}

ordered_draw_matrix <- function(fit, variable) {
  draws <- fit$draws(variables = variable, format = "matrix")
  draws <- as.matrix(draws)
  if (ncol(draws) <= 1L) {
    return(draws)
  }
  pattern <- paste0("^", variable, "\\[([0-9]+)\\]$")
  indices <- suppressWarnings(
    as.integer(sub(pattern, "\\1", colnames(draws)))
  )
  if (all(is.finite(indices))) {
    draws <- draws[, order(indices), drop = FALSE]
  }
  draws
}

shortest_interval_one <- function(x, probability = 0.95) {
  x <- sort(x[is.finite(x)])
  n <- length(x)
  if (n < 2L) {
    return(c(lower = NA_real_, upper = NA_real_))
  }
  k <- floor(probability * n)
  if (k >= n) {
    return(c(lower = x[1L], upper = x[n]))
  }
  widths <- x[(k + 1L):n] - x[seq_len(n - k)]
  start <- which.min(widths)
  c(lower = x[start], upper = x[start + k])
}

vector_recovery_metrics <- function(
    draws,
    truth,
    parameter,
    reference_estimate = NULL,
    probability = 0.95) {
  draws <- as.matrix(draws)
  truth <- as.numeric(truth)
  if (ncol(draws) != length(truth)) {
    stop(
      "Dimension mismatch for ", parameter, ": ",
      ncol(draws), " posterior columns versus ", length(truth), " truths.",
      call. = FALSE
    )
  }
  estimate <- colMeans(draws)
  post_sd <- apply(draws, 2L, stats::sd)
  intervals <- vapply(
    seq_len(ncol(draws)),
    function(j) shortest_interval_one(draws[, j], probability),
    numeric(2L)
  )
  lower <- intervals[1L, ]
  upper <- intervals[2L, ]
  error <- estimate - truth
  calibration <- if (
      length(truth) > 1L &&
      stats::sd(truth) > 0 &&
      all(is.finite(estimate))
    ) {
    stats::coef(stats::lm(estimate ~ truth))
  } else {
    c(NA_real_, NA_real_)
  }
  reference_rmse <- if (is.null(reference_estimate)) {
    NA_real_
  } else {
    sqrt(mean((as.numeric(reference_estimate) - truth)^2))
  }
  data.frame(
    parameter = parameter,
    dimension = length(truth),
    rmse = sqrt(mean(error^2)),
    bias = mean(error),
    coverage = mean(lower <= truth & truth <= upper),
    interval_width = mean(upper - lower),
    posterior_sd = mean(post_sd),
    truth_estimate_correlation = if (length(truth) > 1L) {
      suppressWarnings(stats::cor(truth, estimate))
    } else {
      NA_real_
    },
    calibration_intercept = unname(calibration[1L]),
    calibration_slope = unname(calibration[2L]),
    reference_rmse = reference_rmse,
    rmse_ratio_to_reference = if (is.finite(reference_rmse)) {
      sqrt(mean(error^2)) / reference_rmse
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )
}

scalar_recovery_metrics <- function(
    draws,
    truth,
    parameter,
    probability = 0.95) {
  x <- as.numeric(draws)
  interval <- shortest_interval_one(x, probability)
  estimate <- mean(x)
  data.frame(
    parameter = parameter,
    dimension = 1L,
    rmse = abs(estimate - truth),
    bias = estimate - truth,
    coverage = as.numeric(interval[1L] <= truth && truth <= interval[2L]),
    interval_width = interval[2L] - interval[1L],
    posterior_sd = stats::sd(x),
    truth_estimate_correlation = NA_real_,
    calibration_intercept = NA_real_,
    calibration_slope = NA_real_,
    reference_rmse = NA_real_,
    rmse_ratio_to_reference = NA_real_,
    stringsAsFactors = FALSE
  )
}

parameter_recovery_table <- function(fit, truth, model) {
  model <- match.arg(model, c("sgp", "car", "ind"))
  rows <- list(
    vector_recovery_metrics(
      ordered_draw_matrix(fit, "theta"), truth$theta, "theta"
    ),
    vector_recovery_metrics(
      ordered_draw_matrix(fit, "alpha"),
      truth$alpha,
      "alpha",
      reference_estimate = rep(1, length(truth$alpha))
    ),
    vector_recovery_metrics(
      ordered_draw_matrix(fit, "beta"), truth$beta, "beta"
    ),
    vector_recovery_metrics(
      ordered_draw_matrix(fit, "gamma"), truth$gamma, "gamma"
    ),
    scalar_recovery_metrics(
      ordered_draw_matrix(fit, "sigma_a"), truth$sigma_a, "sigma_a"
    )
  )

  correctly_specified_gp <-
    model == "sgp" &&
    truth$scenario %in% c("gp_aniso_low_nugget", "gp_aniso_high_nugget")
  correctly_specified_car <- model == "car" && truth$scenario == "car_truth"
  correctly_specified_ind <- model == "ind" && truth$scenario == "iid_no_spatial"

  if (correctly_specified_gp) {
    rows <- c(
      rows,
      list(
        vector_recovery_metrics(
          ordered_draw_matrix(fit, "f"), truth$f, "f"
        ),
        vector_recovery_metrics(
          ordered_draw_matrix(fit, "u"), truth$u, "u"
        ),
        scalar_recovery_metrics(
          ordered_draw_matrix(fit, "sigma_f"), truth$sigma_f, "sigma_f"
        ),
        scalar_recovery_metrics(
          ordered_draw_matrix(fit, "sigma_u"), truth$sigma_u, "sigma_u"
        ),
        scalar_recovery_metrics(
          ordered_draw_matrix(fit, "ell")[, 1L], truth$ell[1L], "ell1"
        ),
        scalar_recovery_metrics(
          ordered_draw_matrix(fit, "ell")[, 2L], truth$ell[2L], "ell2"
        ),
        scalar_recovery_metrics(
          ordered_draw_matrix(fit, "nugget_ratio"),
          truth$nugget_ratio,
          "nugget_ratio"
        ),
        scalar_recovery_metrics(
          ordered_draw_matrix(fit, "anisotropy_ratio"),
          truth$anisotropy_ratio,
          "anisotropy_ratio"
        )
      )
    )
  }

  if (correctly_specified_car) {
    rows <- c(
      rows,
      list(
        vector_recovery_metrics(
          ordered_draw_matrix(fit, "f"), truth$f, "f"
        ),
        vector_recovery_metrics(
          ordered_draw_matrix(fit, "u"), truth$u, "u"
        ),
        scalar_recovery_metrics(
          ordered_draw_matrix(fit, "sigma_f"), truth$sigma_f, "sigma_f"
        ),
        scalar_recovery_metrics(
          ordered_draw_matrix(fit, "sigma_u"), truth$sigma_u, "sigma_u"
        ),
        scalar_recovery_metrics(
          ordered_draw_matrix(fit, "rho"), truth$rho, "rho"
        )
      )
    )
  }

  if (correctly_specified_ind) {
    rows <- c(
      rows,
      list(
        scalar_recovery_metrics(
          ordered_draw_matrix(fit, "sigma_beta"),
          truth$sigma_beta,
          "sigma_beta"
        )
      )
    )
  }

  do.call(rbind, rows)
}

log_mean_exp <- function(x) {
  xmax <- max(x)
  xmax + log(mean(exp(x - xmax)))
}

rank_auc <- function(y, score) {
  y <- as.integer(y)
  n1 <- sum(y == 1L)
  n0 <- sum(y == 0L)
  if (n1 == 0L || n0 == 0L) {
    return(NA_real_)
  }
  ranks <- rank(score, ties.method = "average")
  (sum(ranks[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

predictive_metrics <- function(fit, truth) {
  p_draws <- ordered_draw_matrix(fit, "p_test")
  log_lik_draws <- ordered_draw_matrix(fit, "log_lik_test")
  y <- truth$y[truth$test_idx]
  p_true <- truth$prob[truth$test_idx]
  p_bar <- colMeans(p_draws)
  lpd <- apply(log_lik_draws, 2L, log_mean_exp)
  data.frame(
    n_test = length(y),
    neg_log_predictive_density = -mean(lpd),
    brier = mean((p_bar - y)^2),
    auc = rank_auc(y, p_bar),
    accuracy = mean(as.integer(p_bar >= 0.5) == y),
    probability_rmse = sqrt(mean((p_bar - p_true)^2)),
    stringsAsFactors = FALSE
  )
}

first_existing_column <- function(x, candidates) {
  hit <- candidates[candidates %in% names(x)]
  if (length(hit) == 0L) {
    return(NULL)
  }
  x[[hit[1L]]]
}

fit_diagnostics <- function(fit, model, ess_threshold = 400) {
  monitor <- switch(
    model,
    sgp = c(
      "theta", "alpha", "beta", "gamma", "sigma_a",
      "sigma_f", "sigma_u", "ell", "nugget_ratio"
    ),
    car = c(
      "theta", "alpha", "beta", "gamma", "sigma_a",
      "sigma_f", "sigma_u", "rho"
    ),
    ind = c(
      "theta", "alpha", "beta", "gamma", "sigma_a", "sigma_beta"
    )
  )
  summary <- fit$summary(variables = monitor)
  max_rhat <- if (all(is.na(summary$rhat))) {
    NA_real_
  } else {
    max(summary$rhat, na.rm = TRUE)
  }
  min_ess_bulk <- if (all(is.na(summary$ess_bulk))) {
    NA_real_
  } else {
    min(summary$ess_bulk, na.rm = TRUE)
  }
  min_ess_tail <- if (all(is.na(summary$ess_tail))) {
    NA_real_
  } else {
    min(summary$ess_tail, na.rm = TRUE)
  }

  sampler <- tryCatch(
    fit$diagnostic_summary(),
    error = function(e) NULL
  )
  if (is.null(sampler)) {
    divergences <- NA_integer_
    treedepth_hits <- NA_integer_
    min_ebfmi <- NA_real_
  } else {
    div <- first_existing_column(
      sampler,
      c("num_divergent", "num_divergences")
    )
    td <- first_existing_column(
      sampler,
      c("num_max_treedepth", "num_max_treedepth_exceeded")
    )
    bfmi <- first_existing_column(sampler, c("ebfmi", "e_bfmi"))
    divergences <- if (is.null(div)) NA_integer_ else sum(div)
    treedepth_hits <- if (is.null(td)) NA_integer_ else sum(td)
    min_ebfmi <- if (is.null(bfmi)) NA_real_ else min(bfmi)
  }

  diagnostic_pass <-
    is.finite(max_rhat) &&
    max_rhat <= 1.01 &&
    is.finite(min_ess_bulk) &&
    min_ess_bulk >= ess_threshold &&
    is.finite(min_ess_tail) &&
    min_ess_tail >= ess_threshold &&
    !is.na(divergences) &&
    divergences == 0L &&
    !is.na(treedepth_hits) &&
    treedepth_hits == 0L &&
    !is.na(min_ebfmi) &&
    min_ebfmi >= 0.30

  data.frame(
    fit_completed = TRUE,
    max_rhat = max_rhat,
    min_ess_bulk = min_ess_bulk,
    min_ess_tail = min_ess_tail,
    divergences = divergences,
    treedepth_hits = treedepth_hits,
    min_ebfmi = min_ebfmi,
    diagnostic_pass = diagnostic_pass,
    error_message = NA_character_,
    stringsAsFactors = FALSE
  )
}

add_run_identifiers <- function(x, scenario, model, replication, seed) {
  cbind(
    data.frame(
      scenario = scenario,
      model = model,
      replication = replication,
      seed = seed,
      stringsAsFactors = FALSE
    ),
    x
  )
}

write_csv_atomic <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = basename(path), tmpdir = dirname(path))
  utils::write.csv(x, tmp, row.names = FALSE)
  if (!file.rename(tmp, path)) {
    unlink(tmp)
    stop("Could not atomically write ", path, call. = FALSE)
  }
  invisible(path)
}

run_one_model <- function(
    sim,
    model,
    model_object,
    scenario,
    replication,
    seed,
    result_dir,
    chains = 4L,
    parallel_chains = chains,
    threads_per_chain = 1L,
    iter_warmup = 1000L,
    iter_sampling = 1000L,
    adapt_delta = 0.99,
    max_treedepth = 15L,
    a_A = 2,
    b_A = 0.1,
    sigma_gamma_prior = 1,
    save_fit = FALSE) {
  model_dir <- file.path(
    result_dir,
    scenario,
    sprintf("rep_%03d", replication)
  )
  dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

  stan_data <- make_stan_data(
    sim,
    model,
    a_A = a_A,
    b_A = b_A,
    sigma_gamma_prior = sigma_gamma_prior
  )
  result <- tryCatch(
    fit_matched_model(
      model_object = model_object,
      stan_data = stan_data,
      seed = seed,
      chains = chains,
      parallel_chains = parallel_chains,
      threads_per_chain = threads_per_chain,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth
    ),
    error = function(e) e
  )

  if (inherits(result, "error")) {
    diagnostics <- data.frame(
      fit_completed = FALSE,
      max_rhat = NA_real_,
      min_ess_bulk = NA_real_,
      min_ess_tail = NA_real_,
      divergences = NA_integer_,
      treedepth_hits = NA_integer_,
      min_ebfmi = NA_real_,
      diagnostic_pass = FALSE,
      error_message = conditionMessage(result),
      elapsed_seconds = NA_real_,
      stringsAsFactors = FALSE
    )
    diagnostics <- add_run_identifiers(
      diagnostics, scenario, model, replication, seed
    )
    write_csv_atomic(
      diagnostics,
      file.path(model_dir, paste0(model, "_diagnostics.csv"))
    )
    return(invisible(NULL))
  }

  fit <- result$fit
  postprocessed <- tryCatch(
    list(
      parameter_metrics = parameter_recovery_table(fit, sim$truth, model),
      prediction = predictive_metrics(fit, sim$truth),
      diagnostics = fit_diagnostics(fit, model)
    ),
    error = function(e) e
  )
  if (inherits(postprocessed, "error")) {
    diagnostics <- data.frame(
      fit_completed = TRUE,
      max_rhat = NA_real_,
      min_ess_bulk = NA_real_,
      min_ess_tail = NA_real_,
      divergences = NA_integer_,
      treedepth_hits = NA_integer_,
      min_ebfmi = NA_real_,
      diagnostic_pass = FALSE,
      error_message = paste(
        "Post-processing failed:",
        conditionMessage(postprocessed)
      ),
      elapsed_seconds = result$elapsed_seconds,
      stringsAsFactors = FALSE
    )
    diagnostics <- add_run_identifiers(
      diagnostics, scenario, model, replication, seed
    )
    write_csv_atomic(
      diagnostics,
      file.path(model_dir, paste0(model, "_diagnostics.csv"))
    )
    return(invisible(NULL))
  }

  parameter_metrics <- postprocessed$parameter_metrics
  prediction <- postprocessed$prediction
  diagnostics <- postprocessed$diagnostics
  diagnostics$elapsed_seconds <- result$elapsed_seconds

  parameter_metrics <- add_run_identifiers(
    parameter_metrics, scenario, model, replication, seed
  )
  prediction <- add_run_identifiers(
    prediction, scenario, model, replication, seed
  )
  diagnostics <- add_run_identifiers(
    diagnostics, scenario, model, replication, seed
  )

  write_csv_atomic(
    parameter_metrics,
    file.path(model_dir, paste0(model, "_parameter_metrics.csv"))
  )
  write_csv_atomic(
    prediction,
    file.path(model_dir, paste0(model, "_predictive_metrics.csv"))
  )
  write_csv_atomic(
    diagnostics,
    file.path(model_dir, paste0(model, "_diagnostics.csv"))
  )

  if (isTRUE(save_fit)) {
    fit$save_object(
      file = file.path(model_dir, paste0(model, "_fit.rds"))
    )
    saveRDS(
      sim,
      file = file.path(model_dir, "simulation_truth.rds")
    )
  }
  invisible(
    list(
      parameter_metrics = parameter_metrics,
      prediction = prediction,
      diagnostics = diagnostics
    )
  )
}
