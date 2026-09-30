########################## Smoke testing and setting up the PLUS models ###########################################
###################################################################################################################

exact_plus_find_simulation_root <- function() {
  source_file <- tryCatch(
    normalizePath(sys.frame(1L)$ofile, winslash = "/", mustWork = TRUE),
    error = function(e) NA_character_
  )
  candidates <- unique(c(
    if (is.character(source_file) && !is.na(source_file)) {
      dirname(dirname(source_file))
    } else {
      character(0)
    },
    normalizePath(getwd(), winslash = "/", mustWork = TRUE),
    normalizePath(
      file.path(getwd(), "revision_simulation"),
      winslash = "/",
      mustWork = FALSE
    )
  ))
  valid <- candidates[
    file.exists(file.path(candidates, "R", "simulation_framework.R"))
  ]
  if (length(valid) == 0L) {
    stop(
      "Could not locate revision_simulation/R/simulation_framework.R.",
      call. = FALSE
    )
  }
  valid[1L]
}

EXACT_PLUS_SIMULATION_ROOT <- exact_plus_find_simulation_root()

if (!exists("simulate_irt_dgp", mode = "function") ||
    !exists("shortest_interval_one", mode = "function") ||
    !exists("vector_recovery_metrics", mode = "function")) {
  source(
    file.path(
      EXACT_PLUS_SIMULATION_ROOT,
      "R",
      "simulation_framework.R"
    )
  )
}

EXACT_PLUS_MODELS <- c(
  exact_1plus = "1PLUS",
  exact_2plus = "2PLUS",
  exact_3plus = "3PLUS"
)

stopifnot(
  identical(unname(EXACT_PLUS_MODELS), c("1PLUS", "2PLUS", "3PLUS"))
)

exact_plus_require_packages <- function() {
  base_required <- c("deldir", "rstan")
  missing <- base_required[
    !vapply(base_required, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing) > 0L) {
    stop(
      "Install the following packages before running the smoke tests: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

exact_plus_source_candidates <- function() {
  project_parent <- dirname(EXACT_PLUS_SIMULATION_ROOT)
  unique(c(
    file.path(
      project_parent,
      "art_revision",
      "staircase-source",
      "staircase-main",
      "R",
      "functions.R"
    ),
    file.path(
      project_parent,
      "art_revision",
      "staircase-main",
      "R",
      "functions.R"
    ),
    file.path(
      project_parent,
      "staircase-source",
      "staircase-main",
      "R",
      "functions.R"
    ),
    file.path(project_parent, "staircase-main", "R", "functions.R"),
    file.path(
      getwd(),
      "staircase-source",
      "staircase-main",
      "R",
      "functions.R"
    ),
    file.path(getwd(), "staircase-main", "R", "functions.R")
  ))
}

load_exact_staircase_ir_spat <- function() {
  if (requireNamespace("staircase", quietly = TRUE)) {
    return(list(
      fun = getExportedValue("staircase", "ir_spat"),
      version = as.character(utils::packageVersion("staircase")),
      implementation = "installed staircase::ir_spat"
    ))
  }
  
  source_candidates <- exact_plus_source_candidates()
  source_path <- source_candidates[file.exists(source_candidates)][1L]
  if (length(source_path) != 1L || is.na(source_path)) {
    stop(
      paste0(
        "The staircase package is unavailable and the exact published ",
        "R/functions.R source was not found. Expected one of: ",
        paste(source_candidates, collapse = "; ")
      ),
      call. = FALSE
    )
  }
  
  direct_source_packages <- c("dismo", "dplyr", "plyr", "rgeos", "rstan")
  missing <- direct_source_packages[
    !vapply(
      direct_source_packages,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )
  ]
  if (length(missing) > 0L) {
    stop(
      paste0(
        "Direct loading of the published staircase source requires: ",
        paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
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
  
  if (!exists("ir_spat", envir = source_environment, inherits = FALSE)) {
    stop(
      "The published staircase functions.R file did not define ir_spat().",
      call. = FALSE
    )
  }
  
  description_path <- file.path(dirname(dirname(source_path)), "DESCRIPTION")
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

exact_plus_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = FALSE)
  invisible(path)
}

exact_plus_covariate_group <- function(
    x,
    grouping = c("three_equal_width", "single_group")) {
  grouping <- match.arg(grouping)
  if (grouping == "single_group") {
    return(factor(rep(1L, length(x)), levels = 1L))
  }
  
  # X is generated on [-1, 1]. Fixed cut points avoid using information from
  # held-out cells and mirror the grouped-covariate adaptation used for ART.
  group <- cut(
    x,
    breaks = c(-Inf, -1 / 3, 1 / 3, Inf),
    labels = FALSE,
    include.lowest = TRUE
  )
  factor(group, levels = 1:3)
}

make_exact_plus_simulation_data <- function(
    sim,
    indices,
    covariate_grouping = c("three_equal_width", "single_group"),
    require_complete_design = FALSE) {
  covariate_grouping <- match.arg(covariate_grouping)
  truth <- sim$truth
  if (truth$p != 1L || ncol(truth$X) != 1L) {
    stop(
      "The exact PLUS smoke runner currently requires p = 1.",
      call. = FALSE
    )
  }
  
  indices <- as.integer(indices)
  item_id <- truth$item_id[indices]
  result <- data.frame(
    source_index = indices,
    subject_id = as.integer(truth$subj_id[indices]),
    item_id = as.integer(item_id),
    response = as.integer(truth$y[indices]),
    x = as.numeric(truth$X[indices, 1L]),
    true_probability = as.numeric(truth$prob[indices]),
    lon = as.numeric(truth$coords[item_id, 1L]),
    lat = as.numeric(truth$coords[item_id, 2L]),
    stringsAsFactors = FALSE
  )
  result$True_Species <- exact_plus_covariate_group(
    result$x,
    grouping = covariate_grouping
  )
  
 
  result <- result[order(result$item_id, result$subject_id), , drop = FALSE]
  rownames(result) <- NULL
  
  if (isTRUE(require_complete_design) &&
      !identical(sort(unique(result$item_id)), seq_len(truth$J))) {
    stop("Training data do not contain every item.", call. = FALSE)
  }
  if (isTRUE(require_complete_design) &&
      !identical(sort(unique(result$subject_id)), seq_len(truth$I))) {
    stop("Training data do not contain every respondent.", call. = FALSE)
  }
  if (covariate_grouping == "three_equal_width" &&
      length(unique(as.integer(result$True_Species))) != 3L) {
    stop("All three covariate groups must occur in the training data.", call. = FALSE)
  }
  result
}

exact_plus_as_draw_matrix <- function(x, expected_columns, parameter) {
  if (is.null(dim(x))) {
    x <- matrix(x, ncol = expected_columns)
  } else {
    x <- as.matrix(x)
  }
  if (ncol(x) != expected_columns) {
    stop(
      parameter,
      " draws have ",
      ncol(x),
      " columns; expected ",
      expected_columns,
      ".",
      call. = FALSE
    )
  }
  x
}

extract_exact_plus_simulation_draws <- function(
    fit,
    itemtype,
    I,
    J,
    n_groups,
    max_draws = NULL,
    seed = NULL) {
  extracted <- rstan::extract(
    fit,
    pars = c("abil", "difficulty", "diff_species", "alpha", "eta"),
    permuted = TRUE,
    inc_warmup = FALSE
  )
  
  ability <- exact_plus_as_draw_matrix(extracted$abil, I, "abil")
  difficulty <- exact_plus_as_draw_matrix(
    extracted$difficulty,
    J,
    "difficulty"
  )
  diff_species <- exact_plus_as_draw_matrix(
    extracted$diff_species,
    n_groups,
    "diff_species"
  )
  alpha <- exact_plus_as_draw_matrix(extracted$alpha, J, "alpha")
  eta <- exact_plus_as_draw_matrix(extracted$eta, n_groups, "eta")
  
  n_available <- nrow(ability)
  if (any(c(
    nrow(difficulty),
    nrow(diff_species),
    nrow(alpha),
    nrow(eta)
  ) != n_available)) {
    stop("The extracted PLUS draw arrays have inconsistent sizes.", call. = FALSE)
  }
  
  keep <- if (is.null(max_draws) || !is.finite(max_draws) ||
              n_available <= max_draws) {
    seq_len(n_available)
  } else {
    if (!is.null(seed)) set.seed(seed)
    sort(sample.int(n_available, max_draws, replace = FALSE))
  }
  
  if (itemtype == "1PLUS") {
    alpha[,] <- 1
  }
  if (itemtype != "3PLUS") {
    eta[,] <- 0
  }
  
  list(
    ability = ability[keep, , drop = FALSE],
    difficulty = difficulty[keep, , drop = FALSE],
    diff_species = diff_species[keep, , drop = FALSE],
    alpha = alpha[keep, , drop = FALSE],
    eta = eta[keep, , drop = FALSE]
  )
}

subset_exact_plus_draws <- function(draws, max_draws, seed) {
  n_available <- nrow(draws$ability)
  set.seed(seed)
  keep <- if (n_available <= max_draws) {
    seq_len(n_available)
  } else {
    sort(sample.int(n_available, max_draws, replace = FALSE))
  }
  lapply(draws, function(x) x[keep, , drop = FALSE])
}


identify_exact_plus_draws <- function(draws, itemtype) {
  ability_location <- rowMeans(draws$ability)
  centered_ability <- sweep(
    draws$ability,
    1L,
    ability_location,
    FUN = "-"
  )
  ability_scale <- sqrt(rowMeans(centered_ability^2))
  if (any(!is.finite(ability_scale) | ability_scale <= 0)) {
    stop("Invalid posterior ability scale in PLUS identification.", call. = FALSE)
  }
  
  theta <- sweep(centered_ability, 1L, ability_scale, FUN = "/")
  alpha <- sweep(draws$alpha, 1L, ability_scale, FUN = "*")
  mean_group_difficulty <- rowMeans(draws$diff_species)
  centered_item_location <- sweep(
    draws$difficulty,
    1L,
    mean_group_difficulty,
    FUN = "+"
  )
  centered_item_location <- sweep(
    centered_item_location,
    1L,
    ability_location,
    FUN = "-"
  )
  beta <- draws$alpha * centered_item_location
  
  list(
    theta = theta,
    alpha = alpha,
    beta = beta,
    ability_location = ability_location,
    ability_scale = ability_scale,
    mean_group_difficulty = mean_group_difficulty,
    alpha_item_specific = itemtype != "1PLUS"
  )
}

summarize_exact_plus_parameter <- function(
    draws,
    truth,
    parameter,
    item_specific_estimable = TRUE,
    probability = 0.95) {
  draws <- as.matrix(draws)
  truth <- as.numeric(truth)
  if (ncol(draws) != length(truth)) {
    stop(
      "Dimension mismatch for ", parameter, ": ", ncol(draws),
      " posterior columns versus ", length(truth), " truths.",
      call. = FALSE
    )
  }
  
  intervals <- vapply(
    seq_len(ncol(draws)),
    function(j) shortest_interval_one(draws[, j], probability),
    numeric(2L)
  )
  estimate <- colMeans(draws)
  lower <- intervals[1L, ]
  upper <- intervals[2L, ]
  recovery_is_applicable <- parameter != "alpha" || item_specific_estimable
  
  data.frame(
    parameter = parameter,
    index = seq_along(truth),
    truth = truth,
    posterior_mean = estimate,
    posterior_sd = apply(draws, 2L, stats::sd),
    posterior_median = apply(draws, 2L, stats::median),
    lower_95 = lower,
    upper_95 = upper,
    bias = if (recovery_is_applicable) estimate - truth else NA_real_,
    squared_error = if (recovery_is_applicable) {
      (estimate - truth)^2
    } else {
      NA_real_
    },
    covered_95 = if (recovery_is_applicable) {
      as.numeric(lower <= truth & truth <= upper)
    } else {
      NA_real_
    },
    interval_width_95 = upper - lower,
    item_specific_estimable = item_specific_estimable,
    estimate_type = if (parameter == "alpha" && !item_specific_estimable) {
      "derived common discrimination; native 1PLUS alpha_j fixed to 1"
    } else {
      "posterior parameter"
    },
    identification = "draw-wise mean(theta)=0 and mean(theta^2)=1",
    stringsAsFactors = FALSE
  )
}

exact_plus_parameter_outputs <- function(
    identified_draws,
    truth,
    itemtype,
    probability = 0.95) {
  alpha_item_specific <- isTRUE(identified_draws$alpha_item_specific)
  summaries <- rbind(
    summarize_exact_plus_parameter(
      identified_draws$theta,
      truth$theta,
      parameter = "theta",
      probability = probability
    ),
    summarize_exact_plus_parameter(
      identified_draws$alpha,
      truth$alpha,
      parameter = "alpha",
      item_specific_estimable = alpha_item_specific,
      probability = probability
    ),
    summarize_exact_plus_parameter(
      identified_draws$beta,
      truth$beta,
      parameter = "beta",
      probability = probability
    )
  )
  
  recovery <- rbind(
    vector_recovery_metrics(
      identified_draws$theta,
      truth$theta,
      parameter = "theta",
      probability = probability
    ),
    if (alpha_item_specific) {
      vector_recovery_metrics(
        identified_draws$alpha,
        truth$alpha,
        parameter = "alpha",
        probability = probability
      )
    } else {
      data.frame(
        parameter = "alpha",
        dimension = length(truth$alpha),
        rmse = NA_real_,
        bias = NA_real_,
        coverage = NA_real_,
        interval_width = NA_real_,
        posterior_sd = NA_real_,
        truth_estimate_correlation = NA_real_,
        calibration_intercept = NA_real_,
        calibration_slope = NA_real_,
        reference_rmse = NA_real_,
        rmse_ratio_to_reference = NA_real_,
        stringsAsFactors = FALSE
      )
    },
    vector_recovery_metrics(
      identified_draws$beta,
      truth$beta,
      parameter = "beta",
      probability = probability
    )
  )
  recovery$item_specific_estimable <- c(TRUE, alpha_item_specific, TRUE)
  recovery$note <- c(
    "identified ability",
    if (alpha_item_specific) {
      "identified item discrimination"
    } else {
      "not applicable: 1PLUS fixes item-specific alpha_j"
    },
    "identified item intercept beta, not raw PLUS difficulty"
  )
  
  transform_summary <- data.frame(
    quantity = c(
      "native_ability_location",
      "native_ability_scale",
      "mean_group_difficulty"
    ),
    posterior_mean = c(
      mean(identified_draws$ability_location),
      mean(identified_draws$ability_scale),
      mean(identified_draws$mean_group_difficulty)
    ),
    posterior_sd = c(
      stats::sd(identified_draws$ability_location),
      stats::sd(identified_draws$ability_scale),
      stats::sd(identified_draws$mean_group_difficulty)
    ),
    stringsAsFactors = FALSE
  )
  
  list(
    summaries = summaries,
    recovery = recovery,
    transform_summary = transform_summary
  )
}

score_exact_plus_simulation_test <- function(test_data, draws, itemtype) {
  n_draws <- nrow(draws$ability)
  n_test <- nrow(test_data)
  subject_id <- test_data$subject_id
  item_id <- test_data$item_id
  group_id <- as.integer(test_data$True_Species)
  probabilities <- matrix(NA_real_, nrow = n_test, ncol = n_draws)
  
  for (s in seq_len(n_draws)) {
    centered <- draws$ability[s, subject_id] -
      draws$difficulty[s, item_id] -
      draws$diff_species[s, group_id]
    base_probability <- stats::plogis(
      draws$alpha[s, item_id] * centered
    )
    if (itemtype == "3PLUS") {
      guessing <- draws$eta[s, group_id]
      probabilities[, s] <- guessing + (1 - guessing) * base_probability
    } else {
      probabilities[, s] <- base_probability
    }
  }
  
  p_mean <- rowMeans(probabilities)
  p_mean <- pmin(1 - 1e-12, pmax(1e-12, p_mean))
  data.frame(
    source_index = test_data$source_index,
    subject_id = subject_id,
    item_id = item_id,
    covariate_group = group_id,
    y = test_data$response,
    true_probability = test_data$true_probability,
    p_mean = p_mean,
    stringsAsFactors = FALSE
  )
}

exact_plus_auc <- function(y, probability) {
  y <- as.integer(y)
  n_positive <- sum(y == 1L)
  n_negative <- sum(y == 0L)
  if (n_positive == 0L || n_negative == 0L) {
    return(NA_real_)
  }
  ranks <- rank(probability, ties.method = "average")
  (
    sum(ranks[y == 1L]) - n_positive * (n_positive + 1) / 2
  ) / (n_positive * n_negative)
}

exact_plus_prediction_metrics <- function(predictions) {
  y <- predictions$y
  p <- predictions$p_mean
  data.frame(
    n_test = nrow(predictions),
    mean_log_score = mean(y * log(p) + (1 - y) * log1p(-p)),
    brier = mean((y - p)^2),
    auc = exact_plus_auc(y, p),
    accuracy = mean(as.integer(p >= 0.5) == y),
    probability_rmse = sqrt(mean(
      (p - predictions$true_probability)^2
    ))
  )
}

exact_plus_finite_summary <- function(x, fun) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else fun(x)
}

exact_plus_fit_diagnostics <- function(
    fit,
    chains,
    iter_warmup,
    iter_sampling,
    wall_seconds,
    warning_count) {
  fit_summary <- summary(fit)$summary
  sampler <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  
  divergences <- sum(vapply(
    sampler,
    function(x) {
      if ("divergent__" %in% colnames(x)) sum(x[, "divergent__"]) else 0
    },
    numeric(1)
  ))
  treedepth_hits <- sum(vapply(
    sampler,
    function(x) {
      if ("treedepth__" %in% colnames(x)) {
        sum(x[, "treedepth__"] >= 10L)
      } else {
        0
      }
    },
    numeric(1)
  ))
  ebfmi <- vapply(
    sampler,
    function(x) {
      if (!("energy__" %in% colnames(x))) return(NA_real_)
      energy <- x[, "energy__"]
      if (length(energy) < 2L || !is.finite(stats::var(energy)) ||
          stats::var(energy) <= 0) {
        return(NA_real_)
      }
      mean(diff(energy)^2) / stats::var(energy)
    },
    numeric(1)
  )
  
  data.frame(
    max_rhat = if ("Rhat" %in% colnames(fit_summary)) {
      exact_plus_finite_summary(fit_summary[, "Rhat"], max)
    } else {
      NA_real_
    },
    min_n_eff = if ("n_eff" %in% colnames(fit_summary)) {
      exact_plus_finite_summary(fit_summary[, "n_eff"], min)
    } else {
      NA_real_
    },
    divergences = divergences,
    treedepth_hits = treedepth_hits,
    min_ebfmi = exact_plus_finite_summary(ebfmi, min),
    wall_seconds = wall_seconds,
    warning_count = warning_count,
    chains = chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    adapt_delta_setting = 0.8,
    max_treedepth_setting = 10L
  )
}

fit_one_exact_plus_smoke_test <- function(
    exact_ir_spat,
    sim,
    scenario,
    replication,
    model,
    itemtype,
    seed,
    output_dir,
    covariate_grouping,
    chains,
    iter_warmup,
    iter_sampling,
    max_scoring_draws,
    refresh,
    save_fit) {
  model_dir <- file.path(
    output_dir,
    scenario,
    sprintf("rep_%03d", replication),
    model
  )
  dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
  
  train_data <- make_exact_plus_simulation_data(
    sim,
    indices = sim$truth$train_idx,
    covariate_grouping = covariate_grouping,
    require_complete_design = TRUE
  )
  test_data <- make_exact_plus_simulation_data(
    sim,
    indices = sim$truth$test_idx,
    covariate_grouping = covariate_grouping,
    require_complete_design = FALSE
  )
  n_groups <- nlevels(train_data$True_Species)
  
  warning_messages <- character(0)
  started <- Sys.time()
  fit <- withCallingHandlers(
    tryCatch(
      exact_ir_spat(
        data = train_data,
        spat_model = "car",
        itemtype = itemtype,
        abil = "subject_id",
        diff = "item_id",
        y = "response",
        coords = c("lon", "lat"),
        iter = iter_warmup + iter_sampling,
        warmup = iter_warmup,
        chains = chains,
        loglik = TRUE,
        refresh = refresh,
        seed = seed
      ),
      error = function(e) e
    ),
    warning = function(w) {
      warning_messages <<- c(warning_messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  wall_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  warning_messages <- unique(warning_messages)
  writeLines(warning_messages, file.path(model_dir, "warnings.txt"))
  
  if (inherits(fit, "error") || !inherits(fit, "stanfit")) {
    error_message <- if (inherits(fit, "condition")) {
      conditionMessage(fit)
    } else {
      paste(capture.output(str(fit)), collapse = " ")
    }
    writeLines(error_message, file.path(model_dir, "error.txt"))
    failure <- data.frame(
      scenario = scenario,
      replication = replication,
      model = model,
      itemtype = itemtype,
      success = FALSE,
      wall_seconds = wall_seconds,
      warning_count = length(warning_messages),
      error = error_message,
      stringsAsFactors = FALSE
    )
    exact_plus_write_csv(failure, file.path(model_dir, "run_status.csv"))
    return(failure)
  }
  
  if (isTRUE(save_fit)) {
    saveRDS(fit, file.path(model_dir, "fit.rds"))
  }
  
  all_draws <- extract_exact_plus_simulation_draws(
    fit = fit,
    itemtype = itemtype,
    I = sim$truth$I,
    J = sim$truth$J,
    n_groups = n_groups,
    max_draws = NULL
  )
  identified_draws <- identify_exact_plus_draws(
    draws = all_draws,
    itemtype = itemtype
  )
  parameter_outputs <- exact_plus_parameter_outputs(
    identified_draws = identified_draws,
    truth = sim$truth,
    itemtype = itemtype
  )
  scoring_draws <- subset_exact_plus_draws(
    draws = all_draws,
    max_draws = max_scoring_draws,
    seed = seed + 1L
  )
  predictions <- score_exact_plus_simulation_test(
    test_data = test_data,
    draws = scoring_draws,
    itemtype = itemtype
  )
  metrics <- exact_plus_prediction_metrics(predictions)
  metrics$scenario <- scenario
  metrics$replication <- replication
  metrics$model <- model
  metrics$itemtype <- itemtype
  metrics$success <- TRUE
  metrics$wall_seconds <- wall_seconds
  
  diagnostics <- exact_plus_fit_diagnostics(
    fit = fit,
    chains = chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    wall_seconds = wall_seconds,
    warning_count = length(warning_messages)
  )
  diagnostics$scenario <- scenario
  diagnostics$replication <- replication
  diagnostics$model <- model
  diagnostics$itemtype <- itemtype
  
  parameter_outputs$summaries$scenario <- scenario
  parameter_outputs$summaries$replication <- replication
  parameter_outputs$summaries$model <- model
  parameter_outputs$summaries$itemtype <- itemtype
  parameter_outputs$recovery$scenario <- scenario
  parameter_outputs$recovery$replication <- replication
  parameter_outputs$recovery$model <- model
  parameter_outputs$recovery$itemtype <- itemtype
  parameter_outputs$transform_summary$scenario <- scenario
  parameter_outputs$transform_summary$replication <- replication
  parameter_outputs$transform_summary$model <- model
  parameter_outputs$transform_summary$itemtype <- itemtype
  
  run_status <- data.frame(
    scenario = scenario,
    replication = replication,
    model = model,
    itemtype = itemtype,
    success = TRUE,
    wall_seconds = wall_seconds,
    warning_count = length(warning_messages),
    error = "",
    stringsAsFactors = FALSE
  )
  
  exact_plus_write_csv(predictions, file.path(model_dir, "predictions.csv"))
  exact_plus_write_csv(metrics, file.path(model_dir, "prediction_metrics.csv"))
  exact_plus_write_csv(diagnostics, file.path(model_dir, "mcmc_diagnostics.csv"))
  exact_plus_write_csv(
    parameter_outputs$summaries,
    file.path(model_dir, "posterior_parameter_summaries.csv")
  )
  exact_plus_write_csv(
    parameter_outputs$recovery,
    file.path(model_dir, "parameter_recovery_metrics.csv")
  )
  exact_plus_write_csv(
    parameter_outputs$transform_summary,
    file.path(model_dir, "identification_transform_summary.csv")
  )
  exact_plus_write_csv(run_status, file.path(model_dir, "run_status.csv"))
  run_status
}

read_exact_plus_smoke_statuses <- function(output_dir) {
  paths <- list.files(
    output_dir,
    pattern = "^run_status[.]csv$",
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(paths) == 0L) return(data.frame())
  do.call(
    rbind,
    lapply(paths, utils::read.csv, stringsAsFactors = FALSE)
  )
}

read_exact_plus_nested_csv <- function(output_dir, filename) {
  paths <- list.files(
    output_dir,
    pattern = paste0("^", gsub("[.]", "[.]", filename), "$"),
    recursive = TRUE,
    full.names = TRUE
  )
  root_copy <- normalizePath(
    file.path(output_dir, filename),
    winslash = "/",
    mustWork = FALSE
  )
  normalized_paths <- normalizePath(
    paths,
    winslash = "/",
    mustWork = FALSE
  )
  paths <- paths[normalized_paths != root_copy]
  if (length(paths) == 0L) return(data.frame())
  do.call(
    rbind,
    lapply(paths, utils::read.csv, stringsAsFactors = FALSE)
  )
}

run_exact_plus_simulation_smoke_tests <- function(
    scenarios = scenario_grid()$scenario,
    models = names(EXACT_PLUS_MODELS),
    n_rep = 1L,
    I = 200L,
    J = 50L,
    p = 1L,
    sigma_a = 0.35,
    gamma = 0.50,
    test_fraction = 0.20,
    covariate_grouping = c("three_equal_width", "single_group"),
    base_seed = 202600L,
    output_dir = file.path(
      EXACT_PLUS_SIMULATION_ROOT,
      "results",
      "exact_plus_smoke"
    ),
    chains = 1L,
    iter_warmup = 150L,
    iter_sampling = 150L,
    max_scoring_draws = 100L,
    refresh = 50L,
    save_fits = FALSE,
    skip_completed = TRUE) {
  exact_plus_require_packages()
  scenarios <- match.arg(
    scenarios,
    scenario_grid()$scenario,
    several.ok = TRUE
  )
  models <- match.arg(
    models,
    names(EXACT_PLUS_MODELS),
    several.ok = TRUE
  )
  covariate_grouping <- match.arg(covariate_grouping)
  
  stopifnot(
    n_rep >= 1L,
    I >= 10L,
    J >= 4L,
    p == 1L,
    chains >= 1L,
    iter_warmup >= 1L,
    iter_sampling >= 1L,
    max_scoring_draws >= 1L
  )
  
  staircase_runtime <- load_exact_staircase_ir_spat()
  exact_ir_spat <- staircase_runtime$fun
  message("Exact PLUS implementation: ", staircase_runtime$implementation)
  message(
    "Smoke-test settings: ",
    n_rep,
    " replication(s), I = ",
    I,
    ", J = ",
    J,
    ", chains = ",
    chains,
    ", warmup = ",
    iter_warmup,
    ", sampling = ",
    iter_sampling
  )
  
  # The legacy RStan implementation can lose PSOCK workers on Windows. Running
  # chains sequentially changes scheduling only and is safest for smoke tests.
  old_cores <- getOption("mc.cores")
  on.exit(options(mc.cores = old_cores), add = TRUE)
  options(mc.cores = 1L)
  rstan::rstan_options(auto_write = TRUE)
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  manifest <- expand.grid(
    scenario = scenarios,
    replication = seq_len(n_rep),
    model = models,
    stringsAsFactors = FALSE
  )
  manifest$itemtype <- unname(EXACT_PLUS_MODELS[manifest$model])
  manifest$I <- I
  manifest$J <- J
  manifest$p <- p
  manifest$sigma_a <- sigma_a
  manifest$gamma <- gamma
  manifest$test_fraction <- test_fraction
  manifest$covariate_grouping <- covariate_grouping
  manifest$chains <- chains
  manifest$iter_warmup <- iter_warmup
  manifest$iter_sampling <- iter_sampling
  manifest$base_seed <- base_seed
  manifest$data_seed <- base_seed +
    10000L * match(manifest$scenario, scenario_grid()$scenario) +
    manifest$replication
  manifest$fit_seed <- manifest$data_seed +
    1000L * match(manifest$model, names(EXACT_PLUS_MODELS))
  manifest$staircase_version <- staircase_runtime$version
  manifest$staircase_implementation <- staircase_runtime$implementation
  exact_plus_write_csv(
    manifest,
    file.path(output_dir, "smoke_test_manifest.csv")
  )
  
 
  for (scenario in scenarios) {
    for (replication in seq_len(n_rep)) {
      scenario_index <- match(scenario, scenario_grid()$scenario)
      data_seed <- base_seed + 10000L * scenario_index + replication
      message(
        "Scenario ",
        scenario,
        ", replication ",
        replication,
        " of ",
        n_rep
      )
      sim <- simulate_irt_dgp(
        scenario = scenario,
        I = I,
        J = J,
        p = p,
        sigma_a = sigma_a,
        gamma = gamma,
        test_fraction = test_fraction,
        seed = data_seed
      )
      scenario_dir <- file.path(
        output_dir,
        scenario,
        sprintf("rep_%03d", replication)
      )
      dir.create(scenario_dir, recursive = TRUE, showWarnings = FALSE)
      saveRDS(sim$truth, file.path(scenario_dir, "truth.rds"))
      
      for (model in models) {
        itemtype <- unname(EXACT_PLUS_MODELS[model])
        model_dir <- file.path(scenario_dir, model)
        status_path <- file.path(model_dir, "run_status.csv")
        parameter_summary_path <- file.path(
          model_dir,
          "posterior_parameter_summaries.csv"
        )
        recovery_path <- file.path(model_dir, "parameter_recovery_metrics.csv")
        if (isTRUE(skip_completed) && file.exists(status_path) &&
            file.exists(parameter_summary_path) && file.exists(recovery_path)) {
          existing <- utils::read.csv(
            status_path,
            stringsAsFactors = FALSE
          )
          if (nrow(existing) == 1L && isTRUE(existing$success[1L])) {
            message("  Reusing completed ", itemtype)
            next
          }
        }
        
        fit_seed <- data_seed +
          1000L * match(model, names(EXACT_PLUS_MODELS))
        message("  Fitting exact ", itemtype, " (CAR)")
        tryCatch(
          fit_one_exact_plus_smoke_test(
            exact_ir_spat = exact_ir_spat,
            sim = sim,
            scenario = scenario,
            replication = replication,
            model = model,
            itemtype = itemtype,
            seed = fit_seed,
            output_dir = output_dir,
            covariate_grouping = covariate_grouping,
            chains = chains,
            iter_warmup = iter_warmup,
            iter_sampling = iter_sampling,
            max_scoring_draws = max_scoring_draws,
            refresh = refresh,
            save_fit = save_fits
          ),
          error = function(e) {
            dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
            error_message <- conditionMessage(e)
            writeLines(error_message, file.path(model_dir, "error.txt"))
            failure <- data.frame(
              scenario = scenario,
              replication = replication,
              model = model,
              itemtype = itemtype,
              success = FALSE,
              wall_seconds = NA_real_,
              warning_count = NA_integer_,
              error = error_message,
              stringsAsFactors = FALSE
            )
            exact_plus_write_csv(
              failure,
              file.path(model_dir, "run_status.csv")
            )
            message("    FAILED: ", error_message)
            failure
          }
        )
      }
    }
  }
  
  summary <- read_exact_plus_smoke_statuses(output_dir)
  if (nrow(summary) > 0L) {
    summary <- summary[
      order(summary$scenario, summary$replication, summary$model),
      ,
      drop = FALSE
    ]
    exact_plus_write_csv(
      summary,
      file.path(output_dir, "smoke_test_summary.csv")
    )
    message(
      "Completed ",
      sum(summary$success),
      " of ",
      nrow(summary),
      " exact PLUS smoke-test fits successfully."
    )
  }
  
  all_metrics <- read_exact_plus_nested_csv(
    output_dir,
    "prediction_metrics.csv"
  )
  if (nrow(all_metrics) > 0L) {
    all_metrics <- all_metrics[
      order(
        all_metrics$scenario,
        all_metrics$replication,
        all_metrics$model
      ),
      ,
      drop = FALSE
    ]
    exact_plus_write_csv(
      all_metrics,
      file.path(output_dir, "all_prediction_metrics.csv")
    )
  }
  
  all_diagnostics <- read_exact_plus_nested_csv(
    output_dir,
    "mcmc_diagnostics.csv"
  )
  if (nrow(all_diagnostics) > 0L) {
    all_diagnostics <- all_diagnostics[
      order(
        all_diagnostics$scenario,
        all_diagnostics$replication,
        all_diagnostics$model
      ),
      ,
      drop = FALSE
    ]
    exact_plus_write_csv(
      all_diagnostics,
      file.path(output_dir, "all_mcmc_diagnostics.csv")
    )
  }
  
  all_parameter_summaries <- read_exact_plus_nested_csv(
    output_dir,
    "posterior_parameter_summaries.csv"
  )
  if (nrow(all_parameter_summaries) > 0L) {
    all_parameter_summaries <- all_parameter_summaries[
      order(
        all_parameter_summaries$scenario,
        all_parameter_summaries$replication,
        all_parameter_summaries$model,
        all_parameter_summaries$parameter,
        all_parameter_summaries$index
      ),
      ,
      drop = FALSE
    ]
    exact_plus_write_csv(
      all_parameter_summaries,
      file.path(output_dir, "all_posterior_parameter_summaries.csv")
    )
  }
  
  all_recovery <- read_exact_plus_nested_csv(
    output_dir,
    "parameter_recovery_metrics.csv"
  )
  if (nrow(all_recovery) > 0L) {
    all_recovery <- all_recovery[
      order(
        all_recovery$scenario,
        all_recovery$replication,
        all_recovery$model,
        all_recovery$parameter
      ),
      ,
      drop = FALSE
    ]
    exact_plus_write_csv(
      all_recovery,
      file.path(output_dir, "all_parameter_recovery_metrics.csv")
    )
  }
  
  invisible(summary)
}
