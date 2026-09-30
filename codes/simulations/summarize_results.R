# Monte Carlo summaries and uncertainty-map helpers.

if (!exists("ordered_draw_matrix", mode = "function")) {
  source(file.path("R", "simulation_framework.R"))
}

read_result_type <- function(result_dir, suffix) {
  files <- list.files(
    result_dir,
    pattern = paste0(suffix, "[.]csv$"),
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(files) == 0L) {
    return(data.frame())
  }
  do.call(
    rbind,
    lapply(files, utils::read.csv, stringsAsFactors = FALSE)
  )
}

summarize_numeric_columns <- function(data, group_columns, value_columns) {
  if (nrow(data) == 0L) {
    return(data.frame())
  }
  group_key <- do.call(
    interaction,
    c(
      data[group_columns],
      list(drop = TRUE, lex.order = TRUE)
    )
  )
  pieces <- split(data, group_key)
  rows <- lapply(
    pieces,
    function(piece) {
      identifiers <- piece[1L, group_columns, drop = FALSE]
      summaries <- lapply(
        value_columns,
        function(variable) {
          x <- piece[[variable]]
          x <- x[is.finite(x)]
          n <- length(x)
          c(
            n = n,
            mean = if (n > 0L) mean(x) else NA_real_,
            sd = if (n > 1L) stats::sd(x) else NA_real_,
            mcse = if (n > 1L) stats::sd(x) / sqrt(n) else NA_real_
          )
        }
      )
      names(summaries) <- value_columns
      wide <- unlist(summaries)
      names(wide) <- paste(
        rep(value_columns, each = 4L),
        rep(c("n", "mean", "sd", "mcse"), times = length(value_columns)),
        sep = "_"
      )
      cbind(identifiers, as.data.frame(as.list(wide)))
    }
  )
  rownames_out <- NULL
  out <- do.call(rbind, rows)
  rownames(out) <- rownames_out
  out
}

summarize_simulation_directory <- function(
    result_dir,
    output_dir = file.path(result_dir, "summary")) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  parameter <- read_result_type(result_dir, "_parameter_metrics")
  prediction <- read_result_type(result_dir, "_predictive_metrics")
  diagnostics <- read_result_type(result_dir, "_diagnostics")
  
  parameter_summary <- summarize_numeric_columns(
    parameter,
    group_columns = c("scenario", "model", "parameter"),
    value_columns = c(
      "rmse", "bias", "coverage", "interval_width", "posterior_sd",
      "truth_estimate_correlation", "calibration_intercept",
      "calibration_slope", "reference_rmse", "rmse_ratio_to_reference"
    )
  )
  prediction_summary <- summarize_numeric_columns(
    prediction,
    group_columns = c("scenario", "model"),
    value_columns = c(
      "neg_log_predictive_density", "brier", "auc",
      "accuracy", "probability_rmse"
    )
  )
  
  if (nrow(diagnostics) > 0L) {
    diagnostics$fit_failure <- as.numeric(!diagnostics$fit_completed)
    diagnostics$diagnostic_failure <- as.numeric(!diagnostics$diagnostic_pass)
  }
  diagnostic_summary <- summarize_numeric_columns(
    diagnostics,
    group_columns = c("scenario", "model"),
    value_columns = c(
      "fit_failure", "diagnostic_failure", "max_rhat",
      "min_ess_bulk", "min_ess_tail", "divergences",
      "treedepth_hits", "min_ebfmi", "elapsed_seconds"
    )
  )
  
  write_csv_atomic(
    parameter_summary,
    file.path(output_dir, "parameter_recovery_summary.csv")
  )
  write_csv_atomic(
    prediction_summary,
    file.path(output_dir, "predictive_summary.csv")
  )
  write_csv_atomic(
    diagnostic_summary,
    file.path(output_dir, "diagnostic_runtime_summary.csv")
  )
  
  invisible(
    list(
      parameter = parameter_summary,
      prediction = prediction_summary,
      diagnostics = diagnostic_summary
    )
  )
}

# Construct Appendix Table B.1 from replication-level SGP-IRT output.
#
# For a scalar parameter, each replication-level `bias` is the posterior-mean
# error. The Monte Carlo RMSE is therefore sqrt(mean(bias^2)); averaging the
# replication-level absolute errors would instead produce an MAE.
summarize_table_B1 <- function(
    result_dir,
    output_file = file.path(result_dir, "summary", "table_B1_sgp.csv"),
    expected_replications = 5L,
    require_diagnostic_pass = FALSE) {
  parameter <- read_result_type(result_dir, "_parameter_metrics")
  diagnostics <- read_result_type(result_dir, "_diagnostics")
  
  if (nrow(parameter) == 0L) {
    stop("No parameter-metric files were found under ", result_dir, call. = FALSE)
  }
  
  keep <- parameter$model == "sgp" &
    parameter$parameter %in% c("gamma", "sigma_a")
  parameter <- parameter[keep, , drop = FALSE]
  if (nrow(parameter) == 0L) {
    stop("No SGP-IRT gamma or sigma_a metrics were found.", call. = FALSE)
  }
  
  duplicate_key <- duplicated(
    parameter[c("scenario", "replication", "parameter")]
  )
  if (any(duplicate_key)) {
    stop(
      "Duplicate scenario/replication/parameter rows were found in the results.",
      call. = FALSE
    )
  }
  
  if (isTRUE(require_diagnostic_pass)) {
    if (nrow(diagnostics) == 0L) {
      stop(
        "Diagnostic filtering was requested, but no diagnostics were found.",
        call. = FALSE
      )
    }
    diagnostics <- diagnostics[
      diagnostics$model == "sgp",
      c("scenario", "replication", "diagnostic_pass"),
      drop = FALSE
    ]
    parameter <- merge(
      parameter,
      diagnostics,
      by = c("scenario", "replication"),
      all.x = TRUE,
      sort = FALSE
    )
    parameter <- parameter[
      !is.na(parameter$diagnostic_pass) & parameter$diagnostic_pass,
      ,
      drop = FALSE
    ]
  }
  
  scenario_order <- scenario_grid()$scenario
  required_grid <- expand.grid(
    scenario = scenario_order,
    parameter = c("gamma", "sigma_a"),
    stringsAsFactors = FALSE
  )
  count_table <- stats::aggregate(
    replication ~ scenario + parameter,
    data = parameter,
    FUN = function(x) length(unique(x))
  )
  names(count_table)[names(count_table) == "replication"] <- "n_rep"
  count_check <- merge(
    required_grid,
    count_table,
    by = c("scenario", "parameter"),
    all.x = TRUE,
    sort = FALSE
  )
  count_check$n_rep[is.na(count_check$n_rep)] <- 0L
  bad_count <- count_check$n_rep != expected_replications
  if (any(bad_count)) {
    details <- paste(
      paste0(
        count_check$scenario[bad_count], "/",
        count_check$parameter[bad_count], "=",
        count_check$n_rep[bad_count]
      ),
      collapse = ", "
    )
    stop(
      "Table B.1 is incomplete; expected ", expected_replications,
      " replications per cell but found: ", details,
      call. = FALSE
    )
  }
  
  split_key <- interaction(
    factor(parameter$scenario, levels = scenario_order),
    factor(parameter$parameter, levels = c("gamma", "sigma_a")),
    drop = TRUE,
    lex.order = TRUE
  )
  pieces <- split(parameter, split_key)
  rows <- lapply(
    pieces,
    function(piece) {
      data.frame(
        scenario = piece$scenario[1L],
        parameter = piece$parameter[1L],
        truth = if (piece$parameter[1L] == "gamma") 0.50 else 0.35,
        n_rep = length(unique(piece$replication)),
        bias = mean(piece$bias),
        rmse = sqrt(mean(piece$bias^2)),
        coverage = mean(piece$coverage),
        interval_width = mean(piece$interval_width),
        posterior_sd = mean(piece$posterior_sd),
        stringsAsFactors = FALSE
      )
    }
  )
  out <- do.call(rbind, rows)
  out <- out[
    order(
      match(out$scenario, scenario_order),
      match(out$parameter, c("gamma", "sigma_a"))
    ),
    ,
    drop = FALSE
  ]
  rownames(out) <- NULL
  
  write_csv_atomic(out, output_file)
  invisible(out)
}

color_scale <- function(x, palette = "Viridis") {
  n_colors <- 100L
  breaks <- seq(min(x), max(x), length.out = n_colors + 1L)
  if (diff(range(x)) == 0) {
    index <- rep(50L, length(x))
  } else {
    index <- findInterval(x, breaks, all.inside = TRUE)
  }
  grDevices::hcl.colors(n_colors, palette = palette)[index]
}

plot_uncertainty_maps <- function(
    fit,
    sim,
    variable = c("beta", "f"),
    output_file = NULL,
    width = 1800,
    height = 650,
    resolution = 150) {
  variable <- match.arg(variable)
  draws <- ordered_draw_matrix(fit, variable)
  posterior_mean <- colMeans(draws)
  posterior_sd <- apply(draws, 2L, stats::sd)
  truth <- sim$truth[[variable]]
  coords <- sim$truth$coords
  
  if (!is.null(output_file)) {
    grDevices::png(
      filename = output_file,
      width = width,
      height = height,
      res = resolution
    )
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  
  old_par <- graphics::par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::plot(
    coords,
    pch = 21,
    bg = color_scale(truth),
    cex = 1.7,
    xlab = "coordinate 1",
    ylab = "coordinate 2",
    main = paste("True", variable)
  )
  graphics::plot(
    coords,
    pch = 21,
    bg = color_scale(posterior_mean),
    cex = 1.7,
    xlab = "coordinate 1",
    ylab = "coordinate 2",
    main = "Posterior mean"
  )
  graphics::plot(
    coords,
    pch = 21,
    bg = color_scale(posterior_sd, palette = "Inferno"),
    cex = 1.7,
    xlab = "coordinate 1",
    ylab = "coordinate 2",
    main = "Posterior SD"
  )
  invisible(
    data.frame(
      item = seq_len(nrow(coords)),
      s1 = coords[, 1L],
      s2 = coords[, 2L],
      truth = truth,
      posterior_mean = posterior_mean,
      posterior_sd = posterior_sd
    )
  )
}
