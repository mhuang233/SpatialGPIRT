`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

assert_true <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
  invisible(TRUE)
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

write_csv_atomic <- function(x, path) {
  ensure_dir(dirname(path))
  tmp <- tempfile(pattern = "write-", tmpdir = dirname(path), fileext = ".csv")
  utils::write.csv(x, tmp, row.names = FALSE, na = "")
  if (file.exists(path)) {
    unlink(path)
  }
  if (!file.rename(tmp, path)) {
    unlink(tmp)
    stop("Could not write ", path, call. = FALSE)
  }
  invisible(path)
}

save_rds_atomic <- function(x, path) {
  ensure_dir(dirname(path))
  tmp <- tempfile(pattern = "write-", tmpdir = dirname(path), fileext = ".rds")
  saveRDS(x, tmp)
  if (file.exists(path)) {
    unlink(path)
  }
  if (!file.rename(tmp, path)) {
    unlink(tmp)
    stop("Could not write ", path, call. = FALSE)
  }
  invisible(path)
}

clip_probability <- function(p, eps = 1e-10) {
  assert_true(
    is.numeric(p) && length(eps) == 1L && is.finite(eps) &&
      eps > 0 && eps < 0.5,
    "`p` must be numeric and `eps` must be a finite scalar in (0, 0.5)."
  )
  result <- p
  result[!is.na(result) & result < eps] <- eps
  result[!is.na(result) & result > 1 - eps] <- 1 - eps
  result
}

log_sum_exp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

shortest_interval <- function(x, probability = 0.95) {
  x <- sort(as.numeric(x[is.finite(x)]))
  n <- length(x)
  if (n == 0L) {
    return(c(lower = NA_real_, upper = NA_real_))
  }
  if (n == 1L) {
    return(c(lower = x, upper = x))
  }
  m <- max(1L, floor(probability * n))
  if (m >= n) {
    return(c(lower = x[1L], upper = x[n]))
  }
  lower_index <- seq_len(n - m)
  upper_index <- lower_index + m
  k <- which.min(x[upper_index] - x[lower_index])
  c(lower = x[lower_index[k]], upper = x[upper_index[k]])
}

summarize_vector <- function(x, probability = 0.95) {
  hdi <- shortest_interval(x, probability)
  data.frame(
    mean = mean(x),
    sd = stats::sd(x),
    median = stats::median(x),
    lower = unname(hdi["lower"]),
    upper = unname(hdi["upper"])
  )
}

rank_auc <- function(y, p) {
  keep <- is.finite(y) & is.finite(p)
  y <- as.integer(y[keep])
  p <- p[keep]
  n1 <- sum(y == 1L)
  n0 <- sum(y == 0L)
  if (n1 == 0L || n0 == 0L) {
    return(NA_real_)
  }
  ranks <- rank(p, ties.method = "average")
  (sum(ranks[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

expected_calibration_error <- function(y, p, n_bins = 10L) {
  breaks <- unique(stats::quantile(
    p,
    probs = seq(0, 1, length.out = n_bins + 1L),
    na.rm = TRUE,
    names = FALSE
  ))
  if (length(breaks) < 3L) {
    return(NA_real_)
  }
  bin <- cut(p, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  pieces <- split(seq_along(y), bin)
  sum(vapply(
    pieces,
    function(ii) length(ii) / length(y) * abs(mean(y[ii]) - mean(p[ii])),
    numeric(1)
  ))
}

calibration_coefficients <- function(y, p) {
  lp <- stats::qlogis(clip_probability(p))
  if (length(unique(y)) < 2L || stats::sd(lp) < 1e-12) {
    return(c(intercept = NA_real_, slope = NA_real_))
  }
  fit <- suppressWarnings(stats::glm(y ~ lp, family = stats::binomial()))
  cf <- stats::coef(fit)
  c(intercept = unname(cf[1L]), slope = unname(cf[2L]))
}

prediction_metrics <- function(predictions) {
  y <- predictions$y
  p <- clip_probability(predictions$p_mean)
  cal <- calibration_coefficients(y, p)
  data.frame(
    n_target = length(y),
    prevalence = mean(y),
    mean_log_score = mean(y * log(p) + (1 - y) * log1p(-p)),
    brier = mean((y - p)^2),
    auc = rank_auc(y, p),
    accuracy = mean((p >= 0.5) == y),
    calibration_intercept = unname(cal["intercept"]),
    calibration_slope = unname(cal["slope"]),
    ece_10 = expected_calibration_error(y, p, n_bins = 10L)
  )
}

stratified_folds <- function(strata, k, seed) {
  set.seed(seed)
  strata <- as.factor(strata)
  fold <- integer(length(strata))
  for (level in levels(strata)) {
    ii <- which(strata == level)
    fold[ii] <- sample(rep(seq_len(k), length.out = length(ii)))
  }
  fold
}

make_balance_strata <- function(vocabulary, total_score = NULL, n_bins = 4L) {
  safe_cut <- function(x) {
    breaks <- unique(stats::quantile(
      x,
      probs = seq(0, 1, length.out = n_bins + 1L),
      na.rm = TRUE,
      names = FALSE
    ))
    if (length(breaks) < 3L) {
      return(factor(rep(1L, length(x))))
    }
    cut(x, breaks = breaks, include.lowest = TRUE, ordered_result = TRUE)
  }
  if (is.null(total_score)) {
    return(safe_cut(vocabulary))
  }
  interaction(safe_cut(vocabulary), safe_cut(total_score), drop = TRUE)
}

source_art_files <- function(root = ".") {
  files <- c(
    "R/00_utils.R",
    "R/01_data.R",
    "R/02_stan.R",
    "R/03_geometry.R",
    "R/04_prediction.R",
    "R/05_diagnostics.R",
    "R/06_summaries_plots.R",
    "R/07_workflow.R"
  )
  invisible(lapply(file.path(root, files), source, local = .GlobalEnv))
}

write_software_manifest <- function(output_dir, config) {
  ensure_dir(output_dir)
  packages <- c(
    "cmdstanr", "posterior", "ggplot2", "patchwork",
    "dplyr", "tidyr", "readr"
  )
  versions <- vapply(packages, function(package) {
    if (requireNamespace(package, quietly = TRUE)) {
      as.character(utils::packageVersion(package))
    } else {
      NA_character_
    }
  }, character(1))
  manifest <- data.frame(
    component = c("R", "CmdStan", packages),
    version = c(
      R.version.string,
      tryCatch(
        as.character(cmdstanr::cmdstan_version()),
        error = function(e) NA_character_
      ),
      versions
    )
  )
  write_csv_atomic(manifest, file.path(output_dir, "software_versions.csv"))
  capture <- utils::capture.output(str(config, max.level = 3L))
  writeLines(capture, file.path(output_dir, "analysis_config.txt"))
  invisible(manifest)
}
