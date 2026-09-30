make_anchor_target_split <- function(
    Y_test,
    anchor_fraction = 0.5,
    seed = 1L) {
  Y_test <- as.matrix(Y_test)
  I <- nrow(Y_test)
  J <- ncol(Y_test)
  n_anchor <- max(1L, min(J - 1L, round(anchor_fraction * J)))
  set.seed(seed)
  anchors <- vector("list", I)
  targets <- vector("list", I)

  # A cyclic offset keeps target exposure approximately balanced across items.
  base_order <- sample.int(J)
  for (i in seq_len(I)) {
    shifted <- base_order[((seq_len(J) + i - 2L) %% J) + 1L]
    anchors[[i]] <- sort(shifted[seq_len(n_anchor)])
    targets[[i]] <- sort(setdiff(seq_len(J), anchors[[i]]))
  }

  target_counts <- tabulate(unlist(targets), nbins = J)
  assert_true(all(target_counts > 0L), "At least one item never appears in the target set.")
  list(
    anchors = anchors,
    targets = targets,
    n_anchor = n_anchor,
    target_counts = target_counts
  )
}

extract_scoring_draws <- function(
    fit,
    J,
    max_draws = 500L,
    seed = 1L) {
  alpha <- indexed_draws(fit, "alpha", J)
  beta <- indexed_draws(fit, "beta", J)
  gamma <- indexed_draws(fit, "gamma", 1L)
  guessing <- indexed_draws(fit, "guessing", J)
  rows <- thin_draw_rows(nrow(alpha), max_draws, seed)
  list(
    alpha = alpha[rows, , drop = FALSE],
    beta = beta[rows, , drop = FALSE],
    gamma = gamma[rows, 1L],
    guessing = guessing[rows, , drop = FALSE],
    draw_rows = rows
  )
}

integrated_new_person_predictions <- function(
    y,
    vocabulary_z,
    anchor_items,
    target_items,
    scoring_draws,
    theta_grid) {
  S <- nrow(scoring_draws$alpha)
  Q <- length(theta_grid)
  T <- length(target_items)
  prior_log_weight <- stats::dnorm(theta_grid, log = TRUE)
  p_by_draw <- matrix(NA_real_, nrow = S, ncol = T)

  y <- as.numeric(y)
  vocabulary_z <- as.numeric(vocabulary_z)
  theta_grid <- as.numeric(theta_grid)
  anchor_items <- as.integer(anchor_items)
  target_items <- as.integer(target_items)
  assert_true(length(vocabulary_z) == 1L && is.finite(vocabulary_z),
              "The test respondent's standardized vocabulary must be one finite number.")
  assert_true(length(anchor_items) > 0L, "The anchor-item set is empty.")
  assert_true(T > 0L, "The target-item set is empty.")
  assert_true(all(anchor_items >= 1L & anchor_items <= length(y)),
              "Anchor-item indices are outside the response vector.")
  assert_true(all(target_items >= 1L & target_items <= length(y)),
              "Target-item indices are outside the response vector.")

  for (s in seq_len(S)) {
    # A one-row subset of posterior::draws_matrix can retain dimensions.
    # outer() would then create a higher-dimensional array, so explicitly
    # reduce every selected posterior row to an ordinary numeric vector.
    alpha_anchor <- as.numeric(scoring_draws$alpha[s, anchor_items])
    beta_anchor <- as.numeric(scoring_draws$beta[s, anchor_items])
    guess_anchor <- as.numeric(scoring_draws$guessing[s, anchor_items])
    gamma_s <- as.numeric(scoring_draws$gamma[s])
    covariate_term <- vocabulary_z * gamma_s
    assert_true(
      length(alpha_anchor) == length(anchor_items) &&
        length(beta_anchor) == length(anchor_items) &&
        length(guess_anchor) == length(anchor_items),
      "Anchor posterior-draw dimensions do not match the anchor-item set."
    )

    eta_anchor <- outer(theta_grid, alpha_anchor) -
      matrix(beta_anchor, nrow = Q, ncol = length(anchor_items), byrow = TRUE) +
      covariate_term
    probability_anchor <- sweep(
      stats::plogis(eta_anchor),
      2L,
      1 - guess_anchor,
      "*"
    )
    probability_anchor <- sweep(
      probability_anchor,
      2L,
      guess_anchor,
      "+"
    )
    probability_anchor <- clip_probability(probability_anchor)
    y_anchor <- y[anchor_items]
    log_likelihood <- rowSums(
      sweep(log(probability_anchor), 2L, y_anchor, "*") +
        sweep(log1p(-probability_anchor), 2L, 1 - y_anchor, "*")
    )
    log_weight <- prior_log_weight + log_likelihood
    weight <- exp(log_weight - log_sum_exp(log_weight))

    alpha_target <- as.numeric(scoring_draws$alpha[s, target_items])
    beta_target <- as.numeric(scoring_draws$beta[s, target_items])
    guess_target <- as.numeric(scoring_draws$guessing[s, target_items])
    assert_true(
      length(alpha_target) == T &&
        length(beta_target) == T &&
        length(guess_target) == T,
      "Target posterior-draw dimensions do not match the target-item set."
    )
    eta_target <- outer(theta_grid, alpha_target) -
      matrix(beta_target, nrow = Q, ncol = T, byrow = TRUE) +
      covariate_term
    probability_target <- sweep(
      stats::plogis(eta_target),
      2L,
      1 - guess_target,
      "*"
    )
    probability_target <- sweep(
      probability_target,
      2L,
      guess_target,
      "+"
    )
    p_by_draw[s, ] <- colSums(probability_target * weight)
  }
  p_by_draw
}

score_test_respondents <- function(
    fit,
    Y_test,
    vocabulary_test_z,
    subject_source_rows,
    anchor_target,
    theta_grid,
    max_draws = 500L,
    seed = 1L) {
  Y_test <- as.matrix(Y_test)
  I <- nrow(Y_test)
  J <- ncol(Y_test)
  assert_true(length(vocabulary_test_z) == I, "Test vocabulary length mismatch.")
  assert_true(length(subject_source_rows) == I, "Test subject ID length mismatch.")
  scoring_draws <- extract_scoring_draws(
    fit = fit,
    J = J,
    max_draws = max_draws,
    seed = seed
  )

  output <- vector("list", I)
  for (i in seq_len(I)) {
    target_items <- anchor_target$targets[[i]]
    p_by_draw <- integrated_new_person_predictions(
      y = Y_test[i, ],
      vocabulary_z = vocabulary_test_z[i],
      anchor_items = anchor_target$anchors[[i]],
      target_items = target_items,
      scoring_draws = scoring_draws,
      theta_grid = theta_grid
    )
    intervals <- t(apply(p_by_draw, 2L, shortest_interval))
    output[[i]] <- data.frame(
      subject_source_row = subject_source_rows[i],
      test_subject_index = i,
      item_index = target_items,
      item = colnames(Y_test)[target_items],
      y = as.integer(Y_test[i, target_items]),
      p_mean = colMeans(p_by_draw),
      p_median = apply(p_by_draw, 2L, stats::median),
      p_lower = intervals[, "lower"],
      p_upper = intervals[, "upper"],
      n_anchor = length(anchor_target$anchors[[i]])
    )
  }
  do.call(rbind, output)
}

item_level_prediction_metrics <- function(predictions) {
  pieces <- split(predictions, predictions$item)
  do.call(rbind, lapply(pieces, function(piece) {
    result <- prediction_metrics(piece)
    result$item <- piece$item[1L]
    result$item_index <- piece$item_index[1L]
    result
  }))
}

calibration_curve_data <- function(predictions, n_bins = 10L) {
  p <- predictions$p_mean
  breaks <- unique(stats::quantile(
    p,
    probs = seq(0, 1, length.out = n_bins + 1L),
    na.rm = TRUE,
    names = FALSE
  ))
  if (length(breaks) < 3L) {
    return(data.frame())
  }
  bin <- cut(p, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  pieces <- split(seq_len(nrow(predictions)), bin)
  do.call(rbind, lapply(seq_along(pieces), function(k) {
    ii <- pieces[[k]]
    data.frame(
      bin = k,
      n = length(ii),
      mean_predicted = mean(predictions$p_mean[ii]),
      observed = mean(predictions$y[ii])
    )
  }))
}
