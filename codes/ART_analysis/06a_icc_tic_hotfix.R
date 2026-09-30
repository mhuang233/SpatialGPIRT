# Runtime override for posterior_icc_tic().
#
# Source this file AFTER source_art_files(".").  It replaces any legacy copy
# of posterior_icc_tic() that may remain in R/06_summaries_plots.R.

posterior_icc_tic <- function(
    fit,
    theta_grid = seq(-3, 3, length.out = 201L),
    vocabulary_z = 0,
    max_draws = 1000L,
    seed = 1L) {
  theta_grid <- as.numeric(theta_grid)
  vocabulary_z <- as.numeric(vocabulary_z)
  assert_true(
    length(vocabulary_z) == 1L && is.finite(vocabulary_z),
    "The standardized vocabulary value must be one finite number."
  )

  alpha <- indexed_draws(fit, "alpha")
  J <- ncol(alpha)
  beta <- indexed_draws(fit, "beta", J)
  gamma <- as.numeric(indexed_draws(fit, "gamma", 1L)[, 1L])
  guessing <- indexed_draws(fit, "guessing", J)

  assert_true(J >= 1L, "No item discrimination draws were found.")
  assert_true(
    nrow(alpha) == nrow(beta) &&
      nrow(alpha) == nrow(guessing) &&
      nrow(alpha) == length(gamma),
    "ICC/TIC posterior blocks have inconsistent draw counts."
  )

  rows <- thin_draw_rows(nrow(alpha), max_draws, seed)
  alpha <- alpha[rows, , drop = FALSE]
  beta <- beta[rows, , drop = FALSE]
  gamma <- gamma[rows]
  guessing <- guessing[rows, , drop = FALSE]

  Q <- length(theta_grid)
  S <- length(rows)
  icc_draw <- matrix(NA_real_, nrow = S, ncol = Q)
  tic_draw <- matrix(NA_real_, nrow = S, ncol = Q)

  for (s in seq_len(S)) {
    alpha_s <- as.numeric(alpha[s, , drop = TRUE])
    beta_s <- as.numeric(beta[s, , drop = TRUE])
    gamma_s <- as.numeric(gamma[s])
    guessing_s <- as.numeric(guessing[s, , drop = TRUE])

    assert_true(
      length(alpha_s) == J &&
        length(beta_s) == J &&
        length(guessing_s) == J,
      "ICC/TIC item-vector dimensions do not match."
    )

    eta <- outer(theta_grid, alpha_s) -
      matrix(beta_s, nrow = Q, ncol = J, byrow = TRUE) +
      vocabulary_z * gamma_s

    base_probability <- stats::plogis(eta)
    probability <- sweep(
      base_probability,
      2L,
      1 - guessing_s,
      "*"
    )
    probability <- sweep(
      probability,
      2L,
      guessing_s,
      "+"
    )
    probability <- clip_probability(probability)

    icc_draw[s, ] <- rowMeans(probability)

    # Exact Bernoulli information in theta.  For SGP-IRT guessing_s = 0,
    # reducing this expression to sum_j alpha_j^2 p_j(1-p_j).
    derivative <- sweep(
      base_probability * (1 - base_probability),
      2L,
      (1 - guessing_s) * alpha_s,
      "*"
    )
    tic_draw[s, ] <- rowSums(
      derivative^2 / (probability * (1 - probability))
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

ART_ICC_TIC_HOTFIX_VERSION <- "2026-08-21-explicit-dimensions-v1"

