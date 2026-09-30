ART_GEOMETRY_PATCH_VERSION <- "2026-07-30-diag-v2"

posterior_expected_probabilities <- function(
    fit,
    vocabulary_z,
    max_draws = 500L,
    seed = 1L) {
  vocabulary_z <- as.numeric(vocabulary_z)
  I <- length(vocabulary_z)
  theta <- indexed_draws(fit, "theta", I)
  alpha <- indexed_draws(fit, "alpha")
  beta <- indexed_draws(fit, "beta", ncol(alpha))
  gamma <- indexed_draws(fit, "gamma", 1L)
  rows <- thin_draw_rows(nrow(theta), max_draws, seed)

  theta <- theta[rows, , drop = FALSE]
  alpha <- alpha[rows, , drop = FALSE]
  beta <- beta[rows, , drop = FALSE]
  gamma <- as.numeric(gamma[rows, 1L])
  J <- ncol(alpha)
  assert_true(ncol(theta) == I, "Theta draws do not match the number of respondents.")
  assert_true(ncol(beta) == J, "Beta draws do not match the number of items.")
  assert_true(
    nrow(alpha) == nrow(theta) &&
      nrow(beta) == nrow(theta) &&
      length(gamma) == nrow(theta),
    "Posterior variables do not contain the same number of draws."
  )
  probability <- matrix(0, nrow = I, ncol = J)

  for (s in seq_along(rows)) {
    theta_s <- as.numeric(theta[s, , drop = TRUE])
    alpha_s <- as.numeric(alpha[s, , drop = TRUE])
    beta_s <- as.numeric(beta[s, , drop = TRUE])
    eta <- tcrossprod(theta_s, alpha_s)
    eta <- eta -
      matrix(beta_s, nrow = I, ncol = J, byrow = TRUE) +
      matrix(
        vocabulary_z * gamma[s],
        nrow = I,
        ncol = J,
        byrow = FALSE
      )
    assert_true(
      identical(dim(eta), c(I, J)),
      "The reconstructed predictor is not an I-by-J matrix."
    )
    probability <- probability + stats::plogis(eta)
  }
  probability / length(rows)
}

standardized_irt_residuals <- function(Y, probability, bound = 5) {
  Y <- as.matrix(Y)
  probability <- as.matrix(probability)
  assert_true(
    identical(dim(Y), dim(probability)),
    "Y and posterior probabilities must have identical dimensions."
  )
  probability <- clip_probability(probability, eps = 1e-6)
  residual <- (Y - probability) / sqrt(probability * (1 - probability))
  residual[residual < -bound] <- -bound
  residual[residual > bound] <- bound
  residual
}

residual_association <- function(residual, shrinkage = 0.05) {
  residual <- as.matrix(residual)
  assert_true(
    nrow(residual) >= 2L && ncol(residual) >= 2L,
    "Residual association requires a respondent-by-item matrix with at least two rows and two columns."
  )
  association <- stats::cor(residual, use = "pairwise.complete.obs")
  association <- as.matrix(association)
  J <- ncol(residual)
  assert_true(
    identical(dim(association), c(J, J)),
    paste0(
      "Item association has dimensions ",
      paste(dim(association), collapse = " x "),
      "; expected ", J, " x ", J, "."
    )
  )
  association[!is.finite(association)] <- 0
  association <- (association + t(association)) / 2
  diagonal_index <- cbind(seq_len(J), seq_len(J))
  association[diagonal_index] <- 1
  association <- (1 - shrinkage) * association +
    shrinkage * diag(J)
  association <- as.matrix(association)
  association[association < -0.999] <- -0.999
  association[association > 0.999] <- 0.999
  association[diagonal_index] <- 1
  association
}

correlation_dissimilarity <- function(association) {
  association <- as.matrix(association)
  assert_true(
    nrow(association) == ncol(association),
    "The item-association object must be a square matrix."
  )
  J <- nrow(association)
  squared_distance <- 2 * (1 - association)
  squared_distance[squared_distance < 0] <- 0
  distance <- as.matrix(sqrt(squared_distance))
  assert_true(
    identical(dim(distance), c(J, J)),
    "The correlation dissimilarity did not produce a square matrix."
  )
  diagonal_index <- cbind(seq_len(J), seq_len(J))
  distance[diagonal_index] <- 0
  stats::as.dist(distance)
}

rescale_coordinates <- function(points) {
  points <- as.matrix(points)
  result <- points
  for (d in seq_len(ncol(points))) {
    limits <- range(points[, d])
    if (diff(limits) < 1e-12) {
      result[, d] <- 0.5
    } else {
      result[, d] <- (points[, d] - limits[1L]) / diff(limits)
    }
  }
  colnames(result) <- paste0("MDS", seq_len(ncol(result)))
  result
}

classical_mds_geometry <- function(association, dimensions = 2L) {
  J <- nrow(association)
  dimensions <- min(as.integer(dimensions), J - 1L)
  dissimilarity <- correlation_dissimilarity(association)
  fit <- stats::cmdscale(
    dissimilarity,
    k = dimensions,
    eig = TRUE,
    add = TRUE,
    x.ret = TRUE
  )
  raw_points <- as.matrix(fit$points)
  if (ncol(raw_points) < 2L) {
    raw_points <- cbind(raw_points, 0)
  }
  points <- rescale_coordinates(raw_points[, 1:2, drop = FALSE])
  observed_distance <- as.vector(dissimilarity)
  embedded_distance <- as.vector(stats::dist(raw_points))
  stress <- sqrt(
    sum((observed_distance - embedded_distance)^2) /
      sum(observed_distance^2)
  )
  distance_correlation <- suppressWarnings(stats::cor(
    observed_distance,
    embedded_distance,
    method = "spearman"
  ))
  list(
    coords = points,
    raw_coords = raw_points,
    association = association,
    dissimilarity = as.matrix(dissimilarity),
    eigenvalues = fit$eig,
    additive_constant = fit$ac,
    gof = fit$GOF,
    stress = stress,
    distance_correlation = distance_correlation
  )
}

graph_is_connected <- function(adjacency) {
  J <- nrow(adjacency)
  visited <- rep(FALSE, J)
  queue <- 1L
  visited[1L] <- TRUE
  while (length(queue) > 0L) {
    current <- queue[1L]
    queue <- queue[-1L]
    neighbors <- which(adjacency[current, ] > 0 & !visited)
    if (length(neighbors) > 0L) {
      visited[neighbors] <- TRUE
      queue <- c(queue, neighbors)
    }
  }
  all(visited)
}

build_car_graph <- function(coords, starting_neighbors = 3L) {
  coords <- as.matrix(coords)
  J <- nrow(coords)
  assert_true(
    J >= 2L && ncol(coords) >= 1L,
    "CAR graph construction requires at least two items with matrix-valued coordinates."
  )
  distance <- as.matrix(stats::dist(coords))
  assert_true(
    identical(dim(distance), c(J, J)),
    "The item-distance object is not a J-by-J matrix."
  )
  diagonal_index <- cbind(seq_len(J), seq_len(J))
  distance[diagonal_index] <- Inf
  k <- min(max(1L, as.integer(starting_neighbors)), J - 1L)

  repeat {
    adjacency <- matrix(0, nrow = J, ncol = J)
    for (j in seq_len(J)) {
      neighbors <- order(distance[j, ])[seq_len(k)]
      adjacency[j, neighbors] <- 1
    }
    adjacency <- 1L * ((adjacency + t(adjacency)) > 0)
    adjacency <- as.matrix(adjacency)
    adjacency[diagonal_index] <- 0
    if (graph_is_connected(adjacency) || k == J - 1L) {
      break
    }
    k <- k + 1L
  }

  edges <- which(upper.tri(adjacency) & adjacency > 0, arr.ind = TRUE)
  degree <- rowSums(adjacency)
  normalized <- diag(1 / sqrt(degree)) %*% adjacency %*% diag(1 / sqrt(degree))
  eigenvalues <- eigen(normalized, symmetric = TRUE, only.values = TRUE)$values
  list(
    adjacency = adjacency,
    k = k,
    degree = as.numeric(degree),
    N_edges = nrow(edges),
    node1 = as.integer(edges[, 1L]),
    node2 = as.integer(edges[, 2L]),
    normalized_adjacency_eigenvalues = as.numeric(eigenvalues),
    log_degree_sum = sum(log(degree))
  )
}

fit_geometry_model <- function(
    Y_geometry,
    vocabulary_geometry,
    compiled_ind_model,
    output_dir,
    seed,
    config) {
  scaling <- standardize_from_training(vocabulary_geometry)
  stan_data <- long_irt_data(
    Y = Y_geometry,
    vocabulary_z = scaling$training,
    grainsize = config$grainsize,
    sigma_gamma_prior = config$sigma_gamma_prior
  )
  fit <- sample_art_model(
    model = compiled_ind_model,
    model_name = "ind_2pl",
    stan_data = stan_data,
    output_dir = output_dir,
    seed = seed,
    config = config
  )
  list(fit = fit, scaling = scaling, stan_data = stan_data)
}

construct_item_geometry <- function(
    Y_geometry,
    vocabulary_geometry,
    compiled_ind_model,
    output_dir,
    seed,
    config,
    method = c("residual", "raw")) {
  method <- match.arg(method)
  geometry_fit <- fit_geometry_model(
    Y_geometry = Y_geometry,
    vocabulary_geometry = vocabulary_geometry,
    compiled_ind_model = compiled_ind_model,
    output_dir = file.path(output_dir, "nonspatial_geometry_fit"),
    seed = seed,
    config = config
  )
  probability <- posterior_expected_probabilities(
    fit = geometry_fit$fit,
    vocabulary_z = geometry_fit$scaling$training,
    max_draws = config$max_scoring_draws,
    seed = seed + 1L
  )
  residual <- standardized_irt_residuals(Y_geometry, probability)
  association <- if (method == "residual") {
    residual_association(
      residual,
      shrinkage = config$residual_correlation_shrinkage
    )
  } else {
    residual_association(
      scale(Y_geometry, center = TRUE, scale = FALSE),
      shrinkage = config$residual_correlation_shrinkage
    )
  }
  mds <- classical_mds_geometry(
    association,
    dimensions = config$mds_dimensions
  )
  graph <- build_car_graph(
    mds$coords,
    starting_neighbors = config$car_starting_neighbors
  )
  rownames(mds$coords) <- colnames(Y_geometry)
  colnames(residual) <- colnames(Y_geometry)

  write_csv_atomic(
    data.frame(item = rownames(mds$coords), mds$coords, row.names = NULL),
    file.path(output_dir, "item_coordinates.csv")
  )
  write_csv_atomic(
    data.frame(
      stress = mds$stress,
      distance_correlation = mds$distance_correlation,
      additive_constant = mds$additive_constant,
      gof_1 = mds$gof[1L],
      gof_2 = mds$gof[2L],
      car_neighbors = graph$k,
      car_edges = graph$N_edges,
      geometry_method = method
    ),
    file.path(output_dir, "geometry_diagnostics.csv")
  )
  save_rds_atomic(
    list(
      mds = mds,
      graph = graph,
      residual = residual,
      posterior_probability = probability,
      vocabulary_scaling = geometry_fit$scaling,
      method = method
    ),
    file.path(output_dir, "geometry.rds")
  )
  list(
    coords = mds$coords,
    graph = graph,
    residual = residual,
    probability = probability,
    mds = mds,
    geometry_fit = geometry_fit$fit,
    method = method
  )
}

orthogonal_procrustes <- function(candidate, reference) {
  candidate <- as.matrix(candidate)
  reference <- as.matrix(reference)
  reference_center <- colMeans(reference)
  candidate_centered <- scale(candidate, center = TRUE, scale = FALSE)
  reference_centered <- scale(reference, center = TRUE, scale = FALSE)
  decomposition <- svd(t(candidate_centered) %*% reference_centered)
  rotation <- decomposition$u %*% t(decomposition$v)
  aligned_centered <- candidate_centered %*% rotation
  sweep(aligned_centered, 2L, reference_center, "+")
}

bootstrap_mds_stability <- function(
    residual,
    reference_coords,
    shrinkage = 0.05,
    replicates = 200L,
    seed = 1L) {
  residual <- as.matrix(residual)
  reference_coords <- as.matrix(reference_coords)
  I <- nrow(residual)
  J <- ncol(residual)
  set.seed(seed)
  aligned <- array(NA_real_, dim = c(replicates, J, 2L))
  distance_correlation <- rep(NA_real_, replicates)
  procrustes_rmse <- rep(NA_real_, replicates)
  reference_distance <- as.vector(stats::dist(reference_coords))

  for (b in seq_len(replicates)) {
    index <- sample.int(I, I, replace = TRUE)
    association <- residual_association(residual[index, , drop = FALSE], shrinkage)
    candidate <- classical_mds_geometry(association, dimensions = 2L)$coords
    aligned_b <- orthogonal_procrustes(candidate, reference_coords)
    aligned[b, , ] <- aligned_b
    distance_correlation[b] <- suppressWarnings(stats::cor(
      reference_distance,
      as.vector(stats::dist(candidate)),
      method = "spearman"
    ))
    procrustes_rmse[b] <- sqrt(mean((aligned_b - reference_coords)^2))
  }

  item_stability <- data.frame(
    item = rownames(reference_coords) %||% seq_len(J),
    mds1_sd = apply(aligned[, , 1L, drop = FALSE], 2L, stats::sd),
    mds2_sd = apply(aligned[, , 2L, drop = FALSE], 2L, stats::sd)
  )
  summary <- data.frame(
    statistic = c("distance_correlation", "procrustes_rmse"),
    mean = c(mean(distance_correlation), mean(procrustes_rmse)),
    sd = c(stats::sd(distance_correlation), stats::sd(procrustes_rmse)),
    lower = c(
      stats::quantile(distance_correlation, 0.025),
      stats::quantile(procrustes_rmse, 0.025)
    ),
    upper = c(
      stats::quantile(distance_correlation, 0.975),
      stats::quantile(procrustes_rmse, 0.975)
    )
  )
  list(
    item = item_stability,
    summary = summary,
    distance_correlation = distance_correlation,
    procrustes_rmse = procrustes_rmse,
    aligned = aligned
  )
}
