art_item_columns <- function(data, genre_prefix) {
  columns <- grep(paste0("^", genre_prefix, "[0-9]+$"), names(data), value = TRUE)
  item_number <- as.integer(sub(paste0("^", genre_prefix), "", columns))
  columns[order(item_number)]
}

coerce_binary_column <- function(x, name) {
  y <- suppressWarnings(as.integer(as.character(x)))
  if (anyNA(y) || !all(y %in% c(0L, 1L))) {
    stop("Item column ", name, " is not completely coded 0/1.", call. = FALSE)
  }
  y
}

prepare_art_data <- function(
    path,
    genres,
    vocabulary_column = "Vocab_Sum",
    cohort_rule = c("manuscript_321", "filter_flag", "all_complete"),
    audit_dir = NULL) {
  cohort_rule <- match.arg(cohort_rule)
  assert_true(file.exists(path), paste0("ART data file not found: ", path))

  raw <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  raw$.source_row <- seq_len(nrow(raw))

  required_genre_columns <- lapply(
    unname(genres),
    function(prefix) art_item_columns(raw, prefix)
  )
  names(required_genre_columns) <- names(genres)
  required_items <- unlist(required_genre_columns, use.names = FALSE)
  assert_true(length(required_items) > 0L, "No prespecified ART columns were found.")
  assert_true(
    vocabulary_column %in% names(raw),
    paste0("Missing vocabulary column: ", vocabulary_column)
  )

  for (column in required_items) {
    raw[[column]] <- coerce_binary_column(raw[[column]], column)
  }
  raw[[vocabulary_column]] <- suppressWarnings(as.numeric(raw[[vocabulary_column]]))

  attention_pass <- if ("AttentionCheck_Failed" %in% names(raw)) {
    suppressWarnings(as.numeric(raw$AttentionCheck_Failed)) == 0
  } else {
    rep(TRUE, nrow(raw))
  }
  age_observed <- if ("Age" %in% names(raw)) {
    !is.na(suppressWarnings(as.numeric(raw$Age)))
  } else {
    rep(TRUE, nrow(raw))
  }
  filter_pass <- if ("filter_$" %in% names(raw)) {
    suppressWarnings(as.numeric(raw[["filter_$"]])) == 1
  } else {
    rep(TRUE, nrow(raw))
  }
  complete_required <- stats::complete.cases(
    raw[, c(vocabulary_column, required_items), drop = FALSE]
  )

  keep <- switch(
    cohort_rule,
    manuscript_321 = attention_pass & age_observed & complete_required,
    filter_flag = filter_pass & complete_required,
    all_complete = complete_required
  )

  flow <- data.frame(
    criterion = c(
      "Rows in supplied file",
      "Pass attention check",
      "Age observed",
      "Complete vocabulary and target items",
      paste0("Final cohort: ", cohort_rule)
    ),
    n = c(
      nrow(raw),
      sum(attention_pass),
      sum(age_observed),
      sum(complete_required),
      sum(keep)
    )
  )

  if (cohort_rule == "manuscript_321" && nrow(raw) == 337L) {
    assert_true(
      sum(keep) == 321L,
      "The supplied data no longer reproduce the manuscript's N = 321 cohort."
    )
  }

  data <- raw[keep, , drop = FALSE]
  rownames(data) <- NULL
  data$subject_id <- seq_len(nrow(data))

  audit <- do.call(rbind, lapply(names(genres), function(genre) {
    columns <- required_genre_columns[[genre]]
    rates <- colMeans(data[, columns, drop = FALSE])
    data.frame(
      genre = genre,
      item = columns,
      n_respondents = nrow(data),
      endorsement_rate = unname(rates),
      zero_variance = rates %in% c(0, 1)
    )
  }))
  assert_true(!any(audit$zero_variance), "At least one selected ART item has zero variance.")
  assert_true(all(is.finite(data[[vocabulary_column]])), "Vocabulary has missing/nonfinite values.")
  assert_true(stats::sd(data[[vocabulary_column]]) > 0, "Vocabulary has zero variance.")

  if (!is.null(audit_dir)) {
    ensure_dir(audit_dir)
    write_csv_atomic(flow, file.path(audit_dir, "cohort_flow.csv"))
    write_csv_atomic(audit, file.path(audit_dir, "item_audit.csv"))
    writeLines(
      c(
        paste0("Data file: ", normalizePath(path, winslash = "/", mustWork = TRUE)),
        paste0("Cohort rule: ", cohort_rule),
        paste0("Final N: ", nrow(data)),
        paste0("Vocabulary: ", vocabulary_column),
        paste0("Genres: ", paste(names(genres), collapse = ", "))
      ),
      file.path(audit_dir, "analysis_sample.txt")
    )
  }

  list(
    data = data,
    flow = flow,
    item_audit = audit,
    item_columns = required_genre_columns,
    vocabulary_column = vocabulary_column,
    genres = genres,
    cohort_rule = cohort_rule
  )
}

genre_response_matrix <- function(prepared, genre) {
  assert_true(genre %in% names(prepared$item_columns), paste0("Unknown genre: ", genre))
  columns <- prepared$item_columns[[genre]]
  y <- as.matrix(prepared$data[, columns, drop = FALSE])
  storage.mode(y) <- "integer"
  rownames(y) <- prepared$data$.source_row
  colnames(y) <- columns
  y
}

all_genre_total_score <- function(prepared) {
  all_items <- unlist(prepared$item_columns, use.names = FALSE)
  rowSums(as.matrix(prepared$data[, all_items, drop = FALSE]))
}

standardize_from_training <- function(x_training, x_new = NULL) {
  x_training <- as.numeric(x_training)
  center <- mean(x_training)
  scale <- stats::sd(x_training)
  assert_true(is.finite(scale) && scale > 0, "Training vocabulary SD is not positive.")
  result <- list(
    training = as.numeric((x_training - center) / scale),
    center = center,
    scale = scale
  )
  if (!is.null(x_new)) {
    result$new <- as.numeric((as.numeric(x_new) - center) / scale)
  }
  result
}

make_three_way_splits <- function(
    prepared,
    n_folds = 5L,
    geometry_fraction = 0.5,
    seed = 20260730L) {
  n <- nrow(prepared$data)
  vocabulary <- prepared$data[[prepared$vocabulary_column]]
  # Split assignment uses only the external vocabulary covariate. Test ART
  # responses are not inspected to balance the folds.
  strata <- make_balance_strata(vocabulary)
  outer_fold <- stratified_folds(strata, n_folds, seed)

  splits <- vector("list", n_folds)
  for (fold in seq_len(n_folds)) {
    test <- which(outer_fold == fold)
    outer_training <- setdiff(seq_len(n), test)
    local_strata <- droplevels(strata[outer_training])
    set.seed(seed + 1000L * fold)
    geometry <- integer(0)
    for (level in levels(local_strata)) {
      candidates <- outer_training[which(local_strata == level)]
      n_geometry <- max(
        1L,
        min(length(candidates) - 1L, round(geometry_fraction * length(candidates)))
      )
      if (length(candidates) == 1L) {
        n_geometry <- 0L
      }
      if (n_geometry > 0L) {
        geometry <- c(geometry, sample(candidates, n_geometry))
      }
    }
    calibration <- setdiff(outer_training, geometry)
    assert_true(length(intersect(geometry, calibration)) == 0L, "Split overlap.")
    assert_true(length(intersect(geometry, test)) == 0L, "Split overlap.")
    assert_true(length(intersect(calibration, test)) == 0L, "Split overlap.")
    assert_true(
      length(sort(c(geometry, calibration, test))) == n,
      "Three-way split does not partition all respondents."
    )
    splits[[fold]] <- list(
      fold = fold,
      geometry = sort(geometry),
      calibration = sort(calibration),
      test = sort(test)
    )
  }
  splits
}
