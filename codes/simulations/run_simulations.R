# Production runner for the SGP-IRT simulation study.

if (!exists("simulate_irt_dgp", mode = "function")) {
  source(file.path("R", "simulation_framework.R"))
}

run_simulation_study <- function(
    n_rep = 50L,
    scenarios = scenario_grid()$scenario,
    models = c("sgp", "car", "ind"),
    I = 200L,
    J = 50L,
    p = 1L,
    sigma_a = 0.35,
    gamma = 0.50,
    test_fraction = 0.20,
    base_seed = 202600L,
    stan_dir = "stan",
    output_dir = "results/production",
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
    save_first_fit = TRUE,
    skip_completed = TRUE,
    force_recompile = FALSE) {
  required_simulation_packages()
  scenarios <- match.arg(
    scenarios,
    scenario_grid()$scenario,
    several.ok = TRUE
  )
  models <- match.arg(
    models,
    c("sgp", "car", "ind"),
    several.ok = TRUE
  )
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  model_objects <- compile_matched_models(
    stan_dir = stan_dir,
    models = models,
    force_recompile = force_recompile
  )
  
  manifest <- expand.grid(
    scenario = scenarios,
    replication = seq_len(n_rep),
    stringsAsFactors = FALSE
  )
  manifest$I <- I
  manifest$J <- J
  manifest$p <- p
  manifest$sigma_a <- sigma_a
  manifest$test_fraction <- test_fraction
  manifest$gamma <- paste(gamma, collapse = ";")
  manifest$models <- paste(models, collapse = ";")
  manifest$chains <- chains
  manifest$iter_warmup <- iter_warmup
  manifest$iter_sampling <- iter_sampling
  manifest$adapt_delta <- adapt_delta
  manifest$max_treedepth <- max_treedepth
  manifest$inverse_length_shape <- a_A
  manifest$inverse_length_rate <- b_A
  manifest$sigma_gamma_prior <- sigma_gamma_prior
  manifest$data_seed <- base_seed +
    10000L * match(manifest$scenario, scenario_grid()$scenario) +
    manifest$replication
  write_csv_atomic(manifest, file.path(output_dir, "simulation_manifest.csv"))
  software <- c(
    capture.output(sessionInfo()),
    paste(
      "CmdStan version:",
      paste(cmdstanr::cmdstan_version(), collapse = ".")
    )
  )
  writeLines(software, file.path(output_dir, "software_versions.txt"))
  
  for (row in seq_len(nrow(manifest))) {
    scenario <- manifest$scenario[row]
    replication <- manifest$replication[row]
    data_seed <- manifest$data_seed[row]
    
    message(
      "Scenario ", scenario,
      ", replication ", replication, " of ", n_rep
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
    
    for (model_index in seq_along(models)) {
      model <- models[model_index]
      model_dir <- file.path(
        output_dir,
        scenario,
        sprintf("rep_%03d", replication)
      )
      diagnostic_file <- file.path(
        model_dir,
        paste0(model, "_diagnostics.csv")
      )
      parameter_file <- file.path(
        model_dir,
        paste0(model, "_parameter_metrics.csv")
      )
      predictive_file <- file.path(
        model_dir,
        paste0(model, "_predictive_metrics.csv")
      )
      completed_fit <- FALSE
      if (
        file.exists(diagnostic_file) &&
        file.exists(parameter_file) &&
        file.exists(predictive_file)
      ) {
        completed_fit <- tryCatch(
          {
            diagnostic_record <- utils::read.csv(
              diagnostic_file,
              stringsAsFactors = FALSE
            )
            nrow(diagnostic_record) == 1L &&
              isTRUE(as.logical(diagnostic_record$fit_completed[1L]))
          },
          error = function(e) FALSE
        )
      }
      if (isTRUE(skip_completed) && completed_fit) {
        message("  Skipping completed model ", model)
        next
      }
      
      fit_seed <- data_seed + 1000L * model_index
      message("  Fitting ", toupper(model), "-IRT")
      run_one_model(
        sim = sim,
        model = model,
        model_object = model_objects[[model]],
        scenario = scenario,
        replication = replication,
        seed = fit_seed,
        result_dir = output_dir,
        chains = chains,
        parallel_chains = parallel_chains,
        threads_per_chain = threads_per_chain,
        iter_warmup = iter_warmup,
        iter_sampling = iter_sampling,
        adapt_delta = adapt_delta,
        max_treedepth = max_treedepth,
        a_A = a_A,
        b_A = b_A,
        sigma_gamma_prior = sigma_gamma_prior,
        save_fit = isTRUE(save_first_fit) && replication == 1L
      )
    }
  }
  invisible(manifest)
}

run_single_simulation_job <- function(
    scenario,
    replication,
    model,
    I = 200L,
    J = 50L,
    p = 1L,
    sigma_a = 0.35,
    gamma = 0.50,
    test_fraction = 0.20,
    base_seed = 202600L,
    stan_dir = "stan",
    output_dir = "results/production",
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
  scenario <- match.arg(scenario, scenario_grid()$scenario)
  model <- match.arg(model, c("sgp", "car", "ind"))
  scenario_index <- match(scenario, scenario_grid()$scenario)
  data_seed <- base_seed + 10000L * scenario_index + replication
  fit_seed <- data_seed + 1000L * match(model, c("sgp", "car", "ind"))
  
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
  model_object <- compile_matched_models(
    stan_dir = stan_dir,
    models = model
  )[[model]]
  run_one_model(
    sim = sim,
    model = model,
    model_object = model_object,
    scenario = scenario,
    replication = replication,
    seed = fit_seed,
    result_dir = output_dir,
    chains = chains,
    parallel_chains = parallel_chains,
    threads_per_chain = threads_per_chain,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    a_A = a_A,
    b_A = b_A,
    sigma_gamma_prior = sigma_gamma_prior,
    save_fit = save_fit
  )
}

run_scaling_study <- function(
    J_values = c(25L, 50L, 100L, 200L),
    n_rep = 20L,
    scenarios = c("gp_aniso_low_nugget", "car_truth"),
    models = c("sgp", "car", "ind"),
    output_root = "results/scaling",
    ...) {
  for (J in J_values) {
    run_simulation_study(
      n_rep = n_rep,
      scenarios = scenarios,
      models = models,
      J = J,
      output_dir = file.path(output_root, paste0("J_", J)),
      ...
    )
  }
  invisible(TRUE)
}

run_discrimination_recovery_study <- function(
    I_values = c(100L, 200L, 500L, 1000L),
    n_rep = 20L,
    scenario = "gp_aniso_low_nugget",
    models = c("sgp", "car", "ind"),
    output_root = "results/discrimination_scaling",
    ...) {
  for (I in I_values) {
    run_simulation_study(
      n_rep = n_rep,
      scenarios = scenario,
      models = models,
      I = I,
      output_dir = file.path(output_root, paste0("I_", I)),
      ...
    )
  }
  invisible(TRUE)
}
