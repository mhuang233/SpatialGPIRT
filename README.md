# SpatialGPIRT: Spatial Item Response Theory Models via Gaussian Processes

**SpatialGPIRT** provides a computational framework for simulating and fitting spatial Item Response Theory (IRT) models. Currently, it implements a binary 2-parameter logistic (2PL) model featuring anisotropic Gaussian Process (GP) priors to capture spatial dependencies in item difficulties. However, this can be easily extended to multiple outcomes.

Model estimation is performed via Hamiltonian Monte Carlo (HMC) using [`cmdstanr`](https://www.google.com/search?q=%5Bhttps://mc-stan.org/cmdstanr/%5D(https://mc-stan.org/cmdstanr/)).

## Installation

You can install the development version of SpatialGPIRT from GitHub using the `remotes` package:

```R
# install.packages("remotes")
remotes::install_github("YourUsername/SpatialGPIRT")

```


## Toy Example

Run this to verify your installation.

```R
# Load the package
library(SpatialGPIRT)

# Number of respondents (I) and spatial items (J) 
sim_data <- simulate_sgp_irt_binary(
  I = 50,    
  J = 15,    
  p = 1, # covariates
  seed = 42
)

fit_results <- fit_sgp_irt_binary(
  stan_data = sim_data$stan_data,
  chains = 2,              
  parallel_chains = 2,
  iter_warmup = 200,       
  iter_sampling = 200,     
  seed = 42
)


post_means <- posterior_means_from_fit(fit_results)
recovery_metrics <- compute_recovery_metrics(post_means, sim_data$truth)
print(recovery_metrics)

prob_hat <- posterior_mean_probabilities(fit_results, sim_data$truth, ndraws = 100)
pred_y <- ifelse(prob_hat > 0.5, 1L, 0L)
accuracy <- mean(pred_y == sim_data$truth$y)
cat("Toy model classification accuracy:", round(accuracy, 4), "\n")

```

## Reference

For full methodological details, formulations, and simulation studies, please see the associated pre-print:

> **[Spatial Item Response Theory Models via Gaussian Processes](https://arxiv.org/abs/2507.09824)**
