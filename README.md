# SpatialGPIRT: Spatial Item Response Theory Models via Gaussian Processes

**SpatialGPIRT** fits item response models in which item difficulties are spatially dependent. Standard IRT assumes local item independence; in many assessments, items close together in geography or in meaning behave alike, and that dependence carries information. SpatialGPIRT places an anisotropic Gaussian process (GP) prior on item difficulty, so the data decide how far the dependence reaches and whether it differs by direction.

The current release implements a binary two-parameter logistic (2PL) model with anisotropic GP priors on item difficulty. Estimation runs on Hamiltonian Monte Carlo via [`cmdstanr`](https://mc-stan.org/cmdstanr/).

**Features**

- Anisotropic GP priors (Matérn family) on item difficulty over spatial coordinates
- Simulation utilities for generating spatially dependent response data
- Recovery metrics and posterior predictive accuracy out of the box
- Stan backend: full posterior inference, not point estimates

## Installation

SpatialGPIRT requires a working CmdStan installation. If you have not used `cmdstanr` before, install both:

```R
# install.packages("remotes")
remotes::install_github("mhuang233/SpatialGPIRT")

# one-time CmdStan setup, if needed
# install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
# cmdstanr::install_cmdstan()
```

## Toy Example

Run this to verify your installation (about a minute on a laptop):

```R
# Load the package
library(SpatialGPIRT)

# Number of respondents (I) and spatial items (J)
sim_data <- simulate_sgp_irt_binary(
  I = 50,
  J = 15,
  p = 1,   # covariates
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

Note: the toy settings (2 chains, 200/200 iterations) are for a quick installation check only. For real analyses, use longer chains and confirm convergence diagnostics.

## Roadmap

- Polytomous responses (graded-response extension, as developed in the paper)
- Latent similarity spaces for items without geographic coordinates
- Spatio-temporal extension (ST-IRT) for repeated assessments

## Citation

If you use SpatialGPIRT, please cite the paper:

> Huang, M., & Ghosh, S. (2025). Spatial Dependencies in Item Response Theory: Gaussian Process Priors for Geographic and Cognitive Measurement. arXiv:2507.09824. https://arxiv.org/abs/2507.09824

```bibtex
@misc{huang2025spatialgpirt,
  title  = {Spatial Dependencies in Item Response Theory: Gaussian Process
            Priors for Geographic and Cognitive Measurement},
  author = {Huang, Mingya and Ghosh, Soham},
  year   = {2025},
  eprint = {2507.09824},
  archivePrefix = {arXiv},
  url    = {https://arxiv.org/abs/2507.09824}
}
```

## Authors

[Mingya Huang](https://mhuang233.github.io) (University of Chicago) and Soham Ghosh (University of Wisconsin--Madison).

Issues and contributions are welcome via the [issue tracker](https://github.com/mhuang233/SpatialGPIRT/issues).
