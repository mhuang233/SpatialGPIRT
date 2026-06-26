functions {
  real partial_sum_lpmf(array[] int y_slice,
                        int start,
                        int end,
                        array[] int subj_id,
                        array[] int item_id,
                        matrix X,
                        vector theta,
                        vector alpha,
                        vector beta,
                        vector gamma) {
    real lp = 0;
    int slice_size = end - start + 1;

    for (n in 1:slice_size) {
      int m = start + n - 1;
      int i = subj_id[m];
      int j = item_id[m];
      real eta = alpha[j] * theta[i] - beta[j] + row(X, m) * gamma;
      lp += bernoulli_logit_lupmf(y_slice[n] | eta);
    }
    return lp;
  }
}

data {
  int<lower=1> N;                       // number of observations
  int<lower=1> I;                       // number of respondents
  int<lower=1> J;                       // number of items
  array[N] int<lower=0, upper=1> y;     // binary responses
  array[N] int<lower=1, upper=I> subj_id;
  array[N] int<lower=1, upper=J> item_id;
  int<lower=1> p;                       // number of covariates
  matrix[N, p] X;                       // covariate matrix
  matrix[J, J] sqdist1;                 // squared distances, dimension 1
  matrix[J, J] sqdist2;                 // squared distances, dimension 2
  int<lower=1> grainsize;               // reduce_sum grain size
}

parameters {
  vector[I] theta_raw;
  real<lower=0> sigma_theta;

  vector[J] log_alpha_raw;
  real<lower=0> sigma_alpha;

  vector[J] z_beta;
  real<lower=0> sigma_gp;
  real<lower=0> sigma_nug;
  real<lower=0> ell1;
  real<lower=0> ell2;

  vector[p] gamma;
}

transformed parameters {
  vector[I] theta = sigma_theta * theta_raw;
  vector[J] alpha = exp(sigma_alpha * log_alpha_raw);
  matrix[J, J] K;
  matrix[J, J] L_K;
  vector[J] beta;

  for (j1 in 1:J) {
    for (j2 in j1:J) {
      real exponent_term = -0.5 * (
        sqdist1[j1, j2] / square(ell1) +
        sqdist2[j1, j2] / square(ell2)
      );
      real k_val = square(sigma_gp) * exp(exponent_term);
      K[j1, j2] = k_val;
      K[j2, j1] = k_val;
    }
  }

  for (j in 1:J) {
    K[j, j] += square(sigma_nug) + 1e-8;
  }

  L_K = cholesky_decompose(K);
  beta = L_K * z_beta;
  beta -= mean(beta);                   // soft identification for location
}

model {
  theta_raw ~ std_normal();
  sigma_theta ~ normal(0, 1);

  log_alpha_raw ~ std_normal();
  sigma_alpha ~ normal(0, 0.5);

  z_beta ~ std_normal();
  sigma_gp ~ normal(0, 1);
  sigma_nug ~ normal(0, 0.5);
  ell1 ~ lognormal(log(0.30), 0.50);
  ell2 ~ lognormal(log(0.30), 0.50);

  gamma ~ normal(0, 1);

  target += reduce_sum(partial_sum_lupmf,
                       y,
                       grainsize,
                       subj_id,
                       item_id,
                       X,
                       theta,
                       alpha,
                       beta,
                       gamma);
}
