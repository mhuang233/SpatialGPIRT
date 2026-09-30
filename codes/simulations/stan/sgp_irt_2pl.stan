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
  int<lower=1> N_train;
  int<lower=1> N_test;
  int<lower=2> I;
  int<lower=2> J;
  int<lower=1> p;

  array[N_train] int<lower=0, upper=1> y_train;
  array[N_train] int<lower=1, upper=I> subj_train;
  array[N_train] int<lower=1, upper=J> item_train;
  matrix[N_train, p] X_train;

  array[N_test] int<lower=0, upper=1> y_test;
  array[N_test] int<lower=1, upper=I> subj_test;
  array[N_test] int<lower=1, upper=J> item_test;
  matrix[N_test, p] X_test;

  matrix[J, J] sqdist1;
  matrix[J, J] sqdist2;
  int<lower=1> grainsize;

  real<lower=0> a_A;
  real<lower=0> b_A;
  real<lower=0> sigma_gamma_prior;
}

parameters {
  vector[I] z_theta;

  vector[J] z_log_alpha;
  real<lower=0> sigma_a;

  vector[J] z_f;
  vector[J] z_u;
  real<lower=0> sigma_f;
  real<lower=0> sigma_u;

  // A_pow[d] = A_d^2 for D = 2, where A_d = 1 / ell_d.
  vector<lower=0>[2] A_pow;
  vector[p] gamma;
}

transformed parameters {
  vector[I] theta_centered = z_theta - rep_vector(mean(z_theta), I);
  vector[I] theta =
    theta_centered / sqrt(dot_self(theta_centered) / I);

  vector[J] log_alpha = sigma_a * z_log_alpha;
  vector[J] alpha = exp(log_alpha);

  vector[2] inv_ell = sqrt(A_pow);
  vector[2] ell = 1.0 ./ inv_ell;
  matrix[J, J] R;
  matrix[J, J] L_R;
  vector[J] f;
  vector[J] u = sigma_u * z_u;
  vector[J] beta;
  real<lower=0, upper=1> nugget_ratio =
    square(sigma_u) /
    (square(sigma_f) + square(sigma_u) + 1e-12);

  for (j1 in 1:J) {
    for (j2 in j1:J) {
      real exponent_term = -0.5 * (
        sqdist1[j1, j2] * square(inv_ell[1]) +
        sqdist2[j1, j2] * square(inv_ell[2])
      );
      real r_val = exp(exponent_term);
      R[j1, j2] = r_val;
      R[j2, j1] = r_val;
    }
  }
  for (j in 1:J) {
    R[j, j] += 1e-8;
  }

  L_R = cholesky_decompose(R);
  f = sigma_f * (L_R * z_f);
  beta = f + u;
}

model {
  z_theta ~ std_normal();
  z_log_alpha ~ std_normal();
  sigma_a ~ normal(0, 1);

  z_f ~ std_normal();
  z_u ~ std_normal();
  sigma_f ~ normal(0, 1);
  sigma_u ~ normal(0, 1);
  A_pow ~ gamma(a_A, b_A);
  gamma ~ normal(0, sigma_gamma_prior);

  target += reduce_sum(partial_sum_lupmf,
                       y_train,
                       grainsize,
                       subj_train,
                       item_train,
                       X_train,
                       theta,
                       alpha,
                       beta,
                       gamma);
}

generated quantities {
  vector[N_train] log_lik;
  vector[N_test] log_lik_test;
  vector[N_test] p_test;
  array[N_test] int y_rep_test;
  real anisotropy_ratio = ell[1] / ell[2];

  for (n in 1:N_train) {
    real eta =
      alpha[item_train[n]] * theta[subj_train[n]]
      - beta[item_train[n]]
      + row(X_train, n) * gamma;
    log_lik[n] = bernoulli_logit_lpmf(y_train[n] | eta);
  }

  for (n in 1:N_test) {
    real eta =
      alpha[item_test[n]] * theta[subj_test[n]]
      - beta[item_test[n]]
      + row(X_test, n) * gamma;
    p_test[n] = inv_logit(eta);
    log_lik_test[n] = bernoulli_logit_lpmf(y_test[n] | eta);
    y_rep_test[n] = bernoulli_logit_rng(eta);
  }
}
