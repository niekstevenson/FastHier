data {
  int<lower=1> S;
  int<lower=1> N;
  array[S, N] real<lower=0> y;

  // Hyper-prior structure matching EMC2_cache.R style:
  // mu_p ~ Normal(m0_p, s0_p), sigma2_p ~ InvGamma(a0_p, b0_p)
  vector[3] m0;
  vector<lower=0>[3] s0;
  vector<lower=0>[3] a0;
  vector<lower=0>[3] b0;
}

transformed data {
  vector[S] min_y;
  for (s in 1:S) {
    real m = y[s, 1];
    for (n in 2:N) {
      if (y[s, n] < m) m = y[s, n];
    }
    min_y[s] = m;
  }
}

parameters {
  vector[3] mu;
  vector<lower=0>[3] sigma2;
  matrix[S, 3] z_alpha;
}

transformed parameters {
  vector<lower=0>[3] sigma = sqrt(sigma2);
  matrix[S, 3] alpha;

  for (s in 1:S) {
    for (p in 1:3) {
      alpha[s, p] = mu[p] + sigma[p] * z_alpha[s, p];
    }
  }
}

model {
  mu ~ normal(m0, sqrt(s0));
  sigma2 ~ inv_gamma(a0, b0);
  to_vector(z_alpha) ~ normal(0, 1);

  for (s in 1:S) {
    real shape = exp(alpha[s, 1]);
    real scale = exp(alpha[s, 2]);
    real shift = exp(alpha[s, 3]);

    if (shift >= min_y[s]) {
      target += negative_infinity();
    } else {
      for (n in 1:N) {
        target += gamma_lpdf(y[s, n] - shift | shape, inv(scale));
      }
    }
  }
}

generated quantities {
  matrix[S, 3] alpha_draw = alpha;
  matrix[S, 3] theta_draw;
  for (s in 1:S) {
    theta_draw[s, 1] = exp(alpha[s, 1]);                // shape
    theta_draw[s, 2] = exp(alpha[s, 2]);                // scale
    theta_draw[s, 3] = exp(alpha[s, 3]);                // shift
  }
}
