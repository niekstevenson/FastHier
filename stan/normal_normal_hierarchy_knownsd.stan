data {
  int<lower=1> S;
  int<lower=1> N;
  array[S, N] real y;
  real<lower=0> sigma_y;

  real m0;
  real<lower=0> s0;
  real<lower=0> a0;
  real<lower=0> b0;
}

parameters {
  real mu;
  real<lower=0> tau2;
  vector[S] z_alpha;
}

transformed parameters {
  real<lower=0> tau = sqrt(tau2);
  vector[S] alpha = mu + tau * z_alpha;
}

model {
  mu ~ normal(m0, sqrt(s0));
  tau2 ~ inv_gamma(a0, b0);
  z_alpha ~ normal(0, 1);

  for (s in 1:S) {
    y[s] ~ normal(alpha[s], sigma_y);
  }
}

generated quantities {
  vector[S] alpha_draw = alpha;
}
