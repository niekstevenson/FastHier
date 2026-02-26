data {
  int<lower=1> N;
  array[N] real<lower=0> y;
}

transformed data {
  real min_y = min(y);
}

parameters {
  real eta_shape;
  real eta_scale;
  real eta_shift;
}

transformed parameters {
  real<lower=0> shape = exp(eta_shape) + 1e-9;
  real<lower=0> scale = exp(eta_scale) + 1e-9;
  real<lower=0> shift = exp(eta_shift) + 1e-9;
}

model {
  eta_shape ~ normal(0, 1);
  eta_scale ~ normal(0, 1);
  eta_shift ~ normal(0, 1);

  if (shift >= min_y) {
    target += negative_infinity();
  } else {
    for (n in 1:N) {
      target += gamma_lpdf(y[n] - shift | shape, inv(scale));
    }
  }
}
