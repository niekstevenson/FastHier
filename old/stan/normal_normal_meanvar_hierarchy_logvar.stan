data {
  int<lower=1> S;
  int<lower=1> N;
  array[S, N] real y;

  real prior_mu_mean;
  real<lower=0> prior_sd_mu_mean;

  real prior_mu_log_var;
  real<lower=0> prior_sd_mu_log_var;

  real prior_mean_log_tau2_mean;
  real<lower=0> prior_sd_log_tau2_mean;

  real prior_mean_log_tau2_log_var;
  real<lower=0> prior_sd_log_tau2_log_var;
}

parameters {
  real mu_mean;
  real mu_log_var;
  real log_tau2_mean;
  real log_tau2_log_var;
  vector[S] z_subject_mean;
  vector[S] z_subject_log_var;
}

transformed parameters {
  real<lower=0> tau_mean = exp(0.5 * log_tau2_mean);
  real<lower=0> tau_log_var = exp(0.5 * log_tau2_log_var);
  vector[S] subject_mean = mu_mean + tau_mean * z_subject_mean;
  vector[S] subject_log_var = mu_log_var + tau_log_var * z_subject_log_var;
}

model {
  mu_mean ~ normal(prior_mu_mean, prior_sd_mu_mean);
  mu_log_var ~ normal(prior_mu_log_var, prior_sd_mu_log_var);
  log_tau2_mean ~ normal(prior_mean_log_tau2_mean, prior_sd_log_tau2_mean);
  log_tau2_log_var ~ normal(prior_mean_log_tau2_log_var, prior_sd_log_tau2_log_var);

  z_subject_mean ~ normal(0, 1);
  z_subject_log_var ~ normal(0, 1);

  for (s in 1:S) {
    y[s] ~ normal(subject_mean[s], fmax(exp(0.5 * subject_log_var[s]), 1e-9));
  }
}

generated quantities {
  real<lower=0> tau2_mean = exp(log_tau2_mean);
  real<lower=0> tau2_log_var = exp(log_tau2_log_var);
  vector[S] subject_variance = exp(subject_log_var);
}
