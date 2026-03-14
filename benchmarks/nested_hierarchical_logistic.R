#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(cmdstanr)
})

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[1L]))
} else {
  normalizePath("benchmarks/nested_hierarchical_logistic.R")
}
repo_dir <- dirname(dirname(script_path))
setwd(repo_dir)

source("smc_core.R")
source("SMC_super_fast.R")
source("new_SMC_cache.R")
source("nested_population_SMC.R")
source("make_prior.R")

dir.create("benchmarks", showWarnings = FALSE, recursive = TRUE)
dir.create("samples", showWarnings = FALSE, recursive = TRUE)

set.seed(20260311L)

S <- 20L
N_per_group <- 150L
P <- 3L
theta_names <- c("beta0", "beta1", "beta2")
mu_names <- c("mu_beta0", "mu_beta1", "mu_beta2")
sigma2_names <- c("sigma2_beta0", "sigma2_beta1", "sigma2_beta2")

mu_true <- c(-0.4, 0.9, -0.8)
sigma2_true <- c(0.30, 0.18, 0.16)
sigma_true <- sqrt(sigma2_true)

prior_m0 <- c(0.0, 0.0, 0.0)
prior_s0 <- c(1.0, 1.0, 1.0)
prior_a0 <- c(5.0, 5.0, 5.0)
prior_b0 <- (prior_a0 - 1.0) * sigma2_true

alpha_true <- matrix(0, nrow = S, ncol = P)
for (i in seq_len(S)) {
  alpha_true[i, ] <- rnorm(P, mean = mu_true, sd = sigma_true)
}
colnames(alpha_true) <- theta_names

X_array <- array(0, dim = c(S, N_per_group, P))
y_mat <- matrix(0L, nrow = S, ncol = N_per_group)
data_list <- vector("list", length = S)

for (i in seq_len(S)) {
  x1 <- rnorm(N_per_group)
  x2 <- 0.35 * x1 + sqrt(1 - 0.35^2) * rnorm(N_per_group)
  X_i <- cbind(1.0, x1, x2)
  eta_i <- as.numeric(X_i %*% alpha_true[i, ])
  p_i <- plogis(eta_i)
  y_i <- rbinom(N_per_group, size = 1L, prob = p_i)
  X_array[i, , ] <- X_i
  y_mat[i, ] <- y_i
  data_list[[i]] <- list(X = X_i, y = y_i)
}

loglik_hier_logistic <- function(Theta, data_i) {
  Theta <- as.matrix(Theta)
  X <- as.matrix(data_i$X)
  y <- as.numeric(data_i$y)
  eta <- Theta %*% t(X)
  softplus <- pmax(eta, 0) + log1p(exp(-abs(eta)))
  as.numeric(drop(eta %*% y) - rowSums(softplus))
}

gaussian_map_fn <- phi_to_gaussian_params_diag_factory()

alpha_ref_mean <- setNames(prior_m0, theta_names)
alpha_ref_var <- prior_s0 + prior_b0 / (prior_a0 - 1.0)
Sigma_ref <- diag(alpha_ref_var, P)
colnames(Sigma_ref) <- rownames(Sigma_ref) <- theta_names
Sigma_ref_inv <- chol2inv(chol(Sigma_ref))
Sigma_ref_logdet <- as.numeric(determinant(Sigma_ref, logarithm = TRUE)$modulus)

local_particles <- 800L
local_rounds <- 180L
population_particles <- 1000L
local_bank_particles <- 64L
population_rounds <- 50L
mc_cores <- 1L

stan_results_file <- file.path("samples", "stan_hierarchical_logistic_results.rds")

if (file.exists(stan_results_file)) {
  stan_res <- readRDS(stan_results_file)
} else {
  stan_code <- "
data {
  int<lower=1> S;
  int<lower=1> N_total;
  int<lower=1> K;
  array[N_total] int<lower=1, upper=S> g;
  matrix[N_total, K] X;
  array[N_total] int<lower=0, upper=1> y;
  vector[K] m0;
  vector<lower=0>[K] s0;
  vector<lower=0>[K] a0;
  vector<lower=0>[K] b0;
}

parameters {
  vector[K] mu;
  vector<lower=0>[K] sigma2;
  matrix[S, K] z_alpha;
}

transformed parameters {
  matrix[S, K] alpha;
  for (s in 1:S) {
    alpha[s] = mu' + sqrt(sigma2)' .* z_alpha[s];
  }
}

model {
  for (k in 1:K) {
    mu[k] ~ normal(m0[k], sqrt(s0[k]));
    sigma2[k] ~ inv_gamma(a0[k], b0[k]);
  }
  to_vector(z_alpha) ~ normal(0, 1);

  for (n in 1:N_total) {
    y[n] ~ bernoulli_logit(dot_product(to_vector(X[n]), to_vector(alpha[g[n]])));
  }
}
"

  stan_file <- file.path(tempdir(), "hierarchical_logistic_diag.stan")
  writeLines(stan_code, stan_file)

  X_long <- do.call(rbind, lapply(seq_len(S), function(i) data_list[[i]]$X))
  y_long <- unlist(lapply(data_list, `[[`, "y"), use.names = FALSE)
  g_long <- rep(seq_len(S), each = N_per_group)

  stan_model <- cmdstan_model(
    stan_file,
    exe_file = file.path(tempdir(), "hierarchical_logistic_diag"),
    quiet = TRUE
  )

  stan_fit <- stan_model$sample(
    data = list(
      S = S,
      N_total = length(y_long),
      K = P,
      g = g_long,
      X = X_long,
      y = y_long,
      m0 = prior_m0,
      s0 = prior_s0,
      a0 = prior_a0,
      b0 = prior_b0
    ),
    chains = 4,
    parallel_chains = min(4L, parallel::detectCores(logical = TRUE)),
    iter_warmup = 400,
    iter_sampling = 400,
    seed = 123L,
    adapt_delta = 0.95,
    max_treedepth = 11,
    refresh = 100
  )

  stan_res <- list(
    config = list(
      S = S,
      N_per_group = N_per_group,
      P = P,
      seed_data = 20260311L,
      seed_stan = 123L
    ),
    priors = list(
      m0 = prior_m0,
      s0 = prior_s0,
      a0 = prior_a0,
      b0 = prior_b0
    ),
    truth = list(
      mu = setNames(mu_true, theta_names),
      sigma2 = setNames(sigma2_true, theta_names),
      alpha = alpha_true
    ),
    data = list(
      X = X_array,
      y = y_mat
    ),
    draws = list(
      mu = as.matrix(stan_fit$draws(variables = "mu", format = "draws_matrix")),
      sigma2 = as.matrix(stan_fit$draws(variables = "sigma2", format = "draws_matrix"))
    )
  )

  saveRDS(stan_res, stan_results_file)
}

ll_eval_counter_reset()
local_fits <- parallel::mclapply(
  seq_len(S),
  function(i) {
    out <- enhanced_smc_elite(
      data = data_list[[i]],
      loglik_fn = loglik_hier_logistic,
      mu_ref = alpha_ref_mean,
      Sigma_ref = Sigma_ref,
      M = local_particles,
      resample_threshold = 0.6,
      n_mcmc_moves = 3L,
      max_rounds = local_rounds,
      G_mix = 8L,
      da_enable = TRUE,
      gss_enable = TRUE,
      ll_cache_enable = TRUE,
      deterministic_resampling = FALSE,
      n_cores = 1L,
      seed = as.integer(1123L + 31L * i),
      verbose = FALSE
    )
    out$working_prior <- list(
      mu = as.numeric(alpha_ref_mean),
      Sigma = Sigma_ref,
      Sigma_inv = Sigma_ref_inv,
      logdet = Sigma_ref_logdet
    )
    out
  },
  mc.cores = mc_cores
)
local_prefit_loglik <- ll_eval_counter_get()
ll_eval_counter_disable()

local_objs <- lapply(
  seq_len(S),
  function(i) {
    build_local_exact_object(
      smc_out = local_fits[[i]],
      data = data_list[[i]],
      subj_id = i,
      loglik_fn = loglik_hier_logistic,
      base_seed = as.integer(40000L + 1009L * i)
    )
  }
)

prior <- make_prior_phi_diag(
  m0 = prior_m0,
  s0 = prior_s0,
  a = prior_a0,
  b = prior_b0,
  d = P
)

nested_fit <- nested_population_smc(
  rprior_phi = prior$rprior,
  local_objs = local_objs,
  logprior_phi = prior$lprior,
  gaussian_map_fn = gaussian_map_fn,
  N = population_particles,
  M_local = local_bank_particles,
  rho_res = 0.5,
  rho_local = 0.5,
  n_population_moves = 1L,
  n_local_moves = 0L,
  n_population_refresh_moves = 0L,
  max_rounds = population_rounds,
  max_bank_topups = 3L,
  bank_split_tol = 0.01,
  enable_exact_local_refresh = TRUE,
  max_exact_local_refreshes = 1L,
  exact_refresh_particles = 400L,
  exact_refresh_max_rounds = 120L,
  exact_refresh_n_mcmc_moves = 3L,
  seed = 123L,
  verbose = TRUE
)

stan_mu <- as.matrix(stan_res$draws$mu)
stan_sigma2 <- as.matrix(stan_res$draws$sigma2)
colnames(stan_mu) <- mu_names
colnames(stan_sigma2) <- sigma2_names

nested_phi <- as.matrix(nested_fit$phi)
nested_w <- pmax(as.numeric(nested_fit$w), 0)
nested_w <- nested_w / sum(nested_w)
nested_mu <- nested_phi[, seq_len(P), drop = FALSE]
nested_sigma2 <- exp(nested_phi[, P + seq_len(P), drop = FALSE])
colnames(nested_mu) <- mu_names
colnames(nested_sigma2) <- sigma2_names

set.seed(123L)
draw_idx <- sample.int(nrow(nested_phi), size = nrow(stan_mu), replace = TRUE, prob = nested_w)
nested_mu_draws <- nested_mu[draw_idx, , drop = FALSE]
nested_sigma2_draws <- nested_sigma2[draw_idx, , drop = FALSE]

plot_file <- file.path("benchmarks", "nested_hierarchical_logistic_posteriors.png")
png(plot_file, width = 1200, height = 800, res = 120)
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))

for (j in seq_len(P)) {
  plot(
    density(stan_mu[, j]),
    main = mu_names[j],
    xlab = mu_names[j],
    ylab = "Density",
    lwd = 2,
    col = "black"
  )
  lines(density(nested_mu_draws[, j]), lwd = 2, col = "red3")
  abline(v = mu_true[j], col = "grey60", lty = 3)
  legend("topright", legend = c("Stan", "Nested SMC", "Truth"), col = c("black", "red3", "grey60"), lwd = c(2, 2, 1), lty = c(1, 1, 3), bty = "n")
}

for (j in seq_len(P)) {
  plot(
    density(stan_sigma2[, j]),
    main = sigma2_names[j],
    xlab = sigma2_names[j],
    ylab = "Density",
    lwd = 2,
    col = "black"
  )
  lines(density(nested_sigma2_draws[, j]), lwd = 2, col = "red3")
  abline(v = sigma2_true[j], col = "grey60", lty = 3)
  legend("topright", legend = c("Stan", "Nested SMC", "Truth"), col = c("black", "red3", "grey60"), lwd = c(2, 2, 1), lty = c(1, 1, 3), bty = "n")
}

dev.off()

results <- list(
  stan_results_file = stan_results_file,
  local_fits = local_fits,
  nested_fit = nested_fit,
  plot_file = plot_file,
  settings = list(
    local_particles = local_particles,
    local_rounds = local_rounds,
    population_particles = population_particles,
    local_bank_particles = local_bank_particles,
    population_rounds = population_rounds,
    alpha_ref_mean = unname(alpha_ref_mean),
    alpha_ref_var = unname(alpha_ref_var),
    max_bank_topups = 3L,
    bank_split_tol = 0.01,
    enable_exact_local_refresh = TRUE,
    max_exact_local_refreshes = 1L,
    exact_refresh_particles = 400L,
    exact_refresh_max_rounds = 120L,
    exact_refresh_n_mcmc_moves = 3L,
    mc_cores = mc_cores
  )
)

outfile <- file.path("samples", "nested_hierarchical_logistic_results.rds")
saveRDS(results, outfile)

eval_counts <- nested_fit$meta$local_loglik_evals
local_end_to_end <- as.integer(local_prefit_loglik + eval_counts["total"])

stan_mu_mean <- colMeans(stan_mu)
stan_sigma2_mean <- colMeans(stan_sigma2)
nested_mu_mean <- colSums(nested_mu * nested_w)
nested_sigma2_mean <- colSums(nested_sigma2 * nested_w)

cat("NESTED_HIERARCHICAL_LOGISTIC_OK\n")
cat(sprintf("stan_results=%s\n", stan_results_file))
cat(sprintf("results=%s\n", outfile))
cat(sprintf("plot=%s\n", plot_file))
cat(sprintf("population_final_lambda=%.4f\n", tail(nested_fit$meta$lambda_hist, 1L)))
cat(sprintf("population_rounds=%d\n", nested_fit$meta$rounds))
cat(sprintf("implementation=%s\n", nested_fit$meta$implementation))
cat(sprintf("bank_topup_points=%d\n", nested_fit$meta$bank_topup_points))
cat(sprintf("n_exact_local_refresh=%d\n", nested_fit$meta$n_exact_local_refresh))
cat(sprintf("bank_mean_size_final=%.1f\n", tail(nested_fit$meta$bank_size_mean_hist, 1L)))
cat(sprintf("local_loglik_prefit=%d\n", local_prefit_loglik))
cat(sprintf("local_loglik_nested_total=%d\n", eval_counts["total"]))
cat(sprintf("local_loglik_nested_init=%d\n", eval_counts["initialization"]))
cat(sprintf("local_loglik_nested_enrichment=%d\n", eval_counts["enrichment"]))
cat(sprintf("local_loglik_nested_population_refresh=%d\n", eval_counts["population_refreshes"]))
cat(sprintf("local_loglik_end_to_end=%d\n", local_end_to_end))
for (j in seq_len(P)) {
  cat(sprintf("%s_mean_stan=%.6f %s_mean_nested=%.6f truth=%.6f\n", mu_names[j], stan_mu_mean[j], mu_names[j], nested_mu_mean[j], mu_true[j]))
}
for (j in seq_len(P)) {
  cat(sprintf("%s_mean_stan=%.6f %s_mean_nested=%.6f truth=%.6f\n", sigma2_names[j], stan_sigma2_mean[j], sigma2_names[j], nested_sigma2_mean[j], sigma2_true[j]))
}
