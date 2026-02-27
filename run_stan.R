#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(cmdstanr)
})

set.seed(20260226)

dir.create("samples", showWarnings = FALSE, recursive = TRUE)

S <- 20L
N <- 200L
P <- 3L
param_names <- c("shape", "scale", "shift")

# True hyperparameters on log scale
mu_true <- c(log(4.0), log(0.6), log(0.25))
sigma2_true <- c(0.04, 0.03, 0.02)
sigma_true <- sqrt(sigma2_true)

alpha_true <- matrix(0, nrow = S, ncol = P)
for (s in seq_len(S)) alpha_true[s, ] <- rnorm(P, mean = mu_true, sd = sigma_true)
colnames(alpha_true) <- paste0("eta_", param_names)

shape_true <- exp(alpha_true[, 1L])
scale_true <- exp(alpha_true[, 2L])
shift_true <- exp(alpha_true[, 3L])

y <- matrix(0, nrow = S, ncol = N)
for (s in seq_len(S)) y[s, ] <- shift_true[s] + rgamma(N, shape = shape_true[s], scale = scale_true[s])

# Priors: same structure as EMC2 cache prior setup (mu normal, sigma2 inverse-gamma)
m0 <- mu_true
s0 <- rep(0.25, P)
a0 <- rep(5, P)
b0 <- (a0 - 1) * sigma2_true

min_y <- apply(y, 1L, min)
init_fun <- function(chain_id = 1L) {
  mu_init <- rnorm(P, mean = m0, sd = sqrt(s0) * 0.3)
  sigma2_init <- pmax(0.01, rgamma(P, shape = a0, rate = (a0 - 1) / sigma2_true))
  z_init <- matrix(rnorm(S * P, 0, 0.2), nrow = S, ncol = P)
  alpha_shift_target <- log(pmax(1e-6, 0.75 * min_y))
  z_init[, 3L] <- (alpha_shift_target - mu_init[3L]) / sqrt(sigma2_init[3L])
  list(mu = mu_init, sigma2 = sigma2_init, z_alpha = z_init)
}

stan_data <- list(S = S, N = N, y = y, m0 = m0, s0 = s0, a0 = a0, b0 = b0)
stan_file <- file.path("stan", "hierarchical_shifted_gamma_logparam.stan")
exe_file <- file.path(tempdir(), "hierarchical_shifted_gamma_logparam")

model <- cmdstan_model(stan_file, exe_file = exe_file, quiet = TRUE)
fit <- model$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = min(4, parallel::detectCores(logical = TRUE)),
  iter_warmup = 400,
  iter_sampling = 400,
  seed = 123,
  init = init_fun,
  adapt_delta = 0.95,
  max_treedepth = 11,
  refresh = 100
)

out <- list(
  config = list(
    S = S, N = N, P = P,
    seed_data = 20260226,
    seed_stan = 123,
    chains = 4L,
    iter_warmup = 400L,
    iter_sampling = 400L,
    adapt_delta = 0.95,
    max_treedepth = 11L
  ),
  priors = list(m0 = m0, s0 = s0, a0 = a0, b0 = b0),
  truth = list(
    mu = setNames(mu_true, param_names),
    sigma2 = setNames(sigma2_true, param_names),
    alpha = alpha_true
  ),
  data = list(y = y),
  draws = list(
    mu = as.matrix(fit$draws(variables = "mu", format = "draws_matrix")),
    sigma2 = as.matrix(fit$draws(variables = "sigma2", format = "draws_matrix")),
    alpha = as.matrix(fit$draws(variables = "alpha_draw", format = "draws_matrix")),
    theta = as.matrix(fit$draws(variables = "theta_draw", format = "draws_matrix"))
  ),
  sampler_diagnostics = fit$sampler_diagnostics(format = "draws_matrix"),
  csv_files = fit$output_files()
)

saveRDS(out, file.path("samples", "stan_full_results.rds"))
cat("Saved: samples/stan_full_results.rds\n")
