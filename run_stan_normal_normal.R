#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(cmdstanr))
suppressPackageStartupMessages(library(parallel))

set.seed(20260310L)

S <- 20L
N <- 100L
sigma_y <- 1.0
mu_true <- 1.25
tau2_true <- 0.40

alpha_true <- rnorm(S, mean = mu_true, sd = sqrt(tau2_true))
y <- matrix(rnorm(S * N, mean = rep(alpha_true, each = N), sd = sigma_y), nrow = S, byrow = TRUE)

m0 <- mu_true
s0 <- 1.0
a0 <- 5.0
b0 <- (a0 - 1) * tau2_true

stan_data <- list(
  S = S,
  N = N,
  y = y,
  sigma_y = sigma_y,
  m0 = m0,
  s0 = s0,
  a0 = a0,
  b0 = b0
)

stan_file <- file.path("stan", "normal_normal_hierarchy_knownsd.stan")
stan_exe <- file.path(tempdir(), "normal_normal_hierarchy_knownsd")
stan_model <- cmdstan_model(stan_file, exe_file = stan_exe, quiet = TRUE)
stan_fit <- stan_model$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = min(4L, parallel::detectCores(logical = TRUE)),
  iter_warmup = 1000,
  iter_sampling = 1000,
  seed = 123L,
  refresh = 200
)

stan_mu <- as.numeric(stan_fit$draws(variables = "mu", format = "draws_matrix"))
stan_tau2 <- as.numeric(stan_fit$draws(variables = "tau2", format = "draws_matrix"))
stan_summary <- stan_fit$summary(variables = c("mu", "tau2"))

out <- list(
  simulation = list(
    S = S,
    N = N,
    sigma_y = sigma_y,
    mu_true = mu_true,
    tau2_true = tau2_true,
    alpha_true = alpha_true
  ),
  prior = list(
    m0 = m0,
    s0 = s0,
    a0 = a0,
    b0 = b0
  ),
  data = list(y = y),
  stan = list(
    mu = stan_mu,
    tau2 = stan_tau2,
    summary = stan_summary
  )
)

outfile <- file.path("samples", "normal_normal_stan_results.rds")
saveRDS(out, outfile)

cat("STAN_NORMAL_NORMAL_OK\n")
cat(sprintf("output=%s\n", outfile))
cat(sprintf("mu_mean=%.6f\n", mean(stan_mu)))
cat(sprintf("tau2_mean=%.6f\n", mean(stan_tau2)))
print(stan_summary[, c("variable", "mean", "sd", "ess_bulk", "rhat")], row.names = FALSE)
