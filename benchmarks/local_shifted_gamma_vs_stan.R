#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(cmdstanr)
})

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[1L]))
} else {
  normalizePath("benchmarks/local_shifted_gamma_vs_stan.R")
}
repo_dir <- dirname(dirname(script_path))
setwd(repo_dir)

source("SMC_super_fast.R")

loglik_shifted_gamma <- function(Theta, y) {
  Theta <- as.matrix(Theta)
  eps <- 1e-9
  shape <- exp(Theta[, 1L]) + eps
  scale <- exp(Theta[, 2L]) + eps
  shift <- exp(Theta[, 3L]) + eps
  out <- rep(-1e12, nrow(Theta))
  min_y <- min(y)
  ok <- shift < min_y
  if (!any(ok)) return(out)
  for (i in which(ok)) {
    out[i] <- sum(dgamma(y - shift[i], shape = shape[i], scale = scale[i], log = TRUE))
  }
  out[!is.finite(out)] <- -1e12
  out
}

set.seed(20260311)

dir.create("benchmarks", showWarnings = FALSE, recursive = TRUE)
dir.create("samples", showWarnings = FALSE, recursive = TRUE)

theta_names <- c("eta_shape", "eta_scale", "eta_shift")

theta_true <- c(
  eta_shape = log(4.0),
  eta_scale = log(0.6),
  eta_shift = log(0.25)
)

mu_ref <- c(
  eta_shape = log(3.5),
  eta_scale = log(0.7),
  eta_shift = log(0.20)
)

Sigma_ref <- diag(c(0.35, 0.35, 0.25)^2)
colnames(Sigma_ref) <- rownames(Sigma_ref) <- theta_names
sd_ref <- sqrt(diag(Sigma_ref))

N <- 200L
y <- exp(theta_true[3L]) + rgamma(
  N,
  shape = exp(theta_true[1L]),
  scale = exp(theta_true[2L])
)

smc_fit <- enhanced_smc_elite(
  data = y,
  loglik_fn = loglik_shifted_gamma,
  mu_ref = mu_ref,
  Sigma_ref = Sigma_ref,
  M = 4000L,
  resample_threshold = 0.6,
  n_mcmc_moves = 3L,
  max_rounds = 150L,
  G_mix = 8L,
  da_enable = TRUE,
  gss_enable = TRUE,
  ll_cache_enable = TRUE,
  deterministic_resampling = FALSE,
  n_cores = 1L,
  seed = 20260312L,
  verbose = FALSE
)

stan_code <- "
data {
  int<lower=1> N;
  array[N] real<lower=0> y;
  vector[3] mu_ref;
  vector<lower=0>[3] sd_ref;
}

transformed data {
  real min_y = min(y);
}

parameters {
  real<lower=0> shape;
  real<lower=0> scale;
  real<lower=0, upper=min_y> shift;
}

model {
  shape ~ lognormal(mu_ref[1], sd_ref[1]);
  scale ~ lognormal(mu_ref[2], sd_ref[2]);
  shift ~ lognormal(mu_ref[3], sd_ref[3]);

  for (n in 1:N) {
    target += gamma_lpdf(y[n] - shift | shape, inv(scale));
  }
}

generated quantities {
  vector[3] theta;
  theta[1] = shape;
  theta[2] = scale;
  theta[3] = shift;
}
"

stan_file <- file.path(tempdir(), "local_shifted_gamma.stan")
writeLines(stan_code, stan_file)

stan_model <- cmdstan_model(
  stan_file,
  exe_file = file.path(tempdir(), "local_shifted_gamma"),
  quiet = TRUE
)

stan_fit <- stan_model$sample(
  data = list(
    N = N,
    y = y,
    mu_ref = unname(mu_ref),
    sd_ref = unname(sd_ref)
  ),
  chains = 4,
  parallel_chains = min(4, parallel::detectCores(logical = TRUE)),
  iter_warmup = 500,
  iter_sampling = 500,
  seed = 20260313,
  init = function(chain_id) {
    list(
      shape = exp(rnorm(1L, mean = mu_ref[1L], sd = 0.1)),
      scale = exp(rnorm(1L, mean = mu_ref[2L], sd = 0.1)),
      shift = min(y) * runif(1L, min = 0.2, max = 0.6)
    )
  },
  adapt_delta = 0.99,
  max_treedepth = 12,
  refresh = 100
)

stan_theta <- as.matrix(stan_fit$draws(variables = "theta", format = "draws_matrix"))
colnames(stan_theta) <- c("shape", "scale", "shift")

smc_theta <- exp(as.matrix(smc_fit$Theta))
colnames(smc_theta) <- c("shape", "scale", "shift")
smc_w <- as.numeric(smc_fit$w)
smc_w <- smc_w / sum(smc_w)

set.seed(20260314)
smc_idx <- sample.int(nrow(smc_theta), size = nrow(stan_theta), replace = TRUE, prob = smc_w)
smc_theta_draws <- smc_theta[smc_idx, , drop = FALSE]

png(
  filename = file.path("benchmarks", "local_shifted_gamma_vs_stan.png"),
  width = 1200,
  height = 400,
  res = 120
)
par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))

for (j in seq_len(3L)) {
  d_stan <- density(stan_theta[, j])
  d_smc <- density(smc_theta_draws[, j])
  plot(
    d_stan,
    main = colnames(stan_theta)[j],
    xlab = colnames(stan_theta)[j],
    ylab = "Density",
    lwd = 2,
    col = "black"
  )
  lines(d_smc, col = "red3", lwd = 2)
  abline(v = exp(theta_true[j]), col = "grey60", lty = 3)
  legend(
    "topright",
    legend = c("Stan", "Local SMC", "Truth"),
    col = c("black", "red3", "grey60"),
    lwd = c(2, 2, 1),
    lty = c(1, 1, 3),
    bty = "n"
  )
}

dev.off()

out <- list(
  config = list(
    seed_data = 20260311L,
    seed_smc = 20260312L,
    seed_stan = 20260313L,
    N = N
  ),
  truth = exp(theta_true),
  prior = list(mu_ref = mu_ref, Sigma_ref = Sigma_ref),
  data = list(y = y),
  smc = list(theta = smc_theta, w = smc_w, fit = smc_fit),
  stan = list(theta = stan_theta, summary = stan_fit$summary(variables = "theta"))
)

saveRDS(out, file.path("samples", "local_shifted_gamma_vs_stan.rds"))

cat("Saved: benchmarks/local_shifted_gamma_vs_stan.png\n")
cat("Saved: samples/local_shifted_gamma_vs_stan.rds\n")
cat(sprintf("shape mean: stan=%.4f smc=%.4f\n", mean(stan_theta[, 1L]), sum(smc_theta[, 1L] * smc_w)))
cat(sprintf("scale mean: stan=%.4f smc=%.4f\n", mean(stan_theta[, 2L]), sum(smc_theta[, 2L] * smc_w)))
cat(sprintf("shift mean: stan=%.4f smc=%.4f\n", mean(stan_theta[, 3L]), sum(smc_theta[, 3L] * smc_w)))
