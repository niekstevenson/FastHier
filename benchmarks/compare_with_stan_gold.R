#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(cmdstanr)
  library(bridgesampling)
})

source("outer_SMC.R")

weighted_quantile <- function(x, w, probs = c(0.025, 0.5, 0.975)) {
  o <- order(x)
  x <- x[o]
  w <- w[o]
  w <- w / sum(w)
  cw <- cumsum(w)
  sapply(probs, function(p) x[which(cw >= p)[1L]])
}

weighted_summary <- function(X, w) {
  w <- w / sum(w)
  means <- colSums(X * w)
  q <- vapply(seq_len(ncol(X)), function(j) weighted_quantile(X[, j], w), numeric(3))
  data.frame(
    parameter = colnames(X),
    mean = means,
    q2.5 = q[1, ],
    q50 = q[2, ],
    q97.5 = q[3, ],
    row.names = NULL,
    check.names = FALSE
  )
}

plain_summary <- function(X) {
  data.frame(
    parameter = colnames(X),
    mean = colMeans(X),
    q2.5 = apply(X, 2L, quantile, probs = 0.025),
    q50 = apply(X, 2L, quantile, probs = 0.5),
    q97.5 = apply(X, 2L, quantile, probs = 0.975),
    row.names = NULL,
    check.names = FALSE
  )
}

loglik_shifted_gamma <- function(Theta, y) {
  Theta <- as.matrix(Theta)
  eps <- 1e-9
  shape <- exp(Theta[, 1L]) + eps
  scale <- exp(Theta[, 2L]) + eps
  shift <- exp(Theta[, 3L]) + eps
  out <- rep(-Inf, nrow(Theta))
  min_y <- min(y)
  ok <- shift < min_y
  if (!any(ok)) return(out)
  for (i in which(ok)) {
    out[i] <- sum(dgamma(y - shift[i], shape = shape[i], scale = scale[i], log = TRUE))
  }
  out
}

logprior_eta <- function(Theta) {
  Theta <- as.matrix(Theta)
  rowSums(dnorm(Theta, mean = 0, sd = 1, log = TRUE))
}

rprior_eta <- function(n) {
  out <- matrix(rnorm(n * 3L), nrow = n, ncol = 3L)
  colnames(out) <- c("eta_shape", "eta_scale", "eta_shift")
  out
}

bridge_logpost <- function(s_row, data) {
  row <- matrix(as.numeric(s_row), nrow = 1L)
  logprior_eta(row) + loglik_shifted_gamma(row, data$y)
}

build_stan_gold <- function(gold_path) {
  set.seed(20260226)
  truth <- c(shape = 1.8, scale = 0.9, shift = 0.8)
  N <- 200L
  y <- truth["shift"] + rgamma(N, shape = truth["shape"], scale = truth["scale"])

  model <- cmdstan_model(file.path("stan", "shifted_gamma_logparam.stan"), quiet = TRUE)
  init_shift <- log(min(y) * 0.5)
  fit <- model$sample(
    data = list(N = N, y = y),
    chains = 4,
    parallel_chains = min(4, parallel::detectCores(logical = TRUE)),
    iter_warmup = 1000,
    iter_sampling = 1000,
    seed = 123,
    adapt_delta = 0.98,
    max_treedepth = 12,
    init = function() list(eta_shape = 0, eta_scale = 0, eta_shift = init_shift),
    refresh = 200
  )

  stan_draws <- as.matrix(fit$draws(variables = c("shape", "scale", "shift"), format = "draws_matrix"))
  stan_summary <- plain_summary(stan_draws)

  eta_draws_raw <- fit$draws(
    variables = c("eta_shape", "eta_scale", "eta_shift"),
    format = "matrix"
  )
  eta_draws <- matrix(as.numeric(eta_draws_raw), ncol = 3L)
  colnames(eta_draws) <- c("eta_shape", "eta_scale", "eta_shift")
  lb <- setNames(rep(-Inf, 3L), colnames(eta_draws))
  ub <- setNames(rep(Inf, 3L), colnames(eta_draws))
  bridge <- bridge_sampler(
    samples = eta_draws,
    log_posterior = bridge_logpost,
    data = list(y = y),
    lb = lb,
    ub = ub,
    silent = TRUE
  )

  gold <- list(
    created_utc = format(Sys.time(), tz = "UTC"),
    truth = truth,
    N = N,
    y = y,
    stan_summary = stan_summary,
    bridge_logml = as.numeric(bridge$logml),
    stan_setup = list(
      seed = 123,
      chains = 4,
      iter_warmup = 1000,
      iter_sampling = 1000,
      adapt_delta = 0.98,
      max_treedepth = 12
    )
  )
  saveRDS(gold, gold_path)
  gold
}

gold_path <- file.path("benchmarks", "shifted_gamma_stan_gold.rds")
if (!file.exists(gold_path)) {
  dir.create("benchmarks", showWarnings = FALSE, recursive = TRUE)
  cat("Creating Stan gold benchmark at:", gold_path, "\n")
  gold <- build_stan_gold(gold_path)
} else {
  cat("Loading Stan gold benchmark from:", gold_path, "\n")
  gold <- readRDS(gold_path)
}

smc_out <- outer_smc_direct(
  data = gold$y,
  rprior = rprior_eta,
  logprior = logprior_eta,
  loglik_fn = loglik_shifted_gamma,
  M = 4000L,
  cess_target = 0.95,
  resample_threshold = 0.5,
  n_moves = 2L,
  rw_scale_init = 1.0,
  max_rounds = 200L,
  n_cores_loglik = 1L,
  seed = 123,
  verbose = FALSE
)

smc_nat <- cbind(
  shape = exp(smc_out$theta[, 1L]) + 1e-9,
  scale = exp(smc_out$theta[, 2L]) + 1e-9,
  shift = exp(smc_out$theta[, 3L]) + 1e-9
)
smc_summary <- weighted_summary(smc_nat, smc_out$w)

cmp <- merge(gold$stan_summary, smc_summary, by = "parameter", suffixes = c("_stan", "_smc"))
cmp$truth <- gold$truth[cmp$parameter]
cmp$mean_diff <- cmp$mean_smc - cmp$mean_stan

cat("\nGold benchmark created:", gold$created_utc, "UTC\n")
cat("Gold benchmark path:", gold_path, "\n")
cat("\nPosterior mean comparison (Stan gold vs SMC):\n")
print(cmp[, c("parameter", "truth", "mean_stan", "mean_smc", "mean_diff")], row.names = FALSE)

cat("\nMarginal likelihood:\n")
cat(sprintf("  Stan gold bridge log ML: %.6f\n", gold$bridge_logml))
cat(sprintf("  Outer SMC log ML:        %.6f\n", smc_out$log_evidence))
cat(sprintf("  Difference (SMC-Stan):   %.6f\n", smc_out$log_evidence - gold$bridge_logml))
