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
  sds <- sqrt(colSums((sweep(X, 2L, means, "-")^2) * w))
  q <- vapply(seq_len(ncol(X)), function(j) weighted_quantile(X[, j], w), numeric(3))
  out <- data.frame(
    parameter = colnames(X),
    mean = means,
    sd = sds,
    q2.5 = q[1, ],
    q50 = q[2, ],
    q97.5 = q[3, ],
    row.names = NULL,
    check.names = FALSE
  )
  out
}

plain_summary <- function(X) {
  data.frame(
    parameter = colnames(X),
    mean = colMeans(X),
    sd = apply(X, 2L, sd),
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
  yy <- y
  for (i in which(ok)) {
    out[i] <- sum(dgamma(yy - shift[i], shape = shape[i], scale = scale[i], log = TRUE))
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

run_check <- function() {
  set.seed(20260226)
  truth <- c(shape = 1.8, scale = 0.9, shift = 0.8)
  N <- 200L
  y <- truth["shift"] + rgamma(N, shape = truth["shape"], scale = truth["scale"])

  stan_file <- file.path("stan", "shifted_gamma_logparam.stan")
  model <- cmdstan_model(stan_file, quiet = TRUE)
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

  stan_draws <- as.data.frame(
    fit$draws(variables = c("shape", "scale", "shift"), format = "draws_matrix")
  )
  stan_summary <- plain_summary(as.matrix(stan_draws))

  smc_out <- outer_smc_direct(
    data = y,
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
    verbose = TRUE
  )

  smc_nat <- cbind(
    shape = exp(smc_out$theta[, 1L]) + 1e-9,
    scale = exp(smc_out$theta[, 2L]) + 1e-9,
    shift = exp(smc_out$theta[, 3L]) + 1e-9
  )
  smc_summary <- weighted_summary(smc_nat, smc_out$w)

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

  cmp <- merge(stan_summary, smc_summary, by = "parameter", suffixes = c("_stan", "_smc"))
  cmp$truth <- truth[cmp$parameter]
  cmp$mean_diff <- cmp$mean_smc - cmp$mean_stan

  cat("\nTruth:\n")
  print(truth)
  cat("\nPosterior Summary Comparison (Stan vs Outer SMC):\n")
  print(cmp[, c(
    "parameter", "truth",
    "mean_stan", "mean_smc", "mean_diff",
    "q2.5_stan", "q2.5_smc", "q50_stan", "q50_smc", "q97.5_stan", "q97.5_smc"
  )], row.names = FALSE)

  cat("\nMarginal Likelihood Comparison:\n")
  cat(sprintf("  Stan bridge log ML: %.6f\n", bridge$logml))
  cat(sprintf("  Outer SMC log ML:   %.6f\n", smc_out$log_evidence))
  cat(sprintf("  Difference (SMC-Stan): %.6f\n", smc_out$log_evidence - bridge$logml))
}

tryCatch(
  run_check(),
  error = function(e) {
    cat("\nVerification run failed:\n")
    cat(conditionMessage(e), "\n")
    quit(save = "no", status = 1L)
  }
)
