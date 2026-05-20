#!/usr/bin/env Rscript

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[1L]))
} else {
  normalizePath("tests/test_local_smc_defensive_1d.R")
}
repo_dir <- dirname(dirname(script_path))
setwd(repo_dir)

source("smc_core.R")
source("reference_priors.R")
source("SMC_super_fast.R")

weighted_quantiles <- function(x, w, probs) {
  ord <- order(x)
  x <- x[ord]
  w <- w[ord] / sum(w)
  cw <- cumsum(w)
  keep <- c(TRUE, diff(cw) > 0)
  stats::approx(cw[keep], x[keep], xout = probs, rule = 2, ties = "ordered")$y
}

theta_name <- "theta"
S1 <- matrix(1, 1, 1, dimnames = list(theta_name, theta_name))
S2 <- matrix(25, 1, 1, dimnames = list(theta_name, theta_name))
reference_prior <- make_reference_prior_gaussian_mixture(
  component_means = list(setNames(0, theta_name), setNames(0, theta_name)),
  component_covs = list(S1, S2),
  weights = c(0.9, 0.1),
  label = "defensive_1d"
)

data_obj <- list(y = 1.25, sigma_y = 0.6)
loglik_fn <- function(Theta, data) {
  stats::dnorm(data$y, mean = as.matrix(Theta)[, 1L], sd = data$sigma_y, log = TRUE)
}

grid <- seq(-12, 12, length.out = 120001L)
dx <- grid[2L] - grid[1L]
prior_density <- exp(reference_prior_logpdf(
  reference_prior,
  matrix(grid, ncol = 1L, dimnames = list(NULL, theta_name))
))
lik_density <- stats::dnorm(data_obj$y, mean = grid, sd = data_obj$sigma_y)
posterior_mass <- prior_density * lik_density
posterior_mass <- posterior_mass / sum(posterior_mass)

exact <- c(
  mean = sum(grid * posterior_mass),
  q10 = weighted_quantiles(grid, posterior_mass, 0.1),
  q50 = weighted_quantiles(grid, posterior_mass, 0.5),
  q90 = weighted_quantiles(grid, posterior_mass, 0.9),
  logZ = log(sum(prior_density * lik_density) * dx)
)

fit <- run_tempered_smc(
  bridge_stat_fn = function(Theta) loglik_fn(Theta, data_obj),
  reference_prior = reference_prior,
  M = 3000L,
  n_mcmc_moves = 2L,
  max_rounds = 50L,
  G_mix = 6L,
  n_cores = 1L,
  seed = 42L,
  verbose = FALSE
)

theta <- as.numeric(fit$Theta[, 1L])
w <- fit$w / sum(fit$w)
qs <- weighted_quantiles(theta, w, c(0.1, 0.5, 0.9))
estimate <- c(
  mean = sum(w * theta),
  q10 = qs[1L],
  q50 = qs[2L],
  q90 = qs[3L],
  logZ = fit$log_evidence
)
errors <- estimate - exact

print(round(rbind(exact = exact, local_smc = estimate, error = errors), 6))

if (!identical(fit$transport$meta$method, "sparse_triangular")) {
  stop("Local SMC did not use the retained sparse triangular transport.")
}
if (abs(errors["mean"]) > 0.08) {
  stop("Posterior mean regression failed.")
}
if (max(abs(errors[c("q10", "q50", "q90")])) > 0.16) {
  stop("Posterior quantile regression failed.")
}
if (abs(errors["logZ"]) > 0.15) {
  stop("Log-evidence regression failed.")
}

cat("\nLocal SMC defensive-mixture 1D regression passed.\n")
