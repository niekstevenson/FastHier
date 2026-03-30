#!/usr/bin/env Rscript

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[1L]))
} else {
  normalizePath("tests/test_gss_da_exact_1d.R")
}
repo_dir <- dirname(dirname(script_path))
setwd(repo_dir)

source("smc_core.R")
source("reference_priors.R")
source("SMC_super_fast.R")

weighted_quantiles <- function(x, w, probs) {
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  cw <- cumsum(w)
  keep <- c(TRUE, diff(cw) > 0)
  approx(cw[keep], x[keep], xout = probs, rule = 2, ties = "ordered")$y
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

run_case <- function(label, gss_enable, da_enable, M = 8000L) {
  fit <- enhanced_smc_elite(
    data = data_obj,
    loglik_fn = loglik_fn,
    reference_prior = reference_prior,
    M = M,
    n_mcmc_moves = 3L,
    max_rounds = 60L,
    hist_mix_enable = TRUE,
    gss_enable = gss_enable,
    da_enable = da_enable,
    n_cores = 1L,
    seed = 42L,
    verbose = FALSE
  )
  theta <- as.numeric(fit$Theta[, 1L])
  w <- fit$w / sum(fit$w)
  qs <- weighted_quantiles(theta, w, c(0.1, 0.5, 0.9))
  c(
    mean = sum(w * theta),
    q10 = qs[1L],
    q50 = qs[2L],
    q90 = qs[3L],
    logZ = fit$log_evidence
  )
}

results <- rbind(
  exact = exact,
  baseline = run_case("baseline", FALSE, FALSE),
  gss_only = run_case("gss_only", TRUE, FALSE),
  da_only = run_case("da_only", FALSE, TRUE),
  gss_da = run_case("gss_da", TRUE, TRUE)
)
errors <- sweep(results[-1L, , drop = FALSE], 2L, exact, "-")

print(round(results, 6))
cat("\nErrors vs exact\n")
print(round(errors, 6))

mean_tol <- 0.03
quantile_tol <- 0.06
logz_tol <- 0.08

if (any(abs(errors[, "mean"]) > mean_tol)) {
  stop("Posterior mean regression failed for at least one GSS/DA configuration.")
}
if (any(apply(abs(errors[, c("q10", "q50", "q90"), drop = FALSE]), 1L, max) > quantile_tol)) {
  stop("Posterior quantile regression failed for at least one GSS/DA configuration.")
}
if (any(abs(errors[, "logZ"]) > logz_tol)) {
  stop("Log-evidence regression failed for at least one GSS/DA configuration.")
}

cat("\nGSS/DA exact 1D regression passed.\n")
