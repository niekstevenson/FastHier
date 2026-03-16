#!/usr/bin/env Rscript

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[1L]))
} else {
  normalizePath("benchmarks/simple_collapsed_subject_from_stan.R")
}
repo_dir <- dirname(dirname(script_path))
setwd(repo_dir)

suppressPackageStartupMessages({
  library(parallel)
})

source("smc_core.R")
source("SMC_super_fast.R")
source("hierarchical_smc.R")

stan_results_file <- file.path("samples", "normal_normal_meanvar_hierarchy_compare_results.rds")
output_file <- file.path("samples", "simple_collapsed_subject_from_stan_results.rds")
plot_file <- file.path("benchmarks", "simple_collapsed_subject_from_stan.png")

if (!file.exists(stan_results_file)) {
  stop("Missing Stan benchmark results: ", stan_results_file)
}

dir.create("samples", showWarnings = FALSE, recursive = TRUE)
dir.create("benchmarks", showWarnings = FALSE, recursive = TRUE)

stan_bundle <- readRDS(stan_results_file)

y <- stan_bundle$simulation$y
phi_prior_mean <- stan_bundle$prior$mean
phi_prior_sd <- stan_bundle$prior$sd
stan_draws <- stan_bundle$stan_draws

S <- nrow(y)
N_trials <- ncol(y)
data_list <- lapply(seq_len(S), function(i) y[i, ])

theta_names <- c("subject_mean", "subject_log_var")
phi_names <- c("mu_mean", "mu_log_var", "log_tau2_mean", "log_tau2_log_var")

make_phi_prior_normal <- function(mean, sd) {
  mean <- as.numeric(mean)
  sd <- as.numeric(sd)
  nm <- names(mean)

  rprior_phi <- function(n) {
    out <- matrix(
      rnorm(n * length(mean), mean = rep(mean, each = n), sd = rep(sd, each = n)),
      nrow = n,
      ncol = length(mean),
      byrow = FALSE
    )
    colnames(out) <- nm
    out
  }

  logprior_phi <- function(phi_row) {
    sum(dnorm(as.numeric(phi_row), mean = mean, sd = sd, log = TRUE))
  }

  list(rprior = rprior_phi, lprior = logprior_phi)
}

phi_prior <- make_phi_prior_normal(phi_prior_mean, phi_prior_sd)
gaussian_map_fn <- phi_to_gaussian_params_diag_factory(param_names = theta_names)

loglik_fn <- function(Theta, y_i) {
  Theta <- as.matrix(Theta)
  colnames(Theta) <- theta_names
  vapply(
    seq_len(nrow(Theta)),
    function(i) {
      mu_i <- Theta[i, 1L]
      sd_i <- pmax(exp(0.5 * Theta[i, 2L]), 1e-9)
      sum(dnorm(y_i, mean = mu_i, sd = sd_i, log = TRUE))
    },
    numeric(1L)
  )
}

local_mu_ref <- c(
  subject_mean = unname(phi_prior_mean["mu_mean"]),
  subject_log_var = unname(phi_prior_mean["mu_log_var"])
)
local_Sigma_ref <- structure(
  diag(c(4.0, 2.0), 2L),
  dimnames = list(theta_names, theta_names)
)

detected_cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
if (!is.finite(detected_cores) || detected_cores < 1L) {
  detected_cores <- 1L
}
mc_cores <- as.integer(max(1L, min(4L, detected_cores)))
local_particles <- 512L
local_rounds <- 40L
outer_particles <- 64L
outer_moves <- 1L
subject_rejuvenate_every <- 3L

timed <- function(expr) {
  value <- NULL
  tm <- system.time({
    value <- force(expr)
  })
  list(
    value = value,
    elapsed = unname(tm["elapsed"]),
    user = unname(tm["user.self"]),
    system = unname(tm["sys.self"])
  )
}

set.seed(777L)
subject_order <- sample.int(S)

cat(sprintf("Loaded Stan results from %s\n", stan_results_file))
cat(sprintf("Benchmark data: %d subjects x %d trials\n", S, N_trials))
cat("Running local subject fits...\n")

local_stage <- timed(
  run_local_smc_subjects(
    data_list = data_list,
    loglik_fn = loglik_fn,
    mu_ref = local_mu_ref,
    Sigma_ref = local_Sigma_ref,
    M = local_particles,
    max_rounds = local_rounds,
    n_mcmc_moves = 1L,
    gss_enable = FALSE,
    hist_mix_enable = FALSE,
    da_enable = FALSE,
    n_cores = mc_cores,
    base_seed = 1001L,
    verbose = FALSE
  )
)

cat("Running collapsed subject SMC...\n")

collapsed_stage <- timed(
  collapsed_subject_smc(
    data_list = data_list,
    loglik_fn = loglik_fn,
    prior = phi_prior,
    gaussian_map_fn = gaussian_map_fn,
    d_theta = length(theta_names),
    N = outer_particles,
    local_particles = local_particles,
    resample_threshold = 0.5,
    subject_order = subject_order,
    rejuvenation = "recompute",
    rejuvenate_every = subject_rejuvenate_every,
    rejuvenate_after_resample = TRUE,
    n_rejuvenation_moves = outer_moves,
    auxiliary_proposal = "correlated",
    auxiliary_correlation = 0.95,
    warm_start_mode = "nearest",
    theta_names = theta_names,
    local_n_cores = 1L,
    outer_n_cores = mc_cores,
    base_seed = 2026L,
    assume_unbiased_local_evidence = TRUE,
    evidence_cache = TRUE,
    store_local_fits = FALSE,
    verbose = TRUE,
    local_smc_control = list(
      max_rounds = local_rounds,
      n_mcmc_moves = 1L,
      gss_enable = FALSE,
      hist_mix_enable = FALSE,
      da_enable = FALSE
    )
  )
)

resample_weighted_rows <- function(phi, w, n_draws, seed = 1L) {
  set.seed(as.integer(seed))
  phi <- as.matrix(phi)
  w <- pmax(as.numeric(w), 0)
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) {
    w <- rep(1 / nrow(phi), nrow(phi))
  } else {
    w <- w / sw
  }
  idx <- sample.int(nrow(phi), size = n_draws, replace = TRUE, prob = w)
  phi[idx, , drop = FALSE]
}

phi_draws_to_named_list <- function(phi_draws) {
  phi_draws <- as.matrix(phi_draws)
  colnames(phi_draws) <- phi_names
  list(
    mu_mean = phi_draws[, "mu_mean"],
    mu_log_var = phi_draws[, "mu_log_var"],
    tau2_mean = exp(phi_draws[, "log_tau2_mean"]),
    tau2_log_var = exp(phi_draws[, "log_tau2_log_var"])
  )
}

collapsed_draws <- phi_draws_to_named_list(
  resample_weighted_rows(collapsed_stage$value$phi, collapsed_stage$value$w, n_draws = 4000L, seed = 99L)
)

summary_rows <- function(draws, method) {
  rbind(
    data.frame(parameter = "mu_mean", method = method, mean = mean(draws$mu_mean), sd = sd(draws$mu_mean)),
    data.frame(parameter = "tau2_mean", method = method, mean = mean(draws$tau2_mean), sd = sd(draws$tau2_mean)),
    data.frame(parameter = "mu_log_var", method = method, mean = mean(draws$mu_log_var), sd = sd(draws$mu_log_var)),
    data.frame(parameter = "tau2_log_var", method = method, mean = mean(draws$tau2_log_var), sd = sd(draws$tau2_log_var))
  )
}

comparison_summary <- rbind(
  summary_rows(stan_draws, "Stan"),
  summary_rows(collapsed_draws, "Collapsed subject")
)

plot_param <- function(stan_x, smc_x, main, xlab) {
  d_stan <- density(stan_x)
  d_smc <- density(smc_x)
  xlim <- range(c(d_stan$x, d_smc$x))
  ylim <- c(0, 1.05 * max(d_stan$y, d_smc$y))
  plot(d_stan, lwd = 2, col = "black", main = main, xlab = xlab, ylab = "Density", xlim = xlim, ylim = ylim)
  lines(d_smc, lwd = 2, col = "firebrick3")
  legend("topright", legend = c("Stan", "Collapsed subject"), col = c("black", "firebrick3"), lwd = 2, bty = "n", cex = 0.8)
}

png(plot_file, width = 1200, height = 900)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
plot_param(stan_draws$mu_mean, collapsed_draws$mu_mean, "mu_mean", "mu_mean")
plot_param(stan_draws$tau2_mean, collapsed_draws$tau2_mean, "tau2_mean", "tau2_mean")
plot_param(stan_draws$mu_log_var, collapsed_draws$mu_log_var, "mu_log_var", "mu_log_var")
plot_param(stan_draws$tau2_log_var, collapsed_draws$tau2_log_var, "tau2_log_var", "tau2_log_var")
dev.off()

results <- list(
  stan_source = stan_results_file,
  simulation = list(
    S = S,
    N_trials = N_trials,
    y = y,
    subject_order = subject_order
  ),
  settings = list(
    local_particles = local_particles,
    local_rounds = local_rounds,
    outer_particles = outer_particles,
    outer_moves = outer_moves,
    subject_rejuvenate_every = subject_rejuvenate_every,
    auxiliary_proposal = "correlated",
    auxiliary_correlation = 0.95
  ),
  timings = data.frame(
    stage = c("locals", "collapsed_subject"),
    elapsed = c(local_stage$elapsed, collapsed_stage$elapsed),
    user = c(local_stage$user, collapsed_stage$user),
    system = c(local_stage$system, collapsed_stage$system)
  ),
  stan_draws = stan_draws,
  local_fits = local_stage$value,
  collapsed_fit = collapsed_stage$value,
  collapsed_draws = collapsed_draws,
  comparison_summary = comparison_summary,
  plot_file = plot_file
)

saveRDS(results, output_file)

cat(sprintf("Saved results to %s\n", output_file))
cat(sprintf("Saved plot to %s\n", plot_file))
print(comparison_summary)
print(results$timings)
