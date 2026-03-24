#!/usr/bin/env Rscript

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[1L]))
} else {
  normalizePath("benchmarks/normal_normal_meanvar_hierarchy_compare.R")
}
repo_dir <- dirname(dirname(script_path))
setwd(repo_dir)

suppressPackageStartupMessages({
  library(cmdstanr)
  library(parallel)
})

source("smc_core.R")
source("SMC_super_fast.R")
source("hierarchical_smc.R")

dir.create("benchmarks", showWarnings = FALSE, recursive = TRUE)
dir.create("samples", showWarnings = FALSE, recursive = TRUE)
dir.create("stan", showWarnings = FALSE, recursive = TRUE)

results_file <- file.path("samples", "normal_normal_meanvar_hierarchy_compare_results.rds")
cache_dir <- file.path("samples", "normal_normal_meanvar_hierarchy_compare_cache")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(20260314L)

force_recompute <- FALSE
stage_versions <- list(
  stan = 1L,
  shared_local = 1L,
  shared_outer = 1L,
  collapsed_subject = 2L,
  collapsed_bridge = 2L
)

S <- 15L
N_trials <- 100L
mc_cores <- min(4L, parallel::detectCores(logical = TRUE))

stan_chains <- 4L
stan_iter_warmup <- 1000L
stan_iter_sampling <- 1000L

shared_local_particles <- 1200L
shared_local_rounds <- 100L
shared_outer_particles <- 1200L
shared_bank_particles <- 96L
shared_outer_rounds <- 50L

collapsed_outer_particles <- 64L
collapsed_bridge_outer_particles <- 64L
collapsed_local_particles <- 192L
collapsed_local_rounds <- 30L
collapsed_rejuvenation <- "recompute"
collapsed_rejuvenate_every <- 3L
collapsed_rejuvenate_after_resample <- TRUE
collapsed_n_rejuvenation_moves <- 1L
collapsed_bridge_steps <- 4L
collapsed_bridge_rejuvenation <- "recompute"
collapsed_bridge_rejuvenate_every <- 6L
collapsed_bridge_rejuvenate_after_resample <- TRUE
collapsed_bridge_n_rejuvenation_moves <- 1L

n_plot_draws <- 4000L

theta_names <- c("subject_mean", "subject_log_var")
phi_names <- c("mu_mean", "mu_log_var", "log_tau2_mean", "log_tau2_log_var")

true_phi <- c(
  mu_mean = 1.00,
  mu_log_var = log(0.80^2),
  log_tau2_mean = log(0.50),
  log_tau2_log_var = log(0.20)
)

subject_mean_true <- rnorm(S, mean = true_phi["mu_mean"], sd = sqrt(exp(true_phi["log_tau2_mean"])))
subject_log_var_true <- rnorm(S, mean = true_phi["mu_log_var"], sd = sqrt(exp(true_phi["log_tau2_log_var"])))

y <- matrix(NA_real_, nrow = S, ncol = N_trials)
for (s in seq_len(S)) {
  y[s, ] <- rnorm(
    N_trials,
    mean = subject_mean_true[s],
    sd = sqrt(exp(subject_log_var_true[s]))
  )
}
data_list <- lapply(seq_len(S), function(s) y[s, ])

phi_prior_mean <- c(
  mu_mean = 0.0,
  mu_log_var = log(1.0),
  log_tau2_mean = log(0.5),
  log_tau2_log_var = log(0.25)
)
phi_prior_sd <- c(
  mu_mean = 2.0,
  mu_log_var = 1.0,
  log_tau2_mean = 1.0,
  log_tau2_log_var = 1.0
)

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

timed_eval <- function(expr) {
  value <- NULL
  tm <- system.time({
    value <- force(expr)
  })
  list(
    value = value,
    user = unname(tm["user.self"]),
    system = unname(tm["sys.self"]),
    elapsed = unname(tm["elapsed"])
  )
}

cache_key_from_object <- function(obj) {
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(obj, tmp)
  unname(tools::md5sum(tmp))
}

load_or_run_cached <- function(stage, settings, run_fn, force = FALSE) {
  key <- cache_key_from_object(settings)
  cache_file <- file.path(cache_dir, sprintf("%s_%s.rds", stage, key))
  if (!isTRUE(force) && file.exists(cache_file)) {
    out <- readRDS(cache_file)
    out$cache_file <- cache_file
    return(out)
  }

  timed <- timed_eval(run_fn())
  out <- list(
    value = timed$value,
    user = timed$user,
    system = timed$system,
    elapsed = timed$elapsed,
    cache_file = cache_file,
    settings = settings
  )
  saveRDS(out, cache_file)
  out
}

local_working_prior <- make_working_prior_gaussian(
  mu = c(subject_mean = phi_prior_mean["mu_mean"], subject_log_var = phi_prior_mean["mu_log_var"]),
  Sigma = structure(
    diag(c(4.0, 2.0), 2L),
    dimnames = list(theta_names, theta_names)
  )
)

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

stan_data <- list(
  S = S,
  N = N_trials,
  y = y,
  prior_mu_mean = unname(phi_prior_mean["mu_mean"]),
  prior_sd_mu_mean = unname(phi_prior_sd["mu_mean"]),
  prior_mu_log_var = unname(phi_prior_mean["mu_log_var"]),
  prior_sd_mu_log_var = unname(phi_prior_sd["mu_log_var"]),
  prior_mean_log_tau2_mean = unname(phi_prior_mean["log_tau2_mean"]),
  prior_sd_log_tau2_mean = unname(phi_prior_sd["log_tau2_mean"]),
  prior_mean_log_tau2_log_var = unname(phi_prior_mean["log_tau2_log_var"]),
  prior_sd_log_tau2_log_var = unname(phi_prior_sd["log_tau2_log_var"])
)

stan_file <- file.path("stan", "normal_normal_meanvar_hierarchy_logvar.stan")
stan_model <- cmdstanr::cmdstan_model(stan_file, quiet = TRUE)
code_fingerprint <- list(
  benchmark = unname(tools::md5sum(script_path)),
  hierarchical = unname(tools::md5sum("hierarchical_smc.R")),
  local = unname(tools::md5sum("SMC_super_fast.R")),
  core = unname(tools::md5sum("smc_core.R")),
  stan = unname(tools::md5sum(stan_file))
)

simulation_fingerprint <- list(
  seed = 20260314L,
  S = S,
  N_trials = N_trials,
  true_phi = unname(true_phi),
  phi_prior_mean = unname(phi_prior_mean),
  phi_prior_sd = unname(phi_prior_sd)
)

stan_stage <- load_or_run_cached(
  stage = "stan",
  settings = list(
    stage_version = stage_versions$stan,
    code = code_fingerprint,
    simulation = simulation_fingerprint,
    stan_data = stan_data,
    chains = stan_chains,
    iter_warmup = stan_iter_warmup,
    iter_sampling = stan_iter_sampling
  ),
  run_fn = function() {
    fit <- stan_model$sample(
      data = stan_data,
      chains = stan_chains,
      parallel_chains = min(stan_chains, mc_cores),
      iter_warmup = stan_iter_warmup,
      iter_sampling = stan_iter_sampling,
      seed = 123L,
      refresh = 250L
    )
    list(
      draws = list(
        mu_mean = as.numeric(fit$draws("mu_mean", format = "draws_matrix")),
        mu_log_var = as.numeric(fit$draws("mu_log_var", format = "draws_matrix")),
        tau2_mean = as.numeric(fit$draws("tau2_mean", format = "draws_matrix")),
        tau2_log_var = as.numeric(fit$draws("tau2_log_var", format = "draws_matrix"))
      ),
      summary = fit$summary(variables = c("mu_mean", "mu_log_var", "tau2_mean", "tau2_log_var"))
    )
  },
  force = force_recompute
)

stan_draws <- stan_stage$value$draws
stan_summary <- stan_stage$value$summary

shared_local_stage <- load_or_run_cached(
  stage = "shared_local",
  settings = list(
    stage_version = stage_versions$shared_local,
    code = code_fingerprint,
    simulation = simulation_fingerprint,
    working_prior = local_working_prior,
    M = shared_local_particles,
    max_rounds = shared_local_rounds
  ),
  run_fn = function() {
    run_local_smc_subjects(
      data_list = data_list,
      loglik_fn = loglik_fn,
      mu_ref = local_working_prior$mu,
      Sigma_ref = local_working_prior$Sigma,
      M = shared_local_particles,
      n_cores = mc_cores,
      base_seed = 5000L,
      verbose = FALSE,
      max_rounds = shared_local_rounds
    )
  },
  force = force_recompute
)
local_fits <- shared_local_stage$value

shared_outer_stage <- load_or_run_cached(
  stage = "shared_outer",
  settings = list(
    stage_version = stage_versions$shared_outer,
    code = code_fingerprint,
    simulation = simulation_fingerprint,
    shared_local = list(M = shared_local_particles, max_rounds = shared_local_rounds),
    outer_particles = shared_outer_particles,
    bank_particles = shared_bank_particles,
    outer_rounds = shared_outer_rounds
  ),
  run_fn = function() {
    local_objs <- build_local_exact_objects(
      local_fits = local_fits,
      data_list = data_list,
      loglik_fn = loglik_fn,
      working_prior = local_working_prior,
      base_seed = 9000L
    )
    hierarchical_smc(
      method = "shared_bank",
      local_objs = local_objs,
      prior = phi_prior,
      gaussian_map_fn = gaussian_map_fn,
      N = shared_outer_particles,
      M_local = shared_bank_particles,
      rho_res = 0.5,
      rho_local = 0.5,
      n_population_moves = 1L,
      max_rounds = shared_outer_rounds,
      max_bank_topups = 3L,
      bank_split_tol = 0.01,
      seed = 202L,
      verbose = TRUE
    )
  },
  force = force_recompute
)
shared_fit <- shared_outer_stage$value

set.seed(777L)
subject_order <- sample.int(S)

collapsed_subject_stage <- load_or_run_cached(
  stage = "collapsed_subject",
  settings = list(
    stage_version = stage_versions$collapsed_subject,
    code = code_fingerprint,
    simulation = simulation_fingerprint,
    subject_order = subject_order,
    outer_particles = collapsed_outer_particles,
    local_particles = collapsed_local_particles,
    local_rounds = collapsed_local_rounds,
    rejuvenation = collapsed_rejuvenation,
    rejuvenate_every = collapsed_rejuvenate_every,
    rejuvenate_after_resample = collapsed_rejuvenate_after_resample,
    n_rejuvenation_moves = collapsed_n_rejuvenation_moves
  ),
  run_fn = function() {
    hierarchical_smc(
      method = "collapsed_subject",
      data_list = data_list,
      loglik_fn = loglik_fn,
      prior = phi_prior,
      gaussian_map_fn = gaussian_map_fn,
      d_theta = length(theta_names),
      N = collapsed_outer_particles,
      local_particles = collapsed_local_particles,
      resample_threshold = 0.5,
      subject_order = subject_order,
      rejuvenation = collapsed_rejuvenation,
      rejuvenate_every = collapsed_rejuvenate_every,
      rejuvenate_after_resample = collapsed_rejuvenate_after_resample,
      n_rejuvenation_moves = collapsed_n_rejuvenation_moves,
      theta_names = theta_names,
      local_n_cores = 1L,
      outer_n_cores = mc_cores,
      base_seed = 404L,
      verbose = TRUE,
      local_smc_control = list(
        max_rounds = collapsed_local_rounds,
        n_mcmc_moves = 1L,
        gss_enable = FALSE,
        hist_mix_enable = FALSE,
        da_enable = FALSE
      )
    )
  },
  force = force_recompute
)
collapsed_fit <- collapsed_subject_stage$value

collapsed_bridge_stage <- load_or_run_cached(
  stage = "collapsed_bridge",
  settings = list(
    stage_version = stage_versions$collapsed_bridge,
    code = code_fingerprint,
    simulation = simulation_fingerprint,
    subject_order = subject_order,
    outer_particles = collapsed_bridge_outer_particles,
    local_particles = collapsed_local_particles,
    local_rounds = collapsed_local_rounds,
    bridge_steps = collapsed_bridge_steps,
    rejuvenation = collapsed_bridge_rejuvenation,
    rejuvenate_every = collapsed_bridge_rejuvenate_every,
    rejuvenate_after_resample = collapsed_bridge_rejuvenate_after_resample,
    n_rejuvenation_moves = collapsed_bridge_n_rejuvenation_moves
  ),
  run_fn = function() {
    hierarchical_smc(
      method = "collapsed_bridge",
      data_list = data_list,
      loglik_fn = loglik_fn,
      prior = phi_prior,
      gaussian_map_fn = gaussian_map_fn,
      d_theta = length(theta_names),
      N = collapsed_bridge_outer_particles,
      local_particles = collapsed_local_particles,
      resample_threshold = 0.5,
      subject_order = subject_order,
      n_bridge_steps = collapsed_bridge_steps,
      rejuvenation = collapsed_bridge_rejuvenation,
      rejuvenate_every = collapsed_bridge_rejuvenate_every,
      rejuvenate_after_resample = collapsed_bridge_rejuvenate_after_resample,
      n_rejuvenation_moves = collapsed_bridge_n_rejuvenation_moves,
      theta_names = theta_names,
      local_n_cores = 1L,
      outer_n_cores = mc_cores,
      base_seed = 505L,
      verbose = TRUE,
      local_smc_control = list(
        max_rounds = collapsed_local_rounds,
        n_mcmc_moves = 1L,
        gss_enable = FALSE,
        hist_mix_enable = FALSE,
        da_enable = FALSE
      )
    )
  },
  force = force_recompute
)
collapsed_bridge_fit <- collapsed_bridge_stage$value

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

shared_draws <- phi_draws_to_named_list(
  resample_weighted_rows(shared_fit$phi, shared_fit$w, n_plot_draws, seed = 1001L)
)
collapsed_draws <- phi_draws_to_named_list(
  resample_weighted_rows(collapsed_fit$phi, collapsed_fit$w, n_plot_draws, seed = 1002L)
)
collapsed_bridge_draws <- phi_draws_to_named_list(
  resample_weighted_rows(collapsed_bridge_fit$phi, collapsed_bridge_fit$w, n_plot_draws, seed = 1003L)
)

summary_rows <- function(draws, method) {
  rbind(
    data.frame(parameter = "mu_mean", method = method, mean = mean(draws$mu_mean), sd = stats::sd(draws$mu_mean)),
    data.frame(parameter = "tau2_mean", method = method, mean = mean(draws$tau2_mean), sd = stats::sd(draws$tau2_mean)),
    data.frame(parameter = "mu_log_var", method = method, mean = mean(draws$mu_log_var), sd = stats::sd(draws$mu_log_var)),
    data.frame(parameter = "tau2_log_var", method = method, mean = mean(draws$tau2_log_var), sd = stats::sd(draws$tau2_log_var))
  )
}

comparison_summary <- do.call(
  rbind,
  list(
    summary_rows(stan_draws, "Stan"),
    summary_rows(shared_draws, "SMC shared_bank"),
    summary_rows(collapsed_bridge_draws, "SMC collapsed_bridge"),
    summary_rows(collapsed_draws, "SMC collapsed")
  )
)

timing_summary <- rbind(
  data.frame(method = "Stan", stage = "full", elapsed = stan_stage$elapsed, user = stan_stage$user, system = stan_stage$system),
  data.frame(method = "SMC shared_bank", stage = "local_prefit", elapsed = shared_local_stage$elapsed, user = shared_local_stage$user, system = shared_local_stage$system),
  data.frame(method = "SMC shared_bank", stage = "outer", elapsed = shared_outer_stage$elapsed, user = shared_outer_stage$user, system = shared_outer_stage$system),
  data.frame(method = "SMC shared_bank", stage = "total", elapsed = shared_local_stage$elapsed + shared_outer_stage$elapsed, user = shared_local_stage$user + shared_outer_stage$user, system = shared_local_stage$system + shared_outer_stage$system),
  data.frame(method = "SMC collapsed_bridge", stage = "full", elapsed = collapsed_bridge_stage$elapsed, user = collapsed_bridge_stage$user, system = collapsed_bridge_stage$system),
  data.frame(method = "SMC collapsed_subject", stage = "full", elapsed = collapsed_subject_stage$elapsed, user = collapsed_subject_stage$user, system = collapsed_subject_stage$system)
)

plot_file <- file.path("benchmarks", "normal_normal_meanvar_hierarchy_posteriors.png")
png(plot_file, width = 1100, height = 900)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

make_density_or_spike <- function(x, color, spread_ref) {
  x <- as.numeric(x[is.finite(x)])
  if (!length(x)) stop("Need at least one finite draw.")
  sd_x <- stats::sd(x)
  spike_tol <- max(1e-4, 1e-3 * max(spread_ref, 1e-3))
  if (!is.finite(sd_x) || sd_x < spike_tol) {
    return(list(type = "spike", x = mean(x), color = color, peak = Inf))
  }
  d <- density(x, n = 512)
  list(type = "density", density = d, color = color, peak = max(d$y))
}

plot_one <- function(stan_x, shared_x, bridge_x, collapsed_x, truth, main, xlab) {
  spread_ref <- stats::IQR(c(stan_x, shared_x, bridge_x, collapsed_x), na.rm = TRUE)
  curves <- list(
    make_density_or_spike(stan_x, "black", spread_ref),
    make_density_or_spike(shared_x, "steelblue", spread_ref),
    make_density_or_spike(bridge_x, "darkorange3", spread_ref),
    make_density_or_spike(collapsed_x, "firebrick", spread_ref)
  )
  density_x <- unlist(lapply(curves, function(curve) if (identical(curve$type, "density")) curve$density$x else curve$x))
  density_peaks <- vapply(curves, function(curve) if (identical(curve$type, "density")) curve$peak else NA_real_, numeric(1))
  peak_ref <- max(density_peaks, na.rm = TRUE)
  if (!is.finite(peak_ref) || peak_ref <= 0) peak_ref <- 1
  xlim <- range(c(density_x, truth), na.rm = TRUE)
  ylim <- c(0, 1.05 * peak_ref)

  plot(NA_real_, NA_real_, xlim = xlim, ylim = ylim, main = main, xlab = xlab, ylab = "Density")
  for (curve in curves) {
    if (identical(curve$type, "density")) {
      lines(curve$density, lwd = 2, col = curve$color)
    } else {
      segments(curve$x, 0, curve$x, ylim[2], lwd = 2, col = curve$color)
    }
  }
  abline(v = truth, lwd = 2, col = "darkgreen", lty = 2)
  legend(
    "topright",
    legend = c("Stan", "SMC shared_bank", "SMC collapsed_bridge", "SMC collapsed", "Truth"),
    col = c("black", "steelblue", "darkorange3", "firebrick", "darkgreen"),
    lwd = c(2, 2, 2, 2, 2),
    lty = c(1, 1, 1, 1, 2),
    bty = "n",
    cex = 0.85
  )
}

plot_one(
  stan_x = stan_draws$mu_mean,
  shared_x = shared_draws$mu_mean,
  bridge_x = collapsed_bridge_draws$mu_mean,
  collapsed_x = collapsed_draws$mu_mean,
  truth = true_phi["mu_mean"],
  main = expression(mu[mean]),
  xlab = expression(mu[mean])
)
plot_one(
  stan_x = stan_draws$tau2_mean,
  shared_x = shared_draws$tau2_mean,
  bridge_x = collapsed_bridge_draws$tau2_mean,
  collapsed_x = collapsed_draws$tau2_mean,
  truth = exp(true_phi["log_tau2_mean"]),
  main = expression(tau[mean]^2),
  xlab = expression(tau[mean]^2)
)
plot_one(
  stan_x = stan_draws$mu_log_var,
  shared_x = shared_draws$mu_log_var,
  bridge_x = collapsed_bridge_draws$mu_log_var,
  collapsed_x = collapsed_draws$mu_log_var,
  truth = true_phi["mu_log_var"],
  main = expression(mu[log(var)]),
  xlab = expression(mu[log(var)])
)
plot_one(
  stan_x = stan_draws$tau2_log_var,
  shared_x = shared_draws$tau2_log_var,
  bridge_x = collapsed_bridge_draws$tau2_log_var,
  collapsed_x = collapsed_draws$tau2_log_var,
  truth = exp(true_phi["log_tau2_log_var"]),
  main = expression(tau[log(var)]^2),
  xlab = expression(tau[log(var)]^2)
)
dev.off()

results <- list(
  simulation = list(
    S = S,
    N_trials = N_trials,
    true_phi = true_phi,
    subject_mean_true = subject_mean_true,
    subject_log_var_true = subject_log_var_true,
    y = y
  ),
  prior = list(
    mean = phi_prior_mean,
    sd = phi_prior_sd
  ),
  stan_draws = stan_draws,
  shared_draws = shared_draws,
  collapsed_bridge_draws = collapsed_bridge_draws,
  collapsed_draws = collapsed_draws,
  settings = list(
    stage_versions = stage_versions,
    mc_cores = mc_cores,
    stan_chains = stan_chains,
    stan_iter_warmup = stan_iter_warmup,
    stan_iter_sampling = stan_iter_sampling,
    shared_local_particles = shared_local_particles,
    shared_local_rounds = shared_local_rounds,
    shared_outer_particles = shared_outer_particles,
    shared_bank_particles = shared_bank_particles,
    shared_outer_rounds = shared_outer_rounds,
    collapsed_outer_particles = collapsed_outer_particles,
    collapsed_bridge_outer_particles = collapsed_bridge_outer_particles,
    collapsed_local_particles = collapsed_local_particles,
    collapsed_local_rounds = collapsed_local_rounds,
    collapsed_rejuvenation = collapsed_rejuvenation,
    collapsed_rejuvenate_every = collapsed_rejuvenate_every,
    collapsed_rejuvenate_after_resample = collapsed_rejuvenate_after_resample,
    collapsed_n_rejuvenation_moves = collapsed_n_rejuvenation_moves,
    collapsed_bridge_steps = collapsed_bridge_steps,
    collapsed_bridge_rejuvenation = collapsed_bridge_rejuvenation,
    collapsed_bridge_rejuvenate_every = collapsed_bridge_rejuvenate_every,
    collapsed_bridge_rejuvenate_after_resample = collapsed_bridge_rejuvenate_after_resample,
    collapsed_bridge_n_rejuvenation_moves = collapsed_bridge_n_rejuvenation_moves,
    subject_order = subject_order
  ),
  cache = list(
    force_recompute = force_recompute,
    cache_dir = cache_dir,
    stan = list(file = stan_stage$cache_file),
    shared_local = list(file = shared_local_stage$cache_file),
    shared_outer = list(file = shared_outer_stage$cache_file),
    collapsed_subject = list(file = collapsed_subject_stage$cache_file),
    collapsed_bridge = list(file = collapsed_bridge_stage$cache_file)
  ),
  stan_summary = stan_summary,
  timing_summary = timing_summary,
  shared_fit = shared_fit,
  collapsed_bridge_fit = collapsed_bridge_fit,
  collapsed_fit = collapsed_fit,
  comparison_summary = comparison_summary,
  plot_file = plot_file
)

saveRDS(results, results_file)

cat("NORMAL_NORMAL_MEANVAR_HIERARCHY_OK\n")
cat(sprintf("results=%s\n", results_file))
cat(sprintf("plot=%s\n", plot_file))
cat(sprintf("shared_final_lambda=%.4f\n", tail(shared_fit$meta$lambda_hist, 1L)))
cat(sprintf("collapsed_bridge_stages=%d\n", collapsed_bridge_fit$meta$n_stages))
cat(sprintf("collapsed_bridge_local_runs=%d\n", collapsed_bridge_fit$meta$n_local_runs))
cat(sprintf("collapsed_subjects=%d\n", collapsed_fit$meta$n_subjects))
cat(sprintf("collapsed_local_runs=%d\n", collapsed_fit$meta$n_local_runs))
cat("\nTiming summary (seconds):\n")
print(timing_summary, row.names = FALSE)
cat("\nComparison summary:\n")
print(comparison_summary, row.names = FALSE)
