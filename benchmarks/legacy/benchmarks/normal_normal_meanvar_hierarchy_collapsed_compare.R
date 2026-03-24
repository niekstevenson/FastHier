#!/usr/bin/env Rscript

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[1L]))
} else {
  normalizePath("benchmarks/normal_normal_meanvar_hierarchy_collapsed_compare.R")
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

results_file <- file.path("samples", "normal_normal_meanvar_hierarchy_collapsed_compare_results.rds")
cache_dir <- file.path("samples", "normal_normal_meanvar_hierarchy_collapsed_compare_cache")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

force_recompute <- FALSE
stage_versions <- list(
  stan = 1L,
  collapsed_subject_recompute = 5L,
  collapsed_bridge_fixed = 4L,
  collapsed_bridge_adaptive = 4L
)

set.seed(20260314L)

S <- 15L
N_trials <- 100L
mc_cores <- min(4L, parallel::detectCores(logical = TRUE))

stan_chains <- 4L
stan_iter_warmup <- 1000L
stan_iter_sampling <- 1000L

collapsed_outer_particles <- 64L
collapsed_bridge_particles <- 64L
collapsed_local_particles <- 192L
collapsed_local_rounds <- 30L
collapsed_resample_threshold <- 0.5
collapsed_subject_rejuvenate_every <- 3L
collapsed_subject_moves <- 1L
collapsed_bridge_steps <- 4L
collapsed_bridge_rejuvenate_every <- 6L
collapsed_bridge_moves <- 1L
collapsed_adaptive_ess_target <- 0.55
collapsed_adaptive_lambda_tol <- 0.02
collapsed_adaptive_max_steps <- 12L
collapsed_adaptive_max_bisect <- 6L
include_adaptive_bridge <- FALSE

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
  hierarchical = unname(tools::md5sum("hierarchical_smc.R")),
  collapsed = unname(tools::md5sum("hierarchical_collapsed_smc.R")),
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

collapsed_local_control <- list(
  max_rounds = collapsed_local_rounds,
  n_mcmc_moves = 1L,
  gss_enable = FALSE,
  hist_mix_enable = FALSE,
  da_enable = FALSE
)

set.seed(777L)
subject_order <- sample.int(S)

stan_stage <- load_or_run_cached(
  stage = "stan",
  settings = list(
    stage_version = stage_versions$stan,
    code = list(stan = code_fingerprint$stan),
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

run_collapsed_subject_stage <- function(label, rejuvenation, rejuvenate_every, rejuvenate_after_resample, base_seed) {
  load_or_run_cached(
    stage = label,
    settings = list(
      stage_version = stage_versions[[label]],
      code = code_fingerprint,
      simulation = simulation_fingerprint,
      subject_order = subject_order,
      outer_particles = collapsed_outer_particles,
      local_particles = collapsed_local_particles,
      local_control = collapsed_local_control,
      rejuvenation = rejuvenation,
      rejuvenate_every = rejuvenate_every,
      rejuvenate_after_resample = rejuvenate_after_resample,
      n_rejuvenation_moves = collapsed_subject_moves,
      auxiliary_proposal = "correlated",
      auxiliary_correlation = 0.95,
      warm_start_mode = "nearest"
    ),
    run_fn = function() {
      collapsed_subject_smc(
        data_list = data_list,
        loglik_fn = loglik_fn,
        prior = phi_prior,
        gaussian_map_fn = gaussian_map_fn,
        d_theta = length(theta_names),
        N = collapsed_outer_particles,
        local_particles = collapsed_local_particles,
        resample_threshold = collapsed_resample_threshold,
        subject_order = subject_order,
        rejuvenation = rejuvenation,
        rejuvenate_every = rejuvenate_every,
        rejuvenate_after_resample = rejuvenate_after_resample,
        n_rejuvenation_moves = collapsed_subject_moves,
        auxiliary_proposal = "correlated",
        auxiliary_correlation = 0.95,
        warm_start_mode = "nearest",
        theta_names = theta_names,
        local_n_cores = 1L,
        outer_n_cores = mc_cores,
        base_seed = base_seed,
        assume_unbiased_local_evidence = TRUE,
        evidence_cache = TRUE,
        store_local_fits = FALSE,
        verbose = TRUE,
        local_smc_control = collapsed_local_control
      )
    },
    force = force_recompute
  )
}

run_collapsed_bridge_stage <- function(label,
                                       adaptive_bridge,
                                       n_bridge_steps = collapsed_bridge_steps,
                                       base_seed,
                                       rejuvenate_every = collapsed_bridge_rejuvenate_every) {
  load_or_run_cached(
    stage = label,
    settings = list(
      stage_version = stage_versions[[label]],
      code = code_fingerprint,
      simulation = simulation_fingerprint,
      subject_order = subject_order,
      outer_particles = collapsed_bridge_particles,
      local_particles = collapsed_local_particles,
      local_control = collapsed_local_control,
      adaptive_bridge = adaptive_bridge,
      n_bridge_steps = n_bridge_steps,
      rejuvenation = "recompute",
      rejuvenate_every = rejuvenate_every,
      rejuvenate_after_resample = TRUE,
      n_rejuvenation_moves = collapsed_bridge_moves,
      auxiliary_proposal = "correlated",
      auxiliary_correlation = 0.95,
      warm_start_mode = "nearest",
      adaptive_ess_target = collapsed_adaptive_ess_target,
      adaptive_lambda_tol = collapsed_adaptive_lambda_tol,
      adaptive_max_steps = collapsed_adaptive_max_steps,
      adaptive_max_bisect = collapsed_adaptive_max_bisect
    ),
    run_fn = function() {
      collapsed_bridge_smc(
        data_list = data_list,
        loglik_fn = loglik_fn,
        prior = phi_prior,
        gaussian_map_fn = gaussian_map_fn,
        d_theta = length(theta_names),
        N = collapsed_bridge_particles,
        local_particles = collapsed_local_particles,
        resample_threshold = collapsed_resample_threshold,
        subject_order = subject_order,
        n_bridge_steps = n_bridge_steps,
        rejuvenation = "recompute",
        rejuvenate_every = rejuvenate_every,
        rejuvenate_after_resample = TRUE,
        n_rejuvenation_moves = collapsed_bridge_moves,
        auxiliary_proposal = "correlated",
        auxiliary_correlation = 0.95,
        warm_start_mode = "nearest",
        theta_names = theta_names,
        local_n_cores = 1L,
        outer_n_cores = mc_cores,
        base_seed = base_seed,
        assume_unbiased_local_evidence = TRUE,
        evidence_cache = TRUE,
        store_local_fits = FALSE,
        adaptive_bridge = adaptive_bridge,
        adaptive_ess_target = collapsed_adaptive_ess_target,
        adaptive_lambda_tol = collapsed_adaptive_lambda_tol,
        adaptive_max_steps_per_subject = collapsed_adaptive_max_steps,
        adaptive_max_bisect = collapsed_adaptive_max_bisect,
        verbose = TRUE,
        local_smc_control = collapsed_local_control
      )
    },
    force = force_recompute
  )
}

collapsed_subject_recompute_stage <- run_collapsed_subject_stage(
  label = "collapsed_subject_recompute",
  rejuvenation = "recompute",
  rejuvenate_every = collapsed_subject_rejuvenate_every,
  rejuvenate_after_resample = TRUE,
  base_seed = 404L
)

collapsed_bridge_fixed_stage <- run_collapsed_bridge_stage(
  label = "collapsed_bridge_fixed",
  adaptive_bridge = FALSE,
  n_bridge_steps = collapsed_bridge_steps,
  base_seed = 505L
)

collapsed_bridge_adaptive_stage <- if (isTRUE(include_adaptive_bridge)) {
  run_collapsed_bridge_stage(
    label = "collapsed_bridge_adaptive",
    adaptive_bridge = TRUE,
    n_bridge_steps = collapsed_bridge_steps,
    base_seed = 506L
  )
} else {
  NULL
}

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

method_specs <- list(
  list(
    key = "stan",
    label = "Stan",
    draws = stan_draws,
    stage = stan_stage,
    color = "black"
  ),
  list(
    key = "collapsed_subject_recompute",
    label = "SMC collapsed_subject recompute",
    draws = phi_draws_to_named_list(
      resample_weighted_rows(collapsed_subject_recompute_stage$value$phi, collapsed_subject_recompute_stage$value$w, n_plot_draws, seed = 1001L)
    ),
    stage = collapsed_subject_recompute_stage,
    color = "firebrick3"
  ),
  list(
    key = "collapsed_bridge_fixed",
    label = "SMC collapsed_bridge fixed",
    draws = phi_draws_to_named_list(
      resample_weighted_rows(collapsed_bridge_fixed_stage$value$phi, collapsed_bridge_fixed_stage$value$w, n_plot_draws, seed = 1002L)
    ),
    stage = collapsed_bridge_fixed_stage,
    color = "steelblue4"
  )
)
if (isTRUE(include_adaptive_bridge)) {
  method_specs[[length(method_specs) + 1L]] <- list(
    key = "collapsed_bridge_adaptive",
    label = "SMC collapsed_bridge adaptive",
    draws = phi_draws_to_named_list(
      resample_weighted_rows(collapsed_bridge_adaptive_stage$value$phi, collapsed_bridge_adaptive_stage$value$w, n_plot_draws, seed = 1003L)
    ),
    stage = collapsed_bridge_adaptive_stage,
    color = "seagreen4"
  )
}

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
  lapply(method_specs, function(spec) summary_rows(spec$draws, spec$label))
)

timing_summary <- do.call(
  rbind,
  lapply(method_specs, function(spec) {
    data.frame(
      method = spec$label,
      stage = "full",
      elapsed = spec$stage$elapsed,
      user = spec$stage$user,
      system = spec$stage$system
    )
  })
)

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

plot_param <- function(specs, param, truth, main, xlab) {
  all_x <- unlist(lapply(specs, function(spec) spec$draws[[param]]), use.names = FALSE)
  spread_ref <- stats::IQR(all_x, na.rm = TRUE)
  curves <- lapply(specs, function(spec) {
    make_density_or_spike(spec$draws[[param]], spec$color, spread_ref)
  })
  density_x <- unlist(lapply(curves, function(curve) if (identical(curve$type, "density")) curve$density$x else curve$x))
  density_peaks <- vapply(curves, function(curve) if (identical(curve$type, "density")) curve$peak else NA_real_, numeric(1))
  peak_ref <- max(density_peaks, na.rm = TRUE)
  if (!is.finite(peak_ref) || peak_ref <= 0) peak_ref <- 1
  xlim <- range(c(density_x, truth), na.rm = TRUE)
  ylim <- c(0, 1.05 * peak_ref)

  plot(NA_real_, NA_real_, xlim = xlim, ylim = ylim, main = main, xlab = xlab, ylab = "Density")
  for (k in seq_along(curves)) {
    curve <- curves[[k]]
    if (identical(curve$type, "density")) {
      lines(curve$density, lwd = 2, col = curve$color)
    } else {
      segments(curve$x, 0, curve$x, ylim[2], lwd = 2, col = curve$color)
    }
  }
  abline(v = truth, lwd = 2, col = "darkgreen", lty = 2)
  legend(
    "topright",
    legend = c(vapply(specs, `[[`, character(1), "label"), "Truth"),
    col = c(vapply(specs, `[[`, character(1), "color"), "darkgreen"),
    lwd = c(rep(2, length(specs)), 2),
    lty = c(rep(1, length(specs)), 2),
    bty = "n",
    cex = 0.8
  )
}

plot_file <- file.path("benchmarks", "normal_normal_meanvar_hierarchy_collapsed_posteriors.png")
png(plot_file, width = 1200, height = 900)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
plot_param(method_specs, "mu_mean", true_phi["mu_mean"], expression(mu[mean]), expression(mu[mean]))
plot_param(method_specs, "tau2_mean", exp(true_phi["log_tau2_mean"]), expression(tau[mean]^2), expression(tau[mean]^2))
plot_param(method_specs, "mu_log_var", true_phi["mu_log_var"], expression(mu[log(var)]), expression(mu[log(var)]))
plot_param(method_specs, "tau2_log_var", exp(true_phi["log_tau2_log_var"]), expression(tau[log(var)]^2), expression(tau[log(var)]^2))
dev.off()

results <- list(
  simulation = list(
    S = S,
    N_trials = N_trials,
    true_phi = true_phi,
    subject_mean_true = subject_mean_true,
    subject_log_var_true = subject_log_var_true,
    y = y,
    subject_order = subject_order
  ),
  prior = list(
    mean = phi_prior_mean,
    sd = phi_prior_sd
  ),
  settings = list(
    stage_versions = stage_versions,
    mc_cores = mc_cores,
    stan_chains = stan_chains,
    stan_iter_warmup = stan_iter_warmup,
    stan_iter_sampling = stan_iter_sampling,
    collapsed_outer_particles = collapsed_outer_particles,
    collapsed_bridge_particles = collapsed_bridge_particles,
    collapsed_local_particles = collapsed_local_particles,
    collapsed_local_control = collapsed_local_control,
    collapsed_resample_threshold = collapsed_resample_threshold,
    collapsed_subject_rejuvenate_every = collapsed_subject_rejuvenate_every,
    collapsed_bridge_steps = collapsed_bridge_steps,
    collapsed_bridge_rejuvenate_every = collapsed_bridge_rejuvenate_every,
    collapsed_adaptive_ess_target = collapsed_adaptive_ess_target,
    collapsed_adaptive_lambda_tol = collapsed_adaptive_lambda_tol,
    collapsed_adaptive_max_steps = collapsed_adaptive_max_steps,
    collapsed_adaptive_max_bisect = collapsed_adaptive_max_bisect,
    include_adaptive_bridge = include_adaptive_bridge
  ),
  stan_draws = stan_draws,
  method_draws = setNames(lapply(method_specs, `[[`, "draws"), vapply(method_specs, `[[`, character(1), "key")),
  comparison_summary = comparison_summary,
  timing_summary = timing_summary,
  plot_file = plot_file,
  cache_files = setNames(lapply(method_specs, function(spec) spec$stage$cache_file), vapply(method_specs, `[[`, character(1), "key")),
  stan_summary = stan_summary
)

saveRDS(results, results_file)

cat(sprintf("Saved results to %s\n", results_file))
cat(sprintf("Saved plot to %s\n", plot_file))
print(comparison_summary)
print(timing_summary)
