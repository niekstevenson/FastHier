rm(list = ls())

suppressPackageStartupMessages({
  library(parallel)
})

if (!file.exists("smc_core.R")) {
  stop("Run this script from the FastHierarchical repository root.")
}

detected_cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
if (!is.finite(detected_cores) || detected_cores < 1L) detected_cores <- 1L

base_seed <- 20260521L
cores <- as.integer(min(4L, detected_cores))
rho_ladder <- c(0.0005, 0.001, 0.005, 0.01, 0.02, 0.05, 0.10, 0.20, 0.35, 0.50, 0.75, 1.00)
outer_particles <- 1000L
sketch_starts <- 16L
correction_particles <- 600L
refinement_rounds <- 2L
refinement_points_per_round <- 4L

stan_results_file <- file.path("benchmarks", "samples", "shifted_gamma_hierarchy_stan_results.rds")
results_file <- file.path("benchmarks", "results", "shifted_gamma_sketch_support_results.rds")
plot_file <- file.path("benchmarks", "results", "shifted_gamma_sketch_support_rho_ladder.png")

dir.create(file.path("benchmarks", "results"), showWarnings = FALSE, recursive = TRUE)

source("smc_core.R")
source("population_models.R")
source("outer_population_smc.R")
source("utilities.R")
source("local_likelihood_sketches.R")
source("theta_support_pathfinder.R")
source("theta_smc_corrections.R")
source("theta_correction_surrogate.R")

bundle <- readRDS(stan_results_file)
y <- bundle$data$y
data_list <- lapply(seq_len(nrow(y)), function(i) y[i, ])

alpha_names <- c("eta_shape", "eta_scale", "eta_shift")
m0 <- stats::setNames(as.numeric(bundle$priors$m0), alpha_names)
s0 <- stats::setNames(as.numeric(bundle$priors$s0), alpha_names)
a0 <- stats::setNames(as.numeric(bundle$priors$a0), alpha_names)
b0 <- stats::setNames(as.numeric(bundle$priors$b0), alpha_names)

loglik_shifted_gamma <- function(Theta, y_i) {
  Theta <- as.matrix(Theta)
  colnames(Theta) <- alpha_names

  eps <- 1e-9
  shape <- exp(Theta[, "eta_shape"]) + eps
  scale <- exp(Theta[, "eta_scale"]) + eps
  shift <- exp(Theta[, "eta_shift"]) + eps
  min_y <- min(y_i)

  out <- rep(-1e12, nrow(Theta))
  ok <- shift < min_y
  if (!any(ok)) return(out)

  for (i in which(ok)) {
    out[i] <- sum(stats::dgamma(y_i - shift[i], shape = shape[i], scale = scale[i], log = TRUE))
  }
  out[!is.finite(out)] <- -1e12
  out
}

population_model <- make_population_model_diag_gaussian(
  alpha_names = alpha_names,
  mean_prior_mean = m0,
  mean_prior_var = s0,
  sigma2_prior_shape = a0,
  sigma2_prior_rate = b0,
  label = "shifted_gamma_hierarchy"
)

initial_theta <- matrix(c(m0, log(b0 / (a0 - 1))), nrow = 1L)
colnames(initial_theta) <- population_model$hyper_names

stan_draws <- data.frame(
  mu_shape = bundle$draws$mu[, 1L],
  mu_scale = bundle$draws$mu[, 2L],
  mu_shift = bundle$draws$mu[, 3L],
  sigma2_shape = bundle$draws$sigma2[, 1L],
  sigma2_scale = bundle$draws$sigma2[, 2L],
  sigma2_shift = bundle$draws$sigma2[, 3L],
  check.names = FALSE
)

theta_draws_to_population_draws <- function(fit, seed) {
  parts <- smc_posteriors(
    fit,
    n_draws = nrow(stan_draws),
    seed = seed,
    population_model = population_model
  )
  mu <- as.data.frame(parts$mu, check.names = FALSE)
  sigma2 <- as.data.frame(parts$sigma2, check.names = FALSE)
  colnames(mu) <- alpha_names
  colnames(sigma2) <- alpha_names
  data.frame(
    mu_shape = mu[, "eta_shape"],
    mu_scale = mu[, "eta_scale"],
    mu_shift = mu[, "eta_shift"],
    sigma2_shape = sigma2[, "eta_shape"],
    sigma2_scale = sigma2[, "eta_scale"],
    sigma2_shift = sigma2[, "eta_shift"],
    check.names = FALSE
  )
}

weighted_quant <- function(x, probs) {
  stats::quantile(as.numeric(x), probs = probs, na.rm = TRUE, names = FALSE)
}

support_recall_table <- function(stan, support_draws) {
  params <- colnames(stan)
  rows <- lapply(params, function(param) {
    stan_q <- weighted_quant(stan[, param], c(0.01, 0.05, 0.50, 0.95, 0.99))
    support_q <- weighted_quant(support_draws[, param], c(0.005, 0.995))
    data.frame(
      parameter = param,
      stan_q01 = stan_q[1L],
      stan_q05 = stan_q[2L],
      stan_q50 = stan_q[3L],
      stan_q95 = stan_q[4L],
      stan_q99 = stan_q[5L],
      support_q005 = support_q[1L],
      support_q995 = support_q[2L],
      support_min = min(support_draws[, param], na.rm = TRUE),
      support_max = max(support_draws[, param], na.rm = TRUE),
      covers_stan_q01 = support_q[1L] <= stan_q[1L],
      covers_stan_q99 = support_q[2L] >= stan_q[5L],
      left_gap = max(support_q[1L] - stan_q[1L], 0),
      right_gap = max(stan_q[5L] - support_q[2L], 0),
      check.names = FALSE
    )
  })
  do.call(rbind, rows)
}

rho_quantile_table <- function(rho_draws) {
  do.call(rbind, lapply(names(rho_draws), function(label) {
    x <- rho_draws[[label]]
    do.call(rbind, lapply(colnames(x), function(param) {
      q <- weighted_quant(x[, param], c(0.01, 0.50, 0.99))
      data.frame(
        rho = as.numeric(sub("^rho_", "", label)),
        parameter = param,
        q01 = q[1L],
        q50 = q[2L],
        q99 = q[3L],
        check.names = FALSE
      )
    }))
  }))
}

plot_support_ladder <- function(stan, rho_draws, file) {
  cols <- grDevices::hcl.colors(length(rho_draws), palette = "Dark 3")
  grDevices::png(file, width = 1600, height = 1000)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(2, 3), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
  for (param in colnames(stan)) {
    d_stan <- stats::density(stan[, param])
    d_rho <- lapply(rho_draws, function(x) stats::density(x[, param]))
    xlim <- range(c(d_stan$x, unlist(lapply(d_rho, `[[`, "x"))), finite = TRUE)
    ylim <- c(0, 1.05 * max(c(d_stan$y, unlist(lapply(d_rho, `[[`, "y"))), finite = TRUE))
    graphics::plot(
      d_stan,
      lwd = 2.5,
      col = "black",
      xlim = xlim,
      ylim = ylim,
      main = param,
      xlab = param,
      ylab = "Density"
    )
    for (j in seq_along(d_rho)) {
      graphics::lines(d_rho[[j]], col = cols[j], lwd = 1.4)
    }
    graphics::legend(
      "topright",
      legend = c("Stan", paste0("rho=", format(rho_ladder, trim = TRUE))),
      col = c("black", cols),
      lwd = c(2.5, rep(1.4, length(cols))),
      bty = "n",
      cex = 0.70
    )
  }
  graphics::mtext("Tempered analytic sketch support ladder", outer = TRUE, cex = 1.2)
}

cat(sprintf("Sketch support and correction benchmark | subjects=%d | cores=%d | particles=%d\n",
            length(data_list), cores, outer_particles))
cat("Rho ladder:", paste(rho_ladder, collapse = ", "), "\n")

start_time <- Sys.time()

sketch_set <- fit_local_likelihood_sketches(
  data_list = data_list,
  loglik_fn = loglik_shifted_gamma,
  alpha_names = alpha_names,
  population_model = population_model,
  theta_reference = initial_theta,
  n_starts = sketch_starts,
  start_scale = 4,
  n_jobs = cores,
  seed = base_seed
)

factor_sets <- lapply(rho_ladder, function(rho) {
  build_local_likelihood_sketch_factor_set(
    sketch_set = sketch_set,
    population_model = population_model,
    rho = rho
  )
})
names(factor_sets) <- paste0("rho_", format(rho_ladder, trim = TRUE))

fits <- vector("list", length(rho_ladder))
names(fits) <- names(factor_sets)
timing <- data.frame(
  rho = rho_ladder,
  elapsed_sec = NA_real_,
  rounds = NA_integer_,
  update_ess = NA_real_,
  accept_rate = NA_real_,
  check.names = FALSE
)

first_time <- system.time({
  fits[[1L]] <- outer_population_smc(
    factor_sets[[1L]],
    N = outer_particles,
    resample_threshold = 0.5,
    n_mcmc_moves = 2L,
    min_mcmc_moves = 1L,
    max_rounds = 80L,
    rw_scale_init = 0.8,
    cess_target = 0.60,
    n_cores = cores,
    seed = base_seed,
    verbose = TRUE
  )
})
timing$elapsed_sec[1L] <- first_time[["elapsed"]]
timing$rounds[1L] <- fits[[1L]]$meta$rounds
timing$accept_rate[1L] <- tail(fits[[1L]]$meta$accept_hist, 1L)

for (j in seq_along(rho_ladder)[-1L]) {
  update_time <- system.time({
    fits[[j]] <- update_outer_population_fit(
      fit = fits[[j - 1L]],
      old_factor_set = factor_sets[[j - 1L]],
      new_factor_set = factor_sets[[j]],
      n_mcmc_moves = 2L,
      min_mcmc_moves = 1L,
      resample_threshold = 0.5,
      n_cores = cores,
      seed = base_seed + j,
      verbose = TRUE
    )
  })
  update_info <- tail(fits[[j]]$meta$factor_update, 1L)[[1L]]
  timing$elapsed_sec[j] <- update_time[["elapsed"]]
  timing$rounds[j] <- 0L
  timing$update_ess[j] <- update_info$ess
  timing$accept_rate[j] <- update_info$accept_rate
}

rho_draws <- lapply(seq_along(fits), function(j) {
  theta_draws_to_population_draws(fits[[j]], seed = base_seed + 100L + j)
})
names(rho_draws) <- names(fits)

support_draws <- do.call(rbind, rho_draws)
support_table <- support_recall_table(stan_draws, support_draws)
rho_quantiles <- rho_quantile_table(rho_draws)

support_theta <- do.call(rbind, lapply(fits, `[[`, "theta"))
support_w <- unlist(lapply(fits, function(fit) normalize_weights(fit$w) / length(fits)), use.names = FALSE)

pathfinder <- fit_theta_pathfinder_support(
  factor_set = factor_sets[[length(factor_sets)]],
  support_theta = support_theta,
  support_w = support_w,
  n_starts = 24L,
  coverage_probs = c(0.50, 0.90, 0.99),
  maxit = 600L,
  unique_distance = 0.25,
  max_sigma_directions = population_model$hyper_dim,
  n_jobs = cores,
  seed = base_seed + 500L,
  verbose = TRUE
)
pathfinder$candidates <- theta_pathfinder_candidate_distance(
  pathfinder,
  reference_theta = support_theta,
  reference_w = support_w
)
pathfinder$candidates <- theta_pathfinder_score_candidates(
  pathfinder,
  candidates = pathfinder$candidates,
  n_cores = cores
)

pathfinder_summary <- data.frame(
  starts = nrow(pathfinder$starts),
  finite_paths = sum(is.finite(pathfinder$path_summary$logposterior)),
  converged_paths = sum(pathfinder$path_summary$convergence == 0L),
  unique_paths = length(pathfinder$unique_path_ids),
  usable_paths = length(pathfinder$usable_path_ids),
  hessian_valid_unique = sum(pathfinder$path_summary$hessian_valid[pathfinder$path_summary$unique]),
  components = length(pathfinder$components),
  candidates = nrow(pathfinder$candidates),
  nonduplicate_candidates = sum(pathfinder$candidates$nonduplicate_support),
  support_distance_cutoff = unique(pathfinder$candidates$support_distance_cutoff)[1L],
  max_nearest_support_distance = max(pathfinder$candidates$nearest_support_distance, na.rm = TRUE),
  check.names = FALSE
)

correction_design <- select_theta_correction_points(
  final_fit = fits[[length(fits)]],
  population_model = population_model,
  support_theta = support_theta,
  support_w = support_w,
  pathfinder = pathfinder,
  tail_probs = c(0.05, 0.95)
)

correction_time <- system.time({
  initial_correction_set <- fit_theta_smc_corrections(
    design = correction_design,
    factor_set = factor_sets[[length(factor_sets)]],
    data_list = data_list,
    loglik_fn = loglik_shifted_gamma,
    M = correction_particles,
    n_jobs = cores,
    local_n_cores = 1L,
    base_seed = base_seed + 1000L,
    replicate_count = 2L,
    verbose = TRUE
  )
})

summarize_correction_set <- function(set, elapsed_sec) {
  data.frame(
    points = nrow(set$design),
    point_evaluations = nrow(set$observations),
    local_smc_runs = nrow(set$by_subject),
    particles_per_local_run = correction_particles,
    elapsed_sec = elapsed_sec,
    delta_min = min(set$observations$delta_total),
    delta_median = stats::median(set$observations$delta_total),
    delta_max = max(set$observations$delta_total),
    smc_relative_min = min(set$observations$relative_logposterior_smc, na.rm = TRUE),
    smc_relative_median = stats::median(set$observations$relative_logposterior_smc, na.rm = TRUE),
    max_mcse_total = max(set$observations$mcse_total, na.rm = TRUE),
    max_replicate_delta_sd = max(set$replicate_summary$delta_sd, na.rm = TRUE),
    max_degenerate_subjects = max(set$observations$degenerate_subjects, na.rm = TRUE),
    min_final_ess_frac = min(set$observations$min_final_ess_frac, na.rm = TRUE),
    check.names = FALSE
  )
}

initial_correction_summary <- summarize_correction_set(initial_correction_set, correction_time[["elapsed"]])

initial_correction_surrogate <- fit_theta_correction_surrogate(
  correction_set = initial_correction_set,
  support_theta = support_theta,
  support_w = support_w
)

correction_set <- initial_correction_set
correction_surrogate <- initial_correction_surrogate
refinement_designs <- list()
refinement_correction_sets <- list()
refinement_rows <- list()
refinement_elapsed <- 0
for (round_id in seq_len(refinement_rounds)) {
  design_round <- select_theta_correction_refinement_points(
    correction_surrogate,
    n_points = refinement_points_per_round,
    next_point_id = max(correction_set$design$point_id) + 1L
  )
  if (!nrow(design_round)) break
  design_round$refinement_round <- round_id
  refinement_designs[[round_id]] <- design_round
  time_round <- system.time({
    set_round <- fit_theta_smc_corrections(
      design = design_round,
      factor_set = factor_sets[[length(factor_sets)]],
      data_list = data_list,
      loglik_fn = loglik_shifted_gamma,
      M = correction_particles,
      n_jobs = cores,
      local_n_cores = 1L,
      base_seed = base_seed + 2000L + 100L * round_id,
      replicate_ids = design_round$point_id,
      replicate_count = 2L,
      verbose = TRUE
    )
  })
  refinement_elapsed <- refinement_elapsed + time_round[["elapsed"]]
  refinement_correction_sets[[round_id]] <- set_round
  refinement_rows[[round_id]] <- data.frame(
    refinement_round = round_id,
    points = nrow(design_round),
    point_evaluations = nrow(set_round$observations),
    local_smc_runs = nrow(set_round$by_subject),
    particles_per_local_run = correction_particles,
    elapsed_sec = time_round[["elapsed"]],
    check.names = FALSE
  )
  correction_set <- combine_theta_smc_correction_sets(correction_set, set_round)
  correction_surrogate <- fit_theta_correction_surrogate(
    correction_set = correction_set,
    support_theta = support_theta,
    support_w = support_w
  )
}

refinement_design <- if (length(refinement_designs)) do.call(rbind, refinement_designs) else {
  structure(data.frame(), class = c("theta_correction_design", "data.frame"))
}
refinement_summary <- if (length(refinement_rows)) do.call(rbind, refinement_rows) else {
  data.frame(
    refinement_round = integer(0),
    points = integer(0),
    point_evaluations = integer(0),
    local_smc_runs = integer(0),
    particles_per_local_run = integer(0),
    elapsed_sec = numeric(0)
  )
}
refinement_correction_set <- if (length(refinement_correction_sets)) {
  do.call(combine_theta_smc_correction_sets, refinement_correction_sets)
} else {
  NULL
}

correction_summary <- summarize_correction_set(
  correction_set,
  correction_time[["elapsed"]] + refinement_elapsed
)

summarize_surrogate <- function(surrogate) {
  data.frame(
    training_points = nrow(surrogate$training),
    rejected_points = nrow(surrogate$rejected),
    selected_degree = surrogate$degree,
    selected_lambda = surrogate$lambda,
    selected_loo_nlp = surrogate$model_table$loo_nlp[1L],
    selected_loo_rmse = surrogate$model_table$loo_rmse[1L],
    selected_loo_mae = surrogate$model_table$loo_mae[1L],
    selected_loo_weighted_rmse = surrogate$model_table$loo_weighted_rmse[1L],
    selected_loo_max_abs_error = surrogate$model_table$loo_max_abs_error[1L],
    check.names = FALSE
  )
}

initial_surrogate_summary <- summarize_surrogate(initial_correction_surrogate)
surrogate_summary <- summarize_surrogate(correction_surrogate)

plot_support_ladder(stan_draws, rho_draws, plot_file)

elapsed_sec <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

saveRDS(
  list(
    stan_source = stan_results_file,
    sketch_set = sketch_set,
    fits = fits,
    stan_draws = stan_draws,
    rho_draws = rho_draws,
    support_draws = support_draws,
    support_table = support_table,
    rho_quantiles = rho_quantiles,
    support_theta = support_theta,
    support_w = support_w,
    pathfinder = pathfinder,
    pathfinder_summary = pathfinder_summary,
    correction_design = correction_design,
    initial_correction_set = initial_correction_set,
    initial_correction_summary = initial_correction_summary,
    initial_correction_surrogate = initial_correction_surrogate,
    initial_surrogate_summary = initial_surrogate_summary,
    refinement_design = refinement_design,
    refinement_correction_set = refinement_correction_set,
    refinement_summary = refinement_summary,
    correction_set = correction_set,
    correction_summary = correction_summary,
    correction_surrogate = correction_surrogate,
    surrogate_summary = surrogate_summary,
    timing = timing,
    settings = list(
      seed = base_seed,
      cores = cores,
      rho_ladder = rho_ladder,
      outer_particles = outer_particles,
      sketch_starts = sketch_starts,
      correction_particles = correction_particles,
      refinement_rounds = refinement_rounds,
      refinement_points_per_round = refinement_points_per_round,
      elapsed_sec = elapsed_sec
    ),
    plot_file = plot_file
  ),
  results_file
)

cat("Saved results to:", results_file, "\n")
cat("Saved plot to:", plot_file, "\n")
cat(sprintf("Elapsed: %.1f seconds\n", elapsed_sec))
cat("Timing:\n")
print(timing)
cat("Support recall:\n")
print(support_table)
cat("Pathfinder summary:\n")
print(pathfinder_summary)
cat("Pathfinder unique paths:\n")
print(pathfinder$path_summary[pathfinder$path_summary$unique, , drop = FALSE])
cat("Correction design:\n")
print(correction_design[, c("point_id", "source", "role", "label"), drop = FALSE])
cat("Initial correction summary:\n")
print(initial_correction_summary)
cat("Initial correction surrogate summary:\n")
print(initial_surrogate_summary)
cat("Refinement design:\n")
print(refinement_design[, c("point_id", "source", "role", "label", "parent_point_id", "target_point_id", "parent_loo_residual", "parent_relative_logposterior_smc"), drop = FALSE])
cat("Refinement summary:\n")
print(refinement_summary)
cat("Correction summary:\n")
print(correction_summary)
cat("Correction observations:\n")
print(correction_set$observations[, c("point_id", "replicate_id", "source", "role", "label", "delta_total", "mcse_total", "degenerate_subjects", "min_final_ess_frac"), drop = FALSE])
cat("Correction posterior relevance:\n")
print(correction_set$observations[, c("point_id", "replicate_id", "label", "relative_logposterior_smc", "relative_logposterior_tilde"), drop = FALSE])
cat("Correction replicate summary:\n")
print(correction_set$replicate_summary)
cat("Correction surrogate summary:\n")
print(surrogate_summary)
cat("Correction surrogate model comparison:\n")
print(correction_surrogate$model_table[, c("degree", "lambda", "loo_nlp", "loo_rmse", "loo_mae", "loo_weighted_rmse", "loo_max_abs_error", "training_rmse"), drop = FALSE])
cat("Correction surrogate validation:\n")
print(correction_surrogate$validation[, c("point_id", "label", "delta", "loo_pred", "loo_residual", "noise_sd", "relative_logposterior_smc"), drop = FALSE])
cat("Correction surrogate rejected points:\n")
print(correction_surrogate$rejected)
