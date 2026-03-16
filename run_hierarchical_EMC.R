#!/usr/bin/env Rscript

rm(list = ls())

suppressPackageStartupMessages({
  library(EMC2)
  library(parallel)
})

set.seed(123L)

data_path <- "single_DDM30.RData"
local_results_file <- "local_emc_results.rds"
hierarchical_results_file <- "hierarchical_emc_results.rds"

mc.cores <- 12L
local_particles <- 4000L
outer_particles <- 2000L
local_bank_particles <- 64L
outer_rounds <- 75L
hierarchical_method <- "shared_bank"
collapsed_local_particles <- 1000L
collapsed_rejuvenation <- "none"
collapsed_rejuvenate_every <- Inf
collapsed_rejuvenate_after_resample <- FALSE
collapsed_n_rejuvenation_moves <- 1L
collapsed_bridge_steps <- 4L
collapsed_bridge_schedule <- NULL
base_seed <- 123L
outer_seed <- 2026L
verbose <- TRUE
reuse_local_results <- TRUE

source("SMC_super_fast.R")
source("hierarchical_smc.R")

if (!file.exists(data_path)) {
  stop("Missing data file: ", data_path)
}

load(data_path)
if (!exists("single", inherits = FALSE)) {
  stop("The data file does not define an object named 'single'.")
}

emc <- single[[1]]
mu_ref <- emc[[1]]$prior$theta_mu_mean
Sigma_ref <- emc[[1]]$prior$theta_mu_var
param_names <- names(mu_ref)
data_list <- lapply(single, function(x) x[[1]]$data[[1]])
model_factory <- emc[[1]]$model

loglik_emc2 <- function(Theta, data_i) {
  Theta <- as.matrix(Theta)
  colnames(Theta) <- param_names
  as.numeric(EMC2:::calc_ll_manager(Theta, data_i, model_factory, r_cores = 1L))
}

load_local_fits <- function() {
  if (!reuse_local_results || !file.exists(local_results_file)) return(NULL)
  obj <- readRDS(local_results_file)
  fits <- obj$smc %||% obj
  if (!is.list(fits) || length(fits) != length(data_list)) {
    stop("Existing local results do not match the current data.")
  }
  fits
}

population_prior <- default_population_prior_diag(
  mu_ref = mu_ref,
  Sigma_ref = Sigma_ref,
  mean_var_scale = 1.0,
  sigma2_shape = 3.0,
  sigma2_mean = diag(Sigma_ref)
)

local_fits <- NULL
local_objs <- NULL
gaussian_map_fn <- phi_to_gaussian_params_diag_factory(param_names = param_names)

if (identical(hierarchical_method, "shared_bank")) {
  local_fits <- load_local_fits()
  if (is.null(local_fits)) {
    local_fits <- run_local_smc_subjects(
      data_list = data_list,
      loglik_fn = loglik_emc2,
      mu_ref = mu_ref,
      Sigma_ref = Sigma_ref,
      M = local_particles,
      n_cores = mc.cores,
      base_seed = base_seed,
      verbose = verbose
    )
  }

  working_prior <- make_working_prior_gaussian(mu_ref, Sigma_ref)
  local_objs <- build_local_exact_objects(
    local_fits = local_fits,
    data_list = data_list,
    loglik_fn = loglik_emc2,
    working_prior = working_prior
  )

  hier_fit <- hierarchical_smc(
    method = "shared_bank",
    local_objs = local_objs,
    prior = population_prior,
    gaussian_map_fn = gaussian_map_fn,
    N = outer_particles,
    M_local = local_bank_particles,
    rho_res = 0.5,
    rho_local = 0.5,
    n_population_moves = 1L,
    max_rounds = outer_rounds,
    max_bank_topups = 3L,
    bank_split_tol = 0.01,
    seed = outer_seed,
    verbose = verbose
  )
} else if (identical(hierarchical_method, "collapsed_subject")) {
  hier_fit <- hierarchical_smc(
    method = "collapsed_subject",
    data_list = data_list,
    loglik_fn = loglik_emc2,
    prior = population_prior,
    gaussian_map_fn = gaussian_map_fn,
    d_theta = length(param_names),
    N = outer_particles,
    local_particles = collapsed_local_particles,
    resample_threshold = 0.5,
    rejuvenation = collapsed_rejuvenation,
    rejuvenate_every = collapsed_rejuvenate_every,
    rejuvenate_after_resample = collapsed_rejuvenate_after_resample,
    n_rejuvenation_moves = collapsed_n_rejuvenation_moves,
    theta_names = param_names,
    local_n_cores = 1L,
    outer_n_cores = mc.cores,
    base_seed = outer_seed,
    verbose = verbose,
    local_smc_control = list(
      max_rounds = outer_rounds
    )
  )
} else if (identical(hierarchical_method, "collapsed_bridge")) {
  hier_fit <- hierarchical_smc(
    method = "collapsed_bridge",
    data_list = data_list,
    loglik_fn = loglik_emc2,
    prior = population_prior,
    gaussian_map_fn = gaussian_map_fn,
    d_theta = length(param_names),
    N = outer_particles,
    local_particles = collapsed_local_particles,
    resample_threshold = 0.5,
    bridge_schedule = collapsed_bridge_schedule,
    n_bridge_steps = collapsed_bridge_steps,
    rejuvenation = collapsed_rejuvenation,
    rejuvenate_every = collapsed_rejuvenate_every,
    rejuvenate_after_resample = collapsed_rejuvenate_after_resample,
    n_rejuvenation_moves = collapsed_n_rejuvenation_moves,
    theta_names = param_names,
    local_n_cores = 1L,
    outer_n_cores = mc.cores,
    base_seed = outer_seed,
    verbose = verbose,
    local_smc_control = list(
      max_rounds = outer_rounds
    )
  )
} else {
  stop("Unknown hierarchical method: ", hierarchical_method)
}

posterior_summary <- summarize_phi_diag_posterior(
  phi = hier_fit$phi,
  w = hier_fit$w,
  param_names = param_names
)

saveRDS(
  list(
    local_fits = local_fits,
    local_objects = local_objs,
    hierarchical_fit = hier_fit,
    posterior_summary = posterior_summary,
    settings = list(
      data_path = normalizePath(data_path, winslash = "/", mustWork = TRUE),
      local_results_file = local_results_file,
      hierarchical_results_file = hierarchical_results_file,
      mc.cores = mc.cores,
      hierarchical_method = hierarchical_method,
      local_particles = local_particles,
      outer_particles = outer_particles,
      local_bank_particles = local_bank_particles,
      collapsed_local_particles = collapsed_local_particles,
      collapsed_rejuvenation = collapsed_rejuvenation,
      collapsed_bridge_steps = collapsed_bridge_steps,
      collapsed_bridge_schedule = collapsed_bridge_schedule,
      outer_rounds = outer_rounds,
      base_seed = base_seed,
      outer_seed = outer_seed,
      reuse_local_results = reuse_local_results
    )
  ),
  hierarchical_results_file
)

cat("Saved hierarchical results to:", hierarchical_results_file, "\n")
cat(sprintf("Subjects: %d\n", length(data_list)))
cat(sprintf("Approx log evidence: %.4f\n", hier_fit$log_evidence))
cat(sprintf("Method: %s\n", hierarchical_method))

if (identical(hierarchical_method, "shared_bank")) {
  cat(sprintf("Outer rounds: %d\n", hier_fit$meta$rounds))
  cat(sprintf("Final lambda: %.4f\n", tail(hier_fit$meta$lambda_hist, 1L)))
  cat(sprintf("Mean bank size (final): %.1f\n", tail(hier_fit$meta$bank_size_mean_hist, 1L)))
} else {
  cat(sprintf("Subject stages: %d\n", hier_fit$meta$n_subjects))
  if (identical(hierarchical_method, "collapsed_bridge")) {
    cat(sprintf("Bridge stages: %d\n", hier_fit$meta$n_stages))
  }
  cat(sprintf("Local runs: %d\n", hier_fit$meta$n_local_runs))
}

cat("\nPopulation mean summary:\n")
print(posterior_summary$mu, row.names = FALSE)

cat("\nPopulation variance summary:\n")
print(posterior_summary$sigma2, row.names = FALSE)
