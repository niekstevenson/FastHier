rm(list = ls())
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[1L]))
} else {
  normalizePath("benchmarks/run_population_EMC.R")
}
repo_dir <- dirname(dirname(script_path))
setwd(repo_dir)

suppressPackageStartupMessages({
  library(EMC2)
})

set.seed(123L)

data_path <- file.path("benchmarks", "samples", "single_DDM30.RData")
reference_results_file <- file.path("benchmarks", "results", "reference_local_emc_results.rds")
results_file <- file.path("benchmarks", "results", "population_emc_results.rds")

reuse_reference_results <- TRUE

mc.cores <- 12L
outer_particles <- 2000L
outer_mcmc_moves <- 3L
outer_max_rounds <- 100L
base_seed <- 123L
verbose <- TRUE

dir.create(file.path("benchmarks", "samples"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path("benchmarks", "results"), showWarnings = FALSE, recursive = TRUE)

source("hierarchical_locals.R")
source("population_models.R")
source("outer_population_smc.R")

if (!file.exists(data_path)) {
  stop("Missing data file: ", data_path)
}

load(data_path)
if (!exists("single", inherits = FALSE)) {
  stop("The data file does not define an object named 'single'.")
}

emc <- single[[1]]
base_mu <- emc[[1]]$prior$theta_mu_mean
base_Sigma <- emc[[1]]$prior$theta_mu_var
param_names <- names(base_mu)
data_list <- lapply(single, function(x) x[[1]]$data[[1]])
model_factory <- emc[[1]]$model

loglik_emc2 <- function(Theta, data_i) {
  Theta <- as.matrix(Theta)
  colnames(Theta) <- param_names
  as.numeric(EMC2:::calc_ll_manager(Theta, data_i, model_factory, r_cores = 1L))
}

load_reference_stage <- function() {
  if (isTRUE(reuse_reference_results) && file.exists(reference_results_file)) {
    obj <- readRDS(reference_results_file)
    return(obj$stage %||% obj)
  }

  prepare_reference_local_stage(
    data_list = data_list,
    loglik_fn = loglik_emc2,
    base_mu = base_mu,
    base_Sigma = base_Sigma,
    broad_defensive = TRUE,
    n_jobs = mc.cores
  )
}

stage <- load_reference_stage()

population_model <- default_population_model_diag_gaussian(
  base_mu = base_mu,
  base_Sigma = base_Sigma,
  mean_var_scale = 1.0,
  sigma2_shape = 3.0,
  sigma2_mean = diag(base_Sigma)
)

factor_set <- build_population_factor_set(stage$local_objects, population_model)

fit <- outer_population_smc(
  factor_set = factor_set,
  N = outer_particles,
  n_mcmc_moves = outer_mcmc_moves,
  max_rounds = outer_max_rounds,
  n_cores = mc.cores,
  seed = base_seed,
  verbose = verbose
)

posterior_summary <- summarize_population_posterior_diag(
  theta = fit$theta,
  w = fit$w,
  model = population_model
)

theta_mean <- colSums(fit$theta * fit$w)
theta_mean <- matrix(theta_mean, nrow = 1L, dimnames = list(NULL, colnames(fit$theta)))
local_ess <- population_factor_set_local_ess(factor_set, theta = theta_mean)
local_ess_summary <- c(
  min = min(local_ess),
  q10 = as.numeric(stats::quantile(local_ess, probs = 0.10, names = FALSE)),
  median = as.numeric(stats::quantile(local_ess, probs = 0.50, names = FALSE)),
  mean = mean(local_ess)
)

saveRDS(
  list(
    fit = fit,
    posterior_summary = posterior_summary,
    local_ess = local_ess,
    local_ess_summary = local_ess_summary,
    population_model = population_model,
    settings = list(
      data_path = normalizePath(data_path, winslash = "/", mustWork = TRUE),
      reference_results_file = reference_results_file,
      results_file = results_file,
      reuse_reference_results = reuse_reference_results,
      mc.cores = mc.cores,
      outer_particles = outer_particles,
      outer_mcmc_moves = outer_mcmc_moves,
      outer_max_rounds = outer_max_rounds,
      base_seed = base_seed
    )
  ),
  results_file
)

cat("Saved results to:", results_file, "\n")
cat(sprintf("Locals: %d\n", length(stage$local_objects)))
cat(sprintf("Outer rounds: %d\n", fit$meta$rounds))
cat(sprintf("Final beta: %.4f\n", fit$beta))
cat(sprintf("Log evidence: %.4f +/- %.4f\n", fit$log_evidence, fit$mcse_log_evidence))
cat(sprintf("Local ESS at posterior mean: min=%.1f q10=%.1f median=%.1f mean=%.1f\n",
            local_ess_summary["min"],
            local_ess_summary["q10"],
            local_ess_summary["median"],
            local_ess_summary["mean"]))

cat("\nPopulation mean summary:\n")
print(posterior_summary$mu, row.names = FALSE)

cat("\nPopulation variance summary:\n")
print(posterior_summary$sigma2, row.names = FALSE)
