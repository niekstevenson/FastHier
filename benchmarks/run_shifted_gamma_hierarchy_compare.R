rm(list = ls())
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[1L]))
} else {
  normalizePath("benchmarks/run_shifted_gamma_hierarchy_compare.R")
}
repo_dir <- dirname(dirname(script_path))
setwd(repo_dir)

suppressPackageStartupMessages({
  library(parallel)
})

set.seed(20260324L)

stan_results_file <- file.path("benchmarks", "samples", "shifted_gamma_hierarchy_stan_results.rds")
results_file <- file.path("benchmarks", "results", "shifted_gamma_hierarchy_current_results.rds")
plot_file <- file.path("benchmarks", "results", "shifted_gamma_hierarchy_current_posteriors.png")

detected_cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
if (!is.finite(detected_cores) || detected_cores < 1L) {
  detected_cores <- 1L
}

mc.cores <- as.integer(max(1L, min(4L, detected_cores)))
pilot_size <- 10L
pilot_particles <- 800L
full_particles <- 2000L
outer_particles <- 2000L
outer_mcmc_moves <- 3L
outer_max_rounds <- 80L
base_seed <- 20260324L
verbose <- TRUE

dir.create(file.path("benchmarks", "samples"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path("benchmarks", "results"), showWarnings = FALSE, recursive = TRUE)

source("hierarchical_locals.R")
source("population_models.R")
source("outer_population_smc.R")
source("utilities.R")

if (!file.exists(stan_results_file)) {
  stop("Missing Stan benchmark results: ", stan_results_file)
}

bundle <- readRDS(stan_results_file)

y <- bundle$data$y
m0 <- as.numeric(bundle$priors$m0)
s0 <- as.numeric(bundle$priors$s0)
a0 <- as.numeric(bundle$priors$a0)
b0 <- as.numeric(bundle$priors$b0)

alpha_names <- c("eta_shape", "eta_scale", "eta_shift")

data_list <- lapply(seq_len(nrow(y)), function(i) y[i, ])

names(m0) <- alpha_names
names(s0) <- alpha_names
names(a0) <- alpha_names
names(b0) <- alpha_names

base_mu <- m0
base_var <- s0 + b0 / (a0 - 1)
base_Sigma <- diag(base_var, nrow = length(alpha_names))
dimnames(base_Sigma) <- list(alpha_names, alpha_names)

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
  if (!any(ok)) {
    return(out)
  }

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

stan_draws <- data.frame(
  mu_shape = bundle$draws$mu[, 1L],
  mu_scale = bundle$draws$mu[, 2L],
  mu_shift = bundle$draws$mu[, 3L],
  sigma2_shape = bundle$draws$sigma2[, 1L],
  sigma2_scale = bundle$draws$sigma2[, 2L],
  sigma2_shift = bundle$draws$sigma2[, 3L],
  check.names = FALSE
)

cat(sprintf("Loaded Stan benchmark bundle: %s\n", stan_results_file))
cat(sprintf("Data: %d subjects x %d trials\n", nrow(y), ncol(y)))
cat("Running local-reference stage...\n")

stage <- prepare_reference_local_stage(
  data_list = data_list,
  loglik_fn = loglik_shifted_gamma,
  base_mu = base_mu,
  base_Sigma = base_Sigma,
  pilot_size = min(pilot_size, length(data_list)),
  broad_scale = 1,
  pilot_particles = pilot_particles,
  full_particles = full_particles,
  n_jobs = mc.cores,
  base_seed = base_seed,
  pilot_smc_control = list(
    max_rounds = 40L
  ),
  full_smc_control = list(
    hist_mix_enable = FALSE,
    gss_enable = FALSE,
    da_enable = FALSE
  )
)

cat("Running outer population SMC...\n")

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

workflow_theta <- smc_posteriors(
  fit,
  n_draws = nrow(as.matrix(bundle$draws$mu)),
  seed = base_seed + 1L
)

workflow_draws <- data.frame(
  mu_shape = workflow_theta$mu_eta_shape,
  mu_scale = workflow_theta$mu_eta_scale,
  mu_shift = workflow_theta$mu_eta_shift,
  sigma2_shape = exp(workflow_theta$log_sigma2_eta_shape),
  sigma2_scale = exp(workflow_theta$log_sigma2_eta_scale),
  sigma2_shift = exp(workflow_theta$log_sigma2_eta_shift),
  check.names = FALSE
)

grDevices::png(plot_file, width = 1400, height = 900)
plot_posteriors(
  stan_draws,
  workflow_draws,
  labels = c("Stan", "Current workflow"),
  cols = c("black", "firebrick3"),
  n_cols = 3L
)
grDevices::dev.off()

saveRDS(
  list(
    stan_source = stan_results_file,
    stage = stage,
    fit = fit,
    stan_draws = stan_draws,
    workflow_draws = workflow_draws,
    settings = list(
      mc.cores = mc.cores,
      pilot_size = min(pilot_size, length(data_list)),
      pilot_particles = pilot_particles,
      full_particles = full_particles,
      outer_particles = outer_particles,
      outer_mcmc_moves = outer_mcmc_moves,
      outer_max_rounds = outer_max_rounds,
      base_seed = base_seed
    ),
    plot_file = plot_file
  ),
  results_file
)

cat("Saved results to:", results_file, "\n")
cat("Saved plot to:", plot_file, "\n")
