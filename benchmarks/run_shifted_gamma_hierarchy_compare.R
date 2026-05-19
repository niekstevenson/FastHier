rm(list = ls())

suppressPackageStartupMessages({
  library(parallel)
})

if (!file.exists("smc_core.R")) {
  stop("Run this script from the FastHierarchical repository root.")
}

detected_cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
if (!is.finite(detected_cores) || detected_cores < 1L) detected_cores <- 1L

run_label <- "planned_anchor_audit"
base_seed <- 20260324L
cores <- as.integer(min(4L, detected_cores))
verbose <- FALSE

stan_results_file <- file.path("benchmarks", "samples", "shifted_gamma_hierarchy_stan_results.rds")
results_file <- file.path("benchmarks", "results", sprintf("shifted_gamma_%s_results.rds", run_label))
plot_file <- file.path("benchmarks", "results", sprintf("shifted_gamma_%s_posteriors.png", run_label))

local_control <- list(
  planned_anchor_count = 4L
)
outer_control <- list()
certification_control <- list()

dir.create(file.path("benchmarks", "samples"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path("benchmarks", "results"), showWarnings = FALSE, recursive = TRUE)

source("smc_core.R")
source("reference_priors.R")
source("utilities.R")
source("SMC_super_fast.R")
source("hierarchical_locals.R")
source("population_models.R")
source("outer_population_smc.R")

if (!file.exists(stan_results_file)) {
  stop("Missing Stan benchmark results: ", stan_results_file)
}

set.seed(base_seed)
bundle <- readRDS(stan_results_file)

y <- bundle$data$y
data_list <- lapply(seq_len(nrow(y)), function(i) y[i, ])

alpha_names <- c("eta_shape", "eta_scale", "eta_shift")
m0 <- stats::setNames(as.numeric(bundle$priors$m0), alpha_names)
s0 <- stats::setNames(as.numeric(bundle$priors$s0), alpha_names)
a0 <- stats::setNames(as.numeric(bundle$priors$a0), alpha_names)
b0 <- stats::setNames(as.numeric(bundle$priors$b0), alpha_names)

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
cat(sprintf("Run: label=%s | cores=%d | seed=%d\n", run_label, cores, base_seed))
cat("Running local reference stage...\n")

start_time <- proc.time()[["elapsed"]]

stage <- do.call(
  prepare_reference_local_stage,
  c(
    list(
      data_list = data_list,
      loglik_fn = loglik_shifted_gamma,
      base_mu = base_mu,
      base_Sigma = base_Sigma,
      population_model = population_model,
      n_jobs = cores,
      base_seed = base_seed,
      verbose = verbose
    ),
    local_control
  )
)

cat("Running population stage...\n")

certified_result <- do.call(
  fit_certified_population_model,
  list(
    data_list = data_list,
    loglik_fn = loglik_shifted_gamma,
    local_objects = stage$local_objects,
    population_model = population_model,
    outer_control = outer_control,
    certification_control = certification_control,
    n_cores = cores,
    seed = base_seed + 300000L,
    verbose = verbose
  )
)
fit <- certified_result$fit
elapsed_sec <- proc.time()[["elapsed"]] - start_time

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
  labels = c("Stan", run_label),
  cols = c("black", "firebrick3"),
  n_cols = 3L
)
grDevices::dev.off()

certified_compact <- compact_certified_population_result(certified_result)

saveRDS(
  list(
    stan_source = stan_results_file,
    certified = certified_compact,
    stage = list(
      pilot_indices = stage$pilot$selection$indices,
      local_budget = .reference_local_budget_summary(stage$local_objects)
    ),
    stan_draws = stan_draws,
    workflow_draws = workflow_draws,
    settings = list(
      label = run_label,
      cores = cores,
      seed = base_seed,
      local_control = local_control,
      outer_control = outer_control,
      certification_control = certification_control,
      elapsed_sec = elapsed_sec
    ),
    plot_file = plot_file
  ),
  results_file
)

cat("Saved results to:", results_file, "\n")
cat("Saved plot to:", plot_file, "\n")
cat(sprintf("Elapsed: %.1f seconds\n", elapsed_sec))
