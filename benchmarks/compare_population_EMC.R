rm(list = ls())

suppressPackageStartupMessages({
  library(parallel)
  library(EMC2)
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

data_file <- file.path("benchmarks", "samples", "full_EMC2.RData")
results_file <- file.path("benchmarks", "results", sprintf("population_emc_%s_results.rds", run_label))
plot_file <- file.path("benchmarks", "results", sprintf("population_emc_%s_posteriors.png", run_label))

local_control <- list(
  planned_anchor_count = 4L
)
outer_control <- list()
certification_control <- list()

dir.create(file.path("benchmarks", "results"), showWarnings = FALSE, recursive = TRUE)

source("smc_core.R")
source("reference_priors.R")
source("utilities.R")
source("SMC_super_fast.R")
source("hierarchical_locals.R")
source("population_models.R")
source("outer_population_smc.R")

if (!file.exists(data_file)) {
  stop("Missing EMC2 benchmark data: ", data_file)
}

load(data_file)
if (!exists("ELP_DDM", inherits = FALSE)) {
  stop("The EMC2 data file must define ELP_DDM.")
}

emc <- ELP_DDM[[1L]]
data_list <- emc$data
model_factory <- emc$model
alpha_names <- emc$par_names
alpha_dim <- length(alpha_names)

base_mu <- stats::setNames(rep(0, alpha_dim), alpha_names)
base_Sigma <- diag(1, nrow = alpha_dim, ncol = alpha_dim)
dimnames(base_Sigma) <- list(alpha_names, alpha_names)

loglik_emc2 <- function(Theta, data_i) {
  Theta <- as.matrix(Theta)
  colnames(Theta) <- alpha_names
  out <- as.numeric(EMC2:::calc_ll_manager(Theta, data_i, model_factory, r_cores = 1L))
  bad <- !is.finite(out) | out > 100
  if (any(bad)) {
    finite_good <- out[!bad & is.finite(out)]
    out[bad] <- if (length(finite_good)) min(finite_good) else -1e12
  }
  out
}

population_model <- make_population_model_diag_gaussian(
  alpha_names = alpha_names,
  mean_prior_mean = rep(0, alpha_dim),
  mean_prior_var = rep(1, alpha_dim),
  sigma2_prior_shape = rep(2, alpha_dim),
  sigma2_prior_rate = rep(0.3, alpha_dim),
  label = "emc_normal_gamma"
)

emc_mu <- as.data.frame(parameters(ELP_DDM, selection = "mu"), check.names = FALSE)
emc_sigma2 <- as.data.frame(parameters(ELP_DDM, selection = "sigma2"), check.names = FALSE)
colnames(emc_mu) <- paste0("mu_", alpha_names)
colnames(emc_sigma2) <- paste0("sigma2_", alpha_names)
emc_draws <- data.frame(emc_mu, emc_sigma2, check.names = FALSE)

cat(sprintf("Loaded EMC2 benchmark data: %s\n", data_file))
cat(sprintf("Data: %d subjects | parameters: %d\n", length(data_list), alpha_dim))
cat(sprintf("Run: label=%s | cores=%d | seed=%d\n", run_label, cores, base_seed))
cat("Running local reference stage...\n")

start_time <- proc.time()[["elapsed"]]

stage <- do.call(
  prepare_reference_local_stage,
  c(
    list(
      data_list = data_list,
      loglik_fn = loglik_emc2,
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
    loglik_fn = loglik_emc2,
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

workflow_parts <- smc_posteriors(
  fit,
  n_draws = nrow(emc_mu),
  seed = base_seed + 1L,
  population_model = population_model
)
workflow_mu <- as.data.frame(workflow_parts$mu, check.names = FALSE)
workflow_sigma2 <- as.data.frame(workflow_parts$sigma2, check.names = FALSE)
colnames(workflow_mu) <- paste0("mu_", alpha_names)
colnames(workflow_sigma2) <- paste0("sigma2_", alpha_names)
workflow_draws <- data.frame(workflow_mu, workflow_sigma2, check.names = FALSE)

grDevices::png(plot_file, width = 1800, height = 1400)
plot_posteriors(
  emc_draws,
  workflow_draws,
  labels = c("EMC2", run_label),
  cols = c("black", "firebrick3"),
  n_cols = 4L
)
grDevices::dev.off()

elapsed_sec <- proc.time()[["elapsed"]] - start_time
certified_compact <- compact_certified_population_result(certified_result)

saveRDS(
  list(
    data_source = data_file,
    certified = certified_compact,
    stage = list(
      pilot_indices = stage$pilot$selection$indices,
      local_budget = .reference_local_budget_summary(stage$local_objects)
    ),
    emc_draws = emc_draws,
    workflow_draws = workflow_draws,
    posterior_summary = summarize_population_posterior_diag(
      theta = fit$theta,
      w = fit$w,
      model = population_model
    ),
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
