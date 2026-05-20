rm(list = ls())

suppressPackageStartupMessages({
  library(parallel)
})

if (!file.exists("smc_core.R")) {
  stop("Run this script from the FastHierarchical repository root.")
}

detected_cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
if (!is.finite(detected_cores) || detected_cores < 1L) detected_cores <- 1L

run_label <- "bank_smc_full"
base_seed <- 20260519L
cores <- as.integer(min(4L, detected_cores))

local_particles <- 600L
max_bank_nodes <- 5L
target_local_ess <- 0.30

stan_results_file <- file.path("benchmarks", "samples", "shifted_gamma_hierarchy_stan_results.rds")
results_file <- file.path("benchmarks", "results", sprintf("shifted_gamma_%s_results.rds", run_label))
plot_file <- file.path("benchmarks", "results", sprintf("shifted_gamma_%s_posteriors.png", run_label))

dir.create(file.path("benchmarks", "results"), showWarnings = FALSE, recursive = TRUE)

source("smc_core.R")
source("reference_priors.R")
source("utilities.R")
source("SMC_super_fast.R")
source("population_models.R")
source("outer_population_smc.R")
source("bank_smc.R")

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

cat(sprintf("Loaded Stan benchmark bundle: %s\n", stan_results_file))
cat(sprintf("Data: %d subjects x %d trials\n", nrow(y), ncol(y)))
cat(sprintf("Run: label=%s | cores=%d | seed=%d\n", run_label, cores, base_seed))
cat(sprintf("Bank budget: %d nodes x %d particles per subject\n", max_bank_nodes, local_particles))

start_time <- Sys.time()

bank_result <- fit_bank_smc_population_model(
  data_list = data_list,
  loglik_fn = loglik_shifted_gamma,
  population_model = population_model,
  initial_theta = initial_theta,
  local_control = list(
    M = local_particles,
    bridge_particles = local_particles,
    max_nodes = max_bank_nodes,
    max_particles = max_bank_nodes * local_particles,
    target_ess_frac = target_local_ess
  ),
  n_cores = cores,
  seed = base_seed,
  verbose = TRUE
)

elapsed_sec <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
fit <- bank_result$fit

workflow_parts <- smc_posteriors(
  fit,
  n_draws = nrow(as.matrix(bundle$draws$mu)),
  seed = base_seed + 1L,
  population_model = population_model
)
workflow_mu <- as.data.frame(workflow_parts$mu, check.names = FALSE)
workflow_sigma2 <- as.data.frame(workflow_parts$sigma2, check.names = FALSE)
colnames(workflow_mu) <- alpha_names
colnames(workflow_sigma2) <- alpha_names

workflow_draws <- data.frame(
  mu_shape = workflow_mu[, "eta_shape"],
  mu_scale = workflow_mu[, "eta_scale"],
  mu_shift = workflow_mu[, "eta_shift"],
  sigma2_shape = workflow_sigma2[, "eta_shape"],
  sigma2_scale = workflow_sigma2[, "eta_scale"],
  sigma2_shift = workflow_sigma2[, "eta_shift"],
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

bank_state <- bank_result
bank_state$factor_set <- NULL
bank_state$population_model <- population_model

budget_summary <- bank_smc_bank_budget_summary(bank_result$banks)
audit_summary <- lapply(bank_result$audits, function(x) {
  data.frame(
    failures = sum(!x$audit$covered),
    min_ess_frac = min(x$audit$ess_frac, na.rm = TRUE),
    median_ess_frac = stats::median(x$audit$ess_frac, na.rm = TRUE),
    check.names = FALSE
  )
})
audit_summary <- if (length(audit_summary)) do.call(rbind, audit_summary) else NULL

saveRDS(
  list(
    stan_source = stan_results_file,
    bank_state = bank_state,
    stan_draws = stan_draws,
    workflow_draws = workflow_draws,
    posterior_summary = summarize_population_posterior_diag(
      theta = fit$theta,
      w = fit$w,
      model = population_model
    ),
    bank_budget = budget_summary,
    audit_summary = audit_summary,
    settings = list(
      label = run_label,
      cores = cores,
      seed = base_seed,
      local_particles = local_particles,
      max_bank_nodes = max_bank_nodes,
      target_local_ess = target_local_ess,
      elapsed_sec = elapsed_sec
    ),
    plot_file = plot_file
  ),
  results_file
)

cat("Saved results to:", results_file, "\n")
cat("Saved plot to:", plot_file, "\n")
cat(sprintf("Elapsed: %.1f seconds\n", elapsed_sec))
cat(sprintf(
  "Bank particles per subject: min=%d | median=%d | max=%d\n",
  min(budget_summary$particles),
  as.integer(stats::median(budget_summary$particles)),
  max(budget_summary$particles)
))
if (!is.null(audit_summary)) {
  cat(sprintf("Final audit failures: %d\n", tail(audit_summary$failures, 1L)))
}
