rm(list = ls())

suppressPackageStartupMessages({
  library(parallel)
})

if (!file.exists("smc_core.R")) {
  stop("Run this script from the FastHierarchical repository root.")
}

parse_cli_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (!startsWith(arg, "--")) next
    arg <- sub("^--", "", arg)
    parts <- strsplit(arg, "=", fixed = TRUE)[[1L]]
    key <- gsub("-", "_", parts[1L])
    out[[key]] <- if (length(parts) > 1L) paste(parts[-1L], collapse = "=") else "true"
  }
  out
}

arg_chr <- function(args, key, default) {
  val <- args[[key]]
  if (is.null(val) || !nzchar(val)) default else as.character(val)
}

arg_int <- function(args, key, default) {
  val <- args[[key]]
  if (is.null(val) || !nzchar(val)) return(as.integer(default))
  as.integer(val)
}

arg_num <- function(args, key, default) {
  val <- args[[key]]
  if (is.null(val) || !nzchar(val)) return(as.numeric(default))
  as.numeric(val)
}

arg_lgl <- function(args, key, default = FALSE) {
  val <- args[[key]]
  if (is.null(val) || !nzchar(val)) return(isTRUE(default))
  tolower(as.character(val)) %in% c("1", "true", "t", "yes", "y")
}

arg_num_vec <- function(args, key, default) {
  val <- args[[key]]
  if (is.null(val) || !nzchar(val)) return(as.numeric(default))
  as.numeric(strsplit(val, ",", fixed = TRUE)[[1L]])
}

posterior_quantile_distance <- function(x, y, probs = seq(0.01, 0.99, length.out = 99L)) {
  qx <- stats::quantile(as.numeric(x), probs = probs, names = FALSE, type = 8, na.rm = TRUE)
  qy <- stats::quantile(as.numeric(y), probs = probs, names = FALSE, type = 8, na.rm = TRUE)
  mean(abs(qx - qy))
}

compare_posterior_draws <- function(stan_draws, workflow_draws) {
  common <- intersect(names(stan_draws), names(workflow_draws))
  rows <- lapply(common, function(nm) {
    st <- as.numeric(stan_draws[[nm]])
    wf <- as.numeric(workflow_draws[[nm]])
    q_st <- stats::quantile(st, probs = c(0.05, 0.5, 0.95), names = FALSE, type = 8, na.rm = TRUE)
    q_wf <- stats::quantile(wf, probs = c(0.05, 0.5, 0.95), names = FALSE, type = 8, na.rm = TRUE)
    st_sd <- stats::sd(st, na.rm = TRUE)
    data.frame(
      parameter = nm,
      stan_mean = mean(st, na.rm = TRUE),
      workflow_mean = mean(wf, na.rm = TRUE),
      mean_error = mean(wf, na.rm = TRUE) - mean(st, na.rm = TRUE),
      standardized_mean_error = (mean(wf, na.rm = TRUE) - mean(st, na.rm = TRUE)) / max(st_sd, .Machine$double.eps),
      stan_sd = st_sd,
      workflow_sd = stats::sd(wf, na.rm = TRUE),
      sd_ratio = stats::sd(wf, na.rm = TRUE) / max(st_sd, .Machine$double.eps),
      q05_error = q_wf[1L] - q_st[1L],
      q50_error = q_wf[2L] - q_st[2L],
      q95_error = q_wf[3L] - q_st[3L],
      q_wasserstein = posterior_quantile_distance(wf, st),
      stan_inside_workflow_q05_q95 = mean(st >= q_wf[1L] & st <= q_wf[3L], na.rm = TRUE),
      workflow_inside_stan_q05_q95 = mean(wf >= q_st[1L] & wf <= q_st[3L], na.rm = TRUE),
      check.names = FALSE
    )
  })
  do.call(rbind, rows)
}

detected_cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
if (!is.finite(detected_cores) || detected_cores < 1L) detected_cores <- 1L

cli_args <- parse_cli_args(commandArgs(trailingOnly = TRUE))

run_label <- arg_chr(cli_args, "label", "bank_smc_rho_anchor")
base_seed <- arg_int(cli_args, "base_seed", 20260519L)
cores <- arg_int(cli_args, "mc_cores", arg_int(cli_args, "cores", min(4L, detected_cores)))

local_particles <- arg_int(cli_args, "full_particles", arg_int(cli_args, "local_particles", 600L))
max_bank_nodes <- arg_int(cli_args, "max_bank_nodes", 6L)
target_local_ess <- arg_num(cli_args, "target_local_ess", 0.30)
local_mcmc_moves <- arg_int(cli_args, "local_mcmc_moves", 2L)
local_target_cess <- arg_num(cli_args, "local_target_cess", 0.90)
bridge_min_single_ess <- arg_num(cli_args, "bridge_min_single_ess", Inf)
bridge_max_pareto_k <- arg_num(cli_args, "bridge_max_pareto_k", 0.70)
design_max_points <- arg_int(cli_args, "design_max_points", max(2L, max_bank_nodes - 1L))
paired_effect_profiles <- arg_lgl(cli_args, "paired_effect_profiles", TRUE)
audit_max_points <- arg_int(cli_args, "audit_max_points", max(5L, max_bank_nodes - 1L))
audit_target_ess <- arg_num(cli_args, "audit_target_ess", target_local_ess)
audit_max_repairs <- arg_int(cli_args, "audit_max_repairs", NA_integer_)
audit_refine_rounds <- arg_int(cli_args, "audit_refine_rounds", 3L)
force_refine_points <- arg_int(cli_args, "force_refine_points", 0L)
force_profile_dims <- arg_int(cli_args, "force_profile_dims", 1L)
force_tail_probs <- arg_num_vec(cli_args, "force_tail_probs", c(0.01, 0.05, 0.95, 0.99))
rho_anchor_ladder <- arg_num_vec(cli_args, "rho_anchor_ladder", c(0.001, 0.01, 0.05, 0.15, 0.35, 0.75))
rho_anchor_enabled <- arg_lgl(cli_args, "rho_anchor_enabled", TRUE)
outer_particles <- arg_int(cli_args, "outer_particles", 1200L)
outer_mcmc_moves <- arg_int(cli_args, "outer_mcmc_moves", 3L)
outer_max_rounds <- arg_int(cli_args, "outer_max_rounds", 80L)
verbose <- arg_lgl(cli_args, "verbose", TRUE)

stan_results_file <- arg_chr(cli_args, "stan_results_file", file.path("benchmarks", "samples", "shifted_gamma_hierarchy_stan_results.rds"))
results_file <- arg_chr(cli_args, "results_file", file.path("benchmarks", "results", sprintf("shifted_gamma_%s_results.rds", run_label)))
plot_file <- arg_chr(cli_args, "plot_file", file.path("benchmarks", "results", sprintf("shifted_gamma_%s_posteriors.png", run_label)))
comparison_csv <- arg_chr(cli_args, "comparison_csv", file.path("benchmarks", "results", sprintf("shifted_gamma_%s_posterior_comparison.csv", run_label)))

dir.create(dirname(results_file), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(plot_file), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(comparison_csv), showWarnings = FALSE, recursive = TRUE)

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
if (!is.finite(audit_max_repairs)) {
  audit_max_repairs <- max(20L, nrow(y) * audit_max_points)
}

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
cat(sprintf("Bridge eligibility: single-node ESS %.3f | Pareto k %.2f\n",
            bridge_min_single_ess, bridge_max_pareto_k))
cat(sprintf("Design/audit: %d design points | %d audit points | audit ESS %.3f | refine rounds %d\n",
            design_max_points, audit_max_points, audit_target_ess, audit_refine_rounds))
cat(sprintf("Paired effect profiles: %s\n", if (isTRUE(paired_effect_profiles)) "enabled" else "disabled"))
cat(sprintf("Forced posterior-boundary anchors: %d points | %d dims\n",
            force_refine_points, force_profile_dims))
cat("Rho-sketch anchor ladder:", if (isTRUE(rho_anchor_enabled)) paste(rho_anchor_ladder, collapse = ", ") else "disabled", "\n")

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
    target_ess_frac = target_local_ess,
    target_cess = local_target_cess,
    n_mcmc_moves = local_mcmc_moves,
    bridge_min_single_ess_frac = bridge_min_single_ess,
    bridge_max_pareto_k = bridge_max_pareto_k
  ),
  design_control = list(
    max_points = design_max_points,
    paired_effect_profiles = paired_effect_profiles
  ),
  rho_anchor_control = list(
    enabled = rho_anchor_enabled,
    rho_ladder = rho_anchor_ladder,
    outer_particles = 500L,
    sketch_starts = 8L,
    support_weight_floor = 0.50,
    max_challengers = 1L
  ),
  outer_control = list(
    N = outer_particles,
    n_mcmc_moves = outer_mcmc_moves,
    max_rounds = outer_max_rounds
  ),
  audit_control = list(
    max_points = audit_max_points,
    target_ess_frac = audit_target_ess,
    max_repairs = audit_max_repairs,
    refine_rounds = audit_refine_rounds,
    force_refine_points = force_refine_points,
    force_profile_dims = force_profile_dims,
    force_tail_probs = force_tail_probs
  ),
  n_cores = cores,
  seed = base_seed,
  verbose = verbose
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

posterior_comparison <- compare_posterior_draws(stan_draws, workflow_draws)
utils::write.csv(posterior_comparison, comparison_csv, row.names = FALSE)

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
audit_entries <- Filter(
  function(x) is.list(x) && !is.null(x$audit) && "covered" %in% names(x$audit),
  bank_result$audits
)
audit_summary <- lapply(audit_entries, function(x) {
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
    posterior_comparison = posterior_comparison,
    bank_budget = budget_summary,
    audit_summary = audit_summary,
    settings = list(
      label = run_label,
      cores = cores,
      seed = base_seed,
      local_particles = local_particles,
      max_bank_nodes = max_bank_nodes,
      target_local_ess = target_local_ess,
      local_mcmc_moves = local_mcmc_moves,
      local_target_cess = local_target_cess,
      bridge_min_single_ess = bridge_min_single_ess,
      bridge_max_pareto_k = bridge_max_pareto_k,
      design_max_points = design_max_points,
      paired_effect_profiles = paired_effect_profiles,
      audit_max_points = audit_max_points,
      audit_target_ess = audit_target_ess,
      audit_max_repairs = audit_max_repairs,
      audit_refine_rounds = audit_refine_rounds,
      force_refine_points = force_refine_points,
      force_profile_dims = force_profile_dims,
      force_tail_probs = force_tail_probs,
      rho_anchor_enabled = rho_anchor_enabled,
      rho_anchor_ladder = rho_anchor_ladder,
      outer_particles = outer_particles,
      outer_mcmc_moves = outer_mcmc_moves,
      outer_max_rounds = outer_max_rounds,
      elapsed_sec = elapsed_sec
    ),
    plot_file = plot_file
  ),
  results_file
)

cat("Saved results to:", results_file, "\n")
cat("Saved plot to:", plot_file, "\n")
cat("Saved posterior comparison to:", comparison_csv, "\n")
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
cat("\nPosterior comparison to Stan:\n")
print(posterior_comparison[order(abs(posterior_comparison$standardized_mean_error), decreasing = TRUE), ], row.names = FALSE)
