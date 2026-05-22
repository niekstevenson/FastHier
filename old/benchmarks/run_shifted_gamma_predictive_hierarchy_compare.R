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

arg_num_vec <- function(args, key, default) {
  val <- args[[key]]
  if (is.null(val) || !nzchar(val)) return(as.numeric(default))
  as.numeric(strsplit(val, ",", fixed = TRUE)[[1L]])
}

arg_lgl <- function(args, key, default = FALSE) {
  val <- args[[key]]
  if (is.null(val) || !nzchar(val)) return(isTRUE(default))
  tolower(as.character(val)) %in% c("1", "true", "t", "yes", "y")
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

run_label <- arg_chr(cli_args, "label", "q0_predictive")
base_seed <- arg_int(cli_args, "base_seed", 20260522L)
cores <- arg_int(cli_args, "mc_cores", arg_int(cli_args, "cores", min(4L, detected_cores)))

local_particles <- arg_int(cli_args, "local_particles", 600L)
reference_strata <- arg_int(cli_args, "reference_strata", 1L)
predictive_components <- arg_int(cli_args, "predictive_components", 96L)
outer_particles <- arg_int(cli_args, "outer_particles", 1600L)
outer_mcmc_moves <- arg_int(cli_args, "outer_mcmc_moves", 4L)
outer_max_rounds <- arg_int(cli_args, "outer_max_rounds", 80L)
audit_q0 <- arg_int(cli_args, "audit_q0", 80L)
audit_outer <- arg_int(cli_args, "audit_outer", 80L)
audit_tail <- arg_int(cli_args, "audit_tail", 80L)
audit_q0_axis <- arg_int(cli_args, "audit_q0_axis", 40L)
audit_outer_axis <- arg_int(cli_args, "audit_outer_axis", 40L)
outer_axis_inflate <- arg_num(cli_args, "outer_axis_inflate", 1.5)
audit_ess <- arg_num(cli_args, "audit_ess", 0.05)
compute_psis <- arg_lgl(cli_args, "compute_psis", FALSE)
compute_loo <- arg_lgl(cli_args, "compute_loo", FALSE)
max_pareto_k <- arg_num(cli_args, "max_pareto_k", 0.7)
max_loo_delta <- arg_num(cli_args, "max_loo_delta", Inf)
adaptive_refinement <- arg_lgl(cli_args, "adaptive_refinement", TRUE)
pre_outer_support_refinement <- arg_lgl(cli_args, "pre_outer_support_refinement", FALSE)
pre_outer_support_rounds <- arg_int(cli_args, "pre_outer_support_rounds", 1L)
max_adapt_rounds <- arg_int(cli_args, "max_adapt_rounds", 3L)
max_new_states <- arg_int(cli_args, "max_new_states", 45L)
max_pre_outer_states <- arg_int(cli_args, "max_pre_outer_states", max_new_states)
max_new_states_per_local <- arg_int(cli_args, "max_new_states_per_local", 3L)
theta_anchor_scales <- arg_num_vec(cli_args, "theta_anchor_scales", c(1, 2, 4, 8))
min_refinement_impact_weight <- arg_num(cli_args, "min_refinement_impact_weight", 1e-4)
min_surface_score <- arg_num(cli_args, "min_surface_score", 1e-3)
state_reject_uncertainty_ratio <- arg_num(cli_args, "state_reject_uncertainty_ratio", Inf)
state_reject_ess_ratio <- arg_num(cli_args, "state_reject_ess_ratio", 0)
protected_refinement_fraction <- arg_num(cli_args, "protected_refinement_fraction", 0.35)
min_protected_impact_weight <- if (is.null(cli_args[["min_protected_impact_weight"]])) {
  c(explicit = 0, q0_tail = 1e-8, outer_axis = 1e-6, q0_axis = 1e-6)
} else {
  arg_num(cli_args, "min_protected_impact_weight", 1e-6)
}
protected_source_weights <- c(explicit = 1, q0_tail = 3, outer_axis = 2, q0_axis = 1)
repair_enabled <- arg_lgl(cli_args, "repair_enabled", FALSE)
max_repairs <- arg_int(cli_args, "max_repairs", 0L)
adaptive_outer_proposal <- arg_lgl(cli_args, "adaptive_outer_proposal", FALSE)
outer_proposal_defensive_weight <- arg_num(cli_args, "outer_proposal_defensive_weight", 0.35)
local_sketch_enabled <- arg_lgl(cli_args, "local_sketch_enabled", TRUE)
sketch_starts <- arg_int(cli_args, "sketch_starts", 12L)
sketch_scale <- arg_num(cli_args, "sketch_scale", 6)
mbar_enabled <- arg_lgl(cli_args, "mbar_enabled", TRUE)
mbar_min_state_ess_frac <- arg_num(cli_args, "mbar_min_state_ess_frac", 0.01)
mbar_max_abs_shift <- arg_num(cli_args, "mbar_max_abs_shift", 25)
verbose <- arg_lgl(cli_args, "verbose", TRUE)

stan_results_file <- arg_chr(cli_args, "stan_results_file", file.path("benchmarks", "samples", "shifted_gamma_hierarchy_stan_results.rds"))
q0_results_file <- arg_chr(cli_args, "q0_results_file", file.path("benchmarks", "results", "shifted_gamma_support_diagnostics.rds"))
q0_method <- arg_chr(cli_args, "q0_method", "laplace_t,dmis_broad_laplace")
q0_method_weights <- arg_num_vec(cli_args, "q0_method_weights", c(0.95, 0.05))
results_file <- arg_chr(cli_args, "results_file", file.path("benchmarks", "results", sprintf("shifted_gamma_%s_results.rds", run_label)))
plot_file <- arg_chr(cli_args, "plot_file", file.path("benchmarks", "results", sprintf("shifted_gamma_%s_posteriors.png", run_label)))
comparison_csv <- arg_chr(cli_args, "comparison_csv", file.path("benchmarks", "results", sprintf("shifted_gamma_%s_posterior_comparison.csv", run_label)))
audit_csv <- arg_chr(cli_args, "audit_csv", file.path("benchmarks", "results", sprintf("shifted_gamma_%s_local_audit.csv", run_label)))
refinements_csv <- arg_chr(cli_args, "refinements_csv", file.path("benchmarks", "results", sprintf("shifted_gamma_%s_refinements.csv", run_label)))

dir.create(dirname(results_file), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(plot_file), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(comparison_csv), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(audit_csv), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(refinements_csv), showWarnings = FALSE, recursive = TRUE)

source("smc_core.R")
source("reference_priors.R")
source("utilities.R")
source("SMC_super_fast.R")
source("population_models.R")
source("theta_proposals.R")
source("local_likelihood_sketches.R")
source("local_predictive_smc.R")

if (!file.exists(stan_results_file)) {
  stop("Missing Stan benchmark results: ", stan_results_file)
}
if (!file.exists(q0_results_file)) {
  stop("Missing q0 support results: ", q0_results_file)
}

set.seed(base_seed)
bundle <- readRDS(stan_results_file)
q0_bundle <- readRDS(q0_results_file)
q0_methods <- trimws(strsplit(q0_method, ",", fixed = TRUE)[[1L]])
q0_methods <- q0_methods[nzchar(q0_methods)]
if (!length(q0_methods)) stop("q0_method must name at least one fitted q0 proposal.")
q0_parts <- lapply(q0_methods, function(method) {
  proposal <- q0_bundle$fitted_q0_proposals[[method]]
  if (is.null(proposal)) {
    stop("q0 method not found in support results: ", method)
  }
  proposal
})
names(q0_parts) <- q0_methods
if (length(q0_parts) == 1L) {
  theta_proposal <- q0_parts[[1L]]
} else {
  if (length(q0_method_weights) != length(q0_parts)) {
    stop("q0_method_weights must match the number of q0 methods.")
  }
  theta_proposal <- combine_theta_q0_proposals(
    proposals = q0_parts,
    weights = q0_method_weights,
    label = paste0("combined_", paste(q0_methods, collapse = "_"))
  )
}

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
theta_proposal <- normalize_theta_proposal(theta_proposal, population_model = population_model)

initial_theta <- matrix(c(m0, log(b0 / (a0 - 1))), nrow = 1L)
colnames(initial_theta) <- population_model$hyper_names

local_sketch_reference_fn <- if (isTRUE(local_sketch_enabled)) {
  function(local_id, data_i, population_model, theta_proposal) {
    sketch <- fit_local_likelihood_sketch(
      data_i = data_i,
      loglik_fn = loglik_shifted_gamma,
      alpha_names = alpha_names,
      local_id = local_id,
      population_model = population_model,
      theta_reference = initial_theta,
      n_starts = sketch_starts,
      start_scale = 9,
      reference_scale = 50,
      seed = base_seed + 300000L + local_id
    )
    mu <- as.numeric(sketch$component_means[1L, ])
    names(mu) <- alpha_names
    Sigma <- regularize_cov(
      as.numeric(sketch_scale) * sketch$component_covs[[1L]],
      min_eig = 1e-8,
      cond_cap = 1e8
    )
    dimnames(Sigma) <- list(alpha_names, alpha_names)
    list(make_reference_prior_gaussian(
      mu = mu,
      Sigma = Sigma,
      param_names = alpha_names,
      label = sprintf("local_likelihood_sketch_reference_%d", local_id)
    ))
  }
} else {
  NULL
}

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
cat(sprintf("Loaded q0 support bundle: %s | method=%s\n", q0_results_file, q0_method))
if (length(q0_parts) > 1L) {
  cat(sprintf("q0 mixture weights: %s\n", paste(q0_method_weights, collapse = ", ")))
}
cat(sprintf("Data: %d subjects x %d trials\n", nrow(y), ncol(y)))
cat(sprintf("Run: label=%s | cores=%d | seed=%d\n", run_label, cores, base_seed))
cat(sprintf("q0-predictive locals: %d particles | %d alpha reference components\n",
            local_particles, predictive_components))
cat(sprintf("Reference strata per local: %d\n", reference_strata))
cat(sprintf("Local sketch reference: %s\n", if (isTRUE(local_sketch_enabled)) "enabled" else "disabled"))
cat(sprintf("MBAR calibration: %s\n", if (isTRUE(mbar_enabled)) "enabled" else "disabled"))
cat(sprintf("Adaptive refinement: %s | rounds=%d | max new states=%d\n",
            if (isTRUE(adaptive_refinement)) "enabled" else "disabled",
            max_adapt_rounds,
            max_new_states))
cat(sprintf("Pre-outer support refinement: %s | rounds=%d | max states=%d\n",
            if (isTRUE(pre_outer_support_refinement)) "enabled" else "disabled",
            pre_outer_support_rounds,
            max_pre_outer_states))
cat(sprintf("Refinement impact floor: %.3g | theta anchor scales=%s\n",
            min_refinement_impact_weight,
            paste(theta_anchor_scales, collapse = ", ")))
cat(sprintf("Repair: %s\n", if (isTRUE(repair_enabled)) "enabled" else "disabled"))

start_time <- Sys.time()

workflow <- fit_predictive_reference_population_model(
  data_list = data_list,
  loglik_fn = loglik_shifted_gamma,
  population_model = population_model,
  theta_proposal = theta_proposal,
  local_control = list(
    M = local_particles,
    reference_strata = reference_strata,
    predictive_components = predictive_components,
    extra_reference_prior_fn = local_sketch_reference_fn,
    mbar_enabled = mbar_enabled,
    mbar_min_state_ess_frac = mbar_min_state_ess_frac,
    mbar_max_abs_shift = mbar_max_abs_shift,
    n_mcmc_moves = 2L,
    max_rounds = 100L,
    local_n_cores = 1L
  ),
  outer_control = list(
    N = outer_particles,
    n_mcmc_moves = outer_mcmc_moves,
    max_rounds = outer_max_rounds,
    n_cores = cores,
    seed = base_seed + 100000L,
    verbose = verbose
  ),
  audit_control = list(
    enabled = TRUE,
    n_q0 = audit_q0,
    n_outer = audit_outer,
    n_tail = audit_tail,
    n_q0_axis = audit_q0_axis,
    n_outer_axis = audit_outer_axis,
    outer_axis_inflate = outer_axis_inflate,
    min_ess_frac = audit_ess,
    max_pareto_k = max_pareto_k,
    max_loo_delta = max_loo_delta,
    compute_psis = compute_psis,
    compute_loo = compute_loo,
    adaptive_refinement = adaptive_refinement,
    pre_outer_support_refinement = pre_outer_support_refinement,
    pre_outer_support_rounds = pre_outer_support_rounds,
    max_pre_outer_states = max_pre_outer_states,
    max_adapt_rounds = max_adapt_rounds,
    max_new_states = max_new_states,
    max_new_states_per_local = max_new_states_per_local,
    theta_anchor_scales = theta_anchor_scales,
    min_refinement_impact_weight = min_refinement_impact_weight,
    min_surface_score = min_surface_score,
    state_reject_uncertainty_ratio = state_reject_uncertainty_ratio,
    state_reject_ess_ratio = state_reject_ess_ratio,
    protected_refinement_fraction = protected_refinement_fraction,
    protected_source_weights = protected_source_weights,
    min_protected_impact_weight = min_protected_impact_weight,
    adapt_particles = local_particles,
    n_jobs = cores,
    local_n_cores = 1L,
    repair_enabled = repair_enabled,
    max_repairs = max_repairs,
    repair_particles = local_particles,
    adaptive_outer_proposal = adaptive_outer_proposal,
    outer_proposal_defensive_weight = outer_proposal_defensive_weight
  ),
  n_cores = cores,
  seed = base_seed,
  verbose = verbose
)

elapsed_sec <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
fit <- workflow$fit

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
if (!is.null(workflow$audit)) {
  utils::write.csv(workflow$audit, audit_csv, row.names = FALSE)
}
if (!is.null(workflow$refinement_history) && nrow(workflow$refinement_history)) {
  utils::write.csv(workflow$refinement_history, refinements_csv, row.names = FALSE)
}

grDevices::png(plot_file, width = 1400, height = 900)
plot_posteriors(
  stan_draws,
  workflow_draws,
  labels = c("Stan", run_label),
  cols = c("black", "firebrick3"),
  n_cols = 3L
)
grDevices::dev.off()

audit_summary <- if (!is.null(workflow$audit)) {
  data.frame(
    failures = sum(!workflow$audit$covered),
    total = nrow(workflow$audit),
    min_ess_frac = min(workflow$audit$ess_frac, na.rm = TRUE),
    median_ess_frac = stats::median(workflow$audit$ess_frac, na.rm = TRUE),
    q10_ess_frac = as.numeric(stats::quantile(workflow$audit$ess_frac, 0.10, na.rm = TRUE, names = FALSE)),
    check.names = FALSE
  )
} else {
  NULL
}

workflow_state <- workflow
workflow_state$factor_set <- NULL

saveRDS(
  list(
    stan_source = stan_results_file,
    q0_source = q0_results_file,
    q0_method = q0_method,
    q0_method_weights = q0_method_weights,
    workflow = workflow_state,
    stan_draws = stan_draws,
    workflow_draws = workflow_draws,
    posterior_summary = summarize_population_posterior_diag(
      theta = fit$theta,
      w = fit$w,
      model = population_model
    ),
    posterior_comparison = posterior_comparison,
    audit_summary = audit_summary,
    settings = list(
      label = run_label,
      cores = cores,
      seed = base_seed,
      local_particles = local_particles,
      reference_strata = reference_strata,
      predictive_components = predictive_components,
      outer_particles = outer_particles,
      outer_mcmc_moves = outer_mcmc_moves,
      outer_max_rounds = outer_max_rounds,
      audit_q0 = audit_q0,
      audit_outer = audit_outer,
      audit_tail = audit_tail,
      audit_q0_axis = audit_q0_axis,
      audit_outer_axis = audit_outer_axis,
      outer_axis_inflate = outer_axis_inflate,
      audit_ess = audit_ess,
      compute_psis = compute_psis,
      compute_loo = compute_loo,
      max_pareto_k = max_pareto_k,
      max_loo_delta = max_loo_delta,
      adaptive_refinement = adaptive_refinement,
      pre_outer_support_refinement = pre_outer_support_refinement,
      pre_outer_support_rounds = pre_outer_support_rounds,
      max_pre_outer_states = max_pre_outer_states,
      max_adapt_rounds = max_adapt_rounds,
      max_new_states = max_new_states,
      max_new_states_per_local = max_new_states_per_local,
      theta_anchor_scales = theta_anchor_scales,
      min_refinement_impact_weight = min_refinement_impact_weight,
      min_surface_score = min_surface_score,
      state_reject_uncertainty_ratio = state_reject_uncertainty_ratio,
      state_reject_ess_ratio = state_reject_ess_ratio,
      protected_refinement_fraction = protected_refinement_fraction,
      protected_source_weights = protected_source_weights,
      min_protected_impact_weight = min_protected_impact_weight,
      repair_enabled = repair_enabled,
      max_repairs = max_repairs,
      adaptive_outer_proposal = adaptive_outer_proposal,
      outer_proposal_defensive_weight = outer_proposal_defensive_weight,
      local_sketch_enabled = local_sketch_enabled,
      sketch_starts = sketch_starts,
      sketch_scale = sketch_scale,
      mbar_enabled = mbar_enabled,
      mbar_min_state_ess_frac = mbar_min_state_ess_frac,
      mbar_max_abs_shift = mbar_max_abs_shift,
      elapsed_sec = elapsed_sec
    ),
    plot_file = plot_file
  ),
  results_file
)

cat("Saved results to:", results_file, "\n")
cat("Saved plot to:", plot_file, "\n")
cat("Saved posterior comparison to:", comparison_csv, "\n")
if (!is.null(workflow$audit)) cat("Saved local audit to:", audit_csv, "\n")
if (!is.null(workflow$refinement_history) && nrow(workflow$refinement_history)) {
  cat("Saved local refinements to:", refinements_csv, "\n")
}
cat(sprintf("Elapsed: %.1f seconds\n", elapsed_sec))
cat(sprintf("Outer beta: %.4f | rounds=%d | log evidence=%.4f +/- %.4f\n",
            fit$beta, fit$meta$rounds, fit$log_evidence, fit$mcse_log_evidence))
if (!is.null(audit_summary)) {
  cat(sprintf("Local audit failures: %d/%d | min ESS=%.4f | median ESS=%.4f\n",
              audit_summary$failures, audit_summary$total,
              audit_summary$min_ess_frac, audit_summary$median_ess_frac))
}
cat("\nPosterior comparison to Stan:\n")
print(posterior_comparison[order(abs(posterior_comparison$standardized_mean_error), decreasing = TRUE), ], row.names = FALSE)
