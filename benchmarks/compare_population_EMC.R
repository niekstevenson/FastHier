rm(list = ls())

suppressPackageStartupMessages({
  library(parallel)
  library(EMC2)
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

compare_posterior_draws <- function(reference_draws, workflow_draws) {
  common <- intersect(names(reference_draws), names(workflow_draws))
  rows <- lapply(common, function(nm) {
    ref <- as.numeric(reference_draws[[nm]])
    wf <- as.numeric(workflow_draws[[nm]])
    q_ref <- stats::quantile(ref, probs = c(0.05, 0.5, 0.95), names = FALSE, type = 8, na.rm = TRUE)
    q_wf <- stats::quantile(wf, probs = c(0.05, 0.5, 0.95), names = FALSE, type = 8, na.rm = TRUE)
    ref_sd <- stats::sd(ref, na.rm = TRUE)
    data.frame(
      parameter = nm,
      reference_mean = mean(ref, na.rm = TRUE),
      workflow_mean = mean(wf, na.rm = TRUE),
      mean_error = mean(wf, na.rm = TRUE) - mean(ref, na.rm = TRUE),
      standardized_mean_error = (mean(wf, na.rm = TRUE) - mean(ref, na.rm = TRUE)) / max(ref_sd, .Machine$double.eps),
      reference_sd = ref_sd,
      workflow_sd = stats::sd(wf, na.rm = TRUE),
      sd_ratio = stats::sd(wf, na.rm = TRUE) / max(ref_sd, .Machine$double.eps),
      q05_error = q_wf[1L] - q_ref[1L],
      q50_error = q_wf[2L] - q_ref[2L],
      q95_error = q_wf[3L] - q_ref[3L],
      q_wasserstein = posterior_quantile_distance(wf, ref),
      reference_inside_workflow_q05_q95 = mean(ref >= q_wf[1L] & ref <= q_wf[3L], na.rm = TRUE),
      workflow_inside_reference_q05_q95 = mean(wf >= q_ref[1L] & wf <= q_ref[3L], na.rm = TRUE),
      check.names = FALSE
    )
  })
  do.call(rbind, rows)
}

detected_cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
if (!is.finite(detected_cores) || detected_cores < 1L) detected_cores <- 1L

cli_args <- parse_cli_args(commandArgs(trailingOnly = TRUE))

run_label <- arg_chr(cli_args, "label", "elp_bank_smc")
base_seed <- arg_int(cli_args, "base_seed", 20260522L)
cores <- arg_int(cli_args, "mc_cores", arg_int(cli_args, "cores", min(4L, detected_cores)))

local_particles <- arg_int(cli_args, "full_particles", arg_int(cli_args, "local_particles", 800L))
max_bank_nodes <- arg_int(cli_args, "max_bank_nodes", 10L)
target_local_ess <- arg_num(cli_args, "target_local_ess", 0.30)
local_mcmc_moves <- arg_int(cli_args, "local_mcmc_moves", 3L)
local_target_cess <- arg_num(cli_args, "local_target_cess", 0.90)
local_max_steps <- arg_int(cli_args, "local_max_steps", 128L)
bridge_min_single_ess <- arg_num(cli_args, "bridge_min_single_ess", Inf)
bridge_max_pareto_k <- arg_num(cli_args, "bridge_max_pareto_k", 0.70)
design_max_points <- arg_int(cli_args, "design_max_points", max(2L, max_bank_nodes - 1L))
paired_effect_profiles <- arg_lgl(cli_args, "paired_effect_profiles", TRUE)
design_profile_dims <- arg_int(cli_args, "design_profile_dims", 4L)
design_tail_probs <- arg_num_vec(cli_args, "design_tail_probs", c(0.025, 0.1, 0.9, 0.975))
initial_tail_probs <- arg_num_vec(cli_args, "initial_tail_probs", c(0.01, 0.025, 0.1, 0.9, 0.975, 0.99))
audit_max_points <- arg_int(cli_args, "audit_max_points", max(5L, max_bank_nodes - 1L))
audit_target_ess <- arg_num(cli_args, "audit_target_ess", target_local_ess)
audit_max_repairs <- arg_int(cli_args, "audit_max_repairs", NA_integer_)
audit_refine_rounds <- arg_int(cli_args, "audit_refine_rounds", 2L)
force_refine_points <- arg_int(cli_args, "force_refine_points", 0L)
force_profile_dims <- arg_int(cli_args, "force_profile_dims", 1L)
force_tail_probs <- arg_num_vec(cli_args, "force_tail_probs", c(0.01, 0.05, 0.95, 0.99))
anchor_weight <- arg_num(cli_args, "anchor_weight", 0.70)
anchor_scales <- arg_num_vec(cli_args, "anchor_scales", c(16, 4))
calibrate_normalizers <- arg_lgl(cli_args, "calibrate_normalizers", TRUE)
calibration_max_iter <- arg_int(cli_args, "calibration_max_iter", 3000L)
calibration_tol <- arg_num(cli_args, "calibration_tol", 1e-8)
calibration_anchor <- arg_chr(cli_args, "calibration_anchor", "mean_smc")
graph_certify <- arg_lgl(cli_args, "graph_certify", TRUE)
graph_min_edge_ess <- arg_num(cli_args, "graph_min_edge_ess", 0.05)
graph_max_pareto_k <- arg_num(cli_args, "graph_max_pareto_k", 0.70)
graph_use_psis <- arg_lgl(cli_args, "graph_use_psis", FALSE)
graph_max_rounds <- arg_int(cli_args, "graph_max_rounds", 50L)
graph_require_connected <- arg_lgl(cli_args, "graph_require_connected", TRUE)
surface_refine_points <- arg_int(cli_args, "surface_refine_points", 0L)
surface_refine_rounds <- arg_int(cli_args, "surface_refine_rounds", 1L)
surface_profile_dims <- arg_int(cli_args, "surface_profile_dims", design_profile_dims)
surface_tail_probs <- arg_num_vec(cli_args, "surface_tail_probs", c(0.05, 0.95))
surface_mean_tail_probs <- arg_num_vec(cli_args, "surface_mean_tail_probs", surface_tail_probs)
surface_variance_tail_probs <- arg_num_vec(cli_args, "surface_variance_tail_probs", surface_tail_probs)
surface_inflation_scale <- arg_num(cli_args, "surface_inflation_scale", 2.5)
surface_include_inflated_profile <- arg_lgl(cli_args, "surface_include_inflated_profile", TRUE)
surface_include_posterior_tail <- arg_lgl(cli_args, "surface_include_posterior_tail", TRUE)
surface_paired_effect_profiles <- arg_lgl(cli_args, "surface_paired_effect_profiles", TRUE)
surface_paired_effect_grid <- arg_lgl(cli_args, "surface_paired_effect_grid", TRUE)
surface_candidate_multiplier <- arg_int(cli_args, "surface_candidate_multiplier", 4L)
surface_score_local_count <- arg_int(cli_args, "surface_score_local_count", 20L)
surface_score_log_weight <- arg_num(cli_args, "surface_score_log_weight", 1)
surface_score_score_weight <- arg_num(cli_args, "surface_score_score_weight", 0.25)
surface_score_distance_weight <- arg_num(cli_args, "surface_score_distance_weight", 0.05)
rho_anchor_ladder <- arg_num_vec(cli_args, "rho_anchor_ladder", c(0.001, 0.01, 0.05, 0.15, 0.35, 0.75))
rho_anchor_enabled <- arg_lgl(cli_args, "rho_anchor_enabled", TRUE)
outer_particles <- arg_int(cli_args, "outer_particles", 2000L)
outer_mcmc_moves <- arg_int(cli_args, "outer_mcmc_moves", 4L)
outer_max_rounds <- arg_int(cli_args, "outer_max_rounds", 100L)
verbose <- arg_lgl(cli_args, "verbose", TRUE)

data_file <- arg_chr(cli_args, "data_file", file.path("benchmarks", "samples", "full_EMC2.RData"))
results_file <- arg_chr(cli_args, "results_file", file.path("benchmarks", "results", sprintf("population_emc_%s_results.rds", run_label)))
plot_file <- arg_chr(cli_args, "plot_file", file.path("benchmarks", "results", sprintf("population_emc_%s_posteriors.png", run_label)))
comparison_csv <- arg_chr(cli_args, "comparison_csv", file.path("benchmarks", "results", sprintf("population_emc_%s_posterior_comparison.csv", run_label)))

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

mean_prior_mean <- rep(0, alpha_dim)
mean_prior_var <- rep(1, alpha_dim)
sigma2_prior_shape <- rep(2, alpha_dim)
sigma2_prior_rate <- rep(0.3, alpha_dim)

population_model <- make_population_model_diag_gaussian(
  alpha_names = alpha_names,
  mean_prior_mean = mean_prior_mean,
  mean_prior_var = mean_prior_var,
  sigma2_prior_shape = sigma2_prior_shape,
  sigma2_prior_rate = sigma2_prior_rate,
  label = "emc_normal_gamma"
)

initial_theta <- matrix(c(mean_prior_mean, log(sigma2_prior_rate / (sigma2_prior_shape - 1))), nrow = 1L)
colnames(initial_theta) <- population_model$hyper_names

emc_mu <- as.data.frame(parameters(ELP_DDM, selection = "mu"), check.names = FALSE)
emc_sigma2 <- as.data.frame(parameters(ELP_DDM, selection = "sigma2"), check.names = FALSE)
colnames(emc_mu) <- paste0("mu_", alpha_names)
colnames(emc_sigma2) <- paste0("sigma2_", alpha_names)
emc_draws <- data.frame(emc_mu, emc_sigma2, check.names = FALSE)

if (!is.finite(audit_max_repairs)) {
  audit_max_repairs <- max(50L, length(data_list) * audit_max_points)
}

cat(sprintf("Loaded EMC2 benchmark data: %s\n", data_file))
cat(sprintf("Data: %d subjects | parameters: %d\n", length(data_list), alpha_dim))
cat(sprintf("Run: label=%s | cores=%d | seed=%d\n", run_label, cores, base_seed))
cat(sprintf("Bank budget: %d nodes x %d particles per subject\n", max_bank_nodes, local_particles))
cat(sprintf("Local SMC path: CESS %.3f | moves %d | max steps %d\n",
            local_target_cess, local_mcmc_moves, local_max_steps))
cat(sprintf("Bridge eligibility: single-node ESS %.3f | Pareto k %.2f\n",
            bridge_min_single_ess, bridge_max_pareto_k))
cat(sprintf("Design/audit: %d design points | %d audit points | audit ESS %.3f | refine rounds %d\n",
            design_max_points, audit_max_points, audit_target_ess, audit_refine_rounds))
cat(sprintf("Paired effect profiles: %s\n", if (isTRUE(paired_effect_profiles)) "enabled" else "disabled"))
cat(sprintf("Anchor candidate weight: %.3f\n", anchor_weight))
cat(sprintf("Normalizer calibration: %s\n", if (isTRUE(calibrate_normalizers)) calibration_anchor else "disabled"))
cat(sprintf(
  "Overlap graph certification: %s | edge ESS %.3f | Pareto k %.2f | require=%s\n",
  if (isTRUE(graph_certify)) "enabled" else "disabled",
  graph_min_edge_ess,
  graph_max_pareto_k,
  if (isTRUE(graph_require_connected)) "yes" else "no"
))
cat("Rho-sketch anchor ladder:", if (isTRUE(rho_anchor_enabled)) paste(rho_anchor_ladder, collapse = ", ") else "disabled", "\n")

start_time <- Sys.time()

bank_result <- fit_bank_smc_population_model(
  data_list = data_list,
  loglik_fn = loglik_emc2,
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
    max_steps = local_max_steps,
    bridge_min_single_ess_frac = bridge_min_single_ess,
    bridge_max_pareto_k = bridge_max_pareto_k
  ),
  anchor_control = list(
    anchor_weight = anchor_weight,
    scales = anchor_scales
  ),
  design_control = list(
    max_points = design_max_points,
    profile_dims = design_profile_dims,
    tail_probs = design_tail_probs,
    initial_tail_probs = initial_tail_probs,
    paired_effect_profiles = paired_effect_profiles
  ),
  calibration_control = list(
    enabled = calibrate_normalizers,
    max_iter = calibration_max_iter,
    tol = calibration_tol,
    anchor = calibration_anchor
  ),
  graph_control = list(
    enabled = graph_certify,
    min_edge_ess_frac = graph_min_edge_ess,
    max_pareto_k = graph_max_pareto_k,
    use_psis = graph_use_psis,
    max_rounds = graph_max_rounds,
    require_connected = graph_require_connected
  ),
  surface_control = list(
    enabled = surface_refine_points > 0L,
    max_points = surface_refine_points,
    refine_rounds = surface_refine_rounds,
    profile_dims = surface_profile_dims,
    tail_probs = surface_tail_probs,
    paired_mean_tail_probs = surface_mean_tail_probs,
    paired_variance_tail_probs = surface_variance_tail_probs,
    inflation_scale = surface_inflation_scale,
    include_inflated_profile = surface_include_inflated_profile,
    include_posterior_tail = surface_include_posterior_tail,
    paired_effect_profiles = surface_paired_effect_profiles,
    paired_effect_grid = surface_paired_effect_grid,
    candidate_multiplier = surface_candidate_multiplier,
    score_local_count = surface_score_local_count,
    score_log_weight = surface_score_log_weight,
    score_score_weight = surface_score_score_weight,
    score_distance_weight = surface_score_distance_weight
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
  n_draws = nrow(emc_mu),
  seed = base_seed + 1L,
  population_model = population_model
)
workflow_mu <- as.data.frame(workflow_parts$mu, check.names = FALSE)
workflow_sigma2 <- as.data.frame(workflow_parts$sigma2, check.names = FALSE)
colnames(workflow_mu) <- alpha_names
colnames(workflow_sigma2) <- alpha_names
workflow_draws <- data.frame(
  stats::setNames(workflow_mu, paste0("mu_", alpha_names)),
  stats::setNames(workflow_sigma2, paste0("sigma2_", alpha_names)),
  check.names = FALSE
)

posterior_comparison <- compare_posterior_draws(emc_draws, workflow_draws)
utils::write.csv(posterior_comparison, comparison_csv, row.names = FALSE)

grDevices::png(plot_file, width = 1800, height = 1400)
plot_posteriors(
  emc_draws,
  workflow_draws,
  labels = c("EMC2", run_label),
  cols = c("black", "firebrick3"),
  n_cols = 4L
)
grDevices::dev.off()

bank_state <- bank_result
bank_state$factor_set <- NULL
bank_state$population_model <- population_model

budget_summary <- bank_smc_bank_budget_summary(bank_result$banks)
calibration_summary <- bank_smc_calibration_summary(bank_result$banks)
overlap_graph_summary <- bank_smc_overlap_graph_summary(bank_result$banks)
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
surface_entries <- Filter(
  function(x) is.data.frame(x) && "predicted_anchor_error" %in% names(x),
  bank_result$surface_refinements
)
surface_summary <- lapply(surface_entries, function(x) {
  added <- x[as.logical(x$added), , drop = FALSE]
  finite_error <- added$predicted_anchor_error[is.finite(added$predicted_anchor_error)]
  data.frame(
    added = nrow(added),
    max_abs_prediction_error = if (length(finite_error)) max(abs(finite_error)) else NA_real_,
    median_abs_prediction_error = if (length(finite_error)) stats::median(abs(finite_error)) else NA_real_,
    check.names = FALSE
  )
})
surface_summary <- if (length(surface_summary)) do.call(rbind, surface_summary) else NULL

saveRDS(
  list(
    data_source = data_file,
    bank_state = bank_state,
    emc_draws = emc_draws,
    workflow_draws = workflow_draws,
    posterior_summary = summarize_population_posterior_diag(
      theta = fit$theta,
      w = fit$w,
      model = population_model
    ),
    posterior_comparison = posterior_comparison,
    bank_budget = budget_summary,
    calibration_summary = calibration_summary,
    overlap_graph_summary = overlap_graph_summary,
    audit_summary = audit_summary,
    surface_summary = surface_summary,
    settings = list(
      label = run_label,
      cores = cores,
      seed = base_seed,
      local_particles = local_particles,
      max_bank_nodes = max_bank_nodes,
      target_local_ess = target_local_ess,
      local_mcmc_moves = local_mcmc_moves,
      local_target_cess = local_target_cess,
      local_max_steps = local_max_steps,
      bridge_min_single_ess = bridge_min_single_ess,
      bridge_max_pareto_k = bridge_max_pareto_k,
      design_max_points = design_max_points,
      paired_effect_profiles = paired_effect_profiles,
      design_profile_dims = design_profile_dims,
      design_tail_probs = design_tail_probs,
      initial_tail_probs = initial_tail_probs,
      audit_max_points = audit_max_points,
      audit_target_ess = audit_target_ess,
      audit_max_repairs = audit_max_repairs,
      audit_refine_rounds = audit_refine_rounds,
      force_refine_points = force_refine_points,
      force_profile_dims = force_profile_dims,
      force_tail_probs = force_tail_probs,
      anchor_weight = anchor_weight,
      anchor_scales = anchor_scales,
      calibrate_normalizers = calibrate_normalizers,
      calibration_max_iter = calibration_max_iter,
      calibration_tol = calibration_tol,
      calibration_anchor = calibration_anchor,
      graph_certify = graph_certify,
      graph_min_edge_ess = graph_min_edge_ess,
      graph_max_pareto_k = graph_max_pareto_k,
      graph_use_psis = graph_use_psis,
      graph_max_rounds = graph_max_rounds,
      graph_require_connected = graph_require_connected,
      surface_refine_points = surface_refine_points,
      surface_refine_rounds = surface_refine_rounds,
      surface_profile_dims = surface_profile_dims,
      surface_tail_probs = surface_tail_probs,
      surface_mean_tail_probs = surface_mean_tail_probs,
      surface_variance_tail_probs = surface_variance_tail_probs,
      surface_inflation_scale = surface_inflation_scale,
      surface_include_inflated_profile = surface_include_inflated_profile,
      surface_include_posterior_tail = surface_include_posterior_tail,
      surface_paired_effect_profiles = surface_paired_effect_profiles,
      surface_paired_effect_grid = surface_paired_effect_grid,
      surface_candidate_multiplier = surface_candidate_multiplier,
      surface_score_local_count = surface_score_local_count,
      surface_score_log_weight = surface_score_log_weight,
      surface_score_score_weight = surface_score_score_weight,
      surface_score_distance_weight = surface_score_distance_weight,
      rho_anchor_enabled = rho_anchor_enabled,
      rho_anchor_ladder = rho_anchor_ladder,
      outer_particles = outer_particles,
      outer_mcmc_moves = outer_mcmc_moves,
      outer_max_rounds = outer_max_rounds,
      elapsed_sec = elapsed_sec
    ),
    plot_file = plot_file,
    comparison_csv = comparison_csv
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
if (isTRUE(calibrate_normalizers)) {
  cat(sprintf(
    "Normalizer calibration shifts: max=%.3f | median=%.3f | nonconverged=%d\n",
    max(calibration_summary$max_abs_shift_from_smc, na.rm = TRUE),
    stats::median(calibration_summary$max_abs_shift_from_smc, na.rm = TRUE),
    sum(!calibration_summary$converged, na.rm = TRUE)
  ))
}
if (isTRUE(graph_certify)) {
  cat(sprintf(
    "Overlap graph certification: certified=%d/%d | max components=%d | inserted=%d\n",
    sum(overlap_graph_summary$certified, na.rm = TRUE),
    nrow(overlap_graph_summary),
    max(overlap_graph_summary$n_components, na.rm = TRUE),
    sum(overlap_graph_summary$inserted, na.rm = TRUE)
  ))
}
if (!is.null(surface_summary)) {
  cat(sprintf(
    "Surface anchors: added=%d | max prediction error=%.3f | median prediction error=%.3f\n",
    sum(surface_summary$added, na.rm = TRUE),
    max(surface_summary$max_abs_prediction_error, na.rm = TRUE),
    stats::median(surface_summary$median_abs_prediction_error, na.rm = TRUE)
  ))
}
cat("\nPosterior comparison to EMC2:\n")
print(posterior_comparison[order(abs(posterior_comparison$standardized_mean_error), decreasing = TRUE), ], row.names = FALSE)
