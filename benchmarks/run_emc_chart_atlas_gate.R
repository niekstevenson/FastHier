#!/usr/bin/env Rscript

rm(list = ls())

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[1L]))
} else {
  normalizePath("benchmarks/run_emc_chart_atlas_gate.R")
}
repo_dir <- dirname(dirname(script_path))
setwd(repo_dir)

suppressPackageStartupMessages({
  library(parallel)
  library(EMC2)
})

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

announce_step <- function(step, title, detail = NULL) {
  cat(sprintf("\n=== %s: %s ===\n", step, title))
  if (!is.null(detail) && nzchar(detail)) {
    cat(detail, "\n")
  }
  flush.console()
}

arg_num_vec <- function(args, key, default) {
  val <- args[[key]]
  if (is.null(val) || !nzchar(val)) return(as.numeric(default))
  as.numeric(strsplit(val, ",", fixed = TRUE)[[1L]])
}

arg_chr_vec <- function(args, key, default) {
  val <- args[[key]]
  if (is.null(val) || !nzchar(val)) return(as.character(default))
  trimws(strsplit(val, ",", fixed = TRUE)[[1L]])
}

select_theta_profiles <- function(theta, focus_hyper_names, probs, max_points) {
  theta <- as.matrix(theta)
  focus_idx <- match(focus_hyper_names, colnames(theta))
  focus_idx <- focus_idx[is.finite(focus_idx)]
  rows <- integer(0)
  center <- apply(theta, 2L, stats::median)
  rows <- c(rows, which.min(rowSums(sweep(theta, 2L, center, "-")^2)))
  for (j in focus_idx) {
    qj <- stats::quantile(theta[, j], probs = probs, names = FALSE, type = 8, na.rm = TRUE)
    for (target in qj) {
      rows <- c(rows, which.min(abs(theta[, j] - target)))
    }
  }
  rows <- unique(rows)
  rows <- rows[seq_len(min(length(rows), as.integer(max_points)))]
  theta[rows, , drop = FALSE]
}

detected_cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
if (!is.finite(detected_cores) || detected_cores < 1L) detected_cores <- 1L
cli_args <- parse_cli_args(commandArgs(trailingOnly = TRUE))

run_label <- arg_chr(cli_args, "label", "chart_atlas_gate")
base_seed <- arg_int(cli_args, "base_seed", 20260523L)
cores <- arg_int(cli_args, "mc_cores", arg_int(cli_args, "cores", min(4L, detected_cores)))

data_file <- arg_chr(cli_args, "data_file", file.path("benchmarks", "samples", "full_EMC2.RData"))
baseline_file <- arg_chr(
  cli_args,
  "baseline_file",
  file.path("benchmarks", "results", "population_emc_framework_defaults_fixed_results.rds")
)
results_file <- arg_chr(
  cli_args,
  "results_file",
  file.path("benchmarks", "results", sprintf("population_emc_%s_results.rds", run_label))
)
plot_file <- arg_chr(
  cli_args,
  "plot_file",
  file.path("benchmarks", "results", sprintf("population_emc_%s_posteriors.png", run_label))
)
checkpoint_file <- arg_chr(
  cli_args,
  "checkpoint_file",
  file.path("benchmarks", "results", sprintf("population_emc_%s_checkpoint.rds", run_label))
)
resume_checkpoint <- arg_lgl(cli_args, "resume_checkpoint", FALSE)
comparison_csv <- arg_chr(
  cli_args,
  "comparison_csv",
  file.path("benchmarks", "results", sprintf("population_emc_%s_posterior_comparison.csv", run_label))
)
config_file <- arg_chr(
  cli_args,
  "config_file",
  file.path("benchmarks", "results", sprintf("population_emc_%s_config.rds", run_label))
)
build_history_csv <- arg_chr(
  cli_args,
  "build_history_csv",
  file.path("benchmarks", "results", sprintf("population_emc_%s_atlas_build_history.csv", run_label))
)
calibration_history_csv <- arg_chr(
  cli_args,
  "calibration_history_csv",
  file.path("benchmarks", "results", sprintf("population_emc_%s_calibration_history.csv", run_label))
)
graph_summary_csv <- arg_chr(
  cli_args,
  "graph_summary_csv",
  file.path("benchmarks", "results", sprintf("population_emc_%s_graph_summary.csv", run_label))
)
compression_summary_csv <- arg_chr(
  cli_args,
  "compression_summary_csv",
  file.path("benchmarks", "results", sprintf("population_emc_%s_compression_summary.csv", run_label))
)
smc_state_file <- arg_chr(
  cli_args,
  "smc_state_file",
  file.path("benchmarks", "results", sprintf("population_emc_%s_post_outer_smc.rds", run_label))
)

max_locals <- arg_int(cli_args, "max_locals", NA_integer_)
target_local_particles <- arg_int(cli_args, "target_local_particles", 6000L)
default_local_charts <- 4L
default_particles_per_chart <- max(256L, as.integer(round(target_local_particles / default_local_charts)))
root_particles <- arg_int(cli_args, "root_particles", default_particles_per_chart)
candidate_particles <- arg_int(cli_args, "candidate_particles", default_particles_per_chart)
local_mcmc_moves <- arg_int(cli_args, "local_mcmc_moves", 2L)
local_target_cess <- arg_num(cli_args, "local_target_cess", 0.90)
local_max_steps <- arg_int(cli_args, "local_max_steps", 128L)
max_anchors <- arg_int(cli_args, "max_anchors", 4L)
design_method <- arg_chr(cli_args, "design_method", "adaptive_scout")
if (!design_method %in% c("adaptive_scout", "hybrid", "profiles", "metric_farthest")) {
  stop("design_method must be one of: adaptive_scout, hybrid, profiles, metric_farthest.")
}
scout_min_families <- arg_int(cli_args, "scout_min_families", 1L)
scout_max_families <- arg_int(cli_args, "scout_max_families", 4L)
scout_relative_threshold <- arg_num(cli_args, "scout_relative_threshold", 0.25)
scout_min_anchors_per_local <- arg_int(cli_args, "scout_min_anchors_per_local", 3L)
scout_mean_anchors_per_local <- arg_int(cli_args, "scout_mean_anchors_per_local", 4L)
scout_max_anchors_per_local <- arg_int(cli_args, "scout_max_anchors_per_local", max_anchors)
axis_count <- arg_int(cli_args, "axis_count", 6L)
tail_probs <- arg_num_vec(cli_args, "tail_probs", c(0.05, 0.25, 0.5, 0.75, 0.95))
coverage_prob <- arg_num(cli_args, "coverage_prob", 0.995)
coverage_inflation <- arg_num(cli_args, "coverage_inflation", 1.35)
refine_rounds <- arg_int(cli_args, "refine_rounds", if (identical(design_method, "adaptive_scout")) 0L else 1L)
refine_points <- arg_int(cli_args, "refine_points", 4L)
refine_audit_n <- arg_int(cli_args, "refine_audit_n", arg_int(cli_args, "outer_particles", 2200L))
checkpoint_local_batch_size <- arg_int(cli_args, "checkpoint_local_batch_size", 8L)
strict_design_coverage <- arg_lgl(cli_args, "strict_design_coverage", TRUE)
max_chart_distance <- arg_num(cli_args, "max_chart_distance", NA_real_)
max_prediction_range <- arg_num(cli_args, "max_prediction_range", Inf)
edge_neighbors <- arg_int(cli_args, "edge_neighbors", 2L)
max_intermediates <- arg_int(cli_args, "max_intermediates", 2L)
min_overlap_ess <- arg_num(cli_args, "min_overlap_ess", 0.03)
edge_max_se <- arg_num(cli_args, "edge_max_se", 1.25)
edge_max_gap <- arg_num(cli_args, "edge_max_forward_reverse_gap", 1.25)
edge_max_taylor_gap <- arg_num(cli_args, "edge_max_taylor_gap", 3.0)
require_bar <- arg_lgl(cli_args, "require_bar_converged", TRUE)
use_particle_mis <- arg_lgl(cli_args, "use_particle_mis", TRUE)
require_particle_mis <- arg_lgl(cli_args, "require_particle_mis", TRUE)
min_particle_mis_ess <- arg_num(cli_args, "min_particle_mis_ess", 0.05)
min_particle_mis_ess_abs <- arg_num(cli_args, "min_particle_mis_ess_abs", 50)
max_particle_mis_psis_k <- arg_num(cli_args, "max_particle_mis_psis_k", 0.7)
max_quadratic_particle_gap <- arg_num(cli_args, "max_quadratic_particle_gap", Inf)
sparse_chart_min_covering <- arg_int(cli_args, "sparse_chart_min_covering", 3L)
sparse_chart_max_distance <- arg_num(cli_args, "sparse_chart_max_distance", Inf)
max_leave_chart_out_gap <- arg_num(cli_args, "max_leave_chart_out_gap", Inf)
particle_mis_batch <- arg_lgl(cli_args, "particle_mis_batch", TRUE)
stop_on_uncertified <- arg_lgl(cli_args, "stop_on_uncertified", FALSE)
use_uncertified_estimates <- arg_lgl(cli_args, "use_uncertified_estimates", TRUE)
outer_particles <- arg_int(cli_args, "outer_particles", 1200L)
outer_mcmc_moves <- arg_int(cli_args, "outer_mcmc_moves", 3L)
outer_max_rounds <- arg_int(cli_args, "outer_max_rounds", 90L)
outer_rw_scale_init <- arg_num(cli_args, "outer_rw_scale_init", 0.8)
proposal_core_scale <- arg_num(cli_args, "proposal_core_scale", 1.25)
proposal_tail_scale <- arg_num(cli_args, "proposal_tail_scale", 2.5)
proposal_tail_weight <- arg_num(cli_args, "proposal_tail_weight", 0.05)
calibration_rounds <- arg_int(cli_args, "calibration_rounds", 0L)
calibration_particles <- arg_int(cli_args, "calibration_particles", candidate_particles)
calibration_max_points <- arg_int(cli_args, "calibration_max_points", 7L)
calibration_max_updates <- arg_int(cli_args, "calibration_max_updates", 12L)
calibration_max_fresh_probes <- arg_int(cli_args, "calibration_max_fresh_probes", 64L)
calibration_abs_delta <- arg_num(cli_args, "calibration_abs_delta", 0.75)
calibration_z <- arg_num(cli_args, "calibration_z", 4)
calibration_candidate_pool_multiplier <- arg_int(cli_args, "calibration_candidate_pool_multiplier", 3L)
calibration_confirmation_reps <- arg_int(cli_args, "calibration_confirmation_reps", 0L)
calibration_confirmation_particles <- arg_int(cli_args, "calibration_confirmation_particles", calibration_particles)
calibration_confirmation_max_sd <- arg_num(cli_args, "calibration_confirmation_max_sd", 1.5)
calibration_adaptive_confirmation_reps <- arg_int(cli_args, "calibration_adaptive_confirmation_reps", 2L)
calibration_replicate_bootstrap_B <- arg_int(cli_args, "calibration_replicate_bootstrap_B", 200L)
calibration_max_direct_graph_z <- arg_num(cli_args, "calibration_max_direct_graph_z", 3)
calibration_max_direct_graph_chart_shift <- arg_num(cli_args, "calibration_max_direct_graph_chart_shift", 0.35)
calibration_max_direct_graph_existing_shift <- arg_num(cli_args, "calibration_max_direct_graph_existing_shift", 0.15)
calibration_adaptive_replicate_weight_multiplier <- arg_num(cli_args, "calibration_adaptive_replicate_weight_multiplier", 3)
calibration_adaptive_replicate_min_theta_weight <- arg_num(cli_args, "calibration_adaptive_replicate_min_theta_weight", 0)
calibration_adaptive_replicate_max_graph_z <- arg_num(cli_args, "calibration_adaptive_replicate_max_graph_z", 3)
calibration_adaptive_replicate_graph_shift <- arg_num(cli_args, "calibration_adaptive_replicate_graph_shift", 0.25)
pre_outer_calibration_rounds <- arg_int(cli_args, "pre_outer_calibration_rounds", 1L)
pre_outer_calibration_audit_n <- arg_int(cli_args, "pre_outer_calibration_audit_n", max(outer_particles, refine_audit_n))
pre_outer_calibration_max_points <- arg_int(cli_args, "pre_outer_calibration_max_points", 6L)
pre_outer_calibration_max_updates <- arg_int(cli_args, "pre_outer_calibration_max_updates", 16L)
pre_outer_calibration_max_fresh_probes <- arg_int(cli_args, "pre_outer_calibration_max_fresh_probes", 64L)
initial_certification_rounds <- arg_int(cli_args, "initial_certification_rounds", 1L)
initial_certification_max_updates <- arg_int(
  cli_args,
  "initial_certification_max_updates",
  max(pre_outer_calibration_max_updates, 16L)
)
initial_certification_max_fresh_probes <- arg_int(cli_args, "initial_certification_max_fresh_probes", 64L)
initial_certification_stop_on_uncertified <- arg_lgl(cli_args, "initial_certification_stop_on_uncertified", FALSE)
normalizer_robust_method <- arg_chr(cli_args, "normalizer_robust_method", "student_t")
normalizer_student_t_df <- arg_num(cli_args, "normalizer_student_t_df", 30)
if (!normalizer_robust_method %in% c("student_t", "huber", "none")) {
  stop("normalizer_robust_method must be one of: student_t, huber, none.")
}
verbose <- arg_lgl(cli_args, "verbose", FALSE)
progress_verbose <- arg_lgl(cli_args, "progress_verbose", TRUE)
trace_verbose <- arg_lgl(cli_args, "trace_verbose", verbose)

compress_local_particles <- arg_lgl(cli_args, "compress_local_particles", TRUE)
compression_particles <- arg_int(cli_args, "compression_particles", 64L)
compression_theta_points <- arg_int(cli_args, "compression_theta_points", 96L)
compression_holdout_fraction <- arg_num(cli_args, "compression_holdout_fraction", 0.25)
compression_max_holdout_rmse <- arg_num(cli_args, "compression_max_holdout_rmse", Inf)
compression_stop_on_failure <- arg_lgl(cli_args, "compression_stop_on_failure", FALSE)
compression_moment_weight <- arg_num(cli_args, "compression_moment_weight", 0.05)
compression_chart_weight <- arg_num(cli_args, "compression_chart_weight", 0.05)

focus_parameters <- arg_chr_vec(
  cli_args,
  "focus_parameters",
  c("mu_sv", "sigma2_sv", "sigma2_v_LogFreq", "sigma2_v")
)

source("local_charts.R")
source("utilities.R")

options(
  local_charts.normalizer_robust = normalizer_robust_method != "none",
  local_charts.normalizer_robust_method = normalizer_robust_method,
  local_charts.normalizer_student_t_df = normalizer_student_t_df
)

dir.create(dirname(results_file), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(plot_file), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(checkpoint_file), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(comparison_csv), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(config_file), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(build_history_csv), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(calibration_history_csv), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(graph_summary_csv), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(compression_summary_csv), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(smc_state_file), showWarnings = FALSE, recursive = TRUE)

announce_step("0/6", "Inputs")
if (!file.exists(data_file)) {
  stop("Missing EMC2 benchmark data: ", data_file)
}
if (!file.exists(baseline_file)) {
  stop("Missing baseline result: ", baseline_file)
}

load(data_file)
if (!exists("ELP_DDM", inherits = FALSE)) {
  stop("The EMC2 data file must define ELP_DDM.")
}
baseline <- readRDS(baseline_file)
if (is.null(baseline$workflow_draws)) {
  stop("Baseline result must contain workflow_draws.")
}

emc <- ELP_DDM[[1L]]
data_list <- emc$data
if (is.finite(max_locals)) {
  data_list <- data_list[seq_len(min(length(data_list), as.integer(max_locals)))]
}
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
baseline_draws <- baseline$workflow_draws

theta_cloud <- local_atlas_theta_from_draws(baseline_draws, population_model)
theta_root <- matrix(apply(theta_cloud, 2L, stats::median), nrow = 1L)
colnames(theta_root) <- population_model$hyper_names

focus_hyper_names <- unique(c(
  focus_parameters[startsWith(focus_parameters, "mu_")],
  sub("^sigma2_", "log_sigma2_", focus_parameters[startsWith(focus_parameters, "sigma2_")])
))
focus_hyper_names <- intersect(focus_hyper_names, population_model$hyper_names)
design_focus_hyper_names <- if (identical(design_method, "adaptive_scout")) NULL else focus_hyper_names

estimated_local_particles <- root_particles + candidate_particles * max(scout_mean_anchors_per_local - 1L, 0L)
cat(sprintf("Loaded EMC data: %d locals | %d alpha parameters\n", length(data_list), alpha_dim))
cat(sprintf("Baseline: %s\n", baseline_file))
cat(sprintf(
  "Chart atlas budget: target local particles=%d | estimated local particles=%d | root M=%d | candidate M=%d | mean anchors/local=%d | max anchors=%d | cores=%d\n",
  target_local_particles,
  estimated_local_particles,
  root_particles,
  candidate_particles,
  scout_mean_anchors_per_local,
  max_anchors,
  cores
))
cat(sprintf("Design method: %s\n", design_method))
cat(sprintf(
  "Compression: %s | K=%d | theta points=%d | holdout=%.2f | max holdout RMSE=%s\n",
  if (isTRUE(compress_local_particles)) "enabled" else "disabled",
  compression_particles,
  compression_theta_points,
  compression_holdout_fraction,
  if (is.finite(compression_max_holdout_rmse)) sprintf("%.3f", compression_max_holdout_rmse) else "Inf"
))
cat(sprintf(
  "Evaluation: particle-MIS=%s | min ESS frac=%.3f | min ESS abs=%.1f | max PSIS k=%.2f | calibration rounds=%d | normalizer=%s\n",
  if (isTRUE(require_particle_mis)) "required" else if (isTRUE(use_particle_mis)) "enabled" else "disabled",
  min_particle_mis_ess,
  min_particle_mis_ess_abs,
  max_particle_mis_psis_k,
  calibration_rounds,
  normalizer_robust_method
))
cat(sprintf(
  "Pre-outer certification: rounds=%d | audit theta=%d | max theta=%d | max updates=%d | max fresh=%d\n",
  pre_outer_calibration_rounds,
  pre_outer_calibration_audit_n,
  pre_outer_calibration_max_points,
  pre_outer_calibration_max_updates,
  pre_outer_calibration_max_fresh_probes
))
cat(sprintf(
  "Initial-cloud certification: rounds=%d | max updates=%d | max fresh=%d\n",
  initial_certification_rounds,
  initial_certification_max_updates,
  initial_certification_max_fresh_probes
))
cat(sprintf("Checkpoint: %s | resume=%s\n", checkpoint_file, if (isTRUE(resume_checkpoint)) "TRUE" else "FALSE"))
cat(sprintf(
  "Logging: progress=%s | particle_trace=%s | checkpoint_batch=%d\n",
  if (isTRUE(progress_verbose)) "TRUE" else "FALSE",
  if (isTRUE(trace_verbose)) "TRUE" else "FALSE",
  checkpoint_local_batch_size
))
cat("Focus posterior parameters:", paste(focus_parameters, collapse = ", "), "\n")
cat("Focus hyperparameters:", paste(focus_hyper_names, collapse = ", "), "\n")

saveRDS(
  list(
    label = run_label,
    started_at = Sys.time(),
    data_file = data_file,
    baseline_file = baseline_file,
    checkpoint_file = checkpoint_file,
    results_file = results_file,
    comparison_csv = comparison_csv,
    focus_parameters = focus_parameters,
    focus_hyper_names = focus_hyper_names,
    settings = as.list(cli_args)
  ),
  config_file
)
cat(sprintf("Saved run config: %s\n", config_file))

announce_step(
  "1/6",
  "Build and Certify Local Atlases",
  "Root charts, local scout-selected anchors, bridge/graph calibration, pre-outer certification, then particle compression if enabled."
)
start_time <- Sys.time()
workflow <- fit_chart_atlas_population_model(
  data_list = data_list,
  loglik_fn = loglik_emc2,
  population_model = population_model,
  theta_root = theta_root,
  theta_cloud = theta_cloud,
  local_control = list(
    root_M = root_particles,
    candidate_M = candidate_particles,
    target_cess = local_target_cess,
    n_mcmc_moves = local_mcmc_moves,
    max_steps = local_max_steps,
    root_confirm = "auto"
  ),
  design_control = list(
    max_anchors = max_anchors,
    design_method = design_method,
    scout_only = FALSE,
    stop_after_atlas_build = FALSE,
    axis_count = axis_count,
    tail_probs = tail_probs,
    focus_hyper_names = design_focus_hyper_names,
    coverage_prob = coverage_prob,
    coverage_inflation = coverage_inflation,
    refine_rounds = refine_rounds,
    refine_points = refine_points,
    refine_audit_n = refine_audit_n,
    checkpoint_local_batch_size = checkpoint_local_batch_size,
    strict_design_coverage = strict_design_coverage,
    scout_min_families = scout_min_families,
    scout_max_families = scout_max_families,
    scout_relative_threshold = scout_relative_threshold,
    scout_min_anchors_per_local = scout_min_anchors_per_local,
    scout_mean_anchors_per_local = scout_mean_anchors_per_local,
    scout_max_anchors_per_local = scout_max_anchors_per_local
  ),
  edge_control = list(
    edge_neighbors = edge_neighbors,
    max_intermediates = max_intermediates,
    min_overlap_ess = min_overlap_ess,
    max_se = edge_max_se,
    max_forward_reverse_gap = edge_max_gap,
    max_taylor_gap = edge_max_taylor_gap,
    require_bar_converged = require_bar
  ),
  evaluator_control = list(
    max_chart_distance = if (is.finite(max_chart_distance)) max_chart_distance else NULL,
    min_covering_charts = 3L,
    max_prediction_range = max_prediction_range,
    use_particle_mis = use_particle_mis,
    require_particle_mis = require_particle_mis,
    min_particle_mis_ess = min_particle_mis_ess,
    min_particle_mis_ess_abs = min_particle_mis_ess_abs,
    max_particle_mis_psis_k = max_particle_mis_psis_k,
    max_quadratic_particle_gap = max_quadratic_particle_gap,
    sparse_chart_min_covering = sparse_chart_min_covering,
    sparse_chart_max_distance = sparse_chart_max_distance,
    max_leave_chart_out_gap = max_leave_chart_out_gap,
    particle_mis_role = "estimator",
    particle_mis_batch = particle_mis_batch,
    compress_particle_mis = compress_local_particles,
    compressed_particle_mis_K = compression_particles,
    compressed_particle_mis_n_theta = compression_theta_points,
    compressed_particle_mis_holdout_fraction = compression_holdout_fraction,
    compressed_particle_mis_max_holdout_rmse = compression_max_holdout_rmse,
    compressed_particle_mis_stop_on_failure = compression_stop_on_failure,
    compressed_particle_mis_moment_weight = compression_moment_weight,
    compressed_particle_mis_chart_weight = compression_chart_weight,
    stop_on_uncertified = stop_on_uncertified,
    use_uncertified_estimates = use_uncertified_estimates
  ),
  calibration_control = list(
    rounds = calibration_rounds,
    max_points = calibration_max_points,
    probs = c(0.05, 0.5, 0.95),
    focus_hyper_names = focus_hyper_names,
    local_ids = seq_along(data_list),
    M = calibration_particles,
    target_cess = local_target_cess,
    n_mcmc_moves = local_mcmc_moves,
    max_steps = local_max_steps,
    max_updates = calibration_max_updates,
    max_fresh_probes = calibration_max_fresh_probes,
    abs_delta_threshold = calibration_abs_delta,
    z_threshold = calibration_z,
    candidate_pool_multiplier = calibration_candidate_pool_multiplier,
    confirmation_reps = calibration_confirmation_reps,
    confirmation_M = calibration_confirmation_particles,
    confirmation_abs_delta_threshold = calibration_abs_delta,
    confirmation_z_threshold = calibration_z,
    confirmation_max_sd = calibration_confirmation_max_sd,
    adaptive_confirmation_reps = calibration_adaptive_confirmation_reps,
    replicate_bootstrap_B = calibration_replicate_bootstrap_B,
    max_direct_graph_z = calibration_max_direct_graph_z,
    max_direct_graph_chart_shift = calibration_max_direct_graph_chart_shift,
    max_direct_graph_existing_shift = calibration_max_direct_graph_existing_shift,
    adaptive_replicate_weight_multiplier = calibration_adaptive_replicate_weight_multiplier,
    adaptive_replicate_min_theta_weight = calibration_adaptive_replicate_min_theta_weight,
    adaptive_replicate_max_graph_z = calibration_adaptive_replicate_max_graph_z,
    adaptive_replicate_graph_shift = calibration_adaptive_replicate_graph_shift,
    pre_outer_rounds = pre_outer_calibration_rounds,
    pre_outer_audit_n = pre_outer_calibration_audit_n,
    pre_outer_max_points = pre_outer_calibration_max_points,
    pre_outer_max_updates = pre_outer_calibration_max_updates,
    pre_outer_max_fresh_probes = pre_outer_calibration_max_fresh_probes,
    initial_certification_rounds = initial_certification_rounds,
    initial_certification_max_updates = initial_certification_max_updates,
    initial_certification_max_fresh_probes = initial_certification_max_fresh_probes,
    initial_certification_stop_on_uncertified = initial_certification_stop_on_uncertified
  ),
  proposal_control = list(
    max_components = 4L,
    core_weight = 1 - proposal_tail_weight,
    tail_weight = proposal_tail_weight,
    prior_weight = 0,
    core_scale = proposal_core_scale,
    tail_scale = proposal_tail_scale
  ),
  outer_control = list(
    N = outer_particles,
    n_mcmc_moves = outer_mcmc_moves,
    max_rounds = outer_max_rounds,
    rw_scale_init = outer_rw_scale_init,
    verbose = trace_verbose
  ),
  n_cores = cores,
  seed = base_seed,
  verbose = progress_verbose,
  trace_verbose = trace_verbose,
  checkpoint_file = checkpoint_file,
  resume_checkpoint = resume_checkpoint
)
elapsed_sec <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
if (is.data.frame(workflow$atlas_build_history) && nrow(workflow$atlas_build_history)) {
  utils::write.csv(workflow$atlas_build_history, build_history_csv, row.names = FALSE)
}
if (is.data.frame(workflow$calibration_history) && nrow(workflow$calibration_history)) {
  utils::write.csv(workflow$calibration_history, calibration_history_csv, row.names = FALSE)
}
if (is.data.frame(workflow$graph_summary) && nrow(workflow$graph_summary)) {
  utils::write.csv(workflow$graph_summary, graph_summary_csv, row.names = FALSE)
}
if (is.data.frame(workflow$compression_summary) && nrow(workflow$compression_summary)) {
  utils::write.csv(workflow$compression_summary, compression_summary_csv, row.names = FALSE)
}

announce_step(
  "2/6",
  "Store Post-Outer SMC State",
  "This is a compact inspection checkpoint for the final outer run; detailed stage checkpoints are written next to the main checkpoint file."
)
saveRDS(
  list(
    label = run_label,
    saved_at = Sys.time(),
    fit = workflow$fit,
    factor_set = workflow$factor_set,
    compression_summary = workflow$compression_summary,
    graph_summary = workflow$graph_summary,
    settings = workflow$settings
  ),
  smc_state_file,
  compress = FALSE
)

announce_step("3/6", "Draw Posterior Samples")
workflow_draws <- local_atlas_draws_from_fit(
  workflow$fit,
  population_model = population_model,
  n_draws = nrow(emc_draws),
  seed = base_seed + 1L
)
posterior_comparison <- local_atlas_compare_posterior_draws(emc_draws, workflow_draws)
posterior_summary <- .local_atlas_metric_summary(posterior_comparison)
utils::write.csv(posterior_comparison, comparison_csv, row.names = FALSE)

announce_step("4/6", "Plot Posterior Overlap")
grDevices::png(plot_file, width = 1800, height = 1400)
plot_posteriors(
  emc_draws,
  workflow_draws,
  labels = c("EMC2", run_label),
  cols = c("black", "firebrick3"),
  n_cols = 4L
)
grDevices::dev.off()

announce_step("5/6", "Save Benchmark Result Bundle")
saveRDS(
  list(
    data_source = data_file,
    baseline_file = baseline_file,
    workflow = workflow,
    emc_draws = emc_draws,
    workflow_draws = workflow_draws,
    baseline_draws = baseline_draws,
    posterior_comparison = posterior_comparison,
    graph_summary = workflow$graph_summary,
    design_certification = workflow$design_certification,
    calibration_history = workflow$calibration_history,
    settings = list(
      label = run_label,
      seed = base_seed,
      cores = cores,
      elapsed_sec = elapsed_sec,
      max_locals = max_locals,
      target_local_particles = target_local_particles,
      estimated_local_particles = estimated_local_particles,
      root_particles = root_particles,
      candidate_particles = candidate_particles,
      design_method = design_method,
      scout_min_families = scout_min_families,
      scout_max_families = scout_max_families,
      scout_relative_threshold = scout_relative_threshold,
      scout_min_anchors_per_local = scout_min_anchors_per_local,
      scout_mean_anchors_per_local = scout_mean_anchors_per_local,
      scout_max_anchors_per_local = scout_max_anchors_per_local,
      max_anchors = max_anchors,
      checkpoint_local_batch_size = checkpoint_local_batch_size,
      edge_neighbors = edge_neighbors,
      max_intermediates = max_intermediates,
      outer_particles = outer_particles,
      outer_rw_scale_init = outer_rw_scale_init,
      proposal_core_scale = proposal_core_scale,
      proposal_tail_scale = proposal_tail_scale,
      proposal_tail_weight = proposal_tail_weight,
      calibration_rounds = calibration_rounds,
      calibration_particles = calibration_particles,
      calibration_max_points = calibration_max_points,
      calibration_max_updates = calibration_max_updates,
      calibration_max_fresh_probes = calibration_max_fresh_probes,
      calibration_confirmation_reps = calibration_confirmation_reps,
      calibration_confirmation_particles = calibration_confirmation_particles,
      calibration_adaptive_confirmation_reps = calibration_adaptive_confirmation_reps,
      calibration_replicate_bootstrap_B = calibration_replicate_bootstrap_B,
      calibration_max_direct_graph_z = calibration_max_direct_graph_z,
      calibration_max_direct_graph_chart_shift = calibration_max_direct_graph_chart_shift,
      calibration_max_direct_graph_existing_shift = calibration_max_direct_graph_existing_shift,
      calibration_adaptive_replicate_weight_multiplier = calibration_adaptive_replicate_weight_multiplier,
      calibration_adaptive_replicate_min_theta_weight = calibration_adaptive_replicate_min_theta_weight,
      calibration_adaptive_replicate_max_graph_z = calibration_adaptive_replicate_max_graph_z,
      calibration_adaptive_replicate_graph_shift = calibration_adaptive_replicate_graph_shift,
      pre_outer_calibration_rounds = pre_outer_calibration_rounds,
      pre_outer_calibration_audit_n = pre_outer_calibration_audit_n,
      pre_outer_calibration_max_points = pre_outer_calibration_max_points,
      pre_outer_calibration_max_updates = pre_outer_calibration_max_updates,
      pre_outer_calibration_max_fresh_probes = pre_outer_calibration_max_fresh_probes,
      initial_certification_rounds = initial_certification_rounds,
      initial_certification_max_updates = initial_certification_max_updates,
      initial_certification_max_fresh_probes = initial_certification_max_fresh_probes,
      initial_certification_stop_on_uncertified = initial_certification_stop_on_uncertified,
      normalizer_robust_method = normalizer_robust_method,
      normalizer_student_t_df = normalizer_student_t_df,
      min_particle_mis_ess = min_particle_mis_ess,
      min_particle_mis_ess_abs = min_particle_mis_ess_abs,
      max_particle_mis_psis_k = max_particle_mis_psis_k,
      stop_on_uncertified = stop_on_uncertified,
      use_uncertified_estimates = use_uncertified_estimates,
      max_chart_distance = max_chart_distance,
      checkpoint_file = checkpoint_file,
      resume_checkpoint = resume_checkpoint,
      progress_verbose = progress_verbose,
      trace_verbose = trace_verbose,
      config_file = config_file,
      build_history_csv = build_history_csv,
      calibration_history_csv = calibration_history_csv,
      graph_summary_csv = graph_summary_csv
      ,
      compression_enabled = compress_local_particles,
      compression_particles = compression_particles,
      compression_theta_points = compression_theta_points,
      compression_holdout_fraction = compression_holdout_fraction,
      compression_max_holdout_rmse = compression_max_holdout_rmse,
      compression_stop_on_failure = compression_stop_on_failure,
      compression_summary_csv = compression_summary_csv,
      smc_state_file = smc_state_file
    ),
    plot_file = plot_file,
    comparison_csv = comparison_csv
  ),
  results_file
)

cat(sprintf("Saved results: %s\n", results_file))
cat(sprintf("Saved posterior comparison: %s\n", comparison_csv))
cat(sprintf("Saved posterior plot: %s\n", plot_file))
if (file.exists(compression_summary_csv)) cat(sprintf("Saved compression summary: %s\n", compression_summary_csv))
cat(sprintf("Saved post-outer SMC state: %s\n", smc_state_file))
cat("Validation gate not run by the fitting script; use benchmarks/validate_emc_chart_atlas_gate.R on saved results.\n")
cat("\nPosterior summary:\n")
print(posterior_summary, row.names = FALSE)
announce_step("6/6", "Complete")
