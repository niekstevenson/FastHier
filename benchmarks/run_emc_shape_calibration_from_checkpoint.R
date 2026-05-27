#!/usr/bin/env Rscript

rm(list = ls())

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[1L]))
} else {
  normalizePath("benchmarks/run_emc_shape_calibration_from_checkpoint.R")
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
  if (!is.null(detail) && nzchar(detail)) cat(detail, "\n")
  flush.console()
}

save_stage <- function(object, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(object, file)
  cat(sprintf("Saved: %s\n", file))
  flush.console()
}

write_stage_csv <- function(object, file) {
  if (!is.data.frame(object)) return(invisible(FALSE))
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(object, file, row.names = FALSE)
  cat(sprintf("Saved: %s\n", file))
  flush.console()
  invisible(TRUE)
}

source("local_charts.R")
source("outer_population_smc.R")
source("utilities.R")

cli_args <- parse_cli_args(commandArgs(trailingOnly = TRUE))
detected_cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
if (!is.finite(detected_cores) || detected_cores < 1L) detected_cores <- 1L

label <- arg_chr(cli_args, "label", "strict_repair_realistic_v1_fixed3_shape_phase7")
checkpoint_file <- arg_chr(
  cli_args,
  "checkpoint_file",
  file.path("benchmarks", "results", "population_emc_strict_repair_realistic_v1_fixed3_compressed_k256_checkpoint_post_outer.rds")
)
data_file <- arg_chr(cli_args, "data_file", file.path("benchmarks", "samples", "full_EMC2.RData"))
baseline_file <- arg_chr(
  cli_args,
  "baseline_file",
  file.path("benchmarks", "results", "population_emc_strict_repair_realistic_v1_fixed3_compressed_k256_results.rds")
)
out_dir <- arg_chr(cli_args, "out_dir", file.path("benchmarks", "results"))
seed <- arg_int(cli_args, "seed", 20260527L)
cores <- arg_int(cli_args, "cores", min(4L, detected_cores))

max_theta <- arg_int(cli_args, "max_theta", 96L)
tail_points_per_axis <- arg_int(cli_args, "tail_points_per_axis", 2L)
max_repair_pairs <- arg_int(cli_args, "max_repair_pairs", 12L)
max_holdout_pairs <- arg_int(cli_args, "max_holdout_pairs", 8L)
probe_M <- arg_int(cli_args, "probe_M", 192L)
probe_replicates <- arg_int(cli_args, "probe_replicates", 1L)
repair_M <- arg_int(cli_args, "repair_M", probe_M)
repair_replicates <- arg_int(cli_args, "repair_replicates", 1L)
holdout_M <- arg_int(cli_args, "holdout_M", probe_M)
holdout_replicates <- arg_int(cli_args, "holdout_replicates", probe_replicates)
bootstrap_B <- arg_int(cli_args, "bootstrap_B", 100L)
local_mcmc_moves <- arg_int(cli_args, "local_mcmc_moves", 2L)
local_target_cess <- arg_num(cli_args, "local_target_cess", 0.90)
local_max_steps <- arg_int(cli_args, "local_max_steps", 128L)
run_outer_if_required <- arg_lgl(cli_args, "run_outer_if_required", TRUE)
force_outer_rerun <- arg_lgl(cli_args, "force_outer_rerun", FALSE)
outer_particles <- arg_int(cli_args, "outer_particles", 1200L)
outer_mcmc_moves <- arg_int(cli_args, "outer_mcmc_moves", 3L)
outer_max_rounds <- arg_int(cli_args, "outer_max_rounds", 90L)
n_draws <- arg_int(cli_args, "n_draws", 3000L)
smc_verbose <- arg_lgl(cli_args, "smc_verbose", FALSE)

prefix <- file.path(out_dir, paste0("population_emc_", label))
config_file <- paste0(prefix, "_shape_config.rds")
theta_cloud_file <- paste0(prefix, "_shape_theta_cloud.rds")
raw_file <- paste0(prefix, "_shape_raw_certification.rds")
selection_file <- paste0(prefix, "_shape_selection.rds")
probe_file <- paste0(prefix, "_shape_probe.rds")
repair_file <- paste0(prefix, "_shape_repair.rds")
holdout_file <- paste0(prefix, "_shape_holdout.rds")
gate_file <- paste0(prefix, "_shape_reweight_gate.rds")
outer_file <- paste0(prefix, "_shape_outer_rerun.rds")
checkpoint_out_file <- paste0(prefix, "_checkpoint_post_shape.rds")
results_file <- paste0(prefix, "_results.rds")
comparison_csv <- paste0(prefix, "_posterior_comparison.csv")
raw_csv <- paste0(prefix, "_shape_raw_certification.csv")
selection_csv <- paste0(prefix, "_shape_selection.csv")
probe_csv <- paste0(prefix, "_shape_probe_residuals.csv")
repair_csv <- paste0(prefix, "_shape_repair_selected.csv")
holdout_csv <- paste0(prefix, "_shape_holdout_summary.csv")
gate_csv <- paste0(prefix, "_shape_gate_summary.csv")
plot_file <- paste0(prefix, "_posteriors.png")

select_outer_theta_cloud <- function(theta, w, max_theta, population_model, seed, tail_points_per_axis = 2L) {
  theta <- .as_hyper_matrix(theta, population_model$hyper_names, population_model$hyper_dim)
  w <- .local_chart_normalize_weights(w, nrow(theta))
  max_theta <- min(as.integer(max_theta), nrow(theta))
  if (max_theta < 1L) stop("max_theta must select at least one theta row.")

  set.seed(as.integer(seed))
  top_n <- min(nrow(theta), max(4L, floor(max_theta * 0.25)))
  rows <- head(order(w, decreasing = TRUE), top_n)

  if (tail_points_per_axis > 0L) {
    mu_idx <- grep("^mu_", colnames(theta))
    log_sigma_idx <- grep("^log_sigma2_", colnames(theta))
    focus_idx <- unique(c(mu_idx, log_sigma_idx))
    z_info <- .local_evidence_weighted_center_cov(theta, w, ridge = 1e-8)
    z <- sweep(theta, 2L, z_info$center, "-") %*% z_info$whitening
    axis_risk <- colSums(sweep(abs(z), 1L, w, "*"), na.rm = TRUE)
    axis_order <- focus_idx[order(axis_risk[focus_idx], decreasing = TRUE)]
    axis_order <- head(axis_order, max(1L, floor(max_theta * 0.15)))
    probs <- if (tail_points_per_axis <= 1L) 0.5 else seq(0.05, 0.95, length.out = tail_points_per_axis)
    for (j in axis_order) {
      qj <- as.numeric(.local_atlas_weighted_quantile(theta[, j], w, probs = probs))
      for (target in qj) {
        rows <- c(rows, which.min(abs(theta[, j] - target)))
      }
    }
  }

  remaining <- setdiff(seq_len(nrow(theta)), unique(rows))
  need <- max_theta - length(unique(rows))
  if (need > 0L && length(remaining)) {
    sampled <- sample(remaining, size = min(need, length(remaining)), replace = FALSE, prob = w[remaining])
    rows <- c(rows, sampled)
  }
  rows <- unique(rows)
  if (length(rows) > max_theta) {
    keep_top <- head(order(w[rows], decreasing = TRUE), max_theta)
    rows <- rows[keep_top]
  }
  rows <- rows[order(rows)]
  theta_selected <- theta[rows, , drop = FALSE]
  weights <- .local_chart_normalize_weights(w[rows], length(rows))
  build_local_evidence_certification_cloud(
    theta = theta_selected,
    population_model = population_model,
    theta_weights = weights,
    theta_source = "outer_post_shape"
  )
}

if (!file.exists(checkpoint_file)) stop("Missing checkpoint: ", checkpoint_file)
if (!file.exists(data_file)) stop("Missing EMC2 data: ", data_file)

announce_step("0/8", "Load Inputs", checkpoint_file)
checkpoint <- readRDS(checkpoint_file)
factor_set <- validate_local_atlas_factor_set(checkpoint$factor_set)
old_factor_set <- factor_set
fit <- checkpoint$fit
population_model <- factor_set$population_model
initial_proposal <- checkpoint$initial_proposal

load(data_file)
if (!exists("ELP_DDM", inherits = FALSE)) stop("The EMC2 data file must define ELP_DDM.")
emc <- ELP_DDM[[1L]]
data_list <- emc$data
model_factory <- emc$model
alpha_names <- emc$par_names

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

emc_mu <- as.data.frame(parameters(ELP_DDM, selection = "mu"), check.names = FALSE)
emc_sigma2 <- as.data.frame(parameters(ELP_DDM, selection = "sigma2"), check.names = FALSE)
colnames(emc_mu) <- paste0("mu_", alpha_names)
colnames(emc_sigma2) <- paste0("sigma2_", alpha_names)
emc_draws <- data.frame(emc_mu, emc_sigma2, check.names = FALSE)

baseline_draws <- NULL
if (file.exists(baseline_file)) {
  baseline <- readRDS(baseline_file)
  baseline_draws <- baseline$workflow_draws
}

save_stage(
  list(
    label = label,
    checkpoint_file = checkpoint_file,
    data_file = data_file,
    baseline_file = baseline_file,
    started_at = Sys.time(),
    settings = as.list(cli_args),
    resolved = list(
      seed = seed,
      cores = cores,
      max_theta = max_theta,
      max_repair_pairs = max_repair_pairs,
      max_holdout_pairs = max_holdout_pairs,
      probe_M = probe_M,
      repair_M = repair_M,
      holdout_M = holdout_M,
      outer_particles = outer_particles,
      run_outer_if_required = run_outer_if_required,
      force_outer_rerun = force_outer_rerun,
      smc_verbose = smc_verbose
    )
  ),
  config_file
)

local_control <- list(
  target_cess = local_target_cess,
  n_mcmc_moves = local_mcmc_moves,
  max_steps = local_max_steps
)

announce_step("1/8", "Posterior Theta Audit Cloud")
theta_cloud <- select_outer_theta_cloud(
  theta = fit$theta,
  w = fit$w,
  max_theta = max_theta,
  population_model = population_model,
  seed = seed,
  tail_points_per_axis = tail_points_per_axis
)
save_stage(theta_cloud, theta_cloud_file)
cat(sprintf("Theta audit rows: %d\n", nrow(theta_cloud$theta)))

announce_step("2/8", "Raw PMIS Certification")
raw_certification <- evaluate_raw_local_evidence_certification(
  factor_set = factor_set,
  cloud = theta_cloud,
  include_theta = TRUE,
  n_cores = cores
)
save_stage(raw_certification, raw_file)
write_stage_csv(raw_certification, raw_csv)
raw_summary <- summarize_raw_local_evidence_certification(raw_certification)
print(raw_summary$global, row.names = FALSE)

announce_step("3/8", "Shape Probe Selection")
selection <- select_shape_probe_pairs(
  factor_set = factor_set,
  raw_certification_table = raw_certification,
  graph_summary = checkpoint$graph_summary %||% local_atlas_graph_summary(factor_set),
  compression_summary = checkpoint$compression_summary %||% factor_set$compression_summary,
  control = list(
    max_repair_pairs = max_repair_pairs,
    max_holdout_pairs = max_holdout_pairs,
    max_repair_pairs_per_local = 2L,
    max_repair_pairs_per_theta = 3L,
    max_holdout_pairs_per_local = 2L,
    max_holdout_pairs_per_theta = 3L,
    certified_repair_fraction = 0.45,
    uncertified_repair_fraction = 0.30,
    tail_repair_fraction = 0.25,
    certified_holdout_fraction = 0.60,
    tail_holdout_fraction = 0.25
  )
)
save_stage(selection, selection_file)
write_stage_csv(selection$scored_table, selection_csv)
print(selection$selection_summary, row.names = FALSE)

announce_step("4/8", "Direct Shape Probes")
shape_probe <- run_shape_probe_pairs(
  factor_set = factor_set,
  selection = selection,
  data_list = data_list,
  loglik_fn = loglik_emc2,
  M = probe_M,
  n_replicates = probe_replicates,
  local_control = local_control,
  bootstrap_B = bootstrap_B,
  n_cores = cores,
  seed = seed + 10L,
  verbose = smc_verbose,
  geometry_control = list(
    max_directions = 3L,
    max_local_stencils = max(2L * max_repair_pairs, 8L),
    stencil_step = 0.65,
    local_stencil_step = 0.50,
    min_abs_standardized_residual = 1.0,
    min_abs_residual = 0.05,
    max_stencils_per_local = 2L
  )
)
save_stage(shape_probe, probe_file)
write_stage_csv(shape_probe$residuals, probe_csv)
print(shape_probe$summary, row.names = FALSE)

announce_step("5/8", "Strict Shape Repair")
shape_repair <- repair_shape_residual_geometry(
  factor_set = factor_set,
  shape_probe = shape_probe,
  data_list = data_list,
  loglik_fn = loglik_emc2,
  control = list(
    max_repairs = max_repair_pairs,
    M = repair_M,
    n_mcmc_moves = local_mcmc_moves,
    target_cess = local_target_cess,
    max_steps = local_max_steps,
    max_updates = max_repair_pairs,
    direct_confirmation_reps = repair_replicates,
    direct_confirmation_M = repair_M,
    audit_pre_repair = TRUE,
    audit_post_repair = TRUE,
    audit_scope = "repaired",
    stop_on_empty = FALSE,
    min_abs_standardized_residual = 1.0,
    min_abs_residual = 0.05
  ),
  local_control = local_control,
  edge_control = list(
    edge_neighbors = 2L,
    max_intermediates = 4L,
    min_overlap_ess = 0.03,
    max_se = 1.25,
    max_forward_reverse_gap = 1.25,
    max_taylor_gap = 3.0,
    require_bar_converged = TRUE
  ),
  n_cores = cores,
  seed = seed + 20L,
  verbose = smc_verbose
)
save_stage(shape_repair, repair_file)
write_stage_csv(shape_repair$selected_candidates, repair_csv)
print(shape_repair$summary, row.names = FALSE)

announce_step("6/8", "Holdout Validation")
holdout_validation <- validate_shape_repair_holdout(
  shape_repair = shape_repair,
  shape_probe = shape_probe,
  data_list = data_list,
  loglik_fn = loglik_emc2,
  control = list(
    M = holdout_M,
    n_replicates = holdout_replicates,
    bootstrap_B = bootstrap_B,
    require_holdout_improvement = FALSE,
    max_pair_centered_rmse_ratio = 1.10,
    max_local_centered_rmse_ratio = 1.10,
    max_total_centered_rmse_ratio = 1.10,
    max_abs_total_increase = 0.25,
    max_graph_edge_z_increase = 1.0
  ),
  local_control = local_control,
  n_cores = cores,
  seed = seed + 30L,
  verbose = smc_verbose
)
save_stage(holdout_validation, holdout_file)
write_stage_csv(holdout_validation$summary, holdout_csv)
print(holdout_validation$summary, row.names = FALSE)

announce_step("7/8", "Outer Reweight Gate")
gate <- shape_repair_outer_reweight_gate(
  shape_repair = shape_repair,
  fit = fit,
  old_factor_set = old_factor_set,
  holdout_validation = holdout_validation,
  reference_draws = NULL,
  baseline_draws = NULL,
  control = list(
    min_reweight_ess_fraction = 0.50,
    low_ess_rerun_fraction = 0.25,
    max_psis_k = 0.70,
    n_draws = min(n_draws, nrow(emc_draws)),
    require_holdout_acceptance = TRUE
  ),
  n_cores = cores,
  seed = seed + 40L
)
save_stage(gate, gate_file)
write_stage_csv(gate$summary, gate_csv)
print(gate$summary, row.names = FALSE)

final_fit <- fit
final_factor_set <- old_factor_set
final_source <- "original_post_outer"
outer_rerun <- NULL
if (isTRUE(force_outer_rerun) ||
    (isTRUE(run_outer_if_required) && isTRUE(gate$accepted) &&
       identical(gate$decision, "rerun_outer_required_low_reweight_quality"))) {
  announce_step("8/8", "Frozen Outer Rerun", "Gate says reweighting is not enough; running outer SMC against the repaired factor set.")
  outer_rerun <- local_atlas_frozen_outer_rerun(
    factor_set = shape_repair$factor_set,
    N = outer_particles,
    initial_proposal = initial_proposal,
    n_mcmc_moves = outer_mcmc_moves,
    max_rounds = outer_max_rounds,
    n_cores = cores,
    seed = seed + 50L,
    verbose = smc_verbose
  )
  final_fit <- outer_rerun
  final_factor_set <- shape_repair$factor_set
  final_source <- "frozen_outer_rerun_after_shape_repair"
  save_stage(outer_rerun, outer_file)
} else if (isTRUE(gate$accepted) && !is.null(gate$reweighted_fit)) {
  announce_step("8/8", "Use Reweighted Fit", gate$decision)
  final_fit <- gate$reweighted_fit
  final_factor_set <- shape_repair$factor_set
  final_source <- paste0("reweighted_fit:", gate$decision)
} else {
  announce_step("8/8", "Keep Original Fit", gate$decision)
}

draw_n <- min(as.integer(n_draws), nrow(emc_draws))
workflow_draws <- local_atlas_draws_from_fit(
  final_fit,
  population_model = population_model,
  n_draws = draw_n,
  seed = seed + 60L
)
posterior_comparison <- local_atlas_compare_posterior_draws(emc_draws, workflow_draws)
posterior_summary <- .local_atlas_metric_summary(posterior_comparison)
write_stage_csv(posterior_comparison, comparison_csv)
dir.create(dirname(plot_file), recursive = TRUE, showWarnings = FALSE)
grDevices::png(plot_file, width = 1800, height = 1400)
plot_posteriors(
  emc_draws,
  workflow_draws,
  labels = c("EMC2", final_source),
  cols = c("black", "firebrick3"),
  n_cols = 4L
)
grDevices::dev.off()
cat(sprintf("Saved: %s\n", plot_file))

post_shape_checkpoint <- checkpoint
post_shape_checkpoint$stage <- "post_shape_calibration"
post_shape_checkpoint$factor_set <- final_factor_set
post_shape_checkpoint$fit <- final_fit
post_shape_checkpoint$shape_theta_cloud <- theta_cloud
post_shape_checkpoint$shape_raw_certification_summary <- raw_summary
post_shape_checkpoint$shape_selection_summary <- selection$selection_summary
post_shape_checkpoint$shape_probe_summary <- shape_probe$summary
post_shape_checkpoint$shape_repair_summary <- shape_repair$summary
post_shape_checkpoint$shape_holdout_summary <- holdout_validation$summary
post_shape_checkpoint$shape_gate_summary <- gate$summary
post_shape_checkpoint$shape_final_source <- final_source
post_shape_checkpoint$checkpoint_time <- Sys.time()
save_stage(post_shape_checkpoint, checkpoint_out_file)

result <- list(
  source_checkpoint = checkpoint_file,
  final_source = final_source,
  theta_cloud = theta_cloud,
  raw_certification_summary = raw_summary,
  selection = selection,
  shape_probe = shape_probe,
  shape_repair = shape_repair,
  holdout_validation = holdout_validation,
  gate = gate,
  outer_rerun = outer_rerun,
  fit = final_fit,
  factor_set = final_factor_set,
  population_model = population_model,
  emc_draws = emc_draws,
  workflow_draws = workflow_draws,
  posterior_comparison = posterior_comparison,
  posterior_summary = posterior_summary,
  plot_file = plot_file,
  comparison_csv = comparison_csv,
  checkpoint_file = checkpoint_out_file,
  settings = list(
    label = label,
    seed = seed,
    cores = cores,
    max_theta = max_theta,
    max_repair_pairs = max_repair_pairs,
    max_holdout_pairs = max_holdout_pairs,
    probe_M = probe_M,
    probe_replicates = probe_replicates,
    repair_M = repair_M,
    repair_replicates = repair_replicates,
    holdout_M = holdout_M,
    holdout_replicates = holdout_replicates,
    run_outer_if_required = run_outer_if_required,
    force_outer_rerun = force_outer_rerun,
    outer_particles = outer_particles
  )
)
save_stage(result, results_file)

cat("\nPosterior summary:\n")
print(posterior_summary, row.names = FALSE)
cat(sprintf("\nFinal source: %s\n", final_source))
cat(sprintf("Results: %s\n", results_file))
cat(sprintf("Checkpoint: %s\n", checkpoint_out_file))
