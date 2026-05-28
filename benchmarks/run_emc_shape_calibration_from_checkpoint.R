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

arg_chr_list <- function(args, key, default = character()) {
  val <- args[[key]]
  if (is.null(val) || !nzchar(val)) return(default)
  out <- trimws(strsplit(as.character(val), ",", fixed = TRUE)[[1L]])
  out[nzchar(out)]
}

arg_int_list <- function(args, key, default = integer()) {
  val <- arg_chr_list(args, key, character())
  if (!length(val)) return(default)
  as.integer(val)
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
  invisible(file)
}

write_stage_csv <- function(object, file) {
  if (!is.data.frame(object)) return(invisible(FALSE))
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(object, file, row.names = FALSE)
  cat(sprintf("Saved: %s\n", file))
  flush.console()
  invisible(TRUE)
}

source("hierarchical_framework.R")

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
run_outer_if_required <- arg_lgl(cli_args, "run_outer_if_required", FALSE)
force_outer_rerun <- arg_lgl(cli_args, "force_outer_rerun", FALSE)
outer_particles <- arg_int(cli_args, "outer_particles", 1200L)
outer_mcmc_moves <- arg_int(cli_args, "outer_mcmc_moves", 3L)
outer_max_rounds <- arg_int(cli_args, "outer_max_rounds", 90L)
n_draws <- arg_int(cli_args, "n_draws", 3000L)
smc_verbose <- arg_lgl(cli_args, "smc_verbose", FALSE)
selector_raw <- tolower(arg_chr(cli_args, "selector", "risk"))
selector <- if (selector_raw %in% c("active", "active_feature", "active_feature_residual")) "active" else "risk"
residual_source_files <- arg_chr_list(cli_args, "residual_source_files", character())
active_round_ids <- arg_int_list(cli_args, "active_round_ids", integer())
active_exclude_probe_sets <- arg_chr_list(cli_args, "active_exclude_probe_sets", "holdout")
active_lambda_shape <- arg_num(cli_args, "active_lambda_shape", 1)
active_lambda_evidence <- arg_num(cli_args, "active_lambda_evidence", 0.10)
active_lambda_mean_shape <- arg_num(cli_args, "active_lambda_mean_shape", 0.10)
active_lambda_support <- arg_num(cli_args, "active_lambda_support", 0.10)
stop_after_probe <- arg_lgl(cli_args, "stop_after_probe", FALSE)
repair_executor_raw <- tolower(arg_chr(
  cli_args,
  "repair_executor",
  if (identical(selector, "active")) "selected_exact" else "geometry"
))
repair_executor <- if (repair_executor_raw %in% c("selected_exact", "active_exact", "probe_exact")) {
  "selected_exact"
} else {
  "geometry"
}

prefix <- file.path(out_dir, paste0("population_emc_", label))
config_file <- paste0(prefix, "_shape_config.rds")
theta_cloud_file <- paste0(prefix, "_shape_theta_cloud.rds")
raw_file <- paste0(prefix, "_shape_raw_certification.rds")
selection_file <- paste0(prefix, "_shape_selection.rds")
active_pool_file <- paste0(prefix, "_shape_active_residual_pool.rds")
active_feature_file <- paste0(prefix, "_shape_active_feature_model.rds")
active_feature_csv <- paste0(prefix, "_shape_active_feature_summary.csv")
active_acquisition_csv <- paste0(prefix, "_shape_active_acquisition.csv")
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

if (!file.exists(checkpoint_file)) stop("Missing checkpoint: ", checkpoint_file)
if (!file.exists(data_file)) stop("Missing EMC2 data: ", data_file)

announce_step("0/8", "Load Inputs", checkpoint_file)
checkpoint <- readRDS(checkpoint_file)
factor_set <- validate_local_atlas_factor_set(checkpoint$factor_set)
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

model <- define_hierarchical_model(
  data_list = data_list,
  loglik_fn = loglik_emc2,
  alpha_names = alpha_names,
  population_model = population_model,
  label = "emc_elpd_ddm"
)

emc_mu <- as.data.frame(parameters(ELP_DDM, selection = "mu"), check.names = FALSE)
emc_sigma2 <- as.data.frame(parameters(ELP_DDM, selection = "sigma2"), check.names = FALSE)
colnames(emc_mu) <- paste0("mu_", alpha_names)
colnames(emc_sigma2) <- paste0("sigma2_", alpha_names)
emc_draws <- data.frame(emc_mu, emc_sigma2, check.names = FALSE)

settings <- list(
  label = label,
  checkpoint_file = checkpoint_file,
  data_file = data_file,
  seed = seed,
  cores = cores,
  max_theta = max_theta,
  tail_points_per_axis = tail_points_per_axis,
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
  outer_particles = outer_particles,
  selector = selector,
  selector_raw = selector_raw,
  residual_source_files = residual_source_files,
  active_round_ids = active_round_ids,
  active_exclude_probe_sets = active_exclude_probe_sets,
  repair_executor = repair_executor,
  repair_executor_raw = repair_executor_raw,
  smc_verbose = smc_verbose
)
save_stage(
  list(
    label = label,
    started_at = Sys.time(),
    settings = as.list(cli_args),
    resolved = settings
  ),
  config_file
)

local_control <- list(
  target_cess = local_target_cess,
  n_mcmc_moves = local_mcmc_moves,
  max_steps = local_max_steps
)
selection_control <- list(
  max_repair_pairs = max_repair_pairs,
  max_holdout_pairs = max_holdout_pairs,
  max_repair_pairs_per_local = 2L,
  max_repair_pairs_per_theta = 3L,
  max_holdout_pairs_per_local = 2L,
  max_holdout_pairs_per_theta = 3L,
  lambda_shape = active_lambda_shape,
  lambda_evidence = active_lambda_evidence,
  lambda_mean_shape = active_lambda_mean_shape,
  lambda_support = active_lambda_support,
  certified_repair_fraction = 0.45,
  uncertified_repair_fraction = 0.30,
  tail_repair_fraction = 0.25,
  certified_holdout_fraction = 0.60,
  tail_holdout_fraction = 0.25
)
probe_control <- list(
  M = probe_M,
  n_replicates = probe_replicates,
  bootstrap_B = bootstrap_B,
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
repair_control <- list(
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
)
repair_edge_control <- list(
  edge_neighbors = 2L,
  max_intermediates = 4L,
  min_overlap_ess = 0.03,
  max_se = 1.25,
  max_forward_reverse_gap = 1.25,
  max_taylor_gap = 3.0,
  require_bar_converged = TRUE
)
holdout_control <- list(
  M = holdout_M,
  n_replicates = holdout_replicates,
  bootstrap_B = bootstrap_B,
  require_holdout_improvement = FALSE,
  max_pair_centered_rmse_ratio = 1.10,
  max_local_centered_rmse_ratio = 1.10,
  max_total_centered_rmse_ratio = 1.10,
  max_abs_total_increase = 0.25,
  max_graph_edge_z_increase = 1.0
)
gate_control <- list(
  min_reweight_ess_fraction = 0.50,
  low_ess_rerun_fraction = 0.25,
  max_psis_k = 0.70,
  n_draws = min(n_draws, nrow(emc_draws)),
  require_holdout_acceptance = TRUE,
  run_outer_if_required = run_outer_if_required,
  force_outer_rerun = force_outer_rerun,
  outer_particles = outer_particles,
  outer_mcmc_moves = outer_mcmc_moves,
  outer_max_rounds = outer_max_rounds
)

stage_callback <- function(stage, object) {
  switch(
    stage,
    theta_support = {
      announce_step("1/8", "Posterior Theta Audit Cloud")
      save_stage(object, theta_cloud_file)
      cat(sprintf("Theta audit rows: %d\n", nrow(object$theta)))
    },
    local_certification = {
      announce_step("2/8", "Raw PMIS Certification")
      save_stage(object, raw_file)
      write_stage_csv(object$table, raw_csv)
      print(object$summary$global, row.names = FALSE)
    },
    probe_selection = {
      announce_step("3/8", "Shape Probe Selection")
      save_stage(object, selection_file)
      write_stage_csv(object$selection$scored_table, selection_csv)
      if (!is.null(object$residual_pool)) save_stage(object$residual_pool, active_pool_file)
      if (!is.null(object$feature_model)) {
        save_stage(object$feature_model, active_feature_file)
        write_stage_csv(object$feature_model$summary, active_feature_csv)
        print(object$feature_model$summary, row.names = FALSE)
      }
      write_stage_csv(object$acquisition_table, active_acquisition_csv)
      print(object$selection$selection_summary, row.names = FALSE)
    },
    direct_probes = {
      announce_step("4/8", "Direct Shape Probes")
      save_stage(object, probe_file)
      write_stage_csv(object$residuals, probe_csv)
      print(object$summary, row.names = FALSE)
    },
    repair = {
      announce_step("5/8", "Strict Shape Repair")
      save_stage(object, repair_file)
      write_stage_csv(object$selected_candidates, repair_csv)
      print(object$summary, row.names = FALSE)
    },
    holdout = {
      announce_step("6/8", "Holdout Validation")
      save_stage(object, holdout_file)
      write_stage_csv(object$summary, holdout_csv)
      print(object$summary, row.names = FALSE)
    },
    outer_update = {
      announce_step("7/8", "Outer Reweight Gate")
      save_stage(object$gate, gate_file)
      write_stage_csv(object$gate$summary, gate_csv)
      print(object$gate$summary, row.names = FALSE)
      if (!is.null(object$outer_rerun)) {
        announce_step("8/8", "Frozen Outer Rerun")
        save_stage(object$outer_rerun, outer_file)
      } else if (startsWith(object$final_source, "reweighted_fit:")) {
        announce_step("8/8", "Use Reweighted Fit", object$final_source)
      } else {
        announce_step("8/8", "Keep Original Fit", object$final_source)
      }
    }
  )
  invisible(NULL)
}

calibration <- run_hierarchy_shape_calibration(
  model = model,
  factor_set = factor_set,
  fit = fit,
  initial_proposal = initial_proposal,
  max_theta = max_theta,
  tail_points_per_axis = tail_points_per_axis,
  selector = selector,
  residual_sources = residual_source_files,
  active_round_ids = active_round_ids,
  active_exclude_probe_sets = active_exclude_probe_sets,
  stop_after_probe = stop_after_probe,
  repair_executor = repair_executor,
  local_control = local_control,
  probe_control = probe_control,
  selection_control = selection_control,
  repair_control = repair_control,
  repair_edge_control = repair_edge_control,
  holdout_control = holdout_control,
  gate_control = gate_control,
  graph_summary = checkpoint$graph_summary,
  compression_summary = checkpoint$compression_summary %||% factor_set$compression_summary,
  n_cores = cores,
  seed = seed,
  verbose = smc_verbose,
  stage_callback = stage_callback
)

if (isTRUE(stop_after_probe)) {
  result <- list(
    source_checkpoint = checkpoint_file,
    final_source = calibration$final_source,
    theta_cloud = calibration$theta_support,
    raw_certification_summary = calibration$certification$summary,
    selection = calibration$probe_selection$selection,
    active_residual_pool = calibration$probe_selection$residual_pool,
    active_feature_model = calibration$probe_selection$feature_model,
    shape_probe = calibration$probes,
    population_model = population_model,
    settings = settings
  )
  save_stage(result, results_file)
  cat(sprintf("\nStopped after direct shape probes because --stop_after_probe=true.\nResults: %s\n", results_file))
  quit(save = "no", status = 0L)
}

draw_n <- min(as.integer(n_draws), nrow(emc_draws))
workflow_draws <- local_atlas_draws_from_fit(
  calibration$fit,
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
  labels = c("EMC2", calibration$final_source),
  cols = c("black", "firebrick3"),
  n_cols = 4L
)
grDevices::dev.off()
cat(sprintf("Saved: %s\n", plot_file))

post_shape_checkpoint <- checkpoint
post_shape_checkpoint$stage <- "post_shape_calibration"
post_shape_checkpoint$factor_set <- calibration$factor_set
post_shape_checkpoint$fit <- calibration$fit
post_shape_checkpoint$shape_theta_cloud <- calibration$theta_support
post_shape_checkpoint$shape_raw_certification_summary <- calibration$certification$summary
post_shape_checkpoint$shape_selection_summary <- calibration$probe_selection$selection$selection_summary
post_shape_checkpoint$shape_selector <- selector
post_shape_checkpoint$shape_active_feature_summary <- if (is.null(calibration$probe_selection$feature_model)) {
  NULL
} else {
  calibration$probe_selection$feature_model$summary
}
post_shape_checkpoint$shape_probe_summary <- calibration$probes$summary
post_shape_checkpoint$shape_repair_summary <- calibration$repair$summary
post_shape_checkpoint$shape_holdout_summary <- calibration$holdout$summary
post_shape_checkpoint$shape_gate_summary <- calibration$outer_update$gate$summary
post_shape_checkpoint$shape_final_source <- calibration$final_source
save_hierarchy_checkpoint(post_shape_checkpoint, checkpoint_out_file, stage = "post_shape_calibration")
cat(sprintf("Saved: %s\n", checkpoint_out_file))

result <- list(
  source_checkpoint = checkpoint_file,
  final_source = calibration$final_source,
  theta_cloud = calibration$theta_support,
  raw_certification_summary = calibration$certification$summary,
  selection = calibration$probe_selection$selection,
  active_residual_pool = calibration$probe_selection$residual_pool,
  active_feature_model = calibration$probe_selection$feature_model,
  shape_probe = calibration$probes,
  shape_repair = calibration$repair,
  holdout_validation = calibration$holdout,
  gate = calibration$outer_update$gate,
  outer_rerun = calibration$outer_update$outer_rerun,
  fit = calibration$fit,
  factor_set = calibration$factor_set,
  population_model = population_model,
  emc_draws = emc_draws,
  workflow_draws = workflow_draws,
  posterior_comparison = posterior_comparison,
  posterior_summary = posterior_summary,
  plot_file = plot_file,
  comparison_csv = comparison_csv,
  checkpoint_file = checkpoint_out_file,
  settings = settings
)
save_stage(result, results_file)

cat("\nPosterior summary:\n")
print(posterior_summary, row.names = FALSE)
cat(sprintf("\nFinal source: %s\n", calibration$final_source))
cat(sprintf("Results: %s\n", results_file))
cat(sprintf("Checkpoint: %s\n", checkpoint_out_file))
