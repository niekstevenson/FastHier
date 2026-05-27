#!/usr/bin/env Rscript

rm(list = ls())

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[1L]))
} else {
  normalizePath("benchmarks/run_emc_chart_design_scout.R")
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

detected_cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
if (!is.finite(detected_cores) || detected_cores < 1L) detected_cores <- 1L
cli_args <- parse_cli_args(commandArgs(trailingOnly = TRUE))

label <- arg_chr(cli_args, "label", "chart_design_scout")
stage <- arg_chr(cli_args, "stage", "scout")
if (!stage %in% c("scout", "atlas")) {
  stop("stage must be either 'scout' or 'atlas'.")
}
seed <- arg_int(cli_args, "seed", 20260525L)
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
  file.path("benchmarks", "results", sprintf("population_emc_%s_results.rds", label))
)
checkpoint_file <- arg_chr(
  cli_args,
  "checkpoint_file",
  file.path("benchmarks", "results", sprintf("population_emc_%s_checkpoint.rds", label))
)
family_csv <- arg_chr(
  cli_args,
  "family_csv",
  file.path("benchmarks", "results", sprintf("population_emc_%s_family_scores.csv", label))
)
local_csv <- arg_chr(
  cli_args,
  "local_csv",
  file.path("benchmarks", "results", sprintf("population_emc_%s_local_scores.csv", label))
)
local_family_csv <- arg_chr(
  cli_args,
  "local_family_csv",
  file.path("benchmarks", "results", sprintf("population_emc_%s_local_family_scores.csv", label))
)
budget_csv <- arg_chr(
  cli_args,
  "budget_csv",
  file.path("benchmarks", "results", sprintf("population_emc_%s_local_anchor_budgets.csv", label))
)
build_history_csv <- arg_chr(
  cli_args,
  "build_history_csv",
  file.path("benchmarks", "results", sprintf("population_emc_%s_atlas_build_history.csv", label))
)
design_certification_csv <- arg_chr(
  cli_args,
  "design_certification_csv",
  file.path("benchmarks", "results", sprintf("population_emc_%s_design_certification.csv", label))
)
graph_summary_csv <- arg_chr(
  cli_args,
  "graph_summary_csv",
  file.path("benchmarks", "results", sprintf("population_emc_%s_graph_summary.csv", label))
)

max_locals <- arg_int(cli_args, "max_locals", NA_integer_)
root_particles <- arg_int(cli_args, "root_particles", 256L)
candidate_particles <- arg_int(cli_args, "candidate_particles", root_particles)
local_mcmc_moves <- arg_int(cli_args, "local_mcmc_moves", 2L)
local_max_steps <- arg_int(cli_args, "local_max_steps", 128L)
local_target_cess <- arg_num(cli_args, "local_target_cess", 0.9)
max_anchors <- arg_int(cli_args, "max_anchors", 9L)
edge_neighbors <- arg_int(cli_args, "edge_neighbors", 2L)
max_intermediates <- arg_int(cli_args, "max_intermediates", 2L)
min_overlap_ess <- arg_num(cli_args, "min_overlap_ess", 0.03)
edge_max_se <- arg_num(cli_args, "edge_max_se", 1.25)
edge_max_forward_reverse_gap <- arg_num(cli_args, "edge_max_forward_reverse_gap", 1.25)
edge_max_taylor_gap <- arg_num(cli_args, "edge_max_taylor_gap", 3.0)
require_bar_converged <- arg_lgl(cli_args, "require_bar_converged", TRUE)
min_particle_mis_ess <- arg_num(cli_args, "min_particle_mis_ess", 0.05)
min_particle_mis_ess_abs <- arg_num(cli_args, "min_particle_mis_ess_abs", 50)
max_particle_mis_psis_k <- arg_num(cli_args, "max_particle_mis_psis_k", 0.7)
strict_design_coverage <- arg_lgl(cli_args, "strict_design_coverage", TRUE)
scout_min_families <- arg_int(cli_args, "scout_min_families", 1L)
scout_max_families <- arg_int(cli_args, "scout_max_families", 4L)
scout_relative_threshold <- arg_num(cli_args, "scout_relative_threshold", 0.25)
scout_min_anchors_per_local <- arg_int(cli_args, "scout_min_anchors_per_local", 3L)
scout_mean_anchors_per_local <- arg_int(cli_args, "scout_mean_anchors_per_local", 6L)
scout_max_anchors_per_local <- arg_int(cli_args, "scout_max_anchors_per_local", max_anchors)
verbose <- arg_lgl(cli_args, "verbose", TRUE)
trace_verbose <- arg_lgl(cli_args, "trace_verbose", FALSE)
resume_checkpoint <- arg_lgl(cli_args, "resume_checkpoint", TRUE)

source("local_charts.R")
source("utilities.R")

dir.create(dirname(results_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(checkpoint_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(family_csv), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(local_csv), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(local_family_csv), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(budget_csv), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(build_history_csv), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(design_certification_csv), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(graph_summary_csv), recursive = TRUE, showWarnings = FALSE)

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
alpha_names <- emc$par_names
model_factory <- emc$model
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
theta_cloud <- local_atlas_theta_from_draws(baseline$workflow_draws, population_model)
theta_root <- matrix(apply(theta_cloud, 2L, stats::median), nrow = 1L)
colnames(theta_root) <- population_model$hyper_names

cat(sprintf(
  "Chart design %s: locals=%d root_M=%d candidate_M=%d cores=%d\n",
  stage,
  length(data_list),
  root_particles,
  candidate_particles,
  cores
))
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
    design_method = "adaptive_scout",
    scout_only = identical(stage, "scout"),
    stop_after_atlas_build = identical(stage, "atlas"),
    max_anchors = max_anchors,
    refine_rounds = 0L,
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
    max_forward_reverse_gap = edge_max_forward_reverse_gap,
    max_taylor_gap = edge_max_taylor_gap,
    require_bar_converged = require_bar_converged
  ),
  evaluator_control = list(
    min_particle_mis_ess = min_particle_mis_ess,
    min_particle_mis_ess_abs = min_particle_mis_ess_abs,
    max_particle_mis_psis_k = max_particle_mis_psis_k,
    stop_on_uncertified = TRUE,
    use_uncertified_estimates = FALSE
  ),
  n_cores = cores,
  seed = seed,
  verbose = verbose,
  trace_verbose = trace_verbose,
  checkpoint_file = checkpoint_file,
  resume_checkpoint = resume_checkpoint
)
elapsed_sec <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
plan <- workflow$chart_design_plan
utils::write.csv(plan$family_scores, family_csv, row.names = FALSE)
utils::write.csv(plan$local_scores, local_csv, row.names = FALSE)
utils::write.csv(plan$local_family_scores, local_family_csv, row.names = FALSE)
utils::write.csv(plan$local_anchor_budgets, budget_csv, row.names = FALSE)
if (is.data.frame(workflow$atlas_build_history) && nrow(workflow$atlas_build_history)) {
  utils::write.csv(workflow$atlas_build_history, build_history_csv, row.names = FALSE)
}
if (is.data.frame(workflow$design_certification) && nrow(workflow$design_certification)) {
  utils::write.csv(workflow$design_certification, design_certification_csv, row.names = FALSE)
}
if (is.data.frame(workflow$graph_summary) && nrow(workflow$graph_summary)) {
  utils::write.csv(workflow$graph_summary, graph_summary_csv, row.names = FALSE)
}
saveRDS(
  list(
    workflow = workflow,
    chart_design_plan = plan,
    settings = list(
      label = label,
      stage = stage,
      elapsed_sec = elapsed_sec,
      seed = seed,
      cores = cores,
      max_locals = max_locals,
      root_particles = root_particles,
      candidate_particles = candidate_particles,
      max_anchors = max_anchors,
      edge_neighbors = edge_neighbors,
      max_intermediates = max_intermediates,
      min_overlap_ess = min_overlap_ess,
      edge_max_se = edge_max_se,
      edge_max_forward_reverse_gap = edge_max_forward_reverse_gap,
      edge_max_taylor_gap = edge_max_taylor_gap,
      require_bar_converged = require_bar_converged,
      min_particle_mis_ess = min_particle_mis_ess,
      min_particle_mis_ess_abs = min_particle_mis_ess_abs,
      max_particle_mis_psis_k = max_particle_mis_psis_k,
      strict_design_coverage = strict_design_coverage,
      scout_min_families = scout_min_families,
      scout_max_families = scout_max_families,
      scout_relative_threshold = scout_relative_threshold,
      scout_min_anchors_per_local = scout_min_anchors_per_local,
      scout_mean_anchors_per_local = scout_mean_anchors_per_local,
      scout_max_anchors_per_local = scout_max_anchors_per_local
    )
  ),
  results_file
)

cat(sprintf("Saved %s result: %s\n", stage, results_file))
cat(sprintf("Saved family scores: %s\n", family_csv))
if (identical(stage, "atlas")) {
  cat(sprintf("Saved atlas build history: %s\n", build_history_csv))
  cat(sprintf("Saved design certification: %s\n", design_certification_csv))
  cat(sprintf("Saved graph summary: %s\n", graph_summary_csv))
}
cat("Selected families:", paste(plan$selected_families, collapse = ", "), "\n")
if (!is.null(workflow$factor_set)) {
  graph <- local_atlas_graph_summary(workflow$factor_set)
  cat(sprintf("Frozen factor set: locals=%d active_charts=%d\n",
              workflow$factor_set$n_locals,
              sum(graph$n_active_charts)))
}
cat(sprintf("Elapsed seconds: %.1f\n", elapsed_sec))
