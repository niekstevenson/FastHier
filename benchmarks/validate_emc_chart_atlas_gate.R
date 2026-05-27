#!/usr/bin/env Rscript

rm(list = ls())

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[1L]))
} else {
  normalizePath("benchmarks/validate_emc_chart_atlas_gate.R")
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

arg_chr <- function(args, key, default = NULL) {
  val <- args[[key]]
  if (is.null(val) || !nzchar(val)) default else as.character(val)
}

arg_int <- function(args, key, default) {
  val <- args[[key]]
  if (is.null(val) || !nzchar(val)) return(as.integer(default))
  as.integer(val)
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

cli_args <- parse_cli_args(commandArgs(trailingOnly = TRUE))
results_file <- arg_chr(cli_args, "results_file")
if (is.null(results_file) || !nzchar(results_file)) {
  stop("Provide --results_file=<saved chart atlas result RDS>.")
}
results_file <- normalizePath(results_file)
out_prefix <- sub("_results\\.rds$", "", results_file)
gate_file <- arg_chr(cli_args, "gate_file", paste0(out_prefix, "_validation_gate.rds"))
gate_csv <- arg_chr(cli_args, "gate_csv", paste0(out_prefix, "_gate_failures.csv"))
data_file <- arg_chr(cli_args, "data_file", file.path("benchmarks", "samples", "full_EMC2.RData"))
baseline_file <- arg_chr(
  cli_args,
  "baseline_file",
  file.path("benchmarks", "results", "population_emc_framework_defaults_fixed_results.rds")
)

detected_cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
if (!is.finite(detected_cores) || detected_cores < 1L) detected_cores <- 1L
cores <- arg_int(cli_args, "mc_cores", arg_int(cli_args, "cores", min(4L, detected_cores)))

source("local_charts.R")
source("utilities.R")

result <- readRDS(results_file)
workflow <- result$workflow %||% list(
  factor_set = result$factor_set,
  fit = result$fit,
  initial_proposal = result$initial_proposal,
  population_model = result$population_model
)
if (is.null(workflow) || is.null(workflow$factor_set) || is.null(workflow$fit)) {
  stop("results_file must contain workflow$factor_set/workflow$fit or top-level factor_set/fit.")
}
population_model <- workflow$population_model %||% workflow$factor_set$population_model

load(data_file)
if (!exists("ELP_DDM", inherits = FALSE)) {
  stop("The EMC2 data file must define ELP_DDM.")
}
emc <- ELP_DDM[[1L]]
data_list <- emc$data
alpha_names <- emc$par_names
model_factory <- emc$model
emc_draws <- result$emc_draws
if (is.null(emc_draws)) {
  emc_mu <- as.data.frame(parameters(ELP_DDM, selection = "mu"), check.names = FALSE)
  emc_sigma2 <- as.data.frame(parameters(ELP_DDM, selection = "sigma2"), check.names = FALSE)
  colnames(emc_mu) <- paste0("mu_", alpha_names)
  colnames(emc_sigma2) <- paste0("sigma2_", alpha_names)
  emc_draws <- data.frame(emc_mu, emc_sigma2, check.names = FALSE)
}
workflow_draws <- result$workflow_draws %||%
  local_atlas_draws_from_fit(
    workflow$fit,
    population_model = population_model,
    n_draws = nrow(emc_draws),
    seed = arg_int(cli_args, "draw_seed", 20260524L)
  )
baseline_draws <- result$baseline_draws
if (is.null(baseline_draws)) {
  baseline <- readRDS(baseline_file)
  baseline_draws <- baseline$workflow_draws
}
if (is.null(baseline_draws)) {
  stop("Provide baseline_draws in results_file or a baseline_file containing workflow_draws.")
}

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

settings <- result$settings %||% list()
focus_parameters <- arg_chr_vec(
  cli_args,
  "focus_parameters",
  settings$focus_parameters %||% c("mu_sv", "sigma2_sv", "sigma2_v_LogFreq", "sigma2_v")
)
focus_hyper_names <- unique(c(
  focus_parameters[startsWith(focus_parameters, "mu_")],
  sub("^sigma2_", "log_sigma2_", focus_parameters[startsWith(focus_parameters, "sigma2_")])
))
focus_hyper_names <- intersect(focus_hyper_names, population_model$hyper_names)

fresh_particles <- arg_int(cli_args, "fresh_particles", settings$fresh_particles %||% 900L)
fresh_max_points <- arg_int(cli_args, "fresh_max_points", 9L)
fresh_max_locals <- arg_int(cli_args, "fresh_max_locals", settings$fresh_max_locals %||% 6L)
frozen_outer_particles <- arg_int(cli_args, "frozen_outer_particles", settings$frozen_outer_particles %||% 1600L)
outer_mcmc_moves <- arg_int(cli_args, "outer_mcmc_moves", settings$outer_mcmc_moves %||% 3L)
outer_max_rounds <- arg_int(cli_args, "outer_max_rounds", settings$outer_max_rounds %||% 90L)
local_mcmc_moves <- arg_int(cli_args, "local_mcmc_moves", settings$local_mcmc_moves %||% 2L)
local_max_steps <- arg_int(cli_args, "local_max_steps", settings$local_max_steps %||% 128L)

theta_emc <- local_atlas_theta_from_draws(emc_draws, population_model)
fresh_theta <- .local_atlas_unique_theta(rbind(
  select_theta_profiles(workflow$fit$theta, focus_hyper_names, probs = c(0.05, 0.5, 0.95), max_points = fresh_max_points),
  select_theta_profiles(theta_emc, focus_hyper_names, probs = c(0.05, 0.5, 0.95), max_points = fresh_max_points)
), population_model)
fresh_theta <- fresh_theta[seq_len(min(nrow(fresh_theta), fresh_max_points)), , drop = FALSE]
fresh_local_ids <- seq_len(min(length(data_list), fresh_max_locals))

gate <- local_atlas_benchmark_gate(
  factor_set = workflow$factor_set,
  reference_draws = emc_draws,
  workflow_draws = workflow_draws,
  baseline_draws = baseline_draws,
  focus_parameters = focus_parameters,
  theta_audit = workflow$fit$theta,
  theta_audit_weights = workflow$fit$w,
  frozen_outer_control = list(
    N = frozen_outer_particles,
    initial_proposal = workflow$initial_proposal,
    n_mcmc_moves = outer_mcmc_moves,
    max_rounds = outer_max_rounds,
    verbose = FALSE
  ),
  fresh_probe_control = list(
    theta = fresh_theta,
    data_list = data_list,
    loglik_fn = loglik_emc2,
    local_ids = fresh_local_ids,
    M = fresh_particles,
    n_mcmc_moves = local_mcmc_moves,
    max_steps = local_max_steps
  ),
  thresholds = list(
    max_abs_standardized_mean_error = 0.50,
    mean_abs_standardized_mean_error = 0.20,
    max_uncertified_fraction = 0,
    max_uncertified_weight = 0,
    max_graph_standardized_residual = 3,
    max_fresh_probe_abs_z = 3,
    max_frozen_outer_mean_shift = 0.25,
    min_baseline_mean_error_improvement = 0,
    min_baseline_shape_error_improvement = 0,
    max_focus_parameter_worsening = 0
  ),
  n_cores = cores,
  seed = arg_int(cli_args, "seed", 20260523L)
)

dir.create(dirname(gate_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(gate_csv), recursive = TRUE, showWarnings = FALSE)
saveRDS(
  list(
    source_results_file = results_file,
    gate = gate,
    fresh_theta = fresh_theta,
    fresh_local_ids = fresh_local_ids,
    settings = list(
      focus_parameters = focus_parameters,
      fresh_particles = fresh_particles,
      fresh_max_points = fresh_max_points,
      fresh_max_locals = fresh_max_locals,
      frozen_outer_particles = frozen_outer_particles,
      outer_mcmc_moves = outer_mcmc_moves,
      outer_max_rounds = outer_max_rounds,
      local_mcmc_moves = local_mcmc_moves,
      local_max_steps = local_max_steps,
      cores = cores
    )
  ),
  gate_file
)
utils::write.csv(data.frame(failure = gate$failures, check.names = FALSE), gate_csv, row.names = FALSE)

cat(sprintf("Saved validation gate: %s\n", gate_file))
cat(sprintf("Saved gate failures: %s\n", gate_csv))
cat(sprintf("Gate passed: %s\n", if (isTRUE(gate$passed)) "TRUE" else "FALSE"))
if (!isTRUE(gate$passed)) {
  cat("Gate failures:", paste(gate$failures, collapse = ", "), "\n")
}
