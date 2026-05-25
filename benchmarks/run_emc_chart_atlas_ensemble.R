#!/usr/bin/env Rscript

rm(list = ls())

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[1L]))
} else {
  normalizePath("benchmarks/run_emc_chart_atlas_ensemble.R")
}
repo_dir <- dirname(dirname(script_path))
setwd(repo_dir)

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

arg_num <- function(args, key, default) {
  val <- args[[key]]
  if (is.null(val) || !nzchar(val)) return(as.numeric(default))
  as.numeric(val)
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

cli_args <- parse_cli_args(commandArgs(trailingOnly = TRUE))
run_label <- arg_chr(cli_args, "label", "chart_atlas_validation_weighted_ensemble")
member_files <- arg_chr_vec(cli_args, "members", c(
  file.path("benchmarks", "results", "population_emc_chart_atlas_gate_full_refined_particle_results.rds"),
  file.path("benchmarks", "results", "population_emc_chart_atlas_high_local_nocal_results.rds")
))
weight_power <- arg_num(cli_args, "weight_power", 2)
n_draws <- arg_int(cli_args, "n_draws", NA_integer_)
seed <- arg_int(cli_args, "seed", 20260523L)

results_file <- arg_chr(
  cli_args,
  "results_file",
  file.path("benchmarks", "results", paste0("population_emc_", run_label, "_results.rds"))
)
comparison_csv <- arg_chr(
  cli_args,
  "comparison_csv",
  file.path("benchmarks", "results", paste0("population_emc_", run_label, "_posterior_comparison.csv"))
)
plot_file <- arg_chr(
  cli_args,
  "plot_file",
  file.path("benchmarks", "results", paste0("population_emc_", run_label, "_posteriors.png"))
)

source("local_charts.R")
source("utilities.R")

members <- lapply(member_files, readRDS)
names(members) <- sub("_results[.]rds$", "", basename(member_files))
if (!length(members)) {
  stop("At least one member result is required.")
}
if (any(vapply(members, function(x) is.null(x$workflow_draws) || is.null(x$gate), logical(1)))) {
  stop("Every member result must contain workflow_draws and gate.")
}

weights <- local_atlas_ensemble_weights_from_fresh_probes(
  lapply(members, `[[`, "gate"),
  power = weight_power
)
draw_count <- if (is.finite(n_draws)) as.integer(n_draws) else nrow(members[[1L]]$emc_draws)
workflow_draws <- local_atlas_ensemble_draws(
  lapply(members, `[[`, "workflow_draws"),
  weights = weights$ensemble_weight,
  n_draws = draw_count,
  seed = seed
)

reference_draws <- members[[1L]]$emc_draws
baseline_draws <- members[[1L]]$baseline_draws
posterior_comparison <- local_atlas_compare_posterior_draws(reference_draws, workflow_draws)
baseline_improvement <- local_atlas_compare_to_baseline(
  reference_draws = reference_draws,
  workflow_draws = workflow_draws,
  baseline_draws = baseline_draws
)

dir.create(dirname(results_file), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(comparison_csv), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(plot_file), showWarnings = FALSE, recursive = TRUE)
utils::write.csv(posterior_comparison, comparison_csv, row.names = FALSE)

grDevices::png(plot_file, width = 1800, height = 1400)
plot_posteriors(
  reference_draws,
  workflow_draws,
  labels = c("EMC2", run_label),
  cols = c("black", "firebrick3"),
  n_cols = 4L
)
grDevices::dev.off()

saveRDS(
  list(
    member_files = member_files,
    weights = weights,
    workflow_draws = workflow_draws,
    emc_draws = reference_draws,
    baseline_draws = baseline_draws,
    posterior_comparison = posterior_comparison,
    posterior_summary = .local_atlas_metric_summary(posterior_comparison),
    baseline_improvement = baseline_improvement,
    settings = list(
      label = run_label,
      weight_power = weight_power,
      n_draws = draw_count,
      seed = seed
    ),
    plot_file = plot_file,
    comparison_csv = comparison_csv
  ),
  results_file
)

cat(sprintf("Saved ensemble results: %s\n", results_file))
cat(sprintf("Saved posterior comparison: %s\n", comparison_csv))
cat(sprintf("Saved posterior plot: %s\n", plot_file))
cat("\nEnsemble weights:\n")
print(weights, row.names = FALSE)
cat("\nPosterior summary:\n")
print(.local_atlas_metric_summary(posterior_comparison), row.names = FALSE)
