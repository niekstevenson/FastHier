#!/usr/bin/env Rscript

rm(list = ls())

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[1L]))
} else {
  normalizePath("benchmarks/run_emc_nested_reference_strength.R")
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

arg_chr_vec <- function(args, key, default) {
  val <- args[[key]]
  if (is.null(val) || !nzchar(val)) return(as.character(default))
  trimws(strsplit(val, ",", fixed = TRUE)[[1L]])
}

arg_int_vec <- function(args, key, default) {
  as.integer(arg_chr_vec(args, key, as.character(default)))
}

logmeanexp <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  m <- max(x)
  m + log(mean(exp(x - m)))
}

finite_rmse <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x)) sqrt(mean(x^2)) else NA_real_
}

source("smc_core.R")
source("population_models.R")
source("local_charts.R")

cli_args <- parse_cli_args(commandArgs(trailingOnly = TRUE))
label <- arg_chr(cli_args, "label", "emc_nested_reference_strength")
seed <- arg_int(cli_args, "seed", 20260525L)
cores <- arg_int(cli_args, "cores", 2L)
data_file <- arg_chr(cli_args, "data_file", file.path("benchmarks", "samples", "full_EMC2.RData"))
source_results <- arg_chr(
  cli_args,
  "source_results",
  file.path("benchmarks", "results", "emc_local_surface_hypothesis_hard6_both_results.rds")
)
particles <- arg_int_vec(cli_args, "particles", c(384L, 768L, 1536L))
reps <- arg_int(cli_args, "reps", 3L)
max_pairs <- arg_int(cli_args, "max_pairs", 8L)
pairs_arg <- arg_chr(cli_args, "pairs", "")
target_cess <- arg_num(cli_args, "target_cess", 0.9)
n_mcmc_moves <- arg_int(cli_args, "n_mcmc_moves", 2L)
max_steps <- arg_int(cli_args, "max_steps", 176L)

out_prefix <- file.path("benchmarks", "results", label)
results_file <- arg_chr(cli_args, "results_file", paste0(out_prefix, "_results.rds"))
rows_csv <- arg_chr(cli_args, "rows_csv", paste0(out_prefix, "_rows.csv"))
summary_csv <- arg_chr(cli_args, "summary_csv", paste0(out_prefix, "_summary.csv"))
drift_csv <- arg_chr(cli_args, "drift_csv", paste0(out_prefix, "_drift.csv"))
plot_file <- arg_chr(cli_args, "plot_file", paste0(out_prefix, "_diagnostics.png"))
for (path in c(results_file, rows_csv, summary_csv, drift_csv, plot_file)) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
}

if (!file.exists(data_file)) stop("Missing EMC2 data file: ", data_file)
if (!file.exists(source_results)) stop("Missing source local-surface result: ", source_results)
load(data_file)
if (!exists("ELP_DDM", inherits = FALSE)) stop("The EMC2 data file must define ELP_DDM.")
source_obj <- readRDS(source_results)
if (is.null(source_obj$theta_design) || is.null(source_obj$rows)) {
  stop("source_results must contain theta_design and rows.")
}

emc <- ELP_DDM[[1L]]
data_list <- emc$data
model_factory <- emc$model
alpha_names <- emc$par_names
model <- source_obj$theta_design$population_model
theta <- source_obj$theta_design$theta

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

select_pairs <- function(rows, max_pairs) {
  base <- rows[rows$member == "post_outer", , drop = FALSE]
  if (!nrow(base)) stop("source rows contain no post_outer member.")
  base$score <- abs(base$centered_error_local) +
    pmax(base$fresh_log_m_uncertainty, 0, na.rm = TRUE) +
    ifelse(base$atlas_status != "certified", 0.5, 0) +
    pmax(base$particle_mis_psis_k - 0.7, 0, na.rm = TRUE)
  base <- base[order(-base$score), , drop = FALSE]
  out <- unique(base[, c("local", "theta_row", "theta_label", "theta_source", "score"), drop = FALSE])
  out[seq_len(min(nrow(out), as.integer(max_pairs))), , drop = FALSE]
}

if (nzchar(pairs_arg)) {
  parts <- strsplit(pairs_arg, ",", fixed = TRUE)[[1L]]
  selected_pairs <- do.call(rbind, lapply(parts, function(x) {
    z <- strsplit(trimws(x), ":", fixed = TRUE)[[1L]]
    if (length(z) != 2L) stop("--pairs entries must be local:theta_row.")
    data.frame(
      local = z[1L],
      theta_row = as.integer(z[2L]),
      theta_label = NA_character_,
      theta_source = NA_character_,
      score = NA_real_,
      check.names = FALSE
    )
  }))
} else {
  selected_pairs <- select_pairs(source_obj$rows, max_pairs = max_pairs)
}
selected_pairs <- selected_pairs[
  selected_pairs$local %in% names(data_list) &
    selected_pairs$theta_row >= 1L &
    selected_pairs$theta_row <= nrow(theta),
  ,
  drop = FALSE
]
if (!nrow(selected_pairs)) stop("No valid local/theta pairs selected.")

jobs <- expand.grid(
  pair_row = seq_len(nrow(selected_pairs)),
  particles = particles,
  replicate = seq_len(reps),
  KEEP.OUT.ATTRS = FALSE
)

cat(sprintf(
  "Nested reference strength: %d pairs x %d budgets x %d reps = %d SMC jobs\n",
  nrow(selected_pairs), length(particles), reps, nrow(jobs)
))
cat("Budgets:", paste(particles, collapse = ", "), "\n")
print(selected_pairs, row.names = FALSE)

run_job <- function(k) {
  pair <- selected_pairs[jobs$pair_row[k], , drop = FALSE]
  local_id <- as.character(pair$local)
  theta_row <- as.integer(pair$theta_row)
  M <- as.integer(jobs$particles[k])
  rep_id <- as.integer(jobs$replicate[k])
  run <- .local_chart_run_smc(
    local_id = local_id,
    theta_anchor = theta[theta_row, , drop = FALSE],
    data_i = data_list[[local_id]],
    loglik_fn = loglik_emc2,
    population_model = model,
    M = M,
    target_cess = target_cess,
    resample_threshold = 0.5,
    n_mcmc_moves = n_mcmc_moves,
    rw_scale = 0.75,
    G_mix = 8L,
    da_enable = TRUE,
    refit_every = 2L,
    max_steps = max_steps,
    deterministic_resampling = FALSE,
    n_cores = 1L,
    seed = seed + 1000003L * match(local_id, names(data_list)) +
      9176L * theta_row + 104729L * rep_id + 8191L * M,
    verbose = FALSE,
    source = "emc_nested_reference_strength"
  )
  diag <- run$diagnostics %||% list()
  data.frame(
    local = local_id,
    theta_row = theta_row,
    theta_label = as.character(pair$theta_label),
    theta_source = as.character(pair$theta_source),
    particles = M,
    replicate = rep_id,
    log_marginal = run$logZ,
    logZ_se = run$logZ_se,
    final_ess_frac = as.numeric(diag$final_ess_frac %||% NA_real_),
    min_path_ess_frac = as.numeric(diag$min_path_ess_frac %||% NA_real_),
    rounds = as.integer(diag$rounds %||% NA_integer_),
    check.names = FALSE
  )
}

rows_list <- if (cores <= 1L || nrow(jobs) <= 1L) {
  lapply(seq_len(nrow(jobs)), run_job)
} else {
  parallel::mclapply(seq_len(nrow(jobs)), run_job, mc.cores = min(cores, nrow(jobs)))
}
rows <- do.call(rbind, rows_list)

summary <- do.call(rbind, lapply(split(rows, list(rows$local, rows$theta_row, rows$particles), drop = TRUE), function(df) {
  x <- df$log_marginal[is.finite(df$log_marginal)]
  data.frame(
    local = df$local[1L],
    theta_row = df$theta_row[1L],
    theta_label = df$theta_label[1L],
    theta_source = df$theta_source[1L],
    particles = df$particles[1L],
    n = length(x),
    logZ_logmeanexp = logmeanexp(x),
    logZ_median = if (length(x)) stats::median(x) else NA_real_,
    logZ_mean = if (length(x)) mean(x) else NA_real_,
    replicate_sd = if (length(x) > 1L) stats::sd(x) else NA_real_,
    replicate_range = if (length(x) > 1L) diff(range(x)) else NA_real_,
    mean_path_se = mean(df$logZ_se, na.rm = TRUE),
    min_path_ess_frac = min(df$min_path_ess_frac, na.rm = TRUE),
    check.names = FALSE
  )
}))

max_budget <- max(particles)
baseline <- summary[summary$particles == max_budget, c("local", "theta_row", "logZ_logmeanexp")]
names(baseline)[3L] <- "max_budget_logZ"
drift <- merge(summary, baseline, by = c("local", "theta_row"), all.x = TRUE, sort = FALSE)
drift$delta_to_max_budget <- drift$logZ_logmeanexp - drift$max_budget_logZ
drift <- drift[order(drift$local, drift$theta_row, drift$particles), , drop = FALSE]

utils::write.csv(rows, rows_csv, row.names = FALSE)
utils::write.csv(summary, summary_csv, row.names = FALSE)
utils::write.csv(drift, drift_csv, row.names = FALSE)

budget_cols <- setNames(
  grDevices::hcl.colors(length(particles), "Dark 3"),
  as.character(sort(unique(particles)))
)
pair_id <- paste(selected_pairs$local, selected_pairs$theta_row, sep = ":")

grDevices::png(plot_file, width = 2300, height = 1500, res = 150, bg = "white")
old_par <- graphics::par(mfrow = c(2, 2), mar = c(7, 4, 3, 1), oma = c(0, 0, 1.7, 0))
plot_rows <- merge(rows, data.frame(local = selected_pairs$local, theta_row = selected_pairs$theta_row, pair_id = pair_id), by = c("local", "theta_row"))
xpos <- match(plot_rows$pair_id, pair_id) + (match(plot_rows$particles, sort(unique(particles))) - (length(particles) + 1) / 2) * 0.18
graphics::plot(
  xpos,
  plot_rows$log_marginal,
  pch = 21,
  bg = budget_cols[as.character(plot_rows$particles)],
  col = "white",
  xaxt = "n",
  xlab = "",
  ylab = "replicate log m",
  main = "Nested-SMC replicate cloud by budget"
)
graphics::axis(1, at = seq_along(pair_id), labels = pair_id, las = 2)
graphics::legend("topright", legend = names(budget_cols), fill = budget_cols, title = "particles", bty = "n")

graphics::boxplot(
  replicate_sd ~ particles,
  data = summary,
  xlab = "particles",
  ylab = "replicate SD of log m",
  main = "Replicate variation by budget",
  col = "grey85",
  border = "grey30"
)

drift_nonmax <- drift[drift$particles != max_budget, , drop = FALSE]
if (nrow(drift_nonmax)) {
  graphics::plot(
    drift_nonmax$particles,
    abs(drift_nonmax$delta_to_max_budget),
    log = "x",
    pch = 21,
    bg = budget_cols[as.character(drift_nonmax$particles)],
    col = "white",
    xlab = "particles",
    ylab = "|logmeanexp - max-budget logmeanexp|",
    main = "Cross-budget drift"
  )
  graphics::abline(h = c(0.1, 0.25, 0.5), col = c("grey70", "grey55", "grey40"), lty = 2)
} else {
  graphics::plot.new()
  graphics::title("Cross-budget drift unavailable")
  graphics::text(0.5, 0.5, "single particle budget")
}

source_rows <- source_obj$rows[source_obj$rows$member == "post_outer", , drop = FALSE]
source_rows$pair_id <- paste(source_rows$local, source_rows$theta_row, sep = ":")
source_rows <- source_rows[source_rows$pair_id %in% pair_id, , drop = FALSE]
graphics::plot(
  abs(source_rows$centered_error_local),
  source_rows$fresh_log_m_uncertainty,
  pch = 21,
  bg = "firebrick3",
  col = "white",
  xlab = "|atlas centered error| in source benchmark",
  ylab = "source reference uncertainty",
  main = "Was the old reference precise enough?"
)
graphics::abline(0, 1, lty = 2, col = "grey35")
graphics::text(
  abs(source_rows$centered_error_local),
  source_rows$fresh_log_m_uncertainty,
  labels = source_rows$pair_id,
  pos = 3,
  cex = 0.7
)
graphics::mtext("EMC nested local SMC reference-strength benchmark", outer = TRUE, font = 2)
graphics::par(old_par)
grDevices::dev.off()

saveRDS(
  list(
    rows = rows,
    summary = summary,
    drift = drift,
    selected_pairs = selected_pairs,
    source_results = source_results,
    settings = list(
      label = label,
      particles = particles,
      reps = reps,
      target_cess = target_cess,
      n_mcmc_moves = n_mcmc_moves,
      max_steps = max_steps,
      seed = seed
    )
  ),
  results_file
)

cat("Budget summary:\n")
print(aggregate(replicate_sd ~ particles, summary, function(x) c(median = stats::median(x, na.rm = TRUE), max = max(x, na.rm = TRUE))), row.names = FALSE)
cat("Drift summary to max budget:\n")
print(aggregate(abs(delta_to_max_budget) ~ particles, drift, function(x) c(median = stats::median(x, na.rm = TRUE), max = max(x, na.rm = TRUE))), row.names = FALSE)
cat("Saved rows:", rows_csv, "\n")
cat("Saved summary:", summary_csv, "\n")
cat("Saved drift:", drift_csv, "\n")
cat("Saved plot:", plot_file, "\n")
cat("Saved results:", results_file, "\n")
