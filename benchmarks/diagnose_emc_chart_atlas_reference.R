#!/usr/bin/env Rscript

rm(list = ls())

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[1L]))
} else {
  normalizePath("benchmarks/diagnose_emc_chart_atlas_reference.R")
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

cli_args <- parse_cli_args(commandArgs(trailingOnly = TRUE))
label <- arg_chr(cli_args, "label", "chart_atlas_reference_diagnostic")
member_files <- arg_chr_vec(cli_args, "members", c(
  file.path("benchmarks", "results", "population_emc_chart_atlas_gate_full_refined_particle_results.rds"),
  file.path("benchmarks", "results", "population_emc_chart_atlas_high_local_nocal_results.rds")
))
data_file <- arg_chr(cli_args, "data_file", file.path("benchmarks", "samples", "full_EMC2.RData"))
cores <- arg_int(cli_args, "cores", min(4L, parallel::detectCores(logical = TRUE)))
theta_max_points <- arg_int(cli_args, "theta_max_points", 9L)
theta_probs <- as.numeric(strsplit(arg_chr(cli_args, "theta_probs", "0.05,0.25,0.5,0.75,0.95"), ",", fixed = TRUE)[[1L]])
screen_particles <- arg_int(cli_args, "screen_particles", 180L)
confirm_particles <- arg_int(cli_args, "confirm_particles", 450L)
confirm_reps <- arg_int(cli_args, "confirm_reps", 3L)
confirm_top_pairs <- arg_int(cli_args, "confirm_top_pairs", 24L)
local_mcmc_moves <- arg_int(cli_args, "local_mcmc_moves", 1L)
local_max_steps <- arg_int(cli_args, "local_max_steps", 84L)
focus_parameters <- arg_chr_vec(
  cli_args,
  "focus_parameters",
  c("mu_sv", "sigma2_sv", "sigma2_v_LogFreq", "sigma2_v")
)

results_file <- arg_chr(
  cli_args,
  "results_file",
  file.path("benchmarks", "results", paste0("population_emc_", label, "_results.rds"))
)
rows_csv <- arg_chr(
  cli_args,
  "rows_csv",
  file.path("benchmarks", "results", paste0("population_emc_", label, "_rows.csv"))
)
theta_csv <- arg_chr(
  cli_args,
  "theta_csv",
  file.path("benchmarks", "results", paste0("population_emc_", label, "_theta_summary.csv"))
)
local_csv <- arg_chr(
  cli_args,
  "local_csv",
  file.path("benchmarks", "results", paste0("population_emc_", label, "_local_summary.csv"))
)

source("local_charts.R")

if (!file.exists(data_file)) {
  stop("Missing EMC2 benchmark data: ", data_file)
}
load(data_file)
if (!exists("ELP_DDM", inherits = FALSE)) {
  stop("The EMC2 data file must define ELP_DDM.")
}

members <- lapply(member_files, readRDS)
names(members) <- sub("^population_emc_", "", sub("_results[.]rds$", "", basename(member_files)))
if (any(vapply(members, function(x) is.null(x$workflow) || is.null(x$workflow$factor_set), logical(1)))) {
  stop("Every member result must contain workflow$factor_set.")
}

emc <- ELP_DDM[[1L]]
data_list <- emc$data
model_factory <- emc$model
alpha_names <- emc$par_names
model <- members[[1L]]$workflow$population_model

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
theta_emc <- local_atlas_theta_from_draws(emc_draws, model)

focus_hyper_names <- unique(c(
  focus_parameters[startsWith(focus_parameters, "mu_")],
  sub("^sigma2_", "log_sigma2_", focus_parameters[startsWith(focus_parameters, "sigma2_")])
))
focus_hyper_names <- intersect(focus_hyper_names, model$hyper_names)
theta_design <- build_local_evidence_calibration_design(
  theta = theta_emc,
  population_model = model,
  weights = rep(1 / nrow(theta_emc), nrow(theta_emc)),
  focus_hyper_names = focus_hyper_names,
  probs = theta_probs,
  max_points = theta_max_points,
  source = "emc_reference_benchmark",
  label_prefix = "emc_ref"
)
theta_ref <- theta_design$theta

cat(sprintf(
  "Reference diagnostic: %d members | %d locals | %d EMC theta points | screen M=%d\n",
  length(members), length(data_list), nrow(theta_ref), screen_particles
))

factor_sets <- lapply(members, function(x) x$workflow$factor_set)
audit <- build_local_evidence_audit(
  factor_sets = factor_sets,
  theta = theta_design,
  data_list = data_list,
  loglik_fn = loglik_emc2,
  local_ids = seq_along(data_list),
  n_replicates = 1L,
  M = screen_particles,
  local_control = list(n_mcmc_moves = local_mcmc_moves, max_steps = local_max_steps),
  n_cores = cores,
  seed = 8100001L,
  reference_source = "emc_reference_benchmark"
)
selected_pairs <- select_audit_failures(
  audit,
  max_pairs = confirm_top_pairs,
  min_abs_error = 0,
  min_abs_z = 0
)
if (confirm_top_pairs > 0L && confirm_reps > 0L && nrow(selected_pairs)) {
  audit <- extend_local_evidence_audit(
    audit,
    data_list = data_list,
    loglik_fn = loglik_emc2,
    pairs = selected_pairs,
    n_replicates = confirm_reps,
    M = confirm_particles,
    local_control = list(n_mcmc_moves = local_mcmc_moves, max_steps = local_max_steps),
    n_cores = cores,
    seed = 9100003L
  )
}

rows <- audit$rows
rows$reference_log_marginal <- rows$fresh_log_m_center
rows$reference_se <- rows$fresh_log_m_uncertainty
rows$reference_source <- ifelse(rows$n_replicates > 1L, "confirmed", "screen")
rows$screen_log_marginal <- rows$fresh_log_m_center
rows$screen_se <- rows$fresh_log_m_path_se
rows$confirm_log_marginal <- ifelse(rows$n_replicates > 1L, rows$fresh_log_m_center, NA_real_)
rows$confirm_sd <- ifelse(rows$n_replicates > 1L, rows$fresh_log_m_sd, NA_real_)
rows$error_atlas_minus_reference <- rows$atlas_log_marginal - rows$reference_log_marginal
confirm <- unique(rows[rows$n_replicates > 1L, c(
  "local", "local_pos", "theta_row", "theta_label", "reference_log_marginal",
  "reference_se", "fresh_log_m_sd", "n_replicates"
)])
names(confirm)[names(confirm) == "reference_log_marginal"] <- "confirm_log_marginal"
names(confirm)[names(confirm) == "fresh_log_m_sd"] <- "confirm_sd"

theta_frame <- audit$theta_metadata

theta_summary <- do.call(rbind, lapply(split(rows, rows$member), function(df) {
  out <- aggregate(
    df$error_atlas_minus_reference,
    by = list(member = df$member, theta_row = df$theta_row, theta_label = df$theta_label),
    FUN = .local_evidence_finite_sum
  )
  names(out)[4L] <- "total_log_surface_error"
  stats <- aggregate(
    abs(df$error_atlas_minus_reference),
    by = list(member = df$member, theta_row = df$theta_row, theta_label = df$theta_label),
    FUN = .local_evidence_finite_mean
  )
  names(stats)[4L] <- "mean_abs_local_error"
  uncert <- aggregate(
    df$atlas_status != "certified",
    by = list(member = df$member, theta_row = df$theta_row, theta_label = df$theta_label),
    FUN = sum
  )
  names(uncert)[4L] <- "n_uncertified_locals"
  out <- merge(out, stats, by = c("member", "theta_row", "theta_label"))
  out <- merge(out, uncert, by = c("member", "theta_row", "theta_label"))
  out$centered_total_log_surface_error <- out$total_log_surface_error -
    .local_evidence_finite_mean(out$total_log_surface_error)
  out
}))
theta_summary <- merge(theta_summary, theta_frame, by = c("theta_row", "theta_label"), all.x = TRUE)

local_summary <- aggregate(
  rows$error_atlas_minus_reference,
  by = list(member = rows$member, local = rows$local, local_pos = rows$local_pos),
  FUN = .local_evidence_finite_rmse
)
names(local_summary)[4L] <- "rmse_local_error"
local_mean <- aggregate(
  rows$error_atlas_minus_reference,
  by = list(member = rows$member, local = rows$local, local_pos = rows$local_pos),
  FUN = .local_evidence_finite_mean
)
names(local_mean)[4L] <- "mean_signed_error"
local_max <- aggregate(
  abs(rows$error_atlas_minus_reference),
  by = list(member = rows$member, local = rows$local, local_pos = rows$local_pos),
  FUN = .local_evidence_finite_max
)
names(local_max)[4L] <- "max_abs_error"
local_summary <- Reduce(function(x, y) merge(x, y, by = c("member", "local", "local_pos")),
                        list(local_summary, local_mean, local_max))
local_summary <- local_summary[order(local_summary$member, -local_summary$rmse_local_error, na.last = TRUE), , drop = FALSE]

focus_regression <- do.call(rbind, lapply(split(theta_summary, theta_summary$member), function(df) {
  rows_out <- lapply(focus_hyper_names, function(name) {
    x <- as.numeric(df[[name]])
    y <- as.numeric(df$centered_total_log_surface_error)
    sx <- stats::sd(x)
    sy <- stats::sd(y)
    data.frame(
      member = df$member[1L],
      hyper = name,
      correlation = if (is.finite(sx) && sx > 0 && is.finite(sy) && sy > 0) stats::cor(x, y) else NA_real_,
      standardized_slope = if (is.finite(sx) && sx > 0) as.numeric(coef(stats::lm(y ~ x))[2L] * sx) else NA_real_,
      check.names = FALSE
    )
  })
  do.call(rbind, rows_out)
}))

member_summary <- aggregate(
  rows$error_atlas_minus_reference,
  by = list(member = rows$member),
  FUN = .local_evidence_finite_rmse
)
names(member_summary)[2L] <- "local_error_rmse"
member_mae <- aggregate(
  abs(rows$error_atlas_minus_reference),
  by = list(member = rows$member),
  FUN = .local_evidence_finite_mean
)
names(member_mae)[2L] <- "local_error_mae"
member_theta_range <- aggregate(
  theta_summary$centered_total_log_surface_error,
  by = list(member = theta_summary$member),
  FUN = .local_evidence_finite_range
)
names(member_theta_range)[2L] <- "centered_total_error_range"
member_summary <- Reduce(function(x, y) merge(x, y, by = "member"),
                         list(member_summary, member_mae, member_theta_range))

dir.create(dirname(results_file), showWarnings = FALSE, recursive = TRUE)
utils::write.csv(rows, rows_csv, row.names = FALSE)
utils::write.csv(theta_summary, theta_csv, row.names = FALSE)
utils::write.csv(local_summary, local_csv, row.names = FALSE)
saveRDS(
  list(
    member_files = member_files,
    theta_ref = theta_ref,
    rows = rows,
    theta_summary = theta_summary,
    local_summary = local_summary,
    focus_regression = focus_regression,
    member_summary = member_summary,
    confirmation = confirm,
    settings = list(
      label = label,
      screen_particles = screen_particles,
      confirm_particles = confirm_particles,
      confirm_reps = confirm_reps,
      confirm_top_pairs = confirm_top_pairs,
      theta_probs = theta_probs,
      theta_max_points = theta_max_points,
      focus_hyper_names = focus_hyper_names,
      cores = cores
    ),
    rows_csv = rows_csv,
    theta_csv = theta_csv,
    local_csv = local_csv
  ),
  results_file
)

cat(sprintf("Saved diagnostic: %s\n", results_file))
cat(sprintf("Saved rows: %s\n", rows_csv))
cat(sprintf("Saved theta summary: %s\n", theta_csv))
cat(sprintf("Saved local summary: %s\n", local_csv))
cat("\nMember summary:\n")
print(member_summary, row.names = FALSE)
cat("\nFocus regression:\n")
print(focus_regression, row.names = FALSE)
cat("\nWorst local contributors:\n")
print(utils::head(local_summary, 12L), row.names = FALSE)
cat("\nTheta surface errors:\n")
print(theta_summary[, c("member", "theta_label", "total_log_surface_error", "centered_total_log_surface_error",
                        "mean_abs_local_error", "n_uncertified_locals", focus_hyper_names), drop = FALSE],
      row.names = FALSE)
