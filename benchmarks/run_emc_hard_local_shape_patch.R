#!/usr/bin/env Rscript

rm(list = ls())

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[1L]))
} else {
  normalizePath("benchmarks/run_emc_hard_local_shape_patch.R")
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

arg_lgl <- function(args, key, default) {
  val <- args[[key]]
  if (is.null(val) || !nzchar(val)) return(isTRUE(default))
  val <- tolower(as.character(val))
  if (val %in% c("1", "true", "t", "yes", "y")) return(TRUE)
  if (val %in% c("0", "false", "f", "no", "n")) return(FALSE)
  stop(sprintf("Argument %s must be true or false.", key))
}

arg_chr_vec <- function(args, key, default) {
  val <- args[[key]]
  if (is.null(val) || !nzchar(val)) return(as.character(default))
  trimws(strsplit(val, ",", fixed = TRUE)[[1L]])
}

arg_int_vec <- function(args, key, default) {
  val <- args[[key]]
  if (is.null(val) || !nzchar(val)) return(as.integer(default))
  x <- trimws(strsplit(val, ",", fixed = TRUE)[[1L]])
  if (length(x) == 1L && identical(tolower(x), "all")) return(integer())
  as.integer(x)
}

arg_num_vec <- function(args, key, default) {
  val <- args[[key]]
  if (is.null(val) || !nzchar(val)) return(as.numeric(default))
  as.numeric(trimws(strsplit(val, ",", fixed = TRUE)[[1L]]))
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

finite_mae <- function(x) {
  x <- abs(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x)) mean(x) else NA_real_
}

finite_max_abs <- function(x) {
  x <- abs(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x)) max(x) else NA_real_
}

source("smc_core.R")
source("population_models.R")
source("local_charts.R")

cli_args <- parse_cli_args(commandArgs(trailingOnly = TRUE))
detected_cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
if (!is.finite(detected_cores) || detected_cores < 1L) detected_cores <- 1L

label <- arg_chr(cli_args, "label", "emc_hard_local_shape_patch")
seed <- arg_int(cli_args, "seed", 20260525L)
cores <- arg_int(cli_args, "cores", min(4L, detected_cores))
data_file <- arg_chr(cli_args, "data_file", file.path("benchmarks", "samples", "full_EMC2.RData"))
source_results <- arg_chr(
  cli_args,
  "source_results",
  file.path("benchmarks", "results", "emc_local_surface_hypothesis_hard6_both_results.rds")
)
gold_csv <- arg_chr(
  cli_args,
  "gold_csv",
  file.path("benchmarks", "results", "emc_local_surface_hypothesis_hard6_gold6144_rows.csv")
)
local <- arg_chr(cli_args, "local", "464")
member <- arg_chr(cli_args, "member", "post_outer")
probe_theta_rows_arg <- arg_int_vec(cli_args, "probe_theta_rows", integer())
compare_theta_rows_arg <- arg_int_vec(cli_args, "compare_theta_rows", integer())
kernel_hyper_names_arg <- arg_chr_vec(cli_args, "kernel_hyper_names", character())
posterior_comparison_csv <- arg_chr(
  cli_args,
  "posterior_comparison_csv",
  file.path(
    "benchmarks",
    "results",
    "population_emc_chart_atlas_gate_emc_posterior_loose_run1_post_outer_posterior_comparison.csv"
  )
)
posterior_axis_mean_error_weight <- arg_num(cli_args, "posterior_axis_mean_error_weight", 1)
posterior_axis_wasserstein_weight <- arg_num(cli_args, "posterior_axis_wasserstein_weight", 0.25)
posterior_axis_shape_weight <- arg_num(cli_args, "posterior_axis_shape_weight", 0.1)
posterior_axis_min_abs_standardized_mean_error <- arg_num(
  cli_args,
  "posterior_axis_min_abs_standardized_mean_error",
  0
)
probe_particles <- arg_int(cli_args, "probe_particles", 768L)
probe_reps <- arg_int(cli_args, "probe_reps", 4L)
probe_max_steps <- arg_int(cli_args, "probe_max_steps", 176L)
probe_mcmc_moves <- arg_int(cli_args, "probe_mcmc_moves", 2L)
probe_target_cess <- arg_num(cli_args, "probe_target_cess", 0.9)
kernel_scale <- arg_num(cli_args, "kernel_scale", 0.35)
scale_floor <- arg_num(cli_args, "scale_floor", 2)
shrink <- arg_num(cli_args, "shrink", 4)
kernel_scale_grid <- arg_num_vec(cli_args, "kernel_scale_grid", kernel_scale)
scale_floor_grid <- arg_num_vec(cli_args, "scale_floor_grid", scale_floor)
shrink_grid <- arg_num_vec(cli_args, "shrink_grid", shrink)
kernel_selection_score <- arg_chr(cli_args, "kernel_selection_score", "loo")
min_train <- arg_int(cli_args, "min_train", 4L)
probe_sd_floor <- arg_num(cli_args, "probe_sd_floor", 0.05)
reuse_probe_cache <- arg_lgl(cli_args, "reuse_probe_cache", TRUE)
refresh_probe_cache <- arg_lgl(cli_args, "refresh_probe_cache", FALSE)

out_prefix <- file.path("benchmarks", "results", label)
results_file <- arg_chr(cli_args, "results_file", paste0(out_prefix, "_results.rds"))
rows_csv <- arg_chr(cli_args, "rows_csv", paste0(out_prefix, "_rows.csv"))
summary_csv <- arg_chr(cli_args, "summary_csv", paste0(out_prefix, "_summary.csv"))
probe_csv <- arg_chr(cli_args, "probe_csv", paste0(out_prefix, "_probes.csv"))
selection_csv <- arg_chr(cli_args, "selection_csv", paste0(out_prefix, "_selection.csv"))
plot_file <- arg_chr(cli_args, "plot_file", paste0(out_prefix, "_diagnostics.png"))
probe_cache_file <- arg_chr(cli_args, "probe_cache_file", paste0(out_prefix, "_probe_cache.rds"))
for (path in c(results_file, rows_csv, summary_csv, probe_csv, selection_csv, plot_file, probe_cache_file)) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
}

if (!file.exists(data_file)) stop("Missing EMC2 data file: ", data_file)
if (!file.exists(source_results)) stop("Missing source results: ", source_results)
if (!file.exists(gold_csv)) stop("Missing gold comparison CSV: ", gold_csv)
if (!length(kernel_hyper_names_arg) ||
    identical(tolower(kernel_hyper_names_arg[1L]), "auto")) {
  if (!file.exists(posterior_comparison_csv)) {
    stop("Automatic kernel-axis selection requires posterior_comparison_csv: ", posterior_comparison_csv)
  }
}

load(data_file)
if (!exists("ELP_DDM", inherits = FALSE)) stop("The EMC2 data file must define ELP_DDM.")
source_obj <- readRDS(source_results)
gold_rows <- utils::read.csv(gold_csv, check.names = FALSE)
if (is.null(source_obj$theta_design) || is.null(source_obj$rows)) {
  stop("source_results must contain theta_design and rows.")
}

emc <- ELP_DDM[[1L]]
data_list <- emc$data
model_factory <- emc$model
alpha_names <- emc$par_names
population_model <- source_obj$theta_design$population_model
theta <- source_obj$theta_design$theta
theta_metadata <- source_obj$theta_design$metadata
if (!local %in% names(data_list)) {
  stop("Selected local is not present in the EMC data: ", local)
}
if (!member %in% unique(source_obj$rows$member)) {
  stop("Selected member is not present in source rows: ", member)
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

base_rows <- source_obj$rows[
  source_obj$rows$member == member & source_obj$rows$local == local,
  ,
  drop = FALSE
]
base_rows <- base_rows[order(base_rows$theta_row), , drop = FALSE]
base_log_m <- rep(NA_real_, nrow(theta))
base_log_m[base_rows$theta_row] <- base_rows$atlas_log_marginal

gold_local <- unique(gold_rows[gold_rows$local == local, c(
  "local", "theta_row", "theta_label", "gold_log_marginal",
  "gold_replicate_sd", "gold_replicate_range"
), drop = FALSE])
if (length(compare_theta_rows_arg)) {
  gold_local <- gold_local[gold_local$theta_row %in% compare_theta_rows_arg, , drop = FALSE]
}
if (!nrow(gold_local)) {
  stop("No gold rows remain for local ", local, ".")
}

probe_theta_rows <- if (length(probe_theta_rows_arg)) {
  sort(unique(probe_theta_rows_arg))
} else {
  seq_len(nrow(theta))
}
probe_theta_rows <- probe_theta_rows[
  probe_theta_rows >= 1L &
    probe_theta_rows <= nrow(theta) &
    is.finite(base_log_m[probe_theta_rows])
]
if (!length(probe_theta_rows)) {
  stop("No valid probe theta rows.")
}

auto_kernel <- !length(kernel_hyper_names_arg) ||
  identical(tolower(kernel_hyper_names_arg[1L]), "auto")
posterior_comparison <- NULL
posterior_axis <- NULL
axis_selection_table <- NULL
kernel_hyper_names <- if (isTRUE(auto_kernel)) {
  posterior_comparison <- utils::read.csv(posterior_comparison_csv, check.names = FALSE)
  posterior_axis <- local_evidence_detect_kernel_axis(
    posterior_comparison = posterior_comparison,
    theta_names = population_model$hyper_names,
    mean_error_weight = posterior_axis_mean_error_weight,
    wasserstein_weight = posterior_axis_wasserstein_weight,
    shape_weight = posterior_axis_shape_weight,
    min_abs_standardized_mean_error = posterior_axis_min_abs_standardized_mean_error
  )
  axis_selection_table <- posterior_axis$table
  posterior_axis$active_hyper_names
} else {
  intersect(kernel_hyper_names_arg, population_model$hyper_names)
}
if (!length(kernel_hyper_names)) stop("No valid kernel hyperparameters were selected.")

probe_settings <- list(
  data_file = normalizePath(data_file),
  source_results = normalizePath(source_results),
  local = local,
  theta = theta,
  probe_theta_rows = probe_theta_rows,
  probe_particles = as.integer(probe_particles),
  probe_reps = as.integer(probe_reps),
  probe_max_steps = as.integer(probe_max_steps),
  probe_mcmc_moves = as.integer(probe_mcmc_moves),
  probe_target_cess = as.numeric(probe_target_cess),
  seed = as.integer(seed)
)
probe_cache_compatible <- function(cache) {
  if (!is.list(cache) || is.null(cache$probes) || is.null(cache$settings)) return(FALSE)
  s <- cache$settings
  identical(s$data_file, probe_settings$data_file) &&
    identical(s$source_results, probe_settings$source_results) &&
    identical(as.character(s$local), as.character(probe_settings$local)) &&
    identical(as.integer(s$probe_theta_rows), as.integer(probe_settings$probe_theta_rows)) &&
    identical(as.integer(s$probe_particles), as.integer(probe_settings$probe_particles)) &&
    identical(as.integer(s$probe_reps), as.integer(probe_settings$probe_reps)) &&
    identical(as.integer(s$probe_max_steps), as.integer(probe_settings$probe_max_steps)) &&
    identical(as.integer(s$probe_mcmc_moves), as.integer(probe_settings$probe_mcmc_moves)) &&
    isTRUE(all.equal(as.numeric(s$probe_target_cess), probe_settings$probe_target_cess, tolerance = 1e-12)) &&
    isTRUE(all.equal(as.matrix(s$theta), as.matrix(probe_settings$theta), tolerance = 1e-10, check.attributes = FALSE))
}

run_probe <- function(theta_row, rep_id) {
  run <- .local_chart_run_smc(
    local_id = local,
    theta_anchor = theta[theta_row, , drop = FALSE],
    data_i = data_list[[local]],
    loglik_fn = loglik_emc2,
    population_model = population_model,
    M = probe_particles,
    target_cess = probe_target_cess,
    resample_threshold = 0.5,
    n_mcmc_moves = probe_mcmc_moves,
    rw_scale = 0.75,
    G_mix = 8L,
    da_enable = TRUE,
    refit_every = 2L,
    max_steps = probe_max_steps,
    deterministic_resampling = FALSE,
    n_cores = 1L,
    seed = seed + 1000003L * match(local, names(data_list)) +
      9176L * theta_row + 104729L * rep_id + 8191L * probe_particles,
    verbose = FALSE,
    source = "emc_hard_local_shape_patch_probe"
  )
  diag <- run$diagnostics %||% list()
  data.frame(
    local = local,
    theta_row = theta_row,
    theta_label = theta_metadata$theta_label[theta_row],
    replicate = rep_id,
    log_marginal = run$logZ,
    path_se = run$logZ_se,
    final_ess_frac = as.numeric(diag$final_ess_frac %||% NA_real_),
    min_path_ess_frac = as.numeric(diag$min_path_ess_frac %||% NA_real_),
    check.names = FALSE
  )
}

cat(sprintf(
  "EMC hard-local shape patch: local=%s member=%s | %d probe theta x %d reps x %d particles\n",
  local, member, length(probe_theta_rows), probe_reps, probe_particles
))
cat("Probe theta rows:", paste(probe_theta_rows, collapse = ", "), "\n")
cat("Gold theta rows:", paste(gold_local$theta_row, collapse = ", "), "\n")
if (isTRUE(auto_kernel)) {
  cat("Posterior comparison:", posterior_comparison_csv, "\n")
  cat("Kernel axis selected from posterior diagnostics:", paste(kernel_hyper_names, collapse = "+"), "\n")
  print(axis_selection_table[seq_len(min(6L, nrow(axis_selection_table))), ], row.names = FALSE)
} else {
  cat("Kernel hyperparameters:", paste(kernel_hyper_names, collapse = ", "), "\n")
}

if (isTRUE(reuse_probe_cache) && !isTRUE(refresh_probe_cache) && file.exists(probe_cache_file)) {
  cache <- readRDS(probe_cache_file)
  if (!probe_cache_compatible(cache)) {
    stop("Probe cache is incompatible. Use --refresh_probe_cache=true or a different --probe_cache_file.")
  }
  probes <- cache$probes
  cat("Loaded probe cache:", probe_cache_file, "\n")
} else {
  probe_jobs <- expand.grid(
    theta_row = probe_theta_rows,
    replicate = seq_len(probe_reps),
    KEEP.OUT.ATTRS = FALSE
  )
  probe_parts <- if (cores <= 1L || nrow(probe_jobs) <= 1L) {
    lapply(seq_len(nrow(probe_jobs)), function(k) {
      run_probe(probe_jobs$theta_row[k], probe_jobs$replicate[k])
    })
  } else {
    parallel::mclapply(seq_len(nrow(probe_jobs)), function(k) {
      run_probe(probe_jobs$theta_row[k], probe_jobs$replicate[k])
    }, mc.cores = min(cores, nrow(probe_jobs)))
  }
  probes <- do.call(rbind, probe_parts)
  saveRDS(list(probes = probes, settings = probe_settings), probe_cache_file)
  cat("Saved probe cache:", probe_cache_file, "\n")
}

probe_summary <- local_evidence_summarize_probe_replicates(
  probes,
  theta_row_col = "theta_row",
  log_m_col = "log_marginal",
  se_col = "path_se"
)
probe_summary$local <- local
probe_summary$theta_label <- theta_metadata$theta_label[probe_summary$theta_row]
probe_summary$probe_mean_path_se <- vapply(split(probes$path_se, probes$theta_row), mean, numeric(1), na.rm = TRUE)
probe_summary <- probe_summary[, c(
  "local", "theta_row", "theta_label", "probe_log_m", "probe_sd",
  "probe_reps", "probe_min", "probe_max", "probe_mean_path_se"
), drop = FALSE]

if (isTRUE(auto_kernel)) {
  patch <- local_evidence_tune_kernel_residual_patch(
    theta = theta,
    base_log_m = base_log_m,
    probe_theta_row = probe_summary$theta_row,
    probe_log_m = probe_summary$probe_log_m,
    probe_sd = probe_summary$probe_sd,
    theta_weights = theta_metadata$theta_weight,
    active_hyper_names = kernel_hyper_names,
    kernel_scale_grid = kernel_scale_grid,
    scale_floor_grid = scale_floor_grid,
    shrink_grid = shrink_grid,
    min_train = min_train,
    probe_sd_floor = probe_sd_floor,
    selection_score = kernel_selection_score,
    probe_replicates = probes,
    replicate_col = "replicate",
    log_m_col = "log_marginal",
    se_col = "path_se"
  )
  selected <- patch$selected
  kernel_hyper_names <- strsplit(selected$active_hyper_names[1L], "+", fixed = TRUE)[[1L]]
  kernel_scale <- selected$kernel_scale[1L]
  scale_floor <- selected$scale_floor[1L]
  shrink <- selected$shrink[1L]
  selection_table <- patch$selection_table
} else {
  patch <- local_evidence_fit_kernel_residual_patch(
    theta = theta,
    base_log_m = base_log_m,
    probe_theta_row = probe_summary$theta_row,
    probe_log_m = probe_summary$probe_log_m,
    probe_sd = probe_summary$probe_sd,
    theta_weights = theta_metadata$theta_weight,
    active_hyper_names = kernel_hyper_names,
    kernel_scale = kernel_scale,
    scale_floor = scale_floor,
    shrink = shrink,
    min_train = min_train,
    probe_sd_floor = probe_sd_floor
  )
  selection_table <- data.frame(
    candidate_id = 1L,
    active_hyper_names = paste(kernel_hyper_names, collapse = "+"),
    dimension = length(kernel_hyper_names),
    kernel_scale = kernel_scale,
    scale_floor = scale_floor,
    shrink = shrink,
    score = NA_real_,
    score_source = "fixed",
    replicate_validation_rmse = NA_real_,
    loo_rmse = patch$diagnostics$loo_rmse[1L],
    training_rmse = patch$diagnostics$training_rmse[1L],
    signal_sd = patch$diagnostics$signal_sd[1L],
    reliability = patch$diagnostics$reliability[1L],
    shape_scale = patch$diagnostics$shape_scale[1L],
    offset = patch$diagnostics$offset[1L],
    status = patch$diagnostics$status[1L],
    check.names = FALSE
  )
}

if (isTRUE(auto_kernel) && !is.null(posterior_axis)) {
  selection_table$axis_source <- "posterior_diagnostic"
  selection_table$axis_family <- posterior_axis$family
  selection_table$posterior_axis_score <- posterior_axis$selected$posterior_axis_score[1L]
  selection_table$posterior_axis_parameters <- posterior_axis$selected$parameters[1L]
} else {
  selection_table$axis_source <- "fixed"
  selection_table$axis_family <- NA_character_
  selection_table$posterior_axis_score <- NA_real_
  selection_table$posterior_axis_parameters <- NA_character_
}

comparison_rows <- merge(
  gold_local,
  data.frame(
    theta_row = seq_len(nrow(theta)),
    base_log_marginal = base_log_m,
    patched_log_marginal = patch$estimate,
    correction = patch$correction,
    check.names = FALSE
  ),
  by = "theta_row",
  all.x = TRUE,
  sort = FALSE
)
comparison_rows <- merge(
  comparison_rows,
  probe_summary[, c("theta_row", "probe_log_m", "probe_sd", "probe_reps"), drop = FALSE],
  by = "theta_row",
  all.x = TRUE,
  sort = FALSE
)
comparison_rows$base_error_vs_gold <- comparison_rows$base_log_marginal - comparison_rows$gold_log_marginal
comparison_rows$patched_error_vs_gold <- comparison_rows$patched_log_marginal - comparison_rows$gold_log_marginal
comparison_rows$probe_error_vs_gold <- comparison_rows$probe_log_m - comparison_rows$gold_log_marginal
comparison_rows$base_centered_error_vs_gold <-
  comparison_rows$base_error_vs_gold - mean(comparison_rows$base_error_vs_gold, na.rm = TRUE)
comparison_rows$patched_centered_error_vs_gold <-
  comparison_rows$patched_error_vs_gold - mean(comparison_rows$patched_error_vs_gold, na.rm = TRUE)
comparison_rows$probe_centered_error_vs_gold <-
  comparison_rows$probe_error_vs_gold - mean(comparison_rows$probe_error_vs_gold, na.rm = TRUE)

metric_row <- function(method, error, centered_error) {
  data.frame(
    method = method,
    n_gold_points = sum(is.finite(error)),
    rmse_vs_gold = finite_rmse(error),
    centered_rmse_vs_gold = finite_rmse(centered_error),
    mae_vs_gold = finite_mae(error),
    max_abs_centered_error_vs_gold = finite_max_abs(centered_error),
    check.names = FALSE
  )
}
summary <- rbind(
  metric_row("base_atlas", comparison_rows$base_error_vs_gold, comparison_rows$base_centered_error_vs_gold),
  metric_row("shape_patch", comparison_rows$patched_error_vs_gold, comparison_rows$patched_centered_error_vs_gold),
  metric_row("probe_only", comparison_rows$probe_error_vs_gold, comparison_rows$probe_centered_error_vs_gold)
)
diag_cols <- patch$diagnostics
summary_meta <- data.frame(
  local = local,
  member = member,
  probe_particles = probe_particles,
  probe_reps = probe_reps,
  kernel_scale = kernel_scale,
  scale_floor = scale_floor,
  shrink = shrink,
  kernel_hyper_names = paste(kernel_hyper_names, collapse = "+"),
  kernel_auto = isTRUE(auto_kernel),
  kernel_axis_source = if (isTRUE(auto_kernel)) "posterior_diagnostic" else "fixed",
  kernel_axis_family = if (!is.null(posterior_axis)) posterior_axis$family else NA_character_,
  posterior_axis_score = if (!is.null(posterior_axis)) posterior_axis$selected$posterior_axis_score[1L] else NA_real_,
  kernel_selection_score = kernel_selection_score,
  patch_status = diag_cols$status[1L],
  patch_shape_scale = diag_cols$shape_scale[1L],
  patch_loo_rmse = diag_cols$loo_rmse[1L],
  check.names = FALSE
)
summary <- cbind(summary, summary_meta[rep(1L, nrow(summary)), , drop = FALSE])

utils::write.csv(comparison_rows, rows_csv, row.names = FALSE)
utils::write.csv(summary, summary_csv, row.names = FALSE)
utils::write.csv(probes, probe_csv, row.names = FALSE)
utils::write.csv(selection_table, selection_csv, row.names = FALSE)

grDevices::png(plot_file, width = 1800, height = 1200, res = 150, bg = "white")
old_par <- graphics::par(mfrow = c(2, 2), mar = c(5, 4, 3, 1), oma = c(0, 0, 1.2, 0))
ylim <- range(
  comparison_rows$gold_log_marginal,
  comparison_rows$base_log_marginal,
  comparison_rows$patched_log_marginal,
  comparison_rows$probe_log_m,
  na.rm = TRUE
)
graphics::plot(
  comparison_rows$theta_row,
  comparison_rows$gold_log_marginal,
  type = "b",
  pch = 18,
  lwd = 2,
  col = "black",
  ylim = ylim,
  xlab = "theta row",
  ylab = "log m_i(theta)",
  main = "Gold vs atlas/patch/probe"
)
graphics::arrows(
  comparison_rows$theta_row,
  comparison_rows$gold_log_marginal - comparison_rows$gold_replicate_sd,
  comparison_rows$theta_row,
  comparison_rows$gold_log_marginal + comparison_rows$gold_replicate_sd,
  angle = 90,
  code = 3,
  length = 0.04,
  col = "grey35"
)
graphics::lines(comparison_rows$theta_row, comparison_rows$base_log_marginal, type = "b", pch = 19, col = "firebrick3", lwd = 2)
graphics::lines(comparison_rows$theta_row, comparison_rows$patched_log_marginal, type = "b", pch = 19, col = "steelblue4", lwd = 2)
graphics::lines(comparison_rows$theta_row, comparison_rows$probe_log_m, type = "b", pch = 17, col = "darkorange3", lwd = 2)
graphics::legend(
  "bottomright",
  legend = c("6144 nested gold +/- SD", "base atlas", "kernel patch", "probe logmeanexp"),
  col = c("black", "firebrick3", "steelblue4", "darkorange3"),
  pch = c(18, 19, 19, 17),
  lwd = 2,
  bty = "n",
  cex = 0.8
)

err_ylim <- range(
  comparison_rows$base_centered_error_vs_gold,
  comparison_rows$patched_centered_error_vs_gold,
  comparison_rows$probe_centered_error_vs_gold,
  na.rm = TRUE
)
graphics::plot(
  comparison_rows$theta_row,
  comparison_rows$base_centered_error_vs_gold,
  type = "b",
  pch = 19,
  col = "firebrick3",
  ylim = err_ylim,
  xlab = "theta row",
  ylab = "centered error vs gold",
  main = "Shape error"
)
graphics::abline(h = 0, col = "grey35")
graphics::lines(comparison_rows$theta_row, comparison_rows$patched_centered_error_vs_gold, type = "b", pch = 19, col = "steelblue4", lwd = 2)
graphics::lines(comparison_rows$theta_row, comparison_rows$probe_centered_error_vs_gold, type = "b", pch = 17, col = "darkorange3", lwd = 2)
graphics::legend("topright", legend = c("base", "patch", "probe"), col = c("firebrick3", "steelblue4", "darkorange3"), pch = c(19, 19, 17), lwd = 2, bty = "n")

train <- patch$train
graphics::plot(
  train$theta_row,
  train$delta,
  type = "b",
  pch = 19,
  col = "grey25",
  xlab = "probe theta row",
  ylab = "probe - atlas",
  main = "Observed residuals"
)
graphics::abline(h = patch$diagnostics$offset[1L], col = "steelblue4", lty = 2, lwd = 2)

summary_plot <- summary
graphics::barplot(
  summary_plot$centered_rmse_vs_gold,
  names.arg = summary_plot$method,
  col = c("firebrick3", "steelblue4", "darkorange3"),
  border = NA,
  las = 2,
  ylab = "centered RMSE vs gold",
  main = "Hard-local metric"
)
graphics::par(old_par)
graphics::mtext(sprintf("EMC local %s shape patch against stored 6144 nested SMC gold", local), outer = TRUE, font = 2)
grDevices::dev.off()

saveRDS(
  list(
    comparison_rows = comparison_rows,
    summary = summary,
    probes = probes,
    probe_summary = probe_summary,
    patch = patch,
    posterior_axis = posterior_axis,
    posterior_axis_selection_table = axis_selection_table,
    source_results = source_results,
    gold_csv = gold_csv,
    settings = list(
      label = label,
      local = local,
      member = member,
      probe_theta_rows = probe_theta_rows,
      compare_theta_rows = gold_local$theta_row,
      probe_particles = probe_particles,
      probe_reps = probe_reps,
      probe_max_steps = probe_max_steps,
      kernel_hyper_names = kernel_hyper_names,
      kernel_auto = isTRUE(auto_kernel),
      posterior_comparison_csv = posterior_comparison_csv,
      posterior_axis_mean_error_weight = posterior_axis_mean_error_weight,
      posterior_axis_wasserstein_weight = posterior_axis_wasserstein_weight,
      posterior_axis_shape_weight = posterior_axis_shape_weight,
      posterior_axis_min_abs_standardized_mean_error =
        posterior_axis_min_abs_standardized_mean_error,
      kernel_scale = kernel_scale,
      scale_floor = scale_floor,
      shrink = shrink,
      kernel_scale_grid = kernel_scale_grid,
      scale_floor_grid = scale_floor_grid,
      shrink_grid = shrink_grid,
      kernel_selection_score = kernel_selection_score,
      min_train = min_train,
      probe_sd_floor = probe_sd_floor
    )
  ),
  results_file
)

cat("Summary:\n")
print(summary, row.names = FALSE)
cat("Saved rows:", rows_csv, "\n")
cat("Saved summary:", summary_csv, "\n")
cat("Saved probes:", probe_csv, "\n")
cat("Saved selection:", selection_csv, "\n")
cat("Saved plot:", plot_file, "\n")
