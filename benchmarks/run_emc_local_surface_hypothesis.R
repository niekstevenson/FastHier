#!/usr/bin/env Rscript

rm(list = ls())

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[1L]))
} else {
  normalizePath("benchmarks/run_emc_local_surface_hypothesis.R")
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

arg_num_vec <- function(args, key, default) {
  val <- args[[key]]
  if (is.null(val) || !nzchar(val)) return(as.numeric(default))
  as.numeric(strsplit(val, ",", fixed = TRUE)[[1L]])
}

finite_rmse <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x)) sqrt(mean(x^2)) else NA_real_
}

finite_mae <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x)) mean(abs(x)) else NA_real_
}

finite_max_abs <- function(x) {
  x <- abs(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x)) max(x) else NA_real_
}

finite_mean <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x)) mean(x) else NA_real_
}

source("smc_core.R")
source("population_models.R")
source("local_charts.R")

cli_args <- parse_cli_args(commandArgs(trailingOnly = TRUE))
detected_cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
if (!is.finite(detected_cores) || detected_cores < 1L) detected_cores <- 1L

label <- arg_chr(cli_args, "label", "emc_local_surface_hypothesis")
seed <- arg_int(cli_args, "seed", 20260525L)
cores <- arg_int(cli_args, "cores", 1L)
data_file <- arg_chr(cli_args, "data_file", file.path("benchmarks", "samples", "full_EMC2.RData"))
post_outer_checkpoint <- arg_chr(
  cli_args,
  "post_outer_checkpoint",
  file.path("benchmarks", "results", "population_emc_chart_atlas_gate_emc_posterior_loose_run1_checkpoint_post_outer.rds")
)
post_cal_checkpoint <- arg_chr(
  cli_args,
  "post_cal_checkpoint",
  file.path("benchmarks", "results", "population_emc_chart_atlas_gate_emc_posterior_loose_run1_checkpoint_posterior_calibration_round_1.rds")
)
local_ids_arg <- arg_chr_vec(cli_args, "locals", c("464", "650", "758"))
theta_max_points <- arg_int(cli_args, "theta_max_points", 7L)
theta_probs <- arg_num_vec(cli_args, "theta_probs", c(0.05, 0.25, 0.50, 0.75, 0.95))
theta_sources <- arg_chr_vec(cli_args, "theta_sources", c("emc", "outer"))
focus_parameters <- arg_chr_vec(
  cli_args,
  "focus_parameters",
  c("mu_sv", "sigma2_sv", "sigma2_v_LogFreq", "sigma2_v")
)
ref_particles <- arg_int(cli_args, "ref_particles", 384L)
ref_reps <- arg_int(cli_args, "ref_reps", 2L)
ref_max_steps <- arg_int(cli_args, "ref_max_steps", 128L)
ref_mcmc_moves <- arg_int(cli_args, "ref_mcmc_moves", 2L)
ref_target_cess <- arg_num(cli_args, "ref_target_cess", 0.9)
repair_particles <- arg_int(cli_args, "repair_particles", ref_particles)
repair_max_updates_per_local <- arg_int(cli_args, "repair_max_updates_per_local", 2L)
repair_confirmation_reps <- arg_int(cli_args, "repair_confirmation_reps", 1L)
repair_max_steps <- arg_int(cli_args, "repair_max_steps", ref_max_steps)
repair_mcmc_moves <- arg_int(cli_args, "repair_mcmc_moves", ref_mcmc_moves)
min_particle_mis_ess <- arg_num(cli_args, "min_particle_mis_ess", 0.05)
max_particle_mis_psis_k <- arg_num(cli_args, "max_particle_mis_psis_k", 0.7)
use_post_cal <- arg_lgl(cli_args, "use_post_cal", TRUE)
reuse_reference_cache <- arg_lgl(cli_args, "reuse_reference_cache", TRUE)
refresh_reference_cache <- arg_lgl(cli_args, "refresh_reference_cache", FALSE)

out_prefix <- file.path("benchmarks", "results", label)
results_file <- arg_chr(cli_args, "results_file", paste0(out_prefix, "_results.rds"))
rows_csv <- arg_chr(cli_args, "rows_csv", paste0(out_prefix, "_rows.csv"))
summary_csv <- arg_chr(cli_args, "summary_csv", paste0(out_prefix, "_summary.csv"))
repair_csv <- arg_chr(cli_args, "repair_csv", paste0(out_prefix, "_repair_rows.csv"))
plot_file <- arg_chr(cli_args, "plot_file", paste0(out_prefix, "_diagnostics.png"))
surface_plot_file <- arg_chr(cli_args, "surface_plot_file", paste0(out_prefix, "_surface_errors.png"))
reference_cache_file <- arg_chr(cli_args, "reference_cache_file", paste0(out_prefix, "_reference_cache.rds"))

for (path in c(results_file, rows_csv, summary_csv, repair_csv, plot_file, surface_plot_file, reference_cache_file)) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
}

if (!file.exists(data_file)) stop("Missing EMC2 data file: ", data_file)
if (!file.exists(post_outer_checkpoint)) stop("Missing post-outer checkpoint: ", post_outer_checkpoint)
if (isTRUE(use_post_cal) && !file.exists(post_cal_checkpoint)) {
  stop("Missing post-calibration checkpoint: ", post_cal_checkpoint)
}

load(data_file)
if (!exists("ELP_DDM", inherits = FALSE)) stop("The EMC2 data file must define ELP_DDM.")
post_outer <- readRDS(post_outer_checkpoint)
post_cal <- if (isTRUE(use_post_cal)) readRDS(post_cal_checkpoint) else NULL

emc <- ELP_DDM[[1L]]
data_list_all <- emc$data
model_factory <- emc$model
alpha_names <- emc$par_names
model <- post_outer$factor_set$population_model

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
colnames(emc_mu) <- paste0("mu_", model$alpha_names)
colnames(emc_sigma2) <- paste0("sigma2_", model$alpha_names)
emc_draws <- data.frame(emc_mu, emc_sigma2, check.names = FALSE)
theta_emc <- local_atlas_theta_from_draws(emc_draws, model)

focus_hyper_names <- unique(c(
  focus_parameters[startsWith(focus_parameters, "mu_")],
  sub("^sigma2_", "log_sigma2_", focus_parameters[startsWith(focus_parameters, "sigma2_")])
))
focus_hyper_names <- intersect(focus_hyper_names, model$hyper_names)
if (!length(focus_hyper_names)) focus_hyper_names <- model$hyper_names

theta_sources <- intersect(tolower(theta_sources), c("emc", "outer"))
if (!length(theta_sources)) stop("--theta_sources must include emc and/or outer.")
build_theta_part <- function(source_name) {
  if (identical(source_name, "outer")) {
    theta <- post_outer$fit$theta
    weights <- post_outer$fit$w
  } else {
    theta <- theta_emc
    weights <- rep(1 / nrow(theta_emc), nrow(theta_emc))
  }
  build_local_evidence_calibration_design(
    theta = theta,
    population_model = model,
    weights = weights,
    focus_hyper_names = focus_hyper_names,
    probs = theta_probs,
    max_points = theta_max_points,
    n_leverage_points = 1L,
    n_uncertainty_points = 0L,
    n_disagreement_points = 0L,
    source = paste0(source_name, "_reference"),
    label_prefix = source_name
  )
}
theta_parts <- lapply(theta_sources, build_theta_part)
theta_all <- do.call(rbind, lapply(theta_parts, `[[`, "theta"))
metadata_all <- do.call(rbind, lapply(theta_parts, `[[`, "metadata"))
theta_key <- apply(round(theta_all, 10), 1L, paste, collapse = "|")
keep_theta <- !duplicated(theta_key)
theta_all <- theta_all[keep_theta, , drop = FALSE]
metadata_all <- metadata_all[keep_theta, , drop = FALSE]
metadata_all$theta_row <- seq_len(nrow(theta_all))
metadata_all$theta_id <- sprintf("theta_%03d", metadata_all$theta_row)
metadata_all$theta_label <- paste0(metadata_all$theta_source, "_", sprintf("%03d", metadata_all$theta_row))
metadata_all$theta_weight <- .local_chart_normalize_weights(metadata_all$theta_weight, nrow(metadata_all))
for (name in model$hyper_names) {
  metadata_all[[name]] <- theta_all[, name]
}
theta_design <- structure(
  list(
    theta = theta_all,
    metadata = metadata_all,
    population_model = model,
    settings = list(theta_sources = theta_sources, theta_max_points = theta_max_points)
  ),
  class = "local_evidence_theta_design"
)

available_locals <- names(post_outer$factor_set$atlases)
local_ids <- intersect(local_ids_arg, available_locals)
if (!length(local_ids)) {
  stop("None of --locals are present in the checkpoint factor set.")
}
missing_locals <- setdiff(local_ids_arg, local_ids)
if (length(missing_locals)) {
  warning("Skipping locals not present in checkpoint: ", paste(missing_locals, collapse = ", "))
}
if (!all(local_ids %in% names(data_list_all))) {
  stop("EMC data is missing selected local ids.")
}
data_list <- data_list_all[local_ids]

strict_factor_set <- function(fs, local_ids) {
  fs <- build_local_atlas_factor_set(
    atlases = fs$atlases[local_ids],
    population_model = fs$population_model,
    max_chart_distance = Inf,
    min_covering_charts = 3L,
    max_prediction_range = Inf,
    distance_scale = Inf,
    se_floor = 1e-6,
    use_particle_mis = TRUE,
    require_particle_mis = TRUE,
    min_particle_mis_ess = min_particle_mis_ess,
    max_particle_mis_psis_k = max_particle_mis_psis_k,
    max_quadratic_particle_gap = Inf,
    sparse_chart_min_covering = 3L,
    sparse_chart_max_distance = Inf,
    max_leave_chart_out_gap = Inf,
    distance_metric = "fisher",
    surface_method = "derivative_ls",
    particle_mis_role = "estimator",
    particle_mis_batch = TRUE,
    stop_on_uncertified = FALSE,
    use_uncertified_estimates = TRUE
  )
  fs
}

factor_sets <- list(post_outer = strict_factor_set(post_outer$factor_set, local_ids))
if (isTRUE(use_post_cal)) {
  factor_sets$post_cal <- strict_factor_set(post_cal$factor_set, local_ids)
}

reference_settings <- list(
  data_file = normalizePath(data_file),
  post_outer_checkpoint = normalizePath(post_outer_checkpoint),
  locals = local_ids,
  theta = theta_design$theta,
  ref_particles = as.integer(ref_particles),
  ref_reps = as.integer(ref_reps),
  ref_max_steps = as.integer(ref_max_steps),
  ref_mcmc_moves = as.integer(ref_mcmc_moves),
  ref_target_cess = as.numeric(ref_target_cess),
  seed = as.integer(seed)
)

reference_cache_compatible <- function(cache) {
  if (!is.list(cache) || is.null(cache$audit) || is.null(cache$settings)) return(FALSE)
  s <- cache$settings
  identical(as.character(s$locals), as.character(reference_settings$locals)) &&
    identical(as.integer(s$ref_particles), as.integer(reference_settings$ref_particles)) &&
    identical(as.integer(s$ref_reps), as.integer(reference_settings$ref_reps)) &&
    identical(as.integer(s$ref_max_steps), as.integer(reference_settings$ref_max_steps)) &&
    identical(as.integer(s$ref_mcmc_moves), as.integer(reference_settings$ref_mcmc_moves)) &&
    isTRUE(all.equal(as.numeric(s$ref_target_cess), as.numeric(reference_settings$ref_target_cess), tolerance = 1e-12)) &&
    isTRUE(all.equal(as.matrix(s$theta), as.matrix(reference_settings$theta), tolerance = 1e-10, check.attributes = FALSE))
}

cat(sprintf(
  "Local-surface hypothesis benchmark: %d locals x %d theta | reference %d reps x %d particles\n",
  length(local_ids), nrow(theta_design$theta), ref_reps, ref_particles
))
cat("Locals:", paste(local_ids, collapse = ", "), "\n")
cat("Focus hyperparameters:", paste(focus_hyper_names, collapse = ", "), "\n")
cat("Theta sources:", paste(theta_sources, collapse = ", "), "\n")

start_time <- Sys.time()
if (isTRUE(reuse_reference_cache) && !isTRUE(refresh_reference_cache) && file.exists(reference_cache_file)) {
  cache <- readRDS(reference_cache_file)
  if (!reference_cache_compatible(cache)) {
    stop("Reference cache is incompatible. Use --refresh_reference_cache=true or another --reference_cache_file.")
  }
  audit <- cache$audit
  cat("Loaded reference cache:", reference_cache_file, "\n")
} else {
  audit <- build_local_evidence_audit(
    factor_sets = factor_sets,
    theta = theta_design,
    data_list = data_list,
    loglik_fn = loglik_emc2,
    local_ids = seq_along(local_ids),
    n_replicates = ref_reps,
    M = ref_particles,
    local_control = list(
      target_cess = ref_target_cess,
      n_mcmc_moves = ref_mcmc_moves,
      max_steps = ref_max_steps
    ),
    reference_source = "emc_local_surface_hypothesis",
    n_cores = cores,
    seed = seed + 1000L,
    verbose = FALSE
  )
  saveRDS(list(audit = audit, settings = reference_settings), reference_cache_file)
  cat("Saved reference cache:", reference_cache_file, "\n")
}

reference_summary <- .local_evidence_replicate_summary(audit$replicates)

post_outer_rows <- audit$rows[audit$rows$member == "post_outer", , drop = FALSE]
post_outer_rows$repair_score <- (
  ifelse(post_outer_rows$atlas_status != "certified", 10, 0) +
    pmax(post_outer_rows$particle_mis_psis_k - max_particle_mis_psis_k, 0, na.rm = TRUE) +
    pmax(min_particle_mis_ess - post_outer_rows$particle_mis_ess_frac, 0, na.rm = TRUE)
) * sqrt(pmax(post_outer_rows$theta_weight, .Machine$double.eps))
eligible <- post_outer_rows[
  post_outer_rows$atlas_status != "certified" |
    (is.finite(post_outer_rows$particle_mis_psis_k) & post_outer_rows$particle_mis_psis_k > max_particle_mis_psis_k) |
    (is.finite(post_outer_rows$particle_mis_ess_frac) & post_outer_rows$particle_mis_ess_frac < min_particle_mis_ess),
  ,
  drop = FALSE
]
eligible <- eligible[order(eligible$local_pos, -eligible$repair_score), , drop = FALSE]
candidate_pairs <- do.call(rbind, lapply(split(eligible, eligible$local_pos), function(df) {
  df[seq_len(min(nrow(df), repair_max_updates_per_local)), c("local_pos", "theta_row"), drop = FALSE]
}))
if (is.null(candidate_pairs) || !nrow(candidate_pairs)) {
  candidate_pairs <- post_outer_rows[order(-abs(post_outer_rows$error_fresh_minus_atlas)), c("local_pos", "theta_row"), drop = FALSE]
  candidate_pairs <- candidate_pairs[seq_len(min(nrow(candidate_pairs), repair_max_updates_per_local * length(local_ids))), , drop = FALSE]
}
candidate_pairs <- unique(candidate_pairs)

repair <- local_atlas_repair_certification_pairs(
  factor_set = factor_sets$post_outer,
  theta = theta_design$theta,
  data_list = data_list,
  loglik_fn = loglik_emc2,
  theta_weights = theta_design$metadata$theta_weight,
  candidate_pairs = candidate_pairs,
  local_ids = seq_along(local_ids),
  M = repair_particles,
  target_cess = ref_target_cess,
  n_mcmc_moves = repair_mcmc_moves,
  max_steps = repair_max_steps,
  max_updates = nrow(candidate_pairs),
  direct_confirmation_reps = repair_confirmation_reps,
  direct_confirmation_M = repair_particles,
  direct_confirmation_max_sd = 1.5,
  replicate_bootstrap_B = 100L,
  max_direct_graph_z = 3,
  max_direct_graph_chart_shift = 0.35,
  max_direct_graph_existing_shift = 0.15,
  local_control = list(
    candidate_M = repair_particles,
    bridge_particles = repair_particles,
    target_cess = ref_target_cess,
    n_mcmc_moves = repair_mcmc_moves,
    max_steps = repair_max_steps
  ),
  edge_control = list(
    max_intermediates = 4L,
    min_overlap_ess = 0.05,
    max_pareto_k = 0.9,
    max_forward_reverse_gap = 2,
    max_taylor_gap = Inf,
    max_cycle_z = 4,
    distance_metric = "fisher"
  ),
  seed = seed + 900000L,
  verbose = FALSE
)
factor_sets$targeted_repair <- repair$factor_set

evaluate_against_reference <- function(factor_sets, audit, reference_summary) {
  atlas_rows <- .local_evidence_atlas_rows(
    factor_sets = factor_sets,
    theta = audit$theta,
    local_metadata = audit$local_metadata,
    n_cores = cores
  )
  rows <- merge(atlas_rows, audit$local_metadata, by = "local_pos", all.x = TRUE, sort = FALSE)
  rows <- merge(rows, audit$theta_metadata, by = "theta_row", all.x = TRUE, sort = FALSE)
  rows <- merge(
    rows,
    reference_summary,
    by = c("local_pos", "theta_row", "local"),
    all.x = TRUE,
    sort = FALSE
  )
  rows$error_atlas_minus_reference <- rows$atlas_log_marginal - rows$fresh_log_m_center
  rows$centered_error_local <- ave(
    rows$error_atlas_minus_reference,
    interaction(rows$member, rows$local_pos, drop = TRUE),
    FUN = function(x) x - finite_mean(x)
  )
  rows$abs_centered_error_local <- abs(rows$centered_error_local)
  rows$certified <- rows$atlas_status == "certified"
  rows[order(rows$member, rows$local_pos, rows$theta_row), , drop = FALSE]
}

rows <- evaluate_against_reference(factor_sets, audit, reference_summary)
theta_totals <- do.call(rbind, lapply(split(rows, list(rows$member, rows$theta_row), drop = TRUE), function(df) {
  err <- df$error_atlas_minus_reference[is.finite(df$error_atlas_minus_reference)]
  data.frame(
    member = df$member[1L],
    theta_row = df$theta_row[1L],
    theta_label = df$theta_label[1L],
    total_error = if (length(err)) sum(err) else NA_real_,
    mean_abs_local_error = if (length(err)) mean(abs(err)) else NA_real_,
    n_uncertified = sum(df$atlas_status != "certified", na.rm = TRUE),
    check.names = FALSE
  )
}))
theta_totals$centered_total_error <- ave(
  theta_totals$total_error,
  theta_totals$member,
  FUN = function(x) x - finite_mean(x)
)

method_summary <- do.call(rbind, lapply(split(rows, rows$member), function(df) {
  theta_df <- theta_totals[theta_totals$member == df$member[1L], , drop = FALSE]
  data.frame(
    method = df$member[1L],
    n_rows = nrow(df),
    certified_fraction = mean(df$certified, na.rm = TRUE),
    finite_fraction = mean(is.finite(df$error_atlas_minus_reference)),
    raw_rmse = finite_rmse(df$error_atlas_minus_reference),
    local_centered_rmse = finite_rmse(df$centered_error_local),
    local_mae = finite_mae(df$error_atlas_minus_reference),
    max_abs_local_centered_error = finite_max_abs(df$centered_error_local),
    total_centered_rmse = finite_rmse(theta_df$centered_total_error),
    max_abs_total_centered_error = finite_max_abs(theta_df$centered_total_error),
    median_particle_mis_ess_frac = stats::median(df$particle_mis_ess_frac, na.rm = TRUE),
    median_particle_mis_psis_k = stats::median(df$particle_mis_psis_k, na.rm = TRUE),
    check.names = FALSE
  )
}))
method_order <- c("post_outer", "post_cal", "targeted_repair")
method_summary <- method_summary[order(match(method_summary$method, method_order)), , drop = FALSE]

utils::write.csv(rows, rows_csv, row.names = FALSE)
utils::write.csv(method_summary, summary_csv, row.names = FALSE)
utils::write.csv(repair$probes, repair_csv, row.names = FALSE)

cols <- c(post_outer = "firebrick3", post_cal = "grey35", targeted_repair = "steelblue4")
present_methods <- intersect(method_order, method_summary$method)

grDevices::png(plot_file, width = 2200, height = 1500, res = 150, bg = "white")
old_par <- graphics::par(mfrow = c(2, 3), mar = c(7, 4, 3, 1), oma = c(0, 0, 1.7, 0))
barplot(
  method_summary$local_centered_rmse[match(present_methods, method_summary$method)],
  names.arg = present_methods,
  col = cols[present_methods],
  border = NA,
  las = 2,
  ylab = "local centered RMSE",
  main = "Shape error per local"
)
barplot(
  method_summary$total_centered_rmse[match(present_methods, method_summary$method)],
  names.arg = present_methods,
  col = cols[present_methods],
  border = NA,
  las = 2,
  ylab = "theta-total centered RMSE",
  main = "Product-surface error"
)
barplot(
  method_summary$certified_fraction[match(present_methods, method_summary$method)],
  names.arg = present_methods,
  col = cols[present_methods],
  border = NA,
  las = 2,
  ylim = c(0, 1),
  ylab = "certified fraction",
  main = "Strict query certification"
)
base_rows <- rows[rows$member == "post_outer", , drop = FALSE]
graphics::plot(
  base_rows$particle_mis_psis_k,
  abs(base_rows$centered_error_local),
  pch = 21,
  bg = grDevices::adjustcolor("firebrick3", 0.55),
  col = "white",
  xlab = "post-outer particle-MIS PSIS k",
  ylab = "|local centered error|",
  main = "Does PSIS locate error?"
)
graphics::abline(v = max_particle_mis_psis_k, col = "grey35", lty = 2)
graphics::plot(
  base_rows$particle_mis_ess_frac,
  abs(base_rows$centered_error_local),
  pch = 21,
  bg = grDevices::adjustcolor("firebrick3", 0.55),
  col = "white",
  xlab = "post-outer particle-MIS ESS fraction",
  ylab = "|local centered error|",
  main = "Does ESS locate error?"
)
graphics::abline(v = min_particle_mis_ess, col = "grey35", lty = 2)
if (all(c("post_outer", "targeted_repair") %in% rows$member)) {
  paired <- merge(
    rows[rows$member == "post_outer", c("local_pos", "theta_row", "centered_error_local")],
    rows[rows$member == "targeted_repair", c("local_pos", "theta_row", "centered_error_local")],
    by = c("local_pos", "theta_row"),
    suffixes = c("_before", "_after")
  )
  lim <- range(abs(paired$centered_error_local_before), abs(paired$centered_error_local_after), na.rm = TRUE)
  graphics::plot(
    abs(paired$centered_error_local_before),
    abs(paired$centered_error_local_after),
    pch = 21,
    bg = grDevices::adjustcolor("steelblue4", 0.65),
    col = "white",
    xlab = "before targeted repair |centered error|",
    ylab = "after targeted repair |centered error|",
    xlim = lim,
    ylim = lim,
    main = "Repair effect"
  )
  graphics::abline(0, 1, lty = 2, col = "grey35")
} else {
  graphics::plot.new()
  graphics::title("Repair effect unavailable")
}
graphics::mtext("EMC local m_i(theta) hypothesis benchmark against nested local SMC", outer = TRUE, font = 2)
graphics::par(old_par)
grDevices::dev.off()

grDevices::png(surface_plot_file, width = 2400, height = 1500, res = 150, bg = "white")
local_levels <- unique(rows$local)
layout_rows <- ceiling(length(local_levels) / 2)
old_par <- graphics::par(mfrow = c(layout_rows, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 1.7, 0))
for (loc in local_levels) {
  graphics::plot(
    NA,
    xlim = range(rows$theta_row),
    ylim = range(rows$centered_error_local[rows$local == loc], na.rm = TRUE),
    xlab = "theta row",
    ylab = "local centered atlas - reference",
    main = paste("local", loc)
  )
  graphics::abline(h = 0, col = "grey35")
  for (method in present_methods) {
    df <- rows[rows$local == loc & rows$member == method, , drop = FALSE]
    df <- df[order(df$theta_row), , drop = FALSE]
    graphics::lines(df$theta_row, df$centered_error_local, type = "b", pch = 19, col = cols[method], lwd = 2)
  }
  graphics::legend("topright", legend = present_methods, col = cols[present_methods], lwd = 2, pch = 19, bty = "n", cex = 0.75)
}
graphics::mtext("Local surface shape error by theta; additive local offsets removed", outer = TRUE, font = 2)
graphics::par(old_par)
grDevices::dev.off()

saveRDS(
  list(
    audit = audit,
    rows = rows,
    method_summary = method_summary,
    theta_totals = theta_totals,
    repair = repair,
    factor_sets = factor_sets,
    theta_design = theta_design,
    local_ids = local_ids,
    settings = list(
      label = label,
      seed = seed,
      ref_particles = ref_particles,
      ref_reps = ref_reps,
      repair_particles = repair_particles,
      repair_max_updates_per_local = repair_max_updates_per_local,
      min_particle_mis_ess = min_particle_mis_ess,
      max_particle_mis_psis_k = max_particle_mis_psis_k,
      runtime_minutes = as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    )
  ),
  results_file
)

cat("Summary:\n")
print(method_summary, row.names = FALSE)
cat("Repair probes selected/activated:", sum(repair$probes$selected %in% TRUE), "/", sum(repair$probes$activation_success %in% TRUE), "\n")
cat("Saved rows:", rows_csv, "\n")
cat("Saved summary:", summary_csv, "\n")
cat("Saved diagnostics:", plot_file, "\n")
cat("Saved surface plot:", surface_plot_file, "\n")
cat("Saved results:", results_file, "\n")
