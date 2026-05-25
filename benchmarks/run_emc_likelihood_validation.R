#!/usr/bin/env Rscript

rm(list = ls())

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[1L]))
} else {
  normalizePath("benchmarks/run_emc_likelihood_validation.R")
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

arg_int_vec <- function(args, key, default) {
  val <- args[[key]]
  if (is.null(val) || !nzchar(val)) return(as.integer(default))
  as.integer(strsplit(val, ",", fixed = TRUE)[[1L]])
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

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
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

finite_centered_rmse <- function(x) {
  x <- as.numeric(x)
  ok <- is.finite(x)
  if (!any(ok)) return(NA_real_)
  y <- x[ok] - mean(x[ok])
  sqrt(mean(y^2))
}

signed_log10p <- function(x) {
  x <- as.numeric(x)
  sign(x) * log10(1 + abs(x))
}

normalize_log_weights <- function(logw) {
  logw <- as.numeric(logw)
  out <- rep(0, length(logw))
  ok <- is.finite(logw)
  if (!any(ok)) return(rep(NA_real_, length(logw)))
  shifted <- exp(logw[ok] - max(logw[ok]))
  total <- sum(shifted)
  if (!is.finite(total) || total <= 0) return(rep(NA_real_, length(logw)))
  out[ok] <- shifted / total
  out
}

clip_signed <- function(x, limit) {
  x <- as.numeric(x)
  pmax(pmin(x, limit), -limit)
}

signed_limit <- function(x, default = 1) {
  x <- abs(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x)) max(x) else default
}

source("smc_core.R")
source("population_models.R")
source("local_charts.R")
source("utilities.R")

cli_args <- parse_cli_args(commandArgs(trailingOnly = TRUE))
detected_cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
if (!is.finite(detected_cores) || detected_cores < 1L) detected_cores <- 1L

label <- arg_chr(cli_args, "label", "emc_likelihood_validation")
seed <- arg_int(cli_args, "seed", 20260524L)
cores <- arg_int(cli_args, "cores", min(4L, detected_cores))
data_file <- arg_chr(cli_args, "data_file", file.path("benchmarks", "samples", "full_EMC2.RData"))
local_positions <- arg_int_vec(cli_args, "local_positions", c(4L, 14L))
atlas_particles <- arg_int(cli_args, "atlas_particles", 512L)
atlas_max_anchors <- arg_int(cli_args, "atlas_max_anchors", 24L)
atlas_max_steps <- arg_int(cli_args, "atlas_max_steps", 176L)
atlas_max_intermediates <- arg_int(cli_args, "atlas_max_intermediates", 4L)
atlas_target_cess <- arg_num(cli_args, "atlas_target_cess", 0.9)
atlas_mcmc_moves <- arg_int(cli_args, "atlas_mcmc_moves", 2L)
ref_particles <- arg_int(cli_args, "ref_particles", 1536L)
ref_reps <- arg_int(cli_args, "ref_reps", 3L)
ref_max_steps <- arg_int(cli_args, "ref_max_steps", 176L)
ref_target_cess <- arg_num(cli_args, "ref_target_cess", 0.9)
ref_mcmc_moves <- arg_int(cli_args, "ref_mcmc_moves", 2L)
repair_particles <- arg_int(cli_args, "repair_particles", atlas_particles)
repair_max_updates <- arg_int(cli_args, "repair_max_updates", 99L)
repair_confirmation_reps <- arg_int(cli_args, "repair_confirmation_reps", 0L)
repair_confirmation_particles <- arg_int(cli_args, "repair_confirmation_particles", repair_particles)
repair_adaptive_confirmation_reps <- arg_int(cli_args, "repair_adaptive_confirmation_reps", 2L)
repair_replicate_bootstrap_B <- arg_int(cli_args, "repair_replicate_bootstrap_B", 200L)
repair_max_direct_graph_z <- arg_num(cli_args, "repair_max_direct_graph_z", 3)
repair_max_direct_graph_chart_shift <- arg_num(cli_args, "repair_max_direct_graph_chart_shift", 0.35)
repair_max_direct_graph_existing_shift <- arg_num(cli_args, "repair_max_direct_graph_existing_shift", 0.15)
repair_adaptive_replicate_weight_multiplier <- arg_num(cli_args, "repair_adaptive_replicate_weight_multiplier", 3)
repair_adaptive_replicate_min_theta_weight <- arg_num(cli_args, "repair_adaptive_replicate_min_theta_weight", 0)
repair_adaptive_replicate_max_graph_z <- arg_num(cli_args, "repair_adaptive_replicate_max_graph_z", 3)
repair_adaptive_replicate_graph_shift <- arg_num(cli_args, "repair_adaptive_replicate_graph_shift", 0.25)
normalizer_robust_method <- arg_chr(cli_args, "normalizer_robust_method", "student_t")
normalizer_student_t_df <- arg_num(cli_args, "normalizer_student_t_df", 30)
if (!normalizer_robust_method %in% c("student_t", "huber", "none")) {
  stop("normalizer_robust_method must be one of: student_t, huber, none.")
}
options(
  local_charts.normalizer_robust = normalizer_robust_method != "none",
  local_charts.normalizer_robust_method = normalizer_robust_method,
  local_charts.normalizer_student_t_df = normalizer_student_t_df
)
audit_max_points <- arg_int(cli_args, "audit_max_points", 7L)
atlas_design_probs <- arg_num_vec(cli_args, "atlas_design_probs", c(0.05, 0.25, 0.5, 0.75, 0.95))
audit_probs <- arg_num_vec(cli_args, "audit_probs", c(0.01, 0.10, 0.33, 0.50, 0.67, 0.90, 0.99))
focus_parameters <- arg_chr_vec(
  cli_args,
  "focus_parameters",
  c("mu_sv", "sigma2_sv", "sigma2_v_LogFreq", "sigma2_v")
)
distance_metric <- arg_chr(cli_args, "distance_metric", "fisher")
max_chart_distance <- arg_num(cli_args, "max_chart_distance", 1.0)
min_covering_charts <- arg_int(cli_args, "min_covering_charts", 3L)
min_particle_mis_ess <- arg_num(cli_args, "min_particle_mis_ess", 0.05)
max_particle_mis_psis_k <- arg_num(cli_args, "max_particle_mis_psis_k", 0.7)
sparse_chart_min_covering <- arg_int(cli_args, "sparse_chart_min_covering", 3L)
sparse_chart_max_distance <- arg_num(cli_args, "sparse_chart_max_distance", 0.05)
refresh_reference_cache <- arg_lgl(cli_args, "refresh_reference_cache", FALSE)

results_file <- arg_chr(cli_args, "results_file", file.path("benchmarks", "results", paste0(label, "_results.rds")))
rows_csv <- arg_chr(cli_args, "rows_csv", file.path("benchmarks", "results", paste0(label, "_rows.csv")))
summary_csv <- arg_chr(cli_args, "summary_csv", file.path("benchmarks", "results", paste0(label, "_summary.csv")))
reference_cache_file <- arg_chr(cli_args, "reference_cache_file", file.path("benchmarks", "results", paste0(label, "_nested_reference_cache.rds")))
error_plot_file <- arg_chr(cli_args, "error_plot_file", file.path("benchmarks", "results", paste0(label, "_nested_error_comparison.png")))
heat_plot_file <- arg_chr(cli_args, "heat_plot_file", file.path("benchmarks", "results", paste0(label, "_theta_error_heat.png")))
repair_plot_file <- arg_chr(cli_args, "repair_plot_file", file.path("benchmarks", "results", paste0(label, "_repair_effect.png")))
critical_plot_file <- arg_chr(cli_args, "critical_plot_file", file.path("benchmarks", "results", paste0(label, "_critical_failures.png")))
promising_summary_csv <- arg_chr(cli_args, "promising_summary_csv", file.path("benchmarks", "results", paste0(label, "_promising_summary.csv")))
posterior_summary_csv <- arg_chr(cli_args, "posterior_summary_csv", file.path("benchmarks", "results", paste0(label, "_posterior_overlap_summary.csv")))
promising_heat_plot_file <- arg_chr(cli_args, "promising_heat_plot_file", file.path("benchmarks", "results", paste0(label, "_promising_error_heat.png")))
posterior_overlap_plot_file <- arg_chr(cli_args, "posterior_overlap_plot_file", file.path("benchmarks", "results", paste0(label, "_posterior_overlap.png")))
posterior_heat_plot_file <- arg_chr(cli_args, "posterior_heat_plot_file", file.path("benchmarks", "results", paste0(label, "_posterior_mass_heat.png")))
mll_error_heatmap_file <- arg_chr(cli_args, "mll_error_heatmap_file", file.path("benchmarks", "results", paste0(label, "_mll_error_heatmap.png")))
mll_error_density_file <- arg_chr(cli_args, "mll_error_density_file", file.path("benchmarks", "results", paste0(label, "_mll_error_density.png")))
posterior_weight_heatmap_file <- arg_chr(cli_args, "posterior_weight_heatmap_file", file.path("benchmarks", "results", paste0(label, "_posterior_weight_heatmap.png")))
posterior_density_1d_file <- arg_chr(cli_args, "posterior_density_1d_file", file.path("benchmarks", "results", paste0(label, "_posterior_density_1d.png")))
posterior_density_2d_file <- arg_chr(cli_args, "posterior_density_2d_file", file.path("benchmarks", "results", paste0(label, "_posterior_density_2d.png")))

for (path in c(
  results_file, rows_csv, summary_csv, reference_cache_file, error_plot_file,
  heat_plot_file, repair_plot_file, critical_plot_file, promising_summary_csv,
  posterior_summary_csv, promising_heat_plot_file, posterior_overlap_plot_file,
  posterior_heat_plot_file, mll_error_heatmap_file, mll_error_density_file,
  posterior_weight_heatmap_file, posterior_density_1d_file, posterior_density_2d_file
)) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
}

if (!file.exists(data_file)) {
  stop("Missing EMC2 benchmark data: ", data_file)
}
load(data_file)
if (!exists("ELP_DDM", inherits = FALSE)) {
  stop("The EMC2 data file must define ELP_DDM.")
}

emc <- ELP_DDM[[1L]]
data_list <- emc$data
model_factory <- emc$model
alpha_names <- emc$par_names
alpha_dim <- length(alpha_names)

local_positions <- local_positions[is.finite(local_positions)]
if (!length(local_positions) || any(local_positions < 1L | local_positions > length(data_list))) {
  stop("local_positions must identify EMC locals.")
}
local_positions <- unique(as.integer(local_positions))
local_ids <- names(data_list)[local_positions] %||% as.character(local_positions)

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
theta_cloud <- local_atlas_theta_from_draws(emc_draws, population_model)
theta_root <- matrix(apply(theta_cloud, 2L, stats::median), nrow = 1L)
colnames(theta_root) <- population_model$hyper_names

focus_hyper_names <- unique(c(
  focus_parameters[startsWith(focus_parameters, "mu_")],
  sub("^sigma2_", "log_sigma2_", focus_parameters[startsWith(focus_parameters, "sigma2_")])
))
focus_hyper_names <- intersect(focus_hyper_names, population_model$hyper_names)
if (!length(focus_hyper_names)) {
  focus_hyper_names <- population_model$hyper_names
}

atlas_design <- build_local_evidence_calibration_design(
  theta = theta_cloud,
  population_model = population_model,
  weights = rep(1 / nrow(theta_cloud), nrow(theta_cloud)),
  focus_hyper_names = focus_hyper_names,
  probs = atlas_design_probs,
  max_points = atlas_max_anchors,
  n_leverage_points = 2L,
  n_uncertainty_points = 0L,
  n_disagreement_points = 0L,
  source = "emc_likelihood_validation_atlas_design",
  label_prefix = "atlas"
)

audit_design <- build_local_evidence_calibration_design(
  theta = theta_cloud,
  population_model = population_model,
  weights = rep(1 / nrow(theta_cloud), nrow(theta_cloud)),
  focus_hyper_names = focus_hyper_names,
  probs = audit_probs,
  max_points = audit_max_points,
  n_leverage_points = 1L,
  n_uncertainty_points = 0L,
  n_disagreement_points = 0L,
  source = "emc_likelihood_validation_theta_audit",
  label_prefix = "theta"
)

cat(sprintf(
  "EMC likelihood validation: %d locals x %d theta profiles\n",
  length(local_positions), nrow(audit_design$theta)
))
cat(sprintf(
  "Locals: %s\n",
  paste(sprintf("%d/%s(n=%d)", local_positions, local_ids, vapply(data_list[local_positions], nrow, integer(1))), collapse = ", ")
))
cat(sprintf(
  "Reference nested SMC: %d reps x %d particles per local/theta\n",
  ref_reps, ref_particles
))

query_native_eval <- function(atlas, theta) {
  evaluate_local_atlas(
    atlas = atlas,
    theta = theta,
    population_model = population_model,
    max_chart_distance = max_chart_distance,
    min_covering_charts = min_covering_charts,
    max_prediction_range = Inf,
    distance_scale = max_chart_distance,
    se_floor = 1e-6,
    use_particle_mis = TRUE,
    require_particle_mis = TRUE,
    min_particle_mis_ess = min_particle_mis_ess,
    max_particle_mis_psis_k = max_particle_mis_psis_k,
    max_quadratic_particle_gap = Inf,
    sparse_chart_min_covering = sparse_chart_min_covering,
    sparse_chart_max_distance = sparse_chart_max_distance,
    max_leave_chart_out_gap = Inf,
    distance_metric = distance_metric,
    surface_method = "derivative_ls",
    min_surface_charts = 2L,
    max_surface_se = Inf,
    surface_value_nugget = 0.05,
    surface_gradient_weight = 1.0,
    surface_curvature_weight = 0.2,
    surface_ridge = 1e-8,
    particle_mis_role = "estimator"
  )
}

eval_rows_for_atlas <- function(atlas, local_pos, local_id, method, eval_fn) {
  out <- lapply(seq_len(nrow(audit_design$theta)), function(theta_row) {
    ev <- eval_fn(atlas, audit_design$theta[theta_row, , drop = FALSE])
    data.frame(
      method = method,
      local_pos = local_pos,
      local = local_id,
      theta_row = theta_row,
      theta_label = audit_design$metadata$theta_label[theta_row],
      estimate = ev$log_marginal,
      se = ev$se,
      status = ev$status,
      reason = ev$reason,
      nearest_charts = paste(ev$nearest_charts, collapse = ","),
      n_covering_charts = {
        rows <- ev$diagnostics$predictions %||% data.frame(within_distance = logical())
        if ("within_distance" %in% names(rows)) sum(rows$within_distance %in% TRUE) else nrow(rows)
      },
      particle_mis_ess_frac = as.numeric(ev$diagnostics$particle_mis$ess_frac %||% NA_real_),
      particle_mis_psis_k = as.numeric(ev$diagnostics$particle_mis$psis_k %||% NA_real_),
      surface_log_marginal = as.numeric(ev$diagnostics$surface_log_marginal %||% NA_real_),
      surface_se = as.numeric(ev$diagnostics$surface_se %||% NA_real_),
      surface_particle_gap = as.numeric(ev$diagnostics$quadratic_particle_gap %||% NA_real_),
      min_covering_distance = as.numeric(ev$diagnostics$min_covering_distance %||% NA_real_),
      check.names = FALSE
    )
  })
  do.call(rbind, out)
}

reference_settings <- list(
  local_positions = local_positions,
  theta = audit_design$theta,
  ref_particles = ref_particles,
  ref_reps = ref_reps,
  ref_max_steps = ref_max_steps,
  ref_target_cess = ref_target_cess,
  ref_mcmc_moves = ref_mcmc_moves,
  alpha_names = alpha_names
)

cache_valid <- FALSE
if (file.exists(reference_cache_file) && !isTRUE(refresh_reference_cache)) {
  cache <- readRDS(reference_cache_file)
  cache_valid <- isTRUE(identical(cache$settings$local_positions, reference_settings$local_positions)) &&
    isTRUE(all.equal(cache$settings$theta, reference_settings$theta, tolerance = 1e-12)) &&
    identical(cache$settings$ref_particles, reference_settings$ref_particles) &&
    identical(cache$settings$ref_reps, reference_settings$ref_reps) &&
    identical(cache$settings$ref_max_steps, reference_settings$ref_max_steps) &&
    identical(cache$settings$alpha_names, reference_settings$alpha_names)
  if (isTRUE(cache_valid)) {
    reference_raw <- cache$reference_raw
    cat(sprintf("Loaded nested-SMC reference cache: %s\n", reference_cache_file))
  }
}

if (!isTRUE(cache_valid)) {
  jobs <- expand.grid(
    local_pos = local_positions,
    theta_row = seq_len(nrow(audit_design$theta)),
    replicate = seq_len(ref_reps),
    KEEP.OUT.ATTRS = FALSE
  )
  run_reference <- function(k) {
    local_pos <- jobs$local_pos[k]
    theta_row <- jobs$theta_row[k]
    rep_id <- jobs$replicate[k]
    local_id <- names(data_list)[local_pos] %||% as.character(local_pos)
    run <- .local_chart_run_smc(
      local_id = local_id,
      theta_anchor = audit_design$theta[theta_row, , drop = FALSE],
      data_i = data_list[[local_pos]],
      loglik_fn = loglik_emc2,
      population_model = population_model,
      M = ref_particles,
      target_cess = ref_target_cess,
      resample_threshold = 0.5,
      n_mcmc_moves = ref_mcmc_moves,
      rw_scale = 0.75,
      G_mix = 8L,
      da_enable = TRUE,
      refit_every = 2L,
      max_steps = ref_max_steps,
      deterministic_resampling = FALSE,
      n_cores = 1L,
      seed = seed + 1000003L * local_pos + 9176L * theta_row + 104729L * rep_id,
      verbose = FALSE,
      source = "emc_likelihood_validation_nested_reference"
    )
    data.frame(
      local_pos = local_pos,
      local = local_id,
      theta_row = theta_row,
      theta_label = audit_design$metadata$theta_label[theta_row],
      replicate = rep_id,
      reference_log_marginal = run$logZ,
      reference_logZ_se = run$logZ_se,
      reference_final_ess_frac = run$diagnostics$final_ess_frac,
      reference_min_path_ess_frac = run$diagnostics$min_path_ess_frac,
      check.names = FALSE
    )
  }
  cat(sprintf("Running %d nested-SMC reference jobs...\n", nrow(jobs)))
  reference_parts <- if (cores <= 1L || nrow(jobs) <= 1L) {
    lapply(seq_len(nrow(jobs)), run_reference)
  } else {
    parallel::mclapply(seq_len(nrow(jobs)), run_reference, mc.cores = min(cores, nrow(jobs)))
  }
  reference_raw <- do.call(rbind, reference_parts)
  saveRDS(
    list(reference_raw = reference_raw, settings = reference_settings, theta_metadata = audit_design$metadata),
    reference_cache_file
  )
  cat(sprintf("Saved nested-SMC reference cache: %s\n", reference_cache_file))
}

reference_summary <- do.call(rbind, lapply(split(reference_raw, interaction(reference_raw$local_pos, reference_raw$theta_row, drop = TRUE)), function(df) {
  data.frame(
    local_pos = df$local_pos[1L],
    local = df$local[1L],
    theta_row = df$theta_row[1L],
    theta_label = df$theta_label[1L],
    reference_log_marginal = mean(df$reference_log_marginal),
    reference_rep_sd = if (nrow(df) > 1L) stats::sd(df$reference_log_marginal) else 0,
    reference_rep_range = diff(range(df$reference_log_marginal)),
    reference_mean_logZ_se = sqrt(mean(pmax(df$reference_logZ_se, 0)^2) / nrow(df)),
    reference_n = nrow(df),
    check.names = FALSE
  )
}))

atlas_objects <- list()
repair_objects <- list()
method_rows <- list()
repair_tables <- list()

for (idx in seq_along(local_positions)) {
  local_pos <- local_positions[idx]
  local_id <- local_ids[idx]
  cat(sprintf("Building atlas for local %d/%s (%d obs)\n", local_pos, local_id, nrow(data_list[[local_pos]])))
  atlas <- build_local_atlas(
    local_id = local_id,
    data_i = data_list[[local_pos]],
    loglik_fn = loglik_emc2,
    population_model = population_model,
    theta_design = atlas_design$theta,
    theta_root = theta_root,
    local_control = list(
      root_M = atlas_particles,
      candidate_M = atlas_particles,
      bridge_particles = atlas_particles,
      target_cess = atlas_target_cess,
      n_mcmc_moves = atlas_mcmc_moves,
      max_steps = atlas_max_steps,
      root_confirm = "auto"
    ),
    edge_control = list(
      max_intermediates = atlas_max_intermediates,
      min_overlap_ess = 0.05,
      max_pareto_k = 0.9,
      max_forward_reverse_gap = 2,
      max_taylor_gap = Inf,
      max_cycle_z = 4,
      distance_metric = distance_metric
    ),
    n_cores = 1L,
    seed = seed + 100003L * local_pos,
    verbose = FALSE
  )
  atlas_objects[[as.character(local_pos)]] <- atlas
  method_rows[[length(method_rows) + 1L]] <- eval_rows_for_atlas(
    atlas,
    local_pos,
    local_id,
    "query_native",
    query_native_eval
  )

  factor_set <- build_local_atlas_factor_set(
    atlases = setNames(list(atlas), local_id),
    population_model = population_model,
    max_chart_distance = max_chart_distance,
    min_covering_charts = min_covering_charts,
    max_prediction_range = Inf,
    distance_scale = max_chart_distance,
    se_floor = 1e-6,
    use_particle_mis = TRUE,
    require_particle_mis = TRUE,
    min_particle_mis_ess = min_particle_mis_ess,
    max_particle_mis_psis_k = max_particle_mis_psis_k,
    max_quadratic_particle_gap = Inf,
    sparse_chart_min_covering = sparse_chart_min_covering,
    sparse_chart_max_distance = sparse_chart_max_distance,
    max_leave_chart_out_gap = Inf,
    distance_metric = distance_metric,
    surface_method = "derivative_ls",
    particle_mis_role = "estimator",
    stop_on_uncertified = FALSE
  )
  repair <- local_atlas_calibrate_posterior_regions(
    factor_set = factor_set,
    theta = audit_design$theta,
    data_list = setNames(list(data_list[[local_pos]]), local_id),
    loglik_fn = loglik_emc2,
    theta_weights = rep(1 / nrow(audit_design$theta), nrow(audit_design$theta)),
    local_ids = 1L,
    M = repair_particles,
    target_cess = atlas_target_cess,
    n_mcmc_moves = atlas_mcmc_moves,
    max_steps = atlas_max_steps,
    max_updates = min(repair_max_updates, nrow(audit_design$theta)),
    abs_delta_threshold = 0.75,
    z_threshold = 4,
    candidate_pool_multiplier = 1L,
    confirmation_reps = repair_confirmation_reps,
    confirmation_M = repair_confirmation_particles,
    confirmation_max_sd = 1.5,
    adaptive_confirmation_reps = repair_adaptive_confirmation_reps,
    replicate_bootstrap_B = repair_replicate_bootstrap_B,
    max_direct_graph_z = repair_max_direct_graph_z,
    max_direct_graph_chart_shift = repair_max_direct_graph_chart_shift,
    max_direct_graph_existing_shift = repair_max_direct_graph_existing_shift,
    adaptive_replicate_weight_multiplier = repair_adaptive_replicate_weight_multiplier,
    adaptive_replicate_min_theta_weight = repair_adaptive_replicate_min_theta_weight,
    adaptive_replicate_max_graph_z = repair_adaptive_replicate_max_graph_z,
    adaptive_replicate_graph_shift = repair_adaptive_replicate_graph_shift,
    local_control = list(
      candidate_M = repair_particles,
      bridge_particles = repair_particles,
      target_cess = atlas_target_cess,
      n_mcmc_moves = atlas_mcmc_moves,
      max_steps = atlas_max_steps
    ),
    edge_control = list(
      max_intermediates = atlas_max_intermediates,
      min_overlap_ess = 0.05,
      max_pareto_k = 0.9,
      max_forward_reverse_gap = 2,
      max_taylor_gap = Inf,
      max_cycle_z = 4,
      distance_metric = distance_metric
    ),
    n_cores = 1L,
    seed = seed + 900001L + 100003L * local_pos,
    verbose = FALSE
  )
  repaired_atlas <- repair$factor_set$atlases[[1L]]
  repair_objects[[as.character(local_pos)]] <- repaired_atlas
  repair_table <- repair$probes
  repair_table$local_pos <- local_pos
  repair_table$local <- local_id
  repair_tables[[length(repair_tables) + 1L]] <- repair_table
  method_rows[[length(method_rows) + 1L]] <- eval_rows_for_atlas(
    repaired_atlas,
    local_pos,
    local_id,
    "query_native_repaired",
    query_native_eval
  )
}

estimate_rows <- do.call(rbind, method_rows)
estimate_rows <- merge(estimate_rows, reference_summary, by = c("local_pos", "local", "theta_row", "theta_label"), all.x = TRUE, sort = FALSE)
estimate_rows$error <- estimate_rows$estimate - estimate_rows$reference_log_marginal
estimate_rows$centered_error <- ave(
  estimate_rows$error,
  interaction(estimate_rows$method, estimate_rows$local_pos, drop = TRUE),
  FUN = function(x) x - mean(x[is.finite(x)], na.rm = TRUE)
)
estimate_rows$abs_error <- abs(estimate_rows$error)
estimate_rows$certified <- estimate_rows$status == "certified"
estimate_rows <- merge(
  estimate_rows,
  audit_design$metadata[, c("theta_label", population_model$hyper_names), drop = FALSE],
  by = "theta_label",
  all.x = TRUE,
  sort = FALSE
)

method_summary <- do.call(rbind, lapply(split(estimate_rows, estimate_rows$method), function(df) {
  data.frame(
    method = df$method[1L],
    n_rows = nrow(df),
    certified_fraction = mean(df$certified, na.rm = TRUE),
    finite_fraction = mean(is.finite(df$error)),
    rmse = finite_rmse(df$error),
    centered_rmse = finite_centered_rmse(df$error),
    mae = finite_mae(df$error),
    max_abs_error = finite_max_abs(df$error),
    median_particle_mis_ess_frac = stats::median(df$particle_mis_ess_frac, na.rm = TRUE),
    median_particle_mis_psis_k = stats::median(df$particle_mis_psis_k, na.rm = TRUE),
    median_surface_particle_gap = stats::median(df$surface_particle_gap, na.rm = TRUE),
    check.names = FALSE
  )
}))
method_order <- c("query_native", "query_native_repaired")
method_summary <- method_summary[order(match(method_summary$method, method_order)), , drop = FALSE]

repair_rows <- if (length(repair_tables)) do.call(rbind, repair_tables) else data.frame()

promising_methods <- c("query_native", "query_native_repaired")
promising_method_summary <- method_summary[method_summary$method %in% promising_methods, , drop = FALSE]

build_posterior_diagnostics <- function() {
  theta_rows <- sort(unique(reference_summary$theta_row))
  theta_meta <- audit_design$metadata
  theta_meta$theta_row <- seq_len(nrow(theta_meta))
  block_defs <- c(
    setNames(as.list(local_positions), paste0("local_", local_positions)),
    list(combined_selected_locals = local_positions)
  )
  detail_parts <- list()
  summary_parts <- list()
  for (block_key in names(block_defs)) {
    block_locals <- as.integer(block_defs[[block_key]])
    block_label <- if (length(block_locals) == 1L) {
      local_id <- unique(reference_summary$local[reference_summary$local_pos == block_locals])
      sprintf("local %d/%s", block_locals, local_id[1L])
    } else {
      sprintf("combined locals %s", paste(block_locals, collapse = "+"))
    }
    ref_logw <- vapply(theta_rows, function(theta_row) {
      vals <- reference_summary$reference_log_marginal[
        reference_summary$local_pos %in% block_locals &
          reference_summary$theta_row == theta_row
      ]
      if (length(vals) != length(block_locals) || any(!is.finite(vals))) return(-Inf)
      sum(vals)
    }, numeric(1))
    ref_w <- normalize_log_weights(ref_logw)
    for (method in promising_methods) {
      method_rows <- estimate_rows[
        estimate_rows$method == method &
          estimate_rows$local_pos %in% block_locals,
        ,
        drop = FALSE
      ]
      est_logw <- vapply(theta_rows, function(theta_row) {
        vals <- method_rows$estimate[method_rows$theta_row == theta_row]
        if (length(vals) != length(block_locals) || any(!is.finite(vals))) return(-Inf)
        sum(vals)
      }, numeric(1))
      certified <- vapply(theta_rows, function(theta_row) {
        vals <- method_rows$certified[method_rows$theta_row == theta_row]
        length(vals) == length(block_locals) && all(vals %in% TRUE)
      }, logical(1))
      est_w <- normalize_log_weights(est_logw)
      finite_weights <- is.finite(ref_w) & is.finite(est_w)
      overlap <- if (any(finite_weights)) sum(pmin(ref_w[finite_weights], est_w[finite_weights])) else NA_real_
      tv <- if (any(finite_weights)) 0.5 * sum(abs(ref_w[finite_weights] - est_w[finite_weights])) else NA_real_
      hellinger <- if (any(finite_weights)) {
        sqrt(0.5 * sum((sqrt(ref_w[finite_weights]) - sqrt(est_w[finite_weights]))^2))
      } else {
        NA_real_
      }
      detail <- data.frame(
        block_key = block_key,
        block_label = block_label,
        block_n_locals = length(block_locals),
        method = method,
        theta_row = theta_rows,
        reference_log_weight = ref_logw,
        method_log_weight = est_logw,
        log_weight_error = est_logw - ref_logw,
        reference_weight = ref_w,
        method_weight = est_w,
        weight_diff = est_w - ref_w,
        certified = certified,
        check.names = FALSE
      )
      detail <- merge(detail, theta_meta, by = "theta_row", all.x = TRUE, sort = FALSE)
      detail_parts[[length(detail_parts) + 1L]] <- detail
      summary_parts[[length(summary_parts) + 1L]] <- data.frame(
        block_key = block_key,
        block_label = block_label,
        block_n_locals = length(block_locals),
        method = method,
        overlap = overlap,
        total_variation = tv,
        hellinger = hellinger,
        reference_mass_uncertified = sum(ref_w[!certified], na.rm = TRUE),
        top_reference_theta = theta_rows[which.max(ref_w)],
        top_method_theta = theta_rows[which.max(est_w)],
        max_abs_weight_diff = max(abs(est_w - ref_w), na.rm = TRUE),
        check.names = FALSE
      )
    }
  }
  list(
    detail = do.call(rbind, detail_parts),
    summary = do.call(rbind, summary_parts)
  )
}

posterior_diagnostics <- build_posterior_diagnostics()

utils::write.csv(estimate_rows, rows_csv, row.names = FALSE)
utils::write.csv(method_summary, summary_csv, row.names = FALSE)
utils::write.csv(promising_method_summary, promising_summary_csv, row.names = FALSE)
utils::write.csv(posterior_diagnostics$summary, posterior_summary_csv, row.names = FALSE)

point_col <- c(
  query_native = "darkorange3",
  query_native_repaired = "steelblue4"
)
method_label <- c(
  query_native = "Query-native particle-MIS",
  query_native_repaired = "Query-native repaired"
)

plot_reference_replicates <- function() {
  ref <- merge(reference_raw, reference_summary, by = c("local_pos", "local", "theta_row", "theta_label"), all.x = TRUE)
  ref$error <- ref$reference_log_marginal.x - ref$reference_log_marginal.y
  locals <- unique(reference_summary$local_pos)
  grDevices::png(error_plot_file, width = 1900, height = 650 * length(locals))
  graphics::par(mfrow = c(length(locals), 1L), mar = c(4.4, 4.8, 3.2, 1.2))
  for (local_pos in locals) {
    ref_l <- ref[ref$local_pos == local_pos, , drop = FALSE]
    est_l <- estimate_rows[estimate_rows$local_pos == local_pos, , drop = FALSE]
    y_lim <- range(c(ref_l$error, est_l$error), na.rm = TRUE)
    graphics::plot(
      NA,
      xlim = range(reference_summary$theta_row),
      ylim = y_lim,
      xlab = "theta profile",
      ylab = "local log m error vs nested-SMC mean",
      main = sprintf("EMC local %d/%s: expensive nested-SMC reference vs atlas methods", local_pos, ref_l$local[1L])
    )
    graphics::abline(h = 0, lty = 2, col = "grey45")
    rep_center <- mean(seq_len(max(ref_l$replicate, na.rm = TRUE)))
    graphics::points(
      ref_l$theta_row + 0.055 * (ref_l$replicate - rep_center),
      ref_l$error,
      pch = 16,
      cex = 0.8,
      col = grDevices::adjustcolor("grey25", 0.55)
    )
    ref_s <- reference_summary[reference_summary$local_pos == local_pos, , drop = FALSE]
    for (i in seq_len(nrow(ref_s))) {
      graphics::segments(
        ref_s$theta_row[i] - 0.24,
        -ref_s$reference_rep_sd[i],
        ref_s$theta_row[i] + 0.24,
        -ref_s$reference_rep_sd[i],
        col = "grey60"
      )
      graphics::segments(
        ref_s$theta_row[i] - 0.24,
        ref_s$reference_rep_sd[i],
        ref_s$theta_row[i] + 0.24,
        ref_s$reference_rep_sd[i],
        col = "grey60"
      )
    }
    offsets <- c(query_native = -0.12, query_native_repaired = 0.12)
    for (method in names(offsets)) {
      df <- est_l[est_l$method == method, , drop = FALSE]
      graphics::points(
        df$theta_row + offsets[[method]],
        df$error,
        pch = ifelse(df$certified, 17, 4),
        cex = 1.35,
        col = point_col[[method]],
        lwd = 2
      )
    }
    graphics::legend(
      "bottomleft",
      legend = c("nested-SMC replicate - mean", "+/- one replicate SD", method_label[names(offsets)]),
      pch = c(16, NA, 17, 17),
      lty = c(NA, 1, NA, NA),
      col = c("grey25", "grey60", point_col[names(offsets)]),
      bty = "n",
      cex = 0.9
    )
  }
  grDevices::dev.off()
}

plot_theta_heat <- function() {
  x_name <- if ("mu_sv" %in% names(estimate_rows)) "mu_sv" else focus_hyper_names[1L]
  y_name <- if ("log_sigma2_sv" %in% names(estimate_rows)) "log_sigma2_sv" else focus_hyper_names[min(2L, length(focus_hyper_names))]
  grDevices::png(heat_plot_file, width = 1800, height = 1200)
  methods <- c("query_native", "query_native_repaired")
  graphics::par(mfrow = c(length(local_positions), length(methods)), mar = c(4.3, 4.5, 3.1, 1.0))
  for (local_pos in local_positions) {
    for (method in methods) {
      df <- estimate_rows[estimate_rows$local_pos == local_pos & estimate_rows$method == method, , drop = FALSE]
      lim <- signed_limit(df$error)
      pal <- grDevices::hcl.colors(101, "RdBu", rev = TRUE)
      col_idx <- floor((pmax(pmin(df$error, lim), -lim) + lim) / (2 * lim) * 100) + 1L
      col_idx[!is.finite(col_idx)] <- NA_integer_
      cols <- pal[col_idx]
      cols[is.na(cols)] <- "white"
      graphics::plot(
        df[[x_name]],
        df[[y_name]],
        xlab = x_name,
        ylab = y_name,
        main = sprintf("local %s: %s\nerror color scale +/- %.2f", unique(df$local), method_label[[method]], lim),
        pch = ifelse(df$certified, 21, 4),
        bg = cols,
        col = ifelse(df$certified, "grey20", "firebrick3"),
        cex = 1.8,
        lwd = 2
      )
      graphics::text(df[[x_name]], df[[y_name]], labels = df$theta_row, pos = 3, cex = 0.75)
      graphics::legend(
        "topright",
        legend = c("certified", "uncertified"),
        pch = c(21, 4),
        pt.bg = c("grey80", NA),
        col = c("grey20", "firebrick3"),
        bty = "n",
        cex = 0.8
      )
    }
  }
  grDevices::dev.off()
}

plot_repair_effect <- function() {
  grDevices::png(repair_plot_file, width = 1800, height = 1250)
  graphics::layout(matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE))
  graphics::par(mar = c(4.5, 4.8, 3.3, 1.2))

  bar_mat <- t(as.matrix(method_summary[, c("rmse", "centered_rmse", "mae"), drop = FALSE]))
  colnames(bar_mat) <- method_summary$method
  graphics::barplot(
    bar_mat,
    beside = TRUE,
    las = 2,
    col = c("grey35", "steelblue3", "tan3"),
    ylab = "error vs nested-SMC mean",
    main = "Method error summary over EMC local/theta validation set"
  )
  graphics::legend("topright", legend = rownames(bar_mat), fill = c("grey35", "steelblue3", "tan3"), bty = "n")

  graphics::barplot(
    method_summary$certified_fraction,
    names.arg = method_summary$method,
    las = 2,
    ylim = c(0, 1),
    col = point_col[method_summary$method],
    ylab = "certified fraction",
    main = "Certification is part of the likelihood estimate"
  )
  graphics::abline(h = 1, lty = 2, col = "grey50")

  paired <- reshape(
    estimate_rows[, c("local_pos", "theta_row", "method", "error"), drop = FALSE],
    idvar = c("local_pos", "theta_row"),
    timevar = "method",
    direction = "wide"
  )
  graphics::plot(
    signed_log10p(paired$error.query_native),
    signed_log10p(paired$error.query_native_repaired),
    xlab = "query-native error, signed log10(1 + abs(error))",
    ylab = "repaired query-native error, signed log10(1 + abs(error))",
    main = "Before/after repair against expensive nested-SMC",
    pch = 16,
    col = grDevices::adjustcolor("black", 0.6)
  )
  graphics::abline(0, 1, lty = 2, col = "grey45")
  graphics::abline(h = 0, v = 0, lty = 3, col = "grey70")
  graphics::legend("topleft", legend = "points below diagonal have smaller repaired signed error", bty = "n", cex = 0.85)

  if (nrow(repair_rows)) {
    shown <- repair_rows[repair_rows$selected %in% TRUE | repair_rows$atlas_status != "certified", , drop = FALSE]
    shown <- shown[is.finite(shown$delta_fresh_minus_atlas), , drop = FALSE]
    if (nrow(shown)) {
      shown$theta_key <- paste0(shown$local, ":t", shown$theta_row)
      graphics::plot(
        seq_len(nrow(shown)),
        shown$delta_fresh_minus_atlas,
        xaxt = "n",
        xlab = "selected/failed local-theta repair probe",
        ylab = "fresh endpoint log m - atlas log m",
        main = "Repair probes target uncertified or high-discrepancy local/theta pairs",
        pch = ifelse(shown$selected, 17, 1),
        col = ifelse(shown$activation_success %in% TRUE, "steelblue4", "firebrick3"),
        cex = 1.35,
        lwd = 2
      )
      graphics::abline(h = 0, lty = 2, col = "grey45")
      graphics::axis(1, at = seq_len(nrow(shown)), labels = shown$theta_key, las = 2, cex.axis = 0.75)
      graphics::legend("topright", legend = c("activated", "not activated"), pch = 17, col = c("steelblue4", "firebrick3"), bty = "n")
    } else {
      graphics::plot.new()
      graphics::title("No finite repair probe deltas")
    }
  } else {
    graphics::plot.new()
    graphics::title("No repair rows")
  }
  grDevices::dev.off()
}

plot_critical_failures <- function() {
  methods <- c("query_native", "query_native_repaired")
  method_offsets <- c(query_native = -0.14, query_native_repaired = 0.14)
  method_pch <- c(query_native = 16, query_native_repaired = 15)
  pair_keys <- unique(estimate_rows[, c("local_pos", "local", "theta_row"), drop = FALSE])
  pair_keys <- pair_keys[order(pair_keys$local_pos, pair_keys$theta_row), , drop = FALSE]
  pair_keys$pair_label <- paste0("L", pair_keys$local, ":t", pair_keys$theta_row)
  estimate_rows$pair_id <- match(
    paste(estimate_rows$local_pos, estimate_rows$theta_row),
    paste(pair_keys$local_pos, pair_keys$theta_row)
  )

  grDevices::png(critical_plot_file, width = 1900, height = 1250)
  graphics::layout(matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE))
  graphics::par(mar = c(6.2, 4.8, 3.4, 1.2))

  signed_limit_raw <- 1.25
  graphics::plot(
    NA,
    xlim = c(0.5, nrow(pair_keys) + 0.5),
    ylim = c(-signed_limit_raw - 0.25, signed_limit_raw + 0.25),
    xaxt = "n",
    xlab = "",
    ylab = "error vs nested-SMC mean, clipped",
    main = "Pointwise EMC local likelihood error: query-native vs targeted repair"
  )
  graphics::abline(h = 0, lty = 2, col = "grey50")
  graphics::abline(h = c(-signed_limit_raw, signed_limit_raw), lty = 3, col = "grey75")
  graphics::axis(1, at = seq_len(nrow(pair_keys)), labels = pair_keys$pair_label, las = 2, cex.axis = 0.7)
  for (method in methods) {
    df <- estimate_rows[estimate_rows$method == method, , drop = FALSE]
    x <- df$pair_id + method_offsets[[method]]
    finite <- is.finite(df$error)
    y <- clip_signed(df$error, signed_limit_raw)
    outlier <- finite & abs(df$error) > signed_limit_raw
    graphics::points(
      x[finite],
      y[finite],
      pch = ifelse(outlier[finite], 24, method_pch[[method]]),
      bg = ifelse(outlier[finite], point_col[[method]], NA),
      col = point_col[[method]],
      cex = 1.25,
      lwd = 2
    )
    if (any(outlier)) {
      graphics::text(
        x[outlier],
        y[outlier] + 0.11 * sign(y[outlier]),
        labels = sprintf("%.1f", df$error[outlier]),
        col = point_col[[method]],
        cex = 0.72
      )
    }
    if (any(!finite)) {
      graphics::points(
        x[!finite],
        rep(-signed_limit_raw - 0.16, sum(!finite)),
        pch = 4,
        col = point_col[[method]],
        cex = 1.15,
        lwd = 2
      )
    }
  }
  graphics::legend(
    "topright",
    legend = c(method_label[methods], "uncertified/missing", "clipped outlier"),
    pch = c(method_pch[methods], 4, 24),
    pt.bg = c(rep(NA, length(methods) + 1L), "grey70"),
    col = c(point_col[methods], "grey20", "grey20"),
    bty = "n",
    cex = 0.85
  )

  graphics::par(mar = c(5.0, 4.8, 3.4, 1.2))
  paired <- reshape(
    estimate_rows[, c("local_pos", "local", "theta_row", "method", "abs_error"), drop = FALSE],
    idvar = c("local_pos", "local", "theta_row"),
    timevar = "method",
    direction = "wide"
  )
  graphics::plot(
    signed_log10p(paired$abs_error.query_native),
    signed_log10p(paired$abs_error.query_native_repaired),
    xlab = "query-native abs error, log10(1 + error)",
    ylab = "repaired abs error, log10(1 + error)",
    main = "Repair should reduce evidence error without hiding uncertified queries",
    pch = 16,
    col = grDevices::adjustcolor("grey20", 0.75),
    cex = 1.2,
    lwd = 2
  )
  graphics::text(
    signed_log10p(paired$abs_error.query_native),
    signed_log10p(paired$abs_error.query_native_repaired),
    labels = paste0("L", paired$local, ":t", paired$theta_row),
    pos = 3,
    cex = 0.68
  )
  graphics::abline(0, 1, lty = 3, col = "grey70")

  cert <- table(
    factor(estimate_rows$method, levels = methods),
    ifelse(estimate_rows$certified, "certified", "uncertified")
  )
  graphics::barplot(
    t(cert),
    beside = FALSE,
    col = c("grey35", "firebrick3"),
    ylim = c(0, max(rowSums(cert))),
    ylab = "local/theta queries",
    main = "Certification changed the estimator: missing support is now explicit",
    las = 2
  )
  graphics::legend("topright", legend = colnames(cert), fill = c("grey35", "firebrick3"), bty = "n")

  ref <- reference_summary
  abs_rows <- estimate_rows[is.finite(estimate_rows$error), , drop = FALSE]
  abs_rows$abs_error <- abs(abs_rows$error)
  graphics::plot(
    ref$reference_rep_sd[match(paste(abs_rows$local_pos, abs_rows$theta_row), paste(ref$local_pos, ref$theta_row))],
    abs_rows$abs_error,
    log = "xy",
    xlab = "nested-SMC replicate SD",
    ylab = "atlas abs error vs nested mean",
    main = "Atlas error should be judged against expensive nested-SMC noise",
    pch = method_pch[abs_rows$method],
    col = point_col[abs_rows$method],
    cex = 1.15,
    lwd = 2
  )
  graphics::abline(0, 1, lty = 2, col = "grey50")
  graphics::abline(0, 2, lty = 3, col = "grey70")
  graphics::legend("topleft", legend = method_label[methods], pch = method_pch[methods], col = point_col[methods], bty = "n", cex = 0.85)

  grDevices::dev.off()
}

add_signed_color_legend <- function(lim, pal, label, digits = 2) {
  usr <- graphics::par("usr")
  x_span <- usr[2] - usr[1]
  y_span <- usr[4] - usr[3]
  x0 <- usr[1] + 0.885 * x_span
  x1 <- usr[1] + 0.915 * x_span
  y0 <- usr[3] + 0.130 * y_span
  y1 <- usr[3] + 0.465 * y_span
  n <- length(pal)
  y <- seq(y0, y1, length.out = n + 1L)
  graphics::rect(x0, y[-length(y)], x1, y[-1L], col = pal, border = NA)
  graphics::rect(x0, y0, x1, y1, border = "grey30")
  graphics::text(x1 + 0.018 * x_span, y1, sprintf(paste0("%.", digits, "f"), lim), adj = 0, cex = 0.68)
  graphics::text(x1 + 0.018 * x_span, (y0 + y1) / 2, "0", adj = 0, cex = 0.68)
  graphics::text(x1 + 0.018 * x_span, y0, sprintf(paste0("%.", digits, "f"), -lim), adj = 0, cex = 0.68)
  graphics::text((x0 + x1) / 2, y1 + 0.055 * y_span, label, cex = 0.68)
}

draw_numeric_heatmap <- function(mat,
                                 row_labels,
                                 col_labels,
                                 title,
                                 legend_label,
                                 diverging = TRUE,
                                 lim = NULL,
                                 digits = 2,
                                 na_label = "UNCERT") {
  mat <- as.matrix(mat)
  n_row <- nrow(mat)
  n_col <- ncol(mat)
  if (is.null(lim)) {
    vals <- mat[is.finite(mat)]
    lim <- if (length(vals)) max(abs(vals)) else 1
    if (!is.finite(lim) || lim <= 0) lim <- 1
  }
  pal <- if (diverging) grDevices::hcl.colors(101, "RdBu", rev = TRUE) else grDevices::hcl.colors(101, "Blues", rev = FALSE)
  graphics::plot(
    NA,
    xlim = c(0.5, n_col + 0.5),
    ylim = c(0.5, n_row + 0.5),
    xaxt = "n",
    yaxt = "n",
    xlab = "",
    ylab = "",
    main = title
  )
  graphics::axis(1, at = seq_len(n_col), labels = col_labels, las = 2, cex.axis = 0.90)
  graphics::axis(2, at = seq_len(n_row), labels = rev(row_labels), las = 2, cex.axis = 0.86)
  for (i in seq_len(n_row)) {
    for (j in seq_len(n_col)) {
      raw_i <- n_row - i + 1L
      value <- mat[raw_i, j]
      if (is.finite(value)) {
        if (diverging) {
          scaled <- floor((clip_signed(value, lim) + lim) / (2 * lim) * 100) + 1L
        } else {
          scaled <- floor(pmax(pmin(value / lim, 1), 0) * 100) + 1L
        }
        fill <- pal[scaled]
        text <- sprintf(paste0("%.", digits, "f"), value)
        text_col <- if (diverging && abs(value) > 0.65 * lim) "white" else "grey10"
      } else {
        fill <- "grey90"
        text <- na_label
        text_col <- "firebrick3"
      }
      graphics::rect(j - 0.5, i - 0.5, j + 0.5, i + 0.5, col = fill, border = "white")
      graphics::text(j, i, text, col = text_col, cex = 0.95, font = ifelse(is.finite(value), 1, 2))
    }
  }
  graphics::box()
  usr <- graphics::par("usr")
  x0 <- usr[2] - 0.055 * (usr[2] - usr[1])
  x1 <- usr[2] - 0.030 * (usr[2] - usr[1])
  y0 <- usr[3] + 0.13 * (usr[4] - usr[3])
  y1 <- usr[3] + 0.47 * (usr[4] - usr[3])
  y <- seq(y0, y1, length.out = length(pal) + 1L)
  graphics::rect(x0, y[-length(y)], x1, y[-1L], col = pal, border = NA, xpd = NA)
  graphics::rect(x0, y0, x1, y1, border = "grey30", xpd = NA)
  if (diverging) {
    graphics::text(x1 + 0.012 * (usr[2] - usr[1]), y1, sprintf(paste0("%.", digits, "f"), lim), adj = 0, cex = 0.74, xpd = NA)
    graphics::text(x1 + 0.012 * (usr[2] - usr[1]), (y0 + y1) / 2, "0", adj = 0, cex = 0.74, xpd = NA)
    graphics::text(x1 + 0.012 * (usr[2] - usr[1]), y0, sprintf(paste0("%.", digits, "f"), -lim), adj = 0, cex = 0.74, xpd = NA)
  } else {
    graphics::text(x1 + 0.012 * (usr[2] - usr[1]), y1, sprintf(paste0("%.", digits, "f"), lim), adj = 0, cex = 0.74, xpd = NA)
    graphics::text(x1 + 0.012 * (usr[2] - usr[1]), y0, "0", adj = 0, cex = 0.74, xpd = NA)
  }
  graphics::text((x0 + x1) / 2, y1 + 0.055 * (usr[4] - usr[3]), legend_label, cex = 0.74, xpd = NA)
}

weighted_density_curve <- function(x, w, from = NULL, to = NULL, adjust = 1.4) {
  x <- as.numeric(x)
  w <- as.numeric(w)
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (sum(ok) < 2L) return(NULL)
  x <- x[ok]
  w <- w[ok] / sum(w[ok])
  if (is.null(from)) from <- min(x) - 0.2 * diff(range(x))
  if (is.null(to)) to <- max(x) + 0.2 * diff(range(x))
  suppressWarnings(stats::density(x, weights = w, from = from, to = to, adjust = adjust, n = 256))
}

weighted_kde2d <- function(x, y, w, xlim, ylim, n = 90L) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  w <- as.numeric(w)
  ok <- is.finite(x) & is.finite(y) & is.finite(w) & w > 0
  if (!any(ok)) return(NULL)
  x <- x[ok]
  y <- y[ok]
  w <- w[ok] / sum(w[ok])
  gx <- seq(xlim[1L], xlim[2L], length.out = n)
  gy <- seq(ylim[1L], ylim[2L], length.out = n)
  bw_x <- max(stats::bw.nrd0(x), diff(xlim) / 20, 1e-4)
  bw_y <- max(stats::bw.nrd0(y), diff(ylim) / 20, 1e-4)
  z <- matrix(0, nrow = n, ncol = n)
  for (i in seq_along(w)) {
    z <- z + w[i] * outer(stats::dnorm(gx, x[i], bw_x), stats::dnorm(gy, y[i], bw_y))
  }
  list(x = gx, y = gy, z = z)
}

plot_promising_error_heat <- function() {
  x_name <- if ("mu_sv" %in% names(estimate_rows)) "mu_sv" else focus_hyper_names[1L]
  y_name <- if ("log_sigma2_sv" %in% names(estimate_rows)) "log_sigma2_sv" else focus_hyper_names[min(2L, length(focus_hyper_names))]
  grDevices::png(promising_heat_plot_file, width = 1500, height = 560 * length(local_positions))
  graphics::par(mfrow = c(length(local_positions), length(promising_methods)), mar = c(4.3, 4.7, 3.3, 1.1))
  pal <- grDevices::hcl.colors(101, "RdBu", rev = TRUE)
  for (local_pos in local_positions) {
    ref_l <- reference_summary[reference_summary$local_pos == local_pos, , drop = FALSE]
    ref_l <- ref_l[order(ref_l$theta_row), , drop = FALSE]
    ref_w <- normalize_log_weights(ref_l$reference_log_marginal)
    cex_by_ref <- 1.0 + 3.2 * sqrt(ref_w / max(ref_w, na.rm = TRUE))
    for (method in promising_methods) {
      df <- estimate_rows[
        estimate_rows$local_pos == local_pos & estimate_rows$method == method,
        ,
        drop = FALSE
      ]
      df <- df[order(df$theta_row), , drop = FALSE]
      lim <- max(0.5, signed_limit(df$error))
      clipped <- clip_signed(df$error, lim)
      col_idx <- floor((clipped + lim) / (2 * lim) * 100) + 1L
      col_idx[!is.finite(col_idx)] <- NA_integer_
      cols <- pal[col_idx]
      cols[is.na(cols)] <- "white"
      graphics::plot(
        df[[x_name]],
        df[[y_name]],
        xlab = x_name,
        ylab = y_name,
        main = sprintf("%s, %s\ncolor = log m error, size = nested posterior mass",
                       unique(df$local), method_label[[method]]),
        pch = ifelse(df$certified, 21, 4),
        bg = cols,
        col = ifelse(df$certified, "grey20", "firebrick3"),
        cex = cex_by_ref,
        lwd = 2
      )
      graphics::text(df[[x_name]], df[[y_name]], labels = df$theta_row, pos = 3, cex = 0.78)
      graphics::legend(
        "topleft",
        legend = c("point size = nested theta mass", "uncertified = red x"),
        bty = "n",
        cex = 0.82
      )
      add_signed_color_legend(lim, pal, "log m error")
    }
  }
  grDevices::dev.off()
}

plot_posterior_overlap <- function() {
  detail <- posterior_diagnostics$detail
  summary <- posterior_diagnostics$summary
  blocks <- unique(detail$block_key)
  grDevices::png(posterior_overlap_plot_file, width = 1600, height = 520 * length(blocks))
  graphics::par(mfrow = c(length(blocks), 1L), mar = c(4.7, 5.0, 3.3, 1.2))
  offsets <- c(query_native = -0.08, query_native_repaired = 0.08)
  for (block in blocks) {
    df_b <- detail[detail$block_key == block, , drop = FALSE]
    ref <- df_b[df_b$method == promising_methods[1L], , drop = FALSE]
    ref <- ref[order(ref$theta_row), , drop = FALSE]
    y_max <- max(c(ref$reference_weight, df_b$method_weight), na.rm = TRUE)
    graphics::plot(
      ref$theta_row,
      ref$reference_weight,
      type = "b",
      pch = 16,
      lwd = 2,
      col = "grey15",
      ylim = c(0, y_max * 1.15),
      xlab = "theta profile",
      ylab = "normalized discrete theta mass",
      main = sprintf(
        "%s: expensive nested-SMC posterior over audit theta profiles",
        ref$block_label[1L]
      )
    )
    for (method in promising_methods) {
      df <- df_b[df_b$method == method, , drop = FALSE]
      df <- df[order(df$theta_row), , drop = FALSE]
      graphics::lines(
        df$theta_row + offsets[[method]],
        df$method_weight,
        type = "b",
        pch = ifelse(df$certified, 17, 4),
        lwd = 2,
        col = point_col[[method]]
      )
    }
    block_summary <- summary[summary$block_key == block, , drop = FALSE]
    legend_text <- c(
      "nested reference",
      sprintf(
        "%s overlap %.3f, TV %.3f, ref mass uncert %.3f",
        method_label[block_summary$method],
        block_summary$overlap,
        block_summary$total_variation,
        block_summary$reference_mass_uncertified
      )
    )
    graphics::legend(
      "topright",
      legend = legend_text,
      pch = c(16, 17, 17),
      lty = 1,
      col = c("grey15", point_col[block_summary$method]),
      bty = "n",
      cex = 0.82
    )
  }
  grDevices::dev.off()
}

plot_posterior_mass_heat <- function() {
  detail <- posterior_diagnostics$detail
  x_name <- if ("mu_sv" %in% names(detail)) "mu_sv" else focus_hyper_names[1L]
  y_name <- if ("log_sigma2_sv" %in% names(detail)) "log_sigma2_sv" else focus_hyper_names[min(2L, length(focus_hyper_names))]
  blocks <- unique(detail$block_key)
  grDevices::png(posterior_heat_plot_file, width = 1500, height = 560 * length(blocks))
  graphics::par(mfrow = c(length(blocks), length(promising_methods)), mar = c(4.3, 4.7, 3.3, 1.1))
  pal <- grDevices::hcl.colors(101, "RdBu", rev = TRUE)
  for (block in blocks) {
    for (method in promising_methods) {
      df <- detail[detail$block_key == block & detail$method == method, , drop = FALSE]
      lim <- max(0.05, signed_limit(df$weight_diff))
      clipped <- clip_signed(df$weight_diff, lim)
      col_idx <- floor((clipped + lim) / (2 * lim) * 100) + 1L
      col_idx[!is.finite(col_idx)] <- NA_integer_
      cols <- pal[col_idx]
      cols[is.na(cols)] <- "white"
      cex_by_ref <- 1.0 + 3.2 * sqrt(df$reference_weight / max(df$reference_weight, na.rm = TRUE))
      graphics::plot(
        df[[x_name]],
        df[[y_name]],
        xlab = x_name,
        ylab = y_name,
        main = sprintf("%s, %s\ncolor = posterior mass diff, size = nested mass",
                       df$block_label[1L], method_label[[method]]),
        pch = ifelse(df$certified, 21, 4),
        bg = cols,
        col = ifelse(df$certified, "grey20", "firebrick3"),
        cex = cex_by_ref,
        lwd = 2
      )
      graphics::text(df[[x_name]], df[[y_name]], labels = df$theta_row, pos = 3, cex = 0.78)
      graphics::legend(
        "topleft",
        legend = c("point size = nested theta mass", "uncertified = red x"),
        bty = "n",
        cex = 0.82
      )
      add_signed_color_legend(lim, pal, "mass diff", digits = 3)
    }
  }
  grDevices::dev.off()
}

plot_mll_error_heatmap <- function() {
  grDevices::png(mll_error_heatmap_file, width = 1900, height = 680 * length(local_positions))
  graphics::par(mfrow = c(length(local_positions), 1L), mar = c(5.8, 11.5, 3.4, 5.0))
  col_labels <- c("particle-MIS", "repaired")
  for (local_pos in local_positions) {
    ref_l <- reference_summary[reference_summary$local_pos == local_pos, , drop = FALSE]
    ref_l <- ref_l[order(ref_l$theta_row), , drop = FALSE]
    ref_w <- normalize_log_weights(ref_l$reference_log_marginal)
    rows <- ref_l$theta_row
    mat <- matrix(NA_real_, nrow = length(rows), ncol = length(promising_methods))
    rownames(mat) <- rows
    colnames(mat) <- promising_methods
    for (j in seq_along(promising_methods)) {
      method <- promising_methods[j]
      df <- estimate_rows[
        estimate_rows$local_pos == local_pos & estimate_rows$method == method,
        ,
        drop = FALSE
      ]
      mat[, j] <- df$error[match(rows, df$theta_row)]
    }
    row_labels <- sprintf(
      "t%02d | mass %.1f%% | sd %.2f",
      rows,
      100 * ref_w,
      ref_l$reference_rep_sd
    )
    lim <- max(0.5, signed_limit(mat))
    draw_numeric_heatmap(
      mat = mat,
      row_labels = row_labels,
      col_labels = col_labels,
      title = sprintf(
        "EMC local %d/%s: MLL error vs expensive nested-SMC mean",
        local_pos,
        unique(ref_l$local)[1L]
      ),
      legend_label = "log m error",
      diverging = TRUE,
      lim = lim,
      digits = 2
    )
  }
  grDevices::dev.off()
}

plot_mll_error_density <- function() {
  ref <- merge(
    reference_raw,
    reference_summary,
    by = c("local_pos", "local", "theta_row", "theta_label"),
    all.x = TRUE,
    suffixes = c("_rep", "_mean")
  )
  ref$replicate_error <- ref$reference_log_marginal_rep - ref$reference_log_marginal_mean
  est <- estimate_rows[estimate_rows$method %in% promising_methods, , drop = FALSE]
  ref_replicates <- reference_raw[, c("local_pos", "local", "theta_row", "theta_label", "replicate", "reference_log_marginal"), drop = FALSE]
  names(ref_replicates)[names(ref_replicates) == "reference_log_marginal"] <- "nested_replicate_log_marginal"
  rep_errors <- merge(
    est,
    ref_replicates,
    by = c("local_pos", "local", "theta_row", "theta_label"),
    all.x = TRUE,
    sort = FALSE
  )
  rep_errors$replicate_error <- rep_errors$estimate - rep_errors$nested_replicate_log_marginal

  grDevices::png(mll_error_density_file, width = 1500, height = 1100)
  graphics::layout(matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE))
  graphics::par(mar = c(4.8, 5.0, 3.2, 1.2))

  finite_all <- c(ref$replicate_error, rep_errors$replicate_error)
  finite_all <- finite_all[is.finite(finite_all)]
  xlim <- range(finite_all)
  xlim <- xlim + c(-0.08, 0.08) * diff(xlim)
  graphics::plot(
    NA,
    xlim = xlim,
    ylim = c(0, 3.5),
    xlab = "log m estimate - nested-SMC replicate",
    ylab = "density",
    main = "MLL error distribution across nested-SMC reference replicates"
  )
  graphics::abline(v = 0, lty = 2, col = "grey45")
  d_ref <- stats::density(ref$replicate_error[is.finite(ref$replicate_error)], from = xlim[1L], to = xlim[2L], adjust = 1.2)
  graphics::lines(d_ref$x, d_ref$y, col = "grey25", lwd = 2)
  for (method in promising_methods) {
    vals <- rep_errors$replicate_error[rep_errors$method == method]
    vals <- vals[is.finite(vals)]
    if (length(vals) > 2L) {
      d <- stats::density(vals, from = xlim[1L], to = xlim[2L], adjust = 1.2)
      graphics::lines(d$x, d$y, col = point_col[[method]], lwd = 2)
    }
  }
  graphics::legend(
    "topright",
    legend = c("nested replicate noise", method_label[promising_methods]),
    col = c("grey25", point_col[promising_methods]),
    lwd = 2,
    bty = "n"
  )

  abs_vals <- list(
    "nested replicate noise" = abs(ref$replicate_error[is.finite(ref$replicate_error)])
  )
  for (method in promising_methods) {
    vals <- abs(rep_errors$replicate_error[rep_errors$method == method])
    abs_vals[[method_label[[method]]]] <- vals[is.finite(vals)]
  }
  graphics::boxplot(
    abs_vals,
    las = 2,
    ylab = "absolute log m error vs nested replicate",
    main = "Absolute MLL error across reference runs",
    col = c("grey80", point_col[promising_methods])
  )
  graphics::abline(h = stats::median(abs_vals[[1L]], na.rm = TRUE), lty = 2, col = "grey45")

  local_method <- rep_errors[is.finite(rep_errors$replicate_error), , drop = FALSE]
  local_method$group <- paste0("L", local_method$local, "\n", method_label[local_method$method])
  graphics::boxplot(
    abs(replicate_error) ~ group,
    data = local_method,
    las = 2,
    ylab = "absolute log m error",
    main = "Which EMC local is still expensive/error-prone?",
    col = rep(point_col[promising_methods], each = length(local_positions))
  )

  summary_for_bar <- do.call(rbind, lapply(promising_methods, function(method) {
    vals <- rep_errors$replicate_error[rep_errors$method == method]
    vals <- vals[is.finite(vals)]
    data.frame(
      method = method,
      rmse = sqrt(mean(vals^2)),
      median_abs = stats::median(abs(vals)),
      q90_abs = as.numeric(stats::quantile(abs(vals), 0.9)),
      check.names = FALSE
    )
  }))
  bar_mat <- t(as.matrix(summary_for_bar[, c("rmse", "median_abs", "q90_abs"), drop = FALSE]))
  colnames(bar_mat) <- method_label[summary_for_bar$method]
  graphics::barplot(
    bar_mat,
    beside = TRUE,
    las = 2,
    col = c("grey35", "steelblue3", "tan3"),
    ylab = "log m error",
    main = "MLL error magnitude summary"
  )
  graphics::legend("topright", legend = rownames(bar_mat), fill = c("grey35", "steelblue3", "tan3"), bty = "n")
  grDevices::dev.off()
}

posterior_long_weights <- function() {
  detail <- posterior_diagnostics$detail
  blocks <- unique(detail$block_key)
  out <- list()
  for (block in blocks) {
    ref <- detail[detail$block_key == block & detail$method == promising_methods[1L], , drop = FALSE]
    ref$method_display <- "nested reference"
    ref$plot_weight <- ref$reference_weight
    ref$plot_certified <- TRUE
    out[[length(out) + 1L]] <- ref
    for (method in promising_methods) {
      df <- detail[detail$block_key == block & detail$method == method, , drop = FALSE]
      df$method_display <- method_label[[method]]
      df$plot_weight <- df$method_weight
      df$plot_certified <- df$certified
      out[[length(out) + 1L]] <- df
    }
  }
  do.call(rbind, out)
}

plot_posterior_weight_heatmap <- function() {
  long <- posterior_long_weights()
  blocks <- unique(long$block_key)
  display_order <- c("nested reference", method_label[promising_methods])
  grDevices::png(posterior_weight_heatmap_file, width = 1500, height = 600 * length(blocks))
  graphics::par(mfrow = c(length(blocks), 1L), mar = c(5.2, 9.2, 3.2, 4.8))
  for (block in blocks) {
    df_b <- long[long$block_key == block, , drop = FALSE]
    theta_rows <- sort(unique(df_b$theta_row))
    mat <- matrix(NA_real_, nrow = length(display_order), ncol = length(theta_rows))
    uncert <- matrix(FALSE, nrow = length(display_order), ncol = length(theta_rows))
    rownames(mat) <- display_order
    colnames(mat) <- theta_rows
    for (i in seq_along(display_order)) {
      df <- df_b[df_b$method_display == display_order[i], , drop = FALSE]
      mat[i, ] <- df$plot_weight[match(theta_rows, df$theta_row)]
      uncert[i, ] <- !(df$plot_certified[match(theta_rows, df$theta_row)] %in% TRUE)
    }
    percent_mat <- 100 * mat
    draw_numeric_heatmap(
      mat = percent_mat,
      row_labels = display_order,
      col_labels = sprintf("theta %02d", theta_rows),
      title = sprintf("%s: posterior mass over audit theta profiles", unique(df_b$block_label)[1L]),
      legend_label = "mass %",
      diverging = FALSE,
      lim = max(percent_mat, na.rm = TRUE),
      digits = 1,
      na_label = "NA"
    )
    for (i in seq_len(nrow(uncert))) {
      for (j in seq_len(ncol(uncert))) {
        if (uncert[i, j]) {
          graphics::points(j, nrow(uncert) - i + 1L, pch = 4, col = "firebrick3", cex = 2.0, lwd = 2)
        }
      }
    }
  }
  grDevices::dev.off()
}

plot_posterior_density_1d <- function() {
  long <- posterior_long_weights()
  vars <- intersect(c("mu_sv", "log_sigma2_sv"), names(long))
  if (length(vars) < 2L) vars <- focus_hyper_names[seq_len(min(2L, length(focus_hyper_names)))]
  blocks <- unique(long$block_key)
  display_order <- c("nested reference", method_label[promising_methods])
  display_col <- c("nested reference" = "grey15", setNames(point_col[promising_methods], method_label[promising_methods]))
  grDevices::png(posterior_density_1d_file, width = 1500, height = 520 * length(blocks))
  graphics::par(mfrow = c(length(blocks), length(vars)), mar = c(4.5, 4.8, 3.2, 1.2))
  for (block in blocks) {
    df_b <- long[long$block_key == block, , drop = FALSE]
    for (var in vars) {
      x_vals <- df_b[[var]]
      x_rng <- range(x_vals, na.rm = TRUE)
      pad <- 0.22 * diff(x_rng)
      from <- x_rng[1L] - pad
      to <- x_rng[2L] + pad
      curves <- lapply(display_order, function(label) {
        df <- df_b[df_b$method_display == label, , drop = FALSE]
        weighted_density_curve(df[[var]], df$plot_weight, from = from, to = to, adjust = 1.5)
      })
      names(curves) <- display_order
      y_max <- max(vapply(curves, function(d) if (is.null(d)) 0 else max(d$y), numeric(1)), na.rm = TRUE)
      graphics::plot(
        NA,
        xlim = c(from, to),
        ylim = c(0, y_max * 1.08),
        xlab = var,
        ylab = "weighted KDE density",
        main = sprintf("%s: posterior density for %s", unique(df_b$block_label)[1L], var)
      )
      for (label in display_order) {
        d <- curves[[label]]
        if (!is.null(d)) graphics::lines(d$x, d$y, col = display_col[[label]], lwd = 2)
      }
      graphics::rug(df_b[[var]][df_b$method_display == "nested reference"], col = "grey50")
      graphics::legend("topright", legend = display_order, col = display_col[display_order], lwd = 2, bty = "n", cex = 0.8)
    }
  }
  grDevices::dev.off()
}

plot_posterior_density_2d <- function() {
  long <- posterior_long_weights()
  x_name <- if ("mu_sv" %in% names(long)) "mu_sv" else focus_hyper_names[1L]
  y_name <- if ("log_sigma2_sv" %in% names(long)) "log_sigma2_sv" else focus_hyper_names[min(2L, length(focus_hyper_names))]
  blocks <- unique(long$block_key)
  display_order <- c("nested reference", method_label[promising_methods])
  grDevices::png(posterior_density_2d_file, width = 1600, height = 520 * length(blocks))
  graphics::par(mfrow = c(length(blocks), length(display_order)), mar = c(4.4, 4.6, 3.2, 1.1))
  x_rng <- range(long[[x_name]], na.rm = TRUE)
  y_rng <- range(long[[y_name]], na.rm = TRUE)
  x_pad <- 0.18 * diff(x_rng)
  y_pad <- 0.18 * diff(y_rng)
  xlim <- x_rng + c(-x_pad, x_pad)
  ylim <- y_rng + c(-y_pad, y_pad)
  for (block in blocks) {
    df_b <- long[long$block_key == block, , drop = FALSE]
    dens <- lapply(display_order, function(label) {
      df <- df_b[df_b$method_display == label, , drop = FALSE]
      weighted_kde2d(df[[x_name]], df[[y_name]], df$plot_weight, xlim = xlim, ylim = ylim, n = 90L)
    })
    names(dens) <- display_order
    z_lim <- max(vapply(dens, function(d) if (is.null(d)) 0 else max(d$z), numeric(1)), na.rm = TRUE)
    for (label in display_order) {
      d <- dens[[label]]
      if (is.null(d)) {
        graphics::plot.new()
        graphics::title(sprintf("%s\n%s", unique(df_b$block_label)[1L], label))
        next
      }
      graphics::image(
        d$x,
        d$y,
        d$z,
        xlim = xlim,
        ylim = ylim,
        zlim = c(0, z_lim),
        col = grDevices::hcl.colors(80, "YlOrRd", rev = FALSE),
        xlab = x_name,
        ylab = y_name,
        main = sprintf("%s\n%s density", unique(df_b$block_label)[1L], label)
      )
      graphics::contour(d$x, d$y, d$z, add = TRUE, drawlabels = FALSE, col = grDevices::adjustcolor("grey20", 0.55))
      df <- df_b[df_b$method_display == label, , drop = FALSE]
      point_cex <- 0.9 + 2.5 * sqrt(df$plot_weight / max(df$plot_weight, na.rm = TRUE))
      graphics::points(df[[x_name]], df[[y_name]], pch = 21, bg = "white", col = "grey15", cex = point_cex, lwd = 1.2)
      graphics::text(df[[x_name]], df[[y_name]], labels = df$theta_row, pos = 3, cex = 0.65)
    }
  }
  grDevices::dev.off()
}

plot_reference_replicates()
plot_theta_heat()
plot_repair_effect()
plot_critical_failures()
plot_promising_error_heat()
plot_posterior_overlap()
plot_posterior_mass_heat()
plot_mll_error_heatmap()
plot_mll_error_density()
plot_posterior_weight_heatmap()
plot_posterior_density_1d()
plot_posterior_density_2d()

saveRDS(
  list(
    estimate_rows = estimate_rows,
    method_summary = method_summary,
    promising_method_summary = promising_method_summary,
    reference_raw = reference_raw,
    reference_summary = reference_summary,
    posterior_diagnostics = posterior_diagnostics,
    repair_rows = repair_rows,
    atlas_objects = atlas_objects,
    repair_objects = repair_objects,
    atlas_design = atlas_design,
    audit_design = audit_design,
    settings = list(
      label = label,
      seed = seed,
      local_positions = local_positions,
      local_ids = local_ids,
      atlas_particles = atlas_particles,
      atlas_max_anchors = atlas_max_anchors,
      ref_particles = ref_particles,
      ref_reps = ref_reps,
      repair_particles = repair_particles,
      repair_confirmation_reps = repair_confirmation_reps,
      repair_adaptive_confirmation_reps = repair_adaptive_confirmation_reps,
      repair_replicate_bootstrap_B = repair_replicate_bootstrap_B,
      repair_max_direct_graph_z = repair_max_direct_graph_z,
      repair_max_direct_graph_chart_shift = repair_max_direct_graph_chart_shift,
      repair_max_direct_graph_existing_shift = repair_max_direct_graph_existing_shift,
      repair_adaptive_replicate_weight_multiplier = repair_adaptive_replicate_weight_multiplier,
      repair_adaptive_replicate_min_theta_weight = repair_adaptive_replicate_min_theta_weight,
      repair_adaptive_replicate_max_graph_z = repair_adaptive_replicate_max_graph_z,
      repair_adaptive_replicate_graph_shift = repair_adaptive_replicate_graph_shift,
      normalizer_robust_method = normalizer_robust_method,
      normalizer_student_t_df = normalizer_student_t_df,
      distance_metric = distance_metric,
      max_chart_distance = max_chart_distance,
      min_covering_charts = min_covering_charts,
      min_particle_mis_ess = min_particle_mis_ess
    ),
    files = list(
      rows_csv = rows_csv,
      summary_csv = summary_csv,
      promising_summary_csv = promising_summary_csv,
      posterior_summary_csv = posterior_summary_csv,
      reference_cache_file = reference_cache_file,
      error_plot_file = error_plot_file,
      heat_plot_file = heat_plot_file,
      repair_plot_file = repair_plot_file,
      critical_plot_file = critical_plot_file,
      promising_heat_plot_file = promising_heat_plot_file,
      posterior_overlap_plot_file = posterior_overlap_plot_file,
      posterior_heat_plot_file = posterior_heat_plot_file,
      mll_error_heatmap_file = mll_error_heatmap_file,
      mll_error_density_file = mll_error_density_file,
      posterior_weight_heatmap_file = posterior_weight_heatmap_file,
      posterior_density_1d_file = posterior_density_1d_file,
      posterior_density_2d_file = posterior_density_2d_file
    )
  ),
  results_file
)

cat(sprintf("Saved results: %s\n", results_file))
cat(sprintf("Saved rows: %s\n", rows_csv))
cat(sprintf("Saved summary: %s\n", summary_csv))
cat(sprintf("Saved promising summary: %s\n", promising_summary_csv))
cat(sprintf("Saved posterior overlap summary: %s\n", posterior_summary_csv))
cat(sprintf("Saved nested error plot: %s\n", error_plot_file))
cat(sprintf("Saved theta error heat plot: %s\n", heat_plot_file))
cat(sprintf("Saved repair effect plot: %s\n", repair_plot_file))
cat(sprintf("Saved critical failure plot: %s\n", critical_plot_file))
cat(sprintf("Saved promising error heat plot: %s\n", promising_heat_plot_file))
cat(sprintf("Saved posterior overlap plot: %s\n", posterior_overlap_plot_file))
cat(sprintf("Saved posterior mass heat plot: %s\n", posterior_heat_plot_file))
cat(sprintf("Saved MLL error heatmap: %s\n", mll_error_heatmap_file))
cat(sprintf("Saved MLL error density plot: %s\n", mll_error_density_file))
cat(sprintf("Saved posterior weight heatmap: %s\n", posterior_weight_heatmap_file))
cat(sprintf("Saved posterior 1D density plot: %s\n", posterior_density_1d_file))
cat(sprintf("Saved posterior 2D density plot: %s\n", posterior_density_2d_file))
cat("\nMethod summary:\n")
print(method_summary, row.names = FALSE)
cat("\nPromising-method posterior overlap summary:\n")
print(posterior_diagnostics$summary, row.names = FALSE)
