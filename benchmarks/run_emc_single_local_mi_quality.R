#!/usr/bin/env Rscript

rm(list = ls())

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[1L]))
} else {
  normalizePath("benchmarks/run_emc_single_local_mi_quality.R")
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

finite_range_width <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x)) diff(range(x)) else NA_real_
}

source("smc_core.R")
source("population_models.R")
source("local_charts.R")
source("utilities.R")

cli_args <- parse_cli_args(commandArgs(trailingOnly = TRUE))
detected_cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
if (!is.finite(detected_cores) || detected_cores < 1L) detected_cores <- 1L

label <- arg_chr(cli_args, "label", "emc_single_local_mi_quality")
seed <- arg_int(cli_args, "seed", 20260524L)
cores <- arg_int(cli_args, "cores", min(4L, detected_cores))
data_file <- arg_chr(cli_args, "data_file", file.path("benchmarks", "samples", "full_EMC2.RData"))
local_pos_arg <- arg_chr(cli_args, "local_pos", "14")
atlas_particles <- arg_int(cli_args, "atlas_particles", 512L)
atlas_max_anchors <- arg_int(cli_args, "atlas_max_anchors", 18L)
atlas_max_steps <- arg_int(cli_args, "atlas_max_steps", 176L)
atlas_max_intermediates <- arg_int(cli_args, "atlas_max_intermediates", 4L)
atlas_target_cess <- arg_num(cli_args, "atlas_target_cess", 0.9)
atlas_mcmc_moves <- arg_int(cli_args, "atlas_mcmc_moves", 2L)
ref_particles <- arg_int(cli_args, "ref_particles", 768L)
ref_reps <- arg_int(cli_args, "ref_reps", 2L)
ref_max_steps <- arg_int(cli_args, "ref_max_steps", 176L)
ref_target_cess <- arg_num(cli_args, "ref_target_cess", 0.9)
ref_mcmc_moves <- arg_int(cli_args, "ref_mcmc_moves", 2L)
audit_max_points <- arg_int(cli_args, "audit_max_points", 13L)
atlas_design_probs <- arg_num_vec(cli_args, "atlas_design_probs", c(0.05, 0.25, 0.5, 0.75, 0.95))
audit_probs <- arg_num_vec(cli_args, "audit_probs", c(0.01, 0.10, 0.33, 0.50, 0.67, 0.90, 0.99))
audit_inflation <- arg_num_vec(cli_args, "audit_inflation", 1)
focus_parameters <- arg_chr_vec(
  cli_args,
  "focus_parameters",
  c("mu_sv", "sigma2_sv", "sigma2_v_LogFreq", "sigma2_v")
)
distance_metric <- arg_chr(cli_args, "distance_metric", "fisher")
surface_method <- arg_chr(cli_args, "surface_method", "derivative_ls")
particle_mis_role <- arg_chr(cli_args, "particle_mis_role", "estimator")
max_chart_distance <- arg_num(cli_args, "max_chart_distance", 1.0)
min_covering_charts <- arg_int(cli_args, "min_covering_charts", 3L)
max_quadratic_particle_gap <- arg_num(cli_args, "max_quadratic_particle_gap", 0.02)
min_particle_mis_ess <- arg_num(cli_args, "min_particle_mis_ess", 0.05)
max_particle_mis_psis_k <- arg_num(cli_args, "max_particle_mis_psis_k", 0.7)
sparse_chart_min_covering <- arg_int(cli_args, "sparse_chart_min_covering", 3L)
sparse_chart_max_distance <- arg_num(cli_args, "sparse_chart_max_distance", 0.05)
skip_reference <- arg_lgl(cli_args, "skip_reference", FALSE)

results_file <- arg_chr(
  cli_args,
  "results_file",
  file.path("benchmarks", "results", paste0(label, "_results.rds"))
)
rows_csv <- arg_chr(
  cli_args,
  "rows_csv",
  file.path("benchmarks", "results", paste0(label, "_rows.csv"))
)
summary_csv <- arg_chr(
  cli_args,
  "summary_csv",
  file.path("benchmarks", "results", paste0(label, "_summary.csv"))
)
plot_file <- arg_chr(
  cli_args,
  "plot_file",
  file.path("benchmarks", "results", paste0(label, "_diagnostics.png"))
)
miss_plot_file <- arg_chr(
  cli_args,
  "miss_plot_file",
  file.path("benchmarks", "results", paste0(label, "_miss_diagnostics.png"))
)
failure_report_file <- arg_chr(
  cli_args,
  "failure_report_file",
  file.path("benchmarks", "results", paste0(label, "_failure_report.png"))
)
comparison_plot_file <- arg_chr(
  cli_args,
  "comparison_plot_file",
  file.path("benchmarks", "results", paste0(label, "_comparison_clear.png"))
)
coverage_plot_file <- arg_chr(
  cli_args,
  "coverage_plot_file",
  file.path("benchmarks", "results", paste0(label, "_coverage_clear.png"))
)
estimator_plot_file <- arg_chr(
  cli_args,
  "estimator_plot_file",
  file.path("benchmarks", "results", paste0(label, "_estimator_clear.png"))
)
reference_plot_file <- arg_chr(
  cli_args,
  "reference_plot_file",
  file.path("benchmarks", "results", paste0(label, "_reference_replicates.png"))
)
miss_csv <- arg_chr(
  cli_args,
  "miss_csv",
  file.path("benchmarks", "results", paste0(label, "_misses.csv"))
)
reference_cache_file <- arg_chr(
  cli_args,
  "reference_cache_file",
  file.path("benchmarks", "results", paste0(label, "_nested_reference_cache.rds"))
)
reuse_reference_cache <- arg_lgl(cli_args, "reuse_reference_cache", TRUE)
refresh_reference_cache <- arg_lgl(cli_args, "refresh_reference_cache", FALSE)

dir.create(dirname(results_file), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(rows_csv), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(summary_csv), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(plot_file), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(miss_plot_file), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(failure_report_file), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(comparison_plot_file), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(coverage_plot_file), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(estimator_plot_file), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(reference_plot_file), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(miss_csv), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(reference_cache_file), showWarnings = FALSE, recursive = TRUE)

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

if (identical(tolower(local_pos_arg), "auto")) {
  n_obs <- vapply(data_list, nrow, integer(1))
  local_pos <- which.max(n_obs)
} else {
  local_pos <- as.integer(local_pos_arg)
}
if (!is.finite(local_pos) || local_pos < 1L || local_pos > length(data_list)) {
  stop("local_pos must be an integer in 1:", length(data_list), " or 'auto'.")
}
local_id <- names(data_list)[local_pos] %||% paste0("local_", local_pos)
data_i <- data_list[[local_pos]]

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
  source = "emc_single_local_atlas_design",
  label_prefix = "atlas"
)

audit_design <- build_local_evidence_calibration_design(
  theta = theta_cloud,
  population_model = population_model,
  weights = rep(1 / nrow(theta_cloud), nrow(theta_cloud)),
  focus_hyper_names = focus_hyper_names,
  probs = audit_probs,
  max_points = audit_max_points,
  inflation = audit_inflation,
  n_leverage_points = 1L,
  n_uncertainty_points = 0L,
  n_disagreement_points = 0L,
  source = "emc_single_local_audit_design",
  label_prefix = "audit"
)

cat(sprintf(
  "EMC single-local m_i(theta): local_pos=%d local=%s n_obs=%d alpha_dim=%d\n",
  local_pos, local_id, nrow(data_i), alpha_dim
))
cat(sprintf(
  "Atlas: %d anchors x %d particles | reference: %d theta x %d reps x %d particles\n",
  nrow(atlas_design$theta), atlas_particles, nrow(audit_design$theta), ref_reps, ref_particles
))
cat("Focus hyperparameters:", paste(focus_hyper_names, collapse = ", "), "\n")

start_time <- Sys.time()
atlas <- build_local_atlas(
  local_id = local_id,
  data_i = data_i,
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
  seed = seed + 1000L,
  verbose = FALSE
)

eval_list <- lapply(seq_len(nrow(audit_design$theta)), function(theta_row) {
  evaluate_local_atlas(
    atlas = atlas,
    theta = audit_design$theta[theta_row, , drop = FALSE],
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
    max_quadratic_particle_gap = max_quadratic_particle_gap,
    sparse_chart_min_covering = sparse_chart_min_covering,
    sparse_chart_max_distance = sparse_chart_max_distance,
    max_leave_chart_out_gap = Inf,
    distance_metric = distance_metric,
    surface_method = surface_method,
    min_surface_charts = 2L,
    max_surface_se = Inf,
    surface_value_nugget = 0.05,
    surface_gradient_weight = 1.0,
    surface_curvature_weight = 0.2,
    surface_ridge = 1e-8,
    particle_mis_role = particle_mis_role
  )
})

atlas_rows <- data.frame(
  theta_row = seq_len(nrow(audit_design$theta)),
  theta_label = audit_design$metadata$theta_label,
  atlas_log_marginal = vapply(eval_list, `[[`, numeric(1), "log_marginal"),
  atlas_se = vapply(eval_list, `[[`, numeric(1), "se"),
  atlas_status = vapply(eval_list, `[[`, character(1), "status"),
  atlas_reason = vapply(eval_list, `[[`, character(1), "reason"),
  nearest_charts = vapply(eval_list, function(x) paste(x$nearest_charts, collapse = ","), character(1)),
  n_covering_charts = vapply(eval_list, function(x) {
    rows <- x$diagnostics$predictions %||% data.frame(within_distance = logical())
    if ("within_distance" %in% names(rows)) sum(rows$within_distance %in% TRUE) else nrow(rows)
  }, integer(1)),
  particle_mis_ess_frac = vapply(eval_list, function(x) {
    as.numeric(x$diagnostics$particle_mis$ess_frac %||% NA_real_)
  }, numeric(1)),
  particle_mis_psis_k = vapply(eval_list, function(x) {
    as.numeric(x$diagnostics$particle_mis$psis_k %||% NA_real_)
  }, numeric(1)),
  quadratic_particle_gap = vapply(eval_list, function(x) {
    as.numeric(x$diagnostics$quadratic_particle_gap %||% NA_real_)
  }, numeric(1)),
  surface_log_marginal = vapply(eval_list, function(x) {
    as.numeric(x$diagnostics$surface_log_marginal %||% NA_real_)
  }, numeric(1)),
  particle_mis_log_marginal = vapply(eval_list, function(x) {
    as.numeric(x$diagnostics$particle_mis$log_marginal %||% NA_real_)
  }, numeric(1)),
  surface_se = vapply(eval_list, function(x) {
    as.numeric(x$diagnostics$surface_se %||% NA_real_)
  }, numeric(1)),
  surface_residual_sd = vapply(eval_list, function(x) {
    as.numeric(x$diagnostics$surface_residual_sd %||% NA_real_)
  }, numeric(1)),
  surface_n_charts = vapply(eval_list, function(x) {
    as.integer(x$diagnostics$surface_n_charts %||% NA_integer_)
  }, integer(1)),
  surface_n_observations = vapply(eval_list, function(x) {
    as.integer(x$diagnostics$surface_n_observations %||% NA_integer_)
  }, integer(1)),
  min_covering_distance = vapply(eval_list, function(x) {
    as.numeric(x$diagnostics$min_covering_distance %||% NA_real_)
  }, numeric(1)),
  check.names = FALSE
)
atlas_rows <- cbind(atlas_rows, audit_design$metadata[, population_model$hyper_names, drop = FALSE])

run_reference_one <- function(theta_row, rep_id) {
  run <- .local_chart_run_smc(
    local_id = local_id,
    theta_anchor = audit_design$theta[theta_row, , drop = FALSE],
    data_i = data_i,
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
    seed = seed + 1000003L * rep_id + 9176L * theta_row,
    verbose = FALSE,
    source = "emc_single_local_reference"
  )
  data.frame(
    theta_row = theta_row,
    theta_label = audit_design$metadata$theta_label[theta_row],
    replicate = rep_id,
    reference_log_marginal = run$logZ,
    reference_logZ_se = run$logZ_se,
    reference_final_ess_frac = run$diagnostics$final_ess_frac %||% NA_real_,
    reference_min_path_ess_frac = run$diagnostics$min_path_ess_frac %||% NA_real_,
    check.names = FALSE
  )
}

reference_cache_settings <- list(
  data_file = normalizePath(data_file),
  local_pos = as.integer(local_pos),
  local = as.character(local_id),
  alpha_names = alpha_names,
  ref_particles = as.integer(ref_particles),
  ref_reps = as.integer(ref_reps),
  ref_max_steps = as.integer(ref_max_steps),
  ref_target_cess = as.numeric(ref_target_cess),
  ref_mcmc_moves = as.integer(ref_mcmc_moves),
  seed = as.integer(seed),
  audit_max_points = as.integer(audit_max_points),
  audit_probs = as.numeric(audit_probs),
  audit_inflation = as.numeric(audit_inflation),
  focus_hyper_names = focus_hyper_names
)

reference_cache_compatible <- function(cache) {
  if (!is.list(cache) ||
      is.null(cache$reference_raw) ||
      is.null(cache$reference_rows) ||
      is.null(cache$audit_theta) ||
      is.null(cache$settings)) {
    return(FALSE)
  }
  settings <- cache$settings
  same_scalar <- function(name) identical(settings[[name]], reference_cache_settings[[name]])
  scalar_names <- c("data_file", "local_pos", "local", "ref_particles", "ref_reps",
                    "ref_max_steps", "ref_mcmc_moves", "seed", "audit_max_points")
  if (!all(vapply(scalar_names, same_scalar, logical(1)))) {
    return(FALSE)
  }
  if (!isTRUE(all.equal(settings$ref_target_cess, reference_cache_settings$ref_target_cess, tolerance = 1e-12))) {
    return(FALSE)
  }
  if (!identical(as.character(settings$focus_hyper_names), as.character(reference_cache_settings$focus_hyper_names))) {
    return(FALSE)
  }
  if (!isTRUE(all.equal(as.numeric(settings$audit_probs), as.numeric(reference_cache_settings$audit_probs), tolerance = 1e-12))) {
    return(FALSE)
  }
  if (!isTRUE(all.equal(as.numeric(settings$audit_inflation), as.numeric(reference_cache_settings$audit_inflation), tolerance = 1e-12))) {
    return(FALSE)
  }
  isTRUE(all.equal(
    as.matrix(cache$audit_theta),
    as.matrix(audit_design$theta),
    tolerance = 1e-10,
    check.attributes = FALSE
  ))
}

if (isTRUE(skip_reference)) {
  reference_raw <- data.frame()
  reference_rows <- data.frame(
    theta_row = atlas_rows$theta_row,
    theta_label = atlas_rows$theta_label,
    reference_log_marginal = NA_real_,
    reference_rep_sd = NA_real_,
    reference_rep_range = NA_real_,
    reference_mean_logZ_se = NA_real_,
    reference_n = 0L,
    check.names = FALSE
  )
} else {
  loaded_reference_cache <- FALSE
  if (isTRUE(reuse_reference_cache) &&
      !isTRUE(refresh_reference_cache) &&
      file.exists(reference_cache_file)) {
    reference_cache <- readRDS(reference_cache_file)
    if (!reference_cache_compatible(reference_cache)) {
      stop("Nested-SMC reference cache is incompatible with this run. ",
           "Use --refresh_reference_cache=true or choose another --reference_cache_file.")
    }
    reference_raw <- reference_cache$reference_raw
    reference_rows <- reference_cache$reference_rows
    loaded_reference_cache <- TRUE
    cat(sprintf("Loaded nested-SMC reference cache: %s\n", reference_cache_file))
  }

  if (!isTRUE(loaded_reference_cache)) {
    reference_jobs <- expand.grid(
      theta_row = seq_len(nrow(audit_design$theta)),
      rep_id = seq_len(ref_reps),
      KEEP.OUT.ATTRS = FALSE
    )
    reference_parts <- if (cores <= 1L || nrow(reference_jobs) <= 1L) {
      lapply(seq_len(nrow(reference_jobs)), function(k) {
        run_reference_one(reference_jobs$theta_row[k], reference_jobs$rep_id[k])
      })
    } else {
      parallel::mclapply(seq_len(nrow(reference_jobs)), function(k) {
        run_reference_one(reference_jobs$theta_row[k], reference_jobs$rep_id[k])
      }, mc.cores = min(cores, nrow(reference_jobs)))
    }
    reference_raw <- do.call(rbind, reference_parts)
    reference_rows <- do.call(rbind, lapply(split(reference_raw, reference_raw$theta_row), function(df) {
      values <- df$reference_log_marginal[is.finite(df$reference_log_marginal)]
      data.frame(
        theta_row = df$theta_row[1L],
        theta_label = df$theta_label[1L],
        reference_log_marginal = if (length(values)) stats::median(values) else NA_real_,
        reference_rep_sd = if (length(values) > 1L) stats::sd(values) else NA_real_,
        reference_rep_range = if (length(values) > 1L) diff(range(values)) else NA_real_,
        reference_mean_logZ_se = mean(df$reference_logZ_se, na.rm = TRUE),
        reference_n = length(values),
        check.names = FALSE
      )
    }))
    saveRDS(
      list(
        reference_raw = reference_raw,
        reference_rows = reference_rows,
        audit_theta = audit_design$theta,
        audit_metadata = audit_design$metadata,
        population_hyper_names = population_model$hyper_names,
        settings = reference_cache_settings
      ),
      reference_cache_file
    )
    cat(sprintf("Saved nested-SMC reference cache: %s\n", reference_cache_file))
  }
}

rows <- merge(atlas_rows, reference_rows, by = c("theta_row", "theta_label"), all.x = TRUE, sort = FALSE)
rows$atlas_error <- rows$atlas_log_marginal - rows$reference_log_marginal
rows$atlas_centered_error <- rows$atlas_error - mean(rows$atlas_error[is.finite(rows$atlas_error)], na.rm = TRUE)
rows$abs_atlas_centered_error <- abs(rows$atlas_centered_error)
rows$reference_noise_flag <- is.finite(rows$reference_rep_sd) &
  rows$reference_rep_sd > pmax(0.5, 2 * rows$reference_mean_logZ_se)
rows$atlas_certified <- rows$atlas_status %in% c("certified", "ok")
rows$atlas_reference_z <- rows$atlas_error / pmax(rows$reference_rep_sd, rows$reference_mean_logZ_se, 1e-12)
rows$surface_reference_error <- rows$surface_log_marginal - rows$reference_log_marginal
rows$particle_reference_error <- rows$particle_mis_log_marginal - rows$reference_log_marginal
rows$diagnostic_class <- ifelse(
  !rows$atlas_certified,
  paste0("uncertified:", rows$atlas_reason),
  ifelse(is.finite(rows$quadratic_particle_gap) & rows$quadratic_particle_gap > max_quadratic_particle_gap,
         "surface_particle_disagreement",
         rows$atlas_reason)
)
rows$failure_mode <- ifelse(
  !rows$atlas_certified,
  paste("uncertified", rows$atlas_reason, sep = ": "),
  ifelse(rows$diagnostic_class == "surface_particle_disagreement",
         "particle-MIS certified; derivative diagnostic disagrees",
         ifelse(rows$diagnostic_class == "active_exact_anchor",
                "exact atlas anchor",
                "particle-MIS estimate certified"))
)

graph_diag <- atlas$normalizer_solution$diagnostics %||% list()
summary <- data.frame(
  label = label,
  local_pos = local_pos,
  local = local_id,
  n_obs = nrow(data_i),
  alpha_dim = alpha_dim,
  theta_audit_points = nrow(audit_design$theta),
  atlas_design_points = nrow(atlas_design$theta),
  atlas_total_charts = length(atlas$charts),
  atlas_active_charts = sum(vapply(atlas$charts, function(chart) identical(chart$status, "active"), logical(1))),
  atlas_active_edges = sum(vapply(atlas$edges, function(edge) identical(edge$status, "active"), logical(1))),
  atlas_particles = atlas_particles,
  ref_particles = ref_particles,
  ref_reps = ref_reps,
  atlas_rmse_vs_reference = finite_rmse(rows$atlas_error),
  atlas_centered_rmse_vs_reference = finite_rmse(rows$atlas_centered_error),
  atlas_mae_vs_reference = finite_mae(rows$atlas_error),
  atlas_max_abs_centered_error = finite_max_abs(rows$atlas_centered_error),
  atlas_centered_error_range = finite_range_width(rows$atlas_centered_error),
  certified_fraction = mean(rows$atlas_certified, na.rm = TRUE),
  finite_reference_fraction = mean(is.finite(rows$reference_log_marginal)),
  median_reference_rep_sd = stats::median(rows$reference_rep_sd, na.rm = TRUE),
  median_reference_mean_logZ_se = stats::median(rows$reference_mean_logZ_se, na.rm = TRUE),
  median_particle_mis_ess_frac = stats::median(rows$particle_mis_ess_frac, na.rm = TRUE),
  median_particle_mis_psis_k = stats::median(rows$particle_mis_psis_k, na.rm = TRUE),
  median_surface_se = stats::median(rows$surface_se, na.rm = TRUE),
  max_graph_standardized_edge_residual = as.numeric(graph_diag$max_abs_standardized_edge_residual %||% NA_real_),
  max_graph_standardized_direct_residual = as.numeric(graph_diag$max_abs_standardized_direct_residual %||% NA_real_),
  runtime_minutes = as.numeric(difftime(Sys.time(), start_time, units = "mins")),
  check.names = FALSE
)

focus_for_plot <- intersect(c("mu_sv", "log_sigma2_sv", "mu_v_LogFreq", "log_sigma2_v_LogFreq"), names(rows))
if (!length(focus_for_plot)) {
  focus_for_plot <- intersect(population_model$hyper_names, names(rows))
  focus_for_plot <- focus_for_plot[seq_len(min(length(focus_for_plot), 2L))]
}

plot_or_empty <- function(x, y, xlab, ylab, main, ...) {
  ok <- is.finite(x) & is.finite(y)
  if (!any(ok)) {
    graphics::plot.new()
    graphics::title(main)
    return(invisible(FALSE))
  }
  graphics::plot(x[ok], y[ok], xlab = xlab, ylab = ylab, main = main, ...)
  invisible(TRUE)
}

point_col <- function(class) {
  palette <- c(
    active_exact_anchor = "black",
    particle_mis_active_chart_coverage = "darkorange3",
    surface_particle_disagreement = "purple4",
    no_active_chart_coverage = "firebrick3",
    surface_particle_gap = "purple4"
  )
  out <- unname(palette[class])
  missing <- is.na(out)
  out[missing] <- ifelse(grepl("^uncertified", class[missing]), "firebrick3", "grey35")
  out
}

point_pch <- function(certified) ifelse(certified, 19, 1)

class_label <- function(class) {
  labels <- c(
    active_exact_anchor = "exact atlas anchor",
    particle_mis_active_chart_coverage = "particle-MIS certified",
    surface_particle_disagreement = "particle-MIS certified; surface diagnostic disagrees",
    "uncertified:no_active_chart_coverage" = "uncertified: no chart coverage",
    "uncertified:surface_particle_gap" = "uncertified: surface/particle gap",
    no_active_chart_coverage = "no chart coverage",
    surface_particle_gap = "surface/particle gap"
  )
  out <- unname(labels[class])
  out[is.na(out)] <- class[is.na(out)]
  out
}

draw_class_legend <- function(classes, where = "topright", cex = 0.75) {
  classes <- unique(as.character(classes))
  graphics::legend(
    where,
    legend = class_label(classes),
    col = point_col(classes),
    pch = ifelse(grepl("^uncertified", classes), 1, 19),
    pt.cex = 1.1,
    cex = cex,
    bty = "n"
  )
}

label_points <- function(x, y, labels, cex = 0.75) {
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)) {
    graphics::text(x[ok], y[ok], labels = labels[ok], pos = 3, cex = cex, xpd = NA)
  }
}

active_chart_theta <- function(atlas, population_model) {
  active <- atlas$charts[vapply(atlas$charts, function(chart) identical(chart$status, "active"), logical(1))]
  if (!length(active)) {
    out <- matrix(numeric(), nrow = 0L, ncol = population_model$hyper_dim)
    colnames(out) <- population_model$hyper_names
    return(data.frame(chart_id = character(), out, check.names = FALSE))
  }
  theta <- do.call(rbind, lapply(active, function(chart) chart$theta_anchor))
  data.frame(
    chart_id = vapply(active, `[[`, character(1), "chart_id"),
    theta,
    check.names = FALSE
  )
}

grDevices::png(plot_file, width = 1800, height = 1400)
graphics::layout(matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, byrow = TRUE))
graphics::par(mar = c(4, 4, 3, 1))
ok <- is.finite(rows$reference_log_marginal) & is.finite(rows$atlas_log_marginal)
if (any(ok)) {
  lim <- range(c(rows$reference_log_marginal[ok], rows$atlas_log_marginal[ok]))
  graphics::plot(
    rows$reference_log_marginal[ok],
    rows$atlas_log_marginal[ok],
    xlab = "fresh SMC reference log m",
    ylab = "atlas log m",
    main = "Atlas vs fresh anchored SMC",
    pch = 19,
    col = ifelse(rows$atlas_certified[ok], "steelblue3", "firebrick3"),
    xlim = lim,
    ylim = lim
  )
  graphics::abline(0, 1, col = "grey35", lwd = 2)
} else {
  graphics::plot.new()
  graphics::title("Atlas vs fresh anchored SMC")
}
if (plot_or_empty(
  rows$theta_row,
  rows$atlas_centered_error,
  xlab = "audit theta row",
  ylab = "centered atlas error",
  main = "Held-out local evidence error",
  type = "b",
  pch = 19
)) {
  graphics::abline(h = 0, col = "grey35")
}
if (length(focus_for_plot) >= 1L) {
  if (plot_or_empty(
    rows[[focus_for_plot[1L]]],
    rows$atlas_centered_error,
    xlab = focus_for_plot[1L],
    ylab = "centered atlas error",
    main = "Error by focus coordinate",
    pch = 19
  )) {
    graphics::abline(h = 0, col = "grey35")
  }
} else {
  graphics::plot.new()
}
if (length(focus_for_plot) >= 2L) {
  if (plot_or_empty(
    rows[[focus_for_plot[2L]]],
    rows$atlas_centered_error,
    xlab = focus_for_plot[2L],
    ylab = "centered atlas error",
    main = "Error by variance coordinate",
    pch = 19
  )) {
    graphics::abline(h = 0, col = "grey35")
  }
} else {
  graphics::plot.new()
}
plot_or_empty(
  rows$particle_mis_ess_frac,
  rows$abs_atlas_centered_error,
  xlab = "particle-MIS ESS fraction",
  ylab = "|centered error|",
  main = "Does ESS diagnose error?",
  pch = 19
)
plot_or_empty(
  rows$reference_rep_sd,
  rows$abs_atlas_centered_error,
  xlab = "fresh SMC replicate SD",
  ylab = "|centered error|",
  main = "Reference noise vs atlas error",
  pch = 19
)
grDevices::dev.off()

plot_focus <- intersect(
  c("mu_sv", "log_sigma2_sv", "mu_v_LogFreq", "log_sigma2_v_LogFreq", "log_sigma2_v"),
  names(rows)
)
if (length(plot_focus) < 2L) {
  plot_focus <- intersect(population_model$hyper_names, names(rows))
}
plot_focus <- plot_focus[seq_len(min(length(plot_focus), 5L))]
chart_theta <- active_chart_theta(atlas, population_model)
cols <- point_col(rows$diagnostic_class)
pchs <- point_pch(rows$atlas_certified)

grDevices::png(miss_plot_file, width = 2600, height = 2600)
graphics::layout(matrix(seq_len(16L), nrow = 4L, byrow = TRUE))
graphics::par(mar = c(4, 4, 3, 1))
ok <- is.finite(rows$reference_log_marginal) & is.finite(rows$atlas_log_marginal)
if (any(ok)) {
  lim <- range(c(rows$reference_log_marginal[ok], rows$atlas_log_marginal[ok]))
  graphics::plot(
    rows$reference_log_marginal[ok],
    rows$atlas_log_marginal[ok],
    xlab = "nested SMC reference median",
    ylab = "atlas estimate",
    main = "Atlas vs strong nested SMC",
    pch = pchs[ok],
    col = cols[ok],
    xlim = lim,
    ylim = lim
  )
  graphics::abline(0, 1, col = "grey35", lwd = 2)
  if (any(is.finite(rows$reference_rep_sd[ok]))) {
    graphics::segments(
      rows$reference_log_marginal[ok] - rows$reference_rep_sd[ok],
      rows$atlas_log_marginal[ok],
      rows$reference_log_marginal[ok] + rows$reference_rep_sd[ok],
      rows$atlas_log_marginal[ok],
      col = grDevices::adjustcolor(cols[ok], alpha.f = 0.45)
    )
  }
  label_points(rows$reference_log_marginal, rows$atlas_log_marginal, rows$theta_label)
} else {
  graphics::plot.new()
  graphics::title("Atlas vs strong nested SMC")
}

plot_or_empty(
  rows$theta_row,
  rows$atlas_centered_error,
  xlab = "audit theta row",
  ylab = "centered atlas - nested SMC",
  main = "Signed error by held-out theta",
  type = "h",
  lwd = 3,
  col = cols
)
graphics::points(rows$theta_row, rows$atlas_centered_error, pch = pchs, col = cols)
graphics::abline(h = 0, col = "grey35")
label_points(rows$theta_row, rows$atlas_centered_error, rows$theta_label)

if (nrow(reference_raw)) {
  raw <- merge(reference_raw, rows[, c("theta_row", "theta_label", "atlas_log_marginal", "atlas_certified")],
               by = c("theta_row", "theta_label"), all.x = TRUE, sort = FALSE)
  y_lim <- range(c(raw$reference_log_marginal, rows$atlas_log_marginal), na.rm = TRUE)
  graphics::plot(
    raw$theta_row,
    raw$reference_log_marginal,
    xlab = "audit theta row",
    ylab = "log m",
    main = "Nested SMC replicate cloud vs atlas",
    pch = 16,
    col = grDevices::adjustcolor("grey30", alpha.f = 0.45),
    ylim = y_lim
  )
  graphics::points(rows$theta_row, rows$reference_log_marginal, pch = 18, cex = 1.5, col = "black")
  graphics::points(rows$theta_row, rows$atlas_log_marginal, pch = pchs, cex = 1.2, col = cols)
  graphics::segments(
    rows$theta_row,
    rows$reference_log_marginal - rows$reference_rep_sd,
    rows$theta_row,
    rows$reference_log_marginal + rows$reference_rep_sd,
    col = "grey45"
  )
  graphics::abline(v = rows$theta_row, col = grDevices::adjustcolor("grey70", alpha.f = 0.35))
} else {
  graphics::plot.new()
  graphics::title("Nested SMC replicate cloud vs atlas")
}

if (all(c("mu_sv", "log_sigma2_sv") %in% names(rows))) {
  graphics::plot(
    chart_theta$mu_sv,
    chart_theta$log_sigma2_sv,
    xlab = "mu_sv",
    ylab = "log_sigma2_sv",
    main = "Active charts and held-out audit points",
    pch = 4,
    col = "grey45"
  )
  graphics::points(
    rows$mu_sv,
    rows$log_sigma2_sv,
    pch = pchs,
    cex = 1.4,
    col = cols
  )
  label_points(rows$mu_sv, rows$log_sigma2_sv, rows$theta_label)
} else {
  graphics::plot.new()
  graphics::title("Active charts and held-out audit points")
}

for (coord in plot_focus[seq_len(min(length(plot_focus), 4L))]) {
  plot_or_empty(
    rows[[coord]],
    rows$atlas_centered_error,
    xlab = coord,
    ylab = "centered atlas - nested SMC",
    main = paste("Signed error by", coord),
    pch = pchs,
    col = cols
  )
  graphics::abline(h = 0, col = "grey35")
  label_points(rows[[coord]], rows$atlas_centered_error, rows$theta_label)
}
if (length(plot_focus) < 4L) {
  for (unused in seq_len(4L - length(plot_focus))) {
    graphics::plot.new()
  }
}

diagnostic_panels <- list(
  min_covering_distance = rows$min_covering_distance,
  particle_mis_ess_frac = rows$particle_mis_ess_frac,
  particle_mis_psis_k = rows$particle_mis_psis_k,
  quadratic_particle_gap = rows$quadratic_particle_gap,
  surface_se = rows$surface_se
)
for (nm in names(diagnostic_panels)) {
  plot_or_empty(
    diagnostic_panels[[nm]],
    rows$abs_atlas_centered_error,
    xlab = nm,
    ylab = "|centered error|",
    main = paste("Error vs", nm),
    pch = pchs,
    col = cols
  )
  label_points(diagnostic_panels[[nm]], rows$abs_atlas_centered_error, rows$theta_label)
}
grDevices::dev.off()

grDevices::png(reference_plot_file, width = 1800, height = 1200)
graphics::par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
if (nrow(reference_raw)) {
  raw <- merge(reference_raw, rows[, c("theta_row", "theta_label", "atlas_log_marginal", "reference_log_marginal",
                                       "reference_rep_sd", "atlas_centered_error", "atlas_status", "atlas_reason")],
               by = c("theta_row", "theta_label"), all.x = TRUE, sort = FALSE)
  y_lim <- range(c(raw$reference_log_marginal, rows$atlas_log_marginal), na.rm = TRUE)
  graphics::plot(
    raw$theta_row,
    raw$reference_log_marginal,
    xlab = "audit theta row",
    ylab = "log m",
    main = "Fresh nested-SMC replicates",
    pch = 16,
    col = grDevices::adjustcolor("grey25", alpha.f = 0.45),
    ylim = y_lim
  )
  graphics::points(rows$theta_row, rows$reference_log_marginal, pch = 18, cex = 1.6, col = "black")
  graphics::points(rows$theta_row, rows$atlas_log_marginal, pch = pchs, cex = 1.4, col = cols)
  graphics::segments(
    rows$theta_row,
    rows$reference_log_marginal - rows$reference_rep_sd,
    rows$theta_row,
    rows$reference_log_marginal + rows$reference_rep_sd,
    col = "grey35",
    lwd = 2
  )
  graphics::legend(
    "topleft",
    legend = c("nested replicate", "nested median", "atlas"),
    pch = c(16, 18, 19),
    col = c("grey35", "black", "steelblue3"),
    bty = "n"
  )
  plot_or_empty(
    rows$reference_rep_sd,
    rows$abs_atlas_centered_error,
    xlab = "nested SMC replicate SD",
    ylab = "|centered atlas error|",
    main = "Atlas error relative to reference noise",
    pch = pchs,
    col = cols
  )
  graphics::abline(0, 1, col = "grey50", lty = 2)
  plot_or_empty(
    rows$reference_mean_logZ_se,
    rows$reference_rep_sd,
    xlab = "mean reported nested SMC logZ SE",
    ylab = "empirical nested SMC replicate SD",
    main = "Reference MCSE calibration",
    pch = 19,
    col = "grey30"
  )
  graphics::abline(0, 1, col = "grey50", lty = 2)
  plot_or_empty(
    rows$theta_row,
    rows$atlas_reference_z,
    xlab = "audit theta row",
    ylab = "atlas error / nested replicate SD",
    main = "Error in nested-SMC SD units",
    type = "h",
    lwd = 3,
    col = cols
  )
  graphics::points(rows$theta_row, rows$atlas_reference_z, pch = pchs, col = cols)
  graphics::abline(h = c(-2, 0, 2), col = c("grey60", "grey35", "grey60"), lty = c(2, 1, 2))
} else {
  for (unused in 1:4) {
    graphics::plot.new()
  }
}
grDevices::dev.off()

clip_with_flag <- function(x, limit) {
  list(
    value = pmax(pmin(x, limit), -limit),
    clipped = is.finite(x) & abs(x) > limit
  )
}

grDevices::png(failure_report_file, width = 2600, height = 1800)
graphics::layout(matrix(seq_len(6L), nrow = 2L, byrow = TRUE))
graphics::par(mar = c(4.5, 4.5, 3.5, 1.2), oma = c(0, 0, 2, 0))

if (nrow(reference_raw)) {
  raw <- merge(
    reference_raw,
    rows[, c("theta_row", "theta_label", "atlas_log_marginal", "atlas_certified", "diagnostic_class")],
    by = c("theta_row", "theta_label"),
    all.x = TRUE,
    sort = FALSE
  )
  rep_center <- mean(seq_len(max(raw$replicate, na.rm = TRUE)))
  raw$x_plot <- raw$theta_row + 0.055 * (raw$replicate - rep_center)
  y_lim <- range(c(raw$reference_log_marginal, rows$reference_log_marginal,
                   rows$atlas_log_marginal), na.rm = TRUE)
  graphics::plot(
    raw$x_plot,
    raw$reference_log_marginal,
    xlab = "held-out theta point",
    ylab = "log m_i(theta)",
    main = "1. Strong nested SMC reference vs atlas",
    pch = 16,
    cex = 0.75,
    col = grDevices::adjustcolor("grey35", alpha.f = 0.45),
    xaxt = "n",
    ylim = y_lim
  )
  graphics::axis(1, at = rows$theta_row, labels = rows$theta_label, las = 2, cex.axis = 0.75)
  graphics::segments(
    rows$theta_row,
    rows$reference_log_marginal - rows$reference_rep_sd,
    rows$theta_row,
    rows$reference_log_marginal + rows$reference_rep_sd,
    col = "grey25",
    lwd = 2
  )
  graphics::points(rows$theta_row, rows$reference_log_marginal, pch = 18, cex = 1.5, col = "black")
  graphics::points(rows$theta_row, rows$atlas_log_marginal, pch = pchs, cex = 1.35, col = cols)
  graphics::legend(
    "bottomleft",
    legend = c("nested SMC replicate", "nested SMC median +/- 1 SD", "atlas estimate"),
    pch = c(16, 18, 19),
    lty = c(NA, 1, NA),
    col = c("grey35", "black", "steelblue3"),
    cex = 0.75,
    bty = "n"
  )
  draw_class_legend(rows$diagnostic_class, where = "topright", cex = 0.68)
} else {
  graphics::plot.new()
  graphics::title("1. Strong nested SMC reference vs atlas")
}

err_values <- c(rows$atlas_centered_error, -rows$reference_rep_sd, rows$reference_rep_sd)
err_values <- err_values[is.finite(err_values)]
err_ylim <- if (length(err_values)) range(err_values) else c(-1, 1)
graphics::plot(
  rows$theta_row,
  rows$atlas_centered_error,
  type = "n",
  xaxt = "n",
  ylim = err_ylim,
  xlab = "held-out theta point",
  ylab = "centered atlas error",
  main = "2. Shape error after removing common offset"
)
graphics::axis(1, at = rows$theta_row, labels = rows$theta_label, las = 2, cex.axis = 0.75)
graphics::abline(h = 0, col = "grey30")
ref_sd_median <- stats::median(rows$reference_rep_sd, na.rm = TRUE)
if (is.finite(ref_sd_median)) {
  graphics::abline(h = c(-ref_sd_median, ref_sd_median), col = "grey55", lty = 2)
}
ok_err <- is.finite(rows$atlas_centered_error)
if (any(ok_err)) {
  graphics::segments(rows$theta_row[ok_err], 0, rows$theta_row[ok_err], rows$atlas_centered_error[ok_err],
                     col = cols[ok_err], lwd = 4)
}
graphics::points(rows$theta_row, rows$atlas_centered_error, pch = pchs, col = cols, cex = 1.3)
uncert <- !rows$atlas_certified
if (any(uncert)) {
  graphics::points(rows$theta_row[uncert], rep(0, sum(uncert)), pch = 4, col = "firebrick3", cex = 1.6, lwd = 2)
}
label_points(rows$theta_row, rows$atlas_centered_error, rows$theta_label, cex = 0.7)
graphics::legend(
  "bottomleft",
  legend = c("signed centered error", "dashed: median nested-SMC SD", "red x: uncertified"),
  col = c("grey25", "grey55", "firebrick3"),
  lty = c(1, 2, NA),
  pch = c(19, NA, 4),
  cex = 0.72,
  bty = "n"
)

if (all(c("mu_sv", "log_sigma2_sv") %in% names(rows)) &&
    all(c("mu_sv", "log_sigma2_sv") %in% names(chart_theta))) {
  graphics::plot(
    chart_theta$mu_sv,
    chart_theta$log_sigma2_sv,
    xlab = "mu_sv",
    ylab = "log_sigma2_sv",
    main = "3. Coverage in the hard EMC variance plane",
    pch = 4,
    col = "grey55"
  )
  point_cex <- rep(1.25, nrow(rows))
  finite_abs <- rows$abs_atlas_centered_error[is.finite(rows$abs_atlas_centered_error)]
  if (length(finite_abs) && max(finite_abs) > 0) {
    point_cex[is.finite(rows$abs_atlas_centered_error)] <-
      1.1 + 1.5 * rows$abs_atlas_centered_error[is.finite(rows$abs_atlas_centered_error)] / max(finite_abs)
  }
  point_cex[uncert] <- 1.8
  graphics::points(rows$mu_sv, rows$log_sigma2_sv, pch = pchs, cex = point_cex, col = cols, lwd = 2)
  label_points(rows$mu_sv, rows$log_sigma2_sv, rows$theta_label, cex = 0.7)
  graphics::legend(
    "topright",
    legend = c("active atlas chart", "audit point; size = |shape error|", "open red = uncertified"),
    pch = c(4, 19, 1),
    col = c("grey55", "grey25", "firebrick3"),
    cex = 0.72,
    bty = "n"
  )
} else {
  graphics::plot.new()
  graphics::title("3. Coverage in the hard EMC variance plane")
}

clip_limit <- 2
surface_clip <- clip_with_flag(rows$surface_reference_error, clip_limit)
particle_clip <- clip_with_flag(rows$particle_reference_error, clip_limit)
final_clip <- clip_with_flag(rows$atlas_error, clip_limit)
graphics::plot(
  rows$theta_row,
  final_clip$value,
  type = "n",
  xaxt = "n",
  ylim = c(-clip_limit, clip_limit),
  xlab = "held-out theta point",
  ylab = "estimate - nested SMC median",
  main = "4. Which estimator failed? errors clipped at +/-2 log units"
)
graphics::axis(1, at = rows$theta_row, labels = rows$theta_label, las = 2, cex.axis = 0.75)
graphics::abline(h = 0, col = "grey30")
graphics::points(rows$theta_row - 0.18, surface_clip$value, pch = ifelse(surface_clip$clipped, 24, 17),
                 bg = "royalblue3", col = "royalblue3", cex = 1.25)
graphics::points(rows$theta_row, particle_clip$value, pch = ifelse(particle_clip$clipped, 22, 15),
                 bg = "darkorange3", col = "darkorange3", cex = 1.15)
graphics::points(rows$theta_row + 0.18, final_clip$value, pch = ifelse(final_clip$clipped, 21, 19),
                 bg = "black", col = "black", cex = 1.05)
graphics::legend(
  "bottomleft",
  legend = c("derivative surface", "particle-MIS", "final atlas", "triangle/square with fill: clipped"),
  pch = c(17, 15, 19, 24),
  col = c("royalblue3", "darkorange3", "black", "grey35"),
  pt.bg = c("royalblue3", "darkorange3", "black", "grey80"),
  cex = 0.72,
  bty = "n"
)

gap_x <- log10(pmax(rows$quadratic_particle_gap, 1e-6))
if (plot_or_empty(
  gap_x,
  rows$abs_atlas_centered_error,
  xlab = "log10(surface vs particle-MIS gap)",
  ylab = "|centered atlas error|",
  main = "5. Certification gate: surface/particle disagreement",
  pch = pchs,
  col = cols,
  cex = 1.25
)) {
  if (is.finite(max_quadratic_particle_gap) && max_quadratic_particle_gap > 0) {
    graphics::abline(v = log10(max_quadratic_particle_gap), col = "purple4", lty = 2, lwd = 2)
  }
  label_points(gap_x, rows$abs_atlas_centered_error, rows$theta_label, cex = 0.7)
  graphics::legend(
    "topright",
    legend = c("dashed: certification threshold"),
    lty = 2,
    col = "purple4",
    cex = 0.72,
    bty = "n"
  )
}

noise_values <- c(rows$abs_atlas_centered_error, rows$reference_rep_sd)
noise_values <- noise_values[is.finite(noise_values)]
noise_ylim <- if (length(noise_values)) range(noise_values) else c(0, 1)
if (length(noise_values)) {
  graphics::plot(
    rows$theta_row,
    rows$abs_atlas_centered_error,
    type = "b",
    xaxt = "n",
    ylim = c(0, max(noise_ylim)),
    xlab = "held-out theta point",
    ylab = "log units",
    main = "6. Is the miss bigger than nested-SMC noise?",
    pch = 19,
    col = "firebrick3"
  )
  graphics::axis(1, at = rows$theta_row, labels = rows$theta_label, las = 2, cex.axis = 0.75)
  graphics::lines(rows$theta_row, rows$reference_rep_sd, type = "b", pch = 18, col = "grey25")
  graphics::legend(
    "topright",
    legend = c("|centered atlas error|", "nested-SMC replicate SD"),
    pch = c(19, 18),
    col = c("firebrick3", "grey25"),
    cex = 0.72,
    bty = "n"
  )
} else {
  graphics::plot.new()
  graphics::title("6. Is the miss bigger than nested-SMC noise?")
}
graphics::mtext(
  sprintf("EMC local %s: atlas %d x %d, nested reference %d theta x %d reps x %d particles",
          local_id, nrow(atlas_design$theta), atlas_particles, nrow(audit_design$theta), ref_reps, ref_particles),
  outer = TRUE,
  cex = 1.1
)
grDevices::dev.off()

theta_axis <- rows$theta_label

grDevices::png(comparison_plot_file, width = 1900, height = 1350)
graphics::par(mar = c(8, 5.5, 4.5, 2), cex.lab = 1.25, cex.axis = 1.05, cex.main = 1.25)
if (nrow(reference_raw)) {
  raw <- merge(
    reference_raw,
    rows[, c("theta_row", "theta_label", "atlas_log_marginal", "atlas_certified", "diagnostic_class")],
    by = c("theta_row", "theta_label"),
    all.x = TRUE,
    sort = FALSE
  )
  rep_center <- mean(seq_len(max(raw$replicate, na.rm = TRUE)))
  raw$x_plot <- raw$theta_row + 0.07 * (raw$replicate - rep_center)
  y_lim <- range(c(raw$reference_log_marginal, rows$reference_log_marginal,
                   rows$atlas_log_marginal), na.rm = TRUE)
  graphics::plot(
    raw$x_plot,
    raw$reference_log_marginal,
    xaxt = "n",
    xlab = "",
    ylab = "local log evidence  log m_i(theta)",
    main = "Atlas estimate vs strong nested SMC reference",
    pch = 16,
    cex = 1.15,
    col = grDevices::adjustcolor("grey35", alpha.f = 0.45),
    ylim = y_lim
  )
  graphics::axis(1, at = rows$theta_row, labels = theta_axis, las = 2)
  graphics::mtext("held-out theta point", side = 1, line = 6.5, cex = 1.1)
  graphics::segments(
    rows$theta_row,
    rows$reference_log_marginal - rows$reference_rep_sd,
    rows$theta_row,
    rows$reference_log_marginal + rows$reference_rep_sd,
    col = "black",
    lwd = 3
  )
  graphics::points(rows$theta_row, rows$reference_log_marginal, pch = 18, cex = 2.0, col = "black")
  graphics::points(rows$theta_row, rows$atlas_log_marginal, pch = pchs, cex = 2.0, col = cols, lwd = 2)
  if (any(uncert)) {
    graphics::points(rows$theta_row[uncert], rows$reference_log_marginal[uncert],
                     pch = 4, cex = 2.8, col = "firebrick3", lwd = 3)
    graphics::text(rows$theta_row[uncert], rows$reference_log_marginal[uncert],
                   labels = "no atlas estimate", pos = 3, col = "firebrick3", cex = 1.0)
  }
  graphics::legend(
    "bottomleft",
    legend = c("5 independent nested SMC runs", "nested SMC median +/- 1 empirical SD", "atlas estimate", "uncertified atlas point"),
    pch = c(16, 18, 19, 4),
    lty = c(NA, 1, NA, NA),
    col = c("grey35", "black", "darkorange3", "firebrick3"),
    pt.cex = c(1.1, 1.6, 1.8, 2.0),
    lwd = c(NA, 3, NA, NA),
    cex = 1.0,
    bty = "n"
  )
  draw_class_legend(rows$diagnostic_class, where = "topright", cex = 0.95)
  graphics::mtext(
    "If a colored atlas point sits outside the black nested-SMC bar, that theta has local evidence error. Red x means the atlas refused to evaluate.",
    side = 3,
    line = 0.2,
    cex = 0.95
  )
} else {
  graphics::plot.new()
  graphics::title("Atlas estimate vs strong nested SMC reference")
}
grDevices::dev.off()

grDevices::png(coverage_plot_file, width = 1700, height = 1350)
graphics::par(mar = c(5.2, 5.5, 4.5, 2), cex.lab = 1.3, cex.axis = 1.05, cex.main = 1.25)
if (all(c("mu_sv", "log_sigma2_sv") %in% names(rows)) &&
    all(c("mu_sv", "log_sigma2_sv") %in% names(chart_theta))) {
  x_lim <- range(c(chart_theta$mu_sv, rows$mu_sv), na.rm = TRUE)
  y_lim <- range(c(chart_theta$log_sigma2_sv, rows$log_sigma2_sv), na.rm = TRUE)
  graphics::plot(
    chart_theta$mu_sv,
    chart_theta$log_sigma2_sv,
    xlab = "mu_sv",
    ylab = "log_sigma2_sv",
    main = "Where the atlas is missing coverage",
    pch = 4,
    cex = 1.1,
    col = "grey55",
    xlim = x_lim,
    ylim = y_lim
  )
  point_cex <- rep(1.6, nrow(rows))
  finite_abs <- rows$abs_atlas_centered_error[is.finite(rows$abs_atlas_centered_error)]
  if (length(finite_abs) && max(finite_abs) > 0) {
    point_cex[is.finite(rows$abs_atlas_centered_error)] <-
      1.4 + 2.2 * rows$abs_atlas_centered_error[is.finite(rows$abs_atlas_centered_error)] / max(finite_abs)
  }
  point_cex[uncert] <- 2.4
  graphics::points(rows$mu_sv, rows$log_sigma2_sv, pch = pchs, cex = point_cex, col = cols, lwd = 3)
  label_points(rows$mu_sv, rows$log_sigma2_sv, rows$theta_label, cex = 1.0)
  graphics::legend(
    "bottomright",
    legend = c("active atlas chart", "held-out audit point", "point size = |centered error|", "open red point = no active chart coverage"),
    pch = c(4, 19, 19, 1),
    col = c("grey55", "darkorange3", "grey25", "firebrick3"),
    pt.cex = c(1.2, 1.7, 2.5, 2.2),
    cex = 1.0,
    bty = "n"
  )
  if (any(uncert)) {
    graphics::text(rows$mu_sv[uncert], rows$log_sigma2_sv[uncert],
                   labels = "coverage miss", pos = 4, col = "firebrick3", cex = 1.1)
  }
  graphics::mtext(
    "The left-tail audit point is not merely inaccurate; it has no certified atlas coverage.",
    side = 3,
    line = 0.2,
    cex = 1.0
  )
} else {
  graphics::plot.new()
  graphics::title("Where the atlas is missing coverage")
}
grDevices::dev.off()

grDevices::png(estimator_plot_file, width = 1900, height = 1350)
graphics::par(mar = c(8, 5.5, 4.5, 2), cex.lab = 1.25, cex.axis = 1.05, cex.main = 1.25)
clip_limit <- 2
surface_clip <- clip_with_flag(rows$surface_reference_error, clip_limit)
particle_clip <- clip_with_flag(rows$particle_reference_error, clip_limit)
final_clip <- clip_with_flag(rows$atlas_error, clip_limit)
graphics::plot(
  rows$theta_row,
  final_clip$value,
  type = "n",
  xaxt = "n",
  ylim = c(-clip_limit, clip_limit),
  xlab = "",
  ylab = "estimate - nested SMC median, clipped to +/-2",
  main = "How the estimator misses: derivative surface vs particle-MIS vs final atlas"
)
graphics::axis(1, at = rows$theta_row, labels = theta_axis, las = 2)
graphics::mtext("held-out theta point", side = 1, line = 6.5, cex = 1.1)
graphics::abline(h = 0, col = "grey30", lwd = 2)
graphics::abline(h = c(-0.25, 0.25), col = "grey70", lty = 2)
graphics::points(rows$theta_row - 0.20, surface_clip$value, pch = ifelse(surface_clip$clipped, 24, 17),
                 bg = "royalblue3", col = "royalblue3", cex = 1.9, lwd = 2)
graphics::points(rows$theta_row, particle_clip$value, pch = ifelse(particle_clip$clipped, 22, 15),
                 bg = "darkorange3", col = "darkorange3", cex = 1.7, lwd = 2)
graphics::points(rows$theta_row + 0.20, final_clip$value, pch = ifelse(final_clip$clipped, 21, 19),
                 bg = "black", col = "black", cex = 1.55, lwd = 2)
if (any(surface_clip$clipped, na.rm = TRUE)) {
  graphics::text(
    rows$theta_row[surface_clip$clipped] - 0.20,
    surface_clip$value[surface_clip$clipped],
    labels = "clipped",
    pos = 3,
    col = "royalblue3",
    cex = 0.95
  )
}
graphics::legend(
  "bottomleft",
  legend = c("derivative surface", "particle-MIS estimate", "final atlas output", "filled triangle/square at +/-2 = actual error is off-scale"),
  pch = c(17, 15, 19, 24),
  col = c("royalblue3", "darkorange3", "black", "grey35"),
  pt.bg = c("royalblue3", "darkorange3", "black", "grey80"),
  pt.cex = c(1.6, 1.5, 1.5, 1.5),
  cex = 1.0,
  bty = "n"
)
graphics::mtext(
  "Blue off-scale points show derivative-surface diagnostics can be catastrophically wrong; the final non-anchor estimate is particle-MIS.",
  side = 3,
  line = 0.2,
  cex = 0.95
)
grDevices::dev.off()

miss_rows <- rows[order(
  !rows$atlas_certified,
  rows$abs_atlas_centered_error,
  decreasing = TRUE,
  na.last = TRUE
), , drop = FALSE]
utils::write.csv(miss_rows, miss_csv, row.names = FALSE)

utils::write.csv(rows, rows_csv, row.names = FALSE)
utils::write.csv(summary, summary_csv, row.names = FALSE)
saveRDS(
  list(
    atlas = atlas,
    atlas_design = atlas_design,
    audit_design = audit_design,
    rows = rows,
    reference_raw = reference_raw,
    miss_rows = miss_rows,
    summary = summary,
    settings = list(
      label = label,
      seed = seed,
      cores = cores,
      data_file = data_file,
      local_pos = local_pos,
      local = local_id,
      atlas_particles = atlas_particles,
      atlas_max_anchors = atlas_max_anchors,
      atlas_max_steps = atlas_max_steps,
      atlas_max_intermediates = atlas_max_intermediates,
      ref_particles = ref_particles,
      ref_reps = ref_reps,
      ref_max_steps = ref_max_steps,
      audit_max_points = audit_max_points,
      atlas_design_probs = atlas_design_probs,
      audit_probs = audit_probs,
      audit_inflation = audit_inflation,
      focus_hyper_names = focus_hyper_names,
      distance_metric = distance_metric,
      surface_method = surface_method,
      particle_mis_role = particle_mis_role,
      reference_cache_file = reference_cache_file
    )
  ),
  results_file
)

cat(sprintf("Saved results: %s\n", results_file))
cat(sprintf("Saved rows: %s\n", rows_csv))
cat(sprintf("Saved summary: %s\n", summary_csv))
cat(sprintf("Saved diagnostic plot: %s\n\n", plot_file))
cat(sprintf("Saved clear failure report: %s\n", failure_report_file))
cat(sprintf("Saved clear comparison plot: %s\n", comparison_plot_file))
cat(sprintf("Saved clear coverage plot: %s\n", coverage_plot_file))
cat(sprintf("Saved clear estimator plot: %s\n", estimator_plot_file))
cat(sprintf("Saved miss diagnostic plot: %s\n", miss_plot_file))
cat(sprintf("Saved reference replicate plot: %s\n", reference_plot_file))
cat(sprintf("Saved nested-SMC reference cache: %s\n", reference_cache_file))
cat(sprintf("Saved miss table: %s\n\n", miss_csv))
cat("Summary:\n")
print(summary)
cat("\nWorst centered atlas errors:\n")
worst <- rows[order(rows$abs_atlas_centered_error, decreasing = TRUE), , drop = FALSE]
print(worst[seq_len(min(8L, nrow(worst))), c(
  "theta_label", "atlas_status", "atlas_centered_error", "reference_rep_sd",
  "particle_mis_ess_frac", "particle_mis_psis_k", focus_for_plot
), drop = FALSE])
