#!/usr/bin/env Rscript

rm(list = ls())

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[1L]))
} else {
  normalizePath("benchmarks/run_local_mi_quality_toy.R")
}
repo_dir <- dirname(dirname(script_path))
setwd(repo_dir)

suppressPackageStartupMessages({
  library(parallel)
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

weighted_center <- function(x, w) {
  ok <- is.finite(x) & is.finite(w)
  if (!any(ok)) return(NA_real_)
  w <- pmax(w[ok], 0)
  if (sum(w) <= 0) w <- rep(1, length(w))
  sum(x[ok] * w) / sum(w)
}

source("smc_core.R")
source("population_models.R")
source("local_charts.R")
source("utilities.R")

cli_args <- parse_cli_args(commandArgs(trailingOnly = TRUE))
detected_cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
if (!is.finite(detected_cores) || detected_cores < 1L) detected_cores <- 1L

label <- arg_chr(cli_args, "label", "local_mi_quality_toy")
seed <- arg_int(cli_args, "seed", 20260524L)
cores <- arg_int(cli_args, "cores", min(4L, detected_cores))
n_locals <- arg_int(cli_args, "n_locals", 4L)
n_trials <- arg_int(cli_args, "n_trials", 8L)
tau <- arg_num(cli_args, "tau", 0.7)
true_mu <- arg_num(cli_args, "true_mu", 0)
true_sigma <- arg_num(cli_args, "true_sigma", 0.8)
theta_mu <- arg_num_vec(cli_args, "theta_mu", seq(-1.6, 1.6, length.out = 9L))
theta_sigma2 <- arg_num_vec(cli_args, "theta_sigma2", c(0.08, 0.16, 0.32, 0.64, 1.28, 2.56))
nested_particles <- arg_num_vec(cli_args, "nested_particles", c(48, 192))
nested_reps <- arg_int(cli_args, "nested_reps", 2L)
atlas_particles <- arg_int(cli_args, "atlas_particles", 512L)
atlas_max_anchors <- arg_int(cli_args, "atlas_max_anchors", min(18L, length(theta_mu) * length(theta_sigma2)))
atlas_max_steps <- arg_int(cli_args, "atlas_max_steps", 176L)
smc_max_steps <- arg_int(cli_args, "smc_max_steps", 128L)
smc_moves <- arg_int(cli_args, "smc_moves", 2L)
target_cess <- arg_num(cli_args, "target_cess", 0.9)
atlas_max_intermediates <- arg_int(cli_args, "atlas_max_intermediates", 4L)
atlas_edge_min_overlap_ess <- arg_num(cli_args, "atlas_edge_min_overlap_ess", 0.05)
atlas_edge_max_pareto_k <- arg_num(cli_args, "atlas_edge_max_pareto_k", 0.9)
atlas_edge_max_forward_reverse_gap <- arg_num(cli_args, "atlas_edge_max_forward_reverse_gap", 2)
atlas_edge_max_taylor_gap <- arg_num(cli_args, "atlas_edge_max_taylor_gap", Inf)
atlas_eval_max_chart_distance <- arg_num(cli_args, "atlas_eval_max_chart_distance", 1.0)
atlas_eval_min_covering_charts <- arg_int(cli_args, "atlas_eval_min_covering_charts", 3L)
atlas_eval_max_prediction_range <- arg_num(cli_args, "atlas_eval_max_prediction_range", Inf)
atlas_eval_distance_scale <- arg_num(cli_args, "atlas_eval_distance_scale", atlas_eval_max_chart_distance)
atlas_eval_require_particle_mis <- arg_lgl(cli_args, "atlas_eval_require_particle_mis", TRUE)
atlas_eval_min_particle_mis_ess <- arg_num(cli_args, "atlas_eval_min_particle_mis_ess", 0.05)
atlas_eval_max_particle_mis_psis_k <- arg_num(cli_args, "atlas_eval_max_particle_mis_psis_k", 0.7)
atlas_eval_max_quadratic_particle_gap <- arg_num(cli_args, "atlas_eval_max_quadratic_particle_gap", 0.02)
atlas_eval_sparse_chart_min_covering <- arg_int(cli_args, "atlas_eval_sparse_chart_min_covering", 3L)
atlas_eval_sparse_chart_max_distance <- arg_num(cli_args, "atlas_eval_sparse_chart_max_distance", 0.05)
atlas_eval_max_leave_chart_out_gap <- arg_num(cli_args, "atlas_eval_max_leave_chart_out_gap", Inf)
atlas_design_method <- arg_chr(cli_args, "atlas_design_method", "hybrid")
atlas_distance_metric <- arg_chr(cli_args, "atlas_distance_metric", "fisher")
atlas_surface_method <- arg_chr(cli_args, "atlas_surface_method", "derivative_ls")
atlas_particle_mis_role <- arg_chr(cli_args, "atlas_particle_mis_role", "estimator")
atlas_eval_min_surface_charts <- arg_int(cli_args, "atlas_eval_min_surface_charts", 2L)
atlas_eval_max_surface_se <- arg_num(cli_args, "atlas_eval_max_surface_se", Inf)
atlas_surface_value_nugget <- arg_num(cli_args, "atlas_surface_value_nugget", 0.05)
atlas_surface_gradient_weight <- arg_num(cli_args, "atlas_surface_gradient_weight", 1.0)
atlas_surface_curvature_weight <- arg_num(cli_args, "atlas_surface_curvature_weight", 0.2)
atlas_surface_ridge <- arg_num(cli_args, "atlas_surface_ridge", 1e-8)
shape_patch <- arg_lgl(cli_args, "shape_patch", FALSE)
shape_patch_particles <- arg_int(cli_args, "shape_patch_particles", atlas_particles)
shape_patch_reps <- arg_int(cli_args, "shape_patch_reps", 2L)
shape_patch_top_theta <- arg_int(cli_args, "shape_patch_top_theta", 6L)
shape_patch_diag_theta <- arg_int(cli_args, "shape_patch_diag_theta", 4L)
shape_patch_radius <- arg_int(cli_args, "shape_patch_radius", 1L)
shape_patch_max_theta <- arg_int(cli_args, "shape_patch_max_theta", 16L)
shape_patch_kernel_scale <- arg_num(cli_args, "shape_patch_kernel_scale", 0.5)
shape_patch_scale_floor <- arg_num(cli_args, "shape_patch_scale_floor", 0.25)
shape_patch_min_train <- arg_int(cli_args, "shape_patch_min_train", 4L)
shape_patch_shrink <- arg_num(cli_args, "shape_patch_shrink", 2)

results_file <- arg_chr(
  cli_args,
  "results_file",
  file.path("benchmarks", "results", paste0(label, "_results.rds"))
)
summary_csv <- arg_chr(
  cli_args,
  "summary_csv",
  file.path("benchmarks", "results", paste0(label, "_summary.csv"))
)
rows_csv <- arg_chr(
  cli_args,
  "rows_csv",
  file.path("benchmarks", "results", paste0(label, "_rows.csv"))
)
total_csv <- arg_chr(
  cli_args,
  "total_csv",
  file.path("benchmarks", "results", paste0(label, "_total_log_m.csv"))
)
posterior_csv <- arg_chr(
  cli_args,
  "posterior_csv",
  file.path("benchmarks", "results", paste0(label, "_theta_posterior_summary.csv"))
)
high_smc_csv <- arg_chr(
  cli_args,
  "high_smc_csv",
  file.path("benchmarks", "results", paste0(label, "_high_smc_comparison.csv"))
)
shape_patch_csv <- arg_chr(
  cli_args,
  "shape_patch_csv",
  file.path("benchmarks", "results", paste0(label, "_shape_patch_diagnostics.csv"))
)
error_plot_file <- arg_chr(
  cli_args,
  "error_plot_file",
  file.path("benchmarks", "results", paste0(label, "_error_diagnostics.png"))
)
posterior_plot_file <- arg_chr(
  cli_args,
  "posterior_plot_file",
  file.path("benchmarks", "results", paste0(label, "_theta_posterior_diagnostic.png"))
)

dir.create(dirname(results_file), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(summary_csv), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(rows_csv), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(total_csv), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(posterior_csv), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(high_smc_csv), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(shape_patch_csv), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(error_plot_file), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(posterior_plot_file), showWarnings = FALSE, recursive = TRUE)

set.seed(seed)
alpha_true <- stats::rnorm(n_locals, mean = true_mu, sd = true_sigma)
data_list <- lapply(alpha_true, function(alpha_i) {
  stats::rnorm(n_trials, mean = alpha_i, sd = tau)
})
names(data_list) <- sprintf("local_%02d", seq_along(data_list))

population_model <- make_population_model_diag_gaussian(
  alpha_names = "alpha",
  mean_prior_mean = 0,
  mean_prior_var = 9,
  sigma2_prior_shape = 3,
  sigma2_prior_rate = 1,
  label = "toy_normal_normal_local_mi"
)

theta_grid <- as.matrix(expand.grid(
  mu_alpha = theta_mu,
  log_sigma2_alpha = log(theta_sigma2),
  KEEP.OUT.ATTRS = FALSE
))
theta_grid <- theta_grid[, population_model$hyper_names, drop = FALSE]

root_theta <- matrix(
  c(true_mu, log(true_sigma^2)),
  nrow = 1L,
  dimnames = list(NULL, population_model$hyper_names)
)

if (identical(atlas_design_method, "axis_score")) {
  theta_distance <- rowSums(sweep(theta_grid, 2L, as.numeric(root_theta[1L, ]), "-")^2)
  axis_score <- abs(theta_grid[, "mu_alpha"] - true_mu) / max(stats::sd(theta_grid[, "mu_alpha"]), 1e-8) +
    abs(theta_grid[, "log_sigma2_alpha"] - log(true_sigma^2)) /
      max(stats::sd(theta_grid[, "log_sigma2_alpha"]), 1e-8)
  anchor_idx <- unique(c(
    which.min(theta_distance),
    head(order(axis_score, decreasing = TRUE), max(0L, atlas_max_anchors - 1L))
  ))
  anchor_idx <- head(anchor_idx, atlas_max_anchors)
  theta_design <- unique(rbind(root_theta, theta_grid[anchor_idx, , drop = FALSE]))
  theta_design <- theta_design[, population_model$hyper_names, drop = FALSE]
  theta_design <- .local_atlas_unique_theta(theta_design, population_model)
} else {
  theta_design <- .local_atlas_design_from_cloud(
    theta_cloud = theta_grid,
    population_model = population_model,
    theta_root = root_theta,
    max_anchors = atlas_max_anchors,
    axis_count = length(population_model$hyper_names),
    tail_probs = c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1),
    include_axis_profiles = TRUE,
    include_cloud_profiles = TRUE,
    design_method = atlas_design_method,
    distance_metric = atlas_distance_metric
  )
}

loglik_normal_local <- function(Theta, y_i) {
  Theta <- as.matrix(Theta)
  colnames(Theta) <- "alpha"
  alpha <- Theta[, "alpha"]
  vapply(alpha, function(a) {
    sum(stats::dnorm(y_i, mean = a, sd = tau, log = TRUE))
  }, numeric(1))
}

analytic_log_m <- function(y_i, theta) {
  theta <- .as_hyper_matrix(theta, population_model$hyper_names, population_model$hyper_dim)
  n <- length(y_i)
  ybar <- mean(y_i)
  ss <- sum((y_i - ybar)^2)
  mu <- theta[, "mu_alpha"]
  sigma2 <- exp(theta[, "log_sigma2_alpha"])
  -0.5 * (n - 1) * log(2 * pi * tau^2) -
    0.5 * log(n) -
    ss / (2 * tau^2) +
    stats::dnorm(ybar, mean = mu, sd = sqrt(sigma2 + tau^2 / n), log = TRUE)
}

run_nested_one <- function(local_pos, M, rep_id) {
  rows <- lapply(seq_len(nrow(theta_grid)), function(theta_row) {
    run <- .local_chart_run_smc(
      local_id = names(data_list)[local_pos],
      theta_anchor = theta_grid[theta_row, , drop = FALSE],
      data_i = data_list[[local_pos]],
      loglik_fn = loglik_normal_local,
      population_model = population_model,
      M = as.integer(M),
      target_cess = target_cess,
      resample_threshold = 0.5,
      n_mcmc_moves = smc_moves,
      rw_scale = 0.75,
      G_mix = 8L,
      da_enable = TRUE,
      refit_every = 2L,
      max_steps = smc_max_steps,
      deterministic_resampling = FALSE,
      n_cores = 1L,
      seed = seed + 1000003L * rep_id + 1009L * local_pos + 9176L * theta_row + as.integer(M),
      verbose = FALSE,
      source = sprintf("nested_M%s", M)
    )
    data.frame(
      method = sprintf("nested_M%d", as.integer(M)),
      local = names(data_list)[local_pos],
      local_pos = local_pos,
      theta_row = theta_row,
      replicate = rep_id,
      estimate = run$logZ,
      se = run$logZ_se,
      status = "ok",
      reason = NA_character_,
      check.names = FALSE
    )
  })
  do.call(rbind, rows)
}

cat(sprintf("Local m_i(theta) quality toy benchmark: %d locals x %d theta points\n",
            n_locals, nrow(theta_grid)))
cat("Nested particle budgets:", paste(nested_particles, collapse = ", "), "\n")
cat(sprintf("Atlas: %d anchors x %d particles\n", nrow(theta_design), atlas_particles))

reference_rows <- do.call(rbind, lapply(seq_along(data_list), function(local_pos) {
  data.frame(
    method = "analytic",
    local = names(data_list)[local_pos],
    local_pos = local_pos,
    theta_row = seq_len(nrow(theta_grid)),
    replicate = NA_integer_,
    estimate = analytic_log_m(data_list[[local_pos]], theta_grid),
    se = 0,
    status = "ok",
    reason = NA_character_,
    check.names = FALSE
  )
}))

nested_jobs <- expand.grid(
  local_pos = seq_along(data_list),
  M = as.integer(nested_particles),
  rep_id = seq_len(nested_reps),
  KEEP.OUT.ATTRS = FALSE
)
nested_parts <- if (cores <= 1L || nrow(nested_jobs) <= 1L) {
  lapply(seq_len(nrow(nested_jobs)), function(k) {
    run_nested_one(nested_jobs$local_pos[k], nested_jobs$M[k], nested_jobs$rep_id[k])
  })
} else {
  parallel::mclapply(seq_len(nrow(nested_jobs)), function(k) {
    run_nested_one(nested_jobs$local_pos[k], nested_jobs$M[k], nested_jobs$rep_id[k])
  }, mc.cores = min(cores, nrow(nested_jobs)))
}
nested_rows_raw <- do.call(rbind, nested_parts)
nested_rows <- do.call(rbind, lapply(
  split(nested_rows_raw, list(nested_rows_raw$method, nested_rows_raw$local_pos, nested_rows_raw$theta_row), drop = TRUE),
  function(df) {
    data.frame(
      method = df$method[1L],
      local = df$local[1L],
      local_pos = df$local_pos[1L],
      theta_row = df$theta_row[1L],
      replicate = NA_integer_,
      estimate = stats::median(df$estimate),
      se = if (nrow(df) > 1L) stats::sd(df$estimate) else df$se[1L],
      status = if (all(df$status == "ok")) "ok" else "error",
      reason = paste(unique(df$reason[!is.na(df$reason)]), collapse = ","),
      check.names = FALSE
    )
  }
))

atlas_parts <- lapply(seq_along(data_list), function(local_pos) {
  atlas <- build_local_atlas(
    local_id = names(data_list)[local_pos],
    data_i = data_list[[local_pos]],
    loglik_fn = loglik_normal_local,
    population_model = population_model,
    theta_design = theta_design,
    theta_root = root_theta,
    local_control = list(
      root_M = atlas_particles,
      candidate_M = atlas_particles,
      bridge_particles = atlas_particles,
      target_cess = target_cess,
      n_mcmc_moves = smc_moves,
      max_steps = atlas_max_steps,
      root_confirm = "auto"
    ),
    edge_control = list(
      max_intermediates = atlas_max_intermediates,
      min_overlap_ess = atlas_edge_min_overlap_ess,
      max_pareto_k = atlas_edge_max_pareto_k,
      max_forward_reverse_gap = atlas_edge_max_forward_reverse_gap,
      max_taylor_gap = atlas_edge_max_taylor_gap,
      max_cycle_z = 4
    ),
    n_cores = 1L,
    seed = seed + 700001L + 1009L * local_pos,
    verbose = FALSE
  )
  eval_list <- lapply(seq_len(nrow(theta_grid)), function(theta_row) {
    evaluate_local_atlas(
      atlas = atlas,
      theta = theta_grid[theta_row, , drop = FALSE],
      population_model = population_model,
      max_chart_distance = atlas_eval_max_chart_distance,
      min_covering_charts = atlas_eval_min_covering_charts,
      max_prediction_range = atlas_eval_max_prediction_range,
      distance_scale = atlas_eval_distance_scale,
      se_floor = 1e-6,
      use_particle_mis = TRUE,
      require_particle_mis = atlas_eval_require_particle_mis,
      min_particle_mis_ess = atlas_eval_min_particle_mis_ess,
      max_particle_mis_psis_k = atlas_eval_max_particle_mis_psis_k,
      max_quadratic_particle_gap = atlas_eval_max_quadratic_particle_gap,
      sparse_chart_min_covering = atlas_eval_sparse_chart_min_covering,
      sparse_chart_max_distance = atlas_eval_sparse_chart_max_distance,
      max_leave_chart_out_gap = atlas_eval_max_leave_chart_out_gap,
      distance_metric = atlas_distance_metric,
      surface_method = atlas_surface_method,
      min_surface_charts = atlas_eval_min_surface_charts,
      max_surface_se = atlas_eval_max_surface_se,
      surface_value_nugget = atlas_surface_value_nugget,
      surface_gradient_weight = atlas_surface_gradient_weight,
      surface_curvature_weight = atlas_surface_curvature_weight,
      surface_ridge = atlas_surface_ridge,
      particle_mis_role = atlas_particle_mis_role
    )
  })
  eval <- data.frame(
    log_marginal = vapply(eval_list, `[[`, numeric(1), "log_marginal"),
    se = vapply(eval_list, `[[`, numeric(1), "se"),
    status = vapply(eval_list, `[[`, character(1), "status"),
    reason = vapply(eval_list, `[[`, character(1), "reason"),
    nearest_charts = vapply(eval_list, function(x) paste(x$nearest_charts, collapse = ","), character(1)),
    n_covering_charts = vapply(eval_list, function(x) {
      rows <- x$diagnostics$predictions %||% data.frame(within_distance = logical())
      if ("within_distance" %in% names(rows)) sum(rows$within_distance %in% TRUE) else nrow(rows)
    }, integer(1)),
    prediction_range = vapply(eval_list, function(x) {
      as.numeric(x$diagnostics$prediction_range %||% NA_real_)
    }, numeric(1)),
    quadratic_particle_gap = vapply(eval_list, function(x) {
      as.numeric(x$diagnostics$quadratic_particle_gap %||% NA_real_)
    }, numeric(1)),
    particle_mis_ess_frac = vapply(eval_list, function(x) {
      as.numeric(x$diagnostics$particle_mis$ess_frac %||% NA_real_)
    }, numeric(1)),
    particle_mis_psis_k = vapply(eval_list, function(x) {
      as.numeric(x$diagnostics$particle_mis$psis_k %||% NA_real_)
    }, numeric(1)),
    leave_chart_out_gap = vapply(eval_list, function(x) {
      as.numeric(x$diagnostics$leave_chart_out_gap %||% NA_real_)
    }, numeric(1)),
    surface_method = vapply(eval_list, function(x) {
      as.character(x$diagnostics$surface_method %||% NA_character_)
    }, character(1)),
    surface_se = vapply(eval_list, function(x) {
      as.numeric(x$diagnostics$surface_se %||% NA_real_)
    }, numeric(1)),
    surface_residual_sd = vapply(eval_list, function(x) {
      as.numeric(x$diagnostics$surface_residual_sd %||% NA_real_)
    }, numeric(1)),
    min_chart_distance = vapply(eval_list, function(x) {
      rows <- x$diagnostics$predictions %||% data.frame(distance = numeric())
      d <- rows$distance
      d <- d[is.finite(d)]
      if (length(d)) min(d) else NA_real_
    }, numeric(1)),
    max_covering_distance = vapply(eval_list, function(x) {
      rows <- x$diagnostics$predictions %||% data.frame(distance = numeric())
      d <- rows$distance
      d <- d[is.finite(d)]
      if (length(d)) max(d) else NA_real_
    }, numeric(1)),
    check.names = FALSE
  )
  list(
    atlas = atlas,
    rows = data.frame(
      method = "atlas",
      local = names(data_list)[local_pos],
      local_pos = local_pos,
      theta_row = seq_len(nrow(theta_grid)),
      replicate = NA_integer_,
      estimate = eval$log_marginal,
      se = eval$se,
      status = eval$status,
      reason = eval$reason,
      nearest_charts = eval$nearest_charts,
      n_covering_charts = eval$n_covering_charts,
      prediction_range = eval$prediction_range,
      quadratic_particle_gap = eval$quadratic_particle_gap,
      particle_mis_ess_frac = eval$particle_mis_ess_frac,
      particle_mis_psis_k = eval$particle_mis_psis_k,
      leave_chart_out_gap = eval$leave_chart_out_gap,
      surface_method = eval$surface_method,
      surface_se = eval$surface_se,
      surface_residual_sd = eval$surface_residual_sd,
      min_chart_distance = eval$min_chart_distance,
      max_covering_distance = eval$max_covering_distance,
      check.names = FALSE
    )
  )
})
atlas_rows <- do.call(rbind, lapply(atlas_parts, `[[`, "rows"))

theta_neighbors <- function(theta_rows, radius = 1L) {
  radius <- max(0L, as.integer(radius))
  if (!length(theta_rows)) return(integer())
  mu_vals <- sort(unique(theta_grid[, "mu_alpha"]))
  s2_vals <- sort(unique(theta_grid[, "log_sigma2_alpha"]))
  out <- integer()
  for (theta_row in unique(as.integer(theta_rows))) {
    if (!is.finite(theta_row) || theta_row < 1L || theta_row > nrow(theta_grid)) next
    ix <- match(theta_grid[theta_row, "mu_alpha"], mu_vals)
    iy <- match(theta_grid[theta_row, "log_sigma2_alpha"], s2_vals)
    ix_set <- seq.int(max(1L, ix - radius), min(length(mu_vals), ix + radius))
    iy_set <- seq.int(max(1L, iy - radius), min(length(s2_vals), iy + radius))
    keep <- which(theta_grid[, "mu_alpha"] %in% mu_vals[ix_set] &
                    theta_grid[, "log_sigma2_alpha"] %in% s2_vals[iy_set])
    out <- c(out, keep)
  }
  sort(unique(out))
}

run_shape_patch_probe <- function(local_pos, theta_rows, M, rep_id) {
  rows <- lapply(seq_along(theta_rows), function(j) {
    theta_row <- as.integer(theta_rows[j])
    run <- .local_chart_run_smc(
      local_id = names(data_list)[local_pos],
      theta_anchor = theta_grid[theta_row, , drop = FALSE],
      data_i = data_list[[local_pos]],
      loglik_fn = loglik_normal_local,
      population_model = population_model,
      M = as.integer(M),
      target_cess = target_cess,
      resample_threshold = 0.5,
      n_mcmc_moves = smc_moves,
      rw_scale = 0.75,
      G_mix = 8L,
      da_enable = TRUE,
      refit_every = 2L,
      max_steps = smc_max_steps,
      deterministic_resampling = FALSE,
      n_cores = 1L,
      seed = seed + 5000003L * rep_id + 1009L * local_pos + 9176L * theta_row + as.integer(M),
      verbose = FALSE,
      source = "shape_patch_probe"
    )
    data.frame(
      local = names(data_list)[local_pos],
      local_pos = local_pos,
      theta_row = theta_row,
      replicate = rep_id,
      estimate = run$logZ,
      se = run$logZ_se,
      check.names = FALSE
    )
  })
  do.call(rbind, rows)
}

combine_logz_replicates <- function(df) {
  split_rows <- split(df, list(df$local_pos, df$theta_row), drop = TRUE)
  do.call(rbind, lapply(split_rows, function(part) {
    z <- as.numeric(part$estimate)
    z <- z[is.finite(z)]
    se <- as.numeric(part$se)
    se <- se[is.finite(se)]
    center <- if (length(z)) logsumexp(z) - log(length(z)) else NA_real_
    data.frame(
      local = part$local[1L],
      local_pos = part$local_pos[1L],
      theta_row = part$theta_row[1L],
      probe_log_m = center,
      probe_sd = if (length(z) > 1L) stats::sd(z) else if (length(se)) se[1L] else NA_real_,
      probe_reps = length(z),
      probe_min = if (length(z)) min(z) else NA_real_,
      probe_max = if (length(z)) max(z) else NA_real_,
      check.names = FALSE
    )
  }))
}

fit_shape_patch_local <- function(train, base_rows, posterior_weight) {
  theta <- theta_grid[, population_model$hyper_names, drop = FALSE]
  base <- base_rows$estimate[match(seq_len(nrow(theta_grid)), base_rows$theta_row)]
  train <- train[is.finite(train$delta) & is.finite(train$probe_sd), , drop = FALSE]
  if (nrow(train) < as.integer(shape_patch_min_train)) {
    return(list(
      estimate = base,
      diagnostics = data.frame(
        local = if (nrow(train)) train$local[1L] else NA_character_,
        local_pos = if (nrow(train)) train$local_pos[1L] else NA_integer_,
        status = "too_few_patch_points_no_correction",
        n_train = nrow(train),
        offset = NA_real_,
        signal_sd = NA_real_,
        loo_rmse = NA_real_,
        shape_scale = 0,
        check.names = FALSE
      )
    ))
  }
  w_post <- posterior_weight
  w_post[!is.finite(w_post) | w_post < 0] <- 0
  if (sum(w_post) <= 0) w_post <- rep(1 / nrow(theta), nrow(theta)) else w_post <- w_post / sum(w_post)
  center <- colSums(theta * w_post)
  scale <- sqrt(colSums(sweep(theta, 2L, center, "-")^2 * w_post))
  grid_scale <- apply(theta, 2L, stats::sd)
  scale <- pmax(scale, as.numeric(shape_patch_scale_floor) * grid_scale, 1e-8)
  scale[!is.finite(scale) | scale < 1e-8] <- 1
  y <- as.numeric(train$delta)
  local_weight <- posterior_weight[train$theta_row]
  local_weight[!is.finite(local_weight) | local_weight < 0] <- 0
  if (sum(local_weight) <= 0) local_weight <- rep(1, nrow(train))
  noise_var <- pmax(as.numeric(train$probe_sd), 0.05)^2
  w <- local_weight / pmax(noise_var, .Machine$double.eps)
  w <- w / mean(w)

  offset <- weighted_center(y, local_weight)
  y_shape <- y - offset
  z_all <- sweep(theta, 2L, center, "-")
  z_all <- sweep(z_all, 2L, scale, "/")
  z_train <- z_all[train$theta_row, , drop = FALSE]
  bandwidth <- max(as.numeric(shape_patch_kernel_scale), 1e-8)
  smooth_at <- function(z_query, keep = seq_along(y_shape)) {
    if (!length(keep)) return(0)
    d2 <- rowSums(sweep(z_train[keep, , drop = FALSE], 2L, z_query, "-")^2)
    kw <- exp(-0.5 * d2 / bandwidth^2) * w[keep]
    if (!any(is.finite(kw)) || sum(kw) <= 0) return(0)
    sum(kw * y_shape[keep]) / sum(kw)
  }
  pred_train <- vapply(seq_len(nrow(train)), function(j) smooth_at(z_train[j, ]), numeric(1))
  residual <- y_shape - pred_train
  loo <- rep(NA_real_, nrow(train))
  if (nrow(train) > 1L) {
    for (j in seq_len(nrow(train))) {
      keep <- setdiff(seq_len(nrow(train)), j)
      loo[j] <- y_shape[j] - smooth_at(z_train[j, ], keep = keep)
    }
  }
  signal_sd <- sqrt(weighted_center(y_shape^2, local_weight))
  loo_rmse <- finite_rmse(loo)
  if (!is.finite(loo_rmse)) loo_rmse <- sqrt(mean(residual^2))
  reliability <- if (is.finite(signal_sd) && is.finite(loo_rmse)) {
    signal_sd^2 / (signal_sd^2 + loo_rmse^2 + 1e-8)
  } else {
    0
  }
  shape_scale <- max(0, min(1, as.numeric(shape_patch_shrink) * reliability))
  pred_all <- vapply(seq_len(nrow(z_all)), function(j) smooth_at(z_all[j, ]), numeric(1))
  pred_center <- weighted_center(pred_all, w_post)
  dist <- vapply(seq_len(nrow(z_all)), function(i) {
    sqrt(min(rowSums(sweep(z_train, 2L, z_all[i, ], "-")^2)))
  }, numeric(1))
  coverage <- exp(-0.5 * (dist / bandwidth)^2)
  correction <- offset + shape_scale * (pred_all - pred_center)
  list(
    estimate = base + correction,
    diagnostics = data.frame(
      local = train$local[1L],
      local_pos = train$local_pos[1L],
      status = "fitted",
      n_train = nrow(train),
      offset = offset,
      signal_sd = signal_sd,
      loo_rmse = loo_rmse,
      training_rmse = sqrt(mean(residual^2)),
      shape_scale = shape_scale,
      bandwidth = bandwidth,
      scale_floor = as.numeric(shape_patch_scale_floor),
      mean_coverage = mean(coverage),
      min_coverage = min(coverage),
      max_coverage = max(coverage),
      check.names = FALSE
    )
  )
}

shape_patch_rows <- data.frame()
shape_patch_diagnostics <- data.frame()
shape_patch_replicates <- data.frame()
if (isTRUE(shape_patch)) {
  atlas_total <- do.call(rbind, lapply(
    split(atlas_rows, atlas_rows$theta_row),
    function(df) {
      data.frame(
        theta_row = df$theta_row[1L],
        total_log_m = if (all(is.finite(df$estimate)) && nrow(df) == n_locals) sum(df$estimate) else NA_real_,
        check.names = FALSE
      )
    }
  ))
  atlas_total <- atlas_total[match(seq_len(nrow(theta_grid)), atlas_total$theta_row), , drop = FALSE]
  logprior_patch <- population_model_log_hyperprior(population_model, theta_grid)
  logw_patch <- logprior_patch + atlas_total$total_log_m
  logw_patch[!is.finite(logw_patch)] <- -Inf
  posterior_weight <- if (all(!is.finite(logw_patch))) {
    rep(1 / nrow(theta_grid), nrow(theta_grid))
  } else {
    exp(logw_patch - logsumexp(logw_patch))
  }
  top_theta <- head(order(posterior_weight, decreasing = TRUE), max(1L, as.integer(shape_patch_top_theta)))
  diag_score <- rep(0, nrow(theta_grid))
  if (nrow(atlas_rows)) {
    atlas_rows$diag_score <- (
      ifelse(atlas_rows$status != "certified", 5, 0) +
        pmax(as.numeric(atlas_rows$particle_mis_psis_k) - atlas_eval_max_particle_mis_psis_k, 0, na.rm = TRUE) +
        pmax(atlas_eval_min_particle_mis_ess - as.numeric(atlas_rows$particle_mis_ess_frac), 0, na.rm = TRUE) /
          max(atlas_eval_min_particle_mis_ess, 1e-8) +
        pmin(pmax(as.numeric(atlas_rows$se), 0), 2)
    )
    diag_theta <- aggregate(
      diag_score ~ theta_row,
      atlas_rows,
      max,
      na.rm = TRUE
    )
    diag_score[diag_theta$theta_row] <- diag_theta$diag_score
  }
  diag_theta <- head(order(diag_score * sqrt(pmax(posterior_weight, 1e-12)), decreasing = TRUE),
                     max(0L, as.integer(shape_patch_diag_theta)))
  patch_theta_rows <- theta_neighbors(unique(c(top_theta, diag_theta)), radius = shape_patch_radius)
  if (length(patch_theta_rows) > as.integer(shape_patch_max_theta)) {
    rank <- posterior_weight[patch_theta_rows] + 0.1 * diag_score[patch_theta_rows]
    patch_theta_rows <- patch_theta_rows[head(order(rank, decreasing = TRUE), as.integer(shape_patch_max_theta))]
    patch_theta_rows <- sort(unique(patch_theta_rows))
  }
  cat("Shape-patch theta rows:", paste(patch_theta_rows, collapse = ", "), "\n")
  patch_jobs <- expand.grid(
    local_pos = seq_along(data_list),
    rep_id = seq_len(max(1L, as.integer(shape_patch_reps))),
    KEEP.OUT.ATTRS = FALSE
  )
  patch_parts <- if (cores <= 1L || nrow(patch_jobs) <= 1L) {
    lapply(seq_len(nrow(patch_jobs)), function(k) {
      run_shape_patch_probe(
        local_pos = patch_jobs$local_pos[k],
        theta_rows = patch_theta_rows,
        M = shape_patch_particles,
        rep_id = patch_jobs$rep_id[k]
      )
    })
  } else {
    parallel::mclapply(seq_len(nrow(patch_jobs)), function(k) {
      run_shape_patch_probe(
        local_pos = patch_jobs$local_pos[k],
        theta_rows = patch_theta_rows,
        M = shape_patch_particles,
        rep_id = patch_jobs$rep_id[k]
      )
    }, mc.cores = min(cores, nrow(patch_jobs)))
  }
  shape_patch_replicates <- do.call(rbind, patch_parts)
  patch_summary <- combine_logz_replicates(shape_patch_replicates)
  base_patch <- atlas_rows[, c("local", "local_pos", "theta_row", "estimate", "se"), drop = FALSE]
  names(base_patch)[names(base_patch) == "estimate"] <- "atlas_estimate"
  names(base_patch)[names(base_patch) == "se"] <- "atlas_se"
  patch_summary <- merge(patch_summary, base_patch, by = c("local", "local_pos", "theta_row"), all.x = TRUE, sort = FALSE)
  patch_summary$delta <- patch_summary$probe_log_m - patch_summary$atlas_estimate
  patch_summary$theta_weight <- posterior_weight[patch_summary$theta_row]
  patch_fits <- lapply(seq_along(data_list), function(local_pos) {
    train <- patch_summary[patch_summary$local_pos == local_pos, , drop = FALSE]
    base <- atlas_rows[atlas_rows$local_pos == local_pos, , drop = FALSE]
    fit_shape_patch_local(train, base, posterior_weight)
  })
  shape_patch_rows <- do.call(rbind, lapply(seq_along(patch_fits), function(local_pos) {
    est <- patch_fits[[local_pos]]$estimate
    data.frame(
      method = "atlas_shape_patch",
      local = names(data_list)[local_pos],
      local_pos = local_pos,
      theta_row = seq_len(nrow(theta_grid)),
      replicate = NA_integer_,
      estimate = est,
      se = NA_real_,
      status = ifelse(is.finite(est), "ok", "uncertified"),
      reason = ifelse(is.finite(est), "posterior_shape_patch", "shape_patch_failed"),
      check.names = FALSE
    )
  }))
  shape_patch_diagnostics <- do.call(rbind, lapply(patch_fits, `[[`, "diagnostics"))
  design_diag <- data.frame(
    local = "ALL",
    local_pos = NA_integer_,
    status = "design",
    n_train = length(patch_theta_rows),
    offset = NA_real_,
    signal_sd = NA_real_,
    loo_rmse = NA_real_,
    training_rmse = NA_real_,
    shape_scale = NA_real_,
    mean_coverage = NA_real_,
    min_coverage = NA_real_,
    max_coverage = NA_real_,
    patch_theta_rows = paste(patch_theta_rows, collapse = ","),
    check.names = FALSE
  )
  shape_patch_diagnostics$patch_theta_rows <- NA_character_
  diag_cols <- unique(c(names(design_diag), names(shape_patch_diagnostics)))
  for (nm in setdiff(diag_cols, names(design_diag))) design_diag[[nm]] <- NA
  for (nm in setdiff(diag_cols, names(shape_patch_diagnostics))) shape_patch_diagnostics[[nm]] <- NA
  shape_patch_diagnostics <- rbind(
    design_diag[, diag_cols, drop = FALSE],
    shape_patch_diagnostics[, diag_cols, drop = FALSE]
  )
  utils::write.csv(shape_patch_diagnostics, shape_patch_csv, row.names = FALSE)
}

align_rows <- function(parts) {
  cols <- unique(unlist(lapply(parts, names), use.names = FALSE))
  lapply(parts, function(df) {
    missing <- setdiff(cols, names(df))
    for (nm in missing) df[[nm]] <- NA
    df[, cols, drop = FALSE]
  })
}

all_rows <- do.call(rbind, align_rows(list(reference_rows, nested_rows, atlas_rows, shape_patch_rows)))
theta_frame <- data.frame(theta_row = seq_len(nrow(theta_grid)), theta_grid, check.names = FALSE)
truth <- reference_rows[, c("local_pos", "theta_row", "estimate")]
names(truth)[3L] <- "truth"
scored <- merge(all_rows[all_rows$method != "analytic", ], truth, by = c("local_pos", "theta_row"), all.x = TRUE)
scored <- merge(scored, theta_frame, by = "theta_row", all.x = TRUE, sort = FALSE)
scored$error <- scored$estimate - scored$truth
scored$abs_error <- abs(scored$error)

summary_rows <- do.call(rbind, lapply(split(scored, scored$method), function(df) {
  centered_error <- unlist(lapply(split(df, df$local_pos), function(part) {
    part$error - mean(part$error[is.finite(part$error)], na.rm = TRUE)
  }))
  data.frame(
    method = df$method[1L],
    n_rows = nrow(df),
    certified_fraction = mean(df$status == "certified" | df$status == "ok", na.rm = TRUE),
    finite_fraction = mean(is.finite(df$error)),
    rmse = finite_rmse(df$error),
    centered_rmse = finite_rmse(centered_error),
    mae = finite_mae(df$error),
    max_abs_error = finite_max_abs(df$error),
    mean_reported_se = mean(df$se, na.rm = TRUE),
    median_particle_mis_ess_frac = stats::median(df$particle_mis_ess_frac, na.rm = TRUE),
    median_particle_mis_psis_k = stats::median(df$particle_mis_psis_k, na.rm = TRUE),
    median_surface_se = stats::median(df$surface_se, na.rm = TRUE),
    median_surface_residual_sd = stats::median(df$surface_residual_sd, na.rm = TRUE),
    check.names = FALSE
  )
}))
summary_rows <- summary_rows[order(summary_rows$rmse), , drop = FALSE]

estimate_rows <- all_rows[all_rows$method %in% c("analytic", unique(scored$method)), , drop = FALSE]
total_log_m <- do.call(rbind, lapply(
  split(estimate_rows, list(estimate_rows$method, estimate_rows$theta_row), drop = TRUE),
  function(df) {
    complete <- nrow(df) == n_locals && all(is.finite(df$estimate))
    data.frame(
      method = df$method[1L],
      theta_row = df$theta_row[1L],
      total_log_m = if (complete) sum(df$estimate) else NA_real_,
      n_locals = nrow(df),
      n_finite_locals = sum(is.finite(df$estimate)),
      complete = complete,
      check.names = FALSE
    )
  }
))
total_log_m <- merge(total_log_m, theta_frame, by = "theta_row", all.x = TRUE, sort = FALSE)
truth_total <- total_log_m[total_log_m$method == "analytic", c("theta_row", "total_log_m")]
names(truth_total)[2L] <- "analytic_total_log_m"
total_log_m <- merge(total_log_m, truth_total, by = "theta_row", all.x = TRUE, sort = FALSE)
total_log_m$total_error <- total_log_m$total_log_m - total_log_m$analytic_total_log_m
total_log_m <- do.call(rbind, lapply(split(total_log_m, total_log_m$method), function(df) {
  df$centered_total_error <- df$total_error - mean(df$total_error[is.finite(df$total_error)], na.rm = TRUE)
  df
}))
total_summary <- do.call(rbind, lapply(
  split(total_log_m[total_log_m$method != "analytic", ], total_log_m$method[total_log_m$method != "analytic"]),
  function(df) {
    data.frame(
      method = df$method[1L],
      n_theta = nrow(df),
      finite_theta_fraction = mean(is.finite(df$total_error)),
      total_rmse = finite_rmse(df$total_error),
      centered_total_rmse = finite_rmse(df$centered_total_error),
      total_mae = finite_mae(df$total_error),
      max_abs_total_error = finite_max_abs(df$total_error),
      centered_total_range = finite_range_width(df$centered_total_error),
      check.names = FALSE
    )
  }
))
total_summary <- total_summary[order(total_summary$centered_total_rmse), , drop = FALSE]

high_smc_method <- sprintf("nested_M%d", max(as.integer(nested_particles)))
high_ref_local <- all_rows[all_rows$method == high_smc_method, c("local_pos", "theta_row", "estimate"), drop = FALSE]
names(high_ref_local)[3L] <- "high_smc_estimate"
high_cmp_local <- merge(
  all_rows[all_rows$method != high_smc_method, , drop = FALSE],
  high_ref_local,
  by = c("local_pos", "theta_row"),
  all.x = TRUE,
  sort = FALSE
)
high_cmp_local$local_error_vs_high_smc <- high_cmp_local$estimate - high_cmp_local$high_smc_estimate
high_ref_total <- total_log_m[total_log_m$method == high_smc_method, c("theta_row", "total_log_m"), drop = FALSE]
names(high_ref_total)[2L] <- "high_smc_total_log_m"
high_cmp_total <- merge(
  total_log_m[total_log_m$method != high_smc_method, , drop = FALSE],
  high_ref_total,
  by = "theta_row",
  all.x = TRUE,
  sort = FALSE
)
high_cmp_total$total_error_vs_high_smc <- high_cmp_total$total_log_m - high_cmp_total$high_smc_total_log_m
high_cmp_total <- do.call(rbind, lapply(split(high_cmp_total, high_cmp_total$method), function(df) {
  df$centered_total_error_vs_high_smc <- df$total_error_vs_high_smc -
    mean(df$total_error_vs_high_smc[is.finite(df$total_error_vs_high_smc)], na.rm = TRUE)
  df
}))
high_smc_summary <- do.call(rbind, lapply(split(high_cmp_local, high_cmp_local$method), function(df) {
  centered_local <- unlist(lapply(split(df, df$local_pos), function(part) {
    part$local_error_vs_high_smc -
      mean(part$local_error_vs_high_smc[is.finite(part$local_error_vs_high_smc)], na.rm = TRUE)
  }))
  tdf <- high_cmp_total[high_cmp_total$method == df$method[1L], , drop = FALSE]
  data.frame(
    method = df$method[1L],
    high_smc_method = high_smc_method,
    local_rmse_vs_high_smc = finite_rmse(df$local_error_vs_high_smc),
    local_centered_rmse_vs_high_smc = finite_rmse(centered_local),
    local_mae_vs_high_smc = finite_mae(df$local_error_vs_high_smc),
    total_rmse_vs_high_smc = finite_rmse(tdf$total_error_vs_high_smc),
    centered_total_rmse_vs_high_smc = finite_rmse(tdf$centered_total_error_vs_high_smc),
    finite_theta_fraction_vs_high_smc = mean(is.finite(tdf$total_error_vs_high_smc)),
    check.names = FALSE
  )
}))
high_smc_summary <- high_smc_summary[
  order(high_smc_summary$centered_total_rmse_vs_high_smc,
        high_smc_summary$local_centered_rmse_vs_high_smc),
  ,
  drop = FALSE
]

logprior <- population_model_log_hyperprior(population_model, theta_grid)
posterior_grid <- do.call(rbind, lapply(split(total_log_m, total_log_m$method), function(df) {
  df <- df[match(seq_len(nrow(theta_grid)), df$theta_row), , drop = FALSE]
  logp <- logprior + df$total_log_m
  logp[!is.finite(logp)] <- -Inf
  if (all(!is.finite(logp))) {
    weight <- rep(NA_real_, length(logp))
  } else {
    weight <- exp(logp - logsumexp(logp))
  }
  data.frame(
    method = df$method,
    theta_row = df$theta_row,
    weight = weight,
    theta_grid,
    check.names = FALSE
  )
}))
analytic_post <- posterior_grid[posterior_grid$method == "analytic", , drop = FALSE]
analytic_mu_mean <- sum(analytic_post$weight * analytic_post$mu_alpha)
analytic_log_sigma2_mean <- sum(analytic_post$weight * analytic_post$log_sigma2_alpha)
posterior_summary <- do.call(rbind, lapply(split(posterior_grid, posterior_grid$method), function(df) {
  w <- pmax(df$weight, 0)
  if (!any(is.finite(w)) || sum(w, na.rm = TRUE) <= 0) {
    return(data.frame(
      method = df$method[1L],
      mu_mean = NA_real_,
      log_sigma2_mean = NA_real_,
      mu_sd = NA_real_,
      log_sigma2_sd = NA_real_,
      mu_mean_error_vs_analytic = NA_real_,
      log_sigma2_mean_error_vs_analytic = NA_real_,
      max_weight = NA_real_,
      ess_grid = NA_real_,
      check.names = FALSE
    ))
  }
  w[!is.finite(w)] <- 0
  w <- w / sum(w)
  mu_mean <- sum(w * df$mu_alpha)
  log_sigma2_mean <- sum(w * df$log_sigma2_alpha)
  data.frame(
    method = df$method[1L],
    mu_mean = mu_mean,
    log_sigma2_mean = log_sigma2_mean,
    mu_sd = sqrt(sum(w * (df$mu_alpha - mu_mean)^2)),
    log_sigma2_sd = sqrt(sum(w * (df$log_sigma2_alpha - log_sigma2_mean)^2)),
    mu_mean_error_vs_analytic = mu_mean - analytic_mu_mean,
    log_sigma2_mean_error_vs_analytic = log_sigma2_mean - analytic_log_sigma2_mean,
    max_weight = max(w),
    ess_grid = 1 / sum(w^2),
    check.names = FALSE
  )
}))

.grid_matrix <- function(df, value_col, theta_grid) {
  mu_vals <- sort(unique(theta_grid[, "mu_alpha"]))
  log_s2_vals <- sort(unique(theta_grid[, "log_sigma2_alpha"]))
  mat <- matrix(NA_real_, nrow = length(mu_vals), ncol = length(log_s2_vals))
  rownames(mat) <- signif(mu_vals, 4)
  colnames(mat) <- signif(log_s2_vals, 4)
  for (i in seq_len(nrow(df))) {
    ix <- match(df$mu_alpha[i], mu_vals)
    iy <- match(df$log_sigma2_alpha[i], log_s2_vals)
    mat[ix, iy] <- df[[value_col]][i]
  }
  list(x = mu_vals, y = log_s2_vals, z = mat)
}

.plot_grid_image <- function(grid, main, zlim = NULL, palette = grDevices::hcl.colors(81, "RdBu", rev = TRUE)) {
  z <- grid$z
  if (is.null(zlim)) {
    finite <- z[is.finite(z)]
    lim <- if (length(finite)) max(abs(finite)) else 1
    zlim <- c(-lim, lim)
  }
  graphics::image(
    grid$x,
    grid$y,
    z,
    col = palette,
    zlim = zlim,
    xlab = "mu_alpha",
    ylab = "log_sigma2_alpha",
    main = main
  )
  graphics::points(theta_design[, "mu_alpha"], theta_design[, "log_sigma2_alpha"], pch = 4, cex = 0.8)
  if (length(grid$x) > 1L && length(grid$y) > 1L && any(is.finite(z))) {
    graphics::contour(grid$x, grid$y, z, add = TRUE, drawlabels = FALSE, nlevels = 8, col = "grey30")
  }
}

plot_methods <- c("analytic", summary_rows$method)
plot_methods <- unique(plot_methods[plot_methods %in% unique(total_log_m$method)])
error_methods <- summary_rows$method
worst_local <- scored$local_pos[which.max(scored$abs_error)]

grDevices::png(error_plot_file, width = 1800, height = 1400)
graphics::layout(matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, byrow = TRUE))
graphics::par(mar = c(4, 4, 3, 1))
bar_mat <- t(as.matrix(summary_rows[, c("rmse", "centered_rmse", "mae"), drop = FALSE]))
colnames(bar_mat) <- summary_rows$method
graphics::barplot(
  bar_mat,
  beside = TRUE,
  las = 2,
  col = c("grey40", "steelblue3", "tan3"),
  ylab = "local log m error",
  main = "Per-local error vs analytic"
)
graphics::legend("topright", legend = rownames(bar_mat), fill = c("grey40", "steelblue3", "tan3"), bty = "n")
bar_total <- t(as.matrix(total_summary[, c("total_rmse", "centered_total_rmse", "total_mae"), drop = FALSE]))
colnames(bar_total) <- total_summary$method
graphics::barplot(
  bar_total,
  beside = TRUE,
  las = 2,
  col = c("grey40", "steelblue3", "tan3"),
  ylab = "sum_i log m_i error",
  main = "Total local evidence error"
)
graphics::legend("topright", legend = rownames(bar_total), fill = c("grey40", "steelblue3", "tan3"), bty = "n")
for (method in head(error_methods, 4L)) {
  df <- total_log_m[total_log_m$method == method, , drop = FALSE]
  grid <- .grid_matrix(df, "centered_total_error", theta_grid)
  .plot_grid_image(grid, paste(method, "centered total error"))
}
grDevices::dev.off()

grDevices::png(posterior_plot_file, width = 1800, height = 1400)
n_methods <- length(plot_methods)
n_col <- min(3L, n_methods)
n_row <- ceiling(n_methods / n_col)
graphics::par(mfrow = c(n_row, n_col), mar = c(4, 4, 3, 1))
weight_max <- max(posterior_grid$weight[posterior_grid$method %in% plot_methods], na.rm = TRUE)
for (method in plot_methods) {
  df <- posterior_grid[posterior_grid$method == method, , drop = FALSE]
  grid <- .grid_matrix(df, "weight", theta_grid)
  .plot_grid_image(
    grid,
    paste(method, "theta-grid diagnostic posterior"),
    zlim = c(0, weight_max),
    palette = grDevices::hcl.colors(81, "YlOrRd", rev = FALSE)
  )
}
grDevices::dev.off()

utils::write.csv(summary_rows, summary_csv, row.names = FALSE)
utils::write.csv(scored, rows_csv, row.names = FALSE)
utils::write.csv(total_summary, total_csv, row.names = FALSE)
utils::write.csv(posterior_summary, posterior_csv, row.names = FALSE)
utils::write.csv(high_smc_summary, high_smc_csv, row.names = FALSE)
saveRDS(
  list(
    data_list = data_list,
    alpha_true = alpha_true,
    theta_grid = theta_grid,
    theta_design = theta_design,
    reference_rows = reference_rows,
    nested_rows_raw = nested_rows_raw,
    nested_rows = nested_rows,
    atlas_rows = atlas_rows,
    shape_patch_rows = shape_patch_rows,
    shape_patch_replicates = shape_patch_replicates,
    shape_patch_diagnostics = shape_patch_diagnostics,
    atlas_objects = lapply(atlas_parts, `[[`, "atlas"),
    scored = scored,
    summary = summary_rows,
    total_log_m = total_log_m,
    total_summary = total_summary,
    high_smc_summary = high_smc_summary,
    high_smc_comparison_local = high_cmp_local,
    high_smc_comparison_total = high_cmp_total,
    posterior_grid = posterior_grid,
    posterior_summary = posterior_summary,
    settings = list(
      label = label,
      seed = seed,
      cores = cores,
      n_locals = n_locals,
      n_trials = n_trials,
      tau = tau,
      true_mu = true_mu,
      true_sigma = true_sigma,
      nested_particles = nested_particles,
      nested_reps = nested_reps,
      atlas_particles = atlas_particles,
      atlas_max_anchors = atlas_max_anchors,
      atlas_max_intermediates = atlas_max_intermediates,
      atlas_edge_min_overlap_ess = atlas_edge_min_overlap_ess,
      atlas_edge_max_pareto_k = atlas_edge_max_pareto_k,
      atlas_edge_max_forward_reverse_gap = atlas_edge_max_forward_reverse_gap,
      atlas_edge_max_taylor_gap = atlas_edge_max_taylor_gap,
      atlas_design_method = atlas_design_method,
      atlas_distance_metric = atlas_distance_metric,
      atlas_surface_method = atlas_surface_method,
      atlas_particle_mis_role = atlas_particle_mis_role,
      atlas_eval_min_surface_charts = atlas_eval_min_surface_charts,
      atlas_eval_max_surface_se = atlas_eval_max_surface_se,
      atlas_surface_value_nugget = atlas_surface_value_nugget,
      atlas_surface_gradient_weight = atlas_surface_gradient_weight,
      atlas_surface_curvature_weight = atlas_surface_curvature_weight,
      atlas_surface_ridge = atlas_surface_ridge,
      shape_patch = shape_patch,
      shape_patch_particles = shape_patch_particles,
      shape_patch_reps = shape_patch_reps,
      shape_patch_top_theta = shape_patch_top_theta,
      shape_patch_diag_theta = shape_patch_diag_theta,
      shape_patch_radius = shape_patch_radius,
      shape_patch_max_theta = shape_patch_max_theta,
      shape_patch_kernel_scale = shape_patch_kernel_scale,
      shape_patch_scale_floor = shape_patch_scale_floor,
      shape_patch_min_train = shape_patch_min_train,
      shape_patch_shrink = shape_patch_shrink,
      atlas_eval_max_chart_distance = atlas_eval_max_chart_distance,
      atlas_eval_min_covering_charts = atlas_eval_min_covering_charts,
      atlas_eval_max_prediction_range = atlas_eval_max_prediction_range,
      atlas_eval_distance_scale = atlas_eval_distance_scale,
      atlas_eval_require_particle_mis = atlas_eval_require_particle_mis,
      atlas_eval_min_particle_mis_ess = atlas_eval_min_particle_mis_ess,
      atlas_eval_max_particle_mis_psis_k = atlas_eval_max_particle_mis_psis_k,
      atlas_eval_max_quadratic_particle_gap = atlas_eval_max_quadratic_particle_gap,
      atlas_eval_sparse_chart_min_covering = atlas_eval_sparse_chart_min_covering,
      atlas_eval_sparse_chart_max_distance = atlas_eval_sparse_chart_max_distance,
      atlas_eval_max_leave_chart_out_gap = atlas_eval_max_leave_chart_out_gap,
      smc_moves = smc_moves,
      target_cess = target_cess
    ),
    summary_csv = summary_csv,
    rows_csv = rows_csv,
    total_csv = total_csv,
    posterior_csv = posterior_csv,
    high_smc_csv = high_smc_csv,
    shape_patch_csv = shape_patch_csv,
    error_plot_file = error_plot_file,
    posterior_plot_file = posterior_plot_file
  ),
  results_file
)

cat(sprintf("Saved local m_i quality results: %s\n", results_file))
cat(sprintf("Saved summary: %s\n", summary_csv))
cat(sprintf("Saved rows: %s\n", rows_csv))
cat(sprintf("Saved total log m summary: %s\n", total_csv))
cat(sprintf("Saved theta-grid posterior summary: %s\n", posterior_csv))
cat(sprintf("Saved high-SMC comparison: %s\n", high_smc_csv))
if (isTRUE(shape_patch)) {
  cat(sprintf("Saved shape-patch diagnostics: %s\n", shape_patch_csv))
}
cat(sprintf("Saved error diagnostics plot: %s\n", error_plot_file))
cat(sprintf("Saved theta-grid posterior diagnostic plot: %s\n", posterior_plot_file))
cat("\nSummary:\n")
print(summary_rows, row.names = FALSE)
cat("\nTotal log m summary:\n")
print(total_summary, row.names = FALSE)
cat("\nHigh-particle nested SMC comparison:\n")
print(high_smc_summary, row.names = FALSE)
cat("\nTheta-grid diagnostic posterior summary:\n")
print(posterior_summary, row.names = FALSE)
