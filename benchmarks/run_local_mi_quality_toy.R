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

align_rows <- function(parts) {
  cols <- unique(unlist(lapply(parts, names), use.names = FALSE))
  lapply(parts, function(df) {
    missing <- setdiff(cols, names(df))
    for (nm in missing) df[[nm]] <- NA
    df[, cols, drop = FALSE]
  })
}

all_rows <- do.call(rbind, align_rows(list(reference_rows, nested_rows, atlas_rows)))
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
    atlas_objects = lapply(atlas_parts, `[[`, "atlas"),
    scored = scored,
    summary = summary_rows,
    total_log_m = total_log_m,
    total_summary = total_summary,
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
cat(sprintf("Saved error diagnostics plot: %s\n", error_plot_file))
cat(sprintf("Saved theta-grid posterior diagnostic plot: %s\n", posterior_plot_file))
cat("\nSummary:\n")
print(summary_rows, row.names = FALSE)
cat("\nTotal log m summary:\n")
print(total_summary, row.names = FALSE)
cat("\nTheta-grid diagnostic posterior summary:\n")
print(posterior_summary, row.names = FALSE)
