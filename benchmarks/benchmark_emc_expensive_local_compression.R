#!/usr/bin/env Rscript

rm(list = ls())

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[1L]))
} else {
  normalizePath("benchmarks/benchmark_emc_expensive_local_compression.R")
}
repo_dir <- dirname(dirname(script_path))
setwd(repo_dir)

suppressPackageStartupMessages({
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

arg_num <- function(args, key, default) {
  val <- args[[key]]
  if (is.null(val) || !nzchar(val)) return(as.numeric(default))
  as.numeric(val)
}

arg_lgl <- function(args, key, default = FALSE) {
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

finite_rmse <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x)) sqrt(mean(x * x)) else NA_real_
}

center_finite <- function(x) {
  ok <- is.finite(x)
  out <- rep(NA_real_, length(x))
  if (any(ok)) out[ok] <- x[ok] - mean(x[ok])
  out
}

systematic_resample <- function(prob, K, seed = NULL) {
  prob <- pmax(as.numeric(prob), 0)
  prob <- prob / sum(prob)
  if (!is.null(seed)) set.seed(as.integer(seed))
  positions <- (stats::runif(1L) + seq.int(0L, K - 1L)) / K
  idx <- findInterval(positions, cumsum(prob)) + 1L
  idx[idx > length(prob)] <- length(prob)
  tab <- tabulate(idx, nbins = length(prob))
  selected <- which(tab > 0L)
  list(index = selected, weight = as.numeric(tab[selected]) / K)
}

top_weight_select <- function(prob, K) {
  prob <- pmax(as.numeric(prob), 0)
  ord <- order(prob, decreasing = TRUE)
  selected <- ord[seq_len(min(as.integer(K), length(ord)))]
  weight <- prob[selected]
  weight <- weight / sum(weight)
  list(index = selected, weight = weight)
}

fit_simplex_weights <- function(features, target, ridge = 1e-8) {
  features <- as.matrix(features)
  target <- as.numeric(target)
  K <- nrow(features)
  if (K <= 0L) stop("features must have at least one selected row.")
  if (K == 1L) return(1)
  X <- t(features)
  Dmat <- 2 * (crossprod(X) + diag(as.numeric(ridge), K))
  dvec <- as.numeric(2 * crossprod(X, target))
  Amat <- cbind(rep(1, K), diag(K))
  bvec <- c(1, rep(0, K))
  sol <- tryCatch(
    quadprog::solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec, meq = 1L)$solution,
    error = function(e) NULL
  )
  if (is.null(sol) || any(!is.finite(sol))) sol <- rep(1 / K, K)
  sol <- pmax(as.numeric(sol), 0)
  sw <- sum(sol)
  if (!is.finite(sw) || sw <= 0) sol <- rep(1 / K, K) else sol <- sol / sw
  sol
}

sparse_quadrature_select <- function(features, target, K, ridge = 1e-8) {
  features <- as.matrix(features)
  target <- as.numeric(target)
  K <- min(as.integer(K), nrow(features))
  if (K <= 0L) stop("K must be positive.")
  selected <- integer(0)
  active_weight <- numeric(0)
  current <- rep(0, ncol(features))
  available <- rep(TRUE, nrow(features))
  for (k in seq_len(K)) {
    residual <- target - current
    scores <- as.numeric(features %*% residual)
    scores[!available] <- -Inf
    pick <- which.max(scores)
    if (!is.finite(scores[pick])) break
    selected <- c(selected, pick)
    available[pick] <- FALSE
    active_features <- features[selected, , drop = FALSE]
    active_weight <- fit_simplex_weights(active_features, target, ridge = ridge)
    keep <- active_weight > 1e-12
    available[selected[!keep]] <- FALSE
    selected <- selected[keep]
    active_weight <- active_weight[keep]
    active_weight <- active_weight / sum(active_weight)
    current <- as.numeric(crossprod(active_weight, features[selected, , drop = FALSE]))
  }
  ord <- order(selected)
  list(index = selected[ord], weight = active_weight[ord])
}

build_particle_mis_cache <- function(atlas, theta, model) {
  atlas <- validate_local_atlas(atlas)
  model <- normalize_population_model(model)
  theta <- .as_hyper_matrix(theta, model$hyper_names, model$hyper_dim)
  active <- .local_atlas_active_charts(atlas)
  if (!length(active)) stop("Selected atlas has no active charts.")
  chart_ids <- names(active)
  eta <- rep(1 / length(active), length(active))
  log_eta <- log(eta)
  anchors <- do.call(rbind, lapply(active, function(chart) chart$theta_anchor))
  colnames(anchors) <- model$hyper_names
  logZ <- vapply(active, function(chart) as.numeric(chart$logZ_abs), numeric(1))
  alpha <- do.call(rbind, lapply(active, function(chart) chart$alpha_particles))
  chart_index <- rep(seq_along(active), vapply(active, function(chart) nrow(chart$alpha_particles), integer(1)))
  sample_log_weight <- unlist(Map(function(chart, eta_s) {
    log(eta_s) + log(pmax(.local_chart_normalize_weights(chart$weights, nrow(chart$alpha_particles)), .Machine$double.eps))
  }, active, eta), use.names = FALSE)
  sample_weight <- exp(sample_log_weight - logsumexp(sample_log_weight))

  anchor_logp <- population_model_log_alpha_given_theta_many(model, alpha = alpha, theta = anchors)
  den_terms <- sweep(anchor_logp, 1L, log_eta - logZ, "+")
  log_den <- as.numeric(matrixStats::colLogSumExps(den_terms))
  theta_logp <- population_model_log_alpha_given_theta_many(model, alpha = alpha, theta = theta)
  log_h <- sweep(theta_logp, 2L, log_den, "-")
  full_log_terms <- sweep(log_h, 2L, sample_log_weight, "+")
  full_log_m <- .rowLogSumExp(full_log_terms)
  max_terms <- matrixStats::rowMaxs(full_log_terms)
  contribution <- exp(sweep(full_log_terms, 1L, max_terms, "-"))
  contribution <- sweep(contribution, 1L, rowSums(contribution), "/")
  full_ess <- 1 / rowSums(contribution * contribution)
  full_psis_k <- vapply(seq_len(nrow(log_h)), function(j) .local_chart_psis_k(log_h[j, ]), numeric(1))

  list(
    active = active,
    chart_ids = chart_ids,
    alpha = alpha,
    chart_index = chart_index,
    sample_weight = sample_weight,
    sample_log_weight = sample_log_weight,
    log_h = log_h,
    full_log_m = as.numeric(full_log_m),
    full_ess = as.numeric(full_ess),
    full_ess_frac = as.numeric(full_ess / ncol(log_h)),
    full_psis_k = as.numeric(full_psis_k)
  )
}

build_compression_features <- function(cache,
                                       theta_rows,
                                       evidence_weight = 1,
                                       moment_weight = 0.05,
                                       chart_weight = 0.05,
                                       include_moments = TRUE,
                                       include_chart = TRUE) {
  theta_rows <- as.integer(theta_rows)
  evidence <- t(exp(sweep(cache$log_h[theta_rows, , drop = FALSE], 1L, cache$full_log_m[theta_rows], "-")))
  raw <- evidence
  target <- rep(1, length(theta_rows))
  scale <- rep(1, length(theta_rows))
  penalty_weight <- rep(as.numeric(evidence_weight), length(theta_rows))

  if (isTRUE(include_moments)) {
    alpha <- as.matrix(cache$alpha)
    w <- cache$sample_weight
    alpha_mean <- colSums(w * alpha)
    alpha_centered <- sweep(alpha, 2L, alpha_mean, "-")
    alpha_sd <- sqrt(colSums(w * alpha_centered * alpha_centered))
    alpha_sd <- pmax(alpha_sd, 1e-8)
    z <- sweep(alpha_centered, 2L, alpha_sd, "/")
    moment_raw <- cbind(z, z * z)
    moment_target <- colSums(w * moment_raw)
    moment_scale <- sqrt(pmax(colSums(w * sweep(moment_raw, 2L, moment_target, "-")^2), 1e-8))
    raw <- cbind(raw, moment_raw)
    target <- c(target, moment_target)
    scale <- c(scale, moment_scale)
    penalty_weight <- c(penalty_weight, rep(as.numeric(moment_weight), ncol(moment_raw)))
  }

  if (isTRUE(include_chart)) {
    chart <- matrix(0, nrow = nrow(cache$alpha), ncol = length(cache$chart_ids))
    chart[cbind(seq_len(nrow(chart)), cache$chart_index)] <- 1
    chart_target <- colSums(cache$sample_weight * chart)
    chart_scale <- sqrt(pmax(chart_target * (1 - chart_target), 1e-8))
    raw <- cbind(raw, chart)
    target <- c(target, chart_target)
    scale <- c(scale, chart_scale)
    penalty_weight <- c(penalty_weight, rep(as.numeric(chart_weight), ncol(chart)))
  }

  scaled <- sweep(raw, 2L, scale, "/")
  scaled <- sweep(scaled, 2L, sqrt(pmax(penalty_weight, 0)), "*")
  target_scaled <- target / scale * sqrt(pmax(penalty_weight, 0))
  list(features = scaled, target = target_scaled)
}

evaluate_subset <- function(log_h, index, weight) {
  log_weight <- log(pmax(as.numeric(weight), .Machine$double.eps))
  log_terms <- sweep(log_h[, index, drop = FALSE], 2L, log_weight, "+")
  log_m <- .rowLogSumExp(log_terms)
  max_terms <- matrixStats::rowMaxs(log_terms)
  contribution <- exp(sweep(log_terms, 1L, max_terms, "-"))
  contribution <- sweep(contribution, 1L, rowSums(contribution), "/")
  ess <- 1 / rowSums(contribution * contribution)
  psis_k <- vapply(seq_len(nrow(log_h)), function(j) .local_chart_psis_k(log_h[j, index]), numeric(1))
  data.frame(
    compressed_log_marginal = as.numeric(log_m),
    compressed_ess = as.numeric(ess),
    compressed_ess_frac = as.numeric(ess / length(index)),
    compressed_psis_k = as.numeric(psis_k),
    check.names = FALSE
  )
}

exact_anchor_logZ <- function(atlas, theta, model, tolerance = 1e-8) {
  atlas <- validate_local_atlas(atlas)
  theta <- .as_hyper_matrix(theta, model$hyper_names, model$hyper_dim)
  active <- .local_atlas_active_charts(atlas)
  anchors <- do.call(rbind, lapply(active, function(chart) chart$theta_anchor))
  colnames(anchors) <- model$hyper_names
  direct <- rep(NA_real_, nrow(theta))
  direct_se <- rep(NA_real_, nrow(theta))
  graph <- rep(NA_real_, nrow(theta))
  graph_se <- rep(NA_real_, nrow(theta))
  out_chart <- rep(NA_character_, nrow(theta))
  for (j in seq_len(nrow(theta))) {
    distance <- apply(abs(sweep(anchors, 2L, theta[j, ], "-")), 1L, max)
    hit <- which(distance <= tolerance)
    if (length(hit)) {
      chart <- active[[hit[1L]]]
      direct[j] <- as.numeric(chart$diagnostics$pre_graph_logZ_abs %||% chart$logZ_abs)
      direct_se[j] <- as.numeric(chart$diagnostics$pre_graph_logZ_abs_se %||% chart$logZ_abs_se)
      graph[j] <- as.numeric(chart$logZ_abs)
      graph_se[j] <- as.numeric(chart$logZ_abs_se)
      out_chart[j] <- chart$chart_id %||% names(active)[hit[1L]]
    }
  }
  data.frame(
    theta_row = seq_len(nrow(theta)),
    direct_anchor_log_marginal = direct,
    direct_anchor_logZ_se = direct_se,
    graph_anchor_log_marginal = graph,
    graph_anchor_logZ_se = graph_se,
    exact_anchor_log_marginal = direct,
    exact_anchor_logZ_se = direct_se,
    exact_anchor_chart = out_chart,
    check.names = FALSE
  )
}

summarize_method <- function(rows) {
  split_keys <- split(rows, list(rows$method, rows$K, rows$replicate), drop = TRUE)
  do.call(rbind, lapply(split_keys, function(df) {
    gold <- df[df$is_gold, , drop = FALSE]
    data.frame(
      method = df$method[1L],
      K = df$K[1L],
      replicate = df$replicate[1L],
      total_particles = df$total_particles[1L],
      selected_unique_particles = df$selected_unique_particles[1L],
      quadrature_weight_ess = df$quadrature_weight_ess[1L],
      active_compression_ratio = df$total_particles[1L] / df$selected_unique_particles[1L],
      compression_rmse_vs_full = finite_rmse(df$compression_error),
      compression_centered_rmse_vs_full = finite_rmse(center_finite(df$compression_error)),
      gold_rmse_vs_nested = finite_rmse(gold$compressed_error_vs_gold),
      gold_centered_rmse_vs_nested = finite_rmse(center_finite(gold$compressed_error_vs_gold)),
      full_gold_rmse_vs_nested = finite_rmse(gold$full_error_vs_gold),
      full_gold_centered_rmse_vs_nested = finite_rmse(center_finite(gold$full_error_vs_gold)),
      exact_gold_rmse_vs_nested = finite_rmse(gold$exact_error_vs_gold),
      exact_gold_centered_rmse_vs_nested = finite_rmse(center_finite(gold$exact_error_vs_gold)),
      graph_gold_rmse_vs_nested = finite_rmse(gold$graph_error_vs_gold),
      graph_gold_centered_rmse_vs_nested = finite_rmse(center_finite(gold$graph_error_vs_gold)),
      median_compressed_ess_frac = stats::median(df$compressed_ess_frac, na.rm = TRUE),
      median_compressed_psis_k = stats::median(df$compressed_psis_k, na.rm = TRUE),
      check.names = FALSE
    )
  }))
}

plot_results <- function(rows, summary, plot_file, gold_rows) {
  grDevices::png(plot_file, width = 2200, height = 1450, res = 150, bg = "white")
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  }, add = TRUE)
  graphics::par(mfrow = c(2, 3), mar = c(4.3, 4.3, 3, 1))
  base <- rows[rows$method == "full_particle_mis" & rows$replicate == 1L, , drop = FALSE]
  base_gold <- base[base$is_gold, , drop = FALSE]
  if (nrow(base_gold)) {
    y <- rbind(
      data.frame(method = "direct anchor", theta_row = base_gold$theta_row, error = base_gold$exact_error_vs_gold),
      data.frame(method = "graph anchor", theta_row = base_gold$theta_row, error = base_gold$graph_error_vs_gold),
      data.frame(method = "full particle-MIS", theta_row = base_gold$theta_row, error = base_gold$full_error_vs_gold)
    )
    ylim <- range(c(y$error, 0), na.rm = TRUE)
    graphics::plot(
      y$theta_row,
      y$error,
      type = "n",
      xlab = "gold theta row",
      ylab = "log m - 6144 nested gold",
      main = "Expensive local surface vs nested gold",
      ylim = ylim
    )
    graphics::abline(h = 0, col = "grey45")
    graphics::points(y$theta_row[y$method == "direct anchor"], y$error[y$method == "direct anchor"], pch = 19, col = "black", cex = 1.2)
    graphics::points(y$theta_row[y$method == "graph anchor"], y$error[y$method == "graph anchor"], pch = 15, col = "firebrick3", cex = 1.2)
    graphics::points(y$theta_row[y$method == "full particle-MIS"], y$error[y$method == "full particle-MIS"], pch = 17, col = "steelblue4", cex = 1.2)
    graphics::segments(
      gold_rows$theta_row,
      -gold_rows$gold_replicate_sd,
      gold_rows$theta_row,
      gold_rows$gold_replicate_sd,
      col = "grey55"
    )
    graphics::legend("topright", legend = c("direct chart logZ", "graph chart logZ", "full particle-MIS", "gold replicate SD"),
                     col = c("black", "firebrick3", "steelblue4", "grey55"), pch = c(19, 15, 17, NA), lty = c(NA, NA, NA, 1), bty = "n", cex = 0.75)
  } else {
    graphics::plot.new()
    graphics::title("Expensive local surface vs nested gold")
  }

  methods <- setdiff(unique(summary$method), "full_particle_mis")
  cols <- setNames(grDevices::hcl.colors(max(length(methods), 1L), "Dark 3"), methods)
  metric_plot <- function(metric, ylab, main) {
    ok <- summary$method %in% methods & is.finite(summary[[metric]])
    if (!any(ok)) {
      graphics::plot.new()
      graphics::title(main)
      return()
    }
    graphics::plot(NA, xlim = range(summary$K[ok]), ylim = range(summary[[metric]][ok]), log = "x",
                   xlab = "compressed particles K", ylab = ylab, main = main)
    for (method in methods) {
      df <- summary[summary$method == method, , drop = FALSE]
      agg <- aggregate(df[[metric]], list(K = df$K), mean, na.rm = TRUE)
      graphics::lines(agg$K, agg$x, type = "b", pch = 19, lwd = 2, col = cols[method])
    }
    graphics::legend("topright", legend = methods, col = cols[methods], pch = 19, lwd = 2, bty = "n", cex = 0.72)
  }
  metric_plot("compression_rmse_vs_full", "RMSE log m", "Compression error vs full particle-MIS")
  metric_plot("compression_centered_rmse_vs_full", "centered RMSE log m", "Shape error vs full particle-MIS")
  metric_plot("gold_centered_rmse_vs_nested", "centered RMSE log m", "Compressed gold shape vs nested")
  metric_plot("active_compression_ratio", "raw / active particles", "Active compression ratio")

  best <- summary[summary$method == "sparse_quadrature" & is.finite(summary$compression_centered_rmse_vs_full), , drop = FALSE]
  if (nrow(best)) {
    best <- best[order(best$compression_centered_rmse_vs_full, best$K), , drop = FALSE][1L, ]
    df <- rows[rows$method == best$method & rows$K == best$K & rows$replicate == best$replicate, , drop = FALSE]
    ylim <- range(c(df$compression_error, 0), na.rm = TRUE)
    graphics::plot(df$theta_row, df$compression_error, type = "h", lwd = 3,
                   xlab = "theta row", ylab = "compressed - full log m",
                   main = sprintf("Sparse quadrature K=%d", best$K),
                   ylim = ylim)
    graphics::points(df$theta_row, df$compression_error, pch = ifelse(df$is_gold, 21, 19),
                     bg = ifelse(df$is_gold, "gold", "white"))
    graphics::abline(h = 0, col = "grey45")
  } else {
    graphics::plot.new()
    graphics::title("Sparse quadrature")
  }
}

cli_args <- parse_cli_args(commandArgs(trailingOnly = TRUE))
label <- arg_chr(cli_args, "label", "emc_expensive_local_compression_hard464")
seed <- arg_int(cli_args, "seed", 20260526L)
local_id <- arg_chr(cli_args, "local", "464")
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
data_file <- arg_chr(cli_args, "data_file", file.path("benchmarks", "samples", "full_EMC2.RData"))
atlas_particles <- arg_int(cli_args, "atlas_particles", 1024L)
atlas_target_cess <- arg_num(cli_args, "atlas_target_cess", 0.9)
atlas_mcmc_moves <- arg_int(cli_args, "atlas_mcmc_moves", 2L)
atlas_max_steps <- arg_int(cli_args, "atlas_max_steps", 192L)
atlas_max_intermediates <- arg_int(cli_args, "atlas_max_intermediates", 4L)
reuse_atlas_cache <- arg_lgl(cli_args, "reuse_atlas_cache", TRUE)
refresh_atlas_cache <- arg_lgl(cli_args, "refresh_atlas_cache", FALSE)
verbose <- arg_lgl(cli_args, "verbose", FALSE)
k_values <- arg_int_vec(cli_args, "k_values", c(64L, 128L, 256L, 512L))
resample_reps <- arg_int(cli_args, "resample_reps", 4L)
methods <- strsplit(arg_chr(cli_args, "methods", "sparse_quadrature,systematic,top_weight"), ",", fixed = TRUE)[[1L]]
methods <- trimws(methods)
if (length(methods) == 1L && tolower(methods) %in% c("", "none", "full", "full_only")) {
  methods <- character()
}
theta_rows_arg <- arg_chr(cli_args, "theta_rows", "gold")
anchor_theta_rows_arg <- arg_chr(cli_args, "anchor_theta_rows", theta_rows_arg)
out_prefix <- file.path("benchmarks", "results", label)
atlas_cache_file <- arg_chr(cli_args, "atlas_cache_file", paste0(out_prefix, "_atlas_cache.rds"))
rows_csv <- arg_chr(cli_args, "rows_csv", paste0(out_prefix, "_rows.csv"))
summary_csv <- arg_chr(cli_args, "summary_csv", paste0(out_prefix, "_summary.csv"))
plot_file <- arg_chr(cli_args, "plot_file", paste0(out_prefix, "_diagnostics.png"))
results_file <- arg_chr(cli_args, "results_file", paste0(out_prefix, "_results.rds"))

for (path in c(atlas_cache_file, rows_csv, summary_csv, plot_file, results_file)) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}

if (!file.exists(source_results)) stop("Missing source_results: ", source_results)
if (!file.exists(gold_csv)) stop("Missing gold_csv: ", gold_csv)
if (!file.exists(data_file)) stop("Missing data_file: ", data_file)

source("smc_core.R")
source("population_models.R")
source("local_charts.R")
source("utilities.R")
if (any(methods == "sparse_quadrature") && !requireNamespace("quadprog", quietly = TRUE)) {
  stop("quadprog is required for sparse_quadrature.")
}

source_obj <- readRDS(source_results)
theta_design <- source_obj$theta_design
theta_full <- theta_design$theta
metadata_full <- theta_design$metadata
model <- theta_design$population_model
if (is.null(model)) model <- source_obj$factor_sets[[1L]]$population_model
gold_all <- utils::read.csv(gold_csv, check.names = FALSE)
gold_all <- gold_all[as.character(gold_all$local) == as.character(local_id), , drop = FALSE]
parse_theta_rows <- function(value, full_n, gold_rows) {
  if (identical(tolower(value), "gold")) {
    out <- sort(unique(as.integer(gold_rows)))
  } else if (identical(tolower(value), "all")) {
    out <- seq_len(full_n)
  } else {
    out <- as.integer(strsplit(value, ",", fixed = TRUE)[[1L]])
  }
  out <- out[is.finite(out)]
  if (!length(out) || any(out < 1L | out > full_n)) {
    stop("theta row arguments must be 'gold', 'all', or valid row indices in the source theta design.")
  }
  sort(unique(out))
}
selected_theta_rows <- parse_theta_rows(theta_rows_arg, nrow(theta_full), gold_all$theta_row)
selected_anchor_theta_rows <- parse_theta_rows(anchor_theta_rows_arg, nrow(theta_full), gold_all$theta_row)
theta <- theta_full[selected_theta_rows, , drop = FALSE]
metadata <- metadata_full[selected_theta_rows, , drop = FALSE]
metadata$theta_row <- selected_theta_rows
theta_anchor_design <- theta_full[selected_anchor_theta_rows, , drop = FALSE]

load(data_file)
if (!exists("ELP_DDM", inherits = FALSE)) stop("The EMC2 data file must define ELP_DDM.")
emc <- ELP_DDM[[1L]]
data_list <- emc$data
model_factory <- emc$model
alpha_names <- emc$par_names
if (!local_id %in% names(data_list)) stop("Local ", local_id, " not found in EMC data.")
data_i <- data_list[[local_id]]

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

atlas_settings <- list(
  local = local_id,
  source_results = normalizePath(source_results),
  data_file = normalizePath(data_file),
  anchor_theta = theta_anchor_design,
  atlas_particles = as.integer(atlas_particles),
  atlas_target_cess = as.numeric(atlas_target_cess),
  atlas_mcmc_moves = as.integer(atlas_mcmc_moves),
  atlas_max_steps = as.integer(atlas_max_steps),
  atlas_max_intermediates = as.integer(atlas_max_intermediates),
  seed = as.integer(seed)
)

atlas_cache_compatible <- function(cache) {
  if (!is.list(cache) || is.null(cache$atlas) || is.null(cache$settings)) return(FALSE)
  s <- cache$settings
  cached_anchor_theta <- s$anchor_theta %||% s$theta
  identical(as.character(s$local), as.character(atlas_settings$local)) &&
    identical(as.integer(s$atlas_particles), as.integer(atlas_settings$atlas_particles)) &&
    identical(as.integer(s$atlas_mcmc_moves), as.integer(atlas_settings$atlas_mcmc_moves)) &&
    identical(as.integer(s$atlas_max_steps), as.integer(atlas_settings$atlas_max_steps)) &&
    identical(as.integer(s$atlas_max_intermediates), as.integer(atlas_settings$atlas_max_intermediates)) &&
    identical(as.integer(s$seed), as.integer(atlas_settings$seed)) &&
    isTRUE(all.equal(as.numeric(s$atlas_target_cess), as.numeric(atlas_settings$atlas_target_cess), tolerance = 1e-12)) &&
    isTRUE(all.equal(as.matrix(cached_anchor_theta), as.matrix(atlas_settings$anchor_theta), tolerance = 1e-10, check.attributes = FALSE))
}

theta_root <- matrix(apply(theta_anchor_design, 2L, stats::median), nrow = 1L)
colnames(theta_root) <- colnames(theta)

start_time <- Sys.time()
if (isTRUE(reuse_atlas_cache) && !isTRUE(refresh_atlas_cache) && file.exists(atlas_cache_file)) {
  cache <- readRDS(atlas_cache_file)
  if (!atlas_cache_compatible(cache)) {
    stop("Atlas cache is incompatible. Use --refresh_atlas_cache=true or a different --atlas_cache_file.")
  }
  atlas <- cache$atlas
  cat("Loaded expensive local atlas cache:", atlas_cache_file, "\n")
} else {
  cat(sprintf(
    "Fitting expensive local atlas: local=%s | theta anchors=%d | particles=%d\n",
    local_id, nrow(theta_anchor_design), atlas_particles
  ))
  atlas <- build_local_atlas(
    local_id = local_id,
    data_i = data_i,
    loglik_fn = loglik_emc2,
    population_model = model,
    theta_design = theta_anchor_design,
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
      distance_metric = "fisher"
    ),
    n_cores = 1L,
    seed = seed + 1000L,
    verbose = isTRUE(verbose)
  )
  saveRDS(list(atlas = atlas, settings = atlas_settings), atlas_cache_file)
  cat("Saved expensive local atlas cache:", atlas_cache_file, "\n")
}

fit_minutes <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
cache <- build_particle_mis_cache(atlas, theta, model)
exact <- exact_anchor_logZ(atlas, theta, model)
exact$theta_row <- metadata$theta_row

gold <- gold_all
gold_summary <- if (nrow(gold)) {
  aggregate(
    cbind(gold_log_marginal, gold_replicate_sd, gold_replicate_range, gold_mean_path_se) ~
      local + theta_row + theta_label,
    gold,
    function(x) x[which(is.finite(x))[1L]] %||% NA_real_
  )
} else {
  data.frame(local = character(), theta_row = integer(), theta_label = character())
}
gold_rows <- gold_summary

base <- data.frame(
  theta_row = metadata$theta_row,
  theta_label = metadata$theta_label %||% sprintf("theta_%03d", seq_len(nrow(theta))),
  full_log_marginal = cache$full_log_m,
  full_ess_frac = cache$full_ess_frac,
  full_psis_k = cache$full_psis_k,
  check.names = FALSE
)
base <- merge(base, exact, by = "theta_row", all.x = TRUE, sort = FALSE)
gold_key <- paste(gold_summary$theta_row, gold_summary$theta_label, sep = "\r")
base_key <- paste(base$theta_row, base$theta_label, sep = "\r")
gold_match <- match(base_key, gold_key)
for (nm in setdiff(names(gold_summary), c("local", "theta_row", "theta_label"))) {
  base[[nm]] <- gold_summary[[nm]][gold_match]
}
base$is_gold <- is.finite(base$gold_log_marginal)
base$full_error_vs_gold <- base$full_log_marginal - base$gold_log_marginal
base$exact_error_vs_gold <- base$exact_anchor_log_marginal - base$gold_log_marginal
base$graph_error_vs_gold <- base$graph_anchor_log_marginal - base$gold_log_marginal

gold_theta_rows <- base$theta_row[base$is_gold]
train_rows <- which(!base$is_gold)
if (!length(train_rows)) train_rows <- seq_len(nrow(theta))
features <- build_compression_features(cache, train_rows, include_moments = TRUE, include_chart = TRUE)

all_rows <- list()
full_rows <- data.frame(
  method = "full_particle_mis",
  K = nrow(cache$alpha),
  replicate = 1L,
  selected_unique_particles = nrow(cache$alpha),
  quadrature_weight_ess = 1 / sum(cache$sample_weight * cache$sample_weight),
  total_particles = nrow(cache$alpha),
  base,
  compressed_log_marginal = base$full_log_marginal,
  compressed_ess_frac = base$full_ess_frac,
  compressed_psis_k = base$full_psis_k,
  check.names = FALSE
)
full_rows$compression_error <- 0
full_rows$compressed_error_vs_gold <- full_rows$full_error_vs_gold
all_rows[[length(all_rows) + 1L]] <- full_rows

for (K_raw in k_values) {
  K <- min(as.integer(K_raw), nrow(cache$alpha))
  for (method in methods) {
    reps <- if (identical(method, "systematic")) seq_len(resample_reps) else 1L
    for (rep_id in reps) {
      selection <- switch(
        method,
        sparse_quadrature = sparse_quadrature_select(features$features, features$target, K),
        systematic = systematic_resample(cache$sample_weight, K, seed = seed + 1009L * K + rep_id),
        top_weight = top_weight_select(cache$sample_weight, K),
        stop("Unknown compression method: ", method)
      )
      comp <- evaluate_subset(cache$log_h, selection$index, selection$weight)
      rows <- data.frame(
        method = method,
        K = K,
        replicate = rep_id,
        selected_unique_particles = length(selection$index),
        quadrature_weight_ess = 1 / sum(selection$weight * selection$weight),
        total_particles = nrow(cache$alpha),
        base,
        compressed_log_marginal = comp$compressed_log_marginal,
        compressed_ess_frac = comp$compressed_ess_frac,
        compressed_psis_k = comp$compressed_psis_k,
        check.names = FALSE
      )
      rows$compression_error <- rows$compressed_log_marginal - rows$full_log_marginal
      rows$compressed_error_vs_gold <- rows$compressed_log_marginal - rows$gold_log_marginal
      all_rows[[length(all_rows) + 1L]] <- rows
    }
  }
}

rows <- do.call(rbind, all_rows)
summary <- summarize_method(rows)
summary <- summary[order(summary$method, summary$K, summary$replicate), , drop = FALSE]
utils::write.csv(rows, rows_csv, row.names = FALSE)
utils::write.csv(summary, summary_csv, row.names = FALSE)
plot_results(rows, summary, plot_file, gold_rows)
saveRDS(
  list(
    rows = rows,
    summary = summary,
    atlas = atlas,
    settings = list(
      label = label,
      local = local_id,
      source_results = source_results,
      gold_csv = gold_csv,
      atlas_particles = atlas_particles,
      k_values = k_values,
      methods = methods,
      resample_reps = resample_reps,
      theta_rows_arg = theta_rows_arg,
      anchor_theta_rows_arg = anchor_theta_rows_arg,
      selected_theta_rows = selected_theta_rows,
      selected_anchor_theta_rows = selected_anchor_theta_rows,
      train_rows = train_rows,
      gold_theta_rows = gold_theta_rows,
      fit_minutes = fit_minutes
    )
  ),
  results_file
)

cat(sprintf("Saved rows: %s\n", rows_csv))
cat(sprintf("Saved summary: %s\n", summary_csv))
cat(sprintf("Saved plot: %s\n", plot_file))
cat(sprintf("Saved results: %s\n", results_file))
cat("\nExpensive local compression summary:\n")
print(summary[, c(
  "method",
  "K",
  "replicate",
  "total_particles",
  "selected_unique_particles",
  "active_compression_ratio",
  "compression_centered_rmse_vs_full",
  "gold_centered_rmse_vs_nested",
  "full_gold_centered_rmse_vs_nested",
  "exact_gold_centered_rmse_vs_nested",
  "graph_gold_centered_rmse_vs_nested"
)], row.names = FALSE)
