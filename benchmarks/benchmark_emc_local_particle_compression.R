#!/usr/bin/env Rscript

rm(list = ls())

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[1L]))
} else {
  normalizePath("benchmarks/benchmark_emc_local_particle_compression.R")
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

arg_chr <- function(args, key, default = NULL) {
  val <- args[[key]]
  if (is.null(val) || !nzchar(val)) default else as.character(val)
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

center_finite <- function(x) {
  ok <- is.finite(x)
  out <- rep(NA_real_, length(x))
  if (any(ok)) out[ok] <- x[ok] - mean(x[ok])
  out
}

logmeanexp <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  m <- max(x)
  m + log(mean(exp(x - m)))
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

feature_herding_select <- function(features, target, K) {
  features <- as.matrix(features)
  target <- as.numeric(target)
  K <- min(as.integer(K), nrow(features))
  selected <- integer(K)
  available <- rep(TRUE, nrow(features))
  current <- rep(0, ncol(features))
  for (k in seq_len(K)) {
    direction <- target - current
    scores <- as.numeric(features %*% direction)
    scores[!available] <- -Inf
    pick <- which.max(scores)
    selected[k] <- pick
    available[pick] <- FALSE
    current <- current + (features[pick, ] - current) / k
  }
  list(index = selected, weight = rep(1 / K, K))
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
  if (is.null(sol) || any(!is.finite(sol))) {
    sol <- rep(1 / K, K)
  }
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
  for (k in seq_len(K)) {
    residual <- target - current
    scores <- as.numeric(features %*% residual)
    if (length(selected)) {
      scores[selected] <- scores[selected] - 1e-12
    }
    pick <- which.max(scores)
    selected <- unique(c(selected, pick))
    active_features <- features[selected, , drop = FALSE]
    active_weight <- fit_simplex_weights(active_features, target, ridge = ridge)
    keep <- active_weight > 1e-12
    selected <- selected[keep]
    active_weight <- active_weight[keep]
    active_weight <- active_weight / sum(active_weight)
    current <- as.numeric(crossprod(active_weight, features[selected, , drop = FALSE]))
  }
  ord <- order(selected)
  list(index = selected[ord], weight = active_weight[ord])
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
  labels <- paste0("evidence_", theta_rows)

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
    labels <- c(labels, paste0("alpha_moment_", seq_len(ncol(moment_raw))))
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
    labels <- c(labels, paste0("chart_", cache$chart_ids))
  }

  scaled <- sweep(raw, 2L, scale, "/")
  scaled <- sweep(scaled, 2L, sqrt(pmax(penalty_weight, 0)), "*")
  target_scaled <- target / scale * sqrt(pmax(penalty_weight, 0))
  list(features = scaled, target = target_scaled, labels = labels)
}

evaluate_subset <- function(log_h, index, weight) {
  log_weight <- log(pmax(as.numeric(weight), .Machine$double.eps))
  log_terms <- sweep(log_h[, index, drop = FALSE], 2L, log_weight, "+")
  log_m <- .rowLogSumExp(log_terms)
  max_terms <- matrixStats::rowMaxs(log_terms)
  contribution <- exp(sweep(log_terms, 1L, max_terms, "-"))
  contribution <- sweep(contribution, 1L, rowSums(contribution), "/")
  ess <- 1 / rowSums(contribution * contribution)
  psis_k <- vapply(seq_len(nrow(log_h)), function(j) {
    .local_chart_psis_k(log_h[j, index])
  }, numeric(1))
  data.frame(
    compressed_log_marginal = as.numeric(log_m),
    compressed_ess = as.numeric(ess),
    compressed_ess_frac = as.numeric(ess / length(index)),
    compressed_psis_k = as.numeric(psis_k),
    check.names = FALSE
  )
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

summarize_rows <- function(rows) {
  split_keys <- split(rows, list(rows$member, rows$method, rows$K, rows$replicate), drop = TRUE)
  do.call(rbind, lapply(split_keys, function(df) {
    gold <- df[df$is_gold, , drop = FALSE]
    holdout <- df[df$is_holdout, , drop = FALSE]
    data.frame(
      member = df$member[1L],
      method = df$method[1L],
      K = df$K[1L],
      replicate = df$replicate[1L],
      total_particles = df$total_particles[1L],
      selected_unique_particles = df$selected_unique_particles[1L],
      quadrature_weight_ess = df$quadrature_weight_ess[1L],
      compression_ratio = df$total_particles[1L] / df$K[1L],
      active_compression_ratio = df$total_particles[1L] / df$selected_unique_particles[1L],
      train_rmse_vs_full = finite_rmse(df$compression_error[df$is_train]),
      holdout_rmse_vs_full = finite_rmse(holdout$compression_error),
      all_rmse_vs_full = finite_rmse(df$compression_error),
      all_centered_rmse_vs_full = finite_rmse(center_finite(df$compression_error)),
      max_abs_error_vs_full = finite_max_abs(df$compression_error),
      gold_rmse_vs_nested = finite_rmse(gold$compressed_error_vs_gold),
      gold_centered_rmse_vs_nested = finite_rmse(center_finite(gold$compressed_error_vs_gold)),
      full_gold_rmse_vs_nested = finite_rmse(gold$full_error_vs_gold),
      full_gold_centered_rmse_vs_nested = finite_rmse(center_finite(gold$full_error_vs_gold)),
      median_compressed_ess_frac = stats::median(df$compressed_ess_frac, na.rm = TRUE),
      median_compressed_psis_k = stats::median(df$compressed_psis_k, na.rm = TRUE),
      check.names = FALSE
    )
  }))
}

plot_results <- function(summary, rows, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  grDevices::png(file, width = 2200, height = 1500, res = 150, bg = "white")
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    grDevices::dev.off()
  }, add = TRUE)
  par(mfrow = c(2, 3), mar = c(4.2, 4.2, 3, 1))
  methods <- unique(summary$method)
  cols <- setNames(grDevices::hcl.colors(length(methods), "Dark 3"), methods)
  plot_metric <- function(metric, ylab, main) {
    ok <- is.finite(summary[[metric]])
    if (!any(ok)) {
      plot.new()
      title(main)
      text(0.5, 0.5, "not available")
      return()
    }
    ylim <- range(summary[[metric]][ok])
    plot(NA, xlim = range(summary$K[ok]), ylim = ylim, log = "x",
         xlab = "compressed particles K", ylab = ylab, main = main)
    for (method in methods) {
      df <- summary[summary$method == method, , drop = FALSE]
      agg <- aggregate(df[[metric]], list(K = df$K), mean, na.rm = TRUE)
      lines(agg$K, agg$x, type = "b", pch = 19, col = cols[method], lwd = 2)
    }
    legend("topright", legend = methods, col = cols[methods], lwd = 2, pch = 19, bty = "n", cex = 0.72)
  }
  plot_metric("holdout_rmse_vs_full", "RMSE log m", "Held-out compression error vs full particle-MIS")
  plot_metric("all_centered_rmse_vs_full", "centered RMSE log m", "Shape error vs full particle-MIS")
  plot_metric("gold_rmse_vs_nested", "RMSE log m", "Gold nested-SMC points")
  plot_metric("median_compressed_ess_frac", "median ESS fraction", "Compressed particle-MIS ESS")
  plot_metric("median_compressed_psis_k", "median PSIS k", "Compressed particle-MIS PSIS")
  best <- summary[
    summary$method == "sparse_quadrature" & is.finite(summary$holdout_rmse_vs_full),
    ,
    drop = FALSE
  ]
  if (nrow(best)) {
    best <- best[order(best$holdout_rmse_vs_full, best$K), , drop = FALSE][1L, ]
    df <- rows[
      rows$member == best$member &
        rows$method == best$method &
        rows$K == best$K &
        rows$replicate == best$replicate,
      ,
      drop = FALSE
    ]
    ylim <- range(c(df$compression_error, 0), na.rm = TRUE)
    plot(df$theta_row, df$compression_error, type = "h", lwd = 3,
         xlab = "theta row", ylab = "compressed - full log m",
         main = sprintf("Sparse quadrature held-out: K=%d", best$K),
         ylim = ylim)
    points(df$theta_row, df$compression_error, pch = ifelse(df$is_holdout, 21, 19),
           bg = ifelse(df$is_gold, "gold", "white"))
    abline(h = 0, col = "grey35")
    legend("topright", legend = c("filled: gold row", "open: holdout row"), pch = c(21, 21),
           pt.bg = c("gold", "white"), bty = "n", cex = 0.72)
  } else {
    plot.new()
    title("Best held-out")
  }
}

cli_args <- parse_cli_args(commandArgs(trailingOnly = TRUE))
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
label <- arg_chr(cli_args, "label", "emc_local_particle_compression_hard464")
local_id <- arg_chr(cli_args, "local", "464")
members <- arg_chr_vec(cli_args, "members", "post_cal")
k_values <- arg_int_vec(cli_args, "k_values", c(32L, 64L, 128L, 256L, 512L))
methods <- arg_chr_vec(cli_args, "methods", c("sparse_quadrature", "sparse_quadrature_evidence", "systematic", "top_weight"))
resample_reps <- arg_int(cli_args, "resample_reps", 8L)
holdout_rows <- arg_int_vec(cli_args, "holdout_rows", c(3L, 5L, 8L))
seed <- arg_int(cli_args, "seed", 20260526L)
out_prefix <- file.path("benchmarks", "results", label)
rows_csv <- arg_chr(cli_args, "rows_csv", paste0(out_prefix, "_rows.csv"))
summary_csv <- arg_chr(cli_args, "summary_csv", paste0(out_prefix, "_summary.csv"))
plot_file <- arg_chr(cli_args, "plot_file", paste0(out_prefix, "_diagnostics.png"))
results_file <- arg_chr(cli_args, "results_file", paste0(out_prefix, "_results.rds"))

if (!file.exists(source_results)) stop("Missing source_results: ", source_results)
if (!file.exists(gold_csv)) stop("Missing gold_csv: ", gold_csv)

source("local_charts.R")
if (any(methods %in% c("sparse_quadrature", "sparse_quadrature_evidence")) &&
    !requireNamespace("quadprog", quietly = TRUE)) {
  stop("quadprog is required for sparse_quadrature methods.")
}

source_obj <- readRDS(source_results)
if (is.null(source_obj$factor_sets)) stop("source_results must contain factor_sets.")
theta <- source_obj$theta_design$theta
theta_metadata <- source_obj$theta_design$metadata
model <- source_obj$theta_design$population_model
if (is.null(model)) {
  first_member <- names(source_obj$factor_sets)[1L]
  model <- source_obj$factor_sets[[first_member]]$population_model
}
members <- intersect(members, names(source_obj$factor_sets))
if (!length(members)) stop("No requested members are available in source_results.")
gold <- utils::read.csv(gold_csv, check.names = FALSE)
gold <- gold[as.character(gold$local) == as.character(local_id), , drop = FALSE]
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

if (!all(holdout_rows %in% seq_len(nrow(theta)))) {
  stop("holdout_rows must be valid theta row indices.")
}
train_rows <- setdiff(seq_len(nrow(theta)), holdout_rows)

cat(sprintf("Particle-compression benchmark: local=%s | members=%s\n", local_id, paste(members, collapse = ",")))
cat(sprintf("Theta rows: %d total | %d train | %d holdout\n", nrow(theta), length(train_rows), length(holdout_rows)))
cat(sprintf("K values: %s | methods=%s\n", paste(k_values, collapse = ","), paste(methods, collapse = ",")))

all_rows <- list()
for (member in members) {
  fs <- source_obj$factor_sets[[member]]
  if (!local_id %in% names(fs$atlases)) {
    stop("Local ", local_id, " is not present in factor set member ", member, ".")
  }
  atlas <- fs$atlases[[local_id]]
  cache <- build_particle_mis_cache(atlas, theta, model)
  total_particles <- nrow(cache$alpha)
  full_rows <- data.frame(
    theta_row = seq_len(nrow(theta)),
    theta_label = theta_metadata$theta_label %||% sprintf("theta_%03d", seq_len(nrow(theta))),
    full_log_marginal = cache$full_log_m,
    full_ess_frac = cache$full_ess_frac,
    full_psis_k = cache$full_psis_k,
    check.names = FALSE
  )
  gold_key <- paste(gold_summary$theta_row, gold_summary$theta_label, sep = "\r")
  full_key <- paste(full_rows$theta_row, full_rows$theta_label, sep = "\r")
  gold_match <- match(full_key, gold_key)
  for (nm in setdiff(names(gold_summary), c("local", "theta_row", "theta_label"))) {
    full_rows[[nm]] <- gold_summary[[nm]][gold_match]
  }
  full_rows$full_error_vs_gold <- full_rows$full_log_marginal - full_rows$gold_log_marginal
  log_h_train <- cache$log_h[train_rows, , drop = FALSE]
  full_log_m_train <- cache$full_log_m[train_rows]
  feature <- t(exp(sweep(log_h_train, 1L, full_log_m_train, "-")))
  target <- colSums(cache$sample_weight * feature)
  quadrature_features <- build_compression_features(
    cache = cache,
    theta_rows = train_rows,
    include_moments = TRUE,
    include_chart = TRUE
  )
  evidence_quadrature_features <- build_compression_features(
    cache = cache,
    theta_rows = train_rows,
    include_moments = FALSE,
    include_chart = FALSE
  )

  for (K_raw in k_values) {
    K <- min(as.integer(K_raw), total_particles)
    for (method in methods) {
      reps <- if (identical(method, "systematic")) seq_len(resample_reps) else 1L
      for (rep_id in reps) {
        selection <- switch(
          method,
          systematic = systematic_resample(cache$sample_weight, K, seed = seed + 1009L * K + rep_id),
          feature_herding = feature_herding_select(feature, target, K),
          sparse_quadrature = sparse_quadrature_select(
            quadrature_features$features,
            quadrature_features$target,
            K
          ),
          sparse_quadrature_evidence = sparse_quadrature_select(
            evidence_quadrature_features$features,
            evidence_quadrature_features$target,
            K
          ),
          top_weight = top_weight_select(cache$sample_weight, K),
          stop("Unknown compression method: ", method)
        )
        comp <- evaluate_subset(cache$log_h, selection$index, selection$weight)
        rows <- data.frame(
          member = member,
          local = local_id,
          method = method,
          K = K,
          replicate = rep_id,
          selected_unique_particles = length(selection$index),
          quadrature_weight_ess = 1 / sum(selection$weight * selection$weight),
          total_particles = total_particles,
          theta_row = seq_len(nrow(theta)),
          theta_label = full_rows$theta_label,
          is_train = seq_len(nrow(theta)) %in% train_rows,
          is_holdout = seq_len(nrow(theta)) %in% holdout_rows,
          full_log_marginal = full_rows$full_log_marginal,
          full_ess_frac = full_rows$full_ess_frac,
          full_psis_k = full_rows$full_psis_k,
          compressed_log_marginal = comp$compressed_log_marginal,
          compressed_ess_frac = comp$compressed_ess_frac,
          compressed_psis_k = comp$compressed_psis_k,
          gold_log_marginal = full_rows$gold_log_marginal,
          gold_replicate_sd = full_rows$gold_replicate_sd,
          gold_mean_path_se = full_rows$gold_mean_path_se,
          check.names = FALSE
        )
        rows$compression_error <- rows$compressed_log_marginal - rows$full_log_marginal
        rows$full_error_vs_gold <- rows$full_log_marginal - rows$gold_log_marginal
        rows$compressed_error_vs_gold <- rows$compressed_log_marginal - rows$gold_log_marginal
        rows$is_gold <- is.finite(rows$gold_log_marginal)
        all_rows[[length(all_rows) + 1L]] <- rows
      }
    }
  }
}

rows <- do.call(rbind, all_rows)
summary <- summarize_rows(rows)
summary <- summary[order(summary$member, summary$method, summary$K, summary$replicate), , drop = FALSE]
dir.create(dirname(rows_csv), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(rows, rows_csv, row.names = FALSE)
utils::write.csv(summary, summary_csv, row.names = FALSE)
plot_results(summary, rows, plot_file)
saveRDS(
  list(
    rows = rows,
    summary = summary,
    source_results = source_results,
    gold_csv = gold_csv,
    settings = list(
      label = label,
      local = local_id,
      members = members,
      k_values = k_values,
      methods = methods,
      resample_reps = resample_reps,
      train_rows = train_rows,
      holdout_rows = holdout_rows,
      seed = seed
    )
  ),
  results_file
)

cat(sprintf("Saved rows: %s\n", rows_csv))
cat(sprintf("Saved summary: %s\n", summary_csv))
cat(sprintf("Saved plot: %s\n", plot_file))
cat(sprintf("Saved results: %s\n", results_file))
cat("\nCompression summary:\n")
print(summary[, c(
  "member",
  "method",
  "K",
  "replicate",
  "compression_ratio",
  "active_compression_ratio",
  "selected_unique_particles",
  "quadrature_weight_ess",
  "holdout_rmse_vs_full",
  "all_centered_rmse_vs_full",
  "gold_rmse_vs_nested",
  "gold_centered_rmse_vs_nested",
  "median_compressed_ess_frac",
  "median_compressed_psis_k"
)], row.names = FALSE)
