#!/usr/bin/env Rscript

rm(list = ls())

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[1L]))
} else {
  normalizePath("benchmarks/run_emc_total_error_panel_pilot.R")
}
repo_dir <- dirname(dirname(script_path))
setwd(repo_dir)

suppressPackageStartupMessages({
  library(parallel)
  library(EMC2)
})

source("hierarchical_framework.R")
source("utilities.R")

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

arg_lgl <- function(args, key, default = FALSE) {
  val <- args[[key]]
  if (is.null(val) || !nzchar(val)) return(isTRUE(default))
  tolower(as.character(val)) %in% c("1", "true", "t", "yes", "y")
}

announce <- function(title, detail = NULL) {
  cat(sprintf("\n=== %s ===\n", title))
  if (!is.null(detail) && nzchar(detail)) cat(detail, "\n")
  flush.console()
}

save_stage <- function(object, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(object, file)
  cat(sprintf("Saved: %s\n", file))
  flush.console()
  invisible(file)
}

write_csv <- function(object, file) {
  if (!is.data.frame(object)) return(invisible(FALSE))
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(object, file, row.names = FALSE)
  cat(sprintf("Saved: %s\n", file))
  flush.console()
  invisible(TRUE)
}

weighted_mean_matrix <- function(x, w) {
  w <- pmax(as.numeric(w), 0)
  w <- w / sum(w)
  colSums(w * as.matrix(x))
}

nearest_theta_row <- function(theta, target, center, whitening, excluded = integer()) {
  z <- sweep(theta, 2L, center, "-") %*% whitening
  target_z <- matrix(target - center, nrow = 1L) %*% whitening
  d2 <- rowSums(sweep(z, 2L, target_z[1L, ], "-")^2)
  d2[as.integer(excluded)] <- Inf
  which.min(d2)
}

pick_weighted_quantile_row <- function(theta, w, axis, prob, excluded = integer()) {
  q <- .local_atlas_weighted_quantile(theta[, axis], w, probs = prob)
  score <- abs(theta[, axis] - q)
  score[as.integer(excluded)] <- Inf
  which.min(score)
}

select_panel_rows <- function(fit,
                              factor_set,
                              reference_draws,
                              n_train = 8L,
                              n_holdout = 4L,
                              max_axes = 4L,
                              n_draws = 3000L,
                              seed = 1L) {
  model <- factor_set$population_model
  theta <- .as_hyper_matrix(fit$theta, model$hyper_names, model$hyper_dim)
  w <- .local_chart_normalize_weights(fit$w, nrow(theta))
  posterior_geometry <- .local_evidence_weighted_center_cov(theta, w, ridge = 1e-8)
  workflow_draws <- local_atlas_draws_from_fit(
    fit,
    population_model = model,
    n_draws = min(as.integer(n_draws), nrow(theta)),
    seed = as.integer(seed)
  )
  comparison <- local_atlas_compare_posterior_draws(reference_draws, workflow_draws)
  comparison$axis_score <- pmax(abs(comparison$standardized_mean_error), comparison$shape_error, na.rm = TRUE)
  comparison <- comparison[comparison$parameter %in% model$hyper_names, , drop = FALSE]
  comparison <- comparison[order(-comparison$axis_score), , drop = FALSE]
  axes <- head(as.character(comparison$parameter), as.integer(max_axes))
  axis_idx <- match(axes, model$hyper_names)
  axis_idx <- axis_idx[is.finite(axis_idx)]

  selected <- integer()
  roles <- character()
  add <- function(row, role) {
    row <- as.integer(row)
    if (!is.finite(row) || row < 1L || row > nrow(theta) || row %in% selected) return(invisible(FALSE))
    selected <<- c(selected, row)
    roles <<- c(roles, role)
    invisible(TRUE)
  }

  center <- weighted_mean_matrix(theta, w)
  add(nearest_theta_row(theta, center, posterior_geometry$center, posterior_geometry$whitening), "posterior_center")
  for (row in order(-w)) {
    if (length(selected) >= 2L) break
    add(row, "high_weight")
  }
  for (axis in axis_idx) {
    add(pick_weighted_quantile_row(theta, w, axis, 0.05, selected),
        paste0("mismatch_axis_low:", model$hyper_names[axis]))
    add(pick_weighted_quantile_row(theta, w, axis, 0.95, selected),
        paste0("mismatch_axis_high:", model$hyper_names[axis]))
  }
  for (axis in axis_idx) {
    add(pick_weighted_quantile_row(theta, w, axis, 0.15, selected),
        paste0("holdout_axis_low:", model$hyper_names[axis]))
    add(pick_weighted_quantile_row(theta, w, axis, 0.85, selected),
        paste0("holdout_axis_high:", model$hyper_names[axis]))
  }
  for (row in order(-w)) {
    if (length(selected) >= n_train + n_holdout) break
    add(row, "high_weight_fill")
  }

  if (length(selected) < n_train + n_holdout) {
    stop("Could not construct the requested theta panel.")
  }
  selected <- selected[seq_len(n_train + n_holdout)]
  roles <- roles[seq_along(selected)]
  probe_set <- c(rep("repair", n_train), rep("holdout", n_holdout))
  panel_weight <- w[selected]
  panel_weight <- pmax(panel_weight, 0)
  panel_weight <- 0.5 * panel_weight / sum(panel_weight) + 0.5 / length(panel_weight)
  theta_panel <- theta[selected, , drop = FALSE]
  metadata <- data.frame(
    theta_row = seq_along(selected),
    fit_theta_row = selected,
    theta_id = sprintf("panel_theta_%03d", seq_along(selected)),
    theta_source = roles,
    probe_set = probe_set,
    fit_weight = w[selected],
    theta_weight = panel_weight,
    check.names = FALSE
  )
  cbind(metadata, as.data.frame(theta_panel, check.names = FALSE))
}

make_panel_pairs <- function(factor_set, theta_rows, probe_set) {
  grid <- expand.grid(
    local_pos = seq_along(factor_set$atlases),
    theta_row = as.integer(theta_rows),
    KEEP.OUT.ATTRS = FALSE
  )
  grid <- grid[order(grid$theta_row, grid$local_pos), , drop = FALSE]
  grid$probe_set <- probe_set
  grid
}

run_panel_probe <- function(factor_set,
                            panel,
                            pairs,
                            data_list,
                            loglik_fn,
                            M,
                            n_replicates,
                            local_control,
                            bootstrap_B,
                            n_cores,
                            seed,
                            verbose) {
  model <- factor_set$population_model
  theta <- as.matrix(panel[, model$hyper_names, drop = FALSE])
  run_shape_probe_pairs(
    factor_set = factor_set,
    pairs = pairs,
    theta = theta,
    theta_weights = panel$theta_weight,
    data_list = data_list,
    loglik_fn = loglik_fn,
    M = as.integer(M),
    n_replicates = as.integer(n_replicates),
    local_control = local_control,
    bootstrap_B = as.integer(bootstrap_B),
    n_cores = as.integer(n_cores),
    seed = as.integer(seed),
    verbose = verbose,
    geometry_control = list(
      max_directions = 4L,
      max_local_stencils = 48L,
      stencil_step = 0.65,
      local_stencil_step = 0.50,
      min_abs_standardized_residual = 1.0,
      min_abs_residual = 0.05,
      max_stencils_per_local = 3L
    )
  )
}

panel_delta_summary <- function(residuals, n_locals, label) {
  residuals <- as.data.frame(residuals, stringsAsFactors = FALSE, check.names = FALSE)
  if (!nrow(residuals)) return(data.frame())
  split_rows <- split(residuals, residuals$theta_row)
  rows <- lapply(split_rows, function(df) {
    ok <- is.finite(df$residual)
    theta_weight <- unique(df$theta_weight)
    theta_weight <- theta_weight[is.finite(theta_weight)]
    total <- if (any(ok)) sum(df$residual[ok]) else NA_real_
    total_se <- if (any(ok)) sqrt(sum(df$residual_se[ok]^2)) else NA_real_
    data.frame(
      stage = label,
      probe_set = as.character(df$probe_set[1L]),
      theta_row = as.integer(df$theta_row[1L]),
      theta_id = as.character(df$theta_id[1L]),
      theta_source = as.character(df$theta_source[1L]),
      theta_weight = if (length(theta_weight)) theta_weight[1L] else NA_real_,
      n_locals = nrow(df),
      n_ok = sum(ok),
      complete = sum(ok) == as.integer(n_locals),
      total_error = total,
      total_se = total_se,
      total_z = total / pmax(total_se, .Machine$double.eps),
      mean_local_error = if (any(ok)) mean(df$residual[ok]) else NA_real_,
      max_abs_local_error = if (any(ok)) max(abs(df$residual[ok])) else NA_real_,
      certified_biased_locals = sum(
        ok &
          as.character(df$atlas_status) == "certified" &
          abs(df$residual) >= 0.25 &
          abs(df$standardized_residual) >= 2
      ),
      uncertified_locals = sum(as.character(df$atlas_status) != "certified", na.rm = TRUE),
      check.names = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$centered_total_error_posterior_weight <- NA_real_
  out$centered_total_error_equal_weight <- NA_real_
  for (set in unique(out$probe_set)) {
    idx <- out$probe_set == set & out$complete & is.finite(out$total_error)
    if (!any(idx)) next
    w <- pmax(out$theta_weight[idx], 0)
    w <- if (sum(w) > 0) w / sum(w) else rep(1 / sum(idx), sum(idx))
    out$centered_total_error_posterior_weight[idx] <- out$total_error[idx] - sum(w * out$total_error[idx])
    out$centered_total_error_equal_weight[idx] <- out$total_error[idx] - mean(out$total_error[idx])
  }
  out[order(out$probe_set, out$theta_row), , drop = FALSE]
}

panel_metric_summary <- function(delta) {
  if (!nrow(delta)) return(data.frame())
  split_rows <- split(delta, paste(delta$stage, delta$probe_set, sep = "\r"))
  do.call(rbind, lapply(split_rows, function(df) {
    ok <- df$complete & is.finite(df$centered_total_error_equal_weight)
    w_ok <- df$complete & is.finite(df$centered_total_error_posterior_weight) & is.finite(df$theta_weight)
    w <- pmax(df$theta_weight[w_ok], 0)
    w <- if (length(w) && sum(w) > 0) w / sum(w) else numeric()
    data.frame(
      stage = as.character(df$stage[1L]),
      probe_set = as.character(df$probe_set[1L]),
      n_theta = nrow(df),
      complete_theta = sum(df$complete),
      complete_fraction = mean(df$complete),
      equal_centered_delta_rmse = if (any(ok)) sqrt(mean(df$centered_total_error_equal_weight[ok]^2)) else NA_real_,
      posterior_centered_delta_rmse = if (length(w)) {
        sqrt(sum(w * df$centered_total_error_posterior_weight[w_ok]^2))
      } else {
        NA_real_
      },
      max_abs_centered_delta = if (any(ok)) max(abs(df$centered_total_error_equal_weight[ok])) else NA_real_,
      max_abs_total_z = max(abs(df$total_z), na.rm = TRUE),
      total_certified_biased_locals = sum(df$certified_biased_locals, na.rm = TRUE),
      total_uncertified_locals = sum(df$uncertified_locals, na.rm = TRUE),
      check.names = FALSE
    )
  }))
}

plot_delta_panel <- function(delta, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  grDevices::png(file, width = 1800, height = 1100)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  }, add = TRUE)
  graphics::par(mfrow = c(1, 2), mar = c(9, 4, 4, 1))
  for (set in c("repair", "holdout")) {
    df <- delta[delta$probe_set == set, , drop = FALSE]
    if (!nrow(df)) {
      graphics::plot.new()
      graphics::title(set)
      next
    }
    wide <- reshape(
      df[, c("theta_row", "stage", "centered_total_error_equal_weight"), drop = FALSE],
      idvar = "theta_row",
      timevar = "stage",
      direction = "wide"
    )
    ycols <- grep("^centered_total_error_equal_weight\\.", names(wide), value = TRUE)
    y <- as.matrix(wide[, ycols, drop = FALSE])
    colnames(y) <- sub("^centered_total_error_equal_weight\\.", "", ycols)
    graphics::barplot(
      t(y),
      beside = TRUE,
      las = 2,
      col = c("grey55", "firebrick3")[seq_len(ncol(y))],
      main = paste0(set, " theta panel: centered total error"),
      ylab = "Delta_c(theta)",
      names.arg = paste0("theta ", wide$theta_row)
    )
    graphics::abline(h = 0, lty = 2)
    graphics::legend("topright", legend = colnames(y), fill = c("grey55", "firebrick3")[seq_len(ncol(y))], cex = 0.8)
  }
}

plot_residual_heatmap <- function(residuals, file, title) {
  residuals <- as.data.frame(residuals, stringsAsFactors = FALSE, check.names = FALSE)
  locals <- sort(unique(as.integer(residuals$local_pos)))
  theta_rows <- sort(unique(as.integer(residuals$theta_row)))
  mat <- matrix(NA_real_, nrow = length(locals), ncol = length(theta_rows),
                dimnames = list(locals, theta_rows))
  idx <- cbind(match(as.integer(residuals$local_pos), locals), match(as.integer(residuals$theta_row), theta_rows))
  mat[idx] <- residuals$residual
  lim <- max(abs(mat), na.rm = TRUE)
  if (!is.finite(lim) || lim <= 0) lim <- 1
  pal <- grDevices::colorRampPalette(c("navy", "white", "firebrick3"))(101)
  br <- seq(-lim, lim, length.out = length(pal) + 1L)
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  grDevices::png(file, width = 1500, height = 1300)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  }, add = TRUE)
  graphics::par(mar = c(6, 6, 4, 6))
  graphics::image(
    x = seq_along(theta_rows),
    y = seq_along(locals),
    z = t(mat),
    col = pal,
    breaks = br,
    xaxt = "n",
    yaxt = "n",
    xlab = "theta panel row",
    ylab = "local",
    main = title
  )
  graphics::axis(1, at = seq_along(theta_rows), labels = theta_rows, las = 2)
  graphics::axis(2, at = seq_along(locals), labels = locals, las = 2, cex.axis = 0.55)
  graphics::box()
}

combine_probe_residuals <- function(...) {
  out <- .local_atlas_rbind_fill(lapply(list(...), function(x) x$residuals))
  out[order(out$probe_set, out$theta_row, out$local_pos), , drop = FALSE]
}

args <- parse_cli_args(commandArgs(trailingOnly = TRUE))
detected_cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
if (!is.finite(detected_cores) || detected_cores < 1L) detected_cores <- 1L

label <- arg_chr(args, "label", "total_error_panel_pilot_v1")
checkpoint_file <- arg_chr(
  args,
  "checkpoint_file",
  file.path("benchmarks", "results", "population_emc_grouped_active_panel_v1_checkpoint_post_shape.rds")
)
data_file <- arg_chr(args, "data_file", file.path("benchmarks", "samples", "full_EMC2.RData"))
out_dir <- arg_chr(args, "out_dir", file.path("benchmarks", "results"))
seed <- arg_int(args, "seed", 20260529L)
cores <- arg_int(args, "cores", min(4L, detected_cores))
n_train_theta <- arg_int(args, "n_train_theta", 8L)
n_holdout_theta <- arg_int(args, "n_holdout_theta", 4L)
max_axes <- arg_int(args, "max_axes", 4L)
n_draws <- arg_int(args, "n_draws", 3000L)
probe_M <- arg_int(args, "probe_M", 256L)
probe_replicates <- arg_int(args, "probe_replicates", 2L)
bootstrap_B <- arg_int(args, "bootstrap_B", 200L)
repair_M <- arg_int(args, "repair_M", probe_M)
repair_replicates <- arg_int(args, "repair_replicates", probe_replicates)
repair_max_updates <- arg_int(args, "repair_max_updates", 24L)
repair_mode <- tolower(arg_chr(args, "repair_mode", "geometry"))
if (!repair_mode %in% c("geometry", "exact", "none")) {
  stop("repair_mode must be geometry, exact, or none.")
}
local_target_cess <- arg_num(args, "local_target_cess", 0.90)
local_mcmc_moves <- arg_int(args, "local_mcmc_moves", 2L)
local_max_steps <- arg_int(args, "local_max_steps", 128L)
smc_verbose <- arg_lgl(args, "smc_verbose", FALSE)
stop_after_pre_audit <- arg_lgl(args, "stop_after_pre_audit", FALSE)

prefix <- file.path(out_dir, paste0("population_emc_", label))
config_file <- paste0(prefix, "_config.rds")
theta_panel_file <- paste0(prefix, "_theta_panel.csv")
pre_train_file <- paste0(prefix, "_pre_train_probe.rds")
pre_holdout_file <- paste0(prefix, "_pre_holdout_probe.rds")
pre_residual_csv <- paste0(prefix, "_pre_residuals.csv")
repair_file <- paste0(prefix, "_repair.rds")
repair_csv <- paste0(prefix, "_repair_selected.csv")
post_residual_csv <- paste0(prefix, "_post_residuals.csv")
delta_csv <- paste0(prefix, "_delta_panel.csv")
metric_csv <- paste0(prefix, "_delta_metrics.csv")
holdout_file <- paste0(prefix, "_holdout_validation.rds")
holdout_csv <- paste0(prefix, "_holdout_validation_summary.csv")
gate_file <- paste0(prefix, "_reweight_gate.rds")
gate_csv <- paste0(prefix, "_reweight_gate_summary.csv")
posterior_csv <- paste0(prefix, "_posterior_comparison.csv")
results_file <- paste0(prefix, "_results.rds")
delta_plot_file <- paste0(prefix, "_delta_panel.png")
pre_heatmap_file <- paste0(prefix, "_pre_residual_heatmap.png")
post_heatmap_file <- paste0(prefix, "_post_residual_heatmap.png")
posterior_plot_file <- paste0(prefix, "_posteriors.png")

if (!file.exists(checkpoint_file)) stop("Missing checkpoint: ", checkpoint_file)
if (!file.exists(data_file)) stop("Missing EMC2 data: ", data_file)

announce("Load checkpoint and EMC reference", checkpoint_file)
checkpoint <- readRDS(checkpoint_file)
factor_set <- validate_local_atlas_factor_set(checkpoint$factor_set)
fit <- checkpoint$fit

load(data_file)
if (!exists("ELP_DDM", inherits = FALSE)) stop("The EMC2 data file must define ELP_DDM.")
emc <- ELP_DDM[[1L]]
data_list <- emc$data
model_factory <- emc$model
alpha_names <- emc$par_names

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

local_control <- list(
  target_cess = local_target_cess,
  n_mcmc_moves = local_mcmc_moves,
  max_steps = local_max_steps
)

settings <- list(
  label = label,
  checkpoint_file = checkpoint_file,
  data_file = data_file,
  seed = seed,
  cores = cores,
  n_train_theta = n_train_theta,
  n_holdout_theta = n_holdout_theta,
  max_axes = max_axes,
  probe_M = probe_M,
  probe_replicates = probe_replicates,
  repair_M = repair_M,
  repair_replicates = repair_replicates,
  repair_max_updates = repair_max_updates,
  repair_mode = repair_mode,
  local_control = local_control
)
save_stage(settings, config_file)

announce("Select complete theta panels")
panel <- select_panel_rows(
  fit = fit,
  factor_set = factor_set,
  reference_draws = emc_draws,
  n_train = n_train_theta,
  n_holdout = n_holdout_theta,
  max_axes = max_axes,
  n_draws = n_draws,
  seed = seed + 1L
)
write_csv(panel, theta_panel_file)
cat(sprintf("Panel theta: %d repair, %d holdout; locals per theta: %d\n",
            sum(panel$probe_set == "repair"), sum(panel$probe_set == "holdout"), factor_set$n_locals))

train_pairs <- make_panel_pairs(
  factor_set = factor_set,
  theta_rows = panel$theta_row[panel$probe_set == "repair"],
  probe_set = "repair"
)
holdout_pairs <- make_panel_pairs(
  factor_set = factor_set,
  theta_rows = panel$theta_row[panel$probe_set == "holdout"],
  probe_set = "holdout"
)

announce("Run complete repair-theta local evidence panel",
         sprintf("%d local/theta pairs, M=%d, replicates=%d",
                 nrow(train_pairs), probe_M, probe_replicates))
pre_train <- run_panel_probe(
  factor_set = factor_set,
  panel = panel,
  pairs = train_pairs,
  data_list = data_list,
  loglik_fn = loglik_emc2,
  M = probe_M,
  n_replicates = probe_replicates,
  local_control = local_control,
  bootstrap_B = bootstrap_B,
  n_cores = cores,
  seed = seed + 10L,
  verbose = smc_verbose
)
save_stage(pre_train, pre_train_file)

announce("Run complete held-out theta local evidence panel",
         sprintf("%d local/theta pairs, M=%d, replicates=%d",
                 nrow(holdout_pairs), probe_M, probe_replicates))
pre_holdout <- run_panel_probe(
  factor_set = factor_set,
  panel = panel,
  pairs = holdout_pairs,
  data_list = data_list,
  loglik_fn = loglik_emc2,
  M = probe_M,
  n_replicates = probe_replicates,
  local_control = local_control,
  bootstrap_B = bootstrap_B,
  n_cores = cores,
  seed = seed + 20L,
  verbose = smc_verbose
)
save_stage(pre_holdout, pre_holdout_file)

pre_residuals <- combine_probe_residuals(pre_train, pre_holdout)
write_csv(pre_residuals, pre_residual_csv)
pre_delta <- panel_delta_summary(pre_residuals, n_locals = factor_set$n_locals, label = "pre")
write_csv(pre_delta, delta_csv)
plot_residual_heatmap(pre_residuals, pre_heatmap_file, "Pre-repair direct local evidence residuals")

if (isTRUE(stop_after_pre_audit) || identical(repair_mode, "none")) {
  metrics <- panel_metric_summary(pre_delta)
  write_csv(metrics, metric_csv)
  result <- list(
    settings = settings,
    panel = panel,
    pre_train = pre_train,
    pre_holdout = pre_holdout,
    pre_residuals = pre_residuals,
    delta = pre_delta,
    metrics = metrics,
    factor_set = factor_set,
    fit = fit,
    emc_draws = emc_draws
  )
  save_stage(result, results_file)
  quit(save = "no", status = 0L)
}

announce("Repair atlas from complete repair panel", repair_mode)
repair_control <- list(
  max_repairs = repair_max_updates,
  max_updates = repair_max_updates,
  M = repair_M,
  n_mcmc_moves = local_mcmc_moves,
  target_cess = local_target_cess,
  max_steps = local_max_steps,
  direct_confirmation_reps = repair_replicates,
  direct_confirmation_M = repair_M,
  audit_pre_repair = TRUE,
  audit_post_repair = TRUE,
  audit_scope = "repaired",
  stop_on_empty = FALSE,
  min_abs_standardized_residual = 2.0,
  min_abs_residual = 0.25,
  strict_probe_confirmation = TRUE,
  allow_exact_residual_repair = identical(repair_mode, "exact"),
  min_confirmation_abs_residual = 0.25,
  min_confirmation_abs_standardized_residual = 2.0,
  min_confirmed_shape_points = 2L,
  min_confirmed_shape_theta = 2L
)
repair_edge_control <- list(
  edge_neighbors = 2L,
  max_intermediates = 4L,
  min_overlap_ess = 0.03,
  max_se = 1.25,
  max_forward_reverse_gap = 1.25,
  max_taylor_gap = 3.0,
  require_bar_converged = TRUE
)
repair <- if (identical(repair_mode, "exact")) {
  repair_shape_selected_probe_pairs(
    factor_set = factor_set,
    shape_probe = pre_train,
    data_list = data_list,
    loglik_fn = loglik_emc2,
    control = repair_control,
    local_control = local_control,
    edge_control = repair_edge_control,
    n_cores = cores,
    seed = seed + 30L,
    verbose = smc_verbose
  )
} else {
  repair_shape_residual_geometry(
    factor_set = factor_set,
    shape_probe = pre_train,
    data_list = data_list,
    loglik_fn = loglik_emc2,
    control = repair_control,
    local_control = local_control,
    edge_control = repair_edge_control,
    n_cores = cores,
    seed = seed + 30L,
    verbose = smc_verbose
  )
}
save_stage(repair, repair_file)
write_csv(repair$selected_candidates, repair_csv)

announce("Re-audit complete panels against repaired atlas")
post_train_residuals <- .local_shape_residuals_with_new_raw(
  reference_residuals = pre_train$residuals,
  factor_set = repair$factor_set,
  cloud = pre_train$cloud,
  n_cores = cores
)
post_holdout_residuals <- .local_shape_residuals_with_new_raw(
  reference_residuals = pre_holdout$residuals,
  factor_set = repair$factor_set,
  cloud = pre_holdout$cloud,
  n_cores = cores
)
post_residuals <- .local_atlas_rbind_fill(list(post_train_residuals, post_holdout_residuals))
post_residuals <- post_residuals[order(post_residuals$probe_set, post_residuals$theta_row, post_residuals$local_pos), , drop = FALSE]
write_csv(post_residuals, post_residual_csv)
post_delta <- panel_delta_summary(post_residuals, n_locals = factor_set$n_locals, label = "post")
delta <- rbind(pre_delta, post_delta)
write_csv(delta, delta_csv)
metrics <- panel_metric_summary(delta)
write_csv(metrics, metric_csv)
plot_delta_panel(delta, delta_plot_file)
plot_residual_heatmap(post_residuals, post_heatmap_file, "Post-repair direct local evidence residuals")

announce("Holdout validation and outer reweight gate")
holdout_validation <- validate_shape_repair_holdout(
  shape_repair = repair,
  shape_probe = pre_holdout,
  control = list(
    probe_set = "holdout",
    rerun_direct_probes = FALSE,
    require_holdout_improvement = FALSE,
    max_pair_centered_rmse_ratio = 1.10,
    max_local_centered_rmse_ratio = 1.10,
    max_total_centered_rmse_ratio = 1.10,
    max_abs_total_increase = 0.25,
    max_graph_edge_z_increase = 1.0
  ),
  n_cores = cores,
  seed = seed + 40L,
  verbose = smc_verbose
)
save_stage(holdout_validation, holdout_file)
write_csv(holdout_validation$summary, holdout_csv)

gate <- shape_repair_outer_reweight_gate(
  shape_repair = repair,
  fit = fit,
  old_factor_set = factor_set,
  holdout_validation = holdout_validation,
  reference_draws = emc_draws,
  control = list(
    require_holdout_acceptance = FALSE,
    min_reweight_ess_fraction = 0.50,
    low_ess_rerun_fraction = 0.25,
    max_psis_k = 0.70,
    n_draws = min(n_draws, nrow(emc_draws))
  ),
  n_cores = cores,
  seed = seed + 50L
)
save_stage(gate, gate_file)
write_csv(gate$summary, gate_csv)

old_draws <- local_atlas_draws_from_fit(
  fit,
  population_model = factor_set$population_model,
  n_draws = min(n_draws, nrow(emc_draws)),
  seed = seed + 60L
)
new_draws <- local_atlas_draws_from_fit(
  gate$reweighted_fit %||% fit,
  population_model = factor_set$population_model,
  n_draws = min(n_draws, nrow(emc_draws)),
  seed = seed + 61L
)
old_comparison <- local_atlas_compare_posterior_draws(emc_draws, old_draws)
new_comparison <- local_atlas_compare_posterior_draws(emc_draws, new_draws)
old_comparison$stage <- "pre"
new_comparison$stage <- "post_reweight"
posterior_comparison <- rbind(old_comparison, new_comparison)
write_csv(posterior_comparison, posterior_csv)

dir.create(dirname(posterior_plot_file), recursive = TRUE, showWarnings = FALSE)
grDevices::png(posterior_plot_file, width = 1800, height = 1400)
plot_posteriors(
  emc_draws,
  new_draws,
  labels = c("EMC2", "panel-reweighted"),
  cols = c("black", "firebrick3"),
  n_cols = 4L
)
grDevices::dev.off()
cat(sprintf("Saved: %s\n", posterior_plot_file))

result <- list(
  settings = settings,
  panel = panel,
  pre_train = pre_train,
  pre_holdout = pre_holdout,
  pre_residuals = pre_residuals,
  repair = repair,
  post_residuals = post_residuals,
  delta = delta,
  metrics = metrics,
  holdout_validation = holdout_validation,
  gate = gate,
  posterior_comparison = posterior_comparison,
  factor_set = repair$factor_set,
  fit = gate$reweighted_fit %||% fit,
  old_fit = fit,
  emc_draws = emc_draws,
  plot_files = c(delta_panel = delta_plot_file, pre_heatmap = pre_heatmap_file,
                 post_heatmap = post_heatmap_file, posterior = posterior_plot_file)
)
save_stage(result, results_file)

announce("Done")
print(metrics, row.names = FALSE)
print(gate$summary, row.names = FALSE)
