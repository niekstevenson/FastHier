#!/usr/bin/env Rscript

rm(list = ls())

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[1L]))
} else {
  normalizePath("benchmarks/plot_query_native_visuals.R")
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

arg_chr <- function(args, key, default) {
  val <- args[[key]]
  if (is.null(val) || !nzchar(val)) default else as.character(val)
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

finite_limit <- function(x, default = 1) {
  x <- abs(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x)) max(x) else default
}

grid_matrix <- function(df, value_col, x_col = "mu_alpha", y_col = "log_sigma2_alpha") {
  x <- sort(unique(df[[x_col]]))
  y <- sort(unique(df[[y_col]]))
  z <- matrix(NA_real_, nrow = length(x), ncol = length(y), dimnames = list(signif(x, 5), signif(y, 5)))
  for (i in seq_len(nrow(df))) {
    ix <- match(df[[x_col]][i], x)
    iy <- match(df[[y_col]][i], y)
    z[ix, iy] <- df[[value_col]][i]
  }
  list(x = x, y = y, z = z)
}

plot_grid_image <- function(grid,
                            main,
                            xlab = "mu_alpha",
                            ylab = "log_sigma2_alpha",
                            zlim = NULL,
                            palette = grDevices::hcl.colors(101, "RdBu", rev = TRUE),
                            draw_contour = TRUE) {
  if (is.null(zlim)) {
    lim <- finite_limit(grid$z)
    zlim <- c(-lim, lim)
  }
  graphics::image(grid$x, grid$y, grid$z, col = palette, zlim = zlim, xlab = xlab, ylab = ylab, main = main)
  if (isTRUE(draw_contour) && any(is.finite(grid$z))) {
    graphics::contour(grid$x, grid$y, grid$z, add = TRUE, drawlabels = FALSE, nlevels = 8, col = "grey25")
  }
}

normalize_weights <- function(w) {
  w <- as.numeric(w)
  w[!is.finite(w) | w < 0] <- 0
  if (sum(w) <= 0) rep(NA_real_, length(w)) else w / sum(w)
}

method_col <- c(
  analytic = "black",
  atlas = "steelblue4",
  nested_M192 = "darkorange3",
  nested_M64 = "purple4"
)

cli_args <- parse_cli_args(commandArgs(trailingOnly = TRUE))
toy_result <- arg_chr(cli_args, "toy_result", file.path("benchmarks", "results", "query_native_visual_results.rds"))
emc_result <- arg_chr(cli_args, "emc_result", file.path("benchmarks", "results", "emc_single_local_hard14_query_native_roleinfl_24x512_results.rds"))
emc_repair <- arg_chr(cli_args, "emc_repair", file.path("benchmarks", "results", "emc_single_local_hard14_query_native_roleinfl_24x512_repaired_probe.rds"))
label <- arg_chr(cli_args, "label", "query_native_visual")
out_dir <- arg_chr(cli_args, "out_dir", file.path("benchmarks", "results"))

heat_file <- file.path(out_dir, paste0(label, "_heat_plot.png"))
replicate_file <- file.path(out_dir, paste0(label, "_replicate_logm_errors.png"))
posterior_file <- file.path(out_dir, paste0(label, "_posterior_overlap.png"))
summary_file <- file.path(out_dir, paste0(label, "_visual_summary.csv"))

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

toy <- readRDS(toy_result)
emc <- readRDS(emc_result)
repair <- if (file.exists(emc_repair)) readRDS(emc_repair) else NULL

total <- toy$total_log_m
posterior <- toy$posterior_grid
theta_design <- toy$theta_design

heat_methods <- intersect(c("atlas", "nested_M192", "nested_M64"), unique(total$method))
grDevices::png(heat_file, width = 1800, height = 1400)
graphics::layout(matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE))
graphics::par(mar = c(4.4, 4.4, 3.2, 1.2))

analytic_post <- posterior[posterior$method == "analytic", , drop = FALSE]
plot_grid_image(
  grid_matrix(analytic_post, "weight"),
  main = "Analytic theta posterior mass",
  zlim = c(0, max(analytic_post$weight, na.rm = TRUE)),
  palette = grDevices::hcl.colors(101, "YlOrRd", rev = FALSE),
  draw_contour = TRUE
)
graphics::points(theta_design[, "mu_alpha"], theta_design[, "log_sigma2_alpha"], pch = 4, cex = 0.8, lwd = 1.4)
graphics::legend("topright", legend = "atlas chart", pch = 4, bty = "n", cex = 0.9)

for (method in heat_methods) {
  df <- total[total$method == method, , drop = FALSE]
  method_limit <- finite_limit(df$centered_total_error)
  plot_grid_image(
    grid_matrix(df, "centered_total_error"),
    main = sprintf("%s centered total log m error (scale +/- %.2f)", method, method_limit),
    zlim = c(-method_limit, method_limit)
  )
  missing <- df[!is.finite(df$total_log_m), , drop = FALSE]
  if (nrow(missing)) {
    graphics::points(missing$mu_alpha, missing$log_sigma2_alpha, pch = 4, col = "firebrick3", lwd = 2)
  }
  graphics::points(theta_design[, "mu_alpha"], theta_design[, "log_sigma2_alpha"], pch = 3, cex = 0.65, col = "grey20")
}
grDevices::dev.off()

truth_local <- toy$reference_rows[, c("local_pos", "theta_row", "estimate"), drop = FALSE]
names(truth_local)[3L] <- "truth"
nested_raw <- merge(toy$nested_rows_raw, truth_local, by = c("local_pos", "theta_row"), all.x = TRUE)
nested_raw$error <- nested_raw$estimate - nested_raw$truth

truth_total <- total[total$method == "analytic", c("theta_row", "total_log_m"), drop = FALSE]
names(truth_total)[2L] <- "truth_total"
total_raw_parts <- split(nested_raw, interaction(nested_raw$method, nested_raw$replicate, nested_raw$theta_row, drop = TRUE))
total_reps <- do.call(rbind, lapply(total_raw_parts, function(df) {
  data.frame(
    method = df$method[1L],
    replicate = df$replicate[1L],
    theta_row = df$theta_row[1L],
    estimate = if (all(is.finite(df$estimate))) sum(df$estimate) else NA_real_,
    n_locals = nrow(df),
    check.names = FALSE
  )
}))
total_reps <- merge(total_reps, truth_total, by = "theta_row", all.x = TRUE)
total_reps$error <- total_reps$estimate - total_reps$truth_total
atlas_total <- total[total$method == "atlas", c("theta_row", "total_log_m", "analytic_total_log_m", "complete"), drop = FALSE]
atlas_total$error <- atlas_total$total_log_m - atlas_total$analytic_total_log_m

emc_rows <- emc$rows
ref_raw <- merge(
  emc$reference_raw,
  emc_rows[, c("theta_row", "reference_log_marginal", "reference_rep_sd"), drop = FALSE],
  by = "theta_row",
  all.x = TRUE
)
ref_raw$error <- ref_raw$reference_log_marginal.x - ref_raw$reference_log_marginal.y
emc_rows$pre_error <- emc_rows$atlas_log_marginal - emc_rows$reference_log_marginal
repair_rows <- if (!is.null(repair)) repair$rows else data.frame()

grDevices::png(replicate_file, width = 1800, height = 1200)
graphics::layout(matrix(c(1, 2), nrow = 2, byrow = TRUE))
graphics::par(mar = c(4.5, 4.7, 3.4, 1.2))

toy_ylim <- range(c(total_reps$error, atlas_total$error), na.rm = TRUE)
graphics::plot(
  NA,
  xlim = range(total_reps$theta_row, na.rm = TRUE),
  ylim = toy_ylim,
  xlab = "theta grid row",
  ylab = "total log m error vs analytic",
  main = "Toy: replicated nested-SMC total marginal log-likelihood error"
)
graphics::abline(h = 0, lty = 2, col = "grey45")
for (method in intersect(c("nested_M64", "nested_M192"), unique(total_reps$method))) {
  df <- total_reps[total_reps$method == method, , drop = FALSE]
  jitter <- if (identical(method, "nested_M64")) -0.10 else 0.10
  graphics::points(df$theta_row + jitter, df$error, pch = 16, cex = 0.65, col = adjustcolor(method_col[[method]], 0.55))
}
graphics::points(
  atlas_total$theta_row,
  atlas_total$error,
  pch = ifelse(atlas_total$complete, 18, 4),
  cex = 1.25,
  col = method_col[["atlas"]],
  lwd = 2
)
graphics::legend(
  "topleft",
  legend = c("nested_M64 replicates", "nested_M192 replicates", "atlas total estimate", "zero error"),
  pch = c(16, 16, 18, NA),
  lty = c(NA, NA, NA, 2),
  col = c(method_col[["nested_M64"]], method_col[["nested_M192"]], method_col[["atlas"]], "grey45"),
  bty = "n",
  cex = 0.9
)

emc_ylim <- range(c(ref_raw$error, emc_rows$pre_error, repair_rows$error), na.rm = TRUE)
graphics::plot(
  NA,
  xlim = range(emc_rows$theta_row),
  ylim = emc_ylim,
  xlab = "EMC audit theta row",
  ylab = "local log m error vs nested-SMC mean",
  main = "EMC hard local: reference replicates and atlas/repair errors"
)
graphics::abline(h = 0, lty = 2, col = "grey45")
rep_center <- mean(seq_len(max(ref_raw$replicate, na.rm = TRUE)))
graphics::points(
  ref_raw$theta_row + 0.06 * (ref_raw$replicate - rep_center),
  ref_raw$error,
  pch = 16,
  cex = 0.75,
  col = adjustcolor("grey25", 0.55)
)
for (i in seq_len(nrow(emc_rows))) {
  if (is.finite(emc_rows$reference_rep_sd[i])) {
    graphics::segments(
      emc_rows$theta_row[i] - 0.22,
      -emc_rows$reference_rep_sd[i],
      emc_rows$theta_row[i] + 0.22,
      -emc_rows$reference_rep_sd[i],
      col = "grey55"
    )
    graphics::segments(
      emc_rows$theta_row[i] - 0.22,
      emc_rows$reference_rep_sd[i],
      emc_rows$theta_row[i] + 0.22,
      emc_rows$reference_rep_sd[i],
      col = "grey55"
    )
  }
}
graphics::points(
  emc_rows$theta_row - 0.12,
  emc_rows$pre_error,
  pch = ifelse(emc_rows$atlas_certified, 17, 4),
  cex = 1.25,
  col = ifelse(emc_rows$atlas_certified, "darkorange3", "firebrick3"),
  lwd = 2
)
if (nrow(repair_rows)) {
  graphics::points(
    repair_rows$theta_row + 0.12,
    repair_rows$error,
    pch = ifelse(repair_rows$status == "certified", 15, 4),
    cex = 1.15,
    col = ifelse(repair_rows$status == "certified", "steelblue4", "firebrick3"),
    lwd = 2
  )
}
graphics::legend(
  "bottomleft",
  legend = c("nested-SMC replicate - mean", "+/- one replicate SD", "atlas before repair", "repair probe"),
  pch = c(16, NA, 17, 15),
  lty = c(NA, 1, NA, NA),
  col = c("grey25", "grey55", "darkorange3", "steelblue4"),
  bty = "n",
  cex = 0.9
)
grDevices::dev.off()

posterior_methods <- intersect(c("analytic", "atlas", "nested_M192", "nested_M64"), unique(posterior$method))
analytic <- posterior[posterior$method == "analytic", , drop = FALSE]
analytic$weight <- normalize_weights(analytic$weight)
overlap_rows <- do.call(rbind, lapply(setdiff(posterior_methods, "analytic"), function(method) {
  df <- posterior[posterior$method == method, , drop = FALSE]
  df <- df[match(analytic$theta_row, df$theta_row), , drop = FALSE]
  w <- normalize_weights(df$weight)
  data.frame(method = method, overlap = sum(pmin(analytic$weight, w), na.rm = TRUE), check.names = FALSE)
}))

marginal_curve <- function(df, x_col) {
  out <- aggregate(df$weight, by = list(x = df[[x_col]]), sum, na.rm = TRUE)
  names(out)[2L] <- "weight"
  out$weight <- normalize_weights(out$weight)
  out
}

grDevices::png(posterior_file, width = 1800, height = 1350)
graphics::layout(matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE))
graphics::par(mar = c(4.5, 4.7, 3.3, 1.2))

plot_grid_image(
  grid_matrix(analytic, "weight"),
  main = "2D posterior overlap: analytic heat, method contours",
  zlim = c(0, max(analytic$weight, na.rm = TRUE)),
  palette = grDevices::hcl.colors(101, "Grays", rev = FALSE),
  draw_contour = FALSE
)
for (method in posterior_methods) {
  df <- posterior[posterior$method == method, , drop = FALSE]
  grid <- grid_matrix(df, "weight")
  graphics::contour(
    grid$x,
    grid$y,
    grid$z,
    add = TRUE,
    drawlabels = FALSE,
    nlevels = 5,
    col = method_col[[method]] %||% "grey30",
    lwd = ifelse(method == "analytic", 2, 1.5)
  )
}
graphics::legend("topright", legend = posterior_methods, col = method_col[posterior_methods], lwd = c(2, rep(1.5, length(posterior_methods) - 1L)), bty = "n", cex = 0.9)

plot_marginal_overlay <- function(x_col, main, xlab) {
  curves <- lapply(posterior_methods, function(method) {
    marginal_curve(posterior[posterior$method == method, , drop = FALSE], x_col)
  })
  names(curves) <- posterior_methods
  ylim <- range(unlist(lapply(curves, `[[`, "weight")), na.rm = TRUE)
  xlim <- range(unlist(lapply(curves, `[[`, "x")), na.rm = TRUE)
  graphics::plot(NA, xlim = xlim, ylim = ylim, xlab = xlab, ylab = "marginal posterior mass", main = main)
  for (method in posterior_methods) {
    curve <- curves[[method]]
    graphics::lines(curve$x, curve$weight, type = "b", pch = 16, col = method_col[[method]] %||% "grey30", lwd = ifelse(method == "analytic", 2.2, 1.6))
  }
  graphics::legend("topright", legend = posterior_methods, col = method_col[posterior_methods], lwd = 2, pch = 16, bty = "n", cex = 0.85)
}

plot_marginal_overlay("mu_alpha", "mu_alpha marginal overlap", "mu_alpha")
plot_marginal_overlay("log_sigma2_alpha", "log_sigma2_alpha marginal overlap", "log_sigma2_alpha")

graphics::barplot(
  overlap_rows$overlap,
  names.arg = overlap_rows$method,
  ylim = c(0, 1),
  las = 2,
  col = method_col[overlap_rows$method],
  ylab = "discrete overlap with analytic posterior",
  main = "Posterior overlap coefficient"
)
graphics::abline(h = c(0.9, 0.95), lty = 2, col = "grey55")
grDevices::dev.off()

visual_summary <- rbind(
  data.frame(
    diagnostic = "toy_atlas_total",
    value = sprintf(
      "finite_theta=%.3f, centered_total_rmse=%.4f",
      toy$total_summary$finite_theta_fraction[toy$total_summary$method == "atlas"],
      toy$total_summary$centered_total_rmse[toy$total_summary$method == "atlas"]
    ),
    check.names = FALSE
  ),
  data.frame(
    diagnostic = paste0("toy_overlap_", overlap_rows$method),
    value = sprintf("%.4f", overlap_rows$overlap),
    check.names = FALSE
  ),
  data.frame(
    diagnostic = "emc_before_repair",
    value = sprintf(
      "certified=%.3f, centered_rmse=%.4f",
      emc$summary$certified_fraction,
      emc$summary$atlas_centered_rmse_vs_reference
    ),
    check.names = FALSE
  )
)
if (nrow(repair_rows)) {
  finite_repair <- is.finite(repair_rows$error)
  centered <- repair_rows$error - mean(repair_rows$error[finite_repair], na.rm = TRUE)
  visual_summary <- rbind(
    visual_summary,
    data.frame(
      diagnostic = "emc_after_repair_probe",
      value = sprintf(
        "certified=%.3f, rmse=%.4f, centered_rmse=%.4f",
        mean(repair_rows$status == "certified"),
        sqrt(mean(repair_rows$error[finite_repair]^2)),
        sqrt(mean(centered[finite_repair]^2))
      ),
      check.names = FALSE
    )
  )
}
utils::write.csv(visual_summary, summary_file, row.names = FALSE)

cat(sprintf("Saved heat plot: %s\n", heat_file))
cat(sprintf("Saved replicate log-m error plot: %s\n", replicate_file))
cat(sprintf("Saved posterior overlap plot: %s\n", posterior_file))
cat(sprintf("Saved visual summary: %s\n", summary_file))
print(visual_summary, row.names = FALSE)
