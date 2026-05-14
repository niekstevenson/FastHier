rm(list = ls())
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[1L]))
} else {
  normalizePath("benchmarks/run_shifted_gamma_hierarchy_meta_compare.R")
}
repo_dir <- dirname(dirname(script_path))
setwd(repo_dir)

parse_cli_args <- function(args) {
  out <- list()
  if (!length(args)) return(out)
  for (arg in args) {
    if (!startsWith(arg, "--")) next
    arg <- sub("^--", "", arg)
    parts <- strsplit(arg, "=", fixed = TRUE)[[1L]]
    key <- gsub("-", "_", parts[1L])
    value <- if (length(parts) > 1L) paste(parts[-1L], collapse = "=") else "true"
    out[[key]] <- value
  }
  out
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

arg_chr <- function(args, key, default = NULL) {
  val <- args[[key]]
  if (is.null(val) || !nzchar(val)) default else as.character(val)
}

arg_lgl <- function(args, key, default = FALSE) {
  val <- args[[key]]
  if (is.null(val) || !nzchar(val)) return(isTRUE(default))
  tolower(as.character(val)) %in% c("1", "true", "t", "yes", "y")
}

build_config_grid <- function(config_set = "core") {
  config_set <- match.arg(config_set, c("core"))
  data.frame(
    method = c("static", "refresh_dmis"),
    label = c("static", "refresh_dmis"),
    stringsAsFactors = FALSE
  )
}

draw_quantiles <- function(x, probs = c(0.1, 0.5, 0.9)) {
  stats::quantile(as.numeric(x), probs = probs, na.rm = TRUE, names = FALSE)
}

summarize_result <- function(result_file, config, plot_file, log_file) {
  res <- readRDS(result_file)
  workflow_draws <- res$workflow_draws
  stan_draws <- res$stan_draws
  common <- intersect(names(workflow_draws), names(stan_draws))

  row <- data.frame(
    label = config$label,
    method = config$method,
    results_file = normalizePath(result_file, winslash = "/", mustWork = FALSE),
    plot_file = normalizePath(plot_file, winslash = "/", mustWork = FALSE),
    log_file = normalizePath(log_file, winslash = "/", mustWork = FALSE),
    log_evidence = as.numeric(res$fit$log_evidence %||% NA_real_),
    mcse_log_evidence = as.numeric(res$fit$mcse_log_evidence %||% NA_real_),
    outer_rounds = as.integer(res$fit$meta$rounds %||% NA_integer_),
    stringsAsFactors = FALSE
  )

  for (nm in common) {
    wf <- workflow_draws[[nm]]
    st <- stan_draws[[nm]]
    qwf <- draw_quantiles(wf)
    qst <- draw_quantiles(st)
    row[[paste0(nm, "_mean_diff")]] <- mean(wf) - mean(st)
    row[[paste0(nm, "_sd_diff")]] <- stats::sd(wf) - stats::sd(st)
    row[[paste0(nm, "_q10_diff")]] <- qwf[1L] - qst[1L]
    row[[paste0(nm, "_q50_diff")]] <- qwf[2L] - qst[2L]
    row[[paste0(nm, "_q90_diff")]] <- qwf[3L] - qst[3L]
  }
  row
}

plot_evidence_summary <- function(summary_df, file) {
  ok <- is.finite(summary_df$log_evidence)
  if (!any(ok)) return(invisible(NULL))
  x <- summary_df[ok, , drop = FALSE]
  ord <- order(x$log_evidence)
  x <- x[ord, , drop = FALSE]
  y <- seq_len(nrow(x))
  xmin <- min(x$log_evidence - ifelse(is.finite(x$mcse_log_evidence), x$mcse_log_evidence, 0))
  xmax <- max(x$log_evidence + ifelse(is.finite(x$mcse_log_evidence), x$mcse_log_evidence, 0))
  pad <- 0.05 * max(1, xmax - xmin)

  grDevices::png(file, width = 1400, height = max(900, 55 * nrow(x)))
  graphics::par(mar = c(5, 18, 3, 2))
  graphics::plot(
    x$log_evidence, y,
    xlim = c(xmin - pad, xmax + pad),
    ylim = c(0.5, nrow(x) + 0.5),
    yaxt = "n",
    pch = 19,
    col = "firebrick3",
    xlab = "Log evidence",
    ylab = "",
    main = "Shifted-Gamma Evidence By Configuration"
  )
  graphics::axis(2, at = y, labels = x$label, las = 1, cex.axis = 0.9)
  has_mcse <- is.finite(x$mcse_log_evidence)
  if (any(has_mcse)) {
    graphics::segments(
      x0 = x$log_evidence[has_mcse] - x$mcse_log_evidence[has_mcse],
      y0 = y[has_mcse],
      x1 = x$log_evidence[has_mcse] + x$mcse_log_evidence[has_mcse],
      y1 = y[has_mcse],
      lwd = 2,
      col = "gray35"
    )
  }
  graphics::abline(v = max(x$log_evidence), lty = 3, col = "gray60")
  grDevices::dev.off()
}

plot_mu_shift_overlay <- function(result_files, labels, file) {
  results <- lapply(result_files, readRDS)
  stan <- results[[1L]]$stan_draws$mu_shift
  workflows <- lapply(results, function(res) res$workflow_draws$mu_shift)

  dens_all <- c(list(stats::density(stan)), lapply(workflows, stats::density))
  xlim <- range(unlist(lapply(dens_all, `[[`, "x")))
  ylim <- c(0, 1.05 * max(unlist(lapply(dens_all, `[[`, "y"))))
  cols <- c("black", grDevices::hcl.colors(length(workflows), palette = "Dark 3"))

  grDevices::png(file, width = 1400, height = 900)
  graphics::par(mar = c(4, 4, 3, 10), xpd = NA)
  graphics::plot(
    dens_all[[1L]],
    xlim = xlim,
    ylim = ylim,
    lwd = 3,
    col = cols[1L],
    main = "mu_shift Posterior Comparison",
    xlab = "mu_shift",
    ylab = "Density"
  )
  for (i in seq_along(workflows)) {
    graphics::lines(dens_all[[i + 1L]], col = cols[i + 1L], lwd = 2)
  }
  graphics::legend(
    "topright",
    inset = c(-0.28, 0),
    legend = c("Stan", labels),
    col = cols,
    lwd = c(3, rep(2, length(workflows))),
    bty = "n",
    cex = 0.85
  )
  grDevices::dev.off()
}

cli_args <- parse_cli_args(commandArgs(trailingOnly = TRUE))
config_set <- arg_chr(cli_args, "config_set", "core")
results_dir <- arg_chr(cli_args, "results_dir", file.path("benchmarks", "results", "meta_shifted_gamma"))
labels_filter <- arg_chr(cli_args, "labels", NULL)
inner_verbose <- arg_lgl(cli_args, "inner_verbose", FALSE)
fail_fast <- arg_lgl(cli_args, "fail_fast", FALSE)

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
logs_dir <- file.path(results_dir, "logs")
dir.create(logs_dir, showWarnings = FALSE, recursive = TRUE)

configs <- build_config_grid(config_set)
if (!is.null(labels_filter)) {
  keep <- trimws(strsplit(labels_filter, ",", fixed = TRUE)[[1L]])
  configs <- configs[configs$label %in% keep, , drop = FALSE]
}
if (!nrow(configs)) {
  stop("No configurations selected.")
}

forward_keys <- c(
  "mc_cores",
  "pilot_size",
  "pilot_particles",
  "full_particles",
  "outer_particles",
  "outer_mcmc_moves",
  "outer_max_rounds",
  "base_seed"
)
forward_args <- unlist(lapply(
  forward_keys,
  function(key) {
    val <- cli_args[[key]]
    if (is.null(val) || !nzchar(val)) {
      NULL
    } else {
      sprintf("--%s=%s", gsub("_", "-", key), val)
    }
  }
), use.names = FALSE)

summary_rows <- list()

for (i in seq_len(nrow(configs))) {
  cfg <- configs[i, , drop = FALSE]
  label <- cfg$label
  result_file <- file.path(results_dir, sprintf("shifted_gamma_%s_results.rds", label))
  plot_file <- file.path(results_dir, sprintf("shifted_gamma_%s_posteriors.png", label))
  log_file <- file.path(logs_dir, sprintf("shifted_gamma_%s.log", label))

  cat(sprintf("[%d/%d] %s\n", i, nrow(configs), label))

  args <- c(
    "benchmarks/run_shifted_gamma_hierarchy_compare.R",
    sprintf("--label=%s", label),
    sprintf("--method=%s", cfg$method),
    sprintf("--results-file=%s", result_file),
    sprintf("--plot-file=%s", plot_file),
    sprintf("--verbose=%s", tolower(as.character(inner_verbose))),
    forward_args
  )

  status <- tryCatch(
    system2("Rscript", args = args, stdout = log_file, stderr = log_file),
    error = function(e) structure(1L, message = conditionMessage(e))
  )

  if (identical(status, 0L) && file.exists(result_file)) {
    row <- summarize_result(result_file, cfg, plot_file, log_file)
    row$status <- 0L
    row$error_message <- NA_character_
  } else {
    row <- data.frame(
      label = label,
      method = cfg$method,
      results_file = normalizePath(result_file, winslash = "/", mustWork = FALSE),
      plot_file = normalizePath(plot_file, winslash = "/", mustWork = FALSE),
      log_file = normalizePath(log_file, winslash = "/", mustWork = FALSE),
      log_evidence = NA_real_,
      mcse_log_evidence = NA_real_,
      outer_rounds = NA_integer_,
      status = as.integer(status),
      error_message = attr(status, "message") %||% sprintf("Rscript exited with status %s", status),
      stringsAsFactors = FALSE
    )
    if (fail_fast) {
      summary_rows[[length(summary_rows) + 1L]] <- row
      break
    }
  }

  summary_rows[[length(summary_rows) + 1L]] <- row
}

summary_df <- do.call(rbind, summary_rows)
summary_csv <- file.path(results_dir, "shifted_gamma_meta_summary.csv")
summary_rds <- file.path(results_dir, "shifted_gamma_meta_summary.rds")
utils::write.csv(summary_df, summary_csv, row.names = FALSE)
saveRDS(
  list(
    summary = summary_df,
    configs = configs
  ),
  summary_rds
)

ok <- summary_df$status == 0L & file.exists(summary_df$results_file)
if (any(ok)) {
  evidence_plot <- file.path(results_dir, "shifted_gamma_meta_log_evidence.png")
  mu_shift_plot <- file.path(results_dir, "shifted_gamma_meta_mu_shift.png")
  plot_evidence_summary(summary_df[ok, , drop = FALSE], evidence_plot)
  plot_mu_shift_overlay(summary_df$results_file[ok], summary_df$label[ok], mu_shift_plot)
  cat("Saved evidence plot to:", evidence_plot, "\n")
  cat("Saved mu_shift overlay plot to:", mu_shift_plot, "\n")
}

cat("Saved summary CSV to:", summary_csv, "\n")
cat("Saved summary RDS to:", summary_rds, "\n")
