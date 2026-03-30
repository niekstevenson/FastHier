rm(list = ls())

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[1L]))
} else {
  normalizePath("benchmarks/compare_population_EMC.R")
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

arg_int <- function(args, key, default) {
  val <- args[[key]]
  if (is.null(val) || !nzchar(val)) return(as.integer(default))
  as.integer(val)
}

arg_lgl <- function(args, key, default = FALSE) {
  val <- args[[key]]
  if (is.null(val) || !nzchar(val)) return(isTRUE(default))
  tolower(as.character(val)) %in% c("1", "true", "t", "yes", "y")
}

config_label_default <- function(refined_method, transport_method, gss_enable, da_enable) {
  refined_tag <- if (identical(refined_method, "defensive_mixture")) "defmix" else "broad"
  transport_tag <- if (identical(transport_method, "gaussian_copula")) "gcop" else "tri"
  mode_tag <- if (gss_enable && da_enable) {
    "gss_da"
  } else if (gss_enable) {
    "gss"
  } else if (da_enable) {
    "da"
  } else {
    "base"
  }
  paste(refined_tag, transport_tag, mode_tag, sep = "_")
}

build_config_grid <- function(config_set = "full") {
  config_set <- match.arg(config_set, c("full", "core", "gss"))
  refined_methods <- c("defensive_mixture", "broadened_gaussian")
  transport_methods <- c("gaussian_copula", "sparse_triangular")

  configs <- switch(
    config_set,
    full = expand.grid(
      refined_method = refined_methods,
      transport_method = transport_methods,
      gss_enable = c(FALSE, TRUE),
      da_enable = c(FALSE, TRUE),
      stringsAsFactors = FALSE
    ),
    core = expand.grid(
      refined_method = refined_methods,
      transport_method = transport_methods,
      gss_enable = FALSE,
      da_enable = FALSE,
      stringsAsFactors = FALSE
    ),
    gss = expand.grid(
      refined_method = refined_methods,
      transport_method = transport_methods,
      gss_enable = c(FALSE, TRUE),
      da_enable = c(FALSE, TRUE),
      stringsAsFactors = FALSE
    )
  )

  configs$hist_mix_enable <- configs$gss_enable | configs$da_enable
  configs$label <- mapply(
    config_label_default,
    configs$refined_method,
    configs$transport_method,
    configs$gss_enable,
    configs$da_enable,
    USE.NAMES = FALSE
  )

  if (identical(config_set, "gss")) {
    configs <- configs[configs$gss_enable | configs$da_enable, , drop = FALSE]
  }

  configs[order(configs$label), , drop = FALSE]
}

draw_quantiles <- function(x, probs = c(0.1, 0.5, 0.9)) {
  stats::quantile(as.numeric(x), probs = probs, na.rm = TRUE, names = FALSE)
}

run_with_log <- function(log_file, expr) {
  expr <- substitute(expr)
  con <- file(log_file, open = "wt")
  sink(con, split = FALSE)
  sink(con, type = "message")
  on.exit({
    sink(type = "message")
    sink()
    close(con)
  }, add = TRUE)
  eval.parent(expr)
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
    main = "Population EMC Evidence By Configuration"
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

plot_error_summary <- function(summary_df, file) {
  ok <- is.finite(summary_df$mean_abs_posterior_diff)
  if (!any(ok)) return(invisible(NULL))
  x <- summary_df[ok, , drop = FALSE]
  ord <- order(x$mean_abs_posterior_diff)
  x <- x[ord, , drop = FALSE]
  y <- seq_len(nrow(x))

  grDevices::png(file, width = 1400, height = max(900, 55 * nrow(x)))
  graphics::par(mar = c(5, 18, 3, 2))
  graphics::plot(
    x$mean_abs_posterior_diff, y,
    ylim = c(0.5, nrow(x) + 0.5),
    yaxt = "n",
    pch = 19,
    col = "steelblue4",
    xlab = "Mean absolute posterior diff",
    ylab = "",
    main = "Population EMC Posterior Error By Configuration"
  )
  graphics::axis(2, at = y, labels = x$label, las = 1, cex.axis = 0.9)
  grDevices::dev.off()
}

plot_config_posteriors <- function(workflow_draws, emc_draws, label, file) {
  n_params <- ncol(workflow_draws)
  n_cols <- min(4L, n_params)
  n_rows <- ceiling(n_params / n_cols)
  width <- max(1600L, 420L * n_cols)
  height <- max(1000L, 300L * n_rows)

  grDevices::png(file, width = width, height = height)
  plot_posteriors(
    workflow_draws,
    emc_draws,
    labels = c(label, "EMC2"),
    cols = c("firebrick3", "black"),
    n_cols = n_cols
  )
  grDevices::dev.off()
}

summarize_result <- function(result_file, config, plot_file, log_file) {
  res <- readRDS(result_file)
  workflow_draws <- res$workflow_draws
  emc_draws <- res$emc_draws
  common <- intersect(names(workflow_draws), names(emc_draws))

  row <- data.frame(
    label = config$label,
    refined_method = config$refined_method,
    transport_method = config$transport_method,
    hist_mix_enable = config$hist_mix_enable,
    gss_enable = config$gss_enable,
    da_enable = config$da_enable,
    results_file = normalizePath(result_file, winslash = "/", mustWork = FALSE),
    plot_file = normalizePath(plot_file, winslash = "/", mustWork = FALSE),
    log_file = normalizePath(log_file, winslash = "/", mustWork = FALSE),
    elapsed_sec = as.numeric(res$elapsed_sec %||% res$settings$elapsed_sec %||% NA_real_),
    log_evidence = as.numeric(res$fit$log_evidence %||% NA_real_),
    mcse_log_evidence = as.numeric(res$fit$mcse_log_evidence %||% NA_real_),
    outer_rounds = as.integer(res$fit$meta$rounds %||% NA_integer_),
    stringsAsFactors = FALSE
  )

  diff_cols <- character()
  for (nm in common) {
    wf <- workflow_draws[[nm]]
    ref <- emc_draws[[nm]]
    qwf <- draw_quantiles(wf)
    qref <- draw_quantiles(ref)
    mean_col <- paste0(nm, "_mean_diff")
    sd_col <- paste0(nm, "_sd_diff")
    q10_col <- paste0(nm, "_q10_diff")
    q50_col <- paste0(nm, "_q50_diff")
    q90_col <- paste0(nm, "_q90_diff")
    row[[mean_col]] <- mean(wf) - mean(ref)
    row[[sd_col]] <- stats::sd(wf) - stats::sd(ref)
    row[[q10_col]] <- qwf[1L] - qref[1L]
    row[[q50_col]] <- qwf[2L] - qref[2L]
    row[[q90_col]] <- qwf[3L] - qref[3L]
    diff_cols <- c(diff_cols, mean_col, sd_col, q10_col, q50_col, q90_col)
  }

  diff_values <- unlist(row[diff_cols], use.names = FALSE)
  row$mean_abs_posterior_diff <- mean(abs(diff_values))
  row$max_abs_posterior_diff <- max(abs(diff_values))
  row
}

cli_args <- parse_cli_args(commandArgs(trailingOnly = TRUE))
config_set <- arg_chr(cli_args, "config_set", "full")
results_dir <- arg_chr(cli_args, "results_dir", file.path("benchmarks", "results", "meta", "emc"))
labels_filter <- arg_chr(cli_args, "labels", NULL)
inner_verbose <- arg_lgl(cli_args, "inner_verbose", FALSE)
fail_fast <- arg_lgl(cli_args, "fail_fast", FALSE)

suppressPackageStartupMessages({
  library(parallel)
  library(EMC2)
})

detected_cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
if (!is.finite(detected_cores) || detected_cores < 1L) {
  detected_cores <- 1L
}

mc.cores <- as.integer(max(1L, min(arg_int(cli_args, "mc_cores", 4L), detected_cores)))
pilot_size <- arg_int(cli_args, "pilot_size", 30L)
pilot_particles <- arg_int(cli_args, "pilot_particles", 800L)
full_particles <- arg_int(cli_args, "full_particles", 2000L)
outer_particles <- arg_int(cli_args, "outer_particles", 2000L)
outer_mcmc_moves <- arg_int(cli_args, "outer_mcmc_moves", 3L)
outer_max_rounds <- arg_int(cli_args, "outer_max_rounds", 80L)
base_seed <- arg_int(cli_args, "base_seed", 20260324L)
broad_scale <- as.numeric(arg_chr(cli_args, "broad_scale", "1"))

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
logs_dir <- file.path(results_dir, "logs")
dir.create(logs_dir, showWarnings = FALSE, recursive = TRUE)

source("smc_core.R")
source("reference_priors.R")
source("utilities.R")
source("SMC_super_fast.R")
fit_copula_transform_sparse <- fit_copula_transform
source("hierarchical_locals.R")
source("population_models.R")
source("outer_population_smc.R")

data_path <- file.path("benchmarks", "samples", "full_EMC2.RData")
if (!file.exists(data_path)) {
  stop("Missing data file: ", data_path)
}

load(data_path)
if (!exists("ELP_DDM", inherits = FALSE)) {
  stop("The data file does not define an object named 'ELP_DDM'.")
}

emc <- ELP_DDM[[1L]]
data_list <- emc$data
model_factory <- emc$model
param_names <- emc$par_names
d <- length(param_names)

base_mu <- stats::setNames(rep(0, d), param_names)
base_Sigma <- diag(1, nrow = d, ncol = d)
rownames(base_Sigma) <- colnames(base_Sigma) <- param_names

population_model <- make_population_model_diag_gaussian(
  alpha_names = param_names,
  mean_prior_mean = rep(0, d),
  mean_prior_var = rep(1, d),
  sigma2_prior_shape = rep(2, d),
  sigma2_prior_rate = rep(0.3, d),
  label = "emc_normal_gamma"
)

loglik_emc2 <- function(Theta, data_i) {
  Theta <- as.matrix(Theta)
  colnames(Theta) <- param_names
  out <- as.numeric(EMC2:::calc_ll_manager(Theta, data_i, model_factory, r_cores = 1L))
  out[!is.finite(out) | out > 100] <- min(out)
  out
}

reference_mu_draws <- parameters(ELP_DDM, selection = "mu")
reference_sigma2_draws <- parameters(ELP_DDM, selection = "sigma2")
colnames(reference_mu_draws) <- paste0("mu_", param_names)
colnames(reference_sigma2_draws) <- paste0("sigma2_", param_names)
emc_draws <- data.frame(reference_mu_draws, reference_sigma2_draws, check.names = FALSE)

configs <- build_config_grid(config_set)
if (!is.null(labels_filter)) {
  keep <- trimws(strsplit(labels_filter, ",", fixed = TRUE)[[1L]])
  configs <- configs[configs$label %in% keep, , drop = FALSE]
}
if (!nrow(configs)) {
  stop("No configurations selected.")
}

summary_rows <- list()

for (i in seq_len(nrow(configs))) {
  cfg <- configs[i, , drop = FALSE]
  label <- cfg$label
  result_file <- file.path(results_dir, sprintf("population_emc_%s_results.rds", label))
  plot_file <- file.path(results_dir, sprintf("population_emc_%s_posteriors.png", label))
  log_file <- file.path(logs_dir, sprintf("population_emc_%s.log", label))

  cat(sprintf("[%d/%d] %s\n", i, nrow(configs), label))

  if (identical(cfg$transport_method, "gaussian_copula")) {
    fit_copula_transform <- .fit_gaussian_copula_transform
  } else {
    fit_copula_transform <- fit_copula_transform_sparse
  }

  start_time <- proc.time()[["elapsed"]]
  result <- try(
    run_with_log(log_file, {
      cat(sprintf("Configuration: label=%s | refined=%s | transport=%s | hist_mix=%s | gss=%s | da=%s\n",
                  label,
                  cfg$refined_method,
                  cfg$transport_method,
                  cfg$hist_mix_enable,
                  cfg$gss_enable,
                  cfg$da_enable))
      cat(sprintf("Data: %d subjects\n", length(data_list)))
      cat("Running local-reference stage...\n")

      stage <- prepare_reference_local_stage(
        data_list = data_list,
        loglik_fn = loglik_emc2,
        base_mu = base_mu,
        base_Sigma = base_Sigma,
        pilot_size = min(pilot_size, length(data_list)),
        broad_scale = broad_scale,
        pilot_particles = pilot_particles,
        full_particles = full_particles,
        refined_method = cfg$refined_method,
        n_jobs = mc.cores,
        base_seed = base_seed,
        pilot_population_model = population_model,
        verbose = inner_verbose,
        pilot_smc_control = list(
          max_rounds = 40L,
          hist_mix_enable = cfg$hist_mix_enable,
          gss_enable = cfg$gss_enable,
          da_enable = cfg$da_enable
        ),
        full_smc_control = list(
          hist_mix_enable = cfg$hist_mix_enable,
          gss_enable = cfg$gss_enable,
          da_enable = cfg$da_enable
        )
      )

      cat("Running outer population SMC...\n")
      factor_set <- build_population_factor_set(stage$local_objects, population_model)
      fit <- outer_population_smc(
        factor_set = factor_set,
        N = outer_particles,
        n_mcmc_moves = outer_mcmc_moves,
        max_rounds = outer_max_rounds,
        n_cores = mc.cores,
        seed = base_seed,
        verbose = inner_verbose
      )

      workflow_parts <- smc_posteriors(
        fit,
        n_draws = nrow(reference_mu_draws),
        seed = base_seed + 1L,
        population_model = population_model
      )
      workflow_mu <- as.data.frame(workflow_parts$mu, check.names = FALSE)
      workflow_sigma2 <- as.data.frame(workflow_parts$sigma2, check.names = FALSE)
      colnames(workflow_mu) <- paste0("mu_", param_names)
      colnames(workflow_sigma2) <- paste0("sigma2_", param_names)
      workflow_draws <- data.frame(workflow_mu, workflow_sigma2, check.names = FALSE)

      posterior_summary <- summarize_population_posterior_diag(
        theta = fit$theta,
        w = fit$w,
        model = population_model
      )

      plot_config_posteriors(workflow_draws, emc_draws, label, plot_file)

      elapsed_sec <- proc.time()[["elapsed"]] - start_time
      saveRDS(
        list(
          data_source = normalizePath(data_path, winslash = "/", mustWork = TRUE),
          stage = stage,
          fit = fit,
          emc_draws = emc_draws,
          workflow_draws = workflow_draws,
          posterior_summary = posterior_summary,
          settings = list(
            label = label,
            transport_method = cfg$transport_method,
            refined_method = cfg$refined_method,
            hist_mix_enable = cfg$hist_mix_enable,
            gss_enable = cfg$gss_enable,
            da_enable = cfg$da_enable,
            mc.cores = mc.cores,
            pilot_size = min(pilot_size, length(data_list)),
            pilot_particles = pilot_particles,
            full_particles = full_particles,
            outer_particles = outer_particles,
            outer_mcmc_moves = outer_mcmc_moves,
            outer_max_rounds = outer_max_rounds,
            base_seed = base_seed,
            broad_scale = broad_scale,
            elapsed_sec = elapsed_sec
          ),
          elapsed_sec = elapsed_sec,
          plot_file = plot_file
        ),
        result_file
      )

      cat("Saved results to:", result_file, "\n")
      cat("Saved plot to:", plot_file, "\n")
      invisible(NULL)
    }),
    silent = TRUE
  )
  elapsed_sec <- proc.time()[["elapsed"]] - start_time

  if (!inherits(result, "try-error") && file.exists(result_file)) {
    row <- summarize_result(result_file, cfg, plot_file, log_file)
    row$status <- 0L
    row$error_message <- NA_character_
  } else {
    row <- data.frame(
      label = label,
      refined_method = cfg$refined_method,
      transport_method = cfg$transport_method,
      hist_mix_enable = cfg$hist_mix_enable,
      gss_enable = cfg$gss_enable,
      da_enable = cfg$da_enable,
      results_file = normalizePath(result_file, winslash = "/", mustWork = FALSE),
      plot_file = normalizePath(plot_file, winslash = "/", mustWork = FALSE),
      log_file = normalizePath(log_file, winslash = "/", mustWork = FALSE),
      elapsed_sec = elapsed_sec,
      log_evidence = NA_real_,
      mcse_log_evidence = NA_real_,
      outer_rounds = NA_integer_,
      mean_abs_posterior_diff = NA_real_,
      max_abs_posterior_diff = NA_real_,
      status = 1L,
      error_message = as.character(result),
      stringsAsFactors = FALSE
    )
    if (fail_fast) {
      summary_rows[[length(summary_rows) + 1L]] <- row
      break
    }
  }

  summary_rows[[length(summary_rows) + 1L]] <- row
  gc(verbose = FALSE)
}

summary_df <- do.call(function(...) {
  rows <- list(...)
  all_names <- unique(unlist(lapply(rows, names), use.names = FALSE))
  rows <- lapply(rows, function(row) {
    missing <- setdiff(all_names, names(row))
    for (nm in missing) row[[nm]] <- NA
    row[all_names]
  })
  do.call(rbind, rows)
}, summary_rows)

summary_csv <- file.path(results_dir, "population_emc_meta_summary.csv")
summary_rds <- file.path(results_dir, "population_emc_meta_summary.rds")
evidence_plot <- file.path(results_dir, "population_emc_meta_log_evidence.png")
error_plot <- file.path(results_dir, "population_emc_meta_abs_error.png")

utils::write.csv(summary_df, summary_csv, row.names = FALSE)
saveRDS(summary_df, summary_rds)
plot_evidence_summary(summary_df, evidence_plot)
plot_error_summary(summary_df, error_plot)

cat("Saved summary to:", summary_csv, "\n")
cat("Saved summary RDS to:", summary_rds, "\n")
cat("Saved evidence plot to:", evidence_plot, "\n")
cat("Saved posterior error plot to:", error_plot, "\n")
