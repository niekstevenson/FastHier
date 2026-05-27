#!/usr/bin/env Rscript

rm(list = ls())

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[1L]))
} else {
  normalizePath("benchmarks/benchmark_emc_certification_rounds.R")
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
  tolower(as.character(val)) %in% c("1", "true", "t", "yes", "y")
}

arg_int_vec <- function(args, key, default) {
  val <- args[[key]]
  if (is.null(val) || !nzchar(val)) return(as.integer(default))
  as.integer(strsplit(val, ",", fixed = TRUE)[[1L]])
}

active_chart_count <- function(factor_set) {
  sum(vapply(factor_set$atlases, function(atlas) {
    sum(vapply(atlas$charts, function(chart) identical(chart$status, "active"), logical(1)))
  }, integer(1)))
}

certification_state <- function(factor_set, theta, cores) {
  parts <- .local_atlas_factor_set_by_local(factor_set, theta, n_cores = as.integer(cores))
  bad <- .local_atlas_factor_set_uncertified(parts)
  reasons <- if (nrow(bad) && "reason" %in% names(bad)) {
    tab <- sort(table(bad$reason), decreasing = TRUE)
    paste(sprintf("%s:%d", names(tab), as.integer(tab)), collapse = ";")
  } else {
    ""
  }
  list(uncertified = bad, n_uncertified = nrow(bad), reason_counts = reasons)
}

history_stats <- function(history, round_id) {
  if (!length(history)) {
    return(data.frame(
      certification_round = as.integer(round_id),
      fresh_probes = 0L,
      selected_pairs = 0L,
      activated_pairs = 0L,
      n_uncertified_before = NA_integer_,
      check.names = FALSE
    ))
  }
  history <- do.call(rbind, history)
  if (!nrow(history)) {
    return(history_stats(list(), round_id))
  }
  fresh_probed <- if ("fresh_probed" %in% names(history)) {
    !is.na(history$fresh_probed) & history$fresh_probed
  } else {
    rep(FALSE, nrow(history))
  }
  selected <- if ("selected" %in% names(history)) {
    !is.na(history$selected) & history$selected
  } else {
    rep(FALSE, nrow(history))
  }
  activated <- if ("activation_success" %in% names(history)) {
    !is.na(history$activation_success) & history$activation_success
  } else {
    rep(FALSE, nrow(history))
  }
  data.frame(
    certification_round = as.integer(round_id),
    fresh_probes = sum(fresh_probed),
    selected_pairs = sum(selected),
    activated_pairs = sum(activated),
    n_uncertified_before = suppressWarnings(as.integer(history$n_uncertified_pairs[1L])),
    check.names = FALSE
  )
}

metric_summary_row <- function(comparison) {
  out <- .local_atlas_metric_summary(comparison)
  out[1L, , drop = FALSE]
}

write_csv <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, file, row.names = FALSE)
}

plot_summary <- function(summary, file) {
  rows <- summary[is.finite(summary$log_evidence), , drop = FALSE]
  if (!nrow(rows)) return(invisible(FALSE))
  draw_panel <- function(x, y, xlab, ylab, main) {
    ok <- is.finite(x) & is.finite(y)
    if (!any(ok)) {
      plot.new()
      title(main = main)
      text(0.5, 0.5, "not available")
      return(invisible(FALSE))
    }
    plot(x[ok], y[ok], type = "b", pch = 19, xlab = xlab, ylab = ylab, main = main)
  }
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  grDevices::png(file, width = 1600, height = 1200)
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    grDevices::dev.off()
  }, add = TRUE)
  par(mfrow = c(2L, 2L), mar = c(4.2, 4.2, 2.5, 1.0))
  draw_panel(rows$certification_round, rows$mean_shape_error,
             "Theta-cloud certification round", "Mean shape error", "Posterior Shape")
  draw_panel(rows$cumulative_fresh_probes, rows$mean_shape_error,
             "Cumulative fresh probes", "Mean shape error", "Quality per Probe")
  draw_panel(rows$certification_round, rows$log_evidence,
             "Theta-cloud certification round", "Outer log evidence", "Evidence")
  draw_panel(rows$certification_round, rows$n_uncertified_after,
             "Theta-cloud certification round", "Uncertified pairs", "Certification")
  invisible(TRUE)
}

cli_args <- parse_cli_args(commandArgs(trailingOnly = TRUE))

checkpoint_file <- arg_chr(cli_args, "checkpoint_file")
if (is.null(checkpoint_file) || !nzchar(checkpoint_file)) {
  stop("Provide --checkpoint_file=<pre-outer chart atlas checkpoint RDS>.")
}
checkpoint_file <- normalizePath(checkpoint_file)

detected_cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
if (!is.finite(detected_cores) || detected_cores < 1L) detected_cores <- 1L
cores <- arg_int(cli_args, "mc_cores", arg_int(cli_args, "cores", min(4L, detected_cores)))

label <- arg_chr(cli_args, "label", "emc_certification_rounds")
out_dir <- arg_chr(cli_args, "out_dir", file.path("benchmarks", "results"))
out_prefix <- file.path(out_dir, label)
data_file <- arg_chr(cli_args, "data_file", file.path("benchmarks", "samples", "full_EMC2.RData"))
rounds_to_compare <- sort(unique(arg_int_vec(cli_args, "outer_rounds", c(0L, 1L, 2L, 4L, 8L, 12L))))
max_certification_rounds <- arg_int(cli_args, "max_certification_rounds", max(rounds_to_compare))
resume_certification <- arg_lgl(cli_args, "resume_certification", FALSE)
resume_outer <- arg_lgl(cli_args, "resume_outer", FALSE)
compare_emc <- arg_lgl(cli_args, "compare_emc", TRUE)
verbose <- arg_lgl(cli_args, "verbose", FALSE)
trace_verbose <- arg_lgl(cli_args, "trace_verbose", verbose)
seed <- arg_int(cli_args, "seed", 20260526L)
n_draws <- arg_int(cli_args, "n_draws", NA_integer_)

certification_particles <- arg_int(cli_args, "certification_particles", NA_integer_)
certification_mcmc_moves <- arg_int(cli_args, "certification_mcmc_moves", NA_integer_)
certification_max_steps <- arg_int(cli_args, "certification_max_steps", NA_integer_)
certification_max_updates <- arg_int(cli_args, "certification_max_updates", NA_integer_)

outer_particles <- arg_int(cli_args, "outer_particles", NA_integer_)
outer_mcmc_moves <- arg_int(cli_args, "outer_mcmc_moves", NA_integer_)
outer_min_mcmc_moves <- arg_int(cli_args, "outer_min_mcmc_moves", NA_integer_)
outer_max_rounds <- arg_int(cli_args, "outer_max_rounds", NA_integer_)
outer_resample_threshold <- arg_num(cli_args, "outer_resample_threshold", NA_real_)
outer_rw_scale_init <- arg_num(cli_args, "outer_rw_scale_init", NA_real_)
outer_stop_on_uncertified <- arg_lgl(cli_args, "outer_stop_on_uncertified", FALSE)
outer_use_uncertified_estimates <- arg_lgl(cli_args, "outer_use_uncertified_estimates", TRUE)

normalizer_robust_method <- arg_chr(cli_args, "normalizer_robust_method", "student_t")
normalizer_student_t_df <- arg_num(cli_args, "normalizer_student_t_df", 30)
if (!normalizer_robust_method %in% c("student_t", "huber", "none")) {
  stop("normalizer_robust_method must be one of: student_t, huber, none.")
}

source("local_charts.R")
source("outer_population_smc.R")
source("utilities.R")

options(
  local_charts.normalizer_robust = !identical(normalizer_robust_method, "none"),
  local_charts.normalizer_robust_method = normalizer_robust_method,
  local_charts.normalizer_student_t_df = normalizer_student_t_df
)

checkpoint <- readRDS(checkpoint_file)
if (is.null(checkpoint$factor_set)) stop("checkpoint_file does not contain factor_set.")
if (is.null(checkpoint$initial_proposal)) stop("checkpoint_file does not contain initial_proposal.")
factor_set <- validate_local_atlas_factor_set(checkpoint$factor_set)
population_model <- checkpoint$population_model %||% factor_set$population_model
initial_proposal <- normalize_theta_proposal(checkpoint$initial_proposal, population_model = population_model)
settings <- checkpoint$settings %||% list()
local_control <- settings$local_control %||% list()
edge_control <- settings$edge_control %||% list()
calibration_control <- settings$calibration_control %||% list()
outer_control <- settings$outer_control %||% list()

if (is.finite(certification_particles)) {
  calibration_control$M <- as.integer(certification_particles)
  calibration_control$confirmation_M <- as.integer(certification_particles)
}
if (is.finite(certification_mcmc_moves)) {
  calibration_control$n_mcmc_moves <- as.integer(certification_mcmc_moves)
}
if (is.finite(certification_max_steps)) {
  calibration_control$max_steps <- as.integer(certification_max_steps)
}
if (!is.finite(certification_max_updates)) {
  certification_max_updates <- calibration_control$initial_certification_max_updates %||%
    calibration_control$pre_outer_max_updates %||% 128L
}

if (!is.finite(outer_particles)) outer_particles <- outer_control$N %||% 500L
if (!is.finite(outer_mcmc_moves)) outer_mcmc_moves <- outer_control$n_mcmc_moves %||% 2L
if (!is.finite(outer_min_mcmc_moves)) outer_min_mcmc_moves <- outer_control$min_mcmc_moves %||% 1L
if (!is.finite(outer_max_rounds)) outer_max_rounds <- outer_control$max_rounds %||% 80L
if (!is.finite(outer_resample_threshold)) outer_resample_threshold <- outer_control$resample_threshold %||% 0.5
if (!is.finite(outer_rw_scale_init)) outer_rw_scale_init <- outer_control$rw_scale_init %||% 0.8

need_data <- max_certification_rounds > 0L || isTRUE(compare_emc)
emc_draws <- NULL
data_list <- NULL
loglik_emc2 <- NULL
if (isTRUE(need_data)) {
  if (!requireNamespace("EMC2", quietly = TRUE)) {
    stop("EMC2 package is required for EMC certification/comparison.")
  }
  if (!file.exists(data_file)) stop("Missing EMC2 benchmark data: ", data_file)
  load(data_file)
  if (!exists("ELP_DDM", inherits = FALSE)) stop("The EMC2 data file must define ELP_DDM.")
  emc <- ELP_DDM[[1L]]
  data_list <- emc$data[seq_len(factor_set$n_locals)]
  alpha_names <- emc$par_names
  model_factory <- emc$model
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
  if (isTRUE(compare_emc)) {
    emc_mu <- as.data.frame(EMC2::parameters(ELP_DDM, selection = "mu"), check.names = FALSE)
    emc_sigma2 <- as.data.frame(EMC2::parameters(ELP_DDM, selection = "sigma2"), check.names = FALSE)
    colnames(emc_mu) <- paste0("mu_", alpha_names)
    colnames(emc_sigma2) <- paste0("sigma2_", alpha_names)
    emc_draws <- data.frame(emc_mu, emc_sigma2, check.names = FALSE)
  }
}

theta_cloud <- theta_proposal_sample(initial_proposal, n = as.integer(outer_particles), seed = seed)
theta_weights <- rep(1 / nrow(theta_cloud), nrow(theta_cloud))

summary_file <- paste0(out_prefix, "_summary.csv")
history_file <- paste0(out_prefix, "_calibration_history.csv")
summary_plot_file <- paste0(out_prefix, "_summary.png")
combined_results_file <- paste0(out_prefix, "_results.rds")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

round_rows <- list()
history_rows <- list()
round_results <- list()
cumulative_fresh_probes <- 0L

save_round_checkpoint <- function(round_id, factor_set, calibration_history, state, round_stats) {
  file <- sprintf("%s_theta_round_%02d_checkpoint.rds", out_prefix, as.integer(round_id))
  saveRDS(
    list(
      source_checkpoint = checkpoint_file,
      source_stage = checkpoint$stage %||% NA_character_,
      stage = sprintf("theta_cloud_round_%02d", as.integer(round_id)),
      factor_set = factor_set,
      atlases = factor_set$atlases,
      population_model = population_model,
      initial_proposal = initial_proposal,
      theta_cloud = theta_cloud,
      certification_state = state,
      calibration_history = calibration_history,
      round_stats = round_stats,
      settings = list(
        label = label,
        seed = seed,
        certification_particles = calibration_control$M %||% local_control$candidate_M %||% NA_integer_,
        certification_mcmc_moves = calibration_control$n_mcmc_moves %||% local_control$n_mcmc_moves %||% NA_integer_,
        certification_max_steps = calibration_control$max_steps %||% local_control$max_steps %||% NA_integer_,
        certification_max_updates = certification_max_updates,
        outer_particles = outer_particles,
        normalizer_robust_method = normalizer_robust_method,
        normalizer_student_t_df = normalizer_student_t_df
      ),
      checkpoint_time = Sys.time()
    ),
    file
  )
  file
}

run_outer_comparison <- function(round_id, factor_set) {
  result_file <- sprintf("%s_round_%02d_outer_results.rds", out_prefix, as.integer(round_id))
  comparison_file <- sprintf("%s_round_%02d_posterior_comparison.csv", out_prefix, as.integer(round_id))
  plot_file <- sprintf("%s_round_%02d_posteriors.png", out_prefix, as.integer(round_id))
  if (isTRUE(resume_outer) && file.exists(result_file)) {
    return(readRDS(result_file))
  }
  outer_factor_set <- factor_set
  outer_factor_set$evaluator_control <- modifyList(outer_factor_set$evaluator_control, list(
    stop_on_uncertified = isTRUE(outer_stop_on_uncertified),
    use_uncertified_estimates = isTRUE(outer_use_uncertified_estimates)
  ))
  outer_factor_set <- validate_local_atlas_factor_set(outer_factor_set)
  t0 <- Sys.time()
  fit <- outer_population_smc(
    factor_set = outer_factor_set,
    N = as.integer(outer_particles),
    initial_proposal = initial_proposal,
    resample_threshold = outer_resample_threshold,
    n_mcmc_moves = as.integer(outer_mcmc_moves),
    min_mcmc_moves = as.integer(outer_min_mcmc_moves),
    max_rounds = as.integer(outer_max_rounds),
    rw_scale_init = outer_rw_scale_init,
    n_cores = as.integer(cores),
    seed = as.integer(seed) + 7000003L + as.integer(round_id),
    verbose = verbose
  )
  elapsed_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  workflow_draws <- NULL
  posterior_comparison <- NULL
  posterior_summary <- NULL
  if (isTRUE(compare_emc)) {
    draw_n <- if (is.finite(n_draws)) as.integer(n_draws) else nrow(emc_draws)
    workflow_draws <- local_atlas_draws_from_fit(
      fit,
      population_model = population_model,
      n_draws = draw_n,
      seed = as.integer(seed) + 7100003L + as.integer(round_id)
    )
    posterior_comparison <- local_atlas_compare_posterior_draws(emc_draws, workflow_draws)
    posterior_summary <- metric_summary_row(posterior_comparison)
    write_csv(posterior_comparison, comparison_file)
    grDevices::png(plot_file, width = 1800, height = 1400)
    plot_posteriors(
      emc_draws,
      workflow_draws,
      labels = c("EMC2", sprintf("round_%02d", as.integer(round_id))),
      cols = c("black", "firebrick3"),
      n_cols = 4L
    )
    grDevices::dev.off()
  }
  result <- list(
    round = as.integer(round_id),
    fit = fit,
    workflow_draws = workflow_draws,
    posterior_comparison = posterior_comparison,
    posterior_summary = posterior_summary,
    elapsed_sec = elapsed_sec,
    log_evidence = as.numeric(fit$log_evidence %||% NA_real_),
    result_file = result_file,
    comparison_file = if (isTRUE(compare_emc)) comparison_file else NULL,
    plot_file = if (isTRUE(compare_emc)) plot_file else NULL
  )
  saveRDS(result, result_file)
  result
}

record_round <- function(round_id, factor_set, state, stats, checkpoint_path, outer_result = NULL) {
  row <- data.frame(
    certification_round = as.integer(round_id),
    active_charts = active_chart_count(factor_set),
    n_uncertified_after = state$n_uncertified,
    uncertified_reasons_after = state$reason_counts,
    fresh_probes = stats$fresh_probes,
    cumulative_fresh_probes = cumulative_fresh_probes,
    selected_pairs = stats$selected_pairs,
    activated_pairs = stats$activated_pairs,
    n_uncertified_before = stats$n_uncertified_before,
    checkpoint_file = checkpoint_path,
    outer_elapsed_sec = NA_real_,
    log_evidence = NA_real_,
    n_parameters = NA_integer_,
    max_abs_standardized_mean_error = NA_real_,
    mean_abs_standardized_mean_error = NA_real_,
    max_scaled_q_wasserstein = NA_real_,
    mean_scaled_q_wasserstein = NA_real_,
    max_shape_error = NA_real_,
    mean_shape_error = NA_real_,
    posterior_plot_file = NA_character_,
    check.names = FALSE
  )
  if (!is.null(outer_result)) {
    row$outer_elapsed_sec <- outer_result$elapsed_sec
    row$log_evidence <- outer_result$log_evidence
    if (!is.null(outer_result$posterior_summary)) {
      ps <- outer_result$posterior_summary[1L, , drop = FALSE]
      for (nm in names(ps)) row[[nm]] <- ps[[nm]]
    }
    row$posterior_plot_file <- outer_result$plot_file %||% NA_character_
  }
  row
}

cat(sprintf("Certification-round benchmark: %s\n", label))
cat(sprintf("Checkpoint: %s | stage=%s | locals=%d\n", checkpoint_file, checkpoint$stage %||% "", factor_set$n_locals))
cat(sprintf(
  "Theta cloud / outer N=%d | certification rounds=%d | compare rounds=%s\n",
  outer_particles,
  max_certification_rounds,
  paste(rounds_to_compare, collapse = ",")
))
cat(sprintf(
  "Certification M=%s | max_updates=%d | normalizer=%s\n",
  as.character(calibration_control$M %||% local_control$candidate_M %||% NA_integer_),
  certification_max_updates,
  normalizer_robust_method
))

state0 <- certification_state(factor_set, theta_cloud, cores)
stats0 <- history_stats(list(), 0L)
checkpoint0 <- save_round_checkpoint(0L, factor_set, list(), state0, stats0)
outer0 <- NULL
if (0L %in% rounds_to_compare) {
  cat("Running outer comparison for certification round 0\n")
  outer0 <- run_outer_comparison(0L, factor_set)
  round_results[["0"]] <- outer0
}
round_rows[[length(round_rows) + 1L]] <- record_round(0L, factor_set, state0, stats0, checkpoint0, outer0)
write_csv(do.call(rbind, round_rows), summary_file)
plot_summary(do.call(rbind, round_rows), summary_plot_file)

for (round_id in seq_len(as.integer(max_certification_rounds))) {
  round_checkpoint <- sprintf("%s_theta_round_%02d_checkpoint.rds", out_prefix, as.integer(round_id))
  round_history <- list()
  round_stats <- history_stats(list(), round_id)
  if (isTRUE(resume_certification) && file.exists(round_checkpoint)) {
    loaded <- readRDS(round_checkpoint)
    factor_set <- validate_local_atlas_factor_set(loaded$factor_set)
    state <- loaded$certification_state %||% certification_state(factor_set, theta_cloud, cores)
    round_stats <- loaded$round_stats %||% round_stats
    round_history <- loaded$calibration_history %||% list()
    cat(sprintf("Loaded certification round %d checkpoint: %s\n", round_id, round_checkpoint))
  } else {
    before <- certification_state(factor_set, theta_cloud, cores)
    if (before$n_uncertified > 0L) {
      cat(sprintf(
        "Certifying round %d: uncertified=%d | max_updates=%d\n",
        round_id,
        before$n_uncertified,
        certification_max_updates
      ))
      certification <- local_atlas_certify_theta_cloud(
        factor_set = factor_set,
        theta = theta_cloud,
        data_list = data_list,
        loglik_fn = loglik_emc2,
        local_control = local_control,
        edge_control = edge_control,
        calibration_control = calibration_control,
        theta_weights = theta_weights,
        max_rounds = 1L,
        max_updates = as.integer(certification_max_updates),
        n_cores = as.integer(cores),
        seed = as.integer(seed) + 6100003L * as.integer(round_id),
        verbose = verbose,
        trace_verbose = trace_verbose
      )
      factor_set <- validate_local_atlas_factor_set(certification$factor_set)
      round_history <- certification$calibration_history %||% list()
      if (length(round_history)) {
        round_history <- lapply(round_history, function(x) {
          if (is.data.frame(x) && nrow(x)) {
            x$calibration_round <- as.integer(round_id)
          }
          x
        })
      }
      state <- certification_state(factor_set, theta_cloud, cores)
      round_stats <- history_stats(round_history, round_id)
      round_stats$n_uncertified_before <- before$n_uncertified
    } else {
      cat(sprintf("Certification round %d skipped: no uncertified pairs\n", round_id))
      state <- before
    }
    round_checkpoint <- save_round_checkpoint(round_id, factor_set, round_history, state, round_stats)
  }
  cumulative_fresh_probes <- cumulative_fresh_probes + as.integer(round_stats$fresh_probes)
  if (length(round_history)) {
    history_rows <- c(history_rows, round_history)
    write_csv(do.call(rbind, history_rows), history_file)
  }
  outer_result <- NULL
  if (round_id %in% rounds_to_compare) {
    cat(sprintf("Running outer comparison for certification round %d\n", round_id))
    outer_result <- run_outer_comparison(round_id, factor_set)
    round_results[[as.character(round_id)]] <- outer_result
  }
  round_rows[[length(round_rows) + 1L]] <- record_round(
    round_id,
    factor_set,
    state,
    round_stats,
    round_checkpoint,
    outer_result
  )
  current_summary <- do.call(rbind, round_rows)
  write_csv(current_summary, summary_file)
  plot_summary(current_summary, summary_plot_file)
  gc(verbose = FALSE)
}

summary <- do.call(rbind, round_rows)
calibration_history <- if (length(history_rows)) do.call(rbind, history_rows) else data.frame()
saveRDS(
  list(
    source_checkpoint = checkpoint_file,
    source_stage = checkpoint$stage %||% NA_character_,
    summary = summary,
    calibration_history = calibration_history,
    round_results = round_results,
    theta_cloud = theta_cloud,
    settings = list(
      label = label,
      seed = seed,
      outer_rounds = rounds_to_compare,
      max_certification_rounds = max_certification_rounds,
      outer_particles = outer_particles,
      outer_mcmc_moves = outer_mcmc_moves,
      outer_min_mcmc_moves = outer_min_mcmc_moves,
      outer_max_rounds = outer_max_rounds,
      certification_particles = calibration_control$M %||% local_control$candidate_M %||% NA_integer_,
      certification_mcmc_moves = calibration_control$n_mcmc_moves %||% local_control$n_mcmc_moves %||% NA_integer_,
      certification_max_steps = calibration_control$max_steps %||% local_control$max_steps %||% NA_integer_,
      certification_max_updates = certification_max_updates,
      normalizer_robust_method = normalizer_robust_method,
      normalizer_student_t_df = normalizer_student_t_df,
      cores = cores
    )
  ),
  combined_results_file
)

cat(sprintf("Saved summary: %s\n", summary_file))
if (nrow(calibration_history)) cat(sprintf("Saved calibration history: %s\n", history_file))
cat(sprintf("Saved summary plot: %s\n", summary_plot_file))
cat(sprintf("Saved combined results: %s\n", combined_results_file))
cat("\nRound summary:\n")
print(summary[, intersect(c(
  "certification_round",
  "active_charts",
  "n_uncertified_after",
  "fresh_probes",
  "cumulative_fresh_probes",
  "log_evidence",
  "mean_abs_standardized_mean_error",
  "mean_scaled_q_wasserstein",
  "mean_shape_error"
), names(summary)), drop = FALSE], row.names = FALSE)
