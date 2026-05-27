#!/usr/bin/env Rscript

rm(list = ls())

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[1L]))
} else {
  normalizePath("benchmarks/rerun_chart_atlas_outer_from_checkpoint.R")
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

rmse_finite <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x)) sqrt(mean(x * x)) else NA_real_
}

compression_theta_from_checkpoint <- function(checkpoint,
                                              factor_set,
                                              initial_proposal,
                                              n_theta,
                                              seed) {
  model <- factor_set$population_model
  active_anchor_theta <- lapply(factor_set$atlases, function(atlas) {
    active <- Filter(function(chart) identical(chart$status, "active"), atlas$charts)
    if (!length(active)) return(NULL)
    do.call(rbind, lapply(active, function(chart) chart$theta_anchor))
  })
  active_anchor_theta <- active_anchor_theta[!vapply(active_anchor_theta, is.null, logical(1))]
  active_anchor_theta <- if (length(active_anchor_theta)) do.call(rbind, active_anchor_theta) else NULL
  fit_theta <- checkpoint$fit$theta %||% NULL
  if (!is.null(fit_theta) && nrow(as.matrix(fit_theta)) > as.integer(n_theta)) {
    set.seed(as.integer(seed) + 901L)
    fit_theta <- fit_theta[sample.int(nrow(fit_theta), as.integer(n_theta)), , drop = FALSE]
  }
  proposal_n <- max(0L, as.integer(n_theta) - nrow(as.matrix(fit_theta %||% checkpoint$theta_cloud %||% checkpoint$theta_root)))
  proposal_theta <- if (!is.null(initial_proposal) && proposal_n > 0L) {
    theta_proposal_sample(initial_proposal, n = proposal_n, seed = as.integer(seed) + 1009L)
  } else {
    NULL
  }
  pool <- .local_atlas_unique_theta(
    do.call(rbind, Filter(Negate(is.null), list(
      checkpoint$theta_root,
      checkpoint$theta_design,
      checkpoint$theta_cloud,
      active_anchor_theta,
      fit_theta,
      proposal_theta
    ))),
    model
  )
  if (nrow(pool) <= as.integer(n_theta)) return(pool)
  seed_design <- .local_atlas_unique_theta(
    do.call(rbind, Filter(Negate(is.null), list(
      checkpoint$theta_root,
      active_anchor_theta,
      checkpoint$theta_design
    ))),
    model
  )
  seed_keep <- min(nrow(seed_design), max(1L, floor(0.5 * as.integer(n_theta))))
  .local_atlas_metric_farthest_design(
    theta_cloud = pool,
    population_model = model,
    theta_root = checkpoint$theta_root %||% seed_design[1L, , drop = FALSE],
    max_anchors = as.integer(n_theta),
    distance_metric = factor_set$evaluator_control$distance_metric %||% "fisher",
    initial_design = seed_design[seq_len(seed_keep), , drop = FALSE]
  )
}

cli_args <- parse_cli_args(commandArgs(trailingOnly = TRUE))
checkpoint_file <- arg_chr(cli_args, "checkpoint_file")
if (is.null(checkpoint_file) || !nzchar(checkpoint_file)) {
  stop("Provide --checkpoint_file=<chart atlas checkpoint RDS>.")
}
checkpoint_file <- normalizePath(checkpoint_file)
out_prefix <- sub("_checkpoint.*\\.rds$", "", checkpoint_file)
results_file <- arg_chr(cli_args, "results_file", paste0(out_prefix, "_outer_rerun_results.rds"))
post_outer_checkpoint_file <- arg_chr(
  cli_args,
  "post_outer_checkpoint_file",
  paste0(out_prefix, "_checkpoint_post_outer.rds")
)
post_pre_outer_checkpoint_file <- arg_chr(
  cli_args,
  "post_pre_outer_checkpoint_file",
  paste0(out_prefix, "_checkpoint_post_pre_outer.rds")
)
data_file <- arg_chr(cli_args, "data_file", file.path("benchmarks", "samples", "full_EMC2.RData"))
baseline_file <- arg_chr(
  cli_args,
  "baseline_file",
  file.path("benchmarks", "results", "population_emc_framework_defaults_fixed_results.rds")
)
plot_file <- arg_chr(cli_args, "plot_file", paste0(out_prefix, "_outer_rerun_posteriors.png"))
comparison_csv <- arg_chr(cli_args, "comparison_csv", paste0(out_prefix, "_outer_rerun_posterior_comparison.csv"))
compression_summary_csv <- arg_chr(cli_args, "compression_summary_csv", paste0(out_prefix, "_outer_rerun_compression_summary.csv"))
compression_audit_csv <- arg_chr(cli_args, "compression_audit_csv", paste0(out_prefix, "_outer_rerun_compression_audit.csv"))
compare_emc <- arg_lgl(cli_args, "compare_emc", file.exists(data_file))
allow_non_frozen <- arg_lgl(cli_args, "allow_non_frozen", FALSE)

detected_cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
if (!is.finite(detected_cores) || detected_cores < 1L) detected_cores <- 1L
cores <- arg_int(cli_args, "mc_cores", arg_int(cli_args, "cores", min(4L, detected_cores)))

source("local_charts.R")
source("outer_population_smc.R")
source("utilities.R")

checkpoint <- readRDS(checkpoint_file)
checkpoint_stage <- as.character(checkpoint$stage %||% "")
if (!identical(checkpoint_stage, "post_atlas_build") && !isTRUE(allow_non_frozen)) {
  stop(
    "checkpoint_file stage is '", checkpoint_stage,
    "'. Expected frozen pre-outer stage 'post_atlas_build'. ",
    "Use --allow_non_frozen=true only for diagnostics."
  )
}
if (is.null(checkpoint$factor_set)) {
  stop("checkpoint_file does not contain factor_set.")
}
if (is.null(checkpoint$initial_proposal)) {
  stop("checkpoint_file does not contain initial_proposal.")
}
factor_set <- validate_local_atlas_factor_set(checkpoint$factor_set)
override_control <- list()
if (!is.null(cli_args$stop_on_uncertified)) {
  override_control$stop_on_uncertified <- arg_lgl(cli_args, "stop_on_uncertified", factor_set$evaluator_control$stop_on_uncertified)
}
if (!is.null(cli_args$use_uncertified_estimates)) {
  override_control$use_uncertified_estimates <- arg_lgl(cli_args, "use_uncertified_estimates", factor_set$evaluator_control$use_uncertified_estimates)
}
if (!is.null(cli_args$require_particle_mis)) {
  override_control$require_particle_mis <- arg_lgl(cli_args, "require_particle_mis", factor_set$evaluator_control$require_particle_mis)
}
if (!is.null(cli_args$particle_mis_batch)) {
  override_control$particle_mis_batch <- arg_lgl(cli_args, "particle_mis_batch", factor_set$evaluator_control$particle_mis_batch)
}
for (key in c(
  "max_chart_distance",
  "min_particle_mis_ess",
  "min_particle_mis_ess_abs",
  "max_particle_mis_psis_k",
  "max_quadratic_particle_gap",
  "sparse_chart_max_distance",
  "max_leave_chart_out_gap"
)) {
  if (!is.null(cli_args[[key]])) {
    override_control[[key]] <- arg_num(cli_args, key, factor_set$evaluator_control[[key]])
  }
}
if (length(override_control)) {
  factor_set$evaluator_control <- modifyList(factor_set$evaluator_control, override_control)
  factor_set <- validate_local_atlas_factor_set(factor_set)
}
settings <- checkpoint$settings %||% list()
outer_control <- settings$outer_control %||% list()
local_control <- settings$local_control %||% list()
edge_control <- settings$edge_control %||% list()
design_control <- settings$design_control %||% list()
calibration_control <- settings$calibration_control %||% list()

N <- arg_int(cli_args, "outer_particles", outer_control$N %||% 1000L)
n_mcmc_moves <- arg_int(cli_args, "outer_mcmc_moves", outer_control$n_mcmc_moves %||% 3L)
min_mcmc_moves <- arg_int(cli_args, "outer_min_mcmc_moves", outer_control$min_mcmc_moves %||% 1L)
max_rounds <- arg_int(cli_args, "outer_max_rounds", outer_control$max_rounds %||% 80L)
resample_threshold <- arg_num(cli_args, "outer_resample_threshold", outer_control$resample_threshold %||% 0.5)
rw_scale_init <- arg_num(cli_args, "outer_rw_scale_init", outer_control$rw_scale_init %||% 0.8)
seed <- arg_int(cli_args, "seed", (settings$seed %||% 123L) + 9900001L)
verbose <- arg_lgl(cli_args, "verbose", FALSE)
n_draws <- arg_int(cli_args, "n_draws", NA_integer_)
pre_outer_rounds <- arg_int(cli_args, "pre_outer_rounds", calibration_control$pre_outer_rounds %||% 0L)
pre_outer_audit_n <- arg_int(cli_args, "pre_outer_audit_n", calibration_control$pre_outer_audit_n %||% max(N, 100L))
pre_outer_max_points <- arg_int(cli_args, "pre_outer_max_points", calibration_control$pre_outer_max_points %||% 6L)
pre_outer_max_updates <- arg_int(cli_args, "pre_outer_max_updates", calibration_control$pre_outer_max_updates %||% 16L)
pre_outer_particles <- arg_int(cli_args, "pre_outer_particles", calibration_control$M %||% local_control$candidate_M %||% 500L)
pre_outer_mcmc_moves <- arg_int(cli_args, "pre_outer_mcmc_moves", calibration_control$n_mcmc_moves %||% local_control$n_mcmc_moves %||% 2L)
pre_outer_max_steps <- arg_int(cli_args, "pre_outer_max_steps", calibration_control$max_steps %||% local_control$max_steps %||% 128L)
initial_certification_rounds <- arg_int(
  cli_args,
  "initial_certification_rounds",
  calibration_control$initial_certification_rounds %||% if (pre_outer_rounds > 0L) 2L else 1L
)
initial_certification_max_updates <- arg_int(
  cli_args,
  "initial_certification_max_updates",
  calibration_control$initial_certification_max_updates %||% pre_outer_max_updates
)
stop_on_initial_uncertified <- arg_lgl(cli_args, "stop_on_initial_uncertified", TRUE)
compress_local_particles <- arg_lgl(cli_args, "compress_local_particles", FALSE)
compression_particles <- arg_int(cli_args, "compression_particles", 256L)
compression_theta_points <- arg_int(cli_args, "compression_theta_points", 384L)
compression_holdout_fraction <- arg_num(cli_args, "compression_holdout_fraction", 0.33)
compression_max_holdout_rmse <- arg_num(cli_args, "compression_max_holdout_rmse", 0.05)
compression_stop_on_failure <- arg_lgl(cli_args, "compression_stop_on_failure", FALSE)
compression_require_all <- arg_lgl(cli_args, "compression_require_all", FALSE)
compression_audit_n <- arg_int(cli_args, "compression_audit_n", min(N, 384L))
compression_max_total_centered_rmse <- arg_num(cli_args, "compression_max_total_centered_rmse", Inf)

population_model <- checkpoint$population_model %||% factor_set$population_model
initial_proposal <- normalize_theta_proposal(
  checkpoint$initial_proposal,
  population_model = population_model
)
emc_draws <- NULL
baseline_draws <- NULL
data_list <- NULL
loglik_emc2 <- NULL
need_emc_data <- isTRUE(compare_emc) ||
  pre_outer_rounds > 0L ||
  initial_certification_rounds > 0L
if (isTRUE(need_emc_data)) {
  if (!requireNamespace("EMC2", quietly = TRUE)) {
    stop("EMC2 package is required for pre-outer repair or EMC comparison.")
  }
  if (!file.exists(data_file)) {
    stop("Missing EMC2 benchmark data: ", data_file)
  }
  load(data_file)
  if (!exists("ELP_DDM", inherits = FALSE)) {
    stop("The EMC2 data file must define ELP_DDM.")
  }
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
    if (file.exists(baseline_file)) {
      baseline <- readRDS(baseline_file)
      baseline_draws <- baseline$workflow_draws
    }
  }
}

calibration_history <- checkpoint$calibration_history
if (is.data.frame(calibration_history)) {
  calibration_history <- list(calibration_history)
}
if (!is.list(calibration_history)) {
  calibration_history <- list()
}
initial_certification <- NULL
initial_uncertified_error <- NULL
if (pre_outer_rounds > 0L) {
  calibration_control$pre_outer_rounds <- pre_outer_rounds
  calibration_control$pre_outer_audit_n <- pre_outer_audit_n
  calibration_control$pre_outer_max_points <- pre_outer_max_points
  calibration_control$pre_outer_max_updates <- pre_outer_max_updates
  calibration_control$M <- pre_outer_particles
  calibration_control$n_mcmc_moves <- pre_outer_mcmc_moves
  calibration_control$max_steps <- pre_outer_max_steps
  calibration_control$confirmation_M <- calibration_control$confirmation_M %||% pre_outer_particles
  pre_outer <- local_atlas_pre_outer_certify(
    factor_set = factor_set,
    data_list = data_list,
    loglik_fn = loglik_emc2,
    initial_proposal = initial_proposal,
    local_control = local_control,
    edge_control = edge_control,
    calibration_control = calibration_control,
    outer_control = modifyList(outer_control, list(N = N)),
    design_control = design_control,
    calibration_history = calibration_history,
    n_cores = cores,
    seed = seed,
    outer_seed = seed,
    verbose = verbose,
    trace_verbose = verbose
  )
  factor_set <- pre_outer$factor_set
  calibration_history <- pre_outer$calibration_history
}

if (initial_certification_rounds > 0L) {
  calibration_control$M <- pre_outer_particles
  calibration_control$n_mcmc_moves <- pre_outer_mcmc_moves
  calibration_control$max_steps <- pre_outer_max_steps
  calibration_control$confirmation_M <- calibration_control$confirmation_M %||% pre_outer_particles
  initial_theta <- theta_proposal_sample(
    initial_proposal,
    n = as.integer(N),
    seed = as.integer(seed)
  )
  initial_certification <- local_atlas_certify_theta_cloud(
    factor_set = factor_set,
    theta = initial_theta,
    data_list = data_list,
    loglik_fn = loglik_emc2,
    local_control = local_control,
    edge_control = edge_control,
    calibration_control = calibration_control,
    theta_weights = rep(1 / as.integer(N), as.integer(N)),
    max_rounds = as.integer(initial_certification_rounds),
    max_updates = as.integer(initial_certification_max_updates),
    n_cores = cores,
    seed = seed + 1900003L,
    verbose = verbose,
    trace_verbose = verbose
  )
  factor_set <- initial_certification$factor_set
  if (length(initial_certification$calibration_history)) {
    calibration_history <- c(calibration_history %||% list(), initial_certification$calibration_history)
  }
  if (!isTRUE(initial_certification$certified) && isTRUE(stop_on_initial_uncertified)) {
    shown <- utils::head(initial_certification$uncertified, 8L)
    initial_uncertified_error <- paste0(
      "Initial outer theta cloud remains uncertified after exact-cloud repair:\n",
      paste(
        sprintf(
          "local=%s theta_row=%s status=%s reason=%s nearest=%s",
          shown$local,
          shown$theta_row,
          shown$status,
          shown$reason,
          shown$nearest_charts
        ),
        collapse = "\n"
      )
    )
  }
}

if (pre_outer_rounds > 0L || initial_certification_rounds > 0L) {
  post_pre_outer_checkpoint <- checkpoint
  post_pre_outer_checkpoint$stage <- "post_pre_outer"
  post_pre_outer_checkpoint$factor_set <- factor_set
  post_pre_outer_checkpoint$atlases <- factor_set$atlases
  post_pre_outer_checkpoint$initial_proposal <- initial_proposal
  post_pre_outer_checkpoint$calibration_history <- calibration_history
  post_pre_outer_checkpoint$initial_certification <- initial_certification
  post_pre_outer_checkpoint$checkpoint_time <- Sys.time()
  dir.create(dirname(post_pre_outer_checkpoint_file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(post_pre_outer_checkpoint, post_pre_outer_checkpoint_file)
}
if (!is.null(initial_uncertified_error)) {
  stop(initial_uncertified_error)
}

compression_summary <- data.frame()
compression_audit <- data.frame()
compression_audit_summary <- data.frame()
if (isTRUE(compress_local_particles)) {
  compression_theta <- compression_theta_from_checkpoint(
    checkpoint = checkpoint,
    factor_set = factor_set,
    initial_proposal = initial_proposal,
    n_theta = as.integer(compression_theta_points),
    seed = seed + 7100003L
  )
  compressed_factor_set <- compress_local_atlas_factor_set(
    factor_set = factor_set,
    theta = compression_theta,
    K = as.integer(compression_particles),
    holdout_fraction = compression_holdout_fraction,
    max_holdout_rmse = compression_max_holdout_rmse,
    stop_on_failure = compression_stop_on_failure,
    require_compressed_particle_mis = compression_require_all,
    n_cores = cores,
    seed = seed + 7200003L,
    verbose = TRUE
  )
  compression_summary <- local_atlas_compression_summary(compressed_factor_set)
  dir.create(dirname(compression_summary_csv), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(compression_summary, compression_summary_csv, row.names = FALSE)

  audit_theta <- checkpoint$fit$theta %||% compression_theta
  audit_theta <- .as_hyper_matrix(audit_theta, population_model$hyper_names, population_model$hyper_dim)
  if (nrow(audit_theta) > as.integer(compression_audit_n)) {
    set.seed(as.integer(seed) + 7300003L)
    audit_theta <- audit_theta[sample.int(nrow(audit_theta), as.integer(compression_audit_n)), , drop = FALSE]
  }
  raw_factor_set <- factor_set
  raw_factor_set$evaluator_control$use_compressed_particle_mis <- FALSE
  raw_factor_set$evaluator_control$require_compressed_particle_mis <- FALSE
  raw_factor_set <- validate_local_atlas_factor_set(raw_factor_set)
  compressed_loglik <- population_factor_set_loglik(
    compressed_factor_set,
    audit_theta,
    include_constant = FALSE,
    n_cores = cores
  )
  raw_loglik <- population_factor_set_loglik(
    raw_factor_set,
    audit_theta,
    include_constant = FALSE,
    n_cores = cores
  )
  compression_error <- compressed_loglik - raw_loglik
  compression_audit <- data.frame(
    theta_row = seq_len(nrow(audit_theta)),
    raw_loglik = raw_loglik,
    compressed_loglik = compressed_loglik,
    compression_error = compression_error,
    centered_compression_error = compression_error - mean(compression_error[is.finite(compression_error)], na.rm = TRUE),
    check.names = FALSE
  )
  compression_audit_summary <- data.frame(
    n_theta = nrow(compression_audit),
    certified_locals = sum(compression_summary$certified %in% TRUE),
    fallback_locals = sum(!(compression_summary$certified %in% TRUE)),
    median_selected_particles = stats::median(compression_summary$selected_particles, na.rm = TRUE),
    median_raw_particles = stats::median(compression_summary$raw_particles, na.rm = TRUE),
    median_holdout_rmse = stats::median(compression_summary$holdout_rmse, na.rm = TRUE),
    max_holdout_rmse = max(compression_summary$holdout_rmse, na.rm = TRUE),
    total_error_rmse = rmse_finite(compression_audit$compression_error),
    total_centered_error_rmse = rmse_finite(compression_audit$centered_compression_error),
    total_max_abs_centered_error = max(abs(compression_audit$centered_compression_error), na.rm = TRUE),
    check.names = FALSE
  )
  dir.create(dirname(compression_audit_csv), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(compression_audit, compression_audit_csv, row.names = FALSE)
  print(compression_audit_summary, row.names = FALSE)
  if (is.finite(compression_max_total_centered_rmse) &&
      compression_audit_summary$total_centered_error_rmse > compression_max_total_centered_rmse) {
    stop(
      "Compressed factor set failed total-error audit: centered RMSE=",
      signif(compression_audit_summary$total_centered_error_rmse, 4),
      " threshold=", signif(compression_max_total_centered_rmse, 4)
    )
  }
  factor_set <- compressed_factor_set
}

start_time <- Sys.time()
fit <- outer_population_smc(
  factor_set = factor_set,
  N = N,
  initial_proposal = initial_proposal,
  resample_threshold = resample_threshold,
  n_mcmc_moves = n_mcmc_moves,
  min_mcmc_moves = min_mcmc_moves,
  max_rounds = max_rounds,
  rw_scale_init = rw_scale_init,
  n_cores = cores,
  seed = seed,
  verbose = verbose
)
elapsed_sec <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

workflow <- list(
  fit = fit,
  factor_set = factor_set,
  population_model = population_model,
  initial_proposal = initial_proposal,
  theta_root = checkpoint$theta_root,
  theta_design = checkpoint$theta_design,
  theta_cloud = checkpoint$theta_cloud,
  chart_design_plan = checkpoint$chart_design_plan,
  design_certification = checkpoint$design_certification,
  proposal_certification = checkpoint$proposal_certification,
  atlas_build_history = checkpoint$atlas_build_history,
  calibration_history = calibration_history,
  initial_certification = initial_certification,
  graph_summary = local_atlas_graph_summary(factor_set),
  compression_summary = compression_summary,
  compression_audit_summary = compression_audit_summary,
  settings = settings
)

workflow_draws <- NULL
posterior_comparison <- NULL
posterior_summary <- NULL
if (isTRUE(compare_emc)) {
  draw_n <- if (is.finite(n_draws)) as.integer(n_draws) else nrow(emc_draws)
  workflow_draws <- local_atlas_draws_from_fit(
    fit,
    population_model = population_model,
    n_draws = draw_n,
    seed = seed + 1L
  )
  posterior_comparison <- local_atlas_compare_posterior_draws(emc_draws, workflow_draws)
  posterior_summary <- .local_atlas_metric_summary(posterior_comparison)
  dir.create(dirname(comparison_csv), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(posterior_comparison, comparison_csv, row.names = FALSE)
  dir.create(dirname(plot_file), recursive = TRUE, showWarnings = FALSE)
  grDevices::png(plot_file, width = 1800, height = 1400)
  plot_posteriors(
    emc_draws,
    workflow_draws,
    labels = c("EMC2", "explicit_outer"),
    cols = c("black", "firebrick3"),
    n_cols = 4L
  )
  grDevices::dev.off()
}

dir.create(dirname(results_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(post_outer_checkpoint_file), recursive = TRUE, showWarnings = FALSE)
saveRDS(
  list(
    source_checkpoint = checkpoint_file,
    checkpoint_stage = checkpoint_stage,
    workflow = workflow,
    fit = fit,
    factor_set = factor_set,
    population_model = population_model,
    initial_proposal = initial_proposal,
    design_certification = checkpoint$design_certification,
    calibration_history = calibration_history,
    initial_certification = initial_certification,
    compression_summary = compression_summary,
    compression_audit = compression_audit,
    compression_audit_summary = compression_audit_summary,
    emc_draws = emc_draws,
    workflow_draws = workflow_draws,
    baseline_draws = baseline_draws,
    posterior_comparison = posterior_comparison,
    posterior_summary = posterior_summary,
    plot_file = if (isTRUE(compare_emc)) plot_file else NULL,
    comparison_csv = if (isTRUE(compare_emc)) comparison_csv else NULL,
    settings = list(
      outer_particles = N,
      outer_mcmc_moves = n_mcmc_moves,
      outer_min_mcmc_moves = min_mcmc_moves,
      outer_max_rounds = max_rounds,
      outer_resample_threshold = resample_threshold,
      outer_rw_scale_init = rw_scale_init,
      cores = cores,
      seed = seed,
      elapsed_sec = elapsed_sec,
      compare_emc = compare_emc,
      data_file = if (isTRUE(compare_emc)) data_file else NULL,
      baseline_file = if (isTRUE(compare_emc)) baseline_file else NULL,
      post_pre_outer_checkpoint_file = if (pre_outer_rounds > 0L) post_pre_outer_checkpoint_file else NULL,
      post_outer_checkpoint_file = post_outer_checkpoint_file,
      pre_outer_rounds = pre_outer_rounds,
      pre_outer_audit_n = pre_outer_audit_n,
      pre_outer_max_points = pre_outer_max_points,
      pre_outer_max_updates = pre_outer_max_updates,
      pre_outer_particles = pre_outer_particles,
      pre_outer_mcmc_moves = pre_outer_mcmc_moves,
      pre_outer_max_steps = pre_outer_max_steps,
      initial_certification_rounds = initial_certification_rounds,
      initial_certification_max_updates = initial_certification_max_updates,
      stop_on_initial_uncertified = stop_on_initial_uncertified,
      compress_local_particles = compress_local_particles,
      compression_particles = compression_particles,
      compression_theta_points = compression_theta_points,
      compression_holdout_fraction = compression_holdout_fraction,
      compression_max_holdout_rmse = compression_max_holdout_rmse,
      compression_stop_on_failure = compression_stop_on_failure,
      compression_require_all = compression_require_all,
      compression_audit_n = compression_audit_n,
      compression_summary_csv = if (isTRUE(compress_local_particles)) compression_summary_csv else NULL,
      compression_audit_csv = if (isTRUE(compress_local_particles)) compression_audit_csv else NULL,
      evaluator_control_override = override_control
    )
  ),
  results_file
)
post_outer_checkpoint <- checkpoint
post_outer_checkpoint$stage <- "post_outer"
post_outer_checkpoint$fit <- fit
post_outer_checkpoint$factor_set <- factor_set
post_outer_checkpoint$atlases <- factor_set$atlases
post_outer_checkpoint$population_model <- population_model
post_outer_checkpoint$initial_proposal <- initial_proposal
post_outer_checkpoint$calibration_history <- calibration_history
post_outer_checkpoint$graph_summary <- workflow$graph_summary
post_outer_checkpoint$compression_summary <- compression_summary
post_outer_checkpoint$compression_audit_summary <- compression_audit_summary
post_outer_checkpoint$outer_rerun_source_checkpoint <- checkpoint_file
post_outer_checkpoint$outer_rerun_results_file <- results_file
post_outer_checkpoint$checkpoint_time <- Sys.time()
saveRDS(post_outer_checkpoint, post_outer_checkpoint_file)

cat(sprintf("Saved explicit outer rerun: %s\n", results_file))
if (pre_outer_rounds > 0L || initial_certification_rounds > 0L) {
  cat(sprintf("Saved post-pre-outer checkpoint: %s\n", post_pre_outer_checkpoint_file))
}
cat(sprintf("Saved post-outer checkpoint: %s\n", post_outer_checkpoint_file))
cat(sprintf("Source checkpoint: %s stage=%s\n", checkpoint_file, checkpoint_stage))
cat(sprintf("Elapsed seconds: %.1f\n", elapsed_sec))
cat(sprintf("Log evidence: %.6f\n", as.numeric(fit$log_evidence %||% NA_real_)))
if (isTRUE(compress_local_particles)) {
  cat(sprintf("Saved compression summary: %s\n", compression_summary_csv))
  cat(sprintf("Saved compression audit: %s\n", compression_audit_csv))
}
if (isTRUE(compare_emc)) {
  cat(sprintf("Saved posterior comparison: %s\n", comparison_csv))
  cat(sprintf("Saved posterior plot: %s\n", plot_file))
  cat("\nPosterior summary:\n")
  print(posterior_summary, row.names = FALSE)
}
