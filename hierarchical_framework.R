if (!exists("%||%", mode = "function")) {
  source("utilities.R")
}
if (!exists("make_population_model_diag_gaussian", mode = "function")) {
  source("population_models.R")
}
if (!exists("outer_population_smc", mode = "function")) {
  source("outer_population_smc.R")
}
if (!exists("build_local_atlas", mode = "function")) {
  source("local_charts.R")
}

.hierarchy_merge_control <- function(control, defaults) {
  if (is.null(control)) control <- list()
  if (!is.list(control)) stop("control must be a list.")
  utils::modifyList(defaults, control)
}

.hierarchy_as_model <- function(model) {
  if (!inherits(model, "hierarchical_framework_model")) {
    stop("model must be created by define_hierarchical_model().")
  }
  model
}

.hierarchy_call_stage <- function(callback, stage, object) {
  if (is.function(callback)) callback(stage = stage, object = object)
  invisible(object)
}

define_hierarchical_model <- function(data_list,
                                      loglik_fn,
                                      alpha_names = NULL,
                                      population_model = NULL,
                                      local_names = names(data_list),
                                      label = "hierarchical_model") {
  if (!is.list(data_list) || !length(data_list)) {
    stop("data_list must be a non-empty list.")
  }
  if (!is.function(loglik_fn)) {
    stop("loglik_fn must be a function accepting (Theta, data_i).")
  }
  if (is.null(alpha_names)) {
    if (!is.null(population_model$alpha_names)) {
      alpha_names <- population_model$alpha_names
    } else {
      stop("alpha_names are required when population_model does not provide them.")
    }
  }
  alpha_names <- as.character(alpha_names)
  if (!length(alpha_names) || any(!nzchar(alpha_names))) {
    stop("alpha_names must be non-empty names.")
  }
  if (is.null(population_model)) {
    d <- length(alpha_names)
    population_model <- make_population_model_diag_gaussian(
      alpha_names = alpha_names,
      mean_prior_mean = rep(0, d),
      mean_prior_var = rep(1, d),
      sigma2_prior_shape = rep(2, d),
      sigma2_prior_rate = rep(0.3, d),
      label = label
    )
  }
  population_model <- normalize_population_model(population_model)
  if (!identical(as.character(population_model$alpha_names), alpha_names)) {
    stop("population_model alpha_names do not match model alpha_names.")
  }
  if (is.null(local_names)) local_names <- names(data_list)
  if (is.null(local_names)) local_names <- as.character(seq_along(data_list))
  local_names <- as.character(local_names)
  bad <- !nzchar(local_names) | is.na(local_names)
  local_names[bad] <- as.character(which(bad))
  names(data_list) <- local_names
  structure(
    list(
      data_list = data_list,
      loglik_fn = loglik_fn,
      alpha_names = alpha_names,
      population_model = population_model,
      local_names = local_names,
      label = as.character(label)
    ),
    class = "hierarchical_framework_model"
  )
}

build_theta_support <- function(model = NULL,
                                theta = NULL,
                                weights = NULL,
                                fit = NULL,
                                population_model = NULL,
                                max_theta = NULL,
                                tail_points_per_axis = 2L,
                                seed = 123L,
                                source = "theta_support") {
  if (inherits(theta, "local_evidence_certification_cloud")) {
    return(validate_local_evidence_certification_cloud(theta))
  }
  if (!is.null(model)) {
    model <- .hierarchy_as_model(model)
    population_model <- model$population_model
  } else {
    population_model <- normalize_population_model(population_model)
  }
  if (!is.null(fit)) {
    theta <- fit$theta
    weights <- fit$w
    source <- source %||% "outer_posterior"
  }
  if (is.null(theta)) {
    stop("theta or fit is required to build theta support.")
  }
  theta <- .as_hyper_matrix(theta, population_model$hyper_names, population_model$hyper_dim)
  weights <- .local_chart_normalize_weights(weights, nrow(theta))
  if (is.null(max_theta)) max_theta <- nrow(theta)
  max_theta <- min(as.integer(max_theta), nrow(theta))
  if (!is.finite(max_theta) || max_theta < 1L) {
    stop("max_theta must select at least one theta row.")
  }

  set.seed(as.integer(seed))
  top_n <- min(nrow(theta), max(4L, floor(max_theta * 0.25)))
  rows <- head(order(weights, decreasing = TRUE), top_n)

  if (as.integer(tail_points_per_axis) > 0L) {
    mu_idx <- grep("^mu_", colnames(theta))
    log_sigma_idx <- grep("^log_sigma2_", colnames(theta))
    focus_idx <- unique(c(mu_idx, log_sigma_idx))
    if (length(focus_idx)) {
      z_info <- .local_evidence_weighted_center_cov(theta, weights, ridge = 1e-8)
      z <- sweep(theta, 2L, z_info$center, "-") %*% z_info$whitening
      axis_risk <- colSums(sweep(abs(z), 1L, weights, "*"), na.rm = TRUE)
      axis_order <- focus_idx[order(axis_risk[focus_idx], decreasing = TRUE)]
      axis_order <- head(axis_order, max(1L, floor(max_theta * 0.15)))
      probs <- if (as.integer(tail_points_per_axis) <= 1L) {
        0.5
      } else {
        seq(0.05, 0.95, length.out = as.integer(tail_points_per_axis))
      }
      for (j in axis_order) {
        qj <- as.numeric(.local_atlas_weighted_quantile(theta[, j], weights, probs = probs))
        for (target in qj) rows <- c(rows, which.min(abs(theta[, j] - target)))
      }
    }
  }

  rows <- unique(rows)
  need <- max_theta - length(rows)
  if (need > 0L) {
    remaining <- setdiff(seq_len(nrow(theta)), rows)
    if (length(remaining)) {
      sampled <- sample(
        remaining,
        size = min(need, length(remaining)),
        replace = FALSE,
        prob = weights[remaining]
      )
      rows <- c(rows, sampled)
    }
  }
  if (length(rows) > max_theta) {
    rows <- rows[head(order(weights[rows], decreasing = TRUE), max_theta)]
  }
  rows <- rows[order(rows)]
  build_local_evidence_certification_cloud(
    theta = theta[rows, , drop = FALSE],
    population_model = population_model,
    theta_weights = .local_chart_normalize_weights(weights[rows], length(rows)),
    theta_source = source
  )
}

build_local_evidence_atlases <- function(model,
                                         theta_support,
                                         theta_design = NULL,
                                         theta_root = NULL,
                                         initial_proposal = NULL,
                                         local_control = list(),
                                         design_control = list(),
                                         edge_control = list(),
                                         evaluator_control = list(),
                                         proposal_control = list(),
                                         n_cores = 1L,
                                         seed = 123L,
                                         verbose = TRUE,
                                         trace_verbose = FALSE,
                                         checkpoint_file = NULL,
                                         resume_checkpoint = FALSE) {
  model <- .hierarchy_as_model(model)
  theta_support <- build_theta_support(model, theta = theta_support)
  theta_design <- theta_design %||% theta_support$theta
  design_control <- .hierarchy_merge_control(design_control, list(stop_after_atlas_build = TRUE))
  fit_chart_atlas_population_model(
    data_list = model$data_list,
    loglik_fn = model$loglik_fn,
    population_model = model$population_model,
    theta_design = theta_design,
    theta_root = theta_root,
    theta_cloud = theta_support$theta,
    initial_proposal = initial_proposal,
    local_control = local_control,
    design_control = design_control,
    edge_control = edge_control,
    evaluator_control = evaluator_control,
    proposal_control = proposal_control,
    outer_control = list(N = 1L),
    n_cores = n_cores,
    seed = seed,
    verbose = verbose,
    trace_verbose = trace_verbose,
    checkpoint_file = checkpoint_file,
    resume_checkpoint = resume_checkpoint
  )
}

certify_local_evidence <- function(factor_set,
                                   theta_support,
                                   include_theta = TRUE,
                                   n_cores = 1L) {
  factor_set <- validate_local_atlas_factor_set(factor_set)
  theta_support <- build_theta_support(
    population_model = factor_set$population_model,
    theta = theta_support
  )
  table <- evaluate_raw_local_evidence_certification(
    factor_set = factor_set,
    cloud = theta_support,
    include_theta = include_theta,
    n_cores = n_cores
  )
  structure(
    list(
      cloud = theta_support,
      table = table,
      summary = summarize_raw_local_evidence_certification(table)
    ),
    class = "hierarchical_local_evidence_certification"
  )
}

compress_local_evidence <- function(factor_set,
                                    theta_support,
                                    control = list(),
                                    n_cores = 1L,
                                    seed = 123L,
                                    verbose = TRUE) {
  factor_set <- validate_local_atlas_factor_set(factor_set)
  theta_support <- build_theta_support(
    population_model = factor_set$population_model,
    theta = theta_support
  )
  control <- .hierarchy_merge_control(control, list(
    K = 64L,
    holdout_fraction = 0.25,
    ridge = 1e-8,
    evidence_weight = 1,
    moment_weight = 0.05,
    chart_weight = 0.05,
    include_moments = TRUE,
    include_chart = TRUE,
    max_holdout_rmse = Inf,
    stop_on_failure = FALSE,
    require_compressed_particle_mis = FALSE
  ))
  compress_local_atlas_factor_set(
    factor_set = factor_set,
    theta = theta_support$theta,
    K = as.integer(control$K),
    holdout_fraction = control$holdout_fraction,
    ridge = control$ridge,
    evidence_weight = control$evidence_weight,
    moment_weight = control$moment_weight,
    chart_weight = control$chart_weight,
    include_moments = control$include_moments,
    include_chart = control$include_chart,
    max_holdout_rmse = control$max_holdout_rmse,
    stop_on_failure = control$stop_on_failure,
    require_compressed_particle_mis = control$require_compressed_particle_mis,
    n_cores = n_cores,
    seed = seed,
    verbose = verbose
  )
}

run_outer_hierarchy <- function(factor_set,
                                initial_proposal = NULL,
                                control = list(),
                                n_cores = 1L,
                                seed = 123L,
                                verbose = TRUE) {
  factor_set <- validate_local_atlas_factor_set(factor_set)
  control <- .hierarchy_merge_control(control, list(
    N = 1000L,
    n_mcmc_moves = 3L,
    min_mcmc_moves = 1L,
    max_rounds = 80L,
    resample_threshold = 0.5,
    rw_scale_init = 0.8
  ))
  outer_population_smc(
    factor_set = factor_set,
    N = as.integer(control$N),
    initial_proposal = initial_proposal,
    resample_threshold = control$resample_threshold,
    n_mcmc_moves = as.integer(control$n_mcmc_moves),
    min_mcmc_moves = as.integer(control$min_mcmc_moves),
    max_rounds = as.integer(control$max_rounds),
    rw_scale_init = control$rw_scale_init,
    n_cores = n_cores,
    seed = seed,
    verbose = verbose
  )
}

certify_posterior_region <- function(factor_set,
                                     fit,
                                     max_theta = 96L,
                                     tail_points_per_axis = 2L,
                                     seed = 123L,
                                     n_cores = 1L) {
  factor_set <- validate_local_atlas_factor_set(factor_set)
  theta_support <- build_theta_support(
    population_model = factor_set$population_model,
    fit = fit,
    max_theta = max_theta,
    tail_points_per_axis = tail_points_per_axis,
    seed = seed,
    source = "outer_posterior"
  )
  certify_local_evidence(
    factor_set = factor_set,
    theta_support = theta_support,
    include_theta = TRUE,
    n_cores = n_cores
  )
}

select_local_evidence_probes <- function(model,
                                         factor_set,
                                         certification,
                                         selector = c("risk", "active"),
                                         residual_sources = list(),
                                         active_round_ids = integer(),
                                         active_exclude_probe_sets = "holdout",
                                         graph_summary = NULL,
                                         compression_summary = NULL,
                                         control = list(),
                                         n_cores = 1L) {
  model <- .hierarchy_as_model(model)
  factor_set <- validate_local_atlas_factor_set(factor_set)
  selector <- match.arg(tolower(selector), c("risk", "active"))
  raw_table <- if (inherits(certification, "hierarchical_local_evidence_certification")) {
    certification$table
  } else {
    certification
  }
  control <- .hierarchy_merge_control(control, list(
    max_repair_pairs = 12L,
    max_holdout_pairs = 8L,
    max_repair_pairs_per_local = 2L,
    max_repair_pairs_per_theta = 3L,
    max_holdout_pairs_per_local = 2L,
    max_holdout_pairs_per_theta = 3L,
    lambda_shape = 1,
    lambda_evidence = 0.10,
    lambda_mean_shape = 0.10,
    lambda_support = 0.10,
    lambda_raw_weakness = 0.25,
    lambda_graph = 0.20,
    lambda_novelty = 0.25,
    exploration_metric_floor = 0.05,
    exploration_residual_floor = 0.25,
    exploration_repair_fraction = 0.35,
    exploration_holdout_fraction = 0.35
  ))

  residual_pool <- NULL
  feature_model <- NULL
  acquisition_table <- NULL
  if (identical(selector, "active")) {
    if (!length(residual_sources)) {
      stop("active probe selection requires residual_sources.")
    }
    residual_sources <- lapply(residual_sources, function(x) {
      if (is.character(x) && length(x) == 1L) readRDS(x) else x
    })
    source_labels <- names(residual_sources)
    if (is.null(source_labels) || any(!nzchar(source_labels))) {
      source_labels <- paste0("residual_source_", seq_along(residual_sources))
    }
    residual_pool <- build_shape_residual_pool(residual_sources, labels = source_labels)
    feature_model <- fit_shape_residual_feature_model(
      residual_pool,
      round_ids = if (length(active_round_ids)) active_round_ids else NULL,
      exclude_probe_sets = active_exclude_probe_sets
    )
    selection <- select_shape_active_probe_pairs(
      feature_model = feature_model,
      factor_set = factor_set,
      raw_certification_table = raw_table,
      graph_summary = graph_summary %||% local_atlas_graph_summary(factor_set),
      compression_summary = compression_summary %||% factor_set$compression_summary,
      control = control,
      n_cores = n_cores
    )
    acquisition_table <- selection$scored_table
  } else {
    risk_control <- .hierarchy_merge_control(control, list(
      certified_repair_fraction = 0.45,
      uncertified_repair_fraction = 0.30,
      tail_repair_fraction = 0.25,
      certified_holdout_fraction = 0.60,
      tail_holdout_fraction = 0.25
    ))
    selection <- select_shape_probe_pairs(
      factor_set = factor_set,
      raw_certification_table = raw_table,
      graph_summary = graph_summary %||% local_atlas_graph_summary(factor_set),
      compression_summary = compression_summary %||% factor_set$compression_summary,
      control = risk_control
    )
  }
  structure(
    list(
      selection = selection,
      selector = selector,
      residual_pool = residual_pool,
      feature_model = feature_model,
      acquisition_table = acquisition_table
    ),
    class = "hierarchical_probe_selection"
  )
}

run_local_evidence_probes <- function(model,
                                      factor_set,
                                      selection,
                                      control = list(),
                                      local_control = list(),
                                      n_cores = 1L,
                                      seed = 123L,
                                      verbose = FALSE) {
  model <- .hierarchy_as_model(model)
  factor_set <- validate_local_atlas_factor_set(factor_set)
  selection_object <- if (inherits(selection, "hierarchical_probe_selection")) {
    selection$selection
  } else {
    selection
  }
  control <- .hierarchy_merge_control(control, list(
    M = 192L,
    n_replicates = 1L,
    bootstrap_B = 100L,
    geometry_control = list(
      max_directions = 3L,
      max_local_stencils = 24L,
      stencil_step = 0.65,
      local_stencil_step = 0.50,
      min_abs_standardized_residual = 1.0,
      min_abs_residual = 0.05,
      max_stencils_per_local = 2L
    )
  ))
  run_shape_probe_pairs(
    factor_set = factor_set,
    selection = selection_object,
    data_list = model$data_list,
    loglik_fn = model$loglik_fn,
    M = as.integer(control$M),
    n_replicates = as.integer(control$n_replicates),
    local_control = local_control,
    bootstrap_B = as.integer(control$bootstrap_B),
    n_cores = n_cores,
    seed = seed,
    verbose = verbose,
    geometry_control = control$geometry_control
  )
}

repair_local_evidence <- function(model,
                                  factor_set,
                                  probes,
                                  executor = c("geometry", "selected_exact"),
                                  control = list(),
                                  local_control = list(),
                                  edge_control = list(),
                                  n_cores = 1L,
                                  seed = 123L,
                                  verbose = FALSE) {
  model <- .hierarchy_as_model(model)
  factor_set <- validate_local_atlas_factor_set(factor_set)
  executor <- match.arg(tolower(executor), c("selected_exact", "geometry"))
  control <- .hierarchy_merge_control(control, list(
    max_repairs = 12L,
    M = 192L,
    n_mcmc_moves = 2L,
    target_cess = 0.90,
    max_steps = 128L,
    max_updates = 12L,
    direct_confirmation_reps = 1L,
    direct_confirmation_M = 192L,
    audit_pre_repair = TRUE,
    audit_post_repair = TRUE,
    audit_scope = "repaired",
    stop_on_empty = FALSE,
    min_abs_standardized_residual = 1.0,
    min_abs_residual = 0.05
  ))
  edge_control <- .hierarchy_merge_control(edge_control, list(
    edge_neighbors = 2L,
    max_intermediates = 4L,
    min_overlap_ess = 0.03,
    max_se = 1.25,
    max_forward_reverse_gap = 1.25,
    max_taylor_gap = 3.0,
    require_bar_converged = TRUE
  ))
  repair_fn <- if (identical(executor, "selected_exact")) {
    repair_shape_selected_probe_pairs
  } else {
    repair_shape_residual_geometry
  }
  repair_fn(
    factor_set = factor_set,
    shape_probe = probes,
    data_list = model$data_list,
    loglik_fn = model$loglik_fn,
    control = control,
    local_control = local_control,
    edge_control = edge_control,
    n_cores = n_cores,
    seed = seed,
    verbose = verbose
  )
}

validate_local_evidence_repair <- function(model,
                                           repair,
                                           probes,
                                           control = list(),
                                           local_control = list(),
                                           n_cores = 1L,
                                           seed = 123L,
                                           verbose = FALSE) {
  model <- .hierarchy_as_model(model)
  control <- .hierarchy_merge_control(control, list(
    M = 192L,
    n_replicates = 1L,
    bootstrap_B = 100L,
    require_holdout_improvement = FALSE,
    max_pair_centered_rmse_ratio = 1.10,
    max_local_centered_rmse_ratio = 1.10,
    max_total_centered_rmse_ratio = 1.10,
    max_abs_total_increase = 0.25,
    max_graph_edge_z_increase = 1.0
  ))
  validate_shape_repair_holdout(
    shape_repair = repair,
    shape_probe = probes,
    data_list = model$data_list,
    loglik_fn = model$loglik_fn,
    control = control,
    local_control = local_control,
    n_cores = n_cores,
    seed = seed,
    verbose = verbose
  )
}

gate_repair_and_update_outer <- function(repair,
                                         old_fit,
                                         old_factor_set,
                                         holdout_validation,
                                         initial_proposal = NULL,
                                         control = list(),
                                         n_cores = 1L,
                                         seed = 123L,
                                         verbose = FALSE) {
  old_factor_set <- validate_local_atlas_factor_set(old_factor_set)
  control <- .hierarchy_merge_control(control, list(
    min_reweight_ess_fraction = 0.50,
    low_ess_rerun_fraction = 0.25,
    max_psis_k = 0.70,
    n_draws = nrow(old_fit$theta),
    require_holdout_acceptance = TRUE,
    run_outer_if_required = FALSE,
    force_outer_rerun = FALSE,
    outer_particles = nrow(old_fit$theta),
    outer_mcmc_moves = 3L,
    outer_max_rounds = 90L
  ))
  gate <- shape_repair_outer_reweight_gate(
    shape_repair = repair,
    fit = old_fit,
    old_factor_set = old_factor_set,
    holdout_validation = holdout_validation,
    reference_draws = NULL,
    baseline_draws = NULL,
    control = list(
      min_reweight_ess_fraction = control$min_reweight_ess_fraction,
      low_ess_rerun_fraction = control$low_ess_rerun_fraction,
      max_psis_k = control$max_psis_k,
      n_draws = as.integer(control$n_draws),
      require_holdout_acceptance = control$require_holdout_acceptance
    ),
    n_cores = n_cores,
    seed = seed
  )

  final_fit <- old_fit
  final_factor_set <- old_factor_set
  final_source <- "original_post_outer"
  outer_rerun <- NULL
  if (isTRUE(control$force_outer_rerun) ||
      (isTRUE(control$run_outer_if_required) &&
        isTRUE(gate$accepted) &&
        identical(gate$decision, "rerun_outer_required_low_reweight_quality"))) {
    outer_rerun <- run_outer_hierarchy(
      factor_set = repair$factor_set,
      initial_proposal = initial_proposal,
      control = list(
        N = as.integer(control$outer_particles),
        n_mcmc_moves = as.integer(control$outer_mcmc_moves),
        max_rounds = as.integer(control$outer_max_rounds)
      ),
      n_cores = n_cores,
      seed = seed + 10L,
      verbose = verbose
    )
    final_fit <- outer_rerun
    final_factor_set <- repair$factor_set
    final_source <- "frozen_outer_rerun_after_shape_repair"
  } else if (isTRUE(gate$accepted) && !is.null(gate$reweighted_fit)) {
    final_fit <- gate$reweighted_fit
    final_factor_set <- repair$factor_set
    final_source <- paste0("reweighted_fit:", gate$decision)
  }

  structure(
    list(
      gate = gate,
      outer_rerun = outer_rerun,
      fit = final_fit,
      factor_set = final_factor_set,
      final_source = final_source
    ),
    class = "hierarchical_outer_update"
  )
}

save_hierarchy_checkpoint <- function(state, file, stage = NULL) {
  if (!is.null(stage)) state$stage <- stage
  state$checkpoint_time <- Sys.time()
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(state, file)
  invisible(normalizePath(file, winslash = "/", mustWork = FALSE))
}

run_hierarchy_shape_calibration <- function(model,
                                            factor_set,
                                            fit,
                                            initial_proposal = NULL,
                                            max_theta = 96L,
                                            tail_points_per_axis = 2L,
                                            selector = c("risk", "active"),
                                            residual_sources = list(),
                                            active_round_ids = integer(),
                                            active_exclude_probe_sets = "holdout",
                                            stop_after_probe = FALSE,
                                            repair_executor = c("geometry", "selected_exact"),
                                            local_control = list(),
                                            probe_control = list(),
                                            selection_control = list(),
                                            repair_control = list(),
                                            repair_edge_control = list(),
                                            holdout_control = list(),
                                            gate_control = list(),
                                            graph_summary = NULL,
                                            compression_summary = NULL,
                                            n_cores = 1L,
                                            seed = 123L,
                                            verbose = FALSE,
                                            stage_callback = NULL) {
  model <- .hierarchy_as_model(model)
  factor_set <- validate_local_atlas_factor_set(factor_set)
  old_factor_set <- factor_set
  selector <- match.arg(tolower(selector), c("risk", "active"))
  repair_executor <- match.arg(tolower(repair_executor), c("selected_exact", "geometry"))

  theta_support <- build_theta_support(
    model = model,
    fit = fit,
    max_theta = max_theta,
    tail_points_per_axis = tail_points_per_axis,
    seed = seed,
    source = "outer_posterior"
  )
  .hierarchy_call_stage(stage_callback, "theta_support", theta_support)

  certification <- certify_local_evidence(
    factor_set = factor_set,
    theta_support = theta_support,
    include_theta = TRUE,
    n_cores = n_cores
  )
  .hierarchy_call_stage(stage_callback, "local_certification", certification)

  probe_selection <- select_local_evidence_probes(
    model = model,
    factor_set = factor_set,
    certification = certification,
    selector = selector,
    residual_sources = residual_sources,
    active_round_ids = active_round_ids,
    active_exclude_probe_sets = active_exclude_probe_sets,
    graph_summary = graph_summary,
    compression_summary = compression_summary,
    control = selection_control,
    n_cores = n_cores
  )
  .hierarchy_call_stage(stage_callback, "probe_selection", probe_selection)

  probes <- run_local_evidence_probes(
    model = model,
    factor_set = factor_set,
    selection = probe_selection,
    control = probe_control,
    local_control = local_control,
    n_cores = n_cores,
    seed = seed + 10L,
    verbose = verbose
  )
  .hierarchy_call_stage(stage_callback, "direct_probes", probes)

  if (isTRUE(stop_after_probe)) {
    return(structure(
      list(
        final_source = "stopped_after_direct_shape_probe",
        theta_support = theta_support,
        certification = certification,
        probe_selection = probe_selection,
        probes = probes,
        fit = fit,
        factor_set = factor_set,
        population_model = model$population_model
      ),
      class = "hierarchical_shape_calibration"
    ))
  }

  repair <- repair_local_evidence(
    model = model,
    factor_set = factor_set,
    probes = probes,
    executor = repair_executor,
    control = repair_control,
    local_control = local_control,
    edge_control = repair_edge_control,
    n_cores = n_cores,
    seed = seed + 20L,
    verbose = verbose
  )
  .hierarchy_call_stage(stage_callback, "repair", repair)

  holdout <- validate_local_evidence_repair(
    model = model,
    repair = repair,
    probes = probes,
    control = holdout_control,
    local_control = local_control,
    n_cores = n_cores,
    seed = seed + 30L,
    verbose = verbose
  )
  .hierarchy_call_stage(stage_callback, "holdout", holdout)

  update <- gate_repair_and_update_outer(
    repair = repair,
    old_fit = fit,
    old_factor_set = old_factor_set,
    holdout_validation = holdout,
    initial_proposal = initial_proposal,
    control = gate_control,
    n_cores = n_cores,
    seed = seed + 40L,
    verbose = verbose
  )
  .hierarchy_call_stage(stage_callback, "outer_update", update)

  structure(
    list(
      final_source = update$final_source,
      theta_support = theta_support,
      certification = certification,
      probe_selection = probe_selection,
      probes = probes,
      repair = repair,
      holdout = holdout,
      outer_update = update,
      fit = update$fit,
      factor_set = update$factor_set,
      population_model = model$population_model
    ),
    class = "hierarchical_shape_calibration"
  )
}
