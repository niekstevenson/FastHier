#!/usr/bin/env Rscript
# ============================================================================
# Certified local evidence chart primitives
# - Phase 1: chart, edge, and atlas data objects
# - Phase 2: diagonal Gaussian score and curvature for log m_i(theta)
# ============================================================================

if (!exists("%||%", mode = "function")) {
  source("smc_core.R")
}
if (!exists("normalize_population_model", mode = "function") ||
    !exists(".as_hyper_matrix", mode = "function")) {
  source("population_models.R")
}
if (!exists("make_reference_prior_gaussian", mode = "function") ||
    !exists("reference_prior_sample", mode = "function")) {
  source("reference_priors.R")
}
if (!exists("run_tempered_smc", mode = "function")) {
  source("SMC_super_fast.R")
}
if (!exists("fit_theta_q0_proposal", mode = "function") ||
    !exists("theta_proposal_sample", mode = "function")) {
  source("theta_proposals.R")
}
if (!exists("outer_population_smc", mode = "function")) {
  source("outer_population_smc.R")
}

.local_chart_status <- c("candidate", "active", "quarantined", "retired")
.local_edge_status <- c("candidate", "active", "rejected")
.local_edge_methods <- c("BAR", "bridge", "SMC-path", "Taylor-crosscheck")

.local_chart_normalizer_certified <- function(chart) {
  certified <- chart$diagnostics$normalizer_certified
  if (is.null(certified)) {
    return(FALSE)
  }
  isTRUE(certified)
}

.local_chart_set_normalizer_certification <- function(chart,
                                                      certified,
                                                      method,
                                                      reason = NULL,
                                                      details = list()) {
  chart$diagnostics$normalizer_certified <- isTRUE(certified)
  chart$diagnostics$normalizer_certification <- c(
    list(
      certified = isTRUE(certified),
      method = as.character(method),
      reason = as.character(reason %||% "")
    ),
    details
  )
  chart
}

.local_chart_observation_scale <- function(x, key, default = 0) {
  if (is.null(x) || !length(x)) {
    return(as.numeric(default))
  }
  if (!is.null(names(x)) && key %in% names(x)) {
    out <- as.numeric(x[[key]])
  } else {
    out <- as.numeric(x[[1L]])
  }
  if (!is.finite(out) || out < 0) {
    stop("observation overdispersion values must be finite and non-negative.")
  }
  out
}

.local_chart_replicated_logz <- function(runs, bootstrap_B = 200L, seed = NULL) {
  if (!length(runs)) {
    stop("replicated logZ certification requires at least one run.")
  }
  logZ <- vapply(runs, `[[`, numeric(1), "logZ")
  logZ_se <- vapply(runs, `[[`, numeric(1), "logZ_se")
  if (any(!is.finite(logZ))) {
    stop("replicated logZ certification received non-finite logZ values.")
  }
  logZ_se[!is.finite(logZ_se) | logZ_se < 0] <- 0
  R <- length(logZ)
  center <- logsumexp(logZ) - log(R)
  empirical_sd <- if (R > 1L) stats::sd(logZ) else 0

  shifted_z <- exp(logZ - max(logZ))
  mean_shifted_z <- mean(shifted_z)
  delta_var_z <- shifted_z^2 * pmax(exp(logZ_se^2) - 1, 0)
  delta_se <- sqrt(sum(delta_var_z) / (R^2)) / max(mean_shifted_z, .Machine$double.eps)

  boot_se <- 0
  bootstrap_B <- as.integer(bootstrap_B)
  if (R > 1L && bootstrap_B > 1L) {
    if (!is.null(seed)) {
      old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
      } else {
        NULL
      }
      on.exit({
        if (is.null(old_seed)) {
          if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
            rm(".Random.seed", envir = .GlobalEnv)
          }
        } else {
          assign(".Random.seed", old_seed, envir = .GlobalEnv)
        }
      }, add = TRUE)
      set.seed(as.integer(seed))
    }
    boot <- replicate(bootstrap_B, {
      idx <- sample.int(R, R, replace = TRUE)
      logsumexp(logZ[idx]) - log(R)
    })
    boot_se <- stats::sd(boot)
    if (!is.finite(boot_se)) {
      boot_se <- 0
    }
  }

  combined_se <- max(empirical_sd, delta_se, boot_se, .Machine$double.eps)
  closest <- which.min(abs(logZ - center))
  run <- runs[[closest]]
  run$logZ <- as.numeric(center)
  run$logZ_se <- as.numeric(combined_se)
  run$diagnostics$replicated_logZ <- logZ
  run$diagnostics$replicated_logZ_se <- logZ_se
  run$diagnostics$replicated_logZ_center <- as.numeric(center)
  run$diagnostics$replicated_logZ_empirical_sd <- as.numeric(empirical_sd)
  run$diagnostics$replicated_logZ_delta_se <- as.numeric(delta_se)
  run$diagnostics$replicated_logZ_bootstrap_se <- as.numeric(boot_se)
  run$diagnostics$replicated_logZ_combined_se <- as.numeric(combined_se)
  list(
    run = run,
    logZ = logZ,
    logZ_se = logZ_se,
    center = as.numeric(center),
    empirical_sd = as.numeric(empirical_sd),
    delta_se = as.numeric(delta_se),
    bootstrap_se = as.numeric(boot_se),
    combined_se = as.numeric(combined_se)
  )
}

.local_chart_normalize_weights <- function(weights, n) {
  if (is.null(weights)) {
    return(rep(1 / n, n))
  }
  weights <- as.numeric(weights)
  if (length(weights) != n) {
    stop("weights length must match the number of alpha particles.")
  }
  if (any(!is.finite(weights)) || any(weights < 0)) {
    stop("weights must be finite and non-negative.")
  }
  sw <- sum(weights)
  if (!is.finite(sw) || sw <= 0) {
    stop("weights must have positive finite sum.")
  }
  weights / sw
}

.local_chart_align_alpha <- function(alpha, population_model) {
  model <- normalize_population_model(population_model)
  alpha <- as.matrix(alpha)
  if (ncol(alpha) != model$alpha_dim) {
    stop("alpha particle dimension does not match the population model.")
  }
  if (!is.null(colnames(alpha)) && setequal(colnames(alpha), model$alpha_names)) {
    alpha <- alpha[, model$alpha_names, drop = FALSE]
  } else {
    colnames(alpha) <- model$alpha_names
  }
  if (any(!is.finite(alpha))) {
    stop("alpha particles must be finite.")
  }
  alpha
}

.local_chart_align_theta_one <- function(theta, population_model) {
  model <- normalize_population_model(population_model)
  theta <- .as_hyper_matrix(theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  if (nrow(theta) != 1L) {
    stop("theta_anchor must contain exactly one theta row.")
  }
  if (any(!is.finite(theta))) {
    stop("theta_anchor must be finite.")
  }
  theta
}

.local_chart_diag_gaussian_parts <- function(theta, alpha, weights, population_model) {
  model <- normalize_population_model(population_model)
  if (!identical(model$fast_family %||% NA_character_, "gaussian")) {
    stop("local chart derivatives currently support only diagonal Gaussian population models.")
  }
  if (model$hyper_dim != 2L * model$alpha_dim) {
    stop("diagonal Gaussian chart derivatives expect hyper_dim = 2 * alpha_dim.")
  }

  theta <- .local_chart_align_theta_one(theta, model)
  alpha <- .local_chart_align_alpha(alpha, model)
  weights <- .local_chart_normalize_weights(weights, nrow(alpha))

  d <- model$alpha_dim
  mu <- as.numeric(theta[1L, seq_len(d)])
  log_sigma2 <- as.numeric(theta[1L, d + seq_len(d)])
  sigma2 <- pmax(exp(log_sigma2), 1e-12)
  centered <- sweep(alpha, 2L, mu, "-")
  inv_sigma2 <- 1 / sigma2

  list(
    model = model,
    theta = theta,
    alpha = alpha,
    weights = weights,
    centered = centered,
    sigma2 = sigma2,
    inv_sigma2 = inv_sigma2
  )
}

.local_chart_reference_prior_from_theta <- function(population_model,
                                                    theta,
                                                    label = "local_chart_reference") {
  model <- normalize_population_model(population_model)
  theta <- .local_chart_align_theta_one(theta, model)
  components <- population_model_reference_components_from_theta(model, theta)
  make_reference_prior_gaussian(
    mu = components$component_means[[1L]],
    Sigma = components$component_covs[[1L]],
    param_names = model$alpha_names,
    label = label
  )
}

.local_chart_assert_smc_complete <- function(fit, label, tol = 1e-10) {
  final_lambda <- as.numeric(fit$final_lambda %||% NA_real_)
  if (!is.finite(final_lambda) || final_lambda < 1 - tol) {
    stop(label, " SMC stopped before lambda_target; local chart logZ is incomplete.")
  }
  if (!is.finite(fit$log_evidence %||% NA_real_)) {
    stop(label, " SMC did not return a finite local chart logZ.")
  }
  invisible(TRUE)
}

.local_chart_fit_diagnostics <- function(fit,
                                         source,
                                         M,
                                         target_cess,
                                         n_mcmc_moves,
                                         max_steps) {
  w <- .local_chart_normalize_weights(fit$w, length(fit$w))
  ess_frac_hist <- as.numeric(fit$meta$ess_frac_hist %||% NA_real_)
  accept_hist <- as.numeric(fit$meta$accept_rate_hist %||% NA_real_)
  list(
    source = source,
    kernel = "run_tempered_smc",
    M = as.integer(M),
    target_cess = as.numeric(target_cess),
    n_mcmc_moves = as.integer(n_mcmc_moves),
    max_steps = as.integer(max_steps),
    final_lambda = as.numeric(fit$final_lambda %||% NA_real_),
    rounds = as.integer(fit$meta$rounds %||% NA_integer_),
    final_ess = as.numeric(1 / sum(w * w)),
    final_ess_frac = as.numeric(1 / sum(w * w) / length(w)),
    min_path_ess_frac = if (any(is.finite(ess_frac_hist))) min(ess_frac_hist, na.rm = TRUE) else NA_real_,
    mean_accept_rate = if (any(is.finite(accept_hist))) mean(accept_hist, na.rm = TRUE) else NA_real_,
    mcse_logZ = as.numeric(fit$mcse_logZ %||% NA_real_)
  )
}

.local_chart_path_is_weak <- function(diagnostics,
                                      min_final_ess_frac,
                                      min_path_ess_frac,
                                      min_accept_rate,
                                      max_logZ_mcse) {
  weak_final_ess <- is.finite(diagnostics$final_ess_frac) &&
    diagnostics$final_ess_frac < as.numeric(min_final_ess_frac)
  weak_path_ess <- is.finite(diagnostics$min_path_ess_frac) &&
    diagnostics$min_path_ess_frac < as.numeric(min_path_ess_frac)
  weak_accept <- is.finite(diagnostics$mean_accept_rate) &&
    diagnostics$mean_accept_rate < as.numeric(min_accept_rate)
  weak_mcse <- is.finite(diagnostics$mcse_logZ) &&
    diagnostics$mcse_logZ > as.numeric(max_logZ_mcse)
  weak_final_ess || weak_path_ess || weak_accept || weak_mcse
}

.local_chart_run_smc <- function(local_id,
                                 theta_anchor,
                                 data_i,
                                 loglik_fn,
                                 population_model,
                                 M,
                                 target_cess,
                                 resample_threshold,
                                 n_mcmc_moves,
                                 rw_scale,
                                 G_mix,
                                 da_enable,
                                 refit_every,
                                 max_steps,
                                 deterministic_resampling,
                                 n_cores,
                                 seed,
                                 verbose,
                                 source) {
  model <- normalize_population_model(population_model)
  theta_anchor <- .local_chart_align_theta_one(theta_anchor, model)
  M <- as.integer(M)
  if (M <= 0L) {
    stop("M must be positive.")
  }

  reference_prior <- .local_chart_reference_prior_from_theta(
    model,
    theta_anchor,
    label = source
  )
  fit <- run_tempered_smc(
    reference_prior = reference_prior,
    bridge_stat_fn = function(alpha) ll_parallel(
      alpha,
      data_i,
      loglik_fn,
      n_cores = as.integer(n_cores)
    ),
    M = M,
    resample_threshold = as.numeric(resample_threshold),
    n_mcmc_moves = as.integer(n_mcmc_moves),
    max_rounds = as.integer(max_steps),
    lambda_target = 1,
    cess_target = as.numeric(target_cess),
    G_mix = as.integer(G_mix),
    refit_every = as.integer(refit_every),
    rw_scale_init = as.numeric(rw_scale),
    post_adapt_n_mcmc_moves = as.integer(n_mcmc_moves),
    da_enable = isTRUE(da_enable),
    deterministic_resampling = isTRUE(deterministic_resampling),
    n_cores = as.integer(n_cores),
    seed = as.integer(seed),
    verbose = isTRUE(verbose)
  )
  .local_chart_assert_smc_complete(fit, source)

  alpha <- .local_chart_align_alpha(fit$Theta, model)
  weights <- .local_chart_normalize_weights(fit$w, nrow(alpha))
  diagnostics <- .local_chart_fit_diagnostics(
    fit = fit,
    source = source,
    M = M,
    target_cess = target_cess,
    n_mcmc_moves = n_mcmc_moves,
    max_steps = max_steps
  )

  list(
    fit = fit,
    alpha = alpha,
    weights = weights,
    logZ = as.numeric(fit$log_evidence),
    logZ_se = as.numeric(fit$mcse_logZ %||% NA_real_),
    score = local_chart_score(theta_anchor, alpha, weights, model),
    curvature = local_chart_curvature(theta_anchor, alpha, weights, model),
    diagnostics = diagnostics
  )
}

.local_chart_quadratic_geometry <- function(theta_anchor,
                                            score,
                                            curvature,
                                            population_model) {
  model <- normalize_population_model(population_model)
  list(
    type = "quadratic_fisher_louis",
    theta_anchor = .local_chart_align_theta_one(theta_anchor, model),
    score_norm = sqrt(sum(as.numeric(score)^2)),
    curvature_frobenius = sqrt(sum(as.matrix(curvature)^2))
  )
}

.local_chart_confirmation_result <- function(primary,
                                             confirmation,
                                             tolerance) {
  combined_se <- sqrt(
    max(primary$logZ_se, 0, na.rm = TRUE)^2 +
      max(confirmation$logZ_se, 0, na.rm = TRUE)^2
  )
  abs_delta <- abs(confirmation$logZ - primary$logZ)
  list(
    primary_logZ = primary$logZ,
    confirmation_logZ = confirmation$logZ,
    delta = confirmation$logZ - primary$logZ,
    abs_delta = abs_delta,
    combined_se = combined_se,
    tolerance = as.numeric(tolerance),
    passed = is.finite(abs_delta) && abs_delta <= as.numeric(tolerance)
  )
}

.local_chart_logspace_add <- function(a, b) {
  m <- pmax(a, b)
  out <- m + log(exp(a - m) + exp(b - m))
  out[!is.finite(a) & !is.finite(b)] <- -Inf
  out
}

.local_chart_weighted_logmeanexp <- function(log_x, weights) {
  log_x <- as.numeric(log_x)
  weights <- .local_chart_normalize_weights(weights, length(log_x))
  finite <- is.finite(log_x) & weights > 0
  if (!any(finite)) {
    return(list(log_mean = -Inf, se_log = Inf, ess_frac = 0, ess = 0))
  }

  lw <- log(weights[finite]) + log_x[finite]
  log_mean <- logsumexp(lw)
  normalized <- exp(lw - log_mean)
  ess <- 1 / sum(normalized * normalized)
  x_scaled <- exp(log_x[finite] - log_mean)
  mean_scaled <- sum(weights[finite] * x_scaled)
  var_mean_scaled <- sum(weights[finite]^2 * (x_scaled - mean_scaled)^2)
  se_log <- sqrt(max(var_mean_scaled, 0)) / max(mean_scaled, .Machine$double.eps)

  list(
    log_mean = as.numeric(log_mean),
    se_log = as.numeric(se_log),
    ess_frac = as.numeric(ess / length(log_x)),
    ess = as.numeric(ess)
  )
}

.local_chart_psis_k <- function(log_ratios) {
  log_ratios <- as.numeric(log_ratios)
  log_ratios <- log_ratios[is.finite(log_ratios)]
  if (length(log_ratios) < 3L || !requireNamespace("loo", quietly = TRUE)) {
    return(NA_real_)
  }
  if (diff(range(log_ratios)) <= 1e-10) {
    return(-Inf)
  }
  psis <- tryCatch(
    suppressWarnings(loo::psis(matrix(log_ratios - max(log_ratios), ncol = 1L))),
    error = function(e) NULL
  )
  if (is.null(psis)) {
    NA_real_
  } else {
    as.numeric(loo::pareto_k_values(psis)[1L])
  }
}

.local_chart_bar_delta <- function(log_ratio_from,
                                   weights_from,
                                   log_ratio_to,
                                   weights_to,
                                   initial_delta,
                                   max_iter = 200L,
                                   tol = 1e-8) {
  weights_from <- .local_chart_normalize_weights(weights_from, length(log_ratio_from))
  weights_to <- .local_chart_normalize_weights(weights_to, length(log_ratio_to))
  n_from <- 1 / sum(weights_from * weights_from)
  n_to <- 1 / sum(weights_to * weights_to)
  log_n_from <- log(n_from)
  log_n_to <- log(n_to)

  delta <- as.numeric(initial_delta)
  if (!is.finite(delta)) {
    delta <- 0
  }
  converged <- FALSE
  max_abs_step <- Inf
  iter <- 0L

  for (iter in seq_len(as.integer(max_iter))) {
    denom_from <- .local_chart_logspace_add(log_n_from, log_n_to + log_ratio_from - delta)
    denom_to <- .local_chart_logspace_add(log_n_from, log_n_to + log_ratio_to - delta)
    log_num <- logsumexp(log(weights_from) + log_ratio_from - denom_from)
    log_den <- logsumexp(log(weights_to) - denom_to)
    proposal <- as.numeric(log_num - log_den)
    if (!is.finite(proposal)) {
      break
    }
    step <- proposal - delta
    max_abs_step <- abs(step)
    delta <- delta + 0.5 * step
    if (is.finite(max_abs_step) && max_abs_step <= as.numeric(tol)) {
      converged <- TRUE
      break
    }
  }

  list(
    delta = as.numeric(delta),
    converged = isTRUE(converged),
    iterations = as.integer(iter),
    max_abs_step = as.numeric(max_abs_step),
    n_from_eff = as.numeric(n_from),
    n_to_eff = as.numeric(n_to)
  )
}

.local_chart_taylor_delta <- function(from_chart, to_chart, population_model) {
  model <- normalize_population_model(population_model)
  theta_from <- .local_chart_align_theta_one(from_chart$theta_anchor, model)
  theta_to <- .local_chart_align_theta_one(to_chart$theta_anchor, model)
  dtheta <- as.numeric(theta_to - theta_from)

  out <- list(
    from = NA_real_,
    to = NA_real_,
    average = NA_real_,
    disagreement = NA_real_
  )
  if (!is.null(from_chart$score) && !is.null(from_chart$curvature)) {
    g_from <- as.numeric(from_chart$score)
    H_from <- as.matrix(from_chart$curvature)
    out$from <- as.numeric(sum(g_from * dtheta) + 0.5 * drop(t(dtheta) %*% H_from %*% dtheta))
  }
  if (!is.null(to_chart$score) && !is.null(to_chart$curvature)) {
    g_to <- as.numeric(to_chart$score)
    H_to <- as.matrix(to_chart$curvature)
    out$to <- as.numeric(sum(g_to * dtheta) - 0.5 * drop(t(dtheta) %*% H_to %*% dtheta))
  }
  finite <- c(out$from, out$to)
  finite <- finite[is.finite(finite)]
  if (length(finite)) {
    out$average <- mean(finite)
  }
  if (is.finite(out$from) && is.finite(out$to)) {
    out$disagreement <- abs(out$from - out$to)
  }
  out
}

new_local_chart <- function(local_id,
                            chart_id,
                            theta_anchor,
                            population_model,
                            status = "candidate",
                            logZ_abs = NA_real_,
                            logZ_abs_se = NA_real_,
                            alpha_particles = NULL,
                            weights = NULL,
                            score = NULL,
                            curvature = NULL,
                            local_geometry = NULL,
                            source_run_ids = character(),
                            edge_ids = character(),
                            diagnostics = list()) {
  model <- normalize_population_model(population_model)
  status <- match.arg(status, .local_chart_status)
  theta_anchor <- .local_chart_align_theta_one(theta_anchor, model)

  structure(
    list(
      local_id = local_id,
      chart_id = as.character(chart_id),
      theta_anchor = theta_anchor,
      status = status,
      logZ_abs = as.numeric(logZ_abs),
      logZ_abs_se = as.numeric(logZ_abs_se),
      alpha_particles = alpha_particles,
      weights = weights,
      score = score,
      curvature = curvature,
      local_geometry = local_geometry,
      source_run_ids = as.character(source_run_ids),
      edge_ids = as.character(edge_ids),
      diagnostics = diagnostics
    ),
    class = "local_evidence_chart"
  ) |>
    validate_local_chart(population_model = model)
}

validate_local_chart <- function(chart, population_model) {
  model <- normalize_population_model(population_model)
  if (!inherits(chart, "local_evidence_chart")) {
    stop("chart must inherit from 'local_evidence_chart'.")
  }
  if (!chart$status %in% .local_chart_status) {
    stop("chart status is invalid.")
  }
  if (length(chart$chart_id) != 1L || !nzchar(chart$chart_id)) {
    stop("chart_id must be a non-empty scalar string.")
  }
  chart$theta_anchor <- .local_chart_align_theta_one(chart$theta_anchor, model)
  chart$logZ_abs <- as.numeric(chart$logZ_abs)
  chart$logZ_abs_se <- as.numeric(chart$logZ_abs_se)
  if (length(chart$logZ_abs) != 1L || length(chart$logZ_abs_se) != 1L) {
    stop("logZ_abs and logZ_abs_se must be scalar numeric values.")
  }

  has_particles <- !is.null(chart$alpha_particles)
  if (has_particles) {
    chart$alpha_particles <- .local_chart_align_alpha(chart$alpha_particles, model)
    chart$weights <- .local_chart_normalize_weights(chart$weights, nrow(chart$alpha_particles))
  } else if (!is.null(chart$weights)) {
    stop("weights cannot be provided without alpha_particles.")
  }

  if (!is.null(chart$score)) {
    chart$score <- as.numeric(chart$score)
    if (length(chart$score) != model$hyper_dim || any(!is.finite(chart$score))) {
      stop("chart score must be a finite vector with one entry per hyperparameter.")
    }
    names(chart$score) <- model$hyper_names
  }

  if (!is.null(chart$curvature)) {
    chart$curvature <- as.matrix(chart$curvature)
    if (nrow(chart$curvature) != model$hyper_dim ||
        ncol(chart$curvature) != model$hyper_dim ||
        any(!is.finite(chart$curvature))) {
      stop("chart curvature must be a finite hyper_dim x hyper_dim matrix.")
    }
    dimnames(chart$curvature) <- list(model$hyper_names, model$hyper_names)
  }

  if (!is.null(chart$local_geometry) && !is.list(chart$local_geometry)) {
    stop("local_geometry must be NULL or a list.")
  }
  if (!is.list(chart$diagnostics)) {
    stop("diagnostics must be a list.")
  }
  chart$source_run_ids <- as.character(chart$source_run_ids %||% character())
  chart$edge_ids <- as.character(chart$edge_ids %||% character())

  if (identical(chart$status, "active")) {
    if (length(chart$logZ_abs) != 1L || !is.finite(chart$logZ_abs)) {
      stop("active charts require finite logZ_abs.")
    }
    if (length(chart$logZ_abs_se) != 1L || !is.finite(chart$logZ_abs_se) || chart$logZ_abs_se < 0) {
      stop("active charts require finite non-negative logZ_abs_se.")
    }
    if (!has_particles) {
      stop("active charts require alpha_particles and weights.")
    }
    if (is.null(chart$score)) {
      stop("active charts require a score vector.")
    }
    if (is.null(chart$curvature)) {
      stop("active charts require a curvature matrix.")
    }
  }

  chart
}

new_local_edge <- function(local_id,
                           edge_id,
                           from_chart,
                           to_chart,
                           delta = NA_real_,
                           se = NA_real_,
                           method = "bridge",
                           overlap_ess = NA_real_,
                           psis_k = NA_real_,
                           forward_delta = NA_real_,
                           reverse_delta = NA_real_,
                           taylor_delta = NA_real_,
                           cycle_residual = NA_real_,
                           status = "candidate",
                           diagnostics = list()) {
  method <- match.arg(method, .local_edge_methods)
  status <- match.arg(status, .local_edge_status)

  structure(
    list(
      local_id = local_id,
      edge_id = as.character(edge_id),
      from_chart = as.character(from_chart),
      to_chart = as.character(to_chart),
      delta = as.numeric(delta),
      se = as.numeric(se),
      method = method,
      overlap_ess = as.numeric(overlap_ess),
      psis_k = as.numeric(psis_k),
      forward_delta = as.numeric(forward_delta),
      reverse_delta = as.numeric(reverse_delta),
      taylor_delta = as.numeric(taylor_delta),
      cycle_residual = as.numeric(cycle_residual),
      status = status,
      diagnostics = diagnostics
    ),
    class = "local_evidence_edge"
  ) |>
    validate_local_edge()
}

validate_local_edge <- function(edge) {
  if (!inherits(edge, "local_evidence_edge")) {
    stop("edge must inherit from 'local_evidence_edge'.")
  }
  if (!edge$status %in% .local_edge_status) {
    stop("edge status is invalid.")
  }
  if (!edge$method %in% .local_edge_methods) {
    stop("edge method is invalid.")
  }
  for (field in c("edge_id", "from_chart", "to_chart")) {
    value <- edge[[field]]
    if (length(value) != 1L || !nzchar(value)) {
      stop(field, " must be a non-empty scalar string.")
    }
  }
  if (identical(edge$from_chart, edge$to_chart)) {
    stop("edge endpoints must be distinct charts.")
  }
  if (!is.list(edge$diagnostics)) {
    stop("edge diagnostics must be a list.")
  }

  numeric_fields <- c(
    "delta", "se", "overlap_ess", "psis_k", "forward_delta",
    "reverse_delta", "taylor_delta", "cycle_residual"
  )
  for (field in numeric_fields) {
    edge[[field]] <- as.numeric(edge[[field]])
    if (length(edge[[field]]) != 1L) {
      stop(field, " must be a scalar numeric value.")
    }
  }

  if (identical(edge$status, "active")) {
    if (!is.finite(edge$delta)) {
      stop("active edges require finite delta.")
    }
    if (!is.finite(edge$se) || edge$se < 0) {
      stop("active edges require finite non-negative se.")
    }
  }

  edge
}

new_local_atlas <- function(local_id,
                            root_chart_id = NULL,
                            charts = list(),
                            edges = list(),
                            normalizer_solution = NULL,
                            coverage_region = NULL,
                            surface_uncertainty_model = NULL,
                            diagnostics = list()) {
  structure(
    list(
      local_id = local_id,
      root_chart_id = if (is.null(root_chart_id)) NULL else as.character(root_chart_id),
      charts = charts,
      edges = edges,
      normalizer_solution = normalizer_solution,
      coverage_region = coverage_region,
      surface_uncertainty_model = surface_uncertainty_model,
      diagnostics = diagnostics
    ),
    class = "local_evidence_atlas"
  ) |>
    validate_local_atlas()
}

validate_local_atlas <- function(atlas) {
  if (!inherits(atlas, "local_evidence_atlas")) {
    stop("atlas must inherit from 'local_evidence_atlas'.")
  }
  if (!is.list(atlas$charts) || !is.list(atlas$edges)) {
    stop("atlas charts and edges must be lists.")
  }
  if (!is.list(atlas$diagnostics)) {
    stop("atlas diagnostics must be a list.")
  }

  if (length(atlas$charts)) {
    if (any(!vapply(atlas$charts, inherits, logical(1), "local_evidence_chart"))) {
      stop("all atlas charts must inherit from 'local_evidence_chart'.")
    }
    if (any(vapply(atlas$charts, function(chart) {
      !identical(as.character(chart$local_id), as.character(atlas$local_id))
    }, logical(1)))) {
      stop("all atlas charts must have the atlas local_id.")
    }
    chart_ids <- vapply(atlas$charts, `[[`, character(1), "chart_id")
    if (anyDuplicated(chart_ids)) {
      stop("atlas chart ids must be unique.")
    }
    names(atlas$charts) <- chart_ids

    if (!is.null(atlas$root_chart_id)) {
      if (length(atlas$root_chart_id) != 1L || !atlas$root_chart_id %in% chart_ids) {
        stop("root_chart_id must identify an atlas chart.")
      }
      if (!identical(atlas$charts[[atlas$root_chart_id]]$status, "active")) {
        stop("root_chart_id must identify an active chart.")
      }
    }
  } else if (!is.null(atlas$root_chart_id)) {
    stop("root_chart_id cannot be set on an atlas with no charts.")
  }

  if (length(atlas$edges)) {
    if (any(!vapply(atlas$edges, inherits, logical(1), "local_evidence_edge"))) {
      stop("all atlas edges must inherit from 'local_evidence_edge'.")
    }
    if (any(vapply(atlas$edges, function(edge) {
      !identical(as.character(edge$local_id), as.character(atlas$local_id))
    }, logical(1)))) {
      stop("all atlas edges must have the atlas local_id.")
    }
    edge_ids <- vapply(atlas$edges, `[[`, character(1), "edge_id")
    if (anyDuplicated(edge_ids)) {
      stop("atlas edge ids must be unique.")
    }
    names(atlas$edges) <- edge_ids

    if (length(atlas$charts)) {
      chart_ids <- names(atlas$charts)
      for (edge in atlas$edges) {
        if (!edge$from_chart %in% chart_ids || !edge$to_chart %in% chart_ids) {
          stop("atlas edge endpoint is not present in atlas charts.")
        }
        if (identical(edge$status, "active")) {
          if (!identical(atlas$charts[[edge$from_chart]]$status, "active") ||
              !identical(atlas$charts[[edge$to_chart]]$status, "active")) {
            stop("active edges may only connect active charts.")
          }
        }
      }
    }
  }

  atlas
}

build_local_root_chart <- function(local_id,
                                   theta_anchor,
                                   data_i,
                                   loglik_fn,
                                   population_model,
                                   chart_id = "root",
                                   M = 1000L,
                                   target_cess = 0.9,
                                   resample_threshold = 0.5,
                                   n_mcmc_moves = 2L,
                                   rw_scale = 0.75,
                                   G_mix = 8L,
                                   da_enable = TRUE,
                                   refit_every = 2L,
                                   max_steps = 128L,
                                   deterministic_resampling = FALSE,
                                   n_cores = 1L,
                                   seed = 123L,
                                   verbose = FALSE,
                                   confirm = c("auto", "always", "never"),
                                   evidence_leverage = 0,
                                   high_leverage_threshold = 0.05,
                                   confirmation_M = max(100L, ceiling(as.integer(M) / 2L)),
                                   confirmation_abs_tol = 1.0,
                                   min_final_ess_frac = 0.25,
                                   min_path_ess_frac = 0.05,
                                   min_accept_rate = 0.02,
                                   max_logZ_mcse = 0.50) {
  model <- normalize_population_model(population_model)
  theta_anchor <- .local_chart_align_theta_one(theta_anchor, model)
  confirm <- match.arg(confirm)
  if (!is.function(loglik_fn)) {
    stop("loglik_fn must be a function.")
  }
  seed <- as.integer(seed)
  if (length(seed) != 1L || !is.finite(seed)) {
    stop("seed must be a finite scalar integer.")
  }

  primary <- .local_chart_run_smc(
    local_id = local_id,
    theta_anchor = theta_anchor,
    data_i = data_i,
    loglik_fn = loglik_fn,
    population_model = model,
    M = as.integer(M),
    target_cess = target_cess,
    resample_threshold = resample_threshold,
    n_mcmc_moves = n_mcmc_moves,
    rw_scale = rw_scale,
    G_mix = G_mix,
    da_enable = da_enable,
    refit_every = refit_every,
    max_steps = max_steps,
    deterministic_resampling = deterministic_resampling,
    n_cores = n_cores,
    seed = seed,
    verbose = verbose,
    source = "local_chart_root_smc"
  )

  weak_path <- .local_chart_path_is_weak(
    diagnostics = primary$diagnostics,
    min_final_ess_frac = min_final_ess_frac,
    min_path_ess_frac = min_path_ess_frac,
    min_accept_rate = min_accept_rate,
    max_logZ_mcse = max_logZ_mcse
  )
  high_leverage <- is.finite(evidence_leverage) &&
    as.numeric(evidence_leverage) >= as.numeric(high_leverage_threshold)
  needs_confirmation <- identical(confirm, "always") ||
    (identical(confirm, "auto") && (weak_path || high_leverage))

  confirmation_info <- NULL
  if (needs_confirmation) {
    confirmation <- .local_chart_run_smc(
      local_id = local_id,
      theta_anchor = theta_anchor,
      data_i = data_i,
      loglik_fn = loglik_fn,
      population_model = model,
      M = as.integer(confirmation_M),
      target_cess = target_cess,
      resample_threshold = resample_threshold,
      n_mcmc_moves = n_mcmc_moves,
      rw_scale = rw_scale,
      G_mix = G_mix,
      da_enable = da_enable,
      refit_every = refit_every,
      max_steps = max_steps,
      deterministic_resampling = deterministic_resampling,
      n_cores = n_cores,
      seed = seed + 1000003L,
      verbose = verbose,
      source = "local_chart_root_confirmation_smc"
    )
    confirmation_info <- .local_chart_confirmation_result(
      primary = primary,
      confirmation = confirmation,
      tolerance = confirmation_abs_tol
    )
    if (!isTRUE(confirmation_info$passed)) {
      stop(
        "Root chart confirmation failed for local ", local_id,
        ": |delta logZ|=", signif(confirmation_info$abs_delta, 4),
        " exceeds tolerance ", signif(confirmation_info$tolerance, 4), "."
      )
    }
  } else if (identical(confirm, "never") && weak_path) {
    stop(
      "Root chart path diagnostics are weak for local ", local_id,
      " and confirmation is disabled."
    )
  }

  diagnostics <- primary$diagnostics
  diagnostics$status_reason <- "certified_root"
  diagnostics$logZ_role <- "direct_smc_observation"
  diagnostics$direct_observation_role <- "anchor"
  diagnostics$weak_path <- isTRUE(weak_path)
  diagnostics$evidence_leverage <- as.numeric(evidence_leverage)
  diagnostics$high_leverage <- isTRUE(high_leverage)
  diagnostics$confirmation_required <- isTRUE(needs_confirmation)
  diagnostics$confirmation <- confirmation_info
  diagnostics$normalizer_certified <- TRUE
  diagnostics$normalizer_certification <- list(
    certified = TRUE,
    method = if (isTRUE(needs_confirmation)) "confirmed_root_smc" else "root_smc",
    reason = if (isTRUE(needs_confirmation)) {
      "root_smc_path_confirmed"
    } else {
      "root_smc_path"
    }
  )

  new_local_chart(
    local_id = local_id,
    chart_id = chart_id,
    theta_anchor = theta_anchor,
    population_model = model,
    status = "active",
    logZ_abs = primary$logZ,
    logZ_abs_se = primary$logZ_se,
    alpha_particles = primary$alpha,
    weights = primary$weights,
    score = primary$score,
    curvature = primary$curvature,
    local_geometry = .local_chart_quadratic_geometry(theta_anchor, primary$score, primary$curvature, model),
    source_run_ids = chart_id,
    diagnostics = diagnostics
  )
}

build_candidate_chart <- function(local_id,
                                  chart_id,
                                  theta_anchor,
                                  data_i,
                                  loglik_fn,
                                  population_model,
                                  M = 1000L,
                                  target_cess = 0.9,
                                  resample_threshold = 0.5,
                                  n_mcmc_moves = 2L,
                                  rw_scale = 0.75,
                                  G_mix = 8L,
                                  da_enable = TRUE,
                                  refit_every = 2L,
                                  max_steps = 128L,
                                  deterministic_resampling = FALSE,
                                  n_cores = 1L,
                                  seed = 123L,
                                  verbose = FALSE) {
  model <- normalize_population_model(population_model)
  theta_anchor <- .local_chart_align_theta_one(theta_anchor, model)
  if (!is.function(loglik_fn)) {
    stop("loglik_fn must be a function.")
  }
  seed <- as.integer(seed)
  if (length(seed) != 1L || !is.finite(seed)) {
    stop("seed must be a finite scalar integer.")
  }

  primary <- .local_chart_run_smc(
    local_id = local_id,
    theta_anchor = theta_anchor,
    data_i = data_i,
    loglik_fn = loglik_fn,
    population_model = model,
    M = as.integer(M),
    target_cess = target_cess,
    resample_threshold = resample_threshold,
    n_mcmc_moves = n_mcmc_moves,
    rw_scale = rw_scale,
    G_mix = G_mix,
    da_enable = da_enable,
    refit_every = refit_every,
    max_steps = max_steps,
    deterministic_resampling = deterministic_resampling,
    n_cores = n_cores,
    seed = seed,
    verbose = verbose,
    source = "local_chart_candidate_smc"
  )

  diagnostics <- primary$diagnostics
  diagnostics$status_reason <- "uncertified_candidate"
  diagnostics$logZ_role <- "direct_smc_observation"
  diagnostics$direct_observation_role <- "support"
  diagnostics$normalizer_certified <- FALSE
  diagnostics$normalizer_certification <- list(
    certified = FALSE,
    method = "candidate_smc",
    reason = "candidate_requires_certified_relative_edge"
  )

  new_local_chart(
    local_id = local_id,
    chart_id = chart_id,
    theta_anchor = theta_anchor,
    population_model = model,
    status = "candidate",
    logZ_abs = primary$logZ,
    logZ_abs_se = primary$logZ_se,
    alpha_particles = primary$alpha,
    weights = primary$weights,
    score = primary$score,
    curvature = primary$curvature,
    local_geometry = .local_chart_quadratic_geometry(theta_anchor, primary$score, primary$curvature, model),
    source_run_ids = chart_id,
    diagnostics = diagnostics
  )
}

.local_atlas_chart_from_smc_run <- function(run,
                                           local_id,
                                           chart_id,
                                           theta_anchor,
                                           population_model,
                                           status = "candidate",
                                           status_reason = "uncertified_candidate",
                                           source_run_id = chart_id,
                                           direct_observation_role = NULL) {
  model <- normalize_population_model(population_model)
  theta_anchor <- .local_chart_align_theta_one(theta_anchor, model)
  diagnostics <- run$diagnostics
  diagnostics$status_reason <- status_reason
  diagnostics$logZ_role <- "direct_smc_observation"
  diagnostics$direct_observation_role <- direct_observation_role
  new_local_chart(
    local_id = local_id,
    chart_id = chart_id,
    theta_anchor = theta_anchor,
    population_model = model,
    status = status,
    logZ_abs = run$logZ,
    logZ_abs_se = run$logZ_se,
    alpha_particles = run$alpha,
    weights = run$weights,
    score = run$score,
    curvature = run$curvature,
    local_geometry = .local_chart_quadratic_geometry(theta_anchor, run$score, run$curvature, model),
    source_run_ids = source_run_id,
    diagnostics = diagnostics
  )
}

estimate_chart_edge <- function(from_chart,
                                to_chart,
                                population_model,
                                edge_id = NULL,
                                method = c("bridge", "BAR"),
                                max_iter = 200L,
                                tol = 1e-8) {
  model <- normalize_population_model(population_model)
  method <- match.arg(method)
  from_chart <- validate_local_chart(from_chart, model)
  to_chart <- validate_local_chart(to_chart, model)
  if (!identical(as.character(from_chart$local_id), as.character(to_chart$local_id))) {
    stop("chart edge endpoints must have the same local_id.")
  }
  if (is.null(from_chart$alpha_particles) || is.null(from_chart$weights) ||
      is.null(to_chart$alpha_particles) || is.null(to_chart$weights)) {
    stop("chart edge estimation requires particles and weights on both endpoint charts.")
  }
  if (identical(from_chart$chart_id, to_chart$chart_id)) {
    stop("chart edge endpoints must be distinct charts.")
  }

  theta_from <- .local_chart_align_theta_one(from_chart$theta_anchor, model)
  theta_to <- .local_chart_align_theta_one(to_chart$theta_anchor, model)
  log_from_on_from <- population_model_log_alpha_given_theta(
    model,
    alpha = from_chart$alpha_particles,
    theta = theta_from
  )
  log_to_on_from <- population_model_log_alpha_given_theta(
    model,
    alpha = from_chart$alpha_particles,
    theta = theta_to
  )
  log_ratio_from <- log_to_on_from - log_from_on_from

  log_from_on_to <- population_model_log_alpha_given_theta(
    model,
    alpha = to_chart$alpha_particles,
    theta = theta_from
  )
  log_to_on_to <- population_model_log_alpha_given_theta(
    model,
    alpha = to_chart$alpha_particles,
    theta = theta_to
  )
  log_ratio_to <- log_to_on_to - log_from_on_to

  forward <- .local_chart_weighted_logmeanexp(log_ratio_from, from_chart$weights)
  reverse_mean <- .local_chart_weighted_logmeanexp(-log_ratio_to, to_chart$weights)
  forward_delta <- forward$log_mean
  reverse_delta <- -reverse_mean$log_mean
  initial_delta <- mean(c(forward_delta, reverse_delta)[is.finite(c(forward_delta, reverse_delta))])
  if (!is.finite(initial_delta)) {
    initial_delta <- 0
  }

  bar <- .local_chart_bar_delta(
    log_ratio_from = log_ratio_from,
    weights_from = from_chart$weights,
    log_ratio_to = log_ratio_to,
    weights_to = to_chart$weights,
    initial_delta = initial_delta,
    max_iter = max_iter,
    tol = tol
  )
  delta <- if (isTRUE(bar$converged) && is.finite(bar$delta)) {
    bar$delta
  } else {
    initial_delta
  }

  psis_forward <- .local_chart_psis_k(log_ratio_from)
  psis_reverse <- .local_chart_psis_k(-log_ratio_to)
  psis_k <- if (all(is.na(c(psis_forward, psis_reverse)))) {
    NA_real_
  } else {
    max(c(psis_forward, psis_reverse), na.rm = TRUE)
  }
  taylor <- .local_chart_taylor_delta(from_chart, to_chart, model)
  forward_reverse_gap <- abs(forward_delta - reverse_delta)
  taylor_gap <- if (is.finite(taylor$average)) abs(delta - taylor$average) else NA_real_
  se <- sqrt(
    max(c(forward$se_log, reverse_mean$se_log), na.rm = TRUE)^2 +
      (0.5 * forward_reverse_gap)^2
  )

  diagnostics <- list(
    method_requested = method,
    bar = bar,
    forward_se = forward$se_log,
    reverse_se = reverse_mean$se_log,
    forward_ess = forward$ess,
    reverse_ess = reverse_mean$ess,
    forward_ess_frac = forward$ess_frac,
    reverse_ess_frac = reverse_mean$ess_frac,
    forward_reverse_gap = as.numeric(forward_reverse_gap),
    psis_forward = psis_forward,
    psis_reverse = psis_reverse,
    taylor_from = taylor$from,
    taylor_to = taylor$to,
    taylor_gap = as.numeric(taylor_gap),
    taylor_disagreement = taylor$disagreement
  )

  new_local_edge(
    local_id = from_chart$local_id,
    edge_id = edge_id %||% paste(from_chart$chart_id, to_chart$chart_id, sep = "->"),
    from_chart = from_chart$chart_id,
    to_chart = to_chart$chart_id,
    delta = delta,
    se = se,
    method = method,
    overlap_ess = min(forward$ess_frac, reverse_mean$ess_frac),
    psis_k = psis_k,
    forward_delta = forward_delta,
    reverse_delta = reverse_delta,
    taylor_delta = taylor$average,
    cycle_residual = NA_real_,
    status = "candidate",
    diagnostics = diagnostics
  )
}

certify_chart_edge <- function(edge,
                               max_forward_reverse_gap = 1.0,
                               max_se = 1.0,
                               min_overlap_ess = 0.05,
                               max_psis_k = 0.7,
                               max_taylor_gap = 2.0,
                               max_taylor_disagreement = Inf,
                               max_cycle_residual = Inf,
                               require_bar_converged = TRUE) {
  edge <- validate_local_edge(edge)
  failures <- character()
  diagnostics <- edge$diagnostics %||% list()

  if (!is.finite(edge$delta)) {
    failures <- c(failures, "nonfinite_delta")
  }
  if (!is.finite(edge$se) || edge$se > as.numeric(max_se)) {
    failures <- c(failures, "large_or_nonfinite_se")
  }
  if (!is.finite(edge$overlap_ess) || edge$overlap_ess < as.numeric(min_overlap_ess)) {
    failures <- c(failures, "low_overlap_ess")
  }
  if (is.finite(edge$psis_k) && edge$psis_k > as.numeric(max_psis_k)) {
    failures <- c(failures, "high_psis_k")
  }
  gap <- diagnostics$forward_reverse_gap %||% abs(edge$forward_delta - edge$reverse_delta)
  if (!is.finite(gap) || gap > as.numeric(max_forward_reverse_gap)) {
    failures <- c(failures, "forward_reverse_disagreement")
  }
  taylor_gap <- diagnostics$taylor_gap %||% if (is.finite(edge$taylor_delta)) abs(edge$delta - edge$taylor_delta) else NA_real_
  if (is.finite(max_taylor_gap) && (!is.finite(taylor_gap) || taylor_gap > as.numeric(max_taylor_gap))) {
    failures <- c(failures, "taylor_disagreement")
  }
  taylor_disagreement <- diagnostics$taylor_disagreement %||% NA_real_
  if (is.finite(max_taylor_disagreement) &&
      (!is.finite(taylor_disagreement) || taylor_disagreement > as.numeric(max_taylor_disagreement))) {
    failures <- c(failures, "endpoint_taylor_disagreement")
  }
  if (is.finite(max_cycle_residual) &&
      (!is.finite(edge$cycle_residual) || abs(edge$cycle_residual) > as.numeric(max_cycle_residual))) {
    failures <- c(failures, "cycle_residual")
  }
  if (isTRUE(require_bar_converged) && !isTRUE((diagnostics$bar %||% list())$converged)) {
    failures <- c(failures, "bar_not_converged")
  }

  diagnostics$certification <- list(
    passed = !length(failures),
    failures = failures,
    thresholds = list(
      max_forward_reverse_gap = as.numeric(max_forward_reverse_gap),
      max_se = as.numeric(max_se),
      min_overlap_ess = as.numeric(min_overlap_ess),
      max_psis_k = as.numeric(max_psis_k),
      max_taylor_gap = as.numeric(max_taylor_gap),
      max_taylor_disagreement = as.numeric(max_taylor_disagreement),
      max_cycle_residual = as.numeric(max_cycle_residual),
      require_bar_converged = isTRUE(require_bar_converged)
    )
  )
  edge$diagnostics <- diagnostics
  edge$status <- if (length(failures)) "rejected" else "active"
  validate_local_edge(edge)
}

.local_atlas_active_graph <- function(atlas) {
  atlas <- validate_local_atlas(atlas)
  chart_status <- vapply(atlas$charts, `[[`, character(1), "status")
  active_chart_ids <- names(atlas$charts)[chart_status == "active"]
  edge_status <- vapply(atlas$edges, `[[`, character(1), "status")
  active_edge_ids <- names(atlas$edges)[edge_status == "active"]
  active_edges <- atlas$edges[active_edge_ids]
  list(
    atlas = atlas,
    active_chart_ids = active_chart_ids,
    active_edge_ids = active_edge_ids,
    active_edges = active_edges
  )
}

.local_atlas_component_from_root <- function(active_chart_ids, active_edges, root_chart_id) {
  if (!length(active_chart_ids)) {
    return(character())
  }
  seen <- root_chart_id
  frontier <- root_chart_id
  while (length(frontier)) {
    current <- frontier[1L]
    frontier <- frontier[-1L]
    neighbors <- character()
    for (edge in active_edges) {
      if (identical(edge$from_chart, current)) {
        neighbors <- c(neighbors, edge$to_chart)
      }
      if (identical(edge$to_chart, current)) {
        neighbors <- c(neighbors, edge$from_chart)
      }
    }
    neighbors <- setdiff(intersect(neighbors, active_chart_ids), seen)
    if (length(neighbors)) {
      seen <- c(seen, neighbors)
      frontier <- c(frontier, neighbors)
    }
  }
  seen
}

.local_atlas_edge_components <- function(active_chart_ids, active_edges) {
  remaining <- as.character(active_chart_ids)
  components <- list()
  while (length(remaining)) {
    seed <- remaining[1L]
    seen <- seed
    frontier <- seed
    while (length(frontier)) {
      current <- frontier[1L]
      frontier <- frontier[-1L]
      neighbors <- character()
      for (edge in active_edges) {
        if (identical(edge$from_chart, current)) {
          neighbors <- c(neighbors, edge$to_chart)
        }
        if (identical(edge$to_chart, current)) {
          neighbors <- c(neighbors, edge$from_chart)
        }
      }
      neighbors <- setdiff(intersect(neighbors, active_chart_ids), seen)
      if (length(neighbors)) {
        seen <- c(seen, neighbors)
        frontier <- c(frontier, neighbors)
      }
    }
    components[[length(components) + 1L]] <- seen
    remaining <- setdiff(remaining, seen)
  }
  components
}

.local_atlas_weighted_graph_solve <- function(A,
                                              rhs,
                                              row_se,
                                              unknown_ids,
                                              robust = TRUE,
                                              robust_method = c("student_t", "huber", "none"),
                                              student_t_df = 30,
                                              huber_k = 2.5,
                                              max_iter = 8L,
                                              tol = 1e-4) {
  A <- as.matrix(A)
  rhs <- as.numeric(rhs)
  row_se <- pmax(as.numeric(row_se), .Machine$double.eps)
  unknown_ids <- as.character(unknown_ids)
  robust_method <- match.arg(robust_method)
  if (!isTRUE(robust)) {
    robust_method <- "none"
  }
  if (ncol(A) != length(unknown_ids)) {
    stop("weighted graph solve column count does not match unknown_ids.")
  }
  if (!nrow(A) || !ncol(A)) {
    stop("weighted graph solve requires a non-empty design matrix.")
  }
  robust_weight <- rep(1, nrow(A))
  x <- rep(0, ncol(A))
  student_t_df <- as.numeric(student_t_df)
  if (!is.finite(student_t_df) || student_t_df <= 0) {
    stop("student_t_df must be positive and finite.")
  }
  huber_k <- as.numeric(huber_k)
  if (!is.finite(huber_k) || huber_k <= 0) {
    huber_k <- 2.5
  }
  max_iter <- if (!identical(robust_method, "none")) max(1L, as.integer(max_iter)) else 1L
  for (iter in seq_len(max_iter)) {
    Aw <- A * (sqrt(robust_weight) / row_se)
    bw <- rhs * sqrt(robust_weight) / row_se
    qr_A <- qr(Aw)
    if (qr_A$rank < ncol(A)) {
      stop("active chart graph is not identifiable from edge and direct chart observations.")
    }
    x_new <- as.numeric(qr.coef(qr_A, bw))
    if (any(!is.finite(x_new))) {
      stop("graph normalizer solve produced non-finite values.")
    }
    resid_std <- as.numeric((A %*% x_new - rhs) / row_se)
    next_weight <- if (identical(robust_method, "student_t")) {
      (student_t_df + 1) / (student_t_df + resid_std * resid_std)
    } else if (identical(robust_method, "huber")) {
      pmin(1, huber_k / pmax(abs(resid_std), .Machine$double.eps))
    } else {
      rep(1, length(resid_std))
    }
    if (iter > 1L &&
        max(abs(x_new - x), na.rm = TRUE) <= as.numeric(tol) &&
        max(abs(next_weight - robust_weight), na.rm = TRUE) <= as.numeric(tol)) {
      x <- x_new
      robust_weight <- next_weight
      break
    }
    x <- x_new
    robust_weight <- next_weight
  }
  names(x) <- unknown_ids
  standardized_residual <- as.numeric((A %*% x - rhs) / row_se)
  objective <- if (identical(robust_method, "student_t")) {
    sum(0.5 * (student_t_df + 1) * log1p(standardized_residual^2 / student_t_df) + log(row_se))
  } else if (identical(robust_method, "huber")) {
    abs_r <- abs(standardized_residual)
    sum(ifelse(abs_r <= huber_k, 0.5 * abs_r^2, huber_k * (abs_r - 0.5 * huber_k)) + log(row_se))
  } else {
    sum(0.5 * standardized_residual^2 + log(row_se))
  }
  Aw <- A * (sqrt(robust_weight) / row_se)
  information <- crossprod(Aw)
  covariance <- tryCatch(
    chol2inv(chol(information)),
    error = function(e) qr.solve(information)
  )
  dimnames(covariance) <- list(unknown_ids, unknown_ids)
  list(
    x = x,
    covariance = covariance,
    robust_weight = robust_weight,
    standardized_residual = standardized_residual,
    effective_se = row_se / sqrt(pmax(robust_weight, .Machine$double.eps)),
    robust_method = robust_method,
    student_t_df = as.numeric(student_t_df),
    objective = as.numeric(objective)
  )
}

.local_atlas_edge_table <- function(atlas, solution) {
  graph <- .local_atlas_active_graph(atlas)
  if (!length(graph$active_edges)) {
    return(data.frame())
  }
  logZ <- solution$logZ
  rows <- lapply(graph$active_edges, function(edge) {
    predicted <- as.numeric(logZ[[edge$to_chart]] - logZ[[edge$from_chart]])
    residual <- predicted - edge$delta
    se <- max(as.numeric(edge$se), 1e-6, .Machine$double.eps)
    data.frame(
      edge_id = edge$edge_id,
      from_chart = edge$from_chart,
      to_chart = edge$to_chart,
      observed_delta = as.numeric(edge$delta),
      predicted_delta = predicted,
      residual = residual,
      abs_residual = abs(residual),
      se = as.numeric(edge$se),
      standardized_residual = residual / se,
      weight = 1 / (se * se),
      check.names = FALSE
    )
  })
  do.call(rbind, rows)
}

atlas_cycle_diagnostics <- function(atlas, normalizer_solution = NULL) {
  atlas <- validate_local_atlas(atlas)
  solution <- normalizer_solution %||% atlas$normalizer_solution
  if (is.null(solution) || is.null(solution$logZ)) {
    stop("atlas_cycle_diagnostics requires a normalizer solution.")
  }
  graph <- .local_atlas_active_graph(atlas)
  edge_table <- .local_atlas_edge_table(atlas, solution)
  n_active_charts <- length(graph$active_chart_ids)
  n_active_edges <- length(graph$active_edge_ids)
  cycle_df <- max(0L, n_active_edges - n_active_charts + 1L)
  weighted_rss <- if (nrow(edge_table)) {
    sum(edge_table$weight * edge_table$residual^2)
  } else {
    0
  }
  list(
    edge_residuals = edge_table,
    summary = list(
      n_active_charts = as.integer(n_active_charts),
      n_active_edges = as.integer(n_active_edges),
      cycle_degrees_freedom = as.integer(cycle_df),
      weighted_rss = as.numeric(weighted_rss),
      max_abs_edge_residual = if (nrow(edge_table)) max(edge_table$abs_residual) else 0,
      max_abs_standardized_edge_residual = if (nrow(edge_table)) {
        max(abs(edge_table$standardized_residual))
      } else {
        0
      }
    )
  )
}

solve_atlas_normalizers <- function(atlas,
                                    root_chart_id = atlas$root_chart_id,
                                    require_connected = TRUE,
                                    update_charts = TRUE,
                                    update_edges = TRUE,
                                    se_floor = 1e-6,
                                    use_direct_observations = TRUE,
                                    direct_se_inflation = c(anchor = 15, calibration = 3),
                                    direct_overdispersion = c(anchor = 0, calibration = 0.5),
                                    edge_overdispersion = 0,
                                    direct_observation_roles = c("anchor", "calibration"),
                                    root_mode = c("auto", "latent", "fixed"),
                                    robust = getOption("local_charts.normalizer_robust", TRUE),
                                    robust_method = getOption("local_charts.normalizer_robust_method", "student_t"),
                                    student_t_df = getOption("local_charts.normalizer_student_t_df", 30),
                                    robust_huber_k = 2.5,
                                    robust_max_iter = 8L) {
  atlas <- validate_local_atlas(atlas)
  root_mode <- match.arg(root_mode)
  robust_method <- match.arg(robust_method)
  if (!isTRUE(robust)) {
    robust_method <- "none"
  }
  if (is.null(root_chart_id) || length(root_chart_id) != 1L) {
    stop("root_chart_id must be supplied.")
  }
  root_chart_id <- as.character(root_chart_id)
  graph <- .local_atlas_active_graph(atlas)
  active_chart_ids <- graph$active_chart_ids
  active_edge_ids <- graph$active_edge_ids
  active_edges <- graph$active_edges
  if (!root_chart_id %in% active_chart_ids) {
    stop("root_chart_id must identify an active chart.")
  }

  component <- .local_atlas_component_from_root(active_chart_ids, active_edges, root_chart_id)
  disconnected <- setdiff(active_chart_ids, component)
  solve_chart_ids <- if (isTRUE(require_connected)) active_chart_ids else component
  solve_edge_ids <- active_edge_ids[vapply(active_edges, function(edge) {
    edge$from_chart %in% solve_chart_ids && edge$to_chart %in% solve_chart_ids
  }, logical(1))]
  solve_edges <- atlas$edges[solve_edge_ids]

  root_logZ <- as.numeric(atlas$charts[[root_chart_id]]$logZ_abs)
  root_se <- as.numeric(atlas$charts[[root_chart_id]]$logZ_abs_se)
  if (!.local_chart_normalizer_certified(atlas$charts[[root_chart_id]])) {
    stop("root chart normalizer must be explicitly certified before graph solution.")
  }
  if (!is.finite(root_logZ)) {
    stop("root chart must have finite logZ_abs.")
  }
  if (!is.finite(root_se) || root_se < 0) {
    stop("root chart must have finite non-negative logZ_abs_se.")
  }

  direct_ids <- character()
  direct_rhs <- numeric()
  direct_se <- numeric()
  direct_nominal_se <- numeric()
  direct_role <- character()
  direct_inflation <- numeric()
  direct_tau <- numeric()
  direct_se_inflation_for_role <- function(role) {
    x <- direct_se_inflation
    if (is.null(x) || !length(x)) {
      return(15)
    }
    if (!is.null(names(x)) && role %in% names(x)) {
      out <- as.numeric(x[[role]])
    } else {
      out <- as.numeric(x[[1L]])
    }
    if (!is.finite(out) || out <= 0) {
      stop("direct_se_inflation values must be positive finite numbers.")
    }
    out
  }
  direct_overdispersion_for_role <- function(role) {
    .local_chart_observation_scale(direct_overdispersion, role, default = 0)
  }
  edge_tau <- .local_chart_observation_scale(edge_overdispersion, "edge", default = 0)
  if (isTRUE(use_direct_observations) && length(solve_chart_ids)) {
    roles <- as.character(direct_observation_roles)
    for (chart_id in solve_chart_ids) {
      chart <- atlas$charts[[chart_id]]
      role <- as.character(chart$diagnostics$direct_observation_role %||% "")
      if (!role %in% roles || !.local_chart_normalizer_certified(chart)) {
        next
      }
      observed <- as.numeric(chart$diagnostics$pre_graph_logZ_abs %||% chart$logZ_abs)
      observed_se <- as.numeric(chart$diagnostics$pre_graph_logZ_abs_se %||% chart$logZ_abs_se)
      if (is.finite(observed) && is.finite(observed_se) && observed_se >= 0) {
        inflation <- direct_se_inflation_for_role(role)
        tau <- direct_overdispersion_for_role(role)
        nominal_se <- max(observed_se * inflation, as.numeric(se_floor), .Machine$double.eps)
        direct_ids <- c(direct_ids, chart_id)
        direct_rhs <- c(direct_rhs, observed)
        direct_nominal_se <- c(direct_nominal_se, nominal_se)
        direct_se <- c(direct_se, sqrt(nominal_se^2 + tau^2))
        direct_role <- c(direct_role, role)
        direct_inflation <- c(direct_inflation, inflation)
        direct_tau <- c(direct_tau, tau)
      }
    }
  }

  if (length(disconnected) && isTRUE(require_connected)) {
    components <- .local_atlas_edge_components(active_chart_ids, active_edges)
    unanchored <- unlist(Filter(
      function(ids) !any(ids %in% direct_ids),
      components
    ), use.names = FALSE)
    if (length(unanchored)) {
      stop(
        "active chart graph has components without a direct normalizer observation: ",
        paste(unanchored, collapse = ", ")
      )
    }
  }

  use_latent_root <- identical(root_mode, "latent") ||
    (identical(root_mode, "auto") && length(direct_ids) > 0L)

  logZ <- stats::setNames(rep(NA_real_, length(active_chart_ids)), active_chart_ids)
  rel_se <- stats::setNames(rep(NA_real_, length(active_chart_ids)), active_chart_ids)
  abs_se <- rel_se
  covariance <- matrix(0, nrow = length(active_chart_ids), ncol = length(active_chart_ids),
                       dimnames = list(active_chart_ids, active_chart_ids))
  graph_objective <- NA_real_

  if (isTRUE(use_latent_root)) {
    if (!length(direct_ids)) {
      stop("latent root normalizer solve requires at least one direct chart observation.")
    }
    n_rows <- length(solve_edges) + length(direct_ids)
    A <- matrix(0, nrow = n_rows, ncol = length(solve_chart_ids),
                dimnames = list(NULL, solve_chart_ids))
    rhs <- numeric(n_rows)
    row_se <- numeric(n_rows)
    row_type <- character(n_rows)
    row_chart_id <- rep(NA_character_, n_rows)
    row_edge_id <- rep(NA_character_, n_rows)

    for (r in seq_along(solve_edges)) {
      edge <- solve_edges[[r]]
      row_edge_id[r] <- edge$edge_id
      row_type[r] <- "edge"
      A[r, edge$to_chart] <- 1
      A[r, edge$from_chart] <- -1
      rhs[r] <- edge$delta
      edge_se <- max(as.numeric(edge$se), as.numeric(se_floor), .Machine$double.eps)
      row_se[r] <- sqrt(edge_se^2 + edge_tau^2)
    }

    if (length(direct_ids)) {
      offset <- length(solve_edges)
      for (k in seq_along(direct_ids)) {
        r <- offset + k
        chart_id <- direct_ids[k]
        A[r, chart_id] <- 1
        rhs[r] <- direct_rhs[k]
        row_se[r] <- direct_se[k]
        row_type[r] <- "direct"
        row_chart_id[r] <- chart_id
      }
    }

    rownames(A) <- ifelse(
      row_type == "edge",
      paste0("edge:", row_edge_id),
      paste0("direct:", row_chart_id)
    )

    graph_fit <- .local_atlas_weighted_graph_solve(
      A = A,
      rhs = rhs,
      row_se = row_se,
      unknown_ids = solve_chart_ids,
      robust = robust,
      robust_method = robust_method,
      student_t_df = student_t_df,
      huber_k = robust_huber_k,
      max_iter = robust_max_iter
    )
    x <- graph_fit$x
    graph_objective <- as.numeric(graph_fit$objective %||% NA_real_)
    logZ[solve_chart_ids] <- x

    cov_solved <- graph_fit$covariance
    covariance[solve_chart_ids, solve_chart_ids] <- cov_solved
    abs_se[solve_chart_ids] <- sqrt(pmax(diag(cov_solved), 0))
    rel_var <- pmax(
      diag(cov_solved) +
        cov_solved[root_chart_id, root_chart_id] -
        2 * cov_solved[, root_chart_id],
      0
    )
    rel_se[solve_chart_ids] <- sqrt(rel_var)
    rel_se[root_chart_id] <- 0

    direct_residuals <- data.frame(
      chart_id = direct_ids,
      observed_logZ = direct_rhs,
      fitted_logZ = as.numeric(logZ[direct_ids]),
      residual = as.numeric(logZ[direct_ids]) - direct_rhs,
      abs_residual = abs(as.numeric(logZ[direct_ids]) - direct_rhs),
      se = direct_se,
      nominal_se = direct_nominal_se,
      standardized_residual = (as.numeric(logZ[direct_ids]) - direct_rhs) / direct_se,
      role = direct_role,
      se_inflation = direct_inflation,
      overdispersion = direct_tau,
      robust_weight = graph_fit$robust_weight[row_type == "direct"],
      check.names = FALSE
    )
  } else {
    unknown_ids <- setdiff(solve_chart_ids, root_chart_id)
    logZ[root_chart_id] <- root_logZ
    rel_se[root_chart_id] <- 0
    abs_se[root_chart_id] <- root_se

    direct_keep <- direct_ids %in% unknown_ids
    direct_ids <- direct_ids[direct_keep]
    direct_rhs <- direct_rhs[direct_keep]
    direct_se <- direct_se[direct_keep]
    direct_nominal_se <- direct_nominal_se[direct_keep]
    direct_role <- direct_role[direct_keep]
    direct_inflation <- direct_inflation[direct_keep]
    direct_tau <- direct_tau[direct_keep]

    if (length(unknown_ids)) {
    if (!length(solve_edges) && !length(direct_ids)) {
      stop("cannot solve non-root chart normalizers without active edges or direct chart observations.")
    }
    n_rows <- length(solve_edges) + length(direct_ids)
    A <- matrix(0, nrow = n_rows, ncol = length(unknown_ids),
                dimnames = list(NULL, unknown_ids))
    rhs <- numeric(n_rows)
    row_se <- numeric(n_rows)
    row_type <- character(n_rows)
    row_chart_id <- rep(NA_character_, n_rows)
    row_edge_id <- rep(NA_character_, n_rows)

    for (r in seq_along(solve_edges)) {
      edge <- solve_edges[[r]]
      row_edge_id[r] <- edge$edge_id
      row_type[r] <- "edge"
      if (edge$to_chart %in% unknown_ids) {
        A[r, edge$to_chart] <- A[r, edge$to_chart] + 1
      }
      if (edge$from_chart %in% unknown_ids) {
        A[r, edge$from_chart] <- A[r, edge$from_chart] - 1
      }
      fixed <- 0
      if (identical(edge$to_chart, root_chart_id)) {
        fixed <- fixed + root_logZ
      }
      if (identical(edge$from_chart, root_chart_id)) {
        fixed <- fixed - root_logZ
      }
      rhs[r] <- edge$delta - fixed
      edge_se <- max(as.numeric(edge$se), as.numeric(se_floor), .Machine$double.eps)
      row_se[r] <- sqrt(edge_se^2 + edge_tau^2)
    }

    if (length(direct_ids)) {
      offset <- length(solve_edges)
      for (k in seq_along(direct_ids)) {
        r <- offset + k
        chart_id <- direct_ids[k]
        A[r, chart_id] <- 1
        rhs[r] <- direct_rhs[k]
        row_se[r] <- direct_se[k]
        row_type[r] <- "direct"
        row_chart_id[r] <- chart_id
      }
    }

    rownames(A) <- ifelse(
      row_type == "edge",
      paste0("edge:", row_edge_id),
      paste0("direct:", row_chart_id)
    )

    graph_fit <- .local_atlas_weighted_graph_solve(
      A = A,
      rhs = rhs,
      row_se = row_se,
      unknown_ids = unknown_ids,
      robust = robust,
      robust_method = robust_method,
      student_t_df = student_t_df,
      huber_k = robust_huber_k,
      max_iter = robust_max_iter
    )
    x <- graph_fit$x
    graph_objective <- as.numeric(graph_fit$objective %||% NA_real_)
    logZ[unknown_ids] <- x

    cov_unknown <- graph_fit$covariance
    covariance[unknown_ids, unknown_ids] <- cov_unknown
    if (is.finite(root_se) && root_se > 0) {
      covariance[root_chart_id, root_chart_id] <- root_se^2
      covariance[unknown_ids, root_chart_id] <- root_se^2
      covariance[root_chart_id, unknown_ids] <- root_se^2
      covariance[unknown_ids, unknown_ids] <- covariance[unknown_ids, unknown_ids] + root_se^2
    }
    rel_se[unknown_ids] <- sqrt(pmax(diag(cov_unknown), 0))
    abs_se[unknown_ids] <- sqrt(pmax(diag(covariance)[unknown_ids], 0))

    direct_residuals <- if (length(direct_ids)) {
      fitted <- as.numeric(logZ[direct_ids])
      data.frame(
        chart_id = direct_ids,
        observed_logZ = direct_rhs,
        fitted_logZ = fitted,
        residual = fitted - direct_rhs,
        abs_residual = abs(fitted - direct_rhs),
        se = direct_se,
        nominal_se = direct_nominal_se,
        standardized_residual = (fitted - direct_rhs) / direct_se,
        role = direct_role,
        se_inflation = direct_inflation,
        overdispersion = direct_tau,
        robust_weight = graph_fit$robust_weight[row_type == "direct"],
        check.names = FALSE
      )
    } else {
      data.frame()
    }
  } else {
    direct_residuals <- data.frame()
  }
  }

  solution <- list(
    method = if (isTRUE(use_latent_root)) {
      if (!identical(robust_method, "none")) paste0("latent_", robust_method, "_graph") else "latent_weighted_graph_least_squares"
    } else {
      if (!identical(robust_method, "none")) paste0("root_fixed_", robust_method, "_graph") else "root_fixed_weighted_graph_least_squares"
    },
    root_chart_id = root_chart_id,
    root_logZ = as.numeric(logZ[[root_chart_id]]),
    root_logZ_se = as.numeric(abs_se[[root_chart_id]]),
    active_chart_ids = active_chart_ids,
    active_edge_ids = active_edge_ids,
    solved_chart_ids = solve_chart_ids,
    disconnected_chart_ids = disconnected,
    logZ = logZ,
    relative_se = rel_se,
    se = abs_se,
    covariance = covariance,
    root_mode = if (isTRUE(use_latent_root)) "latent" else "fixed",
    direct_se_inflation = direct_se_inflation,
    direct_overdispersion = direct_overdispersion,
    edge_overdispersion = edge_tau,
    direct_observation_roles = as.character(direct_observation_roles),
    connected = !length(disconnected),
    require_connected = isTRUE(require_connected),
    robust = isTRUE(robust),
    robust_method = robust_method,
    student_t_df = as.numeric(student_t_df),
    robust_huber_k = as.numeric(robust_huber_k),
    robust_max_iter = as.integer(robust_max_iter),
    robust_objective = as.numeric(graph_objective)
  )
  cycle <- atlas_cycle_diagnostics(atlas, solution)
  solution$residuals <- cycle$edge_residuals
  solution$direct_residuals <- direct_residuals
  solution$diagnostics <- cycle$summary
  solution$diagnostics$n_direct_observations <- as.integer(nrow(direct_residuals))
  solution$diagnostics$max_abs_direct_residual <- if (nrow(direct_residuals)) {
    max(direct_residuals$abs_residual, na.rm = TRUE)
  } else {
    0
  }
  solution$diagnostics$max_abs_standardized_direct_residual <- if (nrow(direct_residuals)) {
    max(abs(direct_residuals$standardized_residual), na.rm = TRUE)
  } else {
    0
  }

  if (isTRUE(update_charts)) {
    for (chart_id in solve_chart_ids) {
      chart <- atlas$charts[[chart_id]]
      chart$diagnostics$pre_graph_logZ_abs <- chart$diagnostics$pre_graph_logZ_abs %||% chart$logZ_abs
      chart$diagnostics$pre_graph_logZ_abs_se <- chart$diagnostics$pre_graph_logZ_abs_se %||% chart$logZ_abs_se
      chart$logZ_abs <- as.numeric(logZ[[chart_id]])
      chart$logZ_abs_se <- as.numeric(abs_se[[chart_id]])
      chart$diagnostics$graph_normalizer <- list(
        method = solution$method,
        relative_se = as.numeric(rel_se[[chart_id]]),
        absolute_se = as.numeric(abs_se[[chart_id]]),
        root_chart_id = root_chart_id
      )
      if (nrow(direct_residuals) && chart_id %in% direct_residuals$chart_id) {
        chart$diagnostics$graph_direct_observation <- as.list(
          direct_residuals[match(chart_id, direct_residuals$chart_id), , drop = FALSE]
        )
      } else {
        chart$diagnostics$graph_direct_observation <- NULL
      }
      atlas$charts[[chart_id]] <- chart
    }
  }

  if (isTRUE(update_edges) && nrow(solution$residuals)) {
    residuals <- solution$residuals
    for (k in seq_len(nrow(residuals))) {
      edge_id <- residuals$edge_id[k]
      edge <- atlas$edges[[edge_id]]
      edge$cycle_residual <- as.numeric(residuals$residual[k])
      edge$diagnostics$graph_residual <- as.list(residuals[k, , drop = FALSE])
      atlas$edges[[edge_id]] <- validate_local_edge(edge)
    }
  }

  atlas$root_chart_id <- root_chart_id
  atlas$normalizer_solution <- solution
  atlas$diagnostics$normalizer_solution <- solution$diagnostics
  validate_local_atlas(atlas)
}

.local_atlas_active_charts <- function(atlas) {
  atlas <- validate_local_atlas(atlas)
  if (!length(atlas$charts)) {
    return(list())
  }
  active <- atlas$charts[vapply(atlas$charts, function(chart) identical(chart$status, "active"), logical(1))]
  active
}

.local_chart_quadratic_predict <- function(chart,
                                           theta,
                                           population_model,
                                           distance_metric = "euclidean") {
  model <- normalize_population_model(population_model)
  chart <- validate_local_chart(chart, model)
  theta <- .as_hyper_matrix(theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  if (nrow(theta) != 1L) {
    stop("quadratic chart prediction expects one theta row.")
  }
  dtheta <- as.numeric(theta - chart$theta_anchor)
  value <- as.numeric(chart$logZ_abs) +
    sum(as.numeric(chart$score) * dtheta) +
    0.5 * drop(t(dtheta) %*% as.matrix(chart$curvature) %*% dtheta)
  list(
    value = as.numeric(value),
    distance = as.numeric(.local_atlas_theta_distances(
      theta = theta,
      centers = chart$theta_anchor,
      population_model = model,
      distance_metric = distance_metric
    )),
    dtheta = dtheta
  )
}

.local_atlas_prediction_rows <- function(atlas,
                                         theta,
                                         population_model,
                                         active_charts,
                                         max_chart_distance,
                                         distance_metric = "euclidean") {
  rows <- lapply(active_charts, function(chart) {
    pred <- .local_chart_quadratic_predict(
      chart,
      theta,
      population_model,
      distance_metric = distance_metric
    )
    se <- as.numeric(chart$logZ_abs_se)
    if (!is.finite(se) || se < 0) {
      se <- Inf
    }
    data.frame(
      chart_id = chart$chart_id,
      log_marginal = pred$value,
      chart_se = se,
      distance = pred$distance,
      within_distance = is.finite(pred$distance) && pred$distance <= as.numeric(max_chart_distance),
      check.names = FALSE
    )
  })
  do.call(rbind, rows)
}

.local_atlas_combine_predictions <- function(rows,
                                             distance_scale,
                                             se_floor) {
  scale <- max(as.numeric(distance_scale), .Machine$double.eps)
  se <- pmax(as.numeric(rows$chart_se), as.numeric(se_floor), .Machine$double.eps)
  log_precision <- -2 * log(se)
  log_distance_weight <- -0.5 * (rows$distance / scale)^2
  log_w <- log_precision + log_distance_weight
  log_w <- log_w - logsumexp(log_w)
  w <- exp(log_w)
  value <- sum(w * rows$log_marginal)
  between <- sum(w * (rows$log_marginal - value)^2)
  within <- sum(w * se^2)
  list(
    log_marginal = as.numeric(value),
    se = as.numeric(sqrt(max(within + between, 0))),
    weights = w,
    between_sd = sqrt(max(between, 0)),
    within_se = sqrt(max(within, 0))
  )
}

.local_atlas_surface_method <- function(surface_method) {
  surface_method <- as.character(surface_method %||% "derivative_ls")
  if (!length(surface_method) || is.na(surface_method[1L]) || !nzchar(surface_method[1L])) {
    return("derivative_ls")
  }
  surface_method <- surface_method[1L]
  allowed <- c("derivative_ls", "chart_average")
  if (!surface_method %in% allowed) {
    stop("surface_method must be one of: ", paste(allowed, collapse = ", "))
  }
  surface_method
}

.local_atlas_particle_mis_role <- function(particle_mis_role) {
  particle_mis_role <- as.character(particle_mis_role %||% "estimator")
  if (!length(particle_mis_role) || is.na(particle_mis_role[1L]) || !nzchar(particle_mis_role[1L])) {
    return("estimator")
  }
  particle_mis_role <- particle_mis_role[1L]
  allowed <- "estimator"
  if (!particle_mis_role %in% allowed) {
    stop("particle_mis_role must be one of: ", paste(allowed, collapse = ", "))
  }
  particle_mis_role
}

.local_atlas_quadratic_param_pairs <- function(d) {
  pairs <- do.call(rbind, lapply(seq_len(d), function(i) {
    cbind(i = i, j = i:d)
  }))
  storage.mode(pairs) <- "integer"
  pairs
}

.local_atlas_quadratic_value_basis <- function(z, pairs) {
  z <- as.numeric(z)
  d <- length(z)
  out <- numeric(1L + d + nrow(pairs))
  out[1L] <- 1
  out[1L + seq_len(d)] <- z
  offset <- 1L + d
  for (k in seq_len(nrow(pairs))) {
    i <- pairs[k, 1L]
    j <- pairs[k, 2L]
    out[offset + k] <- if (i == j) 0.5 * z[i] * z[j] else z[i] * z[j]
  }
  out
}

.local_atlas_quadratic_gradient_basis <- function(z, component, pairs) {
  z <- as.numeric(z)
  d <- length(z)
  component <- as.integer(component)
  out <- numeric(1L + d + nrow(pairs))
  out[1L + component] <- 1
  offset <- 1L + d
  for (k in seq_len(nrow(pairs))) {
    i <- pairs[k, 1L]
    j <- pairs[k, 2L]
    if (i == component && j == component) {
      out[offset + k] <- z[component]
    } else if (i == component) {
      out[offset + k] <- z[j]
    } else if (j == component) {
      out[offset + k] <- z[i]
    }
  }
  out
}

.local_atlas_quadratic_hessian_basis <- function(d, row, col, pairs) {
  out <- numeric(1L + d + nrow(pairs))
  lo <- min(as.integer(row), as.integer(col))
  hi <- max(as.integer(row), as.integer(col))
  idx <- which(pairs[, 1L] == lo & pairs[, 2L] == hi)
  if (length(idx)) {
    out[1L + d + idx[1L]] <- 1
  }
  out
}

.local_atlas_surface_coordinate_scale <- function(charts, population_model) {
  model <- normalize_population_model(population_model)
  anchors <- do.call(rbind, lapply(charts, function(chart) chart$theta_anchor))
  colnames(anchors) <- model$hyper_names
  scale <- apply(anchors, 2L, stats::sd)
  range_scale <- apply(anchors, 2L, function(x) diff(range(x)))
  scale[!is.finite(scale) | scale <= 1e-8] <- range_scale[!is.finite(scale) | scale <= 1e-8] / 2
  scale[!is.finite(scale) | scale <= 1e-8] <- 1
  names(scale) <- model$hyper_names
  scale
}

.local_atlas_derivative_surface <- function(theta,
                                            charts,
                                            rows,
                                            population_model,
                                            distance_scale,
                                            se_floor = 1e-6,
                                            value_nugget = 0.05,
                                            gradient_weight = 1.0,
                                            curvature_weight = 0.2,
                                            ridge = 1e-8,
                                            min_charts = 2L,
                                            omit_chart_id = NULL) {
  model <- normalize_population_model(population_model)
  theta <- .as_hyper_matrix(theta, model$hyper_names, model$hyper_dim)
  if (nrow(theta) != 1L) {
    stop("derivative surface expects one theta row.")
  }
  rows <- rows[!rows$chart_id %in% as.character(omit_chart_id %||% character()), , drop = FALSE]
  charts <- charts[rows$chart_id]
  charts <- charts[!vapply(charts, is.null, logical(1))]
  if (length(charts) < as.integer(min_charts)) {
    return(list(status = "uncertified", reason = "too_few_surface_charts"))
  }

  d <- model$hyper_dim
  pairs <- .local_atlas_quadratic_param_pairs(d)
  p <- 1L + d + nrow(pairs)
  coord_scale <- .local_atlas_surface_coordinate_scale(charts, model)
  bw <- as.numeric(distance_scale)
  if (!is.finite(bw) || bw <= 0) {
    finite_distance <- as.numeric(rows$distance[is.finite(rows$distance)])
    bw <- if (length(finite_distance)) max(stats::median(finite_distance), 1e-8) else 1
  }
  value_nugget <- max(as.numeric(value_nugget), as.numeric(se_floor), 1e-8)
  gradient_sd <- 1 / sqrt(max(as.numeric(gradient_weight), .Machine$double.eps))
  curvature_sd <- 1 / sqrt(max(as.numeric(curvature_weight), .Machine$double.eps))

  x_rows <- list()
  y <- numeric(0)
  weight <- numeric(0)
  obs_chart <- character()
  obs_type <- character()

  add_obs <- function(basis, response, w, chart_id, type) {
    if (!is.finite(response) || !is.finite(w) || w <= 0) {
      return()
    }
    x_rows[[length(x_rows) + 1L]] <<- basis
    y <<- c(y, as.numeric(response))
    weight <<- c(weight, as.numeric(w))
    obs_chart <<- c(obs_chart, as.character(chart_id))
    obs_type <<- c(obs_type, as.character(type))
  }

  for (chart_id in names(charts)) {
    chart <- validate_local_chart(charts[[chart_id]], model)
    row <- rows[match(chart_id, rows$chart_id), , drop = FALSE]
    if (!nrow(row)) next
    z <- as.numeric((chart$theta_anchor[1L, ] - theta[1L, ]) / coord_scale)
    kernel <- exp(-0.5 * (as.numeric(row$distance[1L]) / bw)^2)
    chart_se <- as.numeric(chart$logZ_abs_se)
    if (!is.finite(chart_se) || chart_se < 0) {
      chart_se <- value_nugget
    }
    value_sd <- sqrt(chart_se^2 + value_nugget^2)
    add_obs(
      basis = .local_atlas_quadratic_value_basis(z, pairs),
      response = as.numeric(chart$logZ_abs),
      w = kernel / value_sd^2,
      chart_id = chart_id,
      type = "value"
    )

    if (!is.null(chart$score)) {
      gradient <- as.numeric(chart$score) * coord_scale
      for (j in seq_len(d)) {
        add_obs(
          basis = .local_atlas_quadratic_gradient_basis(z, j, pairs),
          response = gradient[j],
          w = kernel / gradient_sd^2,
          chart_id = chart_id,
          type = "gradient"
        )
      }
    }
    if (!is.null(chart$curvature) && is.finite(curvature_weight) && curvature_weight > 0) {
      curvature <- as.matrix(chart$curvature) * outer(coord_scale, coord_scale)
      for (k in seq_len(nrow(pairs))) {
        i <- pairs[k, 1L]
        j <- pairs[k, 2L]
        add_obs(
          basis = .local_atlas_quadratic_hessian_basis(d, i, j, pairs),
          response = curvature[i, j],
          w = kernel / curvature_sd^2,
          chart_id = chart_id,
          type = "curvature"
        )
      }
    }
  }

  if (!length(x_rows) || length(unique(obs_chart)) < as.integer(min_charts)) {
    return(list(status = "uncertified", reason = "too_few_surface_observations"))
  }
  X <- do.call(rbind, x_rows)
  w <- pmax(as.numeric(weight), 0)
  finite <- is.finite(y) & rowSums(is.finite(X)) == ncol(X) & is.finite(w) & w > 0
  X <- X[finite, , drop = FALSE]
  y <- y[finite]
  w <- w[finite]
  obs_chart <- obs_chart[finite]
  obs_type <- obs_type[finite]
  if (nrow(X) < 1L || length(unique(obs_chart)) < as.integer(min_charts)) {
    return(list(status = "uncertified", reason = "too_few_surface_observations"))
  }

  sqrt_w <- sqrt(w / max(mean(w), .Machine$double.eps))
  Xw <- X * sqrt_w
  yw <- y * sqrt_w
  penalty <- diag(c(0, rep(as.numeric(ridge), p - 1L)), nrow = p)
  xtx <- crossprod(Xw) + penalty
  rhs <- as.numeric(crossprod(Xw, yw))
  xtx_svd <- tryCatch(
    svd(xtx),
    error = function(e) NULL
  )
  if (is.null(xtx_svd) || !length(xtx_svd$d) || !any(is.finite(xtx_svd$d))) {
    return(list(status = "uncertified", reason = "surface_linear_solve_failed"))
  }
  sv_tol <- max(dim(xtx)) * max(xtx_svd$d, na.rm = TRUE) * .Machine$double.eps
  keep_sv <- xtx_svd$d > sv_tol
  if (!any(keep_sv)) {
    return(list(status = "uncertified", reason = "surface_linear_system_unidentified"))
  }
  beta <- as.numeric(
    xtx_svd$v[, keep_sv, drop = FALSE] %*%
      (as.numeric(crossprod(xtx_svd$u[, keep_sv, drop = FALSE], rhs)) / xtx_svd$d[keep_sv])
  )
  fitted <- as.numeric(X %*% beta)
  residual <- y - fitted
  rss <- sum((sqrt_w * residual)^2)
  rank <- sum(keep_sv)
  df <- max(nrow(X) - rank, 1L)
  residual_sd <- sqrt(rss / df)
  xtx_inv <- xtx_svd$v[, keep_sv, drop = FALSE] %*% (
    diag(1 / xtx_svd$d[keep_sv], nrow = sum(keep_sv)) %*%
      t(xtx_svd$u[, keep_sv, drop = FALSE])
  )
  se_model <- sqrt(max(xtx_inv[1L, 1L], 0)) * residual_sd

  list(
    log_marginal = as.numeric(beta[1L]),
    se = as.numeric(max(se_model, as.numeric(se_floor))),
    status = "certified",
    reason = "derivative_quadratic_surface",
    weights = rows$chart_id,
    beta = beta,
    coord_scale = coord_scale,
    residual_sd = as.numeric(residual_sd),
    n_observations = nrow(X),
    n_charts = length(unique(obs_chart)),
    surface_rank = as.integer(rank),
    surface_parameters = as.integer(p),
    df = as.integer(df),
    obs_type = obs_type,
    obs_chart = obs_chart
  )
}

.local_atlas_surface_leave_one_out <- function(theta,
                                               charts,
                                               rows,
                                               population_model,
                                               distance_scale,
                                               se_floor,
                                               value_nugget,
                                               gradient_weight,
                                               curvature_weight,
                                               ridge,
                                               min_charts) {
  if (nrow(rows) <= as.integer(min_charts)) {
    return(numeric(0))
  }
  vapply(rows$chart_id, function(chart_id) {
    fit <- .local_atlas_derivative_surface(
      theta = theta,
      charts = charts,
      rows = rows,
      population_model = population_model,
      distance_scale = distance_scale,
      se_floor = se_floor,
      value_nugget = value_nugget,
      gradient_weight = gradient_weight,
      curvature_weight = curvature_weight,
      ridge = ridge,
      min_charts = min_charts,
      omit_chart_id = chart_id
    )
    as.numeric(fit$log_marginal %||% NA_real_)
  }, numeric(1))
}

.local_atlas_particle_mis <- function(theta,
                                      charts,
                                      population_model,
                                      eta = NULL,
                                      min_ess_frac = 0.02,
                                      min_ess = 50,
                                      max_psis_k = 0.7) {
  model <- normalize_population_model(population_model)
  theta <- .as_hyper_matrix(theta, model$hyper_names, model$hyper_dim)
  if (!length(charts)) {
    return(list(status = "uncertified", reason = "no_particle_charts"))
  }
  charts <- lapply(charts, validate_local_chart, population_model = model)
  B <- length(charts)
  eta <- pmax(as.numeric(eta %||% rep(1 / B, B)), 0)
  if (length(eta) != B || sum(eta) <= 0 || any(!is.finite(eta))) {
    eta <- rep(1 / B, B)
  } else {
    eta <- eta / sum(eta)
  }
  log_eta <- log(pmax(eta, .Machine$double.eps))
  logZ <- vapply(charts, function(chart) as.numeric(chart$logZ_abs), numeric(1))
  chart_se <- vapply(charts, function(chart) as.numeric(chart$logZ_abs_se), numeric(1))
  if (any(!is.finite(logZ))) {
    return(list(status = "uncertified", reason = "nonfinite_chart_normalizer"))
  }

  log_terms <- numeric(0)
  log_ratios <- numeric(0)
  total_particles <- 0L
  for (s in seq_along(charts)) {
    chart_s <- charts[[s]]
    alpha <- chart_s$alpha_particles
    w <- .local_chart_normalize_weights(chart_s$weights, nrow(alpha))
    log_num <- population_model_log_alpha_given_theta(model, alpha = alpha, theta = theta)
    log_den <- matrix(NA_real_, nrow = nrow(alpha), ncol = B)
    for (b in seq_along(charts)) {
      log_den[, b] <- log_eta[b] +
        population_model_log_alpha_given_theta(model, alpha = alpha, theta = charts[[b]]$theta_anchor) -
        logZ[b]
    }
    log_h <- log_num - .rowLogSumExp(log_den)
    log_terms <- c(log_terms, log_eta[s] + log(pmax(w, .Machine$double.eps)) + log_h)
    log_ratios <- c(log_ratios, log_h)
    total_particles <- total_particles + nrow(alpha)
  }
  finite <- is.finite(log_terms)
  if (!any(finite)) {
    return(list(status = "uncertified", reason = "nonfinite_particle_mis_terms"))
  }
  log_m <- logsumexp(log_terms[finite])
  contribution_w <- exp(log_terms[finite] - log_m)
  ess <- 1 / sum(contribution_w * contribution_w)
  ess_frac <- ess / max(total_particles, 1L)
  psis_k <- .local_chart_psis_k(log_ratios[finite])
  normalizer_se <- sqrt(sum((eta * pmax(chart_se, 0))^2))
  mcse <- 1 / sqrt(max(ess, 1))
  ess_frac_ok <- is.finite(ess_frac) && ess_frac >= as.numeric(min_ess_frac)
  ess_abs_ok <- is.finite(ess) && is.finite(min_ess) && ess >= as.numeric(min_ess)
  ess_ok <- ess_frac_ok || ess_abs_ok
  psis_ok <- !is.finite(max_psis_k) || (!is.na(psis_k) && psis_k <= as.numeric(max_psis_k))
  status <- if (ess_ok && psis_ok) {
    "certified"
  } else {
    "uncertified"
  }
  reason <- if (identical(status, "certified")) {
    "particle_mis_overlap"
  } else if (!ess_ok) {
    "low_particle_mis_ess"
  } else {
    "high_particle_mis_psis"
  }
  list(
    log_marginal = as.numeric(log_m),
    se = as.numeric(sqrt(normalizer_se^2 + mcse^2)),
    status = status,
    reason = reason,
    ess = as.numeric(ess),
    ess_frac = as.numeric(ess_frac),
    psis_k = as.numeric(psis_k),
    eta = eta
  )
}

.local_atlas_particle_mis_many_global <- function(atlas,
                                                  theta,
                                                  population_model,
                                                  max_chart_distance = Inf,
                                                  min_covering_charts = 3L,
                                                  min_ess_frac = 0.05,
                                                  min_ess = 50,
                                                  max_psis_k = 0.7,
                                                  sparse_chart_min_covering = 3L,
                                                  sparse_chart_max_distance = Inf,
                                                  distance_metric = "euclidean",
                                                  se_floor = 1e-6,
                                                  use_compressed_particle_mis = FALSE,
                                                  require_compressed_particle_mis = FALSE,
                                                  use_uncertified_estimates = FALSE) {
  model <- normalize_population_model(population_model)
  atlas <- validate_local_atlas(atlas)
  theta <- .as_hyper_matrix(theta, model$hyper_names, model$hyper_dim)
  active_charts <- .local_atlas_active_charts(atlas)
  if (!length(active_charts) || is.null(atlas$normalizer_solution)) {
    return(data.frame(
      log_marginal = rep(NA_real_, nrow(theta)),
      se = rep(Inf, nrow(theta)),
      status = rep("uncertified", nrow(theta)),
      reason = rep(if (!length(active_charts)) "no_active_charts" else "missing_normalizer_solution", nrow(theta)),
      nearest_charts = rep("", nrow(theta)),
      check.names = FALSE
    ))
  }

  chart_ids <- names(active_charts)
  anchors <- do.call(rbind, lapply(active_charts, function(chart) chart$theta_anchor))
  colnames(anchors) <- model$hyper_names
  distances <- matrix(NA_real_, nrow = nrow(theta), ncol = length(active_charts))
  distances <- .local_atlas_theta_distances(
    theta = theta,
    centers = anchors,
    population_model = model,
    distance_metric = distance_metric
  )
  colnames(distances) <- chart_ids
  covering <- distances <= as.numeric(max_chart_distance)
  covering_count <- rowSums(covering)
  min_distance <- apply(distances, 1L, min)
  nearest <- max.col(-distances, ties.method = "first")
  exact <- min_distance <= 1e-10

  eta <- rep(1 / length(active_charts), length(active_charts))
  log_eta <- log(eta)
  logZ <- vapply(active_charts, function(chart) as.numeric(chart$logZ_abs), numeric(1))
  chart_se <- vapply(active_charts, function(chart) as.numeric(chart$logZ_abs_se), numeric(1))
  compression <- atlas$particle_mis_compression
  compression_ok <- isTRUE(use_compressed_particle_mis) &&
    inherits(compression, "local_atlas_particle_mis_compression") &&
    identical(as.character(compression$chart_ids), as.character(chart_ids))
  if (isTRUE(require_compressed_particle_mis) && !compression_ok) {
    return(data.frame(
      log_marginal = rep(NA_real_, nrow(theta)),
      se = rep(Inf, nrow(theta)),
      status = rep("uncertified", nrow(theta)),
      reason = rep("missing_or_stale_particle_mis_compression", nrow(theta)),
      nearest_charts = rep("", nrow(theta)),
      particle_mis_ess_frac = rep(0, nrow(theta)),
      particle_mis_ess = rep(0, nrow(theta)),
      particle_mis_psis_k = rep(NA_real_, nrow(theta)),
      min_covering_distance = min_distance,
      check.names = FALSE
    ))
  }
  if (compression_ok) {
    alpha <- .local_chart_align_alpha(compression$alpha, model)
    sample_log_weight <- log(pmax(.local_chart_normalize_weights(compression$weights, nrow(alpha)), .Machine$double.eps))
    particle_count_denominator <- nrow(alpha)
  } else {
    alpha <- do.call(rbind, lapply(active_charts, function(chart) chart$alpha_particles))
    sample_log_weight <- unlist(Map(function(chart, eta_s) {
      log(eta_s) + log(pmax(.local_chart_normalize_weights(chart$weights, nrow(chart$alpha_particles)), .Machine$double.eps))
    }, active_charts, eta), use.names = FALSE)
    particle_count_denominator <- nrow(alpha)
  }

  anchor_logp <- population_model_log_alpha_given_theta_many(model, alpha = alpha, theta = anchors)
  den_terms <- sweep(anchor_logp, 1L, log_eta - logZ, "+")
  log_den <- as.numeric(matrixStats::colLogSumExps(den_terms))
  theta_logp <- population_model_log_alpha_given_theta_many(model, alpha = alpha, theta = theta)
  log_ratios <- sweep(theta_logp, 2L, log_den, "-")
  log_terms <- sweep(log_ratios, 2L, sample_log_weight, "+")
  log_marginal <- .rowLogSumExp(log_terms)

  max_terms <- matrixStats::rowMaxs(log_terms)
  contribution <- exp(sweep(log_terms, 1L, max_terms, "-"))
  contribution <- sweep(contribution, 1L, rowSums(contribution), "/")
  ess <- 1 / rowSums(contribution * contribution)
  ess_frac <- ess / max(particle_count_denominator, 1L)
  psis_k <- vapply(seq_len(nrow(log_terms)), function(j) {
    .local_chart_psis_k(log_ratios[j, ])
  }, numeric(1))
  normalizer_se <- sqrt(sum((eta * pmax(chart_se, 0))^2))
  se <- sqrt(normalizer_se^2 + 1 / pmax(ess, 1))
  se <- pmax(se, as.numeric(se_floor))
  if (any(exact)) {
    exact_chart <- nearest[exact]
    exact_idx <- which(exact)
    exact_certified <- vapply(active_charts[exact_chart], .local_chart_normalizer_certified, logical(1))
    if (any(exact_certified)) {
      keep_idx <- exact_idx[exact_certified]
      keep_chart <- exact_chart[exact_certified]
      log_marginal[keep_idx] <- logZ[keep_chart]
      se[keep_idx] <- pmax(chart_se[keep_chart], as.numeric(se_floor))
      ess_frac[keep_idx] <- 1
      psis_k[keep_idx] <- -Inf
    }
  }

  compression_certified <- compression_ok &&
    isTRUE(as.logical(compression$diagnostics$certified[1L] %||% FALSE))
  psis_ok <- if (compression_ok) {
    rep(TRUE, nrow(theta))
  } else {
    !is.finite(max_psis_k) | (!is.na(psis_k) & psis_k <= as.numeric(max_psis_k))
  }
  sparse_ok <- !is.finite(sparse_chart_max_distance) |
    covering_count >= as.integer(sparse_chart_min_covering) |
    min_distance <= as.numeric(sparse_chart_max_distance)
  exact_normalizer_certified <- rep(FALSE, nrow(theta))
  if (any(exact)) {
    exact_chart <- nearest[exact]
    exact_normalizer_certified[exact] <- vapply(active_charts[exact_chart], .local_chart_normalizer_certified, logical(1))
  }
  if (compression_ok) {
    certified <- exact_normalizer_certified | (
      compression_certified &
        covering_count >= as.integer(min_covering_charts) &
        is.finite(log_marginal) &
        sparse_ok
    )
  } else {
    certified <- exact_normalizer_certified | (
      covering_count >= as.integer(min_covering_charts) &
        is.finite(log_marginal) &
        is.finite(ess_frac) &
        (ess_frac >= as.numeric(min_ess_frac) |
           (is.finite(min_ess) & is.finite(ess) & ess >= as.numeric(min_ess))) &
        psis_ok &
        sparse_ok
    )
  }
  reason <- rep("high_particle_mis_psis", nrow(theta))
  reason[!is.finite(ess_frac) | ess_frac < as.numeric(min_ess_frac)] <- "low_particle_mis_ess"
  reason[!sparse_ok] <- "sparse_chart_extrapolation"
  reason[covering_count < as.integer(min_covering_charts)] <- "no_active_chart_coverage"
  reason[exact & !exact_normalizer_certified & covering_count < as.integer(min_covering_charts)] <-
    "exact_anchor_normalizer_uncertified"
  if (compression_ok && !compression_certified) {
    reason[] <- "compressed_particle_mis_not_certified"
  }
  reason[certified] <- "particle_mis_global_active_chart_coverage"
  reason[certified & exact_normalizer_certified] <- "active_exact_anchor"
  nearest_charts <- vapply(seq_len(nrow(theta)), function(i) {
    covered <- chart_ids[covering[i, ]]
    if (length(covered)) paste(covered, collapse = ",") else chart_ids[nearest[i]]
  }, character(1))
  returned <- certified | (isTRUE(use_uncertified_estimates) & is.finite(log_marginal))
  data.frame(
    log_marginal = ifelse(returned, log_marginal, NA_real_),
    se = ifelse(returned, se, Inf),
    status = ifelse(certified, "certified", "uncertified"),
    reason = ifelse(
      certified & compression_ok & reason == "particle_mis_global_active_chart_coverage",
      "compressed_particle_mis_global_active_chart_coverage",
      reason
    ),
    nearest_charts = nearest_charts,
    particle_mis_ess_frac = ess_frac,
    particle_mis_ess = ess,
    particle_mis_psis_k = psis_k,
    min_covering_distance = min_distance,
    check.names = FALSE
  )
}

.local_atlas_particle_mis_cache <- function(atlas,
                                            theta,
                                            population_model) {
  atlas <- validate_local_atlas(atlas)
  model <- normalize_population_model(population_model)
  theta <- .as_hyper_matrix(theta, model$hyper_names, model$hyper_dim)
  active <- .local_atlas_active_charts(atlas)
  if (!length(active)) {
    stop("Cannot compress a local atlas with no active charts.")
  }
  chart_ids <- names(active)
  eta <- rep(1 / length(active), length(active))
  log_eta <- log(eta)
  anchors <- do.call(rbind, lapply(active, function(chart) chart$theta_anchor))
  colnames(anchors) <- model$hyper_names
  logZ <- vapply(active, function(chart) as.numeric(chart$logZ_abs), numeric(1))
  alpha <- do.call(rbind, lapply(active, function(chart) chart$alpha_particles))
  alpha <- .local_chart_align_alpha(alpha, model)
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
    chart_ids = chart_ids,
    anchors = anchors,
    logZ = logZ,
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

.local_atlas_compression_fit_simplex <- function(features,
                                                 target,
                                                 ridge = 1e-8) {
  features <- as.matrix(features)
  target <- as.numeric(target)
  K <- nrow(features)
  if (K <= 0L) {
    stop("Compression simplex fit requires at least one selected particle.")
  }
  if (K == 1L) {
    return(1)
  }
  if (!requireNamespace("quadprog", quietly = TRUE)) {
    stop("quadprog is required for local atlas particle compression.")
  }
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
  total <- sum(sol)
  if (!is.finite(total) || total <= 0) {
    rep(1 / K, K)
  } else {
    sol / total
  }
}

.local_atlas_compression_sparse_quadrature <- function(features,
                                                       target,
                                                       K,
                                                       ridge = 1e-8,
                                                       min_weight = 1e-12) {
  features <- as.matrix(features)
  target <- as.numeric(target)
  K <- min(as.integer(K), nrow(features))
  if (K <= 0L) {
    stop("Compression K must be positive.")
  }
  selected <- integer(0)
  active_weight <- numeric(0)
  current <- rep(0, ncol(features))
  available <- rep(TRUE, nrow(features))
  for (k in seq_len(K)) {
    residual <- target - current
    score <- as.numeric(features %*% residual)
    score[!available] <- -Inf
    pick <- which.max(score)
    if (!length(pick) || !is.finite(score[pick])) {
      break
    }
    selected <- c(selected, pick)
    available[pick] <- FALSE
    active_features <- features[selected, , drop = FALSE]
    active_weight <- .local_atlas_compression_fit_simplex(active_features, target, ridge = ridge)
    keep <- active_weight > as.numeric(min_weight)
    if (!any(keep)) {
      active_weight <- rep(1 / length(selected), length(selected))
      keep <- rep(TRUE, length(selected))
    }
    selected <- selected[keep]
    active_weight <- active_weight[keep]
    active_weight <- active_weight / sum(active_weight)
    current <- as.numeric(crossprod(active_weight, features[selected, , drop = FALSE]))
  }
  ord <- order(selected)
  list(index = selected[ord], weight = active_weight[ord])
}

.local_atlas_compression_features <- function(cache,
                                              train_rows,
                                              evidence_weight = 1,
                                              moment_weight = 0.05,
                                              chart_weight = 0.05,
                                              include_moments = TRUE,
                                              include_chart = TRUE) {
  train_rows <- as.integer(train_rows)
  evidence <- t(exp(sweep(cache$log_h[train_rows, , drop = FALSE], 1L, cache$full_log_m[train_rows], "-")))
  raw <- evidence
  target <- rep(1, length(train_rows))
  scale <- rep(1, length(train_rows))
  penalty_weight <- rep(as.numeric(evidence_weight), length(train_rows))

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

.local_atlas_compression_evaluate <- function(log_h,
                                              index,
                                              weight) {
  log_weight <- log(pmax(as.numeric(weight), .Machine$double.eps))
  log_terms <- sweep(log_h[, index, drop = FALSE], 2L, log_weight, "+")
  log_m <- .rowLogSumExp(log_terms)
  max_terms <- matrixStats::rowMaxs(log_terms)
  contribution <- exp(sweep(log_terms, 1L, max_terms, "-"))
  contribution <- sweep(contribution, 1L, rowSums(contribution), "/")
  ess <- 1 / rowSums(contribution * contribution)
  psis_k <- vapply(seq_len(nrow(log_h)), function(j) .local_chart_psis_k(log_h[j, index]), numeric(1))
  list(
    log_marginal = as.numeric(log_m),
    ess = as.numeric(ess),
    ess_frac = as.numeric(ess / length(index)),
    psis_k = as.numeric(psis_k)
  )
}

compress_local_atlas_particle_mis <- function(atlas,
                                              theta,
                                              population_model,
                                              K = 64L,
                                              holdout_fraction = 0.25,
                                              ridge = 1e-8,
                                              evidence_weight = 1,
                                              moment_weight = 0.05,
                                              chart_weight = 0.05,
                                              include_moments = TRUE,
                                              include_chart = TRUE,
                                              max_holdout_rmse = Inf,
                                              seed = NULL) {
  model <- normalize_population_model(population_model)
  atlas <- validate_local_atlas(atlas)
  theta <- .as_hyper_matrix(theta, model$hyper_names, model$hyper_dim)
  if (nrow(theta) < 2L) {
    stop("Particle-MIS compression needs at least two theta certification points.")
  }
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    } else {
      NULL
    }
    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(as.integer(seed))
  }
  cache <- .local_atlas_particle_mis_cache(atlas, theta, model)
  n_theta <- nrow(theta)
  holdout_n <- floor(n_theta * max(0, min(0.8, as.numeric(holdout_fraction))))
  if (holdout_n > 0L && n_theta - holdout_n >= 2L) {
    holdout_rows <- sort(sample.int(n_theta, holdout_n))
    train_rows <- setdiff(seq_len(n_theta), holdout_rows)
  } else {
    holdout_rows <- integer(0)
    train_rows <- seq_len(n_theta)
  }
  features <- .local_atlas_compression_features(
    cache = cache,
    train_rows = train_rows,
    evidence_weight = evidence_weight,
    moment_weight = moment_weight,
    chart_weight = chart_weight,
    include_moments = include_moments,
    include_chart = include_chart
  )
  selection <- .local_atlas_compression_sparse_quadrature(
    features = features$features,
    target = features$target,
    K = as.integer(K),
    ridge = ridge
  )
  compressed_eval <- .local_atlas_compression_evaluate(cache$log_h, selection$index, selection$weight)
  error <- compressed_eval$log_marginal - cache$full_log_m
  rmse <- function(x) {
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    if (length(x)) sqrt(mean(x * x)) else NA_real_
  }
  centered_rmse <- function(x) {
    ok <- is.finite(x)
    if (!any(ok)) return(NA_real_)
    y <- x[ok] - mean(x[ok])
    sqrt(mean(y * y))
  }
  diagnostics <- data.frame(
    local = as.character(atlas$local_id),
    raw_particles = nrow(cache$alpha),
    selected_particles = length(selection$index),
    quadrature_weight_ess = 1 / sum(selection$weight * selection$weight),
    compression_ratio = nrow(cache$alpha) / max(length(selection$index), 1L),
    theta_points = n_theta,
    train_points = length(train_rows),
    holdout_points = length(holdout_rows),
    train_rmse = rmse(error[train_rows]),
    holdout_rmse = if (length(holdout_rows)) rmse(error[holdout_rows]) else NA_real_,
    all_rmse = rmse(error),
    all_centered_rmse = centered_rmse(error),
    max_abs_error = if (any(is.finite(error))) max(abs(error[is.finite(error)])) else NA_real_,
    median_full_ess_frac = stats::median(cache$full_ess_frac, na.rm = TRUE),
    median_compressed_ess_frac = stats::median(compressed_eval$ess_frac, na.rm = TRUE),
    median_full_psis_k = stats::median(cache$full_psis_k, na.rm = TRUE),
    median_compressed_psis_k = stats::median(compressed_eval$psis_k, na.rm = TRUE),
    certified = if (length(holdout_rows) && is.finite(max_holdout_rmse)) {
      is.finite(rmse(error[holdout_rows])) && rmse(error[holdout_rows]) <= as.numeric(max_holdout_rmse)
    } else {
      TRUE
    },
    check.names = FALSE
  )
  compression <- structure(
    list(
      method = "function_aware_sparse_quadrature",
      chart_ids = cache$chart_ids,
      alpha = cache$alpha[selection$index, , drop = FALSE],
      weights = .local_chart_normalize_weights(selection$weight, length(selection$weight)),
      selected_index = selection$index,
      theta = theta,
      train_rows = train_rows,
      holdout_rows = holdout_rows,
      diagnostics = diagnostics
    ),
    class = "local_atlas_particle_mis_compression"
  )
  atlas$particle_mis_compression <- compression
  atlas
}

compress_local_atlas_factor_set <- function(factor_set,
                                            theta,
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
                                            n_cores = 1L,
                                            seed = 123L,
                                            verbose = TRUE) {
  factor_set <- validate_local_atlas_factor_set(factor_set)
  model <- factor_set$population_model
  theta <- .as_hyper_matrix(theta, model$hyper_names, model$hyper_dim)
  if (!requireNamespace("quadprog", quietly = TRUE)) {
    stop("quadprog is required for local atlas particle compression.")
  }
  compress_one <- function(pos) {
    compress_local_atlas_particle_mis(
      atlas = factor_set$atlases[[pos]],
      theta = theta,
      population_model = model,
      K = K,
      holdout_fraction = holdout_fraction,
      ridge = ridge,
      evidence_weight = evidence_weight,
      moment_weight = moment_weight,
      chart_weight = chart_weight,
      include_moments = include_moments,
      include_chart = include_chart,
      max_holdout_rmse = max_holdout_rmse,
      seed = as.integer(seed) + 7919L * pos
    )
  }
  if (isTRUE(verbose)) {
    .local_atlas_log(
      "compressing local particle-MIS rules: locals=", length(factor_set$atlases),
      " K=", as.integer(K),
      " theta_points=", nrow(theta), "\n",
      verbose = TRUE
    )
  }
  atlases <- if (as.integer(n_cores) <= 1L || length(factor_set$atlases) <= 1L) {
    lapply(seq_along(factor_set$atlases), compress_one)
  } else {
    parallel::mclapply(
      seq_along(factor_set$atlases),
      compress_one,
      mc.cores = as.integer(min(n_cores, length(factor_set$atlases))),
      mc.preschedule = FALSE
    )
  }
  names(atlases) <- names(factor_set$atlases)
  factor_set$atlases <- atlases
  summary <- local_atlas_compression_summary(factor_set)
  factor_set$compression_summary <- summary
  factor_set$evaluator_control$use_compressed_particle_mis <- TRUE
  factor_set$evaluator_control$require_compressed_particle_mis <- TRUE
  certified <- as.logical(summary$certified)
  certified[is.na(certified)] <- FALSE
  if (isTRUE(stop_on_failure) && nrow(summary) && any(!certified)) {
    bad <- summary[!certified, , drop = FALSE]
    bad <- utils::head(bad[order(-bad$holdout_rmse), , drop = FALSE], 8L)
    stop(
      "Local atlas particle compression failed certification:\n",
      paste(sprintf("local=%s holdout_rmse=%s selected=%s raw=%s",
                    bad$local, signif(bad$holdout_rmse, 4), bad$selected_particles, bad$raw_particles),
            collapse = "\n")
    )
  }
  factor_set
}

local_atlas_compression_summary <- function(factor_set) {
  factor_set <- validate_local_atlas_factor_set(factor_set)
  rows <- lapply(factor_set$atlases, function(atlas) {
    compression <- atlas$particle_mis_compression
    if (!inherits(compression, "local_atlas_particle_mis_compression")) {
      return(data.frame(
        local = as.character(atlas$local_id),
        raw_particles = NA_integer_,
        selected_particles = NA_integer_,
        compression_ratio = NA_real_,
        certified = FALSE,
        check.names = FALSE
      ))
    }
    compression$diagnostics
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.local_atlas_compression_theta_design <- function(theta_design,
                                                  theta_cloud,
                                                  initial_proposal,
                                                  population_model,
                                                  theta_root,
                                                  n_theta = 96L,
                                                  distance_metric = "fisher",
                                                  seed = 123L) {
  model <- normalize_population_model(population_model)
  n_theta <- as.integer(n_theta)
  if (!is.finite(n_theta) || n_theta < 2L) {
    stop("compressed_particle_mis_n_theta must be at least 2.")
  }
  theta_design <- .as_hyper_matrix(theta_design, model$hyper_names, model$hyper_dim)
  theta_cloud <- .as_hyper_matrix(theta_cloud, model$hyper_names, model$hyper_dim)
  proposal_theta <- NULL
  proposal_n <- max(0L, n_theta - nrow(theta_design))
  if (!is.null(initial_proposal) && proposal_n > 0L) {
    proposal_theta <- theta_proposal_sample(
      initial_proposal,
      n = proposal_n,
      seed = as.integer(seed) + 1000003L
    )
  }
  pool <- .local_atlas_unique_theta(rbind(theta_design, theta_cloud, proposal_theta), model)
  if (nrow(pool) <= n_theta) {
    return(pool)
  }
  seed_design <- theta_design
  seed_keep <- min(nrow(seed_design), max(1L, floor(0.5 * n_theta)))
  if (nrow(seed_design) > seed_keep) {
    seed_design <- .local_atlas_metric_farthest_design(
      theta_cloud = seed_design,
      population_model = model,
      theta_root = theta_root,
      max_anchors = seed_keep,
      distance_metric = distance_metric
    )
  }
  .local_atlas_metric_farthest_design(
    theta_cloud = pool,
    population_model = model,
    theta_root = theta_root,
    max_anchors = n_theta,
    distance_metric = distance_metric,
    initial_design = seed_design
  )
}

evaluate_local_atlas <- function(atlas,
                                 theta,
                                 population_model,
                                 max_chart_distance = Inf,
                                 min_covering_charts = 3L,
                                 max_prediction_range = Inf,
                                 distance_scale = max_chart_distance,
                                 se_floor = 1e-6,
                                 use_particle_mis = TRUE,
                                 require_particle_mis = TRUE,
                                 min_particle_mis_ess = 0.05,
                                 min_particle_mis_ess_abs = 50,
                                 max_particle_mis_psis_k = 0.7,
                                 max_quadratic_particle_gap = 0.05,
                                 sparse_chart_min_covering = 3L,
                                 sparse_chart_max_distance = Inf,
                                 max_leave_chart_out_gap = Inf,
                                 distance_metric = "euclidean",
                                 surface_method = "derivative_ls",
                                 min_surface_charts = 2L,
                                 max_surface_se = Inf,
                                 surface_value_nugget = 0.05,
                                 surface_gradient_weight = 1.0,
                                 surface_curvature_weight = 0.2,
                                 surface_ridge = 1e-8,
                                 particle_mis_role = "estimator") {
  model <- normalize_population_model(population_model)
  atlas <- validate_local_atlas(atlas)
  theta <- .as_hyper_matrix(theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  surface_method <- .local_atlas_surface_method(surface_method)
  particle_mis_role <- .local_atlas_particle_mis_role(particle_mis_role)
  if (nrow(theta) != 1L) {
    stop("evaluate_local_atlas expects one theta row.")
  }
  active_charts <- .local_atlas_active_charts(atlas)
  if (!length(active_charts)) {
    return(list(
      log_marginal = NA_real_,
      se = Inf,
      status = "uncertified",
      reason = "no_active_charts",
      nearest_charts = character(),
      diagnostics = list()
    ))
  }
  if (is.null(atlas$normalizer_solution)) {
    return(list(
      log_marginal = NA_real_,
      se = Inf,
      status = "uncertified",
      reason = "missing_normalizer_solution",
      nearest_charts = names(active_charts),
      diagnostics = list()
    ))
  }

  rows <- .local_atlas_prediction_rows(
    atlas = atlas,
    theta = theta,
    population_model = model,
    active_charts = active_charts,
    max_chart_distance = max_chart_distance,
    distance_metric = distance_metric
  )
  rows <- rows[is.finite(rows$log_marginal) & is.finite(rows$distance), , drop = FALSE]
  if (!nrow(rows)) {
    return(list(
      log_marginal = NA_real_,
      se = Inf,
      status = "uncertified",
      reason = "no_finite_chart_prediction",
      nearest_charts = character(),
      diagnostics = list(predictions = rows)
    ))
  }

  rows <- rows[order(rows$distance, rows$chart_se), , drop = FALSE]
  exact <- rows[rows$distance <= 1e-10, , drop = FALSE]
  if (nrow(exact)) {
    exact <- exact[order(exact$chart_se), , drop = FALSE]
    exact$normalizer_certified <- vapply(
      active_charts[exact$chart_id],
      .local_chart_normalizer_certified,
      logical(1)
    )
    trusted_exact <- exact[exact$normalizer_certified, , drop = FALSE]
    if (nrow(trusted_exact)) {
      return(list(
        log_marginal = trusted_exact$log_marginal[1L],
        se = max(trusted_exact$chart_se[1L], as.numeric(se_floor)),
        status = "certified",
        reason = "active_exact_anchor",
        nearest_charts = trusted_exact$chart_id[1L],
        diagnostics = list(
          predictions = trusted_exact,
          prediction_range = 0,
          quadratic_particle_gap = 0,
          leave_chart_out_gap = 0,
          min_covering_distance = 0,
          particle_mis = list(
            status = "certified",
            reason = "active_exact_anchor",
            ess = nrow(active_charts[[trusted_exact$chart_id[1L]]]$alpha_particles),
            ess_frac = 1,
            psis_k = -Inf
          ),
          distance_scale = as.numeric(distance_scale)
        )
      ))
    }
  }
  covering <- rows[rows$within_distance, , drop = FALSE]
  if (nrow(covering) < as.integer(min_covering_charts)) {
    return(list(
      log_marginal = NA_real_,
      se = Inf,
      status = "uncertified",
      reason = if (nrow(exact)) "exact_anchor_normalizer_uncertified" else "no_active_chart_coverage",
      nearest_charts = head(rows$chart_id, max(1L, as.integer(min_covering_charts))),
      diagnostics = list(predictions = rows)
    ))
  }
  min_covering_distance <- min(covering$distance[is.finite(covering$distance)])
  if (is.finite(sparse_chart_max_distance) &&
      nrow(covering) < as.integer(sparse_chart_min_covering) &&
      (!is.finite(min_covering_distance) || min_covering_distance > as.numeric(sparse_chart_max_distance))) {
    return(list(
      log_marginal = NA_real_,
      se = Inf,
      status = "uncertified",
      reason = "sparse_chart_extrapolation",
      nearest_charts = covering$chart_id,
      diagnostics = list(
        predictions = covering,
        min_covering_distance = as.numeric(min_covering_distance)
      )
    ))
  }

  prediction_range <- diff(range(covering$log_marginal))
  if (is.finite(max_prediction_range) && prediction_range > as.numeric(max_prediction_range)) {
    return(list(
      log_marginal = NA_real_,
      se = Inf,
      status = "uncertified",
      reason = "chart_prediction_disagreement",
      nearest_charts = covering$chart_id,
      diagnostics = list(
        predictions = covering,
        prediction_range = as.numeric(prediction_range)
      )
    ))
  }

  scale <- if (is.finite(distance_scale)) {
    distance_scale
  } else if (is.finite(max_chart_distance)) {
    max_chart_distance
  } else {
    max(stats::median(covering$distance), 1)
  }
  eta_combined <- .local_atlas_combine_predictions(
    rows = covering,
    distance_scale = scale,
    se_floor = se_floor
  )
  surface <- NULL
  leave_chart_out <- numeric(0)
  if (identical(surface_method, "derivative_ls")) {
    surface <- .local_atlas_derivative_surface(
      theta = theta,
      charts = active_charts,
      rows = covering,
      population_model = model,
      distance_scale = scale,
      se_floor = se_floor,
      value_nugget = surface_value_nugget,
      gradient_weight = surface_gradient_weight,
      curvature_weight = surface_curvature_weight,
      ridge = surface_ridge,
      min_charts = min_surface_charts
    )
    if (identical(surface$status, "certified")) {
      leave_chart_out <- .local_atlas_surface_leave_one_out(
        theta = theta,
        charts = active_charts,
        rows = covering,
        population_model = model,
        distance_scale = scale,
        se_floor = se_floor,
        value_nugget = surface_value_nugget,
        gradient_weight = surface_gradient_weight,
        curvature_weight = surface_curvature_weight,
        ridge = surface_ridge,
        min_charts = min_surface_charts
      )
    }
  } else {
    surface <- eta_combined
    surface$status <- "certified"
    surface$reason <- "chart_average_surface"
    leave_chart_out <- if (nrow(covering) > 1L) {
      vapply(seq_len(nrow(covering)), function(k) {
        .local_atlas_combine_predictions(
          rows = covering[-k, , drop = FALSE],
          distance_scale = scale,
          se_floor = se_floor
        )$log_marginal
      }, numeric(1))
    } else {
      numeric(0)
    }
  }
  leave_chart_out_gap <- if (length(leave_chart_out)) {
    max(abs(leave_chart_out - surface$log_marginal), na.rm = TRUE)
  } else {
    0
  }
  surface_se <- as.numeric(surface$se %||% NA_real_)
  covering$combination_weight <- eta_combined$weights
  if (!isTRUE(use_particle_mis)) {
    return(list(
      log_marginal = NA_real_,
      se = Inf,
      status = "uncertified",
      reason = "particle_mis_required",
      nearest_charts = covering$chart_id,
      diagnostics = list(
        predictions = covering,
        prediction_range = as.numeric(prediction_range),
        surface_method = surface_method,
        surface_log_marginal = surface$log_marginal %||% NA_real_,
        surface_se = surface_se,
        surface = surface,
        leave_chart_out = leave_chart_out,
        leave_chart_out_gap = as.numeric(leave_chart_out_gap),
        min_covering_distance = as.numeric(min_covering_distance),
        particle_mis = NULL,
        distance_scale = as.numeric(scale)
      )
    ))
  }
  particle <- .local_atlas_particle_mis(
    theta = theta,
    charts = active_charts,
    population_model = model,
    eta = NULL,
    min_ess_frac = min_particle_mis_ess,
    min_ess = min_particle_mis_ess_abs,
    max_psis_k = max_particle_mis_psis_k
  )
  gap <- if (is.finite(particle$log_marginal %||% NA_real_) &&
             is.finite(surface$log_marginal %||% NA_real_)) {
    abs(particle$log_marginal - surface$log_marginal)
  } else {
    Inf
  }
  particle$surface_gap <- as.numeric(gap)
  diagnostics <- list(
    predictions = covering,
    prediction_range = as.numeric(prediction_range),
    surface_method = surface_method,
    surface_log_marginal = surface$log_marginal %||% NA_real_,
    surface_se = surface_se,
    surface_residual_sd = as.numeric(surface$residual_sd %||% NA_real_),
    surface_n_observations = as.integer(surface$n_observations %||% NA_integer_),
    surface_n_charts = as.integer(surface$n_charts %||% nrow(covering)),
    surface = surface,
    between_sd = surface$between_sd %||% NA_real_,
    within_se = surface$within_se %||% NA_real_,
    quadratic_particle_gap = as.numeric(gap),
    leave_chart_out = leave_chart_out,
    leave_chart_out_gap = as.numeric(leave_chart_out_gap),
    min_covering_distance = as.numeric(min_covering_distance),
    particle_mis = particle,
    distance_scale = as.numeric(scale)
  )
  if (!identical(particle$status, "certified")) {
    return(list(
      log_marginal = NA_real_,
      se = Inf,
      status = "uncertified",
      reason = paste("particle_mis_failed", particle$reason %||% "unknown", sep = ":"),
      nearest_charts = covering$chart_id,
      diagnostics = diagnostics
    ))
  }
  list(
    log_marginal = particle$log_marginal,
    se = pmax(particle$se, as.numeric(se_floor)),
    status = "certified",
    reason = "particle_mis_active_chart_coverage",
    nearest_charts = covering$chart_id,
    diagnostics = diagnostics
  )
}

evaluate_local_atlas_many <- function(atlas,
                                      theta,
                                      population_model,
                                      max_chart_distance = Inf,
                                      min_covering_charts = 3L,
                                      max_prediction_range = Inf,
                                      distance_scale = max_chart_distance,
                                      se_floor = 1e-6,
                                      use_particle_mis = TRUE,
                                      require_particle_mis = TRUE,
                                      min_particle_mis_ess = 0.05,
                                      min_particle_mis_ess_abs = 50,
                                      max_particle_mis_psis_k = 0.7,
                                      max_quadratic_particle_gap = 0.05,
                                      sparse_chart_min_covering = 3L,
                                      sparse_chart_max_distance = Inf,
                                      max_leave_chart_out_gap = Inf,
                                      distance_metric = "euclidean",
                                      surface_method = "derivative_ls",
                                      min_surface_charts = 2L,
                                      max_surface_se = Inf,
                                      surface_value_nugget = 0.05,
                                      surface_gradient_weight = 1.0,
                                      surface_curvature_weight = 0.2,
                                      surface_ridge = 1e-8,
                                      particle_mis_role = "estimator") {
  model <- normalize_population_model(population_model)
  theta <- .as_hyper_matrix(theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  results <- lapply(seq_len(nrow(theta)), function(i) {
    evaluate_local_atlas(
      atlas = atlas,
      theta = theta[i, , drop = FALSE],
      population_model = model,
      max_chart_distance = max_chart_distance,
      min_covering_charts = min_covering_charts,
      max_prediction_range = max_prediction_range,
      distance_scale = distance_scale,
      se_floor = se_floor,
      use_particle_mis = use_particle_mis,
      require_particle_mis = require_particle_mis,
      min_particle_mis_ess = min_particle_mis_ess,
      min_particle_mis_ess_abs = min_particle_mis_ess_abs,
      max_particle_mis_psis_k = max_particle_mis_psis_k,
      max_quadratic_particle_gap = max_quadratic_particle_gap,
      sparse_chart_min_covering = sparse_chart_min_covering,
      sparse_chart_max_distance = sparse_chart_max_distance,
      max_leave_chart_out_gap = max_leave_chart_out_gap,
      distance_metric = distance_metric,
      surface_method = surface_method,
      min_surface_charts = min_surface_charts,
      max_surface_se = max_surface_se,
      surface_value_nugget = surface_value_nugget,
      surface_gradient_weight = surface_gradient_weight,
      surface_curvature_weight = surface_curvature_weight,
      surface_ridge = surface_ridge,
      particle_mis_role = particle_mis_role
    )
  })
  data.frame(
    log_marginal = vapply(results, `[[`, numeric(1), "log_marginal"),
    se = vapply(results, `[[`, numeric(1), "se"),
    status = vapply(results, `[[`, character(1), "status"),
    reason = vapply(results, `[[`, character(1), "reason"),
    nearest_charts = vapply(
      results,
      function(x) paste(x$nearest_charts, collapse = ","),
      character(1)
    ),
    particle_mis_ess_frac = vapply(
      results,
      function(x) as.numeric(x$diagnostics$particle_mis$ess_frac %||% NA_real_),
      numeric(1)
    ),
    particle_mis_ess = vapply(
      results,
      function(x) as.numeric(x$diagnostics$particle_mis$ess %||% NA_real_),
      numeric(1)
    ),
    particle_mis_psis_k = vapply(
      results,
      function(x) as.numeric(x$diagnostics$particle_mis$psis_k %||% NA_real_),
      numeric(1)
    ),
    quadratic_particle_gap = vapply(
      results,
      function(x) as.numeric(x$diagnostics$quadratic_particle_gap %||% NA_real_),
      numeric(1)
    ),
    surface_method = vapply(
      results,
      function(x) as.character(x$diagnostics$surface_method %||% NA_character_),
      character(1)
    ),
    surface_se = vapply(
      results,
      function(x) as.numeric(x$diagnostics$surface_se %||% NA_real_),
      numeric(1)
    ),
    surface_residual_sd = vapply(
      results,
      function(x) as.numeric(x$diagnostics$surface_residual_sd %||% NA_real_),
      numeric(1)
    ),
    prediction_range = vapply(
      results,
      function(x) as.numeric(x$diagnostics$prediction_range %||% NA_real_),
      numeric(1)
    ),
    leave_chart_out_gap = vapply(
      results,
      function(x) as.numeric(x$diagnostics$leave_chart_out_gap %||% NA_real_),
      numeric(1)
    ),
    min_covering_distance = vapply(
      results,
      function(x) as.numeric(x$diagnostics$min_covering_distance %||% NA_real_),
      numeric(1)
    ),
    check.names = FALSE
  )
}

build_local_atlas_factor_set <- function(atlases,
                                         population_model,
                                         max_chart_distance = Inf,
                                         min_covering_charts = 3L,
                                         max_prediction_range = Inf,
                                         distance_scale = max_chart_distance,
                                         se_floor = 1e-6,
                                         use_particle_mis = TRUE,
                                         require_particle_mis = TRUE,
                                         min_particle_mis_ess = 0.05,
                                         min_particle_mis_ess_abs = 50,
                                         max_particle_mis_psis_k = 0.7,
                                         max_quadratic_particle_gap = 0.05,
                                         sparse_chart_min_covering = 3L,
                                         sparse_chart_max_distance = Inf,
                                         max_leave_chart_out_gap = Inf,
                                         distance_metric = "euclidean",
                                         surface_method = "derivative_ls",
                                         min_surface_charts = 2L,
                                         max_surface_se = Inf,
                                         surface_value_nugget = 0.05,
                                         surface_gradient_weight = 1.0,
                                         surface_curvature_weight = 0.2,
                                         surface_ridge = 1e-8,
                                         particle_mis_role = "estimator",
                                         particle_mis_batch = TRUE,
                                         use_compressed_particle_mis = FALSE,
                                         require_compressed_particle_mis = FALSE,
                                         stop_on_uncertified = TRUE,
                                         use_uncertified_estimates = FALSE) {
  model <- normalize_population_model(population_model)
  if (!is.list(atlases) || !length(atlases)) {
    stop("atlases must be a non-empty list.")
  }
  if (any(!vapply(atlases, inherits, logical(1), "local_evidence_atlas"))) {
    stop("all entries in atlases must inherit from 'local_evidence_atlas'.")
  }
  atlases <- lapply(atlases, validate_local_atlas)
  local_ids <- vapply(atlases, function(atlas) as.character(atlas$local_id), character(1))
  if (anyDuplicated(local_ids)) {
    stop("atlas local_id values must be unique.")
  }
  atlas_names <- names(atlases)
  if (is.null(atlas_names)) {
    atlas_names <- local_ids
  }
  missing_names <- !nzchar(atlas_names)
  atlas_names[missing_names] <- local_ids[missing_names]
  names(atlases) <- atlas_names
  missing_solution <- vapply(atlases, function(atlas) is.null(atlas$normalizer_solution), logical(1))
  if (any(missing_solution)) {
    stop("all atlases must have a graph normalizer solution before outer integration.")
  }

  structure(
    list(
      population_model = model,
      atlases = atlases,
      n_locals = length(atlases),
      log_constant = 0,
      evaluator_control = list(
        max_chart_distance = max_chart_distance,
        min_covering_charts = as.integer(min_covering_charts),
        max_prediction_range = max_prediction_range,
        distance_scale = distance_scale,
        se_floor = as.numeric(se_floor),
        use_particle_mis = isTRUE(use_particle_mis),
        require_particle_mis = isTRUE(require_particle_mis),
        min_particle_mis_ess = as.numeric(min_particle_mis_ess),
        min_particle_mis_ess_abs = as.numeric(min_particle_mis_ess_abs),
        max_particle_mis_psis_k = as.numeric(max_particle_mis_psis_k),
        max_quadratic_particle_gap = max_quadratic_particle_gap,
        sparse_chart_min_covering = as.integer(sparse_chart_min_covering),
        sparse_chart_max_distance = as.numeric(sparse_chart_max_distance),
        max_leave_chart_out_gap = as.numeric(max_leave_chart_out_gap),
        distance_metric = .local_atlas_distance_metric(distance_metric),
        surface_method = .local_atlas_surface_method(surface_method),
        min_surface_charts = as.integer(min_surface_charts),
        max_surface_se = as.numeric(max_surface_se),
        surface_value_nugget = as.numeric(surface_value_nugget),
        surface_gradient_weight = as.numeric(surface_gradient_weight),
        surface_curvature_weight = as.numeric(surface_curvature_weight),
        surface_ridge = as.numeric(surface_ridge),
        particle_mis_role = .local_atlas_particle_mis_role(particle_mis_role),
        particle_mis_batch = isTRUE(particle_mis_batch),
        use_compressed_particle_mis = isTRUE(use_compressed_particle_mis),
        require_compressed_particle_mis = isTRUE(require_compressed_particle_mis),
        stop_on_uncertified = isTRUE(stop_on_uncertified),
        use_uncertified_estimates = isTRUE(use_uncertified_estimates)
      )
    ),
    class = c("local_evidence_atlas_factor_set", "population_factor_set")
  )
}

validate_local_atlas_factor_set <- function(factor_set) {
  if (!inherits(factor_set, "local_evidence_atlas_factor_set")) {
    stop("factor_set must inherit from 'local_evidence_atlas_factor_set'.")
  }
  factor_set$population_model <- normalize_population_model(factor_set$population_model)
  if (!is.list(factor_set$atlases) || !length(factor_set$atlases)) {
    stop("atlas factor set must contain at least one atlas.")
  }
  factor_set$atlases <- lapply(factor_set$atlases, validate_local_atlas)
  factor_set$n_locals <- length(factor_set$atlases)
  factor_set$log_constant <- as.numeric(factor_set$log_constant %||% 0)
  if (length(factor_set$log_constant) != 1L || !is.finite(factor_set$log_constant)) {
    stop("atlas factor set log_constant must be a finite scalar.")
  }
  control <- factor_set$evaluator_control %||% list()
  control$max_chart_distance <- control$max_chart_distance %||% Inf
  control$min_covering_charts <- as.integer(control$min_covering_charts %||% 3L)
  control$max_prediction_range <- control$max_prediction_range %||% Inf
  control$distance_scale <- control$distance_scale %||% control$max_chart_distance
  control$se_floor <- as.numeric(control$se_floor %||% 1e-6)
  control$use_particle_mis <- isTRUE(control$use_particle_mis %||% TRUE)
  control$require_particle_mis <- isTRUE(control$require_particle_mis %||% TRUE)
  control$min_particle_mis_ess <- as.numeric(control$min_particle_mis_ess %||% 0.05)
  control$min_particle_mis_ess_abs <- as.numeric(control$min_particle_mis_ess_abs %||% 50)
  control$max_particle_mis_psis_k <- as.numeric(control$max_particle_mis_psis_k %||% 0.7)
  control$max_quadratic_particle_gap <- control$max_quadratic_particle_gap %||% 0.05
  control$sparse_chart_min_covering <- as.integer(control$sparse_chart_min_covering %||% 3L)
  control$sparse_chart_max_distance <- as.numeric(control$sparse_chart_max_distance %||% Inf)
  control$max_leave_chart_out_gap <- as.numeric(control$max_leave_chart_out_gap %||% Inf)
  control$distance_metric <- .local_atlas_distance_metric(control$distance_metric %||% "euclidean")
  control$surface_method <- .local_atlas_surface_method(control$surface_method %||% "derivative_ls")
  control$min_surface_charts <- as.integer(control$min_surface_charts %||% 2L)
  control$max_surface_se <- as.numeric(control$max_surface_se %||% Inf)
  control$surface_value_nugget <- as.numeric(control$surface_value_nugget %||% 0.05)
  control$surface_gradient_weight <- as.numeric(control$surface_gradient_weight %||% 1.0)
  control$surface_curvature_weight <- as.numeric(control$surface_curvature_weight %||% 0.2)
  control$surface_ridge <- as.numeric(control$surface_ridge %||% 1e-8)
  control$particle_mis_role <- .local_atlas_particle_mis_role(control$particle_mis_role %||% "estimator")
  control$particle_mis_batch <- isTRUE(control$particle_mis_batch %||% TRUE)
  control$use_compressed_particle_mis <- isTRUE(control$use_compressed_particle_mis %||% FALSE)
  control$require_compressed_particle_mis <- isTRUE(control$require_compressed_particle_mis %||% FALSE)
  control$stop_on_uncertified <- isTRUE(control$stop_on_uncertified %||% TRUE)
  control$use_uncertified_estimates <- isTRUE(control$use_uncertified_estimates %||% FALSE)
  factor_set$evaluator_control <- control
  factor_set
}

build_local_evidence_certification_cloud <- function(theta,
                                                     population_model,
                                                     theta_weights = NULL,
                                                     theta_source = "theta_cloud",
                                                     theta_round = 0L,
                                                     theta_id = NULL,
                                                     metadata = NULL,
                                                     deduplicate = TRUE) {
  model <- normalize_population_model(population_model)
  theta <- .as_hyper_matrix(theta, model$hyper_names, model$hyper_dim)
  n_theta <- nrow(theta)
  if (!n_theta) {
    stop("theta must contain at least one row.")
  }
  theta_weights <- .local_chart_normalize_weights(theta_weights, n_theta)
  theta_source <- rep(as.character(theta_source), length.out = n_theta)
  theta_source[is.na(theta_source) | !nzchar(theta_source)] <- "theta_cloud"
  theta_round <- rep(as.integer(theta_round), length.out = n_theta)
  theta_round[!is.finite(theta_round)] <- 0L
  if (is.null(theta_id)) {
    theta_id <- sprintf("cert_theta_%06d", seq_len(n_theta))
  }
  theta_id <- rep(as.character(theta_id), length.out = n_theta)
  theta_id[is.na(theta_id) | !nzchar(theta_id)] <- sprintf(
    "cert_theta_%06d",
    which(is.na(theta_id) | !nzchar(theta_id))
  )

  if (is.null(metadata)) {
    metadata <- data.frame(row.names = seq_len(n_theta))
  } else {
    metadata <- as.data.frame(metadata, stringsAsFactors = FALSE, check.names = FALSE)
    if (nrow(metadata) != n_theta) {
      stop("metadata must have one row per theta row.")
    }
    reserved <- c(
      "theta_row", "theta_id", "theta_source", "theta_weight", "theta_round",
      "multiplicity", "original_theta_ids", model$hyper_names
    )
    metadata <- metadata[, setdiff(names(metadata), reserved), drop = FALSE]
  }

  meta <- data.frame(
    theta_row = seq_len(n_theta),
    theta_id = theta_id,
    theta_source = theta_source,
    theta_weight = theta_weights,
    theta_round = theta_round,
    multiplicity = 1L,
    original_theta_ids = theta_id,
    check.names = FALSE
  )
  if (ncol(metadata)) {
    meta <- cbind(meta, metadata)
  }

  if (isTRUE(deduplicate)) {
    key <- apply(signif(theta, 14L), 1L, paste, collapse = "\r")
    groups <- split(seq_len(n_theta), key, drop = TRUE)
    theta <- do.call(rbind, lapply(groups, function(idx) theta[idx[1L], , drop = FALSE]))
    colnames(theta) <- model$hyper_names
    meta <- do.call(rbind, lapply(seq_along(groups), function(group_id) {
      idx <- groups[[group_id]]
      row <- meta[idx[1L], , drop = FALSE]
      row$theta_id <- if (length(idx) == 1L) {
        meta$theta_id[idx]
      } else {
        sprintf("cert_theta_%06d", group_id)
      }
      row$theta_source <- paste(sort(unique(meta$theta_source[idx])), collapse = ",")
      row$theta_weight <- sum(meta$theta_weight[idx])
      row$theta_round <- min(meta$theta_round[idx])
      row$multiplicity <- sum(meta$multiplicity[idx])
      row$original_theta_ids <- paste(meta$original_theta_ids[idx], collapse = ",")
      row
    }))
    rownames(meta) <- NULL
    meta$theta_row <- seq_len(nrow(theta))
    meta$theta_weight <- .local_chart_normalize_weights(meta$theta_weight, nrow(theta))
  }

  meta <- cbind(meta, as.data.frame(theta, check.names = FALSE))
  structure(
    list(theta = theta, metadata = meta),
    class = "local_evidence_certification_cloud",
    population_model = model
  ) |>
    validate_local_evidence_certification_cloud(population_model = model)
}

validate_local_evidence_certification_cloud <- function(cloud,
                                                        population_model = NULL) {
  if (!inherits(cloud, "local_evidence_certification_cloud")) {
    stop("cloud must inherit from 'local_evidence_certification_cloud'.")
  }
  model <- normalize_population_model(population_model %||% attr(cloud, "population_model"))
  cloud$theta <- .as_hyper_matrix(cloud$theta, model$hyper_names, model$hyper_dim)
  if (!is.data.frame(cloud$metadata) || nrow(cloud$metadata) != nrow(cloud$theta)) {
    stop("cloud metadata must be a data frame with one row per theta row.")
  }
  required <- c("theta_row", "theta_id", "theta_source", "theta_weight", "theta_round")
  missing <- setdiff(required, names(cloud$metadata))
  if (length(missing)) {
    stop("certification cloud metadata is missing: ", paste(missing, collapse = ", "))
  }
  cloud$metadata$theta_row <- as.integer(cloud$metadata$theta_row)
  if (!identical(cloud$metadata$theta_row, seq_len(nrow(cloud$theta)))) {
    cloud$metadata$theta_row <- seq_len(nrow(cloud$theta))
  }
  cloud$metadata$theta_id <- as.character(cloud$metadata$theta_id)
  if (anyDuplicated(cloud$metadata$theta_id)) {
    stop("certification cloud theta_id values must be unique.")
  }
  cloud$metadata$theta_source <- as.character(cloud$metadata$theta_source)
  cloud$metadata$theta_round <- as.integer(cloud$metadata$theta_round)
  cloud$metadata$theta_weight <- .local_chart_normalize_weights(
    cloud$metadata$theta_weight,
    nrow(cloud$theta)
  )
  cloud$metadata[, model$hyper_names] <- as.data.frame(cloud$theta, check.names = FALSE)
  attr(cloud, "population_model") <- model
  cloud
}

.local_evidence_certification_as_cloud <- function(cloud = NULL,
                                                   theta = NULL,
                                                   population_model,
                                                   theta_weights = NULL,
                                                   theta_source = "theta_cloud",
                                                   theta_round = 0L,
                                                   theta_id = NULL,
                                                   metadata = NULL) {
  model <- normalize_population_model(population_model)
  if (!is.null(cloud)) {
    return(validate_local_evidence_certification_cloud(cloud, model))
  }
  if (is.null(theta)) {
    stop("supply either cloud or theta.")
  }
  build_local_evidence_certification_cloud(
    theta = theta,
    population_model = model,
    theta_weights = theta_weights,
    theta_source = theta_source,
    theta_round = theta_round,
    theta_id = theta_id,
    metadata = metadata
  )
}

.local_evidence_eval_column <- function(x, name, default) {
  if (name %in% names(x)) {
    return(x[[name]])
  }
  rep(default, nrow(x))
}

.local_evidence_certification_nearest_geometry <- function(atlas,
                                                          theta,
                                                          population_model,
                                                          distance_metric) {
  model <- normalize_population_model(population_model)
  theta <- .as_hyper_matrix(theta, model$hyper_names, model$hyper_dim)
  atlas <- validate_local_atlas(atlas)
  active <- .local_atlas_active_charts(atlas)
  empty <- data.frame(
    nearest_chart = rep(NA_character_, nrow(theta)),
    local_score_norm = rep(NA_real_, nrow(theta)),
    local_curvature_norm = rep(NA_real_, nrow(theta)),
    check.names = FALSE
  )
  if (!length(active)) {
    return(empty)
  }
  anchors <- do.call(rbind, lapply(active, function(chart) chart$theta_anchor))
  colnames(anchors) <- model$hyper_names
  distances <- .local_atlas_theta_distances(
    theta = theta,
    centers = anchors,
    population_model = model,
    distance_metric = distance_metric
  )
  nearest <- max.col(-distances, ties.method = "first")
  chart_names <- names(active)
  score_norm <- vapply(active, function(chart) {
    if (is.null(chart$score)) NA_real_ else sqrt(sum(as.numeric(chart$score)^2))
  }, numeric(1))
  curvature_norm <- vapply(active, function(chart) {
    if (is.null(chart$curvature)) NA_real_ else sqrt(sum(as.matrix(chart$curvature)^2))
  }, numeric(1))
  data.frame(
    nearest_chart = chart_names[nearest],
    local_score_norm = score_norm[nearest],
    local_curvature_norm = curvature_norm[nearest],
    check.names = FALSE
  )
}

.local_evidence_score_certification_table <- function(table,
                                                      evaluator_control,
                                                      scoring_control = list()) {
  scoring_control <- modifyList(
    list(
      uncertified = 1,
      ess_frac = 1,
      ess_abs = 0.5,
      psis = 1,
      distance = 0.1,
      se = 0.1,
      score_norm = 0.1,
      curvature_norm = 0.02
    ),
    scoring_control
  )
  bad <- is.na(table$status) | table$status != "certified"
  min_ess_frac <- as.numeric(evaluator_control$min_particle_mis_ess %||% 0.05)
  min_ess_abs <- as.numeric(evaluator_control$min_particle_mis_ess_abs %||% 50)
  max_psis_k <- as.numeric(evaluator_control$max_particle_mis_psis_k %||% 0.7)

  ess_frac <- as.numeric(table$particle_mis_ess_frac)
  ess_frac_deficit <- if (is.finite(min_ess_frac) && min_ess_frac > 0) {
    ifelse(is.finite(ess_frac), pmax(min_ess_frac - ess_frac, 0) / min_ess_frac, 1)
  } else {
    rep(0, nrow(table))
  }
  ess_abs <- as.numeric(table$particle_mis_ess)
  ess_abs_deficit <- if (is.finite(min_ess_abs) && min_ess_abs > 0) {
    ifelse(is.finite(ess_abs), pmax(min_ess_abs - ess_abs, 0) / min_ess_abs, 1)
  } else {
    rep(0, nrow(table))
  }
  psis <- as.numeric(table$particle_mis_psis_k)
  psis_excess <- if (is.finite(max_psis_k)) {
    ifelse(is.finite(psis), pmax(psis - max_psis_k, 0), ifelse(bad, 1, 0))
  } else {
    rep(0, nrow(table))
  }
  distance_score <- ifelse(
    is.finite(table$min_covering_distance),
    log1p(pmax(as.numeric(table$min_covering_distance), 0)),
    ifelse(bad, 10, 0)
  )
  se_score <- ifelse(
    is.finite(table$raw_se),
    log1p(pmax(as.numeric(table$raw_se), 0)),
    ifelse(bad, 10, 0)
  )
  score_norm <- ifelse(is.finite(table$local_score_norm), pmax(table$local_score_norm, 0), 0)
  curvature_norm <- ifelse(is.finite(table$local_curvature_norm), pmax(table$local_curvature_norm, 0), 0)
  sensitivity <- 1 +
    as.numeric(scoring_control$score_norm) * log1p(score_norm) +
    as.numeric(scoring_control$curvature_norm) * log1p(curvature_norm)
  severity <-
    as.numeric(scoring_control$uncertified) * as.numeric(bad) +
    as.numeric(scoring_control$ess_frac) * ess_frac_deficit +
    as.numeric(scoring_control$ess_abs) * ess_abs_deficit +
    as.numeric(scoring_control$psis) * psis_excess +
    as.numeric(scoring_control$distance) * distance_score +
    as.numeric(scoring_control$se) * se_score
  severity[!is.finite(severity) | severity < 0] <- 0
  sensitivity[!is.finite(sensitivity) | sensitivity < 1] <- 1
  table$diagnostic_severity <- as.numeric(severity)
  table$sensitivity_score <- as.numeric(sensitivity)
  table$impact_score <- as.numeric(table$theta_weight) * table$diagnostic_severity * table$sensitivity_score
  table
}

evaluate_raw_local_evidence_certification <- function(factor_set,
                                                      cloud = NULL,
                                                      theta = NULL,
                                                      theta_weights = NULL,
                                                      theta_source = "theta_cloud",
                                                      theta_round = 0L,
                                                      theta_id = NULL,
                                                      metadata = NULL,
                                                      local_ids = NULL,
                                                      scoring_control = list(),
                                                      include_theta = FALSE,
                                                      n_cores = 1L) {
  factor_set <- validate_local_atlas_factor_set(factor_set)
  model <- factor_set$population_model
  cloud <- .local_evidence_certification_as_cloud(
    cloud = cloud,
    theta = theta,
    population_model = model,
    theta_weights = theta_weights,
    theta_source = theta_source,
    theta_round = theta_round,
    theta_id = theta_id,
    metadata = metadata
  )
  theta <- cloud$theta
  theta_meta <- cloud$metadata
  control <- factor_set$evaluator_control
  if (is.null(local_ids)) {
    local_pos <- seq_along(factor_set$atlases)
  } else if (is.character(local_ids)) {
    local_pos <- match(local_ids, names(factor_set$atlases))
  } else {
    local_pos <- as.integer(local_ids)
  }
  local_pos <- local_pos[is.finite(local_pos)]
  if (!length(local_pos) || any(local_pos < 1L | local_pos > length(factor_set$atlases))) {
    stop("local_ids must identify atlases in factor_set.")
  }
  local_pos <- unique(local_pos)

  eval_one <- function(pos) {
    atlas <- factor_set$atlases[[pos]]
    raw <- .local_atlas_particle_mis_many_global(
      atlas = atlas,
      theta = theta,
      population_model = model,
      max_chart_distance = control$max_chart_distance,
      min_covering_charts = control$min_covering_charts,
      min_ess_frac = control$min_particle_mis_ess,
      min_ess = control$min_particle_mis_ess_abs,
      max_psis_k = control$max_particle_mis_psis_k,
      sparse_chart_min_covering = control$sparse_chart_min_covering,
      sparse_chart_max_distance = control$sparse_chart_max_distance,
      distance_metric = control$distance_metric,
      se_floor = control$se_floor,
      use_compressed_particle_mis = FALSE,
      require_compressed_particle_mis = FALSE,
      use_uncertified_estimates = TRUE
    )
    geometry <- .local_evidence_certification_nearest_geometry(
      atlas = atlas,
      theta = theta,
      population_model = model,
      distance_metric = control$distance_metric
    )
    out <- data.frame(
      local = names(factor_set$atlases)[pos],
      local_pos = as.integer(pos),
      theta_row = theta_meta$theta_row,
      theta_id = theta_meta$theta_id,
      theta_source = theta_meta$theta_source,
      theta_weight = theta_meta$theta_weight,
      raw_log_m = .local_evidence_eval_column(raw, "log_marginal", NA_real_),
      raw_se = .local_evidence_eval_column(raw, "se", Inf),
      status = .local_evidence_eval_column(raw, "status", "uncertified"),
      reason = .local_evidence_eval_column(raw, "reason", "unknown"),
      particle_mis_ess = .local_evidence_eval_column(raw, "particle_mis_ess", NA_real_),
      particle_mis_ess_frac = .local_evidence_eval_column(raw, "particle_mis_ess_frac", NA_real_),
      particle_mis_psis_k = .local_evidence_eval_column(raw, "particle_mis_psis_k", NA_real_),
      min_covering_distance = .local_evidence_eval_column(raw, "min_covering_distance", NA_real_),
      nearest_charts = .local_evidence_eval_column(raw, "nearest_charts", ""),
      local_score_norm = geometry$local_score_norm,
      local_curvature_norm = geometry$local_curvature_norm,
      selected_for_repair = FALSE,
      repair_status = NA_character_,
      check.names = FALSE
    )
    if (isTRUE(include_theta)) {
      out <- cbind(out, as.data.frame(theta, check.names = FALSE))
    }
    out
  }

  rows <- if (as.integer(n_cores) <= 1L || length(local_pos) <= 1L) {
    lapply(local_pos, eval_one)
  } else {
    parallel::mclapply(
      local_pos,
      eval_one,
      mc.cores = as.integer(min(n_cores, length(local_pos)))
    )
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out <- .local_evidence_score_certification_table(out, control, scoring_control)
  class(out) <- c("local_evidence_certification_table", class(out))
  attr(out, "cloud") <- cloud
  attr(out, "settings") <- list(
    raw_particle_mis = TRUE,
    use_compressed_particle_mis = FALSE,
    evaluator_control = control,
    scoring_control = scoring_control,
    local_ids = names(factor_set$atlases)[local_pos]
  )
  out
}

.local_evidence_weighted_count_quantile <- function(x, weights, prob) {
  q <- .local_atlas_weighted_quantile(x, weights, prob)
  if (!length(q) || !is.finite(q[1L])) NA_real_ else as.numeric(q[1L])
}

summarize_raw_local_evidence_certification <- function(certification_table,
                                                       thresholds = list()) {
  if (!is.data.frame(certification_table) || !nrow(certification_table)) {
    stop("certification_table must be a non-empty data frame.")
  }
  required <- c(
    "local", "local_pos", "theta_row", "theta_id", "theta_source", "theta_weight",
    "status", "reason", "raw_se", "particle_mis_ess_frac", "particle_mis_psis_k",
    "impact_score"
  )
  missing <- setdiff(required, names(certification_table))
  if (length(missing)) {
    stop("certification_table is missing: ", paste(missing, collapse = ", "))
  }
  thresholds <- modifyList(
    list(
      max_weighted_any_uncertified = 0.05,
      max_weighted_mean_uncertified = 0.10,
      max_weighted_q90_uncertified = 0,
      high_weight_quantile = 0.95,
      max_high_weight_uncertified = 1,
      max_worst_local_uncertified_weight = 0.02
    ),
    thresholds
  )
  tab <- certification_table
  tab$uncertified <- is.na(tab$status) | tab$status != "certified"
  n_locals <- length(unique(tab$local_pos))
  n_theta <- length(unique(tab$theta_row))

  by_theta <- do.call(rbind, lapply(split(tab, tab$theta_row, drop = TRUE), function(df) {
    bad <- df$uncertified
    data.frame(
      theta_row = as.integer(df$theta_row[1L]),
      theta_id = as.character(df$theta_id[1L]),
      theta_source = as.character(df$theta_source[1L]),
      theta_weight = as.numeric(df$theta_weight[1L]),
      n_locals = nrow(df),
      uncertified_local_count = sum(bad),
      uncertified_fraction = mean(bad),
      total_impact_score = sum(df$impact_score, na.rm = TRUE),
      uncertified_impact_score = sum(df$impact_score[bad], na.rm = TRUE),
      min_particle_mis_ess_frac = suppressWarnings(min(df$particle_mis_ess_frac, na.rm = TRUE)),
      max_particle_mis_psis_k = suppressWarnings(max(df$particle_mis_psis_k, na.rm = TRUE)),
      max_raw_se = suppressWarnings(max(df$raw_se, na.rm = TRUE)),
      reasons = paste(sort(unique(df$reason[bad])), collapse = ","),
      check.names = FALSE
    )
  }))
  by_theta$min_particle_mis_ess_frac[!is.finite(by_theta$min_particle_mis_ess_frac)] <- NA_real_
  by_theta$max_particle_mis_psis_k[!is.finite(by_theta$max_particle_mis_psis_k)] <- NA_real_
  by_theta$max_raw_se[!is.finite(by_theta$max_raw_se)] <- Inf

  by_local <- do.call(rbind, lapply(split(tab, tab$local_pos, drop = TRUE), function(df) {
    bad <- df$uncertified
    data.frame(
      local = as.character(df$local[1L]),
      local_pos = as.integer(df$local_pos[1L]),
      n_theta = nrow(df),
      certified_fraction = mean(!bad),
      uncertified_fraction = mean(bad),
      uncertified_weight = sum(df$theta_weight[bad]),
      total_impact_score = sum(df$impact_score, na.rm = TRUE),
      uncertified_impact_score = sum(df$impact_score[bad], na.rm = TRUE),
      max_impact_score = suppressWarnings(max(df$impact_score, na.rm = TRUE)),
      min_particle_mis_ess_frac = suppressWarnings(min(df$particle_mis_ess_frac, na.rm = TRUE)),
      max_particle_mis_psis_k = suppressWarnings(max(df$particle_mis_psis_k, na.rm = TRUE)),
      reasons = paste(sort(unique(df$reason[bad])), collapse = ","),
      check.names = FALSE
    )
  }))
  by_local$max_impact_score[!is.finite(by_local$max_impact_score)] <- 0
  by_local$min_particle_mis_ess_frac[!is.finite(by_local$min_particle_mis_ess_frac)] <- NA_real_
  by_local$max_particle_mis_psis_k[!is.finite(by_local$max_particle_mis_psis_k)] <- NA_real_

  theta_weight <- .local_chart_normalize_weights(by_theta$theta_weight, nrow(by_theta))
  any_uncertified <- by_theta$uncertified_local_count > 0
  weighted_any_uncertified <- sum(theta_weight[any_uncertified])
  weighted_mean_uncertified <- sum(theta_weight * by_theta$uncertified_local_count)
  weighted_mean_uncertified_fraction <- weighted_mean_uncertified / max(n_locals, 1L)
  weighted_q90_uncertified <- .local_evidence_weighted_count_quantile(
    by_theta$uncertified_local_count,
    theta_weight,
    0.90
  )
  weight_cut <- as.numeric(stats::quantile(
    theta_weight,
    probs = pmin(pmax(as.numeric(thresholds$high_weight_quantile), 0), 1),
    names = FALSE,
    type = 8
  ))
  high_weight <- theta_weight >= weight_cut
  max_high_weight_uncertified <- if (any(high_weight)) {
    max(by_theta$uncertified_local_count[high_weight], na.rm = TRUE)
  } else {
    0
  }
  worst_local_uncertified_weight <- if (nrow(by_local)) {
    max(by_local$uncertified_weight, na.rm = TRUE)
  } else {
    0
  }
  failures <- character()
  if (weighted_any_uncertified > as.numeric(thresholds$max_weighted_any_uncertified)) {
    failures <- c(failures, "weighted_any_uncertified")
  }
  if (weighted_mean_uncertified > as.numeric(thresholds$max_weighted_mean_uncertified)) {
    failures <- c(failures, "weighted_mean_uncertified")
  }
  if (weighted_q90_uncertified > as.numeric(thresholds$max_weighted_q90_uncertified)) {
    failures <- c(failures, "weighted_q90_uncertified")
  }
  if (max_high_weight_uncertified > as.numeric(thresholds$max_high_weight_uncertified)) {
    failures <- c(failures, "high_weight_uncertified")
  }
  if (worst_local_uncertified_weight > as.numeric(thresholds$max_worst_local_uncertified_weight)) {
    failures <- c(failures, "worst_local_uncertified_weight")
  }

  global <- data.frame(
    n_pairs = nrow(tab),
    n_theta = n_theta,
    n_locals = n_locals,
    certified_fraction = mean(!tab$uncertified),
    uncertified_fraction = mean(tab$uncertified),
    weighted_any_uncertified = weighted_any_uncertified,
    weighted_mean_uncertified_locals = weighted_mean_uncertified,
    weighted_mean_uncertified_fraction = weighted_mean_uncertified_fraction,
    weighted_q90_uncertified_locals = weighted_q90_uncertified,
    max_high_weight_uncertified_locals = max_high_weight_uncertified,
    worst_local_uncertified_weight = worst_local_uncertified_weight,
    total_impact_score = sum(tab$impact_score, na.rm = TRUE),
    uncertified_impact_score = sum(tab$impact_score[tab$uncertified], na.rm = TRUE),
    certified = !length(failures),
    failures = paste(failures, collapse = ","),
    check.names = FALSE
  )
  structure(
    list(
      global = global,
      by_theta = by_theta[order(-by_theta$uncertified_impact_score, -by_theta$total_impact_score), , drop = FALSE],
      by_local = by_local[order(-by_local$uncertified_impact_score, -by_local$uncertified_weight), , drop = FALSE],
      thresholds = thresholds,
      failures = failures
    ),
    class = "local_evidence_certification_summary"
  )
}

.local_evidence_weighted_center_cov <- function(theta,
                                                weights,
                                                ridge = 1e-8) {
  weights <- .local_chart_normalize_weights(weights, nrow(theta))
  center <- colSums(sweep(theta, 1L, weights, "*"))
  centered <- sweep(theta, 2L, center, "-")
  cov <- crossprod(centered, centered * weights)
  cov <- (cov + t(cov)) / 2
  eig <- eigen(cov, symmetric = TRUE)
  positive <- eig$values[is.finite(eig$values) & eig$values > 0]
  floor_value <- max(as.numeric(ridge), if (length(positive)) stats::median(positive) * ridge else ridge)
  values <- pmax(eig$values, floor_value)
  whitening <- sweep(eig$vectors, 2L, sqrt(values), "/")
  unwhitening <- diag(sqrt(values), nrow = length(values)) %*% t(eig$vectors)
  list(
    center = center,
    covariance = cov,
    eigen = eig,
    regularized_values = values,
    whitening = whitening,
    unwhitening = unwhitening
  )
}

learn_certification_failure_geometry <- function(certification_table,
                                                 cloud = NULL,
                                                 max_directions = 6L,
                                                 ridge = 1e-8,
                                                 min_direction_relative_burden = 0.01) {
  if (!is.data.frame(certification_table) || !nrow(certification_table)) {
    stop("certification_table must be a non-empty data frame.")
  }
  if (is.null(cloud)) {
    cloud <- attr(certification_table, "cloud")
  }
  if (is.null(cloud)) {
    stop("cloud is required when certification_table has no cloud attribute.")
  }
  model <- normalize_population_model(attr(cloud, "population_model"))
  cloud <- validate_local_evidence_certification_cloud(cloud, model)
  theta <- cloud$theta
  theta_weight <- cloud$metadata$theta_weight
  if (!"impact_score" %in% names(certification_table)) {
    stop("certification_table must contain impact_score.")
  }
  if (!"theta_row" %in% names(certification_table)) {
    stop("certification_table must contain theta_row.")
  }
  max_directions <- as.integer(max_directions)
  if (!is.finite(max_directions) || max_directions < 1L) {
    stop("max_directions must be a positive integer.")
  }
  burden <- numeric(nrow(theta))
  theta_row <- as.integer(certification_table$theta_row)
  ok <- is.finite(theta_row) & theta_row >= 1L & theta_row <= nrow(theta)
  burden_sum <- rowsum(
    pmax(as.numeric(certification_table$impact_score[ok]), 0),
    theta_row[ok],
    reorder = FALSE
  )
  burden[as.integer(rownames(burden_sum))] <- as.numeric(burden_sum[, 1L])

  whitened <- .local_evidence_weighted_center_cov(theta, theta_weight, ridge = ridge)
  z <- sweep(theta, 2L, whitened$center, "-") %*% whitened$whitening
  colnames(z) <- model$hyper_names
  if (!any(is.finite(burden) & burden > 0)) {
    empty_directions <- data.frame(
      direction = integer(),
      eigenvalue = numeric(),
      relative_burden = numeric(),
      check.names = FALSE
    )
    return(structure(
      list(
        status = "no_failure_burden",
        center = whitened$center,
        covariance = whitened$covariance,
        whitening = whitened$whitening,
        unwhitening = whitened$unwhitening,
        whitened_theta = z,
        theta_burden = data.frame(
          theta_row = seq_len(nrow(theta)),
          theta_id = cloud$metadata$theta_id,
          theta_weight = theta_weight,
          failure_burden = burden,
          normalized_failure_burden = 0,
          check.names = FALSE
        ),
        directions = empty_directions,
        direction_vectors_whitened = matrix(numeric(), nrow = model$hyper_dim, ncol = 0L),
        direction_vectors_theta = matrix(numeric(), nrow = model$hyper_dim, ncol = 0L),
        direction_loadings = data.frame(),
        local_contributions = data.frame()
      ),
      class = "local_evidence_certification_geometry",
      population_model = model,
      cloud = cloud
    ))
  }

  burden[!is.finite(burden) | burden < 0] <- 0
  burden_weight <- burden / sum(burden)
  failure_cov <- crossprod(z, z * burden_weight)
  failure_cov <- (failure_cov + t(failure_cov)) / 2
  eig_fail <- eigen(failure_cov, symmetric = TRUE)
  eig_values <- pmax(eig_fail$values, 0)
  total_eig <- sum(eig_values)
  relative <- if (total_eig > 0) eig_values / total_eig else rep(0, length(eig_values))
  keep <- seq_len(min(max_directions, length(eig_values)))
  keep <- keep[relative[keep] >= as.numeric(min_direction_relative_burden)]
  if (!length(keep)) {
    keep <- seq_len(min(1L, length(eig_values)))
  }
  vectors_z <- eig_fail$vectors[, keep, drop = FALSE]
  vectors_theta <- sapply(seq_along(keep), function(k) {
    whitened$unwhitening %*% vectors_z[, k]
  })
  vectors_theta <- matrix(vectors_theta, nrow = model$hyper_dim)
  rownames(vectors_z) <- model$hyper_names
  rownames(vectors_theta) <- model$hyper_names
  colnames(vectors_z) <- sprintf("direction_%02d", seq_along(keep))
  colnames(vectors_theta) <- colnames(vectors_z)

  direction_rows <- lapply(seq_along(keep), function(k) {
    theta_direction <- as.numeric(vectors_theta[, k])
    theta_unit <- theta_direction / max(sqrt(sum(theta_direction^2)), .Machine$double.eps)
    top <- order(abs(theta_unit), decreasing = TRUE)
    top <- top[seq_len(min(5L, length(top)))]
    data.frame(
      direction = k,
      eigenvalue = eig_values[keep[k]],
      relative_burden = relative[keep[k]],
      top_loadings = paste(
        sprintf("%+.3f %s", theta_unit[top], model$hyper_names[top]),
        collapse = " "
      ),
      check.names = FALSE
    )
  })
  directions <- do.call(rbind, direction_rows)

  loading_rows <- lapply(seq_along(keep), function(k) {
    theta_direction <- as.numeric(vectors_theta[, k])
    theta_unit <- theta_direction / max(sqrt(sum(theta_direction^2)), .Machine$double.eps)
    ord <- order(abs(theta_unit), decreasing = TRUE)
    data.frame(
      direction = k,
      hyper_name = model$hyper_names[ord],
      loading = theta_unit[ord],
      abs_loading = abs(theta_unit[ord]),
      rank = seq_along(ord),
      check.names = FALSE
    )
  })
  direction_loadings <- do.call(rbind, loading_rows)

  local_contributions <- do.call(rbind, lapply(seq_along(keep), function(k) {
    projection2 <- as.numeric(z %*% vectors_z[, k])^2
    weighted <- certification_table
    weighted$direction_weighted_impact <- pmax(weighted$impact_score, 0) * projection2[weighted$theta_row]
    total <- sum(weighted$direction_weighted_impact, na.rm = TRUE)
    rows <- stats::aggregate(
      direction_weighted_impact ~ local + local_pos,
      data = weighted,
      FUN = sum,
      na.rm = TRUE
    )
    rows$direction <- k
    rows$direction_fraction <- if (total > 0) rows$direction_weighted_impact / total else 0
    rows[order(-rows$direction_fraction), c(
      "direction", "local", "local_pos", "direction_weighted_impact", "direction_fraction"
    ), drop = FALSE]
  }))

  theta_burden <- data.frame(
    theta_row = seq_len(nrow(theta)),
    theta_id = cloud$metadata$theta_id,
    theta_source = cloud$metadata$theta_source,
    theta_weight = theta_weight,
    failure_burden = burden,
    normalized_failure_burden = burden_weight,
    check.names = FALSE
  )
  structure(
    list(
      status = "learned",
      center = whitened$center,
      covariance = whitened$covariance,
      whitening = whitened$whitening,
      unwhitening = whitened$unwhitening,
      failure_covariance = failure_cov,
      whitened_theta = z,
      theta_burden = theta_burden[order(-theta_burden$failure_burden), , drop = FALSE],
      directions = directions,
      direction_vectors_whitened = vectors_z,
      direction_vectors_theta = vectors_theta,
      direction_loadings = direction_loadings,
      local_contributions = local_contributions
    ),
    class = "local_evidence_certification_geometry",
    population_model = model,
    cloud = cloud
  )
}

.local_evidence_failure_direction_profiles <- function(cloud,
                                                       geometry,
                                                       probs = c(0.05, 0.95),
                                                       max_directions = ncol(geometry$direction_vectors_whitened),
                                                       include_center = FALSE) {
  cloud <- validate_local_evidence_certification_cloud(
    cloud,
    population_model = attr(cloud, "population_model")
  )
  model <- normalize_population_model(attr(cloud, "population_model"))
  if (!inherits(geometry, "local_evidence_certification_geometry")) {
    stop("geometry must inherit from 'local_evidence_certification_geometry'.")
  }
  dirs <- geometry$direction_vectors_whitened
  if (!is.matrix(dirs) || !ncol(dirs)) {
    return(list(theta = matrix(numeric(), nrow = 0L, ncol = model$hyper_dim), metadata = data.frame()))
  }
  if (nrow(geometry$whitened_theta) != nrow(cloud$theta)) {
    stop("geometry must be learned from the same certification cloud being expanded.")
  }
  max_directions <- as.integer(max_directions)
  if (!is.finite(max_directions) || max_directions < 1L) {
    stop("max_directions must be a positive integer.")
  }
  n_dir <- min(max_directions, ncol(dirs))
  probs <- sort(unique(pmin(pmax(as.numeric(probs), 0), 1)))
  if (!length(probs)) {
    stop("probs must contain at least one probability.")
  }
  z <- geometry$whitened_theta
  weights <- cloud$metadata$theta_weight
  rows <- list()
  meta <- list()
  idx <- 0L
  for (k in seq_len(n_dir)) {
    projection <- as.numeric(z %*% dirs[, k])
    q <- as.numeric(.local_atlas_weighted_quantile(projection, weights, probs))
    for (p_idx in seq_along(q)) {
      value <- q[p_idx]
      if (!is.finite(value)) next
      idx <- idx + 1L
      z_new <- numeric(ncol(z))
      z_new <- value * as.numeric(dirs[, k])
      theta_new <- matrix(geometry$center, nrow = 1L) + matrix(z_new, nrow = 1L) %*% geometry$unwhitening
      colnames(theta_new) <- model$hyper_names
      rows[[idx]] <- theta_new
      meta[[idx]] <- data.frame(
        failure_direction = k,
        projection_probability = probs[p_idx],
        projection_value = value,
        profile_kind = "failure_direction",
        check.names = FALSE
      )
    }
  }
  if (isTRUE(include_center)) {
    idx <- idx + 1L
    center <- matrix(geometry$center, nrow = 1L, dimnames = list(NULL, model$hyper_names))
    rows[[idx]] <- center
    meta[[idx]] <- data.frame(
      failure_direction = NA_integer_,
      projection_probability = NA_real_,
      projection_value = 0,
      profile_kind = "failure_center",
      check.names = FALSE
    )
  }
  if (!length(rows)) {
    return(list(theta = matrix(numeric(), nrow = 0L, ncol = model$hyper_dim), metadata = data.frame()))
  }
  theta <- do.call(rbind, rows)
  colnames(theta) <- model$hyper_names
  profile_meta <- do.call(rbind, meta)
  list(theta = theta, metadata = profile_meta)
}

add_failure_direction_profiles <- function(cloud,
                                           geometry,
                                           probs = c(0.05, 0.95),
                                           max_directions = 3L,
                                           profile_weight_mass = 0.05,
                                           theta_source = "failure_direction",
                                           theta_round = NULL,
                                           include_center = FALSE,
                                           deduplicate = TRUE) {
  cloud <- validate_local_evidence_certification_cloud(
    cloud,
    population_model = attr(cloud, "population_model")
  )
  model <- normalize_population_model(attr(cloud, "population_model"))
  profiles <- .local_evidence_failure_direction_profiles(
    cloud = cloud,
    geometry = geometry,
    probs = probs,
    max_directions = max_directions,
    include_center = include_center
  )
  if (!nrow(profiles$theta)) {
    attr(cloud, "added_profiles") <- profiles
    return(cloud)
  }
  profile_weight_mass <- min(max(as.numeric(profile_weight_mass), 0), 1)
  original_weight <- cloud$metadata$theta_weight * (1 - profile_weight_mass)
  profile_weight <- rep(profile_weight_mass / nrow(profiles$theta), nrow(profiles$theta))
  theta <- rbind(cloud$theta, profiles$theta)
  profile_round <- if (is.null(theta_round)) {
    max(cloud$metadata$theta_round, na.rm = TRUE) + 1L
  } else {
    as.integer(theta_round)
  }
  profile_round_id <- as.integer(profile_round[1L])
  if (!is.finite(profile_round_id)) {
    profile_round_id <- 0L
  }
  profile_direction <- as.integer(profiles$metadata$failure_direction)
  profile_direction[!is.finite(profile_direction)] <- 0L
  metadata <- data.frame(
    failure_direction = c(rep(NA_integer_, nrow(cloud$theta)), profiles$metadata$failure_direction),
    projection_probability = c(rep(NA_real_, nrow(cloud$theta)), profiles$metadata$projection_probability),
    projection_value = c(rep(NA_real_, nrow(cloud$theta)), profiles$metadata$projection_value),
    profile_kind = c(rep("original", nrow(cloud$theta)), profiles$metadata$profile_kind),
    check.names = FALSE
  )
  out <- build_local_evidence_certification_cloud(
    theta = theta,
    population_model = model,
    theta_weights = c(original_weight, profile_weight),
    theta_source = c(cloud$metadata$theta_source, rep(theta_source, nrow(profiles$theta))),
    theta_round = c(cloud$metadata$theta_round, rep(profile_round, nrow(profiles$theta))),
    theta_id = c(
      cloud$metadata$theta_id,
      sprintf("failure_dir_r%02d_d%02d_%03d", profile_round_id, profile_direction, seq_len(nrow(profiles$theta)))
    ),
    metadata = metadata,
    deduplicate = isTRUE(deduplicate)
  )
  attr(out, "added_profiles") <- profiles
  out
}

.local_evidence_repair_selection_limits <- function(max_repairs,
                                                    max_repairs_per_local,
                                                    max_repairs_per_direction,
                                                    max_repairs_per_theta) {
  max_repairs <- as.integer(max_repairs)
  if (!is.finite(max_repairs) || max_repairs < 1L) {
    stop("max_repairs must be a positive integer.")
  }
  clean_limit <- function(x, default) {
    if (is.null(x)) {
      return(as.integer(default))
    }
    if (is.infinite(x)) {
      return(Inf)
    }
    x <- as.integer(x)
    if (!is.finite(x) || x < 1L) {
      stop("repair selection caps must be positive integers or Inf.")
    }
    x
  }
  list(
    max_repairs = max_repairs,
    max_repairs_per_local = clean_limit(max_repairs_per_local, max(1L, ceiling(max_repairs / 4))),
    max_repairs_per_direction = clean_limit(max_repairs_per_direction, max(1L, ceiling(max_repairs / 3))),
    max_repairs_per_theta = clean_limit(max_repairs_per_theta, max(1L, ceiling(max_repairs / 4)))
  )
}

.local_evidence_annotate_repair_candidates <- function(certification_table,
                                                       failure_geometry = NULL,
                                                       only_uncertified = TRUE,
                                                       min_impact_score = 0,
                                                       exclude_selected = TRUE) {
  if (!is.data.frame(certification_table) || !nrow(certification_table)) {
    stop("certification_table must be a non-empty data frame.")
  }
  required <- c(
    "local", "local_pos", "theta_row", "theta_id", "theta_weight",
    "status", "impact_score"
  )
  missing <- setdiff(required, names(certification_table))
  if (length(missing)) {
    stop("certification_table is missing: ", paste(missing, collapse = ", "))
  }
  table <- certification_table
  table$certification_row <- seq_len(nrow(table))
  table$uncertified <- is.na(table$status) | table$status != "certified"
  table$impact_score <- pmax(as.numeric(table$impact_score), 0)
  table$impact_score[!is.finite(table$impact_score)] <- 0

  local_burden <- stats::aggregate(
    impact_score ~ local_pos,
    data = table,
    FUN = sum,
    na.rm = TRUE
  )
  names(local_burden)[names(local_burden) == "impact_score"] <- "local_burden"
  theta_burden <- stats::aggregate(
    impact_score ~ theta_row,
    data = table,
    FUN = sum,
    na.rm = TRUE
  )
  names(theta_burden)[names(theta_burden) == "impact_score"] <- "theta_burden"
  table <- merge(table, local_burden, by = "local_pos", all.x = TRUE, sort = FALSE)
  table <- merge(table, theta_burden, by = "theta_row", all.x = TRUE, sort = FALSE)
  table$local_burden[!is.finite(table$local_burden)] <- 0
  table$theta_burden[!is.finite(table$theta_burden)] <- 0

  table$failure_direction <- NA_integer_
  table$failure_direction_projection <- NA_real_
  table$failure_direction_score <- 0
  table$failure_direction_relative_burden <- NA_real_
  if (!is.null(failure_geometry)) {
    if (!inherits(failure_geometry, "local_evidence_certification_geometry")) {
      stop("failure_geometry must inherit from 'local_evidence_certification_geometry'.")
    }
    dirs <- failure_geometry$direction_vectors_whitened
    z <- failure_geometry$whitened_theta
    if (is.matrix(dirs) && ncol(dirs) &&
        is.matrix(z) && max(table$theta_row, na.rm = TRUE) <= nrow(z)) {
      rel <- rep(1, ncol(dirs))
      if (is.data.frame(failure_geometry$directions) &&
          all(c("direction", "relative_burden") %in% names(failure_geometry$directions))) {
        rel_match <- match(seq_len(ncol(dirs)), failure_geometry$directions$direction)
        rel_from_geometry <- as.numeric(failure_geometry$directions$relative_burden[rel_match])
        rel[is.finite(rel_from_geometry)] <- rel_from_geometry[is.finite(rel_from_geometry)]
      }
      projections <- z[table$theta_row, , drop = FALSE] %*% dirs
      direction_scores <- sweep(projections^2, 2L, pmax(rel, .Machine$double.eps), "*")
      best <- max.col(direction_scores, ties.method = "first")
      best_score <- direction_scores[cbind(seq_len(nrow(direction_scores)), best)]
      has_direction <- is.finite(best_score) & best_score > 0
      table$failure_direction[has_direction] <- best[has_direction]
      table$failure_direction_projection[has_direction] <-
        projections[cbind(which(has_direction), best[has_direction])]
      table$failure_direction_score[has_direction] <-
        table$impact_score[has_direction] * best_score[has_direction]
      table$failure_direction_relative_burden[has_direction] <- rel[best[has_direction]]
    }
  }

  selected_before <- if ("selected_for_repair" %in% names(table)) {
    table$selected_for_repair %in% TRUE
  } else {
    rep(FALSE, nrow(table))
  }
  eligible <- table$impact_score > as.numeric(min_impact_score)
  if (isTRUE(only_uncertified)) {
    eligible <- eligible & table$uncertified
  }
  if (isTRUE(exclude_selected)) {
    eligible <- eligible & !selected_before
  }
  table$repair_eligible <- eligible
  table$repair_rank_score <- table$impact_score *
    (1 + log1p(pmax(table$local_burden, 0)) + log1p(pmax(table$theta_burden, 0))) *
    (1 + log1p(pmax(table$failure_direction_score, 0)))
  table$repair_rank_score[!is.finite(table$repair_rank_score)] <- 0
  table[order(
    -table$repair_rank_score,
    -table$impact_score,
    -table$theta_weight,
    table$local_pos,
    table$theta_row
  ), , drop = FALSE]
}

select_certification_repairs <- function(certification_table,
                                         failure_geometry = NULL,
                                         max_repairs = 16L,
                                         max_repairs_per_local = NULL,
                                         max_repairs_per_direction = NULL,
                                         max_repairs_per_theta = NULL,
                                         high_weight_quantile = 0.90,
                                         min_impact_score = 0,
                                         only_uncertified = TRUE,
                                         exclude_selected = TRUE) {
  limits <- .local_evidence_repair_selection_limits(
    max_repairs = max_repairs,
    max_repairs_per_local = max_repairs_per_local,
    max_repairs_per_direction = max_repairs_per_direction,
    max_repairs_per_theta = max_repairs_per_theta
  )
  table <- .local_evidence_annotate_repair_candidates(
    certification_table = certification_table,
    failure_geometry = failure_geometry,
    only_uncertified = only_uncertified,
    min_impact_score = min_impact_score,
    exclude_selected = exclude_selected
  )
  table$selected_for_repair <- FALSE
  table$repair_selection_order <- NA_integer_
  table$repair_selection_reason <- NA_character_
  table$repair_selection_tags <- NA_character_

  candidates <- table[table$repair_eligible, , drop = FALSE]
  if (!nrow(candidates)) {
    empty <- table[FALSE, , drop = FALSE]
    return(structure(
      list(
        repairs = empty,
        candidate_pairs = data.frame(local_pos = integer(), theta_row = integer()),
        annotated_table = table[order(table$certification_row), , drop = FALSE],
        summary = data.frame(
          max_repairs = limits$max_repairs,
          n_candidates = 0L,
          n_selected = 0L,
          n_selected_locals = 0L,
          n_selected_theta = 0L,
          n_selected_directions = 0L,
          selected_impact_score = 0,
          candidate_impact_score = 0,
          check.names = FALSE
        )
      ),
      class = "local_evidence_certification_repair_selection"
    ))
  }

  selected_rows <- integer()
  selected_keys <- character()
  selection_tags <- rep("", nrow(table))
  local_counts <- integer()
  theta_counts <- integer()
  direction_counts <- integer()
  reasons <- character()
  names(local_counts) <- character()
  names(theta_counts) <- character()
  names(direction_counts) <- character()

  current_count <- function(counts, key) {
    key <- as.character(key)
    if (key %in% names(counts)) counts[[key]] else 0L
  }
  set_count <- function(counts, key, value) {
    key <- as.character(key)
    counts[[key]] <- as.integer(value)
    counts
  }
  within_limit <- function(value, limit) {
    is.infinite(limit) || value < as.integer(limit)
  }
  record_tag <- function(row_index, reason) {
    old <- selection_tags[[row_index]]
    tags <- if (nzchar(old)) strsplit(old, ",", fixed = TRUE)[[1L]] else character()
    tags <- unique(c(tags, reason))
    selection_tags[[row_index]] <<- paste(tags, collapse = ",")
  }
  add_one <- function(row_index, reason) {
    if (length(selected_rows) >= limits$max_repairs) {
      return(FALSE)
    }
    row <- table[row_index, , drop = FALSE]
    key <- paste(row$local_pos, row$theta_row, sep = "\r")
    if (key %in% selected_keys) {
      record_tag(row_index, reason)
      return(FALSE)
    }
    local_key <- as.character(row$local_pos)
    theta_key <- as.character(row$theta_row)
    direction_key <- if (is.finite(row$failure_direction)) {
      as.character(as.integer(row$failure_direction))
    } else {
      NA_character_
    }
    if (!within_limit(current_count(local_counts, local_key), limits$max_repairs_per_local)) {
      return(FALSE)
    }
    if (!within_limit(current_count(theta_counts, theta_key), limits$max_repairs_per_theta)) {
      return(FALSE)
    }
    if (!is.na(direction_key) &&
        !within_limit(current_count(direction_counts, direction_key), limits$max_repairs_per_direction)) {
      return(FALSE)
    }
    selected_rows <<- c(selected_rows, row_index)
    selected_keys <<- c(selected_keys, key)
    record_tag(row_index, reason)
    local_counts <<- set_count(local_counts, local_key, current_count(local_counts, local_key) + 1L)
    theta_counts <<- set_count(theta_counts, theta_key, current_count(theta_counts, theta_key) + 1L)
    if (!is.na(direction_key)) {
      direction_counts <<- set_count(
        direction_counts,
        direction_key,
        current_count(direction_counts, direction_key) + 1L
      )
    }
    reasons <<- c(reasons, reason)
    TRUE
  }
  add_best_by_group <- function(group_col,
                                reason,
                                group_order = NULL,
                                score_col = "repair_rank_score",
                                max_add = Inf) {
    groups <- unique(candidates[[group_col]])
    if (!is.null(group_order)) {
      groups <- group_order[group_order %in% groups]
    }
    added <- 0L
    for (group in groups) {
      if (length(selected_rows) >= limits$max_repairs) break
      if (!is.infinite(max_add) && added >= as.integer(max_add)) break
      idx <- which(table$repair_eligible & table[[group_col]] == group)
      if (!length(idx)) next
      idx <- idx[order(
        -table[[score_col]][idx],
        -table$impact_score[idx],
        -table$theta_weight[idx],
        table$local_pos[idx],
        table$theta_row[idx]
      )]
      for (row_index in idx) {
        if (add_one(row_index, reason)) {
          added <- added + 1L
          break
        }
      }
    }
    invisible(added)
  }
  add_ranked <- function(idx, reason, score_col = "repair_rank_score") {
    idx <- idx[table$repair_eligible[idx]]
    idx <- idx[order(
      -table[[score_col]][idx],
      -table$impact_score[idx],
      -table$theta_weight[idx],
      table$local_pos[idx],
      table$theta_row[idx]
    )]
    for (row_index in idx) {
      if (length(selected_rows) >= limits$max_repairs) break
      add_one(row_index, reason)
    }
  }

  category_quota <- max(1L, floor(limits$max_repairs / 5L))
  local_order <- unique(candidates$local_pos[order(-candidates$local_burden, -candidates$impact_score)])
  add_best_by_group(
    "local_pos",
    "highest_local_burden",
    local_order,
    max_add = min(length(local_order), category_quota)
  )

  theta_order <- unique(candidates$theta_row[order(-candidates$theta_burden, -candidates$impact_score)])
  add_best_by_group(
    "theta_row",
    "highest_theta_burden",
    theta_order,
    max_add = min(length(theta_order), category_quota)
  )

  direction_candidates <- candidates[is.finite(candidates$failure_direction), , drop = FALSE]
  if (nrow(direction_candidates)) {
    direction_order <- unique(direction_candidates$failure_direction[
      order(
        -direction_candidates$failure_direction_relative_burden,
        -direction_candidates$failure_direction_score,
        -direction_candidates$impact_score
      )
    ])
    direction_quota <- max(1L, min(length(direction_order), ceiling(limits$max_repairs / 4L)))
    add_best_by_group("failure_direction", "top_failure_direction", direction_order,
                      score_col = "failure_direction_score",
                      max_add = direction_quota)
  }

  high_weight_quantile <- pmin(pmax(as.numeric(high_weight_quantile), 0), 1)
  weight_cut <- as.numeric(stats::quantile(
    candidates$theta_weight,
    probs = high_weight_quantile,
    names = FALSE,
    type = 8
  ))
  high_weight_theta <- unique(candidates$theta_row[candidates$theta_weight >= weight_cut])
  high_weight_theta <- high_weight_theta[order(
    -vapply(high_weight_theta, function(theta_row) {
      max(candidates$theta_weight[candidates$theta_row == theta_row], na.rm = TRUE)
    }, numeric(1)),
    -vapply(high_weight_theta, function(theta_row) {
      max(candidates$theta_burden[candidates$theta_row == theta_row], na.rm = TRUE)
    }, numeric(1))
  )]
  add_best_by_group(
    "theta_row",
    "high_weight_theta",
    high_weight_theta,
    max_add = min(length(high_weight_theta), category_quota)
  )

  add_ranked(which(table$repair_eligible), "global_impact_fill")

  if (length(selected_rows)) {
    table$selected_for_repair[selected_rows] <- TRUE
    table$repair_selection_order[selected_rows] <- seq_along(selected_rows)
    table$repair_selection_reason[selected_rows] <- reasons
    table$repair_selection_tags[selected_rows] <- selection_tags[selected_rows]
  }
  repairs <- table[selected_rows, , drop = FALSE]
  repairs <- repairs[order(repairs$repair_selection_order), , drop = FALSE]
  candidate_pairs <- unique(repairs[, c("local_pos", "theta_row"), drop = FALSE])
  rownames(candidate_pairs) <- NULL

  summary <- data.frame(
    max_repairs = limits$max_repairs,
    max_repairs_per_local = limits$max_repairs_per_local,
    max_repairs_per_direction = limits$max_repairs_per_direction,
    max_repairs_per_theta = limits$max_repairs_per_theta,
    n_candidates = nrow(candidates),
    n_selected = nrow(repairs),
    n_selected_locals = length(unique(repairs$local_pos)),
    n_selected_theta = length(unique(repairs$theta_row)),
    n_selected_directions = length(unique(repairs$failure_direction[is.finite(repairs$failure_direction)])),
    selected_impact_score = sum(repairs$impact_score, na.rm = TRUE),
    candidate_impact_score = sum(candidates$impact_score, na.rm = TRUE),
    selected_impact_fraction = if (sum(candidates$impact_score, na.rm = TRUE) > 0) {
      sum(repairs$impact_score, na.rm = TRUE) / sum(candidates$impact_score, na.rm = TRUE)
    } else {
      NA_real_
    },
    check.names = FALSE
  )
  structure(
    list(
      repairs = repairs,
      candidate_pairs = candidate_pairs,
      annotated_table = table[order(table$certification_row), , drop = FALSE],
      summary = summary,
      limits = limits
    ),
    class = "local_evidence_certification_repair_selection"
  )
}

validate_raw_atlas_certification <- function(x) {
  if (!inherits(x, "raw_local_evidence_certification")) {
    stop("x must inherit from 'raw_local_evidence_certification'.")
  }
  required <- c(
    "factor_set", "cloud", "certification_table", "summary",
    "geometry", "repair_selection", "rounds", "certified", "status"
  )
  missing <- setdiff(required, names(x))
  if (length(missing)) {
    stop("raw certification object is missing: ", paste(missing, collapse = ", "))
  }
  x$factor_set <- validate_local_atlas_factor_set(x$factor_set)
  x$cloud <- validate_local_evidence_certification_cloud(x$cloud, x$factor_set$population_model)
  if (!is.data.frame(x$certification_table) || !nrow(x$certification_table)) {
    stop("raw certification object must contain a non-empty certification_table.")
  }
  if (!inherits(x$summary, "local_evidence_certification_summary")) {
    stop("raw certification summary has the wrong class.")
  }
  x$certified <- isTRUE(x$certified)
  x$status <- as.character(x$status)
  x
}

.local_evidence_raw_certification_calibration_control <- function(local_control,
                                                                  calibration_control) {
  local_control <- .local_atlas_merge_control(local_control, .local_atlas_default_local_control())
  .local_atlas_merge_control(calibration_control, list(
    M = as.integer(local_control$candidate_M %||% 500L),
    target_cess = local_control$target_cess %||% 0.9,
    n_mcmc_moves = as.integer(local_control$n_mcmc_moves %||% 2L),
    max_steps = as.integer(local_control$max_steps %||% 128L),
    confirmation_reps = 0L,
    confirmation_M = as.integer(local_control$candidate_M %||% 500L),
    confirmation_max_sd = 1.5,
    replicate_bootstrap_B = 200L,
    max_direct_graph_z = 3,
    max_direct_graph_chart_shift = 0.35,
    max_direct_graph_existing_shift = 0.15
  ))
}

.local_evidence_raw_certification_history_row <- function(round_id,
                                                          cloud,
                                                          summary,
                                                          geometry,
                                                          repair_selection,
                                                          calibration = NULL) {
  global <- summary$global
  selection_summary <- if (!is.null(repair_selection)) {
    repair_selection$summary
  } else {
    data.frame(
      n_candidates = NA_integer_,
      n_selected = NA_integer_,
      n_selected_locals = NA_integer_,
      n_selected_theta = NA_integer_,
      n_selected_directions = NA_integer_,
      selected_impact_score = NA_real_,
      candidate_impact_score = NA_real_,
      selected_impact_fraction = NA_real_,
      check.names = FALSE
    )
  }
  top_direction_burden <- if (!is.null(geometry) &&
                              is.data.frame(geometry$directions) &&
                              nrow(geometry$directions)) {
    max(geometry$directions$relative_burden, na.rm = TRUE)
  } else {
    NA_real_
  }
  data.frame(
    round = as.integer(round_id),
    n_theta = nrow(cloud$theta),
    certified = isTRUE(global$certified[1L]),
    failures = as.character(global$failures[1L]),
    uncertified_fraction = as.numeric(global$uncertified_fraction[1L]),
    weighted_any_uncertified = as.numeric(global$weighted_any_uncertified[1L]),
    weighted_mean_uncertified_locals = as.numeric(global$weighted_mean_uncertified_locals[1L]),
    weighted_q90_uncertified_locals = as.numeric(global$weighted_q90_uncertified_locals[1L]),
    total_impact_score = as.numeric(global$total_impact_score[1L]),
    uncertified_impact_score = as.numeric(global$uncertified_impact_score[1L]),
    top_direction_relative_burden = top_direction_burden,
    n_repair_candidates = as.integer(selection_summary$n_candidates[1L]),
    n_selected_repairs = as.integer(selection_summary$n_selected[1L]),
    n_selected_repair_locals = as.integer(selection_summary$n_selected_locals[1L]),
    n_selected_repair_theta = as.integer(selection_summary$n_selected_theta[1L]),
    selected_impact_score = as.numeric(selection_summary$selected_impact_score[1L]),
    selected_impact_fraction = as.numeric(selection_summary$selected_impact_fraction[1L]),
    calibration_probes = if (!is.null(calibration) && is.data.frame(calibration$probes)) nrow(calibration$probes) else NA_integer_,
    calibration_activated = if (!is.null(calibration)) as.integer(calibration$n_activated %||% NA_integer_) else NA_integer_,
    check.names = FALSE
  )
}

certify_and_repair_raw_atlas <- function(factor_set,
                                         cloud = NULL,
                                         theta = NULL,
                                         theta_weights = NULL,
                                         theta_source = "theta_cloud",
                                         theta_round = 0L,
                                         theta_id = NULL,
                                         metadata = NULL,
                                         data_list = NULL,
                                         loglik_fn = NULL,
                                         local_ids = NULL,
                                         local_control = list(),
                                         edge_control = list(),
                                         calibration_control = list(),
                                         certification_thresholds = list(),
                                         scoring_control = list(),
                                         geometry_control = list(),
                                         repair_control = list(),
                                         profile_control = list(),
                                         max_rounds = 2L,
                                         stop_on_uncertified = FALSE,
                                         keep_round_tables = FALSE,
                                         n_cores = 1L,
                                         seed = 123L,
                                         verbose = TRUE,
                                         trace_verbose = FALSE) {
  factor_set <- validate_local_atlas_factor_set(factor_set)
  model <- factor_set$population_model
  cloud <- .local_evidence_certification_as_cloud(
    cloud = cloud,
    theta = theta,
    population_model = model,
    theta_weights = theta_weights,
    theta_source = theta_source,
    theta_round = theta_round,
    theta_id = theta_id,
    metadata = metadata
  )
  max_rounds <- as.integer(max_rounds)
  if (!is.finite(max_rounds) || max_rounds < 0L) {
    stop("max_rounds must be a non-negative integer.")
  }
  if (!is.null(local_ids) && is.character(local_ids)) {
    local_ids <- match(local_ids, names(factor_set$atlases))
  }
  local_ids <- local_ids %||% seq_along(factor_set$atlases)
  local_ids <- as.integer(local_ids)
  if (!length(local_ids) ||
      any(!is.finite(local_ids)) ||
      any(local_ids < 1L | local_ids > length(factor_set$atlases))) {
    stop("local_ids must identify atlases in factor_set.")
  }
  local_ids <- unique(local_ids)
  local_control <- .local_atlas_merge_control(local_control, .local_atlas_default_local_control())
  edge_control <- .local_atlas_merge_control(edge_control, .local_atlas_default_edge_control())
  calibration_control <- .local_evidence_raw_certification_calibration_control(
    local_control = local_control,
    calibration_control = calibration_control
  )
  geometry_control <- modifyList(
    list(
      max_directions = 6L,
      ridge = 1e-8,
      min_direction_relative_burden = 0.01
    ),
    geometry_control
  )
  repair_control <- modifyList(
    list(
      max_repairs = 16L,
      max_repairs_per_local = NULL,
      max_repairs_per_direction = NULL,
      max_repairs_per_theta = NULL,
      high_weight_quantile = 0.90,
      min_impact_score = 0,
      only_uncertified = TRUE,
      exclude_selected = TRUE
    ),
    repair_control
  )
  profile_control <- modifyList(
    list(
      enabled = TRUE,
      probs = c(0.05, 0.95),
      max_directions = 3L,
      profile_weight_mass = 0.05,
      include_center = FALSE
    ),
    profile_control
  )

  rounds <- list()
  history <- list()
  final_table <- NULL
  final_summary <- NULL
  final_geometry <- NULL
  final_selection <- NULL
  status <- "uncertified_budget_exhausted"
  current_cloud <- cloud

  for (round_id in 0:max_rounds) {
    table <- evaluate_raw_local_evidence_certification(
      factor_set = factor_set,
      cloud = current_cloud,
      local_ids = local_ids,
      scoring_control = scoring_control,
      include_theta = FALSE,
      n_cores = as.integer(n_cores)
    )
    summary <- summarize_raw_local_evidence_certification(
      table,
      thresholds = certification_thresholds
    )
    geometry <- learn_certification_failure_geometry(
      table,
      cloud = current_cloud,
      max_directions = as.integer(geometry_control$max_directions),
      ridge = as.numeric(geometry_control$ridge),
      min_direction_relative_burden = as.numeric(geometry_control$min_direction_relative_burden)
    )
    selection <- NULL
    calibration <- NULL
    certified <- isTRUE(summary$global$certified[1L])
    if (isTRUE(certified)) {
      status <- "certified"
      history[[length(history) + 1L]] <- .local_evidence_raw_certification_history_row(
        round_id = round_id,
        cloud = current_cloud,
        summary = summary,
        geometry = geometry,
        repair_selection = NULL,
        calibration = NULL
      )
      rounds[[length(rounds) + 1L]] <- list(
        round = round_id,
        cloud = current_cloud,
        summary = summary,
        geometry = geometry,
        repair_selection = NULL,
        calibration = NULL,
        certification_table = if (isTRUE(keep_round_tables)) table else NULL
      )
      final_table <- table
      final_summary <- summary
      final_geometry <- geometry
      final_selection <- NULL
      break
    }

    if (round_id >= max_rounds) {
      status <- "uncertified_budget_exhausted"
      history[[length(history) + 1L]] <- .local_evidence_raw_certification_history_row(
        round_id = round_id,
        cloud = current_cloud,
        summary = summary,
        geometry = geometry,
        repair_selection = NULL,
        calibration = NULL
      )
      rounds[[length(rounds) + 1L]] <- list(
        round = round_id,
        cloud = current_cloud,
        summary = summary,
        geometry = geometry,
        repair_selection = NULL,
        calibration = NULL,
        certification_table = if (isTRUE(keep_round_tables)) table else NULL
      )
      final_table <- table
      final_summary <- summary
      final_geometry <- geometry
      final_selection <- NULL
      break
    }

    selection <- select_certification_repairs(
      certification_table = table,
      failure_geometry = geometry,
      max_repairs = as.integer(repair_control$max_repairs),
      max_repairs_per_local = repair_control$max_repairs_per_local,
      max_repairs_per_direction = repair_control$max_repairs_per_direction,
      max_repairs_per_theta = repair_control$max_repairs_per_theta,
      high_weight_quantile = as.numeric(repair_control$high_weight_quantile),
      min_impact_score = as.numeric(repair_control$min_impact_score),
      only_uncertified = isTRUE(repair_control$only_uncertified),
      exclude_selected = isTRUE(repair_control$exclude_selected)
    )
    table <- selection$annotated_table
    if (!nrow(selection$repairs)) {
      status <- "uncertified_no_selected_repairs"
      history[[length(history) + 1L]] <- .local_evidence_raw_certification_history_row(
        round_id = round_id,
        cloud = current_cloud,
        summary = summary,
        geometry = geometry,
        repair_selection = selection,
        calibration = NULL
      )
      rounds[[length(rounds) + 1L]] <- list(
        round = round_id,
        cloud = current_cloud,
        summary = summary,
        geometry = geometry,
        repair_selection = selection,
        calibration = NULL,
        certification_table = if (isTRUE(keep_round_tables)) table else NULL
      )
      final_table <- table
      final_summary <- summary
      final_geometry <- geometry
      final_selection <- selection
      break
    }
    if (is.null(data_list) || !is.function(loglik_fn)) {
      stop("Raw atlas repair requires data_list and loglik_fn once certification fails.")
    }
    if (!is.list(data_list) || length(data_list) < max(selection$candidate_pairs$local_pos, na.rm = TRUE)) {
      stop("data_list must contain every local selected for raw atlas repair.")
    }

    .local_atlas_log(
      "raw certification round ", round_id,
      ": uncertified_weight=", sprintf("%.4f", summary$global$weighted_any_uncertified[1L]),
      " selected_repairs=", nrow(selection$repairs),
      " selected_locals=", length(unique(selection$repairs$local_pos)),
      " n_theta=", nrow(current_cloud$theta), "\n",
      verbose = verbose
    )
    calibration <- local_atlas_repair_certification_pairs(
      factor_set = factor_set,
      theta = current_cloud$theta,
      data_list = data_list,
      loglik_fn = loglik_fn,
      theta_weights = current_cloud$metadata$theta_weight,
      candidate_pairs = selection$candidate_pairs,
      local_ids = unique(selection$candidate_pairs$local_pos),
      M = as.integer(calibration_control$M),
      target_cess = calibration_control$target_cess,
      n_mcmc_moves = as.integer(calibration_control$n_mcmc_moves),
      max_steps = as.integer(calibration_control$max_steps),
      max_updates = min(as.integer(calibration_control$max_updates %||% repair_control$max_repairs), nrow(selection$candidate_pairs)),
      direct_confirmation_reps = as.integer(calibration_control$confirmation_reps),
      direct_confirmation_M = as.integer(calibration_control$confirmation_M),
      direct_confirmation_max_sd = calibration_control$confirmation_max_sd,
      replicate_bootstrap_B = as.integer(calibration_control$replicate_bootstrap_B),
      max_direct_graph_z = calibration_control$max_direct_graph_z,
      max_direct_graph_chart_shift = calibration_control$max_direct_graph_chart_shift,
      max_direct_graph_existing_shift = calibration_control$max_direct_graph_existing_shift,
      local_control = local_control,
      edge_control = edge_control,
      seed = as.integer(seed) + 7100003L * (round_id + 1L),
      verbose = isTRUE(trace_verbose)
    )
    factor_set <- calibration$factor_set
    history[[length(history) + 1L]] <- .local_evidence_raw_certification_history_row(
      round_id = round_id,
      cloud = current_cloud,
      summary = summary,
      geometry = geometry,
      repair_selection = selection,
      calibration = calibration
    )
    rounds[[length(rounds) + 1L]] <- list(
      round = round_id,
      cloud = current_cloud,
      summary = summary,
      geometry = geometry,
      repair_selection = selection,
      calibration = calibration,
      certification_table = if (isTRUE(keep_round_tables)) table else NULL
    )

    if (isTRUE(profile_control$enabled)) {
      current_cloud <- add_failure_direction_profiles(
        cloud = current_cloud,
        geometry = geometry,
        probs = profile_control$probs,
        max_directions = as.integer(profile_control$max_directions),
        profile_weight_mass = as.numeric(profile_control$profile_weight_mass),
        theta_round = round_id + 1L,
        include_center = isTRUE(profile_control$include_center),
        deduplicate = TRUE
      )
    }
    final_table <- table
    final_summary <- summary
    final_geometry <- geometry
    final_selection <- selection
    if (!isTRUE(calibration$n_activated > 0L)) {
      status <- "uncertified_no_repairs_activated"
      break
    }
  }

  history <- if (length(history)) do.call(rbind, history) else data.frame()
  out <- structure(
    list(
      factor_set = factor_set,
      atlases = factor_set$atlases,
      cloud = current_cloud,
      certification_table = final_table,
      summary = final_summary,
      geometry = final_geometry,
      repair_selection = final_selection,
      rounds = rounds,
      history = history,
      certified = identical(status, "certified"),
      status = status,
      settings = list(
        max_rounds = max_rounds,
        stop_on_uncertified = isTRUE(stop_on_uncertified),
        keep_round_tables = isTRUE(keep_round_tables),
        local_ids = local_ids,
        certification_thresholds = certification_thresholds,
        scoring_control = scoring_control,
        geometry_control = geometry_control,
        repair_control = repair_control,
        profile_control = profile_control,
        calibration_control = calibration_control,
        seed = as.integer(seed)
      )
    ),
    class = "raw_local_evidence_certification"
  ) |>
    validate_raw_atlas_certification()
  if (!isTRUE(out$certified) && isTRUE(stop_on_uncertified)) {
    failures <- as.character(out$summary$global$failures[1L])
    stop(
      "Raw local evidence certification failed: status=", out$status,
      if (nzchar(failures)) paste0(" failures=", failures) else ""
    )
  }
  out
}

validate_certified_outer_workflow <- function(x) {
  if (!inherits(x, "certified_outer_workflow")) {
    stop("x must inherit from 'certified_outer_workflow'.")
  }
  required <- c(
    "factor_set", "pre_outer_certification", "pilot_fit",
    "post_pilot_certification", "final_fit", "status", "history", "settings"
  )
  missing <- setdiff(required, names(x))
  if (length(missing)) {
    stop("certified outer workflow is missing: ", paste(missing, collapse = ", "))
  }
  x$factor_set <- validate_local_atlas_factor_set(x$factor_set)
  if (!is.null(x$pre_outer_certification)) {
    x$pre_outer_certification <- validate_raw_atlas_certification(x$pre_outer_certification)
  }
  if (!is.null(x$post_pilot_certification)) {
    x$post_pilot_certification <- validate_raw_atlas_certification(x$post_pilot_certification)
  }
  if (!is.data.frame(x$history)) {
    stop("certified outer workflow history must be a data frame.")
  }
  x$status <- as.character(x$status)
  x
}

.local_evidence_outer_workflow_checkpoint <- function(checkpoint_file,
                                                      stage,
                                                      state,
                                                      verbose = TRUE) {
  if (is.null(checkpoint_file) || !nzchar(checkpoint_file)) {
    return(invisible(FALSE))
  }
  .local_atlas_write_checkpoint(
    path = checkpoint_file,
    stage = stage,
    state = state,
    verbose = verbose
  )
  invisible(TRUE)
}

.local_evidence_run_outer <- function(factor_set,
                                      initial_proposal,
                                      outer_control,
                                      n_cores,
                                      seed,
                                      verbose) {
  defaults <- list(
    N = 1000L,
    resample_threshold = 0.5,
    n_mcmc_moves = 3L,
    min_mcmc_moves = 1L,
    max_rounds = 80L,
    rw_scale_init = 0.8,
    verbose = isTRUE(verbose)
  )
  args <- modifyList(defaults, outer_control %||% list())
  args$factor_set <- factor_set
  args$initial_proposal <- initial_proposal
  args$n_cores <- as.integer(n_cores)
  args$seed <- as.integer(seed)
  do.call(outer_population_smc, args)
}

.local_evidence_outer_fit_tail_profiles <- function(theta,
                                                    weights,
                                                    population_model,
                                                    tail_probs = c(0.025, 0.975),
                                                    tail_directions = 4L,
                                                    tail_inflation = 1.25) {
  model <- normalize_population_model(population_model)
  theta <- .as_hyper_matrix(theta, model$hyper_names, model$hyper_dim)
  weights <- .local_chart_normalize_weights(weights, nrow(theta))
  tail_directions <- min(as.integer(tail_directions), model$hyper_dim)
  if (!is.finite(tail_directions) || tail_directions < 1L) {
    return(list(theta = matrix(numeric(), nrow = 0L, ncol = model$hyper_dim), metadata = data.frame()))
  }
  tail_probs <- sort(unique(pmin(pmax(as.numeric(tail_probs), 0), 1)))
  tail_probs <- tail_probs[is.finite(tail_probs)]
  if (!length(tail_probs)) {
    return(list(theta = matrix(numeric(), nrow = 0L, ncol = model$hyper_dim), metadata = data.frame()))
  }
  wcov <- .local_evidence_weighted_center_cov(theta, weights)
  eig <- eigen(wcov$covariance, symmetric = TRUE)
  positive <- which(is.finite(eig$values) & eig$values > 0)
  if (!length(positive)) {
    return(list(theta = matrix(numeric(), nrow = 0L, ncol = model$hyper_dim), metadata = data.frame()))
  }
  keep <- positive[seq_len(min(length(positive), tail_directions))]
  centered <- sweep(theta, 2L, wcov$center, "-")
  rows <- list()
  meta <- list()
  idx <- 0L
  for (direction_id in seq_along(keep)) {
    eig_idx <- keep[direction_id]
    vector <- as.numeric(eig$vectors[, eig_idx])
    projection <- as.numeric(centered %*% vector)
    q <- .local_atlas_weighted_quantile(projection, weights, tail_probs)
    for (prob_id in seq_along(q)) {
      value <- as.numeric(q[prob_id])
      if (!is.finite(value)) next
      idx <- idx + 1L
      theta_new <- matrix(wcov$center + as.numeric(tail_inflation) * value * vector, nrow = 1L)
      colnames(theta_new) <- model$hyper_names
      rows[[idx]] <- theta_new
      meta[[idx]] <- data.frame(
        profile_kind = "outer_posterior_tail",
        tail_direction = direction_id,
        tail_eigenvalue = eig$values[eig_idx],
        tail_probability = tail_probs[prob_id],
        tail_projection = value,
        tail_inflation = as.numeric(tail_inflation),
        check.names = FALSE
      )
    }
  }
  if (!length(rows)) {
    return(list(theta = matrix(numeric(), nrow = 0L, ncol = model$hyper_dim), metadata = data.frame()))
  }
  theta_out <- do.call(rbind, rows)
  colnames(theta_out) <- model$hyper_names
  list(theta = theta_out, metadata = do.call(rbind, meta))
}

.local_evidence_select_outer_fit_rows <- function(theta,
                                                  weights,
                                                  max_points,
                                                  seed = 123L,
                                                  top_weight_fraction = 0.40,
                                                  quantile_probs = c(0.05, 0.5, 0.95),
                                                  max_projection_directions = 4L) {
  weights <- .local_chart_normalize_weights(weights, nrow(theta))
  max_points <- min(as.integer(max_points), nrow(theta))
  if (!is.finite(max_points) || max_points < 1L) {
    stop("max_points must be a positive integer.")
  }
  selected <- integer()
  add_rows <- function(rows) {
    rows <- as.integer(rows)
    rows <- rows[is.finite(rows) & rows >= 1L & rows <= nrow(theta)]
    selected <<- unique(c(selected, rows))
    if (length(selected) > max_points) {
      selected <<- selected[seq_len(max_points)]
    }
  }

  top_n <- max(1L, floor(max_points * as.numeric(top_weight_fraction)))
  add_rows(head(order(weights, decreasing = TRUE), top_n))

  wcov <- .local_evidence_weighted_center_cov(theta, weights)
  eig <- eigen(wcov$covariance, symmetric = TRUE)
  keep <- which(is.finite(eig$values) & eig$values > 0)
  keep <- keep[seq_len(min(length(keep), as.integer(max_projection_directions)))]
  if (length(keep) && length(selected) < max_points) {
    centered <- sweep(theta, 2L, wcov$center, "-")
    quantile_probs <- sort(unique(pmin(pmax(as.numeric(quantile_probs), 0), 1)))
    for (eig_idx in keep) {
      projection <- as.numeric(centered %*% as.numeric(eig$vectors[, eig_idx]))
      q <- .local_atlas_weighted_quantile(projection, weights, quantile_probs)
      add_rows(vapply(q, function(value) which.min(abs(projection - value)), integer(1)))
      if (length(selected) >= max_points) break
    }
  }

  if (length(selected) < max_points) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    } else {
      NULL
    }
    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(as.integer(seed))
    draw_n <- max_points - length(selected)
    draw <- sample.int(nrow(theta), size = max(draw_n * 3L, draw_n), replace = TRUE, prob = weights)
    add_rows(draw)
  }
  if (length(selected) < max_points) {
    add_rows(order(weights, decreasing = TRUE))
  }
  selected[seq_len(min(length(selected), max_points))]
}

build_outer_fit_certification_cloud <- function(fit,
                                                population_model = NULL,
                                                max_posterior_points = 256L,
                                                posterior_weight_mass = 0.90,
                                                include_tail_profiles = TRUE,
                                                tail_probs = c(0.025, 0.975),
                                                tail_directions = 4L,
                                                tail_weight_mass = 1 - posterior_weight_mass,
                                                tail_inflation = 1.25,
                                                seed = 123L,
                                                theta_source = "pilot_outer") {
  model <- normalize_population_model(population_model %||% fit$population_model)
  theta <- .as_hyper_matrix(fit$theta, model$hyper_names, model$hyper_dim)
  weights <- .local_chart_normalize_weights(fit$w, nrow(theta))
  posterior_rows <- .local_evidence_select_outer_fit_rows(
    theta = theta,
    weights = weights,
    max_points = max_posterior_points,
    seed = seed
  )
  posterior_theta <- theta[posterior_rows, , drop = FALSE]
  posterior_weight <- .local_chart_normalize_weights(weights[posterior_rows], length(posterior_rows))

  tail_weight_mass <- if (isTRUE(include_tail_profiles)) {
    min(max(as.numeric(tail_weight_mass), 0), 1)
  } else {
    0
  }
  posterior_weight_mass <- min(max(as.numeric(posterior_weight_mass), 0), 1)
  if (tail_weight_mass + posterior_weight_mass <= 0) {
    posterior_weight_mass <- 1
    tail_weight_mass <- 0
  }
  scale_sum <- posterior_weight_mass + tail_weight_mass
  posterior_weight_mass <- posterior_weight_mass / scale_sum
  tail_weight_mass <- tail_weight_mass / scale_sum

  tail <- if (tail_weight_mass > 0) {
    .local_evidence_outer_fit_tail_profiles(
      theta = theta,
      weights = weights,
      population_model = model,
      tail_probs = tail_probs,
      tail_directions = tail_directions,
      tail_inflation = tail_inflation
    )
  } else {
    list(theta = matrix(numeric(), nrow = 0L, ncol = model$hyper_dim), metadata = data.frame())
  }
  if (!nrow(tail$theta)) {
    tail_weight_mass <- 0
    posterior_weight_mass <- 1
  }

  theta_all <- rbind(posterior_theta, tail$theta)
  weight_all <- c(
    posterior_weight_mass * posterior_weight,
    if (nrow(tail$theta)) rep(tail_weight_mass / nrow(tail$theta), nrow(tail$theta)) else numeric()
  )
  metadata <- data.frame(
    outer_row = c(posterior_rows, rep(NA_integer_, nrow(tail$theta))),
    outer_weight = c(weights[posterior_rows], rep(NA_real_, nrow(tail$theta))),
    profile_kind = c(rep("outer_posterior_particle", nrow(posterior_theta)),
                     as.character(tail$metadata$profile_kind %||% character())),
    tail_direction = c(rep(NA_integer_, nrow(posterior_theta)),
                       as.integer(tail$metadata$tail_direction %||% integer())),
    tail_probability = c(rep(NA_real_, nrow(posterior_theta)),
                         as.numeric(tail$metadata$tail_probability %||% numeric())),
    check.names = FALSE
  )
  build_local_evidence_certification_cloud(
    theta = theta_all,
    population_model = model,
    theta_weights = weight_all,
    theta_source = c(rep(theta_source, nrow(posterior_theta)),
                     rep(paste0(theta_source, "_tail"), nrow(tail$theta))),
    theta_round = 0L,
    theta_id = c(
      sprintf("%s_particle_%06d", theta_source, posterior_rows),
      if (nrow(tail$theta)) sprintf("%s_tail_%03d", theta_source, seq_len(nrow(tail$theta))) else character()
    ),
    metadata = metadata,
    deduplicate = TRUE
  )
}

run_certified_outer_workflow <- function(factor_set,
                                         initial_proposal = NULL,
                                         pre_outer_cloud = NULL,
                                         pre_outer_theta = NULL,
                                         pre_outer_theta_weights = NULL,
                                         data_list = NULL,
                                         loglik_fn = NULL,
                                         local_ids = NULL,
                                         local_control = list(),
                                         edge_control = list(),
                                         pre_outer_control = list(),
                                         pilot_outer_control = list(),
                                         post_pilot_cloud_control = list(),
                                         post_pilot_control = list(),
                                         final_outer_control = list(),
                                         allow_uncertified_pre_outer = FALSE,
                                         allow_uncertified_final_outer = FALSE,
                                         checkpoint_file = NULL,
                                         n_cores = 1L,
                                         seed = 123L,
                                         verbose = TRUE,
                                         trace_verbose = FALSE) {
  factor_set <- validate_local_atlas_factor_set(factor_set)
  model <- factor_set$population_model
  if (!is.null(initial_proposal)) {
    initial_proposal <- normalize_theta_proposal(initial_proposal, population_model = model)
  }
  if (is.null(pre_outer_cloud) && is.null(pre_outer_theta)) {
    if (is.null(initial_proposal)) {
      stop("pre_outer_theta or initial_proposal is required for certified outer workflow.")
    }
    n_pre <- as.integer(pre_outer_control$n_theta %||%
      max(256L, as.integer(pilot_outer_control$N %||% 1000L)))
    pre_outer_theta <- theta_proposal_sample(
      initial_proposal,
      n = n_pre,
      seed = as.integer(seed) + 110001L
    )
    pre_outer_theta_weights <- rep(1 / n_pre, n_pre)
  }
  pre_outer_control <- modifyList(
    list(
      max_rounds = 1L,
      certification_thresholds = list(),
      scoring_control = list(),
      geometry_control = list(),
      repair_control = list(),
      profile_control = list(enabled = TRUE)
    ),
    pre_outer_control
  )
  post_pilot_cloud_control <- modifyList(
    list(
      max_posterior_points = 256L,
      posterior_weight_mass = 0.90,
      include_tail_profiles = TRUE,
      tail_probs = c(0.025, 0.975),
      tail_directions = min(4L, model$hyper_dim),
      tail_weight_mass = 0.10,
      tail_inflation = 1.25
    ),
    post_pilot_cloud_control
  )
  post_pilot_control <- modifyList(
    list(
      max_rounds = 1L,
      certification_thresholds = list(),
      scoring_control = list(),
      geometry_control = list(),
      repair_control = list(),
      profile_control = list(enabled = TRUE)
    ),
    post_pilot_control
  )
  pilot_outer_control <- modifyList(list(N = 500L, max_rounds = 80L), pilot_outer_control)
  final_outer_control <- modifyList(list(N = 2000L, max_rounds = 80L), final_outer_control)
  can_repair <- !is.null(data_list) && is.function(loglik_fn)
  if (!isTRUE(can_repair)) {
    if (isTRUE(allow_uncertified_pre_outer) && as.integer(pre_outer_control$max_rounds) > 0L) {
      pre_outer_control$max_rounds <- 0L
    }
    if (isTRUE(allow_uncertified_final_outer) && as.integer(post_pilot_control$max_rounds) > 0L) {
      post_pilot_control$max_rounds <- 0L
    }
  }

  history <- list()
  status <- "started"

  pre_cert <- certify_and_repair_raw_atlas(
    factor_set = factor_set,
    cloud = pre_outer_cloud,
    theta = pre_outer_theta,
    theta_weights = pre_outer_theta_weights,
    theta_source = "pre_outer",
    data_list = data_list,
    loglik_fn = loglik_fn,
    local_ids = local_ids,
    local_control = local_control,
    edge_control = edge_control,
    calibration_control = pre_outer_control$calibration_control %||% list(),
    certification_thresholds = pre_outer_control$certification_thresholds,
    scoring_control = pre_outer_control$scoring_control,
    geometry_control = pre_outer_control$geometry_control,
    repair_control = pre_outer_control$repair_control,
    profile_control = pre_outer_control$profile_control,
    max_rounds = as.integer(pre_outer_control$max_rounds),
    stop_on_uncertified = FALSE,
    keep_round_tables = isTRUE(pre_outer_control$keep_round_tables),
    n_cores = as.integer(n_cores),
    seed = as.integer(seed) + 210001L,
    verbose = verbose,
    trace_verbose = trace_verbose
  )
  factor_set <- pre_cert$factor_set
  history[[length(history) + 1L]] <- data.frame(
    stage = "pre_outer_certification",
    status = pre_cert$status,
    certified = pre_cert$certified,
    n_theta = nrow(pre_cert$cloud$theta),
    log_evidence = NA_real_,
    mcse_log_evidence = NA_real_,
    check.names = FALSE
  )
  .local_evidence_outer_workflow_checkpoint(
    checkpoint_file,
    "certified_outer_pre_outer",
    list(factor_set = factor_set, pre_outer_certification = pre_cert),
    verbose = verbose
  )
  if (!isTRUE(pre_cert$certified) && !isTRUE(allow_uncertified_pre_outer)) {
    status <- "pre_outer_uncertified"
    out <- structure(
      list(
        factor_set = factor_set,
        pre_outer_certification = pre_cert,
        pilot_fit = NULL,
        post_pilot_cloud = NULL,
        post_pilot_certification = NULL,
        final_fit = NULL,
        status = status,
        history = do.call(rbind, history),
        settings = list(
          allow_uncertified_pre_outer = isTRUE(allow_uncertified_pre_outer),
          allow_uncertified_final_outer = isTRUE(allow_uncertified_final_outer),
          pre_outer_control = pre_outer_control,
          pilot_outer_control = pilot_outer_control,
          post_pilot_cloud_control = post_pilot_cloud_control,
          post_pilot_control = post_pilot_control,
          final_outer_control = final_outer_control,
          seed = as.integer(seed)
        )
      ),
      class = "certified_outer_workflow"
    )
    return(validate_certified_outer_workflow(out))
  }

  .local_atlas_log("starting pilot outer SMC: N=", as.integer(pilot_outer_control$N), "\n", verbose = verbose)
  pilot_fit <- .local_evidence_run_outer(
    factor_set = factor_set,
    initial_proposal = initial_proposal,
    outer_control = pilot_outer_control,
    n_cores = n_cores,
    seed = as.integer(seed) + 310001L,
    verbose = isTRUE(pilot_outer_control$verbose %||% FALSE)
  )
  history[[length(history) + 1L]] <- data.frame(
    stage = "pilot_outer",
    status = if (isTRUE(pilot_fit$beta >= 1 - 1e-12)) "complete" else "incomplete_beta",
    certified = NA,
    n_theta = nrow(pilot_fit$theta),
    log_evidence = as.numeric(pilot_fit$log_evidence),
    mcse_log_evidence = as.numeric(pilot_fit$mcse_log_evidence),
    check.names = FALSE
  )
  .local_evidence_outer_workflow_checkpoint(
    checkpoint_file,
    "certified_outer_pilot",
    list(factor_set = factor_set, pre_outer_certification = pre_cert, pilot_fit = pilot_fit),
    verbose = verbose
  )

  post_cloud <- do.call(
    build_outer_fit_certification_cloud,
    c(list(
      fit = pilot_fit,
      population_model = model,
      seed = as.integer(seed) + 410001L,
      theta_source = "post_pilot"
    ), post_pilot_cloud_control)
  )
  post_cert <- certify_and_repair_raw_atlas(
    factor_set = factor_set,
    cloud = post_cloud,
    data_list = data_list,
    loglik_fn = loglik_fn,
    local_ids = local_ids,
    local_control = local_control,
    edge_control = edge_control,
    calibration_control = post_pilot_control$calibration_control %||% list(),
    certification_thresholds = post_pilot_control$certification_thresholds,
    scoring_control = post_pilot_control$scoring_control,
    geometry_control = post_pilot_control$geometry_control,
    repair_control = post_pilot_control$repair_control,
    profile_control = post_pilot_control$profile_control,
    max_rounds = as.integer(post_pilot_control$max_rounds),
    stop_on_uncertified = FALSE,
    keep_round_tables = isTRUE(post_pilot_control$keep_round_tables),
    n_cores = as.integer(n_cores),
    seed = as.integer(seed) + 510001L,
    verbose = verbose,
    trace_verbose = trace_verbose
  )
  factor_set <- post_cert$factor_set
  history[[length(history) + 1L]] <- data.frame(
    stage = "post_pilot_certification",
    status = post_cert$status,
    certified = post_cert$certified,
    n_theta = nrow(post_cert$cloud$theta),
    log_evidence = NA_real_,
    mcse_log_evidence = NA_real_,
    check.names = FALSE
  )
  .local_evidence_outer_workflow_checkpoint(
    checkpoint_file,
    "certified_outer_post_pilot",
    list(
      factor_set = factor_set,
      pre_outer_certification = pre_cert,
      pilot_fit = pilot_fit,
      post_pilot_cloud = post_cloud,
      post_pilot_certification = post_cert
    ),
    verbose = verbose
  )
  if (!isTRUE(post_cert$certified) && !isTRUE(allow_uncertified_final_outer)) {
    status <- "post_pilot_uncertified"
    out <- structure(
      list(
        factor_set = factor_set,
        pre_outer_certification = pre_cert,
        pilot_fit = pilot_fit,
        post_pilot_cloud = post_cloud,
        post_pilot_certification = post_cert,
        final_fit = NULL,
        status = status,
        history = do.call(rbind, history),
        settings = list(
          allow_uncertified_pre_outer = isTRUE(allow_uncertified_pre_outer),
          allow_uncertified_final_outer = isTRUE(allow_uncertified_final_outer),
          pre_outer_control = pre_outer_control,
          pilot_outer_control = pilot_outer_control,
          post_pilot_cloud_control = post_pilot_cloud_control,
          post_pilot_control = post_pilot_control,
          final_outer_control = final_outer_control,
          seed = as.integer(seed)
        )
      ),
      class = "certified_outer_workflow"
    )
    return(validate_certified_outer_workflow(out))
  }

  .local_atlas_log("starting final outer SMC: N=", as.integer(final_outer_control$N), "\n", verbose = verbose)
  final_fit <- .local_evidence_run_outer(
    factor_set = factor_set,
    initial_proposal = initial_proposal,
    outer_control = final_outer_control,
    n_cores = n_cores,
    seed = as.integer(seed) + 610001L,
    verbose = isTRUE(final_outer_control$verbose %||% FALSE)
  )
  status <- if (isTRUE(final_fit$beta >= 1 - 1e-12)) "complete" else "final_outer_incomplete_beta"
  history[[length(history) + 1L]] <- data.frame(
    stage = "final_outer",
    status = status,
    certified = NA,
    n_theta = nrow(final_fit$theta),
    log_evidence = as.numeric(final_fit$log_evidence),
    mcse_log_evidence = as.numeric(final_fit$mcse_log_evidence),
    check.names = FALSE
  )
  out <- structure(
    list(
      factor_set = factor_set,
      pre_outer_certification = pre_cert,
      pilot_fit = pilot_fit,
      post_pilot_cloud = post_cloud,
      post_pilot_certification = post_cert,
      final_fit = final_fit,
      status = status,
      history = do.call(rbind, history),
      settings = list(
        allow_uncertified_pre_outer = isTRUE(allow_uncertified_pre_outer),
        allow_uncertified_final_outer = isTRUE(allow_uncertified_final_outer),
        pre_outer_control = pre_outer_control,
        pilot_outer_control = pilot_outer_control,
        post_pilot_cloud_control = post_pilot_cloud_control,
        post_pilot_control = post_pilot_control,
        final_outer_control = final_outer_control,
        seed = as.integer(seed)
      )
    ),
    class = "certified_outer_workflow"
  ) |>
    validate_certified_outer_workflow()
  .local_evidence_outer_workflow_checkpoint(
    checkpoint_file,
    "certified_outer_final",
    list(
      factor_set = out$factor_set,
      pre_outer_certification = out$pre_outer_certification,
      pilot_fit = out$pilot_fit,
      post_pilot_cloud = out$post_pilot_cloud,
      post_pilot_certification = out$post_pilot_certification,
      final_fit = out$final_fit,
      workflow = out
    ),
    verbose = verbose
  )
  out
}

.local_atlas_factor_set_by_local <- function(factor_set,
                                             theta,
                                             n_cores = 1L) {
  factor_set <- validate_local_atlas_factor_set(factor_set)
  model <- factor_set$population_model
  theta <- .as_hyper_matrix(theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  control <- factor_set$evaluator_control
  eval_one <- function(atlas) {
    can_batch <- isTRUE(control$use_particle_mis) &&
      isTRUE(control$particle_mis_batch) &&
      !is.finite(control$max_leave_chart_out_gap) &&
      identical(control$particle_mis_role, "estimator")
    if (can_batch) {
      .local_atlas_particle_mis_many_global(
        atlas = atlas,
        theta = theta,
        population_model = model,
        max_chart_distance = control$max_chart_distance,
        min_covering_charts = control$min_covering_charts,
        min_ess_frac = control$min_particle_mis_ess,
        min_ess = control$min_particle_mis_ess_abs,
        max_psis_k = control$max_particle_mis_psis_k,
        sparse_chart_min_covering = control$sparse_chart_min_covering,
        sparse_chart_max_distance = control$sparse_chart_max_distance,
        distance_metric = control$distance_metric,
        se_floor = control$se_floor,
        use_compressed_particle_mis = control$use_compressed_particle_mis,
        require_compressed_particle_mis = control$require_compressed_particle_mis,
        use_uncertified_estimates = control$use_uncertified_estimates
      )
    } else {
      evaluate_local_atlas_many(
        atlas = atlas,
        theta = theta,
        population_model = model,
        max_chart_distance = control$max_chart_distance,
        min_covering_charts = control$min_covering_charts,
        max_prediction_range = control$max_prediction_range,
        distance_scale = control$distance_scale,
        se_floor = control$se_floor,
        use_particle_mis = control$use_particle_mis,
        require_particle_mis = control$require_particle_mis,
        min_particle_mis_ess = control$min_particle_mis_ess,
        min_particle_mis_ess_abs = control$min_particle_mis_ess_abs,
        max_particle_mis_psis_k = control$max_particle_mis_psis_k,
        max_quadratic_particle_gap = control$max_quadratic_particle_gap,
        sparse_chart_min_covering = control$sparse_chart_min_covering,
        sparse_chart_max_distance = control$sparse_chart_max_distance,
        max_leave_chart_out_gap = control$max_leave_chart_out_gap,
        distance_metric = control$distance_metric,
        surface_method = control$surface_method,
        min_surface_charts = control$min_surface_charts,
        max_surface_se = control$max_surface_se,
        surface_value_nugget = control$surface_value_nugget,
        surface_gradient_weight = control$surface_gradient_weight,
        surface_curvature_weight = control$surface_curvature_weight,
        surface_ridge = control$surface_ridge,
        particle_mis_role = control$particle_mis_role
      )
    }
  }
  parts <- if (as.integer(n_cores) <= 1L || length(factor_set$atlases) <= 1L) {
    lapply(factor_set$atlases, eval_one)
  } else {
    parallel::mclapply(
      factor_set$atlases,
      eval_one,
      mc.cores = as.integer(min(n_cores, length(factor_set$atlases)))
    )
  }
  names(parts) <- names(factor_set$atlases)
  parts
}

.local_atlas_factor_set_uncertified <- function(parts) {
  rows <- list()
  for (local_name in names(parts)) {
    part <- parts[[local_name]]
    bad <- which(part$status != "certified")
    if (!length(bad)) next
    rows[[length(rows) + 1L]] <- data.frame(
      local = local_name,
      theta_row = bad,
      status = part$status[bad],
      reason = part$reason[bad],
      nearest_charts = part$nearest_charts[bad],
      check.names = FALSE
    )
  }
  if (length(rows)) do.call(rbind, rows) else data.frame()
}

local_evidence_atlas_factor_set_loglik_by_local <- function(factor_set,
                                                            theta,
                                                            include_constant = TRUE,
                                                            n_cores = 1L) {
  factor_set <- validate_local_atlas_factor_set(factor_set)
  parts <- .local_atlas_factor_set_by_local(factor_set, theta, n_cores = n_cores)
  bad <- .local_atlas_factor_set_uncertified(parts)
  if (nrow(bad) && isTRUE(factor_set$evaluator_control$stop_on_uncertified)) {
    shown <- utils::head(bad, 8L)
    stop(
      "Atlas factor set evaluated uncertified local evidence:\n",
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
  out <- do.call(cbind, lapply(parts, `[[`, "log_marginal"))
  colnames(out) <- names(parts)
  if (isTRUE(include_constant) && ncol(out)) {
    out <- out + factor_set$log_constant / ncol(out)
  }
  out
}

local_evidence_atlas_factor_set_loglik <- function(factor_set,
                                                   theta,
                                                   include_constant = FALSE,
                                                   n_cores = 1L) {
  factor_set <- validate_local_atlas_factor_set(factor_set)
  by_local <- local_evidence_atlas_factor_set_loglik_by_local(
    factor_set = factor_set,
    theta = theta,
    include_constant = FALSE,
    n_cores = n_cores
  )
  out <- rowSums(by_local)
  if (isTRUE(include_constant)) {
    out <- out + factor_set$log_constant
  }
  as.numeric(out)
}

.local_atlas_merge_control <- function(user, defaults) {
  modifyList(defaults, user %||% list())
}

.local_atlas_log <- function(..., verbose = TRUE) {
  if (!isTRUE(verbose)) {
    return(invisible(NULL))
  }
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ...)
  flush.console()
  invisible(NULL)
}

.local_atlas_write_checkpoint <- function(path, stage, state, verbose = TRUE) {
  if (is.null(path) || !nzchar(path)) {
    return(invisible(FALSE))
  }
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  state$stage <- as.character(stage)
  state$checkpoint_time <- Sys.time()
  state$pid <- Sys.getpid()
  tmp <- paste0(path, ".tmp")
  saveRDS(state, tmp)
  ok <- file.rename(tmp, path)
  if (!isTRUE(ok)) {
    file.copy(tmp, path, overwrite = TRUE)
    unlink(tmp)
  }
  stage_path <- sub("\\.rds$", paste0("_", gsub("[^A-Za-z0-9_.-]+", "_", as.character(stage)), ".rds"), path)
  if (!identical(stage_path, path)) {
    stage_tmp <- paste0(stage_path, ".tmp")
    saveRDS(state, stage_tmp)
    stage_ok <- file.rename(stage_tmp, stage_path)
    if (!isTRUE(stage_ok)) {
      file.copy(stage_tmp, stage_path, overwrite = TRUE)
      unlink(stage_tmp)
    }
  }
  .local_atlas_log("checkpoint saved: ", path, " stage=", stage, "\n", verbose = verbose)
  invisible(TRUE)
}

.local_atlas_rbind_fill <- function(rows) {
  if (is.data.frame(rows)) {
    return(rows)
  }
  rows <- Filter(function(x) is.data.frame(x) && nrow(x), rows %||% list())
  if (!length(rows)) {
    return(data.frame())
  }
  cols <- unique(unlist(lapply(rows, names), use.names = FALSE))
  rows <- lapply(rows, function(x) {
    missing <- setdiff(cols, names(x))
    for (nm in missing) {
      x[[nm]] <- NA
    }
    x[, cols, drop = FALSE]
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.local_atlas_default_local_control <- function() {
  list(
    root_M = 1000L,
    candidate_M = 1000L,
    target_cess = 0.9,
    resample_threshold = 0.5,
    n_mcmc_moves = 2L,
    rw_scale = 0.75,
    G_mix = 8L,
    da_enable = TRUE,
    refit_every = 2L,
    max_steps = 128L,
    deterministic_resampling = FALSE,
    n_cores_per_local = 1L,
    root_confirm = "auto",
    root_confirmation_M = NULL,
    root_confirmation_abs_tol = 1.0
  )
}

.local_atlas_default_edge_control <- function() {
  list(
    method = "bridge",
    edge_neighbors = 2L,
    max_intermediates = 2L,
    max_forward_reverse_gap = 1.0,
    max_se = 1.0,
    min_overlap_ess = 0.05,
    max_psis_k = 0.7,
    max_taylor_gap = 2.0,
    max_taylor_disagreement = Inf,
    max_graph_standardized_residual = 3,
    distance_metric = "euclidean",
    require_bar_converged = TRUE,
    bar_max_iter = 200L,
    bar_tol = 1e-8
  )
}

.local_atlas_default_theta <- function(population_model, n_prior = 2000L, seed = 1L) {
  model <- normalize_population_model(population_model)
  if (identical(model$fast_family %||% NA_character_, "gaussian") &&
      !is.null(model$prior_spec)) {
    d <- model$alpha_dim
    spec <- model$prior_spec
    mu <- as.numeric(spec$mean_prior_mean)
    sigma2 <- as.numeric(spec$sigma2_prior_rate) / pmax(as.numeric(spec$sigma2_prior_shape) - 1, 1e-8)
    out <- matrix(c(mu, log(pmax(sigma2, 1e-12))), nrow = 1L)
    colnames(out) <- model$hyper_names
    return(out)
  }
  theta <- population_model_sample_hyper(model, n = as.integer(n_prior), seed = seed)
  out <- matrix(apply(theta, 2L, stats::median), nrow = 1L)
  colnames(out) <- model$hyper_names
  out
}

local_atlas_theta_from_draws <- function(draws, population_model) {
  model <- normalize_population_model(population_model)
  draws <- as.data.frame(draws, check.names = FALSE)
  if (all(model$hyper_names %in% names(draws))) {
    return(.as_hyper_matrix(draws[, model$hyper_names, drop = FALSE], model$hyper_names, model$hyper_dim))
  }
  if (!identical(model$fast_family %||% NA_character_, "gaussian")) {
    stop("draw-to-theta conversion currently supports only diagonal Gaussian population models.")
  }
  mu_names <- paste0("mu_", model$alpha_names)
  sigma2_names <- paste0("sigma2_", model$alpha_names)
  missing <- setdiff(c(mu_names, sigma2_names), names(draws))
  if (length(missing)) {
    stop("draws are missing columns needed for theta conversion: ", paste(missing, collapse = ", "))
  }
  theta <- cbind(
    as.matrix(draws[, mu_names, drop = FALSE]),
    log(pmax(as.matrix(draws[, sigma2_names, drop = FALSE]), 1e-12))
  )
  colnames(theta) <- model$hyper_names
  .as_hyper_matrix(theta, model$hyper_names, model$hyper_dim)
}

.local_atlas_unique_theta <- function(theta, population_model, digits = 10L) {
  model <- normalize_population_model(population_model)
  theta <- .as_hyper_matrix(theta, model$hyper_names, model$hyper_dim)
  key <- as.data.frame(signif(theta, digits = as.integer(digits)), check.names = FALSE)
  theta[!duplicated(key), , drop = FALSE]
}

.local_atlas_weighted_quantile <- function(x, weights, probs) {
  x <- as.numeric(x)
  weights <- .local_chart_normalize_weights(weights, length(x))
  ok <- is.finite(x) & is.finite(weights) & weights > 0
  x <- x[ok]
  weights <- weights[ok]
  if (!length(x)) {
    return(rep(NA_real_, length(probs)))
  }
  ord <- order(x)
  x <- x[ord]
  weights <- weights[ord] / sum(weights[ord])
  cdf <- cumsum(weights)
  vapply(as.numeric(probs), function(p) {
    x[which(cdf >= min(max(p, 0), 1))[1L]]
  }, numeric(1))
}

.local_atlas_distance_metric <- function(distance_metric) {
  distance_metric <- as.character(distance_metric %||% "euclidean")
  if (!length(distance_metric) || is.na(distance_metric[1L]) || !nzchar(distance_metric[1L])) {
    return("euclidean")
  }
  distance_metric <- distance_metric[1L]
  allowed <- c("euclidean", "fisher")
  if (!distance_metric %in% allowed) {
    stop("distance_metric must be one of: ", paste(allowed, collapse = ", "))
  }
  distance_metric
}

.local_atlas_theta_distances <- function(theta,
                                         centers,
                                         population_model,
                                         distance_metric = "euclidean") {
  model <- normalize_population_model(population_model)
  theta <- .as_hyper_matrix(theta, model$hyper_names, model$hyper_dim)
  centers <- .as_hyper_matrix(centers, model$hyper_names, model$hyper_dim)
  distance_metric <- .local_atlas_distance_metric(distance_metric)
  out <- matrix(NA_real_, nrow = nrow(theta), ncol = nrow(centers))
  if (!nrow(theta) || !nrow(centers)) {
    return(out)
  }

  if (identical(distance_metric, "fisher") &&
      identical(model$fast_family %||% NA_character_, "gaussian") &&
      model$hyper_dim == 2L * model$alpha_dim) {
    d <- model$alpha_dim
    mu_idx <- seq_len(d)
    rho_idx <- d + seq_len(d)
    theta_mu <- theta[, mu_idx, drop = FALSE]
    theta_rho <- theta[, rho_idx, drop = FALSE]
    theta_inv_var <- exp(-theta_rho)
    for (j in seq_len(nrow(centers))) {
      center_mu <- matrix(centers[j, mu_idx], nrow = nrow(theta), ncol = d, byrow = TRUE)
      center_rho <- matrix(centers[j, rho_idx], nrow = nrow(theta), ncol = d, byrow = TRUE)
      center_inv_var <- matrix(exp(-centers[j, rho_idx]), nrow = nrow(theta), ncol = d, byrow = TRUE)
      dmu <- theta_mu - center_mu
      drho <- theta_rho - center_rho
      mean_metric <- 0.5 * (theta_inv_var + center_inv_var)
      out[, j] <- sqrt(rowSums(dmu * dmu * mean_metric) + 0.5 * rowSums(drho * drho))
    }
    return(out)
  }

  for (j in seq_len(nrow(centers))) {
    out[, j] <- sqrt(rowSums(sweep(theta, 2L, as.numeric(centers[j, ]), "-")^2))
  }
  out
}

.local_atlas_metric_farthest_design <- function(theta_cloud,
                                                population_model,
                                                theta_root,
                                                max_anchors,
                                                distance_metric = "fisher",
                                                initial_design = NULL) {
  model <- normalize_population_model(population_model)
  theta_cloud <- .local_atlas_unique_theta(theta_cloud, model)
  theta_root <- .local_chart_align_theta_one(theta_root, model)
  max_anchors <- as.integer(max_anchors)
  if (max_anchors <= 1L || nrow(theta_cloud) <= 1L) {
    return(.local_atlas_unique_theta(theta_root, model))
  }
  selected <- .local_atlas_unique_theta(rbind(theta_root, initial_design %||% theta_root), model)
  selected <- selected[seq_len(min(nrow(selected), max_anchors)), , drop = FALSE]

  while (nrow(selected) < min(max_anchors, nrow(theta_cloud))) {
    distances <- .local_atlas_theta_distances(
      theta = theta_cloud,
      centers = selected,
      population_model = model,
      distance_metric = distance_metric
    )
    nearest <- matrixStats::rowMins(distances)
    duplicate <- apply(theta_cloud, 1L, function(x) {
      any(rowSums(sweep(selected, 2L, x, "-")^2) <= 1e-20)
    })
    nearest[duplicate] <- -Inf
    if (!any(is.finite(nearest))) break
    selected <- rbind(selected, theta_cloud[which.max(nearest), , drop = FALSE])
    selected <- .local_atlas_unique_theta(selected, model)
  }
  selected[seq_len(min(nrow(selected), max_anchors)), , drop = FALSE]
}

local_atlas_select_theta_profiles <- function(theta,
                                              population_model,
                                              weights = NULL,
                                              focus_hyper_names = NULL,
                                              probs = c(0.05, 0.5, 0.95),
                                              max_points = 9L,
                                              return_rows = FALSE) {
  model <- normalize_population_model(population_model)
  theta <- .as_hyper_matrix(theta, model$hyper_names, model$hyper_dim)
  weights <- .local_chart_normalize_weights(weights, nrow(theta))
  focus_hyper_names <- focus_hyper_names %||% model$hyper_names
  focus_idx <- match(focus_hyper_names, model$hyper_names)
  focus_idx <- focus_idx[is.finite(focus_idx)]
  center <- vapply(seq_len(ncol(theta)), function(j) {
    .local_atlas_weighted_quantile(theta[, j], weights, 0.5)
  }, numeric(1))
  rows <- which.min(rowSums(sweep(theta, 2L, center, "-")^2))
  for (j in focus_idx) {
    qj <- .local_atlas_weighted_quantile(theta[, j], weights, probs)
    for (target in qj[is.finite(qj)]) {
      rows <- c(rows, which.min(abs(theta[, j] - target)))
    }
  }
  rows <- unique(rows)
  rows <- rows[seq_len(min(length(rows), as.integer(max_points)))]
  out <- .local_atlas_unique_theta(theta[rows, , drop = FALSE], model)
  if (!isTRUE(return_rows)) {
    return(out)
  }
  matched_rows <- vapply(seq_len(nrow(out)), function(i) {
    which.min(rowSums(sweep(theta, 2L, out[i, ], "-")^2))
  }, integer(1))
  list(
    theta = out,
    rows = matched_rows,
    weights = .local_chart_normalize_weights(weights[matched_rows], length(matched_rows))
  )
}

.local_atlas_design_from_cloud <- function(theta_cloud,
                                           population_model,
                                           theta_root = NULL,
                                           max_anchors = 9L,
                                           axis_count = 4L,
                                           tail_probs = c(0.05, 0.5, 0.95),
                                           focus_hyper_names = NULL,
                                           include_axis_profiles = TRUE,
                                           include_cloud_profiles = TRUE,
                                           design_method = "profiles",
                                           distance_metric = "euclidean") {
  model <- normalize_population_model(population_model)
  theta_cloud <- .as_hyper_matrix(theta_cloud, model$hyper_names, model$hyper_dim)
  theta_root <- .as_hyper_matrix(theta_root %||% theta_cloud[1L, , drop = FALSE], model$hyper_names, model$hyper_dim)
  design_method <- match.arg(
    as.character(design_method %||% "profiles"),
    choices = c("profiles", "metric_farthest", "hybrid")
  )
  distance_metric <- .local_atlas_distance_metric(distance_metric)
  center <- as.numeric(theta_root[1L, ])
  names(center) <- model$hyper_names
  spread <- apply(theta_cloud, 2L, stats::sd)
  spread[!is.finite(spread)] <- 0
  focus_idx <- match(focus_hyper_names %||% character(), model$hyper_names)
  focus_idx <- focus_idx[is.finite(focus_idx)]
  ranked <- unique(c(focus_idx, order(spread, decreasing = TRUE)))
  ranked <- ranked[seq_len(min(length(ranked), as.integer(axis_count)))]
  tail_probs <- sort(unique(pmin(pmax(as.numeric(tail_probs), 0), 1)))

  cloud_design <- list()
  if (isTRUE(include_cloud_profiles)) {
    for (j in ranked) {
      qj <- stats::quantile(theta_cloud[, j], probs = tail_probs, names = FALSE, type = 8)
      for (target in qj) {
        idx <- which.min(abs(theta_cloud[, j] - target))
        cloud_design[[length(cloud_design) + 1L]] <- theta_cloud[idx, , drop = FALSE]
      }
    }
  }
  axis_design <- list()
  if (isTRUE(include_axis_profiles)) {
    for (j in ranked) {
      qj <- stats::quantile(theta_cloud[, j], probs = tail_probs, names = FALSE, type = 8)
      for (target in qj) {
        theta <- matrix(center, nrow = 1L, dimnames = list(NULL, model$hyper_names))
        theta[1L, j] <- target
        axis_design[[length(axis_design) + 1L]] <- theta
      }
    }
  }
  profile_parts <- c(list(theta_root), axis_design, cloud_design)
  profile_design <- .local_atlas_unique_theta(do.call(rbind, profile_parts), model)
  if (identical(design_method, "metric_farthest")) {
    out <- .local_atlas_metric_farthest_design(
      theta_cloud = theta_cloud,
      population_model = model,
      theta_root = theta_root,
      max_anchors = as.integer(max_anchors),
      distance_metric = distance_metric
    )
  } else if (identical(design_method, "hybrid")) {
    seed_design <- .local_atlas_unique_theta(do.call(rbind, c(list(theta_root), axis_design)), model)
    seed_n <- min(nrow(seed_design), max(1L, floor(as.integer(max_anchors) / 2L)))
    out <- .local_atlas_metric_farthest_design(
      theta_cloud = theta_cloud,
      population_model = model,
      theta_root = theta_root,
      max_anchors = as.integer(max_anchors),
      distance_metric = distance_metric,
      initial_design = seed_design[seq_len(seed_n), , drop = FALSE]
    )
  } else {
    out <- profile_design
    out <- out[seq_len(min(nrow(out), as.integer(max_anchors))), , drop = FALSE]
  }
  attr(out, "ranked_hyper_names") <- model$hyper_names[ranked]
  attr(out, "design_method") <- design_method
  attr(out, "distance_metric") <- distance_metric
  out
}

local_atlas_hyperparameter_families <- function(population_model) {
  model <- normalize_population_model(population_model)
  if (!identical(model$fast_family %||% NA_character_, "gaussian") ||
      model$hyper_dim != 2L * model$alpha_dim) {
    stop("adaptive chart design currently requires the diagonal Gaussian population model.")
  }
  rows <- lapply(seq_len(model$alpha_dim), function(j) {
    hyper_names <- c(model$hyper_names[j], model$hyper_names[model$alpha_dim + j])
    data.frame(
      family = model$alpha_names[j],
      alpha_name = model$alpha_names[j],
      mu_hyper = hyper_names[1L],
      log_sigma2_hyper = hyper_names[2L],
      hyper_names = paste(hyper_names, collapse = ","),
      check.names = FALSE
    )
  })
  do.call(rbind, rows)
}

.local_atlas_weighted_sd <- function(x, weights) {
  x <- as.numeric(x)
  weights <- .local_chart_normalize_weights(weights, length(x))
  ok <- is.finite(x) & is.finite(weights) & weights > 0
  if (!any(ok)) {
    return(0)
  }
  x <- x[ok]
  weights <- weights[ok] / sum(weights[ok])
  center <- sum(weights * x)
  sqrt(sum(weights * (x - center)^2))
}

.local_atlas_root_scout_atlas <- function(local_id,
                                          data_i,
                                          loglik_fn,
                                          population_model,
                                          theta_root,
                                          local_control,
                                          seed,
                                          verbose = FALSE) {
  model <- normalize_population_model(population_model)
  root <- build_local_root_chart(
    local_id = local_id,
    theta_anchor = theta_root,
    data_i = data_i,
    loglik_fn = loglik_fn,
    population_model = model,
    chart_id = "root",
    M = as.integer(local_control$root_M),
    target_cess = local_control$target_cess,
    resample_threshold = local_control$resample_threshold,
    n_mcmc_moves = as.integer(local_control$n_mcmc_moves),
    rw_scale = local_control$rw_scale,
    G_mix = as.integer(local_control$G_mix),
    da_enable = isTRUE(local_control$da_enable),
    refit_every = as.integer(local_control$refit_every),
    max_steps = as.integer(local_control$max_steps),
    deterministic_resampling = isTRUE(local_control$deterministic_resampling),
    n_cores = as.integer(local_control$n_cores_per_local),
    seed = as.integer(seed),
    verbose = isTRUE(verbose),
    confirm = local_control$root_confirm,
    confirmation_M = as.integer(local_control$root_confirmation_M %||% max(100L, ceiling(as.integer(local_control$root_M) / 2L))),
    confirmation_abs_tol = local_control$root_confirmation_abs_tol
  )
  atlas <- new_local_atlas(
    local_id = local_id,
    root_chart_id = "root",
    charts = list(root),
    edges = list()
  )
  solve_atlas_normalizers(atlas)
}

validate_chart_design_plan <- function(plan) {
  if (!inherits(plan, "chart_design_plan")) {
    stop("plan must inherit from 'chart_design_plan'.")
  }
  required <- c(
    "theta_root",
    "theta_cloud",
    "families",
    "family_scores",
    "local_family_scores",
    "selected_families",
    "selected_hyper_names",
    "local_anchor_budgets",
    "candidate_theta",
    "local_theta_designs",
    "root_atlases",
    "diagnostics"
  )
  missing <- setdiff(required, names(plan))
  if (length(missing)) {
    stop("chart_design_plan is missing: ", paste(missing, collapse = ", "))
  }
  if (!is.data.frame(plan$family_scores) || !nrow(plan$family_scores)) {
    stop("chart_design_plan family_scores must be a non-empty data frame.")
  }
  if (!is.data.frame(plan$local_family_scores) || !nrow(plan$local_family_scores)) {
    stop("chart_design_plan local_family_scores must be a non-empty data frame.")
  }
  if (!is.list(plan$local_theta_designs) || !length(plan$local_theta_designs)) {
    stop("chart_design_plan local_theta_designs must be a non-empty list.")
  }
  plan
}

chart_design_plan_theta_union <- function(plan,
                                          population_model = NULL) {
  plan <- validate_chart_design_plan(plan)
  theta <- do.call(rbind, plan$local_theta_designs)
  if (is.null(population_model)) {
    theta <- .as_hyper_matrix(theta, hyper_names = colnames(plan$theta_root), hyper_dim = ncol(plan$theta_root))
    key <- apply(signif(theta, 14L), 1L, paste, collapse = "\r")
    return(theta[!duplicated(key), , drop = FALSE])
  }
  .local_atlas_unique_theta(theta, normalize_population_model(population_model))
}

new_chart_design_plan <- function(theta_root,
                                  theta_cloud,
                                  theta_weights,
                                  families,
                                  family_scores,
                                  local_scores,
                                  local_family_scores,
                                  selected_families,
                                  selected_hyper_names,
                                  local_anchor_budgets,
                                  candidate_theta,
                                  local_theta_designs,
                                  root_atlases,
                                  diagnostics = list(),
                                  population_model = NULL) {
  structure(
    list(
      theta_root = theta_root,
      theta_cloud = theta_cloud,
      theta_weights = theta_weights,
      families = families,
      family_scores = family_scores,
      local_scores = local_scores,
      local_family_scores = local_family_scores,
      selected_families = selected_families,
      selected_hyper_names = selected_hyper_names,
      local_anchor_budgets = local_anchor_budgets,
      candidate_theta = candidate_theta,
      local_theta_designs = local_theta_designs,
      root_atlases = root_atlases,
      diagnostics = diagnostics
    ),
    class = "chart_design_plan",
    population_model = population_model
  ) |>
    validate_chart_design_plan()
}

build_chart_design_plan <- function(data_list,
                                    loglik_fn,
                                    population_model,
                                    theta_root,
                                    theta_cloud,
                                    theta_weights = NULL,
                                    local_control = list(),
                                    design_control = list(),
                                    n_cores = 1L,
                                    seed = 123L,
                                    verbose = TRUE,
                                    trace_verbose = FALSE) {
  model <- normalize_population_model(population_model)
  if (!is.list(data_list) || !length(data_list)) {
    stop("data_list must be a non-empty list.")
  }
  if (!is.function(loglik_fn)) {
    stop("loglik_fn must be a function.")
  }
  theta_root <- .local_chart_align_theta_one(theta_root, model)
  theta_cloud_raw <- .as_hyper_matrix(theta_cloud, model$hyper_names, model$hyper_dim)
  theta_weights_raw <- theta_weights
  theta_cloud <- .local_atlas_unique_theta(rbind(theta_root, theta_cloud_raw), model)
  if (is.null(theta_weights_raw)) {
    theta_weights <- rep(1 / nrow(theta_cloud), nrow(theta_cloud))
  } else if (length(theta_weights_raw) == nrow(theta_cloud)) {
    theta_weights <- .local_chart_normalize_weights(theta_weights_raw, nrow(theta_cloud))
  } else if (length(theta_weights_raw) == nrow(theta_cloud_raw)) {
    raw_weights <- .local_chart_normalize_weights(theta_weights_raw, nrow(theta_cloud_raw))
    key <- function(x) apply(signif(as.matrix(x), 14L), 1L, paste, collapse = "\r")
    raw_key <- key(theta_cloud_raw)
    cloud_key <- key(theta_cloud)
    theta_weights <- rowsum(raw_weights, raw_key, reorder = FALSE)
    theta_weights <- as.numeric(theta_weights[match(cloud_key, rownames(theta_weights)), , drop = TRUE])
    theta_weights[!is.finite(theta_weights)] <- 0
    if (sum(theta_weights) <= 0) {
      theta_weights <- rep(1 / nrow(theta_cloud), nrow(theta_cloud))
    } else {
      theta_weights <- theta_weights / sum(theta_weights)
    }
  } else {
    stop("theta_weights must match theta_cloud before or after root insertion.")
  }
  local_control <- .local_atlas_merge_control(local_control, .local_atlas_default_local_control())
  local_control$n_cores_per_local <- 1L
  design_control <- .local_atlas_merge_control(design_control, list(
    tail_probs = c(0.05, 0.5, 0.95),
    distance_metric = "fisher",
    max_anchors = 9L,
    scout_min_families = 1L,
    scout_max_families = 4L,
    scout_relative_threshold = 0.25,
    scout_logz_se_weight = 1,
    scout_weak_path_weight = 1,
    scout_min_anchors_per_local = 3L,
    scout_mean_anchors_per_local = 6L,
    scout_max_anchors_per_local = NULL,
    scout_candidate_pool_size = NULL
  ))

  idx <- seq_along(data_list)
  local_names <- names(data_list) %||% as.character(idx)
  .local_atlas_log(
    "building root scout charts: locals=", length(data_list),
    " root_M=", as.integer(local_control$root_M), "\n",
    verbose = verbose
  )
  build_one <- function(i) {
    if (isTRUE(verbose)) {
      cat(sprintf(
        "Building root scout for local %s (%d/%d)\n",
        local_names[i],
        i,
        length(data_list)
      ))
    }
    .local_atlas_root_scout_atlas(
      local_id = local_names[i],
      data_i = data_list[[i]],
      loglik_fn = loglik_fn,
      population_model = model,
      theta_root = theta_root,
      local_control = local_control,
      seed = as.integer(seed) + 100003L * i,
      verbose = isTRUE(trace_verbose)
    )
  }
  root_atlases <- if (as.integer(n_cores) <= 1L || length(idx) <= 1L) {
    lapply(idx, build_one)
  } else {
    parallel::mclapply(idx, build_one, mc.cores = as.integer(min(n_cores, length(idx))))
  }
  names(root_atlases) <- local_names

  families <- local_atlas_hyperparameter_families(model)
  family_hyper <- strsplit(families$hyper_names, ",", fixed = TRUE)
  rows <- list()
  local_rows <- list()
  theta_center <- as.numeric(theta_root[1L, ])
  names(theta_center) <- model$hyper_names
  for (local_pos in seq_along(root_atlases)) {
    atlas <- validate_local_atlas(root_atlases[[local_pos]])
    root <- atlas$charts[[atlas$root_chart_id]]
    g <- as.numeric(root$score)
    names(g) <- model$hyper_names
    H <- as.matrix(root$curvature)
    dimnames(H) <- list(model$hyper_names, model$hyper_names)
    logZ_se <- as.numeric(root$logZ_abs_se)
    weak_path <- isTRUE(root$diagnostics$weak_path)
    reliability <- 1 +
      as.numeric(design_control$scout_logz_se_weight) * max(logZ_se, 0, na.rm = TRUE) +
      as.numeric(design_control$scout_weak_path_weight) * as.numeric(weak_path)
    if (!is.finite(reliability) || reliability <= 0) {
      reliability <- 1
    }
    total_risk <- 0
    for (f in seq_len(nrow(families))) {
      hyper_names <- family_hyper[[f]]
      hyper_idx <- match(hyper_names, model$hyper_names)
      delta <- sweep(theta_cloud[, hyper_idx, drop = FALSE], 2L, theta_center[hyper_idx], "-")
      Hf <- H[hyper_idx, hyper_idx, drop = FALSE]
      local_taylor <- as.numeric(delta %*% g[hyper_idx]) +
        0.5 * rowSums((delta %*% Hf) * delta)
      variation <- .local_atlas_weighted_sd(local_taylor, theta_weights)
      risk <- max(variation, 0) * reliability
      total_risk <- total_risk + risk
      rows[[length(rows) + 1L]] <- data.frame(
        local = local_names[local_pos],
        local_pos = local_pos,
        family = families$family[f],
        alpha_name = families$alpha_name[f],
        hyper_names = paste(hyper_names, collapse = ","),
        taylor_variation = variation,
        reliability = reliability,
        logZ_se = logZ_se,
        weak_path = weak_path,
        risk = risk,
        check.names = FALSE
      )
    }
    local_rows[[length(local_rows) + 1L]] <- data.frame(
      local = local_names[local_pos],
      local_pos = local_pos,
      total_risk = total_risk,
      logZ_se = logZ_se,
      weak_path = weak_path,
      reliability = reliability,
      check.names = FALSE
    )
  }
  local_family_scores <- do.call(rbind, rows)
  local_scores <- do.call(rbind, local_rows)
  family_scores <- stats::aggregate(
    risk ~ family + alpha_name + hyper_names,
    data = local_family_scores,
    FUN = sum
  )
  names(family_scores)[names(family_scores) == "risk"] <- "total_risk"
  family_scores <- family_scores[order(family_scores$total_risk, decreasing = TRUE), , drop = FALSE]
  max_risk <- max(family_scores$total_risk, na.rm = TRUE)
  if (!is.finite(max_risk) || max_risk <= 0) {
    max_risk <- 0
  }
  family_scores$relative_risk <- if (max_risk > 0) family_scores$total_risk / max_risk else 0
  min_families <- max(1L, as.integer(design_control$scout_min_families))
  max_families <- max(min_families, as.integer(design_control$scout_max_families))
  eligible <- family_scores$relative_risk >= as.numeric(design_control$scout_relative_threshold)
  if (sum(eligible) < min_families) {
    eligible[seq_len(min(min_families, nrow(family_scores)))] <- TRUE
  }
  selected_idx <- which(eligible)
  selected_idx <- selected_idx[seq_len(min(length(selected_idx), max_families))]
  family_scores$selected <- seq_len(nrow(family_scores)) %in% selected_idx
  selected_families <- family_scores$family[selected_idx]
  selected_hyper_names <- unlist(strsplit(family_scores$hyper_names[selected_idx], ",", fixed = TRUE), use.names = FALSE)
  selected_hyper_names <- unique(selected_hyper_names)

  local_selected <- local_family_scores[local_family_scores$family %in% selected_families, , drop = FALSE]
  selected_local_risk <- stats::aggregate(risk ~ local + local_pos, data = local_selected, FUN = sum)
  names(selected_local_risk)[names(selected_local_risk) == "risk"] <- "selected_risk"
  local_scores <- merge(local_scores, selected_local_risk, by = c("local", "local_pos"), all.x = TRUE, sort = FALSE)
  local_scores$selected_risk[!is.finite(local_scores$selected_risk)] <- 0

  n_locals <- nrow(local_scores)
  max_anchors <- max(1L, as.integer(design_control$max_anchors))
  min_budget <- max(1L, min(max_anchors, as.integer(design_control$scout_min_anchors_per_local)))
  mean_budget <- max(min_budget, min(max_anchors, as.integer(design_control$scout_mean_anchors_per_local)))
  max_budget <- as.integer(design_control$scout_max_anchors_per_local %||% max_anchors)
  max_budget <- max(min_budget, min(max_anchors, max_budget))
  total_budget <- max(n_locals * min_budget, min(n_locals * max_budget, n_locals * mean_budget))
  extra_total <- max(0L, total_budget - n_locals * min_budget)
  risk_weight <- local_scores$selected_risk
  if (!any(is.finite(risk_weight)) || sum(risk_weight, na.rm = TRUE) <= 0) {
    risk_weight <- rep(1, n_locals)
  }
  risk_weight[!is.finite(risk_weight) | risk_weight < 0] <- 0
  risk_weight <- risk_weight / sum(risk_weight)
  extra <- floor(extra_total * risk_weight)
  remainder <- extra_total - sum(extra)
  if (remainder > 0L) {
    add_order <- order(extra_total * risk_weight - extra, decreasing = TRUE)
    extra[add_order[seq_len(remainder)]] <- extra[add_order[seq_len(remainder)]] + 1L
  }
  budget <- pmin(max_budget, min_budget + extra)
  while (sum(budget) < total_budget && any(budget < max_budget)) {
    add_order <- order(risk_weight, decreasing = TRUE)
    for (j in add_order) {
      if (budget[j] < max_budget) {
        budget[j] <- budget[j] + 1L
        break
      }
    }
  }
  local_anchor_budgets <- data.frame(
    local = local_scores$local,
    local_pos = local_scores$local_pos,
    selected_risk = local_scores$selected_risk,
    anchor_budget = as.integer(budget),
    check.names = FALSE
  )

  candidate_pool_size <- as.integer(design_control$scout_candidate_pool_size %||%
    max(max_anchors, 1L + length(selected_hyper_names) * length(design_control$tail_probs) * 2L))
  candidate_theta <- .local_atlas_design_from_cloud(
    theta_cloud = theta_cloud,
    population_model = model,
    theta_root = theta_root,
    max_anchors = candidate_pool_size,
    axis_count = length(selected_hyper_names),
    tail_probs = design_control$tail_probs,
    focus_hyper_names = selected_hyper_names,
    include_axis_profiles = TRUE,
    include_cloud_profiles = TRUE,
    design_method = "hybrid",
    distance_metric = design_control$distance_metric
  )

  selected_idx_hyper <- match(selected_hyper_names, model$hyper_names)
  local_theta_designs <- vector("list", n_locals)
  names(local_theta_designs) <- local_scores$local
  for (local_pos in seq_len(n_locals)) {
    atlas <- root_atlases[[local_pos]]
    root <- atlas$charts[[atlas$root_chart_id]]
    g <- as.numeric(root$score)
    names(g) <- model$hyper_names
    H <- as.matrix(root$curvature)
    dimnames(H) <- list(model$hyper_names, model$hyper_names)
    delta <- sweep(candidate_theta[, selected_idx_hyper, drop = FALSE], 2L, theta_center[selected_idx_hyper], "-")
    Hs <- H[selected_idx_hyper, selected_idx_hyper, drop = FALSE]
    score <- abs(as.numeric(delta %*% g[selected_idx_hyper]) +
      0.5 * rowSums((delta %*% Hs) * delta))
    is_root <- rowSums(sweep(candidate_theta, 2L, theta_root[1L, ], "-")^2) <= 1e-20
    score[is_root] <- Inf
    order_idx <- order(score, decreasing = TRUE)
    take <- order_idx[seq_len(min(length(order_idx), local_anchor_budgets$anchor_budget[local_pos]))]
    design <- .local_atlas_unique_theta(candidate_theta[take, , drop = FALSE], model)
    if (!any(rowSums(sweep(design, 2L, theta_root[1L, ], "-")^2) <= 1e-20)) {
      design <- .local_atlas_unique_theta(rbind(theta_root, design), model)
      design <- design[seq_len(min(nrow(design), local_anchor_budgets$anchor_budget[local_pos])), , drop = FALSE]
    }
    local_theta_designs[[local_pos]] <- design
  }

  new_chart_design_plan(
    theta_root = theta_root,
    theta_cloud = theta_cloud,
    theta_weights = theta_weights,
    families = families,
    family_scores = family_scores,
    local_scores = local_scores,
    local_family_scores = local_family_scores,
    selected_families = selected_families,
    selected_hyper_names = selected_hyper_names,
    local_anchor_budgets = local_anchor_budgets,
    candidate_theta = candidate_theta,
    local_theta_designs = local_theta_designs,
    root_atlases = root_atlases,
    diagnostics = list(
      design_method = "adaptive_scout",
      distance_metric = design_control$distance_metric,
      tail_probs = design_control$tail_probs,
      min_families = min_families,
      max_families = max_families,
      relative_threshold = as.numeric(design_control$scout_relative_threshold),
      min_anchor_budget = min_budget,
      mean_anchor_budget = mean_budget,
      max_anchor_budget = max_budget,
      candidate_pool_size = nrow(candidate_theta)
    ),
    population_model = model
  )
}

.local_atlas_nearest_design_radius <- function(theta_cloud,
                                               theta_design,
                                               population_model,
                                               coverage_prob = 0.995,
                                               inflation = 1.25,
                                               distance_metric = "euclidean") {
  model <- normalize_population_model(population_model)
  theta_cloud <- .as_hyper_matrix(theta_cloud, model$hyper_names, model$hyper_dim)
  theta_design <- .as_hyper_matrix(theta_design, model$hyper_names, model$hyper_dim)
  dists <- matrixStats::rowMins(.local_atlas_theta_distances(
    theta = theta_cloud,
    centers = theta_design,
    population_model = model,
    distance_metric = distance_metric
  ))
  radius <- as.numeric(stats::quantile(dists, probs = as.numeric(coverage_prob), names = FALSE, type = 8))
  if (!is.finite(radius) || radius <= 0) {
    radius <- max(dists[is.finite(dists)], 1)
  }
  as.numeric(radius * as.numeric(inflation))
}

.local_atlas_unique_id <- function(existing, base) {
  base <- gsub("[^A-Za-z0-9_]+", "_", as.character(base))
  if (!base %in% existing) {
    return(base)
  }
  k <- 2L
  repeat {
    candidate <- paste0(base, "_", k)
    if (!candidate %in% existing) {
      return(candidate)
    }
    k <- k + 1L
  }
}

.local_atlas_chart_distance <- function(chart,
                                        theta,
                                        population_model,
                                        distance_metric = "euclidean") {
  model <- normalize_population_model(population_model)
  theta <- .as_hyper_matrix(theta, model$hyper_names, model$hyper_dim)
  as.numeric(.local_atlas_theta_distances(
    theta = theta,
    centers = chart$theta_anchor,
    population_model = model,
    distance_metric = distance_metric
  ))
}

.local_atlas_nearest_active_ids <- function(atlas,
                                           theta,
                                           population_model,
                                           exclude = character(),
                                           distance_metric = "euclidean") {
  active <- .local_atlas_active_charts(atlas)
  active <- active[setdiff(names(active), exclude)]
  if (!length(active)) {
    return(character())
  }
  d <- vapply(
    active,
    .local_atlas_chart_distance,
    numeric(1),
    theta = theta,
    population_model = population_model,
    distance_metric = distance_metric
  )
  names(sort(d))
}

.local_atlas_quarantine_chart <- function(atlas, chart_id, reason) {
  chart <- atlas$charts[[chart_id]]
  if (!is.null(chart) && !identical(chart$status, "active")) {
    chart$status <- "quarantined"
    chart$diagnostics$status_reason <- as.character(reason)
    atlas$charts[[chart_id]] <- chart
  }
  validate_local_atlas(atlas)
}

.local_atlas_candidate_chart <- function(atlas,
                                         chart_id,
                                         theta_anchor,
                                         data_i,
                                         loglik_fn,
                                         population_model,
                                         local_id,
                                         local_control,
                                         seed,
                                         verbose) {
  chart_id <- .local_atlas_unique_id(names(atlas$charts), chart_id)
  chart <- build_candidate_chart(
    local_id = local_id,
    chart_id = chart_id,
    theta_anchor = theta_anchor,
    data_i = data_i,
    loglik_fn = loglik_fn,
    population_model = population_model,
    M = as.integer(local_control$candidate_M),
    target_cess = local_control$target_cess,
    resample_threshold = local_control$resample_threshold,
    n_mcmc_moves = as.integer(local_control$n_mcmc_moves),
    rw_scale = local_control$rw_scale,
    G_mix = as.integer(local_control$G_mix),
    da_enable = isTRUE(local_control$da_enable),
    refit_every = as.integer(local_control$refit_every),
    max_steps = as.integer(local_control$max_steps),
    deterministic_resampling = isTRUE(local_control$deterministic_resampling),
    n_cores = as.integer(local_control$n_cores_per_local),
    seed = as.integer(seed),
    verbose = isTRUE(verbose)
  )
  atlas$charts[[chart$chart_id]] <- chart
  validate_local_atlas(atlas)
}

.local_atlas_estimate_certified_edge <- function(atlas,
                                                from_id,
                                                to_id,
                                                population_model,
                                                edge_control) {
  edge_id <- .local_atlas_unique_id(
    names(atlas$edges),
    paste0("edge_", from_id, "_", to_id)
  )
  edge <- tryCatch(
    estimate_chart_edge(
      from_chart = atlas$charts[[from_id]],
      to_chart = atlas$charts[[to_id]],
      population_model = population_model,
      edge_id = edge_id,
      method = edge_control$method,
      max_iter = as.integer(edge_control$bar_max_iter),
      tol = edge_control$bar_tol
    ),
    error = function(e) {
      new_local_edge(
        local_id = atlas$local_id,
        edge_id = edge_id,
        from_chart = from_id,
        to_chart = to_id,
        method = "bridge",
        status = "rejected",
        diagnostics = list(certification = list(passed = FALSE, failures = "edge_estimation_error"),
                           error = conditionMessage(e))
      )
    }
  )
  if (!identical(edge$status, "rejected")) {
    edge <- certify_chart_edge(
      edge,
      max_forward_reverse_gap = edge_control$max_forward_reverse_gap,
      max_se = edge_control$max_se,
      min_overlap_ess = edge_control$min_overlap_ess,
      max_psis_k = edge_control$max_psis_k,
      max_taylor_gap = edge_control$max_taylor_gap,
      max_taylor_disagreement = edge_control$max_taylor_disagreement,
      require_bar_converged = isTRUE(edge_control$require_bar_converged)
    )
  }
  edge
}

.local_atlas_try_activate_chart <- function(atlas,
                                           chart_id,
                                           population_model,
                                           edge_control,
                                           from_ids = NULL) {
  atlas <- validate_local_atlas(atlas)
  chart <- atlas$charts[[chart_id]]
  if (is.null(chart)) {
    stop("chart_id is not present in atlas.")
  }
  if (identical(chart$status, "active")) {
    return(list(atlas = atlas, success = TRUE, reason = "already_active", from_id = chart_id))
  }
  if (identical(chart$status, "quarantined")) {
    return(list(atlas = atlas, success = FALSE, reason = "chart_quarantined", from_id = NA_character_))
  }
  if (is.null(from_ids)) {
    from_ids <- .local_atlas_nearest_active_ids(
      atlas,
      chart$theta_anchor,
      population_model,
      exclude = chart_id,
      distance_metric = edge_control$distance_metric %||% "euclidean"
    )
  } else {
    from_ids <- intersect(as.character(from_ids), names(.local_atlas_active_charts(atlas)))
  }
  from_ids <- head(from_ids, as.integer(edge_control$edge_neighbors))
  if (!length(from_ids)) {
    return(list(atlas = atlas, success = FALSE, reason = "no_active_neighbor", from_id = NA_character_))
  }

  edges <- lapply(from_ids, function(from_id) {
    .local_atlas_estimate_certified_edge(
      atlas = atlas,
      from_id = from_id,
      to_id = chart_id,
      population_model = population_model,
      edge_control = edge_control
    )
  })
  active_edges <- edges[vapply(edges, function(edge) identical(edge$status, "active"), logical(1))]
  rejected_edges <- edges[!vapply(edges, function(edge) identical(edge$status, "active"), logical(1))]
  for (edge in rejected_edges) {
    atlas$edges[[edge$edge_id]] <- edge
  }
  if (!length(active_edges)) {
    return(list(atlas = validate_local_atlas(atlas), success = FALSE, reason = "no_certified_edge", from_id = from_ids[1L]))
  }

  active_edges <- active_edges[order(vapply(active_edges, function(edge) edge$se, numeric(1)))]
  chart$status <- "active"
  chart$diagnostics$status_reason <- "certified_by_relative_edge"
  atlas$charts[[chart_id]] <- validate_local_chart(chart, population_model)

  best <- active_edges[[1L]]
  atlas$edges[[best$edge_id]] <- best
  solved <- tryCatch(
    solve_atlas_normalizers(atlas, require_connected = TRUE),
    error = function(e) e
  )
  if (inherits(solved, "error")) {
    best$status <- "rejected"
    best$diagnostics$certification$failures <- c(
      best$diagnostics$certification$failures %||% character(),
      "graph_solve_failed"
    )
    best$diagnostics$graph_error <- conditionMessage(solved)
    atlas$edges[[best$edge_id]] <- best
    atlas <- .local_atlas_quarantine_chart(atlas, chart_id, "graph_solve_failed")
    return(list(atlas = atlas, success = FALSE, reason = "graph_solve_failed", from_id = from_ids[1L]))
  }
  atlas <- solved
  chart <- atlas$charts[[chart_id]]
  if (.local_chart_normalizer_certified(chart)) {
    chart$diagnostics$edge_normalizer_certification <- list(
      certified = TRUE,
      method = "relative_edge",
      reason = "certified_by_relative_edge",
      edge_id = best$edge_id,
      from_chart = best$from_chart
    )
  } else {
    chart <- .local_chart_set_normalizer_certification(
      chart,
      certified = TRUE,
      method = "relative_edge",
      reason = "certified_by_relative_edge",
      details = list(edge_id = best$edge_id, from_chart = best$from_chart)
    )
  }
  atlas$charts[[chart_id]] <- validate_local_chart(chart, population_model)

  if (length(active_edges) > 1L) {
    for (edge in active_edges[-1L]) {
      candidate <- atlas
      candidate$edges[[edge$edge_id]] <- edge
      candidate <- tryCatch(
        solve_atlas_normalizers(candidate, require_connected = TRUE),
        error = function(e) e
      )
      if (inherits(candidate, "error")) {
        edge$status <- "rejected"
        edge$diagnostics$certification$failures <- c(
          edge$diagnostics$certification$failures %||% character(),
          "graph_solve_failed"
        )
        edge$diagnostics$graph_error <- conditionMessage(candidate)
        atlas$edges[[edge$edge_id]] <- edge
        next
      }
      residual <- as.numeric(candidate$diagnostics$normalizer_solution$max_abs_standardized_edge_residual %||% 0)
      if (is.finite(residual) && residual <= as.numeric(edge_control$max_graph_standardized_residual)) {
        atlas <- candidate
      } else {
        edge$status <- "rejected"
        edge$diagnostics$certification$failures <- c(
          edge$diagnostics$certification$failures %||% character(),
          "graph_residual"
        )
        edge$diagnostics$graph_standardized_residual <- residual
        atlas$edges[[edge$edge_id]] <- edge
      }
    }
  }

  list(atlas = validate_local_atlas(atlas), success = TRUE, reason = "certified", from_id = best$from_chart)
}

.local_atlas_activate_existing_chart <- function(atlas,
                                                 target_id,
                                                 data_i,
                                                 loglik_fn,
                                                 population_model,
                                                 local_id,
                                                 local_control,
                                                 edge_control,
                                                 seed,
                                                 verbose = FALSE) {
  model <- normalize_population_model(population_model)
  atlas <- validate_local_atlas(atlas)
  if (is.null(atlas$charts[[target_id]])) {
    stop("target_id is not present in the atlas.")
  }

	  attempt <- .local_atlas_try_activate_chart(atlas, target_id, model, edge_control)
	  atlas <- attempt$atlas
	  last_reason <- attempt$reason
	  if (!isTRUE(attempt$success)) {
	    bridge_from_ids <- .local_atlas_nearest_active_ids(
	      atlas,
	      atlas$charts[[target_id]]$theta_anchor,
	      model,
	      exclude = target_id,
	      distance_metric = edge_control$distance_metric %||% "euclidean"
	    )
	    bridge_from_ids <- unique(c(attempt$from_id, bridge_from_ids))
	    bridge_from_ids <- bridge_from_ids[!is.na(bridge_from_ids) & nzchar(bridge_from_ids)]
	    bridge_from_ids <- head(bridge_from_ids, as.integer(edge_control$edge_neighbors))
	    for (start_id in bridge_from_ids) {
	      if (identical(atlas$charts[[target_id]]$status, "active")) {
	        break
	      }
	      if (is.null(atlas$charts[[start_id]]) ||
	          !identical(atlas$charts[[start_id]]$status, "active")) {
	        next
	      }
	      current_id <- start_id
	      for (depth in seq_len(as.integer(edge_control$max_intermediates))) {
	        if (identical(atlas$charts[[target_id]]$status, "active")) {
	          break
	        }
	        current_theta <- atlas$charts[[current_id]]$theta_anchor
	        target_theta <- atlas$charts[[target_id]]$theta_anchor
	        mid_theta <- matrix(
	          (as.numeric(current_theta[1L, ]) + as.numeric(target_theta[1L, ])) / 2,
	          nrow = 1L,
	          dimnames = list(NULL, model$hyper_names)
	        )
	        mid_id <- sprintf("%s_mid_%s_%02d", target_id, start_id, depth)
	        atlas <- .local_atlas_candidate_chart(
	          atlas = atlas,
	          chart_id = mid_id,
	          theta_anchor = mid_theta,
	          data_i = data_i,
	          loglik_fn = loglik_fn,
	          population_model = model,
	          local_id = local_id,
	          local_control = local_control,
	          seed = as.integer(seed) + 7919L * depth + 104729L * match(start_id, bridge_from_ids),
	          verbose = verbose
	        )
	        mid_id <- tail(names(atlas$charts), 1L)
	        mid_attempt <- .local_atlas_try_activate_chart(
	          atlas,
	          mid_id,
	          model,
	          edge_control,
	          from_ids = current_id
	        )
	        atlas <- mid_attempt$atlas
	        if (!isTRUE(mid_attempt$success)) {
	          atlas <- .local_atlas_quarantine_chart(atlas, mid_id, mid_attempt$reason)
	          last_reason <- mid_attempt$reason
	          break
	        }
	        current_id <- mid_id
	        target_attempt <- .local_atlas_try_activate_chart(
	          atlas,
	          target_id,
	          model,
	          edge_control,
	          from_ids = current_id
	        )
	        atlas <- target_attempt$atlas
	        last_reason <- target_attempt$reason
	      }
	    }
	  }

  success <- identical(atlas$charts[[target_id]]$status, "active")
  if (!isTRUE(success)) {
    atlas <- .local_atlas_quarantine_chart(atlas, target_id, "no_certified_path_to_root")
    last_reason <- last_reason %||% "no_certified_path_to_root"
  }
  atlas <- solve_atlas_normalizers(atlas, require_connected = TRUE)
  list(
    atlas = validate_local_atlas(atlas),
    success = isTRUE(success),
    reason = if (isTRUE(success)) "certified" else as.character(last_reason),
    target_id = target_id
  )
}

.local_atlas_insert_smc_chart <- function(atlas,
                                          run,
                                          chart_id,
                                          theta_anchor,
                                          population_model,
                                          status_reason,
                                          source_run_id = chart_id,
                                          direct_observation_role = NULL,
                                          normalizer_certified = FALSE,
                                          normalizer_certification_method = NULL,
                                          normalizer_certification_reason = NULL) {
  model <- normalize_population_model(population_model)
  atlas <- validate_local_atlas(atlas)
  chart_id <- .local_atlas_unique_id(names(atlas$charts), chart_id)
  chart <- .local_atlas_chart_from_smc_run(
    run = run,
    local_id = atlas$local_id,
    chart_id = chart_id,
    theta_anchor = theta_anchor,
    population_model = model,
    status = "candidate",
    status_reason = status_reason,
    source_run_id = source_run_id,
    direct_observation_role = direct_observation_role
  )
  chart <- .local_chart_set_normalizer_certification(
    chart,
    certified = isTRUE(normalizer_certified),
    method = normalizer_certification_method %||% if (isTRUE(normalizer_certified)) "direct_smc" else "support_probe",
    reason = normalizer_certification_reason %||% if (isTRUE(normalizer_certified)) {
      "direct_smc_normalizer_accepted"
    } else {
      "normalizer_not_directly_certified"
    }
  )
  atlas$charts[[chart_id]] <- chart
  list(atlas = validate_local_atlas(atlas), chart_id = chart_id)
}

.local_atlas_certify_direct_normalizer <- function(atlas,
                                                   chart_id,
                                                   population_model,
                                                   max_direct_z = 3,
                                                   max_chart_shift = 0.35,
                                                   max_existing_shift = 0.15) {
  atlas <- validate_local_atlas(atlas)
  chart <- atlas$charts[[chart_id]]
  if (is.null(chart)) {
    stop("chart_id is not present in atlas.")
  }
  role <- as.character(chart$diagnostics$direct_observation_role %||% "")
  if (!identical(role, "calibration")) {
    return(list(
      atlas = atlas,
      certified = FALSE,
      reason = "not_a_direct_calibration_chart",
      diagnostics = list()
    ))
  }
  full_solution <- atlas$normalizer_solution
  if (is.null(full_solution)) {
    return(list(
      atlas = atlas,
      certified = FALSE,
      reason = "missing_full_graph_solution",
      diagnostics = list()
    ))
  }
  direct <- full_solution$direct_residuals
  direct_row <- if (nrow(direct) && chart_id %in% direct$chart_id) {
    direct[match(chart_id, direct$chart_id), , drop = FALSE]
  } else {
    data.frame()
  }
  direct_z <- if (nrow(direct_row)) {
    abs(as.numeric(direct_row$standardized_residual[1L]))
  } else {
    Inf
  }

  leave_atlas <- atlas
  leave_chart <- leave_atlas$charts[[chart_id]]
  leave_chart$diagnostics$direct_observation_role <- "support"
  leave_chart <- .local_chart_set_normalizer_certification(
    leave_chart,
    certified = TRUE,
    method = "relative_edge",
    reason = "direct_observation_left_out_for_graph_stability",
    details = leave_chart$diagnostics$edge_normalizer_certification %||% list()
  )
  leave_atlas$charts[[chart_id]] <- validate_local_chart(leave_chart, population_model)
  leave_atlas <- tryCatch(
    solve_atlas_normalizers(leave_atlas, require_connected = TRUE),
    error = function(e) e
  )
  if (inherits(leave_atlas, "error")) {
    chart <- .local_chart_set_normalizer_certification(
      chart,
      certified = FALSE,
      method = "replicated_smc_graph_stability",
      reason = "leave_direct_out_graph_failed",
      details = list(error = conditionMessage(leave_atlas))
    )
    atlas$charts[[chart_id]] <- validate_local_chart(chart, population_model)
    return(list(
      atlas = validate_local_atlas(atlas),
      certified = FALSE,
      reason = "leave_direct_out_graph_failed",
      diagnostics = list(error = conditionMessage(leave_atlas), direct_z = direct_z)
    ))
  }

  full_logZ <- full_solution$logZ
  leave_logZ <- leave_atlas$normalizer_solution$logZ
  common <- intersect(names(full_logZ)[is.finite(full_logZ)], names(leave_logZ)[is.finite(leave_logZ)])
  shifts <- abs(as.numeric(full_logZ[common]) - as.numeric(leave_logZ[common]))
  names(shifts) <- common
  chart_shift <- if (chart_id %in% names(shifts)) as.numeric(shifts[[chart_id]]) else Inf
  existing <- setdiff(names(shifts), chart_id)
  existing_shift <- if (length(existing)) max(shifts[existing], na.rm = TRUE) else 0
  if (!is.finite(existing_shift)) {
    existing_shift <- Inf
  }
  diagnostics <- list(
    direct_z = as.numeric(direct_z),
    chart_shift = as.numeric(chart_shift),
    max_existing_shift = as.numeric(existing_shift),
    max_direct_z = as.numeric(max_direct_z),
    max_chart_shift = as.numeric(max_chart_shift),
    max_allowed_existing_shift = as.numeric(max_existing_shift)
  )
  passed <- is.finite(direct_z) &&
    direct_z <= as.numeric(max_direct_z) &&
    is.finite(chart_shift) &&
    chart_shift <= as.numeric(max_chart_shift) &&
    is.finite(existing_shift) &&
    existing_shift <= as.numeric(max_existing_shift)

  if (isTRUE(passed)) {
    chart <- .local_chart_set_normalizer_certification(
      atlas$charts[[chart_id]],
      certified = TRUE,
      method = "replicated_smc_graph_stable",
      reason = "replicated_direct_observation_passed_leave_out_graph_stability",
      details = diagnostics
    )
    atlas$charts[[chart_id]] <- validate_local_chart(chart, population_model)
    return(list(
      atlas = validate_local_atlas(atlas),
      certified = TRUE,
      reason = "replicated_direct_graph_stable",
      diagnostics = diagnostics
    ))
  }

  chart <- leave_atlas$charts[[chart_id]]
  chart <- .local_chart_set_normalizer_certification(
    chart,
    certified = TRUE,
    method = "relative_edge",
    reason = "direct_observation_rejected_by_graph_stability",
    details = diagnostics
  )
  leave_atlas$charts[[chart_id]] <- validate_local_chart(chart, population_model)
  list(
    atlas = validate_local_atlas(leave_atlas),
    certified = FALSE,
    reason = "direct_observation_rejected_by_graph_stability",
    diagnostics = diagnostics
  )
}

.local_atlas_normalizer_shift <- function(before_atlas,
                                          after_atlas,
                                          added_chart_id = NULL) {
  before_solution <- before_atlas$normalizer_solution
  after_solution <- after_atlas$normalizer_solution
  if (is.null(before_solution) || is.null(after_solution)) {
    return(list(
      max_existing_shift = Inf,
      mean_existing_shift = Inf,
      n_common = 0L,
      max_edge_z = Inf,
      max_direct_z = Inf
    ))
  }
  before_logZ <- before_solution$logZ
  after_logZ <- after_solution$logZ
  common <- intersect(names(before_logZ)[is.finite(before_logZ)], names(after_logZ)[is.finite(after_logZ)])
  common <- setdiff(common, as.character(added_chart_id %||% character()))
  shifts <- if (length(common)) abs(as.numeric(after_logZ[common]) - as.numeric(before_logZ[common])) else numeric()
  diag <- after_atlas$diagnostics$normalizer_solution %||% list()
  list(
    max_existing_shift = if (length(shifts)) max(shifts, na.rm = TRUE) else 0,
    mean_existing_shift = if (length(shifts)) mean(shifts, na.rm = TRUE) else 0,
    n_common = as.integer(length(shifts)),
    max_edge_z = as.numeric(diag$max_abs_standardized_edge_residual %||% 0),
    max_direct_z = as.numeric(diag$max_abs_standardized_direct_residual %||% 0)
  )
}

.local_atlas_expand_design <- function(atlas,
                                       local_id,
                                       data_i,
                                       loglik_fn,
                                       population_model,
                                       theta_design,
                                       theta_root,
                                       local_control,
                                       edge_control,
                                       seed,
                                       verbose = FALSE) {
  model <- normalize_population_model(population_model)
  atlas <- validate_local_atlas(atlas)
  theta_design <- .local_atlas_unique_theta(rbind(theta_root, theta_design), model)
  theta_root <- .local_chart_align_theta_one(theta_root, model)
  if (!identical(as.character(atlas$local_id), as.character(local_id))) {
    stop("root atlas local_id does not match requested local_id.")
  }

  target_rows <- seq_len(nrow(theta_design))
  is_root <- apply(theta_design, 1L, function(x) isTRUE(all.equal(as.numeric(x), as.numeric(theta_root[1L, ]))))
  target_rows <- target_rows[!is_root]
  history <- list()
  for (pos in seq_along(target_rows)) {
    design_row <- target_rows[pos]
    target_id <- sprintf("anchor_%03d", design_row)
    atlas <- .local_atlas_candidate_chart(
      atlas = atlas,
      chart_id = target_id,
      theta_anchor = theta_design[design_row, , drop = FALSE],
      data_i = data_i,
      loglik_fn = loglik_fn,
      population_model = model,
      local_id = local_id,
      local_control = local_control,
      seed = as.integer(seed) + 1009L * design_row,
      verbose = verbose
    )
    target_id <- tail(names(atlas$charts), 1L)
    attempt <- .local_atlas_try_activate_chart(atlas, target_id, model, edge_control)
    atlas <- attempt$atlas
    if (!isTRUE(attempt$success)) {
      current_id <- attempt$from_id
      if (!is.na(current_id) && nzchar(current_id)) {
        for (depth in seq_len(as.integer(edge_control$max_intermediates))) {
          if (identical(atlas$charts[[target_id]]$status, "active")) {
            break
          }
          current_theta <- atlas$charts[[current_id]]$theta_anchor
          target_theta <- atlas$charts[[target_id]]$theta_anchor
          mid_theta <- matrix((as.numeric(current_theta[1L, ]) + as.numeric(target_theta[1L, ])) / 2,
                              nrow = 1L, dimnames = list(NULL, model$hyper_names))
          mid_id <- sprintf("%s_mid_%02d", target_id, depth)
          atlas <- .local_atlas_candidate_chart(
            atlas = atlas,
            chart_id = mid_id,
            theta_anchor = mid_theta,
            data_i = data_i,
            loglik_fn = loglik_fn,
            population_model = model,
            local_id = local_id,
            local_control = local_control,
            seed = as.integer(seed) + 1009L * design_row + 7919L * depth,
            verbose = verbose
          )
          mid_id <- tail(names(atlas$charts), 1L)
          mid_attempt <- .local_atlas_try_activate_chart(
            atlas,
            mid_id,
            model,
            edge_control,
            from_ids = current_id
          )
          atlas <- mid_attempt$atlas
          if (!isTRUE(mid_attempt$success)) {
            atlas <- .local_atlas_quarantine_chart(atlas, mid_id, mid_attempt$reason)
            break
          }
          current_id <- mid_id
          target_attempt <- .local_atlas_try_activate_chart(
            atlas,
            target_id,
            model,
            edge_control,
            from_ids = current_id
          )
          atlas <- target_attempt$atlas
        }
      }
    }
    if (!identical(atlas$charts[[target_id]]$status, "active")) {
      atlas <- .local_atlas_quarantine_chart(atlas, target_id, "no_certified_path_to_root")
    }
    atlas <- solve_atlas_normalizers(atlas, require_connected = TRUE)
    history[[length(history) + 1L]] <- data.frame(
      local = as.character(local_id),
      target_chart = target_id,
      target_row = as.integer(design_row),
      status = atlas$charts[[target_id]]$status,
      n_active_charts = sum(vapply(atlas$charts, function(chart) identical(chart$status, "active"), logical(1))),
      n_active_edges = sum(vapply(atlas$edges, function(edge) identical(edge$status, "active"), logical(1))),
      check.names = FALSE
    )
  }

  atlas$diagnostics$build_history <- if (length(history)) do.call(rbind, history) else data.frame()
  validate_local_atlas(atlas)
}

build_local_atlas <- function(local_id,
                              data_i,
                              loglik_fn,
                              population_model,
                              theta_design,
                              theta_root = theta_design[1L, , drop = FALSE],
                              local_control = list(),
                              edge_control = list(),
                              n_cores = 1L,
                              seed = 123L,
                              verbose = FALSE) {
  model <- normalize_population_model(population_model)
  theta_design <- .local_atlas_unique_theta(theta_design, model)
  theta_root <- .local_chart_align_theta_one(theta_root, model)
  local_control <- .local_atlas_merge_control(local_control, .local_atlas_default_local_control())
  edge_control <- .local_atlas_merge_control(edge_control, .local_atlas_default_edge_control())
  local_control$n_cores_per_local <- min(as.integer(local_control$n_cores_per_local), as.integer(n_cores))

  root <- build_local_root_chart(
    local_id = local_id,
    theta_anchor = theta_root,
    data_i = data_i,
    loglik_fn = loglik_fn,
    population_model = model,
    chart_id = "root",
    M = as.integer(local_control$root_M),
    target_cess = local_control$target_cess,
    resample_threshold = local_control$resample_threshold,
    n_mcmc_moves = as.integer(local_control$n_mcmc_moves),
    rw_scale = local_control$rw_scale,
    G_mix = as.integer(local_control$G_mix),
    da_enable = isTRUE(local_control$da_enable),
    refit_every = as.integer(local_control$refit_every),
    max_steps = as.integer(local_control$max_steps),
    deterministic_resampling = isTRUE(local_control$deterministic_resampling),
    n_cores = as.integer(local_control$n_cores_per_local),
    seed = as.integer(seed),
    verbose = isTRUE(verbose),
    confirm = local_control$root_confirm,
    confirmation_M = as.integer(local_control$root_confirmation_M %||% max(100L, ceiling(as.integer(local_control$root_M) / 2L))),
    confirmation_abs_tol = local_control$root_confirmation_abs_tol
  )
  atlas <- solve_atlas_normalizers(new_local_atlas(
    local_id = local_id,
    root_chart_id = "root",
    charts = list(root),
    edges = list()
  ))

  .local_atlas_expand_design(
    atlas = atlas,
    local_id = local_id,
    data_i = data_i,
    loglik_fn = loglik_fn,
    population_model = model,
    theta_design = theta_design,
    theta_root = theta_root,
    local_control = local_control,
    edge_control = edge_control,
    seed = as.integer(seed),
    verbose = isTRUE(verbose)
  )
}

build_local_atlas_from_root <- function(root_atlas,
                                        data_i,
                                        loglik_fn,
                                        population_model,
                                        theta_design,
                                        theta_root = NULL,
                                        local_control = list(),
                                        edge_control = list(),
                                        n_cores = 1L,
                                        seed = 123L,
                                        verbose = FALSE) {
  model <- normalize_population_model(population_model)
  root_atlas <- validate_local_atlas(root_atlas)
  theta_root <- theta_root %||% root_atlas$charts[[root_atlas$root_chart_id]]$theta_anchor
  theta_root <- .local_chart_align_theta_one(theta_root, model)
  local_control <- .local_atlas_merge_control(local_control, .local_atlas_default_local_control())
  edge_control <- .local_atlas_merge_control(edge_control, .local_atlas_default_edge_control())
  local_control$n_cores_per_local <- min(as.integer(local_control$n_cores_per_local), as.integer(n_cores))

  .local_atlas_expand_design(
    atlas = root_atlas,
    local_id = root_atlas$local_id,
    data_i = data_i,
    loglik_fn = loglik_fn,
    population_model = model,
    theta_design = theta_design,
    theta_root = theta_root,
    local_control = local_control,
    edge_control = edge_control,
    seed = as.integer(seed),
    verbose = isTRUE(verbose)
  )
}

local_atlas_pre_outer_certify <- function(factor_set,
                                          data_list,
                                          loglik_fn,
                                          initial_proposal,
                                          local_control = list(),
                                          edge_control = list(),
                                          calibration_control = list(),
                                          outer_control = list(),
                                          design_control = list(),
                                          calibration_history = list(),
                                          n_cores = 1L,
                                          seed = 123L,
                                          outer_seed = as.integer(seed) + 900001L,
                                          verbose = TRUE,
                                          trace_verbose = FALSE,
                                          checkpoint_callback = NULL) {
  factor_set <- validate_local_atlas_factor_set(factor_set)
  if (!is.list(data_list) || length(data_list) < length(factor_set$atlases)) {
    stop("data_list must contain one entry per atlas.")
  }
  if (!is.function(loglik_fn)) {
    stop("loglik_fn must be a function.")
  }
  if (is.null(initial_proposal)) {
    stop("initial_proposal is required for pre-outer certification.")
  }
  local_control <- .local_atlas_merge_control(local_control, .local_atlas_default_local_control())
  edge_control <- .local_atlas_merge_control(edge_control, .local_atlas_default_edge_control())
  calibration_control <- .local_atlas_merge_control(calibration_control, list(
    M = as.integer(local_control$candidate_M %||% 500L),
    target_cess = local_control$target_cess %||% 0.9,
    n_mcmc_moves = as.integer(local_control$n_mcmc_moves %||% 2L),
    max_steps = as.integer(local_control$max_steps %||% 128L),
    max_updates = 12L,
    confirmation_reps = 0L,
    confirmation_M = as.integer(local_control$candidate_M %||% 500L),
    confirmation_max_sd = 1.5,
    replicate_bootstrap_B = 200L,
    max_direct_graph_z = 3,
    max_direct_graph_chart_shift = 0.35,
    max_direct_graph_existing_shift = 0.15,
    pre_outer_rounds = 1L,
    pre_outer_audit_n = NULL,
    pre_outer_max_points = 6L,
    pre_outer_max_updates = 16L
  ))
  outer_control <- .local_atlas_merge_control(outer_control, list(N = 1000L))
  if (is.data.frame(calibration_history)) {
    calibration_history <- list(calibration_history)
  }
  if (!is.list(calibration_history)) {
    calibration_history <- list()
  }

  model <- factor_set$population_model
  max_pre_outer_rounds <- as.integer(calibration_control$pre_outer_rounds %||% 0L)
  for (pre_round in seq_len(max_pre_outer_rounds)) {
    pre_start <- Sys.time()
    pre_audit_n <- as.integer(calibration_control$pre_outer_audit_n %||%
      max(as.integer(outer_control$N), as.integer(design_control$refine_audit_n %||% outer_control$N)))
    .local_atlas_log("pre-outer audit round ", pre_round, ": theta=", pre_audit_n, "\n", verbose = verbose)
    audit_theta <- theta_proposal_sample(
      initial_proposal,
      n = pre_audit_n,
      seed = as.integer(seed) + 39000017L * pre_round
    )
    outer_init_theta <- theta_proposal_sample(
      initial_proposal,
      n = as.integer(outer_control$N),
      seed = as.integer(outer_seed)
    )
    audit_theta <- .local_atlas_unique_theta(rbind(outer_init_theta, audit_theta), model)
    protected_theta_rows <- seq_len(min(nrow(outer_init_theta), nrow(audit_theta)))
    audit_parts <- .local_atlas_factor_set_by_local(factor_set, audit_theta, n_cores = as.integer(n_cores))
    bad <- .local_atlas_factor_set_uncertified(audit_parts)
    if (!nrow(bad)) {
      .local_atlas_log("pre-outer audit round ", pre_round, ": no uncertified pairs\n", verbose = verbose)
      break
    }

    theta_counts <- sort(table(bad$theta_row), decreasing = TRUE)
    local_counts <- table(bad$local)
    bad$local_pos <- match(as.character(bad$local), names(factor_set$atlases))
    bad$theta_bad_count <- as.numeric(theta_counts[match(as.character(bad$theta_row), names(theta_counts))])
    bad$local_bad_count <- as.numeric(local_counts[match(as.character(bad$local), names(local_counts))])
    bad$priority <- bad$theta_bad_count + 0.25 * bad$local_bad_count
    bad <- bad[order(bad$priority, decreasing = TRUE), , drop = FALSE]

    max_pre_theta <- as.integer(calibration_control$pre_outer_max_points)
    max_pre_updates <- as.integer(calibration_control$pre_outer_max_updates)
    counted_rows <- as.integer(names(theta_counts))
    protected_bad_rows <- intersect(protected_theta_rows, counted_rows)
    max_pre_theta <- max(max_pre_theta, length(protected_bad_rows))
    keep_rows <- unique(c(
      protected_bad_rows,
      setdiff(counted_rows, protected_bad_rows)
    ))
    keep_rows <- keep_rows[seq_len(min(length(keep_rows), max_pre_theta))]
    eligible_bad <- bad[bad$theta_row %in% keep_rows, , drop = FALSE]
    if (!nrow(eligible_bad)) {
      eligible_bad <- bad
    }
    eligible_bad$protected_outer_initial <- eligible_bad$theta_row %in% protected_bad_rows
    selected_idx <- integer()
    selected_key <- character()
    add_selected <- function(idx) {
      key <- paste(eligible_bad$local_pos[idx], eligible_bad$theta_row[idx], sep = ":")
      take <- !key %in% selected_key
      if (any(take)) {
        selected_idx <<- c(selected_idx, idx[take])
        selected_key <<- c(selected_key, key[take])
      }
    }
    for (theta_row in protected_bad_rows) {
      theta_idx <- which(eligible_bad$theta_row == theta_row)
      if (!length(theta_idx)) next
      theta_idx <- theta_idx[order(
        eligible_bad$local_bad_count[theta_idx],
        eligible_bad$priority[theta_idx],
        decreasing = TRUE
      )]
      for (idx_one in theta_idx) {
        add_selected(idx_one)
        if (length(selected_idx) >= max_pre_updates) break
      }
      if (length(selected_idx) >= max_pre_updates) break
    }
    for (loc in names(sort(local_counts, decreasing = TRUE))) {
      if (length(selected_idx) >= max_pre_updates) break
      loc_idx <- which(eligible_bad$local == loc)
      if (!length(loc_idx)) next
      loc_idx <- loc_idx[order(
        eligible_bad$theta_bad_count[loc_idx],
        eligible_bad$priority[loc_idx],
        decreasing = TRUE
      )]
      add_selected(loc_idx[1L])
      if (length(selected_idx) >= max_pre_updates) break
    }
    if (length(selected_idx) < max_pre_updates) {
      for (theta_row in keep_rows) {
        theta_idx <- which(eligible_bad$theta_row == theta_row)
        if (!length(theta_idx)) next
        theta_idx <- theta_idx[order(
          eligible_bad$local_bad_count[theta_idx],
          eligible_bad$priority[theta_idx],
          decreasing = TRUE
        )]
        for (idx_one in theta_idx) {
          add_selected(idx_one)
          if (length(selected_idx) >= max_pre_updates) break
        }
        if (length(selected_idx) >= max_pre_updates) break
      }
    }
    if (length(selected_idx) < max_pre_updates) {
      for (idx_one in seq_len(nrow(eligible_bad))) {
        add_selected(idx_one)
        if (length(selected_idx) >= max_pre_updates) break
      }
    }
    selected_bad <- eligible_bad[selected_idx[seq_len(min(length(selected_idx), max_pre_updates))], , drop = FALSE]
    .local_atlas_log(
      "pre-outer audit round ", pre_round,
      ": uncertified_pairs=", nrow(bad),
      " selected_pairs=", nrow(selected_bad),
      " repair_cap=", as.integer(calibration_control$pre_outer_max_updates),
      " selected_theta=", length(unique(selected_bad$theta_row)),
      " selected_locals=", length(unique(selected_bad$local_pos)), "\n",
      verbose = verbose
    )

    repair_weights <- numeric(nrow(audit_theta))
    repair_weights[as.integer(names(theta_counts))] <- as.numeric(theta_counts)
    repair_weights <- .local_chart_normalize_weights(repair_weights, length(repair_weights))
    repair_locals <- unique(selected_bad$local_pos)
    calibration <- local_atlas_repair_certification_pairs(
      factor_set = factor_set,
      theta = audit_theta,
      data_list = data_list,
      loglik_fn = loglik_fn,
      theta_weights = repair_weights,
      candidate_pairs = selected_bad[, c("local_pos", "theta_row"), drop = FALSE],
      local_ids = repair_locals,
      M = as.integer(calibration_control$M),
      target_cess = calibration_control$target_cess,
      n_mcmc_moves = as.integer(calibration_control$n_mcmc_moves),
      max_steps = as.integer(calibration_control$max_steps),
      max_updates = as.integer(calibration_control$pre_outer_max_updates),
      direct_confirmation_reps = as.integer(calibration_control$confirmation_reps),
      direct_confirmation_M = as.integer(calibration_control$confirmation_M),
      direct_confirmation_max_sd = calibration_control$confirmation_max_sd,
      replicate_bootstrap_B = as.integer(calibration_control$replicate_bootstrap_B),
      max_direct_graph_z = calibration_control$max_direct_graph_z,
      max_direct_graph_chart_shift = calibration_control$max_direct_graph_chart_shift,
      max_direct_graph_existing_shift = calibration_control$max_direct_graph_existing_shift,
      local_control = local_control,
      edge_control = edge_control,
      seed = as.integer(seed) + 49000019L * pre_round,
      verbose = isTRUE(trace_verbose)
    )
    calibration_history[[length(calibration_history) + 1L]] <- data.frame(
      calibration_phase = "pre_outer",
      calibration_round = pre_round,
      n_audit_theta = nrow(audit_theta),
      n_uncertified_pairs = nrow(bad),
      calibration$probes,
      check.names = FALSE
    )
    factor_set <- calibration$factor_set
    .local_atlas_log(
      "pre-outer repair round ", pre_round,
      " finished in ", round(as.numeric(difftime(Sys.time(), pre_start, units = "mins")), 3),
      " min | activated=", calibration$n_activated,
      " selected=", calibration$n_selected, "\n",
      verbose = verbose
    )
    if (is.function(checkpoint_callback)) {
      checkpoint_callback(
        stage = paste0("pre_outer_round_", pre_round),
        factor_set = factor_set,
        atlases = factor_set$atlases,
        calibration_history = calibration_history
      )
    }
    if (!isTRUE(calibration$n_activated > 0L)) {
      break
    }
  }

  structure(
    list(
      factor_set = factor_set,
      atlases = factor_set$atlases,
      calibration_history = calibration_history,
      calibration_history_table = .local_atlas_rbind_fill(calibration_history)
    ),
    class = "local_atlas_pre_outer_certification"
  )
}

local_atlas_certify_theta_cloud <- function(factor_set,
                                            theta,
                                            data_list,
                                            loglik_fn,
                                            local_control = list(),
                                            edge_control = list(),
                                            calibration_control = list(),
                                            theta_weights = NULL,
                                            max_rounds = 2L,
                                            max_updates = calibration_control$pre_outer_max_updates %||% calibration_control$max_updates %||% 16L,
                                            n_cores = 1L,
                                            seed = 123L,
                                            verbose = TRUE,
                                            trace_verbose = FALSE) {
  factor_set <- validate_local_atlas_factor_set(factor_set)
  model <- factor_set$population_model
  theta <- .as_hyper_matrix(theta, model$hyper_names, model$hyper_dim)
  theta_weights <- .local_chart_normalize_weights(theta_weights, nrow(theta))
  local_control <- .local_atlas_merge_control(local_control, .local_atlas_default_local_control())
  edge_control <- .local_atlas_merge_control(edge_control, .local_atlas_default_edge_control())
  calibration_control <- .local_atlas_merge_control(calibration_control, list(
    M = as.integer(local_control$candidate_M %||% 500L),
    target_cess = local_control$target_cess %||% 0.9,
    n_mcmc_moves = as.integer(local_control$n_mcmc_moves %||% 2L),
    max_steps = as.integer(local_control$max_steps %||% 128L),
    confirmation_reps = 0L,
    confirmation_M = as.integer(local_control$candidate_M %||% 500L),
    confirmation_max_sd = 1.5,
    replicate_bootstrap_B = 200L,
    max_direct_graph_z = 3,
    max_direct_graph_chart_shift = 0.35,
    max_direct_graph_existing_shift = 0.15
  ))

  history <- list()
  bad <- data.frame()
  for (round_id in seq_len(as.integer(max_rounds))) {
    parts <- .local_atlas_factor_set_by_local(factor_set, theta, n_cores = as.integer(n_cores))
    bad <- .local_atlas_factor_set_uncertified(parts)
    if (!nrow(bad)) {
      .local_atlas_log("theta-cloud certification round ", round_id, ": no uncertified pairs\n", verbose = verbose)
      break
    }
    bad$local_pos <- match(as.character(bad$local), names(factor_set$atlases))
    bad$theta_weight <- theta_weights[bad$theta_row]
    bad <- bad[order(-bad$theta_weight, bad$local_pos, bad$theta_row), , drop = FALSE]
    selected_bad <- bad[seq_len(min(nrow(bad), as.integer(max_updates))), , drop = FALSE]
    .local_atlas_log(
      "theta-cloud certification round ", round_id,
      ": uncertified_pairs=", nrow(bad),
      " selected_pairs=", nrow(selected_bad),
      " repair_cap=", as.integer(max_updates), "\n",
      verbose = verbose
    )
    calibration <- local_atlas_repair_certification_pairs(
      factor_set = factor_set,
      theta = theta,
      data_list = data_list,
      loglik_fn = loglik_fn,
      theta_weights = theta_weights,
      candidate_pairs = selected_bad[, c("local_pos", "theta_row"), drop = FALSE],
      local_ids = unique(selected_bad$local_pos),
      M = as.integer(calibration_control$M),
      target_cess = calibration_control$target_cess,
      n_mcmc_moves = as.integer(calibration_control$n_mcmc_moves),
      max_steps = as.integer(calibration_control$max_steps),
      max_updates = as.integer(max_updates),
      direct_confirmation_reps = as.integer(calibration_control$confirmation_reps),
      direct_confirmation_M = as.integer(calibration_control$confirmation_M),
      direct_confirmation_max_sd = calibration_control$confirmation_max_sd,
      replicate_bootstrap_B = as.integer(calibration_control$replicate_bootstrap_B),
      max_direct_graph_z = calibration_control$max_direct_graph_z,
      max_direct_graph_chart_shift = calibration_control$max_direct_graph_chart_shift,
      max_direct_graph_existing_shift = calibration_control$max_direct_graph_existing_shift,
      local_control = local_control,
      edge_control = edge_control,
      seed = as.integer(seed) + 6100003L * round_id,
      verbose = isTRUE(trace_verbose)
    )
    history[[length(history) + 1L]] <- data.frame(
      calibration_phase = "theta_cloud",
      calibration_round = round_id,
      n_theta = nrow(theta),
      n_uncertified_pairs = nrow(bad),
      calibration$probes,
      check.names = FALSE
    )
    factor_set <- calibration$factor_set
    if (!isTRUE(calibration$n_activated > 0L)) {
      break
    }
  }
  parts <- .local_atlas_factor_set_by_local(factor_set, theta, n_cores = as.integer(n_cores))
  bad <- .local_atlas_factor_set_uncertified(parts)
  structure(
    list(
      factor_set = factor_set,
      atlases = factor_set$atlases,
      calibration_history = history,
      uncertified = bad,
      certified = !nrow(bad)
    ),
    class = "local_atlas_theta_cloud_certification"
  )
}

fit_chart_atlas_population_model <- function(data_list,
                                             loglik_fn,
                                             population_model,
                                             theta_design = NULL,
                                             theta_root = NULL,
                                             theta_cloud = NULL,
                                             initial_proposal = NULL,
                                             local_control = list(),
                                             design_control = list(),
                                             edge_control = list(),
                                             evaluator_control = list(),
                                             calibration_control = list(),
                                             proposal_control = list(),
                                             outer_control = list(),
                                             n_cores = 1L,
                                             seed = 123L,
                                             verbose = TRUE,
                                             trace_verbose = verbose,
                                             checkpoint_file = NULL,
                                             resume_checkpoint = FALSE) {
  model <- normalize_population_model(population_model)
  if (!is.list(data_list) || !length(data_list)) {
    stop("data_list must be a non-empty list.")
  }
  if (!is.function(loglik_fn)) {
    stop("loglik_fn must be a function.")
  }
  design_control_input <- design_control
  evaluator_control_input <- evaluator_control
  design_control <- .local_atlas_merge_control(design_control, list(
    n_cloud = 4000L,
    max_anchors = 9L,
    axis_count = min(4L, model$hyper_dim),
    tail_probs = c(0.05, 0.5, 0.95),
    focus_hyper_names = NULL,
    include_axis_profiles = TRUE,
    include_cloud_profiles = TRUE,
    design_method = "hybrid",
    distance_metric = "fisher",
    coverage_audit_n = NULL,
    coverage_prob = 0.995,
    coverage_inflation = 1.25,
    refine_rounds = 0L,
    refine_points = 4L,
    refine_audit_n = NULL,
    refine_uncertified_weight = 100,
    refine_surface_se_weight = 1,
    refine_leave_gap_weight = 2,
    refine_particle_gap_weight = 2,
    refine_psis_weight = 2,
    refine_surface_se_threshold = Inf,
    refine_leave_gap_threshold = NULL,
    refine_particle_gap_threshold = NULL,
    refine_psis_threshold = NULL,
    refine_min_score = 0,
    checkpoint_local_batch_size = Inf,
    scout_only = FALSE,
    scout_min_families = 1L,
    scout_max_families = 4L,
    scout_relative_threshold = 0.25,
    scout_logz_se_weight = 1,
    scout_weak_path_weight = 1,
    scout_min_anchors_per_local = 3L,
    scout_mean_anchors_per_local = 6L,
    scout_max_anchors_per_local = NULL,
    scout_candidate_pool_size = NULL,
    stop_after_atlas_build = FALSE,
    strict_design_coverage = TRUE
  ))
  proposal_control <- .local_atlas_merge_control(proposal_control, list(
    max_components = 4L,
    core_weight = 0.9,
    tail_weight = 0.1,
    prior_weight = 0,
    df = 7,
    tail_df = 3,
    core_scale = 1.25,
    tail_scale = 4
  ))
  outer_control <- .local_atlas_merge_control(outer_control, list(
    N = 1000L,
    n_mcmc_moves = 3L,
    min_mcmc_moves = 1L,
    max_rounds = 80L,
    resample_threshold = 0.5,
    rw_scale_init = 0.8,
    verbose = verbose
  ))
  calibration_control <- .local_atlas_merge_control(calibration_control, list(
    rounds = 0L,
    max_points = 7L,
    probs = c(0.05, 0.5, 0.95),
    focus_hyper_names = design_control$focus_hyper_names,
    local_ids = NULL,
    M = as.integer(local_control$candidate_M %||% 500L),
    target_cess = local_control$target_cess %||% 0.9,
    n_mcmc_moves = as.integer(local_control$n_mcmc_moves %||% 2L),
    max_steps = as.integer(local_control$max_steps %||% 128L),
    max_updates = 12L,
    confirmation_reps = 0L,
    confirmation_M = as.integer(local_control$candidate_M %||% 500L),
    confirmation_max_sd = 1.5,
    replicate_bootstrap_B = 200L,
    max_direct_graph_z = 3,
    max_direct_graph_chart_shift = 0.35,
    max_direct_graph_existing_shift = 0.15,
    pre_outer_rounds = 1L,
    pre_outer_audit_n = NULL,
    pre_outer_max_points = 6L,
    pre_outer_max_updates = 16L,
    initial_certification_rounds = 2L,
    initial_certification_max_updates = NULL,
    initial_certification_stop_on_uncertified = TRUE
  ))
  calibration_control$initial_certification_max_updates <- as.integer(
    calibration_control$initial_certification_max_updates %||%
      calibration_control$pre_outer_max_updates
  )
  repair_executor <- "local_atlas_repair_certification_pairs"

  theta_root <- .local_chart_align_theta_one(theta_root %||% .local_atlas_default_theta(model, seed = seed), model)
  if (is.null(theta_cloud)) {
    if (!is.null(initial_proposal)) {
      theta_cloud <- theta_proposal_sample(initial_proposal, n = as.integer(design_control$n_cloud), seed = seed + 11L)
    } else {
      theta_cloud <- population_model_sample_hyper(model, n = as.integer(design_control$n_cloud), seed = seed + 11L)
    }
  }
  theta_cloud <- .as_hyper_matrix(theta_cloud, model$hyper_names, model$hyper_dim)
  theta_cloud <- rbind(theta_root, theta_cloud)
  if (is.null(theta_design)) {
    if (identical(as.character(design_control$design_method %||% ""), "adaptive_scout") ||
        isTRUE(design_control$scout_only)) {
      theta_design <- .local_atlas_unique_theta(theta_root, model)
    } else {
      theta_design <- .local_atlas_design_from_cloud(
        theta_cloud = theta_cloud,
        population_model = model,
        theta_root = theta_root,
        max_anchors = as.integer(design_control$max_anchors),
        axis_count = as.integer(design_control$axis_count),
        tail_probs = design_control$tail_probs,
        focus_hyper_names = design_control$focus_hyper_names,
        include_axis_profiles = isTRUE(design_control$include_axis_profiles),
        include_cloud_profiles = isTRUE(design_control$include_cloud_profiles),
        design_method = design_control$design_method %||% "profiles",
        distance_metric = design_control$distance_metric %||% "euclidean"
      )
    }
  } else {
    theta_design <- .local_atlas_unique_theta(rbind(theta_root, theta_design), model)
  }

  if (is.null(initial_proposal)) {
    proposal_theta <- if (nrow(theta_cloud) > nrow(theta_design)) theta_cloud else theta_design
    initial_proposal <- fit_theta_q0_proposal(
      theta = proposal_theta,
      log_weight = rep(0, nrow(proposal_theta)),
      population_model = model,
      max_components = as.integer(proposal_control$max_components),
      core_weight = proposal_control$core_weight,
      tail_weight = proposal_control$tail_weight,
      prior_weight = proposal_control$prior_weight,
      df = proposal_control$df,
      tail_df = proposal_control$tail_df,
      core_scale = proposal_control$core_scale,
      tail_scale = proposal_control$tail_scale,
      seed = seed + 17L,
      label = "chart_atlas_outer_q0"
    )
  }
  initial_proposal <- normalize_theta_proposal(initial_proposal, population_model = model)

  evaluator_control <- .local_atlas_merge_control(evaluator_control, list(
    max_chart_distance = NULL,
    min_covering_charts = 3L,
    max_prediction_range = Inf,
    distance_scale = NULL,
    se_floor = 1e-6,
    use_particle_mis = TRUE,
    require_particle_mis = TRUE,
    min_particle_mis_ess = 0.05,
    min_particle_mis_ess_abs = 50,
    max_particle_mis_psis_k = 0.7,
    max_quadratic_particle_gap = 0.05,
    sparse_chart_min_covering = 3L,
    sparse_chart_max_distance = Inf,
    max_leave_chart_out_gap = Inf,
    distance_metric = design_control$distance_metric %||% "fisher",
    surface_method = "derivative_ls",
    min_surface_charts = 2L,
    max_surface_se = Inf,
    surface_value_nugget = 0.05,
    surface_gradient_weight = 1.0,
    surface_curvature_weight = 0.2,
    surface_ridge = 1e-8,
    particle_mis_role = "estimator",
    particle_mis_batch = TRUE,
    compress_particle_mis = FALSE,
    compressed_particle_mis_K = 64L,
    compressed_particle_mis_n_theta = 96L,
    compressed_particle_mis_holdout_fraction = 0.25,
    compressed_particle_mis_max_holdout_rmse = Inf,
    compressed_particle_mis_stop_on_failure = FALSE,
    compressed_particle_mis_ridge = 1e-8,
    compressed_particle_mis_evidence_weight = 1,
    compressed_particle_mis_moment_weight = 0.05,
    compressed_particle_mis_chart_weight = 0.05,
    compressed_particle_mis_include_moments = TRUE,
    compressed_particle_mis_include_chart = TRUE,
    stop_on_uncertified = TRUE,
    use_uncertified_estimates = FALSE
  ))
  if (is.null(evaluator_control$max_chart_distance)) {
    coverage_audit_n <- as.integer(design_control$coverage_audit_n %||%
      max(nrow(theta_cloud), as.integer(outer_control$N)))
    coverage_theta <- theta_cloud
    if (!is.null(initial_proposal) && coverage_audit_n > 0L) {
      coverage_theta <- rbind(
        coverage_theta,
        theta_proposal_sample(initial_proposal, n = coverage_audit_n, seed = as.integer(seed) + 900001L)
      )
    }
    evaluator_control$max_chart_distance <- .local_atlas_nearest_design_radius(
      theta_cloud = coverage_theta,
      theta_design = theta_design,
      population_model = model,
      coverage_prob = design_control$coverage_prob,
      inflation = design_control$coverage_inflation,
      distance_metric = evaluator_control$distance_metric
    )
  }
  if (is.null(evaluator_control$distance_scale)) {
    evaluator_control$distance_scale <- evaluator_control$max_chart_distance
  }

  checkpoint <- NULL
  resume_stage <- ""
  if (isTRUE(resume_checkpoint) &&
      !is.null(checkpoint_file) &&
      nzchar(checkpoint_file) &&
      file.exists(checkpoint_file)) {
    checkpoint <- readRDS(checkpoint_file)
    resume_stage <- as.character(checkpoint$stage %||% "")
    .local_atlas_log("resuming checkpoint: ", checkpoint_file, " stage=", resume_stage, "\n", verbose = verbose)
  }
  checkpoint_state <- function(stage,
                               factor_set = NULL,
                               atlases = NULL,
                               fit = NULL,
                               atlas_build_history = list(),
                               design_certification = NULL,
                               proposal_certification = NULL,
                               calibration_history = list(),
                               initial_certification = NULL,
                               refine_round = NULL,
                               chart_design_plan = NULL) {
    plan_for_checkpoint <- chart_design_plan
    if (is.null(plan_for_checkpoint) &&
        exists("chart_design_plan", envir = parent.frame(), inherits = FALSE)) {
      plan_for_checkpoint <- get("chart_design_plan", envir = parent.frame(), inherits = FALSE)
    }
    .local_atlas_write_checkpoint(
      checkpoint_file,
      stage = stage,
      state = list(
        factor_set = factor_set,
        atlases = atlases,
        fit = fit,
        population_model = model,
        initial_proposal = initial_proposal,
        theta_root = theta_root,
        theta_design = theta_design,
        theta_cloud = theta_cloud,
        design_certification = design_certification,
        proposal_certification = proposal_certification,
        atlas_build_history = atlas_build_history,
        calibration_history = calibration_history,
        initial_certification = initial_certification,
        refine_round = refine_round,
        chart_design_plan = plan_for_checkpoint,
        settings = list(
          local_control = local_control,
          design_control = design_control,
          edge_control = edge_control,
          evaluator_control = evaluator_control,
          calibration_control = calibration_control,
          proposal_control = proposal_control,
          outer_control = outer_control,
          repair_executor = repair_executor,
          n_cores = as.integer(n_cores),
          seed = as.integer(seed)
        )
      ),
      verbose = verbose
    )
  }

  idx <- seq_along(data_list)
  atlas_build_history <- list()
  design_certification <- NULL
  proposal_certification <- NULL
  factor_set <- NULL
  atlases <- NULL
  fit <- NULL
  calibration_history <- list()
  initial_certification <- NULL
  chart_design_plan <- NULL
  completed_stages <- c("chart_design_plan", "post_atlas_build", "post_pre_outer", "post_compression", "post_outer", "complete")
  partial_atlases <- NULL
  partial_refine_round <- NA_integer_
  if (!is.null(checkpoint) && grepl("^atlas_partial_round_", resume_stage)) {
    theta_design <- checkpoint$theta_design %||% theta_design
    theta_cloud <- checkpoint$theta_cloud %||% theta_cloud
    initial_proposal <- checkpoint$initial_proposal %||% initial_proposal
    atlas_build_history <- checkpoint$atlas_build_history %||% list()
    partial_atlases <- checkpoint$atlases %||% NULL
    partial_refine_round <- as.integer(checkpoint$refine_round %||% NA_integer_)
    .local_atlas_log(
      "resuming partial atlas build: round=", partial_refine_round,
      " completed_locals=", sum(vapply(partial_atlases %||% list(), Negate(is.null), logical(1))),
      "\n",
      verbose = verbose
    )
  }
  if (!is.null(checkpoint) && resume_stage %in% completed_stages) {
    factor_set <- checkpoint$factor_set
    atlases <- checkpoint$atlases
    fit <- checkpoint$fit
    chart_design_plan <- checkpoint$chart_design_plan
    theta_design <- checkpoint$theta_design %||% theta_design
    theta_cloud <- checkpoint$theta_cloud %||% theta_cloud
    initial_proposal <- checkpoint$initial_proposal %||% initial_proposal
    design_certification <- checkpoint$design_certification
    proposal_certification <- checkpoint$proposal_certification
    atlas_build_history <- checkpoint$atlas_build_history %||% list()
    calibration_history <- checkpoint$calibration_history %||% list()
    initial_certification <- checkpoint$initial_certification
    if (!is.null(factor_set)) {
      factor_set$evaluator_control <- modifyList(
        factor_set$evaluator_control %||% list(),
        evaluator_control
      )
      factor_set <- validate_local_atlas_factor_set(factor_set)
    }
  }
  adaptive_scout <- identical(as.character(design_control$design_method %||% ""), "adaptive_scout")
  if (isTRUE(adaptive_scout) && is.null(design_control_input$stop_after_atlas_build)) {
    design_control$stop_after_atlas_build <- TRUE
  }
  if ((isTRUE(adaptive_scout) || isTRUE(design_control$scout_only)) && is.null(chart_design_plan)) {
    chart_design_plan <- build_chart_design_plan(
      data_list = data_list,
      loglik_fn = loglik_fn,
      population_model = model,
      theta_root = theta_root,
      theta_cloud = theta_cloud,
      local_control = local_control,
      design_control = design_control,
      n_cores = as.integer(n_cores),
      seed = as.integer(seed),
      verbose = verbose,
      trace_verbose = trace_verbose
    )
    atlases <- chart_design_plan$root_atlases
    checkpoint_state(
      "chart_design_plan",
      atlases = atlases,
      chart_design_plan = chart_design_plan
    )
  }
  if ((isTRUE(adaptive_scout) || isTRUE(design_control$scout_only)) && !is.null(chart_design_plan)) {
    theta_design <- chart_design_plan_theta_union(chart_design_plan, model)
    if (is.null(evaluator_control_input$max_chart_distance)) {
      coverage_audit_n <- as.integer(design_control$coverage_audit_n %||%
        max(nrow(theta_cloud), as.integer(outer_control$N)))
      coverage_theta <- theta_cloud
      if (!is.null(initial_proposal) && coverage_audit_n > 0L) {
        coverage_theta <- rbind(
          coverage_theta,
          theta_proposal_sample(initial_proposal, n = coverage_audit_n, seed = as.integer(seed) + 900001L)
        )
      }
      evaluator_control$max_chart_distance <- .local_atlas_nearest_design_radius(
        theta_cloud = coverage_theta,
        theta_design = theta_design,
        population_model = model,
        coverage_prob = design_control$coverage_prob,
        inflation = design_control$coverage_inflation,
        distance_metric = evaluator_control$distance_metric
      )
      if (is.null(evaluator_control_input$distance_scale)) {
        evaluator_control$distance_scale <- evaluator_control$max_chart_distance
      }
    }
  }
  if (isTRUE(design_control$scout_only)) {
    return(structure(
      list(
        fit = NULL,
        factor_set = NULL,
        atlases = chart_design_plan$root_atlases,
        population_model = model,
        initial_proposal = initial_proposal,
        theta_root = theta_root,
        theta_design = chart_design_plan_theta_union(chart_design_plan, model),
        theta_cloud = theta_cloud,
        chart_design_plan = chart_design_plan,
        design_certification = NULL,
        proposal_certification = NULL,
        atlas_build_history = data.frame(),
        calibration_history = data.frame(),
        graph_summary = data.frame(),
        settings = list(
          local_control = local_control,
          design_control = design_control,
          edge_control = edge_control,
          evaluator_control = evaluator_control,
          calibration_control = calibration_control,
          proposal_control = proposal_control,
          outer_control = outer_control,
          repair_executor = repair_executor,
          n_cores = as.integer(n_cores),
          seed = as.integer(seed)
        )
      ),
      class = c("chart_design_plan_fit", "chart_atlas_population_fit")
    ))
  }
  max_refine_rounds <- as.integer(design_control$refine_rounds)
  if (isTRUE(adaptive_scout) && max_refine_rounds > 0L) {
    stop("adaptive_scout production expansion uses the fixed chart_design_plan; set refine_rounds = 0 and regenerate the plan instead of running global refinement.")
  }
  if (is.null(factor_set)) {
  local_names <- names(data_list)
  if (is.null(local_names)) {
    local_names <- as.character(idx)
  }
  missing_local_names <- !nzchar(local_names) | is.na(local_names)
  local_names[missing_local_names] <- as.character(idx[missing_local_names])
  .local_atlas_log("building local atlases: locals=", length(data_list),
                   " design_anchors=", if (isTRUE(adaptive_scout)) "local_specific" else nrow(theta_design),
                   " refine_rounds=", max_refine_rounds, "\n", verbose = verbose)
  for (refine_round in seq_len(max_refine_rounds + 1L) - 1L) {
    round_start <- Sys.time()
    build_one <- function(i) {
      if (isTRUE(verbose)) {
        cat(sprintf(
          "Building chart atlas for local %s (%d/%d), design round %d\n",
          local_names[i],
          i,
          length(data_list),
          refine_round
        ))
      }
      if (isTRUE(adaptive_scout)) {
        root_atlas <- chart_design_plan$root_atlases[[local_names[i]]]
        local_design <- chart_design_plan$local_theta_designs[[local_names[i]]]
        if (is.null(root_atlas) || is.null(local_design)) {
          stop("chart_design_plan is missing root atlas or local theta design for local ", local_names[i])
        }
        build_local_atlas_from_root(
          root_atlas = root_atlas,
          data_i = data_list[[i]],
          loglik_fn = loglik_fn,
          population_model = model,
          theta_design = local_design,
          theta_root = theta_root,
          local_control = local_control,
          edge_control = edge_control,
          n_cores = 1L,
          seed = as.integer(seed) + 100003L * i + 10000019L * refine_round,
          verbose = isTRUE(trace_verbose)
        )
      } else {
        build_local_atlas(
          local_id = local_names[i],
          data_i = data_list[[i]],
          loglik_fn = loglik_fn,
          population_model = model,
          theta_design = theta_design,
          theta_root = theta_root,
          local_control = local_control,
          edge_control = edge_control,
          n_cores = 1L,
          seed = as.integer(seed) + 100003L * i + 10000019L * refine_round,
          verbose = isTRUE(trace_verbose)
        )
      }
    }
    batch_size <- suppressWarnings(as.integer(design_control$checkpoint_local_batch_size))
    if (!is.finite(batch_size) || batch_size <= 0L || is.null(checkpoint_file) || !nzchar(checkpoint_file)) {
      batch_size <- length(idx)
    }
    if (!is.null(partial_atlases) && identical(as.integer(partial_refine_round), as.integer(refine_round))) {
      atlases <- partial_atlases
      length(atlases) <- length(data_list)
    } else {
      atlases <- vector("list", length(data_list))
    }
    names(atlases) <- local_names
    missing_idx <- idx[vapply(atlases, is.null, logical(1))]
    if (length(missing_idx)) {
      batches <- split(missing_idx, ceiling(seq_along(missing_idx) / batch_size))
      for (batch_id in seq_along(batches)) {
        batch_idx <- batches[[batch_id]]
        .local_atlas_log(
          "atlas round ", refine_round,
          " local batch ", batch_id, "/", length(batches),
          ": locals ", min(batch_idx), "-", max(batch_idx),
          " (", length(batch_idx), ")\n",
          verbose = verbose
        )
        built <- if (as.integer(n_cores) <= 1L || length(batch_idx) <= 1L) {
          lapply(batch_idx, build_one)
        } else {
          parallel::mclapply(batch_idx, build_one, mc.cores = as.integer(min(n_cores, length(batch_idx))))
        }
        atlases[batch_idx] <- built
        checkpoint_state(
          paste0("atlas_partial_round_", refine_round),
          atlases = atlases,
          atlas_build_history = atlas_build_history,
          design_certification = design_certification,
          proposal_certification = proposal_certification,
          refine_round = refine_round
        )
      }
    }
    partial_atlases <- NULL
    partial_refine_round <- NA_integer_

    factor_set <- build_local_atlas_factor_set(
      atlases = atlases,
      population_model = model,
      max_chart_distance = evaluator_control$max_chart_distance,
      min_covering_charts = evaluator_control$min_covering_charts,
      max_prediction_range = evaluator_control$max_prediction_range,
      distance_scale = evaluator_control$distance_scale,
      se_floor = evaluator_control$se_floor,
      use_particle_mis = evaluator_control$use_particle_mis,
      require_particle_mis = evaluator_control$require_particle_mis,
      min_particle_mis_ess = evaluator_control$min_particle_mis_ess,
      min_particle_mis_ess_abs = evaluator_control$min_particle_mis_ess_abs,
      max_particle_mis_psis_k = evaluator_control$max_particle_mis_psis_k,
      max_quadratic_particle_gap = evaluator_control$max_quadratic_particle_gap,
      sparse_chart_min_covering = evaluator_control$sparse_chart_min_covering,
      sparse_chart_max_distance = evaluator_control$sparse_chart_max_distance,
      max_leave_chart_out_gap = evaluator_control$max_leave_chart_out_gap,
      distance_metric = evaluator_control$distance_metric,
      surface_method = evaluator_control$surface_method,
      min_surface_charts = evaluator_control$min_surface_charts,
      max_surface_se = evaluator_control$max_surface_se,
      surface_value_nugget = evaluator_control$surface_value_nugget,
      surface_gradient_weight = evaluator_control$surface_gradient_weight,
      surface_curvature_weight = evaluator_control$surface_curvature_weight,
      surface_ridge = evaluator_control$surface_ridge,
      particle_mis_role = evaluator_control$particle_mis_role,
      particle_mis_batch = evaluator_control$particle_mis_batch,
      stop_on_uncertified = evaluator_control$stop_on_uncertified,
      use_uncertified_estimates = evaluator_control$use_uncertified_estimates
    )
	    design_certification <- if (isTRUE(adaptive_scout)) {
	      local_atlas_certification_summary_by_local_design(
	        factor_set,
	        local_theta_designs = chart_design_plan$local_theta_designs,
	        n_cores = as.integer(n_cores)
      )
    } else {
      local_atlas_certification_summary(
        factor_set,
        theta = theta_design,
	        n_cores = as.integer(n_cores)
	      )
	    }
	    if (isTRUE(design_control$strict_design_coverage) &&
	        any(design_certification$uncertified_fraction > 0)) {
	      design_repair <- .local_atlas_design_repair_spec(
	        factor_set = factor_set,
	        design_certification = design_certification,
	        population_model = model,
	        theta_design = theta_design,
	        local_theta_designs = if (isTRUE(adaptive_scout)) chart_design_plan$local_theta_designs else NULL
	      )
	      if (nrow(design_repair$pairs)) {
	        .local_atlas_log(
	          "design coverage repair: uncertified_pairs=", nrow(design_repair$pairs),
	          " locals=", length(unique(design_repair$pairs$local_pos)), "\n",
	          verbose = verbose
	        )
	        calibration <- local_atlas_repair_certification_pairs(
	          factor_set = factor_set,
	          theta = design_repair$theta,
	          data_list = data_list,
	          loglik_fn = loglik_fn,
	          theta_weights = rep(1 / nrow(design_repair$theta), nrow(design_repair$theta)),
	          candidate_pairs = design_repair$pairs,
	          local_ids = unique(design_repair$pairs$local_pos),
	          M = as.integer(calibration_control$M),
	          target_cess = calibration_control$target_cess,
	          n_mcmc_moves = as.integer(calibration_control$n_mcmc_moves),
	          max_steps = as.integer(calibration_control$max_steps),
	          max_updates = nrow(design_repair$pairs),
	          direct_confirmation_reps = as.integer(calibration_control$confirmation_reps),
	          direct_confirmation_M = as.integer(calibration_control$confirmation_M),
	          direct_confirmation_max_sd = calibration_control$confirmation_max_sd,
	          replicate_bootstrap_B = as.integer(calibration_control$replicate_bootstrap_B),
	          max_direct_graph_z = calibration_control$max_direct_graph_z,
	          max_direct_graph_chart_shift = calibration_control$max_direct_graph_chart_shift,
	          max_direct_graph_existing_shift = calibration_control$max_direct_graph_existing_shift,
	          local_control = local_control,
	          edge_control = edge_control,
	          seed = as.integer(seed) + 29000017L + 10000019L * refine_round,
	          verbose = isTRUE(trace_verbose)
	        )
	        factor_set <- calibration$factor_set
	        atlases <- factor_set$atlases
	        calibration_history[[length(calibration_history) + 1L]] <- data.frame(
	          calibration_phase = "design",
	          calibration_round = refine_round,
	          n_uncertified_pairs = nrow(design_repair$pairs),
	          calibration$probes,
	          check.names = FALSE
	        )
	        design_certification <- if (isTRUE(adaptive_scout)) {
	          local_atlas_certification_summary_by_local_design(
	            factor_set,
	            local_theta_designs = chart_design_plan$local_theta_designs,
	            n_cores = as.integer(n_cores)
	          )
	        } else {
	          local_atlas_certification_summary(
	            factor_set,
	            theta = theta_design,
	            n_cores = as.integer(n_cores)
	          )
	        }
	        .local_atlas_log(
	          "design coverage repair finished: remaining_uncertified_pairs=",
	          nrow(attr(design_certification, "bad_rows") %||% data.frame()), "\n",
	          verbose = verbose
	        )
	        checkpoint_state(
	          paste0("design_repair_round_", refine_round),
	          factor_set = factor_set,
	          atlases = atlases,
	          atlas_build_history = atlas_build_history,
	          design_certification = design_certification,
	          proposal_certification = proposal_certification,
	          calibration_history = calibration_history,
	          refine_round = refine_round
	        )
	      }
	    }
	    if (isTRUE(design_control$strict_design_coverage) &&
	        any(design_certification$uncertified_fraction > 0)) {
	      bad_table <- attr(design_certification, "bad_rows")
      if (is.null(bad_table) || !is.data.frame(bad_table)) {
        audit_parts <- .local_atlas_factor_set_by_local(factor_set, theta_design, n_cores = as.integer(n_cores))
        bad_rows <- list()
        for (local_name in names(audit_parts)) {
          part <- audit_parts[[local_name]]
          bad <- which(part$status != "certified")
          if (!length(bad)) next
          bad_rows[[length(bad_rows) + 1L]] <- data.frame(
            local = local_name,
            theta_row = bad,
            reason = part$reason[bad],
            nearest_charts = part$nearest_charts[bad],
            check.names = FALSE
          )
        }
        bad_table <- if (length(bad_rows)) do.call(rbind, bad_rows) else data.frame()
      }
      shown <- utils::head(bad_table, 8L)
      stop(
        "chart atlas failed to certify all design anchors:\n",
        paste(
          sprintf(
            "local=%s theta_row=%s reason=%s nearest=%s",
            shown$local,
            shown$theta_row,
            shown$reason,
            shown$nearest_charts
          ),
          collapse = "\n"
        )
      )
    }

    if (refine_round >= max_refine_rounds) {
      proposal_certification <- data.frame(
        refine_round = refine_round,
        n_theta = 0L,
        failed_theta = NA_integer_,
        max_failed_locals = NA_integer_,
        max_refine_score = NA_real_,
        max_surface_excess = NA_real_,
        max_leave_gap_excess = NA_real_,
        max_particle_gap_excess = NA_real_,
        max_psis_excess = NA_real_,
        check.names = FALSE
      )
      atlas_build_history[[length(atlas_build_history) + 1L]] <- proposal_certification
      .local_atlas_log(
        "atlas design round ", refine_round,
        " finished in ", round(as.numeric(difftime(Sys.time(), round_start, units = "mins")), 3),
        " min | refinement audit skipped\n",
        verbose = verbose
      )
      break
    }

    refine_audit_n <- as.integer(design_control$refine_audit_n %||% outer_control$N)
    audit_theta <- theta_proposal_sample(
      initial_proposal,
      n = refine_audit_n,
      seed = as.integer(seed) + 900001L + refine_round
    )
    audit_parts <- .local_atlas_factor_set_by_local(factor_set, audit_theta, n_cores = as.integer(n_cores))
    bad_count <- integer(nrow(audit_theta))
    surface_excess <- numeric(nrow(audit_theta))
    leave_gap_excess <- numeric(nrow(audit_theta))
    particle_gap_excess <- numeric(nrow(audit_theta))
    psis_excess <- numeric(nrow(audit_theta))
    surface_se_threshold <- as.numeric(design_control$refine_surface_se_threshold %||% Inf)
    leave_gap_threshold <- as.numeric(design_control$refine_leave_gap_threshold %||% evaluator_control$max_leave_chart_out_gap)
    particle_gap_threshold <- as.numeric(design_control$refine_particle_gap_threshold %||% evaluator_control$max_quadratic_particle_gap)
    psis_threshold <- as.numeric(design_control$refine_psis_threshold %||% evaluator_control$max_particle_mis_psis_k)
    for (part in audit_parts) {
      bad_count <- bad_count + as.integer(part$status != "certified")
      if ("surface_se" %in% names(part) && is.finite(surface_se_threshold)) {
        x <- as.numeric(part$surface_se)
        x[!is.finite(x)] <- 0
        surface_excess <- surface_excess + pmax(x - surface_se_threshold, 0)
      }
      if ("leave_chart_out_gap" %in% names(part) && is.finite(leave_gap_threshold)) {
        x <- as.numeric(part$leave_chart_out_gap)
        x[!is.finite(x)] <- 0
        leave_gap_excess <- leave_gap_excess + pmax(x - leave_gap_threshold, 0)
      }
      if ("quadratic_particle_gap" %in% names(part) && is.finite(particle_gap_threshold)) {
        x <- as.numeric(part$quadratic_particle_gap)
        x[!is.finite(x)] <- 0
        particle_gap_excess <- particle_gap_excess + pmax(x - particle_gap_threshold, 0)
      }
      if ("particle_mis_psis_k" %in% names(part) && is.finite(psis_threshold)) {
        x <- as.numeric(part$particle_mis_psis_k)
        x[!is.finite(x)] <- -Inf
        psis_excess <- psis_excess + pmax(x - psis_threshold, 0)
      }
    }
    refine_score <- as.numeric(design_control$refine_uncertified_weight) * bad_count +
      as.numeric(design_control$refine_surface_se_weight) * surface_excess +
      as.numeric(design_control$refine_leave_gap_weight) * leave_gap_excess +
      as.numeric(design_control$refine_particle_gap_weight) * particle_gap_excess +
      as.numeric(design_control$refine_psis_weight) * psis_excess
    proposal_certification <- data.frame(
      refine_round = refine_round,
      n_theta = nrow(audit_theta),
      failed_theta = sum(bad_count > 0),
      max_failed_locals = if (length(bad_count)) max(bad_count) else 0L,
      max_refine_score = if (length(refine_score)) max(refine_score) else 0,
      max_surface_excess = if (length(surface_excess)) max(surface_excess) else 0,
      max_leave_gap_excess = if (length(leave_gap_excess)) max(leave_gap_excess) else 0,
      max_particle_gap_excess = if (length(particle_gap_excess)) max(particle_gap_excess) else 0,
      max_psis_excess = if (length(psis_excess)) max(psis_excess) else 0,
      check.names = FALSE
    )
    atlas_build_history[[length(atlas_build_history) + 1L]] <- proposal_certification
    .local_atlas_log(
      "atlas design round ", refine_round,
      " finished in ", round(as.numeric(difftime(Sys.time(), round_start, units = "mins")), 3),
      " min | failed_theta=", proposal_certification$failed_theta,
      " max_failed_locals=", proposal_certification$max_failed_locals, "\n",
      verbose = verbose
    )
    refine_candidates <- which(refine_score > as.numeric(design_control$refine_min_score))
    if (!length(refine_candidates)) {
      break
    }
    ranked_bad <- refine_candidates[order(refine_score[refine_candidates], decreasing = TRUE)]
    add_n <- min(length(ranked_bad), as.integer(design_control$refine_points))
    theta_design <- .local_atlas_unique_theta(
      rbind(theta_design, audit_theta[ranked_bad[seq_len(add_n)], , drop = FALSE]),
      model
    )
  }
  checkpoint_state(
    "post_atlas_build",
    factor_set = factor_set,
    atlases = atlases,
    atlas_build_history = atlas_build_history,
    design_certification = design_certification,
    proposal_certification = proposal_certification
  )
  } else {
    .local_atlas_log("using checkpointed atlas build stage\n", verbose = verbose)
  }
  if (isTRUE(design_control$stop_after_atlas_build)) {
    return(structure(
      list(
        fit = fit,
        factor_set = factor_set,
        atlases = atlases,
        population_model = model,
        initial_proposal = initial_proposal,
        theta_root = theta_root,
        theta_design = theta_design,
        theta_cloud = theta_cloud,
        chart_design_plan = chart_design_plan,
        design_certification = design_certification,
        proposal_certification = proposal_certification,
        atlas_build_history = .local_atlas_rbind_fill(atlas_build_history),
        calibration_history = .local_atlas_rbind_fill(calibration_history),
        initial_certification = initial_certification,
        graph_summary = if (is.null(factor_set)) data.frame() else local_atlas_graph_summary(factor_set),
        compression_summary = if (is.null(factor_set)) data.frame() else factor_set$compression_summary %||% data.frame(),
        settings = list(
          local_control = local_control,
          design_control = design_control,
          edge_control = edge_control,
          evaluator_control = evaluator_control,
          calibration_control = calibration_control,
          proposal_control = proposal_control,
          outer_control = outer_control,
          repair_executor = repair_executor,
          n_cores = as.integer(n_cores),
          seed = as.integer(seed)
        )
      ),
      class = c("chart_atlas_factor_set_fit", "chart_atlas_population_fit")
    ))
  }

  pre_outer_calibration_history <- calibration_history
  if (is.data.frame(pre_outer_calibration_history)) {
    pre_outer_calibration_history <- list(pre_outer_calibration_history)
  }
  if (!is.list(pre_outer_calibration_history)) {
    pre_outer_calibration_history <- list()
  }
  if (resume_stage %in% c("post_pre_outer", "post_compression", "post_outer", "complete")) {
    .local_atlas_log("using checkpointed pre-outer certification stage\n", verbose = verbose)
  } else if (as.integer(calibration_control$pre_outer_rounds %||% 0L) > 0L) {
    pre_outer <- local_atlas_pre_outer_certify(
      factor_set = factor_set,
      data_list = data_list,
      loglik_fn = loglik_fn,
      initial_proposal = initial_proposal,
      local_control = local_control,
      edge_control = edge_control,
      calibration_control = calibration_control,
      outer_control = outer_control,
      design_control = design_control,
      calibration_history = pre_outer_calibration_history,
      n_cores = as.integer(n_cores),
      seed = as.integer(seed),
      verbose = verbose,
      trace_verbose = trace_verbose,
      checkpoint_callback = function(stage, factor_set, atlases, calibration_history) {
        checkpoint_state(
          stage,
          factor_set = factor_set,
          atlases = atlases,
          atlas_build_history = atlas_build_history,
          design_certification = design_certification,
          proposal_certification = proposal_certification,
          calibration_history = calibration_history
        )
      }
    )
    factor_set <- pre_outer$factor_set
    atlases <- pre_outer$atlases
    pre_outer_calibration_history <- pre_outer$calibration_history
    checkpoint_state(
      "post_pre_outer",
      factor_set = factor_set,
      atlases = atlases,
      atlas_build_history = atlas_build_history,
      design_certification = design_certification,
      proposal_certification = proposal_certification,
      calibration_history = pre_outer_calibration_history,
      initial_certification = initial_certification
    )
  }
  if (!resume_stage %in% c("post_pre_outer", "post_compression", "post_outer", "complete") &&
      as.integer(calibration_control$initial_certification_rounds %||% 0L) > 0L) {
    initial_theta <- theta_proposal_sample(
      initial_proposal,
      n = as.integer(outer_control$N),
      seed = as.integer(seed) + 900001L
    )
    initial_certification <- local_atlas_certify_theta_cloud(
      factor_set = factor_set,
      theta = initial_theta,
      data_list = data_list,
      loglik_fn = loglik_fn,
      local_control = local_control,
      edge_control = edge_control,
      calibration_control = calibration_control,
      theta_weights = rep(1 / as.integer(outer_control$N), as.integer(outer_control$N)),
      max_rounds = as.integer(calibration_control$initial_certification_rounds),
      max_updates = as.integer(calibration_control$initial_certification_max_updates),
      n_cores = as.integer(n_cores),
      seed = as.integer(seed) + 1900003L,
      verbose = verbose,
      trace_verbose = trace_verbose
    )
    factor_set <- initial_certification$factor_set
    atlases <- initial_certification$atlases
    if (length(initial_certification$calibration_history)) {
      pre_outer_calibration_history <- c(
        pre_outer_calibration_history %||% list(),
        initial_certification$calibration_history
      )
    }
    checkpoint_state(
      "post_pre_outer",
      factor_set = factor_set,
      atlases = atlases,
      atlas_build_history = atlas_build_history,
      design_certification = design_certification,
      proposal_certification = proposal_certification,
      calibration_history = pre_outer_calibration_history,
      initial_certification = initial_certification
    )
    if (!isTRUE(initial_certification$certified) &&
        isTRUE(calibration_control$initial_certification_stop_on_uncertified)) {
      shown <- utils::head(initial_certification$uncertified, 8L)
      stop(
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

  compress_factor_set_for_outer <- function(stage_label) {
    if (!isTRUE(evaluator_control$compress_particle_mis)) {
      factor_set$evaluator_control$use_compressed_particle_mis <- FALSE
      factor_set$evaluator_control$require_compressed_particle_mis <- FALSE
      return(factor_set)
    }
    compression_theta <- .local_atlas_compression_theta_design(
      theta_design = theta_design,
      theta_cloud = theta_cloud,
      initial_proposal = initial_proposal,
      population_model = model,
      theta_root = theta_root,
      n_theta = as.integer(evaluator_control$compressed_particle_mis_n_theta),
      distance_metric = evaluator_control$distance_metric,
      seed = as.integer(seed) + 7100003L
    )
    out <- compress_local_atlas_factor_set(
      factor_set = factor_set,
      theta = compression_theta,
      K = as.integer(evaluator_control$compressed_particle_mis_K),
      holdout_fraction = evaluator_control$compressed_particle_mis_holdout_fraction,
      ridge = evaluator_control$compressed_particle_mis_ridge,
      evidence_weight = evaluator_control$compressed_particle_mis_evidence_weight,
      moment_weight = evaluator_control$compressed_particle_mis_moment_weight,
      chart_weight = evaluator_control$compressed_particle_mis_chart_weight,
      include_moments = evaluator_control$compressed_particle_mis_include_moments,
      include_chart = evaluator_control$compressed_particle_mis_include_chart,
      max_holdout_rmse = evaluator_control$compressed_particle_mis_max_holdout_rmse,
      stop_on_failure = evaluator_control$compressed_particle_mis_stop_on_failure,
      n_cores = as.integer(n_cores),
      seed = as.integer(seed) + 7200003L,
      verbose = verbose
    )
    summary <- local_atlas_compression_summary(out)
    if (nrow(summary)) {
      .local_atlas_log(
        stage_label,
        " compression: median raw=", stats::median(summary$raw_particles, na.rm = TRUE),
        " median selected=", stats::median(summary$selected_particles, na.rm = TRUE),
        " median ratio=", round(stats::median(summary$compression_ratio, na.rm = TRUE), 1),
        " median holdout RMSE=", round(stats::median(summary$holdout_rmse, na.rm = TRUE), 4),
        "\n",
        verbose = verbose
      )
    }
    out
  }

  if (isTRUE(evaluator_control$compress_particle_mis)) {
    if (resume_stage %in% c("post_compression", "post_outer", "complete")) {
      .local_atlas_log("using checkpointed local particle compression stage\n", verbose = verbose)
      factor_set$evaluator_control$use_compressed_particle_mis <- TRUE
      factor_set$evaluator_control$require_compressed_particle_mis <- TRUE
      factor_set <- validate_local_atlas_factor_set(factor_set)
    } else {
      factor_set <- compress_factor_set_for_outer("pre-outer")
      atlases <- factor_set$atlases
      checkpoint_state(
        "post_compression",
        factor_set = factor_set,
        atlases = atlases,
        atlas_build_history = atlas_build_history,
        design_certification = design_certification,
        proposal_certification = proposal_certification,
        calibration_history = pre_outer_calibration_history,
        initial_certification = initial_certification
      )
    }
  }

  if (is.null(fit) || !resume_stage %in% c("post_outer", "complete")) {
  .local_atlas_log("starting outer SMC: N=", as.integer(outer_control$N), "\n", verbose = verbose)
  outer_start <- Sys.time()
  fit <- outer_population_smc(
    factor_set = factor_set,
    N = as.integer(outer_control$N),
    initial_proposal = initial_proposal,
    resample_threshold = outer_control$resample_threshold,
    n_mcmc_moves = as.integer(outer_control$n_mcmc_moves),
    min_mcmc_moves = as.integer(outer_control$min_mcmc_moves),
    max_rounds = as.integer(outer_control$max_rounds),
    rw_scale_init = outer_control$rw_scale_init,
    n_cores = as.integer(n_cores),
    seed = as.integer(seed) + 900001L,
    verbose = isTRUE(outer_control$verbose)
  )
  .local_atlas_log(
    "outer SMC finished in ", round(as.numeric(difftime(Sys.time(), outer_start, units = "mins")), 3),
    " min\n",
    verbose = verbose
  )
  checkpoint_state(
    "post_outer",
    factor_set = factor_set,
    atlases = atlases,
    fit = fit,
    atlas_build_history = atlas_build_history,
    design_certification = design_certification,
    proposal_certification = proposal_certification,
    calibration_history = pre_outer_calibration_history,
    initial_certification = initial_certification
  )
  } else {
    .local_atlas_log("using checkpointed outer SMC stage\n", verbose = verbose)
  }

  calibration_history <- pre_outer_calibration_history
  if (is.data.frame(calibration_history)) {
    calibration_history <- list(calibration_history)
  }
  max_calibration_rounds <- as.integer(calibration_control$rounds)
  if (resume_stage %in% "complete") {
    .local_atlas_log("using checkpointed posterior calibration stage\n", verbose = verbose)
  } else if (max_calibration_rounds > 0L) {
    for (calibration_round in seq_len(max_calibration_rounds)) {
      posterior_calibration_start <- Sys.time()
      .local_atlas_log("posterior calibration round ", calibration_round, "\n", verbose = verbose)
      calibration_design <- local_atlas_select_theta_profiles(
        theta = fit$theta,
        population_model = model,
        weights = fit$w,
        focus_hyper_names = calibration_control$focus_hyper_names,
        probs = calibration_control$probs,
        max_points = as.integer(calibration_control$max_points),
        return_rows = TRUE
      )
      calibration_theta <- calibration_design$theta
      posterior_repair_local_ids <- calibration_control$local_ids %||% seq_along(data_list)
      posterior_candidate_pairs <- expand.grid(
        local_pos = posterior_repair_local_ids,
        theta_row = seq_len(nrow(calibration_theta)),
        KEEP.OUT.ATTRS = FALSE
      )
      calibration <- local_atlas_repair_certification_pairs(
        factor_set = factor_set,
        theta = calibration_theta,
        data_list = data_list,
        loglik_fn = loglik_fn,
        theta_weights = calibration_design$weights,
        candidate_pairs = posterior_candidate_pairs,
        local_ids = posterior_repair_local_ids,
        M = as.integer(calibration_control$M),
        target_cess = calibration_control$target_cess,
        n_mcmc_moves = as.integer(calibration_control$n_mcmc_moves),
        max_steps = as.integer(calibration_control$max_steps),
        max_updates = as.integer(calibration_control$max_updates),
        direct_confirmation_reps = as.integer(calibration_control$confirmation_reps),
        direct_confirmation_M = as.integer(calibration_control$confirmation_M),
        direct_confirmation_max_sd = calibration_control$confirmation_max_sd,
        replicate_bootstrap_B = as.integer(calibration_control$replicate_bootstrap_B),
        max_direct_graph_z = calibration_control$max_direct_graph_z,
        max_direct_graph_chart_shift = calibration_control$max_direct_graph_chart_shift,
        max_direct_graph_existing_shift = calibration_control$max_direct_graph_existing_shift,
        local_control = local_control,
        edge_control = edge_control,
        seed = as.integer(seed) + 19000019L * calibration_round,
        verbose = isTRUE(trace_verbose)
      )
      calibration_history[[length(calibration_history) + 1L]] <- data.frame(
        calibration_phase = "posterior",
        calibration_round = calibration_round,
        calibration$probes,
        check.names = FALSE
      )
      factor_set <- calibration$factor_set
      atlases <- factor_set$atlases
      if (isTRUE(evaluator_control$compress_particle_mis)) {
        factor_set <- compress_factor_set_for_outer(paste0("posterior calibration round ", calibration_round))
        atlases <- factor_set$atlases
      }
      .local_atlas_log(
        "posterior calibration round ", calibration_round,
        " finished in ", round(as.numeric(difftime(Sys.time(), posterior_calibration_start, units = "mins")), 3),
        " min | activated=", calibration$n_activated,
        " selected=", calibration$n_selected, "\n",
        verbose = verbose
      )
      checkpoint_state(
        paste0("posterior_calibration_round_", calibration_round),
        factor_set = factor_set,
        atlases = atlases,
        fit = fit,
        atlas_build_history = atlas_build_history,
        design_certification = design_certification,
        proposal_certification = proposal_certification,
        calibration_history = calibration_history
      )
      if (isTRUE(calibration$n_activated > 0L)) {
        .local_atlas_log(
          "posterior calibration updated local atlases; outer SMC rerun skipped. Run an explicit outer rerun from the saved checkpoint.\n",
          verbose = verbose
        )
      }
      break
    }
  }

  checkpoint_state(
    "complete",
    factor_set = factor_set,
    atlases = atlases,
    fit = fit,
    atlas_build_history = atlas_build_history,
    design_certification = design_certification,
    proposal_certification = proposal_certification,
    calibration_history = calibration_history,
    initial_certification = initial_certification
  )

  structure(
    list(
      fit = fit,
      factor_set = factor_set,
      atlases = atlases,
      population_model = model,
      initial_proposal = initial_proposal,
      theta_root = theta_root,
      theta_design = theta_design,
      theta_cloud = theta_cloud,
      chart_design_plan = chart_design_plan,
      design_certification = design_certification,
      proposal_certification = proposal_certification,
      atlas_build_history = .local_atlas_rbind_fill(atlas_build_history),
      calibration_history = .local_atlas_rbind_fill(calibration_history),
      initial_certification = initial_certification,
      graph_summary = local_atlas_graph_summary(factor_set),
      compression_summary = factor_set$compression_summary %||% data.frame(),
      settings = list(
        local_control = local_control,
        design_control = design_control,
        edge_control = edge_control,
        evaluator_control = evaluator_control,
        calibration_control = calibration_control,
        proposal_control = proposal_control,
        outer_control = outer_control,
        repair_executor = repair_executor,
        n_cores = as.integer(n_cores),
        seed = as.integer(seed)
      )
    ),
    class = "chart_atlas_population_fit"
  )
}

.local_atlas_posterior_quantile_distance <- function(x,
                                                     y,
                                                     probs = seq(0.01, 0.99, length.out = 99L)) {
  qx <- stats::quantile(as.numeric(x), probs = probs, names = FALSE, type = 8, na.rm = TRUE)
  qy <- stats::quantile(as.numeric(y), probs = probs, names = FALSE, type = 8, na.rm = TRUE)
  mean(abs(qx - qy))
}

local_atlas_compare_posterior_draws <- function(reference_draws, workflow_draws) {
  common <- intersect(names(reference_draws), names(workflow_draws))
  if (!length(common)) {
    stop("reference_draws and workflow_draws have no common columns.")
  }
  rows <- lapply(common, function(nm) {
    ref <- as.numeric(reference_draws[[nm]])
    wf <- as.numeric(workflow_draws[[nm]])
    q_ref <- stats::quantile(ref, probs = c(0.05, 0.5, 0.95), names = FALSE, type = 8, na.rm = TRUE)
    q_wf <- stats::quantile(wf, probs = c(0.05, 0.5, 0.95), names = FALSE, type = 8, na.rm = TRUE)
    ref_sd <- stats::sd(ref, na.rm = TRUE)
    wf_sd <- stats::sd(wf, na.rm = TRUE)
    mean_error <- mean(wf, na.rm = TRUE) - mean(ref, na.rm = TRUE)
    standardized_mean_error <- mean_error / max(ref_sd, .Machine$double.eps)
    q_wasserstein <- .local_atlas_posterior_quantile_distance(wf, ref)
    scaled_q_wasserstein <- q_wasserstein / max(ref_sd, .Machine$double.eps)
    sd_ratio <- wf_sd / max(ref_sd, .Machine$double.eps)
    data.frame(
      parameter = nm,
      reference_mean = mean(ref, na.rm = TRUE),
      workflow_mean = mean(wf, na.rm = TRUE),
      mean_error = mean_error,
      standardized_mean_error = standardized_mean_error,
      reference_sd = ref_sd,
      workflow_sd = wf_sd,
      sd_ratio = sd_ratio,
      q05_error = q_wf[1L] - q_ref[1L],
      q50_error = q_wf[2L] - q_ref[2L],
      q95_error = q_wf[3L] - q_ref[3L],
      q_wasserstein = q_wasserstein,
      scaled_q_wasserstein = scaled_q_wasserstein,
      reference_inside_workflow_q05_q95 = mean(ref >= q_wf[1L] & ref <= q_wf[3L], na.rm = TRUE),
      workflow_inside_reference_q05_q95 = mean(wf >= q_ref[1L] & wf <= q_ref[3L], na.rm = TRUE),
      shape_error = abs(standardized_mean_error) +
        scaled_q_wasserstein +
        abs(log(max(sd_ratio, .Machine$double.eps))),
      check.names = FALSE
    )
  })
  do.call(rbind, rows)
}

.local_atlas_metric_summary <- function(comparison, parameters = NULL) {
  if (!is.null(parameters)) {
    comparison <- comparison[comparison$parameter %in% parameters, , drop = FALSE]
  }
  if (!nrow(comparison)) {
    return(data.frame(
      n_parameters = 0L,
      max_abs_standardized_mean_error = NA_real_,
      mean_abs_standardized_mean_error = NA_real_,
      max_scaled_q_wasserstein = NA_real_,
      mean_scaled_q_wasserstein = NA_real_,
      max_shape_error = NA_real_,
      mean_shape_error = NA_real_,
      check.names = FALSE
    ))
  }
  data.frame(
    n_parameters = nrow(comparison),
    max_abs_standardized_mean_error = max(abs(comparison$standardized_mean_error), na.rm = TRUE),
    mean_abs_standardized_mean_error = mean(abs(comparison$standardized_mean_error), na.rm = TRUE),
    max_scaled_q_wasserstein = max(comparison$scaled_q_wasserstein, na.rm = TRUE),
    mean_scaled_q_wasserstein = mean(comparison$scaled_q_wasserstein, na.rm = TRUE),
    max_shape_error = max(comparison$shape_error, na.rm = TRUE),
    mean_shape_error = mean(comparison$shape_error, na.rm = TRUE),
    check.names = FALSE
  )
}

local_atlas_compare_to_baseline <- function(reference_draws,
                                            workflow_draws,
                                            baseline_draws = NULL,
                                            baseline_comparison = NULL) {
  if (is.null(baseline_comparison)) {
    if (is.null(baseline_draws)) {
      stop("Provide baseline_draws or baseline_comparison.")
    }
    baseline_comparison <- local_atlas_compare_posterior_draws(reference_draws, baseline_draws)
  }
  workflow_comparison <- local_atlas_compare_posterior_draws(reference_draws, workflow_draws)
  common <- intersect(workflow_comparison$parameter, baseline_comparison$parameter)
  if (!length(common)) {
    stop("workflow and baseline comparisons have no common parameters.")
  }
  wf <- workflow_comparison[match(common, workflow_comparison$parameter), , drop = FALSE]
  base <- baseline_comparison[match(common, baseline_comparison$parameter), , drop = FALSE]
  data.frame(
    parameter = common,
    workflow_abs_standardized_mean_error = abs(wf$standardized_mean_error),
    baseline_abs_standardized_mean_error = abs(base$standardized_mean_error),
    standardized_mean_error_improvement =
      abs(base$standardized_mean_error) - abs(wf$standardized_mean_error),
    workflow_scaled_q_wasserstein = wf$scaled_q_wasserstein,
    baseline_scaled_q_wasserstein = base$scaled_q_wasserstein,
    scaled_q_wasserstein_improvement = base$scaled_q_wasserstein - wf$scaled_q_wasserstein,
    workflow_shape_error = wf$shape_error,
    baseline_shape_error = base$shape_error,
    shape_error_improvement = base$shape_error - wf$shape_error,
    check.names = FALSE
  )
}

local_atlas_fresh_probe_rmse <- function(gate_or_probes) {
  probes <- if (is.data.frame(gate_or_probes)) gate_or_probes else gate_or_probes$fresh_probes
  if (is.null(probes) || !nrow(probes) || !"delta_fresh_minus_atlas" %in% names(probes)) {
    return(Inf)
  }
  delta <- as.numeric(probes$delta_fresh_minus_atlas)
  delta <- delta[is.finite(delta)]
  if (!length(delta)) Inf else sqrt(mean(delta^2))
}

local_atlas_ensemble_weights_from_fresh_probes <- function(gates_or_probes,
                                                           power = 2,
                                                           floor = 1e-6) {
  rmse <- vapply(gates_or_probes, local_atlas_fresh_probe_rmse, numeric(1))
  score <- 1 / pmax(rmse, as.numeric(floor))^as.numeric(power)
  if (any(!is.finite(score)) || sum(score) <= 0) {
    score <- rep(1, length(gates_or_probes))
  }
  data.frame(
    member = names(gates_or_probes) %||% paste0("member_", seq_along(gates_or_probes)),
    fresh_probe_rmse = rmse,
    ensemble_weight = score / sum(score),
    check.names = FALSE
  )
}

local_atlas_ensemble_draws <- function(draws_list,
                                       weights,
                                       n_draws = min(vapply(draws_list, nrow, integer(1))),
                                       seed = 123L) {
  if (!is.list(draws_list) || !length(draws_list)) {
    stop("draws_list must be a non-empty list.")
  }
  weights <- .local_chart_normalize_weights(weights, length(draws_list))
  n_draws <- as.integer(n_draws)
  if (n_draws <= 0L) {
    stop("n_draws must be positive.")
  }
  draw_names <- lapply(draws_list, function(x) names(as.data.frame(x, check.names = FALSE)))
  if (any(vapply(draw_names, function(x) !identical(x, draw_names[[1L]]), logical(1)))) {
    stop("All draw sets must have identical columns in identical order.")
  }
  set.seed(as.integer(seed))
  member <- sample(seq_along(draws_list), size = n_draws, replace = TRUE, prob = weights)
  out <- vector("list", length(draws_list))
  for (k in seq_along(draws_list)) {
    take <- sum(member == k)
    if (!take) next
    draws <- as.data.frame(draws_list[[k]], check.names = FALSE)
    out[[k]] <- draws[sample(seq_len(nrow(draws)), size = take, replace = TRUE), , drop = FALSE]
  }
  out <- do.call(rbind, out[!vapply(out, is.null, logical(1))])
  rownames(out) <- NULL
  out
}

local_atlas_draws_from_fit <- function(fit,
                                       population_model = fit$population_model,
                                       n_draws = nrow(fit$theta),
                                       seed = NULL) {
  model <- normalize_population_model(population_model)
  theta <- .as_hyper_matrix(fit$theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  w <- pmax(as.numeric(fit$w %||% rep(1 / nrow(theta), nrow(theta))), 0)
  sw <- sum(w)
  w <- if (!is.finite(sw) || sw <= 0) rep(1 / nrow(theta), nrow(theta)) else w / sw
  if (!is.null(seed)) {
    set.seed(as.integer(seed))
  }
  idx <- sample.int(nrow(theta), size = as.integer(n_draws), replace = TRUE, prob = w)
  theta <- theta[idx, , drop = FALSE]
  d <- model$alpha_dim
  mu <- as.data.frame(theta[, seq_len(d), drop = FALSE], check.names = FALSE)
  sigma2 <- as.data.frame(exp(theta[, d + seq_len(d), drop = FALSE]), check.names = FALSE)
  names(mu) <- paste0("mu_", model$alpha_names)
  names(sigma2) <- paste0("sigma2_", model$alpha_names)
  data.frame(mu, sigma2, check.names = FALSE)
}

local_atlas_graph_summary <- function(factor_set) {
  factor_set <- validate_local_atlas_factor_set(factor_set)
  rows <- lapply(names(factor_set$atlases), function(local_name) {
    atlas <- factor_set$atlases[[local_name]]
    active_charts <- .local_atlas_active_charts(atlas)
    active_edges <- atlas$edges[vapply(atlas$edges, function(edge) identical(edge$status, "active"), logical(1))]
    diag <- if (!is.null(atlas$normalizer_solution)) {
      atlas$normalizer_solution$diagnostics %||% atlas_cycle_diagnostics(atlas)$summary
    } else {
      list(
        max_abs_edge_residual = NA_real_,
        max_abs_standardized_edge_residual = NA_real_,
        n_direct_observations = NA_integer_,
        max_abs_direct_residual = NA_real_,
        max_abs_standardized_direct_residual = NA_real_,
        weighted_rss = NA_real_,
        cycle_degrees_freedom = NA_integer_
      )
    }
    root <- if (!is.null(atlas$root_chart_id) && atlas$root_chart_id %in% names(atlas$charts)) {
      atlas$charts[[atlas$root_chart_id]]
    } else {
      NULL
    }
    active_logZ_se <- vapply(active_charts, function(chart) as.numeric(chart$logZ_abs_se), numeric(1))
    data.frame(
      local = local_name,
      n_active_charts = length(active_charts),
      n_candidate_charts = sum(vapply(atlas$charts, function(chart) identical(chart$status, "candidate"), logical(1))),
      n_quarantined_charts = sum(vapply(atlas$charts, function(chart) identical(chart$status, "quarantined"), logical(1))),
      n_active_edges = length(active_edges),
      root_logZ_se = as.numeric(root$logZ_abs_se %||% NA_real_),
      max_chart_logZ_abs_se = if (length(active_logZ_se)) max(active_logZ_se, na.rm = TRUE) else NA_real_,
      mean_chart_logZ_abs_se = if (length(active_logZ_se)) mean(active_logZ_se, na.rm = TRUE) else NA_real_,
      max_abs_edge_residual = as.numeric(diag$max_abs_edge_residual %||% NA_real_),
      max_abs_standardized_edge_residual = as.numeric(diag$max_abs_standardized_edge_residual %||% NA_real_),
      n_direct_observations = as.integer(diag$n_direct_observations %||% NA_integer_),
      max_abs_direct_residual = as.numeric(diag$max_abs_direct_residual %||% NA_real_),
      max_abs_standardized_direct_residual = as.numeric(diag$max_abs_standardized_direct_residual %||% NA_real_),
      weighted_rss = as.numeric(diag$weighted_rss %||% NA_real_),
      cycle_degrees_freedom = as.integer(diag$cycle_degrees_freedom %||% NA_integer_),
      check.names = FALSE
    )
  })
  do.call(rbind, rows)
}

local_atlas_certification_summary <- function(factor_set,
                                              theta,
                                              theta_weights = NULL,
                                              n_cores = 1L) {
  factor_set <- validate_local_atlas_factor_set(factor_set)
  theta <- .as_hyper_matrix(
    theta,
    hyper_names = factor_set$population_model$hyper_names,
    hyper_dim = factor_set$population_model$hyper_dim
  )
  theta_weights <- .local_chart_normalize_weights(theta_weights, nrow(theta))
  parts <- .local_atlas_factor_set_by_local(factor_set, theta, n_cores = n_cores)
  rows <- lapply(names(parts), function(local_name) {
    part <- parts[[local_name]]
    certified <- part$status == "certified"
    data.frame(
      local = local_name,
      n_theta = nrow(part),
      certified_fraction = mean(certified),
      uncertified_fraction = mean(!certified),
      certified_weight = sum(theta_weights[certified]),
      uncertified_weight = sum(theta_weights[!certified]),
      max_se = suppressWarnings(max(part$se[certified], na.rm = TRUE)),
      median_se = suppressWarnings(stats::median(part$se[certified], na.rm = TRUE)),
      uncertified_reasons = paste(sort(unique(part$reason[!certified])), collapse = ","),
      check.names = FALSE
    )
  })
  bad_rows <- lapply(names(parts), function(local_name) {
    part <- parts[[local_name]]
    bad <- which(part$status != "certified")
    if (!length(bad)) return(data.frame())
    data.frame(
      local = local_name,
      theta_row = bad,
      reason = part$reason[bad],
      nearest_charts = part$nearest_charts[bad],
      check.names = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$max_se[!is.finite(out$max_se)] <- Inf
  out$median_se[!is.finite(out$median_se)] <- Inf
  attr(out, "bad_rows") <- if (length(bad_rows)) do.call(rbind, bad_rows) else data.frame()
  out
}

local_atlas_certification_summary_by_local_design <- function(factor_set,
                                                              local_theta_designs,
                                                              n_cores = 1L) {
  factor_set <- validate_local_atlas_factor_set(factor_set)
  if (!is.list(local_theta_designs) || !length(local_theta_designs)) {
    stop("local_theta_designs must be a non-empty list.")
  }
  atlas_names <- names(factor_set$atlases)
  if (is.null(names(local_theta_designs)) || any(!nzchar(names(local_theta_designs)))) {
    if (length(local_theta_designs) != length(atlas_names)) {
      stop("unnamed local_theta_designs must have one entry per atlas.")
    }
    names(local_theta_designs) <- atlas_names
  }
  missing_design <- setdiff(atlas_names, names(local_theta_designs))
  if (length(missing_design)) {
    stop("local_theta_designs is missing locals: ", paste(missing_design, collapse = ", "))
  }
  control <- factor_set$evaluator_control
  eval_one <- function(local_name) {
    theta <- .as_hyper_matrix(
      local_theta_designs[[local_name]],
      hyper_names = factor_set$population_model$hyper_names,
      hyper_dim = factor_set$population_model$hyper_dim
    )
    part <- evaluate_local_atlas_many(
      atlas = factor_set$atlases[[local_name]],
      theta = theta,
      population_model = factor_set$population_model,
      max_chart_distance = control$max_chart_distance,
      min_covering_charts = control$min_covering_charts,
      max_prediction_range = control$max_prediction_range,
      distance_scale = control$distance_scale,
      se_floor = control$se_floor,
      use_particle_mis = control$use_particle_mis,
      require_particle_mis = control$require_particle_mis,
      min_particle_mis_ess = control$min_particle_mis_ess,
      min_particle_mis_ess_abs = control$min_particle_mis_ess_abs,
      max_particle_mis_psis_k = control$max_particle_mis_psis_k,
      max_quadratic_particle_gap = control$max_quadratic_particle_gap,
      sparse_chart_min_covering = control$sparse_chart_min_covering,
      sparse_chart_max_distance = control$sparse_chart_max_distance,
      max_leave_chart_out_gap = control$max_leave_chart_out_gap,
      distance_metric = control$distance_metric,
      surface_method = control$surface_method,
      min_surface_charts = control$min_surface_charts,
      max_surface_se = control$max_surface_se,
      surface_value_nugget = control$surface_value_nugget,
      surface_gradient_weight = control$surface_gradient_weight,
      surface_curvature_weight = control$surface_curvature_weight,
      surface_ridge = control$surface_ridge,
      particle_mis_role = control$particle_mis_role
    )
    certified <- part$status == "certified"
    summary <- data.frame(
      local = local_name,
      n_theta = nrow(part),
      certified_fraction = mean(certified),
      uncertified_fraction = mean(!certified),
      certified_weight = mean(certified),
      uncertified_weight = mean(!certified),
      max_se = suppressWarnings(max(part$se[certified], na.rm = TRUE)),
      median_se = suppressWarnings(stats::median(part$se[certified], na.rm = TRUE)),
      uncertified_reasons = paste(sort(unique(part$reason[!certified])), collapse = ","),
      check.names = FALSE
    )
    summary$max_se[!is.finite(summary$max_se)] <- Inf
    summary$median_se[!is.finite(summary$median_se)] <- Inf
    bad <- if (any(!certified)) {
      data.frame(
        local = local_name,
        theta_row = which(!certified),
        reason = part$reason[!certified],
        nearest_charts = part$nearest_charts[!certified],
        check.names = FALSE
      )
    } else {
      data.frame()
    }
    list(summary = summary, bad = bad)
  }
  local_results <- if (as.integer(n_cores) <= 1L || length(atlas_names) <= 1L) {
    lapply(atlas_names, eval_one)
  } else {
    parallel::mclapply(atlas_names, eval_one, mc.cores = as.integer(min(n_cores, length(atlas_names))))
  }
  out <- do.call(rbind, lapply(local_results, `[[`, "summary"))
  bad <- do.call(rbind, lapply(local_results, `[[`, "bad"))
  if (is.null(bad)) {
    bad <- data.frame()
  }
  attr(out, "bad_rows") <- bad
  out
}

.local_atlas_design_repair_spec <- function(factor_set,
                                            design_certification,
                                            population_model,
                                            theta_design = NULL,
                                            local_theta_designs = NULL) {
  factor_set <- validate_local_atlas_factor_set(factor_set)
  model <- normalize_population_model(population_model)
  bad <- attr(design_certification, "bad_rows")
  if (is.null(bad) || !is.data.frame(bad) || !nrow(bad)) {
    return(list(theta = matrix(numeric(), nrow = 0L, ncol = model$hyper_dim,
                               dimnames = list(NULL, model$hyper_names)),
                pairs = data.frame(),
                bad_rows = data.frame()))
  }
  atlas_names <- names(factor_set$atlases)
  local_pos <- match(as.character(bad$local), atlas_names)
  theta_row <- as.integer(bad$theta_row)
  keep <- is.finite(local_pos) & is.finite(theta_row) & theta_row >= 1L
  if (!any(keep)) {
    return(list(theta = matrix(numeric(), nrow = 0L, ncol = model$hyper_dim,
                               dimnames = list(NULL, model$hyper_names)),
                pairs = data.frame(),
                bad_rows = data.frame()))
  }
  bad <- bad[keep, , drop = FALSE]
  local_pos <- local_pos[keep]
  theta_row <- theta_row[keep]

  local_design_names <- names(local_theta_designs)
  theta_rows <- vector("list", length(theta_row))
  for (k in seq_along(theta_row)) {
    if (!is.null(local_theta_designs)) {
      local_name <- atlas_names[local_pos[k]]
      design_k <- if (!is.null(local_design_names) && local_name %in% local_design_names) {
        local_theta_designs[[local_name]]
      } else {
        local_theta_designs[[local_pos[k]]]
      }
      design_k <- .as_hyper_matrix(design_k, model$hyper_names, model$hyper_dim)
      if (theta_row[k] > nrow(design_k)) next
      theta_rows[[k]] <- design_k[theta_row[k], , drop = FALSE]
    } else {
      design <- .as_hyper_matrix(theta_design, model$hyper_names, model$hyper_dim)
      if (theta_row[k] > nrow(design)) next
      theta_rows[[k]] <- design[theta_row[k], , drop = FALSE]
    }
  }
  valid <- vapply(theta_rows, Negate(is.null), logical(1))
  if (!any(valid)) {
    return(list(theta = matrix(numeric(), nrow = 0L, ncol = model$hyper_dim,
                               dimnames = list(NULL, model$hyper_names)),
                pairs = data.frame(),
                bad_rows = data.frame()))
  }
  theta <- do.call(rbind, theta_rows[valid])
  theta <- .as_hyper_matrix(theta, model$hyper_names, model$hyper_dim)
  bad <- bad[valid, , drop = FALSE]
  local_pos <- local_pos[valid]
  pairs <- data.frame(
    local_pos = as.integer(local_pos),
    theta_row = seq_len(nrow(theta)),
    check.names = FALSE
  )
  list(theta = theta, pairs = pairs, bad_rows = bad)
}

.local_atlas_prune_inactive <- function(atlas) {
  atlas <- validate_local_atlas(atlas)
  keep_charts <- vapply(atlas$charts, function(chart) identical(chart$status, "active"), logical(1))
  active_ids <- names(atlas$charts)[keep_charts]
  keep_edges <- vapply(atlas$edges, function(edge) {
    identical(edge$status, "active") &&
      edge$from_chart %in% active_ids &&
      edge$to_chart %in% active_ids
  }, logical(1))
  atlas$charts <- atlas$charts[keep_charts]
  atlas$edges <- atlas$edges[keep_edges]
  validate_local_atlas(atlas)
}

local_atlas_quarantine_invariance <- function(factor_set,
                                              theta,
                                              theta_weights = NULL,
                                              n_cores = 1L) {
  factor_set <- validate_local_atlas_factor_set(factor_set)
  theta <- .as_hyper_matrix(
    theta,
    hyper_names = factor_set$population_model$hyper_names,
    hyper_dim = factor_set$population_model$hyper_dim
  )
  theta_weights <- .local_chart_normalize_weights(theta_weights, nrow(theta))
  original <- local_evidence_atlas_factor_set_loglik(
    factor_set,
    theta = theta,
    include_constant = FALSE,
    n_cores = n_cores
  )
  pruned <- factor_set
  pruned$atlases <- lapply(pruned$atlases, .local_atlas_prune_inactive)
  pruned <- validate_local_atlas_factor_set(pruned)
  after <- local_evidence_atlas_factor_set_loglik(
    pruned,
    theta = theta,
    include_constant = FALSE,
    n_cores = n_cores
  )
  delta <- after - original
  data.frame(
    n_theta = length(delta),
    max_abs_delta = max(abs(delta)),
    mean_abs_delta = mean(abs(delta)),
    weighted_mean_abs_delta = sum(theta_weights * abs(delta)),
    check.names = FALSE
  )
}

.local_evidence_finite_values <- function(x) {
  x <- as.numeric(x)
  x[is.finite(x)]
}

.local_evidence_finite_sum <- function(x) {
  x <- .local_evidence_finite_values(x)
  if (length(x)) sum(x) else NA_real_
}

.local_evidence_finite_mean <- function(x) {
  x <- .local_evidence_finite_values(x)
  if (length(x)) mean(x) else NA_real_
}

.local_evidence_finite_rmse <- function(x) {
  x <- .local_evidence_finite_values(x)
  if (length(x)) sqrt(mean(x^2)) else NA_real_
}

.local_evidence_finite_max <- function(x) {
  x <- .local_evidence_finite_values(x)
  if (length(x)) max(x) else NA_real_
}

.local_evidence_finite_range <- function(x) {
  x <- .local_evidence_finite_values(x)
  if (length(x)) diff(range(x)) else NA_real_
}

.local_evidence_representative_rows <- function(weights, n) {
  weights <- .local_chart_normalize_weights(weights, length(weights))
  n <- min(as.integer(n), length(weights))
  if (n <= 0L) {
    return(integer())
  }
  top <- head(order(weights, decreasing = TRUE), min(5L, n))
  cdf <- cumsum(weights)
  probs <- seq(0.5 / n, 1 - 0.5 / n, length.out = n)
  systematic <- vapply(probs, function(p) which(cdf >= p)[1L], integer(1))
  unique(c(top, systematic))[seq_len(min(length(unique(c(top, systematic))), n))]
}

validate_local_evidence_theta_design <- function(design) {
  if (!inherits(design, "local_evidence_theta_design")) {
    stop("design must inherit from 'local_evidence_theta_design'.")
  }
  model <- normalize_population_model(design$population_model)
  design$theta <- .as_hyper_matrix(design$theta, model$hyper_names, model$hyper_dim)
  if (!is.data.frame(design$metadata) || nrow(design$metadata) != nrow(design$theta)) {
    stop("theta design metadata must be a data frame with one row per theta.")
  }
  required <- c("theta_row", "theta_id", "theta_label", "theta_weight", "theta_source")
  missing <- setdiff(required, names(design$metadata))
  if (length(missing)) {
    stop("theta design metadata is missing columns: ", paste(missing, collapse = ", "))
  }
  design$metadata$theta_row <- as.integer(design$metadata$theta_row)
  design$metadata$theta_weight <- .local_chart_normalize_weights(
    design$metadata$theta_weight,
    nrow(design$metadata)
  )
  design$population_model <- model
  design
}

build_local_evidence_calibration_design <- function(theta = NULL,
                                                    population_model = NULL,
                                                    fit = NULL,
                                                    weights = NULL,
                                                    factor_sets = NULL,
                                                    diagnostic_theta = NULL,
                                                    diagnostic_weights = NULL,
                                                    focus_hyper_names = NULL,
                                                    probs = c(0.05, 0.5, 0.95),
                                                    max_points = 9L,
                                                    include_center = TRUE,
                                                    inflation = 1,
                                                    diagnostic_max_points = 256L,
                                                    n_leverage_points = 1L,
                                                    n_uncertainty_points = 2L,
                                                    n_disagreement_points = 2L,
                                                    n_cores = 1L,
                                                    source = "workflow",
                                                    label_prefix = "theta") {
  if (!is.null(fit)) {
    theta <- theta %||% fit$theta
    weights <- weights %||% fit$w
    population_model <- population_model %||% fit$population_model
  }
  if (is.null(population_model)) {
    stop("population_model is required.")
  }
  model <- normalize_population_model(population_model)
  theta <- .as_hyper_matrix(theta, model$hyper_names, model$hyper_dim)
  weights <- .local_chart_normalize_weights(weights, nrow(theta))
  inflation <- as.numeric(inflation)
  inflation <- inflation[is.finite(inflation) & inflation > 0]
  if (!length(inflation)) {
    stop("inflation must contain at least one positive finite value.")
  }
  focus_hyper_names <- focus_hyper_names %||% model$hyper_names
  focus_idx <- match(focus_hyper_names, model$hyper_names)
  focus_idx <- focus_idx[is.finite(focus_idx)]
  if (!length(focus_idx)) {
    focus_idx <- seq_len(model$hyper_dim)
  }

  center <- vapply(seq_len(ncol(theta)), function(j) {
    .local_atlas_weighted_quantile(theta[, j], weights, 0.5)
  }, numeric(1))
  names(center) <- model$hyper_names

  candidates <- list()
  add_candidate <- function(theta_one,
                            row_id,
                            reason,
                            hyper = NA_character_,
                            prob = NA_real_,
                            inflation_value = 1,
                            source_weight = NA_real_,
                            score = 0,
                            mandatory = FALSE) {
    out <- .as_hyper_matrix(theta_one, model$hyper_names, model$hyper_dim)
    if (nrow(out) != 1L) {
      stop("internal calibration design error: candidate theta must have one row.")
    }
    candidates[[length(candidates) + 1L]] <<- list(
      theta = out,
      source_row = as.integer(row_id),
      reason = reason,
      source_hyper = hyper,
      source_prob = as.numeric(prob),
      source_weight = as.numeric(source_weight),
      inflation = as.numeric(inflation_value),
      design_score = as.numeric(score),
      mandatory = isTRUE(mandatory)
    )
  }

  if (isTRUE(include_center)) {
    center_row <- which.min(rowSums(sweep(theta, 2L, center, "-")^2))
    add_candidate(
      theta_one = theta[center_row, , drop = FALSE],
      row_id = center_row,
      reason = "weighted_center",
      source_weight = weights[center_row],
      score = Inf,
      mandatory = TRUE
    )
  }
  for (j in focus_idx) {
    qj <- .local_atlas_weighted_quantile(theta[, j], weights, probs)
    for (k in seq_along(qj)) {
      if (!is.finite(qj[k])) next
      row_id <- which.min(abs(theta[, j] - qj[k]))
      base_theta <- as.numeric(theta[row_id, ])
      tail_score <- abs(as.numeric(stats::qnorm(pmin(pmax(probs[k], 1e-6), 1 - 1e-6))))
      for (inflate in inflation) {
        inflated <- center + inflate * (base_theta - center)
        add_candidate(
          theta_one = matrix(inflated, nrow = 1L, dimnames = list(NULL, model$hyper_names)),
          row_id = row_id,
          reason = if (isTRUE(all.equal(inflate, 1))) "weighted_quantile" else "inflated_weighted_quantile",
          hyper = model$hyper_names[j],
          prob = probs[k],
          inflation_value = inflate,
          source_weight = weights[row_id],
          score = 1 + tail_score + log1p(weights[row_id] * nrow(theta))
        )
      }
    }
  }
  if (as.integer(n_leverage_points) > 0L) {
    leverage_rows <- head(order(weights, decreasing = TRUE), as.integer(n_leverage_points))
    for (row_id in leverage_rows) {
      add_candidate(
        theta_one = theta[row_id, , drop = FALSE],
        row_id = row_id,
        reason = "high_posterior_weight",
        source_weight = weights[row_id],
        score = 3 + log1p(weights[row_id] * nrow(theta))
      )
    }
  }

  diagnostic_rows <- NULL
  if (!is.null(factor_sets) &&
      (as.integer(n_uncertainty_points) > 0L || as.integer(n_disagreement_points) > 0L)) {
    factor_sets <- .local_evidence_factor_set_list(factor_sets)
    diagnostic_theta <- diagnostic_theta %||% theta
    diagnostic_theta <- .as_hyper_matrix(diagnostic_theta, model$hyper_names, model$hyper_dim)
    diagnostic_weights <- .local_chart_normalize_weights(diagnostic_weights, nrow(diagnostic_theta))
    if (nrow(diagnostic_theta) > as.integer(diagnostic_max_points)) {
      take <- .local_evidence_representative_rows(diagnostic_weights, diagnostic_max_points)
      diagnostic_theta <- diagnostic_theta[take, , drop = FALSE]
      diagnostic_weights <- diagnostic_weights[take]
    }
    local_names <- names(factor_sets[[1L]]$atlases)
    local_metadata <- data.frame(
      local_pos = seq_along(factor_sets[[1L]]$atlases),
      local = local_names %||% paste0("local_", seq_along(factor_sets[[1L]]$atlases)),
      check.names = FALSE
    )
    diagnostic_atlas_rows <- .local_evidence_atlas_rows(
      factor_sets = factor_sets,
      theta = diagnostic_theta,
      local_metadata = local_metadata,
      n_cores = n_cores
    )
    theta_parts <- split(diagnostic_atlas_rows, diagnostic_atlas_rows$theta_row)
    diagnostic_rows <- do.call(rbind, lapply(theta_parts, function(df) {
      member_total <- aggregate(
        df$atlas_log_marginal,
        by = list(member = df$member),
        FUN = .local_evidence_finite_sum
      )
      finite_total <- .local_evidence_finite_values(member_total$x)
      data.frame(
        theta_row = df$theta_row[1L],
        theta_weight = diagnostic_weights[df$theta_row[1L]],
        uncertified_locals = sum(df$atlas_status != "certified", na.rm = TRUE),
        total_surface_se = sqrt(sum(pmax(df$atlas_se[is.finite(df$atlas_se)], 0)^2)),
        member_loglik_sd = if (length(finite_total) > 1L) stats::sd(finite_total) else 0,
        check.names = FALSE
      )
    }))
    if (as.integer(n_uncertainty_points) > 0L && nrow(diagnostic_rows)) {
      uncertainty_score <- (
        100 * diagnostic_rows$uncertified_locals +
          diagnostic_rows$total_surface_se
      ) * pmax(diagnostic_rows$theta_weight, .Machine$double.eps)^0.25
      take <- head(order(uncertainty_score, decreasing = TRUE), as.integer(n_uncertainty_points))
      for (idx in take) {
        row_id <- diagnostic_rows$theta_row[idx]
        add_candidate(
          theta_one = diagnostic_theta[row_id, , drop = FALSE],
          row_id = NA_integer_,
          reason = "high_atlas_surface_uncertainty",
          source_weight = diagnostic_rows$theta_weight[idx],
          score = 4 + uncertainty_score[idx]
        )
      }
    }
    if (length(factor_sets) > 1L && as.integer(n_disagreement_points) > 0L && nrow(diagnostic_rows)) {
      disagreement_score <- diagnostic_rows$member_loglik_sd *
        pmax(diagnostic_rows$theta_weight, .Machine$double.eps)^0.25
      take <- head(order(disagreement_score, decreasing = TRUE), as.integer(n_disagreement_points))
      for (idx in take) {
        row_id <- diagnostic_rows$theta_row[idx]
        add_candidate(
          theta_one = diagnostic_theta[row_id, , drop = FALSE],
          row_id = NA_integer_,
          reason = "high_member_disagreement",
          source_weight = diagnostic_rows$theta_weight[idx],
          score = 4 + disagreement_score[idx]
        )
      }
    }
  }
  if (!length(candidates)) {
    stop("calibration design produced no theta points.")
  }

  theta_out <- do.call(rbind, lapply(candidates, `[[`, "theta"))
  key <- as.data.frame(signif(theta_out, digits = 10L), check.names = FALSE)
  keep <- !duplicated(key)
  theta_out <- theta_out[keep, , drop = FALSE]
  candidates <- candidates[keep]
  if (nrow(theta_out) > as.integer(max_points)) {
    mandatory <- vapply(candidates, `[[`, logical(1), "mandatory")
    score <- vapply(candidates, `[[`, numeric(1), "design_score")
    mandatory_idx <- which(mandatory)
    ranked_idx <- setdiff(order(score, decreasing = TRUE), mandatory_idx)
    take <- unique(c(mandatory_idx, ranked_idx))
    take <- take[seq_len(min(length(take), as.integer(max_points)))]
    theta_out <- theta_out[take, , drop = FALSE]
    candidates <- candidates[take]
  }

  raw_weight <- vapply(candidates, `[[`, numeric(1), "source_weight")
  raw_weight[!is.finite(raw_weight) | raw_weight <= 0] <- 1
  metadata <- data.frame(
    theta_row = seq_len(nrow(theta_out)),
    theta_id = sprintf("%s_%03d", label_prefix, seq_len(nrow(theta_out))),
    theta_label = sprintf("%s_%03d", label_prefix, seq_len(nrow(theta_out))),
    theta_weight = raw_weight,
    theta_source = as.character(source),
    source_row = vapply(candidates, `[[`, integer(1), "source_row"),
    source_reason = vapply(candidates, `[[`, character(1), "reason"),
    source_hyper = vapply(candidates, `[[`, character(1), "source_hyper"),
    source_prob = vapply(candidates, `[[`, numeric(1), "source_prob"),
    source_inflation = vapply(candidates, `[[`, numeric(1), "inflation"),
    design_score = vapply(candidates, `[[`, numeric(1), "design_score"),
    mandatory = vapply(candidates, `[[`, logical(1), "mandatory"),
    check.names = FALSE
  )
  metadata <- cbind(metadata, as.data.frame(theta_out, check.names = FALSE))
  metadata$theta_weight <- .local_chart_normalize_weights(metadata$theta_weight, nrow(metadata))

  structure(
    list(
      theta = theta_out,
      metadata = metadata,
      population_model = model,
      settings = list(
        focus_hyper_names = model$hyper_names[focus_idx],
        probs = as.numeric(probs),
        max_points = as.integer(max_points),
        include_center = isTRUE(include_center),
        inflation = as.numeric(inflation),
        diagnostic_max_points = as.integer(diagnostic_max_points),
        n_leverage_points = as.integer(n_leverage_points),
        n_uncertainty_points = as.integer(n_uncertainty_points),
        n_disagreement_points = as.integer(n_disagreement_points),
        source = as.character(source)
      )
    ),
    class = "local_evidence_theta_design"
  ) |>
    validate_local_evidence_theta_design()
}

.local_evidence_factor_set_list <- function(factor_sets) {
  if (inherits(factor_sets, "local_evidence_atlas_factor_set")) {
    factor_sets <- list(member_1 = factor_sets)
  }
  if (!is.list(factor_sets) || !length(factor_sets)) {
    stop("factor_sets must be a local atlas factor set or a non-empty list of them.")
  }
  factor_sets <- lapply(factor_sets, validate_local_atlas_factor_set)
  if (is.null(names(factor_sets)) || any(!nzchar(names(factor_sets)))) {
    names(factor_sets) <- paste0("member_", seq_along(factor_sets))
  }
  model <- factor_sets[[1L]]$population_model
  for (fs in factor_sets[-1L]) {
    if (!identical(fs$population_model$hyper_names, model$hyper_names) ||
        !identical(fs$population_model$alpha_names, model$alpha_names)) {
      stop("All factor sets must use the same population model.")
    }
  }
  factor_sets
}

.local_evidence_theta_inputs <- function(theta,
                                         population_model,
                                         theta_weights = NULL,
                                         theta_metadata = NULL,
                                         reference_source = "workflow") {
  if (inherits(theta, "local_evidence_theta_design") && is.null(population_model)) {
    population_model <- theta$population_model
  }
  model <- normalize_population_model(population_model)
  if (inherits(theta, "local_evidence_theta_design")) {
    design <- validate_local_evidence_theta_design(theta)
    theta <- design$theta
    theta_metadata <- design$metadata
  } else {
    theta <- .as_hyper_matrix(theta, model$hyper_names, model$hyper_dim)
  }
  if (is.null(theta_metadata)) {
    theta_metadata <- data.frame(
      theta_row = seq_len(nrow(theta)),
      theta_id = sprintf("theta_%03d", seq_len(nrow(theta))),
      theta_label = sprintf("theta_%03d", seq_len(nrow(theta))),
      theta_weight = .local_chart_normalize_weights(theta_weights, nrow(theta)),
      theta_source = as.character(reference_source),
      check.names = FALSE
    )
    theta_metadata <- cbind(theta_metadata, as.data.frame(theta, check.names = FALSE))
  } else {
    theta_metadata <- as.data.frame(theta_metadata, check.names = FALSE)
    if (!"theta_row" %in% names(theta_metadata)) {
      theta_metadata$theta_row <- seq_len(nrow(theta_metadata))
    }
    if (!"theta_id" %in% names(theta_metadata)) {
      theta_metadata$theta_id <- sprintf("theta_%03d", seq_len(nrow(theta_metadata)))
    }
    if (!"theta_label" %in% names(theta_metadata)) {
      theta_metadata$theta_label <- theta_metadata$theta_id
    }
    if (!"theta_weight" %in% names(theta_metadata)) {
      theta_metadata$theta_weight <- .local_chart_normalize_weights(theta_weights, nrow(theta))
    }
    if (!"theta_source" %in% names(theta_metadata)) {
      theta_metadata$theta_source <- as.character(reference_source)
    }
    for (name in model$hyper_names) {
      theta_metadata[[name]] <- theta[, name]
    }
  }
  theta_metadata$theta_row <- as.integer(theta_metadata$theta_row)
  theta_metadata$theta_weight <- .local_chart_normalize_weights(theta_metadata$theta_weight, nrow(theta))
  list(theta = theta, metadata = theta_metadata)
}

.local_evidence_resolve_locals <- function(factor_set, data_list, local_ids = seq_along(factor_set$atlases)) {
  factor_set <- validate_local_atlas_factor_set(factor_set)
  atlas_names <- names(factor_set$atlases)
  if (is.character(local_ids)) {
    local_ids <- match(local_ids, atlas_names)
  }
  local_ids <- as.integer(local_ids)
  local_ids <- local_ids[is.finite(local_ids)]
  if (!length(local_ids) || any(local_ids < 1L | local_ids > length(factor_set$atlases))) {
    stop("local_ids must identify atlases in factor_set.")
  }
  if (!is.list(data_list) || !length(data_list)) {
    stop("data_list must contain every selected local.")
  }
  local_names <- atlas_names[local_ids] %||% as.character(local_ids)
  if (!is.null(names(data_list)) && !is.null(atlas_names)) {
    data_pos <- match(local_names, names(data_list))
  } else {
    data_pos <- local_ids
  }
  if (any(!is.finite(data_pos)) || any(data_pos < 1L | data_pos > length(data_list))) {
    stop("data_list must contain every selected local, by name or matching position.")
  }
  data.frame(
    local_pos = local_ids,
    data_pos = as.integer(data_pos),
    local = local_names,
    check.names = FALSE
  )
}

.local_evidence_eval_atlas_one <- function(factor_set, local_pos, theta_one) {
  factor_set <- validate_local_atlas_factor_set(factor_set)
  model <- factor_set$population_model
  atlas <- factor_set$atlases[[local_pos]]
  control <- factor_set$evaluator_control
  can_batch <- isTRUE(control$use_particle_mis) &&
    isTRUE(control$particle_mis_batch) &&
    !is.finite(control$max_leave_chart_out_gap) &&
    identical(control$particle_mis_role, "estimator")
  if (can_batch) {
    out <- .local_atlas_particle_mis_many_global(
      atlas = atlas,
      theta = theta_one,
      population_model = model,
      max_chart_distance = control$max_chart_distance,
      min_covering_charts = control$min_covering_charts,
      min_ess_frac = control$min_particle_mis_ess,
      min_ess = control$min_particle_mis_ess_abs,
      max_psis_k = control$max_particle_mis_psis_k,
      sparse_chart_min_covering = control$sparse_chart_min_covering,
      sparse_chart_max_distance = control$sparse_chart_max_distance,
      distance_metric = control$distance_metric,
      se_floor = control$se_floor,
      use_uncertified_estimates = control$use_uncertified_estimates
    )
    return(out[1L, , drop = FALSE])
  }
  out <- evaluate_local_atlas(
    atlas = atlas,
    theta = theta_one,
    population_model = model,
    max_chart_distance = control$max_chart_distance,
    min_covering_charts = control$min_covering_charts,
    max_prediction_range = control$max_prediction_range,
    distance_scale = control$distance_scale,
    se_floor = control$se_floor,
    use_particle_mis = control$use_particle_mis,
    require_particle_mis = control$require_particle_mis,
    min_particle_mis_ess = control$min_particle_mis_ess,
    min_particle_mis_ess_abs = control$min_particle_mis_ess_abs,
    max_particle_mis_psis_k = control$max_particle_mis_psis_k,
    max_quadratic_particle_gap = control$max_quadratic_particle_gap,
    sparse_chart_min_covering = control$sparse_chart_min_covering,
    sparse_chart_max_distance = control$sparse_chart_max_distance,
    max_leave_chart_out_gap = control$max_leave_chart_out_gap,
    distance_metric = control$distance_metric,
    surface_method = control$surface_method,
    min_surface_charts = control$min_surface_charts,
    max_surface_se = control$max_surface_se,
    surface_value_nugget = control$surface_value_nugget,
    surface_gradient_weight = control$surface_gradient_weight,
    surface_curvature_weight = control$surface_curvature_weight,
    surface_ridge = control$surface_ridge,
    particle_mis_role = control$particle_mis_role
  )
  data.frame(
    log_marginal = out$log_marginal,
    se = out$se,
    status = out$status,
    reason = out$reason,
    nearest_charts = paste(out$nearest_charts, collapse = ","),
    particle_mis_ess_frac = as.numeric(out$diagnostics$particle_mis$ess_frac %||% NA_real_),
    particle_mis_psis_k = as.numeric(out$diagnostics$particle_mis$psis_k %||% NA_real_),
    surface_se = as.numeric(out$diagnostics$surface_se %||% NA_real_),
    surface_residual_sd = as.numeric(out$diagnostics$surface_residual_sd %||% NA_real_),
    leave_chart_out_gap = as.numeric(out$diagnostics$leave_chart_out_gap %||% NA_real_),
    min_covering_distance = as.numeric(out$diagnostics$min_covering_distance %||% NA_real_),
    check.names = FALSE
  )
}

.local_evidence_atlas_rows <- function(factor_sets, theta, local_metadata, n_cores = 1L) {
  jobs <- expand.grid(
    member = names(factor_sets),
    local_pos = local_metadata$local_pos,
    theta_row = seq_len(nrow(theta)),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  eval_job <- function(k) {
    member <- jobs$member[k]
    local_pos <- jobs$local_pos[k]
    theta_row <- jobs$theta_row[k]
    ev <- .local_evidence_eval_atlas_one(
      factor_sets[[member]],
      local_pos = local_pos,
      theta_one = theta[theta_row, , drop = FALSE]
    )
    data.frame(
      member = member,
      local_pos = local_pos,
      theta_row = theta_row,
      atlas_log_marginal = ev$log_marginal,
      atlas_se = ev$se,
      atlas_status = ev$status,
      atlas_reason = ev$reason,
      nearest_charts = ev$nearest_charts,
      particle_mis_ess_frac = ev$particle_mis_ess_frac %||% NA_real_,
      particle_mis_psis_k = ev$particle_mis_psis_k %||% NA_real_,
      surface_se = ev$surface_se %||% NA_real_,
      surface_residual_sd = ev$surface_residual_sd %||% NA_real_,
      leave_chart_out_gap = ev$leave_chart_out_gap %||% NA_real_,
      min_covering_distance = ev$min_covering_distance %||% NA_real_,
      check.names = FALSE
    )
  }
  rows <- if (as.integer(n_cores) <= 1L || nrow(jobs) <= 1L) {
    lapply(seq_len(nrow(jobs)), eval_job)
  } else {
    parallel::mclapply(
      seq_len(nrow(jobs)),
      eval_job,
      mc.cores = as.integer(min(n_cores, nrow(jobs)))
    )
  }
  failed <- vapply(rows, inherits, logical(1), "try-error")
  if (any(failed)) {
    stop(
      "atlas evidence evaluation failed for ",
      sum(failed),
      " job(s); first error: ",
      as.character(rows[[which(failed)[1L]]])
    )
  }
  do.call(rbind, rows)
}

run_local_evidence_replicates <- function(theta,
                                          data_list,
                                          loglik_fn,
                                          population_model,
                                          local_ids = seq_along(data_list),
                                          pairs = NULL,
                                          n_replicates = 1L,
                                          replicate_offset = 0L,
                                          M = 500L,
                                          local_control = list(),
                                          n_cores = 1L,
                                          seed = 123L,
                                          source = "local_evidence_audit",
                                          stop_on_error = TRUE,
                                          verbose = FALSE) {
  model <- normalize_population_model(population_model)
  theta <- .as_hyper_matrix(theta, model$hyper_names, model$hyper_dim)
  if (!is.function(loglik_fn)) {
    stop("loglik_fn must be a function.")
  }
  if (is.null(pairs)) {
    if (is.character(local_ids)) {
      if (is.null(names(data_list))) {
        stop("character local_ids require named data_list.")
      }
      local_ids <- match(local_ids, names(data_list))
    }
    local_ids <- as.integer(local_ids)
    if (!length(local_ids) || any(!is.finite(local_ids)) ||
        any(local_ids < 1L | local_ids > length(data_list))) {
      stop("local_ids must identify entries in data_list.")
    }
    pairs <- expand.grid(
      local_pos = local_ids,
      theta_row = seq_len(nrow(theta)),
      KEEP.OUT.ATTRS = FALSE
    )
    pairs$data_pos <- pairs$local_pos
    pairs$local <- names(data_list)[pairs$data_pos] %||% as.character(pairs$local_pos)
  } else {
    if (!all(c("local_pos", "theta_row") %in% names(pairs))) {
      stop("pairs must contain local_pos and theta_row.")
    }
    keep_cols <- intersect(c("local_pos", "data_pos", "theta_row", "local"), names(pairs))
    pairs <- unique(as.data.frame(pairs[, keep_cols, drop = FALSE]))
    pairs$local_pos <- as.integer(pairs$local_pos)
    if (!"data_pos" %in% names(pairs)) {
      pairs$data_pos <- pairs$local_pos
    }
    pairs$data_pos <- as.integer(pairs$data_pos)
    pairs$theta_row <- as.integer(pairs$theta_row)
    if (!"local" %in% names(pairs)) {
      pairs$local <- names(data_list)[pairs$data_pos] %||% as.character(pairs$local_pos)
    }
    if (any(!is.finite(pairs$local_pos)) || any(!is.finite(pairs$data_pos)) ||
        any(!is.finite(pairs$theta_row)) ||
        any(pairs$data_pos < 1L | pairs$data_pos > length(data_list)) ||
        any(pairs$theta_row < 1L | pairs$theta_row > nrow(theta))) {
      stop("pairs contain invalid data_pos or theta_row.")
    }
  }
  n_replicates <- as.integer(n_replicates)
  if (n_replicates < 0L) {
    stop("n_replicates must be non-negative.")
  }
  if (n_replicates == 0L || !nrow(pairs)) {
    return(data.frame())
  }

  local_control <- .local_atlas_merge_control(local_control, .local_atlas_default_local_control())
  jobs <- merge(
    pairs,
    data.frame(replicate_index = seq_len(n_replicates), check.names = FALSE),
    by = NULL
  )
  jobs$replicate_id <- as.integer(replicate_offset) + jobs$replicate_index

  run_job <- function(k) {
    local_pos <- jobs$local_pos[k]
    data_pos <- jobs$data_pos[k]
    local <- as.character(jobs$local[k])
    theta_row <- jobs$theta_row[k]
    replicate_id <- jobs$replicate_id[k]
    run_seed <- as.integer(seed) + 1000003L * as.integer(replicate_id) +
      1009L * as.integer(local_pos) + 9176L * as.integer(theta_row)
    result <- tryCatch(
      .local_chart_run_smc(
        local_id = local,
        theta_anchor = theta[theta_row, , drop = FALSE],
        data_i = data_list[[data_pos]],
        loglik_fn = loglik_fn,
        population_model = model,
        M = as.integer(M),
        target_cess = local_control$target_cess,
        resample_threshold = local_control$resample_threshold,
        n_mcmc_moves = as.integer(local_control$n_mcmc_moves),
        rw_scale = local_control$rw_scale,
        G_mix = as.integer(local_control$G_mix),
        da_enable = isTRUE(local_control$da_enable),
        refit_every = as.integer(local_control$refit_every),
        max_steps = as.integer(local_control$max_steps),
        deterministic_resampling = isTRUE(local_control$deterministic_resampling),
        n_cores = 1L,
        seed = run_seed,
        verbose = isTRUE(verbose),
        source = source
      ),
      error = function(e) e
    )
    if (inherits(result, "error")) {
      if (isTRUE(stop_on_error)) {
        stop(
          "local evidence replicate failed: local_pos=", local_pos,
          " theta_row=", theta_row,
          " replicate_id=", replicate_id,
          ": ", conditionMessage(result)
        )
      }
      return(data.frame(
        local = local,
        local_pos = local_pos,
        data_pos = data_pos,
        theta_row = theta_row,
        replicate_id = replicate_id,
        M = as.integer(M),
        log_marginal = NA_real_,
        path_se = NA_real_,
        final_ess_frac = NA_real_,
        min_path_ess_frac = NA_real_,
        mean_accept_rate = NA_real_,
        rounds = NA_integer_,
        seed = run_seed,
        status = "error",
        error_message = conditionMessage(result),
        check.names = FALSE
      ))
    }
    diag <- result$diagnostics %||% list()
    data.frame(
      local = local,
      local_pos = local_pos,
      data_pos = data_pos,
      theta_row = theta_row,
      replicate_id = replicate_id,
      M = as.integer(M),
      log_marginal = result$logZ,
      path_se = result$logZ_se,
      final_ess_frac = as.numeric(diag$final_ess_frac %||% NA_real_),
      min_path_ess_frac = as.numeric(diag$min_path_ess_frac %||% NA_real_),
      mean_accept_rate = as.numeric(diag$mean_accept_rate %||% NA_real_),
      rounds = as.integer(diag$rounds %||% NA_integer_),
      seed = run_seed,
      status = "ok",
      error_message = NA_character_,
      check.names = FALSE
    )
  }

  rows <- if (as.integer(n_cores) <= 1L || nrow(jobs) <= 1L) {
    lapply(seq_len(nrow(jobs)), run_job)
  } else {
    parallel::mclapply(
      seq_len(nrow(jobs)),
      run_job,
      mc.cores = as.integer(min(n_cores, nrow(jobs)))
    )
  }
  failed <- vapply(rows, inherits, logical(1), "try-error")
  if (any(failed)) {
    stop(
      "local evidence replicate failed for ",
      sum(failed),
      " job(s); first error: ",
      as.character(rows[[which(failed)[1L]]])
    )
  }
  do.call(rbind, rows)
}

.local_evidence_replicate_summary <- function(replicates, uncertainty_floor = 1e-6) {
  if (is.null(replicates) || !nrow(replicates)) {
    return(data.frame())
  }
  parts <- split(replicates, interaction(replicates$local_pos, replicates$theta_row, drop = TRUE))
  rows <- lapply(parts, function(df) {
    ok <- df$status == "ok" & is.finite(df$log_marginal)
    log_m <- df$log_marginal[ok]
    path_se <- df$path_se[ok]
    n_ok <- length(log_m)
    empirical_sd <- if (n_ok > 1L) stats::sd(log_m) else NA_real_
    path_component <- if (n_ok) sqrt(mean(pmax(path_se, 0)^2, na.rm = TRUE) / n_ok) else NA_real_
    uncertainty <- sqrt(max(
      if (is.finite(empirical_sd)) empirical_sd^2 else -Inf,
      if (is.finite(path_component)) path_component^2 else -Inf,
      as.numeric(uncertainty_floor)^2,
      na.rm = TRUE
    ))
    data.frame(
      local = df$local[1L],
      local_pos = df$local_pos[1L],
      theta_row = df$theta_row[1L],
      n_replicates = nrow(df),
      n_ok = n_ok,
      fresh_log_m_center = if (n_ok) stats::median(log_m) else NA_real_,
      fresh_log_m_mean = if (n_ok) mean(log_m) else NA_real_,
      fresh_log_m_sd = empirical_sd,
      fresh_log_m_path_se = path_component,
      fresh_log_m_uncertainty = uncertainty,
      fresh_min_log_m = if (n_ok) min(log_m) else NA_real_,
      fresh_max_log_m = if (n_ok) max(log_m) else NA_real_,
      check.names = FALSE
    )
  })
  do.call(rbind, rows)
}

.local_evidence_join_rows <- function(atlas_rows,
                                      replicates,
                                      theta_metadata,
                                      local_metadata,
                                      reference_source = "workflow",
                                      uncertainty_floor = 1e-6) {
  summary <- .local_evidence_replicate_summary(replicates, uncertainty_floor = uncertainty_floor)
  rows <- merge(atlas_rows, local_metadata, by = "local_pos", all.x = TRUE)
  rows <- merge(rows, theta_metadata, by = "theta_row", all.x = TRUE)
  if (nrow(summary)) {
    rows <- merge(rows, summary, by = c("local_pos", "theta_row", "local"), all.x = TRUE)
  } else {
    rows$n_replicates <- 0L
    rows$n_ok <- 0L
    rows$fresh_log_m_center <- NA_real_
    rows$fresh_log_m_mean <- NA_real_
    rows$fresh_log_m_sd <- NA_real_
    rows$fresh_log_m_path_se <- NA_real_
    rows$fresh_log_m_uncertainty <- NA_real_
    rows$fresh_min_log_m <- NA_real_
    rows$fresh_max_log_m <- NA_real_
  }
  rows$reference_source <- as.character(reference_source)
  rows$error_fresh_minus_atlas <- rows$fresh_log_m_center - rows$atlas_log_marginal
  rows$abs_error_fresh_minus_atlas <- abs(rows$error_fresh_minus_atlas)
  denom <- sqrt(pmax(rows$atlas_se, 0)^2 + pmax(rows$fresh_log_m_uncertainty, 0)^2)
  rows$standardized_error <- rows$error_fresh_minus_atlas / pmax(denom, .Machine$double.eps)
  rows$abs_standardized_error <- abs(rows$standardized_error)
  rows[order(rows$member, rows$local_pos, rows$theta_row), , drop = FALSE]
}

new_local_evidence_audit <- function(rows,
                                     atlas_rows,
                                     replicates,
                                     theta,
                                     theta_metadata,
                                     local_metadata,
                                     population_model,
                                     factor_set_names,
                                     settings = list()) {
  structure(
    list(
      rows = rows,
      atlas_rows = atlas_rows,
      replicates = replicates,
      theta = theta,
      theta_metadata = theta_metadata,
      local_metadata = local_metadata,
      population_model = population_model,
      factor_set_names = factor_set_names,
      settings = settings
    ),
    class = "local_evidence_audit"
  ) |>
    validate_local_evidence_audit()
}

validate_local_evidence_audit <- function(audit) {
  if (!inherits(audit, "local_evidence_audit")) {
    stop("audit must inherit from 'local_evidence_audit'.")
  }
  model <- normalize_population_model(audit$population_model)
  audit$theta <- .as_hyper_matrix(audit$theta, model$hyper_names, model$hyper_dim)
  for (field in c("rows", "atlas_rows", "replicates", "theta_metadata", "local_metadata")) {
    if (!is.data.frame(audit[[field]])) {
      stop(field, " must be a data frame.")
    }
  }
  row_required <- c(
    "member", "local", "local_pos", "data_pos", "theta_row", "theta_id", "theta_weight",
    "atlas_log_marginal", "atlas_se", "atlas_status", "fresh_log_m_center",
    "fresh_log_m_uncertainty", "error_fresh_minus_atlas"
  )
  missing <- setdiff(row_required, names(audit$rows))
  if (length(missing)) {
    stop("audit rows are missing columns: ", paste(missing, collapse = ", "))
  }
  local_required <- c("local", "local_pos", "data_pos")
  missing <- setdiff(local_required, names(audit$local_metadata))
  if (length(missing)) {
    stop("audit local_metadata is missing columns: ", paste(missing, collapse = ", "))
  }
  audit$population_model <- model
  audit
}

build_local_evidence_audit <- function(factor_sets,
                                       theta,
                                       data_list,
                                       loglik_fn,
                                       local_ids = NULL,
                                       theta_weights = NULL,
                                       theta_metadata = NULL,
                                       n_replicates = 1L,
                                       M = 500L,
                                       local_control = list(),
                                       uncertainty_floor = 1e-6,
                                       reference_source = "workflow",
                                       n_cores = 1L,
                                       seed = 123L,
                                       stop_on_error = TRUE,
                                       verbose = FALSE) {
  factor_sets <- .local_evidence_factor_set_list(factor_sets)
  model <- factor_sets[[1L]]$population_model
  theta_input <- .local_evidence_theta_inputs(
    theta,
    population_model = model,
    theta_weights = theta_weights,
    theta_metadata = theta_metadata,
    reference_source = reference_source
  )
  theta <- theta_input$theta
  theta_metadata <- theta_input$metadata
  local_metadata <- .local_evidence_resolve_locals(
    factor_sets[[1L]],
    data_list = data_list,
    local_ids = local_ids %||% seq_along(factor_sets[[1L]]$atlases)
  )
  atlas_rows <- .local_evidence_atlas_rows(
    factor_sets = factor_sets,
    theta = theta,
    local_metadata = local_metadata,
    n_cores = n_cores
  )
  replicate_pairs <- merge(
    expand.grid(
      local_pos = local_metadata$local_pos,
      theta_row = seq_len(nrow(theta)),
      KEEP.OUT.ATTRS = FALSE
    ),
    local_metadata[, c("local_pos", "data_pos", "local"), drop = FALSE],
    by = "local_pos",
    all.x = TRUE,
    sort = FALSE
  )
  replicates <- run_local_evidence_replicates(
    theta = theta,
    data_list = data_list,
    loglik_fn = loglik_fn,
    population_model = model,
    pairs = replicate_pairs,
    n_replicates = n_replicates,
    M = M,
    local_control = local_control,
    n_cores = n_cores,
    seed = seed,
    source = "local_evidence_audit",
    stop_on_error = stop_on_error,
    verbose = verbose
  )
  rows <- .local_evidence_join_rows(
    atlas_rows = atlas_rows,
    replicates = replicates,
    theta_metadata = theta_metadata,
    local_metadata = local_metadata,
    reference_source = reference_source,
    uncertainty_floor = uncertainty_floor
  )
  new_local_evidence_audit(
    rows = rows,
    atlas_rows = atlas_rows,
    replicates = replicates,
    theta = theta,
    theta_metadata = theta_metadata,
    local_metadata = local_metadata,
    population_model = model,
    factor_set_names = names(factor_sets),
    settings = list(
      n_replicates = as.integer(n_replicates),
      M = as.integer(M),
      local_control = local_control,
      uncertainty_floor = as.numeric(uncertainty_floor),
      reference_source = as.character(reference_source),
      seed = as.integer(seed)
    )
  )
}

extend_local_evidence_audit <- function(audit,
                                        data_list,
                                        loglik_fn,
                                        pairs,
                                        n_replicates = 2L,
                                        M = audit$settings$M %||% 500L,
                                        local_control = audit$settings$local_control %||% list(),
                                        uncertainty_floor = audit$settings$uncertainty_floor %||% 1e-6,
                                        n_cores = 1L,
                                        seed = 123L,
                                        stop_on_error = TRUE,
                                        verbose = FALSE) {
  audit <- validate_local_evidence_audit(audit)
  if (is.null(pairs) || !nrow(pairs)) {
    return(audit)
  }
  pairs <- unique(as.data.frame(pairs[, c("local_pos", "theta_row"), drop = FALSE]))
  pairs <- merge(
    pairs,
    audit$local_metadata[, c("local_pos", "data_pos", "local"), drop = FALSE],
    by = "local_pos",
    all.x = TRUE,
    sort = FALSE
  )
  replicate_offset <- if (nrow(audit$replicates)) {
    max(as.integer(audit$replicates$replicate_id), na.rm = TRUE)
  } else {
    0L
  }
  new_reps <- run_local_evidence_replicates(
    theta = audit$theta,
    data_list = data_list,
    loglik_fn = loglik_fn,
    population_model = audit$population_model,
    pairs = pairs,
    n_replicates = n_replicates,
    replicate_offset = replicate_offset,
    M = M,
    local_control = local_control,
    n_cores = n_cores,
    seed = seed,
    source = "local_evidence_audit_extension",
    stop_on_error = stop_on_error,
    verbose = verbose
  )
  audit$replicates <- rbind(audit$replicates, new_reps)
  audit$rows <- .local_evidence_join_rows(
    atlas_rows = audit$atlas_rows,
    replicates = audit$replicates,
    theta_metadata = audit$theta_metadata,
    local_metadata = audit$local_metadata,
    reference_source = audit$settings$reference_source %||% "workflow",
    uncertainty_floor = uncertainty_floor
  )
  audit$settings$extended <- TRUE
  validate_local_evidence_audit(audit)
}

summarize_local_evidence_audit <- function(audit) {
  audit <- validate_local_evidence_audit(audit)
  rows <- audit$rows
  finite_error <- is.finite(rows$error_fresh_minus_atlas)

  member_summary <- do.call(rbind, lapply(split(rows, rows$member), function(df) {
    err <- df$error_fresh_minus_atlas[is.finite(df$error_fresh_minus_atlas)]
    theta_total <- aggregate(
      df$error_fresh_minus_atlas,
      by = list(theta_row = df$theta_row),
      FUN = .local_evidence_finite_sum
    )
    centered <- theta_total$x - .local_evidence_finite_mean(theta_total$x)
    data.frame(
      member = df$member[1L],
      n_rows = nrow(df),
      n_uncertified = sum(df$atlas_status != "certified", na.rm = TRUE),
      local_error_rmse = if (length(err)) sqrt(mean(err^2)) else NA_real_,
      local_error_mae = if (length(err)) mean(abs(err)) else NA_real_,
      max_abs_local_error = if (length(err)) max(abs(err)) else NA_real_,
      mean_fresh_uncertainty = .local_evidence_finite_mean(df$fresh_log_m_uncertainty),
      centered_total_error_range = .local_evidence_finite_range(centered),
      check.names = FALSE
    )
  }))

  theta_summary <- do.call(rbind, lapply(split(rows, list(rows$member, rows$theta_row), drop = TRUE), function(df) {
    err <- df$error_fresh_minus_atlas[is.finite(df$error_fresh_minus_atlas)]
    data.frame(
      member = df$member[1L],
      theta_row = df$theta_row[1L],
      theta_id = df$theta_id[1L],
      theta_label = df$theta_label[1L],
      theta_weight = df$theta_weight[1L],
      total_log_surface_error = if (length(err)) sum(err) else NA_real_,
      mean_abs_local_error = if (length(err)) mean(abs(err)) else NA_real_,
      rmse_local_error = if (length(err)) sqrt(mean(err^2)) else NA_real_,
      n_uncertified_locals = sum(df$atlas_status != "certified", na.rm = TRUE),
      check.names = FALSE
    )
  }))
  theta_summary <- do.call(rbind, lapply(split(theta_summary, theta_summary$member), function(df) {
    df$centered_total_log_surface_error <- df$total_log_surface_error -
      .local_evidence_finite_mean(df$total_log_surface_error)
    df
  }))
  theta_summary <- merge(
    theta_summary,
    audit$theta_metadata,
    by = c("theta_row", "theta_id", "theta_label", "theta_weight"),
    all.x = TRUE
  )

  local_summary <- do.call(rbind, lapply(split(rows, list(rows$member, rows$local_pos), drop = TRUE), function(df) {
    err <- df$error_fresh_minus_atlas[is.finite(df$error_fresh_minus_atlas)]
    data.frame(
      member = df$member[1L],
      local = df$local[1L],
      local_pos = df$local_pos[1L],
      rmse_local_error = if (length(err)) sqrt(mean(err^2)) else NA_real_,
      mean_signed_error = if (length(err)) mean(err) else NA_real_,
      max_abs_error = if (length(err)) max(abs(err)) else NA_real_,
      n_uncertified = sum(df$atlas_status != "certified", na.rm = TRUE),
      mean_fresh_uncertainty = .local_evidence_finite_mean(df$fresh_log_m_uncertainty),
      check.names = FALSE
    )
  }))
  local_summary <- local_summary[order(local_summary$member, -local_summary$rmse_local_error, na.last = TRUE), , drop = FALSE]

  failure_summary <- data.frame(
    n_rows = nrow(rows),
    n_finite_error = sum(finite_error),
    n_uncertified = sum(rows$atlas_status != "certified", na.rm = TRUE),
    n_abs_error_gt_1 = sum(abs(rows$error_fresh_minus_atlas) > 1, na.rm = TRUE),
    n_abs_z_gt_3 = sum(rows$abs_standardized_error > 3, na.rm = TRUE),
    max_abs_error = .local_evidence_finite_max(abs(rows$error_fresh_minus_atlas)),
    max_abs_z = .local_evidence_finite_max(rows$abs_standardized_error),
    check.names = FALSE
  )

  structure(
    list(
      member_summary = member_summary,
      theta_summary = theta_summary,
      local_summary = local_summary,
      failure_summary = failure_summary
    ),
    class = "local_evidence_audit_summary"
  )
}

select_audit_failures <- function(audit,
                                  member = NULL,
                                  max_pairs = 24L,
                                  min_abs_error = 1,
                                  min_abs_z = 3,
                                  include_uncertified = TRUE,
                                  weight_power = 0.5,
                                  uncertainty_weight = 0.25,
                                  uncertified_bonus = 10) {
  audit <- validate_local_evidence_audit(audit)
  rows <- audit$rows
  if (!is.null(member)) {
    rows <- rows[rows$member %in% member, , drop = FALSE]
  }
  if (!nrow(rows)) {
    return(data.frame())
  }
  bad_uncertified <- !is.na(rows$atlas_status) & rows$atlas_status != "certified"
  bad_error <- is.finite(rows$error_fresh_minus_atlas) &
    abs(rows$error_fresh_minus_atlas) >= as.numeric(min_abs_error)
  bad_z <- is.finite(rows$abs_standardized_error) &
    rows$abs_standardized_error >= as.numeric(min_abs_z)
  eligible <- bad_error | bad_z | (isTRUE(include_uncertified) & bad_uncertified)
  rows <- rows[eligible, , drop = FALSE]
  if (!nrow(rows)) {
    return(data.frame())
  }
  rows$selection_reason <- ifelse(
    bad_uncertified[eligible],
    "uncertified_atlas",
    ifelse(bad_error[eligible] & bad_z[eligible], "large_error_and_z",
           ifelse(bad_error[eligible], "large_error", "large_z"))
  )
  rows$selection_score <- (
    pmax(abs(rows$error_fresh_minus_atlas), 0, na.rm = TRUE) +
      as.numeric(uncertainty_weight) * pmax(rows$fresh_log_m_uncertainty, 0, na.rm = TRUE) +
      ifelse(rows$atlas_status != "certified", as.numeric(uncertified_bonus), 0)
  ) * pmax(rows$theta_weight, .Machine$double.eps)^as.numeric(weight_power)

  parts <- split(rows, interaction(rows$local_pos, rows$theta_row, drop = TRUE))
  selected <- lapply(parts, function(df) {
    best <- df[which.max(df$selection_score), , drop = FALSE]
    data.frame(
      local = best$local,
      local_pos = best$local_pos,
      theta_row = best$theta_row,
      theta_id = best$theta_id,
      theta_label = best$theta_label,
      theta_weight = best$theta_weight,
      score = best$selection_score,
      reason = paste(unique(df$selection_reason), collapse = ","),
      max_abs_error = .local_evidence_finite_max(abs(df$error_fresh_minus_atlas)),
      max_abs_z = .local_evidence_finite_max(df$abs_standardized_error),
      max_fresh_uncertainty = .local_evidence_finite_max(df$fresh_log_m_uncertainty),
      members = paste(unique(df$member), collapse = ","),
      check.names = FALSE
    )
  })
  out <- do.call(rbind, selected)
  out <- out[order(-out$score, -out$max_abs_error), , drop = FALSE]
  utils::head(out, as.integer(max_pairs))
}

.local_evidence_matrix_key <- function(member, local) {
  paste(as.character(member), as.character(local), sep = "\r")
}

.local_evidence_solve <- function(A, B) {
  A <- as.matrix(A)
  B <- as.matrix(B)
  chol_A <- tryCatch(chol(A), error = function(e) NULL)
  if (!is.null(chol_A)) {
    return(backsolve(chol_A, forwardsolve(t(chol_A), B)))
  }
  solve_A <- tryCatch(solve(A, B), error = function(e) NULL)
  if (!is.null(solve_A)) {
    return(solve_A)
  }
  qr.solve(A, B)
}

.local_evidence_default_active_hyper_names <- function(audit,
                                                       max_active_hyper = 6L) {
  audit <- validate_local_evidence_audit(audit)
  model <- audit$population_model
  source_hyper <- character()
  if ("source_hyper" %in% names(audit$theta_metadata)) {
    source_hyper <- as.character(audit$theta_metadata$source_hyper)
    source_hyper <- source_hyper[nzchar(source_hyper) & !is.na(source_hyper)]
  }
  source_hyper <- intersect(unique(source_hyper), model$hyper_names)
  if (length(source_hyper)) {
    return(head(source_hyper, as.integer(max_active_hyper)))
  }
  theta <- .as_hyper_matrix(audit$theta, model$hyper_names, model$hyper_dim)
  w <- .local_chart_normalize_weights(audit$theta_metadata$theta_weight, nrow(theta))
  center <- colSums(theta * w)
  spread <- sqrt(colSums(sweep(theta, 2L, center, "-")^2 * w))
  ranked <- order(spread, decreasing = TRUE)
  model$hyper_names[head(ranked, min(as.integer(max_active_hyper), length(ranked)))]
}

.local_evidence_basis_spec <- function(theta,
                                       population_model,
                                       theta_weights,
                                       active_hyper_names,
                                       include_quadratic = TRUE,
                                       include_interactions = FALSE,
                                       scale_floor = 1e-8) {
  model <- normalize_population_model(population_model)
  theta <- .as_hyper_matrix(theta, model$hyper_names, model$hyper_dim)
  theta_weights <- .local_chart_normalize_weights(theta_weights, nrow(theta))
  active_hyper_names <- intersect(unique(as.character(active_hyper_names)), model$hyper_names)
  if (!length(active_hyper_names)) {
    stop("active_hyper_names must identify at least one population hyperparameter.")
  }
  center <- colSums(theta * theta_weights)
  names(center) <- model$hyper_names
  scale <- sqrt(colSums(sweep(theta, 2L, center, "-")^2 * theta_weights))
  names(scale) <- model$hyper_names
  scale[!is.finite(scale) | scale < as.numeric(scale_floor)] <- 1

  terms <- data.frame(
    term = "intercept",
    kind = "intercept",
    hyper_1 = NA_character_,
    hyper_2 = NA_character_,
    check.names = FALSE
  )
  for (name in active_hyper_names) {
    terms <- rbind(
      terms,
      data.frame(term = paste0("linear:", name), kind = "linear",
                 hyper_1 = name, hyper_2 = NA_character_, check.names = FALSE)
    )
  }
  if (isTRUE(include_quadratic)) {
    for (name in active_hyper_names) {
      terms <- rbind(
        terms,
        data.frame(term = paste0("quadratic:", name), kind = "quadratic",
                   hyper_1 = name, hyper_2 = NA_character_, check.names = FALSE)
      )
    }
  }
  if (isTRUE(include_interactions) && length(active_hyper_names) > 1L) {
    for (pair in utils::combn(active_hyper_names, 2L, simplify = FALSE)) {
      terms <- rbind(
        terms,
        data.frame(term = paste0("interaction:", pair[1L], ":", pair[2L]),
                   kind = "interaction", hyper_1 = pair[1L], hyper_2 = pair[2L],
                   check.names = FALSE)
      )
    }
  }

  structure(
    list(
      population_model = model,
      active_hyper_names = active_hyper_names,
      center = center,
      scale = scale,
      terms = terms,
      include_quadratic = isTRUE(include_quadratic),
      include_interactions = isTRUE(include_interactions)
    ),
    class = "local_evidence_error_basis"
  )
}

.local_evidence_basis_matrix <- function(theta, basis) {
  if (!inherits(basis, "local_evidence_error_basis")) {
    stop("basis must inherit from 'local_evidence_error_basis'.")
  }
  model <- normalize_population_model(basis$population_model)
  theta <- .as_hyper_matrix(theta, model$hyper_names, model$hyper_dim)
  z <- sweep(theta, 2L, basis$center, "-")
  z <- sweep(z, 2L, basis$scale, "/")
  terms <- basis$terms
  X <- matrix(1, nrow = nrow(theta), ncol = nrow(terms))
  colnames(X) <- terms$term
  for (k in seq_len(nrow(terms))) {
    kind <- terms$kind[k]
    if (identical(kind, "intercept")) {
      X[, k] <- 1
    } else if (identical(kind, "linear")) {
      X[, k] <- z[, terms$hyper_1[k]]
    } else if (identical(kind, "quadratic")) {
      X[, k] <- z[, terms$hyper_1[k]]^2 - 1
    } else if (identical(kind, "interaction")) {
      X[, k] <- z[, terms$hyper_1[k]] * z[, terms$hyper_2[k]]
    } else {
      stop("unknown local evidence basis term kind: ", kind)
    }
  }
  X
}

.local_evidence_basis_penalty <- function(basis,
                                          intercept_penalty = 0.25,
                                          linear_penalty = 1,
                                          quadratic_penalty = 2,
                                          interaction_penalty = 4) {
  terms <- basis$terms
  penalty <- numeric(nrow(terms))
  penalty[terms$kind == "intercept"] <- as.numeric(intercept_penalty)
  penalty[terms$kind == "linear"] <- as.numeric(linear_penalty)
  penalty[terms$kind == "quadratic"] <- as.numeric(quadratic_penalty)
  penalty[terms$kind == "interaction"] <- as.numeric(interaction_penalty)
  penalty[!is.finite(penalty) | penalty < 0] <- 0
  penalty
}

.local_evidence_fit_error_one <- function(df,
                                          X_all,
                                          ridge_lambda,
                                          penalty,
                                          residual_floor,
                                          no_data_sd) {
  p <- ncol(X_all)
  zero_cov <- diag(as.numeric(no_data_sd)^2, p)
  dimnames(zero_cov) <- list(colnames(X_all), colnames(X_all))
  if (!nrow(df)) {
    return(list(
      coefficients = setNames(rep(0, p), colnames(X_all)),
      coefficient_cov = zero_cov,
      residual_sd = as.numeric(no_data_sd),
      n_train = 0L,
      n_terms = p,
      rank = 0L,
      status = "no_confirmed_rows",
      diagnostics = list()
    ))
  }

  X <- X_all[df$theta_row, , drop = FALSE]
  y <- as.numeric(df$error_fresh_minus_atlas)
  se <- sqrt(
    pmax(as.numeric(df$fresh_log_m_uncertainty), 0)^2 +
      pmax(as.numeric(df$atlas_se), 0)^2
  )
  se[!is.finite(se) | se <= 0] <- as.numeric(residual_floor)
  precision <- 1 / pmax(se, .Machine$double.eps)^2
  leverage_weight <- .local_chart_normalize_weights(df$theta_weight, nrow(df)) * nrow(df)
  w <- precision * leverage_weight
  finite_positive_w <- w[is.finite(w) & w > 0]
  w[!is.finite(w) | w <= 0] <- if (length(finite_positive_w)) min(finite_positive_w) else 1

  P <- diag(as.numeric(ridge_lambda) * penalty, p)
  A <- crossprod(X, X * w) + P
  b <- crossprod(X, y * w)
  beta <- as.numeric(.local_evidence_solve(A, b))
  names(beta) <- colnames(X_all)
  A_inv <- .local_evidence_solve(A, diag(p))
  dimnames(A_inv) <- list(colnames(X_all), colnames(X_all))
  fitted <- as.numeric(X %*% beta)
  residual <- y - fitted
  z <- residual / se
  edf <- sum(diag(X %*% A_inv %*% t(X * w)))
  df_resid <- max(nrow(X) - edf, 1)
  residual_scale <- sqrt(max(
    sum(w * residual^2) / max(sum(w), .Machine$double.eps),
    as.numeric(residual_floor)^2
  ))
  loo_error <- rep(NA_real_, nrow(df))
  if (nrow(df) > 1L) {
    for (j in seq_len(nrow(df))) {
      keep <- setdiff(seq_len(nrow(df)), j)
      Xj <- X[keep, , drop = FALSE]
      yj <- y[keep]
      wj <- w[keep]
      Aj <- crossprod(Xj, Xj * wj) + P
      bj <- crossprod(Xj, yj * wj)
      betaj <- as.numeric(.local_evidence_solve(Aj, bj))
      loo_error[j] <- y[j] - as.numeric(X[j, , drop = FALSE] %*% betaj)
    }
  }

  list(
    coefficients = beta,
    coefficient_cov = A_inv * max(1, sum(z^2, na.rm = TRUE) / df_resid),
    residual_sd = residual_scale,
    n_train = nrow(df),
    n_terms = p,
    rank = qr(X)$rank,
    status = "fitted",
    diagnostics = list(
      training_rmse = sqrt(mean(residual^2)),
      weighted_rmse = sqrt(sum(leverage_weight * residual^2) / sum(leverage_weight)),
      loo_rmse = .local_evidence_finite_rmse(loo_error),
      max_abs_residual_z = .local_evidence_finite_max(abs(z)),
      raw_error_range = .local_evidence_finite_range(y),
      fitted_error_range = .local_evidence_finite_range(fitted),
      edf = edf,
      df_resid = df_resid,
      residual_mean = mean(residual),
      residual_sd = stats::sd(residual)
    )
  )
}

validate_local_evidence_error_model <- function(error_model) {
  if (!inherits(error_model, "local_evidence_error_model")) {
    stop("error_model must inherit from 'local_evidence_error_model'.")
  }
  error_model$population_model <- normalize_population_model(error_model$population_model)
  if (!inherits(error_model$basis, "local_evidence_error_basis")) {
    stop("error_model basis is invalid.")
  }
  if (!is.list(error_model$models) || !length(error_model$models)) {
    stop("error_model must contain at least one local/member model.")
  }
  error_model
}

validate_total_evidence_error_model <- function(error_model) {
  if (!inherits(error_model, "total_evidence_error_model")) {
    stop("error_model must inherit from 'total_evidence_error_model'.")
  }
  error_model$population_model <- normalize_population_model(error_model$population_model)
  if (!inherits(error_model$basis, "local_evidence_error_basis")) {
    stop("total error model basis is invalid.")
  }
  if (!is.list(error_model$models) || !length(error_model$models)) {
    stop("error_model must contain at least one member-level total model.")
  }
  error_model
}

fit_local_evidence_error_model <- function(audit,
                                           active_hyper_names = NULL,
                                           include_quadratic = TRUE,
                                           include_interactions = FALSE,
                                           max_active_hyper = 6L,
                                           min_replicates = 2L,
                                           min_ok = min_replicates,
                                           allow_screening = FALSE,
                                           ridge_lambda = 1,
                                           intercept_penalty = 0.25,
                                           linear_penalty = 1,
                                           quadratic_penalty = 2,
                                           interaction_penalty = 4,
                                           residual_floor = 0.25,
                                           no_data_sd = 5,
                                           scale_floor = 1e-8) {
  audit <- validate_local_evidence_audit(audit)
  model <- audit$population_model
  active_hyper_names <- active_hyper_names %||%
    .local_evidence_default_active_hyper_names(audit, max_active_hyper = max_active_hyper)
  active_hyper_names <- intersect(unique(as.character(active_hyper_names)), model$hyper_names)
  if (!length(active_hyper_names)) {
    stop("No active hyperparameters available for the local evidence error model.")
  }
  basis <- .local_evidence_basis_spec(
    theta = audit$theta,
    population_model = model,
    theta_weights = audit$theta_metadata$theta_weight,
    active_hyper_names = active_hyper_names,
    include_quadratic = include_quadratic,
    include_interactions = include_interactions,
    scale_floor = scale_floor
  )
  X_all <- .local_evidence_basis_matrix(audit$theta, basis)
  penalty <- .local_evidence_basis_penalty(
    basis,
    intercept_penalty = intercept_penalty,
    linear_penalty = linear_penalty,
    quadratic_penalty = quadratic_penalty,
    interaction_penalty = interaction_penalty
  )
  rows <- audit$rows
  training_eligible <- is.finite(rows$error_fresh_minus_atlas) &
    is.finite(rows$fresh_log_m_uncertainty) &
    rows$n_ok >= as.integer(min_ok) &
    rows$n_replicates >= as.integer(min_replicates)
  if (isTRUE(allow_screening)) {
    training_eligible <- is.finite(rows$error_fresh_minus_atlas) &
      is.finite(rows$fresh_log_m_uncertainty) &
      rows$n_ok > 0L
  }
  rows$training_eligible <- training_eligible
  members <- unique(rows$member)
  local_metadata <- audit$local_metadata[order(audit$local_metadata$local_pos), , drop = FALSE]
  model_entries <- list()
  diag_rows <- list()
  for (member in members) {
    for (j in seq_len(nrow(local_metadata))) {
      local <- local_metadata$local[j]
      local_pos <- local_metadata$local_pos[j]
      df <- rows[rows$member == member &
                   rows$local_pos == local_pos &
                   rows$training_eligible, , drop = FALSE]
      fit <- .local_evidence_fit_error_one(
        df = df,
        X_all = X_all,
        ridge_lambda = ridge_lambda,
        penalty = penalty,
        residual_floor = residual_floor,
        no_data_sd = no_data_sd
      )
      key <- .local_evidence_matrix_key(member, local)
      model_entries[[key]] <- c(
        list(member = member, local = local, local_pos = local_pos),
        fit
      )
      diag <- fit$diagnostics
      diag_rows[[length(diag_rows) + 1L]] <- data.frame(
        member = member,
        local = local,
        local_pos = local_pos,
        status = fit$status,
        n_train = fit$n_train,
        n_terms = fit$n_terms,
        rank = fit$rank,
        training_rmse = as.numeric(diag$training_rmse %||% NA_real_),
        weighted_rmse = as.numeric(diag$weighted_rmse %||% NA_real_),
        loo_rmse = as.numeric(diag$loo_rmse %||% NA_real_),
        max_abs_residual_z = as.numeric(diag$max_abs_residual_z %||% NA_real_),
        raw_error_range = as.numeric(diag$raw_error_range %||% NA_real_),
        fitted_error_range = as.numeric(diag$fitted_error_range %||% NA_real_),
        residual_sd = fit$residual_sd,
        check.names = FALSE
      )
    }
  }
  diagnostics <- if (length(diag_rows)) do.call(rbind, diag_rows) else data.frame()
  training_fit <- rows[rows$training_eligible, , drop = FALSE]
  if (nrow(training_fit)) {
    fitted_delta <- numeric(nrow(training_fit))
    fitted_var <- numeric(nrow(training_fit))
    for (r in seq_len(nrow(training_fit))) {
      key <- .local_evidence_matrix_key(training_fit$member[r], training_fit$local[r])
      entry <- model_entries[[key]]
      xr <- X_all[training_fit$theta_row[r], , drop = FALSE]
      fitted_delta[r] <- as.numeric(xr %*% as.numeric(entry$coefficients))
      fitted_var[r] <- as.numeric(xr %*% as.matrix(entry$coefficient_cov) %*% t(xr)) +
        as.numeric(entry$residual_sd)^2
    }
    training_fit$fitted_delta <- fitted_delta
    training_fit$fitted_sd <- sqrt(pmax(fitted_var, 0))
    training_fit$model_residual <- training_fit$error_fresh_minus_atlas - training_fit$fitted_delta
    training_fit$model_residual_z <- training_fit$model_residual /
      pmax(training_fit$fitted_sd, .Machine$double.eps)
    theta_parts <- split(training_fit, list(training_fit$member, training_fit$theta_row), drop = TRUE)
    theta_summary <- do.call(rbind, lapply(theta_parts, function(df) {
      data.frame(
        member = df$member[1L],
        theta_row = df$theta_row[1L],
        theta_id = df$theta_id[1L],
        theta_label = df$theta_label[1L],
        theta_weight = df$theta_weight[1L],
        observed_total_delta = .local_evidence_finite_sum(df$error_fresh_minus_atlas),
        fitted_total_delta = .local_evidence_finite_sum(df$fitted_delta),
        residual_total_delta = .local_evidence_finite_sum(df$model_residual),
        posterior_weighted_local_rmse = sqrt(sum(df$theta_weight * df$model_residual^2) / sum(df$theta_weight)),
        max_abs_residual_z = .local_evidence_finite_max(abs(df$model_residual_z)),
        check.names = FALSE
      )
    }))
    theta_summary <- do.call(rbind, lapply(split(theta_summary, theta_summary$member), function(df) {
      df$centered_observed_total_delta <- df$observed_total_delta -
        .local_evidence_finite_mean(df$observed_total_delta)
      df$centered_fitted_total_delta <- df$fitted_total_delta -
        .local_evidence_finite_mean(df$fitted_total_delta)
      df$centered_residual_total_delta <- df$residual_total_delta -
        .local_evidence_finite_mean(df$residual_total_delta)
      df
    }))
    calibration_summary <- do.call(rbind, lapply(split(theta_summary, theta_summary$member), function(df) {
      data.frame(
        member = df$member[1L],
        n_training_rows = sum(training_fit$member == df$member[1L]),
        n_training_theta = nrow(df),
        centered_observed_total_range = .local_evidence_finite_range(df$centered_observed_total_delta),
        centered_fitted_total_range = .local_evidence_finite_range(df$centered_fitted_total_delta),
        centered_residual_total_range = .local_evidence_finite_range(df$centered_residual_total_delta),
        posterior_weighted_rmse = sqrt(
          sum(df$theta_weight * df$residual_total_delta^2) / sum(df$theta_weight)
        ),
        max_abs_residual_z = .local_evidence_finite_max(df$max_abs_residual_z),
        check.names = FALSE
      )
    }))
  } else {
    training_fit <- data.frame()
    theta_summary <- data.frame()
    calibration_summary <- data.frame()
  }
  structure(
    list(
      population_model = model,
      basis = basis,
      models = model_entries,
      members = members,
      locals = local_metadata,
      diagnostics = diagnostics,
      training_rows = training_fit,
      theta_diagnostics = theta_summary,
      calibration_summary = calibration_summary,
      settings = list(
        active_hyper_names = active_hyper_names,
        include_quadratic = isTRUE(include_quadratic),
        include_interactions = isTRUE(include_interactions),
        min_replicates = as.integer(min_replicates),
        min_ok = as.integer(min_ok),
        allow_screening = isTRUE(allow_screening),
        ridge_lambda = as.numeric(ridge_lambda),
        penalties = list(
          intercept = as.numeric(intercept_penalty),
          linear = as.numeric(linear_penalty),
          quadratic = as.numeric(quadratic_penalty),
          interaction = as.numeric(interaction_penalty)
        ),
        residual_floor = as.numeric(residual_floor),
        no_data_sd = as.numeric(no_data_sd)
      )
    ),
    class = "local_evidence_error_model"
  ) |>
    validate_local_evidence_error_model()
}

fit_total_evidence_error_model <- function(audit,
                                           active_hyper_names = NULL,
                                           include_quadratic = TRUE,
                                           include_interactions = FALSE,
                                           max_active_hyper = 6L,
                                           min_replicates = 2L,
                                           min_ok = min_replicates,
                                           allow_screening = FALSE,
                                           ridge_lambda = 1,
                                           intercept_penalty = 0.25,
                                           linear_penalty = 1,
                                           quadratic_penalty = 2,
                                           interaction_penalty = 4,
                                           residual_floor = 0.25,
                                           no_data_sd = 10,
                                           scale_floor = 1e-8) {
  audit <- validate_local_evidence_audit(audit)
  model <- audit$population_model
  active_hyper_names <- active_hyper_names %||%
    .local_evidence_default_active_hyper_names(audit, max_active_hyper = max_active_hyper)
  active_hyper_names <- intersect(unique(as.character(active_hyper_names)), model$hyper_names)
  if (!length(active_hyper_names)) {
    stop("No active hyperparameters available for the total evidence error model.")
  }
  basis <- .local_evidence_basis_spec(
    theta = audit$theta,
    population_model = model,
    theta_weights = audit$theta_metadata$theta_weight,
    active_hyper_names = active_hyper_names,
    include_quadratic = include_quadratic,
    include_interactions = include_interactions,
    scale_floor = scale_floor
  )
  X_all <- .local_evidence_basis_matrix(audit$theta, basis)
  penalty <- .local_evidence_basis_penalty(
    basis,
    intercept_penalty = intercept_penalty,
    linear_penalty = linear_penalty,
    quadratic_penalty = quadratic_penalty,
    interaction_penalty = interaction_penalty
  )
  rows <- audit$rows
  eligible <- is.finite(rows$error_fresh_minus_atlas) &
    is.finite(rows$fresh_log_m_uncertainty) &
    rows$n_ok >= as.integer(min_ok) &
    rows$n_replicates >= as.integer(min_replicates)
  if (isTRUE(allow_screening)) {
    eligible <- is.finite(rows$error_fresh_minus_atlas) &
      is.finite(rows$fresh_log_m_uncertainty) &
      rows$n_ok > 0L
  }
  rows <- rows[eligible, , drop = FALSE]
  if (!nrow(rows)) {
    stop("No eligible rows are available for total evidence error fitting.")
  }
  theta_rows <- split(rows, list(rows$member, rows$theta_row), drop = TRUE)
  total_rows <- do.call(rbind, lapply(theta_rows, function(df) {
    err <- df$error_fresh_minus_atlas[is.finite(df$error_fresh_minus_atlas)]
    var <- pmax(df$fresh_log_m_uncertainty, 0)^2 + pmax(df$atlas_se, 0)^2
    var <- var[is.finite(var)]
    data.frame(
      member = df$member[1L],
      theta_row = df$theta_row[1L],
      theta_id = df$theta_id[1L],
      theta_label = df$theta_label[1L],
      theta_weight = df$theta_weight[1L],
      n_locals = nrow(df),
      n_finite = length(err),
      total_error = if (length(err)) sum(err) else NA_real_,
      total_variance = if (length(var)) sum(var) else NA_real_,
      check.names = FALSE
    )
  }))
  total_rows <- total_rows[is.finite(total_rows$total_error), , drop = FALSE]
  members <- unique(total_rows$member)
  models <- list()
  diagnostics <- list()
  for (member in members) {
    df <- total_rows[total_rows$member == member, , drop = FALSE]
    X <- X_all[df$theta_row, , drop = FALSE]
    y <- as.numeric(df$total_error)
    se <- sqrt(pmax(df$total_variance, 0))
    se[!is.finite(se) | se <= 0] <- as.numeric(residual_floor)
    leverage_weight <- .local_chart_normalize_weights(df$theta_weight, nrow(df)) * nrow(df)
    w <- leverage_weight / pmax(se, .Machine$double.eps)^2
    finite_positive_w <- w[is.finite(w) & w > 0]
    w[!is.finite(w) | w <= 0] <- if (length(finite_positive_w)) min(finite_positive_w) else 1
    P <- diag(as.numeric(ridge_lambda) * penalty, ncol(X))
    A <- crossprod(X, X * w) + P
    b <- crossprod(X, y * w)
    beta <- as.numeric(.local_evidence_solve(A, b))
    names(beta) <- colnames(X_all)
    A_inv <- .local_evidence_solve(A, diag(ncol(X)))
    dimnames(A_inv) <- list(colnames(X_all), colnames(X_all))
    fitted <- as.numeric(X %*% beta)
    residual <- y - fitted
    z <- residual / se
    residual_sd <- sqrt(max(
      sum(w * residual^2) / max(sum(w), .Machine$double.eps),
      as.numeric(residual_floor)^2
    ))
    loo_error <- rep(NA_real_, nrow(df))
    if (nrow(df) > 1L) {
      for (j in seq_len(nrow(df))) {
        keep <- setdiff(seq_len(nrow(df)), j)
        Xj <- X[keep, , drop = FALSE]
        yj <- y[keep]
        wj <- w[keep]
        Aj <- crossprod(Xj, Xj * wj) + P
        bj <- crossprod(Xj, yj * wj)
        betaj <- as.numeric(.local_evidence_solve(Aj, bj))
        loo_error[j] <- y[j] - as.numeric(X[j, , drop = FALSE] %*% betaj)
      }
    }
    models[[member]] <- list(
      member = member,
      coefficients = beta,
      coefficient_cov = A_inv * max(1, sum(z^2, na.rm = TRUE) / max(nrow(X) - qr(X)$rank, 1)),
      residual_sd = residual_sd,
      n_train = nrow(df),
      n_terms = ncol(X),
      rank = qr(X)$rank,
      status = "fitted"
    )
    diagnostics[[length(diagnostics) + 1L]] <- data.frame(
      member = member,
      n_train = nrow(df),
      n_terms = ncol(X),
      rank = qr(X)$rank,
      training_rmse = sqrt(mean(residual^2)),
      weighted_rmse = sqrt(sum(leverage_weight * residual^2) / sum(leverage_weight)),
      loo_rmse = .local_evidence_finite_rmse(loo_error),
      max_abs_residual_z = .local_evidence_finite_max(abs(z)),
      observed_total_range = .local_evidence_finite_range(y),
      fitted_total_range = .local_evidence_finite_range(fitted),
      residual_total_range = .local_evidence_finite_range(residual),
      residual_sd = residual_sd,
      check.names = FALSE
    )
  }
  structure(
    list(
      population_model = model,
      basis = basis,
      models = models,
      members = members,
      training_rows = total_rows,
      diagnostics = if (length(diagnostics)) do.call(rbind, diagnostics) else data.frame(),
      settings = list(
        active_hyper_names = active_hyper_names,
        include_quadratic = isTRUE(include_quadratic),
        include_interactions = isTRUE(include_interactions),
        min_replicates = as.integer(min_replicates),
        min_ok = as.integer(min_ok),
        allow_screening = isTRUE(allow_screening),
        ridge_lambda = as.numeric(ridge_lambda),
        residual_floor = as.numeric(residual_floor),
        no_data_sd = as.numeric(no_data_sd)
      )
    ),
    class = "total_evidence_error_model"
  ) |>
    validate_total_evidence_error_model()
}

predict_local_evidence_error <- function(error_model,
                                         theta,
                                         member = NULL,
                                         locals = NULL,
                                         include_residual = TRUE,
                                         covariance = FALSE) {
  error_model <- validate_local_evidence_error_model(error_model)
  model <- error_model$population_model
  theta <- .as_hyper_matrix(theta, model$hyper_names, model$hyper_dim)
  member <- as.character(member %||% error_model$members[1L])
  local_frame <- error_model$locals
  if (!is.null(locals)) {
    local_frame <- local_frame[local_frame$local %in% as.character(locals), , drop = FALSE]
  }
  if (!nrow(local_frame)) {
    stop("No requested locals are present in the local evidence error model.")
  }
  X <- .local_evidence_basis_matrix(theta, error_model$basis)
  mean_mat <- matrix(0, nrow = nrow(theta), ncol = nrow(local_frame))
  var_mat <- matrix(0, nrow = nrow(theta), ncol = nrow(local_frame))
  colnames(mean_mat) <- local_frame$local
  colnames(var_mat) <- local_frame$local
  cov_list <- if (isTRUE(covariance)) vector("list", nrow(local_frame)) else NULL
  if (isTRUE(covariance)) names(cov_list) <- local_frame$local
  row_diag <- list()
  for (j in seq_len(nrow(local_frame))) {
    local <- local_frame$local[j]
    key <- .local_evidence_matrix_key(member, local)
    entry <- error_model$models[[key]]
    if (is.null(entry)) {
      stop("No local evidence error model for member=", member, " local=", local)
    }
    beta <- as.numeric(entry$coefficients)
    cov_beta <- as.matrix(entry$coefficient_cov)
    pred_mean <- as.numeric(X %*% beta)
    if (isTRUE(covariance)) {
      pred_cov <- X %*% cov_beta %*% t(X)
      pred_var <- pmax(diag(pred_cov), 0)
    } else {
      pred_cov <- NULL
      pred_var <- pmax(rowSums((X %*% cov_beta) * X), 0)
    }
    if (isTRUE(include_residual)) {
      pred_var <- pred_var + as.numeric(entry$residual_sd)^2
      if (isTRUE(covariance)) {
        diag(pred_cov) <- diag(pred_cov) + as.numeric(entry$residual_sd)^2
      }
    }
    mean_mat[, j] <- pred_mean
    var_mat[, j] <- pred_var
    if (isTRUE(covariance)) {
      cov_list[[j]] <- pred_cov
    }
    row_diag[[length(row_diag) + 1L]] <- data.frame(
      member = member,
      local = local,
      local_pos = local_frame$local_pos[j],
      n_train = entry$n_train,
      status = entry$status,
      residual_sd = entry$residual_sd,
      check.names = FALSE
    )
  }
  structure(
    list(
      mean = mean_mat,
      variance = var_mat,
      covariance = cov_list,
      member = member,
      locals = local_frame,
      theta = theta,
      model_diagnostics = if (length(row_diag)) do.call(rbind, row_diag) else data.frame()
    ),
    class = "local_evidence_error_prediction"
  )
}

predict_total_evidence_error <- function(error_model,
                                         theta,
                                         member = NULL,
                                         include_residual = TRUE,
                                         covariance = FALSE) {
  error_model <- validate_total_evidence_error_model(error_model)
  model <- error_model$population_model
  theta <- .as_hyper_matrix(theta, model$hyper_names, model$hyper_dim)
  member <- as.character(member %||% error_model$members[1L])
  entry <- error_model$models[[member]]
  if (is.null(entry)) {
    stop("No total evidence error model for member=", member)
  }
  X <- .local_evidence_basis_matrix(theta, error_model$basis)
  beta <- as.numeric(entry$coefficients)
  cov_beta <- as.matrix(entry$coefficient_cov)
  pred_mean <- as.numeric(X %*% beta)
  if (isTRUE(covariance)) {
    pred_cov <- X %*% cov_beta %*% t(X)
    pred_var <- pmax(diag(pred_cov), 0)
  } else {
    pred_cov <- NULL
    pred_var <- pmax(rowSums((X %*% cov_beta) * X), 0)
  }
  if (isTRUE(include_residual)) {
    pred_var <- pred_var + as.numeric(entry$residual_sd)^2
    if (isTRUE(covariance)) {
      diag(pred_cov) <- diag(pred_cov) + as.numeric(entry$residual_sd)^2
    }
  }
  structure(
    list(
      mean = pred_mean,
      variance = pred_var,
      covariance = pred_cov,
      member = member,
      theta = theta,
      model_diagnostics = data.frame(
        member = member,
        n_train = entry$n_train,
        status = entry$status,
        residual_sd = entry$residual_sd,
        check.names = FALSE
      )
    ),
    class = "total_evidence_error_prediction"
  )
}

build_corrected_local_atlas_factor_set <- function(factor_set,
                                                   error_model,
                                                   member = NULL,
                                                   include_residual_uncertainty = TRUE,
                                                   correction_scale = 1,
                                                   reject_uncertified = TRUE,
                                                   uncertified_loglik = -1e100) {
  factor_set <- validate_local_atlas_factor_set(factor_set)
  if (inherits(error_model, "total_evidence_error_model")) {
    error_model <- validate_total_evidence_error_model(error_model)
    error_model_type <- "total"
  } else {
    error_model <- validate_local_evidence_error_model(error_model)
    error_model_type <- "local"
  }
  if (!identical(factor_set$population_model$hyper_names, error_model$population_model$hyper_names) ||
      !identical(factor_set$population_model$alpha_names, error_model$population_model$alpha_names)) {
    stop("factor_set and error_model population models do not match.")
  }
  member <- as.character(member %||% error_model$members[1L])
  if (!member %in% error_model$members) {
    stop("member is not present in error_model.")
  }
  if (identical(error_model_type, "local")) {
    missing_locals <- setdiff(names(factor_set$atlases), error_model$locals$local)
    if (length(missing_locals)) {
      stop("error_model is missing factor_set locals: ", paste(missing_locals, collapse = ", "))
    }
  }
  structure(
    list(
      population_model = factor_set$population_model,
      base_factor_set = factor_set,
      error_model = error_model,
      error_model_type = error_model_type,
      member = member,
      include_residual_uncertainty = isTRUE(include_residual_uncertainty),
      correction_scale = as.numeric(correction_scale),
      reject_uncertified = isTRUE(reject_uncertified),
      uncertified_loglik = as.numeric(uncertified_loglik),
      n_locals = factor_set$n_locals,
      log_constant = factor_set$log_constant
    ),
    class = c("local_evidence_corrected_factor_set", "population_factor_set")
  )
}

validate_corrected_local_atlas_factor_set <- function(factor_set) {
  if (!inherits(factor_set, "local_evidence_corrected_factor_set")) {
    stop("factor_set must inherit from 'local_evidence_corrected_factor_set'.")
  }
  factor_set$base_factor_set <- validate_local_atlas_factor_set(factor_set$base_factor_set)
  if (inherits(factor_set$error_model, "total_evidence_error_model") ||
      identical(factor_set$error_model_type %||% NULL, "total")) {
    factor_set$error_model <- validate_total_evidence_error_model(factor_set$error_model)
    factor_set$error_model_type <- "total"
  } else {
    factor_set$error_model <- validate_local_evidence_error_model(factor_set$error_model)
    factor_set$error_model_type <- "local"
  }
  factor_set$population_model <- factor_set$base_factor_set$population_model
  factor_set$n_locals <- factor_set$base_factor_set$n_locals
  factor_set$log_constant <- factor_set$base_factor_set$log_constant
  factor_set$correction_scale <- as.numeric(factor_set$correction_scale %||% 1)
  if (length(factor_set$correction_scale) != 1L || !is.finite(factor_set$correction_scale)) {
    stop("corrected factor set correction_scale must be a finite scalar.")
  }
  factor_set$reject_uncertified <- isTRUE(factor_set$reject_uncertified %||% TRUE)
  factor_set$uncertified_loglik <- as.numeric(factor_set$uncertified_loglik %||% -1e100)
  if (length(factor_set$uncertified_loglik) != 1L || !is.finite(factor_set$uncertified_loglik)) {
    stop("corrected factor set uncertified_loglik must be a finite scalar.")
  }
  factor_set
}

corrected_factor_set_loglik_by_local <- function(factor_set,
                                                 theta,
                                                 include_constant = TRUE,
                                                 n_cores = 1L) {
  factor_set <- validate_corrected_local_atlas_factor_set(factor_set)
  model <- factor_set$population_model
  theta <- .as_hyper_matrix(theta, model$hyper_names, model$hyper_dim)
  base_factor_set <- factor_set$base_factor_set
  base_factor_set$evaluator_control$stop_on_uncertified <- FALSE
  parts <- .local_atlas_factor_set_by_local(base_factor_set, theta, n_cores = n_cores)
  bad <- .local_atlas_factor_set_uncertified(parts)
  base <- do.call(cbind, lapply(parts, `[[`, "log_marginal"))
  colnames(base) <- names(parts)
  correction_scale <- as.numeric(factor_set$correction_scale)
  if (identical(factor_set$error_model_type, "total")) {
    pred_total <- predict_total_evidence_error(
      error_model = factor_set$error_model,
      theta = theta,
      member = factor_set$member,
      include_residual = FALSE,
      covariance = FALSE
    )
    out <- base + correction_scale * pred_total$mean / max(ncol(base), 1L)
  } else {
    pred <- predict_local_evidence_error(
      error_model = factor_set$error_model,
      theta = theta,
      member = factor_set$member,
      locals = colnames(base),
      include_residual = FALSE,
      covariance = FALSE
    )
    out <- base + correction_scale * pred$mean[, colnames(base), drop = FALSE]
  }
  if (isTRUE(include_constant) && ncol(out)) {
    out <- out + factor_set$log_constant / ncol(out)
  }
  if (isTRUE(factor_set$reject_uncertified) && nrow(bad)) {
    bad_rows <- unique(as.integer(bad$theta_row))
    bad_rows <- bad_rows[is.finite(bad_rows) & bad_rows >= 1L & bad_rows <= nrow(out)]
    if (length(bad_rows) && ncol(out)) {
      out[bad_rows, ] <- 0
      out[bad_rows, 1L] <- factor_set$uncertified_loglik
    }
  }
  attr(out, "uncertified") <- bad
  out
}

corrected_factor_set_loglik <- function(factor_set,
                                        theta,
                                        include_constant = FALSE,
                                        n_cores = 1L) {
  factor_set <- validate_corrected_local_atlas_factor_set(factor_set)
  by_local <- corrected_factor_set_loglik_by_local(
    factor_set = factor_set,
    theta = theta,
    include_constant = FALSE,
    n_cores = n_cores
  )
  out <- rowSums(by_local)
  if (isTRUE(include_constant)) {
    out <- out + factor_set$log_constant
  }
  as.numeric(out)
}

corrected_factor_set_uncertainty <- function(factor_set,
                                             theta,
                                             by_local = FALSE,
                                             covariance = FALSE) {
  factor_set <- validate_corrected_local_atlas_factor_set(factor_set)
  model <- factor_set$population_model
  theta <- .as_hyper_matrix(theta, model$hyper_names, model$hyper_dim)
  local_names <- names(factor_set$base_factor_set$atlases)
  if (identical(factor_set$error_model_type, "total")) {
    pred <- predict_total_evidence_error(
      error_model = factor_set$error_model,
      theta = theta,
      member = factor_set$member,
      include_residual = isTRUE(factor_set$include_residual_uncertainty),
      covariance = covariance
    )
    total_variance <- as.numeric(pred$variance) * as.numeric(factor_set$correction_scale)^2
    out <- list(
      total_variance = total_variance,
      total_sd = sqrt(pmax(total_variance, 0)),
      member = factor_set$member,
      theta = theta
    )
    if (isTRUE(by_local)) {
      allocated <- matrix(
        total_variance / max(length(local_names), 1L),
        nrow = nrow(theta),
        ncol = length(local_names),
        dimnames = list(NULL, local_names)
      )
      out$by_local_variance <- allocated
      out$by_local_sd <- sqrt(pmax(allocated, 0))
    }
    if (isTRUE(covariance)) {
      out$total_covariance <- pred$covariance * as.numeric(factor_set$correction_scale)^2
      if (isTRUE(by_local)) {
        out$by_local_covariance <- stats::setNames(
          replicate(
            length(local_names),
            pred$covariance * as.numeric(factor_set$correction_scale)^2 /
              max(length(local_names), 1L)^2,
            simplify = FALSE
          ),
          local_names
        )
      }
    }
    return(structure(out, class = "local_evidence_corrected_uncertainty"))
  }
  pred <- predict_local_evidence_error(
    error_model = factor_set$error_model,
    theta = theta,
    member = factor_set$member,
    locals = local_names,
    include_residual = isTRUE(factor_set$include_residual_uncertainty),
    covariance = covariance
  )
  variance <- pred$variance[, local_names, drop = FALSE] * as.numeric(factor_set$correction_scale)^2
  total_variance <- rowSums(variance)
  out <- list(
    total_variance = as.numeric(total_variance),
    total_sd = sqrt(pmax(as.numeric(total_variance), 0)),
    member = factor_set$member,
    theta = theta
  )
  if (isTRUE(by_local)) {
    out$by_local_variance <- variance
    out$by_local_sd <- sqrt(pmax(variance, 0))
  }
  if (isTRUE(covariance)) {
    local_cov <- pred$covariance[local_names]
    out$by_local_covariance <- local_cov
    out$total_covariance <- Reduce(`+`, local_cov)
  }
  structure(out, class = "local_evidence_corrected_uncertainty")
}

.local_evidence_quadratic_form <- function(X, cov_beta, a, residual_sd = 0) {
  xb <- colSums(X * as.numeric(a))
  out <- as.numeric(t(xb) %*% as.matrix(cov_beta) %*% xb)
  out + as.numeric(residual_sd)^2 * sum(as.numeric(a)^2)
}

validate_corrected_outer_uncertainty <- function(uncertainty) {
  if (!inherits(uncertainty, "local_evidence_outer_uncertainty")) {
    stop("uncertainty must inherit from 'local_evidence_outer_uncertainty'.")
  }
  required <- c("evidence", "posterior", "local_contributors", "theta_diagnostics", "settings")
  missing <- setdiff(required, names(uncertainty))
  if (length(missing)) {
    stop("outer uncertainty object is missing: ", paste(missing, collapse = ", "))
  }
  uncertainty
}

summarize_corrected_outer_uncertainty <- function(factor_set,
                                                  fit = NULL,
                                                  theta = NULL,
                                                  weights = NULL,
                                                  parameters = NULL,
                                                  include_residual = factor_set$include_residual_uncertainty,
                                                  max_log_evidence_numerical_se = 1,
                                                  max_standardized_posterior_mean_se = 0.1,
                                                  max_posterior_rms_log_factor_sd = 1,
                                                  keep_theta_diagnostics = TRUE) {
  factor_set <- validate_corrected_local_atlas_factor_set(factor_set)
  model <- factor_set$population_model
  if (!is.null(fit)) {
    theta <- theta %||% fit$theta
    weights <- weights %||% fit$w
  }
  theta <- .as_hyper_matrix(theta, model$hyper_names, model$hyper_dim)
  weights <- .local_chart_normalize_weights(weights, nrow(theta))
  parameters <- parameters %||% model$hyper_names
  parameters <- intersect(unique(as.character(parameters)), model$hyper_names)
  if (!length(parameters)) {
    stop("parameters must identify at least one population hyperparameter.")
  }

  if (identical(factor_set$error_model_type, "total")) {
    error_model <- validate_total_evidence_error_model(factor_set$error_model)
    X <- .local_evidence_basis_matrix(theta, error_model$basis) *
      as.numeric(factor_set$correction_scale)
    entry <- error_model$models[[factor_set$member]]
    if (is.null(entry)) {
      stop("No total evidence error model for member=", factor_set$member)
    }
    cov_beta <- as.matrix(entry$coefficient_cov)
    residual_sd <- if (isTRUE(include_residual)) {
      abs(as.numeric(factor_set$correction_scale)) * as.numeric(entry$residual_sd)
    } else {
      0
    }
    total_point_variance <- pmax(rowSums((X %*% cov_beta) * X), 0) + residual_sd^2
    log_evidence_numerical_var <- .local_evidence_quadratic_form(
      X = X,
      cov_beta = cov_beta,
      a = weights,
      residual_sd = residual_sd
    )
    parameter_rows <- lapply(parameters, function(parameter) {
      value <- theta[, parameter]
      mean_value <- sum(weights * value)
      posterior_sd <- sqrt(sum(weights * (value - mean_value)^2))
      centered_value <- value - mean_value
      numerical_se <- sqrt(pmax(
        .local_evidence_quadratic_form(
          X = X,
          cov_beta = cov_beta,
          a = weights * centered_value,
          residual_sd = residual_sd
        ),
        0
      ))
      data.frame(
        parameter = parameter,
        posterior_mean = mean_value,
        posterior_sd = posterior_sd,
        numerical_mean_se = numerical_se,
        standardized_numerical_mean_se = numerical_se / max(posterior_sd, .Machine$double.eps),
        check.names = FALSE
      )
    })
    parameter_sensitivity <- do.call(rbind, parameter_rows)
    parameter_sensitivity <- parameter_sensitivity[
      order(-parameter_sensitivity$standardized_numerical_mean_se),
      ,
      drop = FALSE
    ]
    worst <- parameter_sensitivity[1L, , drop = FALSE]
    total_sd <- sqrt(pmax(total_point_variance, 0))
    total_sd_q <- .local_atlas_weighted_quantile(total_sd, weights, probs = c(0.5, 0.9, 0.99))
    posterior_rms_log_factor_sd <- sqrt(sum(weights * total_point_variance))
    log_evidence_mcse <- as.numeric(if (!is.null(fit)) fit$mcse_log_evidence %||% NA_real_ else NA_real_)
    evidence <- data.frame(
      log_evidence_mean_surface = as.numeric(if (!is.null(fit)) fit$log_evidence %||% NA_real_ else NA_real_),
      smc_mcse_log_evidence = log_evidence_mcse,
      numerical_se_log_evidence = sqrt(pmax(log_evidence_numerical_var, 0)),
      combined_se_log_evidence = sqrt(
        pmax(log_evidence_mcse, 0)^2 + pmax(log_evidence_numerical_var, 0)
      ),
      check.names = FALSE
    )
    posterior <- data.frame(
      n_theta = nrow(theta),
      posterior_mean_total_variance = sum(weights * total_point_variance),
      posterior_rms_log_factor_sd = posterior_rms_log_factor_sd,
      posterior_mean_log_factor_sd = sum(weights * total_sd),
      posterior_q50_log_factor_sd = total_sd_q[1L],
      posterior_q90_log_factor_sd = total_sd_q[2L],
      posterior_q99_log_factor_sd = total_sd_q[3L],
      worst_parameter = worst$parameter,
      worst_standardized_posterior_mean_se = worst$standardized_numerical_mean_se,
      check.names = FALSE
    )
    failures <- character()
    if (is.finite(evidence$numerical_se_log_evidence) &&
        evidence$numerical_se_log_evidence > as.numeric(max_log_evidence_numerical_se)) {
      failures <- c(failures, "log_evidence_numerical_se")
    }
    if (is.finite(posterior$worst_standardized_posterior_mean_se) &&
        posterior$worst_standardized_posterior_mean_se > as.numeric(max_standardized_posterior_mean_se)) {
      failures <- c(failures, "posterior_mean_sensitivity")
    }
    if (is.finite(posterior$posterior_rms_log_factor_sd) &&
        posterior$posterior_rms_log_factor_sd > as.numeric(max_posterior_rms_log_factor_sd)) {
      failures <- c(failures, "posterior_log_factor_uncertainty")
    }
    theta_diagnostics <- if (isTRUE(keep_theta_diagnostics)) {
      data.frame(
        theta_row = seq_len(nrow(theta)),
        theta_weight = weights,
        total_variance = total_point_variance,
        total_sd = total_sd,
        check.names = FALSE
      )
    } else {
      data.frame()
    }
    local_contributors <- data.frame(
      member = factor_set$member,
      local = "TOTAL",
      local_pos = NA_integer_,
      status = entry$status,
      n_train = entry$n_train,
      pointwise_variance_mean = sum(weights * total_point_variance),
      pointwise_sd_rms = posterior_rms_log_factor_sd,
      pointwise_sd_q50 = total_sd_q[1L],
      pointwise_sd_q90 = total_sd_q[2L],
      pointwise_sd_max = max(total_sd),
      log_evidence_numerical_var = log_evidence_numerical_var,
      log_evidence_numerical_se = sqrt(pmax(log_evidence_numerical_var, 0)),
      residual_sd = residual_sd,
      check.names = FALSE
    )
    return(structure(
      list(
        evidence = evidence,
        posterior = posterior,
        parameter_sensitivity = parameter_sensitivity,
        local_contributors = local_contributors,
        theta_diagnostics = theta_diagnostics,
        certification = data.frame(
          status = if (length(failures)) "uncertain" else "negligible",
          failures = paste(failures, collapse = ","),
          check.names = FALSE
        ),
        settings = list(
          include_residual = isTRUE(include_residual),
          covariance_assumption = "single total evidence error model",
          max_log_evidence_numerical_se = as.numeric(max_log_evidence_numerical_se),
          max_standardized_posterior_mean_se = as.numeric(max_standardized_posterior_mean_se),
          max_posterior_rms_log_factor_sd = as.numeric(max_posterior_rms_log_factor_sd)
        )
      ),
      class = "local_evidence_outer_uncertainty"
    ) |>
      validate_corrected_outer_uncertainty())
  }

  error_model <- validate_local_evidence_error_model(factor_set$error_model)
  X <- .local_evidence_basis_matrix(theta, error_model$basis) *
    as.numeric(factor_set$correction_scale)
  local_names <- names(factor_set$base_factor_set$atlases)
  total_point_variance <- numeric(nrow(theta))
  local_rows <- list()
  parameter_accum <- lapply(parameters, function(name) 0)
  names(parameter_accum) <- parameters

  for (local in local_names) {
    key <- .local_evidence_matrix_key(factor_set$member, local)
    entry <- error_model$models[[key]]
    if (is.null(entry)) {
      stop("No local evidence error model for member=", factor_set$member, " local=", local)
    }
    cov_beta <- as.matrix(entry$coefficient_cov)
    residual_sd <- if (isTRUE(include_residual)) {
      abs(as.numeric(factor_set$correction_scale)) * as.numeric(entry$residual_sd)
    } else {
      0
    }
    diag_var <- pmax(rowSums((X %*% cov_beta) * X), 0) + residual_sd^2
    total_point_variance <- total_point_variance + diag_var
    logZ_var <- .local_evidence_quadratic_form(
      X = X,
      cov_beta = cov_beta,
      a = weights,
      residual_sd = residual_sd
    )
    parameter_vars <- numeric(length(parameters))
    names(parameter_vars) <- parameters
    for (parameter in parameters) {
      value <- theta[, parameter]
      centered_value <- value - sum(weights * value)
      a <- weights * centered_value
      parameter_vars[parameter] <- .local_evidence_quadratic_form(
        X = X,
        cov_beta = cov_beta,
        a = a,
        residual_sd = residual_sd
      )
      parameter_accum[[parameter]] <- parameter_accum[[parameter]] + parameter_vars[parameter]
    }
    local_rows[[length(local_rows) + 1L]] <- data.frame(
      member = factor_set$member,
      local = local,
      local_pos = entry$local_pos,
      status = entry$status,
      n_train = entry$n_train,
      pointwise_variance_mean = sum(weights * diag_var),
      pointwise_sd_rms = sqrt(sum(weights * diag_var)),
      pointwise_sd_q50 = .local_atlas_weighted_quantile(sqrt(pmax(diag_var, 0)), weights, 0.50),
      pointwise_sd_q90 = .local_atlas_weighted_quantile(sqrt(pmax(diag_var, 0)), weights, 0.90),
      pointwise_sd_max = max(sqrt(pmax(diag_var, 0))),
      log_evidence_numerical_var = logZ_var,
      log_evidence_numerical_se = sqrt(pmax(logZ_var, 0)),
      residual_sd = residual_sd,
      check.names = FALSE
    )
  }

  local_contributors <- if (length(local_rows)) do.call(rbind, local_rows) else data.frame()
  if (nrow(local_contributors)) {
    local_contributors <- local_contributors[
      order(-local_contributors$log_evidence_numerical_se,
            -local_contributors$pointwise_sd_rms),
      ,
      drop = FALSE
    ]
  }
  log_evidence_numerical_var <- sum(local_contributors$log_evidence_numerical_var)
  total_sd <- sqrt(pmax(total_point_variance, 0))
  total_sd_q <- .local_atlas_weighted_quantile(total_sd, weights, probs = c(0.5, 0.9, 0.99))
  posterior_rms_log_factor_sd <- sqrt(sum(weights * total_point_variance))

  parameter_rows <- lapply(parameters, function(parameter) {
    value <- theta[, parameter]
    mean_value <- sum(weights * value)
    posterior_sd <- sqrt(sum(weights * (value - mean_value)^2))
    numerical_se <- sqrt(pmax(parameter_accum[[parameter]], 0))
    data.frame(
      parameter = parameter,
      posterior_mean = mean_value,
      posterior_sd = posterior_sd,
      numerical_mean_se = numerical_se,
      standardized_numerical_mean_se = numerical_se / max(posterior_sd, .Machine$double.eps),
      check.names = FALSE
    )
  })
  parameter_sensitivity <- do.call(rbind, parameter_rows)
  parameter_sensitivity <- parameter_sensitivity[
    order(-parameter_sensitivity$standardized_numerical_mean_se),
    ,
    drop = FALSE
  ]
  worst <- parameter_sensitivity[1L, , drop = FALSE]

  log_evidence_mcse <- as.numeric(if (!is.null(fit)) fit$mcse_log_evidence %||% NA_real_ else NA_real_)
  evidence <- data.frame(
    log_evidence_mean_surface = as.numeric(if (!is.null(fit)) fit$log_evidence %||% NA_real_ else NA_real_),
    smc_mcse_log_evidence = log_evidence_mcse,
    numerical_se_log_evidence = sqrt(pmax(log_evidence_numerical_var, 0)),
    combined_se_log_evidence = sqrt(
      pmax(log_evidence_mcse, 0)^2 + pmax(log_evidence_numerical_var, 0)
    ),
    check.names = FALSE
  )
  posterior <- data.frame(
    n_theta = nrow(theta),
    posterior_mean_total_variance = sum(weights * total_point_variance),
    posterior_rms_log_factor_sd = posterior_rms_log_factor_sd,
    posterior_mean_log_factor_sd = sum(weights * total_sd),
    posterior_q50_log_factor_sd = total_sd_q[1L],
    posterior_q90_log_factor_sd = total_sd_q[2L],
    posterior_q99_log_factor_sd = total_sd_q[3L],
    worst_parameter = worst$parameter,
    worst_standardized_posterior_mean_se = worst$standardized_numerical_mean_se,
    check.names = FALSE
  )
  failures <- character()
  if (is.finite(evidence$numerical_se_log_evidence) &&
      evidence$numerical_se_log_evidence > as.numeric(max_log_evidence_numerical_se)) {
    failures <- c(failures, "log_evidence_numerical_se")
  }
  if (is.finite(posterior$worst_standardized_posterior_mean_se) &&
      posterior$worst_standardized_posterior_mean_se > as.numeric(max_standardized_posterior_mean_se)) {
    failures <- c(failures, "posterior_mean_sensitivity")
  }
  if (is.finite(posterior$posterior_rms_log_factor_sd) &&
      posterior$posterior_rms_log_factor_sd > as.numeric(max_posterior_rms_log_factor_sd)) {
    failures <- c(failures, "posterior_log_factor_uncertainty")
  }
  certification <- data.frame(
    status = if (length(failures)) "uncertain" else "negligible",
    failures = paste(failures, collapse = ","),
    check.names = FALSE
  )
  theta_diagnostics <- if (isTRUE(keep_theta_diagnostics)) {
    data.frame(
      theta_row = seq_len(nrow(theta)),
      theta_weight = weights,
      total_variance = total_point_variance,
      total_sd = total_sd,
      check.names = FALSE
    )
  } else {
    data.frame()
  }

  structure(
    list(
      evidence = evidence,
      posterior = posterior,
      parameter_sensitivity = parameter_sensitivity,
      local_contributors = local_contributors,
      theta_diagnostics = theta_diagnostics,
      certification = certification,
      settings = list(
        include_residual = isTRUE(include_residual),
        covariance_assumption = "independent local evidence error models",
        max_log_evidence_numerical_se = as.numeric(max_log_evidence_numerical_se),
        max_standardized_posterior_mean_se = as.numeric(max_standardized_posterior_mean_se),
        max_posterior_rms_log_factor_sd = as.numeric(max_posterior_rms_log_factor_sd)
      )
    ),
    class = "local_evidence_outer_uncertainty"
  ) |>
    validate_corrected_outer_uncertainty()
}

run_corrected_outer_smc_with_uncertainty <- function(factor_set,
                                                     ...,
                                                     uncertainty_control = list()) {
  factor_set <- validate_corrected_local_atlas_factor_set(factor_set)
  fit <- outer_population_smc(factor_set = factor_set, ...)
  fit$numerical_uncertainty <- do.call(
    summarize_corrected_outer_uncertainty,
    c(list(factor_set = factor_set, fit = fit), uncertainty_control)
  )
  fit$log_evidence_numerical_se <- fit$numerical_uncertainty$evidence$numerical_se_log_evidence
  fit$log_evidence_total_se <- fit$numerical_uncertainty$evidence$combined_se_log_evidence
  fit
}

.local_evidence_select_correction_theta <- function(factor_set,
                                                    fit,
                                                    n_theta,
                                                    selection = c("weighted_resample", "top_weight", "uncertainty", "hybrid"),
                                                    seed = 123L) {
  factor_set <- validate_corrected_local_atlas_factor_set(factor_set)
  selection <- match.arg(selection)
  model <- factor_set$population_model
  theta <- .as_hyper_matrix(fit$theta, model$hyper_names, model$hyper_dim)
  w <- .local_chart_normalize_weights(fit$w, nrow(theta))
  n_theta <- min(as.integer(n_theta), nrow(theta))
  if (n_theta <= 0L) {
    stop("n_theta must be positive.")
  }
  set.seed(as.integer(seed))
  selected <- integer()
  multiplicity <- integer()
  if (identical(selection, "weighted_resample")) {
    draw <- sample.int(nrow(theta), size = n_theta, replace = TRUE, prob = w)
    tab <- table(draw)
    selected <- as.integer(names(tab))
    multiplicity <- as.integer(tab)
    theta_weight <- multiplicity / sum(multiplicity)
  } else if (identical(selection, "top_weight")) {
    selected <- head(order(w, decreasing = TRUE), n_theta)
    multiplicity <- rep(1L, length(selected))
    theta_weight <- .local_chart_normalize_weights(w[selected], length(selected))
  } else if (identical(selection, "uncertainty")) {
    uncertainty <- corrected_factor_set_uncertainty(factor_set, theta, by_local = FALSE)
    score <- w^0.5 * uncertainty$total_sd
    selected <- head(order(score, decreasing = TRUE), n_theta)
    multiplicity <- rep(1L, length(selected))
    theta_weight <- .local_chart_normalize_weights(w[selected], length(selected))
  } else {
    n_sample <- max(1L, ceiling(0.7 * n_theta))
    draw <- sample.int(nrow(theta), size = n_sample, replace = TRUE, prob = w)
    uncertainty <- corrected_factor_set_uncertainty(factor_set, theta, by_local = FALSE)
    score <- w^0.5 * uncertainty$total_sd
    deterministic <- head(order(score, decreasing = TRUE), max(0L, n_theta - length(unique(draw))))
    combined <- c(draw, deterministic)
    tab <- table(combined)
    selected <- as.integer(names(tab))
    multiplicity <- as.integer(tab)
    empirical <- multiplicity / sum(multiplicity)
    posterior <- .local_chart_normalize_weights(w[selected], length(selected))
    theta_weight <- .local_chart_normalize_weights(empirical + posterior, length(selected))
  }
  theta_sel <- theta[selected, , drop = FALSE]
  metadata <- data.frame(
    theta_row = seq_len(nrow(theta_sel)),
    theta_id = sprintf("da_theta_%03d", seq_len(nrow(theta_sel))),
    theta_label = sprintf("da_theta_%03d", seq_len(nrow(theta_sel))),
    theta_weight = theta_weight,
    theta_source = paste0("outer_", selection),
    outer_row = selected,
    outer_weight = w[selected],
    multiplicity = multiplicity,
    check.names = FALSE
  )
  metadata <- cbind(metadata, as.data.frame(theta_sel, check.names = FALSE))
  list(theta = theta_sel, metadata = metadata, selection = selection)
}

validate_local_evidence_correction_cache <- function(cache) {
  if (!inherits(cache, "local_evidence_correction_cache")) {
    stop("cache must inherit from 'local_evidence_correction_cache'.")
  }
  cache$factor_set <- validate_corrected_local_atlas_factor_set(cache$factor_set)
  model <- cache$factor_set$population_model
  cache$theta <- .as_hyper_matrix(cache$theta, model$hyper_names, model$hyper_dim)
  for (field in c("rows", "replicates", "theta_metadata", "local_metadata")) {
    if (!is.data.frame(cache[[field]])) {
      stop(field, " must be a data frame.")
    }
  }
  required <- c(
    "local", "local_pos", "theta_row", "theta_id", "theta_weight",
    "corrected_log_marginal", "fresh_log_m_center", "fresh_log_m_uncertainty",
    "delta_fresh_minus_corrected"
  )
  missing <- setdiff(required, names(cache$rows))
  if (length(missing)) {
    stop("correction cache rows are missing columns: ", paste(missing, collapse = ", "))
  }
  cache
}

cache_local_evidence_corrections <- function(factor_set,
                                             theta,
                                             data_list,
                                             loglik_fn,
                                             local_ids = NULL,
                                             theta_weights = NULL,
                                             theta_metadata = NULL,
                                             n_replicates = 2L,
                                             M = 500L,
                                             local_control = list(),
                                             uncertainty_floor = 1e-6,
                                             n_cores = 1L,
                                             seed = 123L,
                                             stop_on_error = TRUE,
                                             verbose = FALSE) {
  factor_set <- validate_corrected_local_atlas_factor_set(factor_set)
  model <- factor_set$population_model
  theta_input <- .local_evidence_theta_inputs(
    theta,
    population_model = model,
    theta_weights = theta_weights,
    theta_metadata = theta_metadata,
    reference_source = "delayed_acceptance"
  )
  theta <- theta_input$theta
  theta_metadata <- theta_input$metadata
  local_metadata <- .local_evidence_resolve_locals(
    factor_set$base_factor_set,
    data_list = data_list,
    local_ids = local_ids %||% seq_along(factor_set$base_factor_set$atlases)
  )
  pairs <- merge(
    expand.grid(
      local_pos = local_metadata$local_pos,
      theta_row = seq_len(nrow(theta)),
      KEEP.OUT.ATTRS = FALSE
    ),
    local_metadata[, c("local_pos", "data_pos", "local"), drop = FALSE],
    by = "local_pos",
    all.x = TRUE,
    sort = FALSE
  )
  replicates <- run_local_evidence_replicates(
    theta = theta,
    data_list = data_list,
    loglik_fn = loglik_fn,
    population_model = model,
    pairs = pairs,
    n_replicates = n_replicates,
    M = M,
    local_control = local_control,
    n_cores = n_cores,
    seed = seed,
    source = "delayed_acceptance_correction",
    stop_on_error = stop_on_error,
    verbose = verbose
  )
  summary <- .local_evidence_replicate_summary(
    replicates,
    uncertainty_floor = uncertainty_floor
  )
  corrected <- corrected_factor_set_loglik_by_local(
    factor_set = factor_set,
    theta = theta,
    include_constant = FALSE,
    n_cores = n_cores
  )
  rows <- merge(
    pairs[, c("local_pos", "data_pos", "local", "theta_row"), drop = FALSE],
    theta_metadata,
    by = "theta_row",
    all.x = TRUE,
    sort = FALSE
  )
  rows <- merge(
    rows,
    summary,
    by = c("local_pos", "local", "theta_row"),
    all.x = TRUE,
    sort = FALSE
  )
  rows$corrected_log_marginal <- mapply(
    function(theta_row, local) corrected[theta_row, local],
    rows$theta_row,
    rows$local
  )
  rows$delta_fresh_minus_corrected <- rows$fresh_log_m_center - rows$corrected_log_marginal
  rows$delta_variance <- pmax(rows$fresh_log_m_uncertainty, 0)^2
  rows$delta_se <- sqrt(pmax(rows$delta_variance, 0))
  rows$complete <- is.finite(rows$delta_fresh_minus_corrected) & rows$n_ok > 0L
  rows <- rows[order(rows$theta_row, rows$local_pos), , drop = FALSE]
  structure(
    list(
      factor_set = factor_set,
      theta = theta,
      theta_metadata = theta_metadata,
      local_metadata = local_metadata,
      rows = rows,
      replicates = replicates,
      settings = list(
        n_replicates = as.integer(n_replicates),
        M = as.integer(M),
        local_control = local_control,
        uncertainty_floor = as.numeric(uncertainty_floor),
        seed = as.integer(seed)
      )
    ),
    class = "local_evidence_correction_cache"
  ) |>
    validate_local_evidence_correction_cache()
}

validate_theta_correction_estimate <- function(correction) {
  if (!inherits(correction, "local_evidence_theta_correction")) {
    stop("correction must inherit from 'local_evidence_theta_correction'.")
  }
  for (field in c("theta_summary", "evidence", "posterior_weights", "local_summary", "settings")) {
    if (is.null(correction[[field]])) {
      stop("theta correction is missing: ", field)
    }
  }
  correction
}

estimate_theta_correction <- function(cache,
                                      theta_weights = NULL,
                                      require_complete = TRUE,
                                      corrected_log_evidence = NULL) {
  cache <- validate_local_evidence_correction_cache(cache)
  rows <- cache$rows
  if (!is.null(theta_weights)) {
    theta_weights <- .local_chart_normalize_weights(theta_weights, nrow(cache$theta))
    rows$theta_weight <- theta_weights[rows$theta_row]
  }
  theta_summary <- do.call(rbind, lapply(split(rows, rows$theta_row), function(df) {
    finite <- is.finite(df$delta_fresh_minus_corrected)
    complete <- all(df$complete) && all(finite)
    data.frame(
      theta_row = df$theta_row[1L],
      theta_id = df$theta_id[1L],
      theta_label = df$theta_label[1L],
      theta_weight = df$theta_weight[1L],
      n_locals = nrow(df),
      n_complete_locals = sum(finite),
      complete = complete,
      total_log_correction = if (any(finite)) sum(df$delta_fresh_minus_corrected[finite]) else NA_real_,
      total_variance = if (any(finite)) sum(df$delta_variance[finite]) else NA_real_,
      max_abs_local_correction = .local_evidence_finite_max(abs(df$delta_fresh_minus_corrected)),
      check.names = FALSE
    )
  }))
  theta_summary$total_se <- sqrt(pmax(theta_summary$total_variance, 0))
  usable <- is.finite(theta_summary$total_log_correction)
  if (isTRUE(require_complete)) {
    usable <- usable & theta_summary$complete
  }
  if (!any(usable)) {
    stop("No usable theta corrections are available.")
  }
  weight <- .local_chart_normalize_weights(theta_summary$theta_weight[usable], sum(usable))
  log_ratio <- theta_summary$total_log_correction[usable]
  log_terms <- log(pmax(weight, .Machine$double.eps)) + log_ratio
  log_norm <- logsumexp(log_terms)
  corrected_weight <- exp(log_terms - log_norm)
  ratio <- exp(log_ratio - max(log_ratio))
  mean_ratio_scaled <- sum(weight * ratio)
  theta_neff <- 1 / sum(weight * weight)
  theta_mcse <- sqrt(
    max(sum(weight * (ratio - mean_ratio_scaled)^2), 0) /
      max(theta_neff, 1) /
      max(mean_ratio_scaled^2, .Machine$double.eps)
  )
  local_probe_var <- sum(corrected_weight^2 * theta_summary$total_variance[usable])
  log_correction_se <- sqrt(theta_mcse^2 + local_probe_var)
  psis_k <- .local_chart_psis_k(log_ratio)
  corrected_log_evidence <- as.numeric(corrected_log_evidence %||% NA_real_)
  evidence <- data.frame(
    corrected_mean_log_evidence = corrected_log_evidence,
    local_log_evidence_correction = log_norm,
    corrected_log_evidence = if (is.finite(corrected_log_evidence)) corrected_log_evidence + log_norm else NA_real_,
    correction_theta_mcse = theta_mcse,
    correction_probe_se = sqrt(pmax(local_probe_var, 0)),
    correction_total_se = log_correction_se,
    correction_ess = ESS(corrected_weight),
    correction_ess_fraction = ESS(corrected_weight) / length(corrected_weight),
    correction_psis_k = psis_k,
    n_theta = nrow(theta_summary),
    n_usable_theta = sum(usable),
    n_corrected_locals = length(unique(rows$local)),
    n_total_locals = cache$factor_set$n_locals,
    local_coverage_fraction = length(unique(rows$local)) / cache$factor_set$n_locals,
    complete_local_correction = length(unique(rows$local)) == cache$factor_set$n_locals,
    check.names = FALSE
  )
  posterior_weights <- theta_summary[, c("theta_row", "theta_id", "theta_label", "theta_weight"), drop = FALSE]
  posterior_weights$usable <- usable
  posterior_weights$stage2_weight <- 0
  posterior_weights$stage2_weight[usable] <- corrected_weight
  local_summary <- do.call(rbind, lapply(split(rows, rows$local), function(df) {
    finite <- is.finite(df$delta_fresh_minus_corrected)
    data.frame(
      local = df$local[1L],
      local_pos = df$local_pos[1L],
      n_theta = length(unique(df$theta_row)),
      n_complete = sum(finite),
      mean_delta = .local_evidence_finite_mean(df$delta_fresh_minus_corrected),
      rmse_delta = .local_evidence_finite_rmse(df$delta_fresh_minus_corrected),
      max_abs_delta = .local_evidence_finite_max(abs(df$delta_fresh_minus_corrected)),
      mean_se = .local_evidence_finite_mean(df$delta_se),
      check.names = FALSE
    )
  }))
  local_summary <- local_summary[order(-local_summary$rmse_delta, -local_summary$max_abs_delta), , drop = FALSE]
  structure(
    list(
      theta_summary = theta_summary,
      evidence = evidence,
      posterior_weights = posterior_weights,
      local_summary = local_summary,
      cache = cache,
      settings = list(require_complete = isTRUE(require_complete))
    ),
    class = "local_evidence_theta_correction"
  ) |>
    validate_theta_correction_estimate()
}

validate_local_evidence_accounting <- function(account) {
  if (!inherits(account, "local_evidence_accounting")) {
    stop("account must inherit from 'local_evidence_accounting'.")
  }
  required <- c("evidence", "posterior_diagnostics", "hard_local_contributors", "settings")
  missing <- setdiff(required, names(account))
  if (length(missing)) {
    stop("evidence account is missing: ", paste(missing, collapse = ", "))
  }
  account
}

build_local_evidence_accounting <- function(fit,
                                            correction = NULL,
                                            uncertainty = fit$numerical_uncertainty %||% NULL,
                                            model_comparison_scale = 1) {
  correction <- if (!is.null(correction)) validate_theta_correction_estimate(correction) else NULL
  uncertainty <- if (!is.null(uncertainty)) validate_corrected_outer_uncertainty(uncertainty) else NULL
  mean_logZ <- as.numeric(fit$log_evidence %||% NA_real_)
  outer_se <- as.numeric(fit$mcse_log_evidence %||% NA_real_)
  local_correction <- as.numeric(if (!is.null(correction)) correction$evidence$local_log_evidence_correction else 0)
  correction_se <- as.numeric(if (!is.null(correction)) correction$evidence$correction_total_se else 0)
  numerical_se <- as.numeric(if (!is.null(uncertainty)) uncertainty$evidence$numerical_se_log_evidence else NA_real_)
  local_numerical_se <- sqrt(pmax(correction_se, 0)^2 + pmax(numerical_se, 0)^2)
  total_se <- sqrt(pmax(outer_se, 0)^2 + pmax(local_numerical_se, 0)^2)
  logZ <- mean_logZ + local_correction
  complete_correction <- is.null(correction) || isTRUE(correction$evidence$complete_local_correction)
  precise <- is.finite(total_se) &&
    total_se <= as.numeric(model_comparison_scale) &&
    isTRUE(complete_correction)
  evidence <- data.frame(
    log_evidence_corrected_mean_surface = mean_logZ,
    local_numerical_correction = local_correction,
    log_evidence_estimate = logZ,
    outer_mcse = outer_se,
    surrogate_numerical_se = numerical_se,
    stage2_correction_se = correction_se,
    local_evidence_numerical_se = local_numerical_se,
    total_log_evidence_se = total_se,
    precise_bayes_factor = precise,
    check.names = FALSE
  )
  posterior_diagnostics <- data.frame(
    correction_ess_fraction = as.numeric(if (!is.null(correction)) correction$evidence$correction_ess_fraction else NA_real_),
    correction_psis_k = as.numeric(if (!is.null(correction)) correction$evidence$correction_psis_k else NA_real_),
    correction_local_coverage_fraction = as.numeric(if (!is.null(correction)) correction$evidence$local_coverage_fraction else NA_real_),
    complete_local_correction = as.logical(if (!is.null(correction)) correction$evidence$complete_local_correction else FALSE),
    posterior_rms_log_factor_sd = as.numeric(if (!is.null(uncertainty)) uncertainty$posterior$posterior_rms_log_factor_sd else NA_real_),
    worst_parameter = as.character(if (!is.null(uncertainty)) uncertainty$posterior$worst_parameter else NA_character_),
    worst_standardized_posterior_mean_se = as.numeric(if (!is.null(uncertainty)) uncertainty$posterior$worst_standardized_posterior_mean_se else NA_real_),
    n_replicated_probes = as.integer(if (!is.null(correction)) nrow(correction$cache$replicates) else 0L),
    n_corrected_theta = as.integer(if (!is.null(correction)) correction$evidence$n_usable_theta else 0L),
    check.names = FALSE
  )
  hard_from_uncertainty <- if (!is.null(uncertainty)) {
    uncertainty$local_contributors[, c("local", "local_pos", "log_evidence_numerical_se", "pointwise_sd_rms"), drop = FALSE]
  } else {
    data.frame()
  }
  hard_from_correction <- if (!is.null(correction)) {
    correction$local_summary[, c("local", "local_pos", "rmse_delta", "max_abs_delta", "mean_se"), drop = FALSE]
  } else {
    data.frame()
  }
  hard_local_contributors <- merge(
    hard_from_uncertainty,
    hard_from_correction,
    by = c("local", "local_pos"),
    all = TRUE,
    sort = FALSE
  )
  if (nrow(hard_local_contributors)) {
    score <- rowSums(
      data.frame(
        pmax(hard_local_contributors$log_evidence_numerical_se, 0, na.rm = TRUE),
        pmax(hard_local_contributors$rmse_delta, 0, na.rm = TRUE)
      ),
      na.rm = TRUE
    )
    hard_local_contributors <- hard_local_contributors[order(-score), , drop = FALSE]
  }
  structure(
    list(
      evidence = evidence,
      posterior_diagnostics = posterior_diagnostics,
      hard_local_contributors = hard_local_contributors,
      correction = correction,
      uncertainty = uncertainty,
      settings = list(model_comparison_scale = as.numeric(model_comparison_scale))
    ),
    class = "local_evidence_accounting"
  ) |>
    validate_local_evidence_accounting()
}

outer_population_delayed_acceptance <- function(factor_set,
                                                data_list,
                                                loglik_fn,
                                                fit = NULL,
                                                n_theta = 64L,
                                                selection = c("weighted_resample", "top_weight", "uncertainty", "hybrid"),
                                                local_ids = NULL,
                                                n_replicates = 2L,
                                                M = 500L,
                                                local_control = list(),
                                                uncertainty_floor = 1e-6,
                                                outer_control = list(),
                                                uncertainty_control = list(),
                                                model_comparison_scale = 1,
                                                n_cores = 1L,
                                                seed = 123L,
                                                stop_on_error = TRUE,
                                                verbose = FALSE) {
  factor_set <- validate_corrected_local_atlas_factor_set(factor_set)
  selection <- match.arg(selection)
  if (is.null(fit)) {
    fit <- do.call(
      outer_population_smc,
      modifyList(
        list(
          factor_set = factor_set,
          n_cores = as.integer(n_cores),
          seed = as.integer(seed),
          verbose = isTRUE(verbose)
        ),
        outer_control
      )
    )
  }
  if (is.null(fit$numerical_uncertainty)) {
    fit$numerical_uncertainty <- do.call(
      summarize_corrected_outer_uncertainty,
      c(list(factor_set = factor_set, fit = fit), uncertainty_control)
    )
  }
  selected <- .local_evidence_select_correction_theta(
    factor_set = factor_set,
    fit = fit,
    n_theta = n_theta,
    selection = selection,
    seed = seed + 880001L
  )
  cache <- cache_local_evidence_corrections(
    factor_set = factor_set,
    theta = selected$theta,
    data_list = data_list,
    loglik_fn = loglik_fn,
    local_ids = local_ids,
    theta_metadata = selected$metadata,
    n_replicates = n_replicates,
    M = M,
    local_control = local_control,
    uncertainty_floor = uncertainty_floor,
    n_cores = n_cores,
    seed = seed + 990001L,
    stop_on_error = stop_on_error,
    verbose = verbose
  )
  correction <- estimate_theta_correction(
    cache,
    corrected_log_evidence = fit$log_evidence,
    require_complete = TRUE
  )
  account <- build_local_evidence_accounting(
    fit = fit,
    correction = correction,
    uncertainty = fit$numerical_uncertainty,
    model_comparison_scale = model_comparison_scale
  )
  structure(
    list(
      fit = fit,
      selected_theta = selected,
      correction_cache = cache,
      correction = correction,
      evidence_account = account,
      settings = list(
        n_theta = as.integer(n_theta),
        selection = selection,
        local_ids = local_ids,
        n_replicates = as.integer(n_replicates),
        M = as.integer(M),
        model_comparison_scale = as.numeric(model_comparison_scale),
        seed = as.integer(seed)
      )
    ),
    class = "local_evidence_delayed_acceptance"
  )
}

local_atlas_frozen_outer_rerun <- function(factor_set,
                                           N = 1000L,
                                           initial_proposal = NULL,
                                           n_mcmc_moves = 3L,
                                           min_mcmc_moves = 1L,
                                           max_rounds = 80L,
                                           n_cores = 1L,
                                           seed = 123L,
                                           verbose = FALSE) {
  if (!exists("outer_population_smc", mode = "function")) {
    stop("outer_population_smc() is unavailable; source outer_population_smc.R.")
  }
  outer_population_smc(
    factor_set = factor_set,
    N = as.integer(N),
    initial_proposal = initial_proposal,
    n_mcmc_moves = as.integer(n_mcmc_moves),
    min_mcmc_moves = as.integer(min_mcmc_moves),
    max_rounds = as.integer(max_rounds),
    n_cores = as.integer(n_cores),
    seed = as.integer(seed),
    verbose = isTRUE(verbose)
  )
}

local_atlas_fresh_endpoint_probe <- function(factor_set,
                                             theta,
                                             data_list,
                                             loglik_fn,
                                             local_ids = seq_along(factor_set$atlases),
                                             M = 500L,
                                             target_cess = 0.9,
                                             n_mcmc_moves = 2L,
                                             max_steps = 128L,
                                             n_cores = 1L,
                                             seed = 123L,
                                             verbose = FALSE) {
  factor_set <- validate_local_atlas_factor_set(factor_set)
  model <- factor_set$population_model
  theta <- .as_hyper_matrix(theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  local_ids <- as.integer(local_ids)
  if (any(local_ids < 1L | local_ids > length(factor_set$atlases))) {
    stop("local_ids are outside the atlas factor set.")
  }
  rows <- list()
  for (local_pos in local_ids) {
    atlas <- factor_set$atlases[[local_pos]]
    data_i <- data_list[[local_pos]]
    for (j in seq_len(nrow(theta))) {
      can_batch <- isTRUE(factor_set$evaluator_control$use_particle_mis) &&
        isTRUE(factor_set$evaluator_control$particle_mis_batch) &&
        !is.finite(factor_set$evaluator_control$max_leave_chart_out_gap) &&
        identical(factor_set$evaluator_control$particle_mis_role, "estimator")
      if (can_batch) {
        atlas_eval_row <- .local_atlas_particle_mis_many_global(
          atlas = atlas,
          theta = theta[j, , drop = FALSE],
          population_model = model,
          max_chart_distance = factor_set$evaluator_control$max_chart_distance,
          min_covering_charts = factor_set$evaluator_control$min_covering_charts,
          min_ess_frac = factor_set$evaluator_control$min_particle_mis_ess,
          min_ess = factor_set$evaluator_control$min_particle_mis_ess_abs,
          max_psis_k = factor_set$evaluator_control$max_particle_mis_psis_k,
          sparse_chart_min_covering = factor_set$evaluator_control$sparse_chart_min_covering,
          sparse_chart_max_distance = factor_set$evaluator_control$sparse_chart_max_distance,
          distance_metric = factor_set$evaluator_control$distance_metric,
          se_floor = factor_set$evaluator_control$se_floor,
          use_uncertified_estimates = factor_set$evaluator_control$use_uncertified_estimates
        )
        atlas_eval <- list(
          log_marginal = atlas_eval_row$log_marginal[1L],
          se = atlas_eval_row$se[1L],
          status = atlas_eval_row$status[1L],
          reason = atlas_eval_row$reason[1L]
        )
      } else {
        atlas_eval <- evaluate_local_atlas(
          atlas,
          theta = theta[j, , drop = FALSE],
          population_model = model,
          max_chart_distance = factor_set$evaluator_control$max_chart_distance,
          min_covering_charts = factor_set$evaluator_control$min_covering_charts,
          max_prediction_range = factor_set$evaluator_control$max_prediction_range,
          distance_scale = factor_set$evaluator_control$distance_scale,
          se_floor = factor_set$evaluator_control$se_floor,
          use_particle_mis = factor_set$evaluator_control$use_particle_mis,
          require_particle_mis = factor_set$evaluator_control$require_particle_mis,
          min_particle_mis_ess = factor_set$evaluator_control$min_particle_mis_ess,
          min_particle_mis_ess_abs = factor_set$evaluator_control$min_particle_mis_ess_abs,
          max_particle_mis_psis_k = factor_set$evaluator_control$max_particle_mis_psis_k,
          max_quadratic_particle_gap = factor_set$evaluator_control$max_quadratic_particle_gap,
          sparse_chart_min_covering = factor_set$evaluator_control$sparse_chart_min_covering,
          sparse_chart_max_distance = factor_set$evaluator_control$sparse_chart_max_distance,
          max_leave_chart_out_gap = factor_set$evaluator_control$max_leave_chart_out_gap,
          distance_metric = factor_set$evaluator_control$distance_metric,
          surface_method = factor_set$evaluator_control$surface_method,
          min_surface_charts = factor_set$evaluator_control$min_surface_charts,
          max_surface_se = factor_set$evaluator_control$max_surface_se,
          surface_value_nugget = factor_set$evaluator_control$surface_value_nugget,
          surface_gradient_weight = factor_set$evaluator_control$surface_gradient_weight,
          surface_curvature_weight = factor_set$evaluator_control$surface_curvature_weight,
          surface_ridge = factor_set$evaluator_control$surface_ridge,
          particle_mis_role = factor_set$evaluator_control$particle_mis_role
        )
      }
      fresh <- .local_chart_run_smc(
        local_id = atlas$local_id,
        theta_anchor = theta[j, , drop = FALSE],
        data_i = data_i,
        loglik_fn = loglik_fn,
        population_model = model,
        M = as.integer(M),
        target_cess = target_cess,
        resample_threshold = 0.5,
        n_mcmc_moves = n_mcmc_moves,
        rw_scale = 0.75,
        G_mix = 8L,
        da_enable = TRUE,
        refit_every = 2L,
        max_steps = max_steps,
        deterministic_resampling = FALSE,
        n_cores = n_cores,
        seed = as.integer(seed) + 1009L * local_pos + 9176L * j,
        verbose = verbose,
        source = "local_atlas_fresh_endpoint_probe"
      )
      rows[[length(rows) + 1L]] <- data.frame(
        local = names(factor_set$atlases)[local_pos],
        theta_row = j,
        atlas_log_marginal = atlas_eval$log_marginal,
        atlas_se = atlas_eval$se,
        atlas_status = atlas_eval$status,
        atlas_reason = atlas_eval$reason,
        fresh_log_marginal = fresh$logZ,
        fresh_se = fresh$logZ_se,
        delta_fresh_minus_atlas = if (is.finite(atlas_eval$log_marginal)) {
          fresh$logZ - atlas_eval$log_marginal
        } else {
          NA_real_
        },
        check.names = FALSE
      )
    }
  }
  if (length(rows)) do.call(rbind, rows) else data.frame()
}

local_atlas_repair_certification_pairs <- function(factor_set,
                                                   theta,
                                                   data_list,
                                                   loglik_fn,
                                                   candidate_pairs,
                                                   theta_weights = NULL,
                                                   local_ids = seq_along(factor_set$atlases),
                                                   M = 500L,
                                                   target_cess = 0.9,
                                                   n_mcmc_moves = 2L,
                                                   max_steps = 128L,
                                                   max_updates = nrow(candidate_pairs),
                                                   direct_confirmation_reps = 0L,
                                                   direct_confirmation_M = M,
                                                   direct_confirmation_max_sd = 1.5,
                                                   replicate_bootstrap_B = 200L,
                                                   max_direct_graph_z = 3,
                                                   max_direct_graph_chart_shift = 0.35,
                                                   max_direct_graph_existing_shift = 0.15,
                                                   local_control = list(),
                                                   edge_control = list(),
                                                   seed = 123L,
                                                   verbose = FALSE) {
  factor_set <- validate_local_atlas_factor_set(factor_set)
  model <- factor_set$population_model
  theta <- .as_hyper_matrix(theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  theta_weights <- .local_chart_normalize_weights(theta_weights, nrow(theta))
  if (!is.list(data_list) || length(data_list) < length(factor_set$atlases)) {
    stop("data_list must contain every local atlas selected for repair.")
  }
  if (!is.function(loglik_fn)) {
    stop("loglik_fn must be a function.")
  }
  if (is.null(candidate_pairs) || !nrow(as.data.frame(candidate_pairs))) {
    return(list(
      factor_set = factor_set,
      probes = data.frame(),
      selected = data.frame(),
      n_selected = 0L,
      n_activated = 0L,
      repair_executor = "local_atlas_repair_certification_pairs",
      graph_summary = local_atlas_graph_summary(factor_set)
    ))
  }
  if (is.character(local_ids)) {
    local_ids <- match(local_ids, names(factor_set$atlases))
  }
  local_ids <- as.integer(local_ids)
  local_ids <- local_ids[is.finite(local_ids)]
  if (!length(local_ids) || any(local_ids < 1L | local_ids > length(factor_set$atlases))) {
    stop("local_ids must identify atlases in factor_set.")
  }
  pairs <- as.data.frame(candidate_pairs, stringsAsFactors = FALSE)
  if (!"local_pos" %in% names(pairs)) {
    if (!"local" %in% names(pairs)) {
      stop("candidate_pairs must contain local_pos or local.")
    }
    pairs$local_pos <- match(as.character(pairs$local), names(factor_set$atlases))
  }
  if (!"theta_row" %in% names(pairs)) {
    stop("candidate_pairs must contain theta_row.")
  }
  pairs <- data.frame(
    local_pos = as.integer(pairs$local_pos),
    theta_row = as.integer(pairs$theta_row),
    check.names = FALSE
  )
  pairs <- pairs[
    is.finite(pairs$local_pos) &
      is.finite(pairs$theta_row) &
      pairs$local_pos %in% local_ids &
      pairs$theta_row >= 1L &
      pairs$theta_row <= nrow(theta),
    ,
    drop = FALSE
  ]
  pairs <- unique(pairs)
  max_updates <- as.integer(max_updates)
  if (!is.finite(max_updates) || max_updates < 0L) {
    stop("max_updates must be a non-negative integer.")
  }
  if (!nrow(pairs) || max_updates == 0L) {
    return(list(
      factor_set = factor_set,
      probes = data.frame(),
      selected = data.frame(),
      n_selected = 0L,
      n_activated = 0L,
      repair_executor = "local_atlas_repair_certification_pairs",
      graph_summary = local_atlas_graph_summary(factor_set)
    ))
  }
  pairs <- pairs[seq_len(min(nrow(pairs), max_updates)), , drop = FALSE]

  local_control <- .local_atlas_merge_control(local_control, .local_atlas_default_local_control())
  edge_control <- .local_atlas_merge_control(edge_control, .local_atlas_default_edge_control())
  local_control$candidate_M <- as.integer(M)
  local_control$n_cores_per_local <- 1L
  direct_confirmation_reps <- max(0L, as.integer(direct_confirmation_reps))
  direct_confirmation_M <- as.integer(direct_confirmation_M)
  if (!is.finite(direct_confirmation_M) || direct_confirmation_M <= 0L) {
    direct_confirmation_M <- as.integer(M)
  }

  run_smc <- function(local_pos, theta_row, rep_id, M_run, source) {
    atlas <- factor_set$atlases[[local_pos]]
    .local_chart_run_smc(
      local_id = atlas$local_id,
      theta_anchor = theta[theta_row, , drop = FALSE],
      data_i = data_list[[local_pos]],
      loglik_fn = loglik_fn,
      population_model = model,
      M = as.integer(M_run),
      target_cess = target_cess,
      resample_threshold = local_control$resample_threshold,
      n_mcmc_moves = as.integer(n_mcmc_moves),
      rw_scale = local_control$rw_scale,
      G_mix = as.integer(local_control$G_mix),
      da_enable = isTRUE(local_control$da_enable),
      refit_every = as.integer(local_control$refit_every),
      max_steps = as.integer(max_steps),
      deterministic_resampling = isTRUE(local_control$deterministic_resampling),
      n_cores = 1L,
      seed = as.integer(seed) + 1009L * local_pos + 9176L * theta_row + 104729L * as.integer(rep_id),
      verbose = verbose,
      source = source
    )
  }

  atlases <- factor_set$atlases
  rows <- vector("list", nrow(pairs))
  for (k in seq_len(nrow(pairs))) {
    local_pos <- pairs$local_pos[k]
    theta_row <- pairs$theta_row[k]
    atlas_before <- atlases[[local_pos]]
    graph_before <- atlas_before$diagnostics$normalizer_solution %||% list()
    primary <- run_smc(
      local_pos = local_pos,
      theta_row = theta_row,
      rep_id = 0L,
      M_run = M,
      source = "local_atlas_certification_repair"
    )
    repaired_run <- primary
    replicated <- NULL
    if (direct_confirmation_reps > 0L) {
      confirmation_runs <- lapply(seq_len(direct_confirmation_reps), function(rep_id) {
        run_smc(
          local_pos = local_pos,
          theta_row = theta_row,
          rep_id = rep_id,
          M_run = direct_confirmation_M,
          source = "local_atlas_certification_repair_confirmation"
        )
      })
      replicated <- .local_chart_replicated_logz(
        c(list(primary), confirmation_runs),
        bootstrap_B = as.integer(replicate_bootstrap_B),
        seed = as.integer(seed) + 900001L + 1009L * local_pos + 9176L * theta_row
      )
      repaired_run <- replicated$run
      rm(confirmation_runs)
    }

    chart_base <- sprintf("cert_t%03d", theta_row)
    inserted <- .local_atlas_insert_smc_chart(
      atlas = atlas_before,
      run = repaired_run,
      chart_id = chart_base,
      theta_anchor = theta[theta_row, , drop = FALSE],
      population_model = model,
      status_reason = "raw_certification_repair",
      source_run_id = sprintf("raw_certification_repair_%s_%03d", names(atlases)[local_pos], theta_row),
      direct_observation_role = "support",
      normalizer_certified = FALSE,
      normalizer_certification_method = "support_probe",
      normalizer_certification_reason = "support_chart_requires_certified_relative_edge"
    )
    attempt <- .local_atlas_activate_existing_chart(
      atlas = inserted$atlas,
      target_id = inserted$chart_id,
      data_i = data_list[[local_pos]],
      loglik_fn = loglik_fn,
      population_model = model,
      local_id = atlas_before$local_id,
      local_control = local_control,
      edge_control = edge_control,
      seed = as.integer(seed) + 7000003L + 7919L * k,
      verbose = verbose
    )

    direct_graph_certified <- NA
    direct_graph_reason <- NA_character_
    direct_graph_z <- NA_real_
    direct_graph_chart_shift <- NA_real_
    direct_graph_existing_shift <- NA_real_
    direct_confirmation_passed <- if (is.null(replicated)) NA else {
      is.finite(replicated$empirical_sd) &&
        replicated$empirical_sd <= as.numeric(direct_confirmation_max_sd)
    }
    if (isTRUE(attempt$success) && isTRUE(direct_confirmation_passed)) {
      direct_candidate <- attempt$atlas
      chart <- .local_atlas_chart_from_smc_run(
        run = repaired_run,
        local_id = direct_candidate$local_id,
        chart_id = inserted$chart_id,
        theta_anchor = theta[theta_row, , drop = FALSE],
        population_model = model,
        status = "active",
        status_reason = "raw_certification_repair",
        source_run_id = sprintf("raw_certification_repair_%s_%03d", names(atlases)[local_pos], theta_row),
        direct_observation_role = "calibration"
      )
      chart <- .local_chart_set_normalizer_certification(
        chart,
        certified = TRUE,
        method = "replicated_smc",
        reason = "replicated_direct_candidate"
      )
      direct_candidate$charts[[inserted$chart_id]] <- validate_local_chart(chart, model)
      direct_candidate <- tryCatch(
        solve_atlas_normalizers(direct_candidate, require_connected = TRUE),
        error = function(e) e
      )
      if (!inherits(direct_candidate, "error")) {
        direct_cert <- .local_atlas_certify_direct_normalizer(
          atlas = direct_candidate,
          chart_id = inserted$chart_id,
          population_model = model,
          max_direct_z = max_direct_graph_z,
          max_chart_shift = max_direct_graph_chart_shift,
          max_existing_shift = max_direct_graph_existing_shift
        )
        attempt$atlas <- direct_cert$atlas
        direct_graph_certified <- isTRUE(direct_cert$certified)
        direct_graph_reason <- direct_cert$reason
        direct_graph_z <- as.numeric(direct_cert$diagnostics$direct_z %||% NA_real_)
        direct_graph_chart_shift <- as.numeric(direct_cert$diagnostics$chart_shift %||% NA_real_)
        direct_graph_existing_shift <- as.numeric(direct_cert$diagnostics$max_existing_shift %||% NA_real_)
      } else {
        direct_graph_certified <- FALSE
        direct_graph_reason <- "replicated_direct_graph_solve_failed"
      }
    }

    atlases[[local_pos]] <- attempt$atlas
    graph_after <- attempt$atlas$diagnostics$normalizer_solution %||% list()
    shift <- .local_atlas_normalizer_shift(
      before_atlas = atlas_before,
      after_atlas = attempt$atlas,
      added_chart_id = inserted$chart_id
    )
    chart_after <- attempt$atlas$charts[[inserted$chart_id]]
    rows[[k]] <- data.frame(
      local = names(atlases)[local_pos],
      local_pos = local_pos,
      theta_row = theta_row,
      theta_weight = theta_weights[theta_row],
      repair_executor = "local_atlas_repair_certification_pairs",
      selected = TRUE,
      fresh_probed = TRUE,
      fresh_log_marginal = primary$logZ,
      fresh_se = primary$logZ_se,
      confirmed = !is.null(replicated),
      confirmation_passed = direct_confirmation_passed,
      confirmed_log_marginal = if (is.null(replicated)) NA_real_ else replicated$center,
      confirmed_se = if (is.null(replicated)) NA_real_ else replicated$combined_se,
      confirmation_sd = if (is.null(replicated)) NA_real_ else replicated$empirical_sd,
      confirmation_delta_se = if (is.null(replicated)) NA_real_ else replicated$delta_se,
      confirmation_bootstrap_se = if (is.null(replicated)) NA_real_ else replicated$bootstrap_se,
      graph_edge_z_before = as.numeric(graph_before$max_abs_standardized_edge_residual %||% 0),
      graph_direct_z_before = as.numeric(graph_before$max_abs_standardized_direct_residual %||% 0),
      graph_existing_shift = as.numeric(shift$max_existing_shift),
      graph_edge_z = as.numeric(graph_after$max_abs_standardized_edge_residual %||% NA_real_),
      graph_direct_z = as.numeric(graph_after$max_abs_standardized_direct_residual %||% NA_real_),
      activation_success = isTRUE(attempt$success),
      activation_reason = attempt$reason,
      chart_id = inserted$chart_id,
      chart_status = as.character(chart_after$status %||% "missing"),
      normalizer_certified = if (is.null(chart_after)) FALSE else .local_chart_normalizer_certified(chart_after),
      normalizer_method = as.character((chart_after$diagnostics$normalizer_certification %||% list())$method %||% NA_character_),
      normalizer_reason = as.character((chart_after$diagnostics$normalizer_certification %||% list())$reason %||% NA_character_),
      direct_graph_certified = direct_graph_certified,
      direct_graph_reason = direct_graph_reason,
      direct_graph_z = direct_graph_z,
      direct_graph_chart_shift = direct_graph_chart_shift,
      direct_graph_existing_shift = direct_graph_existing_shift,
      check.names = FALSE
    )
    rm(primary, repaired_run, replicated)
    invisible(gc(FALSE))
  }

  updated_factor_set <- build_local_atlas_factor_set(
    atlases = atlases,
    population_model = model,
    max_chart_distance = factor_set$evaluator_control$max_chart_distance,
    min_covering_charts = factor_set$evaluator_control$min_covering_charts,
    max_prediction_range = factor_set$evaluator_control$max_prediction_range,
    distance_scale = factor_set$evaluator_control$distance_scale,
    se_floor = factor_set$evaluator_control$se_floor,
    use_particle_mis = factor_set$evaluator_control$use_particle_mis,
    require_particle_mis = factor_set$evaluator_control$require_particle_mis,
    min_particle_mis_ess = factor_set$evaluator_control$min_particle_mis_ess,
    min_particle_mis_ess_abs = factor_set$evaluator_control$min_particle_mis_ess_abs,
    max_particle_mis_psis_k = factor_set$evaluator_control$max_particle_mis_psis_k,
    max_quadratic_particle_gap = factor_set$evaluator_control$max_quadratic_particle_gap,
    sparse_chart_min_covering = factor_set$evaluator_control$sparse_chart_min_covering,
    sparse_chart_max_distance = factor_set$evaluator_control$sparse_chart_max_distance,
    max_leave_chart_out_gap = factor_set$evaluator_control$max_leave_chart_out_gap,
    distance_metric = factor_set$evaluator_control$distance_metric,
    surface_method = factor_set$evaluator_control$surface_method,
    min_surface_charts = factor_set$evaluator_control$min_surface_charts,
    max_surface_se = factor_set$evaluator_control$max_surface_se,
    surface_value_nugget = factor_set$evaluator_control$surface_value_nugget,
    surface_gradient_weight = factor_set$evaluator_control$surface_gradient_weight,
    surface_curvature_weight = factor_set$evaluator_control$surface_curvature_weight,
    surface_ridge = factor_set$evaluator_control$surface_ridge,
    particle_mis_role = factor_set$evaluator_control$particle_mis_role,
    particle_mis_batch = factor_set$evaluator_control$particle_mis_batch,
    stop_on_uncertified = factor_set$evaluator_control$stop_on_uncertified,
    use_uncertified_estimates = factor_set$evaluator_control$use_uncertified_estimates
  )
  probes <- if (length(rows)) do.call(rbind, rows) else data.frame()
  list(
    factor_set = updated_factor_set,
    probes = probes,
    selected = probes,
    n_selected = nrow(probes),
    n_activated = sum(probes$activation_success %in% TRUE),
    repair_executor = "local_atlas_repair_certification_pairs",
    graph_summary = local_atlas_graph_summary(updated_factor_set)
  )
}

local_atlas_benchmark_gate <- function(factor_set,
                                       reference_draws,
                                       workflow_draws = NULL,
                                       workflow_fit = NULL,
                                       baseline_draws = NULL,
                                       baseline_comparison = NULL,
                                       focus_parameters = NULL,
                                       theta_audit = NULL,
                                       theta_audit_weights = NULL,
                                       frozen_outer_control = NULL,
                                       fresh_probe_control = NULL,
                                       thresholds = list(),
                                       n_cores = 1L,
                                       seed = 123L) {
  factor_set <- validate_local_atlas_factor_set(factor_set)
  thresholds <- modifyList(
    list(
      max_abs_standardized_mean_error = 0.50,
      mean_abs_standardized_mean_error = 0.20,
      max_uncertified_fraction = 0,
      max_uncertified_weight = 0,
      max_graph_standardized_residual = 3,
      max_root_logZ_se = Inf,
      max_chart_logZ_abs_se = Inf,
      max_quarantine_invariance_delta = 1e-8,
      max_quarantine_invariance_weighted_delta = 1e-8,
      max_fresh_probe_abs_delta = Inf,
      max_fresh_probe_abs_z = 3,
      max_frozen_outer_mean_shift = 0.25,
      min_baseline_mean_error_improvement = 0,
      min_baseline_shape_error_improvement = 0,
      max_focus_parameter_worsening = 0,
      require_baseline_improvement = TRUE,
      require_theta_audit = TRUE,
      require_frozen_outer = TRUE,
      require_fresh_probes = TRUE
    ),
    thresholds
  )
  failures <- character()

  if (is.null(workflow_draws)) {
    if (is.null(workflow_fit)) {
      stop("local_atlas_benchmark_gate requires workflow_draws or workflow_fit.")
    }
    workflow_draws <- local_atlas_draws_from_fit(
      workflow_fit,
      population_model = factor_set$population_model,
      n_draws = nrow(reference_draws),
      seed = seed
    )
  }
  posterior_comparison <- local_atlas_compare_posterior_draws(reference_draws, workflow_draws)
  posterior_summary <- .local_atlas_metric_summary(posterior_comparison)
  focus_summary <- if (is.null(focus_parameters)) NULL else {
    .local_atlas_metric_summary(posterior_comparison, parameters = focus_parameters)
  }
  if (posterior_summary$max_abs_standardized_mean_error > thresholds$max_abs_standardized_mean_error) {
    failures <- c(failures, "posterior_max_standardized_mean_error")
  }
  if (posterior_summary$mean_abs_standardized_mean_error > thresholds$mean_abs_standardized_mean_error) {
    failures <- c(failures, "posterior_mean_standardized_mean_error")
  }

  baseline_improvement <- NULL
  baseline_summary <- NULL
  focus_baseline_summary <- NULL
  baseline_improvement_summary <- NULL
  if (!is.null(baseline_draws) || !is.null(baseline_comparison)) {
    baseline_improvement <- local_atlas_compare_to_baseline(
      reference_draws = reference_draws,
      workflow_draws = workflow_draws,
      baseline_draws = baseline_draws,
      baseline_comparison = baseline_comparison
    )
    baseline_comparison_resolved <- if (is.null(baseline_comparison)) {
      local_atlas_compare_posterior_draws(reference_draws, baseline_draws)
    } else {
      baseline_comparison
    }
    baseline_summary <- .local_atlas_metric_summary(baseline_comparison_resolved)
    improvement_parameters <- focus_parameters %||% baseline_improvement$parameter
    improvement_rows <- baseline_improvement$parameter %in% improvement_parameters
    if (!any(improvement_rows)) {
      failures <- c(failures, "baseline_improvement_parameters_missing")
    } else {
      if (!is.null(focus_parameters)) {
        focus_baseline_summary <- .local_atlas_metric_summary(
          baseline_comparison_resolved,
          parameters = focus_parameters
        )
      }
      mean_improvement <- mean(
        baseline_improvement$standardized_mean_error_improvement[improvement_rows],
        na.rm = TRUE
      )
      shape_improvement <- mean(
        baseline_improvement$shape_error_improvement[improvement_rows],
        na.rm = TRUE
      )
      worst_worsening <- min(
        baseline_improvement$shape_error_improvement[improvement_rows],
        na.rm = TRUE
      )
      baseline_improvement_summary <- data.frame(
        parameters = paste(baseline_improvement$parameter[improvement_rows], collapse = ","),
        mean_standardized_mean_error_improvement = mean_improvement,
        mean_shape_error_improvement = shape_improvement,
        worst_shape_error_improvement = worst_worsening,
        check.names = FALSE
      )
      if (mean_improvement < thresholds$min_baseline_mean_error_improvement) {
        failures <- c(failures, "baseline_mean_error_not_improved")
      }
      if (shape_improvement < thresholds$min_baseline_shape_error_improvement) {
        failures <- c(failures, "baseline_shape_not_improved")
      }
      if (worst_worsening < -thresholds$max_focus_parameter_worsening) {
        failures <- c(failures, "baseline_parameter_worsened")
      }
    }
  } else if (isTRUE(thresholds$require_baseline_improvement)) {
    failures <- c(failures, "missing_baseline_comparison")
  }

  if (is.null(theta_audit) && !is.null(workflow_fit)) {
    theta_audit <- workflow_fit$theta
  }
  if (is.null(theta_audit_weights) && !is.null(workflow_fit)) {
    theta_audit_weights <- workflow_fit$w
  }
  certification <- NULL
  quarantine_invariance <- NULL
  if (!is.null(theta_audit)) {
    theta_audit <- .as_hyper_matrix(
      theta_audit,
      hyper_names = factor_set$population_model$hyper_names,
      hyper_dim = factor_set$population_model$hyper_dim
    )
    theta_audit_weights <- .local_chart_normalize_weights(theta_audit_weights, nrow(theta_audit))
    certification <- local_atlas_certification_summary(
      factor_set,
      theta_audit,
      theta_weights = theta_audit_weights,
      n_cores = n_cores
    )
    max_uncertified <- max(certification$uncertified_fraction, na.rm = TRUE)
    max_uncertified_weight <- max(certification$uncertified_weight, na.rm = TRUE)
    if (max_uncertified > thresholds$max_uncertified_fraction) {
      failures <- c(failures, "uncertified_theta_coverage")
    }
    if (max_uncertified_weight > thresholds$max_uncertified_weight) {
      failures <- c(failures, "uncertified_theta_weight")
    }
    if (is.finite(max_uncertified) && max_uncertified == 0 &&
        is.finite(max_uncertified_weight) && max_uncertified_weight == 0) {
      quarantine_invariance <- local_atlas_quarantine_invariance(
        factor_set,
        theta_audit,
        theta_weights = theta_audit_weights,
        n_cores = n_cores
      )
      if (quarantine_invariance$max_abs_delta > thresholds$max_quarantine_invariance_delta) {
        failures <- c(failures, "quarantine_invariance")
      }
      if (quarantine_invariance$weighted_mean_abs_delta >
          thresholds$max_quarantine_invariance_weighted_delta) {
        failures <- c(failures, "weighted_quarantine_invariance")
      }
    } else {
      quarantine_invariance <- data.frame(
        n_theta = nrow(theta_audit),
        max_abs_delta = NA_real_,
        mean_abs_delta = NA_real_,
        weighted_mean_abs_delta = NA_real_,
        skipped_reason = "uncertified_theta_coverage",
        check.names = FALSE
      )
    }
  } else if (isTRUE(thresholds$require_theta_audit)) {
    failures <- c(failures, "missing_theta_audit")
  }

  graph <- local_atlas_graph_summary(factor_set)
  if (max(graph$max_abs_standardized_edge_residual, na.rm = TRUE) >
      thresholds$max_graph_standardized_residual) {
    failures <- c(failures, "graph_residual")
  }
  if (max(graph$root_logZ_se, na.rm = TRUE) > thresholds$max_root_logZ_se) {
    failures <- c(failures, "root_logZ_uncertainty")
  }
  if (max(graph$max_chart_logZ_abs_se, na.rm = TRUE) > thresholds$max_chart_logZ_abs_se) {
    failures <- c(failures, "chart_logZ_uncertainty")
  }

  frozen_outer <- NULL
  frozen_outer_comparison <- NULL
  frozen_outer_error <- NULL
  if (!is.null(frozen_outer_control)) {
    frozen_outer <- tryCatch(
      do.call(
        local_atlas_frozen_outer_rerun,
        modifyList(list(factor_set = factor_set, seed = seed, n_cores = n_cores), frozen_outer_control)
      ),
      error = function(e) e
    )
    if (inherits(frozen_outer, "error")) {
      frozen_outer_error <- conditionMessage(frozen_outer)
      frozen_outer <- NULL
      failures <- c(failures, "frozen_outer_error")
    } else {
      frozen_draws <- local_atlas_draws_from_fit(
        frozen_outer,
        population_model = factor_set$population_model,
        n_draws = nrow(workflow_draws),
        seed = seed + 1L
      )
      frozen_outer_comparison <- local_atlas_compare_posterior_draws(workflow_draws, frozen_draws)
      max_shift <- max(abs(frozen_outer_comparison$standardized_mean_error), na.rm = TRUE)
      if (max_shift > thresholds$max_frozen_outer_mean_shift) {
        failures <- c(failures, "frozen_outer_instability")
      }
    }
  } else if (isTRUE(thresholds$require_frozen_outer)) {
    failures <- c(failures, "missing_frozen_outer_rerun")
  }

  fresh_probes <- NULL
  fresh_probe_error <- NULL
  if (!is.null(fresh_probe_control)) {
    fresh_probes <- tryCatch(
      do.call(
        local_atlas_fresh_endpoint_probe,
        modifyList(list(factor_set = factor_set, seed = seed, n_cores = n_cores), fresh_probe_control)
      ),
      error = function(e) e
    )
    if (inherits(fresh_probes, "error")) {
      fresh_probe_error <- conditionMessage(fresh_probes)
      fresh_probes <- NULL
      failures <- c(failures, "fresh_probe_error")
    } else if (!nrow(fresh_probes)) {
      failures <- c(failures, "empty_fresh_endpoint_probes")
      fresh_probes$abs_standardized_delta <- numeric()
    } else {
      max_probe_delta <- max(abs(fresh_probes$delta_fresh_minus_atlas), na.rm = TRUE)
      if (is.finite(max_probe_delta) && max_probe_delta > thresholds$max_fresh_probe_abs_delta) {
        failures <- c(failures, "fresh_probe_delta")
      }
      fresh_probe_z <- abs(fresh_probes$delta_fresh_minus_atlas) /
        pmax(
          sqrt(pmax(fresh_probes$fresh_se, 0)^2 + pmax(fresh_probes$atlas_se, 0)^2),
          .Machine$double.eps
        )
      if (any(fresh_probes$atlas_status != "certified")) {
        failures <- c(failures, "fresh_probe_uncertified_atlas")
      }
      if (max(fresh_probe_z, na.rm = TRUE) > thresholds$max_fresh_probe_abs_z) {
        failures <- c(failures, "fresh_probe_standardized_delta")
      }
      fresh_probes$abs_standardized_delta <- fresh_probe_z
    }
  } else if (isTRUE(thresholds$require_fresh_probes)) {
    failures <- c(failures, "missing_fresh_endpoint_probes")
  }

  structure(
    list(
      passed = !length(failures),
      failures = unique(failures),
      thresholds = thresholds,
      posterior_summary = posterior_summary,
      focus_summary = focus_summary,
      posterior_comparison = posterior_comparison,
      baseline_summary = baseline_summary,
      focus_baseline_summary = focus_baseline_summary,
      baseline_improvement = baseline_improvement,
      baseline_improvement_summary = baseline_improvement_summary,
      certification_summary = certification,
      graph_summary = graph,
      quarantine_invariance = quarantine_invariance,
      frozen_outer_fit = frozen_outer,
      frozen_outer_comparison = frozen_outer_comparison,
      frozen_outer_error = frozen_outer_error,
      fresh_probes = fresh_probes,
      fresh_probe_error = fresh_probe_error
    ),
    class = "local_atlas_benchmark_gate"
  )
}

local_evidence_fit_kernel_residual_patch <- function(theta,
                                                     base_log_m,
                                                     probe_theta_row,
                                                     probe_log_m,
                                                     probe_sd = NULL,
                                                     theta_weights = NULL,
                                                     active_hyper_names = NULL,
                                                     kernel_scale = 0.5,
                                                     scale_floor = 0.25,
                                                     shrink = 2,
                                                     min_train = 4L,
                                                     probe_sd_floor = 0.05) {
  theta <- as.data.frame(theta, check.names = FALSE)
  n_theta <- nrow(theta)
  if (!n_theta) {
    stop("theta must contain at least one row.")
  }
  if (length(base_log_m) != n_theta) {
    stop("base_log_m must have one entry per theta row.")
  }
  if (is.null(active_hyper_names)) {
    active_hyper_names <- names(theta)
  }
  active_hyper_names <- intersect(as.character(active_hyper_names), names(theta))
  if (!length(active_hyper_names)) {
    stop("active_hyper_names must select at least one theta column.")
  }
  probe_theta_row <- as.integer(probe_theta_row)
  probe_log_m <- as.numeric(probe_log_m)
  if (is.null(probe_sd)) {
    probe_sd <- rep(as.numeric(probe_sd_floor), length(probe_log_m))
  }
  probe_sd <- as.numeric(probe_sd)
  if (length(probe_theta_row) != length(probe_log_m) ||
      length(probe_theta_row) != length(probe_sd)) {
    stop("probe_theta_row, probe_log_m, and probe_sd must have the same length.")
  }
  if (is.null(theta_weights)) {
    theta_weights <- rep(1 / n_theta, n_theta)
  }
  theta_weights <- .local_chart_normalize_weights(theta_weights, n_theta)

  ok <- is.finite(probe_theta_row) &
    probe_theta_row >= 1L &
    probe_theta_row <= n_theta &
    is.finite(probe_log_m) &
    is.finite(base_log_m[probe_theta_row])
  probe_theta_row <- probe_theta_row[ok]
  probe_log_m <- probe_log_m[ok]
  probe_sd <- probe_sd[ok]
  if (length(probe_theta_row) < as.integer(min_train)) {
    return(structure(
      list(
        estimate = as.numeric(base_log_m),
        correction = rep(0, n_theta),
        diagnostics = data.frame(
          status = "too_few_patch_points_no_correction",
          n_train = length(probe_theta_row),
          offset = NA_real_,
          signal_sd = NA_real_,
          loo_rmse = NA_real_,
          reliability = 0,
          shape_scale = 0,
          kernel_scale = as.numeric(kernel_scale),
          scale_floor = as.numeric(scale_floor),
          active_hyper_names = paste(active_hyper_names, collapse = ","),
          check.names = FALSE
        )
      ),
      class = "local_evidence_kernel_residual_patch"
    ))
  }

  x <- as.matrix(theta[, active_hyper_names, drop = FALSE])
  center <- colSums(x * theta_weights)
  posterior_scale <- sqrt(colSums(sweep(x, 2L, center, "-")^2 * theta_weights))
  grid_scale <- apply(x, 2L, stats::sd)
  scale <- pmax(posterior_scale, as.numeric(scale_floor) * grid_scale, 1e-8)
  scale[!is.finite(scale) | scale <= 0] <- 1
  z_all <- sweep(sweep(x, 2L, center, "-"), 2L, scale, "/")
  z_train <- z_all[probe_theta_row, , drop = FALSE]

  delta <- probe_log_m - as.numeric(base_log_m[probe_theta_row])
  local_weight <- theta_weights[probe_theta_row]
  if (!any(is.finite(local_weight)) || sum(local_weight, na.rm = TRUE) <= 0) {
    local_weight <- rep(1, length(delta))
  }
  noise_var <- pmax(probe_sd, as.numeric(probe_sd_floor))^2
  smooth_weight <- local_weight / pmax(noise_var, .Machine$double.eps)
  if (!any(is.finite(smooth_weight)) || sum(smooth_weight, na.rm = TRUE) <= 0) {
    smooth_weight <- rep(1, length(delta))
  }
  smooth_weight <- smooth_weight / mean(smooth_weight[is.finite(smooth_weight)])

  weighted_mean <- function(y, w) {
    ok <- is.finite(y) & is.finite(w) & w >= 0
    if (!any(ok)) return(NA_real_)
    sum(y[ok] * w[ok]) / sum(w[ok])
  }
  offset <- weighted_mean(delta, local_weight)
  if (!is.finite(offset)) {
    offset <- mean(delta)
  }
  y_shape <- delta - offset
  bandwidth <- max(as.numeric(kernel_scale), 1e-8)
  smooth_at <- function(z_query, keep = seq_along(y_shape)) {
    if (!length(keep)) {
      return(0)
    }
    d2 <- rowSums(sweep(z_train[keep, , drop = FALSE], 2L, z_query, "-")^2)
    kw <- exp(-0.5 * d2 / bandwidth^2) * smooth_weight[keep]
    if (!any(is.finite(kw)) || sum(kw, na.rm = TRUE) <= 0) {
      return(0)
    }
    sum(kw * y_shape[keep]) / sum(kw)
  }
  pred_train <- vapply(seq_along(y_shape), function(j) smooth_at(z_train[j, ]), numeric(1))
  residual <- y_shape - pred_train
  loo <- rep(NA_real_, length(y_shape))
  if (length(y_shape) > 1L) {
    for (j in seq_along(y_shape)) {
      loo[j] <- y_shape[j] - smooth_at(z_train[j, ], keep = setdiff(seq_along(y_shape), j))
    }
  }
  signal_sd <- sqrt(weighted_mean(y_shape^2, local_weight))
  loo_rmse <- sqrt(mean(loo[is.finite(loo)]^2))
  if (!is.finite(loo_rmse)) {
    loo_rmse <- sqrt(mean(residual^2))
  }
  reliability <- if (is.finite(signal_sd) && is.finite(loo_rmse)) {
    signal_sd^2 / (signal_sd^2 + loo_rmse^2 + 1e-8)
  } else {
    0
  }
  shape_scale <- max(0, min(1, as.numeric(shrink) * reliability))
  pred_all <- vapply(seq_len(n_theta), function(j) smooth_at(z_all[j, ]), numeric(1))
  pred_center <- weighted_mean(pred_all, theta_weights)
  if (!is.finite(pred_center)) {
    pred_center <- mean(pred_all[is.finite(pred_all)])
  }
  dist <- vapply(seq_len(n_theta), function(i) {
    sqrt(min(rowSums(sweep(z_train, 2L, z_all[i, ], "-")^2)))
  }, numeric(1))
  coverage <- exp(-0.5 * dist^2 / bandwidth^2)
  correction <- offset + shape_scale * (pred_all - pred_center)

  structure(
    list(
      estimate = as.numeric(base_log_m) + correction,
      correction = correction,
      train = data.frame(
        theta_row = probe_theta_row,
        base_log_m = as.numeric(base_log_m[probe_theta_row]),
        probe_log_m = probe_log_m,
        probe_sd = probe_sd,
        delta = delta,
        fitted_shape_delta = pred_train,
        loo_error = loo,
        check.names = FALSE
      ),
      diagnostics = data.frame(
        status = "fitted",
        n_train = length(delta),
        offset = offset,
        signal_sd = signal_sd,
        loo_rmse = loo_rmse,
        training_rmse = sqrt(mean(residual^2)),
        reliability = reliability,
        shape_scale = shape_scale,
        kernel_scale = bandwidth,
        scale_floor = as.numeric(scale_floor),
        mean_coverage = mean(coverage),
        min_coverage = min(coverage),
        max_coverage = max(coverage),
        active_hyper_names = paste(active_hyper_names, collapse = ","),
        check.names = FALSE
      )
    ),
    class = "local_evidence_kernel_residual_patch"
  )
}

.local_evidence_posterior_parameter_family <- function(parameter) {
  parameter <- as.character(parameter)
  ifelse(
    grepl("^mu_", parameter),
    sub("^mu_", "", parameter),
    ifelse(
      grepl("^sigma2_", parameter),
      sub("^sigma2_", "", parameter),
      ifelse(
        grepl("^log_sigma2_", parameter),
        sub("^log_sigma2_", "", parameter),
        parameter
      )
    )
  )
}

.local_evidence_family_hyper_group <- function(family, theta_names) {
  family <- as.character(family)
  theta_names <- as.character(theta_names)
  group <- intersect(c(paste0("mu_", family), paste0("log_sigma2_", family)), theta_names)
  if (!length(group) && family %in% theta_names) {
    group <- family
  }
  group
}

local_evidence_detect_kernel_axis <- function(posterior_comparison,
                                              theta_names,
                                              mean_error_weight = 1,
                                              wasserstein_weight = 0.25,
                                              shape_weight = 0.1,
                                              min_abs_standardized_mean_error = 0) {
  if (!is.data.frame(posterior_comparison) || !nrow(posterior_comparison)) {
    stop("posterior_comparison must be a non-empty data frame.")
  }
  if (!"parameter" %in% names(posterior_comparison)) {
    stop("posterior_comparison must contain a parameter column.")
  }
  if (!"standardized_mean_error" %in% names(posterior_comparison)) {
    stop("posterior_comparison must contain a standardized_mean_error column.")
  }
  theta_names <- as.character(theta_names)
  comparison <- posterior_comparison
  comparison$family <- .local_evidence_posterior_parameter_family(comparison$parameter)
  comparison$abs_standardized_mean_error <- abs(as.numeric(comparison$standardized_mean_error))
  comparison$scaled_q_wasserstein <- if ("scaled_q_wasserstein" %in% names(comparison)) {
    as.numeric(comparison$scaled_q_wasserstein)
  } else {
    0
  }
  comparison$shape_error <- if ("shape_error" %in% names(comparison)) {
    as.numeric(comparison$shape_error)
  } else {
    0
  }

  families <- unique(comparison$family)
  rows <- lapply(families, function(family) {
    idx <- comparison$family == family
    group <- .local_evidence_family_hyper_group(family, theta_names)
    if (!length(group)) return(NULL)
    data.frame(
      family = family,
      active_hyper_names = paste(group, collapse = "+"),
      dimension = length(group),
      max_abs_standardized_mean_error = max(comparison$abs_standardized_mean_error[idx], na.rm = TRUE),
      max_scaled_q_wasserstein = max(comparison$scaled_q_wasserstein[idx], na.rm = TRUE),
      max_shape_error = max(comparison$shape_error[idx], na.rm = TRUE),
      parameters = paste(comparison$parameter[idx], collapse = ","),
      check.names = FALSE
    )
  })
  table <- do.call(rbind, Filter(Negate(is.null), rows))
  if (!is.data.frame(table) || !nrow(table)) {
    stop("posterior_comparison did not map to any theta hyperparameter axis.")
  }
  table$max_abs_standardized_mean_error[!is.finite(table$max_abs_standardized_mean_error)] <- 0
  table$max_scaled_q_wasserstein[!is.finite(table$max_scaled_q_wasserstein)] <- 0
  table$max_shape_error[!is.finite(table$max_shape_error)] <- 0
  table$posterior_axis_score <-
    as.numeric(mean_error_weight) * table$max_abs_standardized_mean_error +
    as.numeric(wasserstein_weight) * table$max_scaled_q_wasserstein +
    as.numeric(shape_weight) * table$max_shape_error
  table <- table[order(
    -table$posterior_axis_score,
    -table$max_abs_standardized_mean_error,
    table$family
  ), , drop = FALSE]
  eligible <- table$max_abs_standardized_mean_error >=
    as.numeric(min_abs_standardized_mean_error)
  if (!any(eligible)) {
    eligible[1L] <- TRUE
  }
  selected <- table[which(eligible)[1L], , drop = FALSE]
  active_hyper_names <- strsplit(selected$active_hyper_names[1L], "+", fixed = TRUE)[[1L]]
  structure(
    list(
      active_hyper_names = active_hyper_names,
      family = selected$family[1L],
      selected = selected,
      table = table
    ),
    class = "local_evidence_posterior_kernel_axis"
  )
}

local_evidence_summarize_probe_replicates <- function(probes,
                                                      theta_row_col = "theta_row",
                                                      log_m_col = "log_marginal",
                                                      se_col = "path_se") {
  if (!is.data.frame(probes) || !nrow(probes)) {
    stop("probes must be a non-empty data frame.")
  }
  required <- c(theta_row_col, log_m_col)
  missing <- setdiff(required, names(probes))
  if (length(missing)) {
    stop("probe data is missing columns: ", paste(missing, collapse = ", "))
  }
  parts <- split(probes, probes[[theta_row_col]], drop = TRUE)
  do.call(rbind, lapply(parts, function(df) {
    z <- as.numeric(df[[log_m_col]])
    z <- z[is.finite(z)]
    se <- if (se_col %in% names(df)) as.numeric(df[[se_col]]) else NA_real_
    se <- se[is.finite(se)]
    data.frame(
      theta_row = as.integer(df[[theta_row_col]][1L]),
      probe_log_m = if (length(z)) logsumexp(z) - log(length(z)) else NA_real_,
      probe_sd = if (length(z) > 1L) stats::sd(z) else if (length(se)) mean(se) else NA_real_,
      probe_reps = length(z),
      probe_min = if (length(z)) min(z) else NA_real_,
      probe_max = if (length(z)) max(z) else NA_real_,
      check.names = FALSE
    )
  }))
}

local_evidence_tune_kernel_residual_patch <- function(theta,
                                                      base_log_m,
                                                      probe_theta_row,
                                                      probe_log_m,
                                                      probe_sd = NULL,
                                                      theta_weights = NULL,
                                                      active_hyper_names,
                                                      kernel_scale_grid = c(0.35, 0.5, 0.75),
                                                      scale_floor_grid = c(0.25, 1, 2),
                                                      shrink_grid = c(4),
                                                      min_train = 4L,
                                                      probe_sd_floor = 0.05,
                                                      selection_score = c("loo", "replicate_validation"),
                                                      probe_replicates = NULL,
                                                      replicate_col = "replicate",
                                                      log_m_col = "log_marginal",
                                                      se_col = "path_se") {
  theta <- as.data.frame(theta, check.names = FALSE)
  n_theta <- nrow(theta)
  selection_score <- match.arg(selection_score)
  if (is.null(theta_weights)) {
    theta_weights <- rep(1 / n_theta, n_theta)
  }
  theta_weights <- .local_chart_normalize_weights(theta_weights, n_theta)
  active_hyper_names <- intersect(as.character(active_hyper_names), names(theta))
  if (!length(active_hyper_names)) {
    stop("active_hyper_names must select at least one theta column.")
  }
  group <- active_hyper_names

  replicate_validation_score <- function(group, kernel_scale, scale_floor, shrink) {
    if (is.null(probe_replicates) ||
        !is.data.frame(probe_replicates) ||
        !nrow(probe_replicates) ||
        !replicate_col %in% names(probe_replicates)) {
      return(NA_real_)
    }
    reps <- sort(unique(probe_replicates[[replicate_col]]))
    if (length(reps) < 2L) {
      return(NA_real_)
    }
    err2 <- numeric()
    err_w <- numeric()
    for (rep_id in reps) {
      train_raw <- probe_replicates[probe_replicates[[replicate_col]] != rep_id, , drop = FALSE]
      valid_raw <- probe_replicates[probe_replicates[[replicate_col]] == rep_id, , drop = FALSE]
      train <- local_evidence_summarize_probe_replicates(
        train_raw,
        theta_row_col = "theta_row",
        log_m_col = log_m_col,
        se_col = se_col
      )
      valid <- local_evidence_summarize_probe_replicates(
        valid_raw,
        theta_row_col = "theta_row",
        log_m_col = log_m_col,
        se_col = se_col
      )
      fit <- local_evidence_fit_kernel_residual_patch(
        theta = theta,
        base_log_m = base_log_m,
        probe_theta_row = train$theta_row,
        probe_log_m = train$probe_log_m,
        probe_sd = train$probe_sd,
        theta_weights = theta_weights,
        active_hyper_names = group,
        kernel_scale = kernel_scale,
        scale_floor = scale_floor,
        shrink = shrink,
        min_train = min_train,
        probe_sd_floor = probe_sd_floor
      )
      pred <- fit$correction[valid$theta_row]
      truth <- valid$probe_log_m - as.numeric(base_log_m[valid$theta_row])
      ok <- is.finite(pred) & is.finite(truth)
      if (any(ok)) {
        err2 <- c(err2, (pred[ok] - truth[ok])^2)
        err_w <- c(err_w, theta_weights[valid$theta_row][ok])
      }
    }
    if (!length(err2)) {
      return(NA_real_)
    }
    sqrt(sum(err2 * pmax(err_w, .Machine$double.eps)) / sum(pmax(err_w, .Machine$double.eps)))
  }

  fits <- list()
  rows <- list()
  idx <- 0L
  group_key <- paste(group, collapse = "+")
  for (kernel_scale in as.numeric(kernel_scale_grid)) {
    for (scale_floor in as.numeric(scale_floor_grid)) {
      for (shrink in as.numeric(shrink_grid)) {
        idx <- idx + 1L
        fit <- local_evidence_fit_kernel_residual_patch(
          theta = theta,
          base_log_m = base_log_m,
          probe_theta_row = probe_theta_row,
          probe_log_m = probe_log_m,
          probe_sd = probe_sd,
          theta_weights = theta_weights,
          active_hyper_names = group,
          kernel_scale = kernel_scale,
          scale_floor = scale_floor,
          shrink = shrink,
          min_train = min_train,
          probe_sd_floor = probe_sd_floor
        )
        diag <- fit$diagnostics
        validation <- replicate_validation_score(group, kernel_scale, scale_floor, shrink)
        score <- if (identical(selection_score, "replicate_validation") &&
                     is.finite(validation)) {
          validation
        } else {
          as.numeric(diag$loo_rmse[1L])
        }
        fits[[idx]] <- fit
        rows[[idx]] <- data.frame(
          setting_id = idx,
          active_hyper_names = group_key,
          dimension = length(group),
          kernel_scale = kernel_scale,
          scale_floor = scale_floor,
          shrink = shrink,
          score = score,
          score_source = if (identical(selection_score, "replicate_validation") &&
                             is.finite(validation)) "replicate_validation" else "loo",
          replicate_validation_rmse = validation,
          loo_rmse = as.numeric(diag$loo_rmse[1L]),
          training_rmse = as.numeric(diag$training_rmse[1L]),
          signal_sd = as.numeric(diag$signal_sd[1L]),
          reliability = as.numeric(diag$reliability[1L]),
          shape_scale = as.numeric(diag$shape_scale[1L]),
          offset = as.numeric(diag$offset[1L]),
          status = as.character(diag$status[1L]),
          check.names = FALSE
        )
      }
    }
  }
  table <- do.call(rbind, rows)
  selectable <- is.finite(table$score) & table$status == "fitted"
  if (!any(selectable)) {
    best <- which.min(seq_len(nrow(table)))
  } else {
    candidates <- which(selectable)
    best <- candidates[which.min(table$score[candidates])]
  }
  out <- fits[[best]]
  out$selection_table <- table
  out$selected <- table[best, , drop = FALSE]
  class(out) <- c("local_evidence_tuned_kernel_residual_patch", class(out))
  out
}

local_chart_score <- function(theta, alpha, weights = NULL, population_model) {
  parts <- .local_chart_diag_gaussian_parts(
    theta = theta,
    alpha = alpha,
    weights = weights,
    population_model = population_model
  )
  d <- parts$model$alpha_dim
  mu_score_particles <- sweep(parts$centered, 2L, parts$sigma2, "/")
  rho_score_particles <- -0.5 + sweep(parts$centered^2, 2L, 2 * parts$sigma2, "/")

  score <- c(
    colSums(mu_score_particles * parts$weights),
    colSums(rho_score_particles * parts$weights)
  )
  names(score) <- parts$model$hyper_names
  score
}

local_chart_curvature <- function(theta, alpha, weights = NULL, population_model) {
  parts <- .local_chart_diag_gaussian_parts(
    theta = theta,
    alpha = alpha,
    weights = weights,
    population_model = population_model
  )
  model <- parts$model
  d <- model$alpha_dim
  h <- model$hyper_dim

  mu_score_particles <- sweep(parts$centered, 2L, parts$sigma2, "/")
  rho_score_particles <- -0.5 + sweep(parts$centered^2, 2L, 2 * parts$sigma2, "/")
  complete_score <- cbind(mu_score_particles, rho_score_particles)
  colnames(complete_score) <- model$hyper_names

  score_mean <- colSums(complete_score * parts$weights)
  centered_score <- sweep(complete_score, 2L, score_mean, "-")
  score_cov <- crossprod(centered_score, centered_score * parts$weights)

  expected_complete_hessian <- matrix(0, nrow = h, ncol = h)
  for (j in seq_len(d)) {
    mu_pos <- j
    rho_pos <- d + j
    expected_complete_hessian[mu_pos, mu_pos] <- -parts$inv_sigma2[j]
    cross <- -sum(parts$weights * mu_score_particles[, j])
    expected_complete_hessian[mu_pos, rho_pos] <- cross
    expected_complete_hessian[rho_pos, mu_pos] <- cross
    expected_complete_hessian[rho_pos, rho_pos] <-
      -0.5 * sum(parts$weights * parts$centered[, j]^2 * parts$inv_sigma2[j])
  }

  curvature <- expected_complete_hessian + score_cov
  dimnames(curvature) <- list(model$hyper_names, model$hyper_names)
  curvature
}
