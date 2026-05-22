#!/usr/bin/env Rscript
# ============================================================================
# q0-predictive local SMC components
# - Build alpha references induced by theta proposals
# - Fit and bridge reusable local posterior-reference components
# ============================================================================

if (!exists("%||%", mode = "function") ||
    !exists("ESS", mode = "function") ||
    !exists("ll_parallel", mode = "function")) {
  source("smc_core.R")
}
if (!exists("make_reference_prior_gaussian", mode = "function") ||
    !exists("make_reference_prior_gaussian_mixture", mode = "function") ||
    !exists("combine_reference_priors", mode = "function") ||
    !exists("inflate_reference_prior", mode = "function") ||
    !exists("reference_prior_logpdf", mode = "function")) {
  source("reference_priors.R")
}
if (!exists("normalize_population_model", mode = "function") ||
    !exists("population_model_reference_components_from_theta", mode = "function")) {
  source("population_models.R")
}
if (!exists("run_tempered_smc", mode = "function")) {
  source("SMC_super_fast.R")
}
if (!exists("theta_proposal_sample", mode = "function")) {
  source("theta_proposals.R")
}
if (!exists("build_population_local_factor", mode = "function") ||
    !exists("build_population_factor_set", mode = "function") ||
    !exists("outer_population_smc", mode = "function") ||
    !exists("update_outer_population_fit", mode = "function") ||
    !exists("population_local_factor_ess", mode = "function") ||
    !exists("population_local_factor_log_marginal", mode = "function") ||
    !exists("validate_reference_local_object", mode = "function")) {
  source("outer_population_smc.R")
}

.local_predictive_normalize_weights <- function(w) {
  w <- pmax(as.numeric(w), 0)
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) {
    rep(1 / length(w), length(w))
  } else {
    w / sw
  }
}

.assert_local_predictive_smc_complete <- function(fit, label, tol = 1e-10) {
  final_lambda <- as.numeric(fit$final_lambda %||% NA_real_)
  if (!is.finite(final_lambda) || final_lambda < 1 - tol) {
    stop(label, " SMC stopped before lambda_target; log_evidence is not a completed local normalizer.")
  }
  invisible(TRUE)
}

make_population_predictive_reference_prior <- function(population_model,
                                                       theta = NULL,
                                                       theta_weights = NULL,
                                                       theta_proposal = NULL,
                                                       n_components = 64L,
                                                       component_scale = 1,
                                                       defensive_weight = 0.05,
                                                       defensive_scale = 16,
                                                       seed = NULL,
                                                       label = "q0_predictive_alpha_reference") {
  model <- normalize_population_model(population_model)
  if (is.null(theta)) {
    if (is.null(theta_proposal)) {
      stop("Provide theta or theta_proposal.")
    }
    theta <- theta_proposal_sample(theta_proposal, n = as.integer(n_components), seed = seed)
    theta_weights <- rep(1 / nrow(theta), nrow(theta))
  } else {
    theta <- .as_hyper_matrix(theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
    theta_weights <- pmax(as.numeric(theta_weights %||% rep(1, nrow(theta))), 0)
    sw <- sum(theta_weights)
    if (!is.finite(sw) || sw <= 0) {
      stop("theta_weights must have positive finite mass.")
    }
    theta_weights <- theta_weights / sw
    n_components <- as.integer(min(max(1L, n_components), nrow(theta)))
    if (nrow(theta) > n_components) {
      if (!is.null(seed)) set.seed(as.integer(seed))
      idx <- sample.int(nrow(theta), n_components, replace = FALSE, prob = theta_weights)
      theta <- theta[idx, , drop = FALSE]
      theta_weights <- theta_weights[idx] / sum(theta_weights[idx])
    }
  }

  ref_components <- population_model_reference_components_from_theta(model, theta)
  covs <- lapply(ref_components$component_covs, function(S) as.numeric(component_scale) * S)
  core <- make_reference_prior_gaussian_mixture(
    component_means = ref_components$component_means,
    component_covs = covs,
    weights = theta_weights,
    param_names = model$alpha_names,
    label = label
  )
  if (as.numeric(defensive_weight) <= 0) {
    return(core)
  }
  defensive <- inflate_reference_prior(
    reference_prior = core,
    scale = as.numeric(defensive_scale),
    label = paste0(label, "_defensive")
  )
  combine_reference_priors(
    priors = list(core, defensive),
    weights = c(1 - as.numeric(defensive_weight), as.numeric(defensive_weight)),
    label = label
  )
}

.predictive_theta_strata <- function(theta_proposal,
                                     population_model,
                                     n_strata,
                                     n_components,
                                     seed = NULL) {
  model <- normalize_population_model(population_model)
  n_strata <- as.integer(max(1L, n_strata))
  n_components <- as.integer(max(1L, n_components))
  if (n_strata == 1L) {
    return(list(list(theta = NULL, weight = 1)))
  }

  sample_n <- as.integer(max(n_strata * n_components * 4L, n_strata * 64L))
  theta <- theta_proposal_sample(theta_proposal, n = sample_n, seed = seed)
  theta <- .as_hyper_matrix(theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  ok <- apply(theta, 1L, function(x) all(is.finite(x)))
  theta <- theta[ok, , drop = FALSE]
  if (nrow(theta) < n_strata) {
    return(list(list(theta = theta, weight = 1)))
  }

  center <- colMeans(theta)
  scale <- apply(theta, 2L, stats::sd)
  scale[!is.finite(scale) | scale <= 0] <- 1
  z <- sweep(sweep(theta, 2L, center, "-"), 2L, scale, "/")
  if (!is.null(seed)) set.seed(as.integer(seed + 37L))
  cluster <- tryCatch(
    stats::kmeans(z, centers = n_strata, nstart = 8L, iter.max = 80L)$cluster,
    error = function(e) rep.int(1L, nrow(theta))
  )

  split_idx <- split(seq_len(nrow(theta)), cluster)
  split_idx <- split_idx[order(vapply(split_idx, length, integer(1)), decreasing = TRUE)]
  strata <- lapply(split_idx, function(idx) {
    list(
      theta = theta[idx, , drop = FALSE],
      weight = length(idx) / nrow(theta)
    )
  })
  strata[vapply(strata, function(x) nrow(x$theta) > 0L, logical(1))]
}

make_population_predictive_reference_priors <- function(population_model,
                                                        theta_proposal,
                                                        n_strata = 1L,
                                                        n_components = 64L,
                                                        component_scale = 1,
                                                        defensive_weight = 0.05,
                                                        defensive_scale = 16,
                                                        seed = NULL,
                                                        label = "q0_predictive_alpha_reference") {
  model <- normalize_population_model(population_model)
  strata <- .predictive_theta_strata(
    theta_proposal = theta_proposal,
    population_model = model,
    n_strata = n_strata,
    n_components = n_components,
    seed = seed
  )
  lapply(seq_along(strata), function(k) {
    stratum <- strata[[k]]
    if (is.null(stratum$theta)) {
      make_population_predictive_reference_prior(
        population_model = model,
        theta_proposal = theta_proposal,
        n_components = n_components,
        component_scale = component_scale,
        defensive_weight = defensive_weight,
        defensive_scale = defensive_scale,
        seed = if (is.null(seed)) NULL else as.integer(seed + 997L * k),
        label = sprintf("%s_stratum%d", label, k)
      )
    } else {
      make_population_predictive_reference_prior(
        population_model = model,
        theta = stratum$theta,
        theta_weights = rep(1 / nrow(stratum$theta), nrow(stratum$theta)),
        n_components = n_components,
        component_scale = component_scale,
        defensive_weight = defensive_weight,
        defensive_scale = defensive_scale,
        seed = if (is.null(seed)) NULL else as.integer(seed + 997L * k),
        label = sprintf("%s_stratum%d", label, k)
      )
    }
  })
}

fit_local_reference_smc_component <- function(local_id,
                                              data_i,
                                              loglik_fn,
                                              reference_prior,
                                              M,
                                              resample_threshold = 0.6,
                                              n_mcmc_moves = 2L,
                                              max_rounds = 100L,
                                              cess_target = 0.9,
                                              G_mix = 8L,
                                              rw_scale_init = 0.9,
                                              n_cores = 1L,
                                              seed = 123L,
                                              diagnostics = list(),
                                              verbose = FALSE) {
  fit <- run_tempered_smc(
    reference_prior = reference_prior,
    bridge_stat_fn = function(alpha) ll_parallel(alpha, data_i, loglik_fn, n_cores = as.integer(n_cores)),
    M = as.integer(M),
    resample_threshold = as.numeric(resample_threshold),
    n_mcmc_moves = as.integer(n_mcmc_moves),
    post_adapt_n_mcmc_moves = as.integer(n_mcmc_moves),
    max_rounds = as.integer(max_rounds),
    cess_target = as.numeric(cess_target),
    G_mix = as.integer(G_mix),
    rw_scale_init = as.numeric(rw_scale_init),
    n_cores = as.integer(n_cores),
    seed = as.integer(seed),
    verbose = verbose
  )
  .assert_local_predictive_smc_complete(fit, "q0-predictive local")
  local_reference_component_from_fit(
    local_id = local_id,
    fit = fit,
    data_i = data_i,
    loglik_fn = loglik_fn,
    reference_prior = reference_prior,
    n_cores = n_cores,
    diagnostics = diagnostics
  )
}

local_reference_component_from_fit <- function(local_id,
                                               fit,
                                               data_i,
                                               loglik_fn,
                                               reference_prior,
                                               n_cores = 1L,
                                               diagnostics = list()) {
  alpha <- as.matrix(fit$Theta)
  log_likelihood <- as.numeric(fit$log_likelihood %||% fit$loglik)
  if (length(log_likelihood) != nrow(alpha) || anyNA(log_likelihood)) {
    log_likelihood <- ll_parallel(alpha, data_i, loglik_fn, n_cores = as.integer(n_cores))
  }
  structure(
    list(
      local_id = as.integer(local_id),
      particles = alpha,
      weights = .local_predictive_normalize_weights(fit$w),
      reference_prior = reference_prior,
      log_reference_density = reference_prior_logpdf(reference_prior, alpha),
      proposal_type = "posterior_reference",
      log_likelihood = log_likelihood,
      log_evidence = as.numeric(fit$log_evidence),
      mcse_log_evidence = as.numeric(fit$mcse_logZ %||% NA_real_),
      diagnostics = modifyList(
        list(
          source = "local_reference_smc_component",
          final_lambda = as.numeric(fit$final_lambda %||% NA_real_),
          final_ess_frac = as.numeric(ESS(fit$w) / length(fit$w)),
          rounds = as.integer(fit$meta$rounds %||% NA_integer_)
        ),
        diagnostics %||% list()
      )
    ),
    class = "reference_local_component"
  )
}

bridge_local_reference_smc_component <- function(component,
                                                 data_i,
                                                 loglik_fn,
                                                 new_reference_prior,
                                                 M = NULL,
                                                 resample_threshold = 0.6,
                                                 n_mcmc_moves = 1L,
                                                 max_rounds = 60L,
                                                 cess_target = 0.9,
                                                 G_mix = 8L,
                                                 rw_scale_init = 0.7,
                                                 n_cores = 1L,
                                                 seed = 123L,
                                                 diagnostics = list(),
                                                 verbose = FALSE) {
  old_reference_prior <- component$reference_prior
  old_log_evidence <- as.numeric(component$log_evidence)
  if (!is.finite(old_log_evidence)) {
    stop("Cannot bridge a component without finite log_evidence.")
  }

  base_logpdf_fn <- function(alpha) {
    ll_parallel(alpha, data_i, loglik_fn, n_cores = as.integer(n_cores)) +
      reference_prior_logpdf(old_reference_prior, alpha) -
      old_log_evidence
  }
  bridge_stat_fn <- function(alpha) {
    reference_prior_logpdf(new_reference_prior, alpha) -
      reference_prior_logpdf(old_reference_prior, alpha)
  }

  fit <- run_tempered_smc(
    reference_prior = new_reference_prior,
    bridge_stat_fn = bridge_stat_fn,
    base_logpdf_fn = base_logpdf_fn,
    initial_particles = component$particles,
    initial_weights = component$weights,
    initial_log_normalizer = old_log_evidence,
    M = as.integer(M %||% nrow(component$particles)),
    resample_threshold = as.numeric(resample_threshold),
    n_mcmc_moves = as.integer(n_mcmc_moves),
    post_adapt_n_mcmc_moves = as.integer(n_mcmc_moves),
    max_rounds = as.integer(max_rounds),
    cess_target = as.numeric(cess_target),
    G_mix = as.integer(G_mix),
    rw_scale_init = as.numeric(rw_scale_init),
    n_cores = as.integer(n_cores),
    seed = as.integer(seed),
    verbose = verbose
  )
  .assert_local_predictive_smc_complete(fit, "bridge-repair local")
  fit$log_likelihood <- ll_parallel(fit$Theta, data_i, loglik_fn, n_cores = as.integer(n_cores))
  local_reference_component_from_fit(
    local_id = component$local_id,
    fit = fit,
    data_i = data_i,
    loglik_fn = loglik_fn,
    reference_prior = new_reference_prior,
    n_cores = n_cores,
    diagnostics = modifyList(
      list(
        source = "bridged_local_reference_component",
        bridged_from = component$diagnostics$source %||% "component"
      ),
      diagnostics %||% list()
    )
  )
}

reference_local_object_from_components <- function(local_id,
                                                   components,
                                                   mixture_weights = NULL) {
  if (!length(components)) {
    stop("At least one local component is required.")
  }
  if (is.null(mixture_weights)) {
    mixture_weights <- vapply(components, function(component) nrow(component$particles), numeric(1))
  }
  mixture_weights <- pmax(as.numeric(mixture_weights), 0)
  sw <- sum(mixture_weights)
  if (!is.finite(sw) || sw <= 0) {
    stop("mixture_weights must have positive finite mass.")
  }
  structure(
    list(
      local_id = as.integer(local_id),
      components = components,
      mixture_weights = mixture_weights / sw
    ),
    class = "reference_local_object"
  )
}

calibrate_reference_local_object_mbar <- function(local_object,
                                                  population_model,
                                                  data_i = NULL,
                                                  loglik_fn = NULL,
                                                  local_n_cores = 1L,
                                                  state_counts = c("ess", "particles"),
                                                  anchor = 1L,
                                                  max_iter = 1000L,
                                                  tol = 1e-8,
                                                  offset_drop = 50,
                                                  min_state_ess_frac = 0.01,
                                                  max_abs_shift = 100,
                                                  require_converged = FALSE,
                                                  apply_uncertified = FALSE,
                                                  verbose = FALSE) {
  model <- normalize_population_model(population_model)
  object <- validate_reference_local_object(local_object)
  S <- length(object$components)
  if (S <= 1L) {
    object$diagnostics <- modifyList(
      object$diagnostics %||% list(),
      list(mbar = list(enabled = FALSE, reason = "single_state"))
    )
    return(object)
  }

  state_counts <- match.arg(state_counts)
  anchor <- as.integer(max(1L, min(anchor, S)))
  alpha_list <- lapply(object$components, function(component) {
    alpha <- as.matrix(component$particles)
    if (!is.null(colnames(alpha)) && setequal(colnames(alpha), model$alpha_names)) {
      alpha[, model$alpha_names, drop = FALSE]
    } else {
      colnames(alpha) <- model$alpha_names
      alpha
    }
  })
  alpha <- do.call(rbind, alpha_list)
  source_state <- rep.int(seq_len(S), vapply(alpha_list, nrow, integer(1)))

  log_likelihood <- unlist(lapply(object$components, `[[`, "log_likelihood"), use.names = FALSE)
  if (length(log_likelihood) != nrow(alpha) || anyNA(log_likelihood)) {
    if (is.null(data_i) || is.null(loglik_fn)) {
      stop("MBAR calibration requires data_i and loglik_fn unless every component stores log_likelihood.")
    }
    log_likelihood <- ll_parallel(alpha, data_i, loglik_fn, n_cores = as.integer(local_n_cores))
  }
  log_likelihood <- as.numeric(log_likelihood)
  log_likelihood[!is.finite(log_likelihood)] <- -Inf

  log_ref <- vapply(object$components, function(component) {
    reference_prior_logpdf(component$reference_prior, alpha)
  }, numeric(nrow(alpha)))
  if (!is.matrix(log_ref)) {
    log_ref <- matrix(log_ref, ncol = S)
  }
  logu <- sweep(log_ref, 1L, log_likelihood, "+")

  log_sample_weight <- unlist(lapply(object$components, function(component) {
    log(.local_predictive_normalize_weights(component$weights))
  }), use.names = FALSE)
  log_sample_weight[!is.finite(log_sample_weight)] <- -Inf

  counts <- if (identical(state_counts, "ess")) {
    vapply(object$components, function(component) ESS(component$weights), numeric(1))
  } else {
    vapply(object$components, function(component) nrow(component$particles), numeric(1))
  }
  counts <- pmax(as.numeric(counts), 1)
  log_counts <- log(counts)

  initial_logZ <- vapply(object$components, function(component) {
    as.numeric(component$log_evidence %||% NA_real_)
  }, numeric(1))
  if (!any(is.finite(initial_logZ))) {
    stop("MBAR calibration needs at least one finite component log_evidence for the absolute offset.")
  }
  finite_initial <- is.finite(initial_logZ)
  if (!is.finite(initial_logZ[anchor])) {
    anchor <- which(finite_initial)[1L]
  }
  relative_logZ <- initial_logZ
  relative_logZ[!finite_initial] <- initial_logZ[anchor]
  relative_logZ <- relative_logZ - relative_logZ[anchor]

  log_mass <- log_counts[source_state] + log_sample_weight
  converged <- FALSE
  iter_used <- 0L
  for (iter in seq_len(as.integer(max_iter))) {
    log_den <- .rowLogSumExp(sweep(logu, 2L, log_counts - relative_logZ, "+"))
    next_relative <- vapply(seq_len(S), function(k) {
      logsumexp(log_mass + logu[, k] - log_den)
    }, numeric(1))
    next_relative <- next_relative - next_relative[anchor]
    delta <- max(abs(next_relative - relative_logZ), na.rm = TRUE)
    relative_logZ <- next_relative
    iter_used <- iter
    if (is.finite(delta) && delta <= as.numeric(tol)) {
      converged <- TRUE
      break
    }
  }

  log_den <- .rowLogSumExp(sweep(logu, 2L, log_counts - relative_logZ, "+"))
  state_ess_frac <- vapply(seq_len(S), function(k) {
    log_contrib <- log_mass + logu[, k] - log_den
    ok <- is.finite(log_contrib)
    if (!any(ok)) return(0)
    lw <- log_contrib[ok]
    ww <- exp(lw - logsumexp(lw))
    as.numeric(1 / sum(ww * ww) / length(log_contrib))
  }, numeric(1))

  offsets <- initial_logZ - relative_logZ
  offset_keep <- finite_initial & is.finite(offsets)
  if (any(offset_keep)) {
    max_offset <- max(offsets[offset_keep])
    offset_keep <- offset_keep & offsets >= max_offset - as.numeric(offset_drop)
  }
  if (!any(offset_keep)) {
    offset_keep <- finite_initial
  }
  offset_weights <- counts * pmax(state_ess_frac, 1e-6)
  offset_weights[!offset_keep] <- 0
  offset_weights <- .local_predictive_normalize_weights(offset_weights)
  absolute_offset <- sum(offset_weights * offsets, na.rm = TRUE)
  calibrated_logZ <- relative_logZ + absolute_offset
  shifts <- calibrated_logZ - initial_logZ
  max_shift <- max(abs(shifts), na.rm = TRUE)
  min_overlap <- min(state_ess_frac, na.rm = TRUE)
  certified <- is.finite(max_shift) &&
    is.finite(min_overlap) &&
    min_overlap >= as.numeric(min_state_ess_frac) &&
    max_shift <= as.numeric(max_abs_shift) &&
    (!isTRUE(require_converged) || isTRUE(converged))
  applied <- isTRUE(certified) || isTRUE(apply_uncertified)

  for (k in seq_len(S)) {
    object$components[[k]]$raw_log_evidence <- initial_logZ[k]
    if (isTRUE(applied)) {
      object$components[[k]]$log_evidence <- calibrated_logZ[k]
    }
    object$components[[k]]$diagnostics <- modifyList(
      object$components[[k]]$diagnostics %||% list(),
      list(
        mbar_log_evidence_shift = as.numeric(calibrated_logZ[k] - initial_logZ[k]),
        mbar_state_ess_frac = as.numeric(state_ess_frac[k])
      )
    )
  }
  object$diagnostics <- modifyList(
    object$diagnostics %||% list(),
    list(
      mbar = list(
        enabled = TRUE,
        converged = converged,
        iterations = as.integer(iter_used),
        anchor = as.integer(anchor),
        state_counts = state_counts,
        raw_log_evidence = initial_logZ,
        calibrated_log_evidence = calibrated_logZ,
        log_evidence_shift = shifts,
        absolute_offset = absolute_offset,
        offset_keep = offset_keep,
        state_ess_frac = state_ess_frac,
        max_abs_shift = max_shift,
        min_state_ess_frac = min_overlap,
        min_required_state_ess_frac = as.numeric(min_state_ess_frac),
        max_allowed_abs_shift = as.numeric(max_abs_shift),
        require_converged = isTRUE(require_converged),
        certified = isTRUE(certified),
        applied = isTRUE(applied)
      )
    )
  )
  if (isTRUE(verbose)) {
    cat(sprintf(
      "Local %d MBAR: states=%d | converged=%s | certified=%s | applied=%s | max |shift|=%.3f | min state ESS=%.4f\n",
      as.integer(object$local_id),
      S,
      if (converged) "yes" else "no",
      if (certified) "yes" else "no",
      if (applied) "yes" else "no",
      max_shift,
      min_overlap
    ))
  }
  validate_reference_local_object(object)
}

.reference_local_object_mbar_certified <- function(local_object) {
  mbar <- local_object$diagnostics$mbar %||% NULL
  is.list(mbar) && (!isTRUE(mbar$enabled) || isTRUE(mbar$certified))
}

make_population_theta_reference_prior <- function(population_model,
                                                  theta,
                                                  scale = 1,
                                                  label = "theta_alpha_reference") {
  model <- normalize_population_model(population_model)
  theta <- .as_hyper_matrix(theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  if (nrow(theta) != 1L) {
    stop("make_population_theta_reference_prior expects exactly one theta row.")
  }
  components <- population_model_reference_components_from_theta(model, theta)
  make_reference_prior_gaussian(
    mu = components$component_means[[1L]],
    Sigma = as.numeric(scale) * components$component_covs[[1L]],
    param_names = model$alpha_names,
    label = label
  )
}

reference_local_object_factor <- function(local_object,
                                          population_model,
                                          data_i = NULL,
                                          loglik_fn = NULL,
                                          local_n_cores = 1L) {
  build_population_local_factor(
    local_object = validate_reference_local_object(local_object),
    population_model = population_model,
    data_i = data_i,
    loglik_fn = loglik_fn,
    local_n_cores = local_n_cores
  )
}

reference_local_object_component_factors <- function(local_object,
                                                     population_model,
                                                     data_i = NULL,
                                                     loglik_fn = NULL,
                                                     local_n_cores = 1L) {
  object <- validate_reference_local_object(local_object)
  lapply(seq_along(object$components), function(k) {
    reference_local_object_factor(
      local_object = reference_local_object_from_components(
        local_id = object$local_id,
        components = list(object$components[[k]]),
        mixture_weights = 1
      ),
      population_model = population_model,
      data_i = data_i,
      loglik_fn = loglik_fn,
      local_n_cores = local_n_cores
    )
  })
}

.reference_local_object_loo_delta <- function(local_object,
                                              theta,
                                              population_model,
                                              data_i = NULL,
                                              loglik_fn = NULL,
                                              local_n_cores = 1L,
                                              full_log_marginal = NULL) {
  object <- validate_reference_local_object(local_object)
  if (length(object$components) <= 1L) {
    return(0)
  }
  full_log_marginal <- full_log_marginal %||% population_local_factor_log_marginal(
    reference_local_object_factor(
      local_object = object,
      population_model = population_model,
      data_i = data_i,
      loglik_fn = loglik_fn,
      local_n_cores = local_n_cores
    ),
    theta = theta
  )
  deltas <- vapply(seq_along(object$components), function(k) {
    keep <- setdiff(seq_along(object$components), k)
    kept_weights <- object$mixture_weights[keep]
    loo_object <- reference_local_object_from_components(
      local_id = object$local_id,
      components = object$components[keep],
      mixture_weights = kept_weights
    )
    loo_factor <- reference_local_object_factor(
      local_object = loo_object,
      population_model = population_model,
      data_i = data_i,
      loglik_fn = loglik_fn,
      local_n_cores = local_n_cores
    )
    as.numeric(population_local_factor_log_marginal(loo_factor, theta = theta) - full_log_marginal)
  }, numeric(1))
  finite <- is.finite(deltas)
  if (!any(finite)) Inf else max(abs(deltas[finite]))
}

.reference_local_factor_log_weights <- function(factor, theta) {
  stopifnot(inherits(factor, "population_local_factor"))
  theta_prepared <- population_model_prepare_theta(factor$population_model, theta)
  logp <- population_log_alpha_given_sufficient_stats_many(
    model = factor$population_model,
    sufficient_stats = factor$sufficient_stats,
    theta_prepared = theta_prepared
  )[1L, ]
  lw <- factor$log_base + logp
  keep <- is.finite(lw)
  lse <- if (any(keep)) logsumexp(lw[keep]) else -Inf
  list(
    lw = lw,
    finite = keep,
    log_marginal = as.numeric(lse + factor$log_constant),
    theta_prepared = theta_prepared
  )
}

.reference_log_mcse_from_weights <- function(lw_info) {
  keep <- lw_info$finite
  if (!any(keep) || !is.finite(lw_info$log_marginal)) {
    return(Inf)
  }
  lw <- lw_info$lw[keep]
  lse <- logsumexp(lw)
  w <- exp(lw - lse)
  ess <- 1 / sum(w * w)
  n <- length(w)
  if (!is.finite(ess) || ess <= 0 || n <= 1L) {
    return(Inf)
  }
  sqrt(max(1 / ess - 1 / n, 0))
}

.reference_split_log_marginal_delta <- function(lw_info) {
  keep <- which(lw_info$finite)
  if (length(keep) < 4L) {
    return(0)
  }
  a <- keep[seq.int(1L, length(keep), by = 2L)]
  b <- setdiff(keep, a)
  if (!length(a) || !length(b)) {
    return(0)
  }
  log_a <- logsumexp(lw_info$lw[a]) + log(length(keep) / length(a))
  log_b <- logsumexp(lw_info$lw[b]) + log(length(keep) / length(b))
  abs(as.numeric(log_a - log_b))
}

.reference_gaussian_score_matrix <- function(factor, theta_prepared) {
  model <- normalize_population_model(factor$population_model)
  if (!identical(model$fast_family %||% NULL, "gaussian") ||
      is.null(factor$particles) ||
      !identical(theta_prepared$family, "gaussian") ||
      !identical(theta_prepared$quadratic_kind, "diag")) {
    return(NULL)
  }
  alpha <- as.matrix(factor$particles)
  if (!nrow(alpha)) return(NULL)
  theta <- theta_prepared$theta[1L, , drop = FALSE]
  mu <- as.numeric(theta[1L, seq_len(model$alpha_dim)])
  log_sigma2 <- as.numeric(theta[1L, model$alpha_dim + seq_len(model$alpha_dim)])
  sigma2 <- exp(log_sigma2)
  diff <- sweep(alpha, 2L, mu, "-")
  grad_mu <- sweep(diff, 2L, sigma2, "/")
  grad_log_sigma2 <- -0.5 + 0.5 * sweep(diff * diff, 2L, sigma2, "/")
  out <- cbind(grad_mu, grad_log_sigma2)
  colnames(out) <- model$hyper_names
  out
}

.reference_gradient_split_diagnostic <- function(factor, lw_info) {
  score <- .reference_gaussian_score_matrix(factor, lw_info$theta_prepared)
  if (is.null(score)) {
    return(list(delta = 0, z = 0))
  }
  keep <- which(lw_info$finite)
  if (length(keep) < 4L) {
    return(list(delta = 0, z = 0))
  }
  lse <- logsumexp(lw_info$lw[keep])
  w <- exp(lw_info$lw[keep] - lse)
  a_pos <- seq.int(1L, length(keep), by = 2L)
  b_pos <- setdiff(seq_along(keep), a_pos)
  if (!length(a_pos) || !length(b_pos)) {
    return(list(delta = 0, z = 0))
  }
  score_moments <- function(pos) {
    ww <- w[pos]
    sw <- sum(ww)
    if (!is.finite(sw) || sw <= 0) {
      return(list(mean = rep(NA_real_, ncol(score)), se = rep(NA_real_, ncol(score))))
    }
    ww <- ww / sw
    x <- score[keep[pos], , drop = FALSE]
    mu <- colSums(x * ww)
    centered <- sweep(x, 2L, mu, "-")
    variance <- colSums(centered * centered * ww)
    ess <- 1 / sum(ww * ww)
    se <- sqrt(pmax(variance, 0) / max(ess, 1))
    list(mean = mu, se = se)
  }
  a <- score_moments(a_pos)
  b <- score_moments(b_pos)
  ga <- a$mean
  gb <- b$mean
  delta <- ga - gb
  finite <- is.finite(delta)
  if (!any(finite)) {
    return(list(delta = Inf, z = Inf))
  }
  se <- sqrt(a$se * a$se + b$se * b$se)
  z <- abs(delta) / pmax(se, sqrt(.Machine$double.eps))
  z <- z[is.finite(z)]
  list(
    delta = sqrt(mean(delta[finite] * delta[finite])),
    z = if (length(z)) pmax(sqrt(mean(z * z)) - 1, 0) else Inf
  )
}

.reference_gradient_split_delta <- function(factor, lw_info) {
  .reference_gradient_split_diagnostic(factor, lw_info)$delta
}

.reference_component_masses <- function(factor, lw_info, n_components = NULL) {
  component_id <- as.integer(factor$component_id %||% rep.int(1L, length(lw_info$lw)))
  n_components <- as.integer(n_components %||% max(component_id, na.rm = TRUE))
  out <- rep(0, n_components)
  keep <- which(lw_info$finite)
  if (!length(keep) || !is.finite(logsumexp(lw_info$lw[keep]))) {
    return(out)
  }
  lse <- logsumexp(lw_info$lw[keep])
  mass <- exp(lw_info$lw[keep] - lse)
  id <- component_id[keep]
  ok <- is.finite(id) & id >= 1L & id <= n_components
  if (any(ok)) {
    summed <- rowsum(matrix(mass[ok], ncol = 1L), id[ok], reorder = FALSE)
    full <- rep(0, n_components)
    full[as.integer(rownames(summed))] <- as.numeric(summed[, 1L])
    out <- full
  }
  out
}

.reference_component_mass_summary <- function(component_mass) {
  mass <- as.numeric(component_mass)
  mass <- mass[is.finite(mass) & mass >= 0]
  if (!length(mass) || sum(mass) <= 0) {
    return(list(max_mass = NA_real_, ess_frac = NA_real_))
  }
  mass <- mass / sum(mass)
  ess <- 1 / sum(mass * mass)
  list(
    max_mass = as.numeric(max(mass)),
    ess_frac = as.numeric(ess / length(mass))
  )
}

.reference_component_surface_delta <- function(component_log_marginal,
                                               full_log_marginal,
                                               component_mass = NULL,
                                               max_delta = 10) {
  vals <- as.numeric(component_log_marginal)
  finite <- is.finite(vals)
  if (sum(finite) <= 1L || !is.finite(full_log_marginal)) {
    return(0)
  }
  delta <- pmin(abs(vals - as.numeric(full_log_marginal)), as.numeric(max_delta))
  if (!is.null(component_mass) && length(component_mass) == length(vals)) {
    w <- as.numeric(component_mass)
    w[!is.finite(w) | w < 0] <- 0
    w[!finite] <- 0
    sw <- sum(w)
    if (is.finite(sw) && sw > 0) {
      w <- w / sw
      return(as.numeric(sqrt(sum(w * delta * delta))))
    }
  }
  vals <- vals[finite]
  as.numeric(stats::quantile(
    abs(vals - as.numeric(full_log_marginal)),
    probs = 0.9,
    names = FALSE,
    type = 8
  ))
}

.reference_surface_uncertainty <- function(log_mcse,
                                           pareto_k,
                                           loo_max_abs_delta,
                                           split_delta,
                                           component_delta,
                                           gradient_split_delta,
                                           max_delta = 10) {
  cap <- function(x) {
    x <- as.numeric(x)
    x[is.na(x)] <- 0
    x[is.infinite(x)] <- max_delta
    x[!is.finite(x)] <- 0
    pmin(abs(x), max_delta)
  }
  psis_penalty <- if (is.finite(pareto_k)) {
    pmax(as.numeric(pareto_k) - 0.5, 0)
  } else {
    0
  }
  sqrt(
    cap(log_mcse)^2 +
      cap(loo_max_abs_delta)^2 +
      cap(split_delta)^2 +
      cap(component_delta)^2 +
      cap(gradient_split_delta)^2 +
      cap(psis_penalty)^2
  )
}

audit_reference_local_object_support <- function(local_object,
                                                 theta,
                                                 population_model,
                                                 data_i = NULL,
                                                 loglik_fn = NULL,
                                                 local_n_cores = 1L,
                                                 min_ess_frac = 0.05,
                                                 compute_psis = FALSE,
                                                 compute_loo = FALSE) {
  model <- normalize_population_model(population_model)
  object <- validate_reference_local_object(local_object)
  theta <- .as_hyper_matrix(theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  if (nrow(theta) != 1L) {
    stop("audit_reference_local_object_support expects exactly one theta row.")
  }

  factor <- reference_local_object_factor(
    local_object = object,
    population_model = model,
    data_i = data_i,
    loglik_fn = loglik_fn,
    local_n_cores = local_n_cores
  )
  lw_info <- .reference_local_factor_log_weights(factor, theta)
  ess_frac <- as.numeric(population_local_factor_ess(factor, theta) / factor$n_particles)
  log_marginal <- as.numeric(lw_info$log_marginal)
  pareto_k <- if (isTRUE(compute_psis)) {
    as.numeric(population_local_factor_tail_diagnostic(factor, theta, use_psis = TRUE)$pareto_k)
  } else {
    NA_real_
  }
  loo_max_abs_delta <- if (isTRUE(compute_loo)) {
    .reference_local_object_loo_delta(
      local_object = object,
      theta = theta,
      population_model = model,
      data_i = data_i,
      loglik_fn = loglik_fn,
      local_n_cores = local_n_cores,
      full_log_marginal = log_marginal
    )
  } else {
    NA_real_
  }

  component_factors <- reference_local_object_component_factors(
    local_object = object,
    population_model = model,
    data_i = data_i,
    loglik_fn = loglik_fn,
    local_n_cores = local_n_cores
  )
  component_ess <- vapply(
    component_factors,
    function(component_factor) {
      as.numeric(population_local_factor_ess(component_factor, theta) / component_factor$n_particles)
    },
    numeric(1)
  )
  component_log_marginal <- vapply(
    component_factors,
    population_local_factor_log_marginal,
    numeric(1),
    theta = theta
  )
  best_component <- which.max(component_ess)
  component_mass <- .reference_component_masses(
    factor = factor,
    lw_info = lw_info,
    n_components = length(object$components)
  )
  component_mass_summary <- .reference_component_mass_summary(component_mass)
  log_mcse <- .reference_log_mcse_from_weights(lw_info)
  split_delta <- .reference_split_log_marginal_delta(lw_info)
  gradient_split <- .reference_gradient_split_diagnostic(factor, lw_info)
  gradient_split_delta <- as.numeric(gradient_split$delta)
  component_delta <- .reference_component_surface_delta(
    component_log_marginal,
    log_marginal,
    component_mass = component_mass
  )
  surface_uncertainty <- .reference_surface_uncertainty(
    log_mcse = log_mcse,
    pareto_k = pareto_k,
    loo_max_abs_delta = loo_max_abs_delta,
    split_delta = split_delta,
    component_delta = component_delta,
    gradient_split_delta = gradient_split$z
  )

  list(
    local_id = as.integer(object$local_id),
    covered = is.finite(ess_frac) && ess_frac >= as.numeric(min_ess_frac),
    ess_frac = ess_frac,
    log_marginal = log_marginal,
    pareto_k = pareto_k,
    loo_max_abs_delta = loo_max_abs_delta,
    log_marginal_mcse = log_mcse,
    split_log_marginal_delta = split_delta,
    component_log_marginal_delta = component_delta,
    gradient_split_delta = gradient_split_delta,
    gradient_split_z = as.numeric(gradient_split$z),
    surface_uncertainty = surface_uncertainty,
    max_component_mass = as.numeric(component_mass_summary$max_mass),
    component_mass_ess_frac = as.numeric(component_mass_summary$ess_frac),
    best_component = as.integer(best_component),
    best_component_ess_frac = as.numeric(component_ess[best_component]),
    component_ess_frac = component_ess,
    component_log_marginal = component_log_marginal,
    component_mass = component_mass,
    n_components = as.integer(length(object$components)),
    n_particles = as.integer(sum(vapply(object$components, function(component) nrow(component$particles), integer(1))))
  )
}

audit_reference_local_objects <- function(local_objects,
                                          theta_audit,
                                          population_model,
                                          data_list = NULL,
                                          loglik_fn = NULL,
                                          local_ids = seq_along(local_objects),
                                          min_ess_frac = 0.05,
                                          max_pareto_k = Inf,
                                          max_loo_delta = Inf,
                                          compute_psis = FALSE,
                                          compute_loo = FALSE,
                                          n_jobs = 1L,
                                          local_n_cores = 1L) {
  model <- normalize_population_model(population_model)
  theta_source <- attr(theta_audit, "theta_source")
  theta_source_index <- attr(theta_audit, "theta_source_index")
  theta_protected <- attr(theta_audit, "theta_protected")
  theta_source_rank <- attr(theta_audit, "theta_source_rank")
  theta_log_impact <- attr(theta_audit, "theta_log_impact")
  theta_impact_weight <- attr(theta_audit, "theta_impact_weight")
  theta_log_q0 <- attr(theta_audit, "theta_log_q0")
  theta_audit <- .as_hyper_matrix(theta_audit, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  theta_source <- theta_source %||% rep(NA_character_, nrow(theta_audit))
  theta_source_index <- theta_source_index %||% rep(NA_integer_, nrow(theta_audit))
  theta_protected <- theta_protected %||% rep(FALSE, nrow(theta_audit))
  theta_source_rank <- theta_source_rank %||% rep(100L, nrow(theta_audit))
  theta_log_impact <- theta_log_impact %||% rep(NA_real_, nrow(theta_audit))
  theta_impact_weight <- theta_impact_weight %||% rep(NA_real_, nrow(theta_audit))
  theta_log_q0 <- theta_log_q0 %||% rep(NA_real_, nrow(theta_audit))
  local_ids <- sort(unique(as.integer(local_ids)))
  local_ids <- local_ids[local_ids >= 1L & local_ids <= length(local_objects)]
  if (!length(local_ids)) {
    stop("No valid local ids to audit.")
  }

  parts <- parallel::mclapply(
    local_ids,
    function(i) {
      object <- validate_reference_local_object(local_objects[[i]])
      data_i <- if (is.null(data_list)) NULL else data_list[[i]]
      factor <- reference_local_object_factor(
        local_object = object,
        population_model = model,
        data_i = data_i,
        loglik_fn = loglik_fn,
        local_n_cores = local_n_cores
      )
      component_factors <- reference_local_object_component_factors(
        local_object = object,
        population_model = model,
        data_i = data_i,
        loglik_fn = loglik_fn,
        local_n_cores = local_n_cores
      )
      n_particles <- as.integer(sum(vapply(object$components, function(component) nrow(component$particles), integer(1))))
      rows <- vector("list", nrow(theta_audit))
      for (j in seq_len(nrow(theta_audit))) {
        theta <- theta_audit[j, , drop = FALSE]
        lw_info <- .reference_local_factor_log_weights(factor, theta)
        ess_frac <- as.numeric(population_local_factor_ess(factor, theta) / factor$n_particles)
        log_marginal <- as.numeric(lw_info$log_marginal)
        tail <- if (isTRUE(compute_psis)) {
          population_local_factor_tail_diagnostic(factor, theta, use_psis = TRUE)
        } else {
          list(pareto_k = NA_real_)
        }
        loo_max_abs_delta <- if (isTRUE(compute_loo)) {
          .reference_local_object_loo_delta(
            local_object = object,
            theta = theta,
            population_model = model,
            data_i = data_i,
            loglik_fn = loglik_fn,
            local_n_cores = local_n_cores,
            full_log_marginal = log_marginal
          )
        } else {
          NA_real_
        }
        component_ess <- vapply(
          component_factors,
          function(component_factor) {
            as.numeric(population_local_factor_ess(component_factor, theta) / component_factor$n_particles)
          },
          numeric(1)
        )
        component_log_marginal <- vapply(
          component_factors,
          population_local_factor_log_marginal,
          numeric(1),
          theta = theta
        )
        best_component <- which.max(component_ess)
        component_mass <- .reference_component_masses(
          factor = factor,
          lw_info = lw_info,
          n_components = length(object$components)
        )
        component_mass_summary <- .reference_component_mass_summary(component_mass)
        log_mcse <- .reference_log_mcse_from_weights(lw_info)
        split_delta <- .reference_split_log_marginal_delta(lw_info)
        gradient_split <- .reference_gradient_split_diagnostic(factor, lw_info)
        gradient_split_delta <- as.numeric(gradient_split$delta)
        component_delta <- .reference_component_surface_delta(
          component_log_marginal,
          log_marginal,
          component_mass = component_mass
        )
        surface_uncertainty <- .reference_surface_uncertainty(
          log_mcse = log_mcse,
          pareto_k = as.numeric(tail$pareto_k),
          loo_max_abs_delta = loo_max_abs_delta,
          split_delta = split_delta,
          component_delta = component_delta,
          gradient_split_delta = gradient_split$z
        )
        rows[[j]] <- cbind(
          data.frame(
          local_id = as.integer(i),
          theta_id = as.integer(j),
          covered = is.finite(ess_frac) &&
            ess_frac >= as.numeric(min_ess_frac) &&
            (!isTRUE(compute_psis) || !is.finite(tail$pareto_k) || tail$pareto_k <= as.numeric(max_pareto_k)) &&
            (!isTRUE(compute_loo) || !is.finite(loo_max_abs_delta) || loo_max_abs_delta <= as.numeric(max_loo_delta)),
          ess_frac = ess_frac,
          log_marginal = log_marginal,
          pareto_k = as.numeric(tail$pareto_k),
          loo_max_abs_delta = as.numeric(loo_max_abs_delta),
          log_marginal_mcse = as.numeric(log_mcse),
          split_log_marginal_delta = as.numeric(split_delta),
          component_log_marginal_delta = as.numeric(component_delta),
          gradient_split_delta = as.numeric(gradient_split_delta),
          gradient_split_z = as.numeric(gradient_split$z),
          surface_uncertainty = as.numeric(surface_uncertainty),
          max_component_mass = as.numeric(component_mass_summary$max_mass),
          component_mass_ess_frac = as.numeric(component_mass_summary$ess_frac),
          theta_source = as.character(theta_source[j]),
          theta_source_index = as.integer(theta_source_index[j]),
          theta_protected = isTRUE(theta_protected[j]),
          theta_source_rank = as.integer(theta_source_rank[j]),
          theta_log_impact = as.numeric(theta_log_impact[j]),
          theta_impact_weight = as.numeric(theta_impact_weight[j]),
          theta_log_q0 = as.numeric(theta_log_q0[j]),
          check.names = FALSE
          ),
          as.data.frame(theta, check.names = FALSE),
          data.frame(
          best_component = as.integer(best_component),
          best_component_ess_frac = as.numeric(component_ess[best_component]),
          n_components = as.integer(length(object$components)),
          n_particles = n_particles,
          check.names = FALSE
          )
        )
      }
      do.call(rbind, rows)
    },
    mc.cores = as.integer(max(1L, n_jobs))
  )
  .predictive_reference_add_surface_scores(
    audit = do.call(rbind, parts),
    population_model = model
  )
}

.predictive_reference_add_surface_scores <- function(audit, population_model) {
  if (is.null(audit) || !nrow(audit)) {
    return(audit)
  }
  model <- normalize_population_model(population_model)
  theta_cols <- intersect(model$hyper_names, names(audit))
  surface_uncertainty <- if ("surface_uncertainty" %in% names(audit)) {
    pmax(as.numeric(audit$surface_uncertainty), 0)
  } else {
    rep(0, nrow(audit))
  }
  surface_uncertainty[!is.finite(surface_uncertainty)] <- 10
  impact <- if ("theta_impact_weight" %in% names(audit)) {
    pmax(as.numeric(audit$theta_impact_weight), 0)
  } else {
    rep(0, nrow(audit))
  }
  impact[!is.finite(impact)] <- 0

  posterior_leverage <- rep(1, nrow(audit))
  if (length(theta_cols)) {
    theta_key <- as.character(audit$theta_id)
    first <- !duplicated(theta_key)
    theta_unique <- as.matrix(audit[first, theta_cols, drop = FALSE])
    impact_unique <- impact[first]
    if (!any(impact_unique > 0)) {
      impact_unique <- rep(1, length(impact_unique))
    }
    w <- .local_predictive_normalize_weights(impact_unique)
    center <- colSums(theta_unique * w)
    scale <- sqrt(colSums(sweep(theta_unique, 2L, center, "-")^2 * w))
    scale[!is.finite(scale) | scale <= 0] <- 1
    z <- sweep(sweep(as.matrix(audit[, theta_cols, drop = FALSE]), 2L, center, "-"), 2L, scale, "/")
    z[!is.finite(z)] <- 0
    radial <- sqrt(rowMeans(z * z))

    tail_rank <- matrix(0.5, nrow = nrow(audit), ncol = length(theta_cols))
    for (j in seq_along(theta_cols)) {
      vals <- theta_unique[, j]
      ord <- order(vals)
      ranks <- numeric(length(vals))
      ranks[ord] <- seq_along(vals)
      ranks <- (ranks - 0.5) / max(length(vals), 1L)
      map <- ranks[match(theta_key, theta_key[first])]
      tail_rank[, j] <- map
    }
    tail_depth <- apply(abs(tail_rank - 0.5) * 2, 1L, max)
    posterior_leverage <- 1 + radial + tail_depth
  }

  surface_score <- impact * posterior_leverage * surface_uncertainty
  surface_score[!is.finite(surface_score)] <- 0
  audit$posterior_leverage <- as.numeric(posterior_leverage)
  audit$surface_score <- as.numeric(surface_score)
  audit$needs_refinement <- !as.logical(audit$covered)
  audit
}

bridge_repair_reference_local_objects <- function(local_objects,
                                                  audit,
                                                  theta_audit,
                                                  data_list,
                                                  loglik_fn,
                                                  population_model,
                                                  max_repairs = 20L,
                                                  M = NULL,
                                                  min_ess_frac = 0.05,
                                                  reference_scale = 1,
                                                  resample_threshold = 0.6,
                                                  n_mcmc_moves = 1L,
                                                  max_rounds = 80L,
                                                  cess_target = 0.9,
                                                  G_mix = 8L,
                                                  rw_scale_init = 0.7,
                                                  n_jobs = 1L,
                                                  local_n_cores = 1L,
                                                  seed = 123L,
                                                  verbose = FALSE) {
  model <- normalize_population_model(population_model)
  theta_audit <- .as_hyper_matrix(theta_audit, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  failed <- audit[!audit$covered, , drop = FALSE]
  repair_schema <- data.frame(
    local_id = integer(0),
    theta_id = integer(0),
    old_ess_frac = numeric(0),
    best_component = integer(0),
    repair_log_evidence = numeric(0),
    repair_mcse_log_evidence = numeric(0),
    n_components_after = integer(0),
    n_particles_after = integer(0),
    check.names = FALSE
  )
  if (!nrow(failed) || max_repairs <= 0L) {
    return(list(local_objects = local_objects, repairs = repair_schema))
  }

  failed <- failed[order(failed$ess_frac, -failed$best_component_ess_frac), , drop = FALSE]
  failed <- failed[seq_len(min(nrow(failed), as.integer(max_repairs))), , drop = FALSE]
  by_local <- split(failed, failed$local_id)

  repaired <- parallel::mclapply(
    names(by_local),
    function(local_name) {
      i <- as.integer(local_name)
      object <- validate_reference_local_object(local_objects[[i]])
      rows_i <- by_local[[local_name]]
      repair_rows <- vector("list", 0L)
      for (r in seq_len(nrow(rows_i))) {
        theta_id <- as.integer(rows_i$theta_id[r])
        theta <- theta_audit[theta_id, , drop = FALSE]
        component_index <- as.integer(rows_i$best_component[r])
        component_index <- max(1L, min(component_index, length(object$components)))
        source_component <- object$components[[component_index]]
        repair_ref <- make_population_theta_reference_prior(
          population_model = model,
          theta = theta,
          scale = reference_scale,
          label = sprintf("bridge_repair_local%d_theta%d", i, theta_id)
        )
        repair_component <- bridge_local_reference_smc_component(
          component = source_component,
          data_i = data_list[[i]],
          loglik_fn = loglik_fn,
          new_reference_prior = repair_ref,
          M = as.integer(M %||% nrow(source_component$particles)),
          resample_threshold = resample_threshold,
          n_mcmc_moves = n_mcmc_moves,
          max_rounds = max_rounds,
          cess_target = cess_target,
          G_mix = G_mix,
          rw_scale_init = rw_scale_init,
          n_cores = local_n_cores,
          seed = as.integer(seed + 10000L * i + theta_id),
          diagnostics = list(
            role = "bridge_repair",
            theta_id = theta_id,
            source_component = component_index,
            old_ess_frac = as.numeric(rows_i$ess_frac[r]),
            min_ess_frac = as.numeric(min_ess_frac)
          ),
          verbose = verbose
        )
        object$components <- c(object$components, list(repair_component))
        object$mixture_weights <- .local_predictive_normalize_weights(
          vapply(object$components, function(component) nrow(component$particles), numeric(1))
        )
        object <- validate_reference_local_object(object)
        repair_rows[[length(repair_rows) + 1L]] <- data.frame(
          local_id = as.integer(i),
          theta_id = theta_id,
          old_ess_frac = as.numeric(rows_i$ess_frac[r]),
          best_component = component_index,
          repair_log_evidence = as.numeric(repair_component$log_evidence),
          repair_mcse_log_evidence = as.numeric(repair_component$mcse_log_evidence),
          n_components_after = as.integer(length(object$components)),
          n_particles_after = as.integer(sum(vapply(object$components, function(component) nrow(component$particles), integer(1)))),
          check.names = FALSE
        )
      }
      list(
        local_id = i,
        local_object = object,
        repairs = if (length(repair_rows)) do.call(rbind, repair_rows) else repair_schema
      )
    },
    mc.cores = as.integer(max(1L, n_jobs))
  )

  out <- local_objects
  repair_parts <- vector("list", length(repaired))
  for (k in seq_along(repaired)) {
    out[[repaired[[k]]$local_id]] <- repaired[[k]]$local_object
    repair_parts[[k]] <- repaired[[k]]$repairs
  }
  list(
    local_objects = out,
    repairs = if (length(repair_parts)) do.call(rbind, repair_parts) else repair_schema
  )
}

.predictive_reference_merge_control <- function(user, defaults) {
  user <- user %||% list()
  if (!is.list(user)) stop("control must be a list.")
  modifyList(defaults, user)
}

.predictive_reference_weighted_theta_sample <- function(theta, w, n, seed = NULL) {
  theta <- as.matrix(theta)
  n <- as.integer(n)
  if (n <= 0L || !nrow(theta)) {
    return(theta[0L, , drop = FALSE])
  }
  if (!is.null(seed)) set.seed(as.integer(seed))
  w <- .local_predictive_normalize_weights(w %||% rep(1 / nrow(theta), nrow(theta)))
  idx <- sample.int(nrow(theta), n, replace = TRUE, prob = w)
  theta[idx, , drop = FALSE]
}

.predictive_reference_bind_audit_theta <- function(parts,
                                                   population_model,
                                                   factor_set = NULL,
                                                   theta_proposal = NULL,
                                                   n_cores = 1L) {
  model <- normalize_population_model(population_model)
  parts <- parts[vapply(parts, function(part) !is.null(part$theta) && nrow(part$theta) > 0L, logical(1))]
  if (!length(parts)) {
    stop("Audit requires at least one theta source.")
  }
  theta <- do.call(rbind, lapply(parts, `[[`, "theta"))
  source <- unlist(lapply(parts, `[[`, "source"), use.names = FALSE)
  source_index <- unlist(lapply(parts, `[[`, "source_index"), use.names = FALSE)
  protected <- unlist(lapply(parts, function(part) {
    rep(isTRUE(part$protected), nrow(part$theta))
  }), use.names = FALSE)
  source_rank <- unlist(lapply(parts, function(part) {
    rep(as.integer(part$source_rank %||% 100L), nrow(part$theta))
  }), use.names = FALSE)
  theta <- .as_hyper_matrix(theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  ok <- apply(theta, 1L, function(x) all(is.finite(x)))
  theta <- theta[ok, , drop = FALSE]
  source <- source[ok]
  source_index <- source_index[ok]
  protected <- protected[ok]
  source_rank <- source_rank[ok]
  if (!nrow(theta)) stop("Audit theta set contains no finite rows.")

  log_impact <- rep(NA_real_, nrow(theta))
  if (!is.null(factor_set)) {
    log_impact <- tryCatch(
      population_factor_set_logposterior(
        factor_set = factor_set,
        theta = theta,
        include_constant = FALSE,
        n_cores = n_cores
      ),
      error = function(e) rep(NA_real_, nrow(theta))
    )
  }
  finite <- is.finite(log_impact)
  impact_weight <- rep(0, nrow(theta))
  if (any(finite)) {
    lw <- log_impact[finite] - max(log_impact[finite])
    impact_weight[finite] <- exp(lw)
    sw <- sum(impact_weight)
    if (is.finite(sw) && sw > 0) {
      impact_weight <- impact_weight / sw
    }
  }

  attr(theta, "theta_source") <- source
  attr(theta, "theta_source_index") <- as.integer(source_index)
  attr(theta, "theta_protected") <- as.logical(protected)
  attr(theta, "theta_source_rank") <- as.integer(source_rank)
  attr(theta, "theta_log_impact") <- as.numeric(log_impact)
  attr(theta, "theta_impact_weight") <- as.numeric(impact_weight)
  attr(theta, "theta_log_q0") <- if (!is.null(theta_proposal)) {
    theta_proposal_log_density(theta_proposal, theta)
  } else {
    rep(NA_real_, nrow(theta))
  }
  theta
}

.predictive_reference_axis_theta <- function(theta,
                                             w = NULL,
                                             population_model,
                                             probs = c(0.01, 0.025, 0.05, 0.95, 0.975, 0.99),
                                             inflate = 1,
                                             max_points = 32L,
                                             principal_dims = 2L) {
  model <- normalize_population_model(population_model)
  theta <- .as_hyper_matrix(theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  ok <- apply(theta, 1L, function(x) all(is.finite(x)))
  theta <- theta[ok, , drop = FALSE]
  if (!nrow(theta) || as.integer(max_points) <= 0L) {
    return(matrix(numeric(0), nrow = 0L, ncol = model$hyper_dim, dimnames = list(NULL, model$hyper_names)))
  }
  w <- .local_predictive_normalize_weights(w %||% rep(1, nrow(theta)))
  center <- colSums(theta * w)
  names(center) <- model$hyper_names
  spread <- apply(theta, 2L, stats::sd)
  spread[!is.finite(spread)] <- 0
  dims <- order(spread, decreasing = TRUE)
  probs <- as.numeric(probs)
  probs <- probs[is.finite(probs) & probs > 0 & probs < 1]
  if (!length(probs)) {
    return(matrix(numeric(0), nrow = 0L, ncol = model$hyper_dim, dimnames = list(NULL, model$hyper_names)))
  }

  rows <- list()
  principal_dims <- as.integer(max(0L, min(principal_dims, ncol(theta))))
  if (principal_dims > 0L && nrow(theta) > ncol(theta) + 1L) {
    S <- tryCatch(weighted_cov(theta, w), error = function(e) NULL)
    if (!is.null(S) && all(is.finite(S))) {
      S <- regularize_cov(S, min_eig = 1e-8, cond_cap = 1e8)
      eig <- tryCatch(eigen(S, symmetric = TRUE), error = function(e) NULL)
      if (!is.null(eig)) {
        centered <- sweep(theta, 2L, center, "-")
        for (k in seq_len(principal_dims)) {
          direction <- eig$vectors[, k]
          score <- as.numeric(centered %*% direction)
          vals <- .weighted_quantile(score, w, probs = probs)
          for (val in vals) {
            row <- center + as.numeric(inflate) * as.numeric(val) * direction
            names(row) <- model$hyper_names
            rows[[length(rows) + 1L]] <- row
          }
        }
      }
    }
  }
  for (j in dims) {
    vals <- .weighted_quantile(theta[, j], w, probs = probs)
    for (val in vals) {
      row <- center
      row[j] <- center[j] + as.numeric(inflate) * (as.numeric(val) - center[j])
      rows[[length(rows) + 1L]] <- row
    }
  }
  out <- do.call(rbind, rows)
  colnames(out) <- model$hyper_names
  key <- apply(round(out, 8L), 1L, paste, collapse = "\r")
  out <- out[!duplicated(key), , drop = FALSE]
  out[seq_len(min(nrow(out), as.integer(max_points))), , drop = FALSE]
}

.predictive_reference_tail_theta <- function(theta_proposal,
                                             population_model,
                                             n_tail,
                                             seed = NULL) {
  model <- normalize_population_model(population_model)
  n_tail <- as.integer(n_tail)
  if (n_tail <= 0L) {
    return(matrix(numeric(0), nrow = 0L, ncol = model$hyper_dim, dimnames = list(NULL, model$hyper_names)))
  }
  pool_n <- as.integer(max(1000L, 50L * n_tail))
  pool <- theta_proposal_sample(theta_proposal, n = pool_n, seed = seed)
  pool <- .as_hyper_matrix(pool, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  ok <- apply(pool, 1L, function(x) all(is.finite(x)))
  pool <- pool[ok, , drop = FALSE]
  if (!nrow(pool)) {
    return(pool)
  }

  idx <- integer(0)
  probs <- c(0.005, 0.01, 0.025, 0.05, 0.95, 0.975, 0.99, 0.995)
  for (j in seq_len(ncol(pool))) {
    q <- stats::quantile(pool[, j], probs = probs, names = FALSE, type = 8)
    idx <- c(idx, vapply(q, function(v) which.min(abs(pool[, j] - v)), integer(1)))
  }
  center <- colMeans(pool)
  S <- tryCatch(regularize_cov(stats::cov(pool), min_eig = 1e-8, cond_cap = 1e8), error = function(e) diag(ncol(pool)))
  L <- tryCatch(chol(S), error = function(e) diag(ncol(pool)))
  z <- t(backsolve(L, t(sweep(pool, 2L, center, "-")), transpose = TRUE))
  r2 <- rowSums(z * z)
  idx <- unique(idx)
  extra <- as.integer(n_tail) - length(idx)
  if (extra > 0L) {
    marginal_lo <- apply(pool, 2L, stats::quantile, probs = 0.005, names = FALSE, type = 8)
    marginal_hi <- apply(pool, 2L, stats::quantile, probs = 0.995, names = FALSE, type = 8)
    radial_pool <- apply(
      sweep(pool, 2L, marginal_lo, ">=") &
        sweep(pool, 2L, marginal_hi, "<="),
      1L,
      all
    )
    radial_pool[idx] <- FALSE
    radial_idx <- which(radial_pool & is.finite(r2))
    radial_probs <- seq(0.75, 0.975, length.out = extra)
    if (length(radial_idx)) {
      radial_targets <- stats::quantile(r2[radial_idx], probs = radial_probs, names = FALSE, type = 8)
      idx <- unique(c(idx, vapply(
        radial_targets,
        function(target) radial_idx[which.min(abs(r2[radial_idx] - target))],
        integer(1)
      )))
    }
  }
  idx <- idx[seq_len(min(length(idx), n_tail))]
  out <- pool[idx, , drop = FALSE]
  colnames(out) <- model$hyper_names
  out
}

.predictive_reference_audit_theta <- function(theta_proposal,
                                              outer_fit,
                                              factor_set,
                                              population_model,
                                              n_q0 = 64L,
                                              n_outer = 64L,
                                              n_tail = 32L,
                                              n_q0_axis = 32L,
                                              n_outer_axis = 32L,
                                              axis_tail_probs = c(0.005, 0.01, 0.025, 0.05, 0.95, 0.975, 0.99, 0.995),
                                              outer_axis_inflate = 1.5,
                                              explicit_theta = NULL,
                                              n_cores = 1L,
                                              seed = NULL) {
  model <- normalize_population_model(population_model)
  parts <- list()
  if (!is.null(explicit_theta)) {
    explicit <- .as_hyper_matrix(
      explicit_theta,
      hyper_names = model$hyper_names,
      hyper_dim = model$hyper_dim
    )
    parts[[length(parts) + 1L]] <- list(
      theta = explicit,
      source = rep("explicit", nrow(explicit)),
      source_index = seq_len(nrow(explicit)),
      protected = TRUE,
      source_rank = 1L
    )
  }
  if (as.integer(n_q0) > 0L) {
    q0 <- theta_proposal_sample(
      theta_proposal,
      n = as.integer(n_q0),
      seed = seed
    )
    parts[[length(parts) + 1L]] <- list(
      theta = q0,
      source = rep("q0", nrow(q0)),
      source_index = seq_len(nrow(q0)),
      protected = FALSE,
      source_rank = 50L
    )
    if (as.integer(n_q0_axis) > 0L) {
      q0_axis <- .predictive_reference_axis_theta(
        theta = q0,
        population_model = model,
        probs = axis_tail_probs,
        inflate = 1,
        max_points = as.integer(n_q0_axis)
      )
      parts[[length(parts) + 1L]] <- list(
        theta = q0_axis,
        source = rep("q0_axis", nrow(q0_axis)),
        source_index = seq_len(nrow(q0_axis)),
        protected = TRUE,
        source_rank = 2L
      )
    }
  }
  if (!is.null(outer_fit) && as.integer(n_outer) > 0L) {
    outer <- .predictive_reference_weighted_theta_sample(
      theta = outer_fit$theta,
      w = outer_fit$w,
      n = as.integer(n_outer),
      seed = if (is.null(seed)) NULL else as.integer(seed + 17L)
    )
    parts[[length(parts) + 1L]] <- list(
      theta = outer,
      source = rep("outer", nrow(outer)),
      source_index = seq_len(nrow(outer)),
      protected = FALSE,
      source_rank = 40L
    )
    if (as.integer(n_outer_axis) > 0L) {
      outer_axis <- .predictive_reference_axis_theta(
        theta = outer_fit$theta,
        w = outer_fit$w,
        population_model = model,
        probs = axis_tail_probs,
        inflate = outer_axis_inflate,
        max_points = as.integer(n_outer_axis)
      )
      parts[[length(parts) + 1L]] <- list(
        theta = outer_axis,
        source = rep("outer_axis", nrow(outer_axis)),
        source_index = seq_len(nrow(outer_axis)),
        protected = TRUE,
        source_rank = 3L
      )
    }
  }
  if (as.integer(n_tail) > 0L) {
    tail <- .predictive_reference_tail_theta(
      theta_proposal = theta_proposal,
      population_model = model,
      n_tail = as.integer(n_tail),
      seed = if (is.null(seed)) NULL else as.integer(seed + 31L)
    )
    parts[[length(parts) + 1L]] <- list(
      theta = tail,
      source = rep("q0_tail", nrow(tail)),
      source_index = seq_len(nrow(tail)),
      protected = TRUE,
      source_rank = 4L
    )
  }
  .predictive_reference_bind_audit_theta(
    parts = parts,
    population_model = model,
    factor_set = factor_set,
    theta_proposal = theta_proposal,
    n_cores = n_cores
  )
}

.fit_one_predictive_reference_local <- function(i,
                                                data_list,
                                                loglik_fn,
                                                population_model,
                                                theta_proposal,
                                                local_control,
                                                seed,
                                                verbose) {
  refs <- make_population_predictive_reference_priors(
    population_model = population_model,
    theta_proposal = theta_proposal,
    n_strata = local_control$reference_strata,
    n_components = local_control$predictive_components,
    component_scale = local_control$component_scale,
    defensive_weight = local_control$defensive_weight,
    defensive_scale = local_control$defensive_scale,
    seed = as.integer(seed + 101L * i),
    label = sprintf("q0_predictive_alpha_reference_local%d", i)
  )
  if (is.function(local_control$extra_reference_prior_fn)) {
    extra_refs <- local_control$extra_reference_prior_fn(
      local_id = i,
      data_i = data_list[[i]],
      population_model = population_model,
      theta_proposal = theta_proposal
    )
    if (!is.null(extra_refs)) {
      if (!is.list(extra_refs) || !length(extra_refs)) {
        stop("extra_reference_prior_fn must return NULL or a non-empty list of reference priors.")
      }
      refs <- c(refs, extra_refs)
    }
  }
  component_count <- length(refs)
  components <- lapply(seq_along(refs), function(k) {
    fit_local_reference_smc_component(
      local_id = i,
      data_i = data_list[[i]],
      loglik_fn = loglik_fn,
      reference_prior = refs[[k]],
      M = local_control$M,
      resample_threshold = local_control$resample_threshold,
      n_mcmc_moves = local_control$n_mcmc_moves,
      max_rounds = local_control$max_rounds,
      cess_target = local_control$cess_target,
      G_mix = local_control$G_mix,
      rw_scale_init = local_control$rw_scale_init,
      n_cores = local_control$local_n_cores,
      seed = as.integer(seed + 1009L * i + 100003L * k),
      diagnostics = list(
        role = "q0_predictive_base",
        reference_stratum = as.integer(k),
        reference_strata = as.integer(component_count)
      ),
      verbose = verbose
    )
  })
  object <- reference_local_object_from_components(
    local_id = i,
    components = components
  )
  if (isTRUE(local_control$mbar_enabled)) {
    object <- calibrate_reference_local_object_mbar(
      local_object = object,
      population_model = population_model,
      data_i = data_list[[i]],
      loglik_fn = loglik_fn,
      local_n_cores = local_control$local_n_cores,
      state_counts = local_control$mbar_state_counts,
      max_iter = local_control$mbar_max_iter,
      tol = local_control$mbar_tol,
      offset_drop = local_control$mbar_offset_drop,
      min_state_ess_frac = local_control$mbar_min_state_ess_frac,
      max_abs_shift = local_control$mbar_max_abs_shift,
      require_converged = local_control$mbar_require_converged,
      apply_uncertified = local_control$mbar_apply_uncertified,
      verbose = verbose
    )
    if (!.reference_local_object_mbar_certified(object)) {
      stop(sprintf("Local %d base state bank is not MBAR-certified.", i))
    }
  }
  object
}

.select_predictive_refinement_rows <- function(audit,
                                               max_new_states,
                                               max_new_states_per_local,
                                               mode = c("impact", "protected"),
                                               min_impact_weight = 0,
                                               protected_sources = c("explicit", "q0_axis", "outer_axis", "q0_tail"),
                                               min_protected_impact_weight = 0,
                                               min_surface_score = 0,
                                               source_filter = NULL,
                                               exclude_keys = NULL) {
  if (is.null(audit) || !nrow(audit)) {
    return(audit[0L, , drop = FALSE])
  }
  mode <- match.arg(mode)
  surface_score_all <- if ("surface_score" %in% names(audit)) {
    pmax(as.numeric(audit$surface_score), 0)
  } else {
    rep(0, nrow(audit))
  }
  surface_score_all[!is.finite(surface_score_all)] <- 0
  needs_refinement <- if ("needs_refinement" %in% names(audit)) {
    as.logical(audit$needs_refinement)
  } else {
    !as.logical(audit$covered) | surface_score_all > as.numeric(min_surface_score)
  }
  needs_refinement[is.na(needs_refinement)] <- FALSE
  needs_refinement <- needs_refinement | surface_score_all > as.numeric(min_surface_score)
  failed <- audit[needs_refinement, , drop = FALSE]
  if (!nrow(failed) || as.integer(max_new_states) <= 0L) {
    return(failed[0L, , drop = FALSE])
  }
  if (!is.null(exclude_keys) && length(exclude_keys)) {
    key <- paste(failed$local_id, failed$theta_id, sep = "\r")
    failed <- failed[!(key %in% exclude_keys), , drop = FALSE]
    if (!nrow(failed)) return(failed)
  }

  impact <- if ("theta_impact_weight" %in% names(failed)) {
    pmax(as.numeric(failed$theta_impact_weight), 0)
  } else {
    rep(NA_real_, nrow(failed))
  }
  impact[!is.finite(impact)] <- 0
  surface_score <- if ("surface_score" %in% names(failed)) {
    pmax(as.numeric(failed$surface_score), 0)
  } else {
    impact
  }
  surface_score[!is.finite(surface_score)] <- 0

  source <- as.character(failed$theta_source %||% rep(NA_character_, nrow(failed)))
  protected <- if ("theta_protected" %in% names(failed)) {
    as.logical(failed$theta_protected)
  } else {
    source %in% protected_sources
  }
  protected[is.na(protected)] <- FALSE

  source_rank <- if ("theta_source_rank" %in% names(failed)) {
    as.integer(failed$theta_source_rank)
  } else {
    match(source, protected_sources)
  }
  source_rank[!is.finite(source_rank)] <- 100L

  log_q0 <- if ("theta_log_q0" %in% names(failed)) {
    as.numeric(failed$theta_log_q0)
  } else {
    rep(NA_real_, nrow(failed))
  }
  log_q0[!is.finite(log_q0)] <- -Inf
  posterior_leverage <- if ("posterior_leverage" %in% names(failed)) {
    pmax(as.numeric(failed$posterior_leverage), 0)
  } else {
    rep(1, nrow(failed))
  }
  posterior_leverage[!is.finite(posterior_leverage)] <- 1

  if (identical(mode, "protected")) {
    failed <- failed[protected, , drop = FALSE]
    impact <- impact[protected]
    surface_score <- surface_score[protected]
    source_rank <- source_rank[protected]
    log_q0 <- log_q0[protected]
    posterior_leverage <- posterior_leverage[protected]
    if (!nrow(failed)) return(failed)
    if (!is.null(source_filter) && length(source_filter)) {
      keep <- as.character(failed$theta_source) %in% as.character(source_filter)
      failed <- failed[keep, , drop = FALSE]
      impact <- impact[keep]
      surface_score <- surface_score[keep]
      source_rank <- source_rank[keep]
      log_q0 <- log_q0[keep]
      posterior_leverage <- posterior_leverage[keep]
      if (!nrow(failed)) return(failed)
    }
    floor <- .predictive_reference_source_values(
      min_protected_impact_weight,
      source = as.character(failed$theta_source),
      default = 0
    )
    impact_supported <- impact >= floor |
      surface_score > as.numeric(min_surface_score)
    source_current <- as.character(failed$theta_source)
    coverage_failure_current <- !as.logical(failed$covered)
    coverage_failure_current[is.na(coverage_failure_current)] <- TRUE
    direct_coverage_source <- source_current %in% c("explicit", "outer", "outer_axis")
    keep <- (coverage_failure_current & (direct_coverage_source | impact_supported)) |
      impact_supported
    failed <- failed[keep, , drop = FALSE]
    impact <- impact[keep]
    surface_score <- surface_score[keep]
    source_rank <- source_rank[keep]
    log_q0 <- log_q0[keep]
    posterior_leverage <- posterior_leverage[keep]
    impact_supported <- impact_supported[keep]
    if (!nrow(failed)) return(failed)
  } else {
    impact_supported <- rep(TRUE, nrow(failed))
  }

  if (any(is.finite(impact))) {
    actionable <- if (identical(mode, "impact")) {
      impact >= as.numeric(min_impact_weight) | surface_score > as.numeric(min_surface_score)
    } else {
      impact >= 0 | surface_score > as.numeric(min_surface_score)
    }
    if (!any(actionable)) {
      return(failed[0L, , drop = FALSE])
    }
    failed <- failed[actionable, , drop = FALSE]
    impact <- impact[actionable]
    surface_score <- surface_score[actionable]
    posterior_leverage <- posterior_leverage[actionable]
    impact_supported <- impact_supported[actionable]
  } else {
    impact <- rep(0, nrow(failed))
    surface_score <- rep(0, nrow(failed))
    posterior_leverage <- rep(1, nrow(failed))
    impact_supported <- rep(FALSE, nrow(failed))
  }
  coverage_failure <- !as.logical(failed$covered)
  coverage_failure[is.na(coverage_failure)] <- TRUE
  theta_key <- paste(as.character(failed$theta_source), failed$theta_id, sep = "\r")
  theta_fail_count <- ave(rep.int(1L, nrow(failed)), theta_key, FUN = length)
  theta_min_ess <- ave(as.numeric(failed$ess_frac), theta_key, FUN = function(x) min(x, na.rm = TRUE))
  theta_min_ess[!is.finite(theta_min_ess)] <- Inf
  pareto_rank <- if ("pareto_k" %in% names(failed)) {
    ifelse(is.finite(failed$pareto_k), failed$pareto_k, -Inf)
  } else {
    rep(-Inf, nrow(failed))
  }
  loo_rank <- if ("loo_max_abs_delta" %in% names(failed)) {
    ifelse(is.finite(failed$loo_max_abs_delta), failed$loo_max_abs_delta, -Inf)
  } else {
    rep(-Inf, nrow(failed))
  }
  if (identical(mode, "protected")) {
    failed <- failed[order(-coverage_failure, -impact_supported, -theta_fail_count, theta_min_ess, -impact, -posterior_leverage, -log_q0, -surface_score, source_rank, failed$ess_frac, -pareto_rank, -loo_rank), , drop = FALSE]
  } else {
    failed <- failed[order(-surface_score, -impact, -coverage_failure, -theta_fail_count, theta_min_ess, failed$ess_frac, -pareto_rank, -loo_rank), , drop = FALSE]
  }

  keep <- logical(nrow(failed))
  per_local <- list()
  for (r in seq_len(nrow(failed))) {
    i <- as.character(failed$local_id[r])
    current <- per_local[[i]] %||% 0L
    if (current >= as.integer(max_new_states_per_local)) next
    keep[r] <- TRUE
    per_local[[i]] <- current + 1L
    if (sum(keep) >= as.integer(max_new_states)) break
  }
  failed[keep, , drop = FALSE]
}

.predictive_reference_source_values <- function(x, source, default = 0) {
  source <- as.character(source)
  if (is.null(x) || !length(x)) {
    return(rep(as.numeric(default), length(source)))
  }
  values <- as.numeric(x)
  if (!length(values)) {
    return(rep(as.numeric(default), length(source)))
  }
  x_names <- names(x)
  has_names <- !is.null(x_names) && any(!is.na(x_names) & nzchar(x_names))
  if (!has_names) {
    out <- rep(values[1L], length(source))
  } else {
    out <- values[match(source, x_names)]
    out[is.na(out)] <- as.numeric(default)
  }
  out[!is.finite(out)] <- as.numeric(default)
  out
}

.predictive_reference_protected_source_budgets <- function(audit,
                                                           protected_sources,
                                                           protected_source_weights,
                                                           max_new_states,
                                                           min_protected_impact_weight,
                                                           min_surface_score = 0) {
  budget <- as.integer(max_new_states)
  if (is.null(audit) || !nrow(audit) || budget <= 0L || !length(protected_sources)) {
    return(integer(0))
  }
  source <- as.character(audit$theta_source %||% rep(NA_character_, nrow(audit)))
  covered <- as.logical(audit$covered)
  covered[is.na(covered)] <- FALSE
  protected <- if ("theta_protected" %in% names(audit)) {
    as.logical(audit$theta_protected)
  } else {
    source %in% protected_sources
  }
  protected[is.na(protected)] <- FALSE
  impact <- if ("theta_impact_weight" %in% names(audit)) {
    pmax(as.numeric(audit$theta_impact_weight), 0)
  } else {
    rep(0, nrow(audit))
  }
  impact[!is.finite(impact)] <- 0
  surface_score <- if ("surface_score" %in% names(audit)) {
    pmax(as.numeric(audit$surface_score), 0)
  } else {
    rep(0, nrow(audit))
  }
  surface_score[!is.finite(surface_score)] <- 0
  needs_refinement <- if ("needs_refinement" %in% names(audit)) {
    as.logical(audit$needs_refinement)
  } else {
    !covered | surface_score > as.numeric(min_surface_score)
  }
  needs_refinement[is.na(needs_refinement)] <- FALSE
  needs_refinement <- needs_refinement | surface_score > as.numeric(min_surface_score)
  floor <- .predictive_reference_source_values(
    min_protected_impact_weight,
    source = source,
    default = 0
  )
  active <- protected &
    (
      (!covered & source %in% c("explicit", "outer", "outer_axis")) |
        (needs_refinement &
           (impact >= floor | surface_score > as.numeric(min_surface_score)))
    )
  active_sources <- protected_sources[protected_sources %in% source[active]]
  if (!length(active_sources)) {
    return(integer(0))
  }

  weights <- as.numeric(protected_source_weights)
  names(weights) <- names(protected_source_weights)
  if (!length(weights) || any(!is.finite(weights)) || all(weights <= 0)) {
    weights <- rep(1, length(protected_sources))
    names(weights) <- protected_sources
  }
  source_weight <- weights[active_sources]
  source_weight[!is.finite(source_weight) | source_weight <= 0] <- NA_real_
  source_weight[is.na(source_weight)] <- min(weights[is.finite(weights) & weights > 0], na.rm = TRUE)
  source_weight[!is.finite(source_weight) | source_weight <= 0] <- 1

  raw <- budget * source_weight / sum(source_weight)
  quotas <- pmax(1L, floor(raw))
  names(quotas) <- active_sources
  while (sum(quotas) > budget) {
    idx <- order(quotas, source_weight, decreasing = FALSE)[1L]
    quotas[idx] <- quotas[idx] - 1L
    if (quotas[idx] <= 0L) quotas <- quotas[-idx]
  }
  leftover <- budget - sum(quotas)
  if (leftover > 0L && length(quotas)) {
    frac <- raw[names(quotas)] - floor(raw[names(quotas)])
    for (nm in names(sort(frac, decreasing = TRUE))) {
      if (leftover <= 0L) break
      quotas[nm] <- quotas[nm] + 1L
      leftover <- leftover - 1L
    }
  }
  quotas
}

.select_predictive_protected_refinement_rows <- function(audit,
                                                         max_new_states,
                                                         max_new_states_per_local,
                                                         protected_sources,
                                                         protected_source_weights,
                                                         min_protected_impact_weight,
                                                         min_surface_score = 0,
                                                         exclude_keys = NULL) {
  quotas <- .predictive_reference_protected_source_budgets(
    audit = audit,
    protected_sources = protected_sources,
    protected_source_weights = protected_source_weights,
    max_new_states = max_new_states,
    min_protected_impact_weight = min_protected_impact_weight,
    min_surface_score = min_surface_score
  )
  selected <- data.frame()
  selected_key <- exclude_keys %||% character(0)
  for (src in names(quotas)) {
    remaining <- as.integer(max_new_states) - nrow(selected)
    if (remaining <= 0L) break
    picked <- .select_predictive_refinement_rows(
      audit = audit,
      max_new_states = min(as.integer(quotas[[src]]), remaining),
      max_new_states_per_local = max_new_states_per_local,
      mode = "protected",
      protected_sources = protected_sources,
      min_protected_impact_weight = min_protected_impact_weight,
      min_surface_score = min_surface_score,
      source_filter = src,
      exclude_keys = selected_key
    )
    selected <- .predictive_reference_bind_selected_rows(selected, picked)
    selected_key <- if (nrow(selected)) {
      unique(c(selected_key, paste(selected$local_id, selected$theta_id, sep = "\r")))
    } else {
      selected_key
    }
  }
  remaining <- as.integer(max_new_states) - nrow(selected)
  if (remaining > 0L) {
    fill <- .select_predictive_refinement_rows(
      audit = audit,
      max_new_states = remaining,
      max_new_states_per_local = max_new_states_per_local,
      mode = "protected",
      protected_sources = protected_sources,
      min_protected_impact_weight = min_protected_impact_weight,
      min_surface_score = min_surface_score,
      exclude_keys = selected_key
    )
    selected <- .predictive_reference_bind_selected_rows(selected, fill)
  }
  selected
}

.predictive_reference_state_acceptance <- function(old_row,
                                                   new_audit,
                                                   audit_control) {
  old_covered <- isTRUE(old_row$covered[1L])
  new_covered <- isTRUE(new_audit$covered)
  old_ess <- as.numeric(old_row$ess_frac[1L] %||% NA_real_)
  new_ess <- as.numeric(new_audit$ess_frac %||% NA_real_)
  old_uncertainty <- as.numeric(old_row$surface_uncertainty[1L] %||% NA_real_)
  new_uncertainty <- as.numeric(new_audit$surface_uncertainty %||% NA_real_)
  min_ess <- as.numeric(audit_control$min_ess_frac %||% 0.05)
  reject_uncertainty_ratio <- as.numeric(audit_control$state_reject_uncertainty_ratio %||% 2)
  reject_ess_ratio <- as.numeric(audit_control$state_reject_ess_ratio %||% 0.8)

  if (!old_covered && new_covered) {
    return(list(accepted = TRUE, reason = "covered_trigger"))
  }
  if (!old_covered && is.finite(new_ess) && new_ess >= min_ess &&
      (!is.finite(old_uncertainty) || !is.finite(new_uncertainty) || new_uncertainty <= old_uncertainty)) {
    return(list(accepted = TRUE, reason = "ess_reached_floor"))
  }
  if (is.finite(old_uncertainty) && is.finite(new_uncertainty) &&
      is.finite(old_ess) && is.finite(new_ess) && old_ess > 0 &&
      new_uncertainty > reject_uncertainty_ratio * old_uncertainty &&
      new_ess < reject_ess_ratio * old_ess) {
    return(list(accepted = FALSE, reason = "state_degraded_trigger"))
  }
  list(accepted = TRUE, reason = "mbar_certified")
}

.refine_one_predictive_reference_local <- function(i,
                                                   local_object,
                                                   rows_i,
                                                   theta_audit,
                                                   data_list,
                                                   loglik_fn,
                                                   population_model,
                                                   local_control,
                                                   audit_control,
                                                   seed,
                                                   verbose) {
  object <- validate_reference_local_object(local_object)
  added_rows <- vector("list", 0L)
  anchor_scales <- as.numeric(audit_control$theta_anchor_scales %||%
    audit_control$theta_anchor_scale %||% 1)
  anchor_scales <- anchor_scales[is.finite(anchor_scales) & anchor_scales > 0]
  if (!length(anchor_scales)) anchor_scales <- 1

  rejection_reason <- function(candidate) {
    mbar <- candidate$diagnostics$mbar %||% NULL
    if (!isTRUE(local_control$mbar_enabled) || is.null(mbar)) {
      return("accepted")
    }
    if (isTRUE(mbar$certified)) {
      return("accepted")
    }
    sprintf(
      "mbar_uncertified:min_state_ess=%.4g,max_shift=%.4g,converged=%s",
      as.numeric(mbar$min_state_ess_frac %||% NA_real_),
      as.numeric(mbar$max_abs_shift %||% NA_real_),
      if (isTRUE(mbar$converged)) "true" else "false"
    )
  }

  for (r in seq_len(nrow(rows_i))) {
    theta_id <- as.integer(rows_i$theta_id[r])
    theta <- theta_audit[theta_id, , drop = FALSE]
    for (scale in anchor_scales) {
      ref <- make_population_theta_reference_prior(
        population_model = population_model,
        theta = theta,
        scale = scale,
        label = sprintf("theta_anchor_local%d_theta%d_scale%.3g", i, theta_id, scale)
      )
      component <- fit_local_reference_smc_component(
        local_id = i,
        data_i = data_list[[i]],
        loglik_fn = loglik_fn,
        reference_prior = ref,
        M = audit_control$adapt_particles,
        resample_threshold = local_control$resample_threshold,
        n_mcmc_moves = local_control$n_mcmc_moves,
        max_rounds = local_control$max_rounds,
        cess_target = local_control$cess_target,
        G_mix = local_control$G_mix,
        rw_scale_init = local_control$rw_scale_init,
        n_cores = local_control$local_n_cores,
        seed = as.integer(seed + 1000003L * i + 1009L * theta_id + r + 7919L * round(scale)),
        diagnostics = list(
          role = "theta_anchor_refinement",
          theta_id = theta_id,
          theta_anchor_scale = as.numeric(scale),
          trigger_ess_frac = as.numeric(rows_i$ess_frac[r]),
          trigger_pareto_k = as.numeric(rows_i$pareto_k[r] %||% NA_real_),
          trigger_loo_max_abs_delta = as.numeric(rows_i$loo_max_abs_delta[r] %||% NA_real_)
        ),
        verbose = verbose
      )

      candidate <- object
      candidate$components <- c(candidate$components, list(component))
      candidate$mixture_weights <- .local_predictive_normalize_weights(
        vapply(candidate$components, function(component) nrow(component$particles), numeric(1))
      )
      candidate <- validate_reference_local_object(candidate)
      if (isTRUE(local_control$mbar_enabled)) {
        candidate <- calibrate_reference_local_object_mbar(
          local_object = candidate,
          population_model = population_model,
          data_i = data_list[[i]],
          loglik_fn = loglik_fn,
          local_n_cores = local_control$local_n_cores,
          state_counts = local_control$mbar_state_counts,
          max_iter = local_control$mbar_max_iter,
          tol = local_control$mbar_tol,
          offset_drop = local_control$mbar_offset_drop,
          min_state_ess_frac = local_control$mbar_min_state_ess_frac,
          max_abs_shift = local_control$mbar_max_abs_shift,
          require_converged = local_control$mbar_require_converged,
          apply_uncertified = local_control$mbar_apply_uncertified,
          verbose = verbose
        )
      }

      mbar_accepted <- !isTRUE(local_control$mbar_enabled) ||
        .reference_local_object_mbar_certified(candidate)
      candidate_audit <- NULL
      state_acceptance <- list(accepted = FALSE, reason = rejection_reason(candidate))
      if (isTRUE(mbar_accepted)) {
        candidate_audit <- audit_reference_local_object_support(
          local_object = candidate,
          theta = theta,
          population_model = population_model,
          data_i = data_list[[i]],
          loglik_fn = loglik_fn,
          local_n_cores = local_control$local_n_cores,
          min_ess_frac = audit_control$min_ess_frac,
          compute_psis = audit_control$compute_psis,
          compute_loo = audit_control$compute_loo
        )
        state_acceptance <- .predictive_reference_state_acceptance(
          old_row = rows_i[r, , drop = FALSE],
          new_audit = candidate_audit,
          audit_control = audit_control
        )
      }
      accepted <- isTRUE(mbar_accepted) && isTRUE(state_acceptance$accepted)
      mbar <- candidate$diagnostics$mbar %||% list()
      added_rows[[length(added_rows) + 1L]] <- cbind(
        data.frame(
          local_id = as.integer(i),
          theta_id = theta_id,
          theta_source = as.character(rows_i$theta_source[r] %||% NA_character_),
          theta_protected = isTRUE(rows_i$theta_protected[r] %||% FALSE),
          theta_impact_weight = as.numeric(rows_i$theta_impact_weight[r] %||% NA_real_),
          theta_log_q0 = as.numeric(rows_i$theta_log_q0[r] %||% NA_real_),
          old_ess_frac = as.numeric(rows_i$ess_frac[r]),
          old_pareto_k = as.numeric(rows_i$pareto_k[r] %||% NA_real_),
          old_loo_max_abs_delta = as.numeric(rows_i$loo_max_abs_delta[r] %||% NA_real_),
          old_log_marginal_mcse = as.numeric(rows_i$log_marginal_mcse[r] %||% NA_real_),
          old_split_log_marginal_delta = as.numeric(rows_i$split_log_marginal_delta[r] %||% NA_real_),
          old_component_log_marginal_delta = as.numeric(rows_i$component_log_marginal_delta[r] %||% NA_real_),
          old_gradient_split_delta = as.numeric(rows_i$gradient_split_delta[r] %||% NA_real_),
          old_surface_uncertainty = as.numeric(rows_i$surface_uncertainty[r] %||% NA_real_),
          old_posterior_leverage = as.numeric(rows_i$posterior_leverage[r] %||% NA_real_),
          old_surface_score = as.numeric(rows_i$surface_score[r] %||% NA_real_),
          check.names = FALSE
        ),
        as.data.frame(theta, check.names = FALSE),
        data.frame(
          theta_anchor_scale = as.numeric(scale),
          candidate_component = as.integer(length(candidate$components)),
          candidate_raw_log_evidence = as.numeric(component$raw_log_evidence %||% component$log_evidence),
          candidate_log_evidence = as.numeric(component$log_evidence),
          accepted = isTRUE(accepted),
          rejection_reason = if (isTRUE(accepted)) "accepted" else as.character(state_acceptance$reason),
          candidate_ess_frac = as.numeric(candidate_audit$ess_frac %||% NA_real_),
          candidate_surface_uncertainty = as.numeric(candidate_audit$surface_uncertainty %||% NA_real_),
          candidate_log_marginal_mcse = as.numeric(candidate_audit$log_marginal_mcse %||% NA_real_),
          candidate_split_log_marginal_delta = as.numeric(candidate_audit$split_log_marginal_delta %||% NA_real_),
          candidate_component_log_marginal_delta = as.numeric(candidate_audit$component_log_marginal_delta %||% NA_real_),
          candidate_gradient_split_z = as.numeric(candidate_audit$gradient_split_z %||% NA_real_),
          mbar_converged = isTRUE(mbar$converged),
          mbar_min_state_ess_frac = as.numeric(mbar$min_state_ess_frac %||% NA_real_),
          mbar_max_abs_shift = as.numeric(mbar$max_abs_shift %||% NA_real_),
          new_component = if (isTRUE(accepted)) as.integer(length(candidate$components)) else NA_integer_,
          check.names = FALSE
        )
      )
      if (isTRUE(accepted)) {
        object <- candidate
        break
      }
    }
  }
  object <- validate_reference_local_object(object)
  if (isTRUE(local_control$mbar_enabled) &&
      !.reference_local_object_mbar_certified(object)) {
    stop(sprintf("Local %d state bank is not MBAR-certified after refinement.", i))
  }
  list(
    local_id = as.integer(i),
    local_object = object,
    refinements = if (length(added_rows)) do.call(rbind, added_rows) else data.frame()
  )
}

.refine_predictive_reference_local_objects <- function(local_objects,
                                                       selected,
                                                       theta_audit,
                                                       data_list,
                                                       loglik_fn,
                                                       population_model,
                                                       local_control,
                                                       audit_control,
                                                       n_jobs,
                                                       seed,
                                                       verbose) {
  if (is.null(selected) || !nrow(selected)) {
    return(list(local_objects = local_objects, refinements = data.frame()))
  }
  by_local <- split(selected, selected$local_id)
  parts <- parallel::mclapply(
    names(by_local),
    function(local_name) {
      i <- as.integer(local_name)
      .refine_one_predictive_reference_local(
        i = i,
        local_object = local_objects[[i]],
        rows_i = by_local[[local_name]],
        theta_audit = theta_audit,
        data_list = data_list,
        loglik_fn = loglik_fn,
        population_model = population_model,
        local_control = local_control,
        audit_control = audit_control,
        seed = seed,
        verbose = verbose
      )
    },
    mc.cores = as.integer(max(1L, n_jobs))
  )
  out <- local_objects
  refinement_parts <- vector("list", length(parts))
  for (k in seq_along(parts)) {
    out[[parts[[k]]$local_id]] <- parts[[k]]$local_object
    refinement_parts[[k]] <- parts[[k]]$refinements
  }
  list(
    local_objects = out,
    refinements = if (length(refinement_parts)) do.call(rbind, refinement_parts) else data.frame()
  )
}

.predictive_reference_bind_selected_rows <- function(...) {
  parts <- list(...)
  parts <- parts[vapply(parts, function(x) !is.null(x) && nrow(x) > 0L, logical(1))]
  if (!length(parts)) {
    return(data.frame())
  }
  out <- do.call(rbind, parts)
  key <- paste(out$local_id, out$theta_id, sep = "\r")
  out[!duplicated(key), , drop = FALSE]
}

.predictive_reference_adaptive_outer_proposal <- function(base_proposal,
                                                          previous_fit,
                                                          population_model,
                                                          audit_control,
                                                          seed,
                                                          label = "adaptive_outer_initial") {
  if (!isTRUE(audit_control$adaptive_outer_proposal) || is.null(previous_fit)) {
    return(base_proposal)
  }
  model <- normalize_population_model(population_model)
  theta <- .as_hyper_matrix(previous_fit$theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  w <- pmax(as.numeric(previous_fit$w), 0)
  ok <- apply(theta, 1L, function(x) all(is.finite(x))) & is.finite(w) & w > 0
  if (sum(ok) <= model$hyper_dim + 1L) {
    return(base_proposal)
  }
  w <- w[ok] / sum(w[ok])
  posterior_proposal <- fit_theta_q0_proposal(
    theta = theta[ok, , drop = FALSE],
    log_weight = log(w),
    population_model = model,
    max_components = audit_control$outer_proposal_max_components,
    core_weight = audit_control$outer_proposal_core_weight,
    tail_weight = audit_control$outer_proposal_tail_weight,
    prior_weight = 0,
    df = audit_control$outer_proposal_df,
    tail_df = audit_control$outer_proposal_tail_df,
    core_scale = audit_control$outer_proposal_core_scale,
    tail_scale = audit_control$outer_proposal_tail_scale,
    cluster_sample_size = audit_control$outer_proposal_cluster_sample_size,
    min_component_weight = audit_control$outer_proposal_min_component_weight,
    seed = seed,
    label = paste0(label, "_posterior")
  )
  defensive_weight <- min(max(as.numeric(audit_control$outer_proposal_defensive_weight), 0), 1)
  combine_theta_q0_proposals(
    proposals = list(defensive = base_proposal, posterior = posterior_proposal),
    weights = c(defensive = defensive_weight, posterior = 1 - defensive_weight),
    population_model = model,
    label = label
  )
}

fit_predictive_reference_population_model <- function(data_list,
                                                      loglik_fn,
                                                      population_model,
                                                      theta_proposal,
                                                      local_control = list(),
                                                      outer_control = list(),
                                                      audit_control = list(),
                                                      n_cores = 1L,
                                                      seed = 123L,
                                                      verbose = TRUE) {
  model <- normalize_population_model(population_model)
  theta_proposal <- normalize_theta_proposal(theta_proposal, population_model = model)

  local_defaults <- list(
    M = 600L,
    reference_strata = 1L,
    predictive_components = 64L,
    component_scale = 1,
    defensive_weight = 0.05,
    defensive_scale = 16,
    resample_threshold = 0.6,
    n_mcmc_moves = 2L,
    max_rounds = 100L,
    cess_target = 0.9,
    G_mix = 8L,
    rw_scale_init = 0.9,
    local_n_cores = 1L,
    extra_reference_prior_fn = NULL,
    mbar_enabled = TRUE,
    mbar_state_counts = "ess",
    mbar_max_iter = 1000L,
    mbar_tol = 1e-8,
    mbar_offset_drop = 50,
    mbar_min_state_ess_frac = 0.01,
    mbar_max_abs_shift = 25,
    mbar_require_converged = FALSE,
    mbar_apply_uncertified = FALSE
  )
  outer_defaults <- list(
    N = 1200L,
    resample_threshold = 0.5,
    n_mcmc_moves = 3L,
    min_mcmc_moves = 1L,
    max_rounds = 80L,
    beta_target = 1,
    rw_scale_init = 0.8,
    n_cores = n_cores,
    seed = seed + 700000L,
    verbose = verbose
  )
  audit_defaults <- list(
    enabled = TRUE,
    theta = NULL,
    n_q0 = 64L,
    n_outer = 64L,
    n_tail = 32L,
    n_q0_axis = 32L,
    n_outer_axis = 32L,
    axis_tail_probs = c(0.005, 0.01, 0.025, 0.05, 0.95, 0.975, 0.99, 0.995),
    outer_axis_inflate = 1.5,
    min_ess_frac = 0.05,
    max_pareto_k = 0.7,
    max_loo_delta = Inf,
    compute_psis = TRUE,
    compute_loo = FALSE,
    n_jobs = n_cores,
    local_n_cores = 1L,
    adaptive_refinement = TRUE,
    pre_outer_support_refinement = FALSE,
    pre_outer_support_rounds = 1L,
    max_pre_outer_states = NULL,
    max_pre_outer_states_per_local = NULL,
    max_adapt_rounds = 3L,
    max_new_states = 45L,
    max_new_states_per_local = 3L,
    protected_refinement_fraction = 0.35,
    protected_sources = c("explicit", "q0_tail", "outer_axis", "q0_axis"),
    protected_source_weights = c(explicit = 1, q0_tail = 3, outer_axis = 2, q0_axis = 1),
    min_protected_impact_weight = c(explicit = 0, q0_tail = 1e-8, outer_axis = 1e-6, q0_axis = 1e-6),
    adapt_particles = local_control$M %||% local_defaults$M,
    theta_anchor_scale = c(1, 4, 16),
    theta_anchor_scales = NULL,
    min_refinement_impact_weight = 1e-4,
    min_surface_score = 1e-3,
    state_reject_uncertainty_ratio = Inf,
    state_reject_ess_ratio = 0,
    repair_enabled = FALSE,
    max_repairs = 0L,
    repair_particles = local_control$M %||% local_defaults$M,
    repair_reference_scale = 1,
    repair_n_mcmc_moves = 1L,
    repair_max_rounds = 80L,
    adaptive_outer_proposal = FALSE,
    outer_proposal_defensive_weight = 0.35,
    outer_proposal_max_components = 4L,
    outer_proposal_core_weight = 0.65,
    outer_proposal_tail_weight = 0.35,
    outer_proposal_df = 7,
    outer_proposal_tail_df = 3,
    outer_proposal_core_scale = 1.5,
    outer_proposal_tail_scale = 16,
    outer_proposal_cluster_sample_size = 2000L,
    outer_proposal_min_component_weight = 0.02
  )

  local_control <- .predictive_reference_merge_control(local_control, local_defaults)
  outer_control <- .predictive_reference_merge_control(outer_control, outer_defaults)
  audit_control <- .predictive_reference_merge_control(audit_control, audit_defaults)
  audit_control$max_pre_outer_states <- as.integer(
    audit_control$max_pre_outer_states %||% audit_control$max_new_states
  )
  audit_control$max_pre_outer_states_per_local <- as.integer(
    audit_control$max_pre_outer_states_per_local %||% audit_control$max_new_states_per_local
  )

  if (isTRUE(verbose)) {
    cat(sprintf(
      "Fitting q0-predictive locals: %d locals x %d particles | q0 components=%d\n",
      length(data_list),
      as.integer(local_control$M),
      as.integer(local_control$predictive_components)
    ))
  }

  local_objects <- parallel::mclapply(
    seq_along(data_list),
    function(i) {
      .fit_one_predictive_reference_local(
        i = i,
        data_list = data_list,
        loglik_fn = loglik_fn,
        population_model = model,
        theta_proposal = theta_proposal,
        local_control = local_control,
        seed = as.integer(seed),
        verbose = FALSE
      )
    },
    mc.cores = as.integer(max(1L, n_cores))
  )
  names(local_objects) <- names(data_list) %||% paste0("local_", seq_along(data_list))

  build_factor_set_from_local_objects <- function(local_objects_current) {
    build_population_factor_set(
      local_objects = local_objects_current,
      population_model = model,
      data_list = data_list,
      loglik_fn = loglik_fn,
      local_n_cores = audit_control$local_n_cores
    )
  }

  run_outer_from_local_objects <- function(local_objects_current,
                                           seed_offset = 0L,
                                           initial_proposal = theta_proposal) {
    factor_set_current <- build_factor_set_from_local_objects(local_objects_current)
    outer_args <- modifyList(
      list(
        factor_set = factor_set_current,
        initial_proposal = initial_proposal
      ),
      outer_control
    )
    outer_args$seed <- as.integer(outer_control$seed + seed_offset)
    list(
      factor_set = factor_set_current,
      fit = do.call(outer_population_smc, outer_args)
    )
  }

  theta_audit <- NULL
  audit <- NULL
  audit_history <- list()
  refinement_history <- data.frame()
  repair_result <- NULL
  repair_audit <- NULL

  factor_set <- build_factor_set_from_local_objects(local_objects)
  if (isTRUE(audit_control$enabled) &&
      isTRUE(audit_control$pre_outer_support_refinement) &&
      as.integer(audit_control$pre_outer_support_rounds) > 0L) {
    for (pre_round in seq_len(as.integer(audit_control$pre_outer_support_rounds))) {
      theta_audit <- .predictive_reference_audit_theta(
        theta_proposal = theta_proposal,
        outer_fit = NULL,
        factor_set = factor_set,
        population_model = model,
        n_q0 = audit_control$n_q0,
        n_outer = 0L,
        n_tail = audit_control$n_tail,
        n_q0_axis = audit_control$n_q0_axis,
        n_outer_axis = 0L,
        axis_tail_probs = audit_control$axis_tail_probs,
        outer_axis_inflate = audit_control$outer_axis_inflate,
        explicit_theta = audit_control$theta,
        n_cores = audit_control$n_jobs,
        seed = seed + 760000L + 1009L * pre_round
      )
      audit <- audit_reference_local_objects(
        local_objects = local_objects,
        theta_audit = theta_audit,
        population_model = model,
        data_list = data_list,
        loglik_fn = loglik_fn,
        min_ess_frac = audit_control$min_ess_frac,
        max_pareto_k = audit_control$max_pareto_k,
        max_loo_delta = audit_control$max_loo_delta,
        compute_psis = audit_control$compute_psis,
        compute_loo = audit_control$compute_loo,
        n_jobs = audit_control$n_jobs,
        local_n_cores = audit_control$local_n_cores
      )
      audit$adapt_round <- -as.integer(pre_round)
      audit$adapt_phase <- "pre_outer_support"
      audit_history[[length(audit_history) + 1L]] <- audit
      selected <- .select_predictive_protected_refinement_rows(
        audit = audit,
        max_new_states = audit_control$max_pre_outer_states,
        max_new_states_per_local = audit_control$max_pre_outer_states_per_local,
        protected_sources = audit_control$protected_sources,
        protected_source_weights = audit_control$protected_source_weights,
        min_protected_impact_weight = audit_control$min_protected_impact_weight,
        min_surface_score = audit_control$min_surface_score
      )
      if (isTRUE(verbose)) {
        cat(sprintf(
          "Pre-outer support audit %d: failures=%d/%d | selected protected states=%d | min ESS=%.4f | median ESS=%.4f\n",
          pre_round,
          sum(!audit$covered),
          nrow(audit),
          nrow(selected),
          min(audit$ess_frac, na.rm = TRUE),
          stats::median(audit$ess_frac, na.rm = TRUE)
        ))
      }
      if (!nrow(selected)) break
      refinement <- .refine_predictive_reference_local_objects(
        local_objects = local_objects,
        selected = selected,
        theta_audit = theta_audit,
        data_list = data_list,
        loglik_fn = loglik_fn,
        population_model = model,
        local_control = local_control,
        audit_control = audit_control,
        n_jobs = audit_control$n_jobs,
        seed = seed + 780000L + 10007L * pre_round,
        verbose = FALSE
      )
      local_objects <- refinement$local_objects
      if (nrow(refinement$refinements)) {
        refinement$refinements$adapt_round <- -as.integer(pre_round)
        refinement$refinements$adapt_phase <- "pre_outer_support"
        refinement_history <- rbind(refinement_history, refinement$refinements)
      }
      factor_set <- build_factor_set_from_local_objects(local_objects)
    }
  }

  outer_run <- run_outer_from_local_objects(local_objects, seed_offset = 0L)
  factor_set <- outer_run$factor_set
  fit <- outer_run$fit

  if (isTRUE(audit_control$enabled)) {
    max_rounds <- if (isTRUE(audit_control$adaptive_refinement)) {
      as.integer(max(0L, audit_control$max_adapt_rounds))
    } else {
      0L
    }
    for (adapt_round in seq.int(0L, max_rounds)) {
      theta_audit <- .predictive_reference_audit_theta(
        theta_proposal = theta_proposal,
        outer_fit = fit,
        factor_set = factor_set,
        population_model = model,
        n_q0 = audit_control$n_q0,
        n_outer = audit_control$n_outer,
        n_tail = audit_control$n_tail,
        n_q0_axis = audit_control$n_q0_axis,
        n_outer_axis = audit_control$n_outer_axis,
        axis_tail_probs = audit_control$axis_tail_probs,
        outer_axis_inflate = audit_control$outer_axis_inflate,
        explicit_theta = audit_control$theta,
        n_cores = audit_control$n_jobs,
        seed = seed + 800000L + 1009L * adapt_round
      )
      audit <- audit_reference_local_objects(
        local_objects = local_objects,
        theta_audit = theta_audit,
        population_model = model,
        data_list = data_list,
        loglik_fn = loglik_fn,
        min_ess_frac = audit_control$min_ess_frac,
        max_pareto_k = audit_control$max_pareto_k,
        max_loo_delta = audit_control$max_loo_delta,
        compute_psis = audit_control$compute_psis,
        compute_loo = audit_control$compute_loo,
        n_jobs = audit_control$n_jobs,
        local_n_cores = audit_control$local_n_cores
      )
      audit$adapt_round <- as.integer(adapt_round)
      audit$adapt_phase <- "posterior"
      audit_history[[length(audit_history) + 1L]] <- audit
      surface_needs_refinement <- (audit$surface_score %||% 0) >
        as.numeric(audit_control$min_surface_score)
      needs_refinement <- (audit$needs_refinement %||% !audit$covered) |
        surface_needs_refinement
      if (isTRUE(verbose)) {
        cat(sprintf(
          "Local state-bank audit round %d: failures=%d/%d | surface candidates=%d | max surface score=%.4g | min ESS=%.4f | median ESS=%.4f\n",
          adapt_round,
          sum(!audit$covered),
          nrow(audit),
          sum(needs_refinement),
          max(audit$surface_score %||% 0, na.rm = TRUE),
          min(audit$ess_frac, na.rm = TRUE),
          stats::median(audit$ess_frac, na.rm = TRUE)
        ))
      }
      if (adapt_round >= max_rounds || !any(needs_refinement)) {
        break
      }
      protected_fraction <- min(max(as.numeric(audit_control$protected_refinement_fraction), 0), 1)
      protected_budget <- as.integer(min(
        audit_control$max_new_states,
        ceiling(protected_fraction * as.integer(audit_control$max_new_states))
      ))
      protected_selected <- .select_predictive_protected_refinement_rows(
        audit = audit,
        max_new_states = protected_budget,
        max_new_states_per_local = audit_control$max_new_states_per_local,
        protected_sources = audit_control$protected_sources,
        protected_source_weights = audit_control$protected_source_weights,
        min_protected_impact_weight = audit_control$min_protected_impact_weight,
        min_surface_score = audit_control$min_surface_score
      )
      selected_key <- if (nrow(protected_selected)) {
        paste(protected_selected$local_id, protected_selected$theta_id, sep = "\r")
      } else {
        character(0)
      }
      impact_budget <- as.integer(audit_control$max_new_states) - nrow(protected_selected)
      impact_selected <- .select_predictive_refinement_rows(
        audit = audit,
        max_new_states = impact_budget,
        max_new_states_per_local = audit_control$max_new_states_per_local,
        mode = "impact",
        min_impact_weight = audit_control$min_refinement_impact_weight,
        min_surface_score = audit_control$min_surface_score,
        exclude_keys = selected_key
      )
      selected <- .predictive_reference_bind_selected_rows(protected_selected, impact_selected)
      if (!nrow(selected)) break
      refinement <- .refine_predictive_reference_local_objects(
        local_objects = local_objects,
        selected = selected,
        theta_audit = theta_audit,
        data_list = data_list,
        loglik_fn = loglik_fn,
        population_model = model,
        local_control = local_control,
        audit_control = audit_control,
        n_jobs = audit_control$n_jobs,
        seed = seed + 850000L + 10007L * adapt_round,
        verbose = FALSE
      )
      local_objects <- refinement$local_objects
      if (nrow(refinement$refinements)) {
        refinement$refinements$adapt_round <- as.integer(adapt_round)
        refinement$refinements$adapt_phase <- "posterior"
        refinement_history <- rbind(refinement_history, refinement$refinements)
      }
      adaptive_initial_proposal <- .predictive_reference_adaptive_outer_proposal(
        base_proposal = theta_proposal,
        previous_fit = fit,
        population_model = model,
        audit_control = audit_control,
        seed = seed + 830000L + 1009L * adapt_round,
        label = sprintf("adaptive_outer_initial_round%d", adapt_round)
      )
      outer_run <- run_outer_from_local_objects(
        local_objects,
        seed_offset = as.integer(100000L + 1009L * adapt_round),
        initial_proposal = adaptive_initial_proposal
      )
      factor_set <- outer_run$factor_set
      fit <- outer_run$fit
    }

    if (isTRUE(audit_control$repair_enabled) &&
        as.integer(audit_control$max_repairs) > 0L &&
        any(!audit$covered)) {
      repair_result <- bridge_repair_reference_local_objects(
        local_objects = local_objects,
        audit = audit,
        theta_audit = theta_audit,
        data_list = data_list,
        loglik_fn = loglik_fn,
        population_model = model,
        max_repairs = audit_control$max_repairs,
        M = audit_control$repair_particles,
        min_ess_frac = audit_control$min_ess_frac,
        reference_scale = audit_control$repair_reference_scale,
        n_mcmc_moves = audit_control$repair_n_mcmc_moves,
        max_rounds = audit_control$repair_max_rounds,
        n_jobs = audit_control$n_jobs,
        local_n_cores = audit_control$local_n_cores,
        seed = seed + 900000L,
        verbose = FALSE
      )
      local_objects <- repair_result$local_objects
      repaired_ids <- sort(unique(repair_result$repairs$local_id))
      if (length(repaired_ids) && isTRUE(local_control$mbar_enabled)) {
        for (i in repaired_ids) {
          local_objects[[i]] <- calibrate_reference_local_object_mbar(
            local_object = local_objects[[i]],
            population_model = model,
            data_i = data_list[[i]],
            loglik_fn = loglik_fn,
            local_n_cores = local_control$local_n_cores,
            state_counts = local_control$mbar_state_counts,
            max_iter = local_control$mbar_max_iter,
            tol = local_control$mbar_tol,
            offset_drop = local_control$mbar_offset_drop,
            min_state_ess_frac = local_control$mbar_min_state_ess_frac,
            max_abs_shift = local_control$mbar_max_abs_shift,
            require_converged = local_control$mbar_require_converged,
            apply_uncertified = local_control$mbar_apply_uncertified,
            verbose = FALSE
          )
        }
      }
      new_factor_set <- build_population_factor_set(
        local_objects = local_objects,
        population_model = model,
        data_list = data_list,
        loglik_fn = loglik_fn,
        local_n_cores = audit_control$local_n_cores
      )
      fit <- update_outer_population_fit(
        fit = fit,
        old_factor_set = factor_set,
        new_factor_set = new_factor_set,
        n_mcmc_moves = outer_control$n_mcmc_moves,
        min_mcmc_moves = outer_control$min_mcmc_moves,
        resample_threshold = outer_control$resample_threshold,
        n_cores = outer_control$n_cores,
        seed = seed + 950000L,
        verbose = verbose
      )
      factor_set <- new_factor_set
      repair_audit <- audit_reference_local_objects(
        local_objects = local_objects,
        theta_audit = theta_audit,
        population_model = model,
        data_list = data_list,
        loglik_fn = loglik_fn,
        local_ids = repaired_ids,
        min_ess_frac = audit_control$min_ess_frac,
        max_pareto_k = audit_control$max_pareto_k,
        max_loo_delta = audit_control$max_loo_delta,
        compute_psis = audit_control$compute_psis,
        compute_loo = audit_control$compute_loo,
        n_jobs = audit_control$n_jobs,
        local_n_cores = audit_control$local_n_cores
      )
      if (isTRUE(verbose)) {
        cat(sprintf(
          "Bridge repair applied to %d local/theta failures; repaired-local remaining failures=%d\n",
          nrow(repair_result$repairs),
          sum(!repair_audit$covered)
        ))
      }
    }
  }

  structure(
    list(
      fit = fit,
      factor_set = factor_set,
      local_objects = local_objects,
      theta_proposal = theta_proposal,
      theta_audit = theta_audit,
      audit = audit,
      audit_history = audit_history,
      refinement_history = refinement_history,
      repair_result = repair_result,
      repair_audit = repair_audit,
      population_model = model,
      settings = list(
        local_control = local_control,
        outer_control = outer_control,
        audit_control = audit_control,
        n_cores = n_cores,
        seed = seed
      )
    ),
    class = "predictive_reference_population_fit"
  )
}
