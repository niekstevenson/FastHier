#!/usr/bin/env Rscript
# ============================================================================
# Outer population SMC
# - Stacked local marginal-likelihood factors only
# - Tempered outer SMC for the diagonal Gaussian hierarchy
# ============================================================================

if (!exists("%||%", mode = "function") ||
    !exists("weighted_cov", mode = "function") ||
    !exists(".rowLogSumExp", mode = "function") ||
    !exists("next_lambda_via_rCESS", mode = "function") ||
    !exists("cess_target_at_lambda", mode = "function") ||
    !exists("stratified_resample_sorted", mode = "function")) {
  source("smc_core.R")
}
if (!exists("normalize_reference_prior", mode = "function") ||
    !exists("reference_prior_logpdf", mode = "function")) {
  source("reference_priors.R")
}
if (!exists("normalize_population_model", mode = "function") ||
    !exists("population_model_prepare_theta", mode = "function") ||
    !exists("population_model_log_alpha_given_prepared_theta_many", mode = "function")) {
  source("population_models.R")
}

suppressPackageStartupMessages({
  library(parallel)
})

normalize_reference_local_weights <- function(w) {
  w <- pmax(as.numeric(w), 0)
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) {
    rep(1 / length(w), length(w))
  } else {
    w / sw
  }
}

validate_reference_local_component <- function(component) {
  if (!is.list(component)) stop("reference local component must be a list.")
  required <- c("particles", "weights", "reference_prior")
  missing <- setdiff(required, names(component))
  if (length(missing)) {
    stop("reference local component is missing: ", paste(missing, collapse = ", "))
  }

  particles <- as.matrix(component$particles)
  weights <- normalize_reference_local_weights(component$weights)
  reference_prior <- normalize_reference_prior(reference_prior = component$reference_prior)
  log_reference_density <- as.numeric(
    component$log_reference_density %||% reference_prior_logpdf(reference_prior, particles)
  )
  proposal_type <- as.character(component$proposal_type %||% "posterior_reference")
  if (!identical(length(proposal_type), 1L) ||
      !(proposal_type %in% c("posterior_reference", "prior_reference"))) {
    stop("reference local component proposal_type must be 'posterior_reference' or 'prior_reference'.")
  }
  log_likelihood <- component$log_likelihood
  if (is.null(log_likelihood)) {
    log_likelihood <- rep(NA_real_, nrow(particles))
  }
  log_likelihood <- as.numeric(log_likelihood)
  if (length(log_likelihood) != nrow(particles)) {
    stop("reference local component log_likelihood must match the number of particles.")
  }
  if (length(log_reference_density) != nrow(particles)) {
    stop("reference local component log_reference_density must match the number of particles.")
  }

  component$particles <- particles
  component$weights <- weights
  component$reference_prior <- reference_prior
  component$log_reference_density <- log_reference_density
  component$proposal_type <- proposal_type
  component$log_likelihood <- log_likelihood
  component$log_evidence <- as.numeric(component$log_evidence %||% NA_real_)
  component$mcse_log_evidence <- as.numeric(component$mcse_log_evidence %||% NA_real_)
  component
}

validate_reference_local_object <- function(local_object) {
  if (!is.list(local_object)) stop("local_object must be a list.")
  local_object$local_id <- as.integer(local_object$local_id %||% NA_integer_)

  if (!is.null(local_object$components)) {
    if (!length(local_object$components)) {
      stop("local_object components must be non-empty.")
    }
    local_object$components <- lapply(local_object$components, validate_reference_local_component)
    mixture_weights <- local_object$mixture_weights %||%
      rep(1 / length(local_object$components), length(local_object$components))
    mixture_weights <- normalize_reference_local_weights(mixture_weights)
    if (length(mixture_weights) != length(local_object$components)) {
      stop("local_object mixture_weights must match the number of components.")
    }
    local_object$mixture_weights <- mixture_weights
    return(local_object)
  }

  required <- c("particles", "weights", "reference_prior")
  missing <- setdiff(required, names(local_object))
  if (length(missing)) {
    stop("local_object is missing: ", paste(missing, collapse = ", "))
  }

  local_object$components <- list(
    validate_reference_local_component(
      list(
        particles = local_object$particles,
        weights = local_object$weights,
        log_evidence = local_object$log_evidence,
        mcse_log_evidence = local_object$mcse_log_evidence,
        reference_prior = local_object$reference_prior,
        log_reference_density = local_object$log_reference_density,
        proposal_type = local_object$proposal_type %||% "posterior_reference",
        log_likelihood = local_object$log_likelihood,
        diagnostics = local_object$diagnostics %||% list()
      )
    )
  )
  local_object$mixture_weights <- 1
  local_object
}

.is_compressed_population_local_object <- function(local_object) {
  is.list(local_object) &&
    (!is.null(local_object$factor_particles) || !is.null(local_object$factor_sufficient_stats)) &&
    !is.null(local_object$factor_log_base)
}

population_sufficient_stats_from_alpha <- function(population_model, alpha) {
  model <- normalize_population_model(population_model)
  alpha <- as.matrix(alpha)
  if (ncol(alpha) != model$alpha_dim) {
    stop("Sufficient-statistic alpha dimension does not match the population model.")
  }
  if (!is.null(colnames(alpha)) && setequal(colnames(alpha), model$alpha_names)) {
    alpha <- alpha[, model$alpha_names, drop = FALSE]
  } else {
    colnames(alpha) <- model$alpha_names
  }

  if (identical(model$fast_family %||% NULL, "gaussian")) {
    return(structure(
      list(
        family = "gaussian_diag",
        linear = alpha,
        quadratic_diag = alpha * alpha,
        n = as.integer(nrow(alpha)),
        alpha_names = model$alpha_names
      ),
      class = "population_sufficient_stats"
    ))
  }

  structure(
    list(
      family = "raw_alpha",
      alpha = alpha,
      n = as.integer(nrow(alpha)),
      alpha_names = model$alpha_names
    ),
    class = "population_sufficient_stats"
  )
}

population_sufficient_stats_n <- function(sufficient_stats) {
  as.integer(sufficient_stats$n %||%
    if (!is.null(sufficient_stats$linear)) nrow(sufficient_stats$linear) else nrow(sufficient_stats$alpha))
}

validate_population_sufficient_stats <- function(sufficient_stats, population_model) {
  model <- normalize_population_model(population_model)
  if (!is.list(sufficient_stats) || is.null(sufficient_stats$family)) {
    stop("population sufficient statistics must be a list with a family.")
  }
  family <- as.character(sufficient_stats$family)

  if (identical(family, "gaussian_diag")) {
    linear <- as.matrix(sufficient_stats$linear)
    quadratic_diag <- as.matrix(sufficient_stats$quadratic_diag)
    if (ncol(linear) != model$alpha_dim ||
        ncol(quadratic_diag) != model$alpha_dim ||
        nrow(linear) != nrow(quadratic_diag)) {
      stop("Gaussian sufficient statistics have incompatible dimensions.")
    }
    colnames(linear) <- model$alpha_names
    colnames(quadratic_diag) <- model$alpha_names
    sufficient_stats$linear <- linear
    sufficient_stats$quadratic_diag <- quadratic_diag
    sufficient_stats$n <- as.integer(nrow(linear))
    sufficient_stats$alpha_names <- model$alpha_names
    return(structure(sufficient_stats, class = "population_sufficient_stats"))
  }

  if (identical(family, "raw_alpha")) {
    return(population_sufficient_stats_from_alpha(model, sufficient_stats$alpha))
  }

  stop("Unknown population sufficient-statistic family: ", family)
}

population_sufficient_stats_subset <- function(sufficient_stats, idx) {
  idx <- as.integer(idx)
  family <- as.character(sufficient_stats$family)
  if (identical(family, "gaussian_diag")) {
    out <- sufficient_stats
    out$linear <- sufficient_stats$linear[idx, , drop = FALSE]
    out$quadratic_diag <- sufficient_stats$quadratic_diag[idx, , drop = FALSE]
    out$n <- as.integer(length(idx))
    return(out)
  }
  if (identical(family, "raw_alpha")) {
    out <- sufficient_stats
    out$alpha <- sufficient_stats$alpha[idx, , drop = FALSE]
    out$n <- as.integer(length(idx))
    return(out)
  }
  stop("Unknown population sufficient-statistic family: ", family)
}

population_sufficient_stats_bind <- function(population_model, stats_list) {
  model <- normalize_population_model(population_model)
  if (!length(stats_list)) {
    alpha <- matrix(numeric(0), nrow = 0L, ncol = model$alpha_dim)
    colnames(alpha) <- model$alpha_names
    return(population_sufficient_stats_from_alpha(model, alpha))
  }
  stats_list <- lapply(stats_list, validate_population_sufficient_stats, population_model = model)
  families <- vapply(stats_list, `[[`, character(1), "family")
  if (length(unique(families)) != 1L) {
    stop("Cannot bind mixed sufficient-statistic families.")
  }

  if (identical(unname(families[[1L]]), "gaussian_diag")) {
    linear <- do.call(rbind, lapply(stats_list, `[[`, "linear"))
    quadratic_diag <- do.call(rbind, lapply(stats_list, `[[`, "quadratic_diag"))
    colnames(linear) <- model$alpha_names
    colnames(quadratic_diag) <- model$alpha_names
    return(structure(
      list(
        family = "gaussian_diag",
        linear = linear,
        quadratic_diag = quadratic_diag,
        n = as.integer(nrow(linear)),
        alpha_names = model$alpha_names
      ),
      class = "population_sufficient_stats"
    ))
  }

  alpha <- do.call(rbind, lapply(stats_list, `[[`, "alpha"))
  population_sufficient_stats_from_alpha(model, alpha)
}

population_log_alpha_given_sufficient_stats_many <- function(model,
                                                             sufficient_stats,
                                                             theta = NULL,
                                                             theta_prepared = NULL) {
  model <- normalize_population_model(model)
  sufficient_stats <- validate_population_sufficient_stats(sufficient_stats, model)
  if (is.null(theta_prepared)) {
    theta_prepared <- population_model_prepare_theta(model, theta)
  }

  if (identical(sufficient_stats$family, "gaussian_diag") &&
      identical(theta_prepared$family, "gaussian")) {
    if (!identical(theta_prepared$quadratic_kind, "diag")) {
      stop("Only diagonal Gaussian sufficient-statistic evaluation is supported.")
    }
    linear <- theta_prepared$eta %*% t(sufficient_stats$linear)
    quad <- theta_prepared$quadratic_coef %*% t(sufficient_stats$quadratic_diag)
    return(sweep(linear - 0.5 * quad, 1L, theta_prepared$log_kernel_constant, "+"))
  }

  if (!identical(sufficient_stats$family, "raw_alpha")) {
    stop("Sufficient-statistic family is incompatible with the population model.")
  }
  population_model_log_alpha_given_prepared_theta_many(
    model,
    alpha = sufficient_stats$alpha,
    theta_prepared = theta_prepared
  )
}

build_compressed_population_local_factor <- function(local_object,
                                                     population_model) {
  population_model <- normalize_population_model(population_model)
  alpha <- NULL
  if (!is.null(local_object$factor_particles)) {
    alpha <- as.matrix(local_object$factor_particles)
    if (ncol(alpha) != population_model$alpha_dim) {
      stop("Compressed local factor particle dimension does not match the population model.")
    }
    if (!is.null(colnames(alpha)) && setequal(colnames(alpha), population_model$alpha_names)) {
      alpha <- alpha[, population_model$alpha_names, drop = FALSE]
    } else {
      colnames(alpha) <- population_model$alpha_names
    }
  }

  sufficient_stats <- if (!is.null(local_object$factor_sufficient_stats)) {
    validate_population_sufficient_stats(local_object$factor_sufficient_stats, population_model)
  } else {
    population_sufficient_stats_from_alpha(population_model, alpha)
  }
  n_particles <- population_sufficient_stats_n(sufficient_stats)
  if (!is.null(alpha) && nrow(alpha) != n_particles) {
    stop("Compressed local factor particles and sufficient statistics have different lengths.")
  }

  log_base <- as.numeric(local_object$factor_log_base)
  if (length(log_base) != n_particles) {
    stop("Compressed local factor log_base must match the number of particles.")
  }

  structure(
    list(
      local_id = as.integer(local_object$local_id %||% NA_integer_),
      population_model = population_model,
      particles = alpha,
      sufficient_stats = sufficient_stats,
      n_particles = n_particles,
      log_weights = rep(NA_real_, n_particles),
      log_reference_density = rep(NA_real_, n_particles),
      log_base = log_base,
      log_constant = as.numeric(local_object$factor_log_constant %||% 0),
      reference_prior = NULL,
      component_id = rep.int(1L, n_particles)
    ),
    class = "population_local_factor"
  )
}

build_population_local_factor <- function(local_object,
                                          population_model,
                                          data_i = NULL,
                                          loglik_fn = NULL,
                                          local_n_cores = 1L) {
  if (.is_compressed_population_local_object(local_object)) {
    return(build_compressed_population_local_factor(local_object, population_model))
  }

  local_object <- validate_reference_local_object(local_object)
  population_model <- normalize_population_model(population_model)
  components <- local_object$components

  align_alpha <- function(alpha) {
    alpha <- as.matrix(alpha)
    if (ncol(alpha) != population_model$alpha_dim) {
      stop("Local particle dimension does not match the population model.")
    }
    if (!is.null(colnames(alpha)) && setequal(colnames(alpha), population_model$alpha_names)) {
      alpha[, population_model$alpha_names, drop = FALSE]
    } else {
      colnames(alpha) <- population_model$alpha_names
      alpha
    }
  }

  proposal_types <- vapply(components, `[[`, character(1), "proposal_type")

  if (length(components) == 1L && identical(proposal_types, "posterior_reference")) {
    component <- components[[1L]]
    alpha <- align_alpha(component$particles)
    log_weights <- log(component$weights)
    log_reference_density <- as.numeric(
      component$log_reference_density %||% reference_prior_logpdf(component$reference_prior, alpha)
    )
    log_base <- log_weights - log_reference_density
    log_constant <- as.numeric(component$log_evidence %||% 0)
    component_id <- rep.int(1L, nrow(alpha))
    reference_prior <- component$reference_prior
  } else {
    eta <- normalize_reference_local_weights(local_object$mixture_weights)
    log_eta <- ifelse(eta > 0, log(eta), -Inf)
    alpha_list <- lapply(components, function(component) align_alpha(component$particles))
    alpha <- do.call(rbind, alpha_list)
    component_lengths <- vapply(alpha_list, nrow, integer(1))
    component_id <- rep.int(seq_along(alpha_list), component_lengths)

    cached_loglik <- unlist(lapply(components, `[[`, "log_likelihood"), use.names = FALSE)
    if (length(cached_loglik) != nrow(alpha) || anyNA(cached_loglik)) {
      if (is.null(data_i) || is.null(loglik_fn)) {
        stop("Multi-component local factors require data_i and loglik_fn unless every component stores log_likelihood.")
      }
      cached_loglik <- ll_parallel(
        alpha,
        data_i,
        loglik_fn,
        n_cores = as.integer(max(1L, local_n_cores))
      )
    }
    log_likelihood <- as.numeric(cached_loglik)
    log_likelihood[is.na(log_likelihood)] <- -Inf

    logq_mat <- vapply(seq_along(components), function(j) {
      component <- components[[j]]
      as.numeric(reference_prior_logpdf(component$reference_prior, alpha))
    }, numeric(nrow(alpha)))
    if (!is.matrix(logq_mat)) {
      logq_mat <- matrix(logq_mat, ncol = length(components))
    }

    component_logZ <- vapply(components, function(component) {
      as.numeric(component$log_evidence %||% NA_real_)
    }, numeric(1))
    needs_logZ <- proposal_types == "posterior_reference"
    if (any(needs_logZ & !is.finite(component_logZ))) {
      stop("Posterior-reference DMIS components require finite log_evidence.")
    }

    logh_mat <- logq_mat
    if (any(needs_logZ)) {
      logh_mat[, needs_logZ] <- sweep(
        logh_mat[, needs_logZ, drop = FALSE],
        1L,
        log_likelihood,
        "+"
      )
      logh_mat[, needs_logZ] <- sweep(
        logh_mat[, needs_logZ, drop = FALSE],
        2L,
        component_logZ[needs_logZ],
        "-"
      )
    }
    log_hmix <- .rowLogSumExp(sweep(logh_mat, 2L, log_eta, "+"))
    log_sample_weights <- unlist(lapply(components, function(component) log(component$weights)), use.names = FALSE)
    log_weights <- log_sample_weights + log_eta[component_id]
    log_base <- log_weights + log_likelihood - log_hmix
    log_base[!is.finite(log_base)] <- -Inf
    log_constant <- 0
    log_reference_density <- rep(NA_real_, nrow(alpha))
    reference_prior <- NULL
  }

  structure(
    list(
      local_id = as.integer(local_object$local_id),
      population_model = population_model,
      particles = alpha,
      sufficient_stats = population_sufficient_stats_from_alpha(population_model, alpha),
      n_particles = as.integer(nrow(alpha)),
      log_weights = log_weights,
      log_reference_density = log_reference_density,
      log_base = log_base,
      log_constant = log_constant,
      reference_prior = reference_prior,
      component_id = component_id
    ),
    class = "population_local_factor"
  )
}

.gaussian_log_alpha_given_prepared_theta_many <- function(alpha, theta_prepared) {
  alpha <- as.matrix(alpha)
  if (!identical(theta_prepared$quadratic_kind, "diag")) {
    stop("Only diagonal Gaussian population models are supported.")
  }
  linear <- theta_prepared$eta %*% t(alpha)
  quad_t <- t(alpha * alpha)
  quad <- theta_prepared$quadratic_coef %*% quad_t
  sweep(linear - 0.5 * quad, 1L, theta_prepared$log_kernel_constant, "+")
}

.population_log_alpha_given_theta_many <- function(model,
                                                   alpha,
                                                   theta = NULL,
                                                   theta_prepared = NULL) {
  if (is.null(theta_prepared)) {
    theta_prepared <- population_model_prepare_theta(model, theta)
  }
  if (identical(theta_prepared$family, "gaussian")) {
    return(.gaussian_log_alpha_given_prepared_theta_many(alpha, theta_prepared))
  }
  population_model_log_alpha_given_prepared_theta_many(model, alpha = alpha, theta_prepared = theta_prepared)
}

.local_log_marginal_blocked <- function(model,
                                        sufficient_stats,
                                        log_base,
                                        theta,
                                        theta_prepared = NULL,
                                        block_size = 1024L) {
  theta <- .as_hyper_matrix(theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  theta_prepared <- theta_prepared %||% population_model_prepare_theta(model, theta)
  sufficient_stats <- validate_population_sufficient_stats(sufficient_stats, model)
  block_size <- as.integer(max(1L, block_size))
  accum <- rep.int(-Inf, nrow(theta))
  n_alpha <- population_sufficient_stats_n(sufficient_stats)

  for (start in seq.int(1L, n_alpha, by = block_size)) {
    idx <- seq.int(start, min(start + block_size - 1L, n_alpha))
    log_terms <- population_log_alpha_given_sufficient_stats_many(
      model = model,
      sufficient_stats = population_sufficient_stats_subset(sufficient_stats, idx),
      theta_prepared = theta_prepared
    )
    block_lse <- .rowLogSumExp(sweep(log_terms, 2L, log_base[idx], "+"))
    accum <- rlogsumexp2(accum, block_lse)
  }
  as.numeric(accum)
}

population_local_factor_log_marginal_many <- function(factor,
                                                      theta,
                                                      include_constant = TRUE,
                                                      block_size = 1024L) {
  stopifnot(inherits(factor, "population_local_factor"))
  theta_prepared <- population_model_prepare_theta(factor$population_model, theta)
  out <- .local_log_marginal_blocked(
    model = factor$population_model,
    sufficient_stats = factor$sufficient_stats,
    log_base = factor$log_base,
    theta = theta_prepared$theta,
    theta_prepared = theta_prepared,
    block_size = block_size
  )
  if (isTRUE(include_constant)) out <- out + factor$log_constant
  as.numeric(out)
}

population_local_factor_log_marginal <- function(factor, theta, include_constant = TRUE) {
  population_local_factor_log_marginal_many(factor, theta = theta, include_constant = include_constant)[1L]
}

population_local_factor_ess <- function(factor, theta) {
  stopifnot(inherits(factor, "population_local_factor"))
  theta_prepared <- population_model_prepare_theta(factor$population_model, theta)
  logp <- population_log_alpha_given_sufficient_stats_many(
    model = factor$population_model,
    sufficient_stats = factor$sufficient_stats,
    theta_prepared = theta_prepared
  )[1L, ]
  lw <- factor$log_base + logp
  lse <- logsumexp(lw)
  if (!is.finite(lse)) return(0)
  w <- exp(lw - lse)
  as.numeric(1 / sum(w * w))
}

population_local_factor_reweighted_particles <- function(factor, theta) {
  stopifnot(inherits(factor, "population_local_factor"))
  if (is.null(factor$particles)) {
    stop("Raw factor particles are not stored for this sufficient-statistic factor.")
  }
  theta_prepared <- population_model_prepare_theta(factor$population_model, theta)
  logp <- population_log_alpha_given_sufficient_stats_many(
    model = factor$population_model,
    sufficient_stats = factor$sufficient_stats,
    theta_prepared = theta_prepared
  )[1L, ]
  lw <- factor$log_base + logp
  lse <- logsumexp(lw)
  weights <- if (is.finite(lse)) {
    exp(lw - lse)
  } else {
    rep(1 / nrow(factor$particles), nrow(factor$particles))
  }
  list(particles = factor$particles, weights = weights)
}

population_local_factor_tail_diagnostic <- function(factor, theta, use_psis = TRUE) {
  stopifnot(inherits(factor, "population_local_factor"))
  theta_prepared <- population_model_prepare_theta(factor$population_model, theta)
  logp <- population_log_alpha_given_sufficient_stats_many(
    model = factor$population_model,
    sufficient_stats = factor$sufficient_stats,
    theta_prepared = theta_prepared
  )[1L, ]
  lw <- factor$log_base + logp
  keep <- is.finite(lw)
  if (!any(keep)) {
    return(list(ess = 0, pareto_k = NA_real_, max_weight = NA_real_, q99_weight = NA_real_, n_finite = 0L))
  }

  lw <- lw[keep]
  lse <- logsumexp(lw)
  if (!is.finite(lse)) {
    return(list(ess = 0, pareto_k = NA_real_, max_weight = NA_real_, q99_weight = NA_real_, n_finite = length(lw)))
  }

  w <- exp(lw - lse)
  pareto_k <- NA_real_
  if (isTRUE(use_psis) && length(lw) >= 5L && requireNamespace("loo", quietly = TRUE)) {
    psis_obj <- tryCatch(
      suppressWarnings(loo::psis(matrix(lw - max(lw), ncol = 1L))),
      error = function(e) NULL
    )
    if (!is.null(psis_obj)) {
      pareto_k <- as.numeric(loo::pareto_k_values(psis_obj)[1L])
    }
  }

  list(
    ess = as.numeric(1 / sum(w * w)),
    pareto_k = pareto_k,
    max_weight = as.numeric(max(w)),
    q99_weight = as.numeric(stats::quantile(w, probs = 0.99, na.rm = TRUE, names = FALSE)),
    n_finite = length(lw)
  )
}

.population_factor_set_from_factors <- function(factors,
                                                population_model,
                                                particle_block_size = 1024L) {
  population_model <- normalize_population_model(population_model)
  particle_block_size <- as.integer(max(1L, particle_block_size))

  if (length(factors)) {
    local_lengths <- vapply(factors, `[[`, integer(1), "n_particles")
    sufficient_stats <- population_sufficient_stats_bind(
      population_model,
      lapply(factors, `[[`, "sufficient_stats")
    )
    has_all_particles <- all(vapply(factors, function(factor) !is.null(factor$particles), logical(1)))
    alpha <- if (isTRUE(has_all_particles)) {
      do.call(rbind, lapply(factors, `[[`, "particles"))
    } else {
      NULL
    }
    log_base <- unlist(lapply(factors, `[[`, "log_base"), use.names = FALSE)
    local_index <- rep.int(seq_along(factors), local_lengths)
    block_starts <- seq.int(1L, population_sufficient_stats_n(sufficient_stats), by = particle_block_size)
    blocks <- lapply(block_starts, function(start) {
      idx <- seq.int(start, min(start + particle_block_size - 1L, population_sufficient_stats_n(sufficient_stats)))
      rr <- rle(local_index[idx])
      ends <- cumsum(rr$lengths)
      starts <- c(1L, head(ends, -1L) + 1L)
      list(
        idx = idx,
        local = as.integer(rr$values),
        starts = as.integer(starts),
        ends = as.integer(ends)
      )
    })
  } else {
    alpha <- matrix(numeric(0), nrow = 0L, ncol = population_model$alpha_dim)
    colnames(alpha) <- population_model$alpha_names
    sufficient_stats <- population_sufficient_stats_from_alpha(population_model, alpha)
    log_base <- numeric(0)
    local_index <- integer(0)
    blocks <- list()
  }

  structure(
    list(
      population_model = population_model,
      factors = factors,
      n_locals = length(factors),
      log_constant = sum(vapply(factors, `[[`, numeric(1), "log_constant")),
      stack = list(
        alpha = alpha,
        sufficient_stats = sufficient_stats,
        log_base = log_base,
        local_index = local_index,
        blocks = blocks,
        particle_block_size = particle_block_size
      )
    ),
    class = "population_factor_set"
  )
}

.data_for_local_object <- function(local_object, pos, data_list, local_names = NULL) {
  if (is.null(data_list)) return(NULL)
  local_id <- as.integer(local_object$local_id %||% NA_integer_)
  if (is.finite(local_id) && local_id >= 1L && local_id <= length(data_list)) {
    return(data_list[[local_id]])
  }
  if (!is.null(local_names) &&
      !is.na(local_names[pos]) &&
      nzchar(local_names[pos]) &&
      local_names[pos] %in% names(data_list)) {
    return(data_list[[local_names[pos]]])
  }
  data_list[[pos]]
}

build_population_factor_set <- function(local_objects,
                                        population_model,
                                        particle_block_size = 1024L,
                                        data_list = NULL,
                                        loglik_fn = NULL,
                                        local_n_cores = 1L) {
  population_model <- normalize_population_model(population_model)
  local_names <- names(local_objects)
  factors <- lapply(seq_along(local_objects), function(i) {
    build_population_local_factor(
      local_objects[[i]],
      population_model = population_model,
      data_i = .data_for_local_object(local_objects[[i]], i, data_list, local_names),
      loglik_fn = loglik_fn,
      local_n_cores = local_n_cores
    )
  })
  names(factors) <- names(local_objects)
  .population_factor_set_from_factors(
    factors = factors,
    population_model = population_model,
    particle_block_size = particle_block_size
  )
}

update_population_factor_set_locals <- function(factor_set,
                                                local_objects,
                                                local_ids,
                                                data_list = NULL,
                                                loglik_fn = NULL,
                                                local_n_cores = 1L) {
  stopifnot(inherits(factor_set, "population_factor_set"))
  local_ids <- sort(unique(as.integer(local_ids)))
  local_ids <- local_ids[local_ids >= 1L & local_ids <= length(local_objects)]
  if (!length(local_ids)) return(factor_set)

  factors <- factor_set$factors
  local_names <- names(local_objects)
  for (i in local_ids) {
    factors[[i]] <- build_population_local_factor(
      local_objects[[i]],
      population_model = factor_set$population_model,
      data_i = .data_for_local_object(local_objects[[i]], i, data_list, local_names),
      loglik_fn = loglik_fn,
      local_n_cores = local_n_cores
    )
  }
  names(factors) <- names(local_objects)
  .population_factor_set_from_factors(
    factors = factors,
    population_model = factor_set$population_model,
    particle_block_size = factor_set$stack$particle_block_size
  )
}

.evaluate_factor_set_stacked <- function(factor_set, theta, include_constant = FALSE) {
  stopifnot(inherits(factor_set, "population_factor_set"))
  model <- factor_set$population_model
  theta <- .as_hyper_matrix(theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  if (!length(factor_set$factors)) {
    out <- rep.int(0, nrow(theta))
    if (isTRUE(include_constant)) out <- out + factor_set$log_constant
    return(out)
  }

  theta_prepared <- population_model_prepare_theta(model, theta)
  accum <- matrix(-Inf, nrow = nrow(theta), ncol = factor_set$n_locals)

  for (block in factor_set$stack$blocks) {
    log_terms <- population_log_alpha_given_sufficient_stats_many(
      model = model,
      sufficient_stats = population_sufficient_stats_subset(
        factor_set$stack$sufficient_stats,
        block$idx
      ),
      theta_prepared = theta_prepared
    )
    log_terms <- sweep(log_terms, 2L, factor_set$stack$log_base[block$idx], "+")

    for (seg_id in seq_along(block$local)) {
      cols <- seq.int(block$starts[seg_id], block$ends[seg_id])
      seg_lse <- if (length(cols) == 1L) {
        log_terms[, cols]
      } else {
        .rowLogSumExp(log_terms[, cols, drop = FALSE])
      }
      accum[, block$local[seg_id]] <- rlogsumexp2(accum[, block$local[seg_id]], as.numeric(seg_lse))
    }
  }

  out <- rowSums(accum)
  if (isTRUE(include_constant)) out <- out + factor_set$log_constant
  as.numeric(out)
}

population_factor_set_loglik <- function(factor_set, theta, include_constant = FALSE, n_cores = 1L) {
  stopifnot(inherits(factor_set, "population_factor_set"))
  model <- factor_set$population_model
  theta <- .as_hyper_matrix(theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)

  if (!length(factor_set$factors)) {
    out <- rep.int(0, nrow(theta))
    if (isTRUE(include_constant)) out <- out + factor_set$log_constant
    return(out)
  }

  if (as.integer(n_cores) <= 1L || nrow(theta) <= 1L) {
    return(.evaluate_factor_set_stacked(factor_set, theta, include_constant = include_constant))
  }

  split_idx <- cut(seq_len(nrow(theta)), breaks = as.integer(min(n_cores, nrow(theta))), labels = FALSE)
  batches <- split(seq_len(nrow(theta)), split_idx)
  parts <- parallel::mclapply(
    batches,
    function(idx) {
      .evaluate_factor_set_stacked(
        factor_set = factor_set,
        theta = theta[idx, , drop = FALSE],
        include_constant = FALSE
      )
    },
    mc.cores = as.integer(min(n_cores, length(batches)))
  )

  out <- numeric(nrow(theta))
  for (i in seq_along(batches)) {
    out[batches[[i]]] <- parts[[i]]
  }
  if (isTRUE(include_constant)) out <- out + factor_set$log_constant
  out
}

population_factor_set_loglik_by_local <- function(factor_set,
                                                  theta,
                                                  include_constant = TRUE,
                                                  n_cores = 1L) {
  stopifnot(inherits(factor_set, "population_factor_set"))
  model <- factor_set$population_model
  theta <- .as_hyper_matrix(theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  if (!length(factor_set$factors)) {
    return(matrix(numeric(0), nrow = nrow(theta), ncol = 0L))
  }

  eval_one <- function(factor) {
    population_local_factor_log_marginal_many(
      factor,
      theta = theta,
      include_constant = include_constant
    )
  }

  parts <- if (as.integer(n_cores) <= 1L || length(factor_set$factors) <= 1L) {
    lapply(factor_set$factors, eval_one)
  } else {
    parallel::mclapply(
      factor_set$factors,
      eval_one,
      mc.cores = as.integer(min(n_cores, length(factor_set$factors)))
    )
  }
  out <- do.call(cbind, parts)
  colnames(out) <- names(factor_set$factors) %||% paste0("local_", seq_along(factor_set$factors))
  out
}

population_factor_set_local_ess <- function(factor_set, theta) {
  stopifnot(inherits(factor_set, "population_factor_set"))
  vapply(factor_set$factors, population_local_factor_ess, numeric(1), theta = theta)
}

population_factor_set_logposterior <- function(factor_set, theta, include_constant = TRUE, n_cores = 1L) {
  stopifnot(inherits(factor_set, "population_factor_set"))
  model <- factor_set$population_model
  theta <- .as_hyper_matrix(theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  population_model_log_hyperprior(model, theta) +
    population_factor_set_loglik(factor_set, theta, include_constant = include_constant, n_cores = n_cores)
}

.outer_cheap_sort_order <- function(Z) {
  Z <- as.matrix(Z)
  n <- nrow(Z)
  if (n <= 1L) return(seq_len(n))
  if (ncol(Z) <= 1L) return(order(Z[, 1L]))
  order(rowSums(Z))
}

.outer_whiten_theta <- function(theta, w = NULL) {
  theta <- as.matrix(theta)
  n <- nrow(theta)
  d <- ncol(theta)
  if (n <= 1L || d <= 0L) return(theta)

  if (is.null(w)) {
    w <- rep(1 / n, n)
  } else {
    w <- pmax(as.numeric(w), 0)
    sw <- sum(w)
    w <- if (!is.finite(sw) || sw <= 0) rep(1 / n, n) else w / sw
  }

  mu <- colSums(theta * w)
  S <- tryCatch(weighted_cov(theta, w), error = function(e) NULL)
  if (is.null(S) || any(!is.finite(S))) {
    S <- tryCatch(stats::cov(theta), error = function(e) NULL)
  }
  if (is.null(S) || any(!is.finite(S))) {
    S <- diag(d)
  } else {
    S <- regularize_cov(S, min_eig = 1e-8, cond_cap = 1e8)
  }

  L <- tryCatch(chol(S), error = function(e) diag(d))
  theta_centered <- sweep(theta, 2L, mu, "-")
  Z <- t(backsolve(L, t(theta_centered), transpose = TRUE))
  Z[!is.finite(Z)] <- 0
  Z
}

.outer_resample_sort_order <- function(theta,
                                       w = NULL,
                                       ess_frac = 1,
                                       mode = c("adaptive", "hilbert", "cheap1d", "none"),
                                       hilbert_hard_ess = 0.25,
                                       hilbert_max_dim = 8L) {
  mode <- match.arg(mode)
  theta <- as.matrix(theta)
  if (nrow(theta) <= 1L) {
    return(list(order = seq_len(nrow(theta)), sort_used = "none"))
  }
  if (mode == "none") {
    return(list(order = seq_len(nrow(theta)), sort_used = "none"))
  }

  Z <- .outer_whiten_theta(theta, w = w)
  want_hilbert <- identical(mode, "hilbert") ||
    (identical(mode, "adaptive") &&
       ess_frac <= as.numeric(hilbert_hard_ess) &&
       ncol(Z) <= as.integer(hilbert_max_dim))

  if (isTRUE(want_hilbert)) {
    ord <- tryCatch(hilbert_sort_order(Z, bits = 16L), error = function(e) NULL)
    if (!is.null(ord)) {
      return(list(order = ord, sort_used = "hilbert"))
    }
  }

  list(order = .outer_cheap_sort_order(Z), sort_used = "cheap1d")
}

.outer_rejuvenate <- function(theta,
                              logprior,
                              loglik_dynamic,
                              beta,
                              factor_set,
                              w = NULL,
                              min_n_moves = 1L,
                              max_n_moves = 3L,
                              target_accepted_moves = NULL,
                              target_accepted_moves_resampled = 1.5,
                              target_accepted_moves_not_resampled = 0.75,
                              rw_scale = 0.8,
                              target_accept = 0.234,
                              rm_gain_move = 0.15,
                              n_cores = 1L,
                              resampled = FALSE,
                              seed = NULL) {
  if (!is.null(seed)) set.seed(as.integer(seed))
  theta <- as.matrix(theta)
  N <- nrow(theta)
  d <- ncol(theta)
  model <- factor_set$population_model
  if (N <= 1L) {
    return(list(
      theta = theta,
      logprior = logprior,
      loglik_dynamic = loglik_dynamic,
      accept_rate = 0,
      move_accept = numeric(0),
      n_moves_used = 0L,
      cumulative_accept = 0,
      rw_scale = as.numeric(rw_scale)
    ))
  }

  if (is.null(w)) {
    w_cov <- rep(1 / N, N)
  } else {
    w_cov <- pmax(as.numeric(w), 0)
    sw <- sum(w_cov)
    w_cov <- if (!is.finite(sw) || sw <= 0) rep(1 / N, N) else w_cov / sw
  }

  S <- tryCatch(weighted_cov(theta, w_cov), error = function(e) NULL)
  if (is.null(S) || any(!is.finite(S))) {
    S <- tryCatch(stats::cov(theta), error = function(e) NULL)
  }
  if (is.null(S) || any(!is.finite(S))) {
    S <- diag(d)
  } else {
    S <- regularize_cov(S, min_eig = 1e-8, cond_cap = 1e8)
  }
  Lrw <- tryCatch(chol(S), error = function(e) diag(d))
  log_rw_scale <- log(pmax(as.numeric(rw_scale), 1e-6))
  min_n_moves <- as.integer(max(1L, min_n_moves))
  max_n_moves <- as.integer(max(min_n_moves, max_n_moves))
  target_accepted_moves <- as.numeric(
    target_accepted_moves %||%
      if (isTRUE(resampled)) target_accepted_moves_resampled else target_accepted_moves_not_resampled
  )

  acc_total <- 0L
  prop_total <- 0L
  logpost <- logprior + beta * loglik_dynamic
  move_accept <- numeric(0)
  cumulative_accept <- 0

  for (move_id in seq_len(max_n_moves)) {
    step_scale <- exp(log_rw_scale) / sqrt(max(d, 1))
    theta_prop <- theta + step_scale * (matrix(rnorm(N * d), nrow = N, ncol = d) %*% Lrw)
    colnames(theta_prop) <- model$hyper_names

    logprior_prop <- population_model_log_hyperprior(model, theta_prop)
    ok <- is.finite(logprior_prop)
    loglik_prop <- rep(-Inf, N)
    if (any(ok)) {
      loglik_prop[ok] <- population_factor_set_loglik(
        factor_set,
        theta = theta_prop[ok, , drop = FALSE],
        include_constant = FALSE,
        n_cores = n_cores
      )
    }
    logpost_prop <- logprior_prop + beta * loglik_prop
    log_alpha <- logpost_prop - logpost
    accept <- which(log(runif(N)) < pmin(0, log_alpha))
    prop_total <- prop_total + N
    accept_rate_move <- length(accept) / N
    move_accept <- c(move_accept, accept_rate_move)
    cumulative_accept <- cumulative_accept + accept_rate_move
    log_rw_scale <- .clamp(
      log_rw_scale + as.numeric(rm_gain_move) * (accept_rate_move - as.numeric(target_accept)),
      log(0.05),
      log(4.0)
    )

    if (length(accept)) {
      theta[accept, ] <- theta_prop[accept, , drop = FALSE]
      logprior[accept] <- logprior_prop[accept]
      loglik_dynamic[accept] <- loglik_prop[accept]
      logpost[accept] <- logpost_prop[accept]
      acc_total <- acc_total + length(accept)
    }

    if (move_id >= min_n_moves && cumulative_accept >= target_accepted_moves) {
      break
    }
  }

  list(
    theta = theta,
    logprior = logprior,
    loglik_dynamic = loglik_dynamic,
    accept_rate = if (prop_total > 0L) acc_total / prop_total else 0,
    move_accept = move_accept,
    n_moves_used = length(move_accept),
    cumulative_accept = cumulative_accept,
    rw_scale = exp(log_rw_scale)
  )
}

update_outer_population_fit <- function(fit,
                                        old_factor_set,
                                        new_factor_set,
                                        n_mcmc_moves = 3L,
                                        min_mcmc_moves = 1L,
                                        resample_threshold = 0.5,
                                        rw_scale = NULL,
                                        n_cores = 1L,
                                        seed = 123L,
                                        verbose = TRUE) {
  stopifnot(inherits(old_factor_set, "population_factor_set"))
  stopifnot(inherits(new_factor_set, "population_factor_set"))
  theta <- as.matrix(fit$theta)
  w_old <- normalize_reference_local_weights(fit$w)
  old_loglik <- as.numeric(fit$loglik_dynamic %||%
    population_factor_set_loglik(old_factor_set, theta, include_constant = FALSE, n_cores = n_cores))
  new_loglik <- population_factor_set_loglik(
    new_factor_set,
    theta = theta,
    include_constant = FALSE,
    n_cores = n_cores
  )

  log_ratio <- new_loglik - old_loglik
  ok <- is.finite(log_ratio)
  if (!any(ok)) stop("Population factor update produced no finite correction weights.")
  log_ratio[!ok] <- min(log_ratio[ok])
  logw <- log(pmax(w_old, .Machine$double.eps)) + log_ratio
  lse <- logsumexp(logw)
  w <- exp(logw - lse)
  w <- w / sum(w)
  ess_frac <- ESS(w) / length(w)
  resampled <- FALSE

  model <- new_factor_set$population_model
  logprior <- as.numeric(fit$logprior %||% population_model_log_hyperprior(model, theta))
  if (ess_frac < as.numeric(resample_threshold)) {
    sort_info <- .outer_resample_sort_order(theta = theta, w = w, ess_frac = ess_frac)
    ord <- sort_info$order
    idx <- ord[stratified_resample_sorted(w[ord])]
    theta <- theta[idx, , drop = FALSE]
    logprior <- logprior[idx]
    new_loglik <- new_loglik[idx]
    w <- rep(1 / nrow(theta), nrow(theta))
    resampled <- TRUE
  }

  rejuvenated <- .outer_rejuvenate(
    theta = theta,
    logprior = logprior,
    loglik_dynamic = new_loglik,
    beta = 1,
    factor_set = new_factor_set,
    w = w,
    min_n_moves = min_mcmc_moves,
    max_n_moves = n_mcmc_moves,
    rw_scale = as.numeric(rw_scale %||% tail(fit$meta$rw_scale_hist %||% 0.8, 1L)),
    n_cores = n_cores,
    resampled = resampled,
    seed = seed
  )
  if (isTRUE(verbose)) {
    cat(sprintf(
      "Population factor update: ESS=%.3f | rejuvenation accept=%.3f%s\n",
      ess_frac,
      rejuvenated$accept_rate,
      if (resampled) " | resampled" else ""
    ))
  }

  fit$theta <- rejuvenated$theta
  fit$w <- if (isTRUE(resampled)) {
    rep(1 / nrow(rejuvenated$theta), nrow(rejuvenated$theta))
  } else {
    w
  }
  fit$logprior <- rejuvenated$logprior
  fit$loglik_dynamic <- rejuvenated$loglik_dynamic
  fit$population_model <- model
  fit$log_evidence_dynamic <- as.numeric(fit$log_evidence_dynamic %||% 0) + lse
  fit$log_evidence_constant <- new_factor_set$log_constant
  fit$log_evidence <- fit$log_evidence_dynamic + new_factor_set$log_constant
  fit$meta$factor_update <- c(
    fit$meta$factor_update %||% list(),
    list(list(
      ess = ess_frac,
      resampled = resampled,
      accept_rate = rejuvenated$accept_rate,
      moves = rejuvenated$n_moves_used
    ))
  )
  fit
}

outer_population_smc <- function(factor_set,
                                 N = 2000L,
                                 resample_threshold = 0.5,
                                 n_mcmc_moves = 3L,
                                 min_mcmc_moves = 1L,
                                 max_rounds = 80L,
                                 beta_target = 1.0,
                                 rw_scale_init = 0.8,
                                 target_accept = 0.234,
                                 rm_gain = 0.15,
                                 target_accepted_moves = NULL,
                                 target_accepted_moves_resampled = 1.5,
                                 target_accepted_moves_not_resampled = 0.75,
                                 deterministic_resampling = FALSE,
                                 resample_sort_mode = c("adaptive", "hilbert", "cheap1d", "none"),
                                 hilbert_hard_ess = 0.25,
                                 hilbert_max_dim = 8L,
                                 n_cores = 1L,
                                 seed = 123L,
                                 verbose = TRUE) {
  stopifnot(inherits(factor_set, "population_factor_set"))
  model <- factor_set$population_model
  vcat <- function(...) if (isTRUE(verbose)) base::cat(...)

  set.seed(as.integer(seed))
  beta_target <- as.numeric(beta_target)
  if (!is.finite(beta_target) || beta_target <= 0 || beta_target > 1) {
    stop("beta_target must be in (0, 1].")
  }
  resample_sort_mode <- match.arg(resample_sort_mode)
  min_mcmc_moves <- as.integer(max(1L, min_mcmc_moves))
  n_mcmc_moves <- as.integer(max(min_mcmc_moves, n_mcmc_moves))

  theta <- population_model_sample_hyper(model, n = as.integer(N))
  colnames(theta) <- model$hyper_names
  logprior <- population_model_log_hyperprior(model, theta)
  loglik_dynamic <- population_factor_set_loglik(factor_set, theta, include_constant = FALSE, n_cores = n_cores)

  w <- rep(1 / N, N)
  beta <- 0.0
  round <- 0L
  log_evidence_dynamic <- 0.0
  mcse_var_accum <- 0.0
  log_rw_scale <- log(rw_scale_init)

  beta_hist <- beta
  ess_hist <- 1.0
  accept_hist <- numeric(0)
  rw_scale_hist <- exp(log_rw_scale)
  resampled_hist <- logical(0)
  moves_hist <- integer(0)
  sort_hist <- character(0)

  while (beta < beta_target - 1e-12 && round < as.integer(max_rounds)) {
    round <- round + 1L
    target_cess <- cess_target_at_lambda(beta)
    beta_new <- next_lambda_via_rCESS(
      w = w,
      loglik = loglik_dynamic,
      lambda = beta,
      target = target_cess,
      lambda_target = beta_target
    )
    if (beta_new <= beta) {
      beta_new <- min(beta_target, beta + min(1e-4, beta_target - beta))
    }

    delta <- beta_new - beta
    x <- as.numeric(loglik_dynamic)
    ok <- is.finite(x)
    if (!any(ok)) stop("All outer log-likelihood values are non-finite.")
    x[!ok] <- min(x[ok])
    mx <- max(x)
    u <- exp(delta * (x - mx))
    mu1 <- sum(w * u)
    mu2 <- sum(w * u * u)
    neff <- 1 / sum(w * w)
    mcse_var_accum <- mcse_var_accum +
      max((mu2 - mu1^2) / (max(neff, 1) * max(mu1^2, .Machine$double.eps)), 0)

    logw_raw <- log(pmax(w, .Machine$double.eps)) + delta * (x - mx)
    lse <- logsumexp(logw_raw)
    log_evidence_dynamic <- log_evidence_dynamic + delta * mx + lse
    w <- exp(logw_raw - lse)
    w <- w / sum(w)
    ess_now <- ESS(w) / length(w)
    resampled <- FALSE

    vcat(sprintf(
      "\nOuter round %d: beta %.3f -> %.3f | ESS=%.3f | logZ_dyn=%.4f\n",
      round, beta, beta_new, ess_now, log_evidence_dynamic
    ))

    beta <- beta_new
    sort_used <- "not_used"
    if (ess_now < resample_threshold) {
      sort_info <- .outer_resample_sort_order(
        theta = theta,
        w = w,
        ess_frac = ess_now,
        mode = resample_sort_mode,
        hilbert_hard_ess = hilbert_hard_ess,
        hilbert_max_dim = hilbert_max_dim
      )
      ord <- sort_info$order
      sort_used <- sort_info$sort_used
      idx_sorted <- stratified_resample_sorted(w[ord], deterministic = deterministic_resampling)
      idx <- ord[idx_sorted]
      theta <- theta[idx, , drop = FALSE]
      logprior <- logprior[idx]
      loglik_dynamic <- loglik_dynamic[idx]
      w <- rep(1 / N, N)
      resampled <- TRUE
    }

    rejuvenated <- .outer_rejuvenate(
      theta = theta,
      logprior = logprior,
      loglik_dynamic = loglik_dynamic,
      beta = beta,
      factor_set = factor_set,
      w = w,
      min_n_moves = min_mcmc_moves,
      max_n_moves = n_mcmc_moves,
      target_accepted_moves = target_accepted_moves,
      target_accepted_moves_resampled = target_accepted_moves_resampled,
      target_accepted_moves_not_resampled = target_accepted_moves_not_resampled,
      rw_scale = exp(log_rw_scale),
      target_accept = target_accept,
      rm_gain_move = rm_gain,
      n_cores = n_cores,
      resampled = resampled,
      seed = seed + 1009L * round
    )
    theta <- rejuvenated$theta
    logprior <- rejuvenated$logprior
    loglik_dynamic <- rejuvenated$loglik_dynamic
    accept_rate <- rejuvenated$accept_rate
    log_rw_scale <- log(pmax(rejuvenated$rw_scale, 1e-6))

    beta_hist <- c(beta_hist, beta)
    ess_hist <- c(ess_hist, ess_now)
    accept_hist <- c(accept_hist, accept_rate)
    rw_scale_hist <- c(rw_scale_hist, exp(log_rw_scale))
    resampled_hist <- c(resampled_hist, resampled)
    moves_hist <- c(moves_hist, rejuvenated$n_moves_used)
    sort_hist <- c(sort_hist, sort_used)

    vcat(sprintf(
      "  Rejuvenation accept=%.3f | moves=%d | rw_scale=%.3f%s%s\n",
      accept_rate,
      rejuvenated$n_moves_used,
      exp(log_rw_scale),
      if (resampled) " | resampled" else "",
      if (identical(sort_used, "not_used")) "" else paste0(" | sort=", sort_used)
    ))
  }

  if (beta < beta_target - 1e-12) {
    warning("outer_population_smc hit max_rounds before reaching beta_target.")
  }

  list(
    theta = theta,
    w = w,
    logprior = logprior,
    loglik_dynamic = loglik_dynamic,
    beta = beta,
    log_evidence_dynamic = log_evidence_dynamic,
    log_evidence_constant = factor_set$log_constant,
    log_evidence = log_evidence_dynamic + factor_set$log_constant,
    mcse_log_evidence = sqrt(mcse_var_accum),
    population_model = model,
    meta = list(
      rounds = round,
      beta_hist = beta_hist,
      ess_hist = ess_hist,
      accept_hist = accept_hist,
      rw_scale_hist = rw_scale_hist,
      resampled_hist = resampled_hist,
      moves_hist = moves_hist,
      sort_hist = sort_hist
    )
  )
}
