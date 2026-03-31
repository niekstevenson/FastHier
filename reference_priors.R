#!/usr/bin/env Rscript
# ============================================================================
# Reference prior helpers
# - Gaussian and Gaussian-mixture reference priors for reusable local SMC
# - Shared evaluators and constructors used by the local and hierarchical code
# ============================================================================

if (!exists("%||%", mode = "function") ||
    !exists("dmvnorm_chol_log", mode = "function") ||
    !exists(".rowLogSumExp", mode = "function") ||
    !exists("regularize_cov", mode = "function")) {
  source("smc_core.R")
}

suppressPackageStartupMessages({
  library(mvtnorm)
})

.reference_param_names <- function(mu, Sigma, param_names = NULL) {
  param_names %||% names(mu) %||% colnames(Sigma) %||% paste0("theta", seq_along(mu))
}

.reference_matrix <- function(Theta, d, param_names = NULL) {
  Theta <- as.matrix(Theta)
  if (ncol(Theta) != d) {
    if (nrow(Theta) == d && ncol(Theta) == 1L) {
      Theta <- t(Theta)
    } else if (length(Theta) == d) {
      Theta <- matrix(as.numeric(Theta), nrow = 1L, ncol = d)
    } else {
      stop("Theta has incompatible dimensions for the reference prior.")
    }
  }
  if (!is.null(param_names)) {
    if (!is.null(colnames(Theta)) && setequal(colnames(Theta), param_names)) {
      Theta <- Theta[, param_names, drop = FALSE]
    } else {
      colnames(Theta) <- param_names
    }
  }
  Theta
}

.normalize_reference_component <- function(mu, Sigma, weight = 1, param_names = NULL) {
  mu <- as.numeric(mu)
  Sigma <- regularize_cov(as.matrix(Sigma), min_eig = 1e-8, cond_cap = 1e8)
  d <- length(mu)
  if (nrow(Sigma) != d || ncol(Sigma) != d) {
    stop("Reference prior covariance dimensions do not match the mean.")
  }
  nm <- .reference_param_names(mu, Sigma, param_names = param_names)
  names(mu) <- nm
  dimnames(Sigma) <- list(nm, nm)
  L <- tryCatch(chol(Sigma), error = function(e) chol(diag(d) * 1e-6 + Sigma))
  list(
    mean = mu,
    cov = Sigma,
    chol = L,
    weight = as.numeric(weight)
  )
}

.reference_mixture_moments <- function(components, weights = NULL) {
  if (!length(components)) stop("At least one reference-prior component is required.")
  weights <- as.numeric(weights %||% vapply(components, `[[`, numeric(1), "weight"))
  weights <- pmax(weights, 0)
  sw <- sum(weights)
  if (!is.finite(sw) || sw <= 0) stop("Reference prior weights must sum to a positive value.")
  weights <- weights / sw

  means <- lapply(components, `[[`, "mean")
  covs <- lapply(components, `[[`, "cov")
  mu <- Reduce(`+`, Map(function(w, m) w * m, weights, means))
  Sigma <- Reduce(
    `+`,
    Map(
      function(w, m, S) {
        dm <- as.numeric(m - mu)
        w * (S + tcrossprod(dm))
      },
      weights,
      means,
      covs
    )
  )
  Sigma <- regularize_cov(Sigma, min_eig = 1e-8, cond_cap = 1e8)
  list(mu = mu, Sigma = Sigma, weights = weights)
}

.build_reference_prior <- function(components, family, label = NULL) {
  if (!length(components)) stop("Reference prior must contain at least one component.")
  raw_weights <- pmax(vapply(components, `[[`, numeric(1), "weight"), 0)
  sw <- sum(raw_weights)
  if (!is.finite(sw) || sw <= 0) stop("Reference prior component weights must be positive.")
  raw_weights <- raw_weights / sw
  for (k in seq_along(components)) {
    components[[k]]$weight <- raw_weights[k]
    components[[k]]$log_weight <- log(raw_weights[k])
  }

  moments <- .reference_mixture_moments(components, raw_weights)
  nm <- names(components[[1L]]$mean)
  names(moments$mu) <- nm
  dimnames(moments$Sigma) <- list(nm, nm)

  structure(
    list(
      family = family,
      label = label %||% family,
      d = length(moments$mu),
      param_names = nm,
      mean = moments$mu,
      cov = moments$Sigma,
      chol = chol(moments$Sigma),
      n_components = length(components),
      components = components
    ),
    class = "reference_prior"
  )
}

make_reference_prior_gaussian <- function(mu, Sigma, scale = 1, param_names = NULL, label = "gaussian") {
  nm <- .reference_param_names(mu, Sigma, param_names = param_names)
  mu <- stats::setNames(as.numeric(mu), nm)
  Sigma <- as.matrix(Sigma)
  dimnames(Sigma) <- list(nm, nm)
  Sigma <- Sigma * as.numeric(scale)
  component <- .normalize_reference_component(mu, Sigma, weight = 1, param_names = nm)
  .build_reference_prior(list(component), family = "gaussian", label = label)
}

make_reference_prior_gaussian_mixture <- function(component_means,
                                                  component_covs,
                                                  weights = NULL,
                                                  param_names = NULL,
                                                  label = "gaussian_mixture") {
  if (length(component_means) != length(component_covs)) {
    stop("component_means and component_covs must have the same length.")
  }
  if (!length(component_means)) {
    stop("At least one component is required.")
  }
  weights <- as.numeric(weights %||% rep(1 / length(component_means), length(component_means)))
  components <- lapply(
    seq_along(component_means),
    function(k) {
      .normalize_reference_component(
        mu = component_means[[k]],
        Sigma = component_covs[[k]],
        weight = weights[k],
        param_names = param_names
      )
    }
  )
  .build_reference_prior(components, family = "gaussian_mixture", label = label)
}

combine_reference_priors <- function(priors, weights = NULL, label = "gaussian_mixture") {
  if (!length(priors)) stop("At least one reference prior is required.")
  weights <- as.numeric(weights %||% rep(1 / length(priors), length(priors)))
  if (length(weights) != length(priors)) {
    stop("weights must have the same length as priors.")
  }
  flat_components <- list()
  ptr <- 1L
  for (i in seq_along(priors)) {
    prior_i <- normalize_reference_prior(reference_prior = priors[[i]])
    for (comp in prior_i$components) {
      flat_components[[ptr]] <- .normalize_reference_component(
        mu = comp$mean,
        Sigma = comp$cov,
        weight = weights[i] * comp$weight,
        param_names = prior_i$param_names
      )
      ptr <- ptr + 1L
    }
  }
  .build_reference_prior(flat_components, family = "gaussian_mixture", label = label)
}

inflate_reference_prior <- function(reference_prior, scale = 1.5, label = NULL) {
  prior <- normalize_reference_prior(reference_prior = reference_prior)
  components <- lapply(
    prior$components,
    function(comp) {
      .normalize_reference_component(
        mu = comp$mean,
        Sigma = as.numeric(scale) * comp$cov,
        weight = comp$weight,
        param_names = prior$param_names
      )
    }
  )
  .build_reference_prior(
    components,
    family = if (length(components) == 1L && identical(prior$family, "gaussian")) "gaussian" else "gaussian_mixture",
    label = label %||% prior$label
  )
}

make_broad_reference_prior <- function(mu,
                                       Sigma,
                                       scale = 4,
                                       defensive = FALSE,
                                       defensive_scale = 16,
                                       defensive_weight = 0.10,
                                       param_names = NULL,
                                       label = NULL) {
  broad_core <- make_reference_prior_gaussian(
    mu = mu,
    Sigma = Sigma,
    scale = scale,
    param_names = param_names,
    label = label %||% "broad_gaussian"
  )
  if (!isTRUE(defensive)) {
    return(broad_core)
  }
  defensive_component <- make_reference_prior_gaussian(
    mu = mu,
    Sigma = Sigma,
    scale = defensive_scale,
    param_names = broad_core$param_names,
    label = "defensive_component"
  )
  combine_reference_priors(
    priors = list(broad_core, defensive_component),
    weights = c(1 - defensive_weight, defensive_weight),
    label = label %||% "broad_defensive_mixture"
  )
}

normalize_reference_prior <- function(reference_prior = NULL, mu = NULL, Sigma = NULL) {
  if (!is.null(reference_prior)) {
    if (!inherits(reference_prior, "reference_prior")) {
      stop("reference_prior must inherit from 'reference_prior'.")
    }
    return(reference_prior)
  }
  if (is.null(mu) || is.null(Sigma)) {
    stop("Provide either reference_prior or both mu and Sigma.")
  }
  make_reference_prior_gaussian(mu = mu, Sigma = Sigma)
}

reference_prior_is_gaussian <- function(reference_prior) {
  prior <- normalize_reference_prior(reference_prior = reference_prior)
  identical(prior$family, "gaussian") && identical(prior$n_components, 1L)
}

reference_prior_geometry <- function(reference_prior) {
  prior <- normalize_reference_prior(reference_prior = reference_prior)
  list(
    mean = prior$mean,
    cov = prior$cov,
    chol = prior$chol,
    param_names = prior$param_names
  )
}

reference_prior_logpdf <- function(reference_prior, Theta) {
  prior <- normalize_reference_prior(reference_prior = reference_prior)
  Theta <- .reference_matrix(Theta, d = prior$d, param_names = prior$param_names)
  if (prior$n_components == 1L) {
    comp <- prior$components[[1L]]
    return(as.numeric(dmvnorm_chol_log(Theta, comp$mean, comp$chol)))
  }
  logs <- matrix(-Inf, nrow = nrow(Theta), ncol = prior$n_components)
  for (k in seq_len(prior$n_components)) {
    comp <- prior$components[[k]]
    logs[, k] <- comp$log_weight + dmvnorm_chol_log(Theta, comp$mean, comp$chol)
  }
  as.numeric(.rowLogSumExp(logs))
}

reference_prior_sample <- function(reference_prior, n, seed = NULL) {
  prior <- normalize_reference_prior(reference_prior = reference_prior)
  n <- as.integer(n)
  if (n <= 0L) {
    return(matrix(0, nrow = 0L, ncol = prior$d, dimnames = list(NULL, prior$param_names)))
  }
  if (!is.null(seed)) set.seed(as.integer(seed))

  out <- matrix(NA_real_, nrow = n, ncol = prior$d)
  comp_idx <- sample.int(prior$n_components, size = n, replace = TRUE, prob = vapply(prior$components, `[[`, numeric(1), "weight"))
  for (k in seq_len(prior$n_components)) {
    take <- which(comp_idx == k)
    if (!length(take)) next
    comp <- prior$components[[k]]
    out[take, ] <- mvtnorm::rmvnorm(length(take), mean = comp$mean, sigma = comp$cov)
  }
  colnames(out) <- prior$param_names
  out
}

reference_prior_summary <- function(reference_prior) {
  prior <- normalize_reference_prior(reference_prior = reference_prior)
  list(
    family = prior$family,
    label = prior$label,
    d = prior$d,
    n_components = prior$n_components,
    param_names = prior$param_names,
    mean = prior$mean,
    cov = prior$cov,
    weights = vapply(prior$components, `[[`, numeric(1), "weight")
  )
}
