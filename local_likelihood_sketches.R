#!/usr/bin/env Rscript
# ============================================================================
# Local likelihood sketches
# - Cheap theta-independent approximations to L_i(alpha)
# - Phase 1 support-engine building block
# ============================================================================

if (!exists("%||%", mode = "function") ||
    !exists("logsumexp", mode = "function") ||
    !exists(".rowLogSumExp", mode = "function") ||
    !exists("dmvnorm_chol_log", mode = "function") ||
    !exists("regularize_cov", mode = "function") ||
    !exists("safe_logdet", mode = "function")) {
  source("smc_core.R")
}

if (!exists("normalize_population_model", mode = "function") ||
    !exists("population_model_reference_components_from_theta", mode = "function")) {
  source("population_models.R")
}

suppressPackageStartupMessages({
  library(parallel)
})

.sketch_as_alpha_matrix <- function(alpha, alpha_names = NULL) {
  alpha <- as.matrix(alpha)
  if (!is.null(alpha_names)) {
    alpha_names <- as.character(alpha_names)
    if (ncol(alpha) != length(alpha_names)) {
      if (length(alpha) == length(alpha_names)) {
        alpha <- matrix(as.numeric(alpha), nrow = 1L)
      } else {
        stop("alpha has incompatible dimension.")
      }
    }
    if (!is.null(colnames(alpha)) && setequal(colnames(alpha), alpha_names)) {
      alpha <- alpha[, alpha_names, drop = FALSE]
    } else {
      colnames(alpha) <- alpha_names
    }
  }
  alpha
}

.sketch_loglik <- function(alpha, data_i, loglik_fn, alpha_names) {
  alpha <- .sketch_as_alpha_matrix(alpha, alpha_names)
  out <- as.numeric(loglik_fn(alpha, data_i))
  if (length(out) != nrow(alpha)) {
    stop("loglik_fn must return one value per alpha row.")
  }
  out
}

.sketch_reference_start <- function(population_model,
                                    theta_reference,
                                    alpha_names) {
  if (is.null(population_model) || is.null(theta_reference)) return(NULL)
  model <- normalize_population_model(population_model)
  ref <- population_model_reference_components_from_theta(model, theta_reference)
  list(
    mean = as.numeric(ref$component_means[[1L]][alpha_names]),
    cov = ref$component_covs[[1L]][alpha_names, alpha_names, drop = FALSE]
  )
}

.sketch_start_matrix <- function(alpha_names,
                                 starts = NULL,
                                 initial_alpha = NULL,
                                 start_cov = NULL,
                                 population_model = NULL,
                                 theta_reference = NULL,
                                 n_starts = 8L,
                                 start_scale = 1,
                                 seed = NULL) {
  alpha_names <- as.character(alpha_names)
  d <- length(alpha_names)
  ref <- .sketch_reference_start(population_model, theta_reference, alpha_names)

  start_rows <- list()
  if (!is.null(starts)) {
    starts <- .sketch_as_alpha_matrix(starts, alpha_names)
    for (i in seq_len(nrow(starts))) start_rows[[length(start_rows) + 1L]] <- starts[i, ]
  }
  if (!is.null(initial_alpha)) {
    initial_alpha <- .sketch_as_alpha_matrix(initial_alpha, alpha_names)
    for (i in seq_len(nrow(initial_alpha))) {
      start_rows[[length(start_rows) + 1L]] <- initial_alpha[i, ]
    }
  }
  if (!is.null(ref)) {
    start_rows[[length(start_rows) + 1L]] <- stats::setNames(ref$mean, alpha_names)
    if (is.null(start_cov)) start_cov <- ref$cov
  }

  if (!length(start_rows)) {
    stop("Provide starts, initial_alpha, or population_model plus theta_reference.")
  }

  n_starts <- as.integer(max(1L, n_starts))
  base <- as.numeric(start_rows[[1L]])
  cov <- as.matrix(start_cov %||% diag(1, d))
  if (nrow(cov) != d || ncol(cov) != d) {
    stop("start_cov has incompatible dimension.")
  }
  cov <- regularize_cov(cov * as.numeric(start_scale), min_eig = 1e-8, cond_cap = 1e8)

  if (!is.null(seed)) set.seed(as.integer(seed))
  while (length(start_rows) < n_starts) {
    z <- as.numeric(stats::rnorm(d) %*% chol(cov))
    start_rows[[length(start_rows) + 1L]] <- stats::setNames(base + z, alpha_names)
  }

  out <- do.call(rbind, start_rows)
  out <- out[seq_len(min(nrow(out), n_starts)), , drop = FALSE]
  colnames(out) <- alpha_names
  out
}

.sketch_gaussian_reference <- function(population_model,
                                       theta_reference,
                                       alpha_names,
                                       reference_mean = NULL,
                                       reference_cov = NULL,
                                       reference_scale = 25) {
  alpha_names <- as.character(alpha_names)
  d <- length(alpha_names)
  if (is.null(reference_mean) || is.null(reference_cov)) {
    ref <- .sketch_reference_start(population_model, theta_reference, alpha_names)
    if (is.null(ref)) {
      stop("Reference-prior site sketches require reference_mean/reference_cov or population_model plus theta_reference.")
    }
    reference_mean <- reference_mean %||% ref$mean
    reference_cov <- reference_cov %||% ref$cov
  }
  reference_mean <- as.numeric(reference_mean)
  if (length(reference_mean) != d) {
    stop("reference_mean has incompatible dimension.")
  }
  names(reference_mean) <- alpha_names
  reference_cov <- as.matrix(reference_cov)
  if (nrow(reference_cov) != d || ncol(reference_cov) != d) {
    stop("reference_cov has incompatible dimension.")
  }
  reference_cov <- regularize_cov(reference_cov * as.numeric(reference_scale), min_eig = 1e-8, cond_cap = 1e10)
  dimnames(reference_cov) <- list(alpha_names, alpha_names)
  reference_precision <- tryCatch(solve(reference_cov), error = function(e) NULL)
  if (is.null(reference_precision) || any(!is.finite(reference_precision))) {
    stop("reference_cov could not be inverted.")
  }
  reference_precision <- (reference_precision + t(reference_precision)) / 2
  list(
    mean = reference_mean,
    cov = reference_cov,
    precision = reference_precision,
    logdet_cov = safe_logdet(reference_cov),
    scale = as.numeric(reference_scale)
  )
}

.sketch_gaussian_logpdf <- function(alpha, mean, cov) {
  alpha <- as.matrix(alpha)
  dmvnorm_chol_log(alpha, mean, chol(cov))
}

.sketch_objective_factory <- function(data_i,
                                      loglik_fn,
                                      alpha_names,
                                      reference = NULL,
                                      invalid_penalty = 1e100) {
  force(data_i)
  force(loglik_fn)
  force(alpha_names)
  force(reference)
  function(par) {
    alpha <- matrix(as.numeric(par), nrow = 1L)
    colnames(alpha) <- alpha_names
    ll <- .sketch_loglik(alpha, data_i, loglik_fn, alpha_names)[1L]
    if (!is.finite(ll)) return(invalid_penalty)
    ref_log <- if (is.null(reference)) {
      0
    } else {
      .sketch_gaussian_logpdf(alpha, reference$mean, reference$cov)[1L]
    }
    target <- ll + ref_log
    if (!is.finite(target)) return(invalid_penalty)
    -target
  }
}

.sketch_optimize_start <- function(start,
                                   objective,
                                   maxit = 1000L,
                                   nm_maxit = 400L) {
  start <- as.numeric(start)
  nm <- tryCatch(
    stats::optim(
      par = start,
      fn = objective,
      method = "Nelder-Mead",
      control = list(maxit = as.integer(nm_maxit), warn.1d.NelderMead = FALSE)
    ),
    error = function(e) list(par = start, value = Inf, convergence = 999L, message = conditionMessage(e))
  )

  bfgs <- tryCatch(
    stats::optim(
      par = as.numeric(nm$par),
      fn = objective,
      method = "BFGS",
      hessian = TRUE,
      control = list(maxit = as.integer(maxit))
    ),
    error = function(e) list(par = as.numeric(nm$par), value = Inf, convergence = 999L, message = conditionMessage(e), hessian = NULL)
  )

  list(
    par = as.numeric(bfgs$par),
    objective = as.numeric(bfgs$value),
    log_likelihood = -as.numeric(bfgs$value),
    convergence = as.integer(bfgs$convergence %||% 999L),
    message = as.character(bfgs$message %||% ""),
    hessian = bfgs$hessian,
    nelder_mead_objective = as.numeric(nm$value),
    nelder_mead_convergence = as.integer(nm$convergence %||% 999L),
    nelder_mead_message = as.character(nm$message %||% "")
  )
}

.sketch_gaussian_site_from_density <- function(log_weight, mean, cov, alpha_names) {
  d <- length(alpha_names)
  cov <- regularize_cov(cov, min_eig = 1e-8, cond_cap = 1e10)
  P <- solve(cov)
  P <- (P + t(P)) / 2
  mean <- as.numeric(mean)
  h <- as.numeric(P %*% mean)
  c0 <- as.numeric(log_weight - 0.5 * (d * log(2 * pi) + safe_logdet(cov) + sum(mean * h)))
  list(constant = c0, linear = h, precision = P)
}

.sketch_repair_site_natural <- function(precision,
                                        linear,
                                        min_positive = 1e-8,
                                        negative_tol = 1e-7) {
  precision <- (as.matrix(precision) + t(as.matrix(precision))) / 2
  ev <- eigen(precision, symmetric = TRUE)
  values <- ev$values
  linear_rot <- as.numeric(t(ev$vectors) %*% as.numeric(linear))
  keep <- is.finite(values) & values > min_positive
  dropped_negative <- values[is.finite(values) & values < -abs(negative_tol)]
  values_repaired <- ifelse(keep, values, 0)
  linear_rot[!keep] <- 0
  P <- ev$vectors %*% diag(values_repaired, length(values_repaired)) %*% t(ev$vectors)
  h <- as.numeric(ev$vectors %*% linear_rot)
  P <- (P + t(P)) / 2
  list(
    precision = P,
    linear = h,
    rank = sum(keep),
    min_raw_eigen = min(values, na.rm = TRUE),
    max_raw_eigen = max(values, na.rm = TRUE),
    dropped_negative_count = length(dropped_negative),
    dropped_negative_min = if (length(dropped_negative)) min(dropped_negative) else 0
  )
}

.sketch_site_from_tilted_laplace <- function(mode,
                                             log_likelihood_at_mode,
                                             tilted_precision,
                                             reference,
                                             alpha_names,
                                             site_cov_cap = 1e6) {
  alpha_names <- as.character(alpha_names)
  d <- length(alpha_names)
  tilted_precision <- (as.matrix(tilted_precision) + t(as.matrix(tilted_precision))) / 2
  reference_precision <- reference$precision
  raw_precision <- tilted_precision - reference_precision
  raw_linear <- as.numeric(tilted_precision %*% as.numeric(mode) -
                             reference_precision %*% reference$mean)
  repaired <- .sketch_repair_site_natural(raw_precision, raw_linear)
  P <- repaired$precision
  h <- repaired$linear
  mode <- as.numeric(mode)
  constant <- as.numeric(log_likelihood_at_mode - sum(h * mode) + 0.5 * sum(mode * as.numeric(P %*% mode)))

  ev <- eigen(P, symmetric = TRUE)
  precision_for_cov <- ev$vectors %*%
    diag(pmax(ev$values, 1 / as.numeric(site_cov_cap)), d) %*%
    t(ev$vectors)
  cov <- solve((precision_for_cov + t(precision_for_cov)) / 2)
  cov <- (cov + t(cov)) / 2
  mean <- mode
  dimnames(cov) <- list(alpha_names, alpha_names)
  names(mean) <- alpha_names

  list(
    constant = constant,
    linear = h,
    precision = P,
    display_mean = mean,
    display_cov = cov,
    rank = repaired$rank,
    min_raw_eigen = repaired$min_raw_eigen,
    max_raw_eigen = repaired$max_raw_eigen,
    dropped_negative_count = repaired$dropped_negative_count,
    dropped_negative_min = repaired$dropped_negative_min
  )
}

.sketch_precision_to_cov <- function(precision,
                                     fallback_cov,
                                     min_eig = 1e-6,
                                     cond_cap = 1e8) {
  d <- nrow(fallback_cov)
  fallback_cov <- regularize_cov(fallback_cov, min_eig = min_eig, cond_cap = cond_cap)
  fallback_precision <- tryCatch(
    solve(fallback_cov),
    error = function(e) diag(1 / pmax(diag(fallback_cov), min_eig), d)
  )
  precision <- as.matrix(precision)
  valid_precision <- nrow(precision) == d &&
    ncol(precision) == d &&
    all(is.finite(precision))

  if (valid_precision) {
    precision <- (precision + t(precision)) / 2
    ev <- eigen(precision, symmetric = TRUE)
    valid_precision <- all(is.finite(ev$values))
    if (valid_precision) {
      fallback_eigen_precision <- diag(t(ev$vectors) %*% fallback_precision %*% ev$vectors)
      repaired_values <- ev$values
      bad <- !is.finite(repaired_values) | repaired_values < min_eig
      repaired_values[bad] <- pmax(fallback_eigen_precision[bad], min_eig)
      precision_reg <- ev$vectors %*% diag(repaired_values, d) %*% t(ev$vectors)
      precision_reg <- (precision_reg + t(precision_reg)) / 2
      precision_reg <- regularize_cov(precision_reg, min_eig = min_eig, cond_cap = cond_cap)
      cov <- tryCatch(solve(precision_reg), error = function(e) NULL)
      if (!is.null(cov) && all(is.finite(cov))) {
        cov <- (cov + t(cov)) / 2
        return(list(cov = cov, precision = precision_reg, hessian_valid = !any(bad)))
      }
    }
  }

  precision_reg <- fallback_precision
  if (all(is.finite(precision_reg))) {
    precision_reg <- regularize_cov(precision_reg, min_eig = min_eig, cond_cap = cond_cap)
    cov <- tryCatch(solve(precision_reg), error = function(e) NULL)
    if (!is.null(cov) && all(is.finite(cov))) {
      cov <- (cov + t(cov)) / 2
      return(list(cov = cov, precision = precision_reg, hessian_valid = FALSE))
    }
  }

  list(cov = fallback_cov, precision = fallback_precision, hessian_valid = FALSE)
}

new_local_likelihood_sketch <- function(local_id,
                                        alpha_names,
                                        component_log_weights,
                                        component_means,
                                        component_covs,
                                        component_log_likelihood,
                                        component_site_constant = NULL,
                                        component_site_linear = NULL,
                                        component_site_precision = NULL,
                                        diagnostics = list()) {
  alpha_names <- as.character(alpha_names)
  component_means <- as.matrix(component_means)
  if (ncol(component_means) != length(alpha_names)) {
    stop("component_means dimension does not match alpha_names.")
  }
  colnames(component_means) <- alpha_names
  component_covs <- lapply(component_covs, function(S) {
    S <- as.matrix(S)
    if (nrow(S) != length(alpha_names) || ncol(S) != length(alpha_names)) {
      stop("component covariance has incompatible dimension.")
    }
    S <- regularize_cov(S, min_eig = 1e-8, cond_cap = 1e8)
    dimnames(S) <- list(alpha_names, alpha_names)
    S
  })
  component_log_weights <- as.numeric(component_log_weights)
  component_log_likelihood <- as.numeric(component_log_likelihood)
  k <- nrow(component_means)
  if (length(component_covs) != k ||
      length(component_log_weights) != k ||
      length(component_log_likelihood) != k) {
    stop("component fields must have matching lengths.")
  }
  if (any(!is.finite(component_log_weights)) ||
      any(!is.finite(component_log_likelihood))) {
    stop("component log weights and log likelihoods must be finite.")
  }
  if (is.null(component_site_constant) ||
      is.null(component_site_linear) ||
      is.null(component_site_precision)) {
    site_parts <- lapply(seq_len(k), function(j) {
      .sketch_gaussian_site_from_density(
        log_weight = component_log_weights[j],
        mean = component_means[j, ],
        cov = component_covs[[j]],
        alpha_names = alpha_names
      )
    })
    component_site_constant <- vapply(site_parts, `[[`, numeric(1), "constant")
    component_site_linear <- do.call(rbind, lapply(site_parts, `[[`, "linear"))
    component_site_precision <- lapply(site_parts, `[[`, "precision")
  }
  component_site_constant <- as.numeric(component_site_constant)
  component_site_linear <- as.matrix(component_site_linear)
  if (nrow(component_site_linear) != k || ncol(component_site_linear) != length(alpha_names)) {
    stop("component_site_linear has incompatible dimension.")
  }
  colnames(component_site_linear) <- alpha_names
  if (length(component_site_constant) != k || length(component_site_precision) != k) {
    stop("component site fields must match the number of components.")
  }
  component_site_precision <- lapply(component_site_precision, function(P) {
    P <- as.matrix(P)
    if (nrow(P) != length(alpha_names) || ncol(P) != length(alpha_names)) {
      stop("component_site_precision has incompatible dimension.")
    }
    P <- (P + t(P)) / 2
    dimnames(P) <- list(alpha_names, alpha_names)
    P
  })
  if (any(!is.finite(component_site_constant)) ||
      any(!is.finite(component_site_linear)) ||
      any(vapply(component_site_precision, function(P) any(!is.finite(P)), logical(1)))) {
    stop("component site fields must be finite.")
  }

  structure(
    list(
      local_id = as.integer(local_id %||% NA_integer_),
      alpha_names = alpha_names,
      alpha_dim = length(alpha_names),
      n_components = as.integer(k),
      component_log_weights = component_log_weights,
      component_means = component_means,
      component_covs = component_covs,
      component_log_likelihood = component_log_likelihood,
      component_site_constant = component_site_constant,
      component_site_linear = component_site_linear,
      component_site_precision = component_site_precision,
      diagnostics = diagnostics %||% list()
    ),
    class = "local_likelihood_sketch"
  )
}

validate_local_likelihood_sketch <- function(sketch) {
  if (!inherits(sketch, "local_likelihood_sketch")) {
    stop("sketch must inherit from 'local_likelihood_sketch'.")
  }
  new_local_likelihood_sketch(
    local_id = sketch$local_id,
    alpha_names = sketch$alpha_names,
    component_log_weights = sketch$component_log_weights,
    component_means = sketch$component_means,
    component_covs = sketch$component_covs,
    component_log_likelihood = sketch$component_log_likelihood,
    component_site_constant = sketch$component_site_constant,
    component_site_linear = sketch$component_site_linear,
    component_site_precision = sketch$component_site_precision,
    diagnostics = sketch$diagnostics
  )
}

.sketch_site_log_eval <- function(alpha, constant, linear, precision) {
  alpha <- as.matrix(alpha)
  quad <- rowSums((alpha %*% precision) * alpha)
  as.numeric(constant + alpha %*% as.numeric(linear) - 0.5 * quad)
}

local_likelihood_sketch_log_approx <- function(sketch, alpha) {
  sketch <- validate_local_likelihood_sketch(sketch)
  alpha <- .sketch_as_alpha_matrix(alpha, sketch$alpha_names)
  logs <- matrix(NA_real_, nrow = nrow(alpha), ncol = sketch$n_components)
  for (k in seq_len(sketch$n_components)) {
    logs[, k] <- .sketch_site_log_eval(
      alpha = alpha,
      constant = sketch$component_site_constant[k],
      linear = sketch$component_site_linear[k, ],
      precision = sketch$component_site_precision[[k]]
    )
  }
  as.numeric(.rowLogSumExp(logs))
}

local_likelihood_sketch_tempered_components <- function(sketch, rho = 1) {
  sketch <- validate_local_likelihood_sketch(sketch)
  rho <- as.numeric(rho)
  if (length(rho) != 1L || !is.finite(rho) || rho <= 0) {
    stop("rho must be a finite positive scalar.")
  }
  list(
    rho = rho,
    component_site_constant = rho * sketch$component_site_constant,
    component_site_linear = rho * sketch$component_site_linear,
    component_site_precision = lapply(sketch$component_site_precision, function(P) rho * P),
    component_means = sketch$component_means,
    component_power_tempering = TRUE
  )
}

.sketch_site_gaussian_integral_log <- function(site_constant,
                                               site_linear,
                                               site_precision,
                                               prior_mean,
                                               prior_cov) {
  d <- length(prior_mean)
  prior_precision <- solve(prior_cov)
  prior_precision <- (prior_precision + t(prior_precision)) / 2
  prior_linear <- as.numeric(prior_precision %*% as.numeric(prior_mean))
  prior_constant <- as.numeric(
    -0.5 * (
      d * log(2 * pi) +
        safe_logdet(prior_cov) +
        sum(as.numeric(prior_mean) * prior_linear)
    )
  )
  A <- site_precision + prior_precision
  A <- (A + t(A)) / 2
  L <- tryCatch(chol(A), error = function(e) chol(regularize_cov(A, min_eig = 1e-10, cond_cap = 1e12)))
  b <- as.numeric(site_linear) + prior_linear
  sol <- backsolve(L, b, transpose = TRUE)
  quad <- sum(sol^2)
  as.numeric(site_constant + prior_constant + 0.5 * d * log(2 * pi) - sum(log(diag(L))) + 0.5 * quad)
}

local_likelihood_sketch_log_marginal <- function(sketch,
                                                 population_model,
                                                 theta,
                                                 rho = 1) {
  sketch <- validate_local_likelihood_sketch(sketch)
  model <- normalize_population_model(population_model)
  if (!setequal(sketch$alpha_names, model$alpha_names)) {
    stop("sketch alpha_names do not match the population model alpha_names.")
  }
  theta <- .as_hyper_matrix(theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  ref <- population_model_reference_components_from_theta(model, theta)
  temp <- local_likelihood_sketch_tempered_components(sketch, rho = rho)

  terms <- matrix(-Inf, nrow = nrow(theta), ncol = sketch$n_components)
  for (j in seq_len(nrow(theta))) {
    pop_mean <- as.numeric(ref$component_means[[j]][sketch$alpha_names])
    pop_cov <- ref$component_covs[[j]][sketch$alpha_names, sketch$alpha_names, drop = FALSE]
    for (k in seq_len(sketch$n_components)) {
      terms[j, k] <- .sketch_site_gaussian_integral_log(
        site_constant = temp$component_site_constant[k],
        site_linear = temp$component_site_linear[k, ],
        site_precision = temp$component_site_precision[[k]],
        prior_mean = pop_mean,
        prior_cov = pop_cov
      )
    }
  }

  if (sketch$n_components == 1L) as.numeric(terms[, 1L]) else as.numeric(.rowLogSumExp(terms))
}

local_likelihood_sketch_set_log_marginal_matrix <- function(sketch_set,
                                                            population_model,
                                                            theta,
                                                            rho = 1,
                                                            n_jobs = 1L) {
  sketch_set <- validate_local_likelihood_sketch_set(sketch_set)
  model <- normalize_population_model(population_model)
  theta <- .as_hyper_matrix(theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  ids <- seq_along(sketch_set$sketches)
  parts <- parallel::mclapply(
    ids,
    function(i) {
      local_likelihood_sketch_log_marginal(
        sketch = sketch_set$sketches[[i]],
        population_model = model,
        theta = theta,
        rho = rho
      )
    },
    mc.cores = as.integer(max(1L, n_jobs))
  )
  out <- do.call(cbind, parts)
  colnames(out) <- names(sketch_set$sketches) %||% paste0("local_", ids)
  rownames(out) <- rownames(theta)
  out
}

local_likelihood_sketch_set_loglik <- function(sketch_set,
                                               population_model,
                                               theta,
                                               rho = 1,
                                               n_jobs = 1L) {
  mat <- local_likelihood_sketch_set_log_marginal_matrix(
    sketch_set = sketch_set,
    population_model = population_model,
    theta = theta,
    rho = rho,
    n_jobs = n_jobs
  )
  as.numeric(rowSums(mat))
}

local_likelihood_sketch_set_logposterior <- function(sketch_set,
                                                     population_model,
                                                     theta,
                                                     rho = 1,
                                                     n_jobs = 1L) {
  model <- normalize_population_model(population_model)
  theta <- .as_hyper_matrix(theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  population_model_log_hyperprior(model, theta) +
    local_likelihood_sketch_set_loglik(
      sketch_set = sketch_set,
      population_model = model,
      theta = theta,
      rho = rho,
      n_jobs = n_jobs
    )
}

build_local_likelihood_sketch_factor_set <- function(sketch_set,
                                                     population_model,
                                                     rho = 1) {
  sketch_set <- validate_local_likelihood_sketch_set(sketch_set)
  model <- normalize_population_model(population_model)
  rho <- as.numeric(rho)
  if (length(rho) != 1L || !is.finite(rho) || rho <= 0) {
    stop("rho must be a finite positive scalar.")
  }
  if (!setequal(sketch_set$alpha_names, model$alpha_names)) {
    stop("sketch_set alpha_names do not match the population model alpha_names.")
  }
  structure(
    list(
      sketch_set = sketch_set,
      population_model = model,
      rho = rho,
      n_locals = sketch_set$n_locals,
      log_constant = 0
    ),
    class = c("local_likelihood_sketch_factor_set", "population_factor_set")
  )
}

fit_local_likelihood_sketch <- function(data_i,
                                        loglik_fn,
                                        alpha_names,
                                        local_id = NULL,
                                        starts = NULL,
                                        initial_alpha = NULL,
                                        start_cov = NULL,
                                        reference_mean = NULL,
                                        reference_cov = NULL,
                                        reference_scale = 25,
                                        population_model = NULL,
                                        theta_reference = NULL,
                                        n_starts = 8L,
                                        start_scale = 1,
                                        maxit = 1000L,
                                        nm_maxit = 400L,
                                        seed = NULL) {
  alpha_names <- as.character(alpha_names)
  d <- length(alpha_names)
  reference <- .sketch_gaussian_reference(
    population_model = population_model,
    theta_reference = theta_reference,
    alpha_names = alpha_names,
    reference_mean = reference_mean,
    reference_cov = reference_cov,
    reference_scale = reference_scale
  )
  if (is.null(start_cov)) start_cov <- reference$cov
  starts <- .sketch_start_matrix(
    alpha_names = alpha_names,
    starts = starts,
    initial_alpha = initial_alpha,
    start_cov = start_cov,
    population_model = population_model,
    theta_reference = theta_reference,
    n_starts = n_starts,
    start_scale = start_scale,
    seed = seed
  )
  objective <- .sketch_objective_factory(data_i, loglik_fn, alpha_names, reference = reference)
  fits <- lapply(seq_len(nrow(starts)), function(i) {
    .sketch_optimize_start(
      start = starts[i, ],
      objective = objective,
      maxit = maxit,
      nm_maxit = nm_maxit
    )
  })
  log_target <- vapply(fits, `[[`, numeric(1), "log_likelihood")
  if (!any(is.finite(log_target))) {
    stop("No finite local reference-site optimum was found.")
  }
  best <- which.max(log_target)
  mode <- fits[[best]]$par
  names(mode) <- alpha_names
  mode_mat <- matrix(mode, nrow = 1L, dimnames = list(NULL, alpha_names))
  mode_log_likelihood <- .sketch_loglik(mode_mat, data_i, loglik_fn, alpha_names)[1L]
  mode_log_reference <- .sketch_gaussian_logpdf(mode_mat, reference$mean, reference$cov)[1L]
  mode_log_target <- log_target[best]

  fallback_cov <- reference$cov
  fallback_cov <- regularize_cov(as.matrix(fallback_cov), min_eig = 1e-6, cond_cap = 1e8)

  curv <- .sketch_precision_to_cov(
    precision = fits[[best]]$hessian,
    fallback_cov = fallback_cov,
    min_eig = 1e-6,
    cond_cap = 1e8
  )
  site <- .sketch_site_from_tilted_laplace(
    mode = mode,
    log_likelihood_at_mode = mode_log_likelihood,
    tilted_precision = curv$precision,
    reference = reference,
    alpha_names = alpha_names
  )
  cov <- site$display_cov
  log_c <- .sketch_site_log_eval(
    alpha = matrix(site$display_mean, nrow = 1L, dimnames = list(NULL, alpha_names)),
    constant = site$constant,
    linear = site$linear,
    precision = site$precision
  )[1L]

  opt_df <- data.frame(
    start = seq_len(nrow(starts)),
    log_target = log_target,
    log_likelihood = vapply(fits, function(fit) {
      alpha <- matrix(fit$par, nrow = 1L, dimnames = list(NULL, alpha_names))
      .sketch_loglik(alpha, data_i, loglik_fn, alpha_names)[1L]
    }, numeric(1)),
    convergence = vapply(fits, `[[`, integer(1), "convergence"),
    nelder_mead_log_target = -vapply(fits, `[[`, numeric(1), "nelder_mead_objective"),
    nelder_mead_convergence = vapply(fits, `[[`, integer(1), "nelder_mead_convergence"),
    check.names = FALSE
  )

  hessian_values <- tryCatch(
    eigen((fits[[best]]$hessian + t(fits[[best]]$hessian)) / 2, symmetric = TRUE)$values,
    error = function(e) rep(NA_real_, d)
  )
  new_local_likelihood_sketch(
    local_id = local_id,
    alpha_names = alpha_names,
    component_log_weights = log_c,
    component_means = matrix(site$display_mean, nrow = 1L, dimnames = list(NULL, alpha_names)),
    component_covs = list(cov),
    component_log_likelihood = mode_log_likelihood,
    component_site_constant = site$constant,
    component_site_linear = matrix(site$linear, nrow = 1L, dimnames = list(NULL, alpha_names)),
    component_site_precision = list(site$precision),
    diagnostics = list(
      method = "reference_prior_quadratic_site",
      optimizer = "Nelder-Mead+BFGS",
      starts = starts,
      optimization = opt_df,
      best_start = as.integer(best),
      mode_log_likelihood = as.numeric(mode_log_likelihood),
      mode_log_reference = as.numeric(mode_log_reference),
      mode_log_target = as.numeric(mode_log_target),
      hessian_valid = isTRUE(curv$hessian_valid),
      hessian_eigen_min = min(hessian_values, na.rm = TRUE),
      hessian_eigen_max = max(hessian_values, na.rm = TRUE),
      covariance_condition = kappa(cov),
      regularized = !isTRUE(curv$hessian_valid),
      reference_mean = reference$mean,
      reference_cov = reference$cov,
      reference_scale = reference$scale,
      site_rank = as.integer(site$rank),
      site_min_raw_eigen = as.numeric(site$min_raw_eigen),
      site_max_raw_eigen = as.numeric(site$max_raw_eigen),
      site_dropped_negative_count = as.integer(site$dropped_negative_count),
      site_dropped_negative_min = as.numeric(site$dropped_negative_min)
    )
  )
}

fit_local_likelihood_sketches <- function(data_list,
                                          loglik_fn,
                                          alpha_names,
                                          starts = NULL,
                                          initial_alpha = NULL,
                                          start_cov = NULL,
                                          reference_mean = NULL,
                                          reference_cov = NULL,
                                          reference_scale = 25,
                                          population_model = NULL,
                                          theta_reference = NULL,
                                          n_starts = 8L,
                                          start_scale = 1,
                                          n_jobs = 1L,
                                          seed = 123L) {
  ids <- seq_along(data_list)
  sketches <- parallel::mclapply(
    ids,
    function(i) {
      fit_local_likelihood_sketch(
        data_i = data_list[[i]],
        loglik_fn = loglik_fn,
        alpha_names = alpha_names,
        local_id = i,
        starts = if (is.list(starts)) starts[[i]] else starts,
        initial_alpha = if (is.list(initial_alpha)) initial_alpha[[i]] else initial_alpha,
        start_cov = start_cov,
        reference_mean = reference_mean,
        reference_cov = reference_cov,
        reference_scale = reference_scale,
        population_model = population_model,
        theta_reference = theta_reference,
        n_starts = n_starts,
        start_scale = start_scale,
        seed = as.integer(seed + i)
      )
    },
    mc.cores = as.integer(max(1L, n_jobs))
  )
  names(sketches) <- names(data_list) %||% paste0("local_", ids)
  structure(
    list(
      sketches = sketches,
      alpha_names = as.character(alpha_names),
      n_locals = length(sketches),
      settings = list(
        n_starts = as.integer(n_starts),
        start_scale = as.numeric(start_scale),
        reference_scale = as.numeric(reference_scale),
        seed = as.integer(seed)
      )
    ),
    class = "local_likelihood_sketch_set"
  )
}

validate_local_likelihood_sketch_set <- function(sketch_set) {
  if (!inherits(sketch_set, "local_likelihood_sketch_set")) {
    stop("sketch_set must inherit from 'local_likelihood_sketch_set'.")
  }
  sketch_set$sketches <- lapply(sketch_set$sketches, validate_local_likelihood_sketch)
  sketch_set$n_locals <- length(sketch_set$sketches)
  sketch_set
}

local_likelihood_sketch_profile_check <- function(sketch,
                                                  data_i,
                                                  loglik_fn,
                                                  radii = c(-2, -1, 0, 1, 2),
                                                  max_directions = NULL,
                                                  invalid_loglik_threshold = NULL) {
  sketch <- validate_local_likelihood_sketch(sketch)
  if (sketch$n_components != 1L) {
    stop("Profile checks currently expect a single-component sketch.")
  }
  radii <- as.numeric(radii)
  center <- as.numeric(sketch$component_means[1L, ])
  cov <- sketch$component_covs[[1L]]
  ev <- eigen(cov, symmetric = TRUE)
  d <- sketch$alpha_dim
  dirs <- seq_len(d)
  if (!is.null(max_directions)) dirs <- head(dirs, as.integer(max_directions))

  rows <- vector("list", length(dirs) * length(radii))
  ptr <- 0L
  for (j in dirs) {
    step <- sqrt(max(ev$values[j], 1e-12)) * ev$vectors[, j]
    for (r in radii) {
      alpha <- center + r * step
      alpha_mat <- matrix(alpha, nrow = 1L)
      colnames(alpha_mat) <- sketch$alpha_names
      true <- .sketch_loglik(alpha_mat, data_i, loglik_fn, sketch$alpha_names)[1L]
      approx <- local_likelihood_sketch_log_approx(sketch, alpha_mat)[1L]
      ptr <- ptr + 1L
      rows[[ptr]] <- data.frame(
        local_id = sketch$local_id,
        direction = j,
        radius = r,
        log_likelihood = true,
        sketch_log_likelihood = approx,
        error = approx - true,
        check.names = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  finite <- is.finite(out$log_likelihood) & is.finite(out$sketch_log_likelihood)
  if (!is.null(invalid_loglik_threshold)) {
    finite <- finite & out$log_likelihood > as.numeric(invalid_loglik_threshold)
  }
  out$finite <- finite
  out
}

local_likelihood_sketch_set_profile_check <- function(sketch_set,
                                                      data_list,
                                                      loglik_fn,
                                                      radii = c(-2, -1, 0, 1, 2),
                                                      max_directions = NULL,
                                                      invalid_loglik_threshold = NULL,
                                                      n_jobs = 1L) {
  sketch_set <- validate_local_likelihood_sketch_set(sketch_set)
  if (length(data_list) != length(sketch_set$sketches)) {
    stop("data_list and sketch_set lengths do not match.")
  }
  ids <- seq_along(data_list)
  parts <- parallel::mclapply(
    ids,
    function(i) {
      local_likelihood_sketch_profile_check(
        sketch = sketch_set$sketches[[i]],
        data_i = data_list[[i]],
        loglik_fn = loglik_fn,
        radii = radii,
        max_directions = max_directions,
        invalid_loglik_threshold = invalid_loglik_threshold
      )
    },
    mc.cores = as.integer(max(1L, n_jobs))
  )
  do.call(rbind, parts)
}

local_likelihood_sketch_summary <- function(sketch_set) {
  sketch_set <- validate_local_likelihood_sketch_set(sketch_set)
  do.call(rbind, lapply(sketch_set$sketches, function(sketch) {
    data.frame(
      local_id = sketch$local_id,
      n_components = sketch$n_components,
      mode_log_likelihood = sketch$diagnostics$mode_log_likelihood,
      hessian_valid = isTRUE(sketch$diagnostics$hessian_valid),
      hessian_eigen_min = as.numeric(sketch$diagnostics$hessian_eigen_min),
      hessian_eigen_max = as.numeric(sketch$diagnostics$hessian_eigen_max),
      covariance_condition = as.numeric(sketch$diagnostics$covariance_condition),
      regularized = isTRUE(sketch$diagnostics$regularized),
      check.names = FALSE
    )
  }))
}
