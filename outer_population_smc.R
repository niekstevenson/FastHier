#!/usr/bin/env Rscript
# ============================================================================
# Outer population SMC
# - Reusable local marginal-likelihood factors built from reference local fits
# - Tempered outer SMC for hyperposterior and hyper-evidence estimation
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
if (!exists("normalize_population_model", mode = "function")) {
  source("population_models.R")
}

suppressPackageStartupMessages({
  library(parallel)
})

validate_reference_local_object <- function(local_object) {
  if (!is.list(local_object)) stop("local_object must be a list.")
  required <- c("local_id", "particles", "weights", "reference_prior")
  missing <- setdiff(required, names(local_object))
  if (length(missing)) {
    stop("local_object is missing: ", paste(missing, collapse = ", "))
  }
  particles <- as.matrix(local_object$particles)
  weights <- pmax(as.numeric(local_object$weights), 0)
  sw <- sum(weights)
  if (!is.finite(sw) || sw <= 0) {
    stop("local_object weights are invalid.")
  }
  local_object$particles <- particles
  local_object$weights <- weights / sw
  local_object$reference_prior <- normalize_reference_prior(reference_prior = local_object$reference_prior)
  local_object
}

build_population_local_factor <- function(local_object, population_model) {
  local_object <- validate_reference_local_object(local_object)
  population_model <- normalize_population_model(population_model)

  alpha <- as.matrix(local_object$particles)
  if (ncol(alpha) != population_model$alpha_dim) {
    stop("Local particle dimension does not match the population model.")
  }
  if (!is.null(colnames(alpha)) && setequal(colnames(alpha), population_model$alpha_names)) {
    alpha <- alpha[, population_model$alpha_names, drop = FALSE]
  } else {
    colnames(alpha) <- population_model$alpha_names
  }

  logw <- log(local_object$weights)
  logq <- as.numeric(local_object$log_reference_density %||% reference_prior_logpdf(local_object$reference_prior, alpha))
  log_base <- logw - logq

  structure(
    list(
      local_id = as.integer(local_object$local_id),
      population_model = population_model,
      particles = alpha,
      log_weights = logw,
      log_reference_density = logq,
      log_base = log_base,
      log_constant = as.numeric(local_object$log_evidence %||% 0),
      reference_prior = local_object$reference_prior
    ),
    class = "population_local_factor"
  )
}

population_local_factor_log_marginal_many <- function(factor, theta, include_constant = TRUE) {
  stopifnot(inherits(factor, "population_local_factor"))
  logp <- population_model_log_alpha_given_theta_many(
    factor$population_model,
    alpha = factor$particles,
    theta = theta
  )
  out <- .rowLogSumExp(sweep(logp, 2L, factor$log_base, "+"))
  if (isTRUE(include_constant)) out <- out + factor$log_constant
  as.numeric(out)
}

population_local_factor_log_marginal <- function(factor, theta, include_constant = TRUE) {
  population_local_factor_log_marginal_many(factor, theta = theta, include_constant = include_constant)[1L]
}

population_local_factor_ess <- function(factor, theta) {
  stopifnot(inherits(factor, "population_local_factor"))
  logp <- population_model_log_alpha_given_theta(
    factor$population_model,
    alpha = factor$particles,
    theta = theta
  )
  lw <- factor$log_base + logp
  lse <- logsumexp(lw)
  if (!is.finite(lse)) return(0)
  w <- exp(lw - lse)
  as.numeric(1 / sum(w * w))
}

population_local_factor_reweighted_particles <- function(factor, theta) {
  stopifnot(inherits(factor, "population_local_factor"))
  logp <- population_model_log_alpha_given_theta(
    factor$population_model,
    alpha = factor$particles,
    theta = theta
  )
  lw <- factor$log_base + logp
  lse <- logsumexp(lw)
  if (!is.finite(lse)) {
    weights <- rep(1 / nrow(factor$particles), nrow(factor$particles))
  } else {
    weights <- exp(lw - lse)
  }
  list(
    particles = factor$particles,
    weights = weights
  )
}

build_population_factor_set <- function(local_objects, population_model) {
  population_model <- normalize_population_model(population_model)
  factors <- lapply(local_objects, build_population_local_factor, population_model = population_model)
  names(factors) <- names(local_objects)
  structure(
    list(
      population_model = population_model,
      factors = factors,
      n_locals = length(factors),
      log_constant = sum(vapply(factors, `[[`, numeric(1), "log_constant"))
    ),
    class = "population_factor_set"
  )
}

.evaluate_factor_batch <- function(factors, theta, include_constant = FALSE) {
  out <- rep.int(0, nrow(as.matrix(theta)))
  for (factor in factors) {
    out <- out + population_local_factor_log_marginal_many(factor, theta = theta, include_constant = include_constant)
  }
  out
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

  if (as.integer(n_cores) <= 1L || length(factor_set$factors) == 1L) {
    return(.evaluate_factor_batch(factor_set$factors, theta, include_constant = include_constant))
  }

  split_idx <- cut(seq_along(factor_set$factors), breaks = as.integer(min(n_cores, length(factor_set$factors))), labels = FALSE)
  batches <- split(factor_set$factors, split_idx)
  parts <- parallel::mclapply(
    batches,
    .evaluate_factor_batch,
    theta = theta,
    include_constant = FALSE,
    mc.cores = as.integer(min(n_cores, length(batches)))
  )
  out <- Reduce(`+`, parts)
  if (isTRUE(include_constant)) out <- out + factor_set$log_constant
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

.outer_sort_order <- function(theta) {
  theta <- as.matrix(theta)
  if (nrow(theta) <= 1L) return(seq_len(nrow(theta)))
  if (ncol(theta) == 1L) return(order(theta[, 1L]))
  theta_sc <- scale(theta)
  theta_sc[!is.finite(theta_sc)] <- 0
  order(rowSums(theta_sc))
}

.outer_rejuvenate <- function(theta,
                              logprior,
                              loglik_dynamic,
                              beta,
                              factor_set,
                              n_moves = 3L,
                              rw_scale = 0.8,
                              n_cores = 1L,
                              seed = NULL) {
  if (!is.null(seed)) set.seed(as.integer(seed))
  theta <- as.matrix(theta)
  N <- nrow(theta)
  d <- ncol(theta)
  model <- factor_set$population_model

  S <- tryCatch(stats::cov(theta), error = function(e) NULL)
  if (is.null(S) || any(!is.finite(S))) {
    S <- diag(d)
  } else {
    S <- regularize_cov(S, min_eig = 1e-8, cond_cap = 1e8)
  }
  Lrw <- tryCatch(chol(S), error = function(e) diag(d))
  step_scale <- as.numeric(rw_scale) / sqrt(max(d, 1))

  acc_total <- 0L
  prop_total <- 0L
  logpost <- logprior + beta * loglik_dynamic

  for (move_id in seq_len(max(1L, as.integer(n_moves)))) {
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

    if (length(accept)) {
      theta[accept, ] <- theta_prop[accept, , drop = FALSE]
      logprior[accept] <- logprior_prop[accept]
      loglik_dynamic[accept] <- loglik_prop[accept]
      logpost[accept] <- logpost_prop[accept]
      acc_total <- acc_total + length(accept)
    }
  }

  list(
    theta = theta,
    logprior = logprior,
    loglik_dynamic = loglik_dynamic,
    accept_rate = if (prop_total > 0L) acc_total / prop_total else 0
  )
}

outer_population_smc <- function(factor_set,
                                 N = 2000L,
                                 resample_threshold = 0.5,
                                 n_mcmc_moves = 3L,
                                 max_rounds = 100L,
                                 beta_target = 1.0,
                                 rw_scale_init = 0.8,
                                 target_accept = 0.234,
                                 rm_gain = 0.05,
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
    mcse_var_accum <- mcse_var_accum + max((mu2 - mu1^2) / (max(neff, 1) * max(mu1^2, .Machine$double.eps)), 0)

    logw_raw <- log(pmax(w, .Machine$double.eps)) + delta * (x - mx)
    lse <- logsumexp(logw_raw)
    log_evidence_dynamic <- log_evidence_dynamic + delta * mx + lse
    w <- exp(logw_raw - lse)
    w <- w / sum(w)
    ess_now <- ESS(w) / length(w)
    resampled <- FALSE

    vcat(sprintf("\nOuter round %d: beta %.3f -> %.3f | ESS=%.3f | logZ_dyn=%.4f\n",
                 round, beta, beta_new, ess_now, log_evidence_dynamic))

    beta <- beta_new
    if (ess_now < resample_threshold) {
      ord <- .outer_sort_order(theta)
      idx_sorted <- stratified_resample_sorted(w[ord], deterministic = FALSE)
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
      n_moves = n_mcmc_moves,
      rw_scale = exp(log_rw_scale),
      n_cores = n_cores,
      seed = seed + 1009L * round
    )
    theta <- rejuvenated$theta
    logprior <- rejuvenated$logprior
    loglik_dynamic <- rejuvenated$loglik_dynamic
    accept_rate <- rejuvenated$accept_rate

    log_rw_scale <- .clamp(log_rw_scale + rm_gain * (accept_rate - target_accept), log(0.05), log(2.0))

    beta_hist <- c(beta_hist, beta)
    ess_hist <- c(ess_hist, ess_now)
    accept_hist <- c(accept_hist, accept_rate)
    rw_scale_hist <- c(rw_scale_hist, exp(log_rw_scale))
    resampled_hist <- c(resampled_hist, resampled)

    vcat(sprintf("  Rejuvenation accept=%.3f | rw_scale=%.3f%s\n",
                 accept_rate, exp(log_rw_scale), if (resampled) " | resampled" else ""))
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
      resampled_hist = resampled_hist
    )
  )
}
