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
if (!exists("normalize_population_model", mode = "function") ||
    !exists("population_model_prepare_theta", mode = "function") ||
    !exists("population_model_log_alpha_given_prepared_theta_many", mode = "function")) {
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
      reference_prior = local_object$reference_prior,
      fast_family = as.character(population_model$fast_family %||% ""),
      fast_cache = NULL
    ),
    class = "population_local_factor"
  )
}

.pack_alpha_quadratic_terms_t <- function(alpha) {
  alpha <- as.matrix(alpha)
  d <- ncol(alpha)
  out <- matrix(0, nrow = d * (d + 1L) / 2L, ncol = nrow(alpha))
  pos <- 1L
  for (j in seq_len(d)) {
    aj <- alpha[, j]
    out[pos, ] <- aj * aj
    pos <- pos + 1L
    if (j < d) {
      for (k in (j + 1L):d) {
        out[pos, ] <- 2 * aj * alpha[, k]
        pos <- pos + 1L
      }
    }
  }
  out
}

.gaussian_log_alpha_given_prepared_theta_many <- function(alpha, theta_prepared) {
  alpha <- as.matrix(alpha)
  linear <- theta_prepared$eta %*% t(alpha)
  quad_t <- if (identical(theta_prepared$quadratic_kind, "diag")) {
    t(alpha * alpha)
  } else {
    .pack_alpha_quadratic_terms_t(alpha)
  }
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
                                        alpha,
                                        log_base,
                                        theta,
                                        theta_prepared = NULL,
                                        block_size = 1024L) {
  theta <- .as_hyper_matrix(theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  theta_prepared <- theta_prepared %||% population_model_prepare_theta(model, theta)
  block_size <- as.integer(max(1L, block_size))
  accum <- rep.int(-Inf, nrow(theta))
  n_alpha <- nrow(alpha)

  for (start in seq.int(1L, n_alpha, by = block_size)) {
    idx <- seq.int(start, min(start + block_size - 1L, n_alpha))
    log_terms <- .population_log_alpha_given_theta_many(
      model = model,
      alpha = alpha[idx, , drop = FALSE],
      theta_prepared = theta_prepared
    )
    block_lse <- .rowLogSumExp(sweep(log_terms, 2L, log_base[idx], "+"))
    accum <- rlogsumexp2(accum, block_lse)
  }
  as.numeric(accum)
}

population_local_factor_log_marginal_many <- function(factor, theta, include_constant = TRUE, block_size = 1024L) {
  stopifnot(inherits(factor, "population_local_factor"))
  theta_prepared <- population_model_prepare_theta(factor$population_model, theta)
  out <- .local_log_marginal_blocked(
    model = factor$population_model,
    alpha = factor$particles,
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
  logp <- .population_log_alpha_given_theta_many(
    model = factor$population_model,
    alpha = factor$particles,
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
  theta_prepared <- population_model_prepare_theta(factor$population_model, theta)
  logp <- .population_log_alpha_given_theta_many(
    model = factor$population_model,
    alpha = factor$particles,
    theta_prepared = theta_prepared
  )[1L, ]
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

build_population_factor_set <- function(local_objects,
                                        population_model,
                                        particle_block_size = 1024L) {
  population_model <- normalize_population_model(population_model)
  factors <- lapply(local_objects, build_population_local_factor, population_model = population_model)
  names(factors) <- names(local_objects)
  particle_block_size <- as.integer(max(1L, particle_block_size))
  if (length(factors)) {
    local_lengths <- vapply(factors, function(factor) nrow(factor$particles), integer(1))
    alpha <- do.call(rbind, lapply(factors, `[[`, "particles"))
    log_base <- unlist(lapply(factors, `[[`, "log_base"), use.names = FALSE)
    local_index <- rep.int(seq_along(factors), local_lengths)
    block_starts <- seq.int(1L, nrow(alpha), by = particle_block_size)
    blocks <- lapply(block_starts, function(start) {
      idx <- seq.int(start, min(start + particle_block_size - 1L, nrow(alpha)))
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
        log_base = log_base,
        local_index = local_index,
        blocks = blocks,
        particle_block_size = particle_block_size
      )
    ),
    class = "population_factor_set"
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
    log_terms <- .population_log_alpha_given_theta_many(
      model = model,
      alpha = factor_set$stack$alpha[block$idx, , drop = FALSE],
      theta_prepared = theta_prepared
    )
    log_terms <- sweep(log_terms, 2L, factor_set$stack$log_base[block$idx], "+")

    for (seg_id in seq_along(block$local)) {
      cols <- seq.int(block$starts[seg_id], block$ends[seg_id])
      seg_lse <- if (length(cols) == 1L) log_terms[, cols] else .rowLogSumExp(log_terms[, cols, drop = FALSE])
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

outer_population_smc <- function(factor_set,
                                 N = 2000L,
                                 resample_threshold = 0.5,
                                 n_mcmc_moves = 3L,
                                 min_mcmc_moves = 1L,
                                 max_rounds = 100L,
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

    vcat(sprintf("  Rejuvenation accept=%.3f | moves=%d | rw_scale=%.3f%s%s\n",
                 accept_rate,
                 rejuvenated$n_moves_used,
                 exp(log_rw_scale),
                 if (resampled) " | resampled" else "",
                 if (identical(sort_used, "not_used")) "" else paste0(" | sort=", sort_used)))
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
