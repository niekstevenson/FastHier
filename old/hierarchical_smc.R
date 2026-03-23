#!/usr/bin/env Rscript
# ============================================================================
# Hierarchical SMC built on top of the maintained local sampler.
#
# This is an initial pragmatic outer layer:
# - Each subject is summarized by a reusable bank of local posterior draws.
# - The outer sampler targets population hyperparameters phi.
# - Subject contributions are approximated by cached local-bank integrals.
#
# The implementation is deliberately simpler than an exact collapsed
# pseudo-marginal SMC. It is intended as a first hierarchical extension of the
# existing local SMC workflow, not as the final exact algorithm.
# ============================================================================

if (!exists("%||%", mode = "function") ||
    !exists("weighted_cov", mode = "function") ||
    !exists(".rowLogSumExp", mode = "function") ||
    !exists("cess_target_at_lambda", mode = "function")) {
  source("smc_core.R")
}

suppressPackageStartupMessages({
  library(Matrix)
  library(mvtnorm)
  library(parallel)
})

dnorm_log <- function(x, m, v) {
  -0.5 * (log(2 * pi * v) + (x - m)^2 / v)
}

dinvgamma_log <- function(x, a, b) {
  ifelse(x > 0, a * log(b) - lgamma(a) - (a + 1) * log(x) - b / x, -Inf)
}

.normalize_logweights <- function(logw) {
  logw <- as.numeric(logw)
  m <- max(logw)
  if (!is.finite(m)) {
    w <- rep(1 / length(logw), length(logw))
    return(list(log_norm = -Inf, w = w))
  }
  ww <- exp(logw - m)
  s <- sum(ww)
  if (!is.finite(s) || s <= 0) {
    w <- rep(1 / length(logw), length(logw))
    return(list(log_norm = -Inf, w = w))
  }
  list(log_norm = m + log(s), w = ww / s)
}

.ess_from_weights <- function(w) {
  w <- as.numeric(w)
  s <- sum(w)
  if (!is.finite(s) || s <= 0) return(0)
  w <- w / s
  1 / sum(w * w)
}

.resample_indices <- function(w) {
  sample.int(length(w), size = length(w), replace = TRUE, prob = w)
}

.regularize_cov_safe <- function(S, jitter = 1e-8) {
  S <- as.matrix((S + t(S)) / 2)
  S <- tryCatch(as.matrix(Matrix::nearPD(S, conv.tol = 1e-7)$mat), error = function(e) S)
  d <- nrow(S)
  S + diag(jitter, d)
}

.weighted_mean_vec <- function(X, w) {
  X <- as.matrix(X)
  w <- pmax(as.numeric(w), 0)
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) {
    w <- rep(1 / nrow(X), nrow(X))
  } else {
    w <- w / sw
  }
  colSums(X * w)
}

.weighted_cov_mat <- function(X, w) {
  X <- as.matrix(X)
  w <- pmax(as.numeric(w), 0)
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) {
    w <- rep(1 / nrow(X), nrow(X))
  } else {
    w <- w / sw
  }
  mu <- .weighted_mean_vec(X, w)
  Xc <- sweep(X, 2L, mu, `-`)
  S <- crossprod(sqrt(w) * Xc)
  .regularize_cov_safe(S)
}

make_working_prior_gaussian <- function(mu, Sigma, phi_anchor = NULL) {
  Sigma <- .regularize_cov_safe(Sigma)
  Sigma_inv <- chol2inv(chol(Sigma))
  list(
    mu = as.numeric(mu),
    Sigma = Sigma,
    Sigma_inv = Sigma_inv,
    logdet = as.numeric(determinant(Sigma, logarithm = TRUE)$modulus),
    phi_anchor = phi_anchor
  )
}

.normalize_working_prior_gaussian <- function(working_prior, d_theta = NULL) {
  if (is.null(working_prior)) stop("working_prior must not be NULL.")
  mu <- as.numeric(working_prior$mu %||% working_prior$mean)
  if (!length(mu)) stop("working_prior must contain 'mu' or 'mean'.")
  if (!is.null(d_theta) && length(mu) != d_theta) {
    stop("working_prior mean length does not match d_theta.")
  }

  Sigma <- working_prior$Sigma
  Sigma_inv <- working_prior$Sigma_inv
  if (is.null(Sigma) && is.null(Sigma_inv)) {
    stop("working_prior must contain either 'Sigma' or 'Sigma_inv'.")
  }

  if (is.null(Sigma)) {
    Sigma_inv <- as.matrix(Sigma_inv)
    L <- tryCatch(chol(Sigma_inv), error = function(e) NULL)
    if (is.null(L)) stop("working_prior$Sigma_inv is not SPD.")
    Sigma <- chol2inv(L)
  } else {
    Sigma <- .regularize_cov_safe(Sigma)
  }

  if (is.null(Sigma_inv)) {
    L <- tryCatch(chol(Sigma), error = function(e) NULL)
    if (is.null(L)) stop("working_prior$Sigma is not SPD.")
    Sigma_inv <- chol2inv(L)
  } else {
    Sigma_inv <- as.matrix(Sigma_inv)
  }

  logdet <- as.numeric(working_prior$logdet %||% determinant(Sigma, logarithm = TRUE)$modulus)
  structure(
    list(
      mu = mu,
      Sigma = Sigma,
      Sigma_inv = Sigma_inv,
      logdet = logdet,
      phi_anchor = working_prior$phi_anchor %||% NULL
    ),
    class = "working_prior_gaussian"
  )
}

.log_prior_gauss_mat <- function(Theta, pars) {
  Theta <- as.matrix(Theta)
  mu <- as.numeric(pars$mu)
  K <- as.matrix(pars$Sigma_inv)
  if (ncol(Theta) == 1L && length(mu) == 1L && nrow(K) == 1L && ncol(K) == 1L) {
    lt1 <- Theta[, 1L] - mu[1L]
    q <- as.numeric(K[1L, 1L]) * lt1 * lt1
    out <- -0.5 * (q + pars$logdet + log(2 * pi))
    if (!is.null(pars$const)) out <- out + pars$const
    return(as.numeric(out))
  }
  lt <- sweep(Theta, 2L, mu, `-`)
  q <- rowSums((lt %*% K) * lt)
  out <- -0.5 * (q + pars$logdet + ncol(Theta) * log(2 * pi))
  if (!is.null(pars$const)) out <- out + pars$const
  as.numeric(out)
}

log_working_prior_gaussian_mat <- function(Theta, working_prior) {
  wp <- .normalize_working_prior_gaussian(working_prior, d_theta = ncol(as.matrix(Theta)))
  .log_prior_gauss_mat(Theta, wp)
}

make_prior_phi_diag <- function(m0, s0, a, b, d, param_names = NULL) {
  m0 <- rep_len(as.numeric(m0), d)
  s0 <- rep_len(as.numeric(s0), d)
  a <- rep_len(as.numeric(a), d)
  b <- rep_len(as.numeric(b), d)
  if (any(s0 <= 0)) stop("s0 must be strictly positive.")
  if (any(a <= 1)) stop("Inverse-gamma shape must exceed 1 for a finite mean.")
  if (is.null(param_names)) param_names <- paste0("theta", seq_len(d))
  mu_names <- paste0("mu_", param_names)
  log_sigma2_names <- paste0("log_sigma2_", param_names)

  rprior_phi <- function(n) {
    mu <- matrix(
      rnorm(n * d, mean = rep(m0, each = n), sd = sqrt(rep(s0, each = n))),
      nrow = n,
      ncol = d,
      byrow = FALSE
    )
    sigma2 <- matrix(
      1 / rgamma(n * d, shape = rep(a, each = n), rate = rep(b, each = n)),
      nrow = n,
      ncol = d,
      byrow = FALSE
    )
    out <- cbind(mu, log(sigma2))
    colnames(out) <- c(mu_names, log_sigma2_names)
    out
  }

  logprior_phi <- function(phi_row) {
    phi_row <- as.numeric(phi_row)
    mu <- phi_row[seq_len(d)]
    ell <- phi_row[d + seq_len(d)]
    sigma2 <- exp(ell)
    sum(dnorm_log(mu, m0, s0)) +
      sum(dinvgamma_log(sigma2, a, b)) +
      sum(ell)
  }

  list(
    rprior = rprior_phi,
    lprior = logprior_phi,
    d = d,
    param_names = param_names
  )
}

default_population_prior_diag <- function(mu_ref,
                                          Sigma_ref,
                                          mean_var_scale = 1.0,
                                          sigma2_shape = 3.0,
                                          sigma2_mean = NULL) {
  mu_ref <- as.numeric(mu_ref)
  Sigma_ref <- as.matrix(Sigma_ref)
  d <- length(mu_ref)
  if (nrow(Sigma_ref) != d || ncol(Sigma_ref) != d) {
    stop("Sigma_ref dimensions do not match mu_ref.")
  }
  var_diag <- pmax(diag(Sigma_ref), 1e-8)
  sigma2_mean <- sigma2_mean %||% var_diag
  sigma2_mean <- rep_len(as.numeric(sigma2_mean), d)
  shape <- rep_len(as.numeric(sigma2_shape), d)
  if (any(shape <= 1)) stop("sigma2_shape must exceed 1.")
  rate <- sigma2_mean * (shape - 1)
  make_prior_phi_diag(
    m0 = mu_ref,
    s0 = mean_var_scale * var_diag,
    a = shape,
    b = rate,
    d = d,
    param_names = names(mu_ref) %||% colnames(Sigma_ref)
  )
}

phi_to_gaussian_params_diag_factory <- function(param_names = NULL) {
  function(phi, d) {
    stopifnot(length(phi) >= 2 * d)
    mu <- as.numeric(phi[seq_len(d)])
    ell <- as.numeric(phi[d + seq_len(d)])
    sigma2 <- pmax(exp(ell), 1e-10)
    mu_names <- param_names
    if (is.null(mu_names) && !is.null(names(phi))) {
      mu_names <- sub("^mu_", "", names(phi)[seq_len(d)])
    }
    if (!is.null(mu_names) && length(mu_names) == d) {
      names(mu) <- mu_names
    }
    list(
      mu = mu,
      Sigma_inv = diag(1 / sigma2, d),
      logdet = sum(log(sigma2)),
      const = NULL
    )
  }
}

.gaussian_pars_to_mean_cov <- function(pars, theta_names = NULL) {
  Sigma_inv <- as.matrix(pars$Sigma_inv)
  L <- tryCatch(chol(Sigma_inv), error = function(e) NULL)
  if (is.null(L)) stop("gaussian_map_fn returned non-SPD precision.")
  Sigma <- chol2inv(L)
  Sigma <- .regularize_cov_safe(Sigma)
  mu <- as.numeric(pars$mu)
  nm <- theta_names %||% names(pars$mu)
  if (!is.null(nm) && length(nm) == length(mu)) {
    names(mu) <- nm
    dimnames(Sigma) <- list(nm, nm)
  }
  list(mu = mu, Sigma = Sigma)
}

.call_local_log_evidence_fn <- function(local_log_evidence_fn,
                                        phi,
                                        data_i,
                                        lambda,
                                        checkpoint_lambdas = NULL,
                                        seed_plan = NULL,
                                        auxiliary_state = NULL) {
  fn_formals <- names(formals(local_log_evidence_fn) %||% pairlist())
  accepts_lambda <- "lambda" %in% fn_formals || "..." %in% fn_formals
  if (!accepts_lambda && abs(as.numeric(lambda) - 1.0) > 1e-12) {
    stop("local_log_evidence_fn must accept a 'lambda' argument (or ...) for bridge evaluations.")
  }
  args <- list(phi = phi, data_i = data_i)
  if (accepts_lambda) {
    args$lambda <- as.numeric(lambda)
  }
  if ("checkpoint_lambdas" %in% fn_formals || "..." %in% fn_formals) {
    args$checkpoint_lambdas <- checkpoint_lambdas
  }
  if ("seed_plan" %in% fn_formals || "..." %in% fn_formals) {
    args$seed_plan <- seed_plan
  }
  if ("auxiliary_state" %in% fn_formals || "..." %in% fn_formals) {
    args$auxiliary_state <- auxiliary_state
  }
  do.call(local_log_evidence_fn, args)
}

.estimate_local_log_evidence_at_phi <- function(phi,
                                                data_i,
                                                loglik_fn,
                                                gaussian_map_fn,
                                                d_theta,
                                                local_particles,
                                                local_log_evidence_fn = NULL,
                                                lambda = 1.0,
                                                theta_names = NULL,
                                                local_n_cores = 1L,
                                                seed = NULL,
                                                seed_plan = NULL,
                                                checkpoint_lambdas = NULL,
                                                auxiliary_state = NULL,
                                                verbose = FALSE,
                                                warm_start_fit = NULL,
                                                local_smc_control = list()) {
  if (is.function(local_log_evidence_fn)) {
    raw <- .call_local_log_evidence_fn(
        local_log_evidence_fn = local_log_evidence_fn,
        phi = phi,
        data_i = data_i,
        lambda = lambda,
        checkpoint_lambdas = checkpoint_lambdas,
        seed_plan = seed_plan,
        auxiliary_state = auxiliary_state
      )
    if (is.list(raw) && !is.null(raw$log_evidence)) {
      return(list(
        log_evidence = as.numeric(raw$log_evidence),
        local_fit = raw$local_fit %||% NULL,
        checkpoint_lambdas = as.numeric(raw$checkpoint_lambdas %||% checkpoint_lambdas %||% numeric(0)),
        checkpoint_log_evidence = as.numeric(raw$checkpoint_log_evidence %||% numeric(0)),
        checkpoint_log_increments = as.numeric(raw$checkpoint_log_increments %||% numeric(0)),
        auxiliary_state = raw$auxiliary_state %||% auxiliary_state
      ))
    }
    return(list(
      log_evidence = as.numeric(raw),
      local_fit = NULL,
      checkpoint_lambdas = as.numeric(checkpoint_lambdas %||% numeric(0)),
      checkpoint_log_evidence = numeric(0),
      checkpoint_log_increments = numeric(0),
      auxiliary_state = auxiliary_state
    ))
  }
  if (!exists("enhanced_smc_elite", mode = "function")) {
    source("SMC_super_fast.R")
  }
  if (!is.null(seed)) set.seed(as.integer(seed))
  pars <- gaussian_map_fn(phi, d_theta)
  moments <- .gaussian_pars_to_mean_cov(pars, theta_names = theta_names)
  # Nested local evidence estimation is more stable with a conservative inner SMC.
  nested_defaults <- list(
    gss_enable = FALSE,
    hist_mix_enable = FALSE,
    da_enable = FALSE
  )
  call_args <- modifyList(
    list(
      data = data_i,
      loglik_fn = loglik_fn,
      mu_ref = moments$mu,
      Sigma_ref = moments$Sigma,
      M = as.integer(local_particles),
      lambda_target = as.numeric(lambda),
      checkpoint_lambdas = checkpoint_lambdas,
      n_cores = as.integer(local_n_cores),
      seed = as.integer(seed %||% 123L),
      seed_plan = seed_plan,
      verbose = verbose,
      warm_start_fit = warm_start_fit
    ),
    modifyList(nested_defaults, local_smc_control)
  )
  smc_out <- do.call(enhanced_smc_elite, call_args)
  list(
    log_evidence = as.numeric(smc_out$log_evidence),
    local_fit = smc_out,
    checkpoint_lambdas = as.numeric(smc_out$checkpoint_lambdas %||% numeric(0)),
    checkpoint_log_evidence = as.numeric(smc_out$checkpoint_log_evidence %||% numeric(0)),
    checkpoint_log_increments = as.numeric(smc_out$checkpoint_log_increments %||% numeric(0)),
    auxiliary_state = list(
      seed_plan = smc_out$seed_plan %||% seed_plan
    )
  )
}

run_local_smc_subjects <- function(data_list,
                                   loglik_fn,
                                   mu_ref,
                                   Sigma_ref,
                                   M = 4000L,
                                   n_cores = 1L,
                                   base_seed = 123L,
                                   verbose = TRUE,
                                   ...) {
  if (!exists("enhanced_smc_elite", mode = "function")) {
    source("SMC_super_fast.R")
  }
  parallel::mclapply(
    seq_along(data_list),
    function(i) {
      enhanced_smc_elite(
        data = data_list[[i]],
        loglik_fn = loglik_fn,
        mu_ref = mu_ref,
        Sigma_ref = Sigma_ref,
        M = as.integer(M),
        n_cores = 1L,
        seed = as.integer(base_seed + i - 1L),
        verbose = verbose,
        ...
      )
    },
    mc.cores = as.integer(max(1L, n_cores))
  )
}

build_local_exact_object <- function(smc_out,
                                     data = NULL,
                                     subj_id = NA_integer_,
                                     loglik_fn = NULL,
                                     working_prior = NULL,
                                     base_seed = NULL) {
  if (is.null(smc_out$Theta) || is.null(smc_out$w)) {
    stop("smc_out must contain Theta and w.")
  }
  Theta <- as.matrix(smc_out$Theta)
  d_theta <- ncol(Theta)
  wp <- .normalize_working_prior_gaussian(working_prior %||% smc_out$working_prior, d_theta = d_theta)
  w <- pmax(as.numeric(smc_out$w), 0)
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) {
    w <- rep(1 / nrow(Theta), nrow(Theta))
  } else {
    w <- w / sw
  }
  structure(
    list(
      subj_id = as.integer(subj_id),
      data = data,
      loglik_fn = loglik_fn,
      working_prior = wp,
      base_seed = as.integer(base_seed %||% (10007L + 1009L * as.integer(subj_id))),
      d_theta = d_theta,
      theta_bank = Theta,
      theta_bank_w = w,
      theta_bank_loglik = if (!is.null(smc_out$loglik)) as.numeric(smc_out$loglik) else NULL,
      local_fit_diag = list(
        rounds = smc_out$meta$rounds %||% length(smc_out$meta$lambda_hist %||% numeric(0)),
        log_evidence = smc_out$log_evidence %||% NA_real_,
        mcse_log_evidence = smc_out$mcse_logZ %||% NA_real_,
        final_lambda = smc_out$final_lambda %||% 1.0
      )
    ),
    class = "local_exact_object"
  )
}

build_local_exact_objects <- function(local_fits,
                                      data_list,
                                      loglik_fn,
                                      working_prior,
                                      base_seed = 40000L) {
  if (length(local_fits) != length(data_list)) {
    stop("local_fits and data_list must have the same length.")
  }
  lapply(
    seq_along(local_fits),
    function(i) {
      build_local_exact_object(
        smc_out = local_fits[[i]],
        data = data_list[[i]],
        subj_id = i,
        loglik_fn = loglik_fn,
        working_prior = working_prior,
        base_seed = as.integer(base_seed + 1009L * i)
      )
    }
  )
}

.new_local_bank <- function(alpha, w, logg, local_obj, loglik = NULL,
                            source_idx = NULL, proposal_tag = "prefit_gamma0") {
  alpha <- as.matrix(alpha)
  bank <- structure(
    list(
      alpha = alpha,
      w = as.numeric(w),
      logg = as.numeric(logg),
      loglik = if (is.null(loglik)) NULL else as.numeric(loglik),
      source_idx = if (is.null(source_idx)) NULL else as.integer(source_idx),
      d_theta = ncol(alpha),
      local_obj = local_obj,
      proposal_tag = as.character(proposal_tag)
    ),
    class = "local_bank"
  )
  .validate_local_bank(bank)
}

.validate_local_bank <- function(bank) {
  stopifnot(inherits(bank, "local_bank"))
  A <- as.matrix(bank$alpha)
  n <- nrow(A)
  if (length(bank$w) != n) stop("local_bank has inconsistent weight length.")
  if (length(bank$logg) != n) stop("local_bank has inconsistent logg length.")
  if (!is.null(bank$loglik) && length(bank$loglik) != n) {
    stop("local_bank has inconsistent loglik length.")
  }
  if (!is.null(bank$source_idx) && length(bank$source_idx) != n) {
    stop("local_bank has inconsistent source_idx length.")
  }
  sw <- sum(bank$w)
  if (!is.finite(sw) || sw <= 0) stop("local_bank weights are invalid.")
  bank$alpha <- A
  bank$w <- bank$w / sw
  bank$logw <- log(pmax(bank$w, .Machine$double.eps))
  bank
}

build_local_bank <- function(local_obj, M_local = NULL, seed = 123L, proposal_tag = "prefit_gamma0") {
  stopifnot(inherits(local_obj, "local_exact_object"))
  if (!is.null(seed)) set.seed(as.integer(seed))

  Theta_all <- as.matrix(local_obj$theta_bank)
  n_all <- nrow(Theta_all)
  if (n_all <= 0L) stop("local_obj$theta_bank is empty.")

  w_all <- pmax(as.numeric(local_obj$theta_bank_w), 0)
  sw <- sum(w_all)
  if (!is.finite(sw) || sw <= 0) {
    w_all <- rep(1 / n_all, n_all)
  } else {
    w_all <- w_all / sw
  }

  if (is.null(M_local) || as.integer(M_local) >= n_all) {
    idx <- seq_len(n_all)
    Theta <- Theta_all
    w_bank <- w_all
  } else {
    idx <- sample.int(n_all, size = as.integer(M_local), replace = TRUE, prob = w_all)
    Theta <- Theta_all[idx, , drop = FALSE]
    w_bank <- rep(1 / length(idx), length(idx))
  }

  logg <- log_working_prior_gaussian_mat(Theta, local_obj$working_prior)
  loglik <- NULL
  if (!is.null(local_obj$theta_bank_loglik)) {
    loglik <- as.numeric(local_obj$theta_bank_loglik[idx])
  }

  .new_local_bank(
    alpha = Theta,
    w = w_bank,
    logg = logg,
    loglik = loglik,
    source_idx = idx,
    local_obj = local_obj,
    proposal_tag = proposal_tag
  )
}

.local_bank_delta <- function(bank, phi, gaussian_map_fn, validate = TRUE) {
  if (isTRUE(validate)) bank <- .validate_local_bank(bank)
  logp <- .log_prior_gauss_mat(bank$alpha, gaussian_map_fn(as.numeric(phi), bank$d_theta))
  as.numeric(logp - bank$logg)
}

.local_bank_weight_stats <- function(bank, delta, lambda) {
  logw_lambda <- bank$logw + as.numeric(lambda) * as.numeric(delta)
  norm <- .normalize_logweights(logw_lambda)
  list(
    logZ = as.numeric(norm$log_norm),
    ess = .ess_from_weights(norm$w)
  )
}

.estimate_logZ_validated <- function(bank, phi, lambda, gaussian_map_fn) {
  delta <- .local_bank_delta(bank, phi = phi, gaussian_map_fn = gaussian_map_fn, validate = FALSE)
  .local_bank_weight_stats(bank, delta = delta, lambda = lambda)$logZ
}

.estimate_logZ_bank_list <- function(banks, phi, lambda, gaussian_map_fn) {
  vapply(
    banks,
    function(bank) .estimate_logZ_validated(bank, phi = phi, lambda = lambda, gaussian_map_fn = gaussian_map_fn),
    numeric(1L)
  )
}

.prepare_bank_delta_matrix <- function(bank, phi, gaussian_map_fn) {
  bank <- .validate_local_bank(bank)
  phi <- as.matrix(phi)
  N <- nrow(phi)
  M <- nrow(bank$alpha)
  delta_mat <- matrix(0, nrow = N, ncol = M)
  for (n in seq_len(N)) {
    delta_mat[n, ] <- .local_bank_delta(
      bank,
      phi = phi[n, , drop = TRUE],
      gaussian_map_fn = gaussian_map_fn,
      validate = FALSE
    )
  }
  list(bank = bank, delta = delta_mat)
}

.prepare_population_delta_cache <- function(banks, phi, gaussian_map_fn, lambda_ref = NULL) {
  phi <- as.matrix(phi)
  N <- nrow(phi)
  lapply(
    banks,
    function(bank) {
      bank <- .validate_local_bank(bank)
      M <- nrow(bank$alpha)
      delta_mat <- matrix(0, nrow = N, ncol = M)
      for (n in seq_len(N)) {
        delta_mat[n, ] <- .local_bank_delta(
          bank,
          phi = phi[n, , drop = TRUE],
          gaussian_map_fn = gaussian_map_fn,
          validate = FALSE
        )
      }
      cache <- list(
        bank = bank,
        delta = delta_mat,
        logw_rep = matrix(rep(bank$logw, each = N), nrow = N),
        lambda_ref = if (is.null(lambda_ref)) NULL else as.numeric(lambda_ref)
      )
      if (!is.null(lambda_ref)) {
        cache$logZ_ref <- .rowLogSumExp(as.numeric(lambda_ref) * delta_mat + cache$logw_rep)
      }
      cache
    }
  )
}

.logZ_matrix_from_cache <- function(cache, lambda) {
  lambda <- as.numeric(lambda)
  do.call(
    cbind,
    lapply(
      cache,
      function(cache_i) {
        if (!is.null(cache_i$logZ_ref) && !is.null(cache_i$lambda_ref) &&
            abs(lambda - cache_i$lambda_ref) < 1e-15) {
          return(as.numeric(cache_i$logZ_ref))
        }
        .rowLogSumExp(lambda * cache_i$delta + cache_i$logw_rep)
      }
    )
  )
}

.predict_population_log_increment_from_cache <- function(cache, lambda_old, lambda_new) {
  lambda_old <- as.numeric(lambda_old)
  lambda_new <- as.numeric(lambda_new)
  if (abs(lambda_new - lambda_old) < 1e-15) {
    return(rep.int(0, nrow(cache[[1L]]$delta)))
  }
  N <- nrow(cache[[1L]]$delta)
  out <- numeric(N)
  for (i in seq_along(cache)) {
    cache_i <- cache[[i]]
    logZ_old <- cache_i$logZ_ref
    if (is.null(logZ_old)) {
      logZ_old <- .rowLogSumExp(lambda_old * cache_i$delta + cache_i$logw_rep)
    }
    logZ_new <- .rowLogSumExp(lambda_new * cache_i$delta + cache_i$logw_rep)
    out <- out + (logZ_new - logZ_old)
  }
  out
}

.rcess_from_log_increment <- function(w, log_inc) {
  w <- as.numeric(w)
  w[!is.finite(w) | w < 0] <- 0
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) return(0)
  w <- w / sw

  log_inc <- as.numeric(log_inc)
  ok <- is.finite(log_inc)
  if (!any(ok)) return(0)
  log_inc[!ok] <- min(log_inc[ok])

  a1 <- logsumexp(log(w) + log_inc)
  a2 <- logsumexp(log(w) + 2 * log_inc)
  val <- exp(2 * a1 - a2)
  if (!is.finite(val)) return(0)
  pmin(pmax(val, 0), 1)
}

choose_next_lambda_nested <- function(phi,
                                      banks,
                                      w,
                                      lambda_curr,
                                      rho_step,
                                      gaussian_map_fn,
                                      step_target_fn = cess_target_at_lambda,
                                      increment_cache = NULL,
                                      tol = 1e-4,
                                      max_iter = 30L) {
  if (is.null(increment_cache)) {
    increment_cache <- .prepare_population_delta_cache(
      banks = banks,
      phi = phi,
      gaussian_map_fn = gaussian_map_fn,
      lambda_ref = lambda_curr
    )
  }
  rcess_at <- function(lambda_new) {
    log_inc <- .predict_population_log_increment_from_cache(
      cache = increment_cache,
      lambda_old = lambda_curr,
      lambda_new = lambda_new
    )
    .rcess_from_log_increment(w, log_inc)
  }

  target_frac <- if (is.null(rho_step)) {
    as.numeric(step_target_fn(lambda_curr))
  } else {
    as.numeric(rho_step)
  }
  target_frac <- pmin(pmax(target_frac, 0), 0.999)

  rcess_one <- rcess_at(1.0)
  if (rcess_one >= target_frac) return(1.0)

  lo <- lambda_curr
  hi <- 1.0
  for (iter in seq_len(max_iter)) {
    mid <- 0.5 * (lo + hi)
    rcess_mid <- rcess_at(mid)
    if (rcess_mid >= target_frac) {
      lo <- mid
    } else {
      hi <- mid
    }
    if ((hi - lo) < tol) break
  }

  min(1.0, max(lo, lambda_curr + min(1e-4, 1.0 - lambda_curr)))
}

.bank_ess_values_from_delta <- function(bank, delta_mat, lambda) {
  logw_lambda <- sweep(as.numeric(lambda) * delta_mat, 2L, bank$logw, "+")
  log_norm <- .rowLogSumExp(logw_lambda)
  log_sum_sq <- .rowLogSumExp(2 * logw_lambda)
  as.numeric(exp(2 * log_norm - log_sum_sq))
}

.bank_split_diff_from_delta <- function(bank, delta_mat, lambda, idx_odd, idx_even) {
  if (length(idx_odd) < 4L || length(idx_even) < 4L) {
    return(rep.int(Inf, nrow(delta_mat)))
  }

  w_odd <- bank$w[idx_odd]
  sw_odd <- sum(w_odd)
  w_even <- bank$w[idx_even]
  sw_even <- sum(w_even)
  if (!is.finite(sw_odd) || sw_odd <= 0 || !is.finite(sw_even) || sw_even <= 0) {
    return(rep.int(Inf, nrow(delta_mat)))
  }

  logw_odd <- log(pmax(w_odd / sw_odd, .Machine$double.eps))
  logw_even <- log(pmax(w_even / sw_even, .Machine$double.eps))
  logZ_odd <- .rowLogSumExp(sweep(as.numeric(lambda) * delta_mat[, idx_odd, drop = FALSE], 2L, logw_odd, "+"))
  logZ_even <- .rowLogSumExp(sweep(as.numeric(lambda) * delta_mat[, idx_even, drop = FALSE], 2L, logw_even, "+"))
  as.numeric(abs(logZ_odd - logZ_even))
}

.bank_group_diagnostics <- function(bank,
                                    phi,
                                    lambda,
                                    gaussian_map_fn,
                                    ess_prob = 0.1,
                                    split_prob = 0.9) {
  prep <- .prepare_bank_delta_matrix(bank, phi = phi, gaussian_map_fn = gaussian_map_fn)
  idx_odd <- seq.int(1L, nrow(bank$alpha), by = 2L)
  idx_even <- seq.int(2L, nrow(bank$alpha), by = 2L)

  ess_vals <- .bank_ess_values_from_delta(prep$bank, delta_mat = prep$delta, lambda = lambda)
  split_diff <- .bank_split_diff_from_delta(
    prep$bank,
    delta_mat = prep$delta,
    lambda = lambda,
    idx_odd = idx_odd,
    idx_even = idx_even
  )

  list(
    M = nrow(bank$alpha),
    ess_min = min(ess_vals, na.rm = TRUE),
    ess_q = as.numeric(stats::quantile(ess_vals, probs = ess_prob, names = FALSE, na.rm = TRUE)),
    split_q = as.numeric(stats::quantile(split_diff, probs = split_prob, names = FALSE, na.rm = TRUE))
  )
}

.bank_ess_summary <- function(banks,
                              phi,
                              lambda,
                              gaussian_map_fn,
                              probs = c(0.1, 0.5),
                              increment_cache = NULL) {
  cache <- increment_cache
  if (is.null(cache)) {
    cache <- .prepare_population_delta_cache(
      banks = banks,
      phi = phi,
      gaussian_map_fn = gaussian_map_fn,
      lambda_ref = NULL
    )
  }
  if (!length(cache)) {
    return(list(min = NA_real_, quantiles = setNames(rep(NA_real_, length(probs)), paste0("q", probs * 100))))
  }

  N <- nrow(cache[[1L]]$delta)
  ess_mat <- vapply(
    cache,
    function(cache_i) .bank_ess_values_from_delta(cache_i$bank, delta_mat = cache_i$delta, lambda = lambda),
    numeric(N)
  )
  ess_vals <- as.numeric(ess_mat)
  qq <- stats::quantile(ess_vals, probs = probs, names = FALSE, na.rm = TRUE)
  list(
    min = min(ess_vals, na.rm = TRUE),
    quantiles = setNames(as.numeric(qq), paste0("q", probs * 100))
  )
}

.topup_local_bank_from_prefit <- function(bank, n_new, seed = NULL) {
  bank <- .validate_local_bank(bank)
  if (!is.null(seed)) set.seed(as.integer(seed))
  local_obj <- bank$local_obj

  Theta_all <- as.matrix(local_obj$theta_bank)
  n_all <- nrow(Theta_all)
  used <- unique(bank$source_idx %||% integer(0))
  remaining <- setdiff(seq_len(n_all), used)
  if (!length(remaining) || as.integer(n_new) <= 0L) {
    return(list(bank = bank, n_new = 0L))
  }

  n_new <- min(as.integer(n_new), length(remaining))
  w_rem <- pmax(as.numeric(local_obj$theta_bank_w[remaining]), 0)
  sw <- sum(w_rem)
  if (!is.finite(sw) || sw <= 0) {
    w_rem <- rep(1 / length(remaining), length(remaining))
  } else {
    w_rem <- w_rem / sw
  }
  pick <- sample.int(length(remaining), size = n_new, replace = FALSE, prob = w_rem)
  idx_new <- remaining[pick]

  alpha_new <- Theta_all[idx_new, , drop = FALSE]
  logg_new <- log_working_prior_gaussian_mat(alpha_new, local_obj$working_prior)
  loglik_new <- NULL
  if (!is.null(local_obj$theta_bank_loglik)) {
    loglik_new <- as.numeric(local_obj$theta_bank_loglik[idx_new])
  }

  alpha_all <- rbind(bank$alpha, alpha_new)
  logg_all <- c(bank$logg, logg_new)
  loglik_all <- if (is.null(bank$loglik) && is.null(loglik_new)) NULL else c(bank$loglik, loglik_new)
  source_all <- c(bank$source_idx %||% integer(0), idx_new)
  w_all <- rep(1 / nrow(alpha_all), nrow(alpha_all))

  list(
    bank = .new_local_bank(
      alpha = alpha_all,
      w = w_all,
      logg = logg_all,
      loglik = loglik_all,
      source_idx = source_all,
      local_obj = local_obj,
      proposal_tag = "prefit_gamma0_topped_up"
    ),
    n_new = as.integer(n_new)
  )
}

enrich_local_bank <- function(bank,
                              phi,
                              lambda,
                              gaussian_map_fn,
                              rho_local = 0.5,
                              split_tol = 0.03,
                              max_topups = 2L,
                              topup_batch = NULL,
                              seed = NULL) {
  bank <- .validate_local_bank(bank)
  topup_points <- 0L
  enriched <- FALSE

  max_topups <- as.integer(max(0L, max_topups))
  if (max_topups <= 0L) {
    return(list(bank = bank, n_new = topup_points, enriched = enriched))
  }

  for (step in seq_len(max_topups)) {
    diag_now <- .bank_group_diagnostics(
      bank,
      phi = phi,
      lambda = lambda,
      gaussian_map_fn = gaussian_map_fn
    )
    need_ess <- is.finite(diag_now$ess_q) && diag_now$ess_q < rho_local * diag_now$M
    need_split <- !is.finite(diag_now$split_q) || diag_now$split_q > split_tol
    if (!(need_ess || need_split)) break

    batch <- topup_batch
    if (is.null(batch)) {
      batch <- max(8L, min(32L, ceiling(0.5 * nrow(bank$alpha))))
    }
    top <- .topup_local_bank_from_prefit(
      bank = bank,
      n_new = batch,
      seed = if (is.null(seed)) NULL else as.integer(seed + step)
    )
    if (top$n_new <= 0L) break
    bank <- top$bank
    topup_points <- topup_points + top$n_new
    enriched <- TRUE
  }

  list(bank = bank, n_new = as.integer(topup_points), enriched = enriched)
}

mutate_population_particle_nested <- function(phi_row,
                                              logprior_row,
                                              lambda,
                                              banks,
                                              logprior_phi,
                                              gaussian_map_fn,
                                              rw_cov,
                                              current_logZ = NULL,
                                              seed = NULL) {
  if (!is.null(seed)) set.seed(as.integer(seed))
  prop <- as.numeric(phi_row + mvtnorm::rmvnorm(1L, sigma = rw_cov))
  lp_prop <- logprior_phi(prop)
  if (!is.finite(lp_prop)) {
    return(list(phi = phi_row, logprior = logprior_row, accepted = FALSE, logZ = current_logZ))
  }

  log_ratio <- lp_prop - logprior_row
  prop_logZ <- .estimate_logZ_bank_list(banks, phi = prop, lambda = lambda, gaussian_map_fn = gaussian_map_fn)
  curr_logZ <- current_logZ
  if (is.null(curr_logZ)) {
    curr_logZ <- .estimate_logZ_bank_list(banks, phi = phi_row, lambda = lambda, gaussian_map_fn = gaussian_map_fn)
  }
  log_ratio <- log_ratio + sum(prop_logZ - curr_logZ)

  accepted <- is.finite(log_ratio) && (log(runif(1L)) < min(0, log_ratio))
  if (!accepted) {
    return(list(phi = phi_row, logprior = logprior_row, accepted = FALSE, logZ = curr_logZ))
  }

  list(phi = prop, logprior = lp_prop, accepted = TRUE, logZ = prop_logZ)
}

nested_population_smc <- function(rprior_phi,
                                  local_objs,
                                  logprior_phi,
                                  gaussian_map_fn,
                                  N,
                                  M_local,
                                  rho_step = NULL,
                                  rho_res = 0.5,
                                  rho_local = 0.5,
                                  n_population_moves = 1L,
                                  max_rounds = 100L,
                                  max_bank_topups = 2L,
                                  bank_split_tol = 0.03,
                                  bank_topup_batch = NULL,
                                  step_target_fn = cess_target_at_lambda,
                                  seed = 123L,
                                  verbose = TRUE) {
  vcat <- function(...) if (isTRUE(verbose)) cat(...)
  set.seed(as.integer(seed))

  N <- as.integer(N)
  phi <- as.matrix(rprior_phi(N))
  if (nrow(phi) != N) stop("rprior_phi returned wrong number of particles.")
  w <- rep(1 / N, N)
  logprior <- apply(phi, 1L, logprior_phi)

  banks <- lapply(
    seq_along(local_objs),
    function(i) build_local_bank(local_objs[[i]], M_local = M_local, seed = as.integer(seed + 1009L * i))
  )

  lambda <- 0
  round <- 0L
  log_evidence <- 0
  lambda_hist <- c(0)
  ess_hist <- c(N)
  bank_ess_min_hist <- c(M_local)
  bank_ess_q10_hist <- c(M_local)
  bank_size_mean_hist <- c(mean(vapply(banks, function(b) nrow(b$alpha), numeric(1))))
  bank_topup_points <- 0L
  n_population_accept <- 0L

  while (lambda < 1 - 1e-12 && round < max_rounds) {
    round <- round + 1L
    increment_cache <- .prepare_population_delta_cache(
      banks = banks,
      phi = phi,
      gaussian_map_fn = gaussian_map_fn,
      lambda_ref = lambda
    )
    lambda_next <- choose_next_lambda_nested(
      phi = phi,
      banks = banks,
      w = w,
      lambda_curr = lambda,
      rho_step = rho_step,
      gaussian_map_fn = gaussian_map_fn,
      step_target_fn = step_target_fn,
      increment_cache = increment_cache
    )
    lambda_next <- max(lambda_next, lambda)
    vcat(sprintf("Round %d: lambda %.4f -> %.4f\n", round, lambda, lambda_next))

    log_w_inc <- .predict_population_log_increment_from_cache(
      cache = increment_cache,
      lambda_old = lambda,
      lambda_new = lambda_next
    )
    norm <- .normalize_logweights(log(pmax(w, .Machine$double.eps)) + log_w_inc)
    w <- norm$w
    log_evidence <- log_evidence + norm$log_norm
    lambda <- lambda_next
    lambda_hist <- c(lambda_hist, lambda)
    ess_now <- .ess_from_weights(w)
    ess_hist <- c(ess_hist, ess_now)

    bank_diag <- .bank_ess_summary(
      banks,
      phi = phi,
      lambda = lambda,
      gaussian_map_fn = gaussian_map_fn,
      increment_cache = increment_cache
    )
    bank_ess_min_hist <- c(bank_ess_min_hist, bank_diag$min)
    bank_ess_q10_hist <- c(bank_ess_q10_hist, unname(bank_diag$quantiles["q10"]))
    banks_mut <- lapply(increment_cache, `[[`, "bank")
    current_logZ_mat <- .logZ_matrix_from_cache(increment_cache, lambda = lambda)

    if (ess_now < rho_res * N) {
      idx <- .resample_indices(w)
      phi <- phi[idx, , drop = FALSE]
      w <- rep(1 / nrow(phi), nrow(phi))
      logprior <- logprior[idx]
      current_logZ_mat <- current_logZ_mat[idx, , drop = FALSE]
      ess_now <- .ess_from_weights(w)
      ess_hist[length(ess_hist)] <- ess_now
      vcat(sprintf("  population resample at ESS=%.1f\n", ess_now))
    }

    S <- tryCatch(weighted_cov(phi, w), error = function(e) stats::cov(phi))
    S <- as.matrix(Matrix::nearPD(S, conv.tol = 1e-7)$mat)
    diag(S) <- pmax(diag(S), 1e-8)
    rw_cov <- (0.25^2 / max(1L, ncol(phi))) * S + diag(1e-10, ncol(phi))

    for (mv in seq_len(as.integer(n_population_moves))) {
      for (n in seq_len(N)) {
        out <- mutate_population_particle_nested(
          phi_row = phi[n, ],
          logprior_row = logprior[n],
          lambda = lambda,
          banks = banks_mut,
          logprior_phi = logprior_phi,
          gaussian_map_fn = gaussian_map_fn,
          rw_cov = rw_cov,
          current_logZ = current_logZ_mat[n, ],
          seed = as.integer(seed + 500000L * round + 1009L * mv + n)
        )
        phi[n, ] <- out$phi
        logprior[n] <- out$logprior
        current_logZ_mat[n, ] <- out$logZ
        n_population_accept <- n_population_accept + as.integer(out$accepted)
      }
    }

    round_topups <- 0L
    for (i in seq_along(banks)) {
      up <- enrich_local_bank(
        banks[[i]],
        phi = phi,
        lambda = lambda,
        gaussian_map_fn = gaussian_map_fn,
        rho_local = rho_local,
        split_tol = bank_split_tol,
        max_topups = max_bank_topups,
        topup_batch = bank_topup_batch,
        seed = as.integer(seed + 900000L * round + i)
      )
      banks[[i]] <- up$bank
      round_topups <- round_topups + up$n_new
    }
    bank_topup_points <- bank_topup_points + round_topups
    bank_diag_post <- .bank_ess_summary(banks, phi = phi, lambda = lambda, gaussian_map_fn = gaussian_map_fn)
    bank_ess_min_hist[length(bank_ess_min_hist)] <- bank_diag_post$min
    bank_ess_q10_hist[length(bank_ess_q10_hist)] <- unname(bank_diag_post$quantiles["q10"])
    bank_size_mean_hist <- c(bank_size_mean_hist, mean(vapply(banks, function(b) nrow(b$alpha), numeric(1))))
    if (round_topups > 0L) {
      vcat(sprintf("  bank top-ups added %d cached local points\n", round_topups))
    }
    vcat(sprintf("  ESS=%.1f | logZ≈%.4f\n", ess_now, log_evidence))

    if (lambda >= 1 - 1e-12) break
  }

  list(
    phi = phi,
    w = w,
    local_banks = banks,
    log_evidence = log_evidence,
    meta = list(
      rounds = round,
      lambda_hist = lambda_hist,
      ess_hist = ess_hist,
      bank_ess_min_hist = bank_ess_min_hist,
      bank_ess_q10_hist = bank_ess_q10_hist,
      bank_size_mean_hist = bank_size_mean_hist,
      bank_topup_points = as.integer(bank_topup_points),
      n_population_accept = as.integer(n_population_accept),
      implementation = "shared_local_banks_baseline"
    )
  )
}

source("hierarchical_collapsed_smc.R")

hierarchical_smc <- function(method = c("shared_bank", "collapsed_subject", "collapsed_bridge"), ...) {
  method <- match.arg(method)
  if (identical(method, "shared_bank")) {
    hierarchical_smc_from_local_objects(...)
  } else if (identical(method, "collapsed_subject")) {
    collapsed_subject_smc(...)
  } else {
    collapsed_bridge_smc(...)
  }
}

hierarchical_smc_from_local_objects <- function(local_objs,
                                                prior,
                                                gaussian_map_fn,
                                                N = 2000L,
                                                M_local = 64L,
                                                rho_step = NULL,
                                                rho_res = 0.5,
                                                rho_local = 0.5,
                                                n_population_moves = 1L,
                                                max_rounds = 100L,
                                                max_bank_topups = 2L,
                                                bank_split_tol = 0.03,
                                                bank_topup_batch = NULL,
                                                seed = 123L,
                                                verbose = TRUE) {
  nested_population_smc(
    rprior_phi = prior$rprior,
    local_objs = local_objs,
    logprior_phi = prior$lprior,
    gaussian_map_fn = gaussian_map_fn,
    N = as.integer(N),
    M_local = as.integer(M_local),
    rho_step = rho_step,
    rho_res = rho_res,
    rho_local = rho_local,
    n_population_moves = as.integer(n_population_moves),
    max_rounds = as.integer(max_rounds),
    max_bank_topups = as.integer(max_bank_topups),
    bank_split_tol = bank_split_tol,
    bank_topup_batch = bank_topup_batch,
    seed = as.integer(seed),
    verbose = verbose
  )
}

summarize_phi_diag_posterior <- function(phi, w, param_names = NULL) {
  phi <- as.matrix(phi)
  d <- ncol(phi) / 2L
  if (is.null(param_names)) {
    param_names <- paste0("theta", seq_len(d))
  }
  w <- pmax(as.numeric(w), 0)
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) {
    w <- rep(1 / nrow(phi), nrow(phi))
  } else {
    w <- w / sw
  }

  weighted_quantile <- function(x, w, probs = c(0.025, 0.5, 0.975)) {
    ord <- order(x)
    x_ord <- x[ord]
    w_ord <- w[ord]
    w_ord <- w_ord / sum(w_ord)
    cw <- cumsum(w_ord)
    as.numeric(stats::approx(cw, x_ord, xout = probs, rule = 2)$y)
  }

  mu_part <- phi[, seq_len(d), drop = FALSE]
  sigma2_part <- exp(phi[, d + seq_len(d), drop = FALSE])
  colnames(mu_part) <- param_names
  colnames(sigma2_part) <- param_names

  mu_q <- t(vapply(seq_len(d), function(j) weighted_quantile(mu_part[, j], w), numeric(3)))
  sigma2_q <- t(vapply(seq_len(d), function(j) weighted_quantile(sigma2_part[, j], w), numeric(3)))

  list(
    mu = data.frame(
      parameter = param_names,
      mean = colSums(mu_part * w),
      median = mu_q[, 2L],
      q025 = mu_q[, 1L],
      q975 = mu_q[, 3L],
      row.names = NULL,
      check.names = FALSE
    ),
    sigma2 = data.frame(
      parameter = param_names,
      mean = colSums(sigma2_part * w),
      median = sigma2_q[, 2L],
      q025 = sigma2_q[, 1L],
      q975 = sigma2_q[, 3L],
      row.names = NULL,
      check.names = FALSE
    )
  )
}
