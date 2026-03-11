# ========================================================================
# Maintained local / population helpers
# - local surrogate objects built from per-group SMC output
# - Gaussian working-prior helpers
# - exact init / refresh proposal builders for population-conditioned local work
# ========================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(mvtnorm)
})

sc_path <- file.path(getwd(), "smc_core.R")
if (file.exists(sc_path)) source(sc_path)

.ensure_mix_cache <- function(mix) {
  if (is.null(mix$cache)) {
    mix$cache <- prep_mix_cache(mix$meansZ, mix$covsZ, mix$wZ)
  }
  mix
}

.regularize_cov_safe <- function(S, jitter = 1e-8) {
  S <- as.matrix((S + t(S)) / 2)
  S <- tryCatch(as.matrix(Matrix::nearPD(S, conv.tol = 1e-7)$mat), error = function(e) S)
  d <- nrow(S)
  S + diag(jitter, d)
}

.weighted_mean_vec <- function(X, w) {
  w <- pmax(as.numeric(w), 0)
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) {
    w <- rep(1 / nrow(as.matrix(X)), nrow(as.matrix(X)))
  } else {
    w <- w / sw
  }
  colSums(as.matrix(X) * w)
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

.log_prior_gauss_mat <- function(Theta, pars) {
  Theta <- as.matrix(Theta)
  mu <- as.numeric(pars$mu)
  K <- as.matrix(pars$Sigma_inv)
  lt <- sweep(Theta, 2L, mu, `-`)
  q <- rowSums((lt %*% K) * lt)
  out <- -0.5 * (q + pars$logdet + ncol(Theta) * log(2 * pi))
  if (!is.null(pars$const)) out <- out + pars$const
  as.numeric(out)
}

.draw_local_theta_given_population_callbacks <- function(phi, n, callbacks, d_theta) {
  if (!is.null(callbacks$rtheta_given_phi)) {
    th <- tryCatch(
      callbacks$rtheta_given_phi(phi, n, aux = NULL),
      error = function(e) callbacks$rtheta_given_phi(phi, n)
    )
    th <- as.matrix(th)
    if (nrow(th) != n) stop("rtheta_given_phi returned incorrect number of draws.")
    return(th)
  }
  if (is.null(callbacks$gaussian_map_fn)) {
    stop("Need either callbacks$rtheta_given_phi or callbacks$gaussian_map_fn.")
  }
  pars <- callbacks$gaussian_map_fn(phi, d_theta)
  K <- as.matrix(pars$Sigma_inv)
  L <- tryCatch(chol(K), error = function(e) NULL)
  if (is.null(L)) stop("gaussian_map_fn returned non-SPD precision.")
  Sig <- chol2inv(L)
  Sig <- .regularize_cov_safe(Sig)
  mvtnorm::rmvnorm(as.integer(n), mean = as.numeric(pars$mu), sigma = Sig)
}

build_local_surrogate <- function(
    smc_out,
    data = NULL,
    subj_id = NA_integer_,
    loglik_fn = NULL,
    base_seed = NULL
) {
  stopifnot(!is.null(smc_out$Theta), !is.null(smc_out$Z), !is.null(smc_out$w), !is.null(smc_out$transport))

  d_theta <- ncol(smc_out$Theta)
  w <- pmax(as.numeric(smc_out$w), 0)
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) {
    w <- rep(1 / nrow(smc_out$Theta), nrow(smc_out$Theta))
  } else {
    w <- w / sw
  }

  muZ <- .weighted_mean_vec(smc_out$Z, w)
  SigZ <- .weighted_cov_mat(smc_out$Z, w)
  mix_local <- list(meansZ = list(as.numeric(muZ)), covsZ = list(SigZ), wZ = 1)
  mix_local$cache <- prep_mix_cache(mix_local$meansZ, mix_local$covsZ, mix_local$wZ)

  theta_mean <- as.numeric(.weighted_mean_vec(smc_out$Theta, w))
  theta_cov <- .weighted_cov_mat(smc_out$Theta, w)

  structure(
    list(
      subj_id = as.integer(subj_id),
      transport = smc_out$transport,
      mix_local = mix_local,
      scale_local = theta_cov,
      theta_mean = theta_mean,
      d_theta = d_theta,
      base_seed = as.integer(base_seed %||% (10007L + 1009L * as.integer(subj_id))),
      local_fit_diag = list(
        rounds = smc_out$meta$rounds %||% length(smc_out$meta$lambda_hist %||% numeric(0)),
        log_evidence = smc_out$log_evidence %||% NA_real_,
        mcse_log_evidence = smc_out$mcse_logZ %||% NA_real_,
        final_lambda = smc_out$final_lambda %||% 1.0
      ),
      data = data,
      loglik_fn = loglik_fn
    ),
    class = "local_surrogate"
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

.draw_from_working_prior_gaussian <- function(n, working_prior) {
  wp <- .normalize_working_prior_gaussian(working_prior)
  mvtnorm::rmvnorm(as.integer(n), mean = wp$mu, sigma = wp$Sigma)
}

log_working_prior_gaussian_mat <- function(Theta, working_prior) {
  wp <- .normalize_working_prior_gaussian(working_prior, d_theta = ncol(as.matrix(Theta)))
  .log_prior_gauss_mat(Theta, wp)
}

build_local_exact_object <- function(
    smc_out,
    data = NULL,
    subj_id = NA_integer_,
    loglik_fn = NULL,
    working_prior = NULL,
    base_seed = NULL,
    ...
) {
  surrogate <- build_local_surrogate(
    smc_out = smc_out,
    data = data,
    subj_id = subj_id,
    loglik_fn = loglik_fn,
    base_seed = base_seed
  )

  working_prior_use <- working_prior %||% smc_out$working_prior
  if (is.null(working_prior_use)) {
    stop("build_local_exact_object requires a working_prior or smc_out$working_prior.")
  }
  wp <- .normalize_working_prior_gaussian(working_prior_use, d_theta = surrogate$d_theta)

  structure(
    list(
      subj_id = as.integer(subj_id),
      surrogate = surrogate,
      working_prior = wp,
      data = data %||% surrogate$data,
      loglik_fn = loglik_fn %||% surrogate$loglik_fn,
      base_seed = as.integer(base_seed %||% surrogate$base_seed),
      d_theta = surrogate$d_theta,
      theta_bank = as.matrix(smc_out$Theta),
      theta_bank_loglik = if (!is.null(smc_out$loglik)) as.numeric(smc_out$loglik) else NULL,
      theta_bank_w = {
        w <- pmax(as.numeric(smc_out$w), 0)
        sw <- sum(w)
        if (!is.finite(sw) || sw <= 0) rep(1 / nrow(as.matrix(smc_out$Theta)), nrow(as.matrix(smc_out$Theta))) else w / sw
      }
    ),
    class = "local_exact_object"
  )
}

make_local_theta_proposal_init <- function(local_obj, init_ctl = list()) {
  stopifnot(inherits(local_obj, "local_exact_object"))
  ctl <- modifyList(
    list(
      w_local = 0.90,
      w_phi_anchor = 0.00,
      w_defensive = 0.10,
      defensive_df = 3L,
      anchor_n = 32L,
      anchor_cov_inflation = 1.0
    ),
    init_ctl
  )
  callbacks <- list(
    rtheta_given_phi = function(phi, n, aux = NULL) {
      .draw_from_working_prior_gaussian(n, local_obj$working_prior)
    }
  )
  phi_stub <- local_obj$working_prior$phi_anchor %||% local_obj$working_prior$mu
  make_local_theta_proposal(
    phi = phi_stub,
    surrogate = local_obj$surrogate,
    proposal_control = ctl,
    callbacks = callbacks
  )
}

make_local_theta_proposal_refresh <- function(phi, local_obj, refresh_ctl = list(), callbacks = list()) {
  stopifnot(inherits(local_obj, "local_exact_object"))
  ctl <- modifyList(
    list(
      w_local = 0.60,
      w_phi_anchor = 0.30,
      w_defensive = 0.10,
      defensive_df = 3L,
      anchor_n = 96L,
      anchor_cov_inflation = 1.5
    ),
    refresh_ctl
  )
  cb <- modifyList(
    list(
      rtheta_given_phi = function(phi, n, aux = NULL) {
        .draw_from_working_prior_gaussian(n, local_obj$working_prior)
      }
    ),
    callbacks
  )
  make_local_theta_proposal(
    phi = phi,
    surrogate = local_obj$surrogate,
    proposal_control = ctl,
    callbacks = cb
  )
}

make_local_theta_proposal <- function(phi, surrogate, proposal_control = list(), callbacks = list()) {
  stopifnot(inherits(surrogate, "local_surrogate"))
  ctl <- modifyList(
    list(
      w_local = 0.60,
      w_phi_anchor = 0.30,
      w_defensive = 0.10,
      defensive_df = 3L,
      anchor_n = 96L,
      anchor_cov_inflation = 1.5
    ),
    proposal_control
  )

  w_comp <- c(max(0, ctl$w_local), max(0, ctl$w_phi_anchor), max(0, ctl$w_defensive))
  if (sum(w_comp) <= 0) w_comp <- c(0.6, 0.3, 0.1)
  w_comp <- w_comp / sum(w_comp)

  d <- surrogate$d_theta
  Tmap <- surrogate$transport
  mix_local <- .ensure_mix_cache(surrogate$mix_local)

  mix_anchor <- NULL
  Theta_anchor <- tryCatch(
    .draw_local_theta_given_population_callbacks(phi, as.integer(max(16L, ctl$anchor_n)), callbacks, d),
    error = function(e) NULL
  )
  if (!is.null(Theta_anchor) && nrow(Theta_anchor) > 1L) {
    Z_anchor <- Tmap$fwd(Theta_anchor)
    mu_a <- colMeans(Z_anchor)
    Zc <- sweep(Z_anchor, 2L, mu_a, `-`)
    Sig_a <- (t(Zc) %*% Zc) / max(1L, nrow(Zc))
    Sig_a <- .regularize_cov_safe(Sig_a * (as.numeric(ctl$anchor_cov_inflation)^2))
    mix_anchor <- list(meansZ = list(as.numeric(mu_a)), covsZ = list(Sig_a), wZ = 1)
    mix_anchor$cache <- prep_mix_cache(mix_anchor$meansZ, mix_anchor$covsZ, mix_anchor$wZ)
  }
  if (is.null(mix_anchor)) mix_anchor <- mix_local
  mix_anchor <- .ensure_mix_cache(mix_anchor)
  mix_def <- mix_anchor

  draw <- function(n, seed = NULL) {
    n <- as.integer(n)
    if (n <= 0L) return(matrix(numeric(0), 0, d))
    if (!is.null(seed)) set.seed(as.integer(seed))
    cnt <- as.integer(rmultinom(1L, size = n, prob = w_comp)[, 1L])

    Z_list <- list()
    if (cnt[1L] > 0L) {
      Z_list[[length(Z_list) + 1L]] <- sample_gmm_Z_qmc(
        cnt[1L], mix_local$meansZ, mix_local$cache,
        seed = if (is.null(seed)) NULL else seed + 17L
      )
    }
    if (cnt[2L] > 0L) {
      Z_list[[length(Z_list) + 1L]] <- sample_gmm_Z_qmc(
        cnt[2L], mix_anchor$meansZ, mix_anchor$cache,
        seed = if (is.null(seed)) NULL else seed + 31L
      )
    }
    if (cnt[3L] > 0L) {
      Z_list[[length(Z_list) + 1L]] <- rmvt_mixture_Z_qmc(
        cnt[3L], mix_def$meansZ, mix_def$cache,
        nu = as.integer(max(2L, ctl$defensive_df)),
        seed = if (is.null(seed)) NULL else seed + 47L
      )
    }

    Z <- do.call(rbind, Z_list)
    if (nrow(Z) != n) stop("Adaptive proposal generated wrong number of draws.")
    if (nrow(Z) > 1L) Z <- Z[sample.int(nrow(Z)), , drop = FALSE]
    as.matrix(Tmap$inv(Z))
  }

  log_q <- function(Theta) {
    Theta <- as.matrix(Theta)
    Z <- Tmap$fwd(Theta)
    lq <- rep(-Inf, nrow(Z))
    if (w_comp[1L] > 0) {
      l1 <- log(w_comp[1L]) + gmm_logpdf_Z_vec(Z, mix_local$meansZ, mix_local$cache)
      lq <- l1
    }
    if (w_comp[2L] > 0) {
      l2 <- log(w_comp[2L]) + gmm_logpdf_Z_vec(Z, mix_anchor$meansZ, mix_anchor$cache)
      lq <- if (all(is.infinite(lq))) l2 else rlogsumexp2(lq, l2)
    }
    if (w_comp[3L] > 0) {
      l3 <- log(w_comp[3L]) + dmvt_mixture_logpdf_Z_vec(
        Z, mix_def$meansZ, mix_def$cache,
        nu = as.integer(max(2L, ctl$defensive_df))
      )
      lq <- if (all(is.infinite(lq))) l3 else rlogsumexp2(lq, l3)
    }
    as.numeric(lq + Tmap$log_jac(Theta))
  }

  list(
    draw = draw,
    log_q = log_q,
    w_comp = w_comp,
    mix_local = mix_local,
    mix_anchor = mix_anchor,
    mix_def = mix_def
  )
}

# Legacy benchmark-facing names retained as thin aliases.
build_subject_surrogate_from_outer <- build_local_exact_object
build_subject_exact_local_from_outer <- build_local_exact_object
