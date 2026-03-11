# ========================================================================
# Helpers for the explicit population/local state used by the future nested plan
#
# These are not the maintained stage-2 driver. They are just building blocks:
# - explicit local-state bookkeeping
# - exact initialization from local proposal objects
# - exact population and local refresh kernels for an extended-state method
# ========================================================================

source("new_SMC_cache.R")

sc_path <- file.path(getwd(), "smc_core.R")
if (file.exists(sc_path)) source(sc_path)

suppressPackageStartupMessages({
  library(Matrix)
})

materialize_local_state_from_indices <- function(alpha_bank, index_mat) {
  alpha_bank <- as.matrix(alpha_bank)
  index_mat <- as.matrix(index_mat)
  N <- nrow(index_mat)
  S <- ncol(index_mat)
  d <- ncol(alpha_bank)
  out <- array(NA_real_, dim = c(N, S, d))
  for (i in seq_len(S)) {
    out[, i, ] <- alpha_bank[index_mat[, i], , drop = FALSE]
  }
  out
}

refresh_local_indices_exact <- function(index_mat, n_bank, frac = 0.10, seed = NULL) {
  if (!is.null(seed)) set.seed(as.integer(seed))
  index_mat <- as.matrix(index_mat)
  frac <- .clamp(as.numeric(frac), 0, 1)
  if (frac <= 0) return(index_mat)
  N <- nrow(index_mat)
  S <- ncol(index_mat)
  R <- max(1L, ceiling(frac * S))
  for (n in seq_len(N)) {
    idx <- sample.int(S, R)
    index_mat[n, idx] <- sample.int(as.integer(n_bank), length(idx), replace = TRUE)
  }
  index_mat
}

compute_local_suff_stats <- function(alpha_state, X_list = NULL) {
  alpha_state <- array(alpha_state, dim = dim(alpha_state))
  d <- dim(alpha_state)[3L]
  if (!is.null(X_list)) {
    stop("compute_local_suff_stats only supports intercept-only form in the current skeleton.")
  }
  sum_alpha <- apply(alpha_state, c(1L, 3L), sum)
  if (!is.matrix(sum_alpha)) sum_alpha <- matrix(sum_alpha, ncol = d)
  N <- dim(alpha_state)[1L]
  sum_cross <- array(0, dim = c(N, d, d))
  for (n in seq_len(N)) {
    A <- alpha_state[n, , , drop = FALSE]
    A <- matrix(A, ncol = d)
    sum_cross[n, , ] <- crossprod(A)
  }
  list(sum_alpha = sum_alpha, sum_cross = sum_cross, d = d, n_particles = N)
}

update_local_suff_stats_one <- function(stats_obj, alpha_old, alpha_new, x_i = NULL) {
  if (!is.null(x_i)) {
    stop("update_local_suff_stats_one only supports intercept-only form in the current skeleton.")
  }
  alpha_old <- as.numeric(alpha_old)
  alpha_new <- as.numeric(alpha_new)
  if (!all(c("sum_alpha", "sum_cross") %in% names(stats_obj))) {
    stop("stats_obj must contain sum_alpha and sum_cross.")
  }
  stats_obj$sum_alpha <- stats_obj$sum_alpha + (alpha_new - alpha_old)
  stats_obj$sum_cross <- stats_obj$sum_cross +
    tcrossprod(alpha_new) - tcrossprod(alpha_old)
  stats_obj
}

log_tempered_population_contrib <- function(phi_row, alpha_mat, beta, gaussian_map_fn, working_priors) {
  alpha_mat <- as.matrix(alpha_mat)
  S <- nrow(alpha_mat)
  if (length(working_priors) != S) {
    stop("length(working_priors) must match nrow(alpha_mat).")
  }
  pars <- gaussian_map_fn(phi_row, ncol(alpha_mat))
  lp_phi <- .log_prior_gauss_mat(alpha_mat, pars)
  lp_g <- numeric(S)
  for (i in seq_len(S)) {
    lp_g[i] <- log_working_prior_gaussian_mat(alpha_mat[i, , drop = FALSE], working_priors[[i]])
  }
  sum((1 - beta) * lp_g + beta * lp_phi)
}

.log_local_given_population_sum_mat <- function(Phi, alpha_state, gaussian_map_fn) {
  Phi <- as.matrix(Phi)
  N <- nrow(Phi)
  d_theta <- dim(alpha_state)[3L]
  out <- numeric(N)
  for (n in seq_len(N)) {
    pars <- gaussian_map_fn(Phi[n, , drop = TRUE], d_theta)
    A <- matrix(alpha_state[n, , , drop = FALSE], ncol = d_theta)
    out[n] <- sum(.log_prior_gauss_mat(A, pars))
  }
  out
}

init_population_particles_exact <- function(local_objs, phi_particles, phi_anchor = NULL, init_ctl = list()) {
  phi_particles <- as.matrix(phi_particles)
  N <- nrow(phi_particles)
  S <- length(local_objs)
  if (S <= 0L) stop("local_objs must not be empty.")
  d_theta <- local_objs[[1L]]$d_theta
  ctl <- modifyList(
    list(
      seed = NULL,
      mode = "surrogate_is",
      proposal_control = list()
    ),
    init_ctl
  )
  if (!is.null(ctl$seed)) set.seed(as.integer(ctl$seed))
  init_mode <- match.arg(as.character(ctl$mode), c("surrogate_is", "bank_resample"))
  use_init_is_weights <- identical(init_mode, "surrogate_is")

  alpha_state <- array(NA_real_, dim = c(N, S, d_theta))
  subject_loglik <- matrix(NA_real_, nrow = N, ncol = S)
  subject_logg <- matrix(NA_real_, nrow = N, ncol = S)
  subject_logq <- matrix(NA_real_, nrow = N, ncol = S)

  for (j in seq_len(S)) {
    loc <- local_objs[[j]]
    seed_j <- as.integer((ctl$seed %||% loc$base_seed) + 1049L * j)
    prop <- NULL
    if (identical(init_mode, "surrogate_is")) {
      prop <- make_local_theta_proposal_init(loc, init_ctl = ctl$proposal_control)
      Theta_j <- prop$draw(N, seed = seed_j)
      Theta_j <- as.matrix(Theta_j)
      if (nrow(Theta_j) != N) stop("Initialization proposal returned wrong number of draws.")
    } else {
      set.seed(seed_j)
      idx <- sample.int(nrow(loc$theta_bank), size = N, replace = TRUE, prob = loc$theta_bank_w)
      Theta_j <- as.matrix(loc$theta_bank[idx, , drop = FALSE])
    }
    alpha_state[, j, ] <- Theta_j
    subject_loglik[, j] <- ll_parallel(Theta_j, loc$data, loc$loglik_fn, n_cores = 1L)
    subject_logg[, j] <- log_working_prior_gaussian_mat(Theta_j, loc$working_prior)
    subject_logq[, j] <- if (identical(init_mode, "surrogate_is")) prop$log_q(Theta_j) else 0
  }

  if (use_init_is_weights) {
    logw_init <- rowSums(subject_loglik + subject_logg - subject_logq)
    lw_shift <- logw_init - max(logw_init)
    w <- exp(lw_shift)
    w <- w / sum(w)
  } else {
    logw_init <- rep(0, N)
    w <- rep(1 / N, N)
  }

  list(
    alpha_state = alpha_state,
    subject_loglik = subject_loglik,
    subject_logg = subject_logg,
    subject_logq = subject_logq,
    log_weight = logw_init,
    w = w,
    init_mode = init_mode,
    ess_init = as.numeric(1 / sum(w * w)),
    ess_init_frac = as.numeric((1 / sum(w * w)) / max(N, 1L)),
    max_weight_init = as.numeric(max(w)),
    suff_stats = compute_local_suff_stats(alpha_state),
    n_exact_loglik_init = as.integer(N * S),
    phi_anchor = phi_anchor
  )
}

rejuvenate_population_given_local_exact <- function(phi, alpha_state, logprior, w, lambda,
                                             logprior_phi, gaussian_map_fn,
                                             rw_scale = 1.0,
                                             n_moves = 2L,
                                             log_alpha_given_phi_sum = NULL,
                                             rng_seed = NULL) {
  if (!is.null(rng_seed)) set.seed(as.integer(rng_seed))
  phi <- as.matrix(phi)
  N <- nrow(phi)
  d <- ncol(phi)
  if (is.null(log_alpha_given_phi_sum)) {
    log_alpha_given_phi_sum <- .log_local_given_population_sum_mat(phi, alpha_state, gaussian_map_fn)
  }

  S <- tryCatch(weighted_cov(phi, w), error = function(e) stats::cov(phi))
  S <- as.matrix(Matrix::nearPD(S, conv.tol = 1e-7)$mat)
  diag(S) <- pmax(diag(S), 1e-8)
  L <- tryCatch(
    chol((rw_scale^2 / max(d, 1L)) * S + diag(1e-10, d)),
    error = function(e) chol((rw_scale^2 / max(d, 1L)) * (S + diag(1e-6, d)))
  )

  accepted <- 0L
  for (mv in seq_len(n_moves)) {
    prop <- phi + matrix(rnorm(N * d), N, d) %*% L
    lp_prop <- apply(prop, 1L, logprior_phi)
    loga_prop <- .log_local_given_population_sum_mat(prop, alpha_state, gaussian_map_fn)
    loga <- (lp_prop - logprior) + lambda * (loga_prop - log_alpha_given_phi_sum)
    loga[!is.finite(loga)] <- -Inf
    acc <- log(runif(N)) < pmin(0, loga)
    if (any(acc)) {
      phi[acc, ] <- prop[acc, , drop = FALSE]
      logprior[acc] <- lp_prop[acc]
      log_alpha_given_phi_sum[acc] <- loga_prop[acc]
      accepted <- accepted + sum(acc)
    }
  }

  list(
    phi = phi,
    logprior = logprior,
    log_alpha_given_phi_sum = log_alpha_given_phi_sum,
    acc_rate = accepted / max(1L, N * n_moves)
  )
}

refresh_local_block_exact <- function(phi, alpha_state, local_objs, beta,
                                      refresh_frac = 0.05,
                                      delayed_accept = TRUE,
                                      refresh_ctl = list(),
                                      subject_loglik = NULL,
                                      subject_logg = NULL,
                                      gaussian_map_fn = NULL,
                                      seed = NULL) {
  if (!is.null(seed)) set.seed(as.integer(seed))
  phi <- as.matrix(phi)
  N <- nrow(phi)
  S <- length(local_objs)
  d_theta <- dim(alpha_state)[3L]
  R <- max(1L, ceiling(.clamp(as.numeric(refresh_frac), 0, 1) * S))

  if (is.null(subject_loglik)) subject_loglik <- matrix(NA_real_, nrow = N, ncol = S)
  if (is.null(subject_logg)) subject_logg <- matrix(NA_real_, nrow = N, ncol = S)

  n_screen_pass <- 0L
  n_final_acc <- 0L
  n_exact_loglik <- 0L

  for (n in seq_len(N)) {
    subj_idx <- sample.int(S, R)
    pars_n <- if (!is.null(gaussian_map_fn)) gaussian_map_fn(phi[n, , drop = TRUE], d_theta) else NULL
    for (j in subj_idx) {
      loc <- local_objs[[j]]
      curr <- matrix(alpha_state[n, j, ], nrow = 1L)
      if (!is.finite(subject_loglik[n, j])) {
        subject_loglik[n, j] <- ll_parallel(curr, loc$data, loc$loglik_fn, n_cores = 1L)
      }
      if (!is.finite(subject_logg[n, j])) {
        subject_logg[n, j] <- log_working_prior_gaussian_mat(curr, loc$working_prior)
      }

      prop <- make_local_theta_proposal_refresh(
        phi = phi[n, , drop = TRUE],
        local_obj = loc,
        refresh_ctl = refresh_ctl,
        callbacks = list(gaussian_map_fn = gaussian_map_fn)
      )
      cand <- prop$draw(1L, seed = as.integer((seed %||% loc$base_seed) + 1009L * n + 31L * j))
      cand <- as.matrix(cand)
      q_curr <- prop$log_q(curr)
      q_cand <- prop$log_q(cand)
      logg_curr <- subject_logg[n, j]
      logg_cand <- log_working_prior_gaussian_mat(cand, loc$working_prior)

      if (is.null(pars_n)) {
        stop("refresh_local_block_exact requires gaussian_map_fn in the current skeleton.")
      }
      logp_curr <- .log_prior_gauss_mat(curr, pars_n)
      logp_cand <- .log_prior_gauss_mat(cand, pars_n)

      screen_loga <- (1 - beta) * (logg_cand - logg_curr) +
        beta * (logp_cand - logp_curr) +
        (q_curr - q_cand)

      if (!delayed_accept) {
        loglik_cand <- ll_parallel(cand, loc$data, loc$loglik_fn, n_cores = 1L)
        n_exact_loglik <- n_exact_loglik + 1L
        full_loga <- screen_loga + (loglik_cand - subject_loglik[n, j])
        if (is.finite(full_loga) && log(runif(1L)) < min(0, full_loga)) {
          alpha_state[n, j, ] <- cand[1L, ]
          subject_loglik[n, j] <- loglik_cand
          subject_logg[n, j] <- logg_cand
          n_final_acc <- n_final_acc + 1L
        }
        next
      }

      if (!is.finite(screen_loga) || log(runif(1L)) >= min(0, screen_loga)) next
      n_screen_pass <- n_screen_pass + 1L

      loglik_cand <- ll_parallel(cand, loc$data, loc$loglik_fn, n_cores = 1L)
      n_exact_loglik <- n_exact_loglik + 1L
      final_loga <- loglik_cand - subject_loglik[n, j]
      if (is.finite(final_loga) && log(runif(1L)) < min(0, final_loga)) {
        alpha_state[n, j, ] <- cand[1L, ]
        subject_loglik[n, j] <- loglik_cand
        subject_logg[n, j] <- logg_cand
        n_final_acc <- n_final_acc + 1L
      }
    }
  }

  list(
    alpha_state = alpha_state,
    subject_loglik = subject_loglik,
    subject_logg = subject_logg,
    n_refresh = as.integer(N * R),
    n_screen_pass = as.integer(n_screen_pass),
    n_final_accept = as.integer(n_final_acc),
    n_exact_loglik = as.integer(n_exact_loglik)
  )
}
