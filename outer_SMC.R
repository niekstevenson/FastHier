# ========================================================================
# Outer SMC over φ with Option-B inner caches (K RQMC batches)
# - Pseudo-marginal exactness via batch-index auxiliary variables
# - Correlated PM: keep b_idx fixed within moves, refresh sparsely
# - Adaptive tempering via rCESS; resampling (systematic | multinomial)
# - Stable numerics; finite checks; reproducible seeding
# - Diagnostics: per-subject ESS_i/M & PSIS k from active batches
# - Records per-round log-evidence increments
# ========================================================================

source("new_SMC_cache.R")

# Use shared core utilities (numerics, resampling, rCESS, etc.)
sc_path <- file.path(getwd(), "smc_core.R")
if (file.exists(sc_path)) source(sc_path)

suppressPackageStartupMessages({
  library(Matrix)
  library(matrixStats)
})

# ------------------------------ numerics ---------------------------------
# (moved to smc_core: .clamp, logsumexp, logsumexp_w, `%||%`)

# ------------------------------ rCESS step -------------------------------
# (moved to smc_core: rCESS_stat, next_lambda_via_rCESS_stat)

# ------------------------------ resampling -------------------------------
.resample_multinomial <- function(w) {
  # unbiased for evidence; more variance
  N <- length(w)
  sample.int(N, size = N, replace = TRUE, prob = w)
}

# Use smc_core::stratified_resample_sorted for low-variance systematic resampling

# ---------------------------- φ utilities --------------------------------

.choose_probe_phi <- function(phi, w, logprior, loglik, lambda, mode = c("wmean","best")) {
  mode <- match.arg(mode)
  if (mode == "wmean") {
    ww <- pmax(w, 0); ww <- ww / sum(ww)
    drop(t(phi) %*% ww)
  } else {
    idx <- which.max(logprior + lambda * loglik)
    phi[idx, , drop = TRUE]
  }
}

# --------------------------- Inner log-likelihood ------------------------
# Batch-indexed pseudo-marginal log-lik: sum_i log \hat Z_i(φ | b_i)
.outer_loglik_hat_batch <- function(phi_row, caches, b_idx,
                                    log_prior_theta_given_phi_mat = NULL,
                                    gaussian_map_fn = NULL) {
  S <- length(caches)
  s <- 0.0
  if (!is.null(gaussian_map_fn)) {
    for (i in seq_len(S)) {
      s <- s + log_marginal_unbiased_gaussian_batch(caches[[i]], phi_row, b_idx[i])
    }
  } else {
    stopifnot(!is.null(log_prior_theta_given_phi_mat))
    for (i in seq_len(S)) {
      s <- s + log_marginal_unbiased_batch(caches[[i]], phi_row, b_idx[i],
                                           log_prior_theta_given_phi_mat)
    }
  }
  s
}

.outer_loglik_hat_particles <- function(Phi, caches, b_state,
                                        log_prior_theta_given_phi_mat = NULL,
                                        gaussian_map_fn = NULL,
                                        pm_mode = c("strict", "fast")) {
  pm_mode <- match.arg(pm_mode)
  Phi <- as.matrix(Phi)
  N <- nrow(Phi)
  if (N <= 0L) return(numeric(0))

  if (pm_mode == "fast") {
    stopifnot(is.vector(b_state), length(b_state) == length(caches))
    return(apply(Phi, 1L, function(p)
      .outer_loglik_hat_batch(p, caches, b_state,
                              log_prior_theta_given_phi_mat, gaussian_map_fn)))
  }

  stopifnot(is.matrix(b_state), nrow(b_state) == N, ncol(b_state) == length(caches))
  out <- numeric(N)
  for (i in seq_len(N)) {
    out[i] <- .outer_loglik_hat_batch(Phi[i, , drop = TRUE], caches, b_state[i, ],
                                      log_prior_theta_given_phi_mat, gaussian_map_fn)
  }
  out
}

.init_aux_batch_state <- function(caches, M, pm_mode = c("strict", "fast")) {
  pm_mode <- match.arg(pm_mode)
  S <- length(caches)
  if (pm_mode == "fast") {
    return(vapply(caches, function(ci) sample.int(ci$K, 1L), integer(1)))
  }
  B <- matrix(0L, nrow = M, ncol = S)
  for (j in seq_len(S)) {
    B[, j] <- sample.int(caches[[j]]$K, M, replace = TRUE)
  }
  B
}

.probe_aux_row <- function(phi, w, logprior, loglik, lambda, b_state) {
  idx <- which.max(logprior + lambda * loglik)
  if (is.matrix(b_state)) b_state[idx, ] else b_state
}

.refresh_aux_state_block_mh <- function(phi, loglik, b_state, lambda,
                                        caches,
                                        log_prior_theta_given_phi_mat = NULL,
                                        gaussian_map_fn = NULL,
                                        frac = 0.10,
                                        rng_seed = NULL) {
  if (!is.matrix(b_state)) {
    return(list(b_state = b_state, loglik = loglik, acc_rate = NA_real_))
  }
  if (!is.null(rng_seed)) set.seed(rng_seed)
  frac <- .clamp(frac, 0, 1)
  if (frac <= 0) return(list(b_state = b_state, loglik = loglik, acc_rate = 0))

  N <- nrow(phi)
  S <- ncol(b_state)
  R <- max(1L, ceiling(frac * S))
  nacc <- 0L
  for (i in seq_len(N)) {
    idx <- sample.int(S, R)
    b_prop <- b_state[i, ]
    for (j in idx) b_prop[j] <- sample.int(caches[[j]]$K, 1L)
    ll_prop <- .outer_loglik_hat_batch(phi[i, , drop = TRUE], caches, b_prop,
                                       log_prior_theta_given_phi_mat, gaussian_map_fn)
    loga <- lambda * (ll_prop - loglik[i])
    if (!is.finite(loga)) loga <- -Inf
    if (log(runif(1)) < min(0, loga)) {
      b_state[i, ] <- b_prop
      loglik[i] <- ll_prop
      nacc <- nacc + 1L
    }
  }
  list(b_state = b_state, loglik = loglik, acc_rate = nacc / max(N, 1L))
}

# ------------------------ Adaptive PM aux-state API ----------------------
init_pm_aux_state <- function(M_particles, S_subjects, seed = 123L, M_default = 128L,
                              pm_mode = c("strict", "fast")) {
  pm_mode <- match.arg(pm_mode)
  set.seed(as.integer(seed))
  if (length(M_default) == 1L) M_default <- rep(as.integer(M_default), S_subjects)
  M_default <- as.integer(M_default)
  if (length(M_default) != S_subjects) stop("M_default must be scalar or length S_subjects.")

  if (pm_mode == "strict") {
    list(
      seed = matrix(sample.int(.Machine$integer.max, M_particles * S_subjects, replace = TRUE),
                    nrow = M_particles, ncol = S_subjects),
      M_alloc = matrix(rep(M_default, each = M_particles), nrow = M_particles, ncol = S_subjects),
      mode = pm_mode
    )
  } else {
    list(
      seed = as.integer(sample.int(.Machine$integer.max, S_subjects, replace = TRUE)),
      M_alloc = as.integer(M_default),
      mode = pm_mode
    )
  }
}

refresh_pm_aux_state <- function(aux_state, frac = 0.05, seed = NULL) {
  frac <- .clamp(as.numeric(frac), 0, 1)
  if (frac <= 0) return(aux_state)
  if (!is.null(seed)) set.seed(as.integer(seed))

  if (is.matrix(aux_state$seed)) {
    n <- length(aux_state$seed)
    R <- as.integer(max(1L, ceiling(frac * n)))
    idx <- sample.int(n, R)
    aux_state$seed[idx] <- sample.int(.Machine$integer.max, R, replace = TRUE)
  } else {
    n <- length(aux_state$seed)
    R <- as.integer(max(1L, ceiling(frac * n)))
    idx <- sample.int(n, R)
    aux_state$seed[idx] <- sample.int(.Machine$integer.max, R, replace = TRUE)
  }
  aux_state
}

.pm_aux_subset <- function(aux_state, idx, pm_mode = c("strict", "fast")) {
  pm_mode <- match.arg(pm_mode)
  if (pm_mode == "strict") {
    list(
      seed = aux_state$seed[idx, , drop = FALSE],
      M_alloc = aux_state$M_alloc[idx, , drop = FALSE],
      mode = "strict"
    )
  } else {
    aux_state
  }
}

.pm_aux_probe_row <- function(phi, w, logprior, loglik, lambda, aux_state, pm_mode = c("strict", "fast")) {
  pm_mode <- match.arg(pm_mode)
  if (pm_mode == "fast") return(aux_state)
  idx <- which.max(logprior + lambda * loglik)
  list(
    seed = aux_state$seed[idx, ],
    M_alloc = aux_state$M_alloc[idx, ],
    mode = "fast"
  )
}

.pm_aux_apply_subject_M <- function(aux_state, subject_M, pm_mode = c("strict", "fast")) {
  pm_mode <- match.arg(pm_mode)
  subject_M <- as.integer(subject_M)
  if (pm_mode == "strict") {
    aux_state$M_alloc <- matrix(rep(subject_M, each = nrow(aux_state$M_alloc)),
                                nrow = nrow(aux_state$M_alloc), ncol = ncol(aux_state$M_alloc))
  } else {
    aux_state$M_alloc <- subject_M
  }
  aux_state
}

.adaptive_pm_callbacks <- function(log_prior_theta_given_phi_mat = NULL, gaussian_map_fn = NULL,
                                   rtheta_given_phi = NULL, data_list = NULL, loglik_fn = NULL) {
  list(
    log_prior_theta_given_phi_mat = log_prior_theta_given_phi_mat,
    gaussian_map_fn = gaussian_map_fn,
    rtheta_given_phi = rtheta_given_phi,
    data_list = data_list,
    loglik_fn = loglik_fn
  )
}

.outer_loglik_hat_particles_adaptive_pm <- function(Phi, surrogates, aux_state, callbacks,
                                                    control = list(),
                                                    pm_mode = c("strict", "fast"),
                                                    M_override = NULL) {
  pm_mode <- match.arg(pm_mode)
  Phi <- as.matrix(Phi)
  N <- nrow(Phi)
  S <- length(surrogates)
  out <- numeric(N)
  if (N <= 0L) return(out)

  for (i in seq_len(N)) {
    s <- 0.0
    for (j in seq_len(S)) {
      aux_ij <- if (pm_mode == "strict") {
        list(seed = aux_state$seed[i, j], M = aux_state$M_alloc[i, j])
      } else {
        list(seed = aux_state$seed[j], M = aux_state$M_alloc[j])
      }
      cb <- callbacks
      cb$data_i <- callbacks$data_list[[j]]
      est <- estimate_log_marginal_subject_pm(
        phi = Phi[i, , drop = TRUE],
        surrogate = surrogates[[j]],
        callbacks = cb,
        aux_state_i = aux_ij,
        control = modifyList(control, list(M_override = M_override))
      )
      s <- s + est$log_mhat
    }
    out[i] <- s
  }
  out
}

.adaptive_pm_probe_subject_diag <- function(phi_row, surrogates, aux_probe, callbacks, control = list()) {
  S <- length(surrogates)
  rows <- vector("list", S)
  for (j in seq_len(S)) {
    cb <- callbacks
    cb$data_i <- callbacks$data_list[[j]]
    aux_j <- list(seed = aux_probe$seed[j], M = aux_probe$M_alloc[j])
    est <- estimate_log_marginal_subject_pm(
      phi = phi_row,
      surrogate = surrogates[[j]],
      callbacks = cb,
      aux_state_i = aux_j,
      control = control
    )
    rows[[j]] <- data.frame(
      i = j,
      log_mhat = est$log_mhat,
      ess_norm = est$ess_norm,
      khat = est$khat,
      var_proxy = est$var_proxy,
      M_used = est$M_used,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

.adaptive_pm_update_subject_M <- function(subject_M, subject_diag, control = list()) {
  ctl <- modifyList(
    list(
      M_levels = c(64L, 128L, 256L, 512L),
      ess_norm_threshold = 0.15,
      khat_threshold = 0.80,
      var_log_target_low = 0.30,
      var_log_target_high = 1.00,
      var_relax_ess = 0.35,
      decrease_enabled = TRUE
    ),
    control
  )
  levels <- sort(unique(as.integer(ctl$M_levels)))
  cur <- as.integer(subject_M)
  out <- cur
  for (j in seq_along(cur)) {
    row <- subject_diag[subject_diag$i == j, , drop = FALSE]
    if (!nrow(row)) next
    ess_j <- as.numeric(row$ess_norm[1L])
    kh_j <- as.numeric(row$khat[1L])
    var_j <- as.numeric(row$var_proxy[1L])

    pos <- which(levels <= cur[j])
    pos <- if (length(pos)) max(pos) else 1L

    bad <- (!is.finite(var_j) || var_j > ctl$var_log_target_high) ||
      (is.finite(ess_j) && ess_j < ctl$ess_norm_threshold) ||
      (is.finite(kh_j) && kh_j > ctl$khat_threshold)
    if (bad) {
      jump <- if (!is.finite(var_j) || (is.finite(var_j) && var_j > 2 * ctl$var_log_target_high)) 2L else 1L
      pos_new <- min(length(levels), pos + jump)
      out[j] <- levels[pos_new]
      next
    }

    if (isTRUE(ctl$decrease_enabled)) {
      low_var <- is.finite(var_j) && var_j < ctl$var_log_target_low
      strong_ess <- is.finite(ess_j) && ess_j > ctl$var_relax_ess
      ok_k <- !is.finite(kh_j) || kh_j < 0.60
      if (low_var && strong_ess && ok_k) {
        pos_new <- max(1L, pos - 1L)
        out[j] <- levels[pos_new]
      }
    }
  }
  out
}

.rejuvenate_rw_adaptive_pm <- function(phi, loglik, logprior, w, lambda,
                                       surrogates, aux_state, logprior_phi, callbacks,
                                       control = list(),
                                       pm_mode = c("strict", "fast"),
                                       n_moves = 2L,
                                       rw_scale = 1.0,
                                       rng_seed = NULL) {
  pm_mode <- match.arg(pm_mode)
  if (!is.null(rng_seed)) set.seed(as.integer(rng_seed))
  N <- nrow(phi)
  d <- ncol(phi)

  ctl <- modifyList(
    list(
      da_enable = FALSE,
      da_M = 32L
    ),
    control
  )

  S <- weighted_cov(phi, w)
  S <- as.matrix(Matrix::nearPD(S, conv.tol = 1e-7)$mat)
  diag(S) <- pmax(diag(S), 1e-8)
  ev <- eigen(S, symmetric = TRUE)
  lamv <- pmax(ev$values, 1e-8)
  Sprop <- ev$vectors %*% diag(lamv, d) %*% t(ev$vectors)
  L <- tryCatch(chol(Sprop + diag(1e-8, d)), error = function(e) chol(Sprop + diag(1e-6, d)))
  step <- (rw_scale / sqrt(max(d, 1))) * L

  accepted <- 0L
  nprop <- 0L

  for (m in seq_len(n_moves)) {
    Z <- matrix(rnorm(N * d), N, d)
    prop <- phi + Z %*% t(step)
    lp_prop <- apply(prop, 1L, logprior_phi)

    if (isTRUE(ctl$da_enable)) {
      ll1_curr <- .outer_loglik_hat_particles_adaptive_pm(
        Phi = phi,
        surrogates = surrogates,
        aux_state = aux_state,
        callbacks = callbacks,
        control = control,
        pm_mode = pm_mode,
        M_override = as.integer(ctl$da_M)
      )
      ll1_prop <- .outer_loglik_hat_particles_adaptive_pm(
        Phi = prop,
        surrogates = surrogates,
        aux_state = aux_state,
        callbacks = callbacks,
        control = control,
        pm_mode = pm_mode,
        M_override = as.integer(ctl$da_M)
      )
      loga1 <- (lp_prop - logprior) + lambda * (ll1_prop - ll1_curr)
      u1 <- log(runif(N))
      pass1 <- (u1 < pmin(0, loga1))
      if (any(pass1)) {
        idx <- which(pass1)
        ll_prop_full <- .outer_loglik_hat_particles_adaptive_pm(
          Phi = prop[idx, , drop = FALSE],
          surrogates = surrogates,
          aux_state = .pm_aux_subset(aux_state, idx, pm_mode = pm_mode),
          callbacks = callbacks,
          control = control,
          pm_mode = pm_mode,
          M_override = NULL
        )
        loga2 <- lambda * ((ll_prop_full - loglik[idx]) - (ll1_prop[idx] - ll1_curr[idx]))
        u2 <- log(runif(length(idx)))
        acc_idx <- idx[which(u2 < pmin(0, loga2))]
        if (length(acc_idx)) {
          phi[acc_idx, ] <- prop[acc_idx, , drop = FALSE]
          logprior[acc_idx] <- lp_prop[acc_idx]
          loglik[acc_idx] <- ll_prop_full[match(acc_idx, idx)]
          accepted <- accepted + length(acc_idx)
        }
      }
      nprop <- nprop + N
    } else {
      ll_prop <- .outer_loglik_hat_particles_adaptive_pm(
        Phi = prop,
        surrogates = surrogates,
        aux_state = aux_state,
        callbacks = callbacks,
        control = control,
        pm_mode = pm_mode
      )
      loga <- (lp_prop - logprior) + lambda * (ll_prop - loglik)
      u <- log(runif(N))
      acc <- (u < pmin(0, loga))
      if (any(acc)) {
        phi[acc, ] <- prop[acc, , drop = FALSE]
        logprior[acc] <- lp_prop[acc]
        loglik[acc] <- ll_prop[acc]
        accepted <- accepted + sum(acc)
      }
      nprop <- nprop + N
    }
  }

  list(
    phi = phi,
    loglik = loglik,
    logprior = logprior,
    aux_state = aux_state,
    acc_rate = accepted / max(1L, nprop),
    rw_accept_rate = accepted / max(1L, nprop),
    indep_accept_rate = NA_real_
  )
}

# --------------- Per-subject IS diagnostics at a given φ ----------------
.cache_logw_at_phi_batch <- function(cache, phi_row, b,
                                     log_prior_theta_given_phi_mat = NULL) {
  B <- cache$batches[[b]]
  if (!is.null(cache$gaussian_map) && is.null(log_prior_theta_given_phi_mat)) {
    pars <- cache$gaussian_map(phi_row)
    lprior <- .log_prior_gauss_mat(B$Theta, pars)
  } else {
    stopifnot(!is.null(log_prior_theta_given_phi_mat))
    lprior <- log_prior_theta_given_phi_mat(B$Theta, phi_row, aux = NULL)
  }
  B$log_py + lprior - B$log_r
}

.is_stats_one_batch <- function(cache, phi_row, b,
                                log_prior_theta_given_phi_mat = NULL) {
  lw <- .cache_logw_at_phi_batch(cache, phi_row, b, log_prior_theta_given_phi_mat)
  a  <- max(lw)
  w  <- exp(lw - a)
  logZ_hat <- a + log(mean(w))
  w <- w / sum(w)
  ess <- 1 / sum(w * w)
  khat <- NA_real_
  if (requireNamespace("loo", quietly = TRUE)) {
    ps <- tryCatch(loo::psis(lw), error = function(e) NULL)
    if (!is.null(ps)) {
      khat <- tryCatch(mean(loo::pareto_k_values(ps)), error = function(e) NA_real_)
    }
  }
  list(logZ_hat = logZ_hat, ess = ess, M = cache$batches[[b]]$M, k = khat)
}

.is_diagnostics_batch <- function(caches, phi_row, b_idx,
                                  log_prior_theta_given_phi_mat = NULL) {
  S <- length(caches)
  out <- vector("list", S)
  for (i in seq_len(S)) {
    out[[i]] <- .is_stats_one_batch(caches[[i]], phi_row, b_idx[i],
                                    log_prior_theta_given_phi_mat)
  }
  data.frame(
    i = seq_len(S),
    logZ_hat = vapply(out, `[[`, numeric(1), "logZ_hat"),
    ess      = vapply(out, `[[`, numeric(1), "ess"),
    M        = vapply(out, `[[`, numeric(1), "M"),
    psis_k   = vapply(out, `[[`, numeric(1), "k")
  )
}

# -------- Replicate SD of log Z_i(φ) over batches (no model calls) -------
.rep_logZi_sd_over_batches <- function(caches, phi_row, B_reps = 2L,
                                       log_prior_theta_given_phi_mat = NULL,
                                       gaussian_map_fn = NULL,
                                       rng_seed = NULL) {
  if (!is.null(rng_seed)) set.seed(rng_seed)
  S <- length(caches)
  out_sd <- numeric(S)
  for (i in seq_len(S)) {
    K <- caches[[i]]$K
    if (B_reps <= 1L || K <= 1L) {
      out_sd[i] <- NA_real_
    } else {
      vals <- numeric(B_reps)
      for (b in seq_len(B_reps)) {
        bb <- sample.int(K, 1L)
        if (!is.null(gaussian_map_fn)) {
          vals[b] <- log_marginal_unbiased_gaussian_batch(caches[[i]], phi_row, bb)
        } else {
          vals[b] <- log_marginal_unbiased_batch(caches[[i]], phi_row, bb,
                                                 log_prior_theta_given_phi_mat)
        }
      }
      out_sd[i] <- stats::sd(vals)
    }
  }
  out_sd
}

# --------- batch-index refresh (cheap correlated PM; no model calls) -----
# If the cache-side helper isn't available, define a local fallback.
if (!exists("refresh_batch_indices", mode = "function")) {
  refresh_batch_indices <- function(b_vec, K, frac = 0.2) {
    n <- length(b_vec); R <- ceiling(n * .clamp(frac, 0, 1))
    if (R <= 0L) return(b_vec)
    idx <- sample.int(n, R)
    b_vec[idx] <- sample.int(K, R, replace = TRUE)
    b_vec
  }
}

# ------------------------------ Rejuvenation -----------------------------
# One sweep of adaptive MH on φ with RW + independence mixture (correlated PM; b_idx fixed)
.rejuvenate_rw_batch <- function(phi, loglik, logprior, w, lambda,
                                 caches, b_state, logprior_phi,
                                 log_prior_theta_given_phi_mat = NULL,
                                 gaussian_map_fn = NULL,
                                 pm_mode = c("strict", "fast"),
                                 n_moves = 2L,
                                 rw_scale = 1.0,
                                 weak_dim_idx = integer(0),
                                 rw_expand_factor = 2.5,
                                 elite_mix = NULL,
                                 hist_mix = NULL,
                                 indep_prob = 0.30,
                                 indep_t_df = 4L,
                                 indep_t_prob = 0.85,
                                 rng_seed = NULL) {
  pm_mode <- match.arg(pm_mode)
  if (!is.null(rng_seed)) set.seed(rng_seed)
  N <- nrow(phi); d <- ncol(phi)

  S <- weighted_cov(phi, w)
  S <- as.matrix(Matrix::nearPD(S, conv.tol = 1e-7)$mat)
  diag(S) <- pmax(diag(S), 1e-8)
  ev <- eigen(S, symmetric = TRUE)
  lamv <- pmax(ev$values, 1e-8)
  if (length(weak_dim_idx)) lamv[weak_dim_idx] <- lamv[weak_dim_idx] * rw_expand_factor
  Sprop <- ev$vectors %*% diag(lamv, d) %*% t(ev$vectors)
  L <- tryCatch(chol(Sprop + diag(1e-8, d)), error = function(e) chol(Sprop + diag(1e-6, d)))
  step <- (rw_scale / sqrt(max(d, 1))) * L

  rlogsumexp2 <- function(a, b) { m <- pmax(a, b); m + log(exp(a - m) + exp(b - m)) }
  log_q_mixture <- function(Phi, elite, hist, indep_t_prob, indep_t_df) {
    if (.is_empty_mix(elite)) elite <- .default_std_normal_mix(ncol(Phi))
    lq_el_n <- gmm_logpdf_Z_vec(Phi, elite$meansZ, elite$cache)
    lq_el_t <- dmvt_mixture_logpdf_Z_vec(Phi, elite$meansZ, elite$cache, indep_t_df)
    lq_el   <- rlogsumexp2(lq_el_n + log1p(-indep_t_prob), lq_el_t + log(indep_t_prob))
    if (is.null(hist) || is.null(hist$mix) || .is_empty_mix(hist$mix)) return(lq_el)
    hm <- hist$mix
    lq_hi_n <- gmm_logpdf_Z_vec(Phi, hm$meansZ, hm$cache)
    lq_hi_t <- dmvt_mixture_logpdf_Z_vec(Phi, hm$meansZ, hm$cache, indep_t_df)
    lq_hi   <- rlogsumexp2(lq_hi_n + log1p(-indep_t_prob), lq_hi_t + log(indep_t_prob))
    rlogsumexp2(lq_el + log1p(-hist$prob), lq_hi + log(hist$prob))
  }
  sample_from_q_phi <- function(M, elite, hist, indep_t_prob, indep_t_df) {
    if (.is_empty_mix(elite)) elite <- .default_std_normal_mix(d)
    use_hist <- !is.null(hist) && !is.null(hist$mix) && !.is_empty_mix(hist$mix)
    hprob <- if (use_hist) hist$prob else 0.0
    take_hist <- runif(M) < hprob
    t_flags   <- runif(M) < indep_t_prob
    out <- matrix(NA_real_, M, d)
    em <- elite; hm <- if (use_hist) hist$mix else NULL
    n_hi_t <- sum(take_hist & t_flags)
    n_hi_n <- sum(take_hist & !t_flags)
    n_el_t <- sum(!take_hist & t_flags)
    n_el_n <- sum(!take_hist & !t_flags)
    if (n_hi_n) out[take_hist & !t_flags, ] <- sample_gmm_Z_qmc(n_hi_n, hm$meansZ, hm$cache)
    if (n_hi_t) out[take_hist &  t_flags, ] <- rmvt_mixture_Z_qmc(n_hi_t, hm$meansZ, hm$cache, nu = indep_t_df)
    if (n_el_n) out[!take_hist & !t_flags, ] <- sample_gmm_Z_qmc(n_el_n, em$meansZ, em$cache)
    if (n_el_t) out[!take_hist &  t_flags, ] <- rmvt_mixture_Z_qmc(n_el_t, em$meansZ, em$cache, nu = indep_t_df)
    colnames(out) <- colnames(phi)
    out
  }

  accepted <- 0L
  nacc_rw <- 0L; nacc_id <- 0L
  nprop_rw <- 0L; nprop_id <- 0L

  eval_ll <- function(Phi, b_sub) {
    .outer_loglik_hat_particles(
      Phi = Phi,
      caches = caches,
      b_state = b_sub,
      log_prior_theta_given_phi_mat = log_prior_theta_given_phi_mat,
      gaussian_map_fn = gaussian_map_fn,
      pm_mode = pm_mode
    )
  }

  for (m in seq_len(n_moves)) {
    choose_id <- runif(N) < indep_prob
    idx_id <- which(choose_id)
    idx_rw <- which(!choose_id)

    # Independence proposals (independent of current state; keep b_idx fixed for PM)
    if (length(idx_id)) {
      prop_id <- sample_from_q_phi(length(idx_id), elite_mix, hist_mix, indep_t_prob, indep_t_df)
      lp_prop_id <- apply(prop_id, 1L, logprior_phi)
      b_id <- if (pm_mode == "strict") b_state[idx_id, , drop = FALSE] else b_state
      ll_prop_id <- eval_ll(prop_id, b_id)
      lq_curr <- log_q_mixture(phi[idx_id, , drop = FALSE], elite_mix, hist_mix, indep_t_prob, indep_t_df)
      lq_prop <- log_q_mixture(prop_id, elite_mix, hist_mix, indep_t_prob, indep_t_df)
      a_id <- (lp_prop_id - logprior[idx_id]) + lambda * (ll_prop_id - loglik[idx_id]) + (lq_curr - lq_prop)
      u_id <- log(runif(length(idx_id)))
      acc_id_idx <- which(u_id < pmin(0, a_id))
      if (length(acc_id_idx)) {
        which_acc <- idx_id[acc_id_idx]
        phi[which_acc, ]    <- prop_id[acc_id_idx, , drop = FALSE]
        loglik[which_acc]   <- ll_prop_id[acc_id_idx]
        logprior[which_acc] <- lp_prop_id[acc_id_idx]
        accepted <- accepted + length(acc_id_idx)
        nacc_id <- nacc_id + length(acc_id_idx)
      }
      nprop_id <- nprop_id + length(idx_id)
    }

    # Random-walk proposals
    if (length(idx_rw)) {
      Z <- matrix(rnorm(length(idx_rw) * d), length(idx_rw), d)
      prop <- phi[idx_rw, , drop = FALSE] + Z %*% t(step)
      lp_prop <- apply(prop, 1L, logprior_phi)
      b_rw <- if (pm_mode == "strict") b_state[idx_rw, , drop = FALSE] else b_state
      ll_prop <- eval_ll(prop, b_rw)
      a <- (lp_prop - logprior[idx_rw]) + lambda * (ll_prop - loglik[idx_rw])
      u <- log(runif(length(idx_rw)))
      acc <- (u < pmin(0, a))
      if (any(acc)) {
        which_acc <- idx_rw[acc]
        phi[which_acc, ]    <- prop[acc, , drop = FALSE]
        loglik[which_acc]   <- ll_prop[acc]
        logprior[which_acc] <- lp_prop[acc]
        accepted <- accepted + sum(acc)
        nacc_rw <- nacc_rw + sum(acc)
      }
      nprop_rw <- nprop_rw + length(idx_rw)
    }
  }
  total_props <- N * n_moves
  list(phi = phi, loglik = loglik, logprior = logprior, b_state = b_state,
       acc_rate = accepted / max(total_props, 1L),
       rw_accept_rate = if (nprop_rw>0) nacc_rw/nprop_rw else NA_real_,
       indep_accept_rate = if (nprop_id>0) nacc_id/nprop_id else NA_real_)
}

.outer_smc_phi_batch_adaptive_pm <- function(
    surrogates,
    rprior_phi,
    logprior_phi,
    log_prior_theta_given_phi_mat = NULL,
    gaussian_map_fn = NULL,
    rtheta_given_phi = NULL,
    data_list = NULL,
    loglik_fn = NULL,
    M = 2000L,
    cess_target = 0.95,
    resample_threshold = 0.5,
    resampling = c("systematic", "multinomial"),
    n_moves = 2L,
    rw_scale_init = 1.2,
    pm_mode = c("strict", "fast"),
    pm_aux_mode = c("rng_stream", "batch_idx"),
    adaptive_pm_control = list(),
    max_rounds = 200L,
    refresh_batches_after_resample = FALSE,
    refresh_batches_each_round = FALSE,
    refresh_batches_frac = 0.02,
    block_refresh_every = 5L,
    block_refresh_frac = 0.10,
    diag_enable = TRUE,
    diag_probe = c("wmean", "best"),
    collect_round_diagnostics = FALSE,
    seed = 123,
    verbose = TRUE
) {
  resampling <- match.arg(resampling)
  diag_probe <- match.arg(diag_probe)
  pm_mode <- match.arg(pm_mode)
  pm_aux_mode <- match.arg(pm_aux_mode)
  if (pm_aux_mode != "rng_stream" && verbose) {
    message("[adaptive_pm] pm_aux_mode='", pm_aux_mode, "' requested; using RNG streams in adaptive PM mode.")
  }
  set.seed(as.integer(seed))
  t_start <- proc.time()[3]

  S <- length(surrogates)
  if (S <= 0) stop("Empty 'surrogates' list.")
  if (is.null(data_list) || is.null(loglik_fn)) {
    stop("adaptive_pm mode requires data_list and loglik_fn.")
  }
  if (length(data_list) != S) {
    stop("length(data_list) must equal number of surrogates S = ", S)
  }
  for (i in seq_len(S)) {
    if (!inherits(surrogates[[i]], "subject_surrogate_pm")) {
      stop("adaptive_pm mode expects surrogates of class 'subject_surrogate_pm'.")
    }
  }
  if (is.null(gaussian_map_fn) && is.null(log_prior_theta_given_phi_mat)) {
    stop("Provide either gaussian_map_fn or log_prior_theta_given_phi_mat.")
  }
  if (is.null(gaussian_map_fn) && is.null(rtheta_given_phi)) {
    stop("adaptive_pm strict mode requires rtheta_given_phi or gaussian_map_fn.")
  }

  ctl <- modifyList(
    list(
      M_default = 128L,
      initial_subject_M = NULL,
      M_levels = c(64L, 128L, 256L, 512L),
      ess_norm_threshold = 0.15,
      khat_threshold = 0.80,
      var_log_target_low = 0.30,
      var_log_target_high = 1.00,
      var_relax_ess = 0.35,
      decrease_enabled = TRUE,
      M_update_every = 1L,
      adapt_until_round = Inf,
      freeze_after_round = Inf,
      da_enable = (pm_mode == "strict"),
      da_M = 32L,
      max_retries = 2L,
      min_log_mhat = -1e12,
      completion_mode = "default",  # "default" uses max_rounds; "lambda1" ignores max_rounds
      max_wall_time_sec = Inf,      # watchdog (seconds), active for both completion modes
      checkpoint_enable = FALSE,
      checkpoint_path = NULL,
      checkpoint_every_rounds = 5L,
      checkpoint_include_state = FALSE
    ),
    adaptive_pm_control
  )
  completion_mode <- match.arg(as.character(ctl$completion_mode), c("default", "lambda1"))
  max_wall_time_sec <- as.numeric(ctl$max_wall_time_sec)
  if (!is.finite(max_wall_time_sec) || max_wall_time_sec <= 0) max_wall_time_sec <- Inf

  adapt_until_round <- suppressWarnings(as.integer(ctl$adapt_until_round))
  if (!is.finite(adapt_until_round) || is.na(adapt_until_round) || adapt_until_round < 1L) {
    adapt_until_round <- .Machine$integer.max
  }
  freeze_after_round <- suppressWarnings(as.integer(ctl$freeze_after_round))
  if (is.finite(freeze_after_round) && !is.na(freeze_after_round) && freeze_after_round >= 1L) {
    adapt_until_round <- min(adapt_until_round, freeze_after_round)
  }

  if (!is.null(ctl$initial_subject_M)) {
    init_M <- as.integer(ctl$initial_subject_M)
    if (length(init_M) != S) {
      stop("adaptive_pm_control$initial_subject_M must be length S (", S, ").")
    }
    subject_M <- pmax(16L, init_M)
  } else {
    subject_M <- rep(as.integer(ctl$M_default), S)
  }
  callbacks <- .adaptive_pm_callbacks(
    log_prior_theta_given_phi_mat = log_prior_theta_given_phi_mat,
    gaussian_map_fn = gaussian_map_fn,
    rtheta_given_phi = rtheta_given_phi,
    data_list = data_list,
    loglik_fn = loglik_fn
  )

  aux_state <- init_pm_aux_state(
    M_particles = as.integer(M),
    S_subjects = S,
    seed = as.integer(seed + 4049L),
    M_default = subject_M,
    pm_mode = pm_mode
  )

  phi <- rprior_phi(M)
  if (!is.matrix(phi)) phi <- matrix(phi, nrow = M)
  logprior <- apply(phi, 1L, logprior_phi)
  loglik <- .outer_loglik_hat_particles_adaptive_pm(
    Phi = phi,
    surrogates = surrogates,
    aux_state = aux_state,
    callbacks = callbacks,
    control = ctl,
    pm_mode = pm_mode
  )
  if (!all(is.finite(logprior))) stop("Non-finite log prior at initialization (adaptive_pm).")
  if (!all(is.finite(loglik))) stop("Non-finite inner likelihood at initialization (adaptive_pm).")

  lambda <- 0.0
  w <- rep(1 / M, M)
  logZ <- 0.0
  logZ_increments <- numeric()
  mcse_var_accum <- 0.0

  target_acc <- 0.234
  log_rw_scale <- log(rw_scale_init)
  rm_gain <- 0.05
  lambda_hist <- lambda
  acc_hist <- numeric()
  round <- 0L
  round_diag <- list()
  termination_reason <- NA_character_
  lambda_target_reached <- FALSE

  .adaptive_pm_checkpoint <- function(tag = "round") {
    if (!isTRUE(ctl$checkpoint_enable)) return(invisible(NULL))
    path <- as.character(ctl$checkpoint_path %||% "")
    if (!nzchar(path)) return(invisible(NULL))
    ck <- list(
      tag = tag,
      timestamp = Sys.time(),
      seed = as.integer(seed),
      pm_mode = pm_mode,
      completion_mode = completion_mode,
      max_wall_time_sec = max_wall_time_sec,
      round = as.integer(round),
      lambda = as.numeric(lambda),
      log_evidence = as.numeric(logZ),
      mcse_log_evidence = sqrt(mcse_var_accum),
      subject_M = as.integer(subject_M)
    )
    if (isTRUE(ctl$checkpoint_include_state)) {
      ck$state <- list(
        phi = phi,
        w = w,
        loglik = loglik,
        logprior = logprior,
        aux_state = aux_state,
        log_rw_scale = as.numeric(log_rw_scale),
        lambda_hist = lambda_hist,
        acc_hist = acc_hist
      )
    }
    tryCatch(
      saveRDS(ck, path),
      error = function(e) {
        if (isTRUE(verbose)) message("[adaptive_pm] checkpoint write failed: ", conditionMessage(e))
      }
    )
    invisible(NULL)
  }

  while (lambda < 1 - 1e-12 &&
         (completion_mode == "lambda1" || round < max_rounds)) {
    if ((proc.time()[3] - t_start) >= max_wall_time_sec) {
      termination_reason <- "max_wall_time"
      if (isTRUE(verbose)) cat("\n[adaptive_pm] stopping at watchdog limit (max_wall_time_sec).\n")
      break
    }
    round <- round + 1L
    t_round_start <- proc.time()[3]
    if (verbose) cat(sprintf("\n[Round %d] λ=%.3f  ", round, lambda))

    h <- loglik
    target_cess <- cess_target_at_lambda(lambda)
    delta <- next_lambda_via_rCESS_stat(w, h, lambda, target = target_cess) - lambda
    if (!is.finite(delta)) delta <- min(1 - lambda, 1e-4)
    if (delta <= 1e-8) delta <- min(1 - lambda, 1e-4)
    lambda_new <- min(1.0, lambda + delta)

    m <- max(h)
    log_u <- delta * (h - m)
    logZ_inc <- logsumexp_w(log_u, w) + delta * m
    logZ <- logZ + logZ_inc
    logZ_increments <- c(logZ_increments, logZ_inc)

    u <- exp(log_u)
    mu1 <- sum(w * u)
    mu2 <- sum(w * u * u)
    Neff <- 1 / sum(w * w)
    var_logZ_inc <- (mu2 - mu1^2) / (max(Neff, 1) * max(mu1^2, .Machine$double.eps))
    mcse_var_accum <- mcse_var_accum + max(var_logZ_inc, 0)

    w <- w * u
    w <- w / sum(w)
    ess_frac <- 1 / sum(w * w) / M
    if (verbose) cat(sprintf("Δ=%.4f | ESS/N=%.3f | logZ+=%.4f\n", delta, ess_frac, logZ_inc))

    if (ess_frac < resample_threshold) {
      idx <- if (resampling == "multinomial") .resample_multinomial(w) else {
        set.seed(seed + 991 * round)
        mu_phi <- colSums(phi * w)
        Xc <- sweep(phi, 2L, mu_phi, `-`)
        Sphi <- weighted_cov(phi, w)
        Sphi <- as.matrix(Matrix::nearPD(Sphi, conv.tol = 1e-7)$mat)
        diag(Sphi) <- pmax(diag(Sphi), 1e-10)
        R <- tryCatch(chol(Sphi), error = function(e) chol(Sphi + diag(1e-8, ncol(Sphi))))
        Zphi <- t(backsolve(R, t(Xc), transpose = TRUE))
        ord_h <- hilbert_sort_order(Zphi, bits = 16L)
        ws <- w[ord_h]
        sel_sorted <- stratified_resample_sorted(ws, deterministic = FALSE)
        ord_h[sel_sorted]
      }
      phi <- phi[idx, , drop = FALSE]
      loglik <- loglik[idx]
      logprior <- logprior[idx]
      if (pm_mode == "strict") {
        aux_state$seed <- aux_state$seed[idx, , drop = FALSE]
        aux_state$M_alloc <- aux_state$M_alloc[idx, , drop = FALSE]
      }
      w <- rep(1 / M, M)
      if (verbose) cat(sprintf("  Resampled (%s).\n", resampling))

      if (refresh_batches_after_resample) {
        aux_state <- refresh_pm_aux_state(
          aux_state,
          frac = refresh_batches_frac,
          seed = as.integer(seed + 1229L * round)
        )
        loglik <- .outer_loglik_hat_particles_adaptive_pm(
          Phi = phi,
          surrogates = surrogates,
          aux_state = aux_state,
          callbacks = callbacks,
          control = ctl,
          pm_mode = pm_mode
        )
      }
    }

    diag_median_ess_norm <- NA_real_
    diag_max_psis_k <- NA_real_
    diag_median_var_proxy <- NA_real_
    diag_sum_sd2 <- NA_real_
    diag_n_offenders <- NA_integer_
    need_m_update <- (round <= adapt_until_round) &&
      (round %% as.integer(max(1L, ctl$M_update_every)) == 0L)
    need_probe <- isTRUE(diag_enable) || need_m_update
    sdiag <- NULL
    if (need_probe) {
      probe_phi <- .choose_probe_phi(phi, w, logprior, loglik, lambda_new, diag_probe)
      probe_aux <- .pm_aux_probe_row(phi, w, logprior, loglik, lambda_new, aux_state, pm_mode)
      sdiag <- .adaptive_pm_probe_subject_diag(probe_phi, surrogates, probe_aux, callbacks, ctl)
      if (is.data.frame(sdiag) && nrow(sdiag) > 0L) {
        diag_median_ess_norm <- suppressWarnings(stats::median(sdiag$ess_norm, na.rm = TRUE))
        diag_max_psis_k <- suppressWarnings(max(sdiag$khat, na.rm = TRUE))
        diag_median_var_proxy <- suppressWarnings(
          stats::median(sdiag$var_proxy[is.finite(sdiag$var_proxy)], na.rm = TRUE)
        )
        diag_sum_sd2 <- suppressWarnings(sum(sdiag$var_proxy[is.finite(sdiag$var_proxy)]))
        bad <- which((is.finite(sdiag$ess_norm) & sdiag$ess_norm < ctl$ess_norm_threshold) |
                       (is.finite(sdiag$khat) & sdiag$khat > ctl$khat_threshold) |
                       (!is.finite(sdiag$var_proxy) | sdiag$var_proxy > ctl$var_log_target_high))
        diag_n_offenders <- as.integer(length(bad))
      }
      if (isTRUE(diag_enable) && verbose) {
        cat(sprintf("  [AdaptivePM] median ESS_i/M=%.2f | median var=%.2f | max k=%.2f | offenders=%d\n",
                    diag_median_ess_norm, diag_median_var_proxy, diag_max_psis_k, diag_n_offenders))
      }
      if (need_m_update && is.data.frame(sdiag) && nrow(sdiag) > 0L) {
        subject_M_new <- as.integer(.adaptive_pm_update_subject_M(subject_M, sdiag, ctl))
        if (length(subject_M_new) == length(subject_M) && any(subject_M_new != subject_M)) {
          subject_M <- subject_M_new
          aux_state <- .pm_aux_apply_subject_M(aux_state, subject_M, pm_mode = pm_mode)
          # Keep MH target consistent after changing the PM estimator definition.
          loglik <- .outer_loglik_hat_particles_adaptive_pm(
            Phi = phi,
            surrogates = surrogates,
            aux_state = aux_state,
            callbacks = callbacks,
            control = ctl,
            pm_mode = pm_mode
          )
        } else {
          subject_M <- subject_M_new
        }
      }
    }

    move <- .rejuvenate_rw_adaptive_pm(
      phi = phi,
      loglik = loglik,
      logprior = logprior,
      w = w,
      lambda = lambda_new,
      surrogates = surrogates,
      aux_state = aux_state,
      logprior_phi = logprior_phi,
      callbacks = callbacks,
      control = ctl,
      pm_mode = pm_mode,
      n_moves = n_moves,
      rw_scale = exp(log_rw_scale),
      rng_seed = as.integer(seed + 97L * round)
    )
    phi <- move$phi
    loglik <- move$loglik
    logprior <- move$logprior
    aux_state <- move$aux_state
    acc_hist <- c(acc_hist, move$acc_rate)
    log_rw_scale <- .clamp(log_rw_scale + rm_gain * (move$acc_rate - target_acc), log(0.05), log(2.5))
    if (verbose) {
      cat(sprintf("  MH acc=%.3f | rw_scale=%.3f\n", move$acc_rate, exp(log_rw_scale)))
    }

    if (refresh_batches_each_round) {
      aux_state <- refresh_pm_aux_state(
        aux_state,
        frac = refresh_batches_frac,
        seed = as.integer(seed + 9929L * round)
      )
      loglik <- .outer_loglik_hat_particles_adaptive_pm(
        Phi = phi,
        surrogates = surrogates,
        aux_state = aux_state,
        callbacks = callbacks,
        control = ctl,
        pm_mode = pm_mode
      )
    }

    if (pm_mode == "strict" && block_refresh_every > 0L &&
        (round %% as.integer(max(1L, block_refresh_every)) == 0L)) {
      aux_state <- refresh_pm_aux_state(
        aux_state,
        frac = block_refresh_frac,
        seed = as.integer(seed + 7001L * round)
      )
      loglik <- .outer_loglik_hat_particles_adaptive_pm(
        Phi = phi,
        surrogates = surrogates,
        aux_state = aux_state,
        callbacks = callbacks,
        control = ctl,
        pm_mode = pm_mode
      )
    }

    lambda <- lambda_new
    lambda_hist <- c(lambda_hist, lambda)

    if (isTRUE(ctl$checkpoint_enable) &&
        isTRUE(round %% as.integer(max(1L, ctl$checkpoint_every_rounds)) == 0L)) {
      .adaptive_pm_checkpoint(tag = "round")
    }

    if (isTRUE(collect_round_diagnostics)) {
      round_diag[[length(round_diag) + 1L]] <- list(
        round = round,
        elapsed_sec = as.numeric(proc.time()[3] - t_round_start),
        lambda = lambda,
        delta = delta,
        ess_frac = ess_frac,
        logZ_inc = logZ_inc,
        mh_acc = move$acc_rate,
        diag_median_ess_norm = diag_median_ess_norm,
        diag_median_var_proxy = diag_median_var_proxy,
        diag_max_psis_k = diag_max_psis_k,
        diag_sum_sd2 = diag_sum_sd2,
        diag_n_offenders = diag_n_offenders
      )
    }
  }

  if (lambda >= 1 - 1e-12) {
    lambda_target_reached <- TRUE
    if (is.na(termination_reason)) termination_reason <- "lambda_reached"
  } else if (is.na(termination_reason)) {
    termination_reason <- if (completion_mode == "default" && round >= max_rounds) "max_rounds" else "stopped_early"
  }
  .adaptive_pm_checkpoint(tag = paste0("final_", termination_reason))

  list(
    phi = phi,
    w = w,
    loglik = loglik,
    logprior = logprior,
    log_evidence = logZ,
    mcse_log_evidence = sqrt(mcse_var_accum),
    logZ_increments = logZ_increments,
    b_idx = NULL,
    B_idx = NULL,
    aux_state = aux_state,
    meta = list(
      pm_mode = pm_mode,
      pm_aux_mode = "rng_stream",
      inner_ll_mode = "adaptive_pm",
      lambda_hist = lambda_hist,
      acc_hist = acc_hist,
      rw_scale_final = exp(log_rw_scale),
      resampling = resampling,
      refresh_batches_after_resample = refresh_batches_after_resample,
      refresh_batches_each_round = refresh_batches_each_round,
      refresh_batches_frac = refresh_batches_frac,
      block_refresh_every = block_refresh_every,
      block_refresh_frac = block_refresh_frac,
      subject_M_final = subject_M,
      adaptive_pm_control = ctl,
      adapt_until_round = as.integer(adapt_until_round),
      completion_mode = completion_mode,
      max_wall_time_sec = max_wall_time_sec,
      completion_reason = termination_reason,
      lambda_target_reached = lambda_target_reached,
      rounds_completed = as.integer(round),
      checkpoint_enable = isTRUE(ctl$checkpoint_enable),
      checkpoint_path = ctl$checkpoint_path %||% NULL,
      elapsed_sec = as.numeric(proc.time()[3] - t_start),
      round_diagnostics = round_diag
    )
  )
}

materialize_alpha_from_indices <- function(alpha_bank, index_mat) {
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

refresh_alpha_indices_exact <- function(index_mat, n_bank, frac = 0.10, seed = NULL) {
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

compute_alpha_suff_stats <- function(alpha_state, X_list = NULL) {
  alpha_state <- array(alpha_state, dim = dim(alpha_state))
  d <- dim(alpha_state)[3L]
  if (!is.null(X_list)) {
    stop("compute_alpha_suff_stats only supports intercept-only form in the current skeleton.")
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

update_alpha_suff_stats_one <- function(stats_obj, alpha_old, alpha_new, x_i = NULL) {
  if (!is.null(x_i)) {
    stop("update_alpha_suff_stats_one only supports intercept-only form in the current skeleton.")
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

log_tempered_hier_contrib <- function(phi_row, alpha_mat, beta, gaussian_map_fn, working_priors) {
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

.log_alpha_given_phi_sum_mat <- function(Phi, alpha_state, gaussian_map_fn) {
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

init_extended_exact_particles <- function(local_objs, phi_particles, phi_anchor = NULL, init_ctl = list()) {
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
  init_mode <- match.arg(
    as.character(ctl$mode),
    c("surrogate_is", "outer_resample")
  )
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
      prop <- make_theta_proposal_exact_init(loc, init_ctl = ctl$proposal_control)
      Theta_j <- prop$draw(N, seed = seed_j)
      Theta_j <- as.matrix(Theta_j)
      if (nrow(Theta_j) != N) stop("Initialization proposal returned wrong number of draws.")
    } else {
      set.seed(seed_j)
      idx <- sample.int(nrow(loc$outer_theta), size = N, replace = TRUE, prob = loc$outer_w)
      Theta_j <- as.matrix(loc$outer_theta[idx, , drop = FALSE])
    }
    alpha_state[, j, ] <- Theta_j
    subject_loglik[, j] <- ll_parallel(Theta_j, loc$data, loc$loglik_fn, n_cores = 1L)
    subject_logg[, j] <- log_working_prior_gaussian_mat(Theta_j, loc$working_prior)
    if (identical(init_mode, "surrogate_is")) {
      subject_logq[, j] <- prop$log_q(Theta_j)
    } else {
      subject_logq[, j] <- 0
    }
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
  ess_init <- 1 / sum(w * w)
  ess_init_frac <- ess_init / max(N, 1L)
  max_weight_init <- max(w)

  list(
    alpha_state = alpha_state,
    subject_loglik = subject_loglik,
    subject_logg = subject_logg,
    subject_logq = subject_logq,
    log_weight = logw_init,
    w = w,
    init_mode = init_mode,
    ess_init = as.numeric(ess_init),
    ess_init_frac = as.numeric(ess_init_frac),
    max_weight_init = as.numeric(max_weight_init),
    suff_stats = compute_alpha_suff_stats(alpha_state),
    n_exact_loglik_init = as.integer(N * S),
    phi_anchor = phi_anchor
  )
}

rejuvenate_phi_given_alpha_exact <- function(phi, alpha_state, logprior, w, lambda,
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
    log_alpha_given_phi_sum <- .log_alpha_given_phi_sum_mat(phi, alpha_state, gaussian_map_fn)
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
    loga_prop <- .log_alpha_given_phi_sum_mat(prop, alpha_state, gaussian_map_fn)
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

refresh_alpha_block_exact <- function(phi, alpha_state, local_objs, beta,
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

      prop <- make_theta_proposal_exact_refresh(
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
        stop("refresh_alpha_block_exact requires gaussian_map_fn in the current skeleton.")
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

outer_smc_phi_extended_exact <- function(
    local_objs,
    rprior_phi,
    logprior_phi,
    gaussian_map_fn,
    data_list = NULL,
    loglik_fn = NULL,
    M = 1200L,
    cess_target = 0.95,
    resample_threshold = 0.5,
    resampling = c("systematic", "multinomial"),
    n_moves = 2L,
    rw_scale_init = 1.2,
    max_rounds = 200L,
    refresh_frac = 0.05,
    delayed_accept = TRUE,
    init_control = list(),
    seed = 123L,
    verbose = TRUE
) {
  resampling <- match.arg(resampling)
  set.seed(as.integer(seed))
  if (!length(local_objs)) stop("local_objs must not be empty.")
  if (!all(vapply(local_objs, inherits, logical(1), what = "subject_exact_local"))) {
    stop("outer_smc_phi_extended_exact expects local_objs of class 'subject_exact_local'.")
  }

  init_ctl <- modifyList(
    list(
      mode = "surrogate_is",
      phi_mode = "prior",
      phi_center = NULL,
      phi_sd = NULL,
      proposal_control = list()
    ),
    init_control
  )
  phi_mode <- match.arg(as.character(init_ctl$phi_mode), c("prior", "centered"))
  logq_phi <- NULL
  if (identical(phi_mode, "centered")) {
    center <- as.numeric(init_ctl$phi_center)
    if (!length(center)) stop("init_control$phi_center must be provided when phi_mode='centered'.")
    d_phi <- length(center)
    sd_phi <- as.numeric(init_ctl$phi_sd %||% rep(0.10, d_phi))
    if (length(sd_phi) == 1L) sd_phi <- rep(sd_phi, d_phi)
    phi <- matrix(rnorm(M * d_phi), nrow = M, ncol = d_phi)
    phi <- sweep(phi, 2L, sd_phi, `*`)
    phi <- sweep(phi, 2L, center, `+`)
    logq_phi <- Reduce(
      `+`,
      lapply(seq_len(d_phi), function(j) {
        stats::dnorm(phi[, j], mean = center[j], sd = sd_phi[j], log = TRUE)
      })
    )
  } else {
    phi <- rprior_phi(M)
  }
  if (!is.matrix(phi)) phi <- matrix(phi, nrow = M)
  init <- init_extended_exact_particles(
    local_objs = local_objs,
    phi_particles = phi,
    init_ctl = modifyList(
      init_ctl,
      list(seed = as.integer(seed + 101L))
    )
  )

  alpha_state <- init$alpha_state
  subject_loglik <- init$subject_loglik
  subject_logg <- init$subject_logg
  logprior <- apply(phi, 1L, logprior_phi)
  log_alpha_given_phi_sum <- .log_alpha_given_phi_sum_mat(phi, alpha_state, gaussian_map_fn)
  logg_sum <- rowSums(subject_logg)

  logw_init <- init$log_weight
  if (!is.null(logq_phi)) {
    logw_init <- logw_init + (logprior - logq_phi)
  }
  logw <- logw_init - max(logw_init)
  w <- exp(logw)
  w <- w / sum(w)
  ess_init <- 1 / sum(w * w)
  ess_init_frac <- ess_init / max(M, 1L)
  max_weight_init <- max(w)
  logZ <- logsumexp(logw_init) - log(M)
  lambda <- 0.0
  lambda_hist <- lambda
  acc_hist <- numeric()
  ess_hist <- numeric()
  logZ_increments <- numeric()
  log_rw_scale <- log(rw_scale_init)
  round_diag <- list()
    total_exact_loglik <- as.integer(init$n_exact_loglik_init)
  init_diag <- list(
      init_mode = init$init_mode,
      phi_mode = phi_mode,
      ess_init = as.numeric(ess_init),
      ess_init_frac = as.numeric(ess_init_frac),
      max_weight_init = as.numeric(max_weight_init),
      n_exact_loglik_init = as.integer(init$n_exact_loglik_init)
    )

  round <- 0L
  while (lambda < 1 - 1e-12 && round < max_rounds) {
    round <- round + 1L
    target_cess <- cess_target_at_lambda(lambda)

    # Match the maintained outer SMC paths: normalize/repair weights before
    # solving the next tempering step, because rCESS is defined on normalized w.
    w[!is.finite(w) | w < 0] <- 0
    sw <- sum(w)
    if (!is.finite(sw) || sw <= 0) {
      w <- rep(1 / M, M)
    } else {
      w <- w / sw
    }

    h <- log_alpha_given_phi_sum - logg_sum
    if (!all(is.finite(h))) {
      bad <- !is.finite(h)
      finite_h <- h[!bad]
      if (!length(finite_h)) stop("All X1 tempering statistics are non-finite; cannot continue.")
      h[bad] <- min(finite_h)
    }

    lambda_new <- min(
      next_lambda_via_rCESS_stat(w, h, lambda, target = target_cess),
      1.0
    )
    if (!is.finite(lambda_new) || lambda_new <= lambda) {
      delta <- min(1 - lambda, 1e-4)
      lambda_new <- lambda + delta
    } else {
      delta <- lambda_new - lambda
    }

    m <- max(h)
    log_u <- delta * (h - m)
    logZ_inc <- logsumexp_w(log_u, w) + delta * m
    logZ <- logZ + logZ_inc
    logZ_increments <- c(logZ_increments, logZ_inc)

    logw_new_raw <- log(pmax(w, .Machine$double.eps)) + log_u
    lse_new <- logsumexp(logw_new_raw)
    if (!is.finite(lse_new)) {
      w <- rep(1 / M, M)
    } else {
      logw_new <- logw_new_raw - lse_new
      w <- exp(logw_new)
      sw <- sum(w)
      if (!all(is.finite(w)) || !is.finite(sw) || sw <= 0) {
        w <- rep(1 / M, M)
      } else {
        w <- w / sw
      }
    }
    ess_frac <- (1 / sum(w * w)) / M
    ess_hist <- c(ess_hist, ess_frac)
    if (verbose) {
      cat(sprintf("[ExtendedExact Round %d] lambda %.4f -> %.4f | ESS/N=%.3f | logZ+=%.4f\n",
                  round, lambda, lambda_new, ess_frac, logZ_inc))
    }

    if (ess_frac < resample_threshold) {
      idx <- if (resampling == "multinomial") {
        sample.int(M, size = M, replace = TRUE, prob = w)
      } else {
        stratified_resample_sorted(w, deterministic = FALSE)
      }
      phi <- phi[idx, , drop = FALSE]
      alpha_state <- alpha_state[idx, , , drop = FALSE]
      subject_loglik <- subject_loglik[idx, , drop = FALSE]
      subject_logg <- subject_logg[idx, , drop = FALSE]
      logprior <- logprior[idx]
      log_alpha_given_phi_sum <- log_alpha_given_phi_sum[idx]
      logg_sum <- logg_sum[idx]
      w <- rep(1 / M, M)
    }

    move <- rejuvenate_phi_given_alpha_exact(
      phi = phi,
      alpha_state = alpha_state,
      logprior = logprior,
      w = w,
      lambda = lambda_new,
      logprior_phi = logprior_phi,
      gaussian_map_fn = gaussian_map_fn,
      rw_scale = exp(log_rw_scale),
      n_moves = n_moves,
      log_alpha_given_phi_sum = log_alpha_given_phi_sum,
      rng_seed = as.integer(seed + 97L * round)
    )
    phi <- move$phi
    logprior <- move$logprior
    log_alpha_given_phi_sum <- move$log_alpha_given_phi_sum
    acc_hist <- c(acc_hist, move$acc_rate)
    log_rw_scale <- .clamp(log_rw_scale + 0.05 * (move$acc_rate - 0.234), log(0.05), log(2.5))

    refresh <- refresh_alpha_block_exact(
      phi = phi,
      alpha_state = alpha_state,
      local_objs = local_objs,
      beta = lambda_new,
      refresh_frac = refresh_frac,
      delayed_accept = delayed_accept,
      refresh_ctl = list(),
      subject_loglik = subject_loglik,
      subject_logg = subject_logg,
      gaussian_map_fn = gaussian_map_fn,
      seed = as.integer(seed + 701L * round)
    )
    alpha_state <- refresh$alpha_state
    subject_loglik <- refresh$subject_loglik
    subject_logg <- refresh$subject_logg
    logg_sum <- rowSums(subject_logg)
    log_alpha_given_phi_sum <- .log_alpha_given_phi_sum_mat(phi, alpha_state, gaussian_map_fn)
    total_exact_loglik <- total_exact_loglik + refresh$n_exact_loglik

    lambda <- lambda_new
    lambda_hist <- c(lambda_hist, lambda)
    screen_pass_rate <- refresh$n_screen_pass / max(refresh$n_refresh, 1L)
    final_accept_rate <- refresh$n_final_accept / max(refresh$n_refresh, 1L)
    final_accept_given_screen <- refresh$n_final_accept / max(refresh$n_screen_pass, 1L)
    exact_loglik_rate <- refresh$n_exact_loglik / max(refresh$n_refresh, 1L)
    round_diag[[length(round_diag) + 1L]] <- list(
      round = round,
      lambda = lambda,
      delta = delta,
      target_cess = target_cess,
      ess_frac = ess_frac,
      logZ_inc = logZ_inc,
      phi_acc = move$acc_rate,
      n_refresh = refresh$n_refresh,
      n_screen_pass = refresh$n_screen_pass,
      n_final_accept = refresh$n_final_accept,
      n_exact_loglik = refresh$n_exact_loglik,
      da_screen_pass_rate = screen_pass_rate,
      da_final_accept_rate = final_accept_rate,
      da_final_accept_given_screen = final_accept_given_screen,
      exact_loglik_rate = exact_loglik_rate
    )
  }

  da_screen_pass_rate_mean <- if (length(round_diag)) {
    mean(vapply(round_diag, `[[`, numeric(1), "da_screen_pass_rate"))
  } else {
    NA_real_
  }
  da_final_accept_rate_mean <- if (length(round_diag)) {
    mean(vapply(round_diag, `[[`, numeric(1), "da_final_accept_rate"))
  } else {
    NA_real_
  }
  da_final_accept_given_screen_mean <- if (length(round_diag)) {
    mean(vapply(round_diag, `[[`, numeric(1), "da_final_accept_given_screen"))
  } else {
    NA_real_
  }

  list(
    phi = phi,
    w = w,
    alpha_state = alpha_state,
    subject_loglik = subject_loglik,
    logprior = logprior,
    log_alpha_given_phi_sum = log_alpha_given_phi_sum,
    log_evidence = logZ,
    mcse_log_evidence = NA_real_,
    logZ_increments = logZ_increments,
    meta = list(
      inner_ll_mode = "extended_exact",
      lambda_hist = lambda_hist,
      ess_hist = ess_hist,
      acc_hist = acc_hist,
      rw_scale_final = exp(log_rw_scale),
      final_lambda = lambda,
      lambda_final = lambda,
      rounds = round,
      rounds_completed = as.integer(round),
      init_diagnostics = init_diag,
      exact_loglik_calls = total_exact_loglik,
      refresh_frac = refresh_frac,
      delayed_accept = delayed_accept,
      da_screen_pass_rate_mean = da_screen_pass_rate_mean,
      da_final_accept_rate_mean = da_final_accept_rate_mean,
      da_final_accept_given_screen_mean = da_final_accept_given_screen_mean,
      round_diagnostics = round_diag
    )
  )
}

# ------------------------------- Outer SMC -------------------------------
outer_smc_phi_batch <- function(
    caches,                                  # list of subject caches (length S)
    rprior_phi,                              # function(n) -> matrix(n x dphi)
    logprior_phi,                            # function(phi_row) -> scalar log p(φ)
    # Either provide the generic θ|φ evaluator:
    log_prior_theta_given_phi_mat = NULL,    # function(Theta, phi, aux) -> length M_b
    # Or provide the Gaussian fast-path mapper:
    gaussian_map_fn = NULL,                  # function(phi, dθ) -> list(mu, Sigma_inv, logdet, const)
    rtheta_given_phi = NULL,                 # optional draw callback: function(phi, n, aux=NULL)
    # --- NEW: needed for enrichment / cache rebuilds ---
    data_list = NULL,                        # list of per-subject data (length S)
    loglik_fn = NULL,                        # function(Theta_block, data_i) -> vector log p(y_i|θ)
    n_cores_inner = 1L,                      # parallelism for inner loglik recompute
    M = 2000L,
    cess_target = 0.95,
    resample_threshold = 0.5,
    resampling = c("systematic","multinomial"),
    n_moves = 2L,
    rw_scale_init = 1.3,
    pm_mode = c("strict", "fast"),
    inner_ll_mode = c("fixed_cache", "adaptive_pm"),
    pm_aux_mode = c("batch_idx", "rng_stream"),
    adaptive_pm_control = list(),
    # Independence mixture kernel (new)
    indep_prob = 0.30,
    indep_t_df = 4L,
    indep_t_prob = 0.85,
    G_mix = 12L,
    elite_quantile = 0.40,
    cov_inflation = 4.0,
    base_mix_eps = 0.10,
    # History mixture controls
    hist_mix_enable = TRUE,
    hist_mix_lambda_thresh = 0.50,
    hist_mix_prob = 0.20,
    max_rounds = 200L,
    # Correlated PM batch-index refresh knobs
    refresh_batches_after_resample = FALSE,
    refresh_batches_each_round = FALSE,   # ### PATCH: default off to preserve CRNs
    refresh_batches_frac = 0.02,          # ### PATCH: gentler refresh if enabled
    block_refresh_every = 5L,
    block_refresh_frac = 0.10,
    # Diagnostics
    diag_enable = TRUE,
    diag_probe = c("wmean","best"),
    diag_replicates = 3L,                    # B replicates over batches (no model calls)
    offenders_print_cap = 10L,               # how many worst offenders to print
    # Targeted enrichment knobs (Fixes #2–#5)
    auto_enrich_enable = TRUE,
    auto_enrich_every = 10L,
    auto_enrich_lambda_gate = 0.25,
    auto_enrich_ess_thresh = 0.20,
    auto_enrich_k_thresh   = 0.70,
    auto_enrich_max_units  = 2L,     # ### PATCH: safer default
    eps_prior_anchor = 0.03,         # ### PATCH: base; will adapt per-batch
    weak_inflate_factor = 3.0,       # ### PATCH: toned down
    collect_round_diagnostics = FALSE,
    # Seeding
    seed = 123,
    verbose = TRUE
) {
  resampling <- match.arg(resampling)
  diag_probe <- match.arg(diag_probe)
  pm_mode <- match.arg(pm_mode)
  inner_ll_mode <- match.arg(inner_ll_mode)
  pm_aux_mode <- match.arg(pm_aux_mode)

  if (inner_ll_mode == "adaptive_pm") {
    return(.outer_smc_phi_batch_adaptive_pm(
      surrogates = caches,
      rprior_phi = rprior_phi,
      logprior_phi = logprior_phi,
      log_prior_theta_given_phi_mat = log_prior_theta_given_phi_mat,
      gaussian_map_fn = gaussian_map_fn,
      rtheta_given_phi = rtheta_given_phi,
      data_list = data_list,
      loglik_fn = loglik_fn,
      M = M,
      cess_target = cess_target,
      resample_threshold = resample_threshold,
      resampling = resampling,
      n_moves = n_moves,
      rw_scale_init = rw_scale_init,
      pm_mode = pm_mode,
      pm_aux_mode = pm_aux_mode,
      adaptive_pm_control = adaptive_pm_control,
      max_rounds = max_rounds,
      refresh_batches_after_resample = refresh_batches_after_resample,
      refresh_batches_each_round = refresh_batches_each_round,
      refresh_batches_frac = refresh_batches_frac,
      block_refresh_every = block_refresh_every,
      block_refresh_frac = block_refresh_frac,
      diag_enable = diag_enable,
      diag_probe = diag_probe,
      collect_round_diagnostics = collect_round_diagnostics,
      seed = seed,
      verbose = verbose
    ))
  }

  set.seed(seed)
  t_start <- proc.time()[3]

  S <- length(caches)
  if (S <= 0) stop("Empty 'caches' list.")
  if (is.null(gaussian_map_fn) && is.null(log_prior_theta_given_phi_mat)) {
    stop("Provide either gaussian_map_fn or log_prior_theta_given_phi_mat.")
  }
  # --- NEW: enrichment input sanity ---
  if (auto_enrich_enable) {
    if (is.null(data_list) || is.null(loglik_fn)) {
      message("[Enrich] auto_enrich_enable=TRUE but data_list/loglik_fn not provided; disabling enrichment.")
      auto_enrich_enable <- FALSE
    } else if (length(data_list) != S) {
      stop("[Enrich] length(data_list) must equal number of subjects S = ", S)
    }
  }
  # Self-register Gaussian fast-path onto caches if provided
  if (!is.null(gaussian_map_fn)) {
    for (i in seq_along(caches)) {
      if (is.null(caches[[i]]$gaussian_map)) {
        d_i <- ncol(caches[[i]]$batches[[1]]$Theta)
        caches[[i]] <- register_gaussian_prior_map(
          caches[[i]],
          function(phi) gaussian_map_fn(phi, d_i)
        )
      }
    }
  }

  # 0) Initialize auxiliary batch state
  # fast: vector length S, strict: matrix M x S
  b_state <- .init_aux_batch_state(caches, M = M, pm_mode = pm_mode)

  # 1) Initialize particles from prior
  phi <- rprior_phi(M)               # M x dphi
  if (!is.matrix(phi)) phi <- matrix(phi, nrow = M)
  dphi <- ncol(phi)
  logprior <- apply(phi, 1L, logprior_phi)

  # Initial loglik (batch-indexed)
  loglik <- .outer_loglik_hat_particles(
    Phi = phi, caches = caches, b_state = b_state,
    log_prior_theta_given_phi_mat = log_prior_theta_given_phi_mat,
    gaussian_map_fn = gaussian_map_fn,
    pm_mode = pm_mode
  )

  # ### PATCH: sanity message for per-batch size on weakly ID params
  M_per_hint <- vapply(caches, function(ci) ci$batches[[1]]$M, numeric(1))
  if (verbose && any(M_per_hint < 128)) {
    message(sprintf("[Hint] Some batches are small (min M_per=%d). For weakly identified params, prefer M_per >= 128 (increase M or reduce K).",
                    min(M_per_hint)))
  }

  # Finite checks
  if (!all(is.finite(logprior))) stop("Non-finite log prior at initialization.")
  if (!all(is.finite(loglik)))   stop("Non-finite log inner likelihood at initialization.")

  lambda <- 0.0
  w <- rep(1/M, M)
  logZ <- 0.0
  logZ_increments <- numeric()   # record per-round logZ increments
  mcse_var_accum <- 0.0

  # Adaptive RW scale
  target_acc <- 0.234
  log_rw_scale <- log(rw_scale_init)
  rm_gain <- 0.05

  lambda_hist <- lambda
  acc_hist <- numeric()
  # mixture state (φ-space)
  last_elite_mix_phi <- NULL
  elite_history_phi <- list()

  # Diagnostics store
  diag_list <- list()
  aux_refresh_acc_hist <- numeric()

  round <- 0L
  # rolling acceptances for adaptive n_moves
  last_rw_acc <- NA_real_
  last_id_acc <- NA_real_
  moves_cap <- 8L
  while (lambda < 1 - 1e-12 && round < max_rounds) {
    round <- round + 1L
    t_round_start <- proc.time()[3]
    if (verbose) cat(sprintf("\n[Round %d] λ=%.3f  ", round, lambda))

    # 2) Choose δ via rCESS (from smc_core)
    h <- loglik
    target_cess <- cess_target_at_lambda(lambda)
    delta <- next_lambda_via_rCESS_stat(w, h, lambda, target = target_cess) - lambda
    if (delta <= 1e-8) delta <- min(1 - lambda, 1e-4)
    lambda_new <- min(1.0, lambda + delta)

    # 3) Incremental weighting + evidence
    m <- max(h)
    log_u <- delta * (h - m)
    logZ_inc <- logsumexp_w(log_u, w) + delta * m
    logZ <- logZ + logZ_inc
    logZ_increments <- c(logZ_increments, logZ_inc)

    # Delta-method MCSE accumulation for log evidence
    u <- exp(log_u)
    mu1 <- sum(w * u)
    mu2 <- sum(w * u * u)
    Neff <- 1 / sum(w * w)
    var_logZ_inc <- (mu2 - mu1^2) / (max(Neff, 1) * max(mu1^2, .Machine$double.eps))
    mcse_var_accum <- mcse_var_accum + max(var_logZ_inc, 0)

    w <- w * u
    w <- w / sum(w)
    ess_frac <- 1 / sum(w * w) / M
    if (verbose) cat(sprintf("Δ=%.4f | ESS/N=%.3f | logZ+=%.4f\n", delta, ess_frac, logZ_inc))

    # 4) Resample if needed
    resampled <- FALSE
    if (ess_frac < resample_threshold) {
      # 5) Pre-resample light rejuvenation
      if (ess_frac > 0.25 && lambda >= 0.20) {
        if (verbose) cat("  Pre-resample light rejuvenation (1 sweep)...\n")
        # light mixtures for pre-resample
        elite_mix_phi_pre <- tryCatch({
          mix <- fit_elite_mixture_Z(phi, w, elite_quantile = elite_quantile, G = G_mix,
                                     cov_inflation = cov_inflation, housekeeping = TRUE,
                                     min_G_keep = 2, merge_thresh = 0.10,
                                     verbose = FALSE, lambda = lambda_new)
          if (!.is_empty_mix(mix)) blend_mixes(mix, .default_std_normal_mix(ncol(phi)), eps = base_mix_eps) else mix
        }, error = function(e) NULL)
        hist_mix_pre <- NULL
        if (hist_mix_enable) {
          hist_mix_pre <- tryCatch(build_hist_mixture(elite_history_phi, lambda_new,
                                                     hist_mix_lambda_thresh, hist_mix_prob,
                                                     min_G_keep = 2, merge_thresh = 0.10,
                                                     min_eig = 3e-3, force = FALSE), error = function(e) NULL)
          if (!is.null(hist_mix_pre) && !is.null(hist_mix_pre$mix) && !.is_empty_mix(hist_mix_pre$mix)) {
            hist_mix_pre$mix <- blend_mixes(hist_mix_pre$mix, .default_std_normal_mix(ncol(phi)), eps = 0.05)
          }
        }
        # weak dims for RW expansion in this sweep
        weak_idx_pre <- weak_dims_from_mat(phi, w, frac = 0.25, min_keep = 1L)
        pre_move <- .rejuvenate_rw_batch(phi, loglik, logprior, w, lambda_new,
                                         caches, b_state, logprior_phi,
                                         log_prior_theta_given_phi_mat, gaussian_map_fn,
                                         pm_mode = pm_mode,
                                         n_moves = 1L, rw_scale = 0.35,
                                         weak_dim_idx = weak_idx_pre, rw_expand_factor = 1.8,
                                         elite_mix = elite_mix_phi_pre,
                                         hist_mix = hist_mix_pre,
                                         indep_prob = min(0.35, indep_prob),
                                         indep_t_df = indep_t_df, indep_t_prob = indep_t_prob,
                                         rng_seed = seed + 19L * round)
        phi      <- pre_move$phi
        loglik   <- pre_move$loglik
        logprior <- pre_move$logprior
        b_state  <- pre_move$b_state
      }
      idx <- if (resampling == "multinomial") {
        .resample_multinomial(w)
      } else {
        # Hilbert-sorted stratified resampling in whitened φ-space
        set.seed(seed + 991 * round)
        # whiten φ → Zφ using weighted mean/cov
        mu_phi <- colSums(phi * w)
        Xc <- sweep(phi, 2L, mu_phi, `-`)
        Sphi <- weighted_cov(phi, w)
        Sphi <- as.matrix(Matrix::nearPD(Sphi, conv.tol = 1e-7)$mat)
        diag(Sphi) <- pmax(diag(Sphi), 1e-10)
        R <- tryCatch(chol(Sphi), error = function(e) chol(Sphi + diag(1e-8, ncol(Sphi))))
        Zphi <- t(backsolve(R, t(Xc), transpose = TRUE))
        ord_h <- hilbert_sort_order(Zphi, bits = 16L)
        ws <- w[ord_h]
        sel_sorted <- stratified_resample_sorted(ws, deterministic = FALSE)
        ord_inv <- ord_h
        idx <- ord_inv[sel_sorted]
      }
      phi      <- phi[idx, , drop = FALSE]
      loglik   <- loglik[idx]
      logprior <- logprior[idx]
      if (is.matrix(b_state)) b_state <- b_state[idx, , drop = FALSE]
      w <- rep(1/M, M)
      resampled <- TRUE
      if (verbose) cat(sprintf("  Resampled (%s).\n", resampling))

      # Optional: refresh a small fraction of batch indices after resampling
      if (refresh_batches_after_resample) {
        if (is.matrix(b_state)) {
          ref <- .refresh_aux_state_block_mh(
            phi = phi, loglik = loglik, b_state = b_state, lambda = lambda_new,
            caches = caches,
            log_prior_theta_given_phi_mat = log_prior_theta_given_phi_mat,
            gaussian_map_fn = gaussian_map_fn,
            frac = refresh_batches_frac,
            rng_seed = seed + 1229L * round
          )
          b_state <- ref$b_state
          loglik <- ref$loglik
          aux_refresh_acc_hist <- c(aux_refresh_acc_hist, ref$acc_rate)
        } else {
          K_all <- caches[[1]]$K
          b_state <- refresh_batch_indices(b_state, K = K_all, frac = refresh_batches_frac)
          loglik <- .outer_loglik_hat_particles(
            Phi = phi, caches = caches, b_state = b_state,
            log_prior_theta_given_phi_mat = log_prior_theta_given_phi_mat,
            gaussian_map_fn = gaussian_map_fn,
            pm_mode = pm_mode
          )
        }
      }
    }

    # 5) Diagnostics at a probe φ (batch-aware; no model calls)
    if (diag_enable) {
      probe_phi <- .choose_probe_phi(phi, w, logprior, loglik, lambda_new, diag_probe)
      probe_b <- .probe_aux_row(phi, w, logprior, loglik, lambda_new, b_state)
      is_df <- .is_diagnostics_batch(caches, probe_phi, probe_b, log_prior_theta_given_phi_mat)
      ess_norm <- is_df$ess / is_df$M
      psis_k <- is_df$psis_k
      cat(sprintf("  [Diag] median ESS_i/M=%.2f | max k=%.2f | any ESS_i/M<0.10: %s\n",
                  stats::median(ess_norm, na.rm = TRUE),
                  if (all(is.na(psis_k))) NA_real_ else max(psis_k, na.rm = TRUE),
                  any(ess_norm < 0.10, na.rm = TRUE)))
      if (diag_replicates > 0L) {
        rsd <- .rep_logZi_sd_over_batches(caches, probe_phi, B_reps = diag_replicates,
                                          log_prior_theta_given_phi_mat = log_prior_theta_given_phi_mat,
                                          gaussian_map_fn = gaussian_map_fn,
                                          rng_seed = seed + 1777 * round)
        cat(sprintf("  [Diag] median SD_over_batches(logZ_i): %.3f\n",
                    stats::median(rsd, na.rm = TRUE)))
        # --- NEW: variance budget log (Fix #4)
        sum_sd2 <- sum(rsd[is.finite(rsd)]^2)
        cat(sprintf("  [VarBudget] sum_i SD_over_batches(logZ_i)^2 = %.3f (target ≈ 1–2)\n", sum_sd2))
        diag_sum_sd2 <- sum_sd2
      }
      # --- NEW: print top offenders (Fix #4)
      bad_idx <- which( (ess_norm < auto_enrich_ess_thresh) |
                        (!is.na(psis_k) & (psis_k > auto_enrich_k_thresh)) )
      if (length(bad_idx)) {
        ord <- order(ifelse(is.na(psis_k[bad_idx]), -Inf, psis_k[bad_idx]),
                     ess_norm[bad_idx], decreasing = c(TRUE, FALSE))
        show <- bad_idx[ord][seq_len(min(length(bad_idx), offenders_print_cap))]
        msg <- paste(sprintf("#%d ESS/M=%.2f k=%.2f", show, ess_norm[show], psis_k[show]), collapse = " | ")
        cat("  [Diag] Offenders:", msg, "\n")
      }
      diag_median_ess_norm <- stats::median(ess_norm, na.rm = TRUE)
      diag_max_psis_k <- if (all(is.na(psis_k))) NA_real_ else max(psis_k, na.rm = TRUE)
      diag_n_offenders <- length(bad_idx)
    }

    # Build elite mixture on φ for independence proposals (with tiny base mass)
    elite_mix_phi <- tryCatch({
      mix <- fit_elite_mixture_Z(phi, w, elite_quantile = elite_quantile, G = G_mix,
                                 cov_inflation = cov_inflation, housekeeping = TRUE,
                                 min_G_keep = 2, merge_thresh = 0.10,
                                 verbose = FALSE, lambda = lambda_new)
      if (!.is_empty_mix(mix)) blend_mixes(mix, .default_std_normal_mix(ncol(phi)), eps = base_mix_eps) else mix
    }, error = function(e) NULL)
    hist_mix <- NULL
    if (hist_mix_enable) {
      hist_mix <- tryCatch(build_hist_mixture(elite_history_phi, lambda_new,
                                              hist_mix_lambda_thresh, hist_mix_prob,
                                              min_G_keep = 2, merge_thresh = 0.10,
                                              min_eig = 3e-3, force = FALSE), error = function(e) NULL)
      if (!is.null(hist_mix) && !is.null(hist_mix$mix) && !.is_empty_mix(hist_mix$mix)) {
        # light base blend for history as well
        hist_mix$mix <- blend_mixes(hist_mix$mix, .default_std_normal_mix(ncol(phi)), eps = 0.05)
      }
    }

    # Adaptive number of moves (Δλ, difficulty, acceptance-aware)
    pred_rcess_full <- rCESS_stat(w, h, 1 - lambda)
    difficulty <- 1 - pred_rcess_full
    sDelta <- .clamp(delta / 0.05, 0.6, 1.4)
    sDiff  <- 1 + .clamp(1.5 * difficulty, 0, 1.5)
    acc_proxy <- mean(c(last_rw_acc, last_id_acc), na.rm = TRUE)
    if (!is.finite(acc_proxy)) acc_proxy <- 0.25
    sAcc <- if (acc_proxy < 0.10) 0.80 else if (acc_proxy < 0.20) 0.90 else 1.00
    base_moves <- n_moves
    floor_moves <- if (ess_frac < 0.40) 4L else if (ess_frac < 0.70) 3L else 2L
    n_moves_eff <- as.integer(ceiling(base_moves * sDelta * sDiff * sAcc))
    n_moves_eff <- max(floor_moves, min(n_moves_eff, moves_cap))
    if (lambda > 0.90) n_moves_eff <- max(n_moves_eff, 3L)
    if (lambda < 0.20) n_moves_eff <- min(n_moves_eff, 2L)

    # Weak-dimension inflation indices in φ-space
    weak_idx <- weak_dims_from_mat(phi, w, frac = 0.25, min_keep = 1L)

    # 6) Rejuvenate (MH RW + independence; correlated PM via fixed b_idx inside the sweep)
    rw_scale <- exp(log_rw_scale)
    move <- .rejuvenate_rw_batch(phi, loglik, logprior, w, lambda_new,
                                 caches, b_state, logprior_phi,
                                 log_prior_theta_given_phi_mat, gaussian_map_fn,
                                 pm_mode = pm_mode,
                                 n_moves = n_moves_eff,
                                 rw_scale = rw_scale,
                                 weak_dim_idx = weak_idx,
                                 rw_expand_factor = 2.5,
                                 elite_mix = elite_mix_phi,
                                 hist_mix = hist_mix,
                                 indep_prob = indep_prob,
                                 indep_t_df = indep_t_df,
                                 indep_t_prob = indep_t_prob,
                                 rng_seed = seed + 97L * round)
    phi      <- move$phi
    loglik   <- move$loglik
    logprior <- move$logprior
    b_state  <- move$b_state
    acc_hist <- c(acc_hist, move$acc_rate)
    # Adapt RW scale
    log_rw_scale <- .clamp(log_rw_scale + rm_gain * (move$acc_rate - target_acc),
                           log(0.05), log(2.5))
    if (verbose) {
      cat(sprintf("  MH acc=%.3f | rw=%.3f id=%.3f | rw_scale=%.3f\n",
                  move$acc_rate, move$rw_accept_rate %||% NA_real_,
                  move$indep_accept_rate %||% NA_real_, exp(log_rw_scale)))
    }

    last_rw_acc <- move$rw_accept_rate
    last_id_acc <- move$indep_accept_rate

    # Strict-mode auxiliary-state block refresh (MH-correct at fixed λ)
    if (is.matrix(b_state) && block_refresh_every > 0L &&
        (round %% as.integer(max(1L, block_refresh_every)) == 0L)) {
      ref <- .refresh_aux_state_block_mh(
        phi = phi, loglik = loglik, b_state = b_state, lambda = lambda_new,
        caches = caches,
        log_prior_theta_given_phi_mat = log_prior_theta_given_phi_mat,
        gaussian_map_fn = gaussian_map_fn,
        frac = block_refresh_frac,
        rng_seed = seed + 7001L * round
      )
      b_state <- ref$b_state
      loglik <- ref$loglik
      aux_refresh_acc_hist <- c(aux_refresh_acc_hist, ref$acc_rate)
      if (verbose) cat(sprintf("  Aux refresh acc=%.3f\n", ref$acc_rate))
    }

    # update history store (keep last 6)
    if (!is.null(elite_mix_phi) && !.is_empty_mix(elite_mix_phi)) last_elite_mix_phi <- elite_mix_phi
    if (lambda_new > 0.5 && !is.null(last_elite_mix_phi) && !.is_empty_mix(last_elite_mix_phi)) {
      elite_history_phi[[length(elite_history_phi) + 1L]] <- last_elite_mix_phi
      if (length(elite_history_phi) > 6L) elite_history_phi <- tail(elite_history_phi, 6L)
    }

    # 7) NEW: targeted enrichment of worst subjects (Fixes #2–#5)
    if (auto_enrich_enable && lambda_new >= auto_enrich_lambda_gate && (round %% auto_enrich_every == 0)) {
      # use the same probe φ as diagnostics (wmean/best), and build a small panel
      probe_phi <- .choose_probe_phi(phi, w, logprior, loglik, lambda_new, diag_probe)
      probe_phi_wmean <- .choose_probe_phi(phi, w, logprior, loglik, lambda_new, "wmean")
      probe_phi_best  <- .choose_probe_phi(phi, w, logprior, loglik, lambda_new, "best")
      probe_b <- .probe_aux_row(phi, w, logprior, loglik, lambda_new, b_state)
      is_df <- .is_diagnostics_batch(caches, probe_phi, probe_b, log_prior_theta_given_phi_mat)
      ess_norm <- is_df$ess / is_df$M
      psis_k <- is_df$psis_k
      bad_idx <- which( (ess_norm < auto_enrich_ess_thresh) |
                        (!is.na(psis_k) & (psis_k > auto_enrich_k_thresh)) )
      if (length(bad_idx)) {
        ord <- order(ifelse(is.na(psis_k[bad_idx]), -Inf, psis_k[bad_idx]),
                     ess_norm[bad_idx], decreasing = c(TRUE, FALSE))
        select <- bad_idx[ord][seq_len(min(length(bad_idx), auto_enrich_max_units))]
        cat("  [Enrich] subjects:", paste(select, collapse = ","), "\n")
        for (ii in select) {
          if (is.null(data_list[[ii]])) {
            warning(sprintf("[Enrich] data_list[[%d]] is NULL; skipping subject %d", ii, ii))
            next
          }
          # ### PATCH: surgical enrichment on ACTIVE batch only with gating/rollback
          caches[[ii]] <- enrich_subject_cache_with_anchor_surgical(
            cache = caches[[ii]],
            phi_probe = probe_phi,
            b = probe_b[ii],
            log_prior_theta_given_phi_mat = log_prior_theta_given_phi_mat,
            gaussian_map_fn = if (!is.null(gaussian_map_fn)) function(phi, d) gaussian_map_fn(phi, d) else NULL,
            data = data_list[[ii]],
            loglik_fn = loglik_fn,
            eps_prior_base = eps_prior_anchor,
            weak_inflate_factor = weak_inflate_factor,
            elite_q = 0.4, G = 12L, cov_infl = 4.0,
            replace_frac = NULL,                      # adaptive replacement fraction
            sobol_seed = caches[[ii]]$batches[[probe_b[ii]]]$sobol_seed + 1L,  # small shift only
            n_cores = n_cores_inner,
            verbose = verbose,
            subj_id = ii,
            phi_panel = rbind(probe_phi_wmean, probe_phi_best)
          )
        }
        # Aux state unchanged; proposal correlation preserved.
      }
    }

    loglik <- .outer_loglik_hat_particles(
      Phi = phi, caches = caches, b_state = b_state,
      log_prior_theta_given_phi_mat = log_prior_theta_given_phi_mat,
      gaussian_map_fn = gaussian_map_fn,
      pm_mode = pm_mode
    )

    # 8) Optional light batch-index refresh at the end of the round
    if (refresh_batches_each_round) {
      if (is.matrix(b_state)) {
        ref <- .refresh_aux_state_block_mh(
          phi = phi, loglik = loglik, b_state = b_state, lambda = lambda_new,
          caches = caches,
          log_prior_theta_given_phi_mat = log_prior_theta_given_phi_mat,
          gaussian_map_fn = gaussian_map_fn,
          frac = refresh_batches_frac,
          rng_seed = seed + 9929L * round
        )
        b_state <- ref$b_state
        loglik <- ref$loglik
        aux_refresh_acc_hist <- c(aux_refresh_acc_hist, ref$acc_rate)
      } else {
        K_all <- caches[[1]]$K
        b_state <- refresh_batch_indices(b_state, K = K_all, frac = refresh_batches_frac)
        loglik <- .outer_loglik_hat_particles(
          Phi = phi, caches = caches, b_state = b_state,
          log_prior_theta_given_phi_mat = log_prior_theta_given_phi_mat,
          gaussian_map_fn = gaussian_map_fn,
          pm_mode = pm_mode
        )
      }
    }

    # Advance annealing
    lambda <- lambda_new
    lambda_hist <- c(lambda_hist, lambda)

    if (isTRUE(collect_round_diagnostics)) {
      rd <- list(
        round = round,
        elapsed_sec = as.numeric(proc.time()[3] - t_round_start),
        lambda = lambda,
        delta = delta,
        ess_frac = ess_frac,
        logZ_inc = logZ_inc,
        mh_acc = move$acc_rate,
        rw_acc = move$rw_accept_rate %||% NA_real_,
        indep_acc = move$indep_accept_rate %||% NA_real_,
        diag_median_ess_norm = diag_median_ess_norm,
        diag_max_psis_k = diag_max_psis_k,
        diag_sum_sd2 = diag_sum_sd2,
        diag_n_offenders = diag_n_offenders
      )
      diag_list[[length(diag_list) + 1L]] <- rd
    }
  }

  list(
    phi = phi, w = w,
    loglik = loglik, logprior = logprior,
    log_evidence = logZ,
    mcse_log_evidence = sqrt(mcse_var_accum),
    logZ_increments = logZ_increments,
    b_idx = if (!is.matrix(b_state)) b_state else NULL,
    B_idx = if (is.matrix(b_state)) b_state else NULL,
    meta = list(
      pm_mode = pm_mode,
      lambda_hist = lambda_hist,
      acc_hist = acc_hist,
      rw_scale_final = exp(log_rw_scale),
      resampling = resampling,
      refresh_batches_after_resample = refresh_batches_after_resample,
      refresh_batches_each_round = refresh_batches_each_round,
      refresh_batches_frac = refresh_batches_frac,
      block_refresh_every = block_refresh_every,
      block_refresh_frac = block_refresh_frac,
      aux_refresh_acc_hist = aux_refresh_acc_hist,
      n_cores_inner = n_cores_inner,
      elapsed_sec = as.numeric(proc.time()[3] - t_start),
      round_diagnostics = diag_list
    )
  )
}

# ---------------------- Generic outer SMC (direct likelihood) -------------------
.outer_eval_logprior_one <- function(theta_row, logprior_fn) {
  val <- tryCatch(logprior_fn(theta_row), error = function(e) NULL)
  if (is.null(val)) {
    val <- tryCatch(
      logprior_fn(matrix(theta_row, nrow = 1L)),
      error = function(e) stop("logprior failed for both vector and 1-row matrix input.")
    )
  }
  val <- as.numeric(val)
  if (!length(val)) stop("logprior returned empty output for one-row input.")
  val[1L]
}

.outer_eval_logprior_mat <- function(Theta, logprior_fn) {
  out <- tryCatch(logprior_fn(Theta), error = function(e) NULL)
  if (!is.null(out)) {
    out <- as.numeric(out)
    if (length(out) == nrow(Theta)) return(out)
  }
  vapply(seq_len(nrow(Theta)), function(i) .outer_eval_logprior_one(Theta[i, , drop = TRUE], logprior_fn), numeric(1))
}

outer_smc_direct <- function(
    data,
    rprior,                  # function(n) -> matrix(n x d)
    logprior,                # function(theta_matrix_or_row) -> vector/scalar
    loglik_fn,               # function(theta_matrix_or_row, data) -> vector/scalar
    M = 2000L,
    cess_target = 0.95,
    resample_threshold = 0.5,
    n_moves = 2L,
    rw_scale_init = 1.0,
    max_rounds = 200L,
    n_cores_loglik = 1L,
    seed = 123,
    verbose = TRUE
) {
  set.seed(seed)

  Theta <- rprior(M)
  if (!is.matrix(Theta)) Theta <- matrix(Theta, nrow = M)
  d <- ncol(Theta)

  logprior_curr <- .outer_eval_logprior_mat(Theta, logprior)
  loglik_curr <- ll_parallel(Theta, data, loglik_fn, n_cores = n_cores_loglik)
  if (!all(is.finite(logprior_curr))) stop("Non-finite log prior at initialization.")
  if (!any(is.finite(loglik_curr))) stop("All initial log-likelihood values are non-finite.")
  loglik_curr[!is.finite(loglik_curr)] <- -Inf

  lambda <- 0.0
  w <- rep(1 / M, M)
  logZ <- 0.0

  lambda_hist <- lambda
  ess_hist <- numeric()
  acc_hist <- numeric()
  logZ_increments <- numeric()
  log_rw_scale <- log(rw_scale_init)

  round <- 0L
  while (lambda < 1 - 1e-12 && round < max_rounds) {
    round <- round + 1L
    lambda_new <- next_lambda_via_rCESS_stat(w, loglik_curr, lambda, target = cess_target)
    delta <- lambda_new - lambda
    if (!is.finite(delta) || delta <= 1e-8) {
      delta <- min(1 - lambda, 1e-4)
      lambda_new <- lambda + delta
    }

    m <- max(loglik_curr)
    log_u <- delta * (loglik_curr - m)
    logZ_inc <- logsumexp_w(log_u, w) + delta * m
    if (!is.finite(logZ_inc)) stop("Non-finite log-evidence increment.")
    logZ <- logZ + logZ_inc
    logZ_increments <- c(logZ_increments, logZ_inc)

    w <- w * exp(log_u)
    sw <- sum(w)
    if (!is.finite(sw) || sw <= 0) stop("Particle weights became invalid.")
    w <- w / sw
    ess_frac <- (1 / sum(w * w)) / M
    ess_hist <- c(ess_hist, ess_frac)

    if (verbose) {
      cat(sprintf("[Outer Round %d] lambda %.4f -> %.4f | ESS/N=%.3f | logZ+=%.4f\n",
                  round, lambda, lambda_new, ess_frac, logZ_inc))
    }

    if (ess_frac < resample_threshold) {
      idx <- stratified_resample_sorted(w, deterministic = FALSE)
      Theta <- Theta[idx, , drop = FALSE]
      logprior_curr <- logprior_curr[idx]
      loglik_curr <- loglik_curr[idx]
      w <- rep(1 / M, M)
    }

    rw_scale <- exp(log_rw_scale)
    S <- tryCatch(weighted_cov(Theta, w), error = function(e) stats::cov(Theta))
    S <- as.matrix(S)
    diag(S) <- pmax(diag(S), 1e-8)
    L <- tryCatch(
      chol((rw_scale^2 / max(d, 1L)) * S + diag(1e-10, d)),
      error = function(e) chol((rw_scale^2 / max(d, 1L)) * (S + diag(1e-6, d)))
    )

    acc_total <- 0L
    for (mv in seq_len(n_moves)) {
      prop <- Theta + matrix(rnorm(M * d), M, d) %*% L
      lp_prop <- .outer_eval_logprior_mat(prop, logprior)
      ll_prop <- ll_parallel(prop, data, loglik_fn, n_cores = n_cores_loglik)
      loga <- (lp_prop - logprior_curr) + lambda_new * (ll_prop - loglik_curr)
      loga[!is.finite(loga)] <- -Inf
      acc <- log(runif(M)) < pmin(0, loga)
      if (any(acc)) {
        Theta[acc, ] <- prop[acc, , drop = FALSE]
        logprior_curr[acc] <- lp_prop[acc]
        loglik_curr[acc] <- ll_prop[acc]
      }
      acc_total <- acc_total + sum(acc)
    }

    acc_rate <- acc_total / max(1L, M * n_moves)
    acc_hist <- c(acc_hist, acc_rate)
    log_rw_scale <- .clamp(log_rw_scale + 0.05 * (acc_rate - 0.234), log(0.02), log(3.0))

    lambda <- lambda_new
    lambda_hist <- c(lambda_hist, lambda)
  }

  list(
    theta = Theta,
    w = w,
    loglik = loglik_curr,
    logprior = logprior_curr,
    log_evidence = logZ,
    logZ_increments = logZ_increments,
    meta = list(
      lambda_hist = lambda_hist,
      ess_hist = ess_hist,
      acc_hist = acc_hist,
      rw_scale_final = exp(log_rw_scale),
      rounds = round,
      final_lambda = lambda
    )
  )
}
    diag_median_ess_norm <- NA_real_
    diag_max_psis_k <- NA_real_
    diag_sum_sd2 <- NA_real_
    diag_n_offenders <- 0L
