# ======================================================================
# PMMH with Correlated Pseudo-Marginal + Delayed Acceptance + Refresh
# Outer model: per-coordinate Normal(μ_j, τ_j^2), with μ_j ~ N(m0_j, s0_j^2),
# τ_j^2 ~ InvGamma(a0_j, b0_j). Caches provide θ-samples under base prior p0.
# ======================================================================

# ---------- Small utilities ----------
.logsumexp <- function(x) {
  m <- max(x); m + log(sum(exp(x - m)))
}
igamma_logpdf <- function(x, a, b) { # InvGamma(a, b) with "rate" b
  # density: p(x) = b^a / Γ(a) * x^(-a-1) * exp(-b/x), x>0
  out <- a*log(b) - lgamma(a) - (a+1)*log(x) - b/x
  out[x <= 0] <- -Inf
  out
}

# 1) Strictly-increasing cumulative weights (last exactly 1)
.make_cumw <- function(w) {
  w <- pmax(w, 0); cw <- cumsum(w / sum(w))
  eps <- .Machine$double.eps * 8
  cw <- cw + eps * seq_along(cw)
  cw <- cw / cw[length(cw)]; cw[length(cw)] <- 1.0
  cw
}
sample_idx_from_uniforms <- function(cumw, u) {
  K <- length(cumw)
  tiny <- .Machine$double.eps * 8
  u <- pmin(1 - tiny, pmax(u, tiny))
  idx <- 1L + findInterval(u, cumw, rightmost.closed = TRUE)
  idx[idx < 1L] <- 1L; idx[idx > K] <- K
  idx
}

# 3) Robust row-wise normal product (unchanged logic, ensures matrix)
log_prod_normal_rows <- function(Theta, mu, tau2) {
  Theta <- as.matrix(Theta)
  # If names are present, enforce column order to match mu
  if (!is.null(colnames(Theta)) && !is.null(names(mu))) {
    if (!identical(colnames(Theta), names(mu))) {
      if (setequal(colnames(Theta), names(mu))) {
        Theta <- Theta[, names(mu), drop = FALSE]
      } else {
        stop("Name mismatch: cache columns and 'mu' names differ.")
      }
    }
  }
  # Align tau2 by names if provided
  if (!is.null(names(tau2)) && !is.null(names(mu)) && !identical(names(tau2), names(mu))) {
    if (setequal(names(tau2), names(mu))) {
      tau2 <- tau2[names(mu)]
    } else {
      stop("Name mismatch between mu and tau2.")
    }
  }
  d <- ncol(Theta)
  if (length(mu) != d || length(tau2) != d)
    stop(sprintf("Dimension mismatch: ncol(Theta)=%d, length(mu)=%d, length(tau2)=%d", d, length(mu), length(tau2)))
  c0 <- -0.5 * (sum(log(tau2)) + d * log(2*pi))
  quad <- rowSums((sweep(Theta, 2L, mu, `-`)^2) / matrix(tau2, nrow(Theta), d, byrow = TRUE))
  c0 - 0.5 * quad
}
# ---------- Per-subject unbiased Li(η) estimator using a subsample ----------
# Uses: with-replacement subsample indices generated from CPM uniforms
# Returns log Zi_hat(η) and (optionally) the ESS fraction under η using FULL cache (for refresh decision)
logZi_hat_from_cache <- function(cache, mu, tau2, idx_sub, need_full_ess = FALSE) {
  Theta <- as.matrix(cache$theta)
  K <- nrow(Theta)
  # Last-resort guard: clamp indices to [1, K] instead of aborting (can happen if a cache was refreshed)
  if (any(!is.finite(idx_sub) | idx_sub < 1L | idx_sub > K)) {
    idx_sub <- pmin.int(pmax.int(idx_sub, 1L), K)
  }
  d_sub <- log_prod_normal_rows(Theta[idx_sub, , drop = FALSE], mu, tau2) - cache$log_p0[idx_sub]
  lmeanexp <- .logsumexp(d_sub) - log(length(idx_sub))
  logZi <- cache$logZ0 + lmeanexp

  if (!need_full_ess) return(list(logZi = logZi, ess_frac = NA_real_))
  d_all <- log_prod_normal_rows(Theta, mu, tau2) - cache$log_p0
  lw <- log(pmax(cache$w, .Machine$double.eps)) + d_all
  lw <- lw - .logsumexp(lw)
  w_norm <- exp(lw)
  ess <- 1 / sum(w_norm^2)
  list(logZi = logZi, ess_frac = ess / K)
}


# ---------- Correlated seed updater ----------
# Maintain per-subject Gaussian seeds z ~ N(0, I), map via Φ to uniforms,
# build indices via inverse CDF on cache weights.
update_cpm_indices <- function(z_curr, rho) {
  # z':= rho z + sqrt(1-rho^2) ε, ε ~ N(0, I)
  z_prop <- rho * z_curr + sqrt(max(1 - rho^2, 0)) * rnorm(length(z_curr))
  u_prop <- pnorm(z_prop)
  list(z_prop = z_prop, u_prop = u_prop)
}

# ---------- Outer prior on η = (μ, log τ^2) ----------
log_prior_eta <- function(mu, log_tau2, m0, s0, a0, b0) {
  tau2 <- exp(log_tau2)
  # sum_j [ log N(mu_j | m0_j, s0_j^2) + log IG(tau2_j | a0_j, b0_j) ]
  lp_mu <- sum(dnorm(mu, mean = m0, sd = s0, log = TRUE))
  lp_t2 <- sum(igamma_logpdf(tau2, a0, b0) + log_tau2)  # add log|d tau2 / d log_tau2| = log_tau2
  lp_mu + lp_t2
}

# ---------- Cache-based posterior moments for proposal ----------
# For a given cache and outer parameters (mu, tau2), compute
# normalized importance weights w_k ∝ cache$w * N(θ_k | mu, diag(tau2)) / p0(θ_k),
# and return posterior mean and diagonal variance under those weights.
cache_posterior_moments <- function(cache, mu, tau2) {
  Theta <- as.matrix(cache$theta)
  lw <- log(pmax(cache$w, .Machine$double.eps)) +
        (log_prod_normal_rows(Theta, mu, tau2) - cache$log_p0)
  lw <- lw - .logsumexp(lw)
  wnorm <- exp(lw)
  m <- colSums(Theta * wnorm)  # d-vector
  centered <- sweep(Theta, 2L, m, "-")
  vdiag <- colSums(centered^2 * wnorm)  # d-vector (diagonal only, robust)
  list(mean = m, vdiag = vdiag)
}

# Build simple sufficient statistics from list of moments
.build_mom_suff <- function(mom_list) {
  means_mat <- do.call(rbind, lapply(mom_list, `[[`, "mean"))   # N x d
  vdiag_mat <- do.call(rbind, lapply(mom_list, `[[`, "vdiag"))  # N x d
  list(
    means = means_mat,
    vdiag = vdiag_mat,
    sum_means = colSums(means_mat),
    N = nrow(means_mat)
  )
}

# Sample a "Gibbs-ish" independence proposal:
# (1) mu | tau2_cond ~ Normal(m_star, s2_star)
# (2) tau2 | mu_prop ~ InvGamma(a_star, b_star)
.sample_gibbsish <- function(mu_cond, log_tau2_cond, suff, m0, s0, a0, b0) {
  tau2_cond <- exp(log_tau2_cond)
  N <- suff$N
  d <- length(mu_cond)
  s2_star_inv <- (1 / (s0^2)) + (N / tau2_cond)
  s2_star <- 1 / s2_star_inv
  m_star <- s2_star * ((m0 / (s0^2)) + (suff$sum_means / tau2_cond))
  mu_prop <- rnorm(d, mean = m_star, sd = sqrt(s2_star)); names(mu_prop) <- names(mu_cond)

  mu_mat <- matrix(mu_prop, N, d, byrow = TRUE)
  S <- colSums((suff$means - mu_mat)^2 + suff$vdiag)
  a_star <- a0 + N / 2
  b_star <- b0 + 0.5 * S
  tau2_prop <- 1 / rgamma(d, shape = a_star, rate = b_star)
  log_tau2_prop <- log(tau2_prop); names(log_tau2_prop) <- names(mu_cond)
  list(mu_prop = mu_prop, log_tau2_prop = log_tau2_prop)
}

# Log proposal density for the same two-step mechanism, including Jacobian for log_tau2
.logq_gibbsish <- function(mu_draw, log_tau2_draw, mu_cond, log_tau2_cond, suff, m0, s0, a0, b0) {
  tau2_cond <- exp(log_tau2_cond)
  N <- suff$N
  d <- length(mu_draw)
  s2_star_inv <- (1 / (s0^2)) + (N / tau2_cond)
  s2_star <- 1 / s2_star_inv
  m_star <- s2_star * ((m0 / (s0^2)) + (suff$sum_means / tau2_cond))
  lq_mu <- sum(dnorm(mu_draw, mean = m_star, sd = sqrt(s2_star), log = TRUE))
  mu_mat <- matrix(mu_draw, N, d, byrow = TRUE)
  S <- colSums((suff$means - mu_mat)^2 + suff$vdiag)
  a_star <- a0 + N / 2
  b_star <- b0 + 0.5 * S
  tau2_draw <- exp(log_tau2_draw)
  lq_tau <- sum(igamma_logpdf(tau2_draw, a_star, b_star) + log_tau2_draw)  # + log|d tau2 / d log_tau2|
  lq_mu + lq_tau
}

# ---------- Selective refresh policy (ESS-based) ----------
maybe_refresh_subject <- function(i, cache, mu, tau2,
                                  ess_thresh = 0.25,
                                  cooldown = 20L,
                                  last_refresh_iter,
                                  iter,
                                  refresh_fn = NULL,
                                  verbose = FALSE) {
  # Avoid spamming refreshes
  if (is.null(refresh_fn)) return(list(cache = cache, refreshed = FALSE))
  if (!is.null(last_refresh_iter[[i]]) && (iter - last_refresh_iter[[i]] < cooldown)) {
    return(list(cache = cache, refreshed = FALSE))
  }
  ess_info <- logZi_hat_from_cache(cache, mu, tau2, idx_sub = 1L, need_full_ess = TRUE)
  if (is.finite(ess_info$ess_frac) && ess_info$ess_frac < ess_thresh) {
    if (isTRUE(verbose)) cat(sprintf("  [refresh] subject %d: ESS_frac=%.3f < %.2f -> refreshing cache at current η\n",
                                     i, ess_info$ess_frac, ess_thresh))
    # User-supplied: must rebuild cache under prior N(mu, diag(tau2))
    new_cache <- refresh_fn(i = i, mu = mu, tau2 = tau2, old_cache = cache)
    return(list(cache = new_cache, refreshed = TRUE))
  }
  list(cache = cache, refreshed = FALSE)
}

# ---------- Main PMMH with CPM + DA ----------
pmmh_outer_cpm_da <- function(
    caches,                            # list of per-subject caches (from build_inner_cache)
    m0, s0, a0, b0,                    # hyper-prior parameters (scalar or length d)
    n_iter = 5000L,
    init_mu = NULL, init_log_tau2 = NULL, # if NULL, use m0 and log of IG mean if available
    # mixture of proposal kernels
    indep_mix_prob = 0.20,             # probability to use the Gibbs-ish independence proposal
    # DA & CPM controls
    K1 = 16L,                          # per-subject subsample size for Stage-1
    K2 = 128L,                         # per-subject subsample size for Stage-2
    rho = 0.99,                        # CPM correlation in (0,1); 0.98–0.999 works well
    # proposal/adaptation
    target_acc = 0.25,
    adapt_start = 200L,
    cov_eps = 1e-8,
    # refresh policy
    ess_thresh = 0.25,
    ess_check_every = 10L,
    cooldown = 25L,
    refresh_fn = NULL,                 # function(i, mu, tau2, old_cache) -> new_cache
    max_refresh_per_check = 8L,
    # randomness
    seed = 1L,
    verbose = TRUE
) {
  set.seed(seed)
  N <- length(caches)
  d <- ncol(caches[[1]]$theta)
  # broadcast priors if scalars
  if (length(m0)  == 1L) m0  <- rep(m0,  d)
  if (length(s0)  == 1L) s0  <- rep(s0,  d)
  if (length(a0)  == 1L) a0  <- rep(a0,  d)
  if (length(b0)  == 1L) b0  <- rep(b0,  d)

  # init η
  if (is.null(init_mu))        init_mu <- m0
  if (is.null(init_log_tau2))  init_log_tau2 <- log((ifelse(a0>1, b0/(a0-1), 1)))  # IG mean if a0>1, else 1
  theta_dimnames <- colnames(caches[[1]]$theta)
  names(init_mu) <- theta_dimnames

  # ---- Enforce a master parameter order across all caches (by names) ----
  if (is.null(theta_dimnames)) stop("Cache columns must be named.")
  for (i in seq_len(N)) {
    cn <- colnames(caches[[i]]$theta)
    if (is.null(cn)) stop(sprintf("Cache %d has no column names.", i))
    if (!identical(cn, theta_dimnames)) {
      if (!setequal(cn, theta_dimnames)) {
        stop(sprintf("Cache %d names differ from master.\nMaster: %s\nThis:   %s",
                     i, paste(theta_dimnames, collapse=","), paste(cn, collapse=",")))
      }
      # Reorder to master
      caches[[i]]$theta <- caches[[i]]$theta[, theta_dimnames, drop = FALSE]
      # (log_p0 is a total prior log-density per row, independent of column order)
    }
  }

  # Precompute cumulative weights for inverse-CDF sampling
  cumw_list <- lapply(caches, function(ci) .make_cumw(ci$w))

  # Initialize CPM Gaussian seeds (per subject, per stage)
  z1_list <- replicate(N, rnorm(K1), simplify = FALSE)
  z2_list <- replicate(N, rnorm(K2), simplify = FALSE)

  # Helper to build per-subject index sets from seeds
  make_idx_sets <- function(z_list, cumw_list) {
    lapply(seq_along(z_list), function(i) {
      u <- pnorm(z_list[[i]])
      sample_idx_from_uniforms(cumw_list[[i]], u)
    })
  }

  idx1_list <- make_idx_sets(z1_list, cumw_list)
  idx2_list <- make_idx_sets(z2_list, cumw_list)

  # State
  mu_curr <- init_mu
  log_tau2_curr <- init_log_tau2
  tau2_curr <- exp(log_tau2_curr)
  # Stage-1 and Stage-2 log posterior at current state (compute once)
  lp_eta_curr <- log_prior_eta(mu_curr, log_tau2_curr, m0, s0, a0, b0)

  # compute Stage-1 (cheap)
  logLi1_curr <- 0
  for (i in seq_len(N)) {
    logLi1_curr <- logLi1_curr + logZi_hat_from_cache(caches[[i]], mu_curr, tau2_curr, idx1_list[[i]], FALSE)$logZi
  }
  lpost1_curr <- lp_eta_curr + logLi1_curr
  # Stage-2 exact(er)
  logLi2_curr <- 0
  for (i in seq_len(N)) {
    logLi2_curr <- logLi2_curr + logZi_hat_from_cache(caches[[i]], mu_curr, tau2_curr, idx2_list[[i]], FALSE)$logZi
  }
  lpost2_curr <- lp_eta_curr + logLi2_curr

  # Adaptation bookkeeping (Haario + RM scale)
  p_dim <- 2L * d
  x_curr <- c(mu_curr, log_tau2_curr)
  mu_adapt <- x_curr
  C_adapt <- diag(p_dim) * 1e-2
  scale <- (2.38^2 / p_dim)
  log_scale <- log(scale)
  step_t <- function(t) 1 / (t + 10)   # RM gain

  # Storage
  out_mu <- matrix(NA_real_, n_iter, d)
  out_log_tau2 <- matrix(NA_real_, n_iter, d)
  out_acc <- logical(n_iter)
  out_stage1_pass <- logical(n_iter)
  refresh_count <- integer(n_iter)

  # Refresh cooldown track
  last_refresh_iter <- vector("list", N)

  if (isTRUE(verbose)) cat(sprintf("PMMH (CPM+DA): N=%d subjects, d=%d dims, K1=%d, K2=%d, ρ=%.3f\n",
                                   N, d, K1, K2, rho))

  # ------------- main loop -------------
  for (it in seq_len(n_iter)) {
    # Re-map current-state indices from current seeds and (possibly refreshed) cumw.
    # This preserves CPM (seeds unchanged) but guarantees indices match latest cache sizes.
    idx1_list <- lapply(seq_len(N), function(i) sample_idx_from_uniforms(cumw_list[[i]], pnorm(z1_list[[i]])))
    idx2_list <- lapply(seq_len(N), function(i) sample_idx_from_uniforms(cumw_list[[i]], pnorm(z2_list[[i]])))

    # ---------- propose (μ, log τ^2) via adaptive RW ----------
    # Current proposal covariance
    Sig <- exp(log_scale) * (C_adapt + cov_eps * diag(p_dim))
    use_indep <- (runif(1) < indep_mix_prob)
    if (!use_indep) {
      # symmetric random-walk proposal
      prop <- as.numeric(mvtnorm::rmvnorm(1, mean = x_curr, sigma = Sig))
      mu_prop <- prop[seq_len(d)];         names(mu_prop) <- theta_dimnames
      log_tau2_prop <- prop[-seq_len(d)];  names(log_tau2_prop) <- theta_dimnames
      tau2_prop <- exp(log_tau2_prop)
      lp_eta_prop <- log_prior_eta(mu_prop, log_tau2_prop, m0, s0, a0, b0)
    } else {
      # Gibbs-ish independence proposal (two-step, state-dependent via moments)
      moms_curr <- lapply(caches, cache_posterior_moments, mu = mu_curr, tau2 = tau2_curr)
      suff_curr <- .build_mom_suff(moms_curr)
      indep_draw <- .sample_gibbsish(mu_curr, log_tau2_curr, suff_curr, m0, s0, a0, b0)
      mu_prop <- indep_draw$mu_prop; log_tau2_prop <- indep_draw$log_tau2_prop
      tau2_prop <- exp(log_tau2_prop)
      lp_eta_prop <- log_prior_eta(mu_prop, log_tau2_prop, m0, s0, a0, b0)
    }

    # ---------- CPM seeds for proposal ----------
    # Stage-1
    idx1_prop <- vector("list", N)
    z1_prop_list <- vector("list", N)
    for (i in seq_len(N)) {
      up <- update_cpm_indices(z1_list[[i]], rho)
      z1_prop_list[[i]] <- up$z_prop
      idx1_prop[[i]] <- sample_idx_from_uniforms(cumw_list[[i]], up$u_prop)
    }
    # Stage-1 log-posterior
    logLi1_prop <- 0
    for (i in seq_len(N)) {
      logLi1_prop <- logLi1_prop + logZi_hat_from_cache(caches[[i]], mu_prop, tau2_prop, idx1_prop[[i]], FALSE)$logZi
    }
    lpost1_prop <- lp_eta_prop + logLi1_prop

    # ---------- DA Stage-1 decision ----------
    if (!use_indep) {
      # symmetric RW: proposal terms cancel
      log_r1 <- lpost1_prop - lpost1_curr
    } else {
      # include proposal densities q(curr|prop) - q(prop|curr)
      # forward q(prop | curr) based on moments at current state
      lq_prop_given_curr <- .logq_gibbsish(mu_prop, log_tau2_prop,
                                           mu_curr, log_tau2_curr, suff_curr,
                                           m0, s0, a0, b0)
      # reverse q(curr | prop) based on moments at proposed state
      moms_prop <- lapply(caches, cache_posterior_moments, mu = mu_prop, tau2 = tau2_prop)
      suff_prop <- .build_mom_suff(moms_prop)
      lq_curr_given_prop <- .logq_gibbsish(mu_curr, log_tau2_curr,
                                           mu_prop, log_tau2_prop, suff_prop,
                                           m0, s0, a0, b0)
      log_r1 <- (lpost1_prop - lpost1_curr) + (lq_curr_given_prop - lq_prop_given_curr)
    }
    pass1 <- (log(runif(1)) < min(0, log_r1))
    out_stage1_pass[it] <- pass1

    accepted <- FALSE
    if (pass1) {
      # Stage-2 CPM seeds (correlated)
      idx2_prop <- vector("list", N)
      z2_prop_list <- vector("list", N)
      for (i in seq_len(N)) {
        up2 <- update_cpm_indices(z2_list[[i]], rho)
        z2_prop_list[[i]] <- up2$z_prop
        idx2_prop[[i]] <- sample_idx_from_uniforms(cumw_list[[i]], up2$u_prop)
      }
      # Stage-2 log-posterior
      logLi2_prop <- 0
      for (i in seq_len(N)) {
        logLi2_prop <- logLi2_prop + logZi_hat_from_cache(caches[[i]], mu_prop, tau2_prop, idx2_prop[[i]], FALSE)$logZi
      }
      lpost2_prop <- lp_eta_prop + logLi2_prop

      # DA Stage-2 correction
      # Proposal terms cancel in the DA correction ratio, regardless of symmetry.
      # Only the refinement from Stage-1 to Stage-2 enters:
      # r2 = [pi2(prop)/pi1(prop)] / [pi2(curr)/pi1(curr)]
      # hence:
      log_r2 <- (lpost2_prop - lpost2_curr) - (lpost1_prop - lpost1_curr)

      accepted <- (log(runif(1)) < min(0, log_r2))

      if (accepted) {
        # accept full move; switch state & seeds
        mu_curr <- mu_prop
        log_tau2_curr <- log_tau2_prop
        tau2_curr <- tau2_prop
        lp_eta_curr <- lp_eta_prop
        lpost1_curr <- lpost1_prop
        lpost2_curr <- lpost2_prop
        z1_list <- z1_prop_list
        z2_list <- z2_prop_list
        idx1_list <- idx1_prop
        idx2_list <- idx2_prop
        x_curr <- c(mu_curr, log_tau2_curr)
        names(mu_curr) <- theta_dimnames
        names(log_tau2_curr) <- theta_dimnames
      }
    } # end Stage-2

    out_acc[it] <- accepted
    out_mu[it, ] <- mu_curr
    out_log_tau2[it, ] <- log_tau2_curr

    # ---------- Adaptation ----------
    if (it >= adapt_start) {
      # Robbins–Monro scale toward target acceptance
      log_scale <- log_scale + step_t(it - adapt_start + 1) * (as.numeric(accepted) - target_acc)
    }
    # Online empirical covariance (Haario-style, standard recursion without injecting cov_eps)
    mu_old <- mu_adapt
    mu_adapt <- mu_adapt + (x_curr - mu_adapt) / it
    C_adapt <- C_adapt + (tcrossprod(x_curr - mu_old) - C_adapt) / it

    # ---------- Selective refresh (every ess_check_every iters) ----------
    nref_this_check <- 0L
    if (ess_check_every > 0 && (it %% ess_check_every == 0)) {
      ord <- sample.int(N)  # randomize order to avoid bias
      for (ii in ord) {
        chk <- maybe_refresh_subject(
          i = ii, cache = caches[[ii]], mu = mu_curr, tau2 = tau2_curr,
          ess_thresh = ess_thresh, cooldown = cooldown,
          last_refresh_iter = last_refresh_iter, iter = it,
          refresh_fn = refresh_fn, verbose = verbose
        )
        if (isTRUE(chk$refreshed)) {
          caches[[ii]] <- chk$cache
          last_refresh_iter[[ii]] <- it
          # Ensure refreshed cache uses the master column order by names
          cn <- colnames(caches[[ii]]$theta)
          if (is.null(cn) || !setequal(cn, theta_dimnames)) {
            stop(sprintf("Refreshed cache %d has invalid or mismatched names.", ii))
          }
          if (!identical(cn, theta_dimnames)) {
            caches[[ii]]$theta <- caches[[ii]]$theta[, theta_dimnames, drop = FALSE]
          }
          # Rebuild cumw and re-map current-state indices for the refreshed subject
          cumw_list[[ii]] <- .make_cumw(caches[[ii]]$w)
          idx1_list[[ii]] <- sample_idx_from_uniforms(cumw_list[[ii]], pnorm(z1_list[[ii]]))
          idx2_list[[ii]] <- sample_idx_from_uniforms(cumw_list[[ii]], pnorm(z2_list[[ii]]))
          # When we refresh a subject at current η, we should also recompute current log-posteriors
          # at Stage-1 and Stage-2. Recompute whole sums (simple & safe).
          logLi1_curr <- 0; logLi2_curr <- 0
          for (jj in seq_len(N)) {
            logLi1_curr <- logLi1_curr + logZi_hat_from_cache(caches[[jj]], mu_curr, tau2_curr, idx1_list[[jj]], FALSE)$logZi
            logLi2_curr <- logLi2_curr + logZi_hat_from_cache(caches[[jj]], mu_curr, tau2_curr, idx2_list[[jj]], FALSE)$logZi
          }
          lpost1_curr <- lp_eta_curr + logLi1_curr
          lpost2_curr <- lp_eta_curr + logLi2_curr
          nref_this_check <- nref_this_check + 1L
          if (nref_this_check >= max_refresh_per_check) break
        }
      }
    }
    refresh_count[it] <- nref_this_check

    if (isTRUE(verbose) && (it %% max(50L, as.integer(n_iter/20)) == 0L)) {
      acc_rate <- mean(out_acc[seq_len(it)])
      cat(sprintf("iter %d/%d | acc=%.3f | stage1_pass=%.3f | log_scale=%.3f | refreshed=%d\n",
                  it, n_iter, acc_rate, mean(out_stage1_pass[seq_len(it)]), log_scale, nref_this_check))
    }
  } # end main loop

  list(
    mu = out_mu,
    log_tau2 = out_log_tau2,
    accepted = out_acc,
    stage1_pass = out_stage1_pass,
    refresh_count = refresh_count,
    final_state = list(mu = mu_curr, log_tau2 = log_tau2_curr),
    adapt = list(log_scale = log_scale, C = C_adapt),
    meta = list(K1 = K1, K2 = K2, rho = rho, target_acc = target_acc,
                indep_mix_prob = indep_mix_prob)
  )
}
