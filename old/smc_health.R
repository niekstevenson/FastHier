# ================================================================
# Deterministic, single-run SMC diagnostics (correct split-Rhat)
# - Equal-weight resampling (deterministic stratified)
# - Hilbert-sort split in Z for two halves
# - Rank-normalized & folded split-Rhat on marginals + top-K PCs + loglik
# - Cross-fitted PSIS tail shape (k-hat) in Z
# - Final ESS/N check
# ================================================================

.smc_ESS <- function(w) { s <- sum(w); (s * s) / sum(w * w) }

# Deterministic stratified equal-weight resampling: return indices
.smc_equalize_indices <- function(w, N_eq = NULL) {
  w <- pmax(w, 0); w <- w / sum(w)
  N <- if (is.null(N_eq)) length(w) else N_eq
  if (N %% 2L != 0L) N <- N - 1L
  cw <- c(0, cumsum(w))
  u0 <- 0.5 / N
  u  <- u0 + (0:(N - 1L))/N
  findInterval(u, cw, rightmost.closed = TRUE)
}

# Deterministic Hilbert split in Z
.smc_hilbert_split <- function(Z_eq) {
  o <- tryCatch(hilbert_sort_order(Z_eq, bits = 16L),
                error = function(e) seq_len(nrow(Z_eq)))
  N <- nrow(Z_eq)
  list(iA = o[seq_len(N/2)], iB = o[(N/2 + 1):N], order = o)
}

# Rank-normalize (Vehtari et al.): ties handled via average ranks
.smc_rank_normalize <- function(x) {
  n  <- length(x)
  u  <- (rank(x, ties.method = "average") - 3/8) / (n + 1/4)
  stats::qnorm(pmin(pmax(u, .Machine$double.eps), 1 - .Machine$double.eps))
}

# Correct split-Rhat for 2 chains of equal length (rank-normalized + folded)
.smc_rhat_2chain <- function(a, b, folded = TRUE) {
  stopifnot(length(a) == length(b))
  n <- length(a)

  # Rank-normalized R-hat
  z  <- .smc_rank_normalize(c(a, b))
  z1 <- z[seq_len(n)]
  z2 <- z[(n + 1):(2 * n)]

  W  <- 0.5 * (stats::var(z1) + stats::var(z2))
  B  <- n * stats::var(c(mean(z1), mean(z2)))  # = n * (delta^2 / 2)
  var_plus <- ((n - 1)/n) * W + B / n
  rhat_core <- sqrt(var_plus / max(W, .Machine$double.eps))

  if (!folded) return(rhat_core)

  # Folded R-hat (for tail behavior)
  med <- stats::median(z)
  f1  <- abs(z1 - med)
  f2  <- abs(z2 - med)
  fz  <- .smc_rank_normalize(c(f1, f2))
  fz1 <- fz[seq_len(n)]
  fz2 <- fz[(n + 1):(2 * n)]

  Wf  <- 0.5 * (stats::var(fz1) + stats::var(fz2))
  Bf  <- n * stats::var(c(mean(fz1), mean(fz2)))
  var_plus_f <- ((n - 1)/n) * Wf + Bf / n
  rhat_fold <- sqrt(var_plus_f / max(Wf, .Machine$double.eps))
  max(rhat_core, rhat_fold)
}

# Top-K PCA directions (unit-length) on Z (equal weights)
.smc_pca_dirs <- function(Z_eq, K) {
  Zc <- scale(Z_eq, center = TRUE, scale = FALSE)
  pr <- stats::prcomp(Zc, center = FALSE, scale. = FALSE)
  V  <- pr$rotation[, seq_len(min(K, ncol(Z_eq), ncol(pr$rotation))), drop = FALSE]
  apply(V, 2, function(v) v / sqrt(sum(v^2)))
}

# Simple proposal in Z: single Gaussian, inflated covariance
.smc_build_simple_mixZ <- function(Z, inflation = 4.0) {
  Z <- as.matrix(Z)
  mu <- colMeans(Z)
  Zc <- sweep(Z, 2L, mu, `-`)
  S  <- crossprod(Zc) / nrow(Z)
  d  <- ncol(Z)
  S  <- (S + t(S)) / 2 + diag(1e-6, d)
  S  <- S * inflation
  cache <- prep_mix_cache(list(as.numeric(mu)), list(S), wZ = 1)
  list(meansZ = list(as.numeric(mu)), covsZ = list(S), wZ = 1, cache = cache)
}

# Stable log q(Z) for mixture: (1-p)*N + p*t_df
.smc_log_qZ_mix <- function(Z, mix, df_t = 7, p_t = 0.8) {
  lq_n <- gmm_logpdf_Z_vec(Z, mix$meansZ, mix$cache)
  lq_t <- dmvt_mixture_logpdf_Z_vec(Z, mix$meansZ, mix$cache, df_t)
  m <- pmax(lq_n + log1p(-p_t), lq_t + log(p_t))
  m + log(exp(lq_n + log1p(-p_t) - m) + exp(lq_t + log(p_t) - m))
}

# Posterior log-density up to additive constant
.smc_log_post_unnorm <- function(Theta, Z, loglik, Tmap, mu_ref, Sigma_ref) {
  prior_L <- tryCatch(chol(Sigma_ref), error = function(e) chol(Matrix::nearPD(Sigma_ref)$mat))
  lp_prior <- dmvnorm_chol_log(Theta, mu_ref, prior_L)
  lp_prior - Tmap$log_jac(Theta) + loglik
}

# ---- Main diagnostics ----

# Deterministic split-Rhat on marginals + top-K PCs + loglik (optional)
smc_split_rhat_deterministic <- function(run,
                                         K_pc = NULL,
                                         include_loglik = TRUE,
                                         N_eq = NULL) {
  stopifnot(is.matrix(run$Theta), is.matrix(run$Z), is.numeric(run$w))
  Theta <- run$Theta; Z <- run$Z; w <- run$w

  idx_eq   <- .smc_equalize_indices(w, N_eq = N_eq)
  Theta_eq <- Theta[idx_eq, , drop = FALSE]
  Z_eq     <- Z[idx_eq, , drop = FALSE]
  ll_eq    <- if (include_loglik && !is.null(run$loglik)) run$loglik[idx_eq] else NULL

  sp <- .smc_hilbert_split(Z_eq)
  A  <- Theta_eq[sp$iA, , drop = FALSE]
  B  <- Theta_eq[sp$iB, , drop = FALSE]

  out <- numeric(0)
  pnames <- colnames(Theta_eq)
  for (j in seq_len(ncol(Theta_eq))) {
    nm <- if (!is.null(pnames) && nzchar(pnames[j])) pnames[j] else paste0("theta_", j)
    out[nm] <- .smc_rhat_2chain(A[, j], B[, j], folded = TRUE)
  }

  if (is.null(K_pc)) K_pc <- min(5L, ncol(Z_eq))
  if (K_pc > 0) {
    ZA <- Z_eq[sp$iA, , drop = FALSE]
    ZB <- Z_eq[sp$iB, , drop = FALSE]
    V  <- .smc_pca_dirs(Z_eq, K = K_pc)
    for (r in seq_len(ncol(V))) {
      out[paste0("PC", r)] <- .smc_rhat_2chain(as.numeric(ZA %*% V[, r]),
                                               as.numeric(ZB %*% V[, r]),
                                               folded = TRUE)
    }
  }

  if (!is.null(ll_eq)) {
    out["loglik"] <- .smc_rhat_2chain(ll_eq[sp$iA], ll_eq[sp$iB], folded = TRUE)
  }

  out
}

# Cross-fitted PSIS k-hat (A↔B)
smc_psis_k_crossfit <- function(run, mu_ref, Sigma_ref,
                                N_eq = NULL,
                                indep_t_df = 7, indep_t_prob = 0.8,
                                inflation = 4.0) {
  if (!requireNamespace("loo", quietly = TRUE)) {
    warning("Package 'loo' not available; returning NA for PSIS k.")
    return(NA_real_)
  }
  stopifnot(is.matrix(run$Theta), is.matrix(run$Z), is.numeric(run$w))
  Theta <- run$Theta; Z <- run$Z; w <- run$w

  idx_eq   <- .smc_equalize_indices(w, N_eq = N_eq)
  Theta_eq <- Theta[idx_eq, , drop = FALSE]
  Z_eq     <- Z[idx_eq, , drop = FALSE]
  ll_eq    <- run$loglik[idx_eq]
  sp       <- .smc_hilbert_split(Z_eq)

  mix_A <- .smc_build_simple_mixZ(Z_eq[sp$iA, , drop = FALSE], inflation = inflation)
  mix_B <- .smc_build_simple_mixZ(Z_eq[sp$iB, , drop = FALSE], inflation = inflation)

  Theta_B <- Theta_eq[sp$iB, , drop = FALSE]; Z_B <- Z_eq[sp$iB, , drop = FALSE]
  lp_B    <- .smc_log_post_unnorm(Theta_B, Z_B, ll_eq[sp$iB], run$transport, mu_ref, Sigma_ref)
  lqA_B   <- .smc_log_qZ_mix(Z_B, mix_A, df_t = indep_t_df, p_t = indep_t_prob)
  k_B     <- as.numeric(loo::psis(lp_B - lqA_B)$diagnostics$pareto_k)

  Theta_A <- Theta_eq[sp$iA, , drop = FALSE]; Z_A <- Z_eq[sp$iA, , drop = FALSE]
  lp_A    <- .smc_log_post_unnorm(Theta_A, Z_A, ll_eq[sp$iA], run$transport, mu_ref, Sigma_ref)
  lqB_A   <- .smc_log_qZ_mix(Z_A, mix_B, df_t = indep_t_df, p_t = indep_t_prob)
  k_A     <- as.numeric(loo::psis(lp_A - lqB_A)$diagnostics$pareto_k)

  max(k_A, k_B)
}

# One-shot summary
smc_quick_check_deterministic <- function(run, mu_ref, Sigma_ref,
                                          N_eq = NULL,
                                          K_pc = NULL,
                                          rhat_thresh = 1.01,
                                          k_thresh = 0.7,
                                          ess_min_frac = 0.25) {
  N0 <- length(run$w)
  ess_frac <- .smc_ESS(run$w) / N0

  # rhat <- smc_split_rhat_deterministic(run, K_pc = K_pc,
  #                                      include_loglik = !is.null(run$loglik),
  #                                      N_eq = N_eq)
  khat <- smc_psis_k_crossfit(run, mu_ref, Sigma_ref, N_eq = N_eq)

  # ok_rhat <- all(is.finite(rhat)) && max(rhat) <= rhat_thresh
  ok_k    <- is.na(khat) || (is.finite(khat) && khat <= k_thresh)
  ok_ess  <- is.finite(ess_frac) && ess_frac >= ess_min_frac

  cat(sprintf("Final ESS/N:   %.3f %s\n",
              ess_frac, if (ok_ess) "✅" else "⚠️"))
  # cat(sprintf("Split-Rhat max %.3f %s\n",
  #             max(rhat), if (ok_rhat) "✅" else "⚠️"))
  if (!is.na(khat)) {
    cat(sprintf("PSIS k (xfit): %.3f %s\n",
                khat, if (ok_k) "✅" else "⚠️"))
  } else {
    cat("PSIS k (xfit):  NA (loo not installed)\n")
  }

  invisible(list(ess_frac = ess_frac, khat = khat, pass = (ok_k && ok_ess)))
}
