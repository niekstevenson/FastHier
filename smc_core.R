#!/usr/bin/env Rscript
# ============================================================================
# SMC Core Utilities
# - Shared numerics, rCESS scheduling, resampling, mixture helpers, and
#   likelihood batching for the local SMC workflow
# ============================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(matrixStats)
  library(mvtnorm)
  library(qrng)
})

# ------------------------------ small helpers ------------------------------
`%||%` <- function(a, b) if (!is.null(a)) a else b
.clamp <- function(x, a, b) pmin(pmax(x, a), b)

# ----------------------------- log-sum-exp utils ---------------------------
logsumexp <- function(x) { m <- max(x); m + log(sum(exp(x - m))) }
.rowLogSumExp <- function(M) matrixStats::rowLogSumExps(M)
rlogsumexp2 <- function(a, b) { m <- pmax(a, b); m + log(exp(a - m) + exp(b - m)) }

# --------------------------- basic stats helpers ---------------------------
ESS <- function(w) { s <- sum(w); (s * s) / sum(w * w) }
weighted_cov <- function(A, ww) {
  ww <- ww / sum(ww)
  mu <- drop(t(A) %*% ww)
  Z  <- sweep(A, 2L, mu, "-")
  crossprod(sqrt(ww) * Z)
}

# ----------------------------- mvn/mvt logpdf ------------------------------
dmvnorm_chol_log <- function(X, mu, L) {
  Xc <- sweep(X, 2L, mu, `-`)
  sol <- backsolve(L, t(Xc), transpose = TRUE)
  qf  <- colSums(sol^2)
  d   <- ncol(X)
  logdet <- 2 * sum(log(diag(L)))
  -0.5 * (d * log(2 * pi) + logdet + qf)
}

dmvt_chol_log <- function(X, mu, L, df) {
  Xc <- sweep(X, 2L, mu, `-`)
  sol <- backsolve(L, t(Xc), transpose = TRUE)
  qf  <- colSums(sol^2)
  d   <- ncol(X)
  logdet <- 2 * sum(log(diag(L)))
  lgamma((df + d) / 2) - lgamma(df / 2) - 0.5 * logdet - (d / 2) * log(df * pi) -
    ((df + d) / 2) * log1p(qf / df)
}

# ----------------------------- safe logdet ---------------------------------
safe_logdet <- function(S) {
  tryCatch({ L <- chol(S); 2 * sum(log(diag(L))) }, error = function(e) -Inf)
}

# ---------------------------- mixture utilities ---------------------------
prep_mix_cache <- function(meansZ, covsZ, wZ) {
  Ls <- lapply(covsZ, function(cov) {
    tryCatch(chol(cov), error = function(e) {
      d <- nrow(cov)
      ev <- eigen(cov, symmetric = TRUE)
      lam <- pmax(ev$values, 1e-6)
      cov_reg <- ev$vectors %*% diag(lam, d) %*% t(ev$vectors)
      tryCatch(chol(cov_reg), error = function(e2) chol(diag(d) + diag(1e-3, d)))
    })
  })
  logw <- log(wZ / sum(wZ))
  logdets <- vapply(Ls, function(L) 2 * sum(log(diag(L))), numeric(1))
  list(Ls = Ls, logw = logw, logdets = logdets)
}

gmm_logpdf_Z_vec <- function(Zm, meansZ, mix_cache) {
  n <- nrow(Zm); G <- length(mix_cache$logw)
  logs <- matrix(-Inf, n, G)
  for (g in 1:G) {
    logs[, g] <- mix_cache$logw[g] + dmvnorm_chol_log(Zm, meansZ[[g]], mix_cache$Ls[[g]])
  }
  .rowLogSumExp(logs)
}

dmvt_mixture_logpdf_Z_vec <- function(Zm, meansZ, mix_cache, nu) {
  n <- nrow(Zm); G <- length(mix_cache$logw)
  logs <- matrix(-Inf, n, G)
  for (g in 1:G) {
    logs[, g] <- mix_cache$logw[g] + dmvt_chol_log(Zm, meansZ[[g]], mix_cache$Ls[[g]], nu)
  }
  .rowLogSumExp(logs)
}

regularize_cov <- function(Sig, min_eig = 3e-3, cond_cap = 1e4) {
  Sig <- as.matrix((Sig + t(Sig)) / 2)
  ev <- eigen(Sig, symmetric = TRUE)
  vals <- pmax(ev$values, min_eig)
  vmax <- max(vals); vmin <- min(vals)
  if (vmax / vmin > cond_cap) vals <- pmax(vals, vmax / cond_cap)
  S <- ev$vectors %*% diag(vals, length(vals)) %*% t(ev$vectors)
  S + diag(1e-8, nrow(S))
}

merge_components <- function(m1, S1, w1, m2, S2, w2) {
  w <- w1 + w2; alpha <- w1 / w; m <- alpha * m1 + (1 - alpha) * m2
  d1 <- m1 - m; d2 <- m2 - m
  S <- alpha * (S1 + tcrossprod(d1)) + (1 - alpha) * (S2 + tcrossprod(d2))
  list(mean = m, cov = S, weight = w)
}

bhatt_distance <- function(m1, S1, m2, S2) {
  d  <- length(m1)
  S1 <- regularize_cov(S1); S2 <- regularize_cov(S2); S <- (S1 + S2) / 2
  Ls <- tryCatch(chol(Matrix::nearPD(S)$mat), error = function(e) chol(S + diag(1e-8, d)))
  dm <- as.numeric(m2 - m1)
  if (length(dm) != nrow(Ls)) return(NA_real_)
  z  <- backsolve(Ls, dm, transpose = TRUE)
  q  <- sum(z * z)
  0.125 * q + 0.5 * (safe_logdet(S) - 0.5 * (safe_logdet(S1) + safe_logdet(S2)))
}

prune_merge_mixture_Z <- function(meansZ, covsZ, wZ,
                                  w_floor = 0.005, merge_thresh = 0.10, max_G = 32,
                                  min_eig = 3e-3, cond_cap = 1e4, min_G_keep = 2,
                                  verbose = FALSE) {
  G <- length(wZ); if (G == 0) return(NULL)
  wZ <- pmax(wZ, 0)
  keep <- which(wZ / sum(wZ) >= w_floor)
  if (!length(keep)) keep <- which.max(wZ)
  meansZ <- meansZ[keep]; covsZ <- covsZ[keep]; wZ <- wZ[keep]; wZ <- wZ / sum(wZ)
  covsZ <- lapply(covsZ, regularize_cov, min_eig = min_eig, cond_cap = cond_cap)
  if (length(wZ) <= max(2, min_G_keep)) {
    mix_cache <- prep_mix_cache(meansZ, covsZ, wZ)
    return(list(meansZ = meansZ, covsZ = covsZ, wZ = wZ, cache = mix_cache))
  }
  changed <- TRUE; iter_guard <- 0L
  while (changed && iter_guard < 200L) {
    iter_guard <- iter_guard + 1L; changed <- FALSE; G <- length(wZ)
    if (G <= max(2, min_G_keep)) break
    need_merge <- (G > max_G)
    dmin <- Inf; imin <- 1; jmin <- 2
    for (i in 1:(G - 1)) for (j in (i + 1):G) {
      d <- tryCatch(bhatt_distance(meansZ[[i]], covsZ[[i]], meansZ[[j]], covsZ[[j]]), error = function(e) Inf)
      if (!is.finite(d)) d <- Inf
      if (d < dmin) { dmin <- d; imin <- i; jmin <- j }
    }
    if (!is.finite(dmin)) break
    if ((need_merge || dmin < merge_thresh) && (G - 1) >= min_G_keep) {
      m <- merge_components(meansZ[[imin]], covsZ[[imin]], wZ[imin], meansZ[[jmin]], covsZ[[jmin]], wZ[jmin])
      idx <- setdiff(seq_len(G), c(imin, jmin))
      meansZ <- c(meansZ[idx], list(m$mean))
      covsZ  <- c(covsZ[idx],  list(regularize_cov(m$cov, min_eig, cond_cap)))
      wZ     <- c(wZ[idx],      m$weight); wZ <- wZ / sum(wZ); changed <- TRUE
    }
  }
  mix_cache <- prep_mix_cache(meansZ, covsZ, wZ)
  list(meansZ = meansZ, covsZ = covsZ, wZ = wZ, cache = mix_cache)
}

inflate_cov_along_dims <- function(S, idx, factor = 3) {
  if (!length(idx)) return(S)
  S[idx, idx] <- S[idx, idx] * factor
  S
}

inflate_mixture_along_dims <- function(mix, weak_idx, factor = 3, min_eig = 3e-3, cond_cap = 5e3) {
  if (is.null(mix) || is.null(mix$covsZ) || !length(mix$covsZ)) return(mix)
  for (g in seq_along(mix$covsZ)) {
    S <- mix$covsZ[[g]]
    S <- inflate_cov_along_dims(S, weak_idx, factor)
    mix$covsZ[[g]] <- regularize_cov(S, min_eig = min_eig, cond_cap = cond_cap)
  }
  mix$cache <- prep_mix_cache(mix$meansZ, mix$covsZ, mix$wZ)
  mix
}

# ------------------------ weak-dimension identification ------------------------

# ------------------------ mixtures: defaults and blending ------------------------
# Default standard normal mixture in d dimensions
.default_std_normal_mix <- function(d) {
  mu0 <- rep(0, d)
  Sig0 <- diag(d)
  cache0 <- prep_mix_cache(list(mu0), list(Sig0), wZ = 1)
  list(meansZ = list(mu0), covsZ = list(Sig0), wZ = 1, cache = cache0)
}

# Check emptiness of mixture
.is_empty_mix <- function(mix) {
  is.null(mix) || is.null(mix$meansZ) || !length(mix$meansZ)
}

# Blend a base mixture into an existing mixture with weight eps
blend_mixes <- function(mix, base, eps = 0.05) {
  if (.is_empty_mix(base)) return(mix)
  if (.is_empty_mix(mix)) {
    base$cache <- prep_mix_cache(base$meansZ, base$covsZ, base$wZ)
    return(base)
  }
  meansZ <- c(mix$meansZ, base$meansZ)
  covsZ  <- c(mix$covsZ,  base$covsZ)
  wZ     <- c((1 - eps) * mix$wZ, eps * base$wZ)
  wZ     <- wZ / sum(wZ)
  list(meansZ = meansZ, covsZ = covsZ, wZ = wZ, cache = prep_mix_cache(meansZ, covsZ, wZ))
}

# --------------------- elite mixture fitting (whitened EM) ----------------------
weighted_gmm_em_whitened <- function(W, we, G, init = NULL, itmax = 10, tol = 1e-4) {
  n <- nrow(W); d <- ncol(W); we <- we / sum(we)
  if (is.null(init)) {
    means <- matrix(NA_real_, G, d)
    sel <- sample.int(n, 1, prob = we); means[1, ] <- W[sel, ]
    d2 <- rep(Inf, n)
    for (g in 2:G) {
      d2 <- pmin(d2, rowSums((W - matrix(means[g - 1, ], n, d, byrow = TRUE))^2))
      sel <- sample.int(n, 1, prob = we * d2 / sum(we * d2))
      means[g, ] <- W[sel, ]
    }
    Sigmas <- replicate(G, diag(d), simplify = FALSE)
    pis <- rep(1/G, G)
  } else { means <- init$means; Sigmas <- init$covs; pis <- init$weights }
  prev_ll <- -Inf
  for (it in 1:itmax) {
    logR <- matrix(NA_real_, n, G)
    Ls_it <- lapply(Sigmas, function(S) {
      tryCatch(chol(S), error = function(e) {
        dloc <- ncol(S)
        ev <- eigen(S, symmetric = TRUE)
        lam <- pmax(ev$values, 1e-3)
        Sreg <- ev$vectors %*% diag(lam, dloc) %*% t(ev$vectors)
        chol(Sreg)
      })
    })
    for (g in 1:G) logR[, g] <- log(pis[g]) + dmvnorm_chol_log(W, means[g, ], Ls_it[[g]])
    m <- apply(logR, 1L, max)
    R <- exp(logR - m)
    denom <- drop(R %*% rep(1, G))
    ll_curr <- sum(we * (m + log(pmax(denom, .Machine$double.eps))))
    R <- R / pmax(denom, .Machine$double.eps)
    wg <- colSums(R * we); pis <- pmax(wg, 1e-12); pis <- pis / sum(pis)
    for (g in 1:G) {
      rg <- R[, g] * we; sg <- sum(rg); if (sg < 1e-12) next
      mu <- colSums(W * rg) / sg
      Xc <- sweep(W, 2L, mu, `-`)
      Sig <- t(Xc) %*% (Xc * rg) / sg
      ev <- eigen(Sig, symmetric = TRUE)
      lam <- pmax(ev$values, 1e-3)
      Sig <- ev$vectors %*% diag(lam, d) %*% t(ev$vectors)
      means[g, ] <- mu; Sigmas[[g]] <- Sig
    }
    if (is.finite(prev_ll) && (ll_curr - prev_ll) < tol * (1 + abs(prev_ll))) break
    prev_ll <- ll_curr
  }
  list(means = means, covs = Sigmas, weights = pis)
}

fit_elite_mixture_Z <- function(Z, w, elite_quantile = 0.4, G = 12,
                                cov_inflation = 3.0, min_elite = 80,
                                warm_start_mixZ = NULL, em_itmax = 8,
                                housekeeping = TRUE, min_G_keep = 2, merge_thresh = 0.10,
                                verbose = FALSE,
                                lambda = NULL) {
  w <- pmax(w, 0); w <- w / sum(w)
  thr <- stats::quantile(w, 1 - elite_quantile)
  elite_idx <- if (sum(w >= thr) >= min_elite) which(w >= thr)
  else order(w, decreasing = TRUE)[seq_len(min(min_elite, length(w)))]
  Ze <- Z[elite_idx, , drop = FALSE]; we <- w[elite_idx]; we <- we / sum(we)
  d <- ncol(Z)
  muZ <- colSums(t(t(Ze) * we))
  CZ <- matrix(0, d, d)
  for (i in seq_len(nrow(Ze))) { v <- Ze[i, ] - muZ; CZ <- CZ + we[i] * tcrossprod(v, v) }
  L <- chol(CZ + diag(1e-8, d))
  Winv <- solve(L)
  W <- sweep(Ze, 2L, muZ, `-`) %*% Winv
  trim_q <- 0.995
  if (nrow(W) > max(G * 3, 30) && d >= 2) {
    r2 <- rowSums(W^2)
    thr_trim <- stats::qchisq(trim_q, df = d)
    keep <- (r2 <= thr_trim)
    if (sum(keep) >= max(min_elite/2, G * 3)) {
      Ze <- Ze[keep, , drop = FALSE]
      we <- we[keep]; we <- we / sum(we)
      muZ <- colSums(t(t(Ze) * we))
      CZ <- matrix(0, d, d)
      for (i in seq_len(nrow(Ze))) { v <- Ze[i, ] - muZ; CZ <- CZ + we[i] * tcrossprod(v, v) }
      L <- chol(CZ + diag(1e-8, d)); Winv <- solve(L); W <- sweep(Ze, 2L, muZ, `-`) %*% Winv
    }
  }
  init <- NULL
  if (!is.null(warm_start_mixZ)) {
    G0 <- length(warm_start_mixZ$meansZ)
    meansW <- matrix(0, G0, d)
    covsW <- vector("list", G0)
    for (g in 1:G0) {
      meansW[g, ] <- as.numeric((warm_start_mixZ$meansZ[[g]] - muZ) %*% Winv)
      covsW[[g]]  <- Winv %*% warm_start_mixZ$covsZ[[g]] %*% t(Winv)
    }
    init <- list(means = meansW, covs = covsW, weights = warm_start_mixZ$wZ)
    if (G0 != G) {
      init$means <- init$means[seq_len(min(G0, G)), , drop = FALSE]
      init$covs  <- init$covs[seq_len(min(G0, G))]
      init$weights <- init$weights[seq_len(min(G0, G))]
      init$weights <- init$weights / sum(init$weights)
      G <- min(G0, G)
    }
  }
  em <- tryCatch(weighted_gmm_em_whitened(W, we, G, init = init, itmax = em_itmax), error = function(e) NULL)
  if (is.null(em)) {
    meansZ <- list(as.numeric(muZ))
    covsZ <- list(t(L) %*% diag(d) %*% L * cov_inflation)
    wZ <- 1
  } else {
    meansZ <- lapply(1:G, function(g) as.numeric(muZ + (em$means[g, , drop = TRUE] %*% L)))
    covsZ  <- lapply(1:G, function(g) {
      SigW <- em$covs[[g]]
      ev <- eigen(SigW, symmetric = TRUE)
      lam <- pmax(ev$values, 1e-3)
      SigW <- ev$vectors %*% diag(lam, d) %*% t(ev$vectors)
      t(L) %*% SigW %*% L * cov_inflation
    })
    wZ <- as.numeric(em$weights)
  }
  for (g in seq_along(meansZ)) names(meansZ[[g]]) <- colnames(Z)
  min_eig_eff <- if (!is.null(lambda)) max(8e-3, 5e-2 * (1 - lambda)) else 8e-3
  if (housekeeping) {
    prune_merge_mixture_Z(meansZ, covsZ, wZ,
                          w_floor=0.004, merge_thresh=merge_thresh,
                          max_G=32, min_eig=min_eig_eff, cond_cap=5e3,
                          min_G_keep=min_G_keep, verbose=verbose)
  } else {
    mix_cache <- prep_mix_cache(meansZ, covsZ, wZ)
    list(meansZ = meansZ, covsZ = covsZ, wZ = wZ, cache = mix_cache)
  }
}

# ---------------------- history mixture builder ----------------------------
combine_elite_mixtures <- function(elite_mixtures, weights = c(0.6, 0.3, 0.1),
                                   housekeeping = TRUE, min_G_keep = 2, merge_thresh = 0.10,
                                   min_eig = 3e-3,
                                   verbose = FALSE) {
  n_mix <- length(elite_mixtures); if (!n_mix) return(NULL)
  valid_mixtures <- list()
  for (i in seq_len(n_mix)) {
    mix <- elite_mixtures[[i]]
    if (!is.null(mix) && length(mix$meansZ) > 0 && length(mix$covsZ) > 0 && length(mix$wZ) > 0) {
      d_check <- sapply(mix$meansZ, length)
      if (all(d_check == d_check[1]) && d_check[1] > 0) valid_mixtures[[length(valid_mixtures) + 1]] <- mix
    }
  }
  if (!length(valid_mixtures)) return(NULL)
  w_hist <- weights[seq_len(min(length(valid_mixtures), length(weights)))]; w_hist <- w_hist / sum(w_hist)
  all_meansZ <- list(); all_covsZ <- list(); all_wZ <- c()
  for (i in seq_len(length(valid_mixtures))) {
    mix <- valid_mixtures[[i]]; comp_weights <- w_hist[i] * mix$wZ
    all_meansZ <- c(all_meansZ, mix$meansZ); all_covsZ <- c(all_covsZ, mix$covsZ); all_wZ <- c(all_wZ, comp_weights)
  }
  if (!length(all_meansZ)) return(NULL)
  if (housekeeping) {
    prune_merge_mixture_Z(all_meansZ, all_covsZ, all_wZ,
                          w_floor=0.003, merge_thresh=merge_thresh,
                          max_G=48, min_eig=min_eig, cond_cap=1e4,
                          min_G_keep=min_G_keep, verbose=verbose)
  } else {
    all_wZ <- all_wZ / sum(all_wZ); mix_cache <- prep_mix_cache(all_meansZ, all_covsZ, all_wZ)
    list(meansZ = all_meansZ, covsZ = all_covsZ, wZ = all_wZ, cache = mix_cache)
  }
}

build_hist_mixture <- function(elite_history, lambda,
                               hist_mix_lambda_thresh, hist_mix_prob,
                               min_G_keep, merge_thresh,
                               min_eig = 3e-3,
                               force = FALSE) {
  if (!force && (lambda < hist_mix_lambda_thresh || length(elite_history) == 0)) return(NULL)
  ramp <- (lambda - hist_mix_lambda_thresh) / max(1 - hist_mix_lambda_thresh, 1e-8)
  ramp <- pmin(pmax(ramp, 0), 1)
  floor_prob <- min(0.10, 0.5 * hist_mix_prob)
  hist_prob_eff <- floor_prob + (hist_mix_prob - floor_prob) * ramp
  eh <- elite_history
  if (length(eh) > 6) eh <- tail(eh, 6L)
  w_hist <- c(0.40, 0.25, 0.15, 0.10, 0.06, 0.04)
  w_hist <- w_hist[seq_len(min(length(eh), length(w_hist)))]
  combined_mix <- combine_elite_mixtures(
    eh, weights = w_hist,
    housekeeping = TRUE, min_G_keep = min_G_keep,
    merge_thresh = merge_thresh, min_eig = min_eig, verbose = FALSE
  )
  if (!is.null(combined_mix)) list(mix = combined_mix, prob = hist_prob_eff) else NULL
}

# ------------------------------ resampling ---------------------------------

# Hilbert-sort order utilities (Skilling 2004 Axes->Transpose)
.axesto_transpose_uint <- function(x, bits) {
  x <- as.integer(x)
  d <- length(x)
  M <- bitwShiftL(1L, bits - 1L)
  # Inverse-undo
  Q <- M
  while (Q > 1L) {
    P <- Q - 1L
    for (i in seq_len(d)) {
      if (bitwAnd(x[i], Q) != 0L) {
        x[1L] <- bitwXor(x[1L], P)     # invert
      } else {
        t <- bitwAnd(bitwXor(x[1L], x[i]), P)  # exchange
        x[1L] <- bitwXor(x[1L], t)
        x[i]  <- bitwXor(x[i],  t)
      }
    }
    Q <- bitwShiftR(Q, 1L)
  }
  # Gray encode
  if (d >= 2L) for (i in 2L:d) x[i] <- bitwXor(x[i], x[i - 1L])
  t <- 0L
  Q <- M
  while (Q > 1L) {
    if (bitwAnd(x[d], Q) != 0L) t <- bitwXor(t, Q - 1L)
    Q <- bitwShiftR(Q, 1L)
  }
  for (i in seq_len(d)) x[i] <- bitwXor(x[i], t)
  x
}

hilbert_sort_order <- function(Z, bits = 16L) {
  Z <- as.matrix(Z)
  n <- nrow(Z); d <- ncol(Z)
  if (n <= 1L || d == 0L) return(seq_len(n))
  U <- pnorm(Z)
  U[U < 0] <- 0; U[U > 1] <- 1
  scale <- as.integer(bitwShiftL(1L, bits) - 1L)
  A <- pmin(pmax(floor(U * scale + 0.5), 0), scale)
  storage.mode(A) <- "integer"
  Tm <- A
  for (i in seq_len(n)) Tm[i, ] <- .axesto_transpose_uint(A[i, ], bits)
  key <- vector("list", bits * d)
  k <- 1L
  for (p in seq.int(bits - 1L, 0L)) {
    for (j in seq_len(d)) {
      key[[k]] <- as.integer(bitwAnd(bitwShiftR(Tm[, j], p), 1L))
      k <- k + 1L
    }
  }
  do.call(order, key)
}

stratified_resample_sorted <- function(w, deterministic = TRUE) {
  N <- length(w); if (N <= 0L) return(integer(0))
  w <- w / sum(w)
  cw <- c(0, cumsum(w))
  u0 <- if (deterministic) 0.5 / N else runif(1) / N
  u  <- u0 + (0:(N - 1L)) / N
  findInterval(u, cw, rightmost.closed = TRUE)
}

# ----------------------------- rCESS scheduling ---------------------------
rCESS <- function(w, loglik, delta) {
  w <- as.numeric(w)
  w[!is.finite(w) | w < 0] <- 0
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) return(0)
  w <- w / sw

  x <- as.numeric(loglik)
  ok <- is.finite(x)
  if (!any(ok)) return(0)
  x[!ok] <- min(x[ok])
  x <- x - max(x)

  a1 <- logsumexp(log(w) + delta * x)
  a2 <- logsumexp(log(w) + 2 * delta * x)
  val <- exp(2 * a1 - a2)
  if (!is.finite(val)) return(0)
  pmin(pmax(val, 0), 1)
}

rCESS_stat <- function(w, h, delta) {
  w <- as.numeric(w)
  w[!is.finite(w) | w < 0] <- 0
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) return(0)
  w <- w / sw

  x <- as.numeric(h)
  ok <- is.finite(x)
  if (!any(ok)) return(0)
  x[!ok] <- min(x[ok])
  x <- x - max(x)

  a1 <- logsumexp(log(w) + delta * x)
  a2 <- logsumexp(log(w) + 2 * delta * x)
  val <- exp(2 * a1 - a2)
  if (!is.finite(val)) return(0)
  pmin(pmax(val, 0), 1)
}

cess_target_at_lambda <- function(lambda) {
  low <- 0.75; high <- 0.975
  tgt <- low + (high - low) * lambda
  pmin(pmax(tgt, 0), 0.999)
}

next_lambda_via_rCESS <- function(w, loglik, lambda,
                                  target = 0.92,
                                  eps_stop = 1e-6,
                                  itmax = 30,
                                  lambda_target = 1.0) {
  lambda_target <- as.numeric(lambda_target)
  if (lambda >= lambda_target - eps_stop) return(lambda_target)
  rem <- lambda_target - lambda
  r_full <- rCESS(w, loglik, rem)
  if (is.finite(r_full) && r_full >= target) return(lambda_target)
  lo <- 0.0; hi <- rem
  for (it in 1:itmax) {
    mid <- 0.5 * (lo + hi)
    r_mid <- rCESS(w, loglik, mid)
    if (!is.finite(r_mid)) r_mid <- 0
    if (r_mid >= target) lo <- mid else hi <- mid
    if ((hi - lo) <= max(1e-8, 0.02 * lo)) break
  }
  lambda + max(lo, min(hi, rem))
}

next_lambda_via_rCESS_stat <- function(w, h, lambda,
                                       target = 0.92,
                                       eps_stop = 1e-6,
                                       itmax = 30,
                                       lambda_target = 1.0) {
  lambda_target <- as.numeric(lambda_target)
  if (lambda >= lambda_target - eps_stop) return(lambda_target)
  rem <- lambda_target - lambda
  r_full <- rCESS_stat(w, h, rem)
  if (is.finite(r_full) && r_full >= target) return(lambda_target)
  lo <- 0.0; hi <- rem
  for (it in 1:itmax) {
    mid <- 0.5 * (lo + hi)
    r_mid <- rCESS_stat(w, h, mid)
    if (!is.finite(r_mid)) r_mid <- 0
    if (r_mid >= target) lo <- mid else hi <- mid
    if ((hi - lo) <= max(1e-8, 0.02 * lo)) break
  }
  lambda + max(lo, min(hi, rem))
}

# --------------------------- likelihood batching --------------------------
.ll_eval_one <- function(theta_row, data, loglik_fn) {
  val <- tryCatch(loglik_fn(theta_row, data), error = function(e) NULL)
  if (is.null(val)) {
    val <- tryCatch(
      loglik_fn(matrix(theta_row, nrow = 1L), data),
      error = function(e) stop("loglik_fn failed for both vector and 1-row matrix input.")
    )
  }
  val <- as.numeric(val)
  if (!length(val)) stop("loglik_fn returned empty output for one-row input.")
  val[1L]
}

.ll_eval_counter_state <- local({
  env <- new.env(parent = emptyenv())
  env$enabled <- FALSE
  env$total <- 0L
  env
})

.ll_eval_counter_add <- function(n) {
  if (!isTRUE(.ll_eval_counter_state$enabled)) return(invisible(NULL))
  .ll_eval_counter_state$total <- as.integer(.ll_eval_counter_state$total + as.integer(n))
  invisible(NULL)
}

.ll_eval_block <- function(Theta, data, loglik_fn) {
  out <- tryCatch(loglik_fn(Theta, data), error = function(e) NULL)
  if (!is.null(out)) {
    out <- as.numeric(out)
    if (length(out) == nrow(Theta)) {
      .ll_eval_counter_add(nrow(Theta))
      return(out)
    }
  }
  vals <- vapply(seq_len(nrow(Theta)), function(i) .ll_eval_one(Theta[i, , drop = TRUE], data, loglik_fn), numeric(1))
  .ll_eval_counter_add(nrow(Theta))
  vals
}

ll_parallel <- function(Theta, data, loglik_fn, n_cores = 1) {
  if (!is.matrix(Theta)) stop("Theta must be a matrix.")
  if (!is.function(loglik_fn)) stop("loglik_fn must be an R function.")
  n <- nrow(Theta)
  if (n == 0L) return(numeric(0))
  if (n_cores <= 1L || n < 2L) return(.ll_eval_block(Theta, data, loglik_fn))
  n_cores <- min(n_cores, parallel::detectCores(logical = TRUE), n)
  idx_list <- parallel::splitIndices(n, n_cores)
  parts <- parallel::mclapply(
    idx_list,
    function(ii) .ll_eval_block(Theta[ii, , drop = FALSE], data, loglik_fn),
    mc.cores = n_cores
  )
  out <- numeric(n)
  out[unlist(idx_list, use.names = FALSE)] <- unlist(parts, use.names = FALSE)
  out
}

# ----------------------------- QMC sampling -------------------------------
sample_gmm_Z_qmc <- function(M, meansZ, mix_cache, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  G <- length(meansZ); d <- length(meansZ[[1]])
  U <- qrng::sobol(n = M, d = d + 1, randomize = TRUE)
  cw <- c(0, cumsum(exp(mix_cache$logw)))
  comp <- findInterval(U[, 1], cw, rightmost.closed = TRUE)
  Zstd <- qnorm(U[, 2:(d + 1), drop = FALSE])
  Z <- matrix(NA_real_, M, d)
  for (g in seq_len(G)) {
    idx <- which(comp == g); if (!length(idx)) next
    L <- mix_cache$Ls[[g]]
    Z[idx, ] <- Zstd[idx, , drop = FALSE] %*% t(L)
    Z[idx, ] <- sweep(Z[idx, , drop = FALSE], 2L, meansZ[[g]], `+`)
  }
  colnames(Z) <- names(meansZ[[1]])
  Z
}

rmvt_mixture_Z_qmc <- function(M, meansZ, mix_cache, nu, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  G <- length(meansZ); d <- length(meansZ[[1]])
  U <- qrng::sobol(n = M, d = d + 2, randomize = TRUE)
  cw <- c(0, cumsum(exp(mix_cache$logw)))
  comp <- findInterval(U[, 1], cw, rightmost.closed = TRUE)
  Zstd <- qnorm(U[, 2:(d + 1), drop = FALSE])
  s <- sqrt(qchisq(U[, d + 2], df = nu) / nu)
  Z <- matrix(NA_real_, M, d)
  for (g in seq_len(G)) {
    idx <- which(comp == g); if (!length(idx)) next
    L <- mix_cache$Ls[[g]]; Y <- Zstd[idx, , drop = FALSE] %*% t(L)
    Z[idx, ] <- sweep(Y / s[idx], 2L, meansZ[[g]], `+`)
  }
  colnames(Z) <- names(meansZ[[1]])
  Z
}

# ============================================================================
# End of SMC core utilities
# ============================================================================
