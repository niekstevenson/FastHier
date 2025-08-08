# ========================================================================
# ENHANCED SMC SAMPLER WITH TRANSPORT MAPS AND ELITE/HISTORY MIXTURE PROPOSALS
# Patched version (performance fix):
#   • Keep Tier A(1–3), 3A/3B/3C, 4a/4b
#   • CRITICAL FIX: Re-vectorize independence kernel
#       - For mtm_K == 1: fully vectorized (single batched loglik call)
#       - For mtm_K > 1: batched over K (small K), no per-particle loops
#   • Defaults keep simple cases fast (mtm_K = 1); extra work only for hard cases
# ========================================================================

suppressPackageStartupMessages({
  library(mvtnorm)
  library(qrng)
  library(Matrix)
})

# ------------- utilities (unchanged) -------------
logsumexp <- function(x, w = NULL) {
  if (is.null(w)) { m <- max(x); m + log(sum(exp(x - m))) }
  else            { m <- max(x); m + log(sum(w * exp(x - m))) }
}
ESS <- function(w) { s <- sum(w); (s*s) / sum(w*w) }
weighted_mean <- function(A, ww) { ww <- ww / sum(ww); drop(t(A) %*% ww) }
weighted_cov  <- function(A, ww) {
  ww <- ww / sum(ww); mu <- drop(t(A) %*% ww)
  Z  <- sweep(A, 2L, mu, "-"); crossprod(sqrt(ww) * Z)
}
dmvnorm_chol_log <- function(X, mu, L) {
  Xc <- sweep(X, 2L, mu, `-`)
  sol <- backsolve(L, t(Xc), transpose = TRUE)
  qf  <- colSums(sol^2); d <- ncol(X); logdet <- sum(log(diag(L))) * 2
  -0.5 * (d * log(2*pi) + logdet + qf)
}
dmvt_chol_log <- function(X, mu, L, df) {
  Xc <- sweep(X, 2L, mu, `-`)
  sol <- backsolve(L, t(Xc), transpose = TRUE)
  qf  <- colSums(sol^2); d <- ncol(X); logdet <- sum(log(diag(L))) * 2
  lgamma((df + d)/2) - lgamma(df/2) - 0.5*logdet - (d/2)*log(df*pi) -
    ((df + d)/2) * log1p(qf/df)
}

# ------------- tail-safe marginal transform (unchanged) -------------
.make_tail_safe_cdf_bundle <- function(ux, cw, Fraw, slope_raw) {
  stopifnot(length(ux) >= 2L)
  xmin <- ux[1]; xmax <- ux[length(ux)]
  uL <- cw[1]; uR <- cw[length(cw)]
  uL <- pmin(pmax(uL, 1e-12), 1 - 1e-12)
  uR <- pmin(pmax(uR, 1e-12), 1 - 1e-12)
  vL <- qnorm(uL); vR <- qnorm(uR)
  sL <- max(slope_raw(xmin), 1e-12, na.rm = TRUE)
  sR <- max(slope_raw(xmax), 1e-12, na.rm = TRUE)
  gL <- sL / dnorm(vL); gR <- sR / dnorm(vR)
  Fhat <- function(x) {
    x <- as.numeric(x); out <- numeric(length(x))
    lt <- x < xmin; rt <- x > xmax; md <- !(lt | rt)
    if (any(lt)) { v <- vL + gL * (x[lt] - xmin); out[lt] <- pnorm(v) }
    if (any(rt)) { v <- vR + gR * (x[rt] - xmax); out[rt] <- pnorm(v) }
    if (any(md)) out[md] <- Fraw(x[md])
    pmin(pmax(out, 1e-12), 1 - 1e-12)
  }
  slope_safe <- function(x) {
    x <- as.numeric(x); out <- numeric(length(x))
    lt <- x < xmin; rt <- x > xmax; md <- !(lt | rt)
    if (any(lt)) { v <- vL + gL * (x[lt] - xmin); out[lt] <- dnorm(v) * gL }
    if (any(rt)) { v <- vR + gR * (x[rt] - xmax); out[rt] <- dnorm(v) * gR }
    if (any(md)) { val <- slope_raw(x[md]); val[!is.finite(val)] <- 0; out[md] <- pmax(val, 1e-12) }
    out
  }
  Qhat_safe <- function(u) {
    u <- pmin(pmax(u, 1e-12), 1 - 1e-12); out <- numeric(length(u))
    lt <- u < uL; rt <- u > uR; md <- !(lt | rt)
    if (any(md)) out[md] <- splinefun(cw, ux, method = "monoH.FC")(u[md])
    if (any(lt)) { v <- qnorm(u[lt]); out[lt] <- xmin + (v - vL)/gL }
    if (any(rt)) { v <- qnorm(u[rt]); out[rt] <- xmax + (v - vR)/gR }
    out
  }
  list(Fhat = Fhat, slope = slope_safe, Qhat = Qhat_safe,
       xmin = xmin, xmax = xmax, uL = uL, uR = uR, vL = vL, vR = vR, gL = gL, gR = gR)
}

fit_copula_transform <- function(X, w,
                                 ngrid = 400,
                                 tail = 1e-3,
                                 slope_floor = 1e-6,
                                 slope_cap = Inf,
                                 corr_shrink = 0.75,
                                 use_weights_in_corr = TRUE,
                                 disable_jacobian = FALSE) {
  X <- as.matrix(X); d <- ncol(X); w <- pmax(w, 0); w <- w/sum(w)
  pgrid <- seq(tail, 1 - tail, length.out = ngrid)
  wtd_quantile <- function(x, w, p) {
    o <- order(x); x <- x[o]; w <- w[o]
    cw <- cumsum(w); cw <- cw/cw[length(cw)]
    approx(cw, x, xout = p, ties = "ordered", rule = 2)$y
  }
  per_dim <- vector("list", d)
  for (j in seq_len(d)) {
    xj <- X[, j]
    xq <- wtd_quantile(xj, w, pgrid)
    xq <- cummax(xq)
    # ensure at least two distinct points; if flat -> degenerate dim
    dx <- diff(xq)
    if (all(abs(dx) < 1e-12)) {
      # dimension is essentially constant
      x_const <- xq[1]
      per_dim[[j]] <- list(
        Fhat = function(x) rep(0.5, length(x)),           # put all mass at centre
        Qhat = function(u) rep(x_const, length(u)),       # always returns the constant
        pdf  = function(x) rep(1e-12, length(x))          # tiny but finite density
      )
      next                                                     # go to next dimension
    }

    # jitter any ties so splinefun sees strictly increasing x
    eps  <- max(1e-12, 1e-8 * abs(xq[1]))
    for (k in seq_along(xq)[-1]) {
      if (xq[k] <= xq[k - 1]) xq[k] <- xq[k - 1] + eps
    }
    # ---------------------------------------------------------------

    # these two calls are unchanged *except* that the Fraw one is now
    # inside suppressWarnings(...) to silence the duplicate-x message
    Qraw <- splinefun(pgrid, xq, method = "monoH.FC")
    Fraw <- suppressWarnings(splinefun(xq, pgrid, method = "monoH.FC"))
    slope_raw <- function(x) Fraw(x, deriv = 1)
    bundle <- .make_tail_safe_cdf_bundle(xq, pgrid, Fraw, slope_raw)
    per_dim[[j]] <- list(
      Fhat = function(x) bundle$Fhat(x),
      Qhat = function(u) bundle$Qhat(u),
      pdf  = function(x) { out <- bundle$slope(x); out[!is.finite(out)] <- slope_floor
      pmin(pmax(out, slope_floor), slope_cap) }
    )
  }
  fwd_y <- function(Xnew) {
    Xnew <- as.matrix(Xnew); Y <- Xnew
    for (j in seq_len(ncol(Xnew))) {
      u <- per_dim[[j]]$Fhat(Xnew[, j]); Y[, j] <- qnorm(u)
    }
    colnames(Y) <- colnames(Xnew); Y
  }
  Y0 <- fwd_y(X)
  if (use_weights_in_corr) {
    mu <- colSums(Y0 * w); Yc <- sweep(Y0, 2L, mu, `-`)
    S  <- t(Yc) %*% (Yc * w) / sum(w)
    sdv <- sqrt(pmax(diag(S), 1e-12)); Rhat <- S / (sdv %o% sdv)
  } else Rhat <- stats::cor(Y0)
  Rhat[!is.finite(Rhat)] <- 0; diag(Rhat) <- 1
  R <- (1 - corr_shrink) * diag(d) + corr_shrink * Rhat
  R <- as.matrix(Matrix::nearPD(R, conv.tol = 1e-6)$mat)
  U <- chol(R)
  fwd <- function(Xnew) { Y <- fwd_y(Xnew); Zt <- backsolve(U, t(Y), transpose = FALSE); Z <- t(Zt); colnames(Z) <- colnames(Xnew); Z }
  inv <- function(Znew) {
    Znew <- as.matrix(Znew); Y <- Znew %*% U; Xr <- Y
    for (j in seq_len(ncol(Y))) { u <- pnorm(Y[, j]); Xr[, j] <- per_dim[[j]]$Qhat(u) }
    colnames(Xr) <- colnames(Znew); Xr
  }
  logdetU <- sum(log(diag(U)))
  log_jac <- function(Theta) {
    if (disable_jacobian) return(rep(0, nrow(as.matrix(Theta))))
    Theta <- as.matrix(Theta); Y <- fwd_y(Theta); s <- rep(0.0, nrow(Theta))
    for (j in seq_len(ncol(Theta))) s <- s + (log(per_dim[[j]]$pdf(Theta[, j])) - dnorm(Y[, j], log = TRUE))
    s - logdetU
  }
  list(fwd = fwd, inv = inv, log_jac = log_jac, meta = list(R = R, U = U))
}

# ------------- mixture helpers (unchanged) -------------
prep_mix_cache <- function(meansZ, covsZ, wZ) {
  Ls <- lapply(covsZ, function(cov) {
    tryCatch(chol(cov), error = function(e) {
      d <- nrow(cov); ev <- eigen(cov, symmetric = TRUE)
      lam <- pmax(ev$values, 1e-6)
      cov_reg <- ev$vectors %*% diag(lam, d) %*% t(ev$vectors)
      tryCatch(chol(cov_reg), error = function(e2) chol(diag(d) + diag(1e-3, d)))
    })
  })
  logw <- log(wZ / sum(wZ))
  logdets <- vapply(Ls, function(L) 2*sum(log(diag(L))), numeric(1))
  list(Ls = Ls, logw = logw, logdets = logdets)
}
gmm_logpdf_Z_vec <- function(Zm, meansZ, mix_cache) {
  n <- nrow(Zm); G <- length(mix_cache$logw)
  logs <- matrix(-Inf, n, G)
  for (g in seq_len(G)) logs[, g] <- mix_cache$logw[g] + dmvnorm_chol_log(Zm, meansZ[[g]], mix_cache$Ls[[g]])
  apply(logs, 1L, logsumexp)
}
dmvt_mixture_logpdf_Z_vec <- function(Zm, meansZ, mix_cache, nu) {
  n <- nrow(Zm); G <- length(mix_cache$logw)
  logs <- matrix(-Inf, n, G)
  for (g in seq_len(G)) logs[, g] <- mix_cache$logw[g] + dmvt_chol_log(Zm, meansZ[[g]], mix_cache$Ls[[g]], nu)
  apply(logs, 1L, logsumexp)
}

# --- housekeeping (unchanged except min_G_keep) ---
regularize_cov <- function(Sig, min_eig = 3e-3, cond_cap = 1e4) {
  Sig <- as.matrix((Sig + t(Sig)) / 2)
  ev <- eigen(Sig, symmetric = TRUE); vals <- pmax(ev$values, min_eig)
  vmax <- max(vals); vmin <- min(vals)
  if (vmax / vmin > cond_cap) vals <- pmax(vals, vmax / cond_cap)
  S <- ev$vectors %*% diag(vals, length(vals)) %*% t(ev$vectors)
  S + diag(1e-8, nrow(S))
}
safe_logdet <- function(A, floor_eig = 1e-10) {
  A <- as.matrix((A + t(A)) / 2)
  ev <- tryCatch(eigen(A, symmetric = TRUE, only.values = TRUE)$values, error = function(e) NULL)
  if (!is.null(ev)) return(sum(log(pmax(ev, floor_eig))))
  L <- tryCatch(chol(Matrix::nearPD(A)$mat), error = function(e) chol(A + diag(floor_eig, nrow(A))))
  2 * sum(log(diag(L)))
}
bhatt_distance <- function(m1, S1, m2, S2) {
  d <- length(m1); S1 <- regularize_cov(S1); S2 <- regularize_cov(S2); S <- (S1 + S2)/2
  Ls <- tryCatch(chol(Matrix::nearPD(S)$mat), error = function(e) chol(S + diag(1e-8, d)))
  dm  <- m2 - m1; sol <- backsolve(Ls, t(dm), transpose = TRUE); q <- sum(sol^2)
  0.125 * q + 0.5 * (safe_logdet(S) - 0.5*(safe_logdet(S1) + safe_logdet(S2)))
}
merge_components <- function(m1,S1,w1, m2,S2,w2) {
  w <- w1 + w2; alpha <- w1 / w; m <- alpha*m1 + (1 - alpha)*m2
  d1 <- m1 - m; d2 <- m2 - m
  S <- alpha*(S1 + tcrossprod(d1)) + (1 - alpha)*(S2 + tcrossprod(d2))
  list(mean = m, cov = S, weight = w)
}
prune_merge_mixture_Z <- function(meansZ, covsZ, wZ,
                                  w_floor = 0.005, merge_thresh = 0.10, max_G = 32,
                                  min_eig = 3e-3, cond_cap = 1e4, min_G_keep = 2,
                                  verbose = FALSE) {
  G <- length(wZ); if (G == 0) return(NULL)
  wZ <- pmax(wZ, 0); keep <- which(wZ / sum(wZ) >= w_floor); if (!length(keep)) keep <- which.max(wZ)
  meansZ <- meansZ[keep]; covsZ <- covsZ[keep]; wZ <- wZ[keep]; wZ <- wZ / sum(wZ)
  covsZ <- lapply(covsZ, regularize_cov, min_eig=min_eig, cond_cap=cond_cap)
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
    for (i in 1:(G-1)) for (j in (i+1):G) {
      d <- tryCatch(bhatt_distance(meansZ[[i]], covsZ[[i]], meansZ[[j]], covsZ[[j]]), error = function(e) Inf)
      if (!is.finite(d)) d <- Inf; if (d < dmin) { dmin <- d; imin <- i; jmin <- j }
    }
    if (!is.finite(dmin)) break
    if ((need_merge || dmin < merge_thresh) && (G - 1) >= min_G_keep) {
      if (verbose) cat(sprintf("    [mix] merging (%d,%d), dB=%.3f\n", imin, jmin, dmin))
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

# mixture sampling (unchanged)
sample_gmm_Z <- function(M, meansZ, mix_cache, seed = NULL, Zstd = NULL, comp_u = NULL) {
  if (!is.null(seed)) set.seed(seed)
  G <- length(meansZ); d <- length(meansZ[[1]])
  if (is.null(comp_u)) comp_u <- runif(M)
  cw <- c(0, cumsum(exp(mix_cache$logw))); cw[length(cw)] <- 1
  comp <- pmin(pmax(findInterval(comp_u, cw, rightmost.closed = TRUE), 1L), G)
  if (is.null(Zstd)) Zstd <- matrix(rnorm(M * d), M, d)
  Z <- matrix(NA_real_, M, d)
  for (g in seq_len(G)) {
    idx <- which(comp == g); if (!length(idx)) next
    L <- mix_cache$Ls[[g]]
    Z[idx, ] <- Zstd[idx, ] %*% t(L)
    Z[idx, ] <- sweep(Z[idx, ], 2L, meansZ[[g]], `+`)
  }
  colnames(Z) <- names(meansZ[[1]]); Z
}
sample_gmm_Z_qmc <- function(M, meansZ, mix_cache, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  G <- length(meansZ); d <- length(meansZ[[1]])
  U <- qrng::sobol(n = M, d = d + 1, randomize = TRUE)
  cw <- c(0, cumsum(exp(mix_cache$logw)))
  comp <- findInterval(U[,1], cw, rightmost.closed = TRUE)
  Zstd <- qnorm(U[, 2:(d + 1), drop = FALSE])
  Z <- matrix(NA_real_, M, d)
  for (g in seq_len(G)) {
    idx <- which(comp == g); if (!length(idx)) next
    L <- mix_cache$Ls[[g]]
    Z[idx, ] <- Zstd[idx, , drop = FALSE] %*% t(L)
    Z[idx, ] <- sweep(Z[idx, , drop = FALSE], 2L, meansZ[[g]], `+`)
  }
  colnames(Z) <- names(meansZ[[1]]); Z
}
rmvt_mixture_Z_qmc <- function(M, meansZ, mix_cache, nu, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  G <- length(meansZ); d <- length(meansZ[[1]])
  U <- qrng::sobol(n = M, d = d + 2, randomize = TRUE)
  cw <- c(0, cumsum(exp(mix_cache$logw)))
  comp <- findInterval(U[,1], cw, rightmost.closed = TRUE)
  Zstd <- qnorm(U[, 2:(d + 1), drop = FALSE])
  s <- sqrt(qchisq(U[, d + 2], df = nu) / nu)
  Z <- matrix(NA_real_, M, d)
  for (g in seq_len(G)) {
    idx <- which(comp == g); if (!length(idx)) next
    L <- mix_cache$Ls[[g]]; Y <- Zstd[idx, , drop = FALSE] %*% t(L)
    Z[idx, ] <- sweep(Y / s[idx], 2L, meansZ[[g]], `+`)
  }
  colnames(Z) <- names(meansZ[[1]]); Z
}

# ------------- difficulty / adaptation (unchanged) -------------
.clamp <- function(x, a, b) pmin(pmax(x, a), b)
.lerp  <- function(a, b, t) a + (b - a) * t
ema_update <- function(prev, new, beta = 0.3) { if (is.null(prev) || !is.finite(prev)) return(new); (1 - beta) * prev + beta * new }
eff_frac <- function(w) { w <- pmax(w, 0); w <- w / sum(w); n <- length(w); neff <- (sum(w)^2) / sum(w^2); neff / n }

difficulty_metrics <- function(loglik, w, acc_rate, ess_frac_post,
                               lambda, delta_lambda_planned, round,
                               config = list(ess_good=0.60, ess_bad=0.20, acc_good=0.30, acc_bad=0.10, ll_ref=80, slow_ref=0.03)) {
  n <- length(w)
  p_ess <- .clamp((config$ess_good - ess_frac_post) / (config$ess_good - config$ess_bad), 0, 1)
  p_acc <- .clamp((config$acc_good - acc_rate)     / (config$acc_good - config$acc_bad), 0, 1)
  ll_range <- max(loglik) - min(loglik); p_ll <- .clamp(ll_range / config$ll_ref, 0, 1)
  mw <- max(w / sum(w)); mw_min <- 1/n; mw_max <- 0.5; p_w <- .clamp((mw - mw_min) / (mw_max - mw_min), 0, 1)
  p_slow <- 0; if (lambda > 0.20 && round > 4) p_slow <- .clamp((config$slow_ref - delta_lambda_planned) / config$slow_ref, 0, 1)
  wts <- c(ess=.35, acc=.30, ll=.15, w=.10, slow=.10)
  idx <- sum(wts * c(p_ess, p_acc, p_ll, p_w, p_slow)) / sum(wts)
  list(index=idx, components=list(ess=p_ess, acc=p_acc, ll=p_ll, weight=p_w, slow=p_slow),
       raw=list(ess_frac_post=ess_frac_post, acc=acc_rate, ll_range=ll_range, max_w=mw))
}
knob_targets_from_difficulty <- function(base, d) {
  list(
    dalpha_cap         = .lerp(base$dalpha_cap, 0.04, d),
    n_mcmc_moves       = .lerp(base$n_mcmc_moves, base$n_mcmc_moves + 3, d),
    refit_every        = round(.lerp(max(1L, base$refit_every), 1L, d)),
    G_mix              = round(.clamp(.lerp(base$G_mix, min(32L, base$G_mix + 10L), d), 2L, 32L)),
    cov_inflation      = .clamp(.lerp(base$cov_inflation, min(10, base$cov_inflation * 2), d), 1.5, 10),
    rw_scale           = .clamp(.lerp(base$rw_scale, 0.06, d), 0.05, base$rw_scale),
    resample_threshold = .clamp(.lerp(base$resample_threshold, 0.30, d), 0.25, 0.80),
    indep_t_prob       = .clamp(.lerp(0.40, 0.98, d), 0.30, 0.98),
    indep_t_df         = round(.clamp(.lerp(9, 3, d), 3, 15))
  )
}
smooth_update_knobs <- function(current, targets, eta = 0.25) {
  out <- current
  for (nm in names(targets)) {
    v  <- current[[nm]]; vt <- targets[[nm]]
    if (is.null(v) || !is.finite(v)) { out[[nm]] <- vt; next }
    val <- v + eta * (vt - v)
    if (nm %in% c("n_mcmc_moves", "refit_every", "G_mix", "indep_t_df")) val <- round(val)
    out[[nm]] <- val
  }
  out
}
suggest_particles_smooth <- function(d, M, M_cap = 8000, d_start = 0.7) {
  if (d <= d_start || M >= M_cap) return(M)
  s <- (d - d_start) / (1 - d_start)
  min(M_cap, round(M * (1 + 0.5 * s)))
}

# ------------- EM / mixtures (unchanged) -------------
weighted_gmm_em_whitened <- function(W, we, G, init = NULL, itmax = 10) {
  n <- nrow(W); d <- ncol(W); we <- we / sum(we)
  if (is.null(init)) {
    means <- matrix(NA_real_, G, d); sel <- sample.int(n, 1, prob = we); means[1, ] <- W[sel, ]
    d2 <- rep(Inf, n)
    for (g in 2:G) { d2 <- pmin(d2, rowSums((W - matrix(means[g - 1, ], n, d, byrow = TRUE))^2))
    sel <- sample.int(n, 1, prob = we * d2 / sum(we * d2)); means[g, ] <- W[sel, ] }
    Sigmas <- replicate(G, diag(d), simplify = FALSE); pis <- rep(1/G, G)
  } else { means <- init$means; Sigmas <- init$covs; pis <- init$weights }
  for (it in 1:itmax) {
    logR <- matrix(NA_real_, n, G)
    for (g in 1:G) logR[, g] <- log(pis[g]) + mvtnorm::dmvnorm(W, mean = means[g, ], sigma = Sigmas[[g]], log = TRUE)
    m <- apply(logR, 1L, max); R <- exp(logR - m); denom <- drop(R %*% rep(1, G)); R <- R / pmax(denom, .Machine$double.eps)
    wg <- colSums(R * we); pis <- pmax(wg, 1e-12); pis <- pis / sum(pis)
    for (g in 1:G) {
      rg <- R[, g] * we; sg <- sum(rg); if (sg < 1e-12) next
      mu <- colSums(W * rg) / sg; Xc <- sweep(W, 2L, mu, `-`); Sig <- t(Xc) %*% (Xc * rg) / sg
      ev <- eigen(Sig, symmetric = TRUE); lam <- pmax(ev$values, 1e-3)
      Sig <- ev$vectors %*% diag(lam, d) %*% t(ev$vectors)
      means[g, ] <- mu; Sigmas[[g]] <- Sig
    }
  }
  list(means = means, covs = Sigmas, weights = pis)
}
combine_elite_mixtures <- function(elite_mixtures, weights = c(0.6, 0.3, 0.1),
                                   housekeeping = TRUE, min_G_keep = 2, verbose = FALSE) {
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
    prune_merge_mixture_Z(all_meansZ, all_covsZ, all_wZ, w_floor=0.003, merge_thresh=0.08,
                          max_G=48, min_eig=3e-3, cond_cap=1e4, min_G_keep=min_G_keep, verbose=verbose)
  } else {
    all_wZ <- all_wZ / sum(all_wZ); mix_cache <- prep_mix_cache(all_meansZ, all_covsZ, all_wZ)
    list(meansZ = all_meansZ, covsZ = all_covsZ, wZ = all_wZ, cache = mix_cache)
  }
}
fit_elite_mixture_Z <- function(Z, w, elite_quantile = 0.4, G = 12,
                                cov_inflation = 3.0, min_elite = 80,
                                warm_start_mixZ = NULL, em_itmax = 8,
                                housekeeping = TRUE, min_G_keep = 2, verbose = FALSE) {
  w <- pmax(w, 0); w <- w / sum(w)
  thr <- stats::quantile(w, 1 - elite_quantile)
  elite_idx <- if (sum(w >= thr) >= min_elite) which(w >= thr)
  else order(w, decreasing = TRUE)[seq_len(min(min_elite, length(w)))]
  Ze <- Z[elite_idx, , drop = FALSE]; we <- w[elite_idx]; we <- we / sum(we)
  d <- ncol(Z); muZ <- colSums(t(t(Ze) * we)); CZ <- matrix(0, d, d)
  for (i in seq_len(nrow(Ze))) { v <- Ze[i, ] - muZ; CZ <- CZ + we[i] * tcrossprod(v, v) }
  L <- chol(CZ + diag(1e-8, d)); Winv <- solve(L); W <- sweep(Ze, 2L, muZ, `-`) %*% Winv
  init <- NULL
  if (!is.null(warm_start_mixZ)) {
    G0 <- length(warm_start_mixZ$meansZ); meansW <- matrix(0, G0, d); covsW <- vector("list", G0)
    for (g in 1:G0) { meansW[g, ] <- as.numeric((warm_start_mixZ$meansZ[[g]] - muZ) %*% Winv)
    covsW[[g]]  <- Winv %*% warm_start_mixZ$covsZ[[g]] %*% t(Winv) }
    init <- list(means = meansW, covs = covsW, weights = warm_start_mixZ$wZ)
    if (G0 != G) { init$means <- init$means[seq_len(min(G0, G)), , drop = FALSE]
    init$covs  <- init$covs[seq_len(min(G0, G))]
    init$weights <- init$weights[seq_len(min(G0, G))]; init$weights <- init$weights / sum(init$weights); G <- min(G0, G) }
  }
  em <- tryCatch(weighted_gmm_em_whitened(W, we, G, init = init, itmax = em_itmax), error = function(e) NULL)
  if (is.null(em)) { meansZ <- list(as.numeric(muZ)); covsZ <- list(t(L) %*% diag(d) %*% L * cov_inflation); wZ <- 1 }
  else {
    meansZ <- lapply(1:G, function(g) as.numeric(muZ + (em$means[g, , drop = TRUE] %*% L)))
    covsZ  <- lapply(1:G, function(g) { SigW <- em$covs[[g]]
    ev <- eigen(SigW, symmetric = TRUE); lam <- pmax(ev$values, 1e-3)
    SigW <- ev$vectors %*% diag(lam, d) %*% t(ev$vectors); t(L) %*% SigW %*% L * cov_inflation })
    wZ <- as.numeric(em$weights)
  }
  for (g in seq_along(meansZ)) names(meansZ[[g]]) <- colnames(Z)
  if (housekeeping) prune_merge_mixture_Z(meansZ, covsZ, wZ, w_floor=0.004, merge_thresh=0.10,
                                          max_G=32, min_eig=3e-3, cond_cap=5e3, min_G_keep=min_G_keep, verbose=verbose)
  else { mix_cache <- prep_mix_cache(meansZ, covsZ, wZ); list(meansZ = meansZ, covsZ = covsZ, wZ = wZ, cache = mix_cache) }
}

# ------------- resampling -------------
systematic_resample <- function(w) {
  N <- length(w); w <- w / sum(w); u0 <- runif(1) / N; cw <- cumsum(w)
  idx <- integer(N); j <- 1L
  for (i in 1:N) { u <- u0 + (i - 1)/N; while (j < N && cw[j] < u) j <- j + 1L; idx[i] <- j }
  idx
}

# ------------- MCMC moves in Z: RW + Independence (vectorized) -------------
mcmc_moves_z_mix_batched <- function(Z, loglik, lpz, Tmap, lambda,
                                     mu_ref, prior_L, w,
                                     elite_mix, hist_mix = NULL,
                                     data, loglik_fn, n_moves = 1,
                                     rw_prob = 0.35, rw_scale = 0.9,
                                     indep_t_df = 7, indep_t_prob = 0.75,
                                     mtm_K = 1L,
                                     seed = NULL,
                                     param_names = NULL) {
  N <- nrow(Z); d <- ncol(Z)
  if (!is.null(seed)) set.seed(seed)
  if (is.null(param_names)) param_names <- colnames(Z)

  S <- weighted_cov(Z, w); if (any(!is.finite(S))) S <- diag(d)
  Lrw <- tryCatch(chol(S + diag(1e-8, d)), error = function(e) diag(d))
  step_scale <- rw_scale / sqrt(max(d, 1))

  acc_rw <- 0L; acc_id <- 0L
  prop_rw <- 0L; prop_id <- 0L
  esjd_rw_sum <- 0.0; esjd_id_sum <- 0.0

  rlogsumexp2 <- function(a, b) { m <- pmax(a, b); m + log(exp(a - m) + exp(b - m)) }
  log_q_mixture <- function(Zmat, elite, hist, indep_t_prob, hprob) {
    lq_el_n <- gmm_logpdf_Z_vec(Zmat, elite$meansZ, elite$cache)
    lq_el_t <- dmvt_mixture_logpdf_Z_vec(Zmat, elite$meansZ, elite$cache, indep_t_df)
    lq_el   <- rlogsumexp2(lq_el_n + log1p(-indep_t_prob), lq_el_t + log(indep_t_prob))
    if (is.null(hist)) return(lq_el)
    hm <- hist$mix
    lq_hi_n <- gmm_logpdf_Z_vec(Zmat, hm$meansZ, hm$cache)
    lq_hi_t <- dmvt_mixture_logpdf_Z_vec(Zmat, hm$meansZ, hm$cache, indep_t_df)
    lq_hi   <- rlogsumexp2(lq_hi_n + log1p(-indep_t_prob), lq_hi_t + log(indep_t_prob))
    rlogsumexp2(lq_el + log1p(-hprob), lq_hi + log(hprob))
  }

  # draw M samples from q (batched)
  sample_from_q <- function(M) {
    use_hist <- !is.null(hist_mix)
    hprob <- if (use_hist) hist_mix$prob else 0.0
    take_hist <- runif(M) < hprob
    t_flags   <- runif(M) < indep_t_prob
    Zout <- matrix(NA_real_, M, d)

    n_hi_t <- sum(take_hist & t_flags)
    n_hi_n <- sum(take_hist & !t_flags)
    n_el_t <- sum(!take_hist & t_flags)
    n_el_n <- sum(!take_hist & !t_flags)

    if (n_hi_n) { hm <- hist_mix$mix; Zout[take_hist & !t_flags, ] <- sample_gmm_Z_qmc(n_hi_n, hm$meansZ, hm$cache) }
    if (n_hi_t) { hm <- hist_mix$mix; Zout[take_hist &  t_flags, ] <- rmvt_mixture_Z_qmc(n_hi_t, hm$meansZ, hm$cache, nu = indep_t_df) }
    if (n_el_n) { Zout[!take_hist & !t_flags, ] <- sample_gmm_Z_qmc(n_el_n, elite_mix$meansZ, elite_mix$cache) }
    if (n_el_t) { Zout[!take_hist &  t_flags, ] <- rmvt_mixture_Z_qmc(n_el_t, elite_mix$meansZ, elite_mix$cache, nu = indep_t_df) }
    colnames(Zout) <- colnames(Z)
    list(Z = Zout, hprob = hprob)
  }

  # Random-walk (vectorized)
  do_rw <- function(idx) {
    if (!length(idx)) return(invisible(NULL))
    k <- length(idx); prop_rw <<- prop_rw + k
    Zc <- Z[idx, , drop = FALSE]
    eps <- matrix(rnorm(k * d), k, d) %*% t(Lrw)
    Zp <- Zc + eps * step_scale; colnames(Zp) <- param_names
    Theta_p <- Tmap$inv(Zp); colnames(Theta_p) <- param_names
    ll_p <- loglik_fn(Theta_p, data)
    lp_prior_p <- dmvnorm_chol_log(Theta_p, mu_ref, prior_L)
    lpz_p <- as.numeric(lp_prior_p - Tmap$log_jac(Theta_p))
    lt_p <- lpz_p + lambda * ll_p; lt_c <- lpz[idx] + lambda * loglik[idx]
    a <- lt_p - lt_c; u <- log(runif(k)); acc_idx <- which(u < a)
    if (length(acc_idx)) {
      dZ <- Zp[acc_idx, , drop=FALSE] - Zc[acc_idx, , drop=FALSE]
      esjd_rw_sum <<- esjd_rw_sum + sum(rowSums(dZ^2))
      Z[idx[acc_idx], ] <<- Zp[acc_idx, , drop = FALSE]
      loglik[idx[acc_idx]] <<- ll_p[acc_idx]
      lpz[idx[acc_idx]] <<- lpz_p[acc_idx]
      acc_rw <<- acc_rw + length(acc_idx)
    }
    invisible(NULL)
  }

  # Independence move
  do_indep <- function(idx) {
    if (!length(idx)) return(invisible(NULL))
    k <- length(idx); prop_id <<- prop_id + k
    Zc <- Z[idx, , drop = FALSE]
    lq_c <- log_q_mixture(Zc, elite_mix, hist_mix, indep_t_prob,
                          if (!is.null(hist_mix)) hist_mix$prob else 0.0)
    lt_c <- lpz[idx] + lambda * loglik[idx]

    if (mtm_K <= 1L) {
      # ---- FAST PATH (fully vectorized) ----
      smp <- sample_from_q(k); Zp <- smp$Z
      lq_p <- log_q_mixture(Zp, elite_mix, hist_mix, indep_t_prob, smp$hprob)
      Theta_p <- Tmap$inv(Zp); colnames(Theta_p) <- param_names
      ll_p <- loglik_fn(Theta_p, data)
      lpz_p <- as.numeric(dmvnorm_chol_log(Theta_p, mu_ref, prior_L) - Tmap$log_jac(Theta_p))
      lt_p <- lpz_p + lambda * ll_p

      a <- (lt_p - lq_p) - (lt_c - lq_c)
      u <- log(runif(k)); acc_idx <- which(u < pmin(0, -Inf + a) + a*0 + pmin(0, a))  # numerically stable min(1, exp(a))
      # simpler: accept if log u < a
      acc_idx <- which(u < a)
      if (length(acc_idx)) {
        dZ <- Zp[acc_idx, , drop=FALSE] - Zc[acc_idx, , drop=FALSE]
        esjd_id_sum <<- esjd_id_sum + sum(rowSums(dZ^2))
        Z[idx[acc_idx], ] <<- Zp[acc_idx, , drop = FALSE]
        loglik[idx[acc_idx]] <<- ll_p[acc_idx]
        lpz[idx[acc_idx]] <<- lpz_p[acc_idx]
        acc_id <<- acc_id + length(acc_idx)
      }
      return(invisible(NULL))
    }

    # ---- MTM-K (batched over K, still vectorized across particles) ----
    K <- mtm_K
    # forward candidates: make a big block (k*K x d) and then split
    Zcand_list <- vector("list", K)
    lq_fwd_list <- vector("list", K)
    Theta_fwd_all <- NULL
    for (j in 1:K) {
      smp <- sample_from_q(k)
      Zcand_list[[j]] <- smp$Z
      lq_fwd_list[[j]] <- log_q_mixture(smp$Z, elite_mix, hist_mix, indep_t_prob, smp$hprob)
      Theta_fwd_all <- rbind(Theta_fwd_all, Tmap$inv(smp$Z))
    }
    colnames(Theta_fwd_all) <- param_names
    ll_fwd_all <- loglik_fn(Theta_fwd_all, data)
    # split back per-try
    ll_fwd_list <- vector("list", K)
    lpz_fwd_list <- vector("list", K)
    offset <- 0L
    for (j in 1:K) {
      i1 <- offset + 1L; i2 <- offset + k; offset <- i2
      Theta_j <- Theta_fwd_all[i1:i2, , drop = FALSE]
      ll_fwd_list[[j]] <- ll_fwd_all[i1:i2]
      lpz_fwd_list[[j]] <- as.numeric(dmvnorm_chol_log(Theta_j, mu_ref, prior_L) - Tmap$log_jac(Theta_j))
    }
    # weights and selection
    wj <- matrix(0.0, k, K)
    for (j in 1:K) {
      lt_fwd_j <- lpz_fwd_list[[j]] + lambda * ll_fwd_list[[j]]
      wj[, j] <- exp(lt_fwd_j - lq_fwd_list[[j]])
    }
    sw <- rowSums(wj); sw[!is.finite(sw)] <- 0
    # select per row
    u <- runif(k) * sw
    cs <- t(apply(wj, 1L, cumsum))
    choose_j <- 1 + rowSums(sweep(cs, 1L, u, `<`))
    choose_j[!is.finite(choose_j)] <- 1L
    # chosen proposals
    Zp <- matrix(NA_real_, k, d)
    lt_p <- numeric(k); lq_p <- numeric(k)
    for (j in 1:K) {
      idxj <- which(choose_j == j)
      if (!length(idxj)) next
      Zp[idxj, ] <- Zcand_list[[j]][idxj, , drop = FALSE]
      lq_p[idxj] <- lq_fwd_list[[j]][idxj]
      lt_p[idxj] <- lpz_fwd_list[[j]][idxj] + lambda * ll_fwd_list[[j]][idxj]
    }

    # backward weights: current point + (K-1) fresh draws
    Kb <- K - 1L
    if (Kb > 0) {
      Zb_all <- vector("list", Kb)
      lq_b_all <- vector("list", Kb)
      Theta_b_big <- NULL
      for (j in 1:Kb) {
        sb <- sample_from_q(k)
        Zb_all[[j]] <- sb$Z
        lq_b_all[[j]] <- log_q_mixture(sb$Z, elite_mix, hist_mix, indep_t_prob, sb$hprob)
        Theta_b_big <- rbind(Theta_b_big, Tmap$inv(sb$Z))
      }
      colnames(Theta_b_big) <- param_names
      ll_b_big <- loglik_fn(Theta_b_big, data)

      # evaluate backward weights
      wb_sum <- numeric(k)
      # current point contribution
      wb_curr <- exp(lt_c - lq_c)
      wb_sum <- wb_sum + wb_curr
      # K-1 proposals
      off <- 0L
      for (j in 1:Kb) {
        i1 <- off + 1L; i2 <- off + k; off <- i2
        Theta_bj <- Theta_b_big[i1:i2, , drop = FALSE]
        lpz_bj <- as.numeric(dmvnorm_chol_log(Theta_bj, mu_ref, prior_L) - Tmap$log_jac(Theta_bj))
        lt_bj <- lpz_bj + lambda * ll_b_big[i1:i2]
        wb_sum <- wb_sum + exp(lt_bj - lq_b_all[[j]])
      }
      alpha <- pmin(1, sw / pmax(wb_sum, .Machine$double.eps))
    } else {
      # K=1 path would have been handled above; keep for completeness
      alpha <- pmin(1, exp((lt_p - lq_p) - (lt_c - lq_c)))
    }

    acc <- runif(k) < alpha
    if (any(acc)) {
      ia <- which(acc)
      Theta_p <- Tmap$inv(Zp[ia, , drop=FALSE]); colnames(Theta_p) <- param_names
      ll_p <- loglik_fn(Theta_p, data)
      lpz_p <- as.numeric(dmvnorm_chol_log(Theta_p, mu_ref, prior_L) - Tmap$log_jac(Theta_p))
      dZ <- Zp[ia, , drop=FALSE] - Zc[ia, , drop=FALSE]
      esjd_id_sum <<- esjd_id_sum + sum(rowSums(dZ^2))
      Z[idx[ia], ] <<- Zp[ia, , drop = FALSE]
      loglik[idx[ia]] <<- ll_p
      lpz[idx[ia]] <<- lpz_p
      acc_id <<- acc_id + length(ia)
    }
    invisible(NULL)
  }

  for (m in 1:n_moves) {
    rw_flags <- runif(N) < rw_prob
    idx_rw <- which(rw_flags)
    idx_id <- which(!rw_flags)
    do_rw(idx_rw)
    do_indep(idx_id)
  }

  acc_total <- (acc_rw + acc_id) / (N * n_moves)
  list(
    Z = Z, loglik = loglik, lpz = lpz,
    accept_rate = acc_total,
    rw_accept_rate   = if (prop_rw>0) acc_rw/prop_rw else 0,
    indep_accept_rate= if (prop_id>0) acc_id/prop_id else 0,
    esjd_rw = esjd_rw_sum,
    esjd_indep = esjd_id_sum,
    prop_rw = prop_rw,
    prop_indep = prop_id
  )
}

# ------------- convergence helpers, CESS (with early trim), adaptive M (unchanged) -------------
split_rhat_weighted <- function(Theta, w) {
  n <- nrow(Theta); if (n < 4) return(Inf)
  o <- order(runif(n)); i1 <- o[seq(1, n, by = 2)]; i2 <- o[seq(2, n, by = 2)]
  W1 <- w[i1]/sum(w[i1]); W2 <- w[i2]/sum(w[i2])
  mu1 <- drop(t(Theta[i1, , drop=FALSE]) %*% W1)
  mu2 <- drop(t(Theta[i2, , drop=FALSE]) %*% W2)
  s21 <- diag(weighted_cov(Theta[i1,,drop=FALSE], W1))
  s22 <- diag(weighted_cov(Theta[i2,,drop=FALSE], W2))
  sqrt(1 + ((mu1 - mu2)^2) / pmax(0.5*(s21 + s22), 1e-16)) |> max(na.rm=TRUE)
}
cess_fraction <- function(w, loglik, delta, trim = 0) {
  w <- pmax(w, 0); w <- w / sum(w)
  if (any(!is.finite(loglik)) || length(unique(loglik)) == 1) return(0.5)
  x <- loglik - max(loglik)
  if (trim > 0 && length(x) > 20) {
    qs <- quantile(x, probs = c(trim, 1 - trim), names = FALSE, type = 7)
    keep <- (x >= qs[1] & x <= qs[2]); x <- x[keep]; w <- w[keep]; w <- w / sum(w)
  }
  dx <- delta * x; if (any(dx < -500) || any(dx > 500)) return(0.1)
  m1 <- max(dx); a1 <- m1 + log(sum(w * exp(dx - m1)))
  m2 <- max(2*dx); a2 <- m2 + log(sum(w * exp(2*dx - m2)))
  out <- exp(2*a1 - a2); if (!is.finite(out) || out <= 0 || out > 1) return(0.1); out
}
cess_target_at_lambda <- function(lambda) 0.99 - 0.09 * exp(-3 * lambda)
next_lambda_via_cess <- function(w, loglik, lambda, target = 0.9,
                                 dalpha_floor_factor = 0.02, dalpha_cap = 0.25, tol = 1e-4,
                                 trim_early = 0.01, early_thr = 0.20) {
  if (lambda >= 1) return(1)
  rem <- 1 - lambda
  ll_range <- max(loglik) - min(loglik)
  ll_has_extreme <- ll_range > 100 || any(!is.finite(loglik))
  if (ll_has_extreme) {
    safe_delta <- min(dalpha_cap * 0.5, max(dalpha_floor_factor * rem, tol))
    cat(sprintf("    [SAFEGUARD] Extreme likelihoods detected (range=%.1f) -> conservative step\n", ll_range))
    return(min(1.0, lambda + safe_delta))
  }
  trim <- if (lambda < early_thr) trim_early else 0
  f <- function(d) cess_fraction(w, loglik, d, trim = trim) - target
  fr <- f(rem); if (!is.finite(fr) || fr >= 0) return(1.0)
  f0 <- f(0); if (!is.finite(f0) || f0 < 0) return(min(1.0, lambda + min(dalpha_cap, max(dalpha_floor_factor * rem, tol))))
  sol <- tryCatch(uniroot(f, c(0, rem), tol = tol), error = function(e) NULL)
  if (is.null(sol)) return(min(1.0, lambda + min(dalpha_cap, max(dalpha_floor_factor * rem, tol))))
  delta <- max(sol$root, dalpha_floor_factor * rem); delta <- min(delta, dalpha_cap)
  max_safe_jump <- if (lambda < 0.1) 0.15 else if (lambda < 0.5) 0.25 else 0.5
  if (delta > max_safe_jump) { cat(sprintf("    [SAFEGUARD] Capping huge jump: %.4f -> %.4f\n", sol$root, delta)); delta <- max_safe_jump }
  lambda + delta
}
dalpha_floor_eff <- function(lambda, base = 0.02) { if (lambda < 0.30) 0.04 else if (lambda < 0.70) 0.03 else base }

grow_particles_keep <- function(Z, w, lambda, Tmap, mu_ref, prior_L,
                                loglik_fn, data, mix_hist, mix_elite,
                                growth = 0.25, M_cap = 8000L, param_names = NULL,
                                indep_t_df = 5) {
  M <- nrow(Z); if (M >= M_cap) return(list(Z=Z, w=w, loglik=NULL, lpz=NULL, Theta=NULL, grew=FALSE))
  M_add <- min(max(1L, round(M * growth)), M_cap - M); if (M_add <= 0) return(list(Z=Z, w=w, loglik=NULL, lpz=NULL, Theta=NULL, grew=FALSE))
  if (!is.null(mix_hist)) {
    Znew <- sample_gmm_Z_qmc(M_add, mix_hist$meansZ, mix_hist$cache); lq <- gmm_logpdf_Z_vec(Znew, mix_hist$meansZ, mix_hist$cache)
  } else if (!is.null(mix_elite)) {
    Znew <- rmvt_mixture_Z_qmc(M_add, mix_elite$meansZ, mix_elite$cache, nu = indep_t_df)
    lq_n <- gmm_logpdf_Z_vec(Znew, mix_elite$meansZ, mix_elite$cache)
    lq_t <- dmvt_mixture_logpdf_Z_vec(Znew, mix_elite$meansZ, mix_elite$cache, nu = indep_t_df)
    m <- pmax(lq_n, lq_t); lq <- log(0.5*exp(lq_n - m) + 0.5*exp(lq_t - m)) + m
  } else return(list(Z=Z, w=w, loglik=NULL, lpz=NULL, Theta=NULL, grew=FALSE))
  Thetan <- Tmap$inv(Znew); if (!is.null(param_names)) colnames(Thetan) <- param_names
  lln <- loglik_fn(Thetan, data)
  lp_pr <- dmvnorm_chol_log(Thetan, mu_ref, prior_L)
  lpzn <- as.numeric(lp_pr - Tmap$log_jac(Thetan))
  lw_new <- lpzn + lambda * lln - lq
  lw_all <- c(log(w), lw_new); lw_all <- lw_all - logsumexp(lw_all); w_all <- exp(lw_all)
  list(Z = rbind(Z, Znew), w = w_all, loglik = lln, lpz = lpzn, Theta = Thetan, grew = TRUE)
}

# ------------- main SMC (unchanged except print labels) -------------
enhanced_smc_elite <- function(data, loglik_fn, mu_ref, Sigma_ref,
                               M = 2500L,
                               elite_quantile = 0.4,
                               cov_inflation = 3.0,
                               resample_threshold = 0.5,
                               n_mcmc_moves = 2,
                               max_rounds = 200L,
                               G_mix = 12L,
                               gamma_sharp = 1.0,
                               mix_difficulty_threshold = 0.2, # skip elite mixture when difficulty_ema below this
                               dalpha_floor_factor = 0.02,
                               dalpha_cap = 0.20,
                               refit_every = 4L,
                               rw_prob = 0.35,
                               rw_scale_init = 0.9,
                               indep_t_df = 7,
                               indep_t_prob = 0.75,
                               hist_mix_enable = TRUE,
                               hist_mix_lambda_thresh = 0.6,
                               hist_mix_prob = 0.20,
                               max_backtrack_retries = 2L,
                               ess_frac_good = 0.40,
                               ess_frac_min = 0.20,
                               ess_frac_repopulate = 0.08,
                               accept_min = 0.05,
                               split_rhat_max = 1.20,
                               B_logz_boot = 50L,
                               allow_adaptive_M = TRUE,
                               M_cap = 8000L,
                               growth_frac = 0.25,
                               mtm_K = 1L,
                               seed = 123) {
  set.seed(seed); param_names <- names(mu_ref)
  prior_L <- tryCatch(chol(Sigma_ref), error = function(e) chol(Matrix::nearPD(Sigma_ref)$mat))

  cat("Stage 1: Initial sampling from prior...\n")
  Theta <- mvtnorm::rmvnorm(M, mu_ref, Sigma_ref); colnames(Theta) <- param_names
  loglik <- loglik_fn(Theta, data)
  cat(sprintf("  Initial loglik range: [%.2f, %.2f]\n", min(loglik), max(loglik)))

  cat("Building initial transport map (uniform weights)...\n")
  Tmap <- fit_copula_transform(Theta, rep(1/M, M), corr_shrink = 0.75, slope_floor = 1e-5, slope_cap = 50, tail = 1e-3)
  Z <- Tmap$fwd(Theta)
  lp_prior <- dmvnorm_chol_log(Theta, mu_ref, prior_L)
  lpz <- as.numeric(lp_prior - Tmap$log_jac(Theta))

  w <- rep(1/M, M); lambda <- 0; round <- 0L
  acc_hist <- c(); lambda_hist <- c(0); rw_scale <- rw_scale_init
  log_evidence <- 0.0; logZ_se_hist <- numeric(0)
  elite_history <- list(); last_elite_mix <- NULL

  adaptive_settings <- list(dalpha_cap=dalpha_cap, n_mcmc_moves=n_mcmc_moves, refit_every=refit_every,
                            G_mix=G_mix, cov_inflation=cov_inflation, rw_scale=rw_scale_init,
                            resample_threshold=resample_threshold, indep_t_prob=indep_t_prob,
                            indep_t_df=indep_t_df, mtm_K=mtm_K)
  base_settings <- adaptive_settings
  difficulty_ema <- NULL; difficulty_history <- list()
  esjd_rw_ema <- NULL; esjd_id_ema <- NULL; esjd_total_ema <- NULL; target_esjd_pp <- 0.6
  stuck_count <- 0L; last_lambda <- 0; force_refit_next <- FALSE

  emergency_restart <- function() {
    cat("EMERGENCY RESTART: Algorithm completely stuck - resampling from prior\n")
    Theta_new <- mvtnorm::rmvnorm(M, mu_ref, Sigma_ref); colnames(Theta_new) <- param_names
    loglik_new <- loglik_fn(Theta_new, data)
    ll_med <- median(loglik_new); ll_mad <- mad(loglik_new)
    reasonable <- abs(loglik_new - ll_med) < 5 * ll_mad
    if (sum(reasonable) < M/2) { reasonable <- order(loglik_new, decreasing = TRUE)[1:(M/2)]; reasonable <- seq_len(M) %in% reasonable }
    Theta <<- Theta_new[reasonable, ]; loglik <<- loglik_new[reasonable]
    if (nrow(Theta) < M) {
      n_pad <- M - nrow(Theta); best_idx <- which.max(loglik)
      for (i in 1:n_pad) {
        Theta <<- rbind(Theta, Theta[best_idx, ] + rnorm(ncol(Theta), 0, 0.1))
        loglik <<- c(loglik, loglik_fn(Theta[nrow(Theta), , drop = FALSE], data))
      }
    }
    w <<- rep(1/nrow(Theta), nrow(Theta)); lambda <<- 0.05
    Tmap <<- fit_copula_transform(Theta, w, corr_shrink = 0.75, slope_floor = 1e-5, slope_cap = 50, tail = 1e-3)
    Z <<- Tmap$fwd(Theta); lp_pr <<- dmvnorm_chol_log(Theta, mu_ref, prior_L)
    lpz <<- as.numeric(lp_pr - Tmap$log_jac(Theta)); stuck_count <<- 0L
  }

  repopulate_from_mixture <- function(M_add = ceiling(0.5 * nrow(Z))) {
    if (M_add <= 0) return(invisible(NULL))
    cat(sprintf("  Repopulating with %d new particles from mixtures...\n", M_add))
    use_hist <- length(elite_history) > 0
    if (use_hist) {
      hm <- combine_elite_mixtures(elite_history, weights = c(0.6, 0.3, 0.1), housekeeping = TRUE, min_G_keep = 2, verbose = FALSE)
      if (!is.null(hm)) {
        Znew <- sample_gmm_Z_qmc(M_add, hm$meansZ, hm$cache); lq <- gmm_logpdf_Z_vec(Znew, hm$meansZ, hm$cache)
      } else use_hist <- FALSE
    }
    if (!use_hist) {
      em <- last_elite_mix; if (is.null(em) || !length(em$meansZ)) return(invisible(NULL))
      Znew <- rmvt_mixture_Z_qmc(M_add, em$meansZ, em$cache, nu = 5)
      lq_n <- gmm_logpdf_Z_vec(Znew, em$meansZ, em$cache); lq_t <- dmvt_mixture_logpdf_Z_vec(Znew, em$meansZ, em$cache, nu = 5)
      m <- pmax(lq_n, lq_t); lq <- log(0.5*exp(lq_n - m) + 0.5*exp(lq_t - m)) + m
    }
    Thetan <- Tmap$inv(Znew); colnames(Thetan) <- param_names
    lln <- loglik_fn(Thetan, data); lp_pr <- dmvnorm_chol_log(Thetan, mu_ref, prior_L)
    lpzn <- as.numeric(lp_pr - Tmap$log_jac(Thetan)); lw_new <- lpzn + lambda * lln - lq
    Theta2 <- rbind(Theta, Thetan); Z2 <- rbind(Z, Znew); loglik2 <- c(loglik, lln)
    lw2 <- c(log(w), lw_new); lw2 <- lw2 - logsumexp(lw2); w2 <- exp(lw2)
    idx <- systematic_resample(w2)
    Theta <<- Theta2[idx,,drop=FALSE]; Z <<- Z2[idx,,drop=FALSE]; loglik<<- loglik2[idx]
    w <<- rep(1/nrow(Theta), nrow(Theta))
    lp_pr <- dmvnorm_chol_log(Theta, mu_ref, prior_L)
    lpz <<- as.numeric(lp_pr - Tmap$log_jac(Theta))
  }

  while (lambda < 1 - 1e-12 && round < max_rounds) {
    round <- round + 1L
    cat(sprintf("\nRound %d: lambda %.3f -> ", round, lambda))
    cess_target <- cess_target_at_lambda(lambda)
    dalpha_floor <- dalpha_floor_eff(lambda, base = dalpha_floor_factor)
    next_lambda <- next_lambda_via_cess(w, loglik, lambda, target = cess_target,
                                        dalpha_floor_factor = dalpha_floor,
                                        dalpha_cap = adaptive_settings$dalpha_cap,
                                        trim_early = 0.01, early_thr = 0.20)
    dalpha <- next_lambda - lambda
    cat(sprintf("%.3f (dalpha=%.4f, CESS target=%.3f)\n", next_lambda, dalpha, cess_target))

    success <- FALSE; lambda_try <- next_lambda; snapshot <- NULL
    for (retry in 0:max_backtrack_retries) {
      snapshot <- list(Theta=Theta, Z=Z, w=w, loglik=loglik, lpz=lpz,
                       Tmap=Tmap, rw_scale=rw_scale, last_elite_mix=last_elite_mix,
                       elite_history=elite_history, log_evidence=log_evidence,
                       adaptive_settings=adaptive_settings, difficulty_ema=difficulty_ema,
                       logZ_se_hist=logZ_se_hist, rw_prob=rw_prob,
                       esjd_rw_ema=esjd_rw_ema, esjd_id_ema=esjd_id_ema, esjd_total_ema=esjd_total_ema)

      mll <- max(loglik); lw_inc <- (lambda_try - lambda) * (loglik - mll)
      reweight <- exp(lw_inc); c_norm <- sum(w * reweight)
      w_new <- w * reweight; w_new <- w_new / sum(w_new)

      if (B_logz_boot > 0) {
        prob <- w / sum(w); inc_samples <- numeric(B_logz_boot)
        for (b in 1:B_logz_boot) {
          idxb <- sample.int(length(prob), size = length(prob), replace = TRUE, prob = prob)
          cb <- mean(reweight[idxb]); inc_samples[b] <- (lambda_try - lambda) * mll + log(cb)
        }
        se_inc <- sd(inc_samples)
      } else se_inc <- NA_real_

      ess <- ESS(w_new); ess_frac <- ess/length(w_new)
      log_inc <- (lambda_try - lambda) * mll + log(c_norm)
      log_evidence_try <- log_evidence + log_inc
      cat(sprintf("  ESS(pre-move)=%.3f | logZ += %.4f (SE≈%.4f) -> %.4f\n",
                  ess_frac, log_inc, ifelse(is.na(se_inc), 0, se_inc), log_evidence_try))

      if (ess_frac < adaptive_settings$resample_threshold) {
        cat("  Resampling...\n")
        idx <- systematic_resample(w_new)
        Theta  <- Theta[idx, , drop = FALSE]; Z <- Z[idx, , drop = FALSE]
        loglik <- loglik[idx]; lpz <- lpz[idx]; w_new <- rep(1/length(w_new), length(w_new))
        cat("  Pre-move jitter (RW only)...\n")
        jitter <- mcmc_moves_z_mix_batched(
          Z, loglik, lpz, Tmap, lambda_try,
          mu_ref, prior_L, w_new,
          elite_mix = list(meansZ = list(), covsZ = list(), wZ = 1, cache = list(Ls = list(), logw = 0)),
          hist_mix = NULL,
          data, loglik_fn, n_moves = 1,
          rw_prob = 1.0, rw_scale = rw_scale,
          indep_t_df = adaptive_settings$indep_t_df, indep_t_prob = 0.0,
          mtm_K = 1L,
          seed = seed + round * 31,
          param_names = param_names
        )
        Z <- jitter$Z; loglik <- jitter$loglik; lpz <- jitter$lpz
        Theta <- Tmap$inv(Z); colnames(Theta) <- param_names
      }
      w <- w_new; log_evidence <- log_evidence_try
      if (!is.na(se_inc)) logZ_se_hist <- c(logZ_se_hist, se_inc)

      need_refit_base <- force_refit_next ||
        ((round %% adaptive_settings$refit_every) == 0L) ||
        (length(acc_hist) && tail(acc_hist, 1) < 0.15)
      force_refit_next <- FALSE
      small_step <- ((lambda_try - lambda) < 0.08); concentrated <- (max(w) > 0.10)
      if (need_refit_base && small_step && concentrated) {
        cat("  Rebuilding transport map with CURRENT weights...\n")
        Tmap <- fit_copula_transform(Theta, w, corr_shrink = 0.75, slope_floor = 1e-5, slope_cap = 50, tail = 1e-3)
        Z <- Tmap$fwd(Theta); lp_prior <- dmvnorm_chol_log(Theta, mu_ref, prior_L)
        lpz <- as.numeric(lp_prior - Tmap$log_jac(Theta))
        sub <- sample.int(nrow(Theta), min(128L, nrow(Theta)))
        cyc_err <- max(abs(Tmap$inv(Z[sub,,drop=FALSE]) - Theta[sub,,drop=FALSE]))
        if (!is.finite(cyc_err) || cyc_err > 1e-5) {
          cat(sprintf("    invertibility warn (%.2e) -> refit(uniform)\n", cyc_err))
          Tmap <- fit_copula_transform(Theta, rep(1/nrow(Theta), nrow(Theta)), corr_shrink = 0.75, slope_floor = 1e-5, slope_cap = 50, tail = 1e-3)
          Z <- Tmap$fwd(Theta); lp_prior <- dmvnorm_chol_log(Theta, mu_ref, prior_L)
          lpz <- as.numeric(lp_prior - Tmap$log_jac(Theta))
        }
      } else if (need_refit_base) cat("  Skipping transport refit (coarse step or diffuse weights)...\n")

      G_eff <- max(2L, min(adaptive_settings$G_mix, 2L + floor(6 * lambda)))
      minG_keep <- if (lambda_try < 0.10) 4 else 2
      elite_mix <- if (is.null(last_elite_mix))
        list(meansZ = list(), covsZ = list(), wZ = 1, cache = list(Ls = list(), logw = 0))
      else last_elite_mix
      difficulty_metric <- if (is.null(difficulty_ema)) (1 - ess_frac) else difficulty_ema
      if (difficulty_metric >= mix_difficulty_threshold) {
        cat("  Fitting elite mixture...\n")
        w_fit <- (w^gamma_sharp); w_fit <- w_fit / sum(w_fit)
        elite_q_eff <- max(0.20, elite_quantile * (1 - 0.5 * lambda))
        do_housekeeping <- (lambda_try > 0.06) || (round > 2)
        elite_mix <- fit_elite_mixture_Z(
          Z, w_fit, elite_quantile = elite_q_eff, G = G_eff,
          cov_inflation = adaptive_settings$cov_inflation, warm_start_mixZ = last_elite_mix,
          em_itmax = 8, housekeeping = do_housekeeping, min_G_keep = minG_keep, verbose = FALSE
        )
        last_elite_mix <- elite_mix
      } else {
        cat(sprintf("  Skipping elite mixture fit (difficulty=%.2f < %.2f)\n",
                    difficulty_metric, mix_difficulty_threshold))
      }

      use_hist_mix <- hist_mix_enable && (lambda >= hist_mix_lambda_thresh) && length(elite_history) > 0
      hist_mix <- NULL
      if (use_hist_mix) {
        ramp <- (lambda - hist_mix_lambda_thresh) / max(1 - hist_mix_lambda_thresh, 1e-8)
        ramp <- pmin(pmax(ramp, 0), 1); hist_prob_eff <- (0.5 * hist_mix_prob) + (0.5 * hist_mix_prob) * ramp
        combined_mix <- combine_elite_mixtures(elite_history, weights = c(0.6, 0.3, 0.1),
                                               housekeeping = TRUE, min_G_keep = minG_keep, verbose = FALSE)
        if (!is.null(combined_mix)) hist_mix <- list(mix = combined_mix, prob = hist_prob_eff)
      }

      rw_prob <- .clamp(rw_prob, 0.10, 0.90); rw_prob_eff <- max(0.10, rw_prob * (1 - 0.5 * lambda))
      cat(sprintf("  Rejuvenation: %d move(s) [rw_prob=%.2f, rw_scale=%.3f, t_df=%d, G=%d, MTM-K=%d]...\n",
                  adaptive_settings$n_mcmc_moves, rw_prob_eff, rw_scale, adaptive_settings$indep_t_df, G_eff, adaptive_settings$mtm_K))
      move <- mcmc_moves_z_mix_batched(Z, loglik, lpz, Tmap, lambda_try,
                                       mu_ref, prior_L, w,
                                       elite_mix, hist_mix,
                                       data, loglik_fn, n_moves = adaptive_settings$n_mcmc_moves,
                                       rw_prob = rw_prob_eff, rw_scale = rw_scale,
                                       indep_t_df = adaptive_settings$indep_t_df,
                                       indep_t_prob = adaptive_settings$indep_t_prob,
                                       mtm_K = adaptive_settings$mtm_K,
                                       seed = seed + round * 97,
                                       param_names = param_names)
      Z <- move$Z; loglik <- move$loglik; lpz <- move$lpz
      Theta <- Tmap$inv(Z); colnames(Theta) <- param_names
      acc_hist <- c(acc_hist, move$accept_rate)

      ess_post_tmp <- ESS(w)/length(w)
      if (move$accept_rate < 0.005 && ess_post_tmp > 0.95) {
        cat("  ⚠️ Moves stalled (acc≈0 & ESS≈1) -> forcing transport refit (uniform) and jitter.\n")
        Tmap <- fit_copula_transform(Theta, rep(1/nrow(Theta), nrow(Theta)), corr_shrink = 0.75, slope_floor = 1e-5, slope_cap = 50, tail = 1e-3)
        Z <- Tmap$fwd(Theta); lp_prior <- dmvnorm_chol_log(Theta, mu_ref, prior_L)
        lpz <- as.numeric(lp_prior - Tmap$log_jac(Theta))
        jitter <- mcmc_moves_z_mix_batched(
          Z, loglik, lpz, Tmap, lambda_try,
          mu_ref, prior_L, w,
          elite_mix = list(meansZ = list(), covsZ = list(), wZ = 1, cache = list(Ls = list(), logw = 0)),
          hist_mix = NULL, data, loglik_fn, n_moves = 1,
          rw_prob = 1.0, rw_scale = max(0.1, rw_scale*0.8),
          indep_t_df = adaptive_settings$indep_t_df, indep_t_prob = 0.0, mtm_K = 1L,
          seed = seed + round * 313, param_names = param_names
        )
        Z <- jitter$Z; loglik <- jitter$loglik; lpz <- jitter$lpz
        Theta <- Tmap$inv(Z); colnames(Theta) <- param_names
      }

      esjd_rw_pp    <- if (move$prop_rw>0) move$esjd_rw / move$prop_rw else 0
      esjd_ind_pp   <- if (move$prop_indep>0) move$esjd_indep / move$prop_indep else 0
      esjd_total_pp <- (move$esjd_rw + move$esjd_indep) / nrow(Z)
      esjd_rw_ema    <- ema_update(esjd_rw_ema, esjd_rw_pp, beta = 0.3)
      esjd_id_ema    <- ema_update(esjd_id_ema, esjd_ind_pp, beta = 0.3)
      esjd_total_ema <- ema_update(esjd_total_ema, esjd_total_pp, beta = 0.3)
      cat(sprintf("    accepts: total=%.3f | rw=%.3f | indep=%.3f | ESJD(pp): rw=%.3f, indep=%.3f, total=%.3f\n",
                  move$accept_rate, move$rw_accept_rate, move$indep_accept_rate,
                  esjd_rw_pp, esjd_ind_pp, esjd_total_pp))

      ## sharpen adaptation of RW scale (faster shrink / grow)
      if (move$rw_accept_rate < 0.15) rw_scale <- rw_scale * 0.85
      else if (move$rw_accept_rate > 0.40) rw_scale <- rw_scale * 1.12

      ess_post <- ESS(w)/length(w)
      dm <- difficulty_metrics(loglik, w, move$accept_rate, ess_post, lambda, (lambda_try - lambda), round)
      difficulty_ema <- ema_update(difficulty_ema, dm$index, beta = 0.3)
      ## ------------------------------------------------------------------
      ##  ❐  STALL‑DETECTION / EARLY‑RESCUE  --------------------------------
      ## ------------------------------------------------------------------

      stalled_moves <- (move$accept_rate < 0.01 && esjd_total_pp < 1e-4)
      near_uniform  <- (ESS(w) / length(w) > 0.90)

      if (stalled_moves && near_uniform) {
        ## 1)  make independence kernel more aggressive
        adaptive_settings$mtm_K     <- max(3L, adaptive_settings$mtm_K)
        adaptive_settings$indep_t_prob <- min(0.95, adaptive_settings$indep_t_prob + 0.15)
        adaptive_settings$indep_t_df   <- max(4L,  adaptive_settings$indep_t_df - 1L)

        ## 2)  force a *weighted* transport refit next round
        force_refit_next <- TRUE

        ## 3)  shrink RW step immediately
        rw_scale <- max(0.08, rw_scale * 0.75)
      }

      ## retain original difficulty‑based MTM, but relax threshold
      if (!is.null(difficulty_ema)) {
        if (difficulty_ema > 0.60) adaptive_settings$mtm_K <- max(3L, adaptive_settings$mtm_K)
        else if (difficulty_ema < 0.50 && adaptive_settings$mtm_K > 1L)
          adaptive_settings$mtm_K <- adaptive_settings$mtm_K - 1L
      }
      difficulty_history[[round]] <- list(index = difficulty_ema, components = dm$components, raw = dm$raw)

      targets <- knob_targets_from_difficulty(base_settings, difficulty_ema)
      adaptive_settings <- smooth_update_knobs(adaptive_settings, targets, eta = 0.25)
      denom_esjd <- (esjd_rw_ema + esjd_id_ema + 1e-12)
      p_rw_target <- esjd_rw_ema / denom_esjd
      rw_prob <- .clamp(rw_prob + 0.20 * (p_rw_target - rw_prob), 0.10, 0.90)
      if (!is.null(esjd_total_ema)) {
        if (esjd_total_ema < 0.7 * target_esjd_pp) adaptive_settings$n_mcmc_moves <- min(adaptive_settings$n_mcmc_moves + 1L, 6L)
        else if (esjd_total_ema > 1.5 * target_esjd_pp) adaptive_settings$n_mcmc_moves <- max(adaptive_settings$n_mcmc_moves - 1L, 1L)
      }
      if (!is.null(esjd_id_ema)) {
        if (esjd_id_ema < 0.3) adaptive_settings$G_mix <- min(adaptive_settings$G_mix + 1L, 32L)
        else if (esjd_id_ema > 1.2) adaptive_settings$G_mix <- max(adaptive_settings$G_mix - 1L, 6L)
      }
      if (!is.null(difficulty_ema)) {
        if (difficulty_ema > 0.75 && adaptive_settings$mtm_K < 3) adaptive_settings$mtm_K <- adaptive_settings$mtm_K + 1L
        else if (difficulty_ema < 0.55 && adaptive_settings$mtm_K > 1) adaptive_settings$mtm_K <- adaptive_settings$mtm_K - 1L
      }

      n_mcmc_moves       <- adaptive_settings$n_mcmc_moves
      refit_every        <- adaptive_settings$refit_every
      G_mix              <- adaptive_settings$G_mix
      cov_inflation      <- adaptive_settings$cov_inflation
      resample_threshold <- adaptive_settings$resample_threshold
      dalpha_cap         <- adaptive_settings$dalpha_cap
      rw_scale           <- adaptive_settings$rw_scale
      indep_t_prob       <- adaptive_settings$indep_t_prob
      indep_t_df         <- adaptive_settings$indep_t_df

      if (allow_adaptive_M &&
          (difficulty_ema > 0.75 &&           # hard
          (lambda_try <= .3 &&              # still early
          (move$accept_rate < 0.02 || esjd_total_pp   < 1e-3))) &&
          nrow(Z) < M_cap) {
        cat(sprintf("  ⤴️ Growing particle set by %.0f%% to fight early stagnation\n", 100*growth_frac))
        hist_comb <- if (!is.null(hist_mix)) hist_mix$mix else
          if (length(elite_history)>0) combine_elite_mixtures(elite_history, housekeeping=TRUE, min_G_keep = minG_keep) else NULL
        gp <- grow_particles_keep(Z, w, lambda_try, Tmap, mu_ref, prior_L, loglik_fn, data,
                                  mix_hist = hist_comb, mix_elite = elite_mix,
                                  growth = growth_frac, M_cap = M_cap,
                                  param_names = param_names, indep_t_df = indep_t_df)
        if (gp$grew) {
          Z <- gp$Z; w <- gp$w
          if (!is.null(gp$loglik)) loglik <- c(loglik, gp$loglik)
          if (!is.null(gp$lpz))    lpz    <- c(lpz, gp$lpz)
          if (!is.null(gp$Theta))  Theta  <- rbind(Theta, gp$Theta)
          cat(sprintf("    New M = %d\n", nrow(Z)))
        }
      }

      M_suggest <- suggest_particles_smooth(difficulty_ema, nrow(Z), M_cap = M_cap, d_start = 0.7)
      if (M_suggest > nrow(Z)) cat(sprintf("  💡 Suggest increasing M from %d to %d (difficulty=%.2f)\n", nrow(Z), M_suggest, difficulty_ema))

      splitR <- split_rhat_weighted(Theta, w)
      cat(sprintf("  Health: ESS=%.3f | splitR=%.3f | difficulty=%.2f | rw_prob=%.2f | moves=%d | G=%d\n",
                  ess_post, splitR, difficulty_ema, rw_prob, adaptive_settings$n_mcmc_moves, adaptive_settings$G_mix))

      ## tolerate lower acceptance when tempered likelihood is rough
      accept_min_dynamic <- if (lambda < 0.10) 0.03 else if (lambda < 0.25) 0.02 else 0.01
      bad <- (ess_post < ess_frac_min) ||
        (move$accept_rate < accept_min_dynamic) ||
        (splitR > split_rhat_max)
      if (!bad) {
        success <- TRUE; lambda <- lambda_try; lambda_hist <- c(lambda_hist, lambda)
        if (lambda > 0.5 && !is.null(elite_mix)) {
          elite_history[[length(elite_history) + 1L]] <- elite_mix
          if (length(elite_history) > 3L) elite_history <- tail(elite_history, 3L)
        }
        break
      }
      cat("  Step failed convergence checks -> backtracking...\n")
      Theta <- snapshot$Theta; Z <- snapshot$Z; w <- snapshot$w
      loglik <- snapshot$loglik; lpz <- snapshot$lpz; Tmap <- snapshot$Tmap
      rw_scale <- snapshot$rw_scale; last_elite_mix <- snapshot$last_elite_mix
      elite_history <- snapshot$elite_history; log_evidence <- snapshot$log_evidence
      adaptive_settings <- snapshot$adaptive_settings; difficulty_ema <- snapshot$difficulty_ema
      logZ_se_hist <- snapshot$logZ_se_hist; rw_prob <- snapshot$rw_prob
      esjd_rw_ema <- snapshot$esjd_rw_ema; esjd_id_ema <- snapshot$esjd_id_ema; esjd_total_ema <- snapshot$esjd_total_ema
      if (retry < max_backtrack_retries) { lambda_try <- lambda + 0.5*(lambda_try - lambda); cat(sprintf("  Retrying with lambda=%.4f\n", lambda_try)) }
      else break
    }

    if (!success) {
      difficulty_ema <- ema_update(difficulty_ema, 1.0, beta = 0.5)
      targets <- knob_targets_from_difficulty(base_settings, difficulty_ema)
      adaptive_settings <- smooth_update_knobs(adaptive_settings, targets, eta = 0.5)
      ess_curr <- ESS(w)/length(w)
      if (ess_curr < ess_frac_repopulate) repopulate_from_mixture()
    } else {
      delta_lambda <- lambda - last_lambda; if (delta_lambda < 0.001) stuck_count <- stuck_count + 1L else stuck_count <- 0L
      last_lambda <- lambda; if (stuck_count >= 8) { emergency_restart(); round <- round + 1L; next }
    }
  }

  cat(sprintf("\nCompleted in %d rounds, final lambda = %.3f\n", round, lambda))
  if (length(logZ_se_hist)) {
    cum_se <- sqrt(sum(logZ_se_hist^2, na.rm=TRUE))
    cat(sprintf("Log-evidence SE (approx, bootstrap per step, independence assumed): %.4f\n", cum_se))
  }
  if (length(difficulty_history) > 0) {
    idxs <- sapply(difficulty_history, function(x) x$index)
    cat("\nSmooth difficulty summary:\n")
    cat(sprintf("  mean=%.2f | median=%.2f | max=%.2f\n", mean(idxs, na.rm=TRUE), median(idxs, na.rm=TRUE), max(idxs, na.rm=TRUE)))
  }

  list(
    Theta = Theta, Z = Z, loglik = loglik, w = w, transport = Tmap,
    final_lambda = lambda, log_evidence = log_evidence,
    meta = list(
      rounds = round, ess = ESS(w), acc_hist = acc_hist, lambda_hist = lambda_hist,
      rw_scale_final = rw_scale, adaptive_settings_final = adaptive_settings,
      difficulty_history = difficulty_history, logZ_se_hist = logZ_se_hist,
      esjd = list(rw_ema = esjd_rw_ema, indep_ema = esjd_id_ema, total_ema = esjd_total_ema)
    ),
    elite_history = elite_history
  )
}

# ------------- auto_smc (unchanged) -------------
auto_smc <- function(data, loglik_fn, mu_ref, Sigma_ref,
                     M = NULL, seed = 123, max_rounds = 100L, verbose = TRUE) {
  d <- length(mu_ref)
  N <- if (is.list(data) && "y" %in% names(data)) length(data$y)
  else if (is.list(data) && "N" %in% names(data)) data$N else 200
  M_from_d              <- function(d) round(.clamp(2000 + 400 * d^0.8, 1500, 10000))
  dalpha_cap_from_d     <- function(d) .clamp(0.25 * (4/(d + 4)), 0.06, 0.25)
  n_moves_from_d        <- function(d) round(.clamp(2 + log1p(d)/log(2), 2, 6))
  G_from_d              <- function(d) round(.clamp(6 + 1.2 * d, 6, 24))
  rw_scale_from_d       <- function(d) .clamp(0.55 / sqrt(1 + d/3), 0.12, 0.60)
  refit_every_from_d    <- function(d) max(1L, min(3L, round(3 - 0.6 * log1p(d))))
  cov_inflation_from_d  <- function(d) .clamp(2.3 + 0.05 * d, 2.0, 4.0)
  if (is.null(M)) M <- M_from_d(d)
  settings <- list(
    M = M, max_rounds = max_rounds, dalpha_floor_factor = 0.02,
    dalpha_cap = dalpha_cap_from_d(d), n_mcmc_moves = n_moves_from_d(d),
    G_mix = G_from_d(d), rw_prob = 0.35, rw_scale_init = rw_scale_from_d(d),
    indep_t_prob = 0.75, indep_t_df = 5, resample_threshold = 0.5,
    refit_every = refit_every_from_d(d), cov_inflation = cov_inflation_from_d(d),
    elite_quantile = 0.4, mix_difficulty_threshold = 0.2,
    accept_min = 0.05, ess_frac_min = 0.15,
    max_backtrack_retries = 1L, allow_adaptive_M = TRUE, M_cap = 8000L,
    growth_frac = 0.25, mtm_K = 1L, seed = seed
  )
  if (verbose) {
    cat("=== AUTOMATIC SMC SETUP ===\n")
    cat(sprintf("Problem: %d parameters, ~%d observations\n", d, N))
    cat(sprintf("Dim-based init: M=%d | dalpha_cap=%.3f | moves=%d | G=%d | rw_scale=%.3f | refit_every=%d | cov_infl=%.2f\n",
                settings$M, settings$dalpha_cap, settings$n_mcmc_moves, settings$G_mix,
                settings$rw_scale_init, settings$refit_every, settings$cov_inflation))
    cat("Running adaptive SMC...\n\n")
  }
  result <- tryCatch({
    do.call(enhanced_smc_elite, c(
      list(data = data, loglik_fn = loglik_fn, mu_ref = mu_ref, Sigma_ref = Sigma_ref),
      settings
    ))
  }, error = function(e) {
    if (verbose) {
      cat("\n🚨 SMC FAILED WITH CURRENT SEED\n")
      cat("This suggests an extremely pathological likelihood landscape.\n")
      cat("RECOMMENDATIONS:\n")
      cat("• Try different seeds: c(456, 789, 101112, 131415)\n")
      cat("• Check data generation and likelihood function\n")
      cat("• Consider rescaling/transforming parameters\n")
      cat(sprintf("• Error: %s\n", e$message))
    }
    NULL
  })
  if (is.null(result)) return(NULL)
  if (verbose) {
    rounds <- result$meta$rounds
    if (rounds > 50) {
      cat("\n⚠️ SLOW CONVERGENCE DETECTED\n")
      cat(sprintf("Took %d rounds - this seed may be problematic.\n", rounds))
      cat("RECOMMENDATIONS:\n")
      cat("• Try different seeds for more reliable results\n")
      cat("• Consider increasing M or adjusting priors\n")
    }
    if (length(result$meta$difficulty_history) > 0) {
      idxs <- sapply(result$meta$difficulty_history, function(x) x$index)
      mean_d <- mean(idxs, na.rm=TRUE)
      if (mean_d > 0.6) {
        cat("\n=== POST-RUN RECOMMENDATIONS ===\n")
        cat("This was a challenging run (high mean difficulty index).\n")
        cat(sprintf("• Increase M to %d\n", min(settings$M * 1.5, 10000)))
        cat("• Run multiple chains with different seeds\n")
        cat("• Check model identifiability / priors\n")
      }
    }
  }
  result
}
