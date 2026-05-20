# ========================================================================
# Per-group / local SMC sampler
# - Runs one group at a time and returns the local particle system,
#   transport map, and proposal objects used by the EMC example runner
#
# This implementation provides:
# - Adaptive SMC with temperature schedule based on conditional ESS (CESS)
# - Gaussian-copula transport for parameter space whitening and tail protection
# - Rejuvenation using mixture of elite and historical components
# - Vectorized MCMC moves
# - Robust safeguards including restart, repopulation, and covariance regularization
# ========================================================================

suppressPackageStartupMessages({
  library(mvtnorm)
  library(qrng)
  library(Matrix)
  library(matrixStats)
})

# Source the shared local SMC utilities.
source("smc_core.R")
if (!exists("normalize_reference_prior", mode = "function") ||
    !exists("reference_prior_logpdf", mode = "function") ||
    !exists("reference_prior_sample", mode = "function")) {
  source("reference_priors.R")
}


# ----------------------------- transport utilities -----------------------------

# Create tail-safe CDF transformation bundle for marginal distributions
# Handles degenerate dimensions and provides robust tail behavior
# Args:
#   ux: x-values from quantile grid
#   cw: corresponding CDF values
#   Fraw: raw spline-based CDF function
#   slope_raw: raw derivative function
# Returns: list with Fhat (CDF), slope (derivative), and Qhat (quantile) functions
.make_tail_safe_cdf_bundle <- function(ux, cw) {
  stopifnot(length(ux) >= 2L, length(ux) == length(cw))
  eps <- max(1e-12, 1e-8 * diff(range(ux)))
  for (k in 2:length(ux)) if (ux[k] <= ux[k-1]) ux[k] <- ux[k-1] + eps

  xmin <- ux[1]; xmax <- ux[length(ux)]
  uL <- pmin(pmax(cw[1], 1e-12), 1 - 1e-12)
  uR <- pmin(pmax(cw[length(cw)], 1e-12), 1 - 1e-12)
  vL <- qnorm(uL); vR <- qnorm(uR)

  F_mid <- approxfun(ux, cw, method = "linear", ties = "ordered", rule = 2)
  Q_mid <- approxfun(cw, ux, method = "linear", ties = "ordered", rule = 2)

  du <- diff(cw); dx <- pmax(diff(ux), 1e-12)
  seg_slope <- du / dx
  vgrid <- qnorm(cw); dv <- diff(vgrid); seg_slope_v <- dv / dx

  g_scale <- median(abs(seg_slope_v[is.finite(seg_slope_v)])); if (!is.finite(g_scale) || g_scale<=0) g_scale <- 1
  gL <- seg_slope_v[1]; gR <- seg_slope_v[length(seg_slope_v)]
  gL <- min(max(gL, 1e-6), 1e12 * g_scale); gR <- min(max(gR, 1e-6), 1e12 * g_scale)

  slope_mid <- function(x) {
    k <- pmin(pmax(findInterval(x, ux, all.inside = TRUE), 1L), length(seg_slope))
    seg_slope[k]
  }

  Fhat <- function(x) {
    x <- as.numeric(x)
    lt <- x < xmin; rt <- x > xmax; md <- !(lt | rt)
    out <- numeric(length(x))
    if (any(md)) out[md] <- F_mid(x[md])
    if (any(lt)) out[lt] <- pnorm(vL + gL * (x[lt] - xmin))
    if (any(rt)) out[rt] <- pnorm(vR + gR * (x[rt] - xmax))
    pmin(pmax(out, 1e-12), 1 - 1e-12)
  }

  slope <- function(x) {
    x <- as.numeric(x)
    lt <- x < xmin; rt <- x > xmax; md <- !(lt | rt)
    out <- numeric(length(x))
    if (any(md)) out[md] <- slope_mid(x[md])
    if (any(lt)) { v <- vL + gL * (x[lt] - xmin); out[lt] <- dnorm(v) * gL }
    if (any(rt)) { v <- vR + gR * (x[rt] - xmax); out[rt] <- dnorm(v) * gR }
    pmax(out, 1e-300)     # guard inside log only
  }

  Qhat <- function(u) {
    u <- pmin(pmax(u, 1e-12), 1 - 1e-12)
    lt <- u < uL; rt <- u > uR; md <- !(lt | rt)
    out <- numeric(length(u))
    if (any(md)) out[md] <- Q_mid(u[md])
    if (any(lt)) out[lt] <- xmin + (qnorm(u[lt]) - vL)/gL
    if (any(rt)) out[rt] <- xmax + (qnorm(u[rt]) - vR)/gR
    out
  }

  list(Fhat = Fhat, Qhat = Qhat, slope = slope,
       F_mid = F_mid, Q_mid = Q_mid, xmin = xmin, xmax = xmax, vL = vL, vR = vR, gL = gL, gR = gR)
}

# Weighted quantile helper used by the transport marginals and residual bundles.
.weighted_quantile_linear <- function(x, w, p) {
  o <- order(x)
  x <- x[o]
  w <- pmax(w[o], 0)
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) {
    w <- rep(1 / length(x), length(x))
  } else {
    w <- w / sw
  }
  cw <- cumsum(w)
  cw <- cw / cw[length(cw)]
  approx(cw, x, xout = p, ties = "ordered", rule = 2)$y
}

.weighted_mean1 <- function(x, w) {
  w <- pmax(w, 0)
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) return(mean(x))
  sum(x * w) / sw
}

.weighted_var1 <- function(x, w, mu = NULL) {
  mu <- mu %||% .weighted_mean1(x, w)
  w <- pmax(w, 0)
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) return(stats::var(x))
  sum(w * (x - mu)^2) / sw
}

.weighted_sd1 <- function(x, w, mu = NULL, floor = 1e-6) {
  sqrt(max(.weighted_var1(x, w, mu = mu), floor))
}

.weighted_cor1 <- function(x, y, w) {
  mx <- .weighted_mean1(x, w)
  my <- .weighted_mean1(y, w)
  vx <- .weighted_var1(x, w, mu = mx)
  vy <- .weighted_var1(y, w, mu = my)
  if (!is.finite(vx) || !is.finite(vy) || vx <= 0 || vy <= 0) return(0)
  w <- pmax(w, 0)
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) return(0)
  cov_xy <- sum(w * (x - mx) * (y - my)) / sw
  cov_xy / sqrt(vx * vy)
}

.make_triangular_features <- function(Zprev,
                                      parents,
                                      centers = NULL,
                                      scales = NULL,
                                      interaction_max = 3L) {
  n <- if (is.null(Zprev)) 0L else nrow(Zprev)
  if (!length(parents)) {
    X <- matrix(1, nrow = n, ncol = 1L)
    colnames(X) <- "(Intercept)"
    return(list(X = X, centers = numeric(0), scales = numeric(0)))
  }

  Zsel <- as.matrix(Zprev[, parents, drop = FALSE])
  cols <- list()
  cols[[1L]] <- Zsel
  if (ncol(Zsel) >= 1L) {
    cols[[length(cols) + 1L]] <- Zsel^2 - 1
  }
  if (ncol(Zsel) >= 2L && ncol(Zsel) <= interaction_max) {
    cmb <- utils::combn(seq_len(ncol(Zsel)), 2L)
    for (k in seq_len(ncol(cmb))) {
      cols[[length(cols) + 1L]] <- Zsel[, cmb[1L, k]] * Zsel[, cmb[2L, k]]
    }
  }

  Xraw <- do.call(cbind, cols)
  if (is.null(centers)) {
    centers <- colMeans(Xraw)
  }
  if (is.null(scales)) {
    scales <- apply(Xraw, 2L, stats::sd)
    scales[!is.finite(scales) | scales < 1e-8] <- 1
  }
  Xstd <- sweep(Xraw, 2L, centers, "-")
  Xstd <- sweep(Xstd, 2L, scales, "/")
  X <- cbind(1, Xstd)
  colnames(X) <- c("(Intercept)", paste0("f", seq_len(ncol(Xstd))))
  list(X = X, centers = centers, scales = scales)
}

.fit_triangular_mean_model <- function(y,
                                       Zprev,
                                       w,
                                       max_parents = 4L,
                                       min_parent_cor = 0.05,
                                       ridge = 1e-3,
                                       interaction_max = 3L) {
  if (is.null(Zprev) || !ncol(Zprev)) {
    mu <- .weighted_mean1(y, w)
    return(list(
      parents = integer(0),
      centers = numeric(0),
      scales = numeric(0),
      beta = mu,
      interaction_max = interaction_max
    ))
  }

  dprev <- ncol(Zprev)
  max_keep <- max(0L, min(as.integer(max_parents), dprev))
  if (!max_keep) {
    mu <- .weighted_mean1(y, w)
    return(list(
      parents = integer(0),
      centers = numeric(0),
      scales = numeric(0),
      beta = mu,
      interaction_max = interaction_max
    ))
  }

  cors <- vapply(seq_len(dprev), function(j) abs(.weighted_cor1(y, Zprev[, j], w)), numeric(1))
  ord <- order(cors, decreasing = TRUE)
  keep <- ord[seq_len(max_keep)]
  keep <- keep[cors[keep] >= min_parent_cor]
  if (!length(keep)) {
    keep <- ord[seq_len(min(1L, dprev))]
  }

  feat <- .make_triangular_features(Zprev, keep, interaction_max = interaction_max)
  X <- feat$X
  w <- pmax(w, 0)
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) {
    w <- rep(1 / nrow(X), nrow(X))
  } else {
    w <- w / sw
  }
  WX <- X * w
  pen <- diag(ncol(X))
  pen[1L, 1L] <- 0
  ridge_eff <- as.numeric(ridge) * (1 + ncol(X) / max(nrow(X), 1))
  XtWX <- crossprod(X, WX) + ridge_eff * pen
  XtWy <- drop(crossprod(X, w * y))
  beta <- tryCatch(
    as.numeric(solve(XtWX, XtWy)),
    error = function(e) as.numeric(qr.solve(XtWX, XtWy))
  )
  list(
    parents = keep,
    centers = feat$centers,
    scales = feat$scales,
    beta = beta,
    interaction_max = interaction_max
  )
}

.eval_triangular_mean_model <- function(model, Zprev) {
  if (!length(model$parents)) {
    n <- if (is.null(Zprev)) 1L else nrow(as.matrix(Zprev))
    return(rep(model$beta[1L], n))
  }
  feat <- .make_triangular_features(
    Zprev = Zprev,
    parents = model$parents,
    centers = model$centers,
    scales = model$scales,
    interaction_max = model$interaction_max %||% 3L
  )
  drop(feat$X %*% model$beta)
}

.fit_residual_bundle <- function(r, w, pgrid) {
  rng <- diff(range(r))
  if (!is.finite(rng) || rng < 1e-10) {
    mu <- .weighted_mean1(r, w)
    sig <- .weighted_sd1(r, w, mu = mu, floor = 1e-6)
    return(list(
      Fhat = function(x) pnorm((x - mu) / sig),
      Qhat = function(u) mu + sig * qnorm(pmin(pmax(u, 1e-12), 1 - 1e-12)),
      slope = function(x) dnorm((x - mu) / sig) / sig
    ))
  }
  rq <- .weighted_quantile_linear(r, w, pgrid)
  .make_tail_safe_cdf_bundle(rq, pgrid)
}

# Original Gaussian-copula transport kept as a fallback.
.fit_gaussian_copula_transform <- function(X, w, ngrid = 400, tail = 1e-3,
                                           corr_shrink = 0.5) {
  X <- as.matrix(X); d <- ncol(X)
  w <- pmax(w, 0); w <- w / sum(w)
  pgrid <- seq(tail, 1 - tail, length.out = ngrid)

  # per-dim marginal bundles
  per_dim <- vector("list", d)
  for (j in seq_len(d)) {
    xj <- X[, j]; xq <- .weighted_quantile_linear(xj, w, pgrid)
    rng <- xq[length(xq)] - xq[1]
    if (!is.finite(rng) || rng < 1e-12) {
      muj <- sum(xj * w); vj <- sum(w * (xj - muj)^2); sigj <- sqrt(max(vj, 1e-6))
      per_dim[[j]] <- list(
        Fhat = function(x) pnorm((x - muj)/sigj),
        Qhat = function(u) muj + sigj*qnorm(pmin(pmax(u,1e-12),1-1e-12)),
        slope= function(x) dnorm((x - muj)/sigj)/sigj,
        F_mid=function(x) pnorm((x - muj)/sigj), Q_mid=function(u) muj + sigj*qnorm(pmin(pmax(u,1e-12),1-1e-12)),
        xmin=-Inf, xmax=Inf, vL=-8, vR=8, gL=1/sigj, gR=1/sigj
      )
    } else {
      eps <- max(1e-12, 1e-8 * max(1, abs(rng)))
      for (k in 2:length(xq)) if (xq[k] <= xq[k-1]) xq[k] <- xq[k-1] + eps
      per_dim[[j]] <- .make_tail_safe_cdf_bundle(xq, pgrid)
    }
  }

  # Θ -> Y (probit of marginals with linear tails)
  fwd_y <- function(Xnew) {
    Xnew <- as.matrix(Xnew)
    Y <- matrix(0.0, nrow(Xnew), ncol(Xnew))
    for (j in seq_len(ncol(Xnew))) {
      b <- per_dim[[j]]
      xj <- Xnew[, j]
      lt <- xj < b$xmin; rt <- xj > b$xmax; md <- !(lt | rt)
      yj <- numeric(length(xj))
      if (any(md)) yj[md] <- qnorm(b$F_mid(xj[md]))
      if (any(lt)) yj[lt] <- b$vL + b$gL * (xj[lt] - b$xmin)
      if (any(rt)) yj[rt] <- b$vR + b$gR * (xj[rt] - b$xmax)
      Y[, j] <- yj
    }
    colnames(Y) <- colnames(Xnew)
    Y
  }

  # Weighted mean/covariance in Y, shrink toward diagonal variances (low early, high late set by caller)
  Y0 <- fwd_y(X)
  muY <- colSums(Y0 * w)
  Yc  <- sweep(Y0, 2L, muY, `-`)
  S   <- t(Yc) %*% (Yc * w)  # weighted Cov(Y)
  S[!is.finite(S)] <- 0; diag(S) <- pmax(diag(S), 1e-12)
  # `diag(diag(S))` breaks in 1D because `diag()` on a length-1 numeric
  # is interpreted as an identity-matrix size request. Build the diagonal
  # matrix explicitly so both 1D and higher-dimensional cases work.
  D <- diag(as.numeric(diag(S)), nrow = ncol(S), ncol = ncol(S))
  Tmat <- (1 - corr_shrink) * D + corr_shrink * S
  Tmat <- as.matrix(Matrix::nearPD(Tmat, conv.tol = 1e-6)$mat)
  U <- tryCatch(chol(Tmat), error = function(e) chol(Tmat + diag(1e-8, ncol(S))))

  fwd <- function(Xnew) {
    Y  <- fwd_y(Xnew)
    Yc <- sweep(Y, 2L, muY, `-`)
    Z  <- t(backsolve(U, t(Yc), transpose = TRUE))
    colnames(Z) <- colnames(Xnew)
    Z
  }
  inv <- function(Znew) {
    Znew <- as.matrix(Znew)
    # reconstruct Y = muY + Z %*% U
    Y    <- Znew %*% U
    Y    <- sweep(Y, 2L, muY, `+`)
    Xr   <- matrix(0.0, nrow(Znew), ncol(Znew))
    for (j in seq_len(ncol(Znew))) {
      b <- per_dim[[j]]
      v <- Y[, j]
      v_eps <- 1e-8 + 1e-6 * max(1, abs(b$vL), abs(b$vR))
      lt <- v <= (b$vL + v_eps); rt <- v >= (b$vR - v_eps); md <- !(lt | rt)
      xj <- numeric(length(v))
      if (any(md)) xj[md] <- b$Q_mid(pnorm(v[md]))
      if (any(lt)) xj[lt] <- b$xmin + (v[lt] - b$vL)/b$gL
      if (any(rt)) xj[rt] <- b$xmax + (v[rt] - b$vR)/b$gR
      Xr[, j] <- xj
    }
    colnames(Xr) <- colnames(Znew)
    Xr
  }

  logdetU <- sum(log(diag(U)))
  log_jac <- function(Theta) {
    Theta <- as.matrix(Theta); Y <- fwd_y(Theta)
    out <- rep(0.0, nrow(Theta))
    for (j in seq_len(ncol(Theta))) {
      deriv_raw <- per_dim[[j]]$slope(Theta[, j])
      out <- out + (log(deriv_raw) - dnorm(Y[, j], log = TRUE))
    }
    out - logdetU
  }

  list(fwd = fwd, inv = inv, log_jac = log_jac,
       meta = list(C = Tmat, U = U, muY = muY, method = "gaussian_copula"))
}

# Sparse triangular transport:
# 1. Tail-safe marginal Gaussianization x -> y.
# 2. Sequential sparse autoregression y_j ~ m_j(z_{<j}) with residual Gaussianization.
fit_copula_transform <- function(X, w, ngrid = 400, tail = 1e-3,
                                 corr_shrink = 0.5,
                                 max_parents = NULL,
                                 min_parent_cor = 0.05,
                                 ridge = 1e-3,
                                 interaction_max = 3L) {
  X <- as.matrix(X)
  d <- ncol(X)
  w <- pmax(w, 0)
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) {
    w <- rep(1 / nrow(X), nrow(X))
  } else {
    w <- w / sw
  }
  pgrid <- seq(tail, 1 - tail, length.out = ngrid)
  max_parents <- as.integer(max_parents %||% min(4L, max(1L, floor(sqrt(d)))))

  # Per-dimension outer marginals x -> y.
  marginals <- vector("list", d)
  for (j in seq_len(d)) {
    xj <- X[, j]
    xq <- .weighted_quantile_linear(xj, w, pgrid)
    rng <- xq[length(xq)] - xq[1L]
    if (!is.finite(rng) || rng < 1e-12) {
      muj <- .weighted_mean1(xj, w)
      sigj <- .weighted_sd1(xj, w, mu = muj, floor = 1e-6)
      marginals[[j]] <- list(
        Fhat = function(x) pnorm((x - muj) / sigj),
        Qhat = function(u) muj + sigj * qnorm(pmin(pmax(u, 1e-12), 1 - 1e-12)),
        slope = function(x) dnorm((x - muj) / sigj) / sigj,
        F_mid = function(x) pnorm((x - muj) / sigj),
        Q_mid = function(u) muj + sigj * qnorm(pmin(pmax(u, 1e-12), 1 - 1e-12)),
        xmin = -Inf, xmax = Inf, vL = -8, vR = 8, gL = 1 / sigj, gR = 1 / sigj
      )
    } else {
      eps <- max(1e-12, 1e-8 * max(1, abs(rng)))
      for (k in 2:length(xq)) if (xq[k] <= xq[k - 1L]) xq[k] <- xq[k - 1L] + eps
      marginals[[j]] <- .make_tail_safe_cdf_bundle(xq, pgrid)
    }
  }

  fwd_y <- function(Xnew) {
    Xnew <- as.matrix(Xnew)
    Y <- matrix(0.0, nrow(Xnew), ncol(Xnew))
    for (j in seq_len(ncol(Xnew))) {
      b <- marginals[[j]]
      xj <- Xnew[, j]
      lt <- xj < b$xmin
      rt <- xj > b$xmax
      md <- !(lt | rt)
      yj <- numeric(length(xj))
      if (any(md)) yj[md] <- qnorm(b$F_mid(xj[md]))
      if (any(lt)) yj[lt] <- b$vL + b$gL * (xj[lt] - b$xmin)
      if (any(rt)) yj[rt] <- b$vR + b$gR * (xj[rt] - b$xmax)
      Y[, j] <- yj
    }
    colnames(Y) <- colnames(Xnew)
    Y
  }

  Y0 <- fwd_y(X)
  Z0 <- matrix(0.0, nrow(Y0), ncol(Y0))
  cond_models <- vector("list", d)
  resid_bundles <- vector("list", d)
  resid_sds <- numeric(d)
  for (j in seq_len(d)) {
    Zprev <- if (j > 1L) Z0[, seq_len(j - 1L), drop = FALSE] else NULL
    mean_model <- .fit_triangular_mean_model(
      y = Y0[, j],
      Zprev = Zprev,
      w = w,
      max_parents = max_parents,
      min_parent_cor = min_parent_cor,
      ridge = ridge,
      interaction_max = interaction_max
    )
    mean_fit <- .eval_triangular_mean_model(mean_model, Zprev)
    resid <- Y0[, j] - mean_fit
    resid_bundle <- .fit_residual_bundle(resid, w, pgrid)
    Z0[, j] <- qnorm(resid_bundle$Fhat(resid))
    cond_models[[j]] <- mean_model
    resid_bundles[[j]] <- resid_bundle
    resid_sds[j] <- stats::sd(resid)
  }
  colnames(Z0) <- colnames(X)

  forward_details <- function(Xnew) {
    Xnew <- as.matrix(Xnew)
    Y <- fwd_y(Xnew)
    Z <- matrix(0.0, nrow(Y), ncol(Y), dimnames = dimnames(Y))
    R <- matrix(0.0, nrow(Y), ncol(Y), dimnames = dimnames(Y))
    for (j in seq_len(ncol(Y))) {
      Zprev <- if (j > 1L) Z[, seq_len(j - 1L), drop = FALSE] else NULL
      mean_j <- .eval_triangular_mean_model(cond_models[[j]], Zprev)
      R[, j] <- Y[, j] - mean_j
      Z[, j] <- qnorm(resid_bundles[[j]]$Fhat(R[, j]))
    }
    list(Y = Y, R = R, Z = Z)
  }

  fwd <- function(Xnew) {
    forward_details(Xnew)$Z
  }

  inv <- function(Znew) {
    Znew <- as.matrix(Znew)
    Y <- matrix(0.0, nrow(Znew), ncol(Znew), dimnames = dimnames(Znew))
    Xrec <- matrix(0.0, nrow(Znew), ncol(Znew), dimnames = dimnames(Znew))
    for (j in seq_len(ncol(Znew))) {
      Zprev <- if (j > 1L) Znew[, seq_len(j - 1L), drop = FALSE] else NULL
      mean_j <- .eval_triangular_mean_model(cond_models[[j]], Zprev)
      resid_j <- resid_bundles[[j]]$Qhat(stats::pnorm(Znew[, j]))
      Y[, j] <- mean_j + resid_j
      Xrec[, j] <- marginals[[j]]$Qhat(stats::pnorm(Y[, j]))
    }
    colnames(Xrec) <- colnames(Znew)
    Xrec
  }

  log_jac <- function(Theta) {
    Theta <- as.matrix(Theta)
    det <- forward_details(Theta)
    out <- rep(0.0, nrow(Theta))
    for (j in seq_len(ncol(Theta))) {
      out <- out +
        log(marginals[[j]]$slope(Theta[, j])) - stats::dnorm(det$Y[, j], log = TRUE) +
        log(resid_bundles[[j]]$slope(det$R[, j])) - stats::dnorm(det$Z[, j], log = TRUE)
    }
    out
  }

  list(
    fwd = fwd,
    inv = inv,
    log_jac = log_jac,
    meta = list(
      method = "sparse_triangular",
      max_parents = max_parents,
      min_parent_cor = min_parent_cor,
      ridge = ridge,
      interaction_max = interaction_max,
      marginal_bundles = marginals,
      conditional_models = cond_models,
      residual_sds = resid_sds
    )
  )
}

# ---------- TRANSPORT DIAGNOSTICS (drop-in) ----------

compute_transport_diag <- function(Tmap, Theta, reference_prior, nsample = 2000L) {
  n <- nrow(Theta)
  if (n > nsample) Theta <- Theta[sample.int(n, nsample), , drop = FALSE]
  Z   <- Tmap$fwd(Theta)
  lp  <- reference_prior_logpdf(reference_prior, Theta)
  lqZ <- dmvnorm_chol_log(Z, rep(0, ncol(Z)), chol(diag(ncol(Z))))
  r   <- lp - (lqZ + Tmap$log_jac(Theta))

  S <- stats::cov(Z)
  dS <- diag(S)
  diag_dev <- max(abs(dS - 1))
  off_max  <- if (ncol(Z) > 1) max(abs(S[upper.tri(S)])) else 0

  list(ok = is.finite(mean(r)) && is.finite(stats::sd(r)),
       mean = mean(r), sd = stats::sd(r),
       diag_dev = diag_dev, off_max = off_max)
}

# ---------- SCHEDULE (drop-in) ----------
shrinkage_from_lambda <- function(lambda) {
  # low early (~0.25), high late (~0.95), smooth and monotone
  val <- 0.25 + 0.70 * (lambda^0.8)
  pmin(pmax(val, 0.20), 0.95)
}

# ---------- BUILD (drop-in) ----------
build_transport <- function(Theta, w, mu_ref = NULL, Sigma_ref = NULL, reference_prior = NULL,
                            lambda = 0,
                            ngrid = 400, tail = 1e-3,
                            verbose = TRUE) {
  reference_prior <- normalize_reference_prior(reference_prior = reference_prior, mu = mu_ref, Sigma = Sigma_ref)
  Tmap <- fit_copula_transform(
    Theta, w, ngrid = ngrid, tail = tail,
    corr_shrink = shrinkage_from_lambda(lambda)
  )
  diag <- compute_transport_diag(Tmap, Theta, reference_prior = reference_prior)
  if (isTRUE(verbose)) {
    cat(sprintf("  Transport check: mean=%.4f sd=%.4f | diag_dev=%.3f off_max=%.3f\n",
                diag$mean, diag$sd, diag$diag_dev, diag$off_max))
  }
  attr(Tmap, "diag") <- diag   # stash for relative gates
  Tmap
}

# ---------- KL-BASED REFIT (drop-in, simpler & scalar) ----------
# Weighted mean/cov and Gaussian KL to N(0, I)
.wmean_cov <- function(Z, w) {
  w <- pmax(w, 0); w <- w / sum(w)
  mu <- colSums(Z * w)
  Zc <- sweep(Z, 2L, mu, `-`)
  S  <- t(Zc) %*% (Zc * w)
  list(mu = mu, S = as.matrix(S))
}
.kl_to_standard_normal <- function(mu, S) {
  d <- length(mu)
  S <- Matrix::nearPD((S + t(S))/2, conv.tol = 1e-7)$mat
  ld <- tryCatch(determinant(S, logarithm = TRUE)$modulus[[1]], error = function(e) -Inf)
  tr <- sum(diag(S))
  0.5 * (tr - ld - d + sum(mu * mu))
}

# Refit decision: accept Tcand iff KL_cand <= (1 - kl_improve_frac) * KL_cur
# Defaults tuned for ~5–20D: modest 6% relative improvement target.
maybe_refit_transport <- function(Theta, w, lambda, round, resampled, ess_frac,
                                  Tmap, tr_state,
                                  every_rounds = 2L,
                                  kl_improve_frac = 0.06,
                                  verbose = TRUE) {
  vcat <- function(...) if (isTRUE(verbose)) base::cat(...)
  if (is.null(tr_state)) tr_state <- list(last_refit_round = 0L, last_kl = NA_real_)

  # cadence / degeneracy triggers
  time_based <- (round - tr_state$last_refit_round) >= every_rounds
  degeneracy <- isTRUE(resampled) || (ess_frac < 0.35)
  if (!(time_based || degeneracy)) return(list(refit = FALSE, Tmap = Tmap, tr_state = tr_state))

  # build candidate map with mild weight blending + scheduled shrinkage
  rho     <- max(0.30, 1 - lambda)
  w_blend <- rho * rep(1 / nrow(Theta), nrow(Theta)) + (1 - rho) * (w / sum(w))
  shrink  <- shrinkage_from_lambda(lambda)
  Tcand   <- fit_copula_transform(Theta, w_blend, corr_shrink = shrink)

  # compute KLs with current weights
  Z_cur <- Tmap$fwd(Theta)
  wc    <- .wmean_cov(Z_cur, w)
  KLc   <- .kl_to_standard_normal(wc$mu, wc$S)

  Z_can <- Tcand$fwd(Theta)
  wn    <- .wmean_cov(Z_can, w)
  KLn   <- .kl_to_standard_normal(wn$mu, wn$S)

  accept <- is.finite(KLn) && (KLn <= (1 - kl_improve_frac) * KLc)
  if (accept) {
    tr_state$last_refit_round <- round
    tr_state$last_kl <- KLn
    vcat(sprintf("  Transport refit: KL %.3f → %.3f (shrink=%.2f) [accepted]\n", KLc, KLn, shrink))
    return(list(refit = TRUE, Tmap = Tcand, tr_state = tr_state))
  } else {
    tr_state$last_kl <- KLc
    vcat(sprintf("  Transport refit: KL %.3f → %.3f (shrink=%.2f) [rejected]\n", KLc, KLn, shrink))
    return(list(refit = FALSE, Tmap = Tmap, tr_state = tr_state))
  }
}

# Identify weak dimensions from weighted variance in transformed space
# Args:
#   Z: transformed samples matrix (n x d)
#   w: sample weights (length n)
#   frac: fraction of dimensions to identify as weak
#   min_keep: minimum number of dimensions to keep as weak
# Returns: indices of dimensions with lowest weighted variance
weak_dims_from_Z <- function(Z, w, frac = 0.30, min_keep = 1L) {
  w <- pmax(w, 0); w <- w / sum(w)
  mu <- colSums(Z * w)
  v  <- colSums((sweep(Z, 2L, mu, `-`)^2) * w)
  ord <- order(v, decreasing = FALSE)
  head(ord, max(min_keep, ceiling(length(v) * frac)))
}

# Cheap deterministic ordering used as a fast alternative to Hilbert sorting.
.cheap_sort_order <- function(Z) {
  Z <- as.matrix(Z)
  n <- nrow(Z)
  if (n <= 1L) return(seq_len(n))
  if (ncol(Z) <= 1L) return(order(Z[, 1L]))
  order(rowSums(Z))
}

# ----------------------------- mixture helpers -----------------------------
# dmvt_mixture_logpdf_Z_vec, regularize_cov, merge_components,
# prune_merge_mixture_Z

# ------------- Likelihood cache helpers (tiny, string-keyed) -------------
.ll_cache_make <- function(digits = 8L, cap = 1e5, prune_every = 500L) {
  state <- new.env(parent = emptyenv())
  state$n_added <- 0L
  list(
    env = new.env(hash = TRUE, parent = emptyenv()),
    digits = as.integer(digits),
    cap = as.integer(cap),
    n_added = 0L,
    prune_every = as.integer(prune_every),
    state = state
  )
}
# faster row keys only for a small set of representatives
.ll_cache_key_row <- function(x, digits) {
  # sprintf is faster than formatC here
  paste(sprintf("%.*g", digits, x), collapse = ";")
}
# Evaluate likelihood with **batch de-duplication** and optional env cache.
# expect_dups=TRUE only when we *know* many exact duplicates exist (post-resampling).
.ll_cached_eval <- function(Theta, data, loglik_fn, cache = NULL, expect_dups = FALSE, n_cores = 1) {
  if (is.null(cache) || nrow(Theta) <= 2L) {
    return(ll_parallel(Theta, data, loglik_fn, n_cores))
  }
  # If duplicates unlikely, don't pay any key/lookup overhead.
  if (!isTRUE(expect_dups)) {
    return(ll_parallel(Theta, data, loglik_fn, n_cores))
  }
  n <- nrow(Theta)
  # Factor-based grouping is much cheaper than per-row pasting.
  df <- as.data.frame(Theta, optional = TRUE)
  grp <- do.call(interaction, c(df, drop = TRUE, lex.order = TRUE))
  gidx <- as.integer(grp)
  G <- nlevels(grp)
  if (G == n) {
    return(ll_parallel(Theta, data, loglik_fn, n_cores))
  }
  # Representative index for each group
  rep_idx <- tapply(seq_len(n), gidx, `[`, 1L)
  rep_idx <- as.integer(rep_idx)  # vector of length G
  # Build keys only for representatives
  reps <- Theta[rep_idx, , drop = FALSE]
  digits <- cache$digits
  if (is.null(digits) || !is.finite(digits)) digits <- 8L
  cap <- cache$cap
  if (is.null(cap) || !is.finite(cap)) cap <- 1e5
  prune_every <- cache$prune_every
  if (is.null(prune_every) || !is.finite(prune_every) || prune_every <= 0L) prune_every <- 500L
  ev <- cache$env
  if (is.null(ev) || !is.environment(ev)) {
    ev <- new.env(hash = TRUE, parent = emptyenv())
  }
  if (is.null(cache$state) || !is.environment(cache$state)) {
    cache$state <- new.env(parent = emptyenv())
    cache$state$n_added <- as.integer(ifelse(is.null(cache$n_added), 0L, cache$n_added))
  }
  if (!exists("n_added", envir = cache$state, inherits = FALSE) ||
      is.null(cache$state$n_added) || !is.finite(cache$state$n_added)) {
    cache$state$n_added <- 0L
  }
  keys <- apply(reps, 1L, .ll_cache_key_row, digits = digits)
  hit <- logical(G)
  rep_vals <- numeric(G)
  for (j in seq_len(G)) {
    v <- ev[[ keys[j] ]]
    if (!is.null(v)) { hit[j] <- TRUE; rep_vals[j] <- v }
  }
  if (any(!hit)) {
    miss_idx <- which(!hit)
    vals <- ll_parallel(reps[miss_idx, , drop = FALSE], data, loglik_fn, n_cores)
    rep_vals[miss_idx] <- vals
    # insert into env
    for (k in seq_along(miss_idx)) {
      ev[[ keys[miss_idx[k]] ]] <- vals[k]
    }
    cache$state$n_added <- as.integer(cache$state$n_added + length(miss_idx))
    if (cache$state$n_added >= prune_every && cap > 0L) {
      cache$state$n_added <- 0L
      # crude check: if env grew too large, just reset (cheap and safe)
      # avoids expensive ls() on every call
      keys_now <- ls(ev, all.names = TRUE)
      if (length(keys_now) > cap) {
        rm(list = keys_now, envir = ev)
      }
    }
  }
  # Map rep vals back to all rows
  rep_vals[gidx]
}

# ----------------------------- mixture sampling -----------------------------

# Combine elite mixtures from recent history into single mixture
# Args:
#   elite_mixtures: list of mixture objects from previous rounds
#   weights: relative weights for combining mixtures (most recent first)
#   housekeeping: whether to perform pruning and merging
#   min_G_keep: minimum number of mixture components to retain
#   merge_thresh: threshold for merging similar components
#   min_eig: minimum eigenvalue for covariance regularization
# Returns: combined mixture object or NULL if no valid mixtures

# Build historical mixture from elite mixture history
# Args:
#   elite_history: list of historical elite mixtures
#   lambda: current temperature parameter
#   hist_mix_lambda_thresh: threshold temperature for enabling history mixture
#   hist_mix_prob: target probability for history mixture usage
#   min_G_keep: minimum components to keep in mixture
#   merge_thresh: threshold for merging similar components
#   min_eig: minimum eigenvalue for regularization
#   force: whether to force building even below threshold
# Returns: list with mixture and probability, or NULL

# ---- DA screen mixture and calibration helpers ----
.build_da_screen_mix <- function(elite_mix, hist_mix, d) {
  mix <- elite_mix
  if (!is.null(hist_mix) && !is.null(hist_mix$mix) && !.is_empty_mix(hist_mix$mix)) {
    mix <- blend_mixes(mix, hist_mix$mix, eps = 0.10)
  }
  mix <- blend_mixes(mix, .default_std_normal_mix(d), eps = 0.05)
  mix
}

.calibrate_da_surrogate <- function(Z, Theta, lpz, lambda, ref_mix, screen_mix,
                                    data, loglik_fn, Tmap,
                                    ll_cache = NULL,
                                    n = 20L, seed = NULL,
                                    alpha = 1.0,
                                    lambda_floor = 0.20,
                                    gate_lambda = NULL,
                                    n_cores = 1) {
  # gate calibration early; return explicit "skipped" flag for logging
  gate <- if (is.null(gate_lambda)) lambda_floor else gate_lambda
  if (n <= 0 || lambda < gate) {
    return(list(a = 0.0, b = 1.0, r2 = NA_real_, skipped = TRUE, reason = sprintf("λ<%.2f", gate)))
  }
  if (!is.null(seed)) set.seed(seed)
  N <- nrow(Z); if (N == 0) return(list(a = 0.0, b = 1.0, r2 = NA_real_, skipped = TRUE, reason = "no data"))
  idx <- sample.int(N, min(N, n))
  Zs <- Z[idx, , drop = FALSE]; Thetas <- Theta[idx, , drop = FALSE]
  lref <- if (!is.null(ref_mix)) log_r_theta(Thetas, Zs, Tmap, ref_mix) else rep(0.0, nrow(Zs))
  lqmix <- gmm_logpdf_Z_vec(Zs, screen_mix$meansZ, screen_mix$cache)
  # match the stage-1 surrogate used in moves: include alpha and clamp by lambda_floor
  ltilde <- (alpha * lqmix - lpz[idx] - (1 - lambda) * lref) / max(lambda, lambda_floor, 1e-8)
  # For calibration we *don't* expect duplicates; avoid cache overhead
  ll_exact <- ll_parallel(Thetas, data, loglik_fn, n_cores = n_cores)
  # regress ll_exact ~ a + b*ltilde
  fit <- tryCatch(stats::lm(ll_exact ~ ltilde), error = function(e) NULL)
  if (is.null(fit)) return(list(a = 0.0, b = 1.0, r2 = NA_real_, skipped = TRUE, reason = "lm failed"))
  cf <- stats::coef(fit); a <- as.numeric(cf[1]); b <- as.numeric(cf[2])
  if (!is.finite(a)) a <- 0.0
  if (!is.finite(b)) b <- 1.0
  b <- .clamp(b, 0.5, 1.5)
  r2 <- tryCatch(summary(fit)$r.squared, error = function(e) NA_real_)
  list(a = a, b = b, r2 = r2, skipped = FALSE)
}

# ------------- Difficulty index and light-touch adaptations -------------

# Fit weighted Gaussian mixture model in whitened coordinates
# Uses robust EM algorithm with Cholesky regularization and warm start support
# Args:
#   W: whitened data matrix (n x d)
#   we: sample weights (length n)
#   G: number of mixture components
#   init: optional initialization (means, covs, weights)
#   itmax: maximum EM iterations
#   tol: convergence tolerance
# Returns: list with fitted means, covariances, and weights



# Fit mixture model on elite particles in transformed space
# Captures high-weight modes to drive independence proposals
# Args:
#   Z: transformed sample matrix (n x d)
#   w: sample weights (length n)
#   elite_quantile: quantile threshold for elite particles
#   G: target number of mixture components
#   cov_inflation: inflation factor for covariances
#   min_elite: minimum number of elite particles
#   warm_start_mixZ: optional warm start mixture
#   em_itmax: maximum EM iterations
#   housekeeping: whether to perform pruning and merging
#   min_G_keep: minimum components to keep
#   merge_thresh: threshold for merging components
#   verbose: whether to print verbose output
#   lambda: current temperature for eigenvalue regularization
# Returns: fitted mixture object




# ----------------------------- resampling -----------------------------

# Vectorized MCMC moves in transformed space with mixture proposals
# Combines random-walk and independence proposals from elite/history mixtures
# Args:
#   Z: transformed samples matrix (n x d)
#   loglik: log-likelihood values
#   lpz: log-density in transformed space
#   Tmap: transport map object
#   lambda: current temperature parameter
#   mu_ref: prior mean
#   prior_L: prior Cholesky factor
#   w: particle weights
#   elite_mix: elite mixture for independence proposals
#   hist_mix: historical mixture for independence proposals
#   data: data for likelihood computation
#   loglik_fn: log-likelihood function
#   n_moves: number of MCMC moves per particle
#   rw_prob: probability of random-walk moves
#   rw_scale: scaling for random-walk proposals
#   indep_t_df: degrees of freedom for t-distribution proposals
#   indep_t_prob: probability of using t-distribution vs Gaussian
#   pcn_prob: probability of pCN (preconditioned Crank-Nicolson) moves
#   pcn_beta: pCN parameter (adapted by caller)
#   seed: random seed
#   param_names: parameter names
#   weak_dim_idx: indices of weak dimensions for inflation
#   rw_expand_factor: expansion factor for weak dimensions
#   resampled: flag indicating if called right after resampling
# Returns: list with updated Z, loglik, lpz, and acceptance rates
mcmc_moves_z_mix_batched <- function(Z, loglik, lpz, Tmap, lambda,
                                     mu_ref, prior_L, w,
                                     elite_mix, hist_mix = NULL,
                                     reference_prior,
                                     data, loglik_fn, n_moves = 1,
                                     rw_prob = 0.35, rw_scale = 0.9,
                                     indep_t_df = 7, indep_t_prob = 0.75,
                                     pcn_prob = 0.25, pcn_beta = 0.4,
                                     seed = NULL,
                                     param_names = NULL,
                                     weak_dim_idx = integer(0),
                                     rw_expand_factor = 2.5,
                                     resampled = FALSE,
                                     ref_mix = NULL,
                                     lp_ref_vec = NULL,
                                     # --- Delayed Acceptance (DA) config (snapshot per round) ---
                                     da_enable = FALSE,
                                     da_lambda_floor = 0.20,   # <— NEW: clamp denominator for tiny λ
                                     da_alpha = 1.00,         # mixture power α (tuned per round)
                                     da_calib = c(0.0, 1.0),  # (a, b) affine calibration of ltilde
                                     da_screen_mix = NULL,    # mixture in Z used to build surrogate
                                     # --- Likelihood cache (optional) ---
                                     ll_cache = NULL,
                                     base_logpdf_fn = NULL,
                                     allow_pcn = TRUE,
                                     n_cores = 1
                                     ) {
  N <- nrow(Z); d <- ncol(Z)
  if (!is.null(seed)) set.seed(seed)
  if (is.null(param_names)) param_names <- colnames(Z)
  gss_on <- !is.null(ref_mix)
  allow_pcn <- isTRUE(allow_pcn)
  if (!allow_pcn) pcn_prob <- 0
  lpz_from_theta <- function(Theta_mat) {
    base_log <- if (is.null(base_logpdf_fn)) {
      reference_prior_logpdf(reference_prior, Theta_mat)
    } else {
      base_logpdf_fn(Theta_mat)
    }
    as.numeric(base_log - Tmap$log_jac(Theta_mat))
  }
  standard_normal_logpdf_rows <- function(Zmat) {
    rowSums(stats::dnorm(as.matrix(Zmat), log = TRUE))
  }

  # Prior-blended RW metric in Z-space
  S_emp <- weighted_cov(Z, w); if (any(!is.finite(S_emp))) S_emp <- diag(d)
  eta <- max(0.30, 1 - lambda)
  S_prop <- (1 - eta) * S_emp + eta * diag(d)
  ev <- eigen(S_prop, symmetric = TRUE)
  # Boost eigenvalue floor right after resampling to avoid collapse
  floor0 <- 1e-2 * (1 - lambda) + 1e-3
  if (isTRUE(resampled)) floor0 <- floor0 + 2e-2
  lamv <- pmax(ev$values, floor0)
  if (length(weak_dim_idx)) {
    lamv[weak_dim_idx] <- lamv[weak_dim_idx] * rw_expand_factor
  }
  S_prop <- ev$vectors %*% diag(lamv, d) %*% t(ev$vectors)
  Lrw <- tryCatch(chol(S_prop + diag(1e-8, d)), error = function(e) diag(d))
  step_scale <- rw_scale / sqrt(max(d, 1))

  acc_rw <- 0L; acc_id <- 0L; acc_pcn <- 0L
  prop_rw <- 0L; prop_id <- 0L; prop_pcn <- 0L     # proposals per kernel
  # DA stats (overall + per-kernel)
  da_prop <- 0L; da_pass <- 0L
  da_prop_rw <- 0L; da_pass_rw <- 0L
  da_prop_pcn <- 0L; da_pass_pcn <- 0L
  da_prop_id <- 0L; da_pass_id <- 0L

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

  # ---- DA surrogate helpers (snapshot for this call) ----
  # Turn DA on from λ=0 if enabled and we have a screen mixture.
  da_on <- isTRUE(da_enable) && !is.null(da_screen_mix)
  da_eps <- 1e-8
  # Keep one helper for exact pCN target to avoid DA / non-DA drift.
  pcn_exact_target <- function(lpz_vec, Zmat, loglik_vec, lref_vec = NULL) {
    base <- lpz_vec - standard_normal_logpdf_rows(Zmat)
    if (gss_on) {
      if (is.null(lref_vec)) stop("pcn_exact_target requires lref_vec when GSS is enabled.")
      base + lambda * loglik_vec + (1 - lambda) * lref_vec
    } else {
      base + lambda * loglik_vec
    }
  }
  # vectorized Ltilde for a batch; lpz_vec must correspond to thetas supplied
  Ltilde_batch <- function(Zbat, Theta_bat, lpz_vec, lref_vec = NULL) {
    # log q_mix in Z
    lqmix <- gmm_logpdf_Z_vec(Zbat, da_screen_mix$meansZ, da_screen_mix$cache)
    if (is.null(lref_vec)) lref_vec <- rep(0.0, nrow(Zbat))
    # ltilde_mix as calibrated pseudo-likelihood
    # IMPORTANT: clamp denominator to λ_floor so stage-1 remains well-behaved for λ → 0
    ltilde <- (da_alpha * lqmix - lpz_vec - (1 - lambda) * lref_vec) / max(lambda, da_lambda_floor, da_eps)
    ltilde <- da_calib[1] + da_calib[2] * ltilde  # affine calibration
    # Stage-1 surrogate total
    lpz_vec + lambda * ltilde + (1 - lambda) * lref_vec
  }

  sample_from_q <- function(M) {
    use_hist <- !is.null(hist_mix)
    hprob <- if (use_hist) hist_mix$prob else 0.0
    take_hist <- runif(M) < hprob
    t_flags   <- runif(M) < indep_t_prob
    Zout <- matrix(NA_real_, M, d)
    # ensure non-empty mixtures
    em <- elite_mix; hm <- if (use_hist) hist_mix$mix else NULL
    if (.is_empty_mix(em)) em <- .default_std_normal_mix(d)
    if (use_hist && .is_empty_mix(hm)) hm <- .default_std_normal_mix(d)
    n_hi_t <- sum(take_hist & t_flags)
    n_hi_n <- sum(take_hist & !t_flags)
    n_el_t <- sum(!take_hist & t_flags)
    n_el_n <- sum(!take_hist & !t_flags)
    if (n_hi_n) { Zout[take_hist & !t_flags, ] <- sample_gmm_Z_qmc(n_hi_n, hm$meansZ, hm$cache) }
    if (n_hi_t) { Zout[take_hist &  t_flags, ] <- rmvt_mixture_Z_qmc(n_hi_t, hm$meansZ, hm$cache, nu = indep_t_df) }
    if (n_el_n) { Zout[!take_hist & !t_flags, ] <- sample_gmm_Z_qmc(n_el_n, em$meansZ, em$cache) }
    if (n_el_t) { Zout[!take_hist &  t_flags, ] <- rmvt_mixture_Z_qmc(n_el_t, em$meansZ, em$cache, nu = indep_t_df) }
    colnames(Zout) <- colnames(Z)
    list(Z = Zout, hprob = hprob)
  }

  do_pcn <- function(idx) {
    if (!length(idx)) return(invisible(NULL))
    k <- length(idx); prop_pcn <<- prop_pcn + k   # Track proposal count
    Zc <- Z[idx, , drop = FALSE]
    Theta_c <- Tmap$inv(Zc); colnames(Theta_c) <- param_names
    beta <- min(max(pcn_beta, 1e-6), 1.0)
    Xi <- matrix(rnorm(k * d), nrow = k, ncol = d)
    Zp <- sqrt(1 - beta^2) * Zc + beta * Xi
    Theta_p <- Tmap$inv(Zp)
    colnames(Theta_p) <- param_names
    # ---- DA: Stage-1 (cheap) ----
    if (da_on) {
      lref_c <- if (gss_on) { if (!is.null(lp_ref_vec)) lp_ref_vec[idx] else log_r_theta(Theta_c, Zc, Tmap, ref_mix) } else NULL
      lref_p <- if (gss_on) log_r_theta(Theta_p, Zp, Tmap, ref_mix) else NULL
      lpz_p  <- lpz_from_theta(Theta_p)
      Lc_t   <- Ltilde_batch(Zc, Theta_c, lpz[idx], lref_c) - standard_normal_logpdf_rows(Zc)
      Lp_t   <- Ltilde_batch(Zp, Theta_p, lpz_p, lref_p) - standard_normal_logpdf_rows(Zp)
      a1 <- Lp_t - Lc_t                     # symmetric proposal ⇒ no q terms
      u1 <- log(runif(k))
      pass <- which(u1 < pmin(0, a1))
      da_prop <<- da_prop + k; da_pass <<- da_pass + length(pass)
      da_prop_pcn <<- da_prop_pcn + k; da_pass_pcn <<- da_pass_pcn + length(pass)
      if (!length(pass)) return(invisible(NULL))
      # ---- DA: Stage-2 (exact correction on pass) ----
      ll_p_pass <- .ll_cached_eval(Theta_p[pass, , drop = FALSE], data, loglik_fn, ll_cache,
                                   expect_dups = resampled, n_cores = n_cores)
      Lc_ex <- pcn_exact_target(
        lpz[idx[pass]],
        Zc[pass, , drop = FALSE],
        loglik[idx[pass]],
        if (gss_on) lref_c[pass] else NULL
      )
      Lp_ex <- pcn_exact_target(
        lpz_p[pass],
        Zp[pass, , drop = FALSE],
        ll_p_pass,
        if (gss_on) lref_p[pass] else NULL
      )
      corr  <- (Lp_ex - Lc_ex) - (Lp_t[pass] - Lc_t[pass])
      u2 <- log(runif(length(pass)))
      acc_idx <- pass[which(u2 < pmin(0, corr))]
      if (length(acc_idx)) {
        Z[idx[acc_idx], ]    <<- Zp[acc_idx, , drop = FALSE]
        loglik[idx[acc_idx]] <<- ll_p_pass[match(acc_idx, pass)]
        lpz[idx[acc_idx]]    <<- lpz_p[acc_idx]
        acc_pcn              <<- acc_pcn + length(acc_idx)
      }
    } else {
      # Original exact path
      ll_p <- .ll_cached_eval(Theta_p, data, loglik_fn, ll_cache, expect_dups = resampled, n_cores = n_cores)
      if (gss_on) {
        if (!is.null(lp_ref_vec)) {
          lref_c <- lp_ref_vec[idx]
        } else {
          lref_c <- log_r_theta(Theta_c, Zc, Tmap, ref_mix)
        }
        lref_p <- log_r_theta(Theta_p, Zp, Tmap, ref_mix)
        lpz_p <- lpz_from_theta(Theta_p)
        a <- pcn_exact_target(lpz_p, Zp, ll_p, lref_p) -
          pcn_exact_target(lpz[idx], Zc, loglik[idx], lref_c)
      } else {
        lpz_p <- lpz_from_theta(Theta_p)
        a <- pcn_exact_target(lpz_p, Zp, ll_p) -
          pcn_exact_target(lpz[idx], Zc, loglik[idx])
      }
      u <- log(runif(k)); ia <- which(u < pmin(0, a))
      if (length(ia)) {
        Z[idx[ia], ]       <<- Zp[ia, , drop = FALSE]
        loglik[idx[ia]]    <<- ll_p[ia]
        lpz[idx[ia]]       <<- lpz_p[ia]
        acc_pcn            <<- acc_pcn + length(ia)
      }
    }
    invisible(NULL)
  }

  do_rw <- function(idx) {
    if (!length(idx)) return(invisible(NULL))
    k <- length(idx); prop_rw <<- prop_rw + k
    Zc <- Z[idx, , drop = FALSE]
    eps <- matrix(rnorm(k * d), k, d) %*% t(Lrw)
    Zp <- Zc + eps * step_scale; colnames(Zp) <- param_names
    Theta_p <- Tmap$inv(Zp); colnames(Theta_p) <- param_names
    if (da_on) {
      Theta_c <- Tmap$inv(Zc); colnames(Theta_c) <- param_names
      lref_c <- if (gss_on) { if (!is.null(lp_ref_vec)) lp_ref_vec[idx] else log_r_theta(Theta_c, Zc, Tmap, ref_mix) } else NULL
      lref_p <- if (gss_on) log_r_theta(Theta_p, Zp, Tmap, ref_mix) else NULL
      lpz_p  <- lpz_from_theta(Theta_p)
      Lc_t <- Ltilde_batch(Zc, Theta_c, lpz[idx], lref_c)
      Lp_t <- Ltilde_batch(Zp, Theta_p, lpz_p, lref_p)
      a1 <- Lp_t - Lc_t
      u1 <- log(runif(k)); pass <- which(u1 < pmin(0, a1))
      da_prop <<- da_prop + k; da_pass <<- da_pass + length(pass)
      da_prop_rw <<- da_prop_rw + k; da_pass_rw <<- da_pass_rw + length(pass)
      if (!length(pass)) return(invisible(NULL))
      # exact on pass
      ll_p <- .ll_cached_eval(Theta_p[pass, , drop = FALSE], data, loglik_fn, ll_cache,
                               expect_dups = resampled, n_cores = n_cores)
      Lc_ex <- lpz[idx[pass]] + lambda*loglik[idx[pass]] + if (gss_on) (1-lambda)*lref_c[pass] else 0
      Lp_ex <- lpz_p[pass]     + lambda*ll_p             + if (gss_on) (1-lambda)*lref_p[pass] else 0
      corr  <- (Lp_ex - Lc_ex) - (Lp_t[pass] - Lc_t[pass])
      u2 <- log(runif(length(pass)))
      acc_idx <- pass[which(u2 < pmin(0, corr))]
      if (length(acc_idx)) {
        Z[idx[acc_idx], ]          <<- Zp[acc_idx, , drop = FALSE]
        loglik[idx[acc_idx]]       <<- ll_p[match(acc_idx, pass)]
        lpz[idx[acc_idx]]          <<- lpz_p[acc_idx]
        acc_rw                     <<- acc_rw + length(acc_idx)
      }
    } else {
      ll_p <- .ll_cached_eval(Theta_p, data, loglik_fn, ll_cache, expect_dups = resampled, n_cores = n_cores)
      lpz_p <- lpz_from_theta(Theta_p)
      if (gss_on) {
        if (!is.null(lp_ref_vec)) {
          lref_c <- lp_ref_vec[idx]
        } else {
          Theta_c <- Tmap$inv(Zc); colnames(Theta_c) <- param_names
          lref_c <- log_r_theta(Theta_c, Zc, Tmap, ref_mix)
        }
        lref_p <- log_r_theta(Theta_p, Zp, Tmap, ref_mix)
        lt_p <- lpz_p + lambda * ll_p + (1 - lambda) * lref_p
        lt_c <- lpz[idx] + lambda * loglik[idx] + (1 - lambda) * lref_c
      } else {
        lt_p <- lpz_p + lambda * ll_p; lt_c <- lpz[idx] + lambda * loglik[idx]
      }
      a <- lt_p - lt_c; u <- log(runif(k)); acc_idx <- which(u < pmin(0, a))
      if (length(acc_idx)) {
        Z[idx[acc_idx], ] <<- Zp[acc_idx, , drop = FALSE]
        loglik[idx[acc_idx]] <<- ll_p[acc_idx]
        lpz[idx[acc_idx]] <<- lpz_p[acc_idx]
        acc_rw <<- acc_rw + length(acc_idx)
      }
    }
    invisible(NULL)
  }

  do_indep <- function(idx) {
    if (!length(idx)) return(invisible(NULL))
    k <- length(idx); prop_id <<- prop_id + k
    Zc <- Z[idx, , drop = FALSE]
    # Ensure Theta_c is available when DA is on
    Theta_c <- if (da_on) { tmp <- Tmap$inv(Zc); colnames(tmp) <- param_names; tmp } else NULL
    lq_c <- log_q_mixture(Zc, elite_mix, hist_mix, indep_t_prob,
                          if (!is.null(hist_mix)) hist_mix$prob else 0.0)
    if (gss_on) {
      if (!is.null(lp_ref_vec)) {
        lref_c <- lp_ref_vec[idx]
      } else {
        if (is.null(Theta_c)) { Theta_c <- Tmap$inv(Zc); colnames(Theta_c) <- param_names }
        lref_c <- log_r_theta(Theta_c, Zc, Tmap, ref_mix)
      }
      lt_c_exact <- lpz[idx] + lambda * loglik[idx] + (1 - lambda) * lref_c
    } else {
      lt_c_exact <- lpz[idx] + lambda * loglik[idx]
      lref_c <- NULL
    }

    smp <- sample_from_q(k)
    Zp <- smp$Z
    lq_p <- log_q_mixture(Zp, elite_mix, hist_mix, indep_t_prob, smp$hprob)
    Theta_p <- Tmap$inv(Zp); colnames(Theta_p) <- param_names
    if (da_on) {
      lpz_p <- lpz_from_theta(Theta_p)
      lref_p <- if (gss_on) log_r_theta(Theta_p, Zp, Tmap, ref_mix) else NULL
      Lc_t <- Ltilde_batch(Zc, Theta_c, lpz[idx], lref_c)
      Lp_t <- Ltilde_batch(Zp, Theta_p, lpz_p, lref_p)
      a1 <- (Lp_t - lq_p) - (Lc_t - lq_c)
      u1 <- log(runif(k)); pass <- which(u1 < pmin(0, a1))
      da_prop <<- da_prop + k; da_pass <<- da_pass + length(pass)
      da_prop_id <<- da_prop_id + k; da_pass_id <<- da_pass_id + length(pass)
      if (!length(pass)) return(invisible(NULL))
      # exact only for pass
      ll_p <- .ll_cached_eval(Theta_p[pass, , drop = FALSE], data, loglik_fn, ll_cache,
                               expect_dups = resampled, n_cores = n_cores)
      Lp_ex <- lpz_p[pass] + lambda*ll_p + if (gss_on) (1-lambda)*lref_p[pass] else 0
      a2 <- (Lp_ex - lq_p[pass]) - (lt_c_exact[pass] - lq_c[pass]) - (a1[pass])
      u2 <- log(runif(length(pass)))
      acc_idx <- pass[which(u2 < pmin(0, a2))]
      if (length(acc_idx)) {
        Z[idx[acc_idx], ]       <<- Zp[acc_idx, , drop = FALSE]
        loglik[idx[acc_idx]]    <<- ll_p[match(acc_idx, pass)]
        lpz[idx[acc_idx]]       <<- lpz_p[acc_idx]
        acc_id                  <<- acc_id + length(acc_idx)
      }
    } else {
      ll_p <- .ll_cached_eval(Theta_p, data, loglik_fn, ll_cache, expect_dups = resampled, n_cores = n_cores)
      lpz_p <- lpz_from_theta(Theta_p)
      if (gss_on) {
        lref_p <- log_r_theta(Theta_p, Zp, Tmap, ref_mix)
        lt_p <- lpz_p + lambda * ll_p + (1 - lambda) * lref_p
      } else {
        lt_p <- lpz_p + lambda * ll_p
      }
      a <- (lt_p - lq_p) - (lt_c_exact - lq_c)
      u <- log(runif(k))
      acc_idx <- which(u < pmin(0, a))
      if (length(acc_idx)) {
        Z[idx[acc_idx], ] <<- Zp[acc_idx, , drop = FALSE]
        loglik[idx[acc_idx]] <<- ll_p[acc_idx]
        lpz[idx[acc_idx]] <<- lpz_p[acc_idx]
        acc_id <<- acc_id + length(acc_idx)
      }
    }
    invisible(NULL)
  }

  # No prior-axis refresh path in simplified version

  for (m in 1:n_moves) {
    u <- runif(N)
    th1 <- pcn_prob
    th2 <- pcn_prob + rw_prob
    idx_pcn   <- which(u < th1)
    idx_rw    <- which(u >= th1 & u < th2)
    idx_id    <- which(u >= th2)
    if (length(idx_pcn))   do_pcn(idx_pcn)
    if (length(idx_rw))    do_rw(idx_rw)
    if (length(idx_id))    do_indep(idx_id)
  }

  acc_total <- (acc_rw + acc_id + acc_pcn) / (N * n_moves)

  list(
    Z = Z, loglik = loglik, lpz = lpz,
    accept_rate = acc_total,
    rw_accept_rate   = if (prop_rw>0) acc_rw/prop_rw else 0,
    indep_accept_rate= if (prop_id>0) acc_id/prop_id else 0,
    pcn_accept_rate  = if (prop_pcn>0) acc_pcn/prop_pcn else 0,
    # raw proposal counts (for diagnostics)
    prop_rw = prop_rw, prop_pcn = prop_pcn, prop_id = prop_id,
    # DA diagnostics
    da_prop = da_prop, da_pass = da_pass,
    da_prop_rw = da_prop_rw, da_pass_rw = da_pass_rw,
    da_prop_pcn = da_prop_pcn, da_pass_pcn = da_pass_pcn,
    da_prop_id = da_prop_id, da_pass_id = da_pass_id,
    da_pass_rate = if (da_prop>0) da_pass/da_prop else NA_real_
  )
}


## ----- GSS helpers (on-the-fly reference mix) -----
# Compute log r(theta) when 'ref_mix' is a mixture in Z-space
# r(theta) = q_Z(Z) * |dZ/dtheta|
log_r_theta <- function(Theta, Z, Tmap, ref_mix) {
  if (is.null(ref_mix) || .is_empty_mix(ref_mix)) {
    return(rep(-Inf, nrow(Theta)))  # indicates "no ref" (should not be used in GSS mode)
  }
  lqZ <- gmm_logpdf_Z_vec(Z, ref_mix$meansZ, ref_mix$cache)
  lqZ + Tmap$log_jac(Theta)
}

# Decide when to snapshot/refresh the GSS reference mix
# Returns list(new_ref_mix, switched: TRUE/FALSE)
maybe_update_ref_mix <- function(lambda, elite_mix, hist_mix,
                                 ref_mix, allow_refresh = TRUE,
                                 lambda_ref_snap = 0.30,
                                 lambda_ref_refresh = 0.60,
                                 d = NULL,
                                 snapshot_taken = FALSE,
                                 refresh_taken = FALSE,
                                 single_update_policy = TRUE) {
  cand <- elite_mix
  # If you maintain a history mixture, softly blend it to stabilise the ref
  if (!is.null(hist_mix) && !is.null(hist_mix$mix) && !.is_empty_mix(hist_mix$mix)) {
    # Simple convex blend; small epsilon to ensure broad tails
    cand <- blend_mixes(cand, hist_mix$mix, eps = 0.10)
  }
  # Always add a tiny base mass for tail protection
  if (is.null(d)) {
    # Try to infer d if possible; otherwise skip base blend only if we truly cannot infer
    if (!.is_empty_mix(cand)) d <- length(cand$meansZ[[1]])
  }
  if (!is.null(d)) {
    cand <- blend_mixes(cand, .default_std_normal_mix(d), eps = 0.05)
  }

  if (!single_update_policy) {
    # Legacy behavior
    if (is.null(ref_mix) && lambda >= lambda_ref_snap && !.is_empty_mix(cand)) {
      return(list(ref_mix = cand, switched = TRUE, event = "snapshot",
                  snapshot_taken = snapshot_taken, refresh_taken = refresh_taken))
    }
    if (!is.null(ref_mix) && allow_refresh && lambda >= lambda_ref_refresh && !isTRUE(attr(ref_mix, "refreshed"))) {
      attr(cand, "refreshed") <- TRUE
      return(list(ref_mix = cand, switched = TRUE, event = "refresh",
                  snapshot_taken = snapshot_taken, refresh_taken = refresh_taken))
    }
    return(list(ref_mix = ref_mix, switched = FALSE, event = "none",
                snapshot_taken = snapshot_taken, refresh_taken = refresh_taken))
  }

  # Strict-stability policy: at most one snapshot and one refresh attempt overall.
  if (!snapshot_taken && is.null(ref_mix) && lambda >= lambda_ref_snap && !.is_empty_mix(cand)) {
    return(list(ref_mix = cand, switched = TRUE, event = "snapshot",
                snapshot_taken = TRUE, refresh_taken = refresh_taken))
  }
  if (snapshot_taken && !refresh_taken && !is.null(ref_mix) &&
      allow_refresh && lambda >= lambda_ref_refresh && !.is_empty_mix(cand)) {
    return(list(ref_mix = cand, switched = TRUE, event = "refresh",
                snapshot_taken = snapshot_taken, refresh_taken = TRUE))
  }
  list(ref_mix = ref_mix, switched = FALSE, event = "none",
       snapshot_taken = snapshot_taken, refresh_taken = refresh_taken)
}

# ----------------------------- convergence helpers -----------------------------

# Removed unused helper: grow_particles_keep (was never called)

.normalize_checkpoint_lambdas <- function(checkpoint_lambdas, lambda_target) {
  cp <- sort(unique(as.numeric(checkpoint_lambdas %||% numeric(0))))
  cp <- cp[is.finite(cp) & cp > 0 & cp <= as.numeric(lambda_target)]
  if (!length(cp)) {
    return(numeric(0))
  }
  if (abs(tail(cp, 1L) - as.numeric(lambda_target)) > 1e-12) {
    cp <- c(cp, as.numeric(lambda_target))
  }
  cp
}

.seed_plan_block <- function(seed_plan, block, round = NULL, default = NULL) {
  if (is.null(seed_plan)) {
    return(as.integer(default))
  }
  if (is.null(round)) {
    val <- seed_plan[[block]] %||% default
  } else {
    round_tbl <- seed_plan$round %||% NULL
    if (is.null(round_tbl) || NROW(round_tbl) < as.integer(round) || !(block %in% colnames(round_tbl))) {
      val <- default
    } else {
      val <- round_tbl[as.integer(round), block]
    }
  }
  val <- as.integer(val %||% default)
  if (!is.finite(val)) as.integer(default) else val
}

run_tempered_smc <- function(reference_prior,
                             bridge_stat_fn,
                             base_logpdf_fn = NULL,
                             initial_particles = NULL,
                             initial_weights = NULL,
                             initial_log_normalizer = 0,
                             M = NULL,
                             resample_threshold = 0.6,
                             n_mcmc_moves = 3L,
                             max_rounds = 200L,
                             lambda_target = 1.0,
                             checkpoint_lambdas = NULL,
                             cess_target = NULL,
                             G_mix = 16L,
                             gamma_sharp = 0.7,
                             refit_every = 2L,
                             rw_prob = 0.6,
                             rw_scale_init = 0.9,
                             indep_t_df = 4L,
                             indep_t_prob = 0.85,
                             hist_mix_enable = TRUE,
                             hist_mix_lambda_thresh = 0.50,
                             hist_mix_prob = 0.20,
                             gss_enable = FALSE,
                             lambda_ref_snap = 0.15,
                             allow_ref_refresh = TRUE,
                             lambda_ref_refresh = 0.50,
                             da_enable = TRUE,
                             da_lambda_floor = 0.20,
                             da_target_pass = 0.45,
                             da_calibrate_n = 20L,
                             da_alpha_init = 1.00,
                             da_rm_gain = 0.20,
                             adapt_lambda_max = 0.15,
                             freeze_transport_after_gss = TRUE,
                             single_gss_updates = TRUE,
                             post_adapt_n_mcmc_moves = 3L,
                             pre_resample_lambda_gate = 0,
                             deterministic_resampling = FALSE,
                             resample_sort_mode = c("adaptive", "hilbert", "cheap1d", "none"),
                             hilbert_hard_ess = 0.25,
                             ll_cache_enable = TRUE,
                             ll_cache_digits = 8L,
                             ll_cache_cap = 100000L,
                             seed_plan = NULL,
                             warm_start_fit = NULL,
                             warm_start_use_transport = TRUE,
                             warm_start_use_mixture = TRUE,
                             warm_start_max_particles = 2000L,
                             n_cores = 1,
                             seed = 123,
                             verbose = TRUE) {
  vcat <- function(...) { if (verbose) base::cat(...) }
  cat <- vcat
  set.seed(.seed_plan_block(seed_plan, "global", default = seed))
  resample_sort_mode <- match.arg(resample_sort_mode)
  post_adapt_n_mcmc_moves <- as.integer(max(1L, post_adapt_n_mcmc_moves))
  lambda_target <- as.numeric(lambda_target)
  if (!is.finite(lambda_target) || lambda_target <= 0 || lambda_target > 1) {
    stop("lambda_target must be in (0, 1].")
  }
  checkpoint_lambdas <- .normalize_checkpoint_lambdas(checkpoint_lambdas, lambda_target = lambda_target)
  checkpoint_log_evidence <- numeric(length(checkpoint_lambdas))
  checkpoint_log_increments <- numeric(length(checkpoint_lambdas))
  checkpoint_ptr <- 1L
  if (!is.function(bridge_stat_fn)) {
    stop("bridge_stat_fn must be a function of the particle matrix.")
  }
  reference_prior <- normalize_reference_prior(reference_prior = reference_prior)
  ref_geom <- reference_prior_geometry(reference_prior)
  mu_ref <- ref_geom$mean
  param_names <- ref_geom$param_names
  prior_L <- ref_geom$chol
  allow_pcn <- TRUE
  if (isTRUE(verbose)) {
    cat("  Using transported-space pCN moves.\n")
  }
  base_logpdf <- if (is.null(base_logpdf_fn)) {
    function(Theta_mat) as.numeric(reference_prior_logpdf(reference_prior, Theta_mat))
  } else {
    function(Theta_mat) as.numeric(base_logpdf_fn(Theta_mat))
  }
  bridge_stat_kernel <- function(Theta_mat, data = NULL) {
    as.numeric(bridge_stat_fn(Theta_mat))
  }
  # Build likelihood cache (persists across rounds; valid in θ-space)
  ll_cache <- if (ll_cache_enable) .ll_cache_make(digits = ll_cache_digits, cap = ll_cache_cap) else NULL
  cat("Stage 1: initialize particles & build transport...\n")
  set.seed(.seed_plan_block(seed_plan, "init", default = seed))
  if (is.null(initial_particles)) {
    M <- as.integer(M %||% 5000L)
    Theta <- reference_prior_sample(reference_prior, M)
    w <- rep(1 / M, M)
  } else {
    Theta <- as.matrix(initial_particles)
    if (ncol(Theta) != length(param_names)) {
      stop("initial_particles dimension does not match reference_prior.")
    }
    colnames(Theta) <- param_names
    w_in <- pmax(as.numeric(initial_weights %||% rep(1 / nrow(Theta), nrow(Theta))), 0)
    if (length(w_in) != nrow(Theta)) {
      stop("initial_weights must match nrow(initial_particles).")
    }
    sw_in <- sum(w_in)
    if (!is.finite(sw_in) || sw_in <= 0) {
      stop("initial_weights must sum to a positive finite value.")
    }
    w_in <- w_in / sw_in
    M <- as.integer(M %||% nrow(Theta))
    if (M <= 0L) stop("M must be positive.")
    if (nrow(Theta) != M) {
      cw <- c(0, cumsum(w_in))
      cw[length(cw)] <- 1
      u0 <- if (isTRUE(deterministic_resampling)) 0.5 / M else stats::runif(1) / M
      idx_init <- findInterval(u0 + (0:(M - 1L)) / M, cw, rightmost.closed = TRUE)
      Theta <- Theta[idx_init, , drop = FALSE]
      w <- rep(1 / M, M)
    } else {
      w <- w_in
    }
  }
  colnames(Theta) <- param_names
  loglik <- bridge_stat_kernel(Theta)
  warm_theta <- NULL
  warm_w <- NULL
  warm_mix_seed <- NULL
  if (!is.null(warm_start_fit) &&
      !is.null(warm_start_fit$Theta) &&
      ncol(as.matrix(warm_start_fit$Theta)) == ncol(Theta)) {
    warm_theta <- as.matrix(warm_start_fit$Theta)
    if (nrow(warm_theta) > as.integer(warm_start_max_particles)) {
      warm_idx <- seq_len(nrow(warm_theta))
      warm_w_in <- pmax(as.numeric(warm_start_fit$w %||% rep(1 / nrow(warm_theta), nrow(warm_theta))), 0)
      sw_warm <- sum(warm_w_in)
      if (!is.finite(sw_warm) || sw_warm <= 0) {
        warm_w_in <- rep(1 / nrow(warm_theta), nrow(warm_theta))
      } else {
        warm_w_in <- warm_w_in / sw_warm
      }
      set.seed(.seed_plan_block(seed_plan, "warm_subsample", default = seed + 1L))
      warm_idx <- sample.int(
        nrow(warm_theta),
        size = as.integer(warm_start_max_particles),
        replace = FALSE,
        prob = warm_w_in
      )
      warm_theta <- warm_theta[warm_idx, , drop = FALSE]
      warm_w <- warm_w_in[warm_idx]
    } else {
      warm_w <- pmax(as.numeric(warm_start_fit$w %||% rep(1 / nrow(warm_theta), nrow(warm_theta))), 0)
    }
    sw_warm <- sum(warm_w)
    if (!is.finite(sw_warm) || sw_warm <= 0) {
      warm_w <- rep(1 / nrow(warm_theta), nrow(warm_theta))
    } else {
      warm_w <- warm_w / sw_warm
    }
    warm_mix_seed <- warm_start_fit$elite_mix_final %||% NULL
  }
  # Initial transport with safe defaults + diagnostics. If a nearby local fit is
  # available, reuse its particle cloud only to initialize proposal geometry.
  Tmap <- NULL
  if (isTRUE(warm_start_use_transport) && !is.null(warm_theta)) {
    warm_lambda <- as.numeric(warm_start_fit$final_lambda %||% lambda_target %||% 1.0)
    warm_lambda <- min(max(warm_lambda, 0), 1)
    Tmap <- tryCatch(
      build_transport(
        warm_theta,
        warm_w,
        reference_prior = reference_prior,
        lambda = warm_lambda,
        verbose = FALSE
      ),
      error = function(e) NULL
    )
    if (!is.null(Tmap) && isTRUE(verbose)) {
      cat("  Warm-start transport accepted.\n")
    }
  }
  if (is.null(Tmap)) {
    Tmap <- build_transport(
      Theta, w,
      reference_prior = reference_prior,
      lambda = 0,
      verbose = verbose
    )
  }
  Z <- Tmap$fwd(Theta)
  lpz <- as.numeric(base_logpdf(Theta) - Tmap$log_jac(Theta))
  lambda <- 0; round <- 0L
  log_evidence <- as.numeric(initial_log_normalizer %||% 0)
  lambda_hist <- c(0)
  ess_frac_hist <- numeric(0)
  resampled_hist <- logical(0)
  accept_rate_hist <- numeric(0)
  rw_accept_rate_hist <- numeric(0)
  pcn_accept_rate_hist <- numeric(0)
  indep_accept_rate_hist <- numeric(0)
  da_pass_rate_hist <- numeric(0)
  log_increment_hist <- numeric(0)
  target_acc_rw <- 0.234; log_rw_scale <- log(rw_scale_init); rm_gain <- 0.05
  # MCSE accumulator for log-evidence (sum of per-round variances)
  mcse_var_accum <- 0.0
  # pCN Robbins-Monro adaptation targeting 30% acceptance
  target_acc_pcn <- 0.30
  logit_pcn_beta <- qlogis(0.40)   # Start from beta = 0.4
  rm_gain_pcn    <- 0.05
  last_elite_mix <- NULL
  elite_history <- list()
  if (isTRUE(warm_start_use_mixture) && !is.null(warm_theta)) {
    Zw <- tryCatch(Tmap$fwd(warm_theta), error = function(e) NULL)
    if (!is.null(Zw) && nrow(Zw) >= max(10L, min(50L, ncol(Zw) * 2L))) {
      warm_G <- if (!is.null(warm_mix_seed) && !is.null(warm_mix_seed$meansZ)) {
        length(warm_mix_seed$meansZ)
      } else {
        min(as.integer(G_mix), max(2L, floor(sqrt(nrow(Zw)))))
      }
      set.seed(.seed_plan_block(seed_plan, "warm_elite_fit", default = seed + 2L))
      last_elite_mix <- tryCatch(
        fit_elite_mixture_Z(
          Zw,
          warm_w,
          elite_quantile = 0.4,
          G = as.integer(min(G_mix, warm_G)),
          cov_inflation = 4.0,
          warm_start_mixZ = warm_mix_seed,
          em_itmax = 6,
          housekeeping = TRUE,
          min_G_keep = 2,
          merge_thresh = 0.10,
          verbose = FALSE,
          lambda = min(max(as.numeric(lambda_target), 0), 1)
        ),
        error = function(e) NULL
      )
      if (!is.null(last_elite_mix) && isTRUE(verbose)) {
        cat("  Warm-start elite mixture accepted.\n")
      }
    }
  }
  # --- pCN kill-switch state (rolling) ---
  # Rolling acceptance snapshots (used for adaptive n_moves and routing probs)
  last_rw_acc  <- NA_real_
  last_pcn_acc <- NA_real_
  last_id_acc  <- NA_real_
  pcn_disabled <- !allow_pcn
  pcn_low_streak <- 0
  # --- Adaptive move/jitter knobs ---
  moves_cap <- 8L                 # hard cap on MCMC moves per round
  dup_high_thresh <- 0.20         # if >20% duplicates after resampling, prefer RW jitter first
  dup_low_thresh  <- 0.05         # if <=5% duplicates after first jitter, skip second jitter
  pcn_jitter_first_lambda <- 0.20 # below this λ, try pCN-first; otherwise RW-first
  # GSS state
  ref_mix <- NULL     # frozen working density r
  lp_ref  <- NULL     # current log r(theta) for each particle
  gss_snapshot_taken <- FALSE
  gss_refresh_taken <- FALSE
  transport_frozen <- FALSE
  # DA state across rounds (adapt α between rounds only)
  da_alpha <- da_alpha_init
  last_da_pass_rate <- NA_real_
  da_prop_cum <- 0L; da_pass_cum <- 0L
  # Transport refit state
  tr_state <- NULL  # Will be initialized by maybe_refit_transport
  # Exact normalization correction for target changes induced by GSS snapshot/refresh
  # or by dropping the frozen reference before a transport refit.
  apply_weight_correction <- function(w_cur, log_adjust, clamp_label, reset_msg) {
    w_base <- pmax(as.numeric(w_cur), 0)
    sw_base <- sum(w_base)
    if (!is.finite(sw_base) || sw_base <= 0) {
      w_base <- rep(1 / length(w_cur), length(w_cur))
    } else {
      w_base <- w_base / sw_base
    }

    log_adjust <- as.numeric(log_adjust)
    bad <- !is.finite(log_adjust)
    if (all(bad)) {
      return(list(
        ok = FALSE,
        w = rep(1 / length(w_base), length(w_base)),
        log_norm = NA_real_,
        var_log = 0.0,
        all_bad = TRUE,
        reset_msg = reset_msg
      ))
    }
    if (any(bad)) {
      log_adjust[bad] <- min(log_adjust[!bad]) - 50
      cat(sprintf("  [GSS Guard] Clamped %d non-finite %s values.\n", sum(bad), clamp_label))
    }

    amax <- max(log_adjust)
    u <- exp(log_adjust - amax)
    mu1 <- sum(w_base * u)
    if (!is.finite(mu1) || mu1 <= 0) {
      return(list(
        ok = FALSE,
        w = rep(1 / length(w_base), length(w_base)),
        log_norm = NA_real_,
        var_log = 0.0,
        all_bad = FALSE,
        reset_msg = reset_msg
      ))
    }
    mu2 <- sum(w_base * u * u)
    Neff <- 1 / sum(w_base * w_base)
    var_log <- (mu2 - mu1^2) / (max(Neff, 1) * max(mu1^2, .Machine$double.eps))

    logw_corr <- log(pmax(w_base, .Machine$double.eps)) + log_adjust
    lse_corr <- logsumexp(logw_corr)
    if (!is.finite(lse_corr)) {
      return(list(
        ok = FALSE,
        w = rep(1 / length(w_base), length(w_base)),
        log_norm = NA_real_,
        var_log = 0.0,
        all_bad = FALSE,
        reset_msg = reset_msg
      ))
    }

    w_new <- exp(logw_corr - lse_corr)
    sw_new <- sum(w_new)
    if (!is.finite(sw_new) || sw_new <= 0) {
      return(list(
        ok = FALSE,
        w = rep(1 / length(w_base), length(w_base)),
        log_norm = NA_real_,
        var_log = 0.0,
        all_bad = FALSE,
        reset_msg = reset_msg
      ))
    }
    w_new <- w_new / sw_new

    list(
      ok = TRUE,
      w = w_new,
      log_norm = lse_corr,
      var_log = max(var_log, 0),
      all_bad = FALSE,
      reset_msg = reset_msg
    )
  }

  while (lambda < lambda_target - 1e-12 && round < max_rounds) {
    round <- round + 1L
    cat(sprintf("\nRound %d: λ=%.3f -> ", round, lambda))
    resampled <- FALSE
    cess_target_round <- if (is.null(cess_target)) {
      cess_target_at_lambda(lambda)
    } else {
      as.numeric(cess_target)
    }
    # Guard against non-finite or invalid normalized weights.
    w[!is.finite(w) | w < 0] <- 0
    sw <- sum(w)
    if (!is.finite(sw) || sw <= 0) {
      w <- rep(1 / length(w), length(w))
      cat("  [Guard] Reset invalid weights to uniform.\n")
    } else {
      w <- w / sw
    }

    # --- choose step statistic h for rCESS (GSS-aware if ref is active) ---
    if (gss_enable && !is.null(ref_mix) && is.null(lp_ref)) {
      # safety: if ref_mix exists before snapshot logic below
      lp_ref <- log_r_theta(Theta, Z, Tmap, ref_mix)
    }
    h_step <- if (gss_enable && !is.null(ref_mix) && !is.null(lp_ref)) (loglik - lp_ref) else loglik
    if (!all(is.finite(h_step))) {
      bad <- !is.finite(h_step)
      finite_h <- h_step[!bad]
      if (!length(finite_h)) stop("All h_step values are non-finite; cannot continue tempering.")
      h_step[bad] <- min(finite_h)
      cat(sprintf("  [Guard] Replaced %d non-finite h_step values.\n", sum(bad)))
    }

    next_lambda <- if (gss_enable && !is.null(ref_mix) && !is.null(lp_ref)) {
      next_lambda_via_rCESS_stat(w, h_step, lambda, target = cess_target_round, lambda_target = lambda_target)
    } else {
      next_lambda_via_rCESS(w, h_step, lambda, target = cess_target_round, lambda_target = lambda_target)
    }
    next_checkpoint <- if (checkpoint_ptr <= length(checkpoint_lambdas)) checkpoint_lambdas[checkpoint_ptr] else Inf
    lambda_new  <- min(next_lambda, lambda_target, next_checkpoint)
    if (lambda_new <= lambda) {
      delta_floor <- min(1e-4, lambda_target - lambda)
      lambda_new <- lambda + delta_floor
      cat(sprintf("  [Guard] Applied lambda step floor (Δ=%.4g).\n", delta_floor))
    }
    delta       <- lambda_new - lambda
    lambda      <- lambda_new
    adapt_phase <- (lambda <= adapt_lambda_max)

    cat(sprintf("%.3f (Δ=%.4f, CESS target=%.3f)\n", next_lambda, delta, cess_target_round))

    # Update weights in log-space for numerical stability
    mll   <- max(h_step)
    lw_inc <- delta * (h_step - mll)                   # <= 0 by construction
    # --- MCSE(logZ) per-round contribution (delta-method, weighted) ---
    # u_i = exp(lw_inc_i); use PRE-update weights 'w'
    u   <- exp(lw_inc)
    mu1 <- sum(w * u)                                  # E_w[u]
    mu2 <- sum(w * u * u)                              # E_w[u^2]
    Neff <- 1 / sum(w * w)                             # current effective sample size
    var_logZ_inc <- (mu2 - mu1^2) / (max(Neff, 1) * max(mu1^2, .Machine$double.eps))
    mcse_var_accum <- mcse_var_accum + max(var_logZ_inc, 0)
    logw_new_raw <- log(pmax(w, .Machine$double.eps)) + lw_inc
    lse_new <- logsumexp(logw_new_raw)
    # Compute log evidence increment using log-sum-exp (GSS-aware)
    logZ_inc <- if (is.finite(lse_new)) delta * mll + lse_new else 0.0
    # Normalize weights in log domain
    if (!is.finite(lse_new)) {
      w_new <- rep(1 / length(w), length(w))
      logw_new <- rep(-log(length(w)), length(w))
      cat("  [Guard] Non-finite log-weight normalization; reset to uniform.\n")
    } else {
      logw_new <- logw_new_raw - lse_new
      w_new    <- exp(logw_new)
      if (!all(is.finite(w_new)) || sum(w_new) <= 0) {
        w_new <- rep(1 / length(w), length(w))
        logw_new <- rep(-log(length(w)), length(w))
        cat("  [Guard] Non-finite normalized weights; reset to uniform.\n")
      } else {
        w_new <- w_new / sum(w_new)
      }
    }
    ess <- ESS(w_new)
    if (!is.finite(ess)) ess <- 0
    ess_frac <- ess / length(w_new)
    if (!is.finite(ess_frac)) ess_frac <- 0
    log_evidence <- log_evidence + logZ_inc
    pred_rcess_full <- rCESS_stat(w, h_step, lambda_target - lambda)
    cat(sprintf("  rCESS(remain)=%.3f | target=%.3f\n", pred_rcess_full, cess_target_round))
    cat(sprintf("  ESS(pre)=%.3f | logZ += %.4f -> %.4f\n", ess_frac, logZ_inc, log_evidence))
    if (ess_frac < resample_threshold) {
      # Pre-resample light rejuvenation (SKIP when λ is tiny to avoid wasted ll calls)
      if (ess_frac > 0.25 && lambda >= pre_resample_lambda_gate) {
        cat("  Pre-resample light rejuvenation (pCN+RW)...\n")
        pcn_beta_curr <- plogis(logit_pcn_beta)
        emix_pre <- if (!is.null(last_elite_mix) && !.is_empty_mix(last_elite_mix)) {
          last_elite_mix
        } else {
          .default_std_normal_mix(ncol(Z))
        }
        move_pre <- mcmc_moves_z_mix_batched(
          Z, loglik, lpz, Tmap, lambda,
          mu_ref, prior_L, w_new,
          elite_mix = emix_pre, hist_mix = NULL,
          reference_prior = reference_prior,
          data = NULL, loglik_fn = bridge_stat_kernel,
          n_moves = 1,
          rw_prob = 0.40, rw_scale = 0.35,
          pcn_prob = 0.60, pcn_beta = pcn_beta_curr,
          indep_t_df = indep_t_df, indep_t_prob = indep_t_prob,
          seed = .seed_plan_block(seed_plan, "pre_move", round = round, default = seed + round * 19L),
          param_names = param_names,
          weak_dim_idx = integer(0), rw_expand_factor = 1.0,
          resampled = FALSE,
          base_logpdf_fn = base_logpdf,
          allow_pcn = allow_pcn,
          n_cores = n_cores
        )
        Z <- move_pre$Z; loglik <- move_pre$loglik; lpz <- move_pre$lpz
        Theta <- Tmap$inv(Z); colnames(Theta) <- param_names
        if (gss_enable && !is.null(ref_mix)) {
          lp_ref <- log_r_theta(Theta, Z, Tmap, ref_mix)
        }
      }
      sort_used <- "hilbert"
      if (resample_sort_mode == "none") {
        ord <- seq_len(nrow(Z))
        sort_used <- "none"
      } else if (resample_sort_mode == "cheap1d") {
        ord <- .cheap_sort_order(Z)
        sort_used <- "cheap1d"
      } else if (resample_sort_mode == "hilbert") {
        ord <- hilbert_sort_order(Z, bits = 16L)
        sort_used <- "hilbert"
      } else {
        # Adaptive default: spend Hilbert cost only in hard low-ESS regimes.
        if (ess_frac <= hilbert_hard_ess) {
          ord <- hilbert_sort_order(Z, bits = 16L)
          sort_used <- "hilbert"
        } else {
          ord <- .cheap_sort_order(Z)
          sort_used <- "cheap1d"
        }
      }
      cat(sprintf("  Resampling (sorted stratified, order=%s)...\n", sort_used))
      if (!isTRUE(deterministic_resampling)) {
        set.seed(.seed_plan_block(seed_plan, "resample", round = round, default = seed + round * 23L))
      }
      # Stratified resampling on chosen ordering, then map back to original indices.
      # Stratified resampling on sorted weights, then map back
      idx_sorted <- stratified_resample_sorted(
        w_new[ord],
        deterministic = isTRUE(deterministic_resampling)
      )
      idx <- ord[idx_sorted]
      Theta <- Theta[idx, , drop = FALSE]; Z <- Z[idx, , drop = FALSE]
      loglik <- loglik[idx]; lpz <- lpz[idx]
      w_new <- rep(1/length(w_new), length(w_new))
      resampled <- TRUE
    }
    # Keep the clipped/current lambda in state; `next_lambda` can lie beyond a
    # checkpoint or guard-floor adjustment.
    w <- w_new
    lambda_hist <- c(lambda_hist, lambda)
    ess_frac_hist <- c(ess_frac_hist, ess_frac)
    resampled_hist <- c(resampled_hist, resampled)
    log_increment_hist <- c(log_increment_hist, logZ_inc)

    # Post-resampling jitter now handled by the single adaptive branch below.
    if (resampled) {
      ## --- Decide jitter order after resampling ---
      ## If many exact duplicates or λ is not tiny -> RW-first; else pCN-first.
      dup0 <- 1 - (length(unique(idx)) / length(idx))
      do_rw_first <- (!allow_pcn) || (lambda >= pcn_jitter_first_lambda) || (dup0 > dup_high_thresh) ||
                     (is.finite(last_pcn_acc) && last_pcn_acc < 0.10)

      if (do_rw_first) {
        cat(sprintf("  Pre-move jitter (RW first)... [dup0=%.1f%%]\n", 100*dup0))
        mover1 <- mcmc_moves_z_mix_batched(
          Z, loglik, lpz, Tmap, lambda,
          mu_ref, prior_L, w,
          elite_mix = .default_std_normal_mix(ncol(Z)), hist_mix = NULL,
          reference_prior = reference_prior,
          data = NULL, loglik_fn = bridge_stat_kernel,
          n_moves = 1,
          rw_prob = 1, rw_scale = 0.35,
          pcn_prob = 0, pcn_beta = plogis(logit_pcn_beta),
          indep_t_df = indep_t_df, indep_t_prob = 0,
          seed = .seed_plan_block(seed_plan, "jitter1", round = round, default = seed + round * 37L),
          param_names = param_names,
          weak_dim_idx = integer(0), rw_expand_factor = 1.0,
          resampled = TRUE,
          base_logpdf_fn = base_logpdf,
          allow_pcn = allow_pcn,
          n_cores = n_cores)
        Z <- mover1$Z; loglik <- mover1$loglik; lpz <- mover1$lpz
        Theta <- Tmap$inv(Z); colnames(Theta) <- param_names
        if (gss_enable && !is.null(ref_mix)) lp_ref <- log_r_theta(Theta, Z, Tmap, ref_mix)

        dup1 <- 1 - (nrow(unique(Z)) / nrow(Z))
        need_pcn <- allow_pcn && (dup1 > dup_low_thresh) && (lambda < 0.40) && (mover1$rw_accept_rate < 0.35)
        if (need_pcn) {
          cat(sprintf("  Pre-move jitter (pCN second)... [dup1=%.1f%%]\n", 100*dup1))
          movej <- mcmc_moves_z_mix_batched(
            Z, loglik, lpz, Tmap, lambda,
            mu_ref, prior_L, w,
            elite_mix = .default_std_normal_mix(ncol(Z)), hist_mix = NULL,
            reference_prior = reference_prior,
            data = NULL, loglik_fn = bridge_stat_kernel,
            n_moves = 1,
            rw_prob = 0, rw_scale = 0.5,
            pcn_prob = 1, pcn_beta = plogis(logit_pcn_beta),
            indep_t_df = indep_t_df, indep_t_prob = 0,
            seed = .seed_plan_block(seed_plan, "jitter2", round = round, default = seed + round * 31L),
            param_names = param_names,
            weak_dim_idx = integer(0), rw_expand_factor = 1.0,
            resampled = TRUE,
            base_logpdf_fn = base_logpdf,
            allow_pcn = allow_pcn,
            n_cores = n_cores)
          Z <- movej$Z; loglik <- movej$loglik; lpz <- movej$lpz
          Theta <- Tmap$inv(Z); colnames(Theta) <- param_names
          if (gss_enable && !is.null(ref_mix)) lp_ref <- log_r_theta(Theta, Z, Tmap, ref_mix)
          last_pcn_acc <- movej$pcn_accept_rate
        } else {
          cat(sprintf("  Pre-move jitter (pCN) skipped [dup1=%.1f%%, λ=%.3f]\n", 100*dup1, lambda))
        }
        last_rw_acc <- mover1$rw_accept_rate
      } else {
        cat("  Pre-move jitter (pCN first)...\n")
        movej <- mcmc_moves_z_mix_batched(
          Z, loglik, lpz, Tmap, lambda,
          mu_ref, prior_L, w,
          elite_mix = .default_std_normal_mix(ncol(Z)), hist_mix = NULL,
          reference_prior = reference_prior,
          data = NULL, loglik_fn = bridge_stat_kernel,
          n_moves = 1,
          rw_prob = 0, rw_scale = 0.5,
          pcn_prob = 1, pcn_beta = plogis(logit_pcn_beta),
          indep_t_df = indep_t_df, indep_t_prob = 0,
          seed = .seed_plan_block(seed_plan, "jitter1", round = round, default = seed + round * 31L),
          param_names = param_names,
          weak_dim_idx = integer(0), rw_expand_factor = 1.0,
          resampled = TRUE,
          base_logpdf_fn = base_logpdf,
          allow_pcn = allow_pcn,
          n_cores = n_cores)
        Z <- movej$Z; loglik <- movej$loglik; lpz <- movej$lpz
        Theta <- Tmap$inv(Z); colnames(Theta) <- param_names
        if (gss_enable && !is.null(ref_mix)) lp_ref <- log_r_theta(Theta, Z, Tmap, ref_mix)
        last_pcn_acc <- movej$pcn_accept_rate

        dup1 <- 1 - (nrow(unique(Z)) / nrow(Z))
        do_rw_jitter <- (movej$pcn_accept_rate < 0.15) || (dup1 > dup_low_thresh)
        if (do_rw_jitter) {
          cat(sprintf("  Pre-move jitter (RW second)... [pcn_acc=%.3f, dup1=%.1f%%]\n",
                      movej$pcn_accept_rate, 100*dup1))
          mover <- mcmc_moves_z_mix_batched(
          Z, loglik, lpz, Tmap, lambda,
          mu_ref, prior_L, w,
          elite_mix = .default_std_normal_mix(ncol(Z)), hist_mix = NULL,
          reference_prior = reference_prior,
          data = NULL, loglik_fn = bridge_stat_kernel,
          n_moves = 1,
          rw_prob = 1, rw_scale = 0.35,
          pcn_prob = 0, pcn_beta = plogis(logit_pcn_beta),
          indep_t_df = indep_t_df, indep_t_prob = 0,
          seed = .seed_plan_block(seed_plan, "jitter2", round = round, default = seed + round * 37L),
          param_names = param_names,
          weak_dim_idx = integer(0), rw_expand_factor = 1.0,
          resampled = TRUE,
          base_logpdf_fn = base_logpdf,
          allow_pcn = allow_pcn,
          n_cores = n_cores)
          Z <- mover$Z; loglik <- mover$loglik; lpz <- mover$lpz
          Theta <- Tmap$inv(Z); colnames(Theta) <- param_names
          if (gss_enable && !is.null(ref_mix)) lp_ref <- log_r_theta(Theta, Z, Tmap, ref_mix)
          last_rw_acc <- mover$rw_accept_rate
        } else {
          cat(sprintf("  Pre-move jitter (RW) skipped [pcn_acc=%.3f, dup1=%.1f%%]\n",
                      movej$pcn_accept_rate, 100*dup1))
        }
      }
    }

    # Adaptive transport refit with diagnostics.
    # Stability policy:
    # - only adapt maps in the early phase
    # - freeze map adaptation once GSS has been activated
    allow_transport_refit <- adapt_phase && !(freeze_transport_after_gss && transport_frozen)
    ref <- if (allow_transport_refit) {
      maybe_refit_transport(
        Theta = Theta, w = w, lambda = lambda, round = round,
        resampled = resampled, ess_frac = ess_frac,
        Tmap = Tmap, tr_state = tr_state,
        every_rounds = refit_every,     # advisory cadence; gate is relative
        verbose = verbose
      )
    } else {
      list(refit = FALSE, Tmap = Tmap, tr_state = tr_state)
    }
    if (isTRUE(ref$refit)) {
      # If a GSS reference was active, remove its fixed-lambda contribution
      # before invalidating it due to transport-coordinate change.
      if (gss_enable && !is.null(lp_ref)) {
        lp_ref_old <- as.numeric(lp_ref)
        corr <- apply_weight_correction(
          w_cur = w,
          log_adjust = -(1 - lambda) * lp_ref_old,
          clamp_label = "old lp_ref at refit",
          reset_msg = "  [GSS Guard] Refit deweight normalization non-finite; reset weights.\n"
        )
        if (!isTRUE(corr$ok)) {
          w <- corr$w
          cat(corr$reset_msg)
        } else {
          w <- corr$w
          log_evidence <- log_evidence + corr$log_norm
          mcse_var_accum <- mcse_var_accum + corr$var_log
          cat(sprintf("  [GSS] removed frozen reference before refit (logZ += %.4f)\n", corr$log_norm))
        }
      }
      Tmap <- ref$Tmap
      Z    <- Tmap$fwd(Theta)
      lpz  <- as.numeric(base_logpdf(Theta) - Tmap$log_jac(Theta))
      # Transport changed => all Z-space objects are invalid under old coordinates.
      last_elite_mix <- NULL
      elite_history  <- list()
      ref_mix <- NULL
      lp_ref <- NULL
    }
    tr_state <- ref$tr_state

    # Elite mixture (slightly more inflated early to avoid peaky elites)
    cat("  Fit elite mixture\n")
    w_fit <- (w^gamma_sharp); w_fit <- w_fit / sum(w_fit)
    elite_q_eff <- max(0.20, 0.4 * (1 - 0.5 * lambda))
    G_eff <- max(2L, min(G_mix, 2L + floor(6 * lambda)))
    merge_thresh_eff <- if (lambda >= 0.85) 0.06 else if (lambda >= 0.70) 0.08 else 0.10
    cov_inflation_eff <- max(4.0, 3.0 + 2.0 * (1 - lambda))   # Keep minimum 4.0 late in annealing
    # (2,4) Early rounds: broaden elite/history mixtures (+25% inflation before λ≈0.3)
    if (lambda < 0.30) cov_inflation_eff <- cov_inflation_eff * 1.25
    set.seed(.seed_plan_block(seed_plan, "elite_fit", round = round, default = seed + round * 29L))
    elite_mix <- tryCatch(
      fit_elite_mixture_Z(Z, w_fit, elite_quantile = elite_q_eff, G = G_eff,
                          cov_inflation = cov_inflation_eff, warm_start_mixZ = last_elite_mix,
                          em_itmax = 8, housekeeping = TRUE, min_G_keep = 2,
                          merge_thresh = merge_thresh_eff, verbose = FALSE,
                          lambda = lambda),  # Pass lambda to control min eigenvalue
      error = function(e) NULL
    )
    if (is.null(elite_mix) || .is_empty_mix(elite_mix)) {
      mu <- colSums(Z * w_fit)
      Zc <- sweep(Z, 2L, mu, `-`)
      Sig <- t(Zc) %*% (Zc * w_fit)
      Sig <- as.matrix((Sig + t(Sig)) / 2)
      Sig <- as.matrix(Matrix::nearPD(Sig, conv.tol = 1e-7)$mat)
      elite_mix <- list(meansZ = list(as.numeric(mu)), covsZ = list(Sig), wZ = 1)
      elite_mix$cache <- prep_mix_cache(elite_mix$meansZ, elite_mix$covsZ, elite_mix$wZ)
    }
    last_elite_mix <- elite_mix

    # History mixture
    hist_mix <- NULL
    if (hist_mix_enable) {
      hm_wrap <- build_hist_mixture(elite_history, lambda,
                                    hist_mix_lambda_thresh, hist_mix_prob,
                                min_G_keep = 2, merge_thresh = merge_thresh_eff)
      if (!is.null(hm_wrap)) hist_mix <- hm_wrap
    }

    # --- On-the-fly GSS: snapshot/refresh ref_mix and adjust weights if needed ---
    if (gss_enable) {
      up <- maybe_update_ref_mix(lambda, elite_mix, hist_mix, ref_mix,
                                 allow_refresh = allow_ref_refresh,
                                 lambda_ref_snap = lambda_ref_snap,
                                 lambda_ref_refresh = lambda_ref_refresh,
                                 d = ncol(Z),
                                 snapshot_taken = gss_snapshot_taken,
                                 refresh_taken = gss_refresh_taken,
                                 single_update_policy = single_gss_updates)
      gss_snapshot_taken <- isTRUE(up$snapshot_taken)
      gss_refresh_taken <- isTRUE(up$refresh_taken)
      if (isTRUE(up$switched)) {
        ref_mix_old <- ref_mix
        ref_mix <- up$ref_mix
        # compute new lp_ref at current lambda
        lp_ref_new <- log_r_theta(Theta, Z, Tmap, ref_mix)
        bad_ref <- !is.finite(lp_ref_new)
        if (all(bad_ref)) {
          cat("  [GSS] reference update dropped (all lp_ref non-finite).\n")
          ref_mix <- NULL
          lp_ref <- NULL
        } else {
          if (is.null(lp_ref)) {
            # first snapshot: incorporate (1 - lambda) * log r into current target
            corr <- apply_weight_correction(
              w_cur = w,
              log_adjust = (1 - lambda) * lp_ref_new,
              clamp_label = "lp_ref",
              reset_msg = "  [GSS Guard] Snapshot normalization non-finite; reset weights.\n"
            )
            if (!isTRUE(corr$ok)) {
              w <- corr$w
              cat(corr$reset_msg)
            } else {
              w <- corr$w
              log_evidence <- log_evidence + corr$log_norm
              mcse_var_accum <- mcse_var_accum + corr$var_log
            }
            lp_ref <- lp_ref_new
            cat("  [GSS] reference SNAPSHOT at lambda=", sprintf("%.3f", lambda),
                sprintf(" (weights updated, logZ += %.4f)\n", if (isTRUE(corr$ok)) corr$log_norm else 0))
            if (freeze_transport_after_gss) transport_frozen <- TRUE
          } else {
            # refresh: exact corrective reweighting at fixed lambda
            corr <- apply_weight_correction(
              w_cur = w,
              log_adjust = (1 - lambda) * (lp_ref_new - lp_ref),
              clamp_label = "lp_ref",
              reset_msg = "  [GSS Guard] Refresh normalization non-finite; reset weights.\n"
            )
            if (!isTRUE(corr$ok)) {
              w <- corr$w
              cat(corr$reset_msg)
            } else {
              w <- corr$w
              log_evidence <- log_evidence + corr$log_norm
              mcse_var_accum <- mcse_var_accum + corr$var_log
            }
            lp_ref <- lp_ref_new
            cat("  [GSS] reference REFRESH at lambda=", sprintf("%.3f", lambda),
                sprintf(" (weights corrected, logZ += %.4f)\n", if (isTRUE(corr$ok)) corr$log_norm else 0))
          }
        }
      } else if (!is.null(ref_mix) && is.null(lp_ref)) {
        # safety: if ref_mix was injected externally before entering loop
        lp_ref <- log_r_theta(Theta, Z, Tmap, ref_mix)
        bad_ref <- !is.finite(lp_ref)
        if (all(bad_ref)) {
          cat("  [GSS] safety snapshot dropped (all lp_ref non-finite).\n")
          ref_mix <- NULL
          lp_ref <- NULL
        } else {
          corr <- apply_weight_correction(
            w_cur = w,
            log_adjust = (1 - lambda) * lp_ref,
            clamp_label = "safety lp_ref",
            reset_msg = "  [GSS Guard] Safety normalization non-finite; reset weights.\n"
          )
          if (!isTRUE(corr$ok)) {
            w <- corr$w
            cat(corr$reset_msg)
          } else {
            w <- corr$w
            log_evidence <- log_evidence + corr$log_norm
            mcse_var_accum <- mcse_var_accum + corr$var_log
          }
        }
      }
    }


    indep_t_df_eff <- indep_t_df
    # (2,4) Early rounds: heavier t tails for independence; later tone down
    indep_t_prob_eff <-
      if (lambda < 0.30) max(indep_t_prob, 0.90)
      else if (lambda > 0.70) min(indep_t_prob, 0.75)
      else indep_t_prob
    # Always blend small base mass into mixtures; stronger late blend for tail coverage
    base_mix <- .default_std_normal_mix(ncol(Z))
    eps_elite <- if (lambda < 0.30) 0.12 else 0.10
    eps_hist  <- if (lambda < 0.30) 0.09 else 0.07
    if (!is.null(elite_mix) && !.is_empty_mix(elite_mix)) {
      elite_mix <- blend_mixes(elite_mix, base_mix, eps = eps_elite)
    }
    if (!is.null(hist_mix) && !is.null(hist_mix$mix) && !.is_empty_mix(hist_mix$mix)) {
      hist_mix$mix <- blend_mixes(hist_mix$mix, base_mix, eps = eps_hist)
    }

    # --- DA screen (snapshot for this round) ---
    da_screen_mix <- if (da_enable) .build_da_screen_mix(elite_mix, hist_mix, ncol(Z)) else NULL
    da_calib_active <- da_enable && adapt_phase && !is.null(da_screen_mix)
    # Micro-calibration (tiny exact budget) to tighten stage-1 (early adaptive phase only).
    da_calib_info <- if (da_calib_active) {
      .calibrate_da_surrogate(Z, Theta, lpz, lambda, if (gss_enable) ref_mix else NULL,
                              da_screen_mix, NULL, bridge_stat_kernel, Tmap,
                              ll_cache = ll_cache,
                              n = da_calibrate_n,
                              seed = .seed_plan_block(seed_plan, "da_calib", round = round, default = seed + 7L * round),
                              alpha = da_alpha,
                              lambda_floor = da_lambda_floor,
                              gate_lambda = da_lambda_floor,
                              n_cores = n_cores)
    } else list(
      a = 0.0, b = 1.0, r2 = NA_real_, skipped = TRUE,
      reason = if (!da_enable) "DA off"
               else if (!adapt_phase) "post-adaptation phase"
               else "no screen mix"
    )
    if (da_enable) {
      if (isTRUE(da_calib_info$skipped)) {
        cat(sprintf("  [DA] α=%.3f, calib=skipped (%s)\n", da_alpha, da_calib_info$reason))
      } else {
        cat(sprintf("  [DA] α=%.3f, slope b=%.3f, R²=%.2f\n", da_alpha, da_calib_info$b, da_calib_info$r2))
      }
    }

    # Weak dimension inflation with temperature-aware factor, strong floor late in annealing
    weak_idx <- weak_dims_from_Z(Z, w, frac = 0.25, min_keep = 1L)
    factor_eff <- 4.0 + 1.5 * (1 - lambda)   # Early: 5.5, Late: 4.0
    if (lambda < 0.30) factor_eff <- factor_eff * 1.25  # (4) more breadth early
    if (!is.null(elite_mix) && !.is_empty_mix(elite_mix)) {
      elite_mix <- inflate_mixture_along_dims(elite_mix, weak_idx, factor = factor_eff)
    }
    if (!is.null(hist_mix) && !is.null(hist_mix$mix) && !.is_empty_mix(hist_mix$mix)) {
      hist_mix$mix <- inflate_mixture_along_dims(hist_mix$mix, weak_idx, factor = factor_eff)
    }

    # MCMC moves
    rw_scale <- exp(log_rw_scale)
        # Base channel probabilities (before guards)
    pcn_prob_eff <- if (resampled) {
      max(0.50, 0.90 * (1 - lambda))
    } else {
      .clamp(0.10 + 0.20*(1 - lambda), 0.10, 0.40)
    }
    rw_prob_eff  <- max(0.10, min(0.70, rw_prob * (1 - 0.5 * lambda)))
    # Late, avoid wasting mass on independence if it's weak; early ensure some ID mass.
    id_min <- if (lambda > 0.70) 0.10 else 0.25
    id_prob_eff <- 1 - (pcn_prob_eff + rw_prob_eff)
    if (id_prob_eff < id_min) {
      rw_prob_eff  <- max(0, rw_prob_eff - (id_min - id_prob_eff))
      id_prob_eff  <- 1 - (pcn_prob_eff + rw_prob_eff)
    }
    # (1) pCN kill-switch: if rolling acceptance < 0.05 for 2 rounds, turn pCN off permanently
    if (pcn_disabled) {
      # Split former pCN mass between RW and ID; give ID more if it was accepting
      share_to_rw <- if (is.finite(last_id_acc) && last_id_acc >= 0.10) 0.60 else 0.80
      rw_prob_eff <- min(0.70, rw_prob_eff + pcn_prob_eff * share_to_rw)
      pcn_prob_eff <- 0.0
      id_prob_eff  <- 1 - rw_prob_eff
    }
    ## --- Early adaptive phase vs fixed-kernel phase ---
    if (adapt_phase) {
      # Calmer, capped adaptive number of MCMC moves.
      difficulty <- 1 - pred_rcess_full
      sDelta <- .clamp(delta / 0.05, 0.6, 1.4)             # Δ-driven scaling
      sDiff  <- 1 + .clamp(1.5 * difficulty, 0, 1.5)       # up to 2.5×
      acc_proxy <- mean(c(last_rw_acc, last_pcn_acc), na.rm = TRUE)
      if (!is.finite(acc_proxy)) acc_proxy <- 0.25
      sAcc <- if (acc_proxy < 0.10) 0.80 else              # if proposals are failing, don't throw more moves at it
              if (acc_proxy < 0.20) 0.90 else 1.00
      base_moves <- n_mcmc_moves
      floor_moves <- if (ess_frac < 0.40) 4L else if (ess_frac < 0.70) 3L else 2L
      n_moves_eff <- as.integer(ceiling(base_moves * sDelta * sDiff * sAcc))
      n_moves_eff <- max(floor_moves, min(n_moves_eff, moves_cap))
      if (lambda > 0.90) n_moves_eff <- max(n_moves_eff, 3L)
    } else {
      # Fixed-kernel phase: freeze probabilities/scales/move count for stability.
      pcn_prob_eff <- if (pcn_disabled) 0.0 else 0.20
      rw_prob_eff  <- .clamp(rw_prob, 0.10, 0.75)
      id_prob_eff  <- max(0, 1 - (pcn_prob_eff + rw_prob_eff))
      if (id_prob_eff < 0.05) {
        rw_prob_eff <- max(0.05, rw_prob_eff - (0.05 - id_prob_eff))
      }
      n_moves_eff <- as.integer(max(1L, post_adapt_n_mcmc_moves))
    }
    # Make very-early rounds cheap: cap to 2 moves and remove ID mass.
    if (lambda < da_lambda_floor) {
      n_moves_eff <- min(n_moves_eff, 2L)
      # Reassign independence mass to RW (ID is usually wasted early)
      id_prob_eff <- 0.0
      rw_prob_eff <- min(0.80, 1 - pcn_prob_eff)
    }
    cat(sprintf("  Moves: %d [pcn=%.2f, rw=%.2f, rw_scale=%.3f]\n", n_moves_eff, pcn_prob_eff, rw_prob_eff, rw_scale))

    # Use current adaptive pCN beta in the main moves
    pcn_beta_curr <- plogis(logit_pcn_beta)
    move <- mcmc_moves_z_mix_batched(Z, loglik, lpz, Tmap, lambda,
                                         mu_ref, prior_L, w,
                                         elite_mix, hist_mix,
                                         reference_prior = reference_prior,
                                         data = NULL, loglik_fn = bridge_stat_kernel, n_moves = n_moves_eff,
                                         rw_prob = rw_prob_eff, rw_scale = rw_scale,
                                         pcn_prob = pcn_prob_eff, pcn_beta = pcn_beta_curr,
                                         indep_t_df = indep_t_df_eff, indep_t_prob = indep_t_prob_eff,
                                         seed = .seed_plan_block(seed_plan, "main_move", round = round, default = seed + round * 97),
                                          param_names = param_names,
                                          weak_dim_idx = weak_idx,
                                          rw_expand_factor = 2.5,
                                          ref_mix = ref_mix,
                                          lp_ref_vec = lp_ref,
                                          # DA snapshot params
                                          da_enable = da_enable,
                                          da_lambda_floor = da_lambda_floor,
                                          da_alpha = da_alpha,
                                          da_calib = c(da_calib_info$a, da_calib_info$b),
                                          da_screen_mix = da_screen_mix,
                                          ll_cache = ll_cache,
                                          base_logpdf_fn = base_logpdf,
                                          allow_pcn = allow_pcn,
                                          n_cores = n_cores)
    Z <- move$Z; loglik <- move$loglik; lpz <- move$lpz
    Theta <- Tmap$inv(Z); colnames(Theta) <- param_names
    # keep lp_ref in sync if GSS active and ref_mix frozen
    if (gss_enable && !is.null(ref_mix)) {
      lp_ref <- log_r_theta(Theta, Z, Tmap, ref_mix)
    }
    # Adapt RW scale and pCN beta via Robbins-Monro (early adaptive phase only).
    if (adapt_phase) {
      log_rw_scale <- .clamp(log_rw_scale + rm_gain * (move$rw_accept_rate - target_acc_rw),
                             log(0.05), log(1.5))
      logit_pcn_beta <- .clamp(logit_pcn_beta + rm_gain_pcn * (move$pcn_accept_rate - target_acc_pcn),
                               qlogis(0.05), qlogis(0.95))
    }
    # Update rolling acceptance snapshots for the next round's adaptive n_moves
    last_rw_acc <- move$rw_accept_rate; last_pcn_acc <- move$pcn_accept_rate; last_id_acc <- move$indep_accept_rate
    accept_rate_hist <- c(accept_rate_hist, move$accept_rate)
    rw_accept_rate_hist <- c(rw_accept_rate_hist, move$rw_accept_rate)
    pcn_accept_rate_hist <- c(pcn_accept_rate_hist, move$pcn_accept_rate)
    indep_accept_rate_hist <- c(indep_accept_rate_hist, move$indep_accept_rate)
    da_pass_rate_hist <- c(da_pass_rate_hist, move$da_pass_rate)
    # ---- DA diagnostics & α update ----
    if (da_enable && !is.na(move$da_pass_rate)) {
      last_da_pass_rate <- move$da_pass_rate
      da_prop_cum <- da_prop_cum + move$da_prop
      da_pass_cum <- da_pass_cum + move$da_pass
      # per-kernel pass rates
      pr_rw  <- if (move$da_prop_rw  > 0) move$da_pass_rw  / move$da_prop_rw  else NA_real_
      pr_pcn <- if (move$da_prop_pcn > 0) move$da_pass_pcn / move$da_prop_pcn else NA_real_
      pr_id  <- if (move$da_prop_id  > 0) move$da_pass_id  / move$da_prop_id  else NA_real_
      # conditional acceptance among those that passed
      acc_cnt <- move$rw_accept_rate*move$prop_rw + move$pcn_accept_rate*move$prop_pcn + move$indep_accept_rate*move$prop_id
      cond_acc <- if (move$da_pass > 0) acc_cnt / move$da_pass else NA_real_
      saved_frac <- 1 - last_da_pass_rate
      cat(sprintf("  [DA] pass=%.2f (rw=%.2f, pcn=%.2f, id=%.2f) | cond=%.2f | saved≈%.0f%% | α=%.3f, b=%.3f\n",
                  last_da_pass_rate, pr_rw, pr_pcn, pr_id, cond_acc, 100*saved_frac, da_alpha, da_calib_info$b))
      if (adapt_phase) {
        # RM update on log α (keep α within [0.6, 1.5])
        log_alpha <- log(da_alpha)
        log_alpha <- log_alpha + da_rm_gain * (da_target_pass - last_da_pass_rate)
        da_alpha  <- exp(.clamp(log_alpha, log(0.6), log(1.5)))
        cat(sprintf("  [DA] α tuned → %.3f\n", da_alpha))
      }
    }
    # (1) Update pCN kill-switch streak
    if (!pcn_disabled) {
      if (move$pcn_accept_rate < 0.05) pcn_low_streak <- pcn_low_streak + 1L else pcn_low_streak <- 0L
      if (pcn_low_streak >= 2L) {
        pcn_disabled <- TRUE
        cat("  [pCN] kill-switch engaged: reassigning mass to RW/ID for the remainder.\n")
      }
    }

    # Compute KL divergence proxy for diagnostics
    kl_proxy <- mean(lpz) - mean(base_logpdf(Theta))
    cat(sprintf("  acc[pCN=%.2f, rw=%.2f, id=%.2f] | KLproxy=%.3f\n",
                move$pcn_accept_rate, move$rw_accept_rate,
                move$indep_accept_rate, kl_proxy))

    # Save elite mixture to history
    if (lambda > 0.5 && !is.null(last_elite_mix)) {
      elite_history[[length(elite_history) + 1L]] <- last_elite_mix
      if (length(elite_history) > 6L) elite_history <- tail(elite_history, 6L)
    }
    while (checkpoint_ptr <= length(checkpoint_lambdas) &&
           lambda >= checkpoint_lambdas[checkpoint_ptr] - 1e-12) {
      checkpoint_log_evidence[checkpoint_ptr] <- log_evidence
      checkpoint_log_increments[checkpoint_ptr] <- if (checkpoint_ptr == 1L) {
        checkpoint_log_evidence[checkpoint_ptr]
      } else {
        checkpoint_log_evidence[checkpoint_ptr] - checkpoint_log_evidence[checkpoint_ptr - 1L]
      }
      checkpoint_ptr <- checkpoint_ptr + 1L
    }
  }

  # Final MCSE(logZ) from accumulated per-round variances
  mcse_logZ <- sqrt(mcse_var_accum)
  cat(sprintf("\nDone in %d rounds. Final λ=%.3f/%.3f | logZ≈%.4f ± %.4f (MCSE)\n",
              round, lambda, lambda_target, log_evidence, mcse_logZ))
  list(
    Theta = Theta, Z = Z, loglik = loglik, w = w, transport = Tmap,
    ref_mix = ref_mix, lp_ref = lp_ref,
    final_lambda = lambda, log_evidence = log_evidence,
    mcse_logZ = mcse_logZ, lpz = lpz,
    elite_mix_final = last_elite_mix,
    elite_history = elite_history,
    checkpoint_lambdas = checkpoint_lambdas,
    checkpoint_log_evidence = checkpoint_log_evidence,
    checkpoint_log_increments = checkpoint_log_increments,
    seed_plan = seed_plan,
    meta = list(rounds = round, ess = ESS(w), lambda_hist = lambda_hist,
                ess_frac_hist = ess_frac_hist,
                resampled_hist = resampled_hist,
                accept_rate_hist = accept_rate_hist,
                rw_accept_rate_hist = rw_accept_rate_hist,
                pcn_accept_rate_hist = pcn_accept_rate_hist,
                indep_accept_rate_hist = indep_accept_rate_hist,
                da_pass_rate_hist = da_pass_rate_hist,
                log_increment_hist = log_increment_hist,
                rw_scale_final = exp(log_rw_scale),
                pcn_beta_final = plogis(logit_pcn_beta),
                da_alpha_final = da_alpha,
                da_last_pass = last_da_pass_rate,
                da_pass_cum = if (da_prop_cum>0) da_pass_cum/da_prop_cum else NA_real_,
                ll_cache_size = if (!is.null(ll_cache)) length(ls(ll_cache$env, all.names = TRUE)) else 0L)
  )
}
