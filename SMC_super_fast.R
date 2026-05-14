# ============================================================================
# Per-group / local SMC sampler
# - One tempered target: reference_prior(theta) * likelihood(theta)^lambda
# - One transport: sparse triangular Gaussianizing map
# - Rejuvenation: random-walk, pCN, and current elite-mixture independence moves
# ============================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(matrixStats)
  library(mvtnorm)
  library(qrng)
})

source("smc_core.R")
if (!exists("normalize_reference_prior", mode = "function") ||
    !exists("reference_prior_logpdf", mode = "function") ||
    !exists("reference_prior_sample", mode = "function")) {
  source("reference_priors.R")
}

# ----------------------------- transport utilities -----------------------------

.make_tail_safe_cdf_bundle <- function(ux, cw) {
  stopifnot(length(ux) >= 2L, length(ux) == length(cw))
  eps <- max(1e-12, 1e-8 * diff(range(ux)))
  for (k in 2:length(ux)) {
    if (ux[k] <= ux[k - 1L]) ux[k] <- ux[k - 1L] + eps
  }

  xmin <- ux[1L]
  xmax <- ux[length(ux)]
  uL <- pmin(pmax(cw[1L], 1e-12), 1 - 1e-12)
  uR <- pmin(pmax(cw[length(cw)], 1e-12), 1 - 1e-12)
  vL <- qnorm(uL)
  vR <- qnorm(uR)

  F_mid <- stats::approxfun(ux, cw, method = "linear", ties = "ordered", rule = 2)
  Q_mid <- stats::approxfun(cw, ux, method = "linear", ties = "ordered", rule = 2)

  dx <- pmax(diff(ux), 1e-12)
  seg_slope <- diff(cw) / dx
  vgrid <- qnorm(cw)
  seg_slope_v <- diff(vgrid) / dx
  g_scale <- stats::median(abs(seg_slope_v[is.finite(seg_slope_v)]))
  if (!is.finite(g_scale) || g_scale <= 0) g_scale <- 1
  gL <- min(max(seg_slope_v[1L], 1e-6), 1e12 * g_scale)
  gR <- min(max(seg_slope_v[length(seg_slope_v)], 1e-6), 1e12 * g_scale)

  slope_mid <- function(x) {
    k <- pmin(pmax(findInterval(x, ux, all.inside = TRUE), 1L), length(seg_slope))
    seg_slope[k]
  }

  Fhat <- function(x) {
    x <- as.numeric(x)
    lt <- x < xmin
    rt <- x > xmax
    md <- !(lt | rt)
    out <- numeric(length(x))
    if (any(md)) out[md] <- F_mid(x[md])
    if (any(lt)) out[lt] <- pnorm(vL + gL * (x[lt] - xmin))
    if (any(rt)) out[rt] <- pnorm(vR + gR * (x[rt] - xmax))
    pmin(pmax(out, 1e-12), 1 - 1e-12)
  }

  slope <- function(x) {
    x <- as.numeric(x)
    lt <- x < xmin
    rt <- x > xmax
    md <- !(lt | rt)
    out <- numeric(length(x))
    if (any(md)) out[md] <- slope_mid(x[md])
    if (any(lt)) {
      v <- vL + gL * (x[lt] - xmin)
      out[lt] <- dnorm(v) * gL
    }
    if (any(rt)) {
      v <- vR + gR * (x[rt] - xmax)
      out[rt] <- dnorm(v) * gR
    }
    pmax(out, 1e-300)
  }

  Qhat <- function(u) {
    u <- pmin(pmax(u, 1e-12), 1 - 1e-12)
    lt <- u < uL
    rt <- u > uR
    md <- !(lt | rt)
    out <- numeric(length(u))
    if (any(md)) out[md] <- Q_mid(u[md])
    if (any(lt)) out[lt] <- xmin + (qnorm(u[lt]) - vL) / gL
    if (any(rt)) out[rt] <- xmax + (qnorm(u[rt]) - vR) / gR
    out
  }

  list(
    Fhat = Fhat, Qhat = Qhat, slope = slope,
    F_mid = F_mid, Q_mid = Q_mid,
    xmin = xmin, xmax = xmax, vL = vL, vR = vR, gL = gL, gR = gR
  )
}

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
  stats::approx(cw, x, xout = p, ties = "ordered", rule = 2)$y
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

.normal_bundle <- function(mu, sig) {
  local({
    mu0 <- as.numeric(mu)
    sig0 <- as.numeric(sig)
    list(
      Fhat = function(x) pnorm((x - mu0) / sig0),
      Qhat = function(u) mu0 + sig0 * qnorm(pmin(pmax(u, 1e-12), 1 - 1e-12)),
      slope = function(x) dnorm((x - mu0) / sig0) / sig0,
      F_mid = function(x) pnorm((x - mu0) / sig0),
      Q_mid = function(u) mu0 + sig0 * qnorm(pmin(pmax(u, 1e-12), 1 - 1e-12)),
      xmin = -Inf,
      xmax = Inf,
      vL = -8,
      vR = 8,
      gL = 1 / sig0,
      gR = 1 / sig0
    )
  })
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
  cols <- list(Zsel, Zsel^2 - 1)
  if (ncol(Zsel) >= 2L && ncol(Zsel) <= interaction_max) {
    cmb <- utils::combn(seq_len(ncol(Zsel)), 2L)
    for (k in seq_len(ncol(cmb))) {
      cols[[length(cols) + 1L]] <- Zsel[, cmb[1L, k]] * Zsel[, cmb[2L, k]]
    }
  }

  Xraw <- do.call(cbind, cols)
  if (is.null(centers)) centers <- colMeans(Xraw)
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
  if (!length(keep)) keep <- ord[seq_len(min(1L, dprev))]

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
    return(.normal_bundle(mu, sig))
  }
  rq <- .weighted_quantile_linear(r, w, pgrid)
  .make_tail_safe_cdf_bundle(rq, pgrid)
}

# Sparse triangular transport:
# 1. Tail-safe marginal Gaussianization x -> y.
# 2. Sequential sparse autoregression y_j ~ m_j(z_<j) with residual Gaussianization.
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

  marginals <- vector("list", d)
  for (j in seq_len(d)) {
    xj <- X[, j]
    xq <- .weighted_quantile_linear(xj, w, pgrid)
    rng <- xq[length(xq)] - xq[1L]
    if (!is.finite(rng) || rng < 1e-12) {
      muj <- .weighted_mean1(xj, w)
      sigj <- .weighted_sd1(xj, w, mu = muj, floor = 1e-6)
      marginals[[j]] <- .normal_bundle(muj, sigj)
    } else {
      eps <- max(1e-12, 1e-8 * max(1, abs(rng)))
      for (k in 2:length(xq)) {
        if (xq[k] <= xq[k - 1L]) xq[k] <- xq[k - 1L] + eps
      }
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

compute_transport_diag <- function(Tmap, Theta, reference_prior, nsample = 2000L) {
  n <- nrow(Theta)
  if (n > nsample) Theta <- Theta[sample.int(n, nsample), , drop = FALSE]
  Z <- Tmap$fwd(Theta)
  lp <- reference_prior_logpdf(reference_prior, Theta)
  lqZ <- dmvnorm_chol_log(Z, rep(0, ncol(Z)), chol(diag(ncol(Z))))
  r <- lp - (lqZ + Tmap$log_jac(Theta))

  S <- stats::cov(Z)
  dS <- diag(S)
  diag_dev <- max(abs(dS - 1))
  off_max <- if (ncol(Z) > 1) max(abs(S[upper.tri(S)])) else 0

  list(
    ok = is.finite(mean(r)) && is.finite(stats::sd(r)),
    mean = mean(r),
    sd = stats::sd(r),
    diag_dev = diag_dev,
    off_max = off_max
  )
}

shrinkage_from_lambda <- function(lambda) {
  val <- 0.25 + 0.70 * (lambda^0.8)
  pmin(pmax(val, 0.20), 0.95)
}

build_transport <- function(Theta,
                            w,
                            mu_ref = NULL,
                            Sigma_ref = NULL,
                            reference_prior = NULL,
                            lambda = 0,
                            ngrid = 400,
                            tail = 1e-3,
                            verbose = TRUE) {
  reference_prior <- normalize_reference_prior(reference_prior = reference_prior, mu = mu_ref, Sigma = Sigma_ref)
  Tmap <- fit_copula_transform(
    Theta,
    w,
    ngrid = ngrid,
    tail = tail,
    corr_shrink = shrinkage_from_lambda(lambda)
  )
  diag <- compute_transport_diag(Tmap, Theta, reference_prior = reference_prior)
  if (isTRUE(verbose)) {
    cat(sprintf(
      "  Transport check: mean=%.4f sd=%.4f | diag_dev=%.3f off_max=%.3f\n",
      diag$mean, diag$sd, diag$diag_dev, diag$off_max
    ))
  }
  attr(Tmap, "diag") <- diag
  Tmap
}

.wmean_cov <- function(Z, w) {
  w <- pmax(w, 0)
  sw <- sum(w)
  w <- if (!is.finite(sw) || sw <= 0) rep(1 / nrow(Z), nrow(Z)) else w / sw
  mu <- colSums(Z * w)
  Zc <- sweep(Z, 2L, mu, "-")
  S <- t(Zc) %*% (Zc * w)
  list(mu = mu, S = as.matrix(S))
}

.kl_to_standard_normal <- function(mu, S) {
  d <- length(mu)
  S <- as.matrix(Matrix::nearPD((S + t(S)) / 2, conv.tol = 1e-7)$mat)
  ld <- tryCatch(determinant(S, logarithm = TRUE)$modulus[[1]], error = function(e) -Inf)
  tr <- sum(diag(S))
  0.5 * (tr - ld - d + sum(mu * mu))
}

maybe_refit_transport <- function(Theta,
                                  w,
                                  lambda,
                                  round,
                                  resampled,
                                  ess_frac,
                                  Tmap,
                                  tr_state,
                                  every_rounds = 2L,
                                  kl_improve_frac = 0.06,
                                  verbose = TRUE) {
  vcat <- function(...) if (isTRUE(verbose)) base::cat(...)
  if (is.null(tr_state)) tr_state <- list(last_refit_round = 0L, last_kl = NA_real_)

  time_based <- (round - tr_state$last_refit_round) >= every_rounds
  degeneracy <- isTRUE(resampled) || ess_frac < 0.35
  if (!(time_based || degeneracy)) {
    return(list(refit = FALSE, Tmap = Tmap, tr_state = tr_state))
  }

  rho <- max(0.30, 1 - lambda)
  w_blend <- rho * rep(1 / nrow(Theta), nrow(Theta)) + (1 - rho) * (w / sum(w))
  shrink <- shrinkage_from_lambda(lambda)
  Tcand <- tryCatch(
    fit_copula_transform(Theta, w_blend, corr_shrink = shrink),
    error = function(e) NULL
  )
  if (is.null(Tcand)) {
    return(list(refit = FALSE, Tmap = Tmap, tr_state = tr_state))
  }

  Z_cur <- Tmap$fwd(Theta)
  wc <- .wmean_cov(Z_cur, w)
  KLc <- .kl_to_standard_normal(wc$mu, wc$S)

  Z_can <- Tcand$fwd(Theta)
  wn <- .wmean_cov(Z_can, w)
  KLn <- .kl_to_standard_normal(wn$mu, wn$S)

  accept <- is.finite(KLn) && (!is.finite(KLc) || KLn <= (1 - kl_improve_frac) * KLc)
  if (accept) {
    tr_state$last_refit_round <- round
    tr_state$last_kl <- KLn
    vcat(sprintf("  Transport refit: KL %.3f -> %.3f (shrink=%.2f) [accepted]\n", KLc, KLn, shrink))
    return(list(refit = TRUE, Tmap = Tcand, tr_state = tr_state))
  }

  tr_state$last_kl <- KLc
  vcat(sprintf("  Transport refit: KL %.3f -> %.3f (shrink=%.2f) [rejected]\n", KLc, KLn, shrink))
  list(refit = FALSE, Tmap = Tmap, tr_state = tr_state)
}

.cheap_sort_order <- function(Z) {
  Z <- as.matrix(Z)
  n <- nrow(Z)
  if (n <= 1L) return(seq_len(n))
  if (ncol(Z) <= 1L) return(order(Z[, 1L]))
  order(rowSums(Z))
}

.std_normal_logpdf_Z <- function(Z) {
  rowSums(stats::dnorm(as.matrix(Z), log = TRUE))
}

.normalize_checkpoint_lambdas <- function(checkpoint_lambdas, lambda_target) {
  cp <- sort(unique(as.numeric(checkpoint_lambdas %||% numeric(0))))
  cp <- cp[is.finite(cp) & cp > 0 & cp <= as.numeric(lambda_target)]
  if (!length(cp)) return(numeric(0))
  if (abs(tail(cp, 1L) - as.numeric(lambda_target)) > 1e-12) {
    cp <- c(cp, as.numeric(lambda_target))
  }
  cp
}

.evaluate_z_state <- function(Z, Tmap, data, loglik_fn, reference_prior, n_cores = 1L) {
  Z <- as.matrix(Z)
  Theta <- Tmap$inv(Z)
  colnames(Theta) <- reference_prior$param_names
  ok <- rowSums(!is.finite(Theta)) == 0
  loglik <- rep(-Inf, nrow(Z))
  lpz <- rep(-Inf, nrow(Z))
  if (any(ok)) {
    loglik[ok] <- ll_parallel(Theta[ok, , drop = FALSE], data, loglik_fn, n_cores = n_cores)
    lpz[ok] <- reference_prior_logpdf(reference_prior, Theta[ok, , drop = FALSE]) -
      Tmap$log_jac(Theta[ok, , drop = FALSE])
  }
  loglik[!is.finite(loglik)] <- -Inf
  lpz[!is.finite(lpz)] <- -Inf
  list(Theta = Theta, loglik = loglik, lpz = lpz)
}

mcmc_moves_z_mix_batched <- function(Z,
                                     loglik,
                                     lpz,
                                     Tmap,
                                     lambda,
                                     data,
                                     loglik_fn,
                                     reference_prior,
                                     elite_mix,
                                     batch_size = 4096L,
                                     rw_prob = 0.35,
                                     pcn_prob = 0.25,
                                     rw_scale = 0.8,
                                     pcn_beta = 0.35,
                                     indep_t_df = 7,
                                     n_cores = 1L,
                                     seed = NULL,
                                     allow_pcn = TRUE) {
  if (!is.null(seed)) set.seed(as.integer(seed))
  Z <- as.matrix(Z)
  N <- nrow(Z)
  d <- ncol(Z)
  batch_size <- as.integer(max(1L, batch_size))
  use_indep <- !is.null(elite_mix) && !.is_empty_mix(elite_mix)

  probs <- c(rw = max(rw_prob, 0), pcn = if (isTRUE(allow_pcn)) max(pcn_prob, 0) else 0, indep = if (use_indep) 1 else 0)
  if (sum(probs) <= 0) probs["rw"] <- 1
  probs <- probs / sum(probs)
  kernel <- sample(names(probs), size = N, replace = TRUE, prob = probs)

  accepted <- rep(FALSE, N)
  n_prop <- 0L
  n_accept <- 0L
  by_kernel <- setNames(rep(0L, length(probs)), names(probs))
  acc_by_kernel <- by_kernel

  current_target <- lpz + lambda * loglik
  current_target_pcn <- current_target - .std_normal_logpdf_Z(Z)

  run_block <- function(idx, type) {
    if (!length(idx)) return(NULL)
    k <- length(idx)
    Zc <- Z[idx, , drop = FALSE]
    if (identical(type, "rw")) {
      Zp <- Zc + (rw_scale / sqrt(max(d, 1))) * matrix(rnorm(k * d), nrow = k, ncol = d)
      log_q_ratio <- rep(0, k)
      use_pcn_ratio <- FALSE
    } else if (identical(type, "pcn")) {
      beta <- pmin(pmax(as.numeric(pcn_beta), 1e-6), 0.999)
      Zp <- sqrt(1 - beta^2) * Zc + beta * matrix(rnorm(k * d), nrow = k, ncol = d)
      log_q_ratio <- rep(0, k)
      use_pcn_ratio <- TRUE
    } else {
      Zp <- rmvt_mixture_Z_qmc(
        k,
        elite_mix$meansZ,
        elite_mix$cache,
        nu = indep_t_df,
        seed = sample.int(.Machine$integer.max, 1L)
      )
      log_q_c <- dmvt_mixture_logpdf_Z_vec(Zc, elite_mix$meansZ, elite_mix$cache, indep_t_df)
      log_q_p <- dmvt_mixture_logpdf_Z_vec(Zp, elite_mix$meansZ, elite_mix$cache, indep_t_df)
      log_q_ratio <- log_q_c - log_q_p
      use_pcn_ratio <- FALSE
    }
    colnames(Zp) <- colnames(Z)

    prop <- .evaluate_z_state(
      Zp,
      Tmap = Tmap,
      data = data,
      loglik_fn = loglik_fn,
      reference_prior = reference_prior,
      n_cores = n_cores
    )
    prop_target <- prop$lpz + lambda * prop$loglik
    if (isTRUE(use_pcn_ratio)) {
      log_alpha <- (prop_target - .std_normal_logpdf_Z(Zp)) - current_target_pcn[idx]
    } else {
      log_alpha <- prop_target - current_target[idx] + log_q_ratio
    }
    accept <- is.finite(log_alpha) & (log(runif(k)) < pmin(0, log_alpha))
    list(idx = idx, accept = accept, Zp = Zp, prop = prop)
  }

  for (type in names(probs)) {
    idx_all <- which(kernel == type)
    by_kernel[type] <- length(idx_all)
    if (!length(idx_all)) next
    for (start in seq.int(1L, length(idx_all), by = batch_size)) {
      idx <- idx_all[seq.int(start, min(start + batch_size - 1L, length(idx_all)))]
      res <- run_block(idx, type)
      if (is.null(res)) next
      take <- res$accept
      n_prop <- n_prop + length(idx)
      n_accept <- n_accept + sum(take)
      acc_by_kernel[type] <- acc_by_kernel[type] + sum(take)
      if (any(take)) {
        rows <- idx[take]
        Z[rows, ] <- res$Zp[take, , drop = FALSE]
        loglik[rows] <- res$prop$loglik[take]
        lpz[rows] <- res$prop$lpz[take]
        accepted[rows] <- TRUE
      }
    }
  }

  state <- .evaluate_z_state(
    Z,
    Tmap = Tmap,
    data = data,
    loglik_fn = loglik_fn,
    reference_prior = reference_prior,
    n_cores = n_cores
  )

  list(
    Z = Z,
    Theta = state$Theta,
    loglik = loglik,
    lpz = lpz,
    acc_rate = if (n_prop > 0L) n_accept / n_prop else 0,
    accepted = accepted,
    kernel_counts = by_kernel,
    kernel_accept = acc_by_kernel
  )
}

enhanced_smc_elite <- function(data,
                               loglik_fn,
                               mu_ref = NULL,
                               Sigma_ref = NULL,
                               reference_prior = NULL,
                               M = 4000L,
                               resample_threshold = 0.50,
                               n_mcmc_moves = 3L,
                               max_rounds = 100L,
                               lambda_target = 1.0,
                               G_mix = 12L,
                               elite_quantile = 0.4,
                               gamma_sharp = 1.0,
                               refit_every = 2L,
                               rw_prob = 0.35,
                               rw_scale_init = 0.8,
                               pcn_prob = 0.25,
                               pcn_beta = 0.35,
                               indep_t_df = 7,
                               batch_size = 4096L,
                               n_cores = 1L,
                               seed = 123L,
                               verbose = TRUE,
                               checkpoint_lambdas = NULL,
                               deterministic_resampling = FALSE,
                               resample_sort_mode = c("cheap1d", "none")) {
  reference_prior <- normalize_reference_prior(reference_prior = reference_prior, mu = mu_ref, Sigma = Sigma_ref)
  M <- as.integer(M)
  if (M <= 1L) stop("M must be greater than one.")
  lambda_target <- as.numeric(lambda_target)
  if (!is.finite(lambda_target) || lambda_target <= 0 || lambda_target > 1) {
    stop("lambda_target must be in (0, 1].")
  }
  resample_sort_mode <- match.arg(resample_sort_mode)
  checkpoint_lambdas <- .normalize_checkpoint_lambdas(checkpoint_lambdas, lambda_target = lambda_target)
  checkpoint_log_evidence <- rep(NA_real_, length(checkpoint_lambdas))
  checkpoint_ptr <- 1L
  vcat <- function(...) if (isTRUE(verbose)) base::cat(...)

  set.seed(as.integer(seed))
  Theta <- reference_prior_sample(reference_prior, M)
  loglik <- ll_parallel(Theta, data, loglik_fn, n_cores = n_cores)
  loglik[!is.finite(loglik)] <- -Inf
  w <- rep(1 / M, M)

  Tmap <- build_transport(
    Theta,
    w,
    reference_prior = reference_prior,
    lambda = 0,
    verbose = verbose
  )
  Z <- Tmap$fwd(Theta)
  colnames(Z) <- colnames(Theta)
  lpz <- reference_prior_logpdf(reference_prior, Theta) - Tmap$log_jac(Theta)

  lambda <- 0.0
  round <- 0L
  log_evidence <- 0.0
  mcse_var_accum <- 0.0
  rw_scale <- as.numeric(rw_scale_init)
  elite_mix <- .default_std_normal_mix(ncol(Z))
  tr_state <- list(last_refit_round = 0L, last_kl = NA_real_)

  lambda_hist <- lambda
  ess_hist <- 1.0
  accept_hist <- numeric(0)
  rw_scale_hist <- rw_scale
  resampled_hist <- logical(0)
  move_kernel_hist <- list()

  while (lambda < lambda_target - 1e-12 && round < as.integer(max_rounds)) {
    round <- round + 1L
    next_checkpoint <- if (checkpoint_ptr <= length(checkpoint_lambdas)) {
      checkpoint_lambdas[checkpoint_ptr]
    } else {
      lambda_target
    }
    cess_target <- cess_target_at_lambda(lambda)
    lambda_new <- next_lambda_via_rCESS(
      w = w,
      loglik = loglik,
      lambda = lambda,
      target = cess_target,
      lambda_target = min(lambda_target, next_checkpoint)
    )
    if (lambda_new <= lambda) {
      lambda_new <- min(lambda_target, lambda + min(1e-4, lambda_target - lambda))
    }

    delta <- lambda_new - lambda
    x <- as.numeric(loglik)
    ok <- is.finite(x)
    if (!any(ok)) stop("All local log-likelihood values are non-finite.")
    x[!ok] <- min(x[ok])
    mx <- max(x)
    u <- exp(delta * (x - mx))
    mu1 <- sum(w * u)
    mu2 <- sum(w * u * u)
    neff <- 1 / sum(w * w)
    mcse_var_accum <- mcse_var_accum +
      max((mu2 - mu1^2) / (max(neff, 1) * max(mu1^2, .Machine$double.eps)), 0)

    logw_raw <- log(pmax(w, .Machine$double.eps)) + delta * (x - mx)
    lse <- logsumexp(logw_raw)
    log_evidence <- log_evidence + delta * mx + lse
    w <- exp(logw_raw - lse)
    w <- w / sum(w)
    lambda <- lambda_new
    ess_frac <- ESS(w) / length(w)
    resampled <- FALSE

    vcat(sprintf(
      "\nRound %d: lambda %.3f | ESS=%.3f | logZ=%.4f\n",
      round, lambda, ess_frac, log_evidence
    ))

    if (ess_frac < resample_threshold) {
      ord <- if (identical(resample_sort_mode, "cheap1d")) .cheap_sort_order(Z) else seq_len(M)
      idx_sorted <- stratified_resample_sorted(w[ord], deterministic = deterministic_resampling)
      idx <- ord[idx_sorted]
      Theta <- Theta[idx, , drop = FALSE]
      Z <- Z[idx, , drop = FALSE]
      loglik <- loglik[idx]
      lpz <- lpz[idx]
      w <- rep(1 / M, M)
      resampled <- TRUE
      vcat("  Resampled\n")
    }

    refit <- maybe_refit_transport(
      Theta = Theta,
      w = w,
      lambda = lambda,
      round = round,
      resampled = resampled,
      ess_frac = ess_frac,
      Tmap = Tmap,
      tr_state = tr_state,
      every_rounds = refit_every,
      verbose = verbose
    )
    Tmap <- refit$Tmap
    tr_state <- refit$tr_state
    if (isTRUE(refit$refit)) {
      Z <- Tmap$fwd(Theta)
      colnames(Z) <- colnames(Theta)
      lpz <- reference_prior_logpdf(reference_prior, Theta) - Tmap$log_jac(Theta)
    }

    elite_mix <- tryCatch(
      fit_elite_mixture_Z(
        Z,
        w^gamma_sharp,
        elite_quantile = elite_quantile,
        G = as.integer(G_mix),
        warm_start_mixZ = elite_mix,
        lambda = lambda,
        verbose = FALSE
      ),
      error = function(e) .default_std_normal_mix(ncol(Z))
    )

    move_accept_round <- numeric(0)
    kernel_counts_round <- NULL
    kernel_accept_round <- NULL
    if (as.integer(n_mcmc_moves) > 0L) {
      for (move_id in seq_len(as.integer(n_mcmc_moves))) {
        move <- mcmc_moves_z_mix_batched(
          Z = Z,
          loglik = loglik,
          lpz = lpz,
          Tmap = Tmap,
          lambda = lambda,
          data = data,
          loglik_fn = loglik_fn,
          reference_prior = reference_prior,
          elite_mix = elite_mix,
          batch_size = batch_size,
          rw_prob = rw_prob,
          pcn_prob = pcn_prob,
          rw_scale = rw_scale,
          pcn_beta = pcn_beta,
          indep_t_df = indep_t_df,
          n_cores = n_cores,
          seed = seed + 1009L * round + move_id
        )
        Z <- move$Z
        Theta <- move$Theta
        loglik <- move$loglik
        lpz <- move$lpz
        move_accept_round <- c(move_accept_round, move$acc_rate)
        kernel_counts_round <- (kernel_counts_round %||% setNames(rep(0L, length(move$kernel_counts)), names(move$kernel_counts))) +
          move$kernel_counts
        kernel_accept_round <- (kernel_accept_round %||% setNames(rep(0L, length(move$kernel_accept)), names(move$kernel_accept))) +
          move$kernel_accept
      }

      acc_mean <- mean(move_accept_round)
      rw_scale <- exp(.clamp(log(rw_scale) + 0.20 * (acc_mean - 0.234), log(0.05), log(4.0)))
      vcat(sprintf("  MCMC accept=%.3f | rw_scale=%.3f\n", acc_mean, rw_scale))
    } else {
      acc_mean <- NA_real_
    }

    lambda_hist <- c(lambda_hist, lambda)
    ess_hist <- c(ess_hist, ess_frac)
    accept_hist <- c(accept_hist, acc_mean)
    rw_scale_hist <- c(rw_scale_hist, rw_scale)
    resampled_hist <- c(resampled_hist, resampled)
    move_kernel_hist[[round]] <- list(counts = kernel_counts_round, accept = kernel_accept_round)

    while (checkpoint_ptr <= length(checkpoint_lambdas) &&
           lambda >= checkpoint_lambdas[checkpoint_ptr] - 1e-12) {
      checkpoint_log_evidence[checkpoint_ptr] <- log_evidence
      checkpoint_ptr <- checkpoint_ptr + 1L
    }
  }

  if (lambda < lambda_target - 1e-12) {
    warning("enhanced_smc_elite hit max_rounds before reaching lambda_target.")
  }

  mcse_logZ <- sqrt(mcse_var_accum)
  if (isTRUE(verbose)) {
    cat(sprintf(
      "\nFinished local SMC: rounds=%d | lambda=%.3f/%.3f | logZ=%.4f +/- %.4f\n",
      round, lambda, lambda_target, log_evidence, mcse_logZ
    ))
  }

  list(
    Theta = Theta,
    Z = Z,
    w = w,
    loglik = loglik,
    lpz = lpz,
    final_lambda = lambda,
    log_evidence = log_evidence,
    mcse_logZ = mcse_logZ,
    checkpoint_lambdas = checkpoint_lambdas,
    checkpoint_log_evidence = checkpoint_log_evidence,
    transport = Tmap,
    reference_prior = reference_prior,
    elite_mix_final = elite_mix,
    meta = list(
      rounds = round,
      ess = ESS(w),
      lambda_hist = lambda_hist,
      ess_hist = ess_hist,
      accept_hist = accept_hist,
      rw_scale_hist = rw_scale_hist,
      resampled_hist = resampled_hist,
      move_kernel_hist = move_kernel_hist,
      transport = "sparse_triangular"
    )
  )
}
