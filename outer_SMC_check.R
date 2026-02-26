# --------- normalize + resampling (multinomial/stratified/systematic) ----------
normalize_w <- function(w) w / sum(w)

resample_indices <- function(w, n, method = c("systematic","stratified","multinomial")) {
  w <- normalize_w(w); method <- match.arg(method)
  cw <- cumsum(w); N <- length(w)
  if (method == "multinomial") {
    return(sample.int(N, size = n, replace = TRUE, prob = w))
  }
  if (method == "stratified") {
    u <- (runif(n) + (0:(n-1))) / n
  } else { # systematic
    u0 <- runif(1, 0, 1/n); u <- u0 + (0:(n-1))/n
  }
  findInterval(u, cw) + 1L
}

# ------------------- φ posterior draws from an SMC result ----------------------
# out$phi: M x dphi, out$w: length M
posterior_draws_phi <- function(out, n = 5000L, method = "systematic") {
  w <- normalize_w(out$w)
  idx <- resample_indices(w, n, method)
  out$phi[idx, , drop = FALSE]
}

# --------------- Weighted expectation / quantiles without resampling ----------
post_expectation_phi <- function(out, f) {
  # f: function(phi_row) -> scalar; computed as self-normalized IS
  w <- normalize_w(out$w); phi <- out$phi
  vals <- apply(phi, 1L, f)
  sum(w * vals)
}

wquantile <- function(x, w, probs = c(0.025, 0.5, 0.975)) {
  o <- order(x); x <- x[o]; w <- normalize_w(w[o])
  cw <- cumsum(w)
  sapply(probs, function(p) {
    i <- which(cw >= p)[1]
    if (is.na(i)) return(tail(x, 1L))
    if (i == 1L) return(x[1L])
    # linear interpolation
    x1 <- x[i-1]; x2 <- x[i]; w1 <- cw[i-1]; w2 <- cw[i]
    if (w2 == w1) x2 else x1 + (p - w1) * (x2 - x1) / (w2 - w1)
  })
}

wq_colwise <- function(X, w, probs = c(0.025, 0.5, 0.975)) {
  res <- vapply(seq_len(ncol(X)), function(j) wquantile(X[, j], w, probs),
                numeric(length(probs)))
  rownames(res) <- paste0("q", c("2.5","50","97.5"))
  colnames(res) <- colnames(X) %||% paste0("dim", seq_len(ncol(X)))
  res
}
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ---------------- θ | φ draws from a subject cache (using your cache) -----------
# Works with either the Gaussian fast path (recommended; cache$gaussian_map must be set)
# or a generic log_prior_theta_given_phi_mat(Theta, phi, aux).
.sample_theta_given_phi <- function(cache, phi,
                                    log_prior_theta_given_phi_mat = NULL) {
  # compute log-weights: log p(y|θ) + log p(θ|φ) - log r(θ)
  if (!is.null(cache$gaussian_map) && is.null(log_prior_theta_given_phi_mat)) {
    pars <- cache$gaussian_map(phi)
    lprior <- .log_prior_gauss_mat(cache$Theta, pars)  # from your cache code
  } else {
    stopifnot(!is.null(log_prior_theta_given_phi_mat))
    lprior <- log_prior_theta_given_phi_mat(cache$Theta, phi, aux = NULL)
  }
  lw <- cache$log_py + lprior - cache$log_r
  # sample index with probs ∝ exp(lw)
  m <- max(lw); p <- exp(lw - m); p <- p / sum(p)
  sample.int(cache$M, size = 1L, replace = TRUE, prob = p)
}

# Draw K samples of θ_i | φ for one subject cache
rtheta_given_phi <- function(cache, phi, K = 1L, log_prior_theta_given_phi_mat = NULL) {
  idx <- replicate(K, .sample_theta_given_phi(cache, phi, log_prior_theta_given_phi_mat))
  cache$Theta[idx, , drop = FALSE]
}

# ----------------- Full joint posterior draws (φ, θ_{1:S}) --------------------
# WARNING: This can be large; consider drawing in chunks or only the θ's you need.
posterior_draws_joint <- function(out, caches, n = 1000L,
                                  method = "systematic",
                                  log_prior_theta_given_phi_mat = NULL) {
  phi_samp <- posterior_draws_phi(out, n, method)
  S <- length(caches)
  thetas <- vector("list", S)
  for (i in seq_len(S)) {
    # for each φ draw, sample one θ_i | φ
    idx_i <- vapply(seq_len(n), function(j)
      .sample_theta_given_phi(caches[[i]], phi_samp[j, ],
                              log_prior_theta_given_phi_mat),
      integer(1))
    thetas[[i]] <- caches[[i]]$Theta[idx_i, , drop = FALSE]
  }
  list(phi = phi_samp, theta_list = thetas)
}
