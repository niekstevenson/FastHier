# ----------------------- log-densities (reused) -----------------------
dnorm_log <- function(x, m, v) { -0.5 * (log(2*pi*v) + (x - m)^2 / v) }
dinvgamma_log <- function(x, a, b) {  # shape a, rate b
  ifelse(x > 0, a*log(b) - lgamma(a) - (a+1)*log(x) - b/x, -Inf)
}

# ------------------ Dimension-specific prior over φ -------------------
# φ = (μ_1,...,μ_d,  ℓ_1,...,ℓ_d) with ℓ_j = log σ²_j
# Inputs m0, s0, a, b can be scalars or vectors of length d.
make_prior_phi_diag <- function(m0, s0, a, b, d) {
  m0 <- rep_len(as.numeric(m0), d)
  s0 <- rep_len(as.numeric(s0), d)
  a  <- rep_len(as.numeric(a),  d)
  b  <- rep_len(as.numeric(b),  d)

  rprior_phi <- function(n) {
    # sample μ_j ~ N(m0_j, s0_j) independently over j
    mu <- matrix(
      rnorm(n * d, mean = rep(m0, each = n), sd = sqrt(rep(s0, each = n))),
      nrow = n, ncol = d, byrow = FALSE
    )
    # sample σ²_j ~ InvGamma(a_j, b_j) independently over j
    sigma2 <- matrix(
      1 / rgamma(n * d, shape = rep(a, each = n), rate = rep(b, each = n)),
      nrow = n, ncol = d, byrow = FALSE
    )
    out <- cbind(mu, log_sigma2 = log(sigma2))
    colnames(out) <- c(paste0("mu", seq_len(d)), paste0("log_sigma2_", seq_len(d)))
    out
  }

  logprior_phi <- function(phi_row) {
    mu  <- as.numeric(phi_row[seq_len(d)])
    ell <- as.numeric(phi_row[d + seq_len(d)])
    sigma2 <- exp(ell)
    sum(dnorm_log(mu, m0, s0)) +
      sum(dinvgamma_log(sigma2, a, b)) +
      sum(ell)  # Jacobian for σ² = e^ℓ
  }

  list(rprior = rprior_phi, lprior = logprior_phi, d = d)
}

# --------- θ|φ Gaussian map (fast path), diagonal & per-dim -----------
# Uses the first d entries of μ and the next d entries of log σ².
phi_to_gaussian_params_diag_factory <- function() {
  function(phi, d) {
    stopifnot(length(phi) >= 2*d)
    mu  <- as.numeric(phi[seq_len(d)])
    ell <- as.numeric(phi[d + seq_len(d)])
    sigma2 <- exp(ell)
    list(
      mu = mu,
      Sigma_inv = diag(1/sigma2, d),
      logdet = sum(log(sigma2)),   # log |diag(σ²)| = Σ_j log σ²_j
      const = NULL
    )
  }
}

# Weighted quantile for one vector
wquantile <- function(x, w, probs = c(0.025, 0.5, 0.975)) {
  o <- order(x); x <- x[o]; w <- w[o]; w <- w / sum(w)
  cw <- cumsum(w)
  sapply(probs, function(p) {
    i <- which(cw >= p)[1]
    if (is.na(i)) return(tail(x, 1L))
    if (i == 1L) return(x[1L])
    # linear interpolation between the two surrounding points
    w1 <- cw[i-1]; w2 <- cw[i]; x1 <- x[i-1]; x2 <- x[i]
    if (w2 == w1) x2 else x1 + (p - w1) * (x2 - x1) / (w2 - w1)
  })
}

# Column-wise weighted quantiles (returns a matrix: rows=probs, cols=dimensions)
wq_colwise <- function(X, w, probs = c(0.025, 0.5, 0.975)) {
  res <- vapply(seq_len(ncol(X)), function(j) wquantile(X[, j], w, probs),
                numeric(length(probs)))
  rownames(res) <- paste0(c("q", "q", "q"), c("2.5", "50", "97.5"))
  colnames(res) <- colnames(X) %||% paste0("dim", seq_len(ncol(X)))
  res
}
