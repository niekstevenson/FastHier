# phi_surrogate_setup.R
# φ-aware surrogate posterior setup using base-proxy anchors (Patch 4) — DIAGONAL prior
# Assumes (from your environment):
#   - build_subject_cache_from_smc(...)
#   - log_marginal_proxy_from_base(cache, phi, log_prior_theta_given_phi_mat, ...)
#   - rprior_phi(n) and logprior_phi(phi_row)
# Inputs you already have per subject when calling the builders below:
#   - smc_list, data_list, subjects, loglik_fn

suppressPackageStartupMessages({
  library(Matrix)
  library(stats)
})

sc_path <- file.path(getwd(), "smc_core.R")
if (file.exists(sc_path)) source(sc_path)

source("new_SMC_cache.R")

# If your cache builders live here; adjust/remove if already sourced elsewhere.
# source("new_SMC_cache.R")

# -----------------------------------------------------------------------------
# Utilities
# -----------------------------------------------------------------------------

.spd_fix <- function(S, floor_eig = 1e-8) {
  S <- as.matrix((S + t(S)) / 2)
  ev <- eigen(S, symmetric = TRUE)
  lam <- pmax(ev$values, floor_eig)
  ev$vectors %*% diag(lam, nrow = length(lam)) %*% t(ev$vectors)
}

# -----------------------------------------------------------------------------
# φ → Gaussian prior map (Option: DIAGONAL precision)
# -----------------------------------------------------------------------------
# φ packs: [μ (length dθ), then η (length dθ)] with diagonal precisions τ = exp(η).
# Returns both Sigma_inv (as a diagonal matrix) and tau (vector) for fast paths.
make_gaussian_map_fn_diagprec <- function(d_theta) {
  force(d_theta)
  function(phi, d = d_theta) {
    mu  <- as.numeric(phi[seq_len(d)])
    eta <- as.numeric(phi[d + seq_len(d)])
    tau <- exp(eta)                                 # precisions (diagonal)
    list(
      mu = mu,
      Sigma_inv = diag(tau, d),                     # keep for compatibility
      tau = tau,                                    # fast path
      logdet = -sum(log(tau)),                      # log|Σ| = -log|K| = -sum log τ
      const = 0
    )
  }
}

# -----------------------------------------------------------------------------
# Base-proxy wrapper: Θ (M×d) × φ → log p(Θ|φ) (length M)  (DIAGONAL fast path)
# -----------------------------------------------------------------------------
make_log_prior_theta_given_phi_mat_gauss <- function(gaussian_map_fn) {
  function(Theta, phi, aux = NULL) {
    pars <- gaussian_map_fn(phi, ncol(Theta))
    mu   <- pars$mu
    if (!is.null(pars$tau)) {
      # diagonal fast path
      lt <- sweep(Theta, 2L, mu, `-`)
      q  <- lt^2 %*% pars$tau
      -0.5 * (as.numeric(q) + pars$logdet + ncol(Theta) * log(2*pi))
    } else {
      # generic (shouldn't hit for diag map, but keep as fallback)
      K  <- as.matrix(pars$Sigma_inv)
      lt <- sweep(Theta, 2L, mu, `-`)
      q  <- rowSums((lt %*% K) * lt)
      -0.5 * (q + pars$logdet + ncol(Theta) * log(2*pi))
    }
  }
}

# -----------------------------------------------------------------------------
# Global φ-anchors from base proxies (no new likelihood calls)
# -----------------------------------------------------------------------------
make_phi_anchors_from_base_proxy <- function(
    caches,                 # list of subject caches (length S)
    rprior_phi,             # function(n) -> matrix n x dφ
    logprior_phi,           # function(phi_row) -> scalar
    gaussian_map_fn,        # diag map created above
    K = 5L,
    Npilot = 4000L,
    Nselect = 2000L,
    seed = 1L,
    report = TRUE
) {
  set.seed(seed)
  Phi <- as.matrix(rprior_phi(Npilot))         # Npilot x dφ
  N   <- nrow(Phi); S <- length(caches)
  lp  <- apply(Phi, 1L, logprior_phi)

  lp_theta_fn <- make_log_prior_theta_given_phi_mat_gauss(gaussian_map_fn)

  ll_proxy <- numeric(N)
  for (j in 1:N) {
    phi_j <- Phi[j, , drop = TRUE]
    s <- 0
    for (i in seq_len(S)) {
      s <- s + log_marginal_proxy_from_base(caches[[i]], phi_j, lp_theta_fn)
    }
    ll_proxy[j] <- s
  }

  lw <- lp + ll_proxy
  w  <- exp(lw - max(lw)); w <- w / sum(w)

  Nselect <- min(Nselect, N)
  idx <- sample.int(N, size = Nselect, replace = TRUE, prob = w)
  PhiW <- Phi[idx, , drop = FALSE]

  km   <- kmeans(PhiW, centers = K, iter.max = 100)
  centers <- matrix(km$centers, nrow = K, ncol = ncol(Phi),
                    dimnames = list(NULL, colnames(Phi)))
  w_cent <- as.numeric(km$size); w_cent <- w_cent / sum(w_cent)
  anchors <- list(centers = centers, weights = w_cent)
  if (report) {
    cat(sprintf("[Anchors] base-proxy: K=%d, Npilot=%d, Nselect=%d | weighted by cluster size\n", K, Npilot, Nselect))
  }
  anchors
}

# -----------------------------------------------------------------------------
# Per-subject quadratic *likelihood* site via BLOCK evaluation (no row-wise calls)
# -----------------------------------------------------------------------------
# Diagonal site: log p(y_i|θ) ≈ a + b^T θ - 1/2 ∑_j c_j θ_j^2   (no cross terms)
# We fit (a, b, c_diag) by ridge least squares on a *batch* of design points,
# evaluated in one call to your C++ vectorized loglik_fn(Θ_block, data_i).
fit_quadratic_site_diag_matrix <- function(
    data_i, loglik_fn, gaussian_map_fn,
    phi_star, d_theta,
    n_design = NULL,           # defaults to max(6d, 2d+20)
    radius   = 1.6,            # design spread: multiples of prior sd at φ*
    ridge    = 1e-3,           # ridge on quadratic terms
    c_min    = 1e-6,           # floor for curvature to keep SPD
    seed     = 1L,
    par_names= NULL
) {
  set.seed(seed)
  pars <- gaussian_map_fn(phi_star, d_theta)
  mu0  <- as.numeric(pars$mu)
  tau  <- if (!is.null(pars$tau)) pars$tau else diag(pars$Sigma_inv)
  sd0  <- 1 / sqrt(as.numeric(tau))

  if (is.null(n_design)) n_design <- max(6L * d_theta, 2L * d_theta + 20L)

  # Design: Gaussian cloud around mu0 with axis-wise scales ~ prior sd
  Z <- matrix(rnorm(n_design * d_theta), n_design, d_theta)
  Theta <- sweep(Z, 2L, radius * sd0, `*`)
  Theta <- sweep(Theta, 2L, mu0, `+`)
  colnames(Theta) <- par_names
  # One BLOCK likelihood call (C++): returns vector of length n_design
  y <- ll_parallel(Theta, data_i, loglik_fn)

  # Diagonal quadratic regression: y ≈ a + b^T θ - 1/2 ∑ c_j θ_j^2
  X <- cbind(1, Theta, -0.5 * (Theta^2))
  p  <- ncol(X)
  pen <- diag(c(rep(1e-6, 1 + d_theta), rep(ridge, d_theta)), p, p)
  XtX <- crossprod(X) + pen
  Xty <- crossprod(X, y)
  beta <- as.numeric(solve(XtX, Xty))

  a <- beta[1L]
  b <- beta[1 + seq_len(d_theta)]
  c <- pmax(beta[1 + d_theta + seq_len(d_theta)], c_min)

  J <- diag(c, d_theta)  # site curvature (diagonal)
  h <- b                 # site linear term

  list(J = J, h = h, theta_hat = mu0, converged = TRUE, a = a)
}

# -----------------------------------------------------------------------------
# Create a φ-aware cache for one subject using DIAGONAL prior map
# -----------------------------------------------------------------------------
create_phi_aware_cache <- function(
    smc_out_i,
    data_i,
    subj_id_i,
    loglik_fn,
    gaussian_map_fn,    # from make_gaussian_map_fn_diagprec
    phi_anchors,
    quad_site,          # list(J, h) from fit_quadratic_site_diag_matrix
    M = 512L,
    K_batches = 4L,
    surr_n = 200L,
    surr_c = 1.6,
    surr_weight = 0.35,
    blend_std_norm = 0.15,
    defensive_t_eps = 0.25,
    defensive_t_df = 3L
) {
  build_subject_cache_from_smc(
    smc_out          = smc_out_i,
    data             = data_i,
    loglik_fn        = loglik_fn,
    M                = M,
    K_batches        = K_batches,
    subj_id          = subj_id_i,
    # φ-aware surrogate inputs (Patch 4):
    phi_anchors      = phi_anchors,
    quad_site        = quad_site,
    gaussian_map_fn  = function(phi, d) gaussian_map_fn(phi, d),
    # surrogate mixture tuning
    surr_n           = surr_n,
    surr_c           = surr_c,
    surr_weight      = surr_weight,
    # defensive arms as in your setup
    blend_std_norm   = blend_std_norm,
    defensive_t_eps  = defensive_t_eps,
    defensive_t_df   = defensive_t_df
  )
}

# -----------------------------------------------------------------------------
# Multi-subject orchestrator:
#   1) build SMALL caches (cheap) for all subjects (anchors stage),
#   2) compute global φ-anchors from base proxies,
#   3) fit per-subject diagonal sites via BLOCK likelihood,
#   4) build final φ-aware caches for all subjects.
# -----------------------------------------------------------------------------
build_phi_aware_caches_from_base_proxy <- function(
    smc_list,            # list of smc_out (length S)
    data_list,           # list of data_i (length S)
    subjects,            # vector of subject ids (length S)
    loglik_fn,           # function(Θ_block, data_i) -> vector log p(y_i|θ)
    rprior_phi,          # function(n) -> matrix n x dφ
    logprior_phi,        # function(phi_row) -> scalar log p(φ)
    # anchor & sampling knobs
    K_anchors = 5L,
    Npilot = 4000L,
    Nselect = 2000L,
    seed = 1L,
    # “small” cache (for anchors) and “final” cache sizes
    M_small = 64L,
    K_small = 2L,
    M_final = 512L,
    K_final = 4L,
    # φ-aware surrogate tuning
    surr_n = 200L,
    surr_c = 1.6,
    surr_weight = 0.35,
    # defensive proposal knobs
    blend_std_norm = 0.15,
    defensive_t_eps = 0.25,
    defensive_t_df = 3L,
    # parallelism
    mc.cores = 1L,       # subjects-level parallelism
    n_cores_inner = 1L   # inner likelihood parallelism per subject (used by cache builder)
) {
  stopifnot(length(smc_list) == length(data_list),
            length(data_list) == length(subjects),
            length(smc_list) >= 1L)
  S <- length(smc_list)
  d_theta <- ncol(smc_list[[1]]$Theta)

  gaussian_map_fn <- make_gaussian_map_fn_diagprec(d_theta)

  cat(sprintf("[Stage 1] Building SMALL caches for %d subjects (M=%d, K=%d)...\n",
              S, M_small, K_small))
  caches_small <- parallel::mcmapply(
    smc_list, data_list, subjects,
    FUN = function(smc_out_i, data_i, subj_id_i) {
      build_subject_cache_from_smc(
        smc_out          = smc_out_i,
        data             = data_i,
        loglik_fn        = loglik_fn,
        M                = M_small,
        K_batches        = K_small,
        subj_id          = subj_id_i,
        blend_std_norm   = blend_std_norm,
        defensive_t_eps  = defensive_t_eps,
        defensive_t_df   = defensive_t_df,
        n_cores          = n_cores_inner
      )
    },
    mc.cores = mc.cores,
    SIMPLIFY = FALSE
  )

  cat("[Stage 2] Building global φ-anchors from base proxies...\n")
  phi_anchors <- make_phi_anchors_from_base_proxy(
    caches          = caches_small,
    rprior_phi      = rprior_phi,
    logprior_phi    = logprior_phi,
    gaussian_map_fn = gaussian_map_fn,
    K               = K_anchors,
    Npilot          = Npilot,
    Nselect         = Nselect,
    seed            = seed,
    report          = TRUE
  )
  # weighted mean center for site fitting
  phi_star <- drop(crossprod(phi_anchors$weights, phi_anchors$centers))

  cat("[Stage 3] Fitting per-subject diagonal Laplace sites...\n")
  quad_sites <- parallel::mcmapply(
    data_list, smc_list,
    FUN = function(data_i, smc_out_i) {
      fit_quadratic_site_diag_matrix(
        data_i, loglik_fn, gaussian_map_fn,
        phi_star = phi_star, d_theta = d_theta,
        n_design = max(6L * d_theta, 2L * d_theta + 20L),
        radius   = 1.6, ridge = 1e-3, c_min = 1e-6, seed = seed,
        par_names = colnames(smc_out_i$Theta)
      )
    },
    mc.cores = mc.cores,
    SIMPLIFY = FALSE
  )

  cat("[Stage 4] Building φ-aware caches for all subjects...\n")
  caches_phi_aware <- parallel::mcmapply(
    smc_list, data_list, subjects, quad_sites,
    FUN = function(smc_out_i, data_i, subj_id_i, quad_site_i) {
      create_phi_aware_cache(
        smc_out_i       = smc_out_i,
        data_i          = data_i,
        subj_id_i       = subj_id_i,
        loglik_fn       = loglik_fn,
        gaussian_map_fn = gaussian_map_fn,
        phi_anchors     = phi_anchors,
        quad_site       = quad_site_i,
        M               = M_final,
        K_batches       = K_final,
        surr_n          = surr_n,
        surr_c          = surr_c,
        surr_weight     = surr_weight,
        blend_std_norm  = blend_std_norm,
        defensive_t_eps = defensive_t_eps,
        defensive_t_df  = defensive_t_df
      )
    },
    mc.cores = mc.cores,
    SIMPLIFY = FALSE
  )

  cat("[Done] φ-aware caches built for all subjects.\n")
  return(list(
    caches_small      = caches_small,
    phi_anchors       = phi_anchors,
    gaussian_map_fn   = gaussian_map_fn,
    quad_sites        = quad_sites,
    caches_phi_aware  = caches_phi_aware
  ))
}
