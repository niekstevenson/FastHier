# ========================================================================
# Subject-level caches (Option B: K RQMC batches, batch-index randomization)
# - Proposal r_i(θ): transported GMM in Z with defensive t-arm + std-normal blend
# - K independent RQMC batches (fixed once), batch index = auxiliary U_i
# - Unbiased IS estimator per batch; outer SMC/MCMC uses batch indices for CRNs
# - Stable log-sum-exp; finite checks; nearPD + jitter; parallel likelihood
# - DA proxies from base SMC; PSIS-k diagnostics
# ========================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(mvtnorm)
  library(qrng)
})

# Consolidate shared helpers from smc_core
sc_path <- file.path(getwd(), "smc_core.R")
if (file.exists(sc_path)) source(sc_path)


# === PATCH: DM-MIS helpers =================================================

# Largest-remainder integer rounding to match an exact total.
.fracs_to_counts <- function(target, total) {
  target <- pmax(target, 0)
  base   <- floor(target)
  rem    <- as.integer(round(total - sum(base)))
  if (rem > 0L) {
    frac <- target - base
    ord  <- order(frac, decreasing = TRUE)
    base[ord[seq_len(min(rem, length(base)))]] <- base[ord[seq_len(min(rem, length(base)))] ] + 1L
  } else if (rem < 0L) {
    # remove from smallest fractional parts first (but keep >= 0)
    frac <- target - base
    ord  <- order(frac, decreasing = FALSE)
    j <- 1L
    while (rem < 0L && j <= length(base)) {
      if (base[ord[j]] > 0L) { base[ord[j]] <- base[ord[j]] - 1L; rem <- rem + 1L }
      j <- j + 1L
    }
  }
  as.integer(base)
}

# Decompose n into a sum of powers-of-two (binary decomposition).
.split_to_pow2_counts <- function(n) {
  n <- as.integer(n)
  if (n <= 0L) return(integer(0))
  out <- integer(0)
  rem <- n
  while (rem > 0L) {
    p <- 2L ^ floor(log2(rem))
    out <- c(out, p)
    rem <- rem - p
  }
  out
}


## duplicated numerics/mixture helpers removed; using smc_core versions

.ensure_mix_cache <- function(mix) {
  if (is.null(mix$cache)) {
    mix$cache <- prep_mix_cache(mix$meansZ, mix$covsZ, mix$wZ)
  }
  mix
}

# --- PATCH (4): combine two mixtures with external weights ---------------
.combine_two_mixes <- function(mixA, mixB, wA = 0.5, wB = 0.5) {
  if (.is_empty_mix(mixA)) return(mixB)
  if (.is_empty_mix(mixB)) return(mixA)
  wA <- max(0, wA); wB <- max(0, wB)
  if ((wA + wB) <= 0) { wA <- 0.5; wB <- 0.5 }
  meansZ <- c(mixA$meansZ, mixB$meansZ)
  covsZ  <- c(mixA$covsZ,  mixB$covsZ)
  wZ     <- c(rep(1, length(mixA$meansZ)) * wA / max(length(mixA$meansZ), 1L),
              rep(1, length(mixB$meansZ)) * wB / max(length(mixB$meansZ), 1L))
  wZ     <- wZ / sum(wZ)
  list(meansZ = meansZ, covsZ = covsZ, wZ = wZ,
       cache = prep_mix_cache(meansZ, covsZ, wZ))
}

# --- bridge component from p(theta|phi_anchor), mapped through transport ----
.build_bridge_mixZ <- function(Tmap, phi_anchor, gaussian_map_fn,
                               n_bridge = 300L, cov_inflation = 1.2) {
  if (is.null(phi_anchor) || is.null(gaussian_map_fn)) return(NULL)
  pars <- tryCatch(gaussian_map_fn(phi_anchor), error = function(e) NULL)
  if (is.null(pars) || is.null(pars$mu) || is.null(pars$Sigma_inv)) return(NULL)
  K <- as.matrix(pars$Sigma_inv)
  L <- tryCatch(chol(K), error = function(e) NULL)
  if (is.null(L)) return(NULL)
  Sig <- chol2inv(L)
  Sig <- as.matrix((Sig + t(Sig)) / 2)
  Sig <- Sig * (cov_inflation^2)
  Sig <- as.matrix(Matrix::nearPD(Sig, conv.tol = 1e-7)$mat)
  Theta <- mvtnorm::rmvnorm(as.integer(max(16L, n_bridge)), mean = as.numeric(pars$mu), sigma = Sig)
  Z <- Tmap$fwd(Theta)
  muZ <- colMeans(Z)
  Zc <- sweep(Z, 2L, muZ, `-`)
  SigZ <- (t(Zc) %*% Zc) / max(1L, nrow(Zc))
  SigZ <- as.matrix((SigZ + t(SigZ)) / 2)
  SigZ <- as.matrix(Matrix::nearPD(SigZ, conv.tol = 1e-7)$mat)
  mix <- list(meansZ = list(as.numeric(muZ)), covsZ = list(SigZ), wZ = 1)
  mix$cache <- prep_mix_cache(mix$meansZ, mix$covsZ, mix$wZ)
  mix
}

.as_phi_anchor_matrix <- function(phi_anchors) {
  if (is.null(phi_anchors)) return(NULL)
  if (is.matrix(phi_anchors)) return(phi_anchors)
  if (is.vector(phi_anchors)) return(matrix(phi_anchors, nrow = 1L))
  if (is.list(phi_anchors) && !is.null(phi_anchors$centers)) {
    return(as.matrix(phi_anchors$centers))
  }
  NULL
}

# Build a multi-anchor bridge mixture from a panel of phi anchors.
.build_bridge_panel_mixZ <- function(Tmap, phi_anchors, gaussian_map_fn,
                                     n_bridge_total = 400L, cov_inflation = 1.2) {
  centers <- .as_phi_anchor_matrix(phi_anchors)
  if (is.null(centers) || !nrow(centers) || is.null(gaussian_map_fn)) return(NULL)

  n_anchor <- nrow(centers)
  n_each <- as.integer(max(16L, ceiling(n_bridge_total / n_anchor)))
  meansZ <- list()
  covsZ <- list()

  for (k in seq_len(n_anchor)) {
    pars <- tryCatch(gaussian_map_fn(centers[k, , drop = TRUE]), error = function(e) NULL)
    if (is.null(pars) || is.null(pars$mu) || is.null(pars$Sigma_inv)) next

    K <- as.matrix(pars$Sigma_inv)
    L <- tryCatch(chol(K), error = function(e) NULL)
    if (is.null(L)) next

    Sig <- chol2inv(L)
    Sig <- as.matrix((Sig + t(Sig)) / 2)
    Sig <- Sig * (cov_inflation^2)
    Sig <- as.matrix(Matrix::nearPD(Sig, conv.tol = 1e-7)$mat)

    Theta <- mvtnorm::rmvnorm(n_each, mean = as.numeric(pars$mu), sigma = Sig)
    Z <- Tmap$fwd(Theta)
    if (!is.matrix(Z) || nrow(Z) <= 1L) next

    muZ <- colMeans(Z)
    Zc <- sweep(Z, 2L, muZ, `-`)
    SigZ <- (t(Zc) %*% Zc) / max(1L, nrow(Zc))
    SigZ <- as.matrix((SigZ + t(SigZ)) / 2)
    SigZ <- as.matrix(Matrix::nearPD(SigZ, conv.tol = 1e-7)$mat)

    meansZ[[length(meansZ) + 1L]] <- as.numeric(muZ)
    covsZ[[length(covsZ) + 1L]] <- SigZ
  }

  if (!length(meansZ)) return(NULL)
  wZ <- rep(1 / length(meansZ), length(meansZ))
  list(
    meansZ = meansZ,
    covsZ = covsZ,
    wZ = wZ,
    cache = prep_mix_cache(meansZ, covsZ, wZ)
  )
}

# --- PATCH (4): φ-aware surrogate mixture via quadratic site -------------
.phi_aware_surrogate_mixZ <- function(Tmap, phi_anchors, quad_site,
                                      gaussian_map_fn,
                                      surr_n = 200L, surr_c = 1.6, surr_weight = 0.35) {
  if (is.null(phi_anchors) || is.null(quad_site) || is.null(gaussian_map_fn)) return(NULL)
  stopifnot(!is.null(quad_site$J), !is.null(quad_site$h))
  J <- as.matrix(quad_site$J)
  h <- as.numeric(quad_site$h)
  # Support matrix or list(centers, weights)
  if (is.matrix(phi_anchors)) {
    centers <- phi_anchors
    w_anchor <- rep(1 / nrow(centers), nrow(centers))
  } else if (is.list(phi_anchors) && !is.null(phi_anchors$centers)) {
    centers <- as.matrix(phi_anchors$centers)
    w_anchor <- phi_anchors$weights
    if (is.null(w_anchor)) w_anchor <- rep(1 / nrow(centers), nrow(centers))
    w_anchor <- pmax(w_anchor, 0); w_anchor <- w_anchor / sum(w_anchor)
  } else {
    return(NULL)
  }
  comps_mu <- list(); comps_S <- list()
  for (k in seq_len(nrow(centers))) {
    pars <- gaussian_map_fn(centers[k, , drop = TRUE])
    K    <- as.matrix(pars$Sigma_inv)
    mu0  <- as.numeric(pars$mu)
    # surrogate posterior in θ
    Sig_t <- tryCatch(solve(J + K), error = function(e) MASS::ginv(J + K))
    mu_t  <- Sig_t %*% (h + K %*% mu0)
    # inflate and sample, then push-forward through T
    Sig_i <- as.matrix((Sig_t + t(Sig_t))/2)
    Sig_i <- surr_c^2 * Sig_i
    Theta_s <- mvtnorm::rmvnorm(surr_n, mean = as.numeric(mu_t), sigma = Sig_i)
    Z_s     <- Tmap$fwd(Theta_s)
    muZ <- colMeans(Z_s)
    Zc  <- sweep(Z_s, 2L, muZ, `-`)
    SigZ <- (t(Zc) %*% Zc) / nrow(Zc)
    comps_mu[[length(comps_mu)+1L]] <- as.numeric(muZ)
    comps_S[[length(comps_S)+1L]]   <- as.matrix((SigZ + t(SigZ))/2)
  }
  wZ <- if (exists("w_anchor")) rep(w_anchor, length.out = length(comps_mu)) else rep(1/length(comps_mu), length(comps_mu))
  list(meansZ = comps_mu, covsZ = comps_S, wZ = wZ,
       cache = prep_mix_cache(comps_mu, comps_S, wZ),
       surr_weight = surr_weight)
}

## inflate_mixture_along_dims: use smc_core::inflate_mixture_along_dims

# ---------------------- build r_i(θ) from SMC outputs --------------------
# r(θ) = q_Z(Z; mix_g ⊕ ε_t t) * |J(θ)|, with Z = T(θ)
.build_ri_from_smc <- function(Z, w, Tmap,
                               ref_mix = NULL, use_ref_mix_first = TRUE,
                               G = 12L, elite_quantile = 0.40, cov_inflation = 4.0,
                               lambda = 1.0,
                               blend_std_norm = 0.15,
                               defensive_t_eps = 0.15, defensive_t_df = 3L,
                               # --- WS1/WS2 knobs:
                               cache_bridge_enable = FALSE,
                               cache_bridge_weight = 0.20,
                               cache_defensive_weight = 0.00,
                               cache_defensive_df = 3L,
                               phi_anchor = NULL,
                               # --- PATCH (4) new args:
                               phi_anchors = NULL,
                               quad_site   = NULL,  # list(J=..., h=...)
                               gaussian_map_fn = NULL,
                               surr_n = 200L, surr_c = 1.6, surr_weight = 0.35) {
  d <- ncol(Z)
  mix_g <- NULL

  if (use_ref_mix_first && !is.null(ref_mix) && !.is_empty_mix(ref_mix)) {
    mix_g <- ref_mix
  } else if (exists("fit_elite_mixture_Z", mode = "function")) {
    w_fit <- pmax(w, 0); w_fit <- w_fit / sum(w_fit)
    G_eff <- max(2L, min(G, 2L + floor(6 * lambda)))
    merge_thresh_eff <- if (lambda >= 0.85) 0.06 else if (lambda >= 0.70) 0.08 else 0.10
    mix_g <- tryCatch(
      fit_elite_mixture_Z(Z, w_fit,
                          elite_quantile = max(0.20, elite_quantile),
                          G = G_eff,
                          cov_inflation = cov_inflation,
                          housekeeping = TRUE, min_G_keep = 2,
                          merge_thresh = merge_thresh_eff, verbose = FALSE,
                          lambda = lambda),
      error = function(e) NULL
    )
  }

  if (is.null(mix_g) || .is_empty_mix(mix_g)) {
    w_fit <- pmax(w, 0); w_fit <- w_fit / sum(w_fit)
    mu <- colSums(Z * w_fit)
    Zc <- sweep(Z, 2L, mu, `-`)
    Sig <- t(Zc) %*% (Zc * w_fit)
    Sig <- (Sig + t(Sig)) / 2
    Sig <- as.matrix(Matrix::nearPD(Sig, conv.tol = 1e-7)$mat)
    Sig <- Sig * cov_inflation
    mix_g <- list(meansZ = list(as.numeric(mu)), covsZ = list(Sig), wZ = 1)
    mix_g$cache <- prep_mix_cache(mix_g$meansZ, mix_g$covsZ, mix_g$wZ)
  }

  # --- PATCH (4): φ-aware surrogate components (prepend before base/defense)
  if (!is.null(phi_anchors) && !is.null(quad_site) && !is.null(gaussian_map_fn)) {
    mix_surr <- .phi_aware_surrogate_mixZ(Tmap, phi_anchors, quad_site,
                                          gaussian_map_fn, surr_n = surr_n,
                                          surr_c = surr_c, surr_weight = surr_weight)
    if (!is.null(mix_surr) && !.is_empty_mix(mix_surr)) {
      # combine with target proportion in favor of surrogate components
      mix_g <- .combine_two_mixes(mix_surr, mix_g, wA = mix_surr$surr_weight, wB = 1 - mix_surr$surr_weight)
    }
  }

  bridge_mix <- NULL
  if (isTRUE(cache_bridge_enable) && !is.null(gaussian_map_fn)) {
    if (!is.null(phi_anchors)) {
      bridge_mix <- .build_bridge_panel_mixZ(
        Tmap = Tmap,
        phi_anchors = phi_anchors,
        gaussian_map_fn = gaussian_map_fn,
        n_bridge_total = 400L,
        cov_inflation = 1.2
      )
    }
    if (is.null(bridge_mix) && !is.null(phi_anchor)) {
      bridge_mix <- .build_bridge_mixZ(
        Tmap = Tmap,
        phi_anchor = phi_anchor,
        gaussian_map_fn = gaussian_map_fn,
        n_bridge = 300L,
        cov_inflation = 1.2
      )
    }
    if (!is.null(bridge_mix) && !.is_empty_mix(bridge_mix)) {
      w_bridge <- .clamp(cache_bridge_weight, 0, 0.95)
      # reserve mass for optional defensive arm
      w_local <- max(1e-6, 1 - w_bridge - .clamp(cache_defensive_weight, 0, 0.90))
      mix_g <- .combine_two_mixes(mix_g, bridge_mix, wA = w_local, wB = w_bridge)
    }
  }

  std_mix <- .default_std_normal_mix(d)
  mix_g <- blend_mixes(mix_g, std_mix, eps = blend_std_norm)
  mix_g <- .ensure_mix_cache(mix_g)

  # WS2: optional defensive t-arm centered by bridge if available.
  t_eps <- if (cache_defensive_weight > 0) .clamp(cache_defensive_weight, 0, 0.95) else .clamp(defensive_t_eps, 0, 0.95)
  t_df <- if (cache_defensive_weight > 0) as.integer(max(2L, cache_defensive_df)) else as.integer(max(2L, defensive_t_df))
  mix_t <- NULL
  if (t_eps > 0) {
    mix_t <- if (!is.null(bridge_mix) && !.is_empty_mix(bridge_mix)) bridge_mix else mix_g
  }

  list(
    mix_g = mix_g,
    mix_t = mix_t,
    t_eps = t_eps,
    t_df  = t_df,
    Tmap  = Tmap
  )
}

# ------------------------ log r(θ) under mixture -------------------------
.log_r_theta_mix <- function(Theta, Z, ri) {
  lq_n <- gmm_logpdf_Z_vec(Z, ri$mix_g$meansZ, ri$mix_g$cache)
  if (!is.null(ri$mix_t) && ri$t_eps > 0) {
    lq_t <- dmvt_mixture_logpdf_Z_vec(Z, ri$mix_t$meansZ, ri$mix_t$cache, nu = ri$t_df)
    lq   <- rlogsumexp2(lq_n + log1p(-ri$t_eps), lq_t + log(ri$t_eps))
  } else {
    lq <- lq_n
  }
  lq + ri$Tmap$log_jac(Theta)
}

# ------------------------ unified RQMC sampler (arm+comp) ----------------
# Single Sobol net drives: arm (Normal vs t), component, chi^2 (for t), and d normals.
# === PATCH: deterministic-mixture + stratified + power-of-two RQMC ========
.sample_from_ri_rqmc <- function(M, ri, sobol_seed = NULL, deterministic_counts = FALSE) {
  d <- length(ri$mix_g$meansZ[[1]])
  if (M <= 0L) return(list(Theta = matrix(numeric(0), 0, d), Z = matrix(numeric(0), 0, d)))

  # mixture arm weights
  use_t <- (!is.null(ri$mix_t) && ri$t_eps > 0)
  p_t   <- if (use_t) ri$t_eps else 0
  p_n   <- 1 - p_t

  # component weights
  w_comp <- exp(ri$mix_g$cache$logw); w_comp <- w_comp / sum(w_comp)
  G <- length(w_comp)

  # seeds per stratum (deterministic but subject/batch specific via sobol_seed)
  base_seed <- if (!is.null(sobol_seed)) as.integer(sobol_seed) else 11L

  # Counts for (normal arm, components) and (t arm, components).
  # Default is randomized (strict correctness); deterministic is debug-only.
  target_counts <- c(p_n * w_comp, p_t * w_comp) * M
  probs <- target_counts / sum(target_counts)
  counts <- if (isTRUE(deterministic_counts)) {
    .fracs_to_counts(target_counts, total = M)
  } else {
    set.seed(base_seed + 104729L)
    as.integer(rmultinom(1L, size = M, prob = probs)[, 1L])
  }
  n_ng <- counts[seq_len(G)]
  n_tg <- counts[G + seq_len(G)]

  chunk_id <- 0L

  Z_blocks <- vector("list", 2L * G)  # upper bound; we may store multiple chunks per stratum
  out_ptr  <- 0L

  # Coerce any qrng::sobol() output to an n x d matrix
  .force_mat <- function(U, n, d) {
    if (is.matrix(U) && nrow(U) == n && ncol(U) == d) return(U)
    matrix(as.numeric(U), nrow = n, ncol = d, byrow = FALSE)
  }

  # Clamp to (0,1) to avoid qnorm/ log issues
  .clamp01 <- function(U) pmin(pmax(U, 1e-12), 1 - 1e-12)

  # Generate N(0,1) nets as an n x d matrix
  .sobol_norm <- function(n, d, seed) {
    U <- qrng::sobol(n, d = d, randomize = "digital.shift", seed = seed)
    U <- .force_mat(U, n, d)
    qnorm(.clamp01(U))
  }

  # Generate U(0,1) as a length-n vector (robust when d=1 or n=1)
  .sobol_unif1 <- function(n, seed) {
    U <- qrng::sobol(n, d = 1, randomize = "digital.shift", seed = seed)
    U <- .force_mat(U, n, 1L)
    as.numeric(.clamp01(U[, 1, drop = TRUE]))
  }

  for (g in seq_len(G)) {
    L <- ri$mix_g$cache$Ls[[g]]
    mu <- ri$mix_g$meansZ[[g]]

    # ---- Normal arm, component g
    if (n_ng[g] > 0L) {
      chunks <- .split_to_pow2_counts(n_ng[g])  # each is a power of two
      for (n_chunk in chunks) {
        chunk_id <- chunk_id + 1L
        N   <- .sobol_norm(n_chunk, d, seed = base_seed + 1299709L * chunk_id)
        Zg  <- N %*% t(L)
        Zg  <- sweep(Zg, 2L, mu, `+`)
        out_ptr <- out_ptr + 1L
        Z_blocks[[out_ptr]] <- Zg
      }
    }

    # ---- t arm, component g
    if (use_t && n_tg[g] > 0L) {
      chunks <- .split_to_pow2_counts(n_tg[g])
      for (n_chunk in chunks) {
        chunk_id <- chunk_id + 1L
        N    <- .sobol_norm(n_chunk, d, seed = base_seed + 1299709L * chunk_id)
        Uchi <- .sobol_unif1(n_chunk, seed = base_seed + 1299709L * chunk_id + 7L)
        Gchi <- qchisq(Uchi, df = ri$t_df)
        scale <- sqrt(ri$t_df / Gchi)                     # length n_chunk
        Zg  <- N %*% t(L)
        Zg  <- sweep(Zg, 1L, scale, `*`)                  # row-wise scale
        Zg  <- sweep(Zg, 2L, mu, `+`)
        out_ptr <- out_ptr + 1L
        Z_blocks[[out_ptr]] <- Zg
      }
    }
  }

  # bind and map back to θ
  if (out_ptr == 0L) {
    Z <- matrix(numeric(0), 0, d)
  } else {
    Z <- do.call(rbind, Z_blocks[seq_len(out_ptr)])
  }
  colnames(Z) <- names(ri$mix_g$meansZ[[1]])
  Theta <- ri$Tmap$inv(Z)
  list(
    Theta = Theta,
    Z = Z,
    strata_counts = as.integer(counts),
    strata_probs = as.numeric(probs)
  )
}

# -------------------- Gaussian prior fast path (matrix interface) --------
.register_gauss_map <- function(phi_to_gaussian_params) phi_to_gaussian_params

.log_prior_gauss_mat <- function(Theta, pars) {
  Theta <- as.matrix(Theta)
  mu <- as.numeric(pars$mu)
  K  <- as.matrix(pars$Sigma_inv)
  lt <- sweep(Theta, 2L, mu, `-`)
  q  <- rowSums((lt %*% K) * lt)
  out <- -0.5 * (q + pars$logdet + ncol(Theta) * log(2*pi))
  if (!is.null(pars$const)) out <- out + pars$const
  as.numeric(out)
}

register_gaussian_prior_map <- function(cache, phi_to_gaussian_params) {
  cache$gaussian_map <- .register_gauss_map(phi_to_gaussian_params)
  cache
}

# ------------------- Base-prior proxy from inner SMC ---------------------
.make_base_cache_from_smc <- function(smc_out) {
  stopifnot(!is.null(smc_out$Theta), !is.null(smc_out$w),
            !is.null(smc_out$lpz), !is.null(smc_out$transport),
            !is.null(smc_out$log_evidence))
  # lpz = log p(Theta) - log|J|  => lpz + log|J| = log p(Theta)
  log_p0_vec <- as.numeric(smc_out$lpz + smc_out$transport$log_jac(smc_out$Theta))
  list(
    theta_base = smc_out$Theta,
    w_bar      = as.numeric(pmax(smc_out$w, 0) / sum(pmax(smc_out$w, 0))),
    log_Z0     = as.numeric(smc_out$log_evidence),
    log_p0_vec = log_p0_vec
  )
}

# ------------------- SUBJECT-LEVEL CACHE (K-batch builder) ---------------
# === PATCH: subject-specific seeding =====================================
build_subject_cache_from_smc <- function(
    smc_out,
    data, subj_id, loglik_fn,
    M = 256L,
    K_batches = 4L,
    use_ref_mix_first = TRUE,
    G = 12L,
    elite_quantile = 0.40,
    cov_inflation = 4.0,
    blend_std_norm = 0.15,
    defensive_t_eps = 0.25,
    defensive_t_df = 3L,
    # WS1/WS2 toggles
    cache_bridge_enable = FALSE,
    cache_bridge_weight = 0.20,
    cache_defensive_weight = 0.00,
    cache_defensive_df = 3L,
    phi_anchor = NULL,
    deterministic_counts = FALSE, # debug mode; default randomized for strict correctness
    sobol_seed = NULL,          # <-- PATCH: allow NULL
    n_cores = 1L,
    # --- PATCH (4) new args for φ-aware surrogate:
    phi_anchors = NULL,
    quad_site = NULL,              # list(J=..., h=...)
    gaussian_map_fn = NULL,
    surr_n = 200L,
    surr_c = 1.6,
    surr_weight = 0.35
) {
  stopifnot(!is.null(smc_out$Z), !is.null(smc_out$w), !is.null(smc_out$transport))
  ri <- .build_ri_from_smc(
    Z = smc_out$Z, w = smc_out$w, Tmap = smc_out$transport,
    ref_mix = smc_out$ref_mix, use_ref_mix_first = use_ref_mix_first,
    G = G, elite_quantile = elite_quantile, cov_inflation = cov_inflation,
    lambda = if (!is.null(smc_out$final_lambda)) smc_out$final_lambda else 1.0,
    blend_std_norm = blend_std_norm,
    defensive_t_eps = defensive_t_eps, defensive_t_df = defensive_t_df,
    cache_bridge_enable = cache_bridge_enable,
    cache_bridge_weight = cache_bridge_weight,
    cache_defensive_weight = cache_defensive_weight,
    cache_defensive_df = cache_defensive_df,
    phi_anchor = phi_anchor,
    phi_anchors = phi_anchors, quad_site = quad_site,
    gaussian_map_fn = if (!is.null(gaussian_map_fn)) function(phi) gaussian_map_fn(phi, ncol(smc_out$Theta)) else NULL,
    surr_n = surr_n, surr_c = surr_c, surr_weight = surr_weight
  )

  # --- PATCH: subject-specific base seed (reproducible, independent across subjects)
  base_seed <- if (!is.null(sobol_seed)) as.integer(sobol_seed) else (11L + 104729L * as.integer(subj_id))

  # K independent RQMC batches
  K_batches <- max(1L, as.integer(K_batches))
  M_per <- as.integer(ceiling(M / K_batches))
  seeds <- base_seed + 7919L * seq_len(K_batches)  # subject-specific stride

  batches <- vector("list", K_batches)
  for (k in seq_len(K_batches)) {
    smp <- .sample_from_ri_rqmc(
      M_per, ri, sobol_seed = seeds[k],
      deterministic_counts = deterministic_counts
    )   # <-- uses stratified DM-MIS sampler
    Theta_k <- smp$Theta; Z_k <- smp$Z
    colnames(Theta_k) <- colnames(smc_out$Theta)
    log_py_k <- ll_parallel(Theta_k, data, loglik_fn, n_cores)
    log_r_k  <- .log_r_theta_mix(Theta_k, Z_k, ri)

    if (!all(is.finite(log_py_k))) stop("Non-finite log p(y|θ) in batch ", k)
    if (!all(is.finite(log_r_k)))  stop("Non-finite log r(θ) in batch ", k)

    batches[[k]] <- list(
      Theta = Theta_k,
      Z = Z_k,
      log_py = as.numeric(log_py_k),
      log_r  = as.numeric(log_r_k),
      M = nrow(Theta_k),
      sobol_seed = seeds[k]
    )
  }

  structure(list(
    batches = batches,
    K = K_batches,
    ri = ri,
    base = .make_base_cache_from_smc(smc_out),
    gaussian_map = NULL,
    meta = list(d = ncol(smc_out$Theta),
                subj_id = subj_id,
                phi_anchors = phi_anchors,
                has_quad_site = !is.null(quad_site),
                cache_bridge_enable = isTRUE(cache_bridge_enable),
                cache_bridge_weight = as.numeric(cache_bridge_weight),
                cache_defensive_weight = as.numeric(cache_defensive_weight),
                cache_defensive_df = as.integer(cache_defensive_df),
                has_phi_anchor = !is.null(phi_anchor),
                deterministic_counts = isTRUE(deterministic_counts))
  ), class = "subject_cache_smcK")
}

# ------------------------ Adaptive PM surrogate stack --------------------

.regularize_cov_safe <- function(S, jitter = 1e-8) {
  S <- as.matrix((S + t(S)) / 2)
  S <- tryCatch(as.matrix(Matrix::nearPD(S, conv.tol = 1e-7)$mat), error = function(e) S)
  d <- nrow(S)
  S + diag(jitter, d)
}

.weighted_mean_vec <- function(X, w) {
  w <- pmax(as.numeric(w), 0)
  w <- w / sum(w)
  colSums(as.matrix(X) * w)
}

.weighted_cov_mat <- function(X, w) {
  w <- pmax(as.numeric(w), 0)
  w <- w / sum(w)
  mu <- .weighted_mean_vec(X, w)
  Xc <- sweep(as.matrix(X), 2L, mu, `-`)
  S <- crossprod(sqrt(w) * Xc)
  .regularize_cov_safe(S)
}

.draw_theta_given_phi_callbacks <- function(phi, n, callbacks, d_theta) {
  if (!is.null(callbacks$rtheta_given_phi)) {
    th <- tryCatch(
      callbacks$rtheta_given_phi(phi, n, aux = NULL),
      error = function(e) callbacks$rtheta_given_phi(phi, n)
    )
    th <- as.matrix(th)
    if (nrow(th) != n) stop("rtheta_given_phi returned incorrect number of draws.")
    return(th)
  }
  if (is.null(callbacks$gaussian_map_fn)) {
    stop("Need either callbacks$rtheta_given_phi or callbacks$gaussian_map_fn for adaptive PM.")
  }
  pars <- callbacks$gaussian_map_fn(phi, d_theta)
  K <- as.matrix(pars$Sigma_inv)
  L <- tryCatch(chol(K), error = function(e) NULL)
  if (is.null(L)) stop("gaussian_map_fn returned non-SPD precision.")
  Sig <- chol2inv(L)
  Sig <- .regularize_cov_safe(Sig)
  mvtnorm::rmvnorm(as.integer(n), mean = as.numeric(pars$mu), sigma = Sig)
}

build_subject_surrogate_from_outer <- function(
    smc_out,
    data = NULL,
    subj_id = NA_integer_,
    loglik_fn = NULL,
    G = 10L,
    elite_quantile = 0.40,
    cov_inflation = 2.5,
    blend_std_norm = 0.05,
    base_seed = NULL
) {
  stopifnot(!is.null(smc_out$Theta), !is.null(smc_out$Z), !is.null(smc_out$w), !is.null(smc_out$transport))
  d_theta <- ncol(smc_out$Theta)
  w <- pmax(as.numeric(smc_out$w), 0)
  w <- w / sum(w)

  mix_local <- tryCatch(
    .build_ri_from_smc(
      Z = smc_out$Z,
      w = w,
      Tmap = smc_out$transport,
      ref_mix = smc_out$ref_mix,
      use_ref_mix_first = TRUE,
      G = G,
      elite_quantile = elite_quantile,
      cov_inflation = cov_inflation,
      lambda = smc_out$final_lambda %||% 1.0,
      blend_std_norm = blend_std_norm,
      defensive_t_eps = 0,
      defensive_t_df = 3L,
      cache_bridge_enable = FALSE,
      cache_bridge_weight = 0,
      cache_defensive_weight = 0,
      cache_defensive_df = 3L,
      phi_anchor = NULL,
      phi_anchors = NULL,
      quad_site = NULL,
      gaussian_map_fn = NULL
    )$mix_g,
    error = function(e) NULL
  )
  if (is.null(mix_local) || .is_empty_mix(mix_local)) {
    muZ <- .weighted_mean_vec(smc_out$Z, w)
    SigZ <- .weighted_cov_mat(smc_out$Z, w)
    mix_local <- list(meansZ = list(as.numeric(muZ)), covsZ = list(SigZ), wZ = 1)
    mix_local$cache <- prep_mix_cache(mix_local$meansZ, mix_local$covsZ, mix_local$wZ)
  } else {
    mix_local <- .ensure_mix_cache(mix_local)
  }

  theta_mean <- as.numeric(.weighted_mean_vec(smc_out$Theta, w))
  theta_cov <- .weighted_cov_mat(smc_out$Theta, w)

  structure(
    list(
      subj_id = as.integer(subj_id),
      transport = smc_out$transport,
      mix_local = mix_local,
      scale_local = theta_cov,
      theta_mean = theta_mean,
      d_theta = d_theta,
      base_seed = as.integer(base_seed %||% (10007L + 1009L * as.integer(subj_id))),
      outer_diag = list(
        rounds = smc_out$meta$rounds %||% length(smc_out$meta$lambda_hist %||% numeric(0)),
        log_evidence = smc_out$log_evidence %||% NA_real_,
        mcse_log_evidence = smc_out$mcse_logZ %||% NA_real_,
        final_lambda = smc_out$final_lambda %||% 1.0
      ),
      data = data,
      loglik_fn = loglik_fn
    ),
    class = "subject_surrogate_pm"
  )
}

.normalize_working_prior_gaussian <- function(working_prior, d_theta = NULL) {
  if (is.null(working_prior)) stop("working_prior must not be NULL.")
  mu <- as.numeric(working_prior$mu %||% working_prior$mean)
  if (!length(mu)) stop("working_prior must contain 'mu' or 'mean'.")
  if (!is.null(d_theta) && length(mu) != d_theta) {
    stop("working_prior mean length does not match d_theta.")
  }

  Sigma <- working_prior$Sigma
  Sigma_inv <- working_prior$Sigma_inv
  if (is.null(Sigma) && is.null(Sigma_inv)) {
    stop("working_prior must contain either 'Sigma' or 'Sigma_inv'.")
  }

  if (is.null(Sigma)) {
    Sigma_inv <- as.matrix(Sigma_inv)
    L <- tryCatch(chol(Sigma_inv), error = function(e) NULL)
    if (is.null(L)) stop("working_prior$Sigma_inv is not SPD.")
    Sigma <- chol2inv(L)
  } else {
    Sigma <- .regularize_cov_safe(Sigma)
  }

  if (is.null(Sigma_inv)) {
    L <- tryCatch(chol(Sigma), error = function(e) NULL)
    if (is.null(L)) stop("working_prior$Sigma is not SPD.")
    Sigma_inv <- chol2inv(L)
  } else {
    Sigma_inv <- as.matrix(Sigma_inv)
  }

  logdet <- as.numeric(working_prior$logdet %||% determinant(Sigma, logarithm = TRUE)$modulus)
  structure(
    list(
      mu = mu,
      Sigma = Sigma,
      Sigma_inv = Sigma_inv,
      logdet = logdet,
      phi_anchor = working_prior$phi_anchor %||% NULL
    ),
    class = "working_prior_gaussian"
  )
}

.draw_from_working_prior_gaussian <- function(n, working_prior) {
  wp <- .normalize_working_prior_gaussian(working_prior)
  mvtnorm::rmvnorm(as.integer(n), mean = wp$mu, sigma = wp$Sigma)
}

log_working_prior_gaussian_mat <- function(Theta, working_prior) {
  wp <- .normalize_working_prior_gaussian(working_prior, d_theta = ncol(as.matrix(Theta)))
  .log_prior_gauss_mat(Theta, wp)
}

build_subject_exact_local_from_outer <- function(
    smc_out,
    data = NULL,
    subj_id = NA_integer_,
    loglik_fn = NULL,
    working_prior = NULL,
    base_seed = NULL,
    ...
) {
  surrogate <- build_subject_surrogate_from_outer(
    smc_out = smc_out,
    data = data,
    subj_id = subj_id,
    loglik_fn = loglik_fn,
    base_seed = base_seed,
    ...
  )
  working_prior_use <- working_prior %||% smc_out$working_prior
  if (is.null(working_prior_use)) {
    stop("build_subject_exact_local_from_outer requires a working_prior or smc_out$working_prior.")
  }
  wp <- .normalize_working_prior_gaussian(working_prior_use, d_theta = surrogate$d_theta)

  structure(
    list(
      subj_id = as.integer(subj_id),
      surrogate = surrogate,
      working_prior = wp,
      data = data %||% surrogate$data,
      loglik_fn = loglik_fn %||% surrogate$loglik_fn,
      base_seed = as.integer(base_seed %||% surrogate$base_seed),
      d_theta = surrogate$d_theta,
      outer_theta = as.matrix(smc_out$Theta),
      outer_w = {
        w <- pmax(as.numeric(smc_out$w), 0)
        sw <- sum(w)
        if (!is.finite(sw) || sw <= 0) rep(1 / nrow(as.matrix(smc_out$Theta)), nrow(as.matrix(smc_out$Theta))) else w / sw
      }
    ),
    class = "subject_exact_local"
  )
}

make_theta_proposal_exact_init <- function(local_obj, init_ctl = list()) {
  stopifnot(inherits(local_obj, "subject_exact_local"))
  ctl <- modifyList(
    list(
      w_local = 0.90,
      w_phi_anchor = 0.00,
      w_defensive = 0.10,
      defensive_df = 3L,
      anchor_n = 32L,
      anchor_cov_inflation = 1.0
    ),
    init_ctl
  )
  callbacks <- list(
    rtheta_given_phi = function(phi, n, aux = NULL) {
      .draw_from_working_prior_gaussian(n, local_obj$working_prior)
    }
  )
  phi_stub <- local_obj$working_prior$phi_anchor %||% local_obj$working_prior$mu
  make_theta_proposal_adaptive(
    phi = phi_stub,
    surrogate = local_obj$surrogate,
    adaptive_pm_control = ctl,
    callbacks = callbacks
  )
}

make_theta_proposal_exact_refresh <- function(phi, local_obj, refresh_ctl = list(), callbacks = list()) {
  stopifnot(inherits(local_obj, "subject_exact_local"))
  ctl <- modifyList(
    list(
      w_local = 0.60,
      w_phi_anchor = 0.30,
      w_defensive = 0.10,
      defensive_df = 3L,
      anchor_n = 96L,
      anchor_cov_inflation = 1.5
    ),
    refresh_ctl
  )
  cb <- modifyList(
    list(
      rtheta_given_phi = function(phi, n, aux = NULL) {
        .draw_from_working_prior_gaussian(n, local_obj$working_prior)
      }
    ),
    callbacks
  )
  make_theta_proposal_adaptive(
    phi = phi,
    surrogate = local_obj$surrogate,
    adaptive_pm_control = ctl,
    callbacks = cb
  )
}

make_theta_proposal_adaptive <- function(phi, surrogate, adaptive_pm_control = list(), callbacks = list()) {
  stopifnot(inherits(surrogate, "subject_surrogate_pm"))
  ctl <- modifyList(
    list(
      w_local = 0.60,
      w_phi_anchor = 0.30,
      w_defensive = 0.10,
      defensive_df = 3L,
      anchor_n = 96L,
      anchor_cov_inflation = 1.5
    ),
    adaptive_pm_control
  )
  w_comp <- c(max(0, ctl$w_local), max(0, ctl$w_phi_anchor), max(0, ctl$w_defensive))
  if (sum(w_comp) <= 0) w_comp <- c(0.6, 0.3, 0.1)
  w_comp <- w_comp / sum(w_comp)

  d <- surrogate$d_theta
  Tmap <- surrogate$transport
  mix_local <- .ensure_mix_cache(surrogate$mix_local)

  mix_anchor <- NULL
  Theta_anchor <- tryCatch(
    .draw_theta_given_phi_callbacks(phi, as.integer(max(16L, ctl$anchor_n)), callbacks, d),
    error = function(e) NULL
  )
  if (!is.null(Theta_anchor) && nrow(Theta_anchor) > 1L) {
    Z_anchor <- Tmap$fwd(Theta_anchor)
    mu_a <- colMeans(Z_anchor)
    Zc <- sweep(Z_anchor, 2L, mu_a, `-`)
    Sig_a <- (t(Zc) %*% Zc) / max(1L, nrow(Zc))
    Sig_a <- .regularize_cov_safe(Sig_a * (as.numeric(ctl$anchor_cov_inflation)^2))
    mix_anchor <- list(meansZ = list(as.numeric(mu_a)), covsZ = list(Sig_a), wZ = 1)
    mix_anchor$cache <- prep_mix_cache(mix_anchor$meansZ, mix_anchor$covsZ, mix_anchor$wZ)
  }
  if (is.null(mix_anchor) || .is_empty_mix(mix_anchor)) {
    mix_anchor <- mix_local
  }
  mix_def <- mix_anchor

  draw <- function(n, seed = NULL) {
    n <- as.integer(n)
    if (n <= 0L) return(matrix(numeric(0), 0, d))
    if (!is.null(seed)) set.seed(as.integer(seed))
    cnt <- as.integer(rmultinom(1L, size = n, prob = w_comp)[, 1L])

    Z_list <- list()
    if (cnt[1L] > 0L) {
      Z_list[[length(Z_list) + 1L]] <- sample_gmm_Z_qmc(cnt[1L], mix_local$meansZ, mix_local$cache, seed = if (is.null(seed)) NULL else seed + 17L)
    }
    if (cnt[2L] > 0L) {
      Z_list[[length(Z_list) + 1L]] <- sample_gmm_Z_qmc(cnt[2L], mix_anchor$meansZ, mix_anchor$cache, seed = if (is.null(seed)) NULL else seed + 31L)
    }
    if (cnt[3L] > 0L) {
      Z_list[[length(Z_list) + 1L]] <- rmvt_mixture_Z_qmc(cnt[3L], mix_def$meansZ, mix_def$cache, nu = as.integer(max(2L, ctl$defensive_df)), seed = if (is.null(seed)) NULL else seed + 47L)
    }
    Z <- do.call(rbind, Z_list)
    if (nrow(Z) != n) stop("Adaptive proposal generated wrong number of draws.")
    if (nrow(Z) > 1L) Z <- Z[sample.int(nrow(Z)), , drop = FALSE]
    Theta <- Tmap$inv(Z)
    as.matrix(Theta)
  }

  log_q <- function(Theta) {
    Theta <- as.matrix(Theta)
    Z <- Tmap$fwd(Theta)
    lq <- rep(-Inf, nrow(Z))
    if (w_comp[1L] > 0) {
      l1 <- log(w_comp[1L]) + gmm_logpdf_Z_vec(Z, mix_local$meansZ, mix_local$cache)
      lq <- l1
    }
    if (w_comp[2L] > 0) {
      l2 <- log(w_comp[2L]) + gmm_logpdf_Z_vec(Z, mix_anchor$meansZ, mix_anchor$cache)
      lq <- if (all(is.infinite(lq))) l2 else rlogsumexp2(lq, l2)
    }
    if (w_comp[3L] > 0) {
      l3 <- log(w_comp[3L]) + dmvt_mixture_logpdf_Z_vec(Z, mix_def$meansZ, mix_def$cache, nu = as.integer(max(2L, ctl$defensive_df)))
      lq <- if (all(is.infinite(lq))) l3 else rlogsumexp2(lq, l3)
    }
    as.numeric(lq + Tmap$log_jac(Theta))
  }

  list(
    draw = draw,
    log_q = log_q,
    w_comp = w_comp,
    mix_local = mix_local,
    mix_anchor = mix_anchor,
    mix_def = mix_def
  )
}

estimate_log_marginal_subject_pm <- function(phi, surrogate, callbacks, aux_state_i = list(), control = list()) {
  stopifnot(inherits(surrogate, "subject_surrogate_pm"))
  ctl <- modifyList(
    list(
      M_default = 128L,
      M_override = NULL,
      max_retries = 2L,
      min_log_mhat = -1e12,
      defensive_df = 3L,
      da_M = 32L
    ),
    control
  )
  M_use <- as.integer(ctl$M_override %||% aux_state_i$M %||% ctl$M_default)
  M_use <- as.integer(max(16L, M_use))
  seed0 <- as.integer(aux_state_i$seed %||% surrogate$base_seed %||% 12345L)

  data_i <- callbacks$data_i %||% surrogate$data
  loglik_fn <- callbacks$loglik_fn %||% surrogate$loglik_fn
  if (is.null(loglik_fn) || is.null(data_i)) {
    stop("estimate_log_marginal_subject_pm requires data_i and loglik_fn.")
  }

  .eval_prior <- function(Theta) {
    if (!is.null(callbacks$log_prior_theta_given_phi_mat)) {
      return(as.numeric(tryCatch(
        callbacks$log_prior_theta_given_phi_mat(Theta, phi, aux = NULL),
        error = function(e) callbacks$log_prior_theta_given_phi_mat(Theta, phi)
      )))
    }
    if (!is.null(callbacks$gaussian_map_fn)) {
      pars <- callbacks$gaussian_map_fn(phi, ncol(Theta))
      return(as.numeric(.log_prior_gauss_mat(Theta, pars)))
    }
    stop("Need callbacks$log_prior_theta_given_phi_mat or callbacks$gaussian_map_fn.")
  }

  last_err <- NULL
  for (attempt in seq_len(as.integer(max(1L, ctl$max_retries + 1L)))) {
    seed_use <- as.integer(seed0 + 1009L * (attempt - 1L))
    # Keep proposal construction reproducible under PM auxiliary RNG.
    set.seed(as.integer(seed_use + 577L))
    prop_ctl <- ctl
    prop_ctl$defensive_df <- as.integer(max(2L, ctl$defensive_df + attempt - 1L))
    prop <- tryCatch(make_theta_proposal_adaptive(phi, surrogate, prop_ctl, callbacks), error = function(e) e)
    if (inherits(prop, "error")) {
      last_err <- conditionMessage(prop)
      next
    }

    Theta <- tryCatch(prop$draw(M_use, seed = as.integer(seed_use + 991L)), error = function(e) e)
    if (inherits(Theta, "error")) {
      last_err <- conditionMessage(Theta)
      next
    }
    Theta <- as.matrix(Theta)
    if (nrow(Theta) != M_use) {
      last_err <- "Adaptive proposal returned wrong sample size."
      next
    }

    log_py <- tryCatch(ll_parallel(Theta, data_i, loglik_fn, n_cores = 1L), error = function(e) e)
    if (inherits(log_py, "error")) {
      last_err <- conditionMessage(log_py)
      next
    }
    lprior <- tryCatch(.eval_prior(Theta), error = function(e) e)
    if (inherits(lprior, "error")) {
      last_err <- conditionMessage(lprior)
      next
    }
    lq <- tryCatch(prop$log_q(Theta), error = function(e) e)
    if (inherits(lq, "error")) {
      last_err <- conditionMessage(lq)
      next
    }

    lw <- as.numeric(log_py + lprior - lq)
    if (!all(is.finite(lw))) {
      last_err <- "Non-finite importance weights in adaptive PM estimator."
      next
    }

    a <- max(lw)
    u <- exp(lw - a)
    mu1 <- mean(u)
    mu2 <- mean(u^2)
    log_mhat <- a + log(mu1)
    ww <- u / sum(u)
    ess_abs <- 1 / sum(ww * ww)
    ess_norm <- ess_abs / M_use
    var_proxy <- max(mu2 - mu1^2, 0) / (M_use * max(mu1^2, .Machine$double.eps))
    khat <- NA_real_
    if (requireNamespace("loo", quietly = TRUE)) {
      ps <- suppressWarnings(tryCatch(loo::psis(lw), error = function(e) NULL))
      if (!is.null(ps)) {
        khat <- suppressWarnings(tryCatch(mean(loo::pareto_k_values(ps)), error = function(e) NA_real_))
      }
    }

    aux_new <- aux_state_i
    aux_new$seed <- seed0
    aux_new$M <- M_use
    aux_new$last_ok <- TRUE

    return(list(
      log_mhat = as.numeric(log_mhat),
      ess_is = as.numeric(ess_abs),
      ess_norm = as.numeric(ess_norm),
      var_proxy = as.numeric(var_proxy),
      khat = as.numeric(khat),
      M_used = as.integer(M_use),
      aux_state_i_new = aux_new
    ))
  }

  aux_fail <- aux_state_i
  aux_fail$seed <- seed0
  aux_fail$M <- M_use
  aux_fail$last_ok <- FALSE
  list(
    log_mhat = as.numeric(ctl$min_log_mhat),
    ess_is = 1.0,
    ess_norm = 1.0 / M_use,
    var_proxy = Inf,
    khat = NA_real_,
    M_used = as.integer(M_use),
    aux_state_i_new = aux_fail,
    error = last_err %||% "adaptive PM estimation failed"
  )
}


# --------------------------- unbiased estimators -------------------------
# Batch-indexed unbiased log marginal for a given φ (use in outer sampler).
log_marginal_unbiased_batch <- function(cache, phi, b, log_prior_theta_given_phi_mat, aux = NULL) {
  stopifnot(inherits(cache, "subject_cache_smcK"))
  B <- cache$batches[[b]]
  if (is.null(B)) stop("Invalid batch index: ", b)
  lprior <- log_prior_theta_given_phi_mat(B$Theta, phi, aux)
  lw <- B$log_py + lprior - B$log_r
  logsumexp(lw) - log(B$M)
}

log_marginal_unbiased_gaussian_batch <- function(cache, phi, b) {
  stopifnot(!is.null(cache$gaussian_map))
  B <- cache$batches[[b]]
  pars <- cache$gaussian_map(phi)
  lw <- B$log_py + .log_prior_gauss_mat(B$Theta, pars) - B$log_r
  logsumexp(lw) - log(B$M)
}

# Lower-variance pooled estimator using all fixed cache batches for one subject.
log_marginal_unbiased_gaussian_pooled <- function(cache, phi) {
  stopifnot(!is.null(cache$gaussian_map))
  pars <- cache$gaussian_map(phi)
  lw_all <- vector("list", cache$K)
  total_M <- 0L
  for (b in seq_len(cache$K)) {
    B <- cache$batches[[b]]
    lw_all[[b]] <- B$log_py + .log_prior_gauss_mat(B$Theta, pars) - B$log_r
    total_M <- total_M + as.integer(B$M)
  }
  lw <- unlist(lw_all, use.names = FALSE)
  logsumexp(lw) - log(total_M)
}

# Diagnostics (ESS, PSIS k) for a chosen batch
log_marginal_unbiased_diag_batch <- function(cache, phi, b, log_prior_theta_given_phi_mat, aux = NULL) {
  B <- cache$batches[[b]]
  lprior <- log_prior_theta_given_phi_mat(B$Theta, phi, aux)
  lw <- B$log_py + lprior - B$log_r
  m  <- max(lw); w <- exp(lw - m)
  ess <- (sum(w)^2) / sum(w^2)
  list(logZ = m + log(mean(w)), ESS_IS = ess, M = B$M,
       pareto_k = pareto_k_proxy(lw))
}

# Cache quality diagnostics at a given phi row.
subject_cache_quality_at_phi <- function(cache, phi, log_prior_theta_given_phi_mat = NULL) {
  stopifnot(inherits(cache, "subject_cache_smcK"))
  K <- cache$K
  per_batch <- data.frame(
    b = seq_len(K),
    M = NA_real_,
    ess_abs = NA_real_,
    ess_norm = NA_real_,
    logZ = NA_real_
  )
  for (b in seq_len(K)) {
    B <- cache$batches[[b]]
    if (!is.null(cache$gaussian_map) && is.null(log_prior_theta_given_phi_mat)) {
      pars <- cache$gaussian_map(phi)
      lprior <- .log_prior_gauss_mat(B$Theta, pars)
    } else {
      stopifnot(!is.null(log_prior_theta_given_phi_mat))
      lprior <- log_prior_theta_given_phi_mat(B$Theta, phi, aux = NULL)
    }
    lw <- B$log_py + lprior - B$log_r
    m <- max(lw)
    ww <- exp(lw - m)
    ww <- ww / sum(ww)
    ess_abs <- 1 / sum(ww * ww)
    per_batch$M[b] <- B$M
    per_batch$ess_abs[b] <- ess_abs
    per_batch$ess_norm[b] <- ess_abs / B$M
    per_batch$logZ[b] <- m + log(mean(exp(lw - m)))
  }
  list(
    per_batch = per_batch,
    summary = list(
      median_ess_abs = stats::median(per_batch$ess_abs, na.rm = TRUE),
      median_ess_norm = stats::median(per_batch$ess_norm, na.rm = TRUE),
      sd_logZ_batches = stats::sd(per_batch$logZ, na.rm = TRUE)
    )
  )
}

# ----------------------------- DA proxies --------------------------------
log_marginal_proxy_from_base <- function(cache, phi, log_prior_theta_given_phi_mat, aux = NULL) {
  base <- cache$base
  ratio <- log_prior_theta_given_phi_mat(base$theta_base, phi, aux) - base$log_p0_vec
  base$log_Z0 + logsumexp(log(base$w_bar) + ratio)
}

log_marginal_proxy_downsample_batch <- function(cache, phi, b, log_prior_theta_given_phi_mat, M_DA = 32L, aux = NULL) {
  B <- cache$batches[[b]]
  id <- if (M_DA < B$M) sample.int(B$M, M_DA) else seq_len(B$M)
  lw <- B$log_py[id] + log_prior_theta_given_phi_mat(B$Theta[id, , drop = FALSE], phi, aux) - B$log_r[id]
  logsumexp(lw) - log(length(id))
}

# ------------------------------- PSIS k ----------------------------------
pareto_k_proxy <- function(weights_log) {
  if (!requireNamespace("loo", quietly = TRUE)) return(NA_real_)
  lw <- weights_log - max(weights_log)     # stabilize; invariant to constant shift
  ps <- tryCatch(loo::psis(lw), error = function(e) NULL)
  if (is.null(ps)) return(NA_real_)
  vals <- tryCatch(loo::pareto_k_values(ps), error = function(e) NA_real_)
  if (is.null(vals)) return(NA_real_)
  mean(vals)
}

pareto_k_from_base_proxy <- function(cache, phi, log_prior_theta_given_phi_mat, aux = NULL) {
  base <- cache$base
  logw <- log_prior_theta_given_phi_mat(base$theta_base, phi, aux) - base$log_p0_vec + log(base$w_bar)
  pareto_k_proxy(logw)
}

pareto_k_from_r_proxy_batch <- function(cache, phi, b, log_prior_theta_given_phi_mat, aux = NULL) {
  B <- cache$batches[[b]]
  logw <- B$log_py + log_prior_theta_given_phi_mat(B$Theta, phi, aux) - B$log_r
  pareto_k_proxy(logw)
}

# ----------------- correlated-PM refresh: batch indices only -------------
# In the outer sampler, keep a vector b_i ∈ {1..K} (one per unit).
# Use the SAME b_i for current/proposed φ to get strong correlation.
# Occasionally call this to refresh a fraction cheaply (no model calls).
refresh_batch_indices <- function(b_vec, K, frac = 0.2) {
  n <- length(b_vec); R <- ceiling(n * .clamp(frac, 0, 1))
  if (R <= 0L) return(b_vec)
  idx <- sample.int(n, R)
  b_vec[idx] <- sample.int(K, R, replace = TRUE)
  b_vec
}

# ----------------- heavy adaptation hooks (full recompute) ---------------
enrich_cache_with_mixZ <- function(cache, new_mixZ, data, loglik_fn, sobol_seed = NULL, n_cores = 1L) {
  stopifnot(inherits(cache, "subject_cache_smcK"))
  new_mixZ <- .ensure_mix_cache(new_mixZ)
  cache$ri$mix_g <- new_mixZ
  cache$ri$mix_t <- if (cache$ri$t_eps > 0) new_mixZ else NULL
  if (!is.null(sobol_seed)) {
    base_seed <- sobol_seed
  } else {
    # reuse stored seeds if available; else derive from existing
    base_seed <- cache$batches[[1]]$sobol_seed %||% 11L
  }
  # Rebuild ALL batches with new proposal
  K <- cache$K
  seeds <- base_seed + 7919L * seq_len(K)
  for (k in seq_len(K)) {
    M_k <- cache$batches[[k]]$M
    smp <- .sample_from_ri_rqmc(M_k, cache$ri, sobol_seed = seeds[k])
    Theta_k <- smp$Theta; Z_k <- smp$Z
    log_py_k <- ll_parallel(Theta_k, data, loglik_fn, n_cores)
    log_r_k  <- .log_r_theta_mix(Theta_k, Z_k, cache$ri)
    if (!all(is.finite(log_py_k))) stop("Non-finite log p(y|θ) in enriched batch ", k)
    if (!all(is.finite(log_r_k)))  stop("Non-finite log r(θ) in enriched batch ", k)
    cache$batches[[k]] <- list(
      Theta = Theta_k, Z = Z_k,
      log_py = as.numeric(log_py_k), log_r = as.numeric(log_r_k),
      M = nrow(Theta_k), sobol_seed = seeds[k]
    )
  }
  cache
}

augment_cache_M <- function(cache, data, loglik_fn, M_add = 256L, sobol_seed = NULL, n_cores = 1L) {
  stopifnot(inherits(cache, "subject_cache_smcK"))
  K <- cache$K
  add_per <- as.integer(ceiling(M_add / K))
  base_seed <- (sobol_seed %||% cache$batches[[1]]$sobol_seed %||% 11L) + 104729L  # different stride
  for (k in seq_len(K)) {
    if (add_per <= 0L) next
    seeds_k <- base_seed + 7919L * k
    smp <- .sample_from_ri_rqmc(add_per, cache$ri, sobol_seed = seeds_k)
    Theta_new <- rbind(cache$batches[[k]]$Theta, smp$Theta)
    Z_new     <- rbind(cache$batches[[k]]$Z,     smp$Z)
    log_py_new <- c(cache$batches[[k]]$log_py, ll_parallel(smp$Theta, data, loglik_fn, n_cores))
    log_r_new  <- c(cache$batches[[k]]$log_r,  .log_r_theta_mix(smp$Theta, smp$Z, cache$ri))
    if (!all(is.finite(log_py_new))) stop("Non-finite log p(y|θ) after augment in batch ", k)
    if (!all(is.finite(log_r_new)))  stop("Non-finite log r(θ) after augment in batch ", k)
    cache$batches[[k]]$Theta <- Theta_new
    cache$batches[[k]]$Z     <- Z_new
    cache$batches[[k]]$log_py <- as.numeric(log_py_new)
    cache$batches[[k]]$log_r  <- as.numeric(log_r_new)
    cache$batches[[k]]$M      <- nrow(Theta_new)
  }
  cache
}

## --- NEW: build a prior-anchored component in Z (Fix #2) ----------------
# Requires either gaussian_map_fn (fast path) or a draw_prior_theta_fn
build_prior_anchor_mixZ <- function(Tmap, phi_anchor,
                                    gaussian_map_fn = NULL,
                                    draw_prior_theta_fn = NULL,
                                    n = 200L, eps_prior = 0.04) {
  if (is.null(gaussian_map_fn) && is.null(draw_prior_theta_fn)) {
    warning("No prior sampler provided; skipping prior-anchored component.")
    return(NULL)
  }
  # Draw θ ~ p(θ|φ_anchor)
  if (!is.null(gaussian_map_fn)) {
    pars <- gaussian_map_fn(phi_anchor, d = NULL)  # Expect pars$mu, pars$Sigma_inv
    # ### PATCH: stable SPD covariance via precision Cholesky
    L <- tryCatch(chol(pars$Sigma_inv), error = function(e) NULL)
    if (is.null(L)) stop("gaussian_map_fn returned non-SPD precision.")
    Sigma <- chol2inv(L)
    Theta <- mvtnorm::rmvnorm(n, mean = as.numeric(pars$mu), sigma = Sigma)
  } else {
    Theta <- draw_prior_theta_fn(phi_anchor, n)
  }
  Z <- Tmap$fwd(Theta)
  mu <- colMeans(Z)
  Zc <- sweep(Z, 2L, mu, `-`)
  Sig <- (t(Zc) %*% Zc) / nrow(Zc)
  Sig <- as.matrix((Sig + t(Sig))/2)
  mix <- list(meansZ = list(as.numeric(mu)), covsZ = list(Sig), wZ = 1)
  mix$cache <- prep_mix_cache(mix$meansZ, mix$covsZ, mix$wZ)
  # ### PATCH: return eig diagnostics for logging
  ez <- eigen(Sig, symmetric = TRUE, only.values = TRUE)$values
  list(mix = mix, eps_prior = eps_prior,
       eigZ = ez, condZ = max(ez)/max(min(ez), .Machine$double.eps))
}

## --- NEW: batch-aware enrichment recipe (Fixes #2–#5) -------------------
# ### PATCH: surgical enrichment on ACTIVE batch only, partial replacement + gating
enrich_subject_cache_with_anchor_surgical <- function(cache, phi_probe, b,
                                                      log_prior_theta_given_phi_mat = NULL,
                                                      gaussian_map_fn = NULL,
                                                      data, loglik_fn,
                                                      eps_prior_base = 0.04,
                                                      weak_inflate_factor = 3.0,
                                                      elite_q = 0.4, G = 12L, cov_infl = 4.0,
                                                      replace_frac = NULL,
                                                      sobol_seed = NULL, n_cores = 1L,
                                                      verbose = TRUE,
                                                      subj_id = NA_integer_,
                                                      phi_panel = NULL) {
  B <- cache$batches[[b]]
  if (B$M <= 8L) return(cache)  # pathologically small

  # --- helpers for diagnostics from log-weights
  .diag_from_lw <- function(lw) {
    m <- max(lw); u <- exp(lw - m)
    mu1 <- mean(u); mu2 <- mean(u^2)
    ess <- (sum(u)^2)/sum(u^2)
    mcse <- sqrt( max(mu2 - mu1^2, 0) / (length(u) * max(mu1, .Machine$double.eps)^2) )
    k <- pareto_k_proxy(lw)
    list(ess_norm = ess / length(u), k = k, logZ = m + log(mu1), mcse = mcse)
  }

  # 1) Pre diagnostics at φ_probe
  if (!is.null(cache$gaussian_map) && is.null(log_prior_theta_given_phi_mat)) {
    pars <- cache$gaussian_map(phi_probe)
    lprior <- .log_prior_gauss_mat(B$Theta, pars)
  } else {
    stopifnot(!is.null(log_prior_theta_given_phi_mat))
    lprior <- log_prior_theta_given_phi_mat(B$Theta, phi_probe, aux = NULL)
  }
  lw_pre <- B$log_py + lprior - B$log_r
  pre <- .diag_from_lw(lw_pre)

  # --- adapt replacement fraction (more conservative)
  if (is.null(replace_frac)) {
    if (is.finite(pre$k) && (pre$k > 1.0 || pre$ess_norm < 0.06)) {
      replace_frac <- 0.20
    } else if (is.finite(pre$k) && pre$k > 0.8) {
      replace_frac <- 0.22
    } else {
      replace_frac <- 0.20
    }
  }

  # 2) Fit elite mixture in Z at φ_probe
  a  <- max(lw_pre); w <- exp(lw_pre - a); w <- w / sum(w)
  mix_elite <- NULL
  if (exists("fit_elite_mixture_Z", mode = "function")) {
    mix_elite <- fit_elite_mixture_Z(B$Z, w, elite_quantile = elite_q, G = G,
                                     cov_inflation = cov_infl, housekeeping = TRUE,
                                     min_G_keep = 2, merge_thresh = 0.10,
                                     verbose = FALSE, lambda = 1.0)
  }
  if (is.null(mix_elite) || .is_empty_mix(mix_elite)) {
    mu <- colSums(B$Z * w)
    Zc <- sweep(B$Z, 2L, mu, `-`)
    Sig <- t(Zc) %*% (Zc * w); Sig <- as.matrix((Sig + t(Sig))/2)
    mix_elite <- list(meansZ = list(as.numeric(mu)), covsZ = list(Sig), wZ = 1)
    mix_elite$cache <- prep_mix_cache(mix_elite$meansZ, mix_elite$covsZ, mix_elite$wZ)
  }

  # 3) PC-based weak inflation (based on weighted covariance in Z)
  Zw_mu <- colSums(B$Z * w)
  Zw_c  <- sweep(B$Z, 2L, Zw_mu, `-`)
  S_w   <- (t(Zw_c) %*% (Zw_c * w))
  S_w   <- as.matrix((S_w + t(S_w))/2)
  ev <- eigen(S_w, symmetric = TRUE)
  lam <- pmax(ev$values, 0)
  thr <- stats::quantile(lam, probs = 0.25, names = FALSE)
  weak_mask <- lam <= thr
  # build R = V diag(r) V^T with r = sqrt(f) on weak PCs else 1
  r <- rep(1, length(lam)); r[weak_mask] <- sqrt(max(weak_inflate_factor, 1))
  R <- ev$vectors %*% diag(r, nrow = length(r)) %*% t(ev$vectors)
  # apply to each component covariance: S' = R S R
  mix_elite_pc <- list(
    meansZ = mix_elite$meansZ,
    covsZ  = lapply(mix_elite$covsZ, function(Si) {
      Si <- as.matrix(Si); Si2 <- R %*% Si %*% R; as.matrix((Si2 + t(Si2))/2)
    }),
    wZ = mix_elite$wZ
  )
  mix_elite_pc$cache <- prep_mix_cache(mix_elite_pc$meansZ, mix_elite_pc$covsZ, mix_elite_pc$wZ)

  # 4) Prior-anchored component with ADAPTIVE weight (gentler slope + lower cap)
  eps_prior_adapt <- min(eps_prior_base + 0.05 * max(0, (pre$k %||% 0) - 0.7), 0.08)
  anchor <- build_prior_anchor_mixZ(cache$ri$Tmap, phi_probe,
                                    gaussian_map_fn = if (!is.null(gaussian_map_fn)) function(phi, d=NULL)
                                      gaussian_map_fn(phi, ncol(B$Theta)) else NULL,
                                    draw_prior_theta_fn = NULL,
                                    n = 200L, eps_prior = eps_prior_adapt)
  mix_final <- mix_elite_pc
  if (!is.null(anchor)) {
    # Guard: if Z-cov of anchor is ill-conditioned, try halving eps once; else skip
    if (is.finite(anchor$condZ) && anchor$condZ >= 1e6 && eps_prior_adapt > 0) {
      anchor_half <- build_prior_anchor_mixZ(cache$ri$Tmap, phi_probe,
                          gaussian_map_fn = if (!is.null(gaussian_map_fn)) function(phi, d=NULL)
                            gaussian_map_fn(phi, ncol(B$Theta)) else NULL,
                          draw_prior_theta_fn = NULL,
                          n = 200L, eps_prior = 0.5 * eps_prior_adapt)
      if (!is.null(anchor_half) && is.finite(anchor_half$condZ) && anchor_half$condZ < 1e6) {
        anchor <- anchor_half
      } else {
        if (verbose) cat(sprintf("    [Enrich warn] Anchor skipped due to condZ=%.2e\n", anchor$condZ %||% NA_real_))
        anchor <- NULL
      }
    }
    if (!is.null(anchor)) {
      mix_final <- blend_mixes(mix_elite_pc, anchor$mix, eps = anchor$eps_prior)
    }
  }

  # 5) Propose SURGICAL replacement of worst rows in active batch only
  Rn <- max(1L, floor(B$M * replace_frac))
  ord <- order(lw_pre, decreasing = FALSE)
  idx_rep <- ord[seq_len(Rn)]
  seed_local <- sobol_seed %||% (B$sobol_seed + 1L)
  # draw replacement from improved proposal
  smp <- .sample_from_ri_rqmc(Rn, list(mix_g = mix_final, mix_t = if (cache$ri$t_eps > 0) mix_final else NULL,
                                       t_eps = cache$ri$t_eps, t_df = cache$ri$t_df, Tmap = cache$ri$Tmap),
                              sobol_seed = seed_local)
  Theta_new <- smp$Theta; Z_new <- smp$Z
  colnames(Theta_new) <- colnames(cache$base$theta_base)
  log_py_new <- ll_parallel(Theta_new, data, loglik_fn, n_cores)
  log_r_new  <- .log_r_theta_mix(Theta_new, Z_new, list(mix_g = mix_final, mix_t = if (cache$ri$t_eps > 0) mix_final else NULL,
                                                        t_eps = cache$ri$t_eps, t_df = cache$ri$t_df, Tmap = cache$ri$Tmap))

  # build candidate batch (copy, then patch)
  B_cand <- B
  B_cand$Theta[idx_rep, ] <- Theta_new
  B_cand$Z[idx_rep, ]     <- Z_new
  B_cand$log_py[idx_rep]  <- as.numeric(log_py_new)
  B_cand$log_r[idx_rep]   <- as.numeric(log_r_new)
  # keep the original sobol_seed for batch b (only subset changed)

  # 6) Post diagnostics at φ_probe and optional panel
  diag_at_phi <- function(phi_row, Bobj) {
    if (!is.null(cache$gaussian_map) && is.null(log_prior_theta_given_phi_mat)) {
      pars <- cache$gaussian_map(phi_row)
      lprior_x <- .log_prior_gauss_mat(Bobj$Theta, pars)
    } else {
      lprior_x <- log_prior_theta_given_phi_mat(Bobj$Theta, phi_row, aux = NULL)
    }
    lw_x <- Bobj$log_py + lprior_x - Bobj$log_r
    .diag_from_lw(lw_x)
  }
  post <- diag_at_phi(phi_probe, B_cand)
  # Optional: panel gating with additional φ values
  panel_list <- list(phi_probe)
  if (!is.null(phi_panel)) {
    if (is.matrix(phi_panel)) {
      for (k in seq_len(nrow(phi_panel))) panel_list[[length(panel_list)+1L]] <- phi_panel[k, , drop = TRUE]
    } else if (is.list(phi_panel)) {
      panel_list <- c(panel_list, phi_panel)
    }
  }
  pre_panel <- list(diag_at_phi(phi_probe, B))
  post_panel <- list(post)
  if (length(panel_list) > 1L) {
    for (idx in 2:length(panel_list)) {
      pre_panel[[idx]]  <- diag_at_phi(panel_list[[idx]], B)
      post_panel[[idx]] <- diag_at_phi(panel_list[[idx]], B_cand)
    }
  }
  avg_pre_k   <- mean(vapply(pre_panel, function(x) x$k %||% NA_real_, numeric(1)), na.rm = TRUE)
  avg_post_k  <- mean(vapply(post_panel, function(x) x$k %||% NA_real_, numeric(1)), na.rm = TRUE)
  avg_pre_ess <- mean(vapply(pre_panel, function(x) x$ess_norm, numeric(1)), na.rm = TRUE)
  avg_post_ess<- mean(vapply(post_panel, function(x) x$ess_norm, numeric(1)), na.rm = TRUE)

  # 7) Logs + gating / rollback
  if (verbose) {
    pc_inflated <- sum(weak_mask)
    ez <- anchor$eigZ %||% NA
    sid <- if (is.na(subj_id)) "?" else as.character(subj_id)
    cat(sprintf("  [Enrich try] subj=%s batch=%d\n", sid, b))
    cat(sprintf("    pre:  ESS/M=%.3f  k=%.2f  logZ=%.4f\n", pre$ess_norm, pre$k %||% NA_real_, pre$logZ))
    cat(sprintf("    mix:  eps_prior=%.3f  weak_inflate=%.2f  #PC_inflated=%d\n",
                eps_prior_adapt, weak_inflate_factor, pc_inflated))
    if (is.numeric(ez)) {
      cat(sprintf("         anchor eig(Z): min=%.3e med=%.3e max=%.3e  cond=%.2e\n",
                  min(ez), stats::median(ez), max(ez), anchor$condZ %||% NA_real_))
    }
    dZ  <- post$logZ - pre$logZ
    cat(sprintf("    post: ESS/M=%.3f  k=%.2f  logZ=%.4f  ΔlogZ=%.4f  (MCSE≈%.3f)\n",
                post$ess_norm, post$k %||% NA_real_, post$logZ, dZ, post$mcse))
  }

  # improve_k   <- (is.finite(pre$k) && is.finite(post$k) && (pre$k - post$k) >= 0.15) || (!is.finite(pre$k) && is.finite(post$k))
  # improve_ess <- (post$ess_norm - pre$ess_norm) >= 0.05
  # stable_dZ   <- abs((post$logZ - pre$logZ)) <= (3 * post$mcse + 3e-3)

  # Tightened acceptance using panel averages and guards
  k_drop    <- avg_pre_k - avg_post_k
  ess_gain  <- avg_post_ess - avg_pre_ess
  dz <-  post$logZ - pre$logZ
  safe_dz   <- is.finite(dz) && is.finite(post$mcse) && (abs(dz) <= 2 * post$mcse)
  suspicious_dz <- is.finite(dz) && is.finite(post$mcse) && (abs(dz) > 3 * post$mcse)
  any_k_worse <- FALSE
  for (ii in seq_along(pre_panel)) {
    ki_pre <- pre_panel[[ii]]$k %||% NA_real_
    ki_post<- post_panel[[ii]]$k %||% NA_real_
    if (is.finite(ki_pre) && is.finite(ki_post) && ((ki_post - ki_pre) >= 0.03)) { any_k_worse <- TRUE; break }
  }
  accept <- (!suspicious_dz) && safe_dz && (!any_k_worse) && (k_drop >= 0.08) && (ess_gain >= 0.02)

  decision <- if (accept) "ACCEPT" else "REJECT (rollback)"
  if (verbose) {
    cat(sprintf("    decision: %s\n", decision))
    if (length(panel_list) > 1L) {
      cat(sprintf("    panel avg: ESS/M=%.3f  k=%.2f\n", avg_post_ess, avg_post_k))
    }
  }
  if (!accept) return(cache)

  # Commit: update only batch b; keep proposal ri unchanged globally (surgical)
  cache$batches[[b]] <- B_cand
  cache
}

# ----------------------------- helpers -----------------------------------
# `%||%` is available from smc_core
