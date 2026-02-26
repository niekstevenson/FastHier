# -----------------------------------------------------------------------------
# build_inner_cache(): run enhanced_smc_elite once and package a reusable cache
# -----------------------------------------------------------------------------
# What it stores (per subject):
# - theta:      SMC particles in theta-space (matrix K x d)
# - w:          normalized SMC weights (length K)
# - log_p0:     log base-prior density at each particle (for IS reweighting)
# - logZ0:      inner evidence estimate (log marginal likelihood under base prior)
# - mcse_logZ:  MCSE for logZ0
# - prior:      list(mu, Sigma) to document the base prior used
# - mixZ:       optional elite mixture in Z-space (for later independence moves / DA)
# - transport:  optional light-weight summary of the transport (no closures)
# - meta:       sampler metadata
#
# Notes:
# - Default avoids storing the full transport closures (can be large); set store_transport=TRUE to keep them.
# - You can reduce size via target_size (weighted resample to equal weights).
# - dedup=TRUE collapses exact duplicates (common after resampling).
# -----------------------------------------------------------------------------

build_inner_cache <- function(
    data, loglik_fn,
    mu0, Sigma0,
    M = 5000L,
    seed = 123,
    # size & memory controls
    target_size = NULL,         # if not NULL and < M: systematic resample to equal weights
    dedup = TRUE,               # collapse exact duplicates before optional downsampling
    dedup_digits = 12L,         # key precision for duplicate detection
    store_transport = FALSE,    # keep full Tmap closures (large); otherwise store summary only
    store_mixZ = TRUE,          # fit & store elite mixture in Z (useful later)
    mix_G = 12L,                # components for elite mixture
    mix_elite_q = 0.4,          # elite quantile for mixture fitting
    verbose = TRUE,
    ...                         # forwarded to enhanced_smc_elite()
) {
  if (verbose) cat("Running enhanced_smc_elite to build cache...\n")
  fit <- enhanced_smc_elite(
    data        = data,
    loglik_fn   = loglik_fn,
    mu_ref      = mu0,
    Sigma_ref   = Sigma0,
    M           = M,
    seed        = seed,
    verbose     = verbose,
    ...
  )
  Theta <- fit$Theta
  w     <- as.numeric(fit$w)
  stopifnot(nrow(Theta) == length(w))

  # (1) base prior log-density for reweighting later
  log_p0 <- mvtnorm::dmvnorm(Theta, mean = mu0, sigma = Sigma0, log = TRUE)

  # (2) optional exact-duplicate collapse (fast; exact equality from resampling)
  if (isTRUE(dedup)) {
    if (verbose) cat("Collapsing exact duplicate particles...\n")
    # stable grouping by sprintf string keys (like your ll-cache)
    key <- apply(Theta, 1L, function(x) paste(sprintf("%.*g", dedup_digits, x), collapse = ";"))
    ord <- order(key)
    Theta <- Theta[ord, , drop = FALSE]; w <- w[ord]; log_p0 <- log_p0[ord]; key <- key[ord]
    first <- !duplicated(key)
    if (any(!first)) {
      grp <- cumsum(first)
      w   <- as.numeric(tapply(w, grp, sum))
      Theta <- Theta[first, , drop = FALSE]
      # recompute log_p0 on representatives (robust if keys grouped very close values)
      log_p0 <- mvtnorm::dmvnorm(Theta, mean = mu0, sigma = Sigma0, log = TRUE)
    }
    w <- w / sum(w)
  }

  # (3) optional weighted resample to a smaller cache with equal weights
  sys_resample <- function(w, n_out) {
    # systematic resampling indices
    N <- length(w)
    cs <- cumsum(w)
    u0 <- runif(1, max(0, 1e-12), 1.0 / n_out)
    pts <- u0 + (0:(n_out - 1)) / n_out
    idx <- findInterval(pts, cs) + 1L
    pmax(pmin(idx, N), 1L)
  }
  if (!is.null(target_size) && is.finite(target_size) && target_size > 0L && target_size < nrow(Theta)) {
    if (verbose) cat(sprintf("Downsampling cache to %d particles (equal weights)...\n", target_size))
    set.seed(seed + 1L)
    idx <- sys_resample(w / sum(w), as.integer(target_size))
    Theta <- Theta[idx, , drop = FALSE]
    w     <- rep(1 / length(idx), length(idx))
    log_p0 <- log_p0[idx]
  } else {
    # always normalize
    w <- w / sum(w)
  }

  # (4) optional elite mixture in Z-space for later independence/DA proposals
  mixZ <- NULL
  if (isTRUE(store_mixZ)) {
    if (verbose) cat("Fitting elite mixture in Z-space for cache (surrogate/proposals)...\n")
    # weights for elite fitting (mild sharpening like sampler)
    gamma_sharp <- 0.7
    w_fit <- (w ^ gamma_sharp); w_fit <- w_fit / sum(w_fit)
    # reuse the same fitter you use during SMC
    mixZ <- try(
      fit_elite_mixture_Z(
        Z = fit$Z, w = w_fit,
        elite_quantile = mix_elite_q,
        G = mix_G,
        cov_inflation = 3.5,
        warm_start_mixZ = NULL,
        em_itmax = 8,
        housekeeping = TRUE,
        min_G_keep = 2,
        merge_thresh = 0.10,
        verbose = FALSE,
        lambda = fit$final_lambda
      ),
      silent = TRUE
    )
    if (inherits(mixZ, "try-error")) mixZ <- NULL
  }

  # (5) transport storage: full closures are big; default to a light summary
  transport <- NULL
  transport_summary <- NULL
  if (isTRUE(store_transport)) {
    transport <- fit$transport
  } else {
    # keep only meta (enough for diagnostics / quick recon)
    tm <- fit$transport
    if (!is.null(tm) && !is.null(tm$meta)) {
      transport_summary <- tm$meta   # list(C, U, muY)
    }
  }

  cache <- list(
    theta      = Theta,
    w          = w,
    log_p0     = log_p0,
    logZ0      = as.numeric(fit$log_evidence),
    mcse_logZ  = as.numeric(fit$mcse_logZ),
    prior      = list(mu = as.numeric(mu0), Sigma = as.matrix(Sigma0)),
    param_names = colnames(Theta),
    mixZ       = mixZ,                 # mixture in Z (optional)
    transport  = if (isTRUE(store_transport)) transport else NULL,
    transport_summary = if (!isTRUE(store_transport)) transport_summary else NULL,
    meta = list(
      data          = data,
      loglik_fn     = loglik_fn,
      M_in          = M,
      N_cache       = nrow(Theta),
      seed          = seed,
      rounds        = fit$meta$rounds,
      ess_final     = fit$meta$ess,
      final_lambda  = fit$final_lambda,
      rw_scale_final= fit$meta$rw_scale_final,
      pcn_beta_final= fit$meta$pcn_beta_final,
      built_at      = Sys.time()
    )
  )
  class(cache) <- c("inner_cache", class(cache))
  if (verbose) cat("Cache built.\n")
  cache
}
