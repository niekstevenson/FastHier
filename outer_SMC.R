# ========================================================================
# Outer SMC over φ with Option-B inner caches (K RQMC batches)
# - Pseudo-marginal exactness via batch-index auxiliary variables
# - Correlated PM: keep b_idx fixed within moves, refresh sparsely
# - Adaptive tempering via rCESS; resampling (systematic | multinomial)
# - Stable numerics; finite checks; reproducible seeding
# - Diagnostics: per-subject ESS_i/M & PSIS k from active batches
# - Records per-round log-evidence increments
# ========================================================================

source("new_SMC_cache.R")

# Use shared core utilities (numerics, resampling, rCESS, etc.)
sc_path <- file.path(getwd(), "smc_core.R")
if (file.exists(sc_path)) source(sc_path)

suppressPackageStartupMessages({
  library(Matrix)
  library(matrixStats)
})

# ------------------------------ numerics ---------------------------------
# (moved to smc_core: .clamp, logsumexp, logsumexp_w, `%||%`)

# ------------------------------ rCESS step -------------------------------
# (moved to smc_core: rCESS_stat, next_lambda_via_rCESS_stat)

# ------------------------------ resampling -------------------------------
.resample_multinomial <- function(w) {
  # unbiased for evidence; more variance
  N <- length(w)
  sample.int(N, size = N, replace = TRUE, prob = w)
}

# Use smc_core::stratified_resample_sorted for low-variance systematic resampling

# ---------------------------- φ utilities --------------------------------

.choose_probe_phi <- function(phi, w, logprior, loglik, lambda, mode = c("wmean","best")) {
  mode <- match.arg(mode)
  if (mode == "wmean") {
    ww <- pmax(w, 0); ww <- ww / sum(ww)
    drop(t(phi) %*% ww)
  } else {
    idx <- which.max(logprior + lambda * loglik)
    phi[idx, , drop = TRUE]
  }
}

# --------------------------- Inner log-likelihood ------------------------
# Batch-indexed pseudo-marginal log-lik: sum_i log \hat Z_i(φ | b_i)
.outer_loglik_hat_batch <- function(phi_row, caches, b_idx,
                                    log_prior_theta_given_phi_mat = NULL,
                                    gaussian_map_fn = NULL) {
  S <- length(caches)
  s <- 0.0
  if (!is.null(gaussian_map_fn)) {
    for (i in seq_len(S)) {
      s <- s + log_marginal_unbiased_gaussian_batch(caches[[i]], phi_row, b_idx[i])
    }
  } else {
    stopifnot(!is.null(log_prior_theta_given_phi_mat))
    for (i in seq_len(S)) {
      s <- s + log_marginal_unbiased_batch(caches[[i]], phi_row, b_idx[i],
                                           log_prior_theta_given_phi_mat)
    }
  }
  s
}

.outer_loglik_hat_particles <- function(Phi, caches, b_state,
                                        log_prior_theta_given_phi_mat = NULL,
                                        gaussian_map_fn = NULL,
                                        pm_mode = c("strict", "fast")) {
  pm_mode <- match.arg(pm_mode)
  Phi <- as.matrix(Phi)
  N <- nrow(Phi)
  if (N <= 0L) return(numeric(0))

  if (pm_mode == "fast") {
    stopifnot(is.vector(b_state), length(b_state) == length(caches))
    return(apply(Phi, 1L, function(p)
      .outer_loglik_hat_batch(p, caches, b_state,
                              log_prior_theta_given_phi_mat, gaussian_map_fn)))
  }

  stopifnot(is.matrix(b_state), nrow(b_state) == N, ncol(b_state) == length(caches))
  out <- numeric(N)
  for (i in seq_len(N)) {
    out[i] <- .outer_loglik_hat_batch(Phi[i, , drop = TRUE], caches, b_state[i, ],
                                      log_prior_theta_given_phi_mat, gaussian_map_fn)
  }
  out
}

.init_aux_batch_state <- function(caches, M, pm_mode = c("strict", "fast")) {
  pm_mode <- match.arg(pm_mode)
  S <- length(caches)
  if (pm_mode == "fast") {
    return(vapply(caches, function(ci) sample.int(ci$K, 1L), integer(1)))
  }
  B <- matrix(0L, nrow = M, ncol = S)
  for (j in seq_len(S)) {
    B[, j] <- sample.int(caches[[j]]$K, M, replace = TRUE)
  }
  B
}

.probe_aux_row <- function(phi, w, logprior, loglik, lambda, b_state) {
  idx <- which.max(logprior + lambda * loglik)
  if (is.matrix(b_state)) b_state[idx, ] else b_state
}

.refresh_aux_state_block_mh <- function(phi, loglik, b_state, lambda,
                                        caches,
                                        log_prior_theta_given_phi_mat = NULL,
                                        gaussian_map_fn = NULL,
                                        frac = 0.10,
                                        rng_seed = NULL) {
  if (!is.matrix(b_state)) {
    return(list(b_state = b_state, loglik = loglik, acc_rate = NA_real_))
  }
  if (!is.null(rng_seed)) set.seed(rng_seed)
  frac <- .clamp(frac, 0, 1)
  if (frac <= 0) return(list(b_state = b_state, loglik = loglik, acc_rate = 0))

  N <- nrow(phi)
  S <- ncol(b_state)
  R <- max(1L, ceiling(frac * S))
  nacc <- 0L
  for (i in seq_len(N)) {
    idx <- sample.int(S, R)
    b_prop <- b_state[i, ]
    for (j in idx) b_prop[j] <- sample.int(caches[[j]]$K, 1L)
    ll_prop <- .outer_loglik_hat_batch(phi[i, , drop = TRUE], caches, b_prop,
                                       log_prior_theta_given_phi_mat, gaussian_map_fn)
    loga <- lambda * (ll_prop - loglik[i])
    if (!is.finite(loga)) loga <- -Inf
    if (log(runif(1)) < min(0, loga)) {
      b_state[i, ] <- b_prop
      loglik[i] <- ll_prop
      nacc <- nacc + 1L
    }
  }
  list(b_state = b_state, loglik = loglik, acc_rate = nacc / max(N, 1L))
}

# --------------- Per-subject IS diagnostics at a given φ ----------------
.cache_logw_at_phi_batch <- function(cache, phi_row, b,
                                     log_prior_theta_given_phi_mat = NULL) {
  B <- cache$batches[[b]]
  if (!is.null(cache$gaussian_map) && is.null(log_prior_theta_given_phi_mat)) {
    pars <- cache$gaussian_map(phi_row)
    lprior <- .log_prior_gauss_mat(B$Theta, pars)
  } else {
    stopifnot(!is.null(log_prior_theta_given_phi_mat))
    lprior <- log_prior_theta_given_phi_mat(B$Theta, phi_row, aux = NULL)
  }
  B$log_py + lprior - B$log_r
}

.is_stats_one_batch <- function(cache, phi_row, b,
                                log_prior_theta_given_phi_mat = NULL) {
  lw <- .cache_logw_at_phi_batch(cache, phi_row, b, log_prior_theta_given_phi_mat)
  a  <- max(lw)
  w  <- exp(lw - a)
  logZ_hat <- a + log(mean(w))
  w <- w / sum(w)
  ess <- 1 / sum(w * w)
  khat <- NA_real_
  if (requireNamespace("loo", quietly = TRUE)) {
    ps <- tryCatch(loo::psis(lw), error = function(e) NULL)
    if (!is.null(ps)) {
      khat <- tryCatch(mean(loo::pareto_k_values(ps)), error = function(e) NA_real_)
    }
  }
  list(logZ_hat = logZ_hat, ess = ess, M = cache$batches[[b]]$M, k = khat)
}

.is_diagnostics_batch <- function(caches, phi_row, b_idx,
                                  log_prior_theta_given_phi_mat = NULL) {
  S <- length(caches)
  out <- vector("list", S)
  for (i in seq_len(S)) {
    out[[i]] <- .is_stats_one_batch(caches[[i]], phi_row, b_idx[i],
                                    log_prior_theta_given_phi_mat)
  }
  data.frame(
    i = seq_len(S),
    logZ_hat = vapply(out, `[[`, numeric(1), "logZ_hat"),
    ess      = vapply(out, `[[`, numeric(1), "ess"),
    M        = vapply(out, `[[`, numeric(1), "M"),
    psis_k   = vapply(out, `[[`, numeric(1), "k")
  )
}

# -------- Replicate SD of log Z_i(φ) over batches (no model calls) -------
.rep_logZi_sd_over_batches <- function(caches, phi_row, B_reps = 2L,
                                       log_prior_theta_given_phi_mat = NULL,
                                       gaussian_map_fn = NULL,
                                       rng_seed = NULL) {
  if (!is.null(rng_seed)) set.seed(rng_seed)
  S <- length(caches)
  out_sd <- numeric(S)
  for (i in seq_len(S)) {
    K <- caches[[i]]$K
    if (B_reps <= 1L || K <= 1L) {
      out_sd[i] <- NA_real_
    } else {
      vals <- numeric(B_reps)
      for (b in seq_len(B_reps)) {
        bb <- sample.int(K, 1L)
        if (!is.null(gaussian_map_fn)) {
          vals[b] <- log_marginal_unbiased_gaussian_batch(caches[[i]], phi_row, bb)
        } else {
          vals[b] <- log_marginal_unbiased_batch(caches[[i]], phi_row, bb,
                                                 log_prior_theta_given_phi_mat)
        }
      }
      out_sd[i] <- stats::sd(vals)
    }
  }
  out_sd
}

# --------- batch-index refresh (cheap correlated PM; no model calls) -----
# If the cache-side helper isn't available, define a local fallback.
if (!exists("refresh_batch_indices", mode = "function")) {
  refresh_batch_indices <- function(b_vec, K, frac = 0.2) {
    n <- length(b_vec); R <- ceiling(n * .clamp(frac, 0, 1))
    if (R <= 0L) return(b_vec)
    idx <- sample.int(n, R)
    b_vec[idx] <- sample.int(K, R, replace = TRUE)
    b_vec
  }
}

# ------------------------------ Rejuvenation -----------------------------
# One sweep of adaptive MH on φ with RW + independence mixture (correlated PM; b_idx fixed)
.rejuvenate_rw_batch <- function(phi, loglik, logprior, w, lambda,
                                 caches, b_state, logprior_phi,
                                 log_prior_theta_given_phi_mat = NULL,
                                 gaussian_map_fn = NULL,
                                 pm_mode = c("strict", "fast"),
                                 n_moves = 2L,
                                 rw_scale = 1.0,
                                 weak_dim_idx = integer(0),
                                 rw_expand_factor = 2.5,
                                 elite_mix = NULL,
                                 hist_mix = NULL,
                                 indep_prob = 0.30,
                                 indep_t_df = 4L,
                                 indep_t_prob = 0.85,
                                 rng_seed = NULL) {
  pm_mode <- match.arg(pm_mode)
  if (!is.null(rng_seed)) set.seed(rng_seed)
  N <- nrow(phi); d <- ncol(phi)

  S <- weighted_cov(phi, w)
  S <- as.matrix(Matrix::nearPD(S, conv.tol = 1e-7)$mat)
  diag(S) <- pmax(diag(S), 1e-8)
  ev <- eigen(S, symmetric = TRUE)
  lamv <- pmax(ev$values, 1e-8)
  if (length(weak_dim_idx)) lamv[weak_dim_idx] <- lamv[weak_dim_idx] * rw_expand_factor
  Sprop <- ev$vectors %*% diag(lamv, d) %*% t(ev$vectors)
  L <- tryCatch(chol(Sprop + diag(1e-8, d)), error = function(e) chol(Sprop + diag(1e-6, d)))
  step <- (rw_scale / sqrt(max(d, 1))) * L

  rlogsumexp2 <- function(a, b) { m <- pmax(a, b); m + log(exp(a - m) + exp(b - m)) }
  log_q_mixture <- function(Phi, elite, hist, indep_t_prob, indep_t_df) {
    if (.is_empty_mix(elite)) elite <- .default_std_normal_mix(ncol(Phi))
    lq_el_n <- gmm_logpdf_Z_vec(Phi, elite$meansZ, elite$cache)
    lq_el_t <- dmvt_mixture_logpdf_Z_vec(Phi, elite$meansZ, elite$cache, indep_t_df)
    lq_el   <- rlogsumexp2(lq_el_n + log1p(-indep_t_prob), lq_el_t + log(indep_t_prob))
    if (is.null(hist) || is.null(hist$mix) || .is_empty_mix(hist$mix)) return(lq_el)
    hm <- hist$mix
    lq_hi_n <- gmm_logpdf_Z_vec(Phi, hm$meansZ, hm$cache)
    lq_hi_t <- dmvt_mixture_logpdf_Z_vec(Phi, hm$meansZ, hm$cache, indep_t_df)
    lq_hi   <- rlogsumexp2(lq_hi_n + log1p(-indep_t_prob), lq_hi_t + log(indep_t_prob))
    rlogsumexp2(lq_el + log1p(-hist$prob), lq_hi + log(hist$prob))
  }
  sample_from_q_phi <- function(M, elite, hist, indep_t_prob, indep_t_df) {
    if (.is_empty_mix(elite)) elite <- .default_std_normal_mix(d)
    use_hist <- !is.null(hist) && !is.null(hist$mix) && !.is_empty_mix(hist$mix)
    hprob <- if (use_hist) hist$prob else 0.0
    take_hist <- runif(M) < hprob
    t_flags   <- runif(M) < indep_t_prob
    out <- matrix(NA_real_, M, d)
    em <- elite; hm <- if (use_hist) hist$mix else NULL
    n_hi_t <- sum(take_hist & t_flags)
    n_hi_n <- sum(take_hist & !t_flags)
    n_el_t <- sum(!take_hist & t_flags)
    n_el_n <- sum(!take_hist & !t_flags)
    if (n_hi_n) out[take_hist & !t_flags, ] <- sample_gmm_Z_qmc(n_hi_n, hm$meansZ, hm$cache)
    if (n_hi_t) out[take_hist &  t_flags, ] <- rmvt_mixture_Z_qmc(n_hi_t, hm$meansZ, hm$cache, nu = indep_t_df)
    if (n_el_n) out[!take_hist & !t_flags, ] <- sample_gmm_Z_qmc(n_el_n, em$meansZ, em$cache)
    if (n_el_t) out[!take_hist &  t_flags, ] <- rmvt_mixture_Z_qmc(n_el_t, em$meansZ, em$cache, nu = indep_t_df)
    colnames(out) <- colnames(phi)
    out
  }

  accepted <- 0L
  nacc_rw <- 0L; nacc_id <- 0L
  nprop_rw <- 0L; nprop_id <- 0L

  eval_ll <- function(Phi, b_sub) {
    .outer_loglik_hat_particles(
      Phi = Phi,
      caches = caches,
      b_state = b_sub,
      log_prior_theta_given_phi_mat = log_prior_theta_given_phi_mat,
      gaussian_map_fn = gaussian_map_fn,
      pm_mode = pm_mode
    )
  }

  for (m in seq_len(n_moves)) {
    choose_id <- runif(N) < indep_prob
    idx_id <- which(choose_id)
    idx_rw <- which(!choose_id)

    # Independence proposals (independent of current state; keep b_idx fixed for PM)
    if (length(idx_id)) {
      prop_id <- sample_from_q_phi(length(idx_id), elite_mix, hist_mix, indep_t_prob, indep_t_df)
      lp_prop_id <- apply(prop_id, 1L, logprior_phi)
      b_id <- if (pm_mode == "strict") b_state[idx_id, , drop = FALSE] else b_state
      ll_prop_id <- eval_ll(prop_id, b_id)
      lq_curr <- log_q_mixture(phi[idx_id, , drop = FALSE], elite_mix, hist_mix, indep_t_prob, indep_t_df)
      lq_prop <- log_q_mixture(prop_id, elite_mix, hist_mix, indep_t_prob, indep_t_df)
      a_id <- (lp_prop_id - logprior[idx_id]) + lambda * (ll_prop_id - loglik[idx_id]) + (lq_curr - lq_prop)
      u_id <- log(runif(length(idx_id)))
      acc_id_idx <- which(u_id < pmin(0, a_id))
      if (length(acc_id_idx)) {
        which_acc <- idx_id[acc_id_idx]
        phi[which_acc, ]    <- prop_id[acc_id_idx, , drop = FALSE]
        loglik[which_acc]   <- ll_prop_id[acc_id_idx]
        logprior[which_acc] <- lp_prop_id[acc_id_idx]
        accepted <- accepted + length(acc_id_idx)
        nacc_id <- nacc_id + length(acc_id_idx)
      }
      nprop_id <- nprop_id + length(idx_id)
    }

    # Random-walk proposals
    if (length(idx_rw)) {
      Z <- matrix(rnorm(length(idx_rw) * d), length(idx_rw), d)
      prop <- phi[idx_rw, , drop = FALSE] + Z %*% t(step)
      lp_prop <- apply(prop, 1L, logprior_phi)
      b_rw <- if (pm_mode == "strict") b_state[idx_rw, , drop = FALSE] else b_state
      ll_prop <- eval_ll(prop, b_rw)
      a <- (lp_prop - logprior[idx_rw]) + lambda * (ll_prop - loglik[idx_rw])
      u <- log(runif(length(idx_rw)))
      acc <- (u < pmin(0, a))
      if (any(acc)) {
        which_acc <- idx_rw[acc]
        phi[which_acc, ]    <- prop[acc, , drop = FALSE]
        loglik[which_acc]   <- ll_prop[acc]
        logprior[which_acc] <- lp_prop[acc]
        accepted <- accepted + sum(acc)
        nacc_rw <- nacc_rw + sum(acc)
      }
      nprop_rw <- nprop_rw + length(idx_rw)
    }
  }
  total_props <- N * n_moves
  list(phi = phi, loglik = loglik, logprior = logprior, b_state = b_state,
       acc_rate = accepted / max(total_props, 1L),
       rw_accept_rate = if (nprop_rw>0) nacc_rw/nprop_rw else NA_real_,
       indep_accept_rate = if (nprop_id>0) nacc_id/nprop_id else NA_real_)
}

# ------------------------------- Outer SMC -------------------------------
outer_smc_phi_batch <- function(
    caches,                                  # list of subject caches (length S)
    rprior_phi,                              # function(n) -> matrix(n x dphi)
    logprior_phi,                            # function(phi_row) -> scalar log p(φ)
    # Either provide the generic θ|φ evaluator:
    log_prior_theta_given_phi_mat = NULL,    # function(Theta, phi, aux) -> length M_b
    # Or provide the Gaussian fast-path mapper:
    gaussian_map_fn = NULL,                  # function(phi, dθ) -> list(mu, Sigma_inv, logdet, const)
    # --- NEW: needed for enrichment / cache rebuilds ---
    data_list = NULL,                        # list of per-subject data (length S)
    loglik_fn = NULL,                        # function(Theta_block, data_i) -> vector log p(y_i|θ)
    n_cores_inner = 1L,                      # parallelism for inner loglik recompute
    M = 2000L,
    cess_target = 0.95,
    resample_threshold = 0.5,
    resampling = c("systematic","multinomial"),
    n_moves = 2L,
    rw_scale_init = 1.3,
    pm_mode = c("strict", "fast"),
    # Independence mixture kernel (new)
    indep_prob = 0.30,
    indep_t_df = 4L,
    indep_t_prob = 0.85,
    G_mix = 12L,
    elite_quantile = 0.40,
    cov_inflation = 4.0,
    base_mix_eps = 0.10,
    # History mixture controls
    hist_mix_enable = TRUE,
    hist_mix_lambda_thresh = 0.50,
    hist_mix_prob = 0.20,
    max_rounds = 200L,
    # Correlated PM batch-index refresh knobs
    refresh_batches_after_resample = FALSE,
    refresh_batches_each_round = FALSE,   # ### PATCH: default off to preserve CRNs
    refresh_batches_frac = 0.02,          # ### PATCH: gentler refresh if enabled
    block_refresh_every = 5L,
    block_refresh_frac = 0.10,
    # Diagnostics
    diag_enable = TRUE,
    diag_probe = c("wmean","best"),
    diag_replicates = 3L,                    # B replicates over batches (no model calls)
    offenders_print_cap = 10L,               # how many worst offenders to print
    # Targeted enrichment knobs (Fixes #2–#5)
    auto_enrich_enable = TRUE,
    auto_enrich_every = 10L,
    auto_enrich_lambda_gate = 0.25,
    auto_enrich_ess_thresh = 0.20,
    auto_enrich_k_thresh   = 0.70,
    auto_enrich_max_units  = 2L,     # ### PATCH: safer default
    eps_prior_anchor = 0.03,         # ### PATCH: base; will adapt per-batch
    weak_inflate_factor = 3.0,       # ### PATCH: toned down
    # Seeding
    seed = 123,
    verbose = TRUE
) {
  resampling <- match.arg(resampling)
  diag_probe <- match.arg(diag_probe)
  pm_mode <- match.arg(pm_mode)
  set.seed(seed)

  S <- length(caches)
  if (S <= 0) stop("Empty 'caches' list.")
  if (is.null(gaussian_map_fn) && is.null(log_prior_theta_given_phi_mat)) {
    stop("Provide either gaussian_map_fn or log_prior_theta_given_phi_mat.")
  }
  # --- NEW: enrichment input sanity ---
  if (auto_enrich_enable) {
    if (is.null(data_list) || is.null(loglik_fn)) {
      message("[Enrich] auto_enrich_enable=TRUE but data_list/loglik_fn not provided; disabling enrichment.")
      auto_enrich_enable <- FALSE
    } else if (length(data_list) != S) {
      stop("[Enrich] length(data_list) must equal number of subjects S = ", S)
    }
  }
  # Self-register Gaussian fast-path onto caches if provided
  if (!is.null(gaussian_map_fn)) {
    for (i in seq_along(caches)) {
      if (is.null(caches[[i]]$gaussian_map)) {
        d_i <- ncol(caches[[i]]$batches[[1]]$Theta)
        caches[[i]] <- register_gaussian_prior_map(
          caches[[i]],
          function(phi) gaussian_map_fn(phi, d_i)
        )
      }
    }
  }

  # 0) Initialize auxiliary batch state
  # fast: vector length S, strict: matrix M x S
  b_state <- .init_aux_batch_state(caches, M = M, pm_mode = pm_mode)

  # 1) Initialize particles from prior
  phi <- rprior_phi(M)               # M x dphi
  if (!is.matrix(phi)) phi <- matrix(phi, nrow = M)
  dphi <- ncol(phi)
  logprior <- apply(phi, 1L, logprior_phi)

  # Initial loglik (batch-indexed)
  loglik <- .outer_loglik_hat_particles(
    Phi = phi, caches = caches, b_state = b_state,
    log_prior_theta_given_phi_mat = log_prior_theta_given_phi_mat,
    gaussian_map_fn = gaussian_map_fn,
    pm_mode = pm_mode
  )

  # ### PATCH: sanity message for per-batch size on weakly ID params
  M_per_hint <- vapply(caches, function(ci) ci$batches[[1]]$M, numeric(1))
  if (verbose && any(M_per_hint < 128)) {
    message(sprintf("[Hint] Some batches are small (min M_per=%d). For weakly identified params, prefer M_per >= 128 (increase M or reduce K).",
                    min(M_per_hint)))
  }

  # Finite checks
  if (!all(is.finite(logprior))) stop("Non-finite log prior at initialization.")
  if (!all(is.finite(loglik)))   stop("Non-finite log inner likelihood at initialization.")

  lambda <- 0.0
  w <- rep(1/M, M)
  logZ <- 0.0
  logZ_increments <- numeric()   # record per-round logZ increments
  mcse_var_accum <- 0.0

  # Adaptive RW scale
  target_acc <- 0.234
  log_rw_scale <- log(rw_scale_init)
  rm_gain <- 0.05

  lambda_hist <- lambda
  acc_hist <- numeric()
  # mixture state (φ-space)
  last_elite_mix_phi <- NULL
  elite_history_phi <- list()

  # Diagnostics store
  diag_list <- list()
  aux_refresh_acc_hist <- numeric()

  round <- 0L
  # rolling acceptances for adaptive n_moves
  last_rw_acc <- NA_real_
  last_id_acc <- NA_real_
  moves_cap <- 8L
  while (lambda < 1 - 1e-12 && round < max_rounds) {
    round <- round + 1L
    if (verbose) cat(sprintf("\n[Round %d] λ=%.3f  ", round, lambda))

    # 2) Choose δ via rCESS (from smc_core)
    h <- loglik
    target_cess <- cess_target_at_lambda(lambda)
    delta <- next_lambda_via_rCESS_stat(w, h, lambda, target = target_cess) - lambda
    if (delta <= 1e-8) delta <- min(1 - lambda, 1e-4)
    lambda_new <- min(1.0, lambda + delta)

    # 3) Incremental weighting + evidence
    m <- max(h)
    log_u <- delta * (h - m)
    logZ_inc <- logsumexp_w(log_u, w) + delta * m
    logZ <- logZ + logZ_inc
    logZ_increments <- c(logZ_increments, logZ_inc)

    # Delta-method MCSE accumulation for log evidence
    u <- exp(log_u)
    mu1 <- sum(w * u)
    mu2 <- sum(w * u * u)
    Neff <- 1 / sum(w * w)
    var_logZ_inc <- (mu2 - mu1^2) / (max(Neff, 1) * max(mu1^2, .Machine$double.eps))
    mcse_var_accum <- mcse_var_accum + max(var_logZ_inc, 0)

    w <- w * u
    w <- w / sum(w)
    ess_frac <- 1 / sum(w * w) / M
    if (verbose) cat(sprintf("Δ=%.4f | ESS/N=%.3f | logZ+=%.4f\n", delta, ess_frac, logZ_inc))

    # 4) Resample if needed
    resampled <- FALSE
    if (ess_frac < resample_threshold) {
      # 5) Pre-resample light rejuvenation
      if (ess_frac > 0.25 && lambda >= 0.20) {
        if (verbose) cat("  Pre-resample light rejuvenation (1 sweep)...\n")
        # light mixtures for pre-resample
        elite_mix_phi_pre <- tryCatch({
          mix <- fit_elite_mixture_Z(phi, w, elite_quantile = elite_quantile, G = G_mix,
                                     cov_inflation = cov_inflation, housekeeping = TRUE,
                                     min_G_keep = 2, merge_thresh = 0.10,
                                     verbose = FALSE, lambda = lambda_new)
          if (!.is_empty_mix(mix)) blend_mixes(mix, .default_std_normal_mix(ncol(phi)), eps = base_mix_eps) else mix
        }, error = function(e) NULL)
        hist_mix_pre <- NULL
        if (hist_mix_enable) {
          hist_mix_pre <- tryCatch(build_hist_mixture(elite_history_phi, lambda_new,
                                                     hist_mix_lambda_thresh, hist_mix_prob,
                                                     min_G_keep = 2, merge_thresh = 0.10,
                                                     min_eig = 3e-3, force = FALSE), error = function(e) NULL)
          if (!is.null(hist_mix_pre) && !is.null(hist_mix_pre$mix) && !.is_empty_mix(hist_mix_pre$mix)) {
            hist_mix_pre$mix <- blend_mixes(hist_mix_pre$mix, .default_std_normal_mix(ncol(phi)), eps = 0.05)
          }
        }
        # weak dims for RW expansion in this sweep
        weak_idx_pre <- weak_dims_from_mat(phi, w, frac = 0.25, min_keep = 1L)
        pre_move <- .rejuvenate_rw_batch(phi, loglik, logprior, w, lambda_new,
                                         caches, b_state, logprior_phi,
                                         log_prior_theta_given_phi_mat, gaussian_map_fn,
                                         pm_mode = pm_mode,
                                         n_moves = 1L, rw_scale = 0.35,
                                         weak_dim_idx = weak_idx_pre, rw_expand_factor = 1.8,
                                         elite_mix = elite_mix_phi_pre,
                                         hist_mix = hist_mix_pre,
                                         indep_prob = min(0.35, indep_prob),
                                         indep_t_df = indep_t_df, indep_t_prob = indep_t_prob,
                                         rng_seed = seed + 19L * round)
        phi      <- pre_move$phi
        loglik   <- pre_move$loglik
        logprior <- pre_move$logprior
        b_state  <- pre_move$b_state
      }
      idx <- if (resampling == "multinomial") {
        .resample_multinomial(w)
      } else {
        # Hilbert-sorted stratified resampling in whitened φ-space
        set.seed(seed + 991 * round)
        # whiten φ → Zφ using weighted mean/cov
        mu_phi <- colSums(phi * w)
        Xc <- sweep(phi, 2L, mu_phi, `-`)
        Sphi <- weighted_cov(phi, w)
        Sphi <- as.matrix(Matrix::nearPD(Sphi, conv.tol = 1e-7)$mat)
        diag(Sphi) <- pmax(diag(Sphi), 1e-10)
        R <- tryCatch(chol(Sphi), error = function(e) chol(Sphi + diag(1e-8, ncol(Sphi))))
        Zphi <- t(backsolve(R, t(Xc), transpose = TRUE))
        ord_h <- hilbert_sort_order(Zphi, bits = 16L)
        ws <- w[ord_h]
        sel_sorted <- stratified_resample_sorted(ws, deterministic = FALSE)
        ord_inv <- ord_h
        idx <- ord_inv[sel_sorted]
      }
      phi      <- phi[idx, , drop = FALSE]
      loglik   <- loglik[idx]
      logprior <- logprior[idx]
      if (is.matrix(b_state)) b_state <- b_state[idx, , drop = FALSE]
      w <- rep(1/M, M)
      resampled <- TRUE
      if (verbose) cat(sprintf("  Resampled (%s).\n", resampling))

      # Optional: refresh a small fraction of batch indices after resampling
      if (refresh_batches_after_resample) {
        if (is.matrix(b_state)) {
          ref <- .refresh_aux_state_block_mh(
            phi = phi, loglik = loglik, b_state = b_state, lambda = lambda_new,
            caches = caches,
            log_prior_theta_given_phi_mat = log_prior_theta_given_phi_mat,
            gaussian_map_fn = gaussian_map_fn,
            frac = refresh_batches_frac,
            rng_seed = seed + 1229L * round
          )
          b_state <- ref$b_state
          loglik <- ref$loglik
          aux_refresh_acc_hist <- c(aux_refresh_acc_hist, ref$acc_rate)
        } else {
          K_all <- caches[[1]]$K
          b_state <- refresh_batch_indices(b_state, K = K_all, frac = refresh_batches_frac)
          loglik <- .outer_loglik_hat_particles(
            Phi = phi, caches = caches, b_state = b_state,
            log_prior_theta_given_phi_mat = log_prior_theta_given_phi_mat,
            gaussian_map_fn = gaussian_map_fn,
            pm_mode = pm_mode
          )
        }
      }
    }

    # 5) Diagnostics at a probe φ (batch-aware; no model calls)
    if (diag_enable) {
      probe_phi <- .choose_probe_phi(phi, w, logprior, loglik, lambda_new, diag_probe)
      probe_b <- .probe_aux_row(phi, w, logprior, loglik, lambda_new, b_state)
      is_df <- .is_diagnostics_batch(caches, probe_phi, probe_b, log_prior_theta_given_phi_mat)
      ess_norm <- is_df$ess / is_df$M
      psis_k <- is_df$psis_k
      cat(sprintf("  [Diag] median ESS_i/M=%.2f | max k=%.2f | any ESS_i/M<0.10: %s\n",
                  stats::median(ess_norm, na.rm = TRUE),
                  if (all(is.na(psis_k))) NA_real_ else max(psis_k, na.rm = TRUE),
                  any(ess_norm < 0.10, na.rm = TRUE)))
      if (diag_replicates > 0L) {
        rsd <- .rep_logZi_sd_over_batches(caches, probe_phi, B_reps = diag_replicates,
                                          log_prior_theta_given_phi_mat = log_prior_theta_given_phi_mat,
                                          gaussian_map_fn = gaussian_map_fn,
                                          rng_seed = seed + 1777 * round)
        cat(sprintf("  [Diag] median SD_over_batches(logZ_i): %.3f\n",
                    stats::median(rsd, na.rm = TRUE)))
        # --- NEW: variance budget log (Fix #4)
        sum_sd2 <- sum(rsd[is.finite(rsd)]^2)
        cat(sprintf("  [VarBudget] sum_i SD_over_batches(logZ_i)^2 = %.3f (target ≈ 1–2)\n", sum_sd2))
      }
      # --- NEW: print top offenders (Fix #4)
      bad_idx <- which( (ess_norm < auto_enrich_ess_thresh) |
                        (!is.na(psis_k) & (psis_k > auto_enrich_k_thresh)) )
      if (length(bad_idx)) {
        ord <- order(ifelse(is.na(psis_k[bad_idx]), -Inf, psis_k[bad_idx]),
                     ess_norm[bad_idx], decreasing = c(TRUE, FALSE))
        show <- bad_idx[ord][seq_len(min(length(bad_idx), offenders_print_cap))]
        msg <- paste(sprintf("#%d ESS/M=%.2f k=%.2f", show, ess_norm[show], psis_k[show]), collapse = " | ")
        cat("  [Diag] Offenders:", msg, "\n")
      }
    }

    # Build elite mixture on φ for independence proposals (with tiny base mass)
    elite_mix_phi <- tryCatch({
      mix <- fit_elite_mixture_Z(phi, w, elite_quantile = elite_quantile, G = G_mix,
                                 cov_inflation = cov_inflation, housekeeping = TRUE,
                                 min_G_keep = 2, merge_thresh = 0.10,
                                 verbose = FALSE, lambda = lambda_new)
      if (!.is_empty_mix(mix)) blend_mixes(mix, .default_std_normal_mix(ncol(phi)), eps = base_mix_eps) else mix
    }, error = function(e) NULL)
    hist_mix <- NULL
    if (hist_mix_enable) {
      hist_mix <- tryCatch(build_hist_mixture(elite_history_phi, lambda_new,
                                              hist_mix_lambda_thresh, hist_mix_prob,
                                              min_G_keep = 2, merge_thresh = 0.10,
                                              min_eig = 3e-3, force = FALSE), error = function(e) NULL)
      if (!is.null(hist_mix) && !is.null(hist_mix$mix) && !.is_empty_mix(hist_mix$mix)) {
        # light base blend for history as well
        hist_mix$mix <- blend_mixes(hist_mix$mix, .default_std_normal_mix(ncol(phi)), eps = 0.05)
      }
    }

    # Adaptive number of moves (Δλ, difficulty, acceptance-aware)
    pred_rcess_full <- rCESS_stat(w, h, 1 - lambda)
    difficulty <- 1 - pred_rcess_full
    sDelta <- .clamp(delta / 0.05, 0.6, 1.4)
    sDiff  <- 1 + .clamp(1.5 * difficulty, 0, 1.5)
    acc_proxy <- mean(c(last_rw_acc, last_id_acc), na.rm = TRUE)
    if (!is.finite(acc_proxy)) acc_proxy <- 0.25
    sAcc <- if (acc_proxy < 0.10) 0.80 else if (acc_proxy < 0.20) 0.90 else 1.00
    base_moves <- n_moves
    floor_moves <- if (ess_frac < 0.40) 4L else if (ess_frac < 0.70) 3L else 2L
    n_moves_eff <- as.integer(ceiling(base_moves * sDelta * sDiff * sAcc))
    n_moves_eff <- max(floor_moves, min(n_moves_eff, moves_cap))
    if (lambda > 0.90) n_moves_eff <- max(n_moves_eff, 3L)
    if (lambda < 0.20) n_moves_eff <- min(n_moves_eff, 2L)

    # Weak-dimension inflation indices in φ-space
    weak_idx <- weak_dims_from_mat(phi, w, frac = 0.25, min_keep = 1L)

    # 6) Rejuvenate (MH RW + independence; correlated PM via fixed b_idx inside the sweep)
    rw_scale <- exp(log_rw_scale)
    move <- .rejuvenate_rw_batch(phi, loglik, logprior, w, lambda_new,
                                 caches, b_state, logprior_phi,
                                 log_prior_theta_given_phi_mat, gaussian_map_fn,
                                 pm_mode = pm_mode,
                                 n_moves = n_moves_eff,
                                 rw_scale = rw_scale,
                                 weak_dim_idx = weak_idx,
                                 rw_expand_factor = 2.5,
                                 elite_mix = elite_mix_phi,
                                 hist_mix = hist_mix,
                                 indep_prob = indep_prob,
                                 indep_t_df = indep_t_df,
                                 indep_t_prob = indep_t_prob,
                                 rng_seed = seed + 97L * round)
    phi      <- move$phi
    loglik   <- move$loglik
    logprior <- move$logprior
    b_state  <- move$b_state
    acc_hist <- c(acc_hist, move$acc_rate)
    # Adapt RW scale
    log_rw_scale <- .clamp(log_rw_scale + rm_gain * (move$acc_rate - target_acc),
                           log(0.05), log(2.5))
    if (verbose) {
      cat(sprintf("  MH acc=%.3f | rw=%.3f id=%.3f | rw_scale=%.3f\n",
                  move$acc_rate, move$rw_accept_rate %||% NA_real_,
                  move$indep_accept_rate %||% NA_real_, exp(log_rw_scale)))
    }

    last_rw_acc <- move$rw_accept_rate
    last_id_acc <- move$indep_accept_rate

    # Strict-mode auxiliary-state block refresh (MH-correct at fixed λ)
    if (is.matrix(b_state) && block_refresh_every > 0L &&
        (round %% as.integer(max(1L, block_refresh_every)) == 0L)) {
      ref <- .refresh_aux_state_block_mh(
        phi = phi, loglik = loglik, b_state = b_state, lambda = lambda_new,
        caches = caches,
        log_prior_theta_given_phi_mat = log_prior_theta_given_phi_mat,
        gaussian_map_fn = gaussian_map_fn,
        frac = block_refresh_frac,
        rng_seed = seed + 7001L * round
      )
      b_state <- ref$b_state
      loglik <- ref$loglik
      aux_refresh_acc_hist <- c(aux_refresh_acc_hist, ref$acc_rate)
      if (verbose) cat(sprintf("  Aux refresh acc=%.3f\n", ref$acc_rate))
    }

    # update history store (keep last 6)
    if (!is.null(elite_mix_phi) && !.is_empty_mix(elite_mix_phi)) last_elite_mix_phi <- elite_mix_phi
    if (lambda_new > 0.5 && !is.null(last_elite_mix_phi) && !.is_empty_mix(last_elite_mix_phi)) {
      elite_history_phi[[length(elite_history_phi) + 1L]] <- last_elite_mix_phi
      if (length(elite_history_phi) > 6L) elite_history_phi <- tail(elite_history_phi, 6L)
    }

    # 7) NEW: targeted enrichment of worst subjects (Fixes #2–#5)
    if (auto_enrich_enable && lambda_new >= auto_enrich_lambda_gate && (round %% auto_enrich_every == 0)) {
      # use the same probe φ as diagnostics (wmean/best), and build a small panel
      probe_phi <- .choose_probe_phi(phi, w, logprior, loglik, lambda_new, diag_probe)
      probe_phi_wmean <- .choose_probe_phi(phi, w, logprior, loglik, lambda_new, "wmean")
      probe_phi_best  <- .choose_probe_phi(phi, w, logprior, loglik, lambda_new, "best")
      probe_b <- .probe_aux_row(phi, w, logprior, loglik, lambda_new, b_state)
      is_df <- .is_diagnostics_batch(caches, probe_phi, probe_b, log_prior_theta_given_phi_mat)
      ess_norm <- is_df$ess / is_df$M
      psis_k <- is_df$psis_k
      bad_idx <- which( (ess_norm < auto_enrich_ess_thresh) |
                        (!is.na(psis_k) & (psis_k > auto_enrich_k_thresh)) )
      if (length(bad_idx)) {
        ord <- order(ifelse(is.na(psis_k[bad_idx]), -Inf, psis_k[bad_idx]),
                     ess_norm[bad_idx], decreasing = c(TRUE, FALSE))
        select <- bad_idx[ord][seq_len(min(length(bad_idx), auto_enrich_max_units))]
        cat("  [Enrich] subjects:", paste(select, collapse = ","), "\n")
        for (ii in select) {
          if (is.null(data_list[[ii]])) {
            warning(sprintf("[Enrich] data_list[[%d]] is NULL; skipping subject %d", ii, ii))
            next
          }
          # ### PATCH: surgical enrichment on ACTIVE batch only with gating/rollback
          caches[[ii]] <- enrich_subject_cache_with_anchor_surgical(
            cache = caches[[ii]],
            phi_probe = probe_phi,
            b = probe_b[ii],
            log_prior_theta_given_phi_mat = log_prior_theta_given_phi_mat,
            gaussian_map_fn = if (!is.null(gaussian_map_fn)) function(phi, d) gaussian_map_fn(phi, d) else NULL,
            data = data_list[[ii]],
            loglik_fn = loglik_fn,
            eps_prior_base = eps_prior_anchor,
            weak_inflate_factor = weak_inflate_factor,
            elite_q = 0.4, G = 12L, cov_infl = 4.0,
            replace_frac = NULL,                      # adaptive replacement fraction
            sobol_seed = caches[[ii]]$batches[[probe_b[ii]]]$sobol_seed + 1L,  # small shift only
            n_cores = n_cores_inner,
            verbose = verbose,
            subj_id = ii,
            phi_panel = rbind(probe_phi_wmean, probe_phi_best)
          )
        }
        # Aux state unchanged; proposal correlation preserved.
      }
    }

    loglik <- .outer_loglik_hat_particles(
      Phi = phi, caches = caches, b_state = b_state,
      log_prior_theta_given_phi_mat = log_prior_theta_given_phi_mat,
      gaussian_map_fn = gaussian_map_fn,
      pm_mode = pm_mode
    )

    # 8) Optional light batch-index refresh at the end of the round
    if (refresh_batches_each_round) {
      if (is.matrix(b_state)) {
        ref <- .refresh_aux_state_block_mh(
          phi = phi, loglik = loglik, b_state = b_state, lambda = lambda_new,
          caches = caches,
          log_prior_theta_given_phi_mat = log_prior_theta_given_phi_mat,
          gaussian_map_fn = gaussian_map_fn,
          frac = refresh_batches_frac,
          rng_seed = seed + 9929L * round
        )
        b_state <- ref$b_state
        loglik <- ref$loglik
        aux_refresh_acc_hist <- c(aux_refresh_acc_hist, ref$acc_rate)
      } else {
        K_all <- caches[[1]]$K
        b_state <- refresh_batch_indices(b_state, K = K_all, frac = refresh_batches_frac)
        loglik <- .outer_loglik_hat_particles(
          Phi = phi, caches = caches, b_state = b_state,
          log_prior_theta_given_phi_mat = log_prior_theta_given_phi_mat,
          gaussian_map_fn = gaussian_map_fn,
          pm_mode = pm_mode
        )
      }
    }

    # Advance annealing
    lambda <- lambda_new
    lambda_hist <- c(lambda_hist, lambda)
  }

  list(
    phi = phi, w = w,
    loglik = loglik, logprior = logprior,
    log_evidence = logZ,
    mcse_log_evidence = sqrt(mcse_var_accum),
    logZ_increments = logZ_increments,
    b_idx = if (!is.matrix(b_state)) b_state else NULL,
    B_idx = if (is.matrix(b_state)) b_state else NULL,
    meta = list(
      pm_mode = pm_mode,
      lambda_hist = lambda_hist,
      acc_hist = acc_hist,
      rw_scale_final = exp(log_rw_scale),
      resampling = resampling,
      refresh_batches_after_resample = refresh_batches_after_resample,
      refresh_batches_each_round = refresh_batches_each_round,
      refresh_batches_frac = refresh_batches_frac,
      block_refresh_every = block_refresh_every,
      block_refresh_frac = block_refresh_frac,
      aux_refresh_acc_hist = aux_refresh_acc_hist,
      n_cores_inner = n_cores_inner
    )
  )
}

# ---------------------- Generic outer SMC (direct likelihood) -------------------
.outer_eval_logprior_one <- function(theta_row, logprior_fn) {
  val <- tryCatch(logprior_fn(theta_row), error = function(e) NULL)
  if (is.null(val)) {
    val <- tryCatch(
      logprior_fn(matrix(theta_row, nrow = 1L)),
      error = function(e) stop("logprior failed for both vector and 1-row matrix input.")
    )
  }
  val <- as.numeric(val)
  if (!length(val)) stop("logprior returned empty output for one-row input.")
  val[1L]
}

.outer_eval_logprior_mat <- function(Theta, logprior_fn) {
  out <- tryCatch(logprior_fn(Theta), error = function(e) NULL)
  if (!is.null(out)) {
    out <- as.numeric(out)
    if (length(out) == nrow(Theta)) return(out)
  }
  vapply(seq_len(nrow(Theta)), function(i) .outer_eval_logprior_one(Theta[i, , drop = TRUE], logprior_fn), numeric(1))
}

outer_smc_direct <- function(
    data,
    rprior,                  # function(n) -> matrix(n x d)
    logprior,                # function(theta_matrix_or_row) -> vector/scalar
    loglik_fn,               # function(theta_matrix_or_row, data) -> vector/scalar
    M = 2000L,
    cess_target = 0.95,
    resample_threshold = 0.5,
    n_moves = 2L,
    rw_scale_init = 1.0,
    max_rounds = 200L,
    n_cores_loglik = 1L,
    seed = 123,
    verbose = TRUE
) {
  set.seed(seed)

  Theta <- rprior(M)
  if (!is.matrix(Theta)) Theta <- matrix(Theta, nrow = M)
  d <- ncol(Theta)

  logprior_curr <- .outer_eval_logprior_mat(Theta, logprior)
  loglik_curr <- ll_parallel(Theta, data, loglik_fn, n_cores = n_cores_loglik)
  if (!all(is.finite(logprior_curr))) stop("Non-finite log prior at initialization.")
  if (!any(is.finite(loglik_curr))) stop("All initial log-likelihood values are non-finite.")
  loglik_curr[!is.finite(loglik_curr)] <- -Inf

  lambda <- 0.0
  w <- rep(1 / M, M)
  logZ <- 0.0

  lambda_hist <- lambda
  ess_hist <- numeric()
  acc_hist <- numeric()
  logZ_increments <- numeric()
  log_rw_scale <- log(rw_scale_init)

  round <- 0L
  while (lambda < 1 - 1e-12 && round < max_rounds) {
    round <- round + 1L
    lambda_new <- next_lambda_via_rCESS_stat(w, loglik_curr, lambda, target = cess_target)
    delta <- lambda_new - lambda
    if (!is.finite(delta) || delta <= 1e-8) {
      delta <- min(1 - lambda, 1e-4)
      lambda_new <- lambda + delta
    }

    m <- max(loglik_curr)
    log_u <- delta * (loglik_curr - m)
    logZ_inc <- logsumexp_w(log_u, w) + delta * m
    if (!is.finite(logZ_inc)) stop("Non-finite log-evidence increment.")
    logZ <- logZ + logZ_inc
    logZ_increments <- c(logZ_increments, logZ_inc)

    w <- w * exp(log_u)
    sw <- sum(w)
    if (!is.finite(sw) || sw <= 0) stop("Particle weights became invalid.")
    w <- w / sw
    ess_frac <- (1 / sum(w * w)) / M
    ess_hist <- c(ess_hist, ess_frac)

    if (verbose) {
      cat(sprintf("[Outer Round %d] lambda %.4f -> %.4f | ESS/N=%.3f | logZ+=%.4f\n",
                  round, lambda, lambda_new, ess_frac, logZ_inc))
    }

    if (ess_frac < resample_threshold) {
      idx <- stratified_resample_sorted(w, deterministic = FALSE)
      Theta <- Theta[idx, , drop = FALSE]
      logprior_curr <- logprior_curr[idx]
      loglik_curr <- loglik_curr[idx]
      w <- rep(1 / M, M)
    }

    rw_scale <- exp(log_rw_scale)
    S <- tryCatch(weighted_cov(Theta, w), error = function(e) stats::cov(Theta))
    S <- as.matrix(S)
    diag(S) <- pmax(diag(S), 1e-8)
    L <- tryCatch(
      chol((rw_scale^2 / max(d, 1L)) * S + diag(1e-10, d)),
      error = function(e) chol((rw_scale^2 / max(d, 1L)) * (S + diag(1e-6, d)))
    )

    acc_total <- 0L
    for (mv in seq_len(n_moves)) {
      prop <- Theta + matrix(rnorm(M * d), M, d) %*% L
      lp_prop <- .outer_eval_logprior_mat(prop, logprior)
      ll_prop <- ll_parallel(prop, data, loglik_fn, n_cores = n_cores_loglik)
      loga <- (lp_prop - logprior_curr) + lambda_new * (ll_prop - loglik_curr)
      loga[!is.finite(loga)] <- -Inf
      acc <- log(runif(M)) < pmin(0, loga)
      if (any(acc)) {
        Theta[acc, ] <- prop[acc, , drop = FALSE]
        logprior_curr[acc] <- lp_prop[acc]
        loglik_curr[acc] <- ll_prop[acc]
      }
      acc_total <- acc_total + sum(acc)
    }

    acc_rate <- acc_total / max(1L, M * n_moves)
    acc_hist <- c(acc_hist, acc_rate)
    log_rw_scale <- .clamp(log_rw_scale + 0.05 * (acc_rate - 0.234), log(0.02), log(3.0))

    lambda <- lambda_new
    lambda_hist <- c(lambda_hist, lambda)
  }

  list(
    theta = Theta,
    w = w,
    loglik = loglik_curr,
    logprior = logprior_curr,
    log_evidence = logZ,
    logZ_increments = logZ_increments,
    meta = list(
      lambda_hist = lambda_hist,
      ess_hist = ess_hist,
      acc_hist = acc_hist,
      rw_scale_final = exp(log_rw_scale),
      rounds = round,
      final_lambda = lambda
    )
  )
}
