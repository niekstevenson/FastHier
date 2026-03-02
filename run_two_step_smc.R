#!/usr/bin/env Rscript

# Consolidated maintained modes only: E4 (fixed-cache baseline) and
# H1 (anchor-and-refine strict adaptive PM).
single_experiment_id <- "H1"
publish_experiment_id <- "H1"
publish_seed <- 123L

env_chr <- function(name, default) {
  v <- Sys.getenv(name, unset = "")
  if (!nzchar(v)) default else v
}
env_int <- function(name, default) {
  v <- Sys.getenv(name, unset = "")
  if (!nzchar(v)) default else as.integer(v)
}
env_int_vec <- function(name, default) {
  v <- Sys.getenv(name, unset = "")
  if (!nzchar(v)) return(default)
  parts <- strsplit(v, ",", fixed = TRUE)[[1L]]
  parts <- trimws(parts)
  as.integer(parts[nzchar(parts)])
}

single_experiment_id <- env_chr("FH_SINGLE_EXPERIMENT_ID", single_experiment_id)
publish_experiment_id <- env_chr("FH_PUBLISH_EXPERIMENT_ID", publish_experiment_id)
publish_seed <- env_int("FH_PUBLISH_SEED", publish_seed)
seed_override <- env_int_vec("FH_SEEDS", integer(0))

suppressPackageStartupMessages({
  library(parallel)
})

source("smc_core.R")
source("SMC_super_fast.R")
source("new_SMC_cache.R")
source("outer_SMC.R")
source("make_prior.R")

stan_path <- file.path("samples", "stan_full_results.rds")
if (!file.exists(stan_path)) {
  stop("Missing Stan results at ", stan_path, ". Run run_stan.R first.")
}
stan_res <- readRDS(stan_path)

# ----------------------------- Configuration ------------------------------
base_config <- list(
  seeds = if (length(seed_override)) seed_override else 123:125,
  pilot_inner = list(
    M = 450L,
    max_rounds = 120L,
    n_moves = 1L,
    rw_scale_init = 1.1
  ),
  cache_strategy = list(
    M = 256L,
    K_batches = 4L,
    min_K_batches = 1L,
    min_M_per = 0L,
    preserve_K_with_minM = FALSE,
    hard_cache_multiplier = 1.0,
    deterministic_counts = FALSE,
    cache_bridge_enable = FALSE,
    cache_bridge_weight = 0.20,
    cache_defensive_weight = 0.00,
    cache_defensive_df = 3L,
    defensive_t_eps = 0.25,
    defensive_t_df = 3L
  ),
  quality_gate = list(
    enable = FALSE,
    ess_abs_threshold = 8.0,
    sd_logZ_threshold = 1.0,
    target_sum_sd2 = 15.0
  ),
  enrichment = list(
    max_passes = 2L,
    max_subjects_per_pass = 5L,
    m_add_per_subject = 256L,
    eps_prior_base = 0.04,
    weak_inflate_factor = 3.0,
    min_sum_sd2_improve = 0.25
  ),
  outer_escalation = list(
    enable = FALSE,
    M_subject_default = 300L,
    M_subject_hard = 1200L,
    n_mcmc_moves_default = 2L,
    n_mcmc_moves_hard = 4L,
    max_hard_subjects = 5L,
    near_boundary_ratio = 0.92,
    round_quantile = 0.75,
    mcse_quantile = 0.75,
    shift_sd_quantile = 0.75
  ),
  pilot_strategy = list(
    enable = TRUE,
    subset_enable = TRUE,
    subset_threshold_S = 60L,
    subset_size_cap = 60L,
    subset_frac = 0.20,
    seed = 123L,
    M_subject_pilot = 120L,
    n_mcmc_moves_pilot = 1L,
    max_rounds_pilot = 80L
  ),
  anchor_estimation = list(
    method = "eb_robust",
    winsor_probs = c(0.05, 0.95),
    sigma2_floor = 1e-3,
    sigma2_ceiling = 2.5,
    outer_ref_inflation = 1.5
  ),
  outer_refinement = list(
    max_retries = 2L,
    strict_require_lambda1 = TRUE,
    lambda_subject_threshold = 0.999,
    lambda_strict_threshold = 0.9999,
    max_rounds_default = 220L,
    max_rounds_hard = 350L
  ),
  pm_quality_gate = list(
    M_probe = 512L,
    ess_abs_threshold = 12.0,
    var_proxy_threshold = 0.60,
    var_proxy_q90_threshold = 0.80,
    khat_report_threshold = 2.5,
    M_ladder = c(256L, 512L, 1024L, 2048L),
    max_passes = 2L,
    enable_mixture_fallback = TRUE,
    fallback_w_local = 0.85,
    fallback_w_phi_anchor = 0.10,
    fallback_w_defensive = 0.05,
    rerun_outer_on_fail = TRUE,
    outer_rerun_max_subjects = 8L,
    outer_rerun_M_default = 1200L,
    outer_rerun_M_hard = 2400L,
    outer_rerun_n_moves_default = 4L,
    outer_rerun_n_moves_hard = 6L
  ),
  inner_strict = list(
    adapt_until_round = 15L,
    freeze_after_round = 15L,
    completion_mode = "lambda1",
    strict_require_lambda1 = TRUE,
    lambda_strict_threshold = 0.9999
  ),
  inner = list(
    M = 1200L,
    cess_target = 0.95,
    resample_threshold = 0.5,
    n_moves = 2L,
    rw_scale_init = 1.2,
    block_refresh_every = 5L,
    block_refresh_frac = 0.10,
    max_rounds = 200L,
    inner_ll_mode = "fixed_cache",
    pm_aux_mode = "batch_idx"
  ),
  adaptive_pm = list(
    M_default = 128L,
    M_levels = c(64L, 128L, 256L, 512L),
    ess_norm_threshold = 0.15,
    khat_threshold = 0.80,
    var_log_target_low = 0.30,
    var_log_target_high = 1.00,
    var_relax_ess = 0.35,
    decrease_enabled = TRUE,
    M_update_every = 1L,
    da_enable = TRUE,
    da_M = 32L,
    max_retries = 2L,
    min_log_mhat = -1e12,
    completion_mode = "lambda1",
    max_wall_time_sec = Inf,
    w_local = 0.60,
    w_phi_anchor = 0.30,
    w_defensive = 0.10,
    defensive_df = 3L,
    anchor_n = 96L,
    anchor_cov_inflation = 1.5,
    adapt_until_round = 15L,
    freeze_after_round = 15L
  )
)
experiment_grid <- list(
  E4 = list(mode = "fixed_cache"),
  H1 = list(mode = "adaptive_pm")
)

if (!single_experiment_id %in% names(experiment_grid)) {
  stop("FH_SINGLE_EXPERIMENT_ID must be one of: ", paste(names(experiment_grid), collapse = ", "))
}
experiment_grid <- experiment_grid[single_experiment_id]

# ------------------------------- Helpers ----------------------------------
weighted_mean <- function(X, w) {
  w <- pmax(as.numeric(w), 0)
  w <- w / sum(w)
  colSums(as.matrix(X) * w)
}

weighted_sd <- function(x, w) {
  w <- pmax(as.numeric(w), 0)
  w <- w / sum(w)
  m <- sum(w * x)
  sqrt(sum(w * (x - m)^2))
}

estimate_phi_anchor_from_outer <- function(
    outer_subject,
    method = "eb_robust",
    winsor_probs = c(0.05, 0.95),
    sigma2_floor = 1e-3,
    sigma2_ceiling = 2.5,
    eps_var = 1e-6
) {
  S <- length(outer_subject)
  if (S <= 0L) stop("outer_subject is empty.")
  d <- ncol(as.matrix(outer_subject[[1L]]$Theta))
  subj_mean <- matrix(NA_real_, nrow = S, ncol = d)
  subj_var <- matrix(NA_real_, nrow = S, ncol = d)
  for (i in seq_len(S)) {
    Theta <- as.matrix(outer_subject[[i]]$Theta)
    w <- as.numeric(outer_subject[[i]]$w)
    w <- pmax(w, 0)
    sw <- sum(w)
    if (!is.finite(sw) || sw <= 0) w <- rep(1 / nrow(Theta), nrow(Theta)) else w <- w / sw
    mu_i <- colSums(Theta * w)
    vc_i <- colSums((sweep(Theta, 2L, mu_i, `-`)^2) * w)
    subj_mean[i, ] <- mu_i
    subj_var[i, ] <- pmax(vc_i, 0)
  }

  winsorize_vec <- function(x, probs = c(0.05, 0.95)) {
    q <- as.numeric(stats::quantile(x, probs = probs, names = FALSE, na.rm = TRUE, type = 8))
    if (length(q) != 2L || any(!is.finite(q))) return(x)
    pmin(pmax(x, q[1L]), q[2L])
  }

  if (identical(method, "eb_robust")) {
    subj_mean_w <- apply(subj_mean, 2L, winsorize_vec, probs = winsor_probs)
    if (!is.matrix(subj_mean_w)) subj_mean_w <- matrix(subj_mean_w, ncol = d)
    mu_anchor <- colMeans(subj_mean_w, na.rm = TRUE)
    between <- if (S > 1L) apply(subj_mean_w, 2L, stats::var, na.rm = TRUE) else rep(0, d)
    within <- apply(subj_var, 2L, stats::median, na.rm = TRUE)
  } else {
    mu_anchor <- colMeans(subj_mean, na.rm = TRUE)
    between <- if (S > 1L) apply(subj_mean, 2L, stats::var, na.rm = TRUE) else rep(0, d)
    within <- colMeans(subj_var, na.rm = TRUE)
  }

  sigma2_anchor <- pmax(within + between, eps_var)
  sigma2_anchor <- pmax(as.numeric(sigma2_floor), pmin(as.numeric(sigma2_ceiling), sigma2_anchor))
  c(mu_anchor, log(sigma2_anchor))
}

select_pilot_subjects <- function(data_list, pilot_cfg, seed_run) {
  S <- length(data_list)
  idx_all <- seq_len(S)
  if (!isTRUE(pilot_cfg$enable)) return(idx_all)
  if (!isTRUE(pilot_cfg$subset_enable) || S <= as.integer(pilot_cfg$subset_threshold_S)) return(idx_all)

  n_take <- min(
    as.integer(pilot_cfg$subset_size_cap),
    as.integer(ceiling(as.numeric(pilot_cfg$subset_frac) * S))
  )
  n_take <- max(1L, min(S, n_take))
  if (n_take >= S) return(idx_all)

  min_y <- vapply(data_list, min, numeric(1L))
  rough <- vapply(data_list, function(y) {
    sy <- stats::sd(y)
    my <- abs(mean(y)) + 1e-6
    sy / my
  }, numeric(1L))

  q1 <- unique(as.numeric(stats::quantile(min_y, probs = seq(0, 1, length.out = 5), na.rm = TRUE, type = 8)))
  if (length(q1) < 2L) q1 <- c(min(min_y), max(min_y) + 1e-9)
  q2 <- unique(as.numeric(stats::quantile(rough, probs = seq(0, 1, length.out = 5), na.rm = TRUE, type = 8)))
  if (length(q2) < 2L) q2 <- c(min(rough), max(rough) + 1e-9)

  b1 <- as.integer(cut(min_y, breaks = q1, include.lowest = TRUE, labels = FALSE))
  b2 <- as.integer(cut(rough, breaks = q2, include.lowest = TRUE, labels = FALSE))
  strata <- interaction(b1, b2, drop = TRUE, lex.order = TRUE)
  split_idx <- split(idx_all, strata)

  set.seed(as.integer(pilot_cfg$seed + 1009L * seed_run))
  selected <- integer(0)
  # First pass: one per stratum to guarantee coverage.
  for (g in split_idx) {
    if (length(selected) >= n_take) break
    selected <- c(selected, sample(g, 1L))
  }
  selected <- unique(selected)
  if (length(selected) < n_take) {
    rem <- setdiff(idx_all, selected)
    selected <- c(selected, sample(rem, n_take - length(selected)))
  }
  sort(unique(as.integer(selected)))
}

weighted_quantile <- function(x, w, probs = c(0.025, 0.5, 0.975)) {
  x <- as.numeric(x)
  w <- pmax(as.numeric(w), 0)
  w <- w / sum(w)
  o <- order(x)
  x <- x[o]
  w <- w[o]
  cw <- cumsum(w)
  sapply(probs, function(p) x[which(cw >= p)[1L]])
}

resolve_cache_budget <- function(cache_cfg, is_hard = FALSE) {
  mult <- if (is_hard) as.numeric(cache_cfg$hard_cache_multiplier %||% 1.0) else 1.0
  M_i <- as.integer(max(1L, round(cache_cfg$M * mult)))
  K_i <- as.integer(max(cache_cfg$min_K_batches %||% 1L, cache_cfg$K_batches))
  min_per <- as.integer(cache_cfg$min_M_per %||% 0L)
  preserve_K <- isTRUE(cache_cfg$preserve_K_with_minM)

  if (min_per > 0L) {
    if (preserve_K) {
      M_i <- max(M_i, as.integer(K_i * min_per))
    } else {
      if (M_i < min_per) M_i <- min_per
      if (ceiling(M_i / K_i) < min_per) {
        K_i <- max(cache_cfg$min_K_batches %||% 1L, floor(M_i / min_per))
        if (K_i < 1L) K_i <- 1L
      }
    }
  }

  list(M = M_i, K_batches = K_i)
}

loglik_shifted_gamma <- function(Theta, y) {
  Theta <- as.matrix(Theta)
  eps <- 1e-9
  shape <- exp(Theta[, 1L]) + eps
  scale <- exp(Theta[, 2L]) + eps
  shift <- exp(Theta[, 3L]) + eps
  out <- rep(-1e12, nrow(Theta))
  min_y <- min(y)
  ok <- shift < min_y
  if (!any(ok)) return(out)
  for (i in which(ok)) out[i] <- sum(dgamma(y - shift[i], shape = shape[i], scale = scale[i], log = TRUE))
  out[!is.finite(out)] <- -1e12
  out
}

inner_weighted_summary <- function(inner_fit) {
  phi <- as.matrix(inner_fit$phi)
  w <- as.numeric(inner_fit$w)
  w <- w / sum(w)
  list(
    mean = c(
      mu_shape = sum(phi[, 1L] * w),
      mu_scale = sum(phi[, 2L] * w),
      mu_shift = sum(phi[, 3L] * w),
      sigma2_shape = sum(exp(phi[, 4L]) * w),
      sigma2_scale = sum(exp(phi[, 5L]) * w),
      sigma2_shift = sum(exp(phi[, 6L]) * w)
    ),
    median = c(
      mu_shape = weighted_quantile(phi[, 1L], w, probs = 0.5),
      mu_scale = weighted_quantile(phi[, 2L], w, probs = 0.5),
      mu_shift = weighted_quantile(phi[, 3L], w, probs = 0.5),
      sigma2_shape = weighted_quantile(exp(phi[, 4L]), w, probs = 0.5),
      sigma2_scale = weighted_quantile(exp(phi[, 5L]), w, probs = 0.5),
      sigma2_shift = weighted_quantile(exp(phi[, 6L]), w, probs = 0.5)
    )
  )
}

stan_reference_summary <- function(stan_res) {
  stan_mu <- as.matrix(stan_res$draws$mu)
  stan_sigma2 <- as.matrix(stan_res$draws$sigma2)
  list(
    mean = c(
      mu_shape = mean(stan_mu[, 1L]),
      mu_scale = mean(stan_mu[, 2L]),
      mu_shift = mean(stan_mu[, 3L]),
      sigma2_shape = mean(stan_sigma2[, 1L]),
      sigma2_scale = mean(stan_sigma2[, 2L]),
      sigma2_shift = mean(stan_sigma2[, 3L])
    ),
    q2.5 = c(
      mu_shape = as.numeric(quantile(stan_mu[, 1L], 0.025)),
      mu_scale = as.numeric(quantile(stan_mu[, 2L], 0.025)),
      mu_shift = as.numeric(quantile(stan_mu[, 3L], 0.025)),
      sigma2_shape = as.numeric(quantile(stan_sigma2[, 1L], 0.025)),
      sigma2_scale = as.numeric(quantile(stan_sigma2[, 2L], 0.025)),
      sigma2_shift = as.numeric(quantile(stan_sigma2[, 3L], 0.025))
    ),
    q97.5 = c(
      mu_shape = as.numeric(quantile(stan_mu[, 1L], 0.975)),
      mu_scale = as.numeric(quantile(stan_mu[, 2L], 0.975)),
      mu_shift = as.numeric(quantile(stan_mu[, 3L], 0.975)),
      sigma2_shape = as.numeric(quantile(stan_sigma2[, 1L], 0.975)),
      sigma2_scale = as.numeric(quantile(stan_sigma2[, 2L], 0.975)),
      sigma2_shift = as.numeric(quantile(stan_sigma2[, 3L], 0.975))
    )
  )
}

flag_hard_subjects <- function(outer_subject, data_list, cfg) {
  S <- length(outer_subject)
  shift_sd <- numeric(S)
  near_boundary_mass <- numeric(S)
  rounds <- numeric(S)
  mcse <- numeric(S)

  for (i in seq_len(S)) {
    out <- outer_subject[[i]]
    w <- as.numeric(out$w)
    w <- w / sum(w)
    eta_shift <- as.numeric(out$Theta[, 3L])
    shift <- exp(eta_shift)
    min_y <- min(data_list[[i]])

    shift_sd[i] <- weighted_sd(eta_shift, w)
    near_boundary_mass[i] <- sum(w * (shift > cfg$near_boundary_ratio * min_y))
    rounds[i] <- out$meta$rounds %||% length(out$meta$lambda_hist)
    mcse[i] <- out$mcse_logZ %||% NA_real_
  }

  q_sd <- stats::quantile(shift_sd, cfg$shift_sd_quantile, names = FALSE)
  q_rd <- stats::quantile(rounds, cfg$round_quantile, names = FALSE)
  q_mcse <- stats::quantile(mcse[is.finite(mcse)], cfg$mcse_quantile, names = FALSE)

  flag <- which(
    (shift_sd >= q_sd) |
      (near_boundary_mass > 0.20) |
      (rounds >= q_rd) |
      (is.finite(mcse) & mcse >= q_mcse)
  )

  if (!length(flag)) return(integer(0))

  score <- rank(shift_sd, ties.method = "average") +
    rank(near_boundary_mass, ties.method = "average") +
    rank(rounds, ties.method = "average") +
    rank(ifelse(is.finite(mcse), mcse, min(mcse, na.rm = TRUE)), ties.method = "average")

  flag <- flag[order(score[flag], decreasing = TRUE)]
  as.integer(head(flag, cfg$max_hard_subjects))
}


outer_final_lambda <- function(out) {
  lam <- out$final_lambda
  if (is.null(lam) || !length(lam)) lam <- tail(out$meta$lambda_hist %||% numeric(0), 1L)
  if (is.null(lam) || !length(lam)) return(NA_real_)
  as.numeric(lam[1L])
}

run_outer_subjects <- function(
    data_list,
    seed_run,
    outer_cfg,
    hard_idx = integer(0),
    mu_ref_override = NULL,
    Sigma_ref_override = NULL,
    phase_label = "default",
    M_subject_override = NULL,
    n_mcmc_moves_override = NULL,
    max_rounds_override = NULL
) {
  theta_names <- c("eta_shape", "eta_scale", "eta_shift")
  mu_ref <- if (is.null(mu_ref_override)) setNames(rep(0, 3L), theta_names) else setNames(as.numeric(mu_ref_override), theta_names)
  Sigma_ref <- if (is.null(Sigma_ref_override)) diag(1, 3L) else as.matrix(Sigma_ref_override)
  colnames(Sigma_ref) <- rownames(Sigma_ref) <- theta_names
  Sigma_ref <- as.matrix(Matrix::nearPD(Sigma_ref, conv.tol = 1e-7)$mat)
  diag(Sigma_ref) <- pmax(diag(Sigma_ref), 1e-8)

  M_default <- as.integer(M_subject_override$default %||% outer_cfg$M_subject_default)
  M_hard <- as.integer(M_subject_override$hard %||% outer_cfg$M_subject_hard)
  moves_default <- as.integer(n_mcmc_moves_override$default %||% outer_cfg$n_mcmc_moves_default)
  moves_hard <- as.integer(n_mcmc_moves_override$hard %||% outer_cfg$n_mcmc_moves_hard)
  max_rounds_use <- as.integer(max_rounds_override %||% 150L)
  phase_offset <- sum(utf8ToInt(as.character(phase_label)))

  S <- length(data_list)
  mclapply(
    seq_len(S),
    function(i) {
      is_hard <- i %in% hard_idx
      enhanced_smc_elite(
        data = data_list[[i]],
        loglik_fn = loglik_shifted_gamma,
        mu_ref = mu_ref,
        Sigma_ref = Sigma_ref,
        M = if (is_hard) M_hard else M_default,
        resample_threshold = 0.6,
        n_mcmc_moves = if (is_hard) moves_hard else moves_default,
        max_rounds = max_rounds_use,
        G_mix = 8L,
        da_enable = TRUE,
        gss_enable = TRUE,
        ll_cache_enable = TRUE,
        deterministic_resampling = FALSE,
        n_cores = 1L,
        seed = as.integer(seed_run + 1000L + 31L * phase_offset + i),
        verbose = FALSE
      )
    },
    mc.cores = 1L
  )
}

build_caches <- function(outer_subject, data_list, cache_cfg, phi_anchor, gaussian_map_fn,
                         hard_idx = integer(0)) {
  S <- length(outer_subject)
  budget_idx <- sort(unique(as.integer(hard_idx)))
  caches <- mclapply(
    seq_len(S),
    function(i) {
      budget <- resolve_cache_budget(cache_cfg, is_hard = i %in% budget_idx)
      build_subject_cache_from_smc(
        smc_out = outer_subject[[i]],
        data = data_list[[i]],
        subj_id = i,
        loglik_fn = loglik_shifted_gamma,
        M = budget$M,
        K_batches = budget$K_batches,
        deterministic_counts = cache_cfg$deterministic_counts,
        n_cores = 1L,
        cache_bridge_enable = isTRUE(cache_cfg$cache_bridge_enable),
        cache_bridge_weight = cache_cfg$cache_bridge_weight,
        cache_defensive_weight = cache_cfg$cache_defensive_weight,
        cache_defensive_df = cache_cfg$cache_defensive_df,
        defensive_t_eps = cache_cfg$defensive_t_eps,
        defensive_t_df = cache_cfg$defensive_t_df,
        phi_anchor = phi_anchor,
        phi_anchors = NULL,
        gaussian_map_fn = gaussian_map_fn
      )
    },
    mc.cores = 1L
  )
  for (i in seq_along(caches)) {
    d_i <- ncol(caches[[i]]$batches[[1]]$Theta)
    if (is.null(caches[[i]]$gaussian_map)) {
      caches[[i]] <- register_gaussian_prior_map(
        caches[[i]],
        function(phi) gaussian_map_fn(phi, d_i)
      )
    }
  }
  caches
}

build_surrogates <- function(outer_subject, data_list, loglik_fn, seed_base = 12345L) {
  S <- length(outer_subject)
  mclapply(
    seq_len(S),
    function(i) {
      build_subject_surrogate_from_outer(
        smc_out = outer_subject[[i]],
        data = data_list[[i]],
        subj_id = i,
        loglik_fn = loglik_fn,
        base_seed = as.integer(seed_base + 1009L * i)
      )
    },
    mc.cores = 1L
  )
}

.make_rtheta_given_phi <- function(gaussian_map_fn) {
  function(phi, n, aux = NULL) {
    d <- as.integer(length(phi) / 2L)
    pars <- gaussian_map_fn(phi, d)
    K <- as.matrix(pars$Sigma_inv)
    L <- chol(K)
    Sig <- chol2inv(L)
    mvtnorm::rmvnorm(
      n = as.integer(n),
      mean = as.numeric(pars$mu),
      sigma = as.matrix((Sig + t(Sig)) / 2)
    )
  }
}

probe_pm_quality <- function(
    phi_probe,
    surrogates,
    data_list,
    gaussian_map_fn,
    adaptive_ctl,
    M_by_subject,
    seed_base = 123L
) {
  S <- length(surrogates)
  rtheta_given_phi <- .make_rtheta_given_phi(gaussian_map_fn)
  rows <- vector("list", S)
  for (j in seq_len(S)) {
    cb <- list(
      gaussian_map_fn = gaussian_map_fn,
      rtheta_given_phi = rtheta_given_phi,
      loglik_fn = loglik_shifted_gamma,
      data_i = data_list[[j]]
    )
    est <- estimate_log_marginal_subject_pm(
      phi = as.numeric(phi_probe),
      surrogate = surrogates[[j]],
      callbacks = cb,
      aux_state_i = list(seed = as.integer(seed_base + 971L * j), M = as.integer(M_by_subject[j])),
      control = adaptive_ctl
    )
    rows[[j]] <- data.frame(
      subj = j,
      M = as.integer(M_by_subject[j]),
      ess_norm = as.numeric(est$ess_norm),
      ess_abs = as.numeric(est$ess_is),
      khat = as.numeric(est$khat),
      var_proxy = as.numeric(est$var_proxy),
      has_error = !is.null(est$error),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

pm_quality_gate_surrogates <- function(
    outer_subject,
    surrogates,
    data_list,
    phi_anchor,
    gaussian_map_fn,
    adaptive_ctl,
    gate_cfg,
    seed_run,
    outer_cfg,
    mu_ref_anchor,
    Sigma_ref_anchor,
    outer_ref_cfg
) {
  S <- length(surrogates)
  ladder <- sort(unique(as.integer(gate_cfg$M_ladder)))
  if (!length(ladder)) ladder <- c(256L, 512L, 1024L, 2048L)
  m0 <- as.integer(gate_cfg$M_probe %||% adaptive_ctl$M_default %||% 128L)
  m0 <- ladder[which.min(abs(ladder - m0))]
  M_by_subject <- rep(as.integer(m0), S)

  diag_history <- list()
  mix_fallback <- FALSE
  fail_idx <- integer(0)

  inc_ladder <- function(m) {
    cand <- ladder[ladder > as.integer(m)]
    if (!length(cand)) as.integer(max(ladder)) else as.integer(cand[1L])
  }

  for (pass in seq_len(as.integer(max(1L, gate_cfg$max_passes)))) {
    probe <- probe_pm_quality(
      phi_probe = phi_anchor,
      surrogates = surrogates,
      data_list = data_list,
      gaussian_map_fn = gaussian_map_fn,
      adaptive_ctl = adaptive_ctl,
      M_by_subject = M_by_subject,
      seed_base = as.integer(seed_run + 13000L + 811L * pass)
    )
    median_ess_abs <- suppressWarnings(stats::median(probe$ess_abs, na.rm = TRUE))
    q90_var <- suppressWarnings(as.numeric(stats::quantile(probe$var_proxy[is.finite(probe$var_proxy)], probs = 0.90, names = FALSE, na.rm = TRUE)))
    if (!is.finite(q90_var)) q90_var <- Inf
    q90_khat <- suppressWarnings(as.numeric(stats::quantile(probe$khat[is.finite(probe$khat)], probs = 0.90, names = FALSE, na.rm = TRUE)))
    if (!is.finite(q90_khat)) q90_khat <- Inf
    fail_idx <- which(
      (as.logical(probe$has_error) %||% rep(FALSE, nrow(probe))) |
        !is.finite(probe$ess_abs) | probe$ess_abs < as.numeric(gate_cfg$ess_abs_threshold) |
        !is.finite(probe$var_proxy) | probe$var_proxy > as.numeric(gate_cfg$var_proxy_threshold)
    )
    global_ok <- is.finite(median_ess_abs) &&
      median_ess_abs >= as.numeric(gate_cfg$ess_abs_threshold) &&
      is.finite(q90_var) &&
      q90_var <= as.numeric(gate_cfg$var_proxy_q90_threshold)
    diag_history[[length(diag_history) + 1L]] <- list(
      pass = pass,
      probe = probe,
      median_ess_abs = median_ess_abs,
      q90_var = q90_var,
      q90_khat = q90_khat,
      fail_idx = as.integer(fail_idx),
      global_ok = global_ok,
      mix_fallback = mix_fallback
    )
    cat(sprintf(
      "[PM gate] pass %d: median ESS_abs=%.2f, q90 var=%.3f, q90 khat=%.3f, failing=%d\n",
      pass, median_ess_abs, q90_var, q90_khat, length(fail_idx)
    ))
    if (global_ok) break
    if (length(fail_idx)) {
      for (j in fail_idx) M_by_subject[j] <- inc_ladder(M_by_subject[j])
    } else {
      # Global failure without explicit offenders: bump top-variance subjects.
      ord_var <- order(probe$var_proxy, decreasing = TRUE, na.last = NA)
      if (length(ord_var)) {
        bump <- head(ord_var, min(5L, length(ord_var)))
        for (j in bump) M_by_subject[j] <- inc_ladder(M_by_subject[j])
      }
    }
    if (isTRUE(gate_cfg$enable_mixture_fallback) && !mix_fallback) {
      adaptive_ctl$w_local <- as.numeric(gate_cfg$fallback_w_local)
      adaptive_ctl$w_phi_anchor <- as.numeric(gate_cfg$fallback_w_phi_anchor)
      adaptive_ctl$w_defensive <- as.numeric(gate_cfg$fallback_w_defensive)
      mix_fallback <- TRUE
    }
  }

  # Optional one-shot outer rerun for the worst failing subjects.
  rerun_info <- NULL
  final_probe <- diag_history[[length(diag_history)]]$probe
  final_fail <- as.integer(diag_history[[length(diag_history)]]$fail_idx)
  if (isTRUE(gate_cfg$rerun_outer_on_fail) && length(final_fail)) {
    sev <- rep(0, nrow(final_probe))
    ess_thr <- as.numeric(gate_cfg$ess_abs_threshold)
    var_thr <- as.numeric(gate_cfg$var_proxy_threshold)
    sev <- pmax((ess_thr - final_probe$ess_abs) / pmax(ess_thr, 1e-8), 0) +
      pmax((final_probe$var_proxy - var_thr) / pmax(var_thr, 1e-8), 0)
    sev[!is.finite(sev)] <- 10
    ord_fail <- final_fail[order(sev[final_fail], decreasing = TRUE)]
    max_rerun <- as.integer(gate_cfg$outer_rerun_max_subjects %||% length(ord_fail))
    rerun_subjects <- head(ord_fail, max(1L, min(length(ord_fail), max_rerun)))

    cat("[PM gate] rerunning outer for failing subjects: ", paste(rerun_subjects, collapse = ","), "\n", sep = "")
    outer_retry <- run_outer_subjects(
      data_list = data_list[rerun_subjects],
      seed_run = as.integer(seed_run + 70000L),
      outer_cfg = outer_cfg,
      hard_idx = seq_along(rerun_subjects),
      mu_ref_override = mu_ref_anchor,
      Sigma_ref_override = Sigma_ref_anchor,
      phase_label = "pm_gate_rerun",
      M_subject_override = list(
        default = as.integer(gate_cfg$outer_rerun_M_default %||% outer_cfg$M_subject_hard),
        hard = as.integer(gate_cfg$outer_rerun_M_hard %||% (2L * (outer_cfg$M_subject_hard %||% 1200L)))
      ),
      n_mcmc_moves_override = list(
        default = as.integer(gate_cfg$outer_rerun_n_moves_default %||% max(outer_cfg$n_mcmc_moves_hard %||% 4L, 4L)),
        hard = as.integer(gate_cfg$outer_rerun_n_moves_hard %||% max(outer_cfg$n_mcmc_moves_hard %||% 4L, 6L))
      ),
      max_rounds_override = as.integer(outer_ref_cfg$max_rounds_hard %||% 350L)
    )
    for (k in seq_along(rerun_subjects)) outer_subject[[rerun_subjects[k]]] <- outer_retry[[k]]
    sur_retry <- build_surrogates(
      outer_subject = outer_subject[rerun_subjects],
      data_list = data_list[rerun_subjects],
      loglik_fn = loglik_shifted_gamma,
      seed_base = as.integer(seed_run + 71000L)
    )
    for (k in seq_along(rerun_subjects)) surrogates[[rerun_subjects[k]]] <- sur_retry[[k]]

    probe2 <- probe_pm_quality(
      phi_probe = phi_anchor,
      surrogates = surrogates,
      data_list = data_list,
      gaussian_map_fn = gaussian_map_fn,
      adaptive_ctl = adaptive_ctl,
      M_by_subject = M_by_subject,
      seed_base = as.integer(seed_run + 72000L)
    )
    median_ess2 <- suppressWarnings(stats::median(probe2$ess_abs, na.rm = TRUE))
    q90_var2 <- suppressWarnings(as.numeric(stats::quantile(probe2$var_proxy[is.finite(probe2$var_proxy)], probs = 0.90, names = FALSE, na.rm = TRUE)))
    if (!is.finite(q90_var2)) q90_var2 <- Inf
    q90_khat2 <- suppressWarnings(as.numeric(stats::quantile(probe2$khat[is.finite(probe2$khat)], probs = 0.90, names = FALSE, na.rm = TRUE)))
    if (!is.finite(q90_khat2)) q90_khat2 <- Inf
    final_probe <- probe2
    final_fail <- which(
      (as.logical(probe2$has_error) %||% rep(FALSE, nrow(probe2))) |
        !is.finite(probe2$ess_abs) | probe2$ess_abs < as.numeric(gate_cfg$ess_abs_threshold) |
        !is.finite(probe2$var_proxy) | probe2$var_proxy > as.numeric(gate_cfg$var_proxy_threshold)
    )
    rerun_info <- list(
      subjects = as.integer(rerun_subjects),
      median_ess = median_ess2,
      q90_var = q90_var2,
      q90_khat = q90_khat2
    )
    diag_history[[length(diag_history) + 1L]] <- list(
      pass = "post_rerun",
      probe = probe2,
      median_ess_abs = median_ess2,
      q90_var = q90_var2,
      q90_khat = q90_khat2,
      fail_idx = as.integer(final_fail),
      global_ok = is.finite(median_ess2) &&
        median_ess2 >= as.numeric(gate_cfg$ess_abs_threshold) &&
        is.finite(q90_var2) &&
        q90_var2 <= as.numeric(gate_cfg$var_proxy_q90_threshold),
      mix_fallback = mix_fallback
    )
    cat(sprintf(
      "[PM gate] post-rerun: median ESS_abs=%.2f, q90 var=%.3f, q90 khat=%.3f, failing=%d\n",
      median_ess2, q90_var2, q90_khat2, length(final_fail)
    ))
  }

  final_median_ess <- suppressWarnings(stats::median(final_probe$ess_abs, na.rm = TRUE))
  final_q90_var <- suppressWarnings(as.numeric(stats::quantile(final_probe$var_proxy[is.finite(final_probe$var_proxy)], probs = 0.90, names = FALSE, na.rm = TRUE)))
  if (!is.finite(final_q90_var)) final_q90_var <- Inf
  final_q90_khat <- suppressWarnings(as.numeric(stats::quantile(final_probe$khat[is.finite(final_probe$khat)], probs = 0.90, names = FALSE, na.rm = TRUE)))
  if (!is.finite(final_q90_khat)) final_q90_khat <- Inf
  global_ok <- is.finite(final_median_ess) &&
    final_median_ess >= as.numeric(gate_cfg$ess_abs_threshold) &&
    is.finite(final_q90_var) &&
    final_q90_var <= as.numeric(gate_cfg$var_proxy_q90_threshold)

  list(
    outer_subject = outer_subject,
    surrogates = surrogates,
    adaptive_ctl = adaptive_ctl,
    M_by_subject = as.integer(M_by_subject),
    diag_history = diag_history,
    probe_final = final_probe,
    fail_subjects = as.integer(final_fail),
    global_ok = global_ok,
    median_ess = final_median_ess,
    q90_var = final_q90_var,
    q90_khat = final_q90_khat,
    rerun_info = rerun_info
  )
}

run_inner <- function(caches, prior, gaussian_map_fn, inner_cfg, seed_run,
                      data_list = NULL,
                      inner_ll_mode = "fixed_cache",
                      pm_aux_mode = "batch_idx",
                      adaptive_pm_control = NULL,
                      rtheta_given_phi = NULL) {
  inner_data <- if (identical(inner_ll_mode, "adaptive_pm")) data_list else NULL
  outer_smc_phi_batch(
    caches = caches,
    rprior_phi = prior$rprior,
    logprior_phi = prior$lprior,
    gaussian_map_fn = gaussian_map_fn,
    rtheta_given_phi = rtheta_given_phi,
    data_list = inner_data,
    loglik_fn = if (!is.null(inner_data)) loglik_shifted_gamma else NULL,
    M = inner_cfg$M,
    cess_target = inner_cfg$cess_target,
    resample_threshold = inner_cfg$resample_threshold,
    n_moves = inner_cfg$n_moves,
    rw_scale_init = inner_cfg$rw_scale_init,
    pm_mode = "strict",
    inner_ll_mode = inner_ll_mode,
    pm_aux_mode = pm_aux_mode,
    adaptive_pm_control = adaptive_pm_control %||% list(),
    block_refresh_every = inner_cfg$block_refresh_every,
    block_refresh_frac = inner_cfg$block_refresh_frac,
    auto_enrich_enable = FALSE,
    auto_enrich_every = 10L,
    auto_enrich_lambda_gate = 0.25,
    auto_enrich_ess_thresh = 0.20,
    auto_enrich_k_thresh = 0.70,
    auto_enrich_max_units = 2L,
    diag_enable = FALSE,
    collect_round_diagnostics = TRUE,
    max_rounds = inner_cfg$max_rounds,
    seed = seed_run,
    verbose = FALSE
  )
}

summarize_cache_quality <- function(caches, phi_anchor) {
  phi_anchor <- as.numeric(phi_anchor)
  subject_rows <- vector("list", length(caches))
  per_batch <- vector("list", length(caches))
  for (i in seq_along(caches)) {
    q <- subject_cache_quality_at_phi(caches[[i]], phi_anchor)

    subject_rows[[i]] <- data.frame(
      subj = i,
      median_ess_abs = as.numeric(q$summary$median_ess_abs),
      median_ess_norm = as.numeric(q$summary$median_ess_norm),
      sd_logZ_batches = as.numeric(q$summary$sd_logZ_batches),
      stringsAsFactors = FALSE
    )

    btab <- q$per_batch
    per_batch[[i]] <- data.frame(
      b = btab$b,
      ess_abs = btab$ess_abs,
      ess_norm = btab$ess_norm,
      logZ = btab$logZ,
      stringsAsFactors = FALSE
    )
  }

  subj_df <- do.call(rbind, subject_rows)
  list(
    subject = subj_df,
    per_batch = per_batch,
    global = list(
      median_ess_abs = stats::median(subj_df$median_ess_abs, na.rm = TRUE),
      median_ess_norm = stats::median(subj_df$median_ess_norm, na.rm = TRUE),
      sum_sd2 = sum(subj_df$sd_logZ_batches[is.finite(subj_df$sd_logZ_batches)]^2)
    )
  )
}

quality_gate_enrich <- function(caches, data_list, phi_anchor, gaussian_map_fn,
                                gate_cfg, enrich_cfg, seed_run) {
  diag_history <- list()
  actions <- list()

  q <- summarize_cache_quality(caches, phi_anchor)
  diag_history[[1L]] <- q
  if (!isTRUE(gate_cfg$enable)) {
    return(list(caches = caches, diag_history = diag_history, actions = actions, final = q))
  }

  prev_sum_sd2 <- q$global$sum_sd2
  for (pass in seq_len(enrich_cfg$max_passes)) {
    fail <- which(
      q$subject$median_ess_abs < gate_cfg$ess_abs_threshold |
        q$subject$sd_logZ_batches > gate_cfg$sd_logZ_threshold
    )
    if (!length(fail)) break
    if (q$global$sum_sd2 <= gate_cfg$target_sum_sd2) break

    severity <- numeric(nrow(q$subject))
    severity[fail] <- pmax(q$subject$sd_logZ_batches[fail] / gate_cfg$sd_logZ_threshold, 0) +
      pmax((gate_cfg$ess_abs_threshold - q$subject$median_ess_abs[fail]) / gate_cfg$ess_abs_threshold, 0)
    select <- fail[order(severity[fail], decreasing = TRUE)]
    select <- head(select, enrich_cfg$max_subjects_per_pass)

    pass_action <- list(pass = pass, subjects = as.integer(select), details = list())
    for (idx in seq_along(select)) {
      s <- select[idx]
      btab <- q$per_batch[[s]]
      b_worst <- btab$b[which.min(btab$ess_norm)]

      cache_new <- enrich_subject_cache_with_anchor_surgical(
        cache = caches[[s]],
        phi_probe = phi_anchor,
        b = b_worst,
        log_prior_theta_given_phi_mat = NULL,
        gaussian_map_fn = gaussian_map_fn,
        data = data_list[[s]],
        loglik_fn = loglik_shifted_gamma,
        eps_prior_base = enrich_cfg$eps_prior_base,
        weak_inflate_factor = enrich_cfg$weak_inflate_factor,
        replace_frac = NULL,
        sobol_seed = as.integer(seed_run + 50000L + 97L * s + 11L * pass),
        n_cores = 1L,
        verbose = FALSE,
        subj_id = s
      )

      cache_new <- augment_cache_M(
        cache = cache_new,
        data = data_list[[s]],
        loglik_fn = loglik_shifted_gamma,
        M_add = enrich_cfg$m_add_per_subject,
        sobol_seed = as.integer(seed_run + 60000L + 131L * s + 13L * pass),
        n_cores = 1L
      )

      caches[[s]] <- cache_new
      pass_action$details[[length(pass_action$details) + 1L]] <- list(
        subj = s,
        worst_batch = b_worst,
        m_add = enrich_cfg$m_add_per_subject
      )
    }

    actions[[length(actions) + 1L]] <- pass_action
    q_new <- summarize_cache_quality(caches, phi_anchor)
    diag_history[[length(diag_history) + 1L]] <- q_new

    improve <- prev_sum_sd2 - q_new$global$sum_sd2
    q <- q_new
    if (q$global$sum_sd2 <= gate_cfg$target_sum_sd2) break
    if (improve < enrich_cfg$min_sum_sd2_improve) break
    prev_sum_sd2 <- q$global$sum_sd2
  }

  list(caches = caches, diag_history = diag_history, actions = actions, final = q)
}

run_single_experiment <- function(experiment_id, toggles, seed_run, base_cfg,
                                  data_list, prior, gaussian_map_fn) {
  t0_total <- proc.time()[3]

  cfg <- base_cfg
  cfg$experiment_id <- experiment_id
  cfg$seed <- seed_run
  mode <- as.character(toggles$mode %||% "")
  if (!mode %in% c("fixed_cache", "adaptive_pm")) {
    stop("Unsupported experiment mode: ", mode)
  }
  if (mode == "fixed_cache") {
    # Consolidated E4 baseline settings.
    cfg$cache_strategy$cache_bridge_enable <- TRUE
    cfg$cache_strategy$cache_bridge_weight <- 0.20
    cfg$cache_strategy$cache_defensive_weight <- 0.05
    cfg$cache_strategy$cache_defensive_df <- 3L
    cfg$quality_gate$enable <- TRUE
    cfg$outer_escalation$enable <- TRUE
    cfg$inner$inner_ll_mode <- "fixed_cache"
    cfg$inner$pm_aux_mode <- "batch_idx"
  } else {
    # Consolidated H1 strict-anchor settings.
    cfg$inner$inner_ll_mode <- "adaptive_pm"
    cfg$inner$pm_aux_mode <- "rng_stream"
    cfg$outer_escalation$enable <- TRUE
  }

  cat(sprintf("\n=== %s | seed %d ===\n", experiment_id, seed_run))

  hard_idx <- integer(0)
  outer_subject <- NULL
  elapsed_outer <- 0
  rtheta_given_phi_fn <- .make_rtheta_given_phi(gaussian_map_fn)

  pilot_caches <- NULL
  final_caches <- NULL
  pilot_cache_quality <- list(
    subject = data.frame(subj = integer(0), median_ess_abs = numeric(0), median_ess_norm = numeric(0), sd_logZ_batches = numeric(0)),
    per_batch = list(),
    global = list(median_ess_abs = NA_real_, median_ess_norm = NA_real_, sum_sd2 = NA_real_)
  )
  gate_res <- list(
    caches = NULL,
    diag_history = list(),
    actions = list(),
    final = list(global = list(median_ess_abs = NA_real_, median_ess_norm = NA_real_, sum_sd2 = NA_real_))
  )
  elapsed_cache_pilot <- 0
  elapsed_pilot <- 0
  elapsed_cache_final <- 0
  elapsed_gate <- 0

  if (identical(cfg$inner$inner_ll_mode, "adaptive_pm")) {
    t0_outer <- proc.time()[3]
    pilot_idx <- select_pilot_subjects(data_list, cfg$pilot_strategy, seed_run)
    pilot_data <- data_list[pilot_idx]
    cat(sprintf("Pilot outer subjects: %d/%d\n", length(pilot_idx), length(data_list)))

    outer_pilot <- run_outer_subjects(
      data_list = pilot_data,
      seed_run = as.integer(seed_run + 5000L),
      outer_cfg = cfg$outer_escalation,
      hard_idx = integer(0),
      phase_label = "pilot",
      M_subject_override = list(
        default = cfg$pilot_strategy$M_subject_pilot,
        hard = cfg$pilot_strategy$M_subject_pilot
      ),
      n_mcmc_moves_override = list(
        default = cfg$pilot_strategy$n_mcmc_moves_pilot,
        hard = cfg$pilot_strategy$n_mcmc_moves_pilot
      ),
      max_rounds_override = cfg$pilot_strategy$max_rounds_pilot
    )

    t0_pilot <- proc.time()[3]
    phi_anchor <- as.numeric(estimate_phi_anchor_from_outer(
      outer_pilot,
      method = cfg$anchor_estimation$method %||% "eb_robust",
      winsor_probs = cfg$anchor_estimation$winsor_probs %||% c(0.05, 0.95),
      sigma2_floor = cfg$anchor_estimation$sigma2_floor %||% 1e-3,
      sigma2_ceiling = cfg$anchor_estimation$sigma2_ceiling %||% 2.5
    ))
    pilot_fit <- list(
      phi = matrix(phi_anchor, nrow = 1L),
      w = 1.0,
      meta = list(
        method = "pilot_outer_anchor",
        pilot_subjects = as.integer(pilot_idx)
      )
    )
    elapsed_pilot <- as.numeric(proc.time()[3] - t0_pilot)

    theta_names <- c("eta_shape", "eta_scale", "eta_shift")
    mu_ref_anchor <- setNames(phi_anchor[1:3], theta_names)
    Sigma_ref_anchor <- diag(exp(phi_anchor[4:6]) * as.numeric(cfg$anchor_estimation$outer_ref_inflation %||% 1.5), 3L)
    colnames(Sigma_ref_anchor) <- rownames(Sigma_ref_anchor) <- theta_names

    outer_subject <- run_outer_subjects(
      data_list = data_list,
      seed_run = as.integer(seed_run),
      outer_cfg = cfg$outer_escalation,
      hard_idx = integer(0),
      mu_ref_override = mu_ref_anchor,
      Sigma_ref_override = Sigma_ref_anchor,
      phase_label = "refine",
      max_rounds_override = cfg$outer_refinement$max_rounds_default %||% 220L
    )

    retry <- 0L
    bad_outer <- which(!is.finite(vapply(outer_subject, outer_final_lambda, numeric(1L))) |
      vapply(outer_subject, outer_final_lambda, numeric(1L)) < as.numeric(cfg$outer_refinement$lambda_subject_threshold %||% 0.999))
    while (length(bad_outer) > 0L && retry < as.integer(cfg$outer_refinement$max_retries %||% 2L)) {
      retry <- retry + 1L
      cat("Refining failed outer subjects (retry ", retry, "): ", paste(bad_outer, collapse = ","), "\n", sep = "")
      outer_retry <- run_outer_subjects(
        data_list = data_list[bad_outer],
        seed_run = as.integer(seed_run + 20000L + 3000L * retry),
        outer_cfg = cfg$outer_escalation,
        hard_idx = seq_along(bad_outer),
        mu_ref_override = mu_ref_anchor,
        Sigma_ref_override = Sigma_ref_anchor,
        phase_label = paste0("refine_retry_", retry),
        max_rounds_override = cfg$outer_refinement$max_rounds_hard %||% 350L
      )
      for (k in seq_along(bad_outer)) outer_subject[[bad_outer[k]]] <- outer_retry[[k]]
      bad_outer <- which(!is.finite(vapply(outer_subject, outer_final_lambda, numeric(1L))) |
        vapply(outer_subject, outer_final_lambda, numeric(1L)) < as.numeric(cfg$outer_refinement$lambda_subject_threshold %||% 0.999))
    }
    elapsed_outer <- as.numeric(proc.time()[3] - t0_outer)
    hard_idx <- as.integer(bad_outer)
    if (isTRUE(cfg$outer_refinement$strict_require_lambda1) && length(bad_outer) > 0L) {
      stop("Strict outer gate failed: subjects not converged to lambda gate: ", paste(bad_outer, collapse = ","))
    }

    t0_cache_pilot <- proc.time()[3]
    pilot_caches <- build_surrogates(
      outer_subject = outer_subject,
      data_list = data_list,
      loglik_fn = loglik_shifted_gamma,
      seed_base = as.integer(seed_run + 40000L)
    )
    elapsed_cache_pilot <- as.numeric(proc.time()[3] - t0_cache_pilot)

    t0_cache_final <- proc.time()[3]
    final_caches <- pilot_caches
    gate_res$caches <- final_caches
    gate_res$final <- pilot_cache_quality
    elapsed_cache_final <- as.numeric(proc.time()[3] - t0_cache_final)

    t0_gate <- proc.time()[3]
    pm_gate <- pm_quality_gate_surrogates(
      outer_subject = outer_subject,
      surrogates = final_caches,
      data_list = data_list,
      phi_anchor = phi_anchor,
      gaussian_map_fn = gaussian_map_fn,
      adaptive_ctl = cfg$adaptive_pm,
      gate_cfg = cfg$pm_quality_gate,
      seed_run = seed_run,
      outer_cfg = cfg$outer_escalation,
      mu_ref_anchor = mu_ref_anchor,
      Sigma_ref_anchor = Sigma_ref_anchor,
      outer_ref_cfg = cfg$outer_refinement
    )
    outer_subject <- pm_gate$outer_subject
    final_caches <- pm_gate$surrogates
    elapsed_gate <- as.numeric(proc.time()[3] - t0_gate)

    adaptive_ctl_inner <- pm_gate$adaptive_ctl
    adaptive_ctl_inner$initial_subject_M <- as.integer(pm_gate$M_by_subject)
    adaptive_ctl_inner$completion_mode <- cfg$inner_strict$completion_mode %||% adaptive_ctl_inner$completion_mode
    adaptive_ctl_inner$adapt_until_round <- as.integer(cfg$inner_strict$adapt_until_round %||% adaptive_ctl_inner$adapt_until_round)
    adaptive_ctl_inner$freeze_after_round <- as.integer(cfg$inner_strict$freeze_after_round %||% adaptive_ctl_inner$freeze_after_round)

    if (isTRUE(cfg$outer_refinement$strict_require_lambda1) && !isTRUE(pm_gate$global_ok)) {
      stop(
        sprintf(
          "Strict PM quality gate failed: median ESS_abs=%.2f, q90 var=%.3f (q90 khat=%.3f).",
          pm_gate$median_ess, pm_gate$q90_var, pm_gate$q90_khat
        )
      )
    }

    t0_inner <- proc.time()[3]
    inner_fit <- run_inner(
      caches = final_caches,
      prior = prior,
      gaussian_map_fn = gaussian_map_fn,
      inner_cfg = cfg$inner,
      seed_run = seed_run,
      data_list = data_list,
      inner_ll_mode = "adaptive_pm",
      pm_aux_mode = cfg$inner$pm_aux_mode,
      adaptive_pm_control = adaptive_ctl_inner,
      rtheta_given_phi = rtheta_given_phi_fn
    )
    elapsed_inner <- as.numeric(proc.time()[3] - t0_inner)

    final_lambda_inner <- suppressWarnings(as.numeric(tail(inner_fit$meta$lambda_hist %||% numeric(0), 1L)))
    if (isTRUE(cfg$inner_strict$strict_require_lambda1) &&
        (!is.finite(final_lambda_inner) || final_lambda_inner < as.numeric(cfg$inner_strict$lambda_strict_threshold %||% 0.9999))) {
      stop("Strict inner gate failed: final lambda did not reach strict threshold.")
    }

    gate_res <- list(
      caches = final_caches,
      diag_history = pm_gate$diag_history,
      actions = list(),
      final = list(
        global = list(
          median_ess_abs = pm_gate$median_ess,
          median_ess_norm = NA_real_,
          sum_sd2 = NA_real_
        ),
        pm_probe = pm_gate$probe_final,
        pm_gate_ok = pm_gate$global_ok,
        pm_q90_var = pm_gate$q90_var,
        pm_q90_khat = pm_gate$q90_khat,
        pm_fail_subjects = pm_gate$fail_subjects,
        rerun_info = pm_gate$rerun_info
      )
    )
  } else {
    t0_outer <- proc.time()[3]
    outer_subject <- run_outer_subjects(
      data_list = data_list,
      seed_run = seed_run,
      outer_cfg = cfg$outer_escalation,
      hard_idx = integer(0),
      phase_label = "fixed_outer"
    )

    if (cfg$outer_escalation$enable) {
      hard_idx <- flag_hard_subjects(outer_subject, data_list, cfg$outer_escalation)
      if (length(hard_idx)) {
        cat("Escalating outer for subjects:", paste(hard_idx, collapse = ","), "\n")
        outer_hard <- run_outer_subjects(
          data_list = data_list[hard_idx],
          seed_run = seed_run + 20000L,
          outer_cfg = cfg$outer_escalation,
          hard_idx = seq_along(hard_idx),
          phase_label = "fixed_outer_hard"
        )
        for (k in seq_along(hard_idx)) outer_subject[[hard_idx[k]]] <- outer_hard[[k]]
      }
    }
    elapsed_outer <- as.numeric(proc.time()[3] - t0_outer)

    t0_cache_pilot <- proc.time()[3]
    pilot_cache_cfg <- cfg$cache_strategy
    pilot_cache_cfg$cache_bridge_enable <- FALSE
    pilot_cache_cfg$cache_bridge_weight <- 0
    pilot_cache_cfg$cache_defensive_weight <- 0
    pilot_cache_cfg$hard_cache_multiplier <- 1.0
    pilot_cache_cfg$min_M_per <- 0L
    pilot_caches <- build_caches(
      outer_subject, data_list, pilot_cache_cfg,
      phi_anchor = NULL, gaussian_map_fn = gaussian_map_fn,
      hard_idx = integer(0)
    )
    elapsed_cache_pilot <- as.numeric(proc.time()[3] - t0_cache_pilot)

    t0_pilot <- proc.time()[3]
    pilot_inner_cfg <- cfg$inner
    pilot_inner_cfg$M <- cfg$pilot_inner$M
    pilot_inner_cfg$max_rounds <- cfg$pilot_inner$max_rounds
    pilot_inner_cfg$n_moves <- cfg$pilot_inner$n_moves
    pilot_inner_cfg$rw_scale_init <- cfg$pilot_inner$rw_scale_init
    pilot_fit <- run_inner(
      caches = pilot_caches,
      prior = prior,
      gaussian_map_fn = gaussian_map_fn,
      inner_cfg = pilot_inner_cfg,
      seed_run = as.integer(seed_run + 30000L),
      data_list = NULL,
      inner_ll_mode = "fixed_cache",
      pm_aux_mode = "batch_idx",
      adaptive_pm_control = NULL,
      rtheta_given_phi = NULL
    )
    phi_anchor <- as.numeric(weighted_mean(pilot_fit$phi, pilot_fit$w))
    pilot_cache_quality <- summarize_cache_quality(
      pilot_caches,
      phi_anchor = phi_anchor
    )
    elapsed_pilot <- as.numeric(proc.time()[3] - t0_pilot)

    t0_cache_final <- proc.time()[3]
    final_caches <- build_caches(
      outer_subject = outer_subject,
      data_list = data_list,
      cache_cfg = cfg$cache_strategy,
      phi_anchor = phi_anchor,
      gaussian_map_fn = gaussian_map_fn,
      hard_idx = hard_idx
    )
    elapsed_cache_final <- as.numeric(proc.time()[3] - t0_cache_final)

    t0_gate <- proc.time()[3]
    gate_res <- quality_gate_enrich(
      caches = final_caches,
      data_list = data_list,
      phi_anchor = phi_anchor,
      gaussian_map_fn = gaussian_map_fn,
      gate_cfg = cfg$quality_gate,
      enrich_cfg = cfg$enrichment,
      seed_run = seed_run
    )
    final_caches <- gate_res$caches
    elapsed_gate <- as.numeric(proc.time()[3] - t0_gate)

    t0_inner <- proc.time()[3]
    inner_fit <- run_inner(
      caches = final_caches,
      prior = prior,
      gaussian_map_fn = gaussian_map_fn,
      inner_cfg = cfg$inner,
      seed_run = seed_run,
      data_list = NULL,
      inner_ll_mode = "fixed_cache",
      pm_aux_mode = "batch_idx",
      adaptive_pm_control = NULL,
      rtheta_given_phi = NULL
    )
    elapsed_inner <- as.numeric(proc.time()[3] - t0_inner)
  }

  runtime <- list(
    outer_sec = elapsed_outer,
    cache_pilot_sec = elapsed_cache_pilot,
    pilot_inner_sec = elapsed_pilot,
    cache_final_sec = elapsed_cache_final,
    quality_gate_sec = elapsed_gate,
    inner_sec = elapsed_inner,
    total_sec = as.numeric(proc.time()[3] - t0_total)
  )

  inner_fit$meta$experiment_id <- experiment_id
  inner_fit$meta$seed <- seed_run
  inner_fit$meta$config <- cfg
  inner_fit$meta$runtime <- runtime
  inner_fit$meta$phi_anchor <- phi_anchor
  inner_fit$meta$quality_gate <- list(
    diag_history = gate_res$diag_history,
    actions = gate_res$actions
  )

  list(
    experiment_id = experiment_id,
    seed = seed_run,
    config = cfg,
    outer_subject = outer_subject,
    pilot_caches = pilot_caches,
    final_caches = final_caches,
    pilot_fit = pilot_fit,
    inner_fit = inner_fit,
    hard_subjects = hard_idx,
    pilot_cache_quality = pilot_cache_quality,
    runtime = runtime,
    cache_quality_final = gate_res$final
  )
}

score_run <- function(run_obj, stan_ref) {
  inner_sum <- inner_weighted_summary(run_obj$inner_fit)
  mean_err <- inner_sum$mean - stan_ref$mean
  overlap <- (inner_sum$median >= stan_ref$q2.5) & (inner_sum$median <= stan_ref$q97.5)

  data.frame(
    experiment_id = run_obj$experiment_id,
    seed = run_obj$seed,
    total_sec = run_obj$runtime$total_sec,
    inner_sec = run_obj$runtime$inner_sec,
    outer_sec = run_obj$runtime$outer_sec,
    log_evidence = run_obj$inner_fit$log_evidence,
    mcse_log_evidence = run_obj$inner_fit$mcse_log_evidence,
    final_lambda = tail(run_obj$inner_fit$meta$lambda_hist, 1L),
    mean_err_mu_shape = mean_err["mu_shape"],
    mean_err_mu_scale = mean_err["mu_scale"],
    mean_err_mu_shift = mean_err["mu_shift"],
    mean_err_sigma2_shape = mean_err["sigma2_shape"],
    mean_err_sigma2_scale = mean_err["sigma2_scale"],
    mean_err_sigma2_shift = mean_err["sigma2_shift"],
    overlap_all_medians = all(overlap),
    overlap_mu_shift = overlap["mu_shift"],
    overlap_sigma2_shift = overlap["sigma2_shift"],
    cache_median_ess_abs = run_obj$cache_quality_final$global$median_ess_abs,
    cache_sum_sd2 = run_obj$cache_quality_final$global$sum_sd2,
    phi_neff = 1 / sum(as.numeric(run_obj$inner_fit$w)^2),
    n_hard_subjects = length(run_obj$hard_subjects),
    stringsAsFactors = FALSE
  )
}

strict_publish_ok <- function(run_obj) {
  if (!identical(run_obj$config$inner$inner_ll_mode, "adaptive_pm")) return(TRUE)
  lam_outer <- vapply(run_obj$outer_subject, outer_final_lambda, numeric(1L))
  outer_thr <- as.numeric(run_obj$config$outer_refinement$lambda_strict_threshold %||% 0.9999)
  inner_thr <- as.numeric(run_obj$config$inner_strict$lambda_strict_threshold %||% 0.9999)
  lam_inner <- suppressWarnings(as.numeric(tail(run_obj$inner_fit$meta$lambda_hist %||% numeric(0), 1L)))
  outer_ok <- all(is.finite(lam_outer) & lam_outer >= outer_thr)
  inner_ok <- is.finite(lam_inner) && lam_inner >= inner_thr
  outer_ok && inner_ok
}

# ------------------------------- Data setup --------------------------------
dir.create("samples", showWarnings = FALSE, recursive = TRUE)

y_mat <- stan_res$data$y
S <- nrow(y_mat)
data_list <- lapply(seq_len(S), function(i) as.numeric(y_mat[i, ]))

prior <- make_prior_phi_diag(
  m0 = stan_res$priors$m0,
  s0 = stan_res$priors$s0,
  a = stan_res$priors$a0,
  b = stan_res$priors$b0,
  d = 3L
)
gaussian_map_fn <- phi_to_gaussian_params_diag_factory()
stan_ref <- stan_reference_summary(stan_res)

# ----------------------------- Run experiments -----------------------------
score_rows <- list()

for (exp_id in names(experiment_grid)) {
  toggles <- experiment_grid[[exp_id]]
  for (seed_run in base_config$seeds) {
    run_obj <- run_single_experiment(
      experiment_id = exp_id,
      toggles = toggles,
      seed_run = as.integer(seed_run),
      base_cfg = base_config,
      data_list = data_list,
      prior = prior,
      gaussian_map_fn = gaussian_map_fn
    )

    key <- sprintf("%s_seed%d", exp_id, seed_run)
    score_rows[[length(score_rows) + 1L]] <- score_run(run_obj, stan_ref)

    if (identical(exp_id, publish_experiment_id) && as.integer(seed_run) == as.integer(publish_seed)) {
      if (strict_publish_ok(run_obj)) {
        saveRDS(run_obj$outer_subject, file.path("samples", "outer_smc_results.rds"))
        saveRDS(run_obj$inner_fit, file.path("samples", "inner_smc_results.rds"))
        cat("Published canonical outputs for", key, "\n")
      } else {
        cat("Skipped canonical publish for", key, "(strict gates not satisfied).\n")
      }
    }
  }
}

score_df <- do.call(rbind, score_rows)
score_df <- score_df[order(score_df$experiment_id, score_df$seed), , drop = FALSE]

saveRDS(score_df, file.path("samples", "two_step_experiment_scorecard.rds"))

cat("\nSaved: samples/two_step_experiment_scorecard.rds\n")
cat("\nScorecard preview:\n")
print(score_df[, c(
  "experiment_id", "seed", "total_sec", "log_evidence", "mean_err_mu_shift",
  "mean_err_sigma2_shift", "cache_median_ess_abs", "cache_sum_sd2", "overlap_all_medians"
)], row.names = FALSE)
