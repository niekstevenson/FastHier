#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(mvtnorm)
})

source("smc_core.R")
source("new_SMC_cache.R")

.normalize_logweights <- function(logw) {
  logw <- as.numeric(logw)
  m <- max(logw)
  if (!is.finite(m)) {
    w <- rep(1 / length(logw), length(logw))
    return(list(log_norm = -Inf, w = w))
  }
  ww <- exp(logw - m)
  s <- sum(ww)
  if (!is.finite(s) || s <= 0) {
    w <- rep(1 / length(logw), length(logw))
    return(list(log_norm = -Inf, w = w))
  }
  list(log_norm = m + log(s), w = ww / s)
}

.ess_from_weights <- function(w) {
  w <- as.numeric(w)
  s <- sum(w)
  if (!is.finite(s) || s <= 0) return(0)
  w <- w / s
  1 / sum(w * w)
}

.resample_indices <- function(w) {
  sample.int(length(w), size = length(w), replace = TRUE, prob = w)
}

.zero_local_loglik_evals <- function() {
  setNames(rep.int(0L, 3L), c("initialization", "enrichment", "population_refreshes"))
}

.add_local_loglik_evals <- function(a, b) {
  nms <- names(.zero_local_loglik_evals())
  a <- setNames(as.integer(a[nms]), nms)
  b <- setNames(as.integer(b[nms]), nms)
  a[is.na(a)] <- 0L
  b[is.na(b)] <- 0L
  a + b
}

.finalize_local_loglik_evals <- function(x) {
  x <- setNames(as.integer(x[names(.zero_local_loglik_evals())]), names(.zero_local_loglik_evals()))
  x[is.na(x)] <- 0L
  c(x, total = sum(x))
}

.new_local_bank <- function(alpha, w, logg, local_obj, loglik = NULL,
                            source_idx = NULL,
                            proposal_tag = "prefit_gamma0") {
  alpha <- as.matrix(alpha)
  bank <- structure(
    list(
      alpha = alpha,
      w = as.numeric(w),
      logg = as.numeric(logg),
      loglik = if (is.null(loglik)) NULL else as.numeric(loglik),
      source_idx = if (is.null(source_idx)) NULL else as.integer(source_idx),
      d_theta = ncol(alpha),
      local_obj = local_obj,
      proposal_tag = as.character(proposal_tag)
    ),
    class = "local_bank"
  )
  .validate_local_bank(bank)
}

.validate_local_bank <- function(bank) {
  stopifnot(inherits(bank, "local_bank"))
  A <- bank$alpha
  if (!is.matrix(A)) A <- as.matrix(A)
  n <- nrow(A)
  if (length(bank$w) != n) stop("local_bank has inconsistent weight length.")
  if (length(bank$logg) != n) stop("local_bank has inconsistent logg length.")
  if (!is.null(bank$loglik) && length(bank$loglik) != n) stop("local_bank has inconsistent loglik length.")
  if (!is.null(bank$source_idx) && length(bank$source_idx) != n) stop("local_bank has inconsistent source_idx length.")
  sw <- sum(bank$w)
  if (!is.finite(sw) || sw <= 0) stop("local_bank weights are invalid.")
  bank$alpha <- A
  bank$w <- bank$w / sw
  bank$logw <- log(pmax(bank$w, .Machine$double.eps))
  bank
}

.local_bank_delta <- function(bank, phi, gaussian_map_fn, validate = TRUE) {
  if (isTRUE(validate)) bank <- .validate_local_bank(bank)
  logp <- .log_prior_gauss_mat(bank$alpha, gaussian_map_fn(as.numeric(phi), bank$d_theta))
  as.numeric(logp - bank$logg)
}

.local_bank_weight_stats <- function(bank, delta, lambda) {
  logw_lambda <- bank$logw + as.numeric(lambda) * as.numeric(delta)
  norm <- .normalize_logweights(logw_lambda)
  list(
    logZ = as.numeric(norm$log_norm),
    ess = .ess_from_weights(norm$w)
  )
}

build_local_bank <- function(local_obj, M_local = NULL, seed = 123L, proposal_tag = "prefit_gamma0") {
  stopifnot(inherits(local_obj, "local_exact_object"))
  if (!is.null(seed)) set.seed(as.integer(seed))

  Theta_all <- as.matrix(local_obj$theta_bank)
  n_all <- nrow(Theta_all)
  if (!is.finite(n_all) || n_all <= 0L) stop("local_obj$theta_bank is empty.")

  w_all <- pmax(as.numeric(local_obj$theta_bank_w), 0)
  sw <- sum(w_all)
  if (!is.finite(sw) || sw <= 0) {
    w_all <- rep(1 / n_all, n_all)
  } else {
    w_all <- w_all / sw
  }

  if (is.null(M_local) || as.integer(M_local) >= n_all) {
    idx <- seq_len(n_all)
    Theta <- Theta_all
    w_bank <- w_all
  } else {
    idx <- sample.int(n_all, size = as.integer(M_local), replace = TRUE, prob = w_all)
    Theta <- Theta_all[idx, , drop = FALSE]
    w_bank <- rep(1 / length(idx), length(idx))
  }

  logg <- log_working_prior_gaussian_mat(Theta, local_obj$working_prior)
  loglik <- NULL
  if (!is.null(local_obj$theta_bank_loglik)) {
    loglik <- as.numeric(local_obj$theta_bank_loglik[idx])
  }

  counts <- .zero_local_loglik_evals()
  list(
    bank = .new_local_bank(
      alpha = Theta,
      w = w_bank,
      logg = logg,
      loglik = loglik,
      source_idx = idx,
      local_obj = local_obj,
      proposal_tag = proposal_tag
    ),
    local_loglik_evals = counts
  )
}

print.local_bank <- function(x, ...) {
  cat(sprintf("<local_bank M=%d proposal=%s>\n", nrow(x$alpha), x$proposal_tag))
  invisible(x)
}

estimate_logZ <- function(x, ...) UseMethod("estimate_logZ")
predict_increment <- function(x, ...) UseMethod("predict_increment")
ess_under <- function(x, ...) UseMethod("ess_under")
enrich <- function(x, ...) UseMethod("enrich")

estimate_logZ.local_bank <- function(x, phi, lambda, gaussian_map_fn, ...) {
  x <- .validate_local_bank(x)
  delta <- .local_bank_delta(x, phi = phi, gaussian_map_fn = gaussian_map_fn, validate = FALSE)
  .local_bank_weight_stats(x, delta = delta, lambda = lambda)$logZ
}

predict_increment.local_bank <- function(x, phi, lambda_old, lambda_new, gaussian_map_fn, ...) {
  x <- .validate_local_bank(x)
  lambda_old <- as.numeric(lambda_old)
  lambda_new <- as.numeric(lambda_new)
  if (lambda_new < lambda_old - 1e-12) {
    stop("predict_increment.local_bank does not support decreasing lambda.")
  }
  delta <- .local_bank_delta(x, phi = phi, gaussian_map_fn = gaussian_map_fn, validate = FALSE)
  old_stats <- .local_bank_weight_stats(x, delta = delta, lambda = lambda_old)
  new_stats <- .local_bank_weight_stats(x, delta = delta, lambda = lambda_new)
  list(
    log_u = as.numeric(new_stats$logZ - old_stats$logZ),
    logZ_old = as.numeric(old_stats$logZ),
    logZ_new = as.numeric(new_stats$logZ),
    ess_new = as.numeric(new_stats$ess)
  )
}

ess_under.local_bank <- function(x, phi, lambda, gaussian_map_fn, ...) {
  x <- .validate_local_bank(x)
  delta <- .local_bank_delta(x, phi = phi, gaussian_map_fn = gaussian_map_fn, validate = FALSE)
  .local_bank_weight_stats(x, delta = delta, lambda = lambda)$ess
}

.estimate_logZ_subset <- function(bank, idx, phi, lambda, gaussian_map_fn) {
  idx <- as.integer(idx)
  if (!length(idx)) return(NA_real_)
  w_sub <- bank$w[idx]
  sw <- sum(w_sub)
  if (!is.finite(sw) || sw <= 0) return(NA_real_)
  w_sub <- w_sub / sw
  delta <- .local_bank_delta(bank, phi = phi, gaussian_map_fn = gaussian_map_fn)[idx]
  logsumexp(log(pmax(w_sub, .Machine$double.eps)) + as.numeric(lambda) * delta)
}

.prepare_bank_delta_matrix <- function(bank, phi, gaussian_map_fn) {
  bank <- .validate_local_bank(bank)
  phi <- as.matrix(phi)
  N <- nrow(phi)
  M <- nrow(bank$alpha)
  delta_mat <- matrix(0, nrow = N, ncol = M)
  for (n in seq_len(N)) {
    delta_mat[n, ] <- .local_bank_delta(
      bank,
      phi = phi[n, , drop = TRUE],
      gaussian_map_fn = gaussian_map_fn,
      validate = FALSE
    )
  }
  list(bank = bank, delta = delta_mat)
}

.bank_ess_values_from_delta <- function(bank, delta_mat, lambda) {
  logw_lambda <- sweep(as.numeric(lambda) * delta_mat, 2L, bank$logw, "+")
  log_norm <- .rowLogSumExp(logw_lambda)
  log_sum_sq <- .rowLogSumExp(2 * logw_lambda)
  as.numeric(exp(2 * log_norm - log_sum_sq))
}

.bank_split_diff_from_delta <- function(bank, delta_mat, lambda, idx_odd, idx_even) {
  if (length(idx_odd) < 4L || length(idx_even) < 4L) {
    return(rep.int(Inf, nrow(delta_mat)))
  }

  w_odd <- bank$w[idx_odd]
  sw_odd <- sum(w_odd)
  w_even <- bank$w[idx_even]
  sw_even <- sum(w_even)
  if (!is.finite(sw_odd) || sw_odd <= 0 || !is.finite(sw_even) || sw_even <= 0) {
    return(rep.int(Inf, nrow(delta_mat)))
  }

  logw_odd <- log(pmax(w_odd / sw_odd, .Machine$double.eps))
  logw_even <- log(pmax(w_even / sw_even, .Machine$double.eps))
  logZ_odd <- .rowLogSumExp(sweep(as.numeric(lambda) * delta_mat[, idx_odd, drop = FALSE], 2L, logw_odd, "+"))
  logZ_even <- .rowLogSumExp(sweep(as.numeric(lambda) * delta_mat[, idx_even, drop = FALSE], 2L, logw_even, "+"))
  as.numeric(abs(logZ_odd - logZ_even))
}

.bank_group_diagnostics_from_delta <- function(bank, delta_mat, lambda, ess_prob = 0.1, split_prob = 0.9) {
  idx_odd <- seq.int(1L, nrow(bank$alpha), by = 2L)
  idx_even <- seq.int(2L, nrow(bank$alpha), by = 2L)
  have_split <- length(idx_odd) >= 4L && length(idx_even) >= 4L

  ess_vals <- .bank_ess_values_from_delta(bank, delta_mat = delta_mat, lambda = lambda)
  split_diff <- .bank_split_diff_from_delta(
    bank,
    delta_mat = delta_mat,
    lambda = lambda,
    idx_odd = idx_odd,
    idx_even = idx_even
  )

  list(
    M = nrow(bank$alpha),
    ess_min = min(ess_vals, na.rm = TRUE),
    ess_q = as.numeric(stats::quantile(ess_vals, probs = ess_prob, names = FALSE, na.rm = TRUE)),
    split_q = as.numeric(stats::quantile(split_diff, probs = split_prob, names = FALSE, na.rm = TRUE)),
    have_split = have_split
  )
}

.bank_group_diagnostics <- function(bank, phi, lambda, gaussian_map_fn, ess_prob = 0.1, split_prob = 0.9) {
  prep <- .prepare_bank_delta_matrix(bank, phi = phi, gaussian_map_fn = gaussian_map_fn)
  .bank_group_diagnostics_from_delta(
    prep$bank,
    delta_mat = prep$delta,
    lambda = lambda,
    ess_prob = ess_prob,
    split_prob = split_prob
  )
}

.topup_local_bank_from_prefit <- function(bank, n_new, seed = NULL) {
  bank <- .validate_local_bank(bank)
  if (!is.null(seed)) set.seed(as.integer(seed))
  local_obj <- bank$local_obj

  Theta_all <- as.matrix(local_obj$theta_bank)
  n_all <- nrow(Theta_all)
  used <- unique(bank$source_idx %||% integer(0))
  remaining <- setdiff(seq_len(n_all), used)
  if (!length(remaining) || as.integer(n_new) <= 0L) {
    return(list(bank = bank, n_new = 0L))
  }

  n_new <- min(as.integer(n_new), length(remaining))
  w_rem <- pmax(as.numeric(local_obj$theta_bank_w[remaining]), 0)
  sw <- sum(w_rem)
  if (!is.finite(sw) || sw <= 0) {
    w_rem <- rep(1 / length(remaining), length(remaining))
  } else {
    w_rem <- w_rem / sw
  }
  pick <- sample.int(length(remaining), size = n_new, replace = FALSE, prob = w_rem)
  idx_new <- remaining[pick]

  alpha_new <- Theta_all[idx_new, , drop = FALSE]
  logg_new <- log_working_prior_gaussian_mat(alpha_new, local_obj$working_prior)
  loglik_new <- NULL
  if (!is.null(local_obj$theta_bank_loglik)) {
    loglik_new <- as.numeric(local_obj$theta_bank_loglik[idx_new])
  }

  alpha_all <- rbind(bank$alpha, alpha_new)
  logg_all <- c(bank$logg, logg_new)
  loglik_all <- if (is.null(bank$loglik) && is.null(loglik_new)) NULL else c(bank$loglik, loglik_new)
  source_all <- c(bank$source_idx %||% integer(0), idx_new)
  w_all <- rep(1 / nrow(alpha_all), nrow(alpha_all))

  list(
    bank = .new_local_bank(
      alpha = alpha_all,
      w = w_all,
      logg = logg_all,
      loglik = loglik_all,
      source_idx = source_all,
      local_obj = local_obj,
      proposal_tag = "prefit_gamma0_topped_up"
    ),
    n_new = as.integer(n_new)
  )
}

.refresh_local_bank_exact <- function(bank,
                                      refresh_particles = NULL,
                                      refresh_max_rounds = NULL,
                                      refresh_n_mcmc_moves = 3L,
                                      seed = NULL) {
  bank <- .validate_local_bank(bank)
  local_obj <- bank$local_obj
  theta_bank <- as.matrix(local_obj$theta_bank)
  d <- ncol(theta_bank)
  theta_names <- colnames(theta_bank)
  if (is.null(theta_names)) {
    theta_names <- names(local_obj$working_prior$mu)
    if (is.null(theta_names)) theta_names <- paste0("theta", seq_len(d))
  }

  mu_ref <- setNames(as.numeric(local_obj$working_prior$mu), theta_names)
  Sigma_ref <- as.matrix(local_obj$working_prior$Sigma)
  M_refresh <- as.integer(refresh_particles %||% nrow(theta_bank))
  max_rounds_use <- as.integer(refresh_max_rounds %||% 200L)

  ll_eval_counter_reset(enabled = TRUE)
  smc_out <- enhanced_smc_elite(
    data = local_obj$data,
    loglik_fn = local_obj$loglik_fn,
    mu_ref = mu_ref,
    Sigma_ref = Sigma_ref,
    M = M_refresh,
    resample_threshold = 0.6,
    n_mcmc_moves = as.integer(refresh_n_mcmc_moves),
    max_rounds = max_rounds_use,
    G_mix = 8L,
    da_enable = TRUE,
    gss_enable = TRUE,
    ll_cache_enable = TRUE,
    deterministic_resampling = FALSE,
    n_cores = 1L,
    seed = as.integer(seed %||% local_obj$base_seed),
    verbose = FALSE
  )
  n_ll <- ll_eval_counter_get()
  ll_eval_counter_disable()

  smc_out$working_prior <- local_obj$working_prior
  refreshed_local_obj <- build_local_exact_object(
    smc_out = smc_out,
    data = local_obj$data,
    subj_id = local_obj$subj_id,
    loglik_fn = local_obj$loglik_fn,
    working_prior = local_obj$working_prior,
    base_seed = local_obj$base_seed
  )
  refreshed <- build_local_bank(
    refreshed_local_obj,
    M_local = NULL,
    seed = seed,
    proposal_tag = "exact_refresh_gamma0"
  )

  list(
    bank = refreshed$bank,
    local_loglik_evals = c(initialization = 0L, enrichment = as.integer(n_ll), population_refreshes = 0L)
  )
}

enrich.local_bank <- function(x, phi, lambda, gaussian_map_fn,
                              rho_local = 0.5,
                              split_tol = 0.03,
                              max_topups = 2L,
                              topup_batch = NULL,
                              enable_exact_refresh = FALSE,
                              max_exact_refreshes = 1L,
                              exact_refresh_particles = NULL,
                              exact_refresh_max_rounds = NULL,
                              exact_refresh_n_mcmc_moves = 3L,
                              seed = NULL,
                              ...) {
  bank <- .validate_local_bank(x)
  counts <- .zero_local_loglik_evals()
  topup_points <- 0L
  enriched <- FALSE
  n_exact_refresh <- 0L

  max_topups <- as.integer(max(0L, max_topups))
  max_exact_refreshes <- as.integer(max(0L, max_exact_refreshes))
  if (max_topups <= 0L && max_exact_refreshes <= 0L) {
    return(list(bank = bank, local_loglik_evals = counts, n_new = topup_points, enriched = enriched))
  }

  for (step in seq_len(max_topups)) {
    diag_now <- .bank_group_diagnostics(bank, phi = phi, lambda = lambda, gaussian_map_fn = gaussian_map_fn)
    need_ess <- is.finite(diag_now$ess_q) && diag_now$ess_q < rho_local * diag_now$M
    need_split <- !is.finite(diag_now$split_q) || diag_now$split_q > split_tol
    if (!(need_ess || need_split)) break

    batch <- topup_batch
    if (is.null(batch)) {
      batch <- max(8L, min(32L, ceiling(0.5 * nrow(bank$alpha))))
    }
    top <- .topup_local_bank_from_prefit(
      bank = bank,
      n_new = batch,
      seed = if (is.null(seed)) NULL else as.integer(seed + step)
    )
    if (top$n_new <= 0L) break
    bank <- top$bank
    topup_points <- topup_points + top$n_new
    enriched <- TRUE
  }

  if (isTRUE(enable_exact_refresh) && max_exact_refreshes > 0L) {
    for (rr in seq_len(max_exact_refreshes)) {
      diag_now <- .bank_group_diagnostics(bank, phi = phi, lambda = lambda, gaussian_map_fn = gaussian_map_fn)
      need_ess <- is.finite(diag_now$ess_q) && diag_now$ess_q < rho_local * diag_now$M
      need_split <- !is.finite(diag_now$split_q) || diag_now$split_q > split_tol
      if (!(need_ess || need_split)) break

      refreshed <- .refresh_local_bank_exact(
        bank = bank,
        refresh_particles = exact_refresh_particles,
        refresh_max_rounds = exact_refresh_max_rounds,
        refresh_n_mcmc_moves = exact_refresh_n_mcmc_moves,
        seed = if (is.null(seed)) NULL else as.integer(seed + 100000L * rr)
      )
      bank <- refreshed$bank
      counts <- .add_local_loglik_evals(counts, refreshed$local_loglik_evals)
      enriched <- TRUE
      n_exact_refresh <- n_exact_refresh + 1L
    }
  }

  list(
    bank = bank,
    local_loglik_evals = counts,
    n_new = as.integer(topup_points),
    enriched = enriched,
    n_exact_refresh = as.integer(n_exact_refresh)
  )
}

predict_population_log_increment <- function(banks, phi, lambda_old, lambda_new, gaussian_map_fn,
                                             increment_cache = NULL) {
  cache <- increment_cache
  if (is.null(cache)) {
    cache <- .prepare_population_delta_cache(
      banks = banks,
      phi = phi,
      gaussian_map_fn = gaussian_map_fn,
      lambda_ref = lambda_old
    )
  }
  .predict_population_log_increment_from_cache(
    cache = cache,
    lambda_old = lambda_old,
    lambda_new = lambda_new
  )
}

.prepare_population_delta_cache <- function(banks, phi, gaussian_map_fn, lambda_ref = NULL) {
  phi <- as.matrix(phi)
  N <- nrow(phi)
  lapply(
    banks,
    function(bank) {
      bank <- .validate_local_bank(bank)
      M <- nrow(bank$alpha)
      delta_mat <- matrix(0, nrow = N, ncol = M)
      for (n in seq_len(N)) {
        delta_mat[n, ] <- .local_bank_delta(
          bank,
          phi = phi[n, , drop = TRUE],
          gaussian_map_fn = gaussian_map_fn,
          validate = FALSE
        )
      }
      cache <- list(
        bank = bank,
        delta = delta_mat,
        logw_rep = rep(bank$logw, each = N),
        lambda_ref = if (is.null(lambda_ref)) NULL else as.numeric(lambda_ref)
      )
      if (!is.null(lambda_ref)) {
        cache$logZ_ref <- .rowLogSumExp(as.numeric(lambda_ref) * delta_mat + cache$logw_rep)
      }
      cache
    }
  )
}

.logZ_matrix_from_cache <- function(cache, lambda) {
  lambda <- as.numeric(lambda)
  do.call(
    cbind,
    lapply(
      cache,
      function(cache_i) {
        if (!is.null(cache_i$logZ_ref) && !is.null(cache_i$lambda_ref) &&
            abs(lambda - cache_i$lambda_ref) < 1e-15) {
          return(as.numeric(cache_i$logZ_ref))
        }
        .rowLogSumExp(lambda * cache_i$delta + cache_i$logw_rep)
      }
    )
  )
}

.estimate_logZ_validated <- function(bank, phi, lambda, gaussian_map_fn) {
  delta <- .local_bank_delta(bank, phi = phi, gaussian_map_fn = gaussian_map_fn, validate = FALSE)
  .local_bank_weight_stats(bank, delta = delta, lambda = lambda)$logZ
}

.estimate_logZ_bank_list <- function(banks, phi, lambda, gaussian_map_fn) {
  vapply(
    banks,
    function(bank) .estimate_logZ_validated(bank, phi = phi, lambda = lambda, gaussian_map_fn = gaussian_map_fn),
    numeric(1L)
  )
}

.predict_population_log_increment_from_cache <- function(cache, lambda_old, lambda_new) {
  lambda_old <- as.numeric(lambda_old)
  lambda_new <- as.numeric(lambda_new)
  if (abs(lambda_new - lambda_old) < 1e-15) {
    return(rep.int(0, nrow(cache[[1L]]$delta)))
  }
  N <- nrow(cache[[1L]]$delta)
  out <- numeric(N)
  for (i in seq_along(cache)) {
    cache_i <- cache[[i]]
    logZ_old <- cache_i$logZ_ref
    if (is.null(logZ_old)) {
      logZ_old <- .rowLogSumExp(lambda_old * cache_i$delta + cache_i$logw_rep)
    }
    logZ_new <- .rowLogSumExp(lambda_new * cache_i$delta + cache_i$logw_rep)
    out <- out + (logZ_new - logZ_old)
  }
  out
}

.ess_from_logweights <- function(logw) {
  norm <- .normalize_logweights(logw)
  .ess_from_weights(norm$w)
}

.rcess_from_log_increment <- function(w, log_inc) {
  w <- as.numeric(w)
  w[!is.finite(w) | w < 0] <- 0
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) return(0)
  w <- w / sw

  log_inc <- as.numeric(log_inc)
  ok <- is.finite(log_inc)
  if (!any(ok)) return(0)
  log_inc[!ok] <- min(log_inc[ok])

  a1 <- logsumexp(log(w) + log_inc)
  a2 <- logsumexp(log(w) + 2 * log_inc)
  val <- exp(2 * a1 - a2)
  if (!is.finite(val)) return(0)
  pmin(pmax(val, 0), 1)
}

choose_next_lambda_nested <- function(phi, banks, w, lambda_curr, rho_step,
                                      gaussian_map_fn, step_target_fn = cess_target_at_lambda,
                                      increment_cache = NULL,
                                      tol = 1e-4, max_iter = 30L) {
  if (is.null(increment_cache)) {
    increment_cache <- .prepare_population_delta_cache(
      banks = banks,
      phi = phi,
      gaussian_map_fn = gaussian_map_fn,
      lambda_ref = lambda_curr
    )
  }
  rcess_at <- function(lambda_new) {
    log_inc <- .predict_population_log_increment_from_cache(
      cache = increment_cache,
      lambda_old = lambda_curr,
      lambda_new = lambda_new
    )
    .rcess_from_log_increment(w, log_inc)
  }

  target_frac <- if (is.null(rho_step)) {
    as.numeric(step_target_fn(lambda_curr))
  } else {
    as.numeric(rho_step)
  }
  target_frac <- pmin(pmax(target_frac, 0), 0.999)

  rcess_one <- rcess_at(1.0)
  if (rcess_one >= target_frac) return(1.0)

  lo <- lambda_curr
  hi <- 1.0
  for (iter in seq_len(max_iter)) {
    mid <- 0.5 * (lo + hi)
    rcess_mid <- rcess_at(mid)
    if (rcess_mid >= target_frac) {
      lo <- mid
    } else {
      hi <- mid
    }
    if ((hi - lo) < tol) break
  }

  min(1.0, max(lo, lambda_curr + min(1e-4, 1.0 - lambda_curr)))
}

.bank_ess_summary <- function(banks, phi, lambda, gaussian_map_fn, probs = c(0.1, 0.5),
                              increment_cache = NULL) {
  cache <- increment_cache
  if (is.null(cache)) {
    cache <- .prepare_population_delta_cache(
      banks = banks,
      phi = phi,
      gaussian_map_fn = gaussian_map_fn,
      lambda_ref = NULL
    )
  }
  if (!length(cache)) {
    return(list(min = NA_real_, quantiles = setNames(rep(NA_real_, length(probs)), paste0("q", probs * 100))))
  }

  N <- nrow(cache[[1L]]$delta)
  ess_mat <- vapply(
    cache,
    function(cache_i) .bank_ess_values_from_delta(cache_i$bank, delta_mat = cache_i$delta, lambda = lambda),
    numeric(N)
  )
  ess_vals <- as.numeric(ess_mat)
  if (!length(ess_vals)) {
    return(list(min = NA_real_, quantiles = setNames(rep(NA_real_, length(probs)), paste0("q", probs * 100))))
  }
  qq <- stats::quantile(ess_vals, probs = probs, names = FALSE, na.rm = TRUE)
  list(
    min = min(ess_vals, na.rm = TRUE),
    quantiles = setNames(as.numeric(qq), paste0("q", probs * 100))
  )
}

.resample_population_particles <- function(phi, w, logprior) {
  idx <- .resample_indices(w)
  list(
    phi = phi[idx, , drop = FALSE],
    w = rep(1 / nrow(phi), nrow(phi)),
    logprior = logprior[idx],
    idx = idx
  )
}

mutate_population_particle_nested <- function(phi_row, logprior_row, lambda,
                                              banks, logprior_phi, gaussian_map_fn, rw_cov,
                                              current_logZ = NULL,
                                              seed = NULL) {
  if (!is.null(seed)) set.seed(as.integer(seed))
  prop <- as.numeric(phi_row + mvtnorm::rmvnorm(1L, sigma = rw_cov))
  lp_prop <- logprior_phi(prop)
  if (!is.finite(lp_prop)) {
    return(list(phi = phi_row, logprior = logprior_row, accepted = FALSE, logZ = current_logZ))
  }

  log_ratio <- lp_prop - logprior_row
  prop_logZ <- .estimate_logZ_bank_list(banks, phi = prop, lambda = lambda, gaussian_map_fn = gaussian_map_fn)
  curr_logZ <- current_logZ
  if (is.null(curr_logZ)) {
    curr_logZ <- .estimate_logZ_bank_list(banks, phi = phi_row, lambda = lambda, gaussian_map_fn = gaussian_map_fn)
  }
  log_ratio <- log_ratio + sum(prop_logZ - curr_logZ)

  accepted <- is.finite(log_ratio) && (log(runif(1L)) < min(0, log_ratio))
  if (!accepted) {
    return(list(phi = phi_row, logprior = logprior_row, accepted = FALSE, logZ = curr_logZ))
  }

  list(phi = prop, logprior = lp_prop, accepted = TRUE, logZ = prop_logZ)
}

nested_population_smc <- function(rprior_phi,
                                  local_objs,
                                  logprior_phi,
                                  gaussian_map_fn,
                                  N,
                                  M_local,
                                  rho_step = NULL,
                                  rho_res = 0.5,
                                  rho_local = 0.5,
                                  n_population_moves = 1L,
                                  n_local_moves = 1L,
                                  n_population_refresh_moves = 0L,
                                  max_rounds = 100L,
                                  max_bank_topups = 2L,
                                  bank_split_tol = 0.03,
                                  bank_topup_batch = NULL,
                                  enable_exact_local_refresh = FALSE,
                                  max_exact_local_refreshes = 1L,
                                  exact_refresh_particles = NULL,
                                  exact_refresh_max_rounds = NULL,
                                  exact_refresh_n_mcmc_moves = 3L,
                                  step_target_fn = cess_target_at_lambda,
                                  seed = 123,
                                  verbose = TRUE) {
  vcat <- function(...) if (isTRUE(verbose)) cat(...)
  set.seed(as.integer(seed))

  N <- as.integer(N)
  phi <- as.matrix(rprior_phi(N))
  if (nrow(phi) != N) stop("rprior_phi returned wrong number of particles.")
  w <- rep(1 / N, N)
  logprior <- apply(phi, 1L, logprior_phi)

  bank_init <- lapply(
    seq_along(local_objs),
    function(i) build_local_bank(local_objs[[i]], M_local = M_local, seed = as.integer(seed + 1009L * i))
  )
  banks <- lapply(bank_init, `[[`, "bank")
  local_loglik_evals <- Reduce(
    .add_local_loglik_evals,
    lapply(bank_init, `[[`, "local_loglik_evals"),
    init = .zero_local_loglik_evals()
  )

  lambda <- 0
  round <- 0L
  lambda_hist <- c(0)
  ess_hist <- c(N)
  bank_ess_min_hist <- c(M_local)
  bank_ess_q10_hist <- c(M_local)
  bank_size_mean_hist <- c(mean(vapply(banks, function(b) nrow(b$alpha), numeric(1))))
  bank_topup_points <- 0L
  n_exact_local_refresh <- 0L
  n_population_accept <- 0L

  if (as.integer(n_local_moves) > 0L || as.integer(n_population_refresh_moves) > 0L) {
    vcat("Shared-bank path: per-particle local refresh moves are disabled; adequacy uses adaptive shared-bank top-ups instead.\n")
  }

  while (lambda < 1 - 1e-12 && round < max_rounds) {
    round <- round + 1L
    increment_cache <- .prepare_population_delta_cache(
      banks = banks,
      phi = phi,
      gaussian_map_fn = gaussian_map_fn,
      lambda_ref = lambda
    )
    lambda_next <- choose_next_lambda_nested(
      phi = phi,
      banks = banks,
      w = w,
      lambda_curr = lambda,
      rho_step = rho_step,
      gaussian_map_fn = gaussian_map_fn,
      step_target_fn = step_target_fn,
      increment_cache = increment_cache
    )
    lambda_next <- max(lambda_next, lambda)
    vcat(sprintf("Round %d: lambda %.4f -> %.4f\n", round, lambda, lambda_next))

    log_w_inc <- predict_population_log_increment(
      banks = banks,
      phi = phi,
      lambda_old = lambda,
      lambda_new = lambda_next,
      gaussian_map_fn = gaussian_map_fn,
      increment_cache = increment_cache
    )
    norm <- .normalize_logweights(log(pmax(w, .Machine$double.eps)) + log_w_inc)
    w <- norm$w
    lambda <- lambda_next
    lambda_hist <- c(lambda_hist, lambda)
    ess_now <- .ess_from_weights(w)
    ess_hist <- c(ess_hist, ess_now)

    bank_diag <- .bank_ess_summary(
      banks,
      phi = phi,
      lambda = lambda,
      gaussian_map_fn = gaussian_map_fn,
      increment_cache = increment_cache
    )
    bank_ess_min_hist <- c(bank_ess_min_hist, bank_diag$min)
    bank_ess_q10_hist <- c(bank_ess_q10_hist, unname(bank_diag$quantiles["q10"]))
    banks_mut <- lapply(increment_cache, `[[`, "bank")
    current_logZ_mat <- .logZ_matrix_from_cache(increment_cache, lambda = lambda)

    if (ess_now < rho_res * N) {
      rs <- .resample_population_particles(phi, w, logprior)
      phi <- rs$phi
      w <- rs$w
      logprior <- rs$logprior
      current_logZ_mat <- current_logZ_mat[rs$idx, , drop = FALSE]
      ess_now <- .ess_from_weights(w)
      ess_hist[length(ess_hist)] <- ess_now
      vcat(sprintf("  population resample at ESS=%.1f\n", ess_now))
    }

    S <- tryCatch(weighted_cov(phi, w), error = function(e) stats::cov(phi))
    S <- as.matrix(Matrix::nearPD(S, conv.tol = 1e-7)$mat)
    diag(S) <- pmax(diag(S), 1e-8)
    rw_cov <- (0.25^2 / max(1L, ncol(phi))) * S + diag(1e-10, ncol(phi))

    for (mv in seq_len(as.integer(n_population_moves))) {
      for (n in seq_len(N)) {
        out <- mutate_population_particle_nested(
          phi_row = phi[n, ],
          logprior_row = logprior[n],
          lambda = lambda,
          banks = banks_mut,
          logprior_phi = logprior_phi,
          gaussian_map_fn = gaussian_map_fn,
          rw_cov = rw_cov,
          current_logZ = current_logZ_mat[n, ],
          seed = as.integer(seed + 500000L * round + 1009L * mv + n)
        )
        phi[n, ] <- out$phi
        logprior[n] <- out$logprior
        current_logZ_mat[n, ] <- out$logZ
        n_population_accept <- n_population_accept + as.integer(out$accepted)
      }
    }

    round_topups <- 0L
    for (i in seq_along(banks)) {
      up <- enrich(
        banks[[i]],
        phi = phi,
        lambda = lambda,
        gaussian_map_fn = gaussian_map_fn,
        rho_local = rho_local,
        split_tol = bank_split_tol,
        max_topups = max_bank_topups,
        topup_batch = bank_topup_batch,
        enable_exact_refresh = enable_exact_local_refresh,
        max_exact_refreshes = max_exact_local_refreshes,
        exact_refresh_particles = exact_refresh_particles,
        exact_refresh_max_rounds = exact_refresh_max_rounds,
        exact_refresh_n_mcmc_moves = exact_refresh_n_mcmc_moves,
        seed = as.integer(seed + 900000L * round + i)
      )
      banks[[i]] <- up$bank
      local_loglik_evals <- .add_local_loglik_evals(local_loglik_evals, up$local_loglik_evals)
      round_topups <- round_topups + up$n_new
      n_exact_local_refresh <- n_exact_local_refresh + as.integer(up$n_exact_refresh %||% 0L)
    }
    bank_topup_points <- bank_topup_points + round_topups
    bank_diag_post <- .bank_ess_summary(banks, phi = phi, lambda = lambda, gaussian_map_fn = gaussian_map_fn)
    bank_ess_min_hist[length(bank_ess_min_hist)] <- bank_diag_post$min
    bank_ess_q10_hist[length(bank_ess_q10_hist)] <- unname(bank_diag_post$quantiles["q10"])
    bank_size_mean_hist <- c(bank_size_mean_hist, mean(vapply(banks, function(b) nrow(b$alpha), numeric(1))))
    if (round_topups > 0L) {
      vcat(sprintf("  bank top-ups added %d cached local points\n", round_topups))
    }
    if (isTRUE(enable_exact_local_refresh) && n_exact_local_refresh > 0L) {
      vcat(sprintf("  exact local refreshes so far: %d\n", n_exact_local_refresh))
    }

    if (lambda >= 1 - 1e-12) break
  }

  local_loglik_evals <- .finalize_local_loglik_evals(local_loglik_evals)

  list(
    phi = phi,
    w = w,
    local_banks = banks,
    meta = list(
      rounds = round,
      lambda_hist = lambda_hist,
      ess_hist = ess_hist,
      bank_ess_min_hist = bank_ess_min_hist,
      bank_ess_q10_hist = bank_ess_q10_hist,
      bank_size_mean_hist = bank_size_mean_hist,
      bank_topup_points = as.integer(bank_topup_points),
      n_exact_local_refresh = as.integer(n_exact_local_refresh),
      local_loglik_evals = local_loglik_evals,
      n_exact_loglik = as.integer(local_loglik_evals["total"]),
      n_population_accept = as.integer(n_population_accept),
      n_local_accept = 0L,
      implementation = if (isTRUE(enable_exact_local_refresh)) "shared_local_banks_with_exact_refresh" else "shared_local_banks_baseline",
      ignored_args = list(
        n_local_moves = as.integer(n_local_moves),
        n_population_refresh_moves = as.integer(n_population_refresh_moves)
      )
    )
  )
}
