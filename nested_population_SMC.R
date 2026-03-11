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

.new_local_bank <- function(alpha, w, logg, local_obj, loglik = NULL, proposal_tag = "prefit_gamma0") {
  bank <- structure(
    list(
      alpha = as.matrix(alpha),
      w = as.numeric(w),
      logg = as.numeric(logg),
      loglik = if (is.null(loglik)) NULL else as.numeric(loglik),
      d_theta = ncol(as.matrix(alpha)),
      local_obj = local_obj,
      proposal_tag = as.character(proposal_tag)
    ),
    class = "local_bank"
  )
  .validate_local_bank(bank)
}

.validate_local_bank <- function(bank) {
  stopifnot(inherits(bank, "local_bank"))
  A <- as.matrix(bank$alpha)
  n <- nrow(A)
  if (length(bank$w) != n) stop("local_bank has inconsistent weight length.")
  if (length(bank$logg) != n) stop("local_bank has inconsistent logg length.")
  if (!is.null(bank$loglik) && length(bank$loglik) != n) stop("local_bank has inconsistent loglik length.")
  sw <- sum(bank$w)
  if (!is.finite(sw) || sw <= 0) stop("local_bank weights are invalid.")
  bank$w <- bank$w / sw
  bank
}

.local_bank_delta <- function(bank, phi, gaussian_map_fn) {
  bank <- .validate_local_bank(bank)
  logp <- .log_prior_gauss_mat(bank$alpha, gaussian_map_fn(as.numeric(phi), bank$d_theta))
  as.numeric(logp - bank$logg)
}

build_local_bank <- function(local_obj, M_local = NULL, seed = 123L) {
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
      local_obj = local_obj
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
  delta <- .local_bank_delta(x, phi = phi, gaussian_map_fn = gaussian_map_fn)
  logsumexp(log(pmax(x$w, .Machine$double.eps)) + as.numeric(lambda) * delta)
}

predict_increment.local_bank <- function(x, phi, lambda_old, lambda_new, gaussian_map_fn, ...) {
  lambda_old <- as.numeric(lambda_old)
  lambda_new <- as.numeric(lambda_new)
  if (lambda_new < lambda_old - 1e-12) {
    stop("predict_increment.local_bank does not support decreasing lambda.")
  }
  logZ_old <- estimate_logZ(x, phi = phi, lambda = lambda_old, gaussian_map_fn = gaussian_map_fn)
  logZ_new <- estimate_logZ(x, phi = phi, lambda = lambda_new, gaussian_map_fn = gaussian_map_fn)
  list(
    log_u = as.numeric(logZ_new - logZ_old),
    logZ_old = as.numeric(logZ_old),
    logZ_new = as.numeric(logZ_new),
    ess_new = ess_under(x, phi = phi, lambda = lambda_new, gaussian_map_fn = gaussian_map_fn)
  )
}

ess_under.local_bank <- function(x, phi, lambda, gaussian_map_fn, ...) {
  x <- .validate_local_bank(x)
  delta <- .local_bank_delta(x, phi = phi, gaussian_map_fn = gaussian_map_fn)
  norm <- .normalize_logweights(log(pmax(x$w, .Machine$double.eps)) + as.numeric(lambda) * delta)
  .ess_from_weights(norm$w)
}

enrich.local_bank <- function(x, ...) {
  list(
    bank = .validate_local_bank(x),
    local_loglik_evals = .zero_local_loglik_evals(),
    n_new = 0L,
    enriched = FALSE
  )
}

predict_population_log_increment <- function(banks, phi, lambda_old, lambda_new, gaussian_map_fn) {
  N <- nrow(phi)
  out <- numeric(N)
  for (n in seq_len(N)) {
    acc <- 0
    for (i in seq_along(banks)) {
      acc <- acc + predict_increment(
        banks[[i]],
        phi = phi[n, , drop = TRUE],
        lambda_old = lambda_old,
        lambda_new = lambda_new,
        gaussian_map_fn = gaussian_map_fn
      )$log_u
    }
    out[n] <- acc
  }
  out
}

.ess_from_logweights <- function(logw) {
  norm <- .normalize_logweights(logw)
  .ess_from_weights(norm$w)
}

choose_next_lambda_nested <- function(phi, banks, w, lambda_curr, rho_step,
                                      gaussian_map_fn, tol = 1e-4, max_iter = 30L) {
  target_ess <- rho_step * nrow(phi)
  ess_at <- function(lambda_new) {
    log_inc <- predict_population_log_increment(
      banks = banks,
      phi = phi,
      lambda_old = lambda_curr,
      lambda_new = lambda_new,
      gaussian_map_fn = gaussian_map_fn
    )
    .ess_from_logweights(log(pmax(w, .Machine$double.eps)) + log_inc)
  }

  ess_one <- ess_at(1.0)
  if (ess_one >= target_ess) return(1.0)

  lo <- lambda_curr
  hi <- 1.0
  for (iter in seq_len(max_iter)) {
    mid <- 0.5 * (lo + hi)
    ess_mid <- ess_at(mid)
    if (ess_mid >= target_ess) {
      lo <- mid
    } else {
      hi <- mid
    }
    if ((hi - lo) < tol) break
  }
  lo
}

.bank_ess_summary <- function(banks, phi, lambda, gaussian_map_fn, probs = c(0.1, 0.5)) {
  ess_vals <- numeric(0)
  for (i in seq_along(banks)) {
    for (n in seq_len(nrow(phi))) {
      ess_vals <- c(
        ess_vals,
        ess_under(
          banks[[i]],
          phi = phi[n, , drop = TRUE],
          lambda = lambda,
          gaussian_map_fn = gaussian_map_fn
        )
      )
    }
  }
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
                                              seed = NULL) {
  if (!is.null(seed)) set.seed(as.integer(seed))
  prop <- as.numeric(phi_row + mvtnorm::rmvnorm(1L, sigma = rw_cov))
  lp_prop <- logprior_phi(prop)
  if (!is.finite(lp_prop)) {
    return(list(phi = phi_row, logprior = logprior_row, accepted = FALSE))
  }

  log_ratio <- lp_prop - logprior_row
  for (i in seq_along(banks)) {
    log_ratio <- log_ratio +
      estimate_logZ(banks[[i]], phi = prop, lambda = lambda, gaussian_map_fn = gaussian_map_fn) -
      estimate_logZ(banks[[i]], phi = phi_row, lambda = lambda, gaussian_map_fn = gaussian_map_fn)
  }

  accepted <- is.finite(log_ratio) && (log(runif(1L)) < min(0, log_ratio))
  if (!accepted) {
    return(list(phi = phi_row, logprior = logprior_row, accepted = FALSE))
  }

  list(phi = prop, logprior = lp_prop, accepted = TRUE)
}

nested_population_smc <- function(rprior_phi,
                                  local_objs,
                                  logprior_phi,
                                  gaussian_map_fn,
                                  N,
                                  M_local,
                                  rho_step = 0.8,
                                  rho_res = 0.5,
                                  rho_local = 0.5,
                                  n_population_moves = 1L,
                                  n_local_moves = 1L,
                                  n_population_refresh_moves = 0L,
                                  max_rounds = 100L,
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
  n_population_accept <- 0L

  if (as.integer(n_local_moves) > 0L || as.integer(n_population_refresh_moves) > 0L) {
    vcat("Shared-bank baseline: local refresh moves are currently disabled; reweighting uses cached local banks only.\n")
  }

  while (lambda < 1 - 1e-12 && round < max_rounds) {
    round <- round + 1L
    lambda_next <- choose_next_lambda_nested(
      phi = phi,
      banks = banks,
      w = w,
      lambda_curr = lambda,
      rho_step = rho_step,
      gaussian_map_fn = gaussian_map_fn
    )
    lambda_next <- max(lambda_next, lambda)
    vcat(sprintf("Round %d: lambda %.4f -> %.4f\n", round, lambda, lambda_next))

    log_w_inc <- predict_population_log_increment(
      banks = banks,
      phi = phi,
      lambda_old = lambda,
      lambda_new = lambda_next,
      gaussian_map_fn = gaussian_map_fn
    )
    norm <- .normalize_logweights(log(pmax(w, .Machine$double.eps)) + log_w_inc)
    w <- norm$w
    lambda <- lambda_next
    lambda_hist <- c(lambda_hist, lambda)
    ess_now <- .ess_from_weights(w)
    ess_hist <- c(ess_hist, ess_now)

    bank_diag <- .bank_ess_summary(banks, phi = phi, lambda = lambda, gaussian_map_fn = gaussian_map_fn)
    bank_ess_min_hist <- c(bank_ess_min_hist, bank_diag$min)
    bank_ess_q10_hist <- c(bank_ess_q10_hist, unname(bank_diag$quantiles["q10"]))

    if (ess_now < rho_res * N) {
      rs <- .resample_population_particles(phi, w, logprior)
      phi <- rs$phi
      w <- rs$w
      logprior <- rs$logprior
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
          banks = banks,
          logprior_phi = logprior_phi,
          gaussian_map_fn = gaussian_map_fn,
          rw_cov = rw_cov,
          seed = as.integer(seed + 500000L * round + 1009L * mv + n)
        )
        phi[n, ] <- out$phi
        logprior[n] <- out$logprior
        n_population_accept <- n_population_accept + as.integer(out$accepted)
      }
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
      local_loglik_evals = local_loglik_evals,
      n_exact_loglik = as.integer(local_loglik_evals["total"]),
      n_population_accept = as.integer(n_population_accept),
      n_local_accept = 0L,
      implementation = "shared_local_banks_baseline",
      ignored_args = list(
        rho_local = rho_local,
        n_local_moves = as.integer(n_local_moves),
        n_population_refresh_moves = as.integer(n_population_refresh_moves)
      )
    )
  )
}
