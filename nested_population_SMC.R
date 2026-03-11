#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(mvtnorm)
})

source("smc_core.R")
source("new_SMC_cache.R")

.phi_equal <- function(a, b) {
  if (is.null(a) || is.null(b)) return(FALSE)
  length(a) == length(b) && isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = 0))
}

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
  setNames(
    rep.int(0L, 4L),
    c("initialization", "bridge_updates", "local_rejuvenation", "population_refreshes")
  )
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

.new_local_module <- function(alpha, w, lambda, phi, loglik, logg, logp, local_obj, init_log_norm = NA_real_) {
  module <- structure(
    list(
      alpha = as.matrix(alpha),
      w = as.numeric(w),
      lambda = as.numeric(lambda),
      phi = if (is.null(phi)) NULL else as.numeric(phi),
      loglik = as.numeric(loglik),
      logg = as.numeric(logg),
      logp = if (is.null(logp)) NULL else as.numeric(logp),
      d_theta = ncol(as.matrix(alpha)),
      local_obj = local_obj,
      init_log_norm = as.numeric(init_log_norm)
    ),
    class = "local_module"
  )
  .validate_local_module(module)
}

.validate_local_module <- function(module) {
  stopifnot(inherits(module, "local_module"))
  A <- as.matrix(module$alpha)
  n <- nrow(A)
  if (length(module$w) != n) stop("local_module has inconsistent weight length.")
  if (length(module$loglik) != n) stop("local_module has inconsistent loglik length.")
  if (length(module$logg) != n) stop("local_module has inconsistent logg length.")
  if (!is.null(module$logp) && length(module$logp) != n) stop("local_module has inconsistent logp length.")
  sw <- sum(module$w)
  if (!is.finite(sw) || sw <= 0) stop("local_module weights are invalid.")
  module$w <- module$w / sw
  module
}

.attach_population_state <- function(module, phi, gaussian_map_fn) {
  stopifnot(inherits(module, "local_module"))
  module$phi <- as.numeric(phi)
  module$logp <- .log_prior_gauss_mat(module$alpha, gaussian_map_fn(module$phi, module$d_theta))
  .validate_local_module(module)
}

.local_target_log <- function(module) {
  module$loglik + module$lambda * module$logp + (1 - module$lambda) * module$logg
}

.build_local_module_prototype <- function(local_obj, M_local, seed = NULL) {
  prop <- make_local_theta_proposal_init(local_obj)
  Theta <- as.matrix(prop$draw(as.integer(M_local), seed = seed))
  loglik <- ll_parallel(Theta, local_obj$data, local_obj$loglik_fn, n_cores = 1L)
  logg <- log_working_prior_gaussian_mat(Theta, local_obj$working_prior)
  logq <- prop$log_q(Theta)
  norm <- .normalize_logweights(loglik + logg - logq)
  module <- .new_local_module(
    alpha = Theta,
    w = norm$w,
    lambda = 0,
    phi = NULL,
    loglik = loglik,
    logg = logg,
    logp = NULL,
    local_obj = local_obj,
    init_log_norm = norm$log_norm
  )
  counts <- .zero_local_loglik_evals()
  counts["initialization"] <- as.integer(nrow(Theta))
  list(module = module, local_loglik_evals = counts)
}

build_local_module <- function(local_obj, phi, lambda = 0, M_local, seed = 123, gaussian_map_fn) {
  if (!isTRUE(all.equal(lambda, 0))) {
    stop("Phase 1 build_local_module only supports lambda = 0 initialization.")
  }
  init <- .build_local_module_prototype(local_obj, M_local = M_local, seed = seed)
  .attach_population_state(init$module, phi = phi, gaussian_map_fn = gaussian_map_fn)
}

print.local_module <- function(x, ...) {
  cat(sprintf(
    "<local_module M=%d lambda=%.4f has_phi=%s>\n",
    nrow(x$alpha),
    x$lambda,
    if (is.null(x$phi)) "FALSE" else "TRUE"
  ))
  invisible(x)
}

predict_increment <- function(x, ...) UseMethod("predict_increment")
advance <- function(x, ...) UseMethod("advance")
move_population <- function(x, ...) UseMethod("move_population")
clone <- function(x, ...) UseMethod("clone")

clone.local_module <- function(x, ...) {
  .new_local_module(
    alpha = x$alpha,
    w = x$w,
    lambda = x$lambda,
    phi = x$phi,
    loglik = x$loglik,
    logg = x$logg,
    logp = x$logp,
    local_obj = x$local_obj,
    init_log_norm = x$init_log_norm
  )
}

predict_increment.local_module <- function(x, phi, lambda_new, gaussian_map_fn, ...) {
  x <- .validate_local_module(x)
  lambda_old <- as.numeric(x$lambda)
  lambda_new <- as.numeric(lambda_new)
  if (lambda_new < lambda_old - 1e-12) {
    stop("predict_increment.local_module does not support decreasing lambda.")
  }

  phi <- as.numeric(phi)
  logp_new <- if (.phi_equal(phi, x$phi)) {
    x$logp
  } else {
    .log_prior_gauss_mat(x$alpha, gaussian_map_fn(phi, x$d_theta))
  }

  inc <- lambda_new * logp_new - lambda_old * x$logp - (lambda_new - lambda_old) * x$logg
  norm <- .normalize_logweights(log(pmax(x$w, .Machine$double.eps)) + inc)
  list(
    log_u = norm$log_norm,
    w_new = norm$w,
    ess_new = .ess_from_weights(norm$w),
    logp_new = as.numeric(logp_new),
    inc = as.numeric(inc)
  )
}

.rejuvenate_local_module <- function(module, n_local_moves, gaussian_map_fn, seed = NULL, eval_bucket) {
  counts <- .zero_local_loglik_evals()
  n_accept <- 0L
  if (as.integer(n_local_moves) <= 0L) {
    return(list(module = module, local_loglik_evals = counts, n_accept = n_accept))
  }
  if (!is.null(seed)) set.seed(as.integer(seed))

  n_particles <- nrow(module$alpha)
  prop <- make_local_theta_proposal_refresh(
    phi = module$phi,
    local_obj = module$local_obj,
    callbacks = list(gaussian_map_fn = gaussian_map_fn)
  )

  for (mv in seq_len(as.integer(n_local_moves))) {
    curr <- module$alpha
    cand <- as.matrix(prop$draw(n_particles))
    q_curr <- prop$log_q(curr)
    q_cand <- prop$log_q(cand)
    loglik_cand <- ll_parallel(cand, module$local_obj$data, module$local_obj$loglik_fn, n_cores = 1L)
    logg_cand <- log_working_prior_gaussian_mat(cand, module$local_obj$working_prior)
    logp_cand <- .log_prior_gauss_mat(cand, gaussian_map_fn(module$phi, module$d_theta))
    logtar_curr <- .local_target_log(module)
    logtar_cand <- as.numeric(loglik_cand + module$lambda * logp_cand + (1 - module$lambda) * logg_cand)
    loga <- logtar_cand - logtar_curr + q_curr - q_cand
    acc <- is.finite(loga) & (log(runif(n_particles)) < pmin(0, loga))
    if (any(acc)) {
      module$alpha[acc, ] <- cand[acc, , drop = FALSE]
      module$loglik[acc] <- as.numeric(loglik_cand[acc])
      module$logg[acc] <- as.numeric(logg_cand[acc])
      module$logp[acc] <- as.numeric(logp_cand[acc])
      n_accept <- n_accept + sum(acc)
    }
    counts[eval_bucket] <- counts[eval_bucket] + as.integer(n_particles)
  }

  list(module = .validate_local_module(module), local_loglik_evals = counts, n_accept = n_accept)
}

.advance_local_module_impl <- function(module, phi, lambda_new, rho_local = 0.5,
                                       n_local_moves = 1L,
                                       gaussian_map_fn,
                                       seed = NULL,
                                       eval_bucket = "local_rejuvenation") {
  stopifnot(eval_bucket %in% names(.zero_local_loglik_evals()))
  module <- .validate_local_module(module)
  if (!is.null(seed)) set.seed(as.integer(seed))

  pred <- predict_increment(module, phi = phi, lambda_new = lambda_new, gaussian_map_fn = gaussian_map_fn)
  module$w <- pred$w_new
  module$phi <- as.numeric(phi)
  module$lambda <- as.numeric(lambda_new)
  module$logp <- pred$logp_new

  n_particles <- nrow(module$alpha)
  resampled <- FALSE
  ess_now <- pred$ess_new
  if (ess_now < rho_local * n_particles) {
    idx <- .resample_indices(module$w)
    module$alpha <- module$alpha[idx, , drop = FALSE]
    module$loglik <- module$loglik[idx]
    module$logg <- module$logg[idx]
    module$logp <- module$logp[idx]
    module$w <- rep(1 / n_particles, n_particles)
    ess_now <- n_particles
    resampled <- TRUE
  }

  rej <- .rejuvenate_local_module(
    module = module,
    n_local_moves = n_local_moves,
    gaussian_map_fn = gaussian_map_fn,
    seed = if (is.null(seed)) NULL else seed + 17L,
    eval_bucket = eval_bucket
  )

  list(
    module = rej$module,
    log_u = pred$log_u,
    ess = ess_now,
    resampled = resampled,
    local_loglik_evals = rej$local_loglik_evals,
    n_accept = rej$n_accept
  )
}

advance.local_module <- function(x, phi, lambda_new, rho_local = 0.5,
                                 n_local_moves = 1L,
                                 gaussian_map_fn,
                                 seed = NULL,
                                 ...) {
  .advance_local_module_impl(
    module = x,
    phi = phi,
    lambda_new = lambda_new,
    rho_local = rho_local,
    n_local_moves = n_local_moves,
    gaussian_map_fn = gaussian_map_fn,
    seed = seed,
    eval_bucket = "local_rejuvenation"
  )
}

move_population.local_module <- function(x, phi_new, n_local_moves = 0L,
                                         rho_local = 0.5,
                                         gaussian_map_fn,
                                         seed = NULL,
                                         ...) {
  .advance_local_module_impl(
    module = x,
    phi = phi_new,
    lambda_new = x$lambda,
    rho_local = rho_local,
    n_local_moves = n_local_moves,
    gaussian_map_fn = gaussian_map_fn,
    seed = seed,
    eval_bucket = "population_refreshes"
  )
}

# Compatibility wrappers retained for the current Phase 1 code path.
clone_local_module <- function(module) clone(module)
predict_increment_local_module <- function(module, ...) predict_increment(module, ...)
advance_local_module <- function(module, ...) advance(module, ...)
move_population_local_module <- function(module, ...) move_population(module, ...)

predict_population_log_increment <- function(modules, phi, lambda_new, gaussian_map_fn) {
  N <- nrow(phi)
  out <- numeric(N)
  for (n in seq_len(N)) {
    acc <- 0
    for (i in seq_along(modules[[n]])) {
      acc <- acc + predict_increment(
        modules[[n]][[i]],
        phi = phi[n, , drop = TRUE],
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

choose_next_lambda_nested <- function(phi, modules, w, lambda_curr, rho_step,
                                      gaussian_map_fn, tol = 1e-4, max_iter = 30L) {
  target_ess <- rho_step * nrow(phi)
  ess_at <- function(lambda_new) {
    log_inc <- predict_population_log_increment(modules, phi, lambda_new, gaussian_map_fn)
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

resample_population_modules <- function(phi, modules, w, logprior) {
  idx <- .resample_indices(w)
  list(
    phi = phi[idx, , drop = FALSE],
    modules = lapply(idx, function(k) lapply(modules[[k]], clone)),
    w = rep(1 / nrow(phi), nrow(phi)),
    logprior = logprior[idx],
    idx = idx
  )
}

mutate_population_particle_nested <- function(phi_row, modules_row, logprior_row, lambda,
                                              logprior_phi, gaussian_map_fn, rw_cov,
                                              rho_local = 0.5,
                                              n_population_refresh_moves = 0L,
                                              seed = NULL) {
  if (!is.null(seed)) set.seed(as.integer(seed))
  prop <- as.numeric(phi_row + mvtnorm::rmvnorm(1L, sigma = rw_cov))
  lp_prop <- logprior_phi(prop)
  if (!is.finite(lp_prop)) {
    return(list(
      phi = phi_row,
      modules = modules_row,
      logprior = logprior_row,
      accepted = FALSE,
      local_loglik_evals = .zero_local_loglik_evals(),
      n_local_accept = 0L
    ))
  }

  log_ratio <- lp_prop - logprior_row
  for (i in seq_along(modules_row)) {
    log_ratio <- log_ratio + predict_increment(
      modules_row[[i]],
      phi = prop,
      lambda_new = lambda,
      gaussian_map_fn = gaussian_map_fn
    )$log_u
  }

  if (!is.finite(log_ratio) || log(runif(1L)) >= min(0, log_ratio)) {
    return(list(
      phi = phi_row,
      modules = modules_row,
      logprior = logprior_row,
      accepted = FALSE,
      local_loglik_evals = .zero_local_loglik_evals(),
      n_local_accept = 0L
    ))
  }

  modules_new <- modules_row
  counts <- .zero_local_loglik_evals()
  n_local_accept <- 0L
  for (i in seq_along(modules_new)) {
    mv <- move_population(
      modules_new[[i]],
      phi_new = prop,
      n_local_moves = n_population_refresh_moves,
      rho_local = rho_local,
      gaussian_map_fn = gaussian_map_fn,
      seed = if (is.null(seed)) NULL else seed + i
    )
    modules_new[[i]] <- mv$module
    counts <- .add_local_loglik_evals(counts, mv$local_loglik_evals)
    n_local_accept <- n_local_accept + mv$n_accept
  }

  list(
    phi = prop,
    modules = modules_new,
    logprior = lp_prop,
    accepted = TRUE,
    local_loglik_evals = counts,
    n_local_accept = n_local_accept
  )
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

  local_init <- lapply(
    seq_along(local_objs),
    function(i) .build_local_module_prototype(local_objs[[i]], M_local = M_local, seed = as.integer(seed + 1009L * i))
  )
  local_loglik_evals <- Reduce(
    .add_local_loglik_evals,
    lapply(local_init, `[[`, "local_loglik_evals"),
    init = .zero_local_loglik_evals()
  )

  modules <- vector("list", N)
  for (n in seq_len(N)) {
    modules[[n]] <- lapply(seq_along(local_objs), function(i) {
      .attach_population_state(clone(local_init[[i]]$module), phi = phi[n, , drop = TRUE], gaussian_map_fn = gaussian_map_fn)
    })
  }

  lambda <- 0
  round <- 0L
  lambda_hist <- c(0)
  ess_hist <- c(N)
  n_population_accept <- 0L
  n_local_accept <- 0L

  while (lambda < 1 - 1e-12 && round < max_rounds) {
    round <- round + 1L
    lambda_next <- choose_next_lambda_nested(
      phi = phi,
      modules = modules,
      w = w,
      lambda_curr = lambda,
      rho_step = rho_step,
      gaussian_map_fn = gaussian_map_fn
    )
    lambda_next <- max(lambda_next, lambda)
    vcat(sprintf("Round %d: lambda %.4f -> %.4f\n", round, lambda, lambda_next))

    log_w_inc <- numeric(N)
    for (n in seq_len(N)) {
      for (i in seq_along(modules[[n]])) {
        up <- advance(
          modules[[n]][[i]],
          phi = phi[n, , drop = TRUE],
          lambda_new = lambda_next,
          rho_local = rho_local,
          n_local_moves = n_local_moves,
          gaussian_map_fn = gaussian_map_fn,
          seed = as.integer(seed + 100000L * round + 1009L * n + i)
        )
        modules[[n]][[i]] <- up$module
        log_w_inc[n] <- log_w_inc[n] + up$log_u
        local_loglik_evals <- .add_local_loglik_evals(local_loglik_evals, up$local_loglik_evals)
        n_local_accept <- n_local_accept + up$n_accept
      }
    }

    norm <- .normalize_logweights(log(pmax(w, .Machine$double.eps)) + log_w_inc)
    w <- norm$w
    lambda <- lambda_next
    lambda_hist <- c(lambda_hist, lambda)
    ess_now <- .ess_from_weights(w)
    ess_hist <- c(ess_hist, ess_now)

    if (ess_now < rho_res * N) {
      rs <- resample_population_modules(phi, modules, w, logprior)
      phi <- rs$phi
      modules <- rs$modules
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
          modules_row = modules[[n]],
          logprior_row = logprior[n],
          lambda = lambda,
          logprior_phi = logprior_phi,
          gaussian_map_fn = gaussian_map_fn,
          rw_cov = rw_cov,
          rho_local = rho_local,
          n_population_refresh_moves = n_population_refresh_moves,
          seed = as.integer(seed + 500000L * round + 1009L * mv + n)
        )
        phi[n, ] <- out$phi
        modules[[n]] <- out$modules
        logprior[n] <- out$logprior
        local_loglik_evals <- .add_local_loglik_evals(local_loglik_evals, out$local_loglik_evals)
        n_local_accept <- n_local_accept + out$n_local_accept
        n_population_accept <- n_population_accept + as.integer(out$accepted)
      }
    }

    if (lambda >= 1 - 1e-12) break
  }

  local_loglik_evals <- .finalize_local_loglik_evals(local_loglik_evals)

  list(
    phi = phi,
    w = w,
    modules = modules,
    meta = list(
      rounds = round,
      lambda_hist = lambda_hist,
      ess_hist = ess_hist,
      local_loglik_evals = local_loglik_evals,
      n_exact_loglik = as.integer(local_loglik_evals["total"]),
      n_population_accept = as.integer(n_population_accept),
      n_local_accept = as.integer(n_local_accept)
    )
  )
}
