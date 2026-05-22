#!/usr/bin/env Rscript
# ============================================================================
# Theta support Pathfinder
# - Optimizer/Hessian support proposals for cheap population-level targets
# - Used after conservative rho-SMC support discovery, before local SMC anchors
# - Candidate contours are Gaussian probability ellipsoids; novelty is measured
#   against the empirical spacing of the existing rho-SMC support cloud.
# ============================================================================

if (!exists("%||%", mode = "function") ||
    !exists("weighted_cov", mode = "function") ||
    !exists("regularize_cov", mode = "function")) {
  source("smc_core.R")
}

if (!exists("normalize_population_model", mode = "function") ||
    !exists("population_model_log_hyperprior", mode = "function")) {
  source("population_models.R")
}

if (!exists("population_factor_set_logposterior", mode = "function")) {
  source("outer_population_smc.R")
}

if (!exists("local_likelihood_sketch_set_loglik", mode = "function") &&
    file.exists("local_likelihood_sketches.R")) {
  source("local_likelihood_sketches.R")
}

suppressPackageStartupMessages({
  library(parallel)
})

.theta_pf_as_matrix <- function(theta, model) {
  .as_hyper_matrix(theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
}

.theta_pf_normalize_weights <- function(w, n) {
  if (is.null(w)) return(rep(1 / n, n))
  w <- pmax(as.numeric(w), 0)
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) rep(1 / n, n) else w / sw
}

.theta_pf_weighted_quantile <- function(x, w, probs) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) {
    return(as.numeric(stats::quantile(as.numeric(x), probs = probs, na.rm = TRUE, names = FALSE)))
  }
  x <- x[ok]
  w <- w[ok]
  ord <- order(x)
  x <- as.numeric(x[ord])
  w <- .theta_pf_normalize_weights(w[ord], length(w))
  cw <- c(0, cumsum(w))
  qx <- c(x[1L], x)
  keep <- !duplicated(cw)
  as.numeric(stats::approx(cw[keep], qx[keep], xout = probs, rule = 2)$y)
}

theta_pathfinder_whitener <- function(theta, w = NULL, model = NULL) {
  if (!is.null(model)) {
    model <- normalize_population_model(model)
    theta <- .theta_pf_as_matrix(theta, model)
  } else {
    theta <- as.matrix(theta)
  }
  n <- nrow(theta)
  d <- ncol(theta)
  w <- .theta_pf_normalize_weights(w, n)
  center <- colSums(theta * w)
  S <- if (n > 1L) {
    tryCatch(weighted_cov(theta, w), error = function(e) stats::cov(theta))
  } else {
    diag(d)
  }
  if (is.null(S) || any(!is.finite(S))) S <- diag(d)
  S <- regularize_cov(S, min_eig = 1e-8, cond_cap = 1e8)
  L <- chol(S)
  names(center) <- colnames(theta)
  dimnames(S) <- list(colnames(theta), colnames(theta))
  list(center = center, cov = S, chol = L, names = colnames(theta))
}

theta_pathfinder_to_z <- function(theta, whitener) {
  theta <- as.matrix(theta)
  if (!is.null(whitener$names)) {
    if (!is.null(colnames(theta)) && setequal(colnames(theta), whitener$names)) {
      theta <- theta[, whitener$names, drop = FALSE]
    } else {
      colnames(theta) <- whitener$names
    }
  }
  centered <- sweep(theta, 2L, whitener$center, "-")
  z <- t(backsolve(whitener$chol, t(centered), transpose = TRUE))
  colnames(z) <- whitener$names
  z
}

theta_pathfinder_from_z <- function(z, whitener) {
  z <- as.matrix(z)
  theta <- sweep(z %*% whitener$chol, 2L, whitener$center, "+")
  colnames(theta) <- whitener$names
  theta
}

theta_pathfinder_select_starts <- function(theta,
                                           w = NULL,
                                           n_starts = 24L,
                                           n_pc = 4L,
                                           quantile_probs = c(0.01, 0.05, 0.50, 0.95, 0.99),
                                           seed = NULL,
                                           model = NULL) {
  if (!is.null(model)) {
    model <- normalize_population_model(model)
    theta <- .theta_pf_as_matrix(theta, model)
  } else {
    theta <- as.matrix(theta)
  }
  n <- nrow(theta)
  d <- ncol(theta)
  if (n == 0L) stop("theta must contain at least one row.")
  w <- .theta_pf_normalize_weights(w, n)
  n_starts <- as.integer(max(1L, n_starts))
  if (!is.null(seed)) set.seed(as.integer(seed))

  whitener <- theta_pathfinder_whitener(theta, w = w)
  z <- theta_pathfinder_to_z(theta, whitener)
  center <- matrix(whitener$center, nrow = 1L, dimnames = list(NULL, colnames(theta)))
  rows <- list(center)

  S <- tryCatch(weighted_cov(theta, w), error = function(e) NULL)
  if (!is.null(S) && all(is.finite(S)) && d > 0L) {
    eig <- eigen(regularize_cov(S, min_eig = 1e-8, cond_cap = 1e8), symmetric = TRUE)
    n_pc <- as.integer(min(max(1L, n_pc), d))
    for (j in seq_len(n_pc)) {
      score <- as.numeric(scale(theta, center = whitener$center, scale = FALSE) %*% eig$vectors[, j])
      qs <- .theta_pf_weighted_quantile(score, w, quantile_probs)
      for (q in qs) {
        idx <- which.min(abs(score - q))
        rows[[length(rows) + 1L]] <- theta[idx, , drop = FALSE]
      }
    }
  }

  if (length(rows) < n_starts) {
    idx <- sample.int(n, size = n_starts - length(rows), replace = n < (n_starts - length(rows)), prob = w)
    rows <- c(rows, lapply(idx, function(i) theta[i, , drop = FALSE]))
  }

  starts <- do.call(rbind, rows)
  z_key <- round(theta_pathfinder_to_z(starts, whitener), 4)
  keep <- !duplicated(apply(z_key, 1L, paste, collapse = "\r"))
  starts <- starts[keep, , drop = FALSE]
  if (nrow(starts) > n_starts) starts <- starts[seq_len(n_starts), , drop = FALSE]
  colnames(starts) <- colnames(theta)
  starts
}

.theta_pf_logposterior <- function(factor_set, theta) {
  as.numeric(population_factor_set_logposterior(
    factor_set = factor_set,
    theta = theta,
    include_constant = TRUE,
    n_cores = 1L
  ))
}

.theta_pf_fd_grad <- function(fn, z, eps = 1e-4) {
  z <- as.numeric(z)
  g <- numeric(length(z))
  for (j in seq_along(z)) {
    step <- rep(0, length(z))
    step[j] <- eps * max(1, abs(z[j]))
    fp <- fn(z + step)
    fm <- fn(z - step)
    g[j] <- (fp - fm) / (2 * step[j])
  }
  g
}

.theta_pf_cov_from_hessian <- function(hessian_z,
                                       theta_whitener,
                                       min_eig = 1e-6,
                                       cond_cap = 1e8) {
  d <- ncol(hessian_z)
  H <- as.matrix(hessian_z)
  valid <- nrow(H) == d && ncol(H) == d && all(is.finite(H))
  if (valid) {
    H <- (H + t(H)) / 2
    ev <- eigen(H, symmetric = TRUE)
    valid <- all(is.finite(ev$values)) && min(ev$values) > min_eig
  }
  if (valid) {
    H <- regularize_cov(H, min_eig = min_eig, cond_cap = cond_cap)
    cov_z <- tryCatch(solve(H), error = function(e) NULL)
    valid <- !is.null(cov_z) && all(is.finite(cov_z))
  }
  if (!valid) {
    cov_z <- diag(d)
  }
  cov_z <- regularize_cov(cov_z, min_eig = min_eig, cond_cap = cond_cap)
  L <- theta_whitener$chol
  cov_theta <- t(L) %*% cov_z %*% L
  cov_theta <- regularize_cov(cov_theta, min_eig = min_eig, cond_cap = cond_cap)
  dimnames(cov_theta) <- list(theta_whitener$names, theta_whitener$names)
  list(cov_z = cov_z, cov_theta = cov_theta, hessian_valid = valid)
}

.theta_pf_fit_one <- function(z_start,
                              factor_set,
                              whitener,
                              maxit = 500L) {
  model <- factor_set$population_model
  objective <- function(z) {
    theta <- theta_pathfinder_from_z(matrix(z, nrow = 1L), whitener)
    lp <- .theta_pf_logposterior(factor_set, theta)
    if (!is.finite(lp)) return(1e100)
    -lp
  }

  fit <- tryCatch(
    stats::optim(
      par = as.numeric(z_start),
      fn = objective,
      method = "BFGS",
      control = list(maxit = as.integer(maxit))
    ),
    error = function(e) list(par = as.numeric(z_start), value = Inf, convergence = 999L, message = conditionMessage(e))
  )
  hessian <- tryCatch(stats::optimHess(fit$par, objective), error = function(e) matrix(NA_real_, length(fit$par), length(fit$par)))
  grad <- tryCatch(.theta_pf_fd_grad(objective, fit$par), error = function(e) rep(NA_real_, length(fit$par)))
  theta <- theta_pathfinder_from_z(matrix(fit$par, nrow = 1L), whitener)
  theta <- .theta_pf_as_matrix(theta, model)
  logpost <- .theta_pf_logposterior(factor_set, theta)
  if (is.finite(fit$value)) logpost <- -as.numeric(fit$value)
  cov <- .theta_pf_cov_from_hessian(hessian, whitener)

  list(
    theta = theta,
    z = matrix(fit$par, nrow = 1L, dimnames = list(NULL, whitener$names)),
    logposterior = as.numeric(logpost),
    objective = as.numeric(fit$value),
    convergence = as.integer(fit$convergence %||% 999L),
    message = as.character(fit$message %||% ""),
    gradient_norm = sqrt(sum(grad^2)),
    hessian_valid = isTRUE(cov$hessian_valid),
    cov = cov$cov_theta,
    cov_z = cov$cov_z
  )
}

.theta_pf_unique_paths <- function(paths, distance_threshold = 0.25) {
  finite <- vapply(paths, function(x) is.finite(x$logposterior), logical(1))
  idx <- which(finite)
  if (!length(idx)) return(integer(0))
  idx <- idx[order(vapply(paths[idx], `[[`, numeric(1), "logposterior"), decreasing = TRUE)]
  keep <- integer(0)
  kept_z <- NULL
  for (i in idx) {
    zi <- as.numeric(paths[[i]]$z[1L, ])
    if (is.null(kept_z)) {
      keep <- c(keep, i)
      kept_z <- matrix(zi, nrow = 1L)
      next
    }
    dist <- sqrt(rowSums((sweep(kept_z, 2L, zi, "-"))^2))
    if (min(dist) > as.numeric(distance_threshold)) {
      keep <- c(keep, i)
      kept_z <- rbind(kept_z, zi)
    }
  }
  keep
}

.theta_pf_usable_paths <- function(paths,
                                   unique_idx) {
  if (!length(unique_idx)) return(integer(0))
  unique_idx[vapply(unique_idx, function(i) {
    path <- paths[[i]]
    is.finite(path$logposterior) &&
      identical(path$convergence, 0L) &&
      isTRUE(path$hessian_valid)
  }, logical(1))]
}

theta_pathfinder_component_sigma_points <- function(mean,
                                                    cov,
                                                    coverage_prob = 0.95,
                                                    max_directions = NULL) {
  mean <- as.numeric(mean)
  cov <- regularize_cov(cov, min_eig = 1e-8, cond_cap = 1e8)
  eig <- eigen(cov, symmetric = TRUE)
  d <- length(mean)
  radius <- sqrt(stats::qchisq(as.numeric(coverage_prob), df = d))
  dirs <- seq_len(d)
  if (!is.null(max_directions)) dirs <- head(dirs, as.integer(max_directions))
  rows <- list(mean)
  labels <- "center"
  for (j in dirs) {
    step <- as.numeric(radius) * sqrt(max(eig$values[j], 0)) * eig$vectors[, j]
    rows[[length(rows) + 1L]] <- mean + step
    labels <- c(labels, paste0("pc", j, "_plus"))
    rows[[length(rows) + 1L]] <- mean - step
    labels <- c(labels, paste0("pc", j, "_minus"))
  }
  out <- do.call(rbind, rows)
  list(theta = out, label = labels)
}

fit_theta_pathfinder_support <- function(factor_set,
                                         support_theta,
                                         support_w = NULL,
                                         starts = NULL,
                                         n_starts = 24L,
                                         coverage_probs = c(0.50, 0.90, 0.99),
                                         maxit = 500L,
                                         unique_distance = 0.25,
                                         max_sigma_directions = NULL,
                                         n_jobs = 1L,
                                         seed = 123L,
                                         verbose = TRUE) {
  stopifnot(inherits(factor_set, "population_factor_set"))
  model <- normalize_population_model(factor_set$population_model)
  support_theta <- .theta_pf_as_matrix(support_theta, model)
  support_w <- .theta_pf_normalize_weights(support_w, nrow(support_theta))

  if (is.null(starts)) {
    starts <- theta_pathfinder_select_starts(
      theta = support_theta,
      w = support_w,
      n_starts = n_starts,
      seed = seed,
      model = model
    )
  } else {
    starts <- .theta_pf_as_matrix(starts, model)
    if (nrow(starts) > as.integer(n_starts)) {
      starts <- theta_pathfinder_select_starts(
        theta = starts,
        n_starts = n_starts,
        seed = seed,
        model = model
      )
    }
  }

  whitener <- theta_pathfinder_whitener(support_theta, w = support_w, model = model)
  z_starts <- theta_pathfinder_to_z(starts, whitener)
  ids <- seq_len(nrow(z_starts))
  if (isTRUE(verbose)) {
    cat(sprintf("Theta Pathfinder: %d starts | dim=%d\n", length(ids), model$hyper_dim))
  }
  paths <- parallel::mclapply(
    ids,
    function(i) {
      .theta_pf_fit_one(
        z_start = z_starts[i, ],
        factor_set = factor_set,
        whitener = whitener,
        maxit = maxit
      )
    },
    mc.cores = as.integer(max(1L, n_jobs))
  )
  unique_idx <- .theta_pf_unique_paths(paths, distance_threshold = unique_distance)
  usable_idx <- .theta_pf_usable_paths(paths, unique_idx = unique_idx)

  components <- list()
  candidates <- list()
  comp_id <- 0L
  cand_id <- 0L
  for (path_id in usable_idx) {
    path <- paths[[path_id]]
    for (prob in as.numeric(coverage_probs)) {
      comp_id <- comp_id + 1L
      radius <- sqrt(stats::qchisq(prob, df = model$hyper_dim))
      components[[comp_id]] <- list(
        component_id = comp_id,
        path_id = path_id,
        coverage_prob = prob,
        radius = radius,
        mean = path$theta[1L, ],
        cov = path$cov,
        logposterior = path$logposterior,
        hessian_valid = path$hessian_valid
      )
      sigma <- theta_pathfinder_component_sigma_points(
        mean = path$theta[1L, ],
        cov = path$cov,
        coverage_prob = prob,
        max_directions = max_sigma_directions
      )
      colnames(sigma$theta) <- model$hyper_names
      for (j in seq_len(nrow(sigma$theta))) {
        cand_id <- cand_id + 1L
        candidates[[cand_id]] <- data.frame(
          candidate_id = cand_id,
          component_id = comp_id,
          path_id = path_id,
          coverage_prob = prob,
          radius = radius,
          point = sigma$label[j],
          matrix(sigma$theta[j, ], nrow = 1L, dimnames = list(NULL, model$hyper_names)),
          check.names = FALSE
        )
      }
    }
  }
  candidate_df <- if (length(candidates)) do.call(rbind, candidates) else {
    data.frame(candidate_id = integer(0), component_id = integer(0), path_id = integer(0), coverage_prob = numeric(0), radius = numeric(0), point = character(0))
  }

  path_summary <- do.call(rbind, lapply(seq_along(paths), function(i) {
    data.frame(
      path_id = i,
      logposterior = paths[[i]]$logposterior,
      convergence = paths[[i]]$convergence,
      gradient_norm = paths[[i]]$gradient_norm,
      hessian_valid = paths[[i]]$hessian_valid,
      unique = i %in% unique_idx,
      usable = i %in% usable_idx,
      check.names = FALSE
    )
  }))

  structure(
    list(
      paths = paths,
      path_summary = path_summary,
      unique_path_ids = unique_idx,
      usable_path_ids = usable_idx,
      components = components,
      candidates = candidate_df,
      starts = starts,
      support_theta = support_theta,
      support_w = support_w,
      whitener = whitener,
      factor_set = factor_set,
      settings = list(
        n_starts = as.integer(nrow(starts)),
        coverage_probs = as.numeric(coverage_probs),
        maxit = as.integer(maxit),
        unique_distance = as.numeric(unique_distance),
        seed = as.integer(seed)
      )
    ),
    class = "theta_pathfinder_support"
  )
}

theta_pathfinder_candidate_distance <- function(pathfinder,
                                                reference_theta = pathfinder$support_theta,
                                                reference_w = pathfinder$support_w,
                                                distance_quantile = 0.95,
                                                spacing_sample_size = 2000L) {
  if (!inherits(pathfinder, "theta_pathfinder_support")) {
    stop("pathfinder must inherit from 'theta_pathfinder_support'.")
  }
  model <- normalize_population_model(pathfinder$factor_set$population_model)
  candidates <- pathfinder$candidates
  if (!nrow(candidates)) return(candidates)
  cand_theta <- .theta_pf_as_matrix(candidates[, model$hyper_names, drop = FALSE], model)
  reference_theta <- .theta_pf_as_matrix(reference_theta, model)
  whitener <- theta_pathfinder_whitener(reference_theta, w = reference_w, model = model)
  cand_z <- theta_pathfinder_to_z(cand_theta, whitener)
  ref_z <- theta_pathfinder_to_z(reference_theta, whitener)
  ref_z <- ref_z[!duplicated(as.data.frame(ref_z)), , drop = FALSE]

  nearest <- numeric(nrow(cand_z))
  for (i in seq_len(nrow(cand_z))) {
    dz <- sweep(ref_z, 2L, cand_z[i, ], "-")
    nearest[i] <- sqrt(min(rowSums(dz^2)))
  }
  z_spacing <- ref_z
  if (nrow(z_spacing) > as.integer(spacing_sample_size)) {
    idx <- unique(as.integer(round(seq(1, nrow(z_spacing), length.out = as.integer(spacing_sample_size)))))
    z_spacing <- z_spacing[idx, , drop = FALSE]
  }
  spacing <- numeric(nrow(z_spacing))
  for (i in seq_len(nrow(z_spacing))) {
    dz <- sweep(ref_z, 2L, z_spacing[i, ], "-")
    d2 <- rowSums(dz^2)
    d2[d2 == 0] <- Inf
    spacing[i] <- sqrt(min(d2))
  }
  spacing <- spacing[is.finite(spacing)]
  distance_threshold <- if (length(spacing)) {
    as.numeric(stats::quantile(spacing, probs = as.numeric(distance_quantile), names = FALSE, na.rm = TRUE))
  } else {
    Inf
  }
  candidates$nearest_support_distance <- nearest
  candidates$support_distance_cutoff <- distance_threshold
  candidates$nonduplicate_support <- nearest > distance_threshold
  candidates
}

theta_pathfinder_score_candidates <- function(pathfinder,
                                              candidates = pathfinder$candidates,
                                              n_cores = 1L) {
  if (!inherits(pathfinder, "theta_pathfinder_support")) {
    stop("pathfinder must inherit from 'theta_pathfinder_support'.")
  }
  model <- normalize_population_model(pathfinder$factor_set$population_model)
  if (!nrow(candidates)) return(candidates)
  theta <- .theta_pf_as_matrix(candidates[, model$hyper_names, drop = FALSE], model)
  logposterior <- as.numeric(population_factor_set_logposterior(
    factor_set = pathfinder$factor_set,
    theta = theta,
    include_constant = TRUE,
    n_cores = n_cores
  ))
  candidates$sketch_logposterior <- logposterior
  finite <- is.finite(logposterior)
  candidates$relative_logposterior <- NA_real_
  if (any(finite)) {
    candidates$relative_logposterior[finite] <- logposterior[finite] - max(logposterior[finite])
  }
  candidates
}
