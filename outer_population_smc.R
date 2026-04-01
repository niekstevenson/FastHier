#!/usr/bin/env Rscript
# ============================================================================
# Outer population SMC
# - Reusable local marginal-likelihood factors built from reference local fits
# - Tempered outer SMC for hyperposterior and hyper-evidence estimation
# ============================================================================

if (!exists("%||%", mode = "function") ||
    !exists("weighted_cov", mode = "function") ||
    !exists(".rowLogSumExp", mode = "function") ||
    !exists("next_lambda_via_rCESS", mode = "function") ||
    !exists("cess_target_at_lambda", mode = "function") ||
    !exists("stratified_resample_sorted", mode = "function")) {
  source("smc_core.R")
}
if (!exists("normalize_reference_prior", mode = "function") ||
    !exists("reference_prior_logpdf", mode = "function")) {
  source("reference_priors.R")
}
if (!exists("normalize_population_model", mode = "function") ||
    !exists("population_model_prepare_theta", mode = "function") ||
    !exists("population_model_log_alpha_given_prepared_theta_many", mode = "function")) {
  source("population_models.R")
}

suppressPackageStartupMessages({
  library(parallel)
})

normalize_reference_local_weights <- function(w) {
  w <- pmax(as.numeric(w), 0)
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) {
    rep(1 / length(w), length(w))
  } else {
    w / sw
  }
}

validate_reference_local_component <- function(component) {
  if (!is.list(component)) stop("reference local component must be a list.")
  required <- c("particles", "weights", "reference_prior")
  missing <- setdiff(required, names(component))
  if (length(missing)) {
    stop("reference local component is missing: ", paste(missing, collapse = ", "))
  }
  particles <- as.matrix(component$particles)
  weights <- normalize_reference_local_weights(component$weights)
  reference_prior <- normalize_reference_prior(reference_prior = component$reference_prior)
  log_reference_density <- as.numeric(
    component$log_reference_density %||% reference_prior_logpdf(reference_prior, particles)
  )
  component$particles <- particles
  component$weights <- weights
  component$reference_prior <- reference_prior
  component$log_reference_density <- log_reference_density
  component$log_evidence <- as.numeric(component$log_evidence %||% NA_real_)
  component$mcse_log_evidence <- as.numeric(component$mcse_log_evidence %||% NA_real_)
  component
}

validate_reference_local_object <- function(local_object) {
  if (!is.list(local_object)) stop("local_object must be a list.")
  local_object$local_id <- as.integer(local_object$local_id %||% NA_integer_)

  if (!is.null(local_object$components)) {
    if (!length(local_object$components)) {
      stop("local_object components must be non-empty.")
    }
    local_object$components <- lapply(local_object$components, validate_reference_local_component)
    mixture_weights <- local_object$mixture_weights %||% rep(1 / length(local_object$components), length(local_object$components))
    mixture_weights <- normalize_reference_local_weights(mixture_weights)
    if (length(mixture_weights) != length(local_object$components)) {
      stop("local_object mixture_weights must match the number of components.")
    }
    local_object$mixture_weights <- mixture_weights
    return(local_object)
  }

  required <- c("particles", "weights", "reference_prior")
  missing <- setdiff(required, names(local_object))
  if (length(missing)) {
    stop("local_object is missing: ", paste(missing, collapse = ", "))
  }
  local_object$components <- list(
    validate_reference_local_component(
      list(
        particles = local_object$particles,
        weights = local_object$weights,
        log_evidence = local_object$log_evidence,
        mcse_log_evidence = local_object$mcse_log_evidence,
        reference_prior = local_object$reference_prior,
        log_reference_density = local_object$log_reference_density,
        diagnostics = local_object$diagnostics %||% list()
      )
    )
  )
  local_object$mixture_weights <- 1
  local_object
}

build_population_local_factor <- function(local_object, population_model) {
  local_object <- validate_reference_local_object(local_object)
  population_model <- normalize_population_model(population_model)
  components <- local_object$components

  align_alpha <- function(alpha) {
    alpha <- as.matrix(alpha)
    if (ncol(alpha) != population_model$alpha_dim) {
      stop("Local particle dimension does not match the population model.")
    }
    if (!is.null(colnames(alpha)) && setequal(colnames(alpha), population_model$alpha_names)) {
      alpha[, population_model$alpha_names, drop = FALSE]
    } else {
      colnames(alpha) <- population_model$alpha_names
      alpha
    }
  }

  if (length(components) == 1L) {
    component <- components[[1L]]
    alpha <- align_alpha(component$particles)
    factor_log_weights <- log(component$weights)
    factor_log_reference_density <- as.numeric(component$log_reference_density %||% reference_prior_logpdf(component$reference_prior, alpha))
    log_base <- factor_log_weights - factor_log_reference_density
    log_constant <- as.numeric(component$log_evidence %||% 0)
    component_id <- rep.int(1L, nrow(alpha))
  } else {
    eta <- normalize_reference_local_weights(local_object$mixture_weights)
    log_eta <- ifelse(eta > 0, log(eta), -Inf)
    alpha_list <- lapply(components, function(component) align_alpha(component$particles))
    alpha <- do.call(rbind, alpha_list)
    component_lengths <- vapply(alpha_list, nrow, integer(1))
    component_id <- rep.int(seq_along(alpha_list), component_lengths)

    component_logw <- unlist(
      Map(function(component, eta_r) {
        if (!is.finite(eta_r)) {
          rep(-Inf, nrow(component$particles))
        } else {
          log(component$weights) + eta_r
        }
      }, components, as.list(log_eta)),
      use.names = FALSE
    )

    component_logZ <- vapply(components, function(component) as.numeric(component$log_evidence %||% NA_real_), numeric(1))
    if (any(!is.finite(component_logZ))) {
      stop("DMIS local factors require finite log_evidence for every component.")
    }

    logq_mat <- vapply(
      components,
      function(component) as.numeric(reference_prior_logpdf(component$reference_prior, alpha)),
      numeric(nrow(alpha))
    )
    if (!is.matrix(logq_mat)) {
      logq_mat <- matrix(logq_mat, ncol = length(components))
    }
      log_d <- .rowLogSumExp(sweep(logq_mat, 2L, log_eta - component_logZ, "+"))
    log_base <- component_logw - log_d
    log_constant <- 0
    factor_log_weights <- component_logw
    factor_log_reference_density <- rep(NA_real_, nrow(alpha))
  }

  structure(
    list(
      local_id = as.integer(local_object$local_id),
      population_model = population_model,
      particles = alpha,
      log_weights = factor_log_weights,
      log_reference_density = factor_log_reference_density,
      log_base = log_base,
      log_constant = log_constant,
      reference_prior = if (length(components) == 1L) components[[1L]]$reference_prior else NULL,
      component_id = component_id,
      fast_family = as.character(population_model$fast_family %||% ""),
      fast_cache = NULL
    ),
    class = "population_local_factor"
  )
}

build_population_local_site <- function(local_object, population_model) {
  base_factor <- build_population_local_factor(
    local_object = local_object,
    population_model = population_model
  )
  structure(
    list(
      local_id = as.integer(base_factor$local_id),
      population_model = base_factor$population_model,
      base_factor = base_factor,
      corrections = list(),
      diagnostics = list()
    ),
    class = "population_local_site"
  )
}

population_local_site_log_marginal_many <- function(site, theta, include_constant = TRUE, block_size = 1024L) {
  stopifnot(inherits(site, "population_local_site"))
  out <- population_local_factor_log_marginal_many(
    factor = site$base_factor,
    theta = theta,
    include_constant = include_constant,
    block_size = block_size
  )
  if (!length(site$corrections)) {
    return(out)
  }

  theta_mat <- .as_hyper_matrix(
    theta,
    hyper_names = site$population_model$hyper_names,
    hyper_dim = site$population_model$hyper_dim
  )
  correction_total <- rep.int(0, nrow(theta_mat))
  for (correction in site$corrections) {
    if (is.null(correction$surrogate) || is.null(correction$hyper_names)) {
      next
    }
    block_values <- theta_mat[, correction$hyper_names, drop = FALSE]
    correction_total <- correction_total + predict_site_block_surrogate(correction$surrogate, block_values)
  }
  as.numeric(out + correction_total)
}

population_local_site_log_marginal <- function(site, theta, include_constant = TRUE) {
  population_local_site_log_marginal_many(
    site = site,
    theta = theta,
    include_constant = include_constant
  )[1L]
}

.fit_ridge_regression <- function(X, y, ridge = 1e-4, penalize_intercept = FALSE) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  pen <- diag(ncol(X))
  if (!isTRUE(penalize_intercept) && ncol(X) >= 1L) {
    pen[1L, 1L] <- 0
  }
  XtX <- crossprod(X) + as.numeric(ridge) * pen
  Xty <- crossprod(X, y)
  tryCatch(
    as.numeric(solve(XtX, Xty)),
    error = function(e) as.numeric(qr.solve(XtX, Xty))
  )
}

.make_rbf_features <- function(X, centers, scale) {
  X <- as.matrix(X)
  centers <- as.matrix(centers)
  if (!nrow(X) || !nrow(centers)) {
    return(matrix(0, nrow = nrow(X), ncol = 0L))
  }
  d2 <- outer(
    seq_len(nrow(X)),
    seq_len(nrow(centers)),
    Vectorize(function(i, j) sum((X[i, ] - centers[j, ])^2))
  )
  exp(-d2 / (2 * scale^2))
}

predict_site_block_surrogate <- function(surrogate, block_values) {
  if (is.null(surrogate)) {
    return(rep.int(0, nrow(as.matrix(block_values))))
  }
  block_values <- as.matrix(block_values)
  centers <- as.numeric(surrogate$centers %||% rep(0, ncol(block_values)))
  scales <- as.numeric(surrogate$scales %||% rep(1, ncol(block_values)))
  scales[!is.finite(scales) | scales < 1e-8] <- 1
  Xs <- sweep(sweep(block_values, 2L, centers, "-"), 2L, scales, "/")
  colnames(Xs) <- surrogate$hyper_names %||% colnames(block_values)
  design <- .build_quadratic_design(Xs)
  if (!is.null(surrogate$rbf_centers_scaled) &&
      nrow(as.matrix(surrogate$rbf_centers_scaled)) > 0L &&
      is.finite(surrogate$rbf_scale) &&
      surrogate$rbf_scale > 0) {
    design <- cbind(
      design,
      .make_rbf_features(Xs, surrogate$rbf_centers_scaled, surrogate$rbf_scale)
    )
  }
  pred <- as.numeric(design %*% surrogate$coef)
  clip <- as.numeric(surrogate$clip %||% c(-Inf, Inf))
  pred <- pmin(pmax(pred, clip[1L]), clip[2L])
  pred
}

fit_site_block_surrogate <- function(block_values,
                                     delta,
                                     family = c("quadratic_rbf", "quadratic"),
                                     ridge = 1e-4,
                                     rbf_scale = NULL,
                                     clip_expand = 0.5,
                                     compute_cv = TRUE) {
  family <- match.arg(family)
  block_values <- as.matrix(block_values)
  delta <- as.numeric(delta)
  if (nrow(block_values) != length(delta)) {
    stop("block_values and delta must have matching lengths.")
  }
  if (!nrow(block_values)) {
    stop("At least one anchor point is required.")
  }
  centers <- colMeans(block_values)
  scales <- apply(block_values, 2L, stats::sd)
  scales[!is.finite(scales) | scales < 1e-8] <- 1
  Xs <- sweep(sweep(block_values, 2L, centers, "-"), 2L, scales, "/")
  colnames(Xs) <- colnames(block_values)

  design <- .build_quadratic_design(Xs)
  use_rbf <- identical(family, "quadratic_rbf") && nrow(block_values) >= max(6L, ncol(block_values) + 3L)
  rbf_centers_scaled <- NULL
  rbf_scale_eff <- NA_real_
  if (use_rbf) {
    if (is.null(rbf_scale)) {
      dmat <- as.matrix(dist(Xs))
      dvec <- dmat[upper.tri(dmat)]
      dvec <- dvec[is.finite(dvec) & dvec > 0]
      rbf_scale_eff <- if (length(dvec)) stats::median(dvec) else 1
    } else {
      rbf_scale_eff <- as.numeric(rbf_scale)
    }
    rbf_scale_eff <- max(rbf_scale_eff, 0.25)
    rbf_centers_scaled <- Xs
    design <- cbind(
      design,
      .make_rbf_features(Xs, rbf_centers_scaled, rbf_scale_eff)
    )
  }

  coef <- .fit_ridge_regression(design, delta, ridge = ridge, penalize_intercept = FALSE)
  fitted <- as.numeric(design %*% coef)
  delta_range <- diff(range(delta))
  clip_pad <- max(0.1, as.numeric(clip_expand)) * max(delta_range, 0.05)
  clip <- c(min(delta) - clip_pad, max(delta) + clip_pad)

  cv_rmse <- NA_real_
  if (isTRUE(compute_cv) && nrow(block_values) >= 4L) {
    preds <- numeric(nrow(block_values))
    for (i in seq_len(nrow(block_values))) {
      fit_i <- fit_site_block_surrogate(
        block_values = block_values[-i, , drop = FALSE],
        delta = delta[-i],
        family = if (use_rbf) "quadratic_rbf" else "quadratic",
        ridge = ridge,
        rbf_scale = if (use_rbf) rbf_scale_eff else NULL,
        clip_expand = clip_expand,
        compute_cv = FALSE
      )
      preds[i] <- predict_site_block_surrogate(fit_i, block_values[i, , drop = FALSE])[1L]
    }
    cv_rmse <- sqrt(mean((delta - preds)^2))
  }

  structure(
    list(
      family = if (use_rbf) "quadratic_rbf" else "quadratic",
      hyper_names = colnames(block_values),
      centers = centers,
      scales = scales,
      coef = coef,
      rbf_centers_scaled = rbf_centers_scaled,
      rbf_scale = rbf_scale_eff,
      clip = clip,
      diagnostics = list(
        rmse = sqrt(mean((delta - fitted)^2)),
        cv_rmse = cv_rmse,
        n_anchors = nrow(block_values)
      )
    ),
    class = "site_block_surrogate"
  )
}

.pack_alpha_quadratic_terms_t <- function(alpha) {
  alpha <- as.matrix(alpha)
  d <- ncol(alpha)
  out <- matrix(0, nrow = d * (d + 1L) / 2L, ncol = nrow(alpha))
  pos <- 1L
  for (j in seq_len(d)) {
    aj <- alpha[, j]
    out[pos, ] <- aj * aj
    pos <- pos + 1L
    if (j < d) {
      for (k in (j + 1L):d) {
        out[pos, ] <- 2 * aj * alpha[, k]
        pos <- pos + 1L
      }
    }
  }
  out
}

.gaussian_log_alpha_given_prepared_theta_many <- function(alpha, theta_prepared) {
  alpha <- as.matrix(alpha)
  linear <- theta_prepared$eta %*% t(alpha)
  quad_t <- if (identical(theta_prepared$quadratic_kind, "diag")) {
    t(alpha * alpha)
  } else {
    .pack_alpha_quadratic_terms_t(alpha)
  }
  quad <- theta_prepared$quadratic_coef %*% quad_t
  sweep(linear - 0.5 * quad, 1L, theta_prepared$log_kernel_constant, "+")
}

.population_log_alpha_given_theta_many <- function(model,
                                                   alpha,
                                                   theta = NULL,
                                                   theta_prepared = NULL) {
  if (is.null(theta_prepared)) {
    theta_prepared <- population_model_prepare_theta(model, theta)
  }
  if (identical(theta_prepared$family, "gaussian")) {
    return(.gaussian_log_alpha_given_prepared_theta_many(alpha, theta_prepared))
  }
  population_model_log_alpha_given_prepared_theta_many(model, alpha = alpha, theta_prepared = theta_prepared)
}

.local_log_marginal_blocked <- function(model,
                                        alpha,
                                        log_base,
                                        theta,
                                        theta_prepared = NULL,
                                        block_size = 1024L) {
  theta <- .as_hyper_matrix(theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  theta_prepared <- theta_prepared %||% population_model_prepare_theta(model, theta)
  block_size <- as.integer(max(1L, block_size))
  accum <- rep.int(-Inf, nrow(theta))
  n_alpha <- nrow(alpha)

  for (start in seq.int(1L, n_alpha, by = block_size)) {
    idx <- seq.int(start, min(start + block_size - 1L, n_alpha))
    log_terms <- .population_log_alpha_given_theta_many(
      model = model,
      alpha = alpha[idx, , drop = FALSE],
      theta_prepared = theta_prepared
    )
    block_lse <- .rowLogSumExp(sweep(log_terms, 2L, log_base[idx], "+"))
    accum <- rlogsumexp2(accum, block_lse)
  }
  as.numeric(accum)
}

population_local_factor_log_marginal_many <- function(factor, theta, include_constant = TRUE, block_size = 1024L) {
  stopifnot(inherits(factor, "population_local_factor"))
  theta_prepared <- population_model_prepare_theta(factor$population_model, theta)
  out <- .local_log_marginal_blocked(
    model = factor$population_model,
    alpha = factor$particles,
    log_base = factor$log_base,
    theta = theta_prepared$theta,
    theta_prepared = theta_prepared,
    block_size = block_size
  )
  if (isTRUE(include_constant)) out <- out + factor$log_constant
  as.numeric(out)
}

population_local_factor_log_marginal <- function(factor, theta, include_constant = TRUE) {
  population_local_factor_log_marginal_many(factor, theta = theta, include_constant = include_constant)[1L]
}

population_local_factor_ess <- function(factor, theta) {
  stopifnot(inherits(factor, "population_local_factor"))
  theta_prepared <- population_model_prepare_theta(factor$population_model, theta)
  logp <- .population_log_alpha_given_theta_many(
    model = factor$population_model,
    alpha = factor$particles,
    theta_prepared = theta_prepared
  )[1L, ]
  lw <- factor$log_base + logp
  lse <- logsumexp(lw)
  if (!is.finite(lse)) return(0)
  w <- exp(lw - lse)
  as.numeric(1 / sum(w * w))
}

population_local_factor_reweighted_particles <- function(factor, theta) {
  stopifnot(inherits(factor, "population_local_factor"))
  theta_prepared <- population_model_prepare_theta(factor$population_model, theta)
  logp <- .population_log_alpha_given_theta_many(
    model = factor$population_model,
    alpha = factor$particles,
    theta_prepared = theta_prepared
  )[1L, ]
  lw <- factor$log_base + logp
  lse <- logsumexp(lw)
  if (!is.finite(lse)) {
    weights <- rep(1 / nrow(factor$particles), nrow(factor$particles))
  } else {
    weights <- exp(lw - lse)
  }
  list(
    particles = factor$particles,
    weights = weights
  )
}

population_local_factor_tail_diagnostic <- function(factor, theta, use_psis = TRUE) {
  stopifnot(inherits(factor, "population_local_factor"))
  theta_prepared <- population_model_prepare_theta(factor$population_model, theta)
  logp <- .population_log_alpha_given_theta_many(
    model = factor$population_model,
    alpha = factor$particles,
    theta_prepared = theta_prepared
  )[1L, ]
  lw <- factor$log_base + logp
  keep <- is.finite(lw)
  if (!any(keep)) {
    return(list(
      ess = 0,
      pareto_k = NA_real_,
      max_weight = NA_real_,
      q99_weight = NA_real_,
      n_finite = 0L
    ))
  }

  lw <- lw[keep]
  lse <- logsumexp(lw)
  if (!is.finite(lse)) {
    return(list(
      ess = 0,
      pareto_k = NA_real_,
      max_weight = NA_real_,
      q99_weight = NA_real_,
      n_finite = length(lw)
    ))
  }

  w <- exp(lw - lse)
  pareto_k <- NA_real_
  if (isTRUE(use_psis) &&
      length(lw) >= 5L &&
      requireNamespace("loo", quietly = TRUE)) {
    psis_obj <- tryCatch(
      suppressWarnings(loo::psis(matrix(lw - max(lw), ncol = 1L))),
      error = function(e) NULL
    )
    if (!is.null(psis_obj)) {
      pareto_k <- as.numeric(loo::pareto_k_values(psis_obj)[1L])
    }
  }

  list(
    ess = as.numeric(1 / sum(w * w)),
    pareto_k = pareto_k,
    max_weight = as.numeric(max(w)),
    q99_weight = as.numeric(stats::quantile(w, probs = 0.99, na.rm = TRUE, names = FALSE)),
    n_finite = length(lw)
  )
}

.weighted_mean_vec <- function(x, w) {
  w <- normalize_reference_local_weights(w)
  as.numeric(sum(as.numeric(x) * w))
}

.weighted_quantile_outer <- function(x, w, probs) {
  x <- as.numeric(x)
  w <- normalize_reference_local_weights(w)
  ord <- order(x)
  x_ord <- x[ord]
  w_ord <- w[ord]
  cw <- cumsum(w_ord)
  as.numeric(stats::approx(cw, x_ord, xout = probs, rule = 2)$y)
}

.weighted_marginal_geometry <- function(x, w) {
  mu <- .weighted_mean_vec(x, w)
  xc <- as.numeric(x) - mu
  var_x <- .weighted_mean_vec(xc * xc, w)
  sd_x <- sqrt(max(var_x, 1e-12))
  z <- xc / sd_x
  skew <- .weighted_mean_vec(z^3, w)
  kurt_excess <- .weighted_mean_vec(z^4, w) - 3
  qq <- .weighted_quantile_outer(x, w, probs = c(0.05, 0.25, 0.50, 0.75, 0.95))
  left_tail <- max(qq[3L] - qq[1L], 1e-8)
  right_tail <- max(qq[5L] - qq[3L], 1e-8)
  left_core <- max(qq[3L] - qq[2L], 1e-8)
  right_core <- max(qq[4L] - qq[3L], 1e-8)
  tail_asymmetry <- abs(log(right_tail / left_tail))
  core_asymmetry <- abs(log(right_core / left_core))
  score <- abs(skew) + 0.35 * abs(kurt_excess) + 0.50 * tail_asymmetry + 0.25 * core_asymmetry
  list(
    mean = mu,
    sd = sd_x,
    skew = skew,
    kurt_excess = kurt_excess,
    tail_asymmetry = tail_asymmetry,
    core_asymmetry = core_asymmetry,
    score = score
  )
}

.make_hyper_block_id <- function(hyper_names) {
  paste(as.character(hyper_names), collapse = "__")
}

.build_quadratic_design <- function(X) {
  X <- as.matrix(X)
  d <- ncol(X)
  cols <- list("(Intercept)" = rep.int(1, nrow(X)))
  for (j in seq_len(d)) {
    cols[[colnames(X)[j] %||% paste0("x", j)]] <- X[, j]
  }
  for (j in seq_len(d)) {
    nm <- paste0(colnames(X)[j] %||% paste0("x", j), "^2")
    cols[[nm]] <- X[, j]^2
  }
  if (d > 1L) {
    for (j in seq_len(d - 1L)) {
      for (k in (j + 1L):d) {
        nm <- paste0(colnames(X)[j] %||% paste0("x", j), ":", colnames(X)[k] %||% paste0("x", k))
        cols[[nm]] <- X[, j] * X[, k]
      }
    }
  }
  out <- as.matrix(as.data.frame(cols, check.names = FALSE))
  storage.mode(out) <- "double"
  out
}

.fit_quadratic_shape <- function(block_values, log_values) {
  block_values <- as.matrix(block_values)
  log_values <- as.numeric(log_values)
  if (!length(log_values) || !any(is.finite(log_values))) {
    return(list(scaled_rmse = NA_real_, rank = 0L))
  }
  centers <- colMeans(block_values)
  scales <- apply(block_values, 2L, stats::sd)
  scales[!is.finite(scales) | scales < 1e-8] <- 1
  Xs <- sweep(sweep(block_values, 2L, centers, "-"), 2L, scales, "/")
  colnames(Xs) <- colnames(block_values)
  design <- .build_quadratic_design(Xs)
  if (nrow(design) < ncol(design)) {
    return(list(scaled_rmse = NA_real_, rank = qr(design)$rank))
  }
  coef <- tryCatch(qr.solve(design, log_values), error = function(e) NULL)
  if (is.null(coef)) {
    return(list(scaled_rmse = NA_real_, rank = qr(design)$rank))
  }
  fitted <- as.numeric(design %*% coef)
  value_range <- diff(range(log_values))
  scaled_rmse <- sqrt(mean((log_values - fitted)^2)) / max(value_range, 1e-8)
  list(
    scaled_rmse = as.numeric(scaled_rmse),
    rank = qr(design)$rank
  )
}

.build_block_grid <- function(theta,
                              w,
                              theta_center,
                              block,
                              probs_1d = c(0.05, 0.20, 0.50, 0.80, 0.95),
                              probs_2d = c(0.10, 0.50, 0.90)) {
  theta <- as.matrix(theta)
  theta_center <- .as_hyper_matrix(
    theta_center,
    hyper_names = colnames(theta),
    hyper_dim = ncol(theta)
  )
  idx <- as.integer(block$indices)
  dim_block <- length(idx)
  if (dim_block < 1L || dim_block > 2L) {
    stop("Only 1D and 2D blocks are supported in phase 1.")
  }

  if (dim_block == 1L) {
    vals <- unique(.weighted_quantile_outer(theta[, idx], w, probs = probs_1d))
    block_values <- matrix(vals, ncol = 1L)
    colnames(block_values) <- block$hyper_names
    theta_grid <- matrix(rep(theta_center, each = nrow(block_values)), nrow = nrow(block_values), byrow = FALSE)
    colnames(theta_grid) <- colnames(theta)
    theta_grid[, idx] <- vals
    is_boundary <- seq_len(nrow(block_values)) %in% c(1L, nrow(block_values))
  } else {
    vals_list <- lapply(idx, function(j) unique(.weighted_quantile_outer(theta[, j], w, probs = probs_2d)))
    names(vals_list) <- block$hyper_names
    grid_df <- expand.grid(vals_list, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    block_values <- as.matrix(grid_df)
    colnames(block_values) <- block$hyper_names
    theta_grid <- matrix(rep(theta_center, each = nrow(block_values)), nrow = nrow(block_values), byrow = FALSE)
    colnames(theta_grid) <- colnames(theta)
    theta_grid[, idx] <- block_values
    is_boundary <- rep(FALSE, nrow(block_values))
    for (j in seq_len(ncol(block_values))) {
      is_boundary <- is_boundary |
        block_values[, j] %in% range(block_values[, j], finite = TRUE)
    }
  }

  list(
    theta = theta_grid,
    block_values = block_values,
    is_boundary = is_boundary
  )
}

discover_population_hard_blocks <- function(theta,
                                            w,
                                            population_model = NULL,
                                            max_block_dim = 2L,
                                            max_single_blocks = 4L,
                                            max_pair_blocks = 4L,
                                            single_score_threshold = 0.15,
                                            pair_score_threshold = 0.25,
                                            min_abs_corr = 0.10) {
  theta <- as.matrix(theta)
  hyper_names <- colnames(theta)
  if (is.null(hyper_names) && !is.null(population_model)) {
    hyper_names <- population_model$hyper_names
    colnames(theta) <- hyper_names
  }
  if (is.null(hyper_names)) {
    hyper_names <- paste0("theta_", seq_len(ncol(theta)))
    colnames(theta) <- hyper_names
  }
  w <- normalize_reference_local_weights(w)
  max_block_dim <- as.integer(max(1L, max_block_dim))

  singles <- do.call(
    rbind,
    lapply(seq_len(ncol(theta)), function(j) {
      geom <- .weighted_marginal_geometry(theta[, j], w)
      data.frame(
        block_id = .make_hyper_block_id(hyper_names[j]),
        dim = 1L,
        hyper_names = hyper_names[j],
        indices = j,
        geometry_score = as.numeric(geom$score),
        abs_corr = NA_real_,
        skew = as.numeric(geom$skew),
        kurt_excess = as.numeric(geom$kurt_excess),
        tail_asymmetry = as.numeric(geom$tail_asymmetry),
        core_asymmetry = as.numeric(geom$core_asymmetry),
        stringsAsFactors = FALSE
      )
    })
  )
  singles <- singles[order(-singles$geometry_score, singles$hyper_names), , drop = FALSE]
  keep_single <- singles$geometry_score >= single_score_threshold
  if (sum(keep_single) > max_single_blocks) {
    keep_single <- seq_len(nrow(singles)) %in% head(which(keep_single), max_single_blocks)
  }
  singles_selected <- singles[keep_single, , drop = FALSE]

  pairs_selected <- singles[FALSE, , drop = FALSE]
  if (max_block_dim >= 2L && ncol(theta) >= 2L && max_pair_blocks > 0L) {
    S <- tryCatch(weighted_cov(theta, w), error = function(e) NULL)
    if (!is.null(S)) {
      denom <- sqrt(pmax(diag(S), 1e-12))
      corr <- S / tcrossprod(denom)
      corr[!is.finite(corr)] <- 0
      single_scores <- singles$geometry_score
      names(single_scores) <- singles$hyper_names
      pair_rows <- lapply(utils::combn(seq_len(ncol(theta)), 2L, simplify = FALSE), function(idx) {
        h1 <- hyper_names[idx[1L]]
        h2 <- hyper_names[idx[2L]]
        abs_corr <- abs(corr[idx[1L], idx[2L]])
        pair_score <- abs_corr * (1 + 0.5 * (single_scores[h1] + single_scores[h2]))
        data.frame(
          block_id = .make_hyper_block_id(c(h1, h2)),
          dim = 2L,
          hyper_names = paste(c(h1, h2), collapse = ","),
          indices = paste(idx, collapse = ","),
          geometry_score = as.numeric(pair_score),
          abs_corr = as.numeric(abs_corr),
          skew = NA_real_,
          kurt_excess = NA_real_,
          tail_asymmetry = NA_real_,
          core_asymmetry = NA_real_,
          stringsAsFactors = FALSE
        )
      })
      pairs <- do.call(rbind, pair_rows)
      if (nrow(pairs)) {
        pairs <- pairs[order(-pairs$geometry_score, -pairs$abs_corr, pairs$hyper_names), , drop = FALSE]
        keep_pair <- pairs$geometry_score >= pair_score_threshold & pairs$abs_corr >= min_abs_corr
        if (sum(keep_pair) > max_pair_blocks) {
          keep_pair <- seq_len(nrow(pairs)) %in% head(which(keep_pair), max_pair_blocks)
        }
        pairs_selected <- pairs[keep_pair, , drop = FALSE]
      }
    }
  }

  selected <- rbind(singles_selected, pairs_selected)
  if (!nrow(selected)) {
    return(list(summary = selected, blocks = list()))
  }

  blocks <- lapply(seq_len(nrow(selected)), function(i) {
    idx <- if (selected$dim[i] == 1L) {
      as.integer(selected$indices[i])
    } else {
      as.integer(strsplit(selected$indices[i], ",", fixed = TRUE)[[1L]])
    }
    list(
      id = selected$block_id[i],
      dim = as.integer(selected$dim[i]),
      indices = idx,
      hyper_names = hyper_names[idx],
      geometry_score = as.numeric(selected$geometry_score[i]),
      abs_corr = as.numeric(selected$abs_corr[i])
    )
  })
  names(blocks) <- selected$block_id
  list(summary = selected, blocks = blocks)
}

population_local_site_block_diagnostic <- function(site,
                                                   theta,
                                                   w,
                                                   block,
                                                   theta_center = NULL,
                                                   probs_1d = c(0.05, 0.20, 0.50, 0.80, 0.95),
                                                   probs_2d = c(0.10, 0.50, 0.90),
                                                   use_psis = TRUE) {
  stopifnot(inherits(site, "population_local_site"))
  theta <- .as_hyper_matrix(
    theta,
    hyper_names = site$population_model$hyper_names,
    hyper_dim = site$population_model$hyper_dim
  )
  w <- normalize_reference_local_weights(w)
  theta_center <- theta_center %||% matrix(colSums(theta * w), nrow = 1L)
  colnames(theta_center) <- colnames(theta)

  grid <- .build_block_grid(
    theta = theta,
    w = w,
    theta_center = theta_center,
    block = block,
    probs_1d = probs_1d,
    probs_2d = probs_2d
  )
  log_vals <- population_local_site_log_marginal_many(site, grid$theta, include_constant = TRUE)
  tail_diag <- lapply(seq_len(nrow(grid$theta)), function(i) {
    population_local_factor_tail_diagnostic(
      factor = site$base_factor,
      theta = grid$theta[i, , drop = FALSE],
      use_psis = use_psis
    )
  })
  ess <- vapply(tail_diag, `[[`, numeric(1), "ess")
  pareto_k <- vapply(tail_diag, `[[`, numeric(1), "pareto_k")
  shape <- .fit_quadratic_shape(grid$block_values, log_vals)
  n_particles <- max(1L, nrow(site$base_factor$particles))
  mode_idx <- which.max(log_vals)

  data.frame(
    local_id = as.integer(site$local_id),
    block_id = as.character(block$id),
    dim = as.integer(block$dim),
    hyper_names = paste(block$hyper_names, collapse = ","),
    min_ess = min(ess, na.rm = TRUE),
    min_ess_frac = min(ess, na.rm = TRUE) / n_particles,
    max_pareto_k = if (any(is.finite(pareto_k))) max(pareto_k, na.rm = TRUE) else NA_real_,
    scaled_rmse = as.numeric(shape$scaled_rmse),
    boundary_mode = isTRUE(grid$is_boundary[mode_idx]),
    grid_n = nrow(grid$theta),
    stringsAsFactors = FALSE
  )
}

summarize_population_site_block_diagnostics <- function(local_diag) {
  finite_k <- local_diag$max_pareto_k[is.finite(local_diag$max_pareto_k)]
  finite_rmse <- local_diag$scaled_rmse[is.finite(local_diag$scaled_rmse)]
  data.frame(
    min_ess_frac = min(local_diag$min_ess_frac, na.rm = TRUE),
    q10_ess_frac = as.numeric(stats::quantile(local_diag$min_ess_frac, probs = 0.10, na.rm = TRUE, names = FALSE)),
    mean_ess_frac = mean(local_diag$min_ess_frac, na.rm = TRUE),
    max_pareto_k = if (length(finite_k)) max(finite_k) else NA_real_,
    q90_pareto_k = if (length(finite_k)) as.numeric(stats::quantile(finite_k, probs = 0.90, na.rm = TRUE, names = FALSE)) else NA_real_,
    boundary_rate = mean(local_diag$boundary_mode, na.rm = TRUE),
    mean_scaled_rmse = if (length(finite_rmse)) mean(finite_rmse) else NA_real_,
    q90_scaled_rmse = if (length(finite_rmse)) as.numeric(stats::quantile(finite_rmse, probs = 0.90, na.rm = TRUE, names = FALSE)) else NA_real_,
    stringsAsFactors = FALSE
  )
}

detect_population_hard_blocks <- function(site_set,
                                          theta,
                                          w,
                                          theta_center = NULL,
                                          max_block_dim = 2L,
                                          max_single_blocks = 4L,
                                          max_pair_blocks = 4L,
                                          single_score_threshold = 0.15,
                                          pair_score_threshold = 0.25,
                                          min_abs_corr = 0.10,
                                          probs_1d = c(0.05, 0.20, 0.50, 0.80, 0.95),
                                          probs_2d = c(0.10, 0.50, 0.90),
                                          use_psis = TRUE,
                                          pareto_k_threshold = 0.70,
                                          ess_frac_threshold = 0.05,
                                          rmse_threshold = 0.15,
                                          boundary_threshold = 0.25) {
  stopifnot(inherits(site_set, "population_site_set"))
  theta <- .as_hyper_matrix(
    theta,
    hyper_names = site_set$population_model$hyper_names,
    hyper_dim = site_set$population_model$hyper_dim
  )
  w <- normalize_reference_local_weights(w)
  theta_center <- theta_center %||% matrix(colSums(theta * w), nrow = 1L)
  colnames(theta_center) <- colnames(theta)

  discovered <- discover_population_hard_blocks(
    theta = theta,
    w = w,
    population_model = site_set$population_model,
    max_block_dim = max_block_dim,
    max_single_blocks = max_single_blocks,
    max_pair_blocks = max_pair_blocks,
    single_score_threshold = single_score_threshold,
    pair_score_threshold = pair_score_threshold,
    min_abs_corr = min_abs_corr
  )
  if (!length(discovered$blocks)) {
    return(list(
      theta_center = theta_center,
      candidates = discovered$summary,
      block_diagnostics = list(),
      summary = data.frame(),
      hard_blocks = data.frame()
    ))
  }

  block_diagnostics <- lapply(discovered$blocks, function(block) {
    local_diag <- do.call(
      rbind,
      lapply(
        site_set$sites,
        population_local_site_block_diagnostic,
        theta = theta,
        w = w,
        block = block,
        theta_center = theta_center,
        probs_1d = probs_1d,
        probs_2d = probs_2d,
        use_psis = use_psis
      )
    )
    list(
      block = block,
      local = local_diag,
      summary = summarize_population_site_block_diagnostics(local_diag)
    )
  })
  names(block_diagnostics) <- names(discovered$blocks)

  summary_rows <- lapply(names(block_diagnostics), function(block_id) {
    diag_block <- block_diagnostics[[block_id]]
    geom_row <- discovered$summary[discovered$summary$block_id == block_id, , drop = FALSE]
    summ <- diag_block$summary
    site_score <- max(0, as.numeric(summ$max_pareto_k) - 0.5) +
      4 * max(0, ess_frac_threshold - as.numeric(summ$q10_ess_frac)) +
      2 * max(0, as.numeric(summ$q90_scaled_rmse) - 0.10) +
      max(0, as.numeric(summ$boundary_rate) - 0.10)
    hard <- (is.finite(summ$max_pareto_k) && summ$max_pareto_k > pareto_k_threshold) ||
      (is.finite(summ$q10_ess_frac) && summ$q10_ess_frac < ess_frac_threshold) ||
      (is.finite(summ$q90_scaled_rmse) && summ$q90_scaled_rmse > rmse_threshold) ||
      (
        is.finite(summ$boundary_rate) &&
          summ$boundary_rate > boundary_threshold &&
          (
            as.numeric(geom_row$geometry_score[1L]) > 0.5 ||
              (is.finite(summ$max_pareto_k) && summ$max_pareto_k > 0.7) ||
              (is.finite(summ$q10_ess_frac) && summ$q10_ess_frac < 0.10)
          )
      )
    data.frame(
      block_id = block_id,
      dim = as.integer(diag_block$block$dim),
      hyper_names = paste(diag_block$block$hyper_names, collapse = ","),
      geometry_score = as.numeric(geom_row$geometry_score[1L]),
      abs_corr = as.numeric(geom_row$abs_corr[1L]),
      min_ess_frac = as.numeric(summ$min_ess_frac),
      q10_ess_frac = as.numeric(summ$q10_ess_frac),
      mean_ess_frac = as.numeric(summ$mean_ess_frac),
      max_pareto_k = as.numeric(summ$max_pareto_k),
      q90_pareto_k = as.numeric(summ$q90_pareto_k),
      boundary_rate = as.numeric(summ$boundary_rate),
      mean_scaled_rmse = as.numeric(summ$mean_scaled_rmse),
      q90_scaled_rmse = as.numeric(summ$q90_scaled_rmse),
      site_score = as.numeric(site_score),
      total_score = as.numeric(geom_row$geometry_score[1L] + site_score),
      hard = isTRUE(hard),
      stringsAsFactors = FALSE
    )
  })
  summary_df <- do.call(rbind, summary_rows)
  summary_df <- summary_df[order(-summary_df$total_score, -summary_df$geometry_score), , drop = FALSE]

  list(
    theta_center = theta_center,
    candidates = discovered$summary,
    block_diagnostics = block_diagnostics,
    summary = summary_df,
    hard_blocks = summary_df[summary_df$hard, , drop = FALSE]
  )
}

print_population_hard_blocks <- function(block_diag, title = "Detected hard hyper-blocks:") {
  cat(title, "\n", sep = "")
  if (is.null(block_diag) || !nrow(block_diag$summary)) {
    cat("  No hard candidate blocks were detected.\n")
    return(invisible(NULL))
  }
  print(
    utils::head(
      block_diag$summary[, c(
        "block_id",
        "dim",
        "geometry_score",
        "max_pareto_k",
        "q10_ess_frac",
        "boundary_rate",
        "q90_scaled_rmse",
        "total_score",
        "hard"
      ), drop = FALSE],
      8L
    )
  )
  invisible(NULL)
}

population_factor_set_local_tail_diagnostics <- function(factor_set, theta, use_psis = TRUE) {
  stopifnot(inherits(factor_set, "population_factor_set"))
  model <- factor_set$population_model
  theta <- .as_hyper_matrix(theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  if (nrow(theta) != 1L) {
    stop("population_factor_set_local_tail_diagnostics expects a single theta row.")
  }
  rows <- lapply(seq_along(factor_set$factors), function(i) {
    diag_i <- population_local_factor_tail_diagnostic(factor_set$factors[[i]], theta = theta, use_psis = use_psis)
    data.frame(
      local = names(factor_set$factors)[i] %||% as.character(i),
      ess = as.numeric(diag_i$ess),
      pareto_k = as.numeric(diag_i$pareto_k),
      max_weight = as.numeric(diag_i$max_weight),
      q99_weight = as.numeric(diag_i$q99_weight),
      n_finite = as.integer(diag_i$n_finite),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

population_factor_set_sv_tail_grid <- function(factor_set,
                                               theta_center,
                                               sv_grid,
                                               log_sigma2_name = "log_sigma2_sv",
                                               use_psis = TRUE,
                                               top_n = 5L) {
  stopifnot(inherits(factor_set, "population_factor_set"))
  model <- factor_set$population_model
  theta_center <- .as_hyper_matrix(theta_center, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  if (nrow(theta_center) != 1L) {
    stop("theta_center must contain a single theta row.")
  }
  if (!(log_sigma2_name %in% colnames(theta_center))) {
    stop("log_sigma2_name is not a column in theta_center.")
  }

  summary_rows <- vector("list", length(sv_grid))
  worst_rows <- vector("list", length(sv_grid))

  for (idx in seq_along(sv_grid)) {
    theta <- theta_center
    theta[1L, log_sigma2_name] <- as.numeric(sv_grid[idx])
    local_diag <- population_factor_set_local_tail_diagnostics(
      factor_set = factor_set,
      theta = theta,
      use_psis = use_psis
    )
    finite_k <- local_diag$pareto_k[is.finite(local_diag$pareto_k)]
    summary_rows[[idx]] <- data.frame(
      log_sigma2_sv = as.numeric(sv_grid[idx]),
      sigma2_sv = exp(as.numeric(sv_grid[idx])),
      min_ess = min(local_diag$ess, na.rm = TRUE),
      q10_ess = as.numeric(stats::quantile(local_diag$ess, probs = 0.10, na.rm = TRUE, names = FALSE)),
      mean_ess = mean(local_diag$ess, na.rm = TRUE),
      max_pareto_k = if (length(finite_k)) max(finite_k) else NA_real_,
      q90_pareto_k = if (length(finite_k)) as.numeric(stats::quantile(finite_k, probs = 0.90, na.rm = TRUE, names = FALSE)) else NA_real_,
      n_k_gt_0_7 = sum(local_diag$pareto_k > 0.7, na.rm = TRUE),
      n_k_gt_1_0 = sum(local_diag$pareto_k > 1.0, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    sort_key <- ifelse(is.finite(local_diag$pareto_k), -local_diag$pareto_k, Inf)
    ord <- order(sort_key, local_diag$ess)
    worst <- local_diag[ord, , drop = FALSE]
    worst <- utils::head(worst, as.integer(max(1L, top_n)))
    worst$log_sigma2_sv <- as.numeric(sv_grid[idx])
    worst$sigma2_sv <- exp(as.numeric(sv_grid[idx]))
    worst_rows[[idx]] <- worst
  }

  list(
    summary = do.call(rbind, summary_rows),
    worst = do.call(rbind, worst_rows)
  )
}

build_population_factor_set <- function(local_objects,
                                        population_model,
                                        particle_block_size = 1024L) {
  population_model <- normalize_population_model(population_model)
  factors <- lapply(local_objects, build_population_local_factor, population_model = population_model)
  names(factors) <- names(local_objects)
  particle_block_size <- as.integer(max(1L, particle_block_size))
  if (length(factors)) {
    local_lengths <- vapply(factors, function(factor) nrow(factor$particles), integer(1))
    alpha <- do.call(rbind, lapply(factors, `[[`, "particles"))
    log_base <- unlist(lapply(factors, `[[`, "log_base"), use.names = FALSE)
    local_index <- rep.int(seq_along(factors), local_lengths)
    block_starts <- seq.int(1L, nrow(alpha), by = particle_block_size)
    blocks <- lapply(block_starts, function(start) {
      idx <- seq.int(start, min(start + particle_block_size - 1L, nrow(alpha)))
      rr <- rle(local_index[idx])
      ends <- cumsum(rr$lengths)
      starts <- c(1L, head(ends, -1L) + 1L)
      list(
        idx = idx,
        local = as.integer(rr$values),
        starts = as.integer(starts),
        ends = as.integer(ends)
      )
    })
  } else {
    alpha <- matrix(numeric(0), nrow = 0L, ncol = population_model$alpha_dim)
    colnames(alpha) <- population_model$alpha_names
    log_base <- numeric(0)
    local_index <- integer(0)
    blocks <- list()
  }
  structure(
    list(
      population_model = population_model,
      factors = factors,
      n_locals = length(factors),
      log_constant = sum(vapply(factors, `[[`, numeric(1), "log_constant")),
      stack = list(
        alpha = alpha,
        log_base = log_base,
        local_index = local_index,
        blocks = blocks,
        particle_block_size = particle_block_size
      )
    ),
    class = "population_factor_set"
  )
}

build_population_site_set <- function(local_objects,
                                      population_model,
                                      particle_block_size = 1024L) {
  population_model <- normalize_population_model(population_model)
  sites <- lapply(local_objects, build_population_local_site, population_model = population_model)
  names(sites) <- names(local_objects)
  base_factor_set <- build_population_factor_set(
    local_objects = local_objects,
    population_model = population_model,
    particle_block_size = particle_block_size
  )
  structure(
    list(
      population_model = population_model,
      sites = sites,
      n_locals = length(sites),
      base_factor_set = base_factor_set,
      log_constant = base_factor_set$log_constant,
      n_corrections = 0L
    ),
    class = "population_site_set"
  )
}

population_site_set_add_correction <- function(site_set, local_key, correction) {
  stopifnot(inherits(site_set, "population_site_set"))
  local_key_chr <- as.character(local_key)
  site_idx <- match(local_key_chr, names(site_set$sites))
  if (is.na(site_idx)) {
    site_idx <- match(as.integer(local_key), vapply(site_set$sites, `[[`, integer(1), "local_id"))
  }
  if (is.na(site_idx)) {
    stop("Could not find local site for key: ", local_key_chr)
  }
  block_id <- as.character(correction$block_id %||% .make_hyper_block_id(correction$hyper_names))
  site_set$sites[[site_idx]]$corrections[[block_id]] <- correction
  site_set$n_corrections <- sum(vapply(site_set$sites, function(site) length(site$corrections), integer(1)))
  site_set
}

.evaluate_site_set <- function(site_set, theta, include_constant = FALSE, n_cores = 1L) {
  stopifnot(inherits(site_set, "population_site_set"))
  base <- population_factor_set_loglik(
    factor_set = site_set$base_factor_set,
    theta = theta,
    include_constant = include_constant,
    n_cores = n_cores
  )
  if (!length(site_set$sites) || !site_set$n_corrections) {
    return(base)
  }
  theta_mat <- .as_hyper_matrix(
    theta,
    hyper_names = site_set$population_model$hyper_names,
    hyper_dim = site_set$population_model$hyper_dim
  )
  corr_total <- rep.int(0, nrow(theta_mat))
  for (site in site_set$sites) {
    if (!length(site$corrections)) next
    for (correction in site$corrections) {
      corr_total <- corr_total + predict_site_block_surrogate(
        surrogate = correction$surrogate,
        block_values = theta_mat[, correction$hyper_names, drop = FALSE]
      )
    }
  }
  as.numeric(base + corr_total)
}

.evaluate_factor_set_stacked <- function(factor_set, theta, include_constant = FALSE) {
  stopifnot(inherits(factor_set, "population_factor_set"))
  model <- factor_set$population_model
  theta <- .as_hyper_matrix(theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  if (!length(factor_set$factors)) {
    out <- rep.int(0, nrow(theta))
    if (isTRUE(include_constant)) out <- out + factor_set$log_constant
    return(out)
  }

  theta_prepared <- population_model_prepare_theta(model, theta)
  accum <- matrix(-Inf, nrow = nrow(theta), ncol = factor_set$n_locals)

  for (block in factor_set$stack$blocks) {
    log_terms <- .population_log_alpha_given_theta_many(
      model = model,
      alpha = factor_set$stack$alpha[block$idx, , drop = FALSE],
      theta_prepared = theta_prepared
    )
    log_terms <- sweep(log_terms, 2L, factor_set$stack$log_base[block$idx], "+")

    for (seg_id in seq_along(block$local)) {
      cols <- seq.int(block$starts[seg_id], block$ends[seg_id])
      seg_lse <- if (length(cols) == 1L) log_terms[, cols] else .rowLogSumExp(log_terms[, cols, drop = FALSE])
      accum[, block$local[seg_id]] <- rlogsumexp2(accum[, block$local[seg_id]], as.numeric(seg_lse))
    }
  }

  out <- rowSums(accum)
  if (isTRUE(include_constant)) out <- out + factor_set$log_constant
  as.numeric(out)
}

population_factor_set_loglik <- function(factor_set, theta, include_constant = FALSE, n_cores = 1L) {
  stopifnot(inherits(factor_set, "population_factor_set") || inherits(factor_set, "population_site_set"))
  if (inherits(factor_set, "population_site_set")) {
    return(.evaluate_site_set(
      site_set = factor_set,
      theta = theta,
      include_constant = include_constant,
      n_cores = n_cores
    ))
  }
  model <- factor_set$population_model
  theta <- .as_hyper_matrix(theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  if (!length(factor_set$factors)) {
    out <- rep.int(0, nrow(theta))
    if (isTRUE(include_constant)) out <- out + factor_set$log_constant
    return(out)
  }

  if (as.integer(n_cores) <= 1L || nrow(theta) <= 1L) {
    return(.evaluate_factor_set_stacked(factor_set, theta, include_constant = include_constant))
  }

  split_idx <- cut(seq_len(nrow(theta)), breaks = as.integer(min(n_cores, nrow(theta))), labels = FALSE)
  batches <- split(seq_len(nrow(theta)), split_idx)
  parts <- parallel::mclapply(
    batches,
    function(idx) {
      .evaluate_factor_set_stacked(
        factor_set = factor_set,
        theta = theta[idx, , drop = FALSE],
        include_constant = FALSE
      )
    },
    mc.cores = as.integer(min(n_cores, length(batches)))
  )
  out <- numeric(nrow(theta))
  for (i in seq_along(batches)) {
    out[batches[[i]]] <- parts[[i]]
  }
  if (isTRUE(include_constant)) out <- out + factor_set$log_constant
  out
}

population_factor_set_local_ess <- function(factor_set, theta) {
  stopifnot(inherits(factor_set, "population_factor_set"))
  vapply(factor_set$factors, population_local_factor_ess, numeric(1), theta = theta)
}

population_factor_set_logposterior <- function(factor_set, theta, include_constant = TRUE, n_cores = 1L) {
  stopifnot(inherits(factor_set, "population_factor_set") || inherits(factor_set, "population_site_set"))
  model <- factor_set$population_model
  theta <- .as_hyper_matrix(theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  population_model_log_hyperprior(model, theta) +
    population_factor_set_loglik(factor_set, theta, include_constant = include_constant, n_cores = n_cores)
}

.outer_cheap_sort_order <- function(Z) {
  Z <- as.matrix(Z)
  n <- nrow(Z)
  if (n <= 1L) return(seq_len(n))
  if (ncol(Z) <= 1L) return(order(Z[, 1L]))
  order(rowSums(Z))
}

.outer_whiten_theta <- function(theta, w = NULL) {
  theta <- as.matrix(theta)
  n <- nrow(theta)
  d <- ncol(theta)
  if (n <= 1L || d <= 0L) return(theta)

  if (is.null(w)) {
    w <- rep(1 / n, n)
  } else {
    w <- pmax(as.numeric(w), 0)
    sw <- sum(w)
    w <- if (!is.finite(sw) || sw <= 0) rep(1 / n, n) else w / sw
  }

  mu <- colSums(theta * w)
  S <- tryCatch(weighted_cov(theta, w), error = function(e) NULL)
  if (is.null(S) || any(!is.finite(S))) {
    S <- tryCatch(stats::cov(theta), error = function(e) NULL)
  }
  if (is.null(S) || any(!is.finite(S))) {
    S <- diag(d)
  } else {
    S <- regularize_cov(S, min_eig = 1e-8, cond_cap = 1e8)
  }

  L <- tryCatch(chol(S), error = function(e) diag(d))
  theta_centered <- sweep(theta, 2L, mu, "-")
  Z <- t(backsolve(L, t(theta_centered), transpose = TRUE))
  Z[!is.finite(Z)] <- 0
  Z
}

.outer_resample_sort_order <- function(theta,
                                       w = NULL,
                                       ess_frac = 1,
                                       mode = c("adaptive", "hilbert", "cheap1d", "none"),
                                       hilbert_hard_ess = 0.25,
                                       hilbert_max_dim = 8L) {
  mode <- match.arg(mode)
  theta <- as.matrix(theta)
  if (nrow(theta) <= 1L) {
    return(list(order = seq_len(nrow(theta)), sort_used = "none"))
  }
  if (mode == "none") {
    return(list(order = seq_len(nrow(theta)), sort_used = "none"))
  }

  Z <- .outer_whiten_theta(theta, w = w)
  want_hilbert <- identical(mode, "hilbert") ||
    (identical(mode, "adaptive") &&
       ess_frac <= as.numeric(hilbert_hard_ess) &&
       ncol(Z) <= as.integer(hilbert_max_dim))

  if (isTRUE(want_hilbert)) {
    ord <- tryCatch(hilbert_sort_order(Z, bits = 16L), error = function(e) NULL)
    if (!is.null(ord)) {
      return(list(order = ord, sort_used = "hilbert"))
    }
  }

  list(order = .outer_cheap_sort_order(Z), sort_used = "cheap1d")
}

.outer_rejuvenate <- function(theta,
                              logprior,
                              loglik_dynamic,
                              beta,
                              factor_set,
                              w = NULL,
                              min_n_moves = 1L,
                              max_n_moves = 3L,
                              target_accepted_moves = NULL,
                              target_accepted_moves_resampled = 1.5,
                              target_accepted_moves_not_resampled = 0.75,
                              rw_scale = 0.8,
                              target_accept = 0.234,
                              rm_gain_move = 0.15,
                              n_cores = 1L,
                              resampled = FALSE,
                              seed = NULL) {
  if (!is.null(seed)) set.seed(as.integer(seed))
  theta <- as.matrix(theta)
  N <- nrow(theta)
  d <- ncol(theta)
  model <- factor_set$population_model
  if (N <= 1L) {
    return(list(
      theta = theta,
      logprior = logprior,
      loglik_dynamic = loglik_dynamic,
      accept_rate = 0,
      move_accept = numeric(0),
      n_moves_used = 0L,
      cumulative_accept = 0,
      rw_scale = as.numeric(rw_scale)
    ))
  }

  if (is.null(w)) {
    w_cov <- rep(1 / N, N)
  } else {
    w_cov <- pmax(as.numeric(w), 0)
    sw <- sum(w_cov)
    w_cov <- if (!is.finite(sw) || sw <= 0) rep(1 / N, N) else w_cov / sw
  }

  S <- tryCatch(weighted_cov(theta, w_cov), error = function(e) NULL)
  if (is.null(S) || any(!is.finite(S))) {
    S <- tryCatch(stats::cov(theta), error = function(e) NULL)
  }
  if (is.null(S) || any(!is.finite(S))) {
    S <- diag(d)
  } else {
    S <- regularize_cov(S, min_eig = 1e-8, cond_cap = 1e8)
  }
  Lrw <- tryCatch(chol(S), error = function(e) diag(d))
  log_rw_scale <- log(pmax(as.numeric(rw_scale), 1e-6))
  min_n_moves <- as.integer(max(1L, min_n_moves))
  max_n_moves <- as.integer(max(min_n_moves, max_n_moves))
  target_accepted_moves <- as.numeric(
    target_accepted_moves %||%
      if (isTRUE(resampled)) target_accepted_moves_resampled else target_accepted_moves_not_resampled
  )

  acc_total <- 0L
  prop_total <- 0L
  logpost <- logprior + beta * loglik_dynamic
  move_accept <- numeric(0)
  cumulative_accept <- 0

  for (move_id in seq_len(max_n_moves)) {
    step_scale <- exp(log_rw_scale) / sqrt(max(d, 1))
    theta_prop <- theta + step_scale * (matrix(rnorm(N * d), nrow = N, ncol = d) %*% Lrw)
    colnames(theta_prop) <- model$hyper_names

    logprior_prop <- population_model_log_hyperprior(model, theta_prop)
    ok <- is.finite(logprior_prop)
    loglik_prop <- rep(-Inf, N)
    if (any(ok)) {
      loglik_prop[ok] <- population_factor_set_loglik(
        factor_set,
        theta = theta_prop[ok, , drop = FALSE],
        include_constant = FALSE,
        n_cores = n_cores
      )
    }
    logpost_prop <- logprior_prop + beta * loglik_prop
    log_alpha <- logpost_prop - logpost
    accept <- which(log(runif(N)) < pmin(0, log_alpha))
    prop_total <- prop_total + N
    accept_rate_move <- length(accept) / N
    move_accept <- c(move_accept, accept_rate_move)
    cumulative_accept <- cumulative_accept + accept_rate_move
    log_rw_scale <- .clamp(
      log_rw_scale + as.numeric(rm_gain_move) * (accept_rate_move - as.numeric(target_accept)),
      log(0.05),
      log(4.0)
    )

    if (length(accept)) {
      theta[accept, ] <- theta_prop[accept, , drop = FALSE]
      logprior[accept] <- logprior_prop[accept]
      loglik_dynamic[accept] <- loglik_prop[accept]
      logpost[accept] <- logpost_prop[accept]
      acc_total <- acc_total + length(accept)
    }

    if (move_id >= min_n_moves && cumulative_accept >= target_accepted_moves) {
      break
    }
  }

  list(
    theta = theta,
    logprior = logprior,
    loglik_dynamic = loglik_dynamic,
    accept_rate = if (prop_total > 0L) acc_total / prop_total else 0,
    move_accept = move_accept,
    n_moves_used = length(move_accept),
    cumulative_accept = cumulative_accept,
    rw_scale = exp(log_rw_scale)
  )
}

outer_population_smc <- function(factor_set,
                                 N = 2000L,
                                 resample_threshold = 0.5,
                                 n_mcmc_moves = 3L,
                                 min_mcmc_moves = 1L,
                                 max_rounds = 100L,
                                 beta_target = 1.0,
                                 rw_scale_init = 0.8,
                                 target_accept = 0.234,
                                 rm_gain = 0.15,
                                 target_accepted_moves = NULL,
                                 target_accepted_moves_resampled = 1.5,
                                 target_accepted_moves_not_resampled = 0.75,
                                 deterministic_resampling = FALSE,
                                 resample_sort_mode = c("adaptive", "hilbert", "cheap1d", "none"),
                                 hilbert_hard_ess = 0.25,
                                 hilbert_max_dim = 8L,
                                 n_cores = 1L,
                                 seed = 123L,
                                 verbose = TRUE) {
  stopifnot(inherits(factor_set, "population_factor_set") || inherits(factor_set, "population_site_set"))
  model <- factor_set$population_model
  vcat <- function(...) if (isTRUE(verbose)) base::cat(...)

  set.seed(as.integer(seed))
  beta_target <- as.numeric(beta_target)
  if (!is.finite(beta_target) || beta_target <= 0 || beta_target > 1) {
    stop("beta_target must be in (0, 1].")
  }
  resample_sort_mode <- match.arg(resample_sort_mode)
  min_mcmc_moves <- as.integer(max(1L, min_mcmc_moves))
  n_mcmc_moves <- as.integer(max(min_mcmc_moves, n_mcmc_moves))

  theta <- population_model_sample_hyper(model, n = as.integer(N))
  colnames(theta) <- model$hyper_names
  logprior <- population_model_log_hyperprior(model, theta)
  loglik_dynamic <- population_factor_set_loglik(factor_set, theta, include_constant = FALSE, n_cores = n_cores)

  w <- rep(1 / N, N)
  beta <- 0.0
  round <- 0L
  log_evidence_dynamic <- 0.0
  mcse_var_accum <- 0.0
  log_rw_scale <- log(rw_scale_init)

  beta_hist <- beta
  ess_hist <- 1.0
  accept_hist <- numeric(0)
  rw_scale_hist <- exp(log_rw_scale)
  resampled_hist <- logical(0)
  moves_hist <- integer(0)
  sort_hist <- character(0)

  while (beta < beta_target - 1e-12 && round < as.integer(max_rounds)) {
    round <- round + 1L
    target_cess <- cess_target_at_lambda(beta)
    beta_new <- next_lambda_via_rCESS(
      w = w,
      loglik = loglik_dynamic,
      lambda = beta,
      target = target_cess,
      lambda_target = beta_target
    )
    if (beta_new <= beta) {
      beta_new <- min(beta_target, beta + min(1e-4, beta_target - beta))
    }

    delta <- beta_new - beta
    x <- as.numeric(loglik_dynamic)
    ok <- is.finite(x)
    if (!any(ok)) stop("All outer log-likelihood values are non-finite.")
    x[!ok] <- min(x[ok])
    mx <- max(x)
    u <- exp(delta * (x - mx))
    mu1 <- sum(w * u)
    mu2 <- sum(w * u * u)
    neff <- 1 / sum(w * w)
    mcse_var_accum <- mcse_var_accum + max((mu2 - mu1^2) / (max(neff, 1) * max(mu1^2, .Machine$double.eps)), 0)

    logw_raw <- log(pmax(w, .Machine$double.eps)) + delta * (x - mx)
    lse <- logsumexp(logw_raw)
    log_evidence_dynamic <- log_evidence_dynamic + delta * mx + lse
    w <- exp(logw_raw - lse)
    w <- w / sum(w)
    ess_now <- ESS(w) / length(w)
    resampled <- FALSE

    vcat(sprintf("\nOuter round %d: beta %.3f -> %.3f | ESS=%.3f | logZ_dyn=%.4f\n",
                 round, beta, beta_new, ess_now, log_evidence_dynamic))

    beta <- beta_new
    sort_used <- "not_used"
    if (ess_now < resample_threshold) {
      sort_info <- .outer_resample_sort_order(
        theta = theta,
        w = w,
        ess_frac = ess_now,
        mode = resample_sort_mode,
        hilbert_hard_ess = hilbert_hard_ess,
        hilbert_max_dim = hilbert_max_dim
      )
      ord <- sort_info$order
      sort_used <- sort_info$sort_used
      idx_sorted <- stratified_resample_sorted(w[ord], deterministic = deterministic_resampling)
      idx <- ord[idx_sorted]
      theta <- theta[idx, , drop = FALSE]
      logprior <- logprior[idx]
      loglik_dynamic <- loglik_dynamic[idx]
      w <- rep(1 / N, N)
      resampled <- TRUE
    }

    rejuvenated <- .outer_rejuvenate(
      theta = theta,
      logprior = logprior,
      loglik_dynamic = loglik_dynamic,
      beta = beta,
      factor_set = factor_set,
      w = w,
      min_n_moves = min_mcmc_moves,
      max_n_moves = n_mcmc_moves,
      target_accepted_moves = target_accepted_moves,
      target_accepted_moves_resampled = target_accepted_moves_resampled,
      target_accepted_moves_not_resampled = target_accepted_moves_not_resampled,
      rw_scale = exp(log_rw_scale),
      target_accept = target_accept,
      rm_gain_move = rm_gain,
      n_cores = n_cores,
      resampled = resampled,
      seed = seed + 1009L * round
    )
    theta <- rejuvenated$theta
    logprior <- rejuvenated$logprior
    loglik_dynamic <- rejuvenated$loglik_dynamic
    accept_rate <- rejuvenated$accept_rate
    log_rw_scale <- log(pmax(rejuvenated$rw_scale, 1e-6))

    beta_hist <- c(beta_hist, beta)
    ess_hist <- c(ess_hist, ess_now)
    accept_hist <- c(accept_hist, accept_rate)
    rw_scale_hist <- c(rw_scale_hist, exp(log_rw_scale))
    resampled_hist <- c(resampled_hist, resampled)
    moves_hist <- c(moves_hist, rejuvenated$n_moves_used)
    sort_hist <- c(sort_hist, sort_used)

    vcat(sprintf("  Rejuvenation accept=%.3f | moves=%d | rw_scale=%.3f%s%s\n",
                 accept_rate,
                 rejuvenated$n_moves_used,
                 exp(log_rw_scale),
                 if (resampled) " | resampled" else "",
                 if (identical(sort_used, "not_used")) "" else paste0(" | sort=", sort_used)))
  }

  if (beta < beta_target - 1e-12) {
    warning("outer_population_smc hit max_rounds before reaching beta_target.")
  }

  list(
    theta = theta,
    w = w,
    logprior = logprior,
    loglik_dynamic = loglik_dynamic,
    beta = beta,
    log_evidence_dynamic = log_evidence_dynamic,
    log_evidence_constant = factor_set$log_constant,
    log_evidence = log_evidence_dynamic + factor_set$log_constant,
    mcse_log_evidence = sqrt(mcse_var_accum),
    population_model = model,
    meta = list(
      rounds = round,
      beta_hist = beta_hist,
      ess_hist = ess_hist,
      accept_hist = accept_hist,
      rw_scale_hist = rw_scale_hist,
      resampled_hist = resampled_hist,
      moves_hist = moves_hist,
      sort_hist = sort_hist
    )
  )
}
