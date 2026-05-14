#!/usr/bin/env Rscript
# ============================================================================
# Population model helpers
# - Simple population-model interface for outer SMC
# - Default diagonal Gaussian random-effects model with hyperpriors
# ============================================================================

if (!exists("%||%", mode = "function")) {
  source("smc_core.R")
}

.dnorm_log_scalar <- function(x, mean, var) {
  -0.5 * (log(2 * pi * var) + (x - mean)^2 / var)
}

.dinvgamma_log_scalar <- function(x, shape, rate) {
  ifelse(x > 0, shape * log(rate) - lgamma(shape) - (shape + 1) * log(x) - rate / x, -Inf)
}

.weighted_quantile <- function(x, w, probs = c(0.025, 0.5, 0.975)) {
  ord <- order(x)
  x_ord <- x[ord]
  w_ord <- pmax(as.numeric(w[ord]), 0)
  sw <- sum(w_ord)
  if (!is.finite(sw) || sw <= 0) {
    w_ord <- rep(1 / length(w_ord), length(w_ord))
  } else {
    w_ord <- w_ord / sw
  }
  cw <- cumsum(w_ord)
  as.numeric(stats::approx(cw, x_ord, xout = probs, rule = 2)$y)
}

.as_hyper_matrix <- function(theta, hyper_names = NULL, hyper_dim = NULL) {
  theta <- as.matrix(theta)
  if (nrow(theta) == 1L && ncol(theta) == 1L && !is.null(hyper_dim) && hyper_dim > 1L && length(theta) == hyper_dim) {
    theta <- matrix(as.numeric(theta), nrow = 1L, ncol = hyper_dim)
  }
  if (!is.null(hyper_dim) && ncol(theta) != hyper_dim) {
    if (length(theta) == hyper_dim) {
      theta <- matrix(as.numeric(theta), nrow = 1L, ncol = hyper_dim)
    } else {
      stop("Theta has incompatible dimension for the population model.")
    }
  }
  if (!is.null(hyper_names)) {
    if (!is.null(colnames(theta)) && setequal(colnames(theta), hyper_names)) {
      theta <- theta[, hyper_names, drop = FALSE]
    } else {
      colnames(theta) <- hyper_names
    }
  }
  theta
}

prepare_gaussian_theta_diag <- function(mean,
                                        sigma2,
                                        alpha_names = colnames(mean)) {
  mean <- as.matrix(mean)
  sigma2 <- as.matrix(sigma2)
  n_theta <- nrow(mean)
  d <- ncol(mean)

  if (ncol(sigma2) != d) {
    stop("sigma2 dimensions do not match mean.")
  }
  if (nrow(sigma2) == 1L && n_theta > 1L) {
    sigma2 <- sigma2[rep.int(1L, n_theta), , drop = FALSE]
  }
  if (nrow(sigma2) != n_theta) {
    stop("sigma2 must provide one row per theta row.")
  }

  sigma2 <- pmax(sigma2, 1e-12)
  inv_sigma2 <- 1 / sigma2
  eta <- mean * inv_sigma2

  out <- list(
    family = "gaussian",
    theta = mean,
    mean = mean,
    eta = eta,
    quadratic_kind = "diag",
    quadratic_coef = inv_sigma2,
    log_kernel_constant = -0.5 * (
      d * log(2 * pi) +
        rowSums(log(sigma2)) +
        rowSums(mean * eta)
    ),
    alpha_dim = d,
    alpha_names = alpha_names
  )
  colnames(out$eta) <- alpha_names
  colnames(out$quadratic_coef) <- alpha_names
  out
}

normalize_population_model <- function(model) {
  if (!is.list(model)) stop("population model must be a list.")
  required <- c("name", "alpha_dim", "hyper_dim", "alpha_names", "hyper_names", "sample_hyper", "log_hyperprior", "log_alpha_given_theta")
  missing <- setdiff(required, names(model))
  if (length(missing)) {
    stop("population model is missing: ", paste(missing, collapse = ", "))
  }
  if (!is.function(model$sample_hyper) || !is.function(model$log_hyperprior) || !is.function(model$log_alpha_given_theta)) {
    stop("population model functions are invalid.")
  }
  model$alpha_dim <- as.integer(model$alpha_dim)
  model$hyper_dim <- as.integer(model$hyper_dim)
  if (length(model$alpha_names) != model$alpha_dim) stop("alpha_names length does not match alpha_dim.")
  if (length(model$hyper_names) != model$hyper_dim) stop("hyper_names length does not match hyper_dim.")
  structure(model, class = c("population_model", class(model)))
}

population_model_sample_hyper <- function(model, n, seed = NULL) {
  model <- normalize_population_model(model)
  if (!is.null(seed)) set.seed(as.integer(seed))
  theta <- as.matrix(model$sample_hyper(as.integer(n)))
  .as_hyper_matrix(theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
}

population_model_log_hyperprior <- function(model, theta) {
  model <- normalize_population_model(model)
  theta <- .as_hyper_matrix(theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  out <- as.numeric(model$log_hyperprior(theta))
  if (length(out) != nrow(theta)) {
    stop("log_hyperprior must return one value per theta row.")
  }
  out
}

population_model_log_alpha_given_theta <- function(model, alpha, theta) {
  model <- normalize_population_model(model)
  alpha <- as.matrix(alpha)
  if (ncol(alpha) != model$alpha_dim) stop("alpha dimension does not match the population model.")
  if (!is.null(colnames(alpha)) && setequal(colnames(alpha), model$alpha_names)) {
    alpha <- alpha[, model$alpha_names, drop = FALSE]
  } else {
    colnames(alpha) <- model$alpha_names
  }
  theta <- .as_hyper_matrix(theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  if (nrow(theta) != 1L) stop("population_model_log_alpha_given_theta expects a single theta row.")
  out <- as.numeric(model$log_alpha_given_theta(alpha, theta[1L, , drop = FALSE]))
  if (length(out) != nrow(alpha)) {
    stop("log_alpha_given_theta must return one value per alpha row.")
  }
  out
}

population_model_log_alpha_given_theta_many <- function(model, alpha, theta) {
  model <- normalize_population_model(model)
  alpha <- as.matrix(alpha)
  theta <- .as_hyper_matrix(theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  if (ncol(alpha) != model$alpha_dim) stop("alpha dimension does not match the population model.")
  if (!is.null(colnames(alpha)) && setequal(colnames(alpha), model$alpha_names)) {
    alpha <- alpha[, model$alpha_names, drop = FALSE]
  } else {
    colnames(alpha) <- model$alpha_names
  }

  if (is.function(model$log_alpha_given_theta_many)) {
    out <- model$log_alpha_given_theta_many(alpha, theta)
    out <- as.matrix(out)
    if (nrow(out) != nrow(theta) || ncol(out) != nrow(alpha)) {
      stop("log_alpha_given_theta_many must return an n_theta x n_alpha matrix.")
    }
    return(out)
  }

  out <- t(vapply(
    seq_len(nrow(theta)),
    function(i) population_model_log_alpha_given_theta(model, alpha, theta[i, , drop = FALSE]),
    numeric(nrow(alpha))
  ))
  rownames(out) <- rownames(theta)
  out
}

population_model_prepare_theta <- function(model, theta) {
  model <- normalize_population_model(model)
  theta <- .as_hyper_matrix(theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)

  if (!is.function(model$prepare_theta)) {
    return(list(family = "generic", theta = theta))
  }

  out <- model$prepare_theta(theta)
  if (!is.list(out)) {
    stop("prepare_theta() must return a list.")
  }

  out$family <- as.character(out$family %||% "generic")
  out$theta <- .as_hyper_matrix(out$theta %||% theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)

  if (identical(out$family, "gaussian")) {
    out$eta <- as.matrix(out$eta)
    out$quadratic_kind <- as.character(out$quadratic_kind %||% "diag")
    if (!identical(out$quadratic_kind, "diag")) {
      stop("Only diagonal Gaussian population models are supported in the fast prepared path.")
    }
    out$quadratic_coef <- as.matrix(out$quadratic_coef)
    out$log_kernel_constant <- as.numeric(out$log_kernel_constant)
    if (nrow(out$eta) != nrow(theta) ||
        nrow(out$quadratic_coef) != nrow(theta) ||
        ncol(out$eta) != model$alpha_dim ||
        ncol(out$quadratic_coef) != model$alpha_dim ||
        length(out$log_kernel_constant) != nrow(theta)) {
      stop("Gaussian prepare_theta() outputs must provide one row per theta row.")
    }
  }

  out
}

population_model_log_alpha_given_prepared_theta_many <- function(model, alpha, theta_prepared) {
  model <- normalize_population_model(model)
  alpha <- as.matrix(alpha)
  if (is.function(model$log_alpha_given_prepared_theta_many)) {
    out <- model$log_alpha_given_prepared_theta_many(alpha, theta_prepared)
    out <- as.matrix(out)
    if (nrow(out) != nrow(theta_prepared$theta) || ncol(out) != nrow(alpha)) {
      stop("log_alpha_given_prepared_theta_many must return an n_theta x n_alpha matrix.")
    }
    return(out)
  }
  population_model_log_alpha_given_theta_many(model, alpha = alpha, theta = theta_prepared$theta)
}

population_model_reference_components_from_theta <- function(model, theta) {
  model <- normalize_population_model(model)
  theta <- .as_hyper_matrix(theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  if (!is.function(model$reference_components_from_theta)) {
    stop("population model does not define reference_components_from_theta().")
  }

  out <- model$reference_components_from_theta(theta)
  if (!is.list(out) ||
      is.null(out$component_means) ||
      is.null(out$component_covs)) {
    stop("reference_components_from_theta() must return a list with component_means and component_covs.")
  }
  if (length(out$component_means) != nrow(theta) || length(out$component_covs) != nrow(theta)) {
    stop("reference_components_from_theta() must return one alpha-space component per theta row.")
  }

  component_means <- lapply(out$component_means, function(mu) {
    mu <- as.numeric(mu)
    if (length(mu) != model$alpha_dim) {
      stop("reference component mean has incompatible dimension.")
    }
    names(mu) <- model$alpha_names
    mu
  })
  component_covs <- lapply(out$component_covs, function(S) {
    S <- as.matrix(S)
    if (nrow(S) != model$alpha_dim || ncol(S) != model$alpha_dim) {
      stop("reference component covariance has incompatible dimension.")
    }
    dimnames(S) <- list(model$alpha_names, model$alpha_names)
    S
  })

  list(
    component_means = component_means,
    component_covs = component_covs
  )
}

make_population_model_diag_gaussian <- function(alpha_names,
                                                mean_prior_mean,
                                                mean_prior_var,
                                                sigma2_prior_shape,
                                                sigma2_prior_rate,
                                                label = "diag_gaussian") {
  alpha_names <- as.character(alpha_names)
  d <- length(alpha_names)
  mean_prior_mean <- rep_len(as.numeric(mean_prior_mean), d)
  mean_prior_var <- rep_len(as.numeric(mean_prior_var), d)
  sigma2_prior_shape <- rep_len(as.numeric(sigma2_prior_shape), d)
  sigma2_prior_rate <- rep_len(as.numeric(sigma2_prior_rate), d)

  if (any(mean_prior_var <= 0)) stop("mean_prior_var must be strictly positive.")
  if (any(sigma2_prior_shape <= 1)) stop("sigma2_prior_shape must exceed 1.")
  if (any(sigma2_prior_rate <= 0)) stop("sigma2_prior_rate must be strictly positive.")

  mu_names <- paste0("mu_", alpha_names)
  log_sigma2_names <- paste0("log_sigma2_", alpha_names)
  hyper_names <- c(mu_names, log_sigma2_names)

  sample_hyper <- function(n) {
    mu <- matrix(
      rnorm(n * d, mean = rep(mean_prior_mean, each = n), sd = sqrt(rep(mean_prior_var, each = n))),
      nrow = n,
      ncol = d,
      byrow = FALSE
    )
    sigma2 <- matrix(
      1 / rgamma(n * d, shape = rep(sigma2_prior_shape, each = n), rate = rep(sigma2_prior_rate, each = n)),
      nrow = n,
      ncol = d,
      byrow = FALSE
    )
    out <- cbind(mu, log(sigma2))
    colnames(out) <- hyper_names
    out
  }

  log_hyperprior <- function(theta) {
    theta <- .as_hyper_matrix(theta, hyper_names = hyper_names, hyper_dim = 2L * d)
    mu <- theta[, seq_len(d), drop = FALSE]
    ell <- theta[, d + seq_len(d), drop = FALSE]
    sigma2 <- exp(ell)
    mean_log <- vapply(
      seq_len(d),
      function(j) .dnorm_log_scalar(mu[, j], mean_prior_mean[j], mean_prior_var[j]),
      numeric(nrow(theta))
    )
    sigma_log <- vapply(
      seq_len(d),
      function(j) .dinvgamma_log_scalar(sigma2[, j], sigma2_prior_shape[j], sigma2_prior_rate[j]),
      numeric(nrow(theta))
    )
    out <- rowSums(mean_log) + rowSums(sigma_log) + rowSums(ell)
    as.numeric(out)
  }

  log_alpha_given_theta <- function(alpha, theta_row) {
    alpha <- as.matrix(alpha)
    theta_row <- .as_hyper_matrix(theta_row, hyper_names = hyper_names, hyper_dim = 2L * d)
    mu <- as.numeric(theta_row[1L, seq_len(d)])
    sigma2 <- pmax(exp(as.numeric(theta_row[1L, d + seq_len(d)])), 1e-12)
    centered <- sweep(alpha, 2L, mu, "-")
    as.numeric(-0.5 * (d * log(2 * pi) + sum(log(sigma2)) + rowSums(sweep(centered^2, 2L, sigma2, "/"))))
  }

  log_alpha_given_theta_many <- function(alpha, theta) {
    alpha <- as.matrix(alpha)
    theta <- .as_hyper_matrix(theta, hyper_names = hyper_names, hyper_dim = 2L * d)
    n_theta <- nrow(theta)
    n_alpha <- nrow(alpha)
    mu <- theta[, seq_len(d), drop = FALSE]
    sigma2 <- pmax(exp(theta[, d + seq_len(d), drop = FALSE]), 1e-12)

    out <- matrix(
      -0.5 * (d * log(2 * pi) + rowSums(log(sigma2))),
      nrow = n_theta,
      ncol = n_alpha
    )
    for (j in seq_len(d)) {
      diff <- outer(mu[, j], alpha[, j], "-")
      out <- out - 0.5 * sweep(diff^2, 1L, sigma2[, j], "/")
    }
    out
  }

  reference_components_from_theta <- function(theta) {
    theta <- .as_hyper_matrix(theta, hyper_names = hyper_names, hyper_dim = 2L * d)
    mu <- theta[, seq_len(d), drop = FALSE]
    sigma2 <- pmax(exp(theta[, d + seq_len(d), drop = FALSE]), 1e-12)

    list(
      component_means = lapply(seq_len(nrow(theta)), function(i) {
        out <- as.numeric(mu[i, ])
        names(out) <- alpha_names
        out
      }),
      component_covs = lapply(seq_len(nrow(theta)), function(i) {
        out <- diag(as.numeric(sigma2[i, ]), nrow = d, ncol = d)
        dimnames(out) <- list(alpha_names, alpha_names)
        out
      })
    )
  }

  prepare_theta <- function(theta) {
    theta <- .as_hyper_matrix(theta, hyper_names = hyper_names, hyper_dim = 2L * d)
    mu <- theta[, seq_len(d), drop = FALSE]
    sigma2 <- pmax(exp(theta[, d + seq_len(d), drop = FALSE]), 1e-12)
    prepared <- prepare_gaussian_theta_diag(
      mean = mu,
      sigma2 = sigma2,
      alpha_names = alpha_names
    )
    prepared$theta <- theta
    prepared
  }

  normalize_population_model(
    list(
      name = label,
      fast_family = "gaussian",
      alpha_dim = d,
      hyper_dim = 2L * d,
      alpha_names = alpha_names,
      hyper_names = hyper_names,
      sample_hyper = sample_hyper,
      log_hyperprior = log_hyperprior,
      log_alpha_given_theta = log_alpha_given_theta,
      log_alpha_given_theta_many = log_alpha_given_theta_many,
      prepare_theta = prepare_theta,
      reference_components_from_theta = reference_components_from_theta,
      prior_spec = list(
        mean_prior_mean = mean_prior_mean,
        mean_prior_var = mean_prior_var,
        sigma2_prior_shape = sigma2_prior_shape,
        sigma2_prior_rate = sigma2_prior_rate
      ),
      theta_to_list = function(theta) {
        theta <- .as_hyper_matrix(theta, hyper_names = hyper_names, hyper_dim = 2L * d)
        list(
          mu = theta[, seq_len(d), drop = FALSE],
          sigma2 = exp(theta[, d + seq_len(d), drop = FALSE])
        )
      }
    )
  )
}

default_population_model_diag_gaussian <- function(base_mu,
                                                   base_Sigma,
                                                   mean_var_scale = 1.0,
                                                   sigma2_shape = 3.0,
                                                   sigma2_mean = NULL,
                                                   label = "diag_gaussian") {
  base_mu <- as.numeric(base_mu)
  base_Sigma <- as.matrix(base_Sigma)
  d <- length(base_mu)
  if (nrow(base_Sigma) != d || ncol(base_Sigma) != d) {
    stop("base_Sigma dimensions do not match base_mu.")
  }
  var_diag <- pmax(diag(base_Sigma), 1e-8)
  sigma2_mean <- rep_len(as.numeric(sigma2_mean %||% var_diag), d)
  sigma2_shape <- rep_len(as.numeric(sigma2_shape), d)
  sigma2_rate <- sigma2_mean * (sigma2_shape - 1)

  make_population_model_diag_gaussian(
    alpha_names = names(base_mu) %||% colnames(base_Sigma) %||% paste0("theta", seq_len(d)),
    mean_prior_mean = base_mu,
    mean_prior_var = mean_var_scale * var_diag,
    sigma2_prior_shape = sigma2_shape,
    sigma2_prior_rate = sigma2_rate,
    label = label
  )
}

summarize_population_posterior_diag <- function(theta, w, model) {
  model <- normalize_population_model(model)
  theta <- .as_hyper_matrix(theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  w <- pmax(as.numeric(w), 0)
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) {
    w <- rep(1 / nrow(theta), nrow(theta))
  } else {
    w <- w / sw
  }

  d <- model$alpha_dim
  mu_part <- theta[, seq_len(d), drop = FALSE]
  sigma2_part <- exp(theta[, d + seq_len(d), drop = FALSE])
  colnames(mu_part) <- model$alpha_names
  colnames(sigma2_part) <- model$alpha_names

  mu_q <- t(vapply(seq_len(d), function(j) .weighted_quantile(mu_part[, j], w), numeric(3)))
  sigma2_q <- t(vapply(seq_len(d), function(j) .weighted_quantile(sigma2_part[, j], w), numeric(3)))

  list(
    mu = data.frame(
      parameter = model$alpha_names,
      mean = colSums(mu_part * w),
      median = mu_q[, 2L],
      q025 = mu_q[, 1L],
      q975 = mu_q[, 3L],
      row.names = NULL,
      check.names = FALSE
    ),
    sigma2 = data.frame(
      parameter = model$alpha_names,
      mean = colSums(sigma2_part * w),
      median = sigma2_q[, 2L],
      q025 = sigma2_q[, 1L],
      q975 = sigma2_q[, 3L],
      row.names = NULL,
      check.names = FALSE
    )
  )
}
