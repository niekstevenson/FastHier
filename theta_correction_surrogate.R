#!/usr/bin/env Rscript
# ============================================================================
# Theta-level correction surrogate
# - Phase 6 support-engine building block
# - Fits and validates Delta(theta) from exact local-SMC correction points
# ============================================================================

if (!exists("%||%", mode = "function") ||
    !exists("weighted_cov", mode = "function") ||
    !exists("regularize_cov", mode = "function")) {
  source("smc_core.R")
}

if (!exists("normalize_population_model", mode = "function")) {
  source("population_models.R")
}

if (!exists("theta_pathfinder_whitener", mode = "function") &&
    file.exists("theta_support_pathfinder.R")) {
  source("theta_support_pathfinder.R")
}

.theta_sur_as_matrix <- function(theta, model) {
  .as_hyper_matrix(theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
}

.theta_sur_normalize_weights <- function(w, n) {
  if (is.null(w)) return(rep(1 / n, n))
  if (length(w) != n) stop("weights length does not match theta rows.")
  w <- pmax(as.numeric(w), 0)
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) rep(1 / n, n) else w / sw
}

.theta_sur_support_whitener <- function(theta, w, model) {
  theta <- .theta_sur_as_matrix(theta, model)
  w <- .theta_sur_normalize_weights(w, nrow(theta))
  if (exists("theta_pathfinder_whitener", mode = "function")) {
    return(theta_pathfinder_whitener(theta, w = w, model = model))
  }
  center <- colSums(theta * w)
  S <- tryCatch(weighted_cov(theta, w), error = function(e) stats::cov(theta))
  S <- regularize_cov(S, min_eig = 1e-8, cond_cap = 1e8)
  list(center = center, cov = S, chol = chol(S), names = model$hyper_names)
}

.theta_sur_to_z <- function(theta, whitener, model) {
  theta <- .theta_sur_as_matrix(theta, model)
  centered <- sweep(theta[, whitener$names, drop = FALSE], 2L, whitener$center, "-")
  z <- t(backsolve(whitener$chol, t(centered), transpose = TRUE))
  colnames(z) <- whitener$names
  z
}

.theta_sur_from_z <- function(z, whitener, model) {
  z <- as.matrix(z)
  theta <- sweep(z %*% whitener$chol, 2L, whitener$center, "+")
  colnames(theta) <- whitener$names
  .theta_sur_as_matrix(theta, model)
}

.theta_sur_design_matrix <- function(z, degree) {
  z <- as.matrix(z)
  degree <- as.integer(degree)
  if (!degree %in% 0:2) stop("degree must be 0, 1, or 2.")
  n <- nrow(z)
  d <- ncol(z)
  X <- matrix(1, nrow = n, ncol = 1L)
  colnames(X) <- "intercept"
  if (degree >= 1L) {
    linear <- z
    colnames(linear) <- paste0("z_", colnames(z) %||% seq_len(d))
    X <- cbind(X, linear)
  }
  if (degree >= 2L) {
    quad <- vector("list", 0L)
    names_quad <- character(0)
    for (j in seq_len(d)) {
      quad[[length(quad) + 1L]] <- z[, j]^2
      names_quad <- c(names_quad, paste0("z_", colnames(z)[j], "^2"))
    }
    if (d > 1L) {
      for (a in seq_len(d - 1L)) {
        for (b in (a + 1L):d) {
          quad[[length(quad) + 1L]] <- z[, a] * z[, b]
          names_quad <- c(names_quad, paste0("z_", colnames(z)[a], ":z_", colnames(z)[b]))
        }
      }
    }
    Q <- do.call(cbind, quad)
    colnames(Q) <- names_quad
    X <- cbind(X, Q)
  }
  X
}

.theta_sur_ridge_fit <- function(X, y, noise_sd, lambda) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  noise_sd <- pmax(as.numeric(noise_sd), .Machine$double.eps)
  w <- 1 / noise_sd
  Xw <- X * w
  yw <- y * w
  P <- diag(ncol(X))
  P[1L, 1L] <- 0
  A <- crossprod(Xw) + as.numeric(lambda) * P
  b <- crossprod(Xw, yw)
  coef <- tryCatch(
    as.numeric(solve(A, b)),
    error = function(e) {
      sv <- svd(A)
      tol <- max(dim(A)) * .Machine$double.eps * max(sv$d)
      keep <- sv$d > tol
      if (!any(keep)) return(rep(0, ncol(A)))
      as.numeric(sv$v[, keep, drop = FALSE] %*% ((t(sv$u[, keep, drop = FALSE]) %*% b) / sv$d[keep]))
    }
  )
  names(coef) <- colnames(X)
  fitted <- as.numeric(X %*% coef)
  list(coef = coef, fitted = fitted)
}

.theta_sur_loo_score <- function(X, y, noise_sd, lambda) {
  n <- nrow(X)
  pred <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    train <- setdiff(seq_len(n), i)
    fit <- .theta_sur_ridge_fit(X[train, , drop = FALSE], y[train], noise_sd[train], lambda)
    pred[i] <- as.numeric(X[i, , drop = FALSE] %*% fit$coef)
  }
  err <- y - pred
  noise_sd <- pmax(noise_sd, .Machine$double.eps)
  data.frame(
    loo_nlp = mean(0.5 * (err / noise_sd)^2 + log(noise_sd)),
    loo_rmse = sqrt(mean(err^2)),
    loo_mae = mean(abs(err)),
    loo_weighted_rmse = sqrt(mean((err / noise_sd)^2)),
    loo_max_abs_error = max(abs(err)),
    check.names = FALSE
  )
}

.theta_sur_point_training <- function(correction_set) {
  if (!inherits(correction_set, "theta_smc_correction_set")) {
    stop("correction_set must inherit from 'theta_smc_correction_set'.")
  }
  model <- normalize_population_model(correction_set$factor_set$population_model)
  obs <- correction_set$observations
  obs$calibration_usable <- is.finite(obs$delta_total) &
    is.finite(obs$mcse_total) &
    is.finite(obs$logposterior_smc) &
    is.finite(obs$logposterior_tilde) &
    obs$degenerate_subjects == 0

  by_point <- split(obs, obs$point_id)
  rows <- list()
  rejected <- list()
  for (point_name in names(by_point)) {
    x_all <- by_point[[point_name]]
    x <- x_all[x_all$calibration_usable, , drop = FALSE]
    if (!nrow(x)) {
      rejected[[length(rejected) + 1L]] <- data.frame(
        point_id = x_all$point_id[1L],
        label = x_all$label[1L],
        role = x_all$role[1L],
        source = x_all$source[1L],
        reason = "not usable for smooth calibration",
        degenerate_subjects_max = max(x_all$degenerate_subjects, na.rm = TRUE),
        check.names = FALSE
      )
      next
    }
    delta <- as.numeric(x$delta_total)
    mcse <- pmax(as.numeric(x$mcse_total), .Machine$double.eps)
    n_rep <- length(delta)
    mcse_mean_sd <- sqrt(sum(mcse^2)) / n_rep
    replicate_mean_sd <- if (n_rep > 1L) stats::sd(delta) / sqrt(n_rep) else NA_real_
    noise_sd <- max(mcse_mean_sd, replicate_mean_sd %||% 0, .Machine$double.eps, na.rm = TRUE)
    theta_row <- as.numeric(x[1L, model$hyper_names])
    names(theta_row) <- model$hyper_names
    rows[[length(rows) + 1L]] <- data.frame(
      point_id = x$point_id[1L],
      label = x$label[1L],
      role = x$role[1L],
      source = x$source[1L],
      matrix(theta_row, nrow = 1L, dimnames = list(NULL, model$hyper_names)),
      delta = mean(delta),
      noise_sd = noise_sd,
      n_replicates = n_rep,
      replicate_sd = if (n_rep > 1L) stats::sd(delta) else NA_real_,
      mcse_mean = mean(mcse),
      relative_logposterior_smc = mean(x$relative_logposterior_smc),
      relative_logposterior_tilde = mean(x$relative_logposterior_tilde),
      check.names = FALSE
    )
  }
  training <- if (length(rows)) do.call(rbind, rows) else data.frame()
  rejected <- if (length(rejected)) do.call(rbind, rejected) else data.frame()
  rownames(training) <- NULL
  rownames(rejected) <- NULL
  list(training = training, rejected = rejected)
}

fit_theta_correction_surrogate <- function(correction_set,
                                           support_theta,
                                           support_w = NULL,
                                           degrees = 0:2,
                                           lambda_grid = c(0, 10^seq(-4, 4, length.out = 17))) {
  if (!inherits(correction_set, "theta_smc_correction_set")) {
    stop("correction_set must inherit from 'theta_smc_correction_set'.")
  }
  model <- normalize_population_model(correction_set$factor_set$population_model)
  split <- .theta_sur_point_training(correction_set)
  training <- split$training
  if (nrow(training) < 2L) {
    stop("At least two usable correction points are required to fit a surrogate.")
  }
  theta <- .theta_sur_as_matrix(training[, model$hyper_names, drop = FALSE], model)
  support_theta <- .theta_sur_as_matrix(support_theta, model)
  support_w <- .theta_sur_normalize_weights(support_w, nrow(support_theta))
  whitener <- .theta_sur_support_whitener(support_theta, support_w, model)
  z <- .theta_sur_to_z(theta, whitener, model)
  y <- as.numeric(training$delta)
  noise_sd <- pmax(as.numeric(training$noise_sd), .Machine$double.eps)

  candidates <- list()
  ptr <- 1L
  for (degree in as.integer(degrees)) {
    X <- .theta_sur_design_matrix(z, degree)
    for (lambda in as.numeric(lambda_grid)) {
      score <- .theta_sur_loo_score(X, y, noise_sd, lambda)
      fit <- .theta_sur_ridge_fit(X, y, noise_sd, lambda)
      candidates[[ptr]] <- cbind(
        data.frame(
          candidate_id = ptr,
          degree = degree,
          lambda = lambda,
          n_features = ncol(X),
          training_points = nrow(X),
          check.names = FALSE
        ),
        score
      )
      candidates[[ptr]]$training_rmse <- sqrt(mean((y - fit$fitted)^2))
      ptr <- ptr + 1L
    }
  }
  model_table <- do.call(rbind, candidates)
  ord <- order(model_table$loo_nlp, model_table$degree, model_table$lambda)
  model_table <- model_table[ord, , drop = FALSE]
  rownames(model_table) <- NULL
  selected <- model_table[1L, , drop = FALSE]

  X_selected <- .theta_sur_design_matrix(z, selected$degree[1L])
  fit <- .theta_sur_ridge_fit(X_selected, y, noise_sd, selected$lambda[1L])
  loo_pred <- rep(NA_real_, nrow(X_selected))
  for (i in seq_len(nrow(X_selected))) {
    train <- setdiff(seq_len(nrow(X_selected)), i)
    fit_i <- .theta_sur_ridge_fit(
      X_selected[train, , drop = FALSE],
      y[train],
      noise_sd[train],
      selected$lambda[1L]
    )
    loo_pred[i] <- as.numeric(X_selected[i, , drop = FALSE] %*% fit_i$coef)
  }
  validation <- data.frame(
    point_id = training$point_id,
    label = training$label,
    role = training$role,
    source = training$source,
    delta = y,
    fitted = fit$fitted,
    residual = y - fit$fitted,
    loo_pred = loo_pred,
    loo_residual = y - loo_pred,
    noise_sd = noise_sd,
    relative_logposterior_smc = training$relative_logposterior_smc,
    check.names = FALSE
  )

  structure(
    list(
      model = model,
      training = training,
      rejected = split$rejected,
      whitener = whitener,
      degree = as.integer(selected$degree[1L]),
      lambda = as.numeric(selected$lambda[1L]),
      coef = fit$coef,
      feature_names = names(fit$coef),
      model_table = model_table,
      validation = validation,
      settings = list(
        degrees = as.integer(degrees),
        lambda_grid = as.numeric(lambda_grid)
      )
    ),
    class = "theta_correction_surrogate"
  )
}

predict_theta_correction_surrogate <- function(surrogate, theta) {
  if (!inherits(surrogate, "theta_correction_surrogate")) {
    stop("surrogate must inherit from 'theta_correction_surrogate'.")
  }
  theta <- .theta_sur_as_matrix(theta, surrogate$model)
  z <- .theta_sur_to_z(theta, surrogate$whitener, surrogate$model)
  X <- .theta_sur_design_matrix(z, surrogate$degree)
  as.numeric(X[, surrogate$feature_names, drop = FALSE] %*% surrogate$coef)
}

select_theta_correction_refinement_points <- function(surrogate,
                                                      n_points = 4L,
                                                      next_point_id = NULL,
                                                      relevance_temperature = 10,
                                                      min_distance = 0.15) {
  if (!inherits(surrogate, "theta_correction_surrogate")) {
    stop("surrogate must inherit from 'theta_correction_surrogate'.")
  }
  model <- surrogate$model
  training <- surrogate$training
  validation <- surrogate$validation
  if (!nrow(training) || !nrow(validation)) {
    stop("surrogate must contain training and validation rows.")
  }

  theta <- .theta_sur_as_matrix(training[, model$hyper_names, drop = FALSE], model)
  z <- .theta_sur_to_z(theta, surrogate$whitener, model)
  rownames(z) <- as.character(training$point_id)
  validation$abs_loo_residual <- abs(validation$loo_residual)
  validation$relevance_weight <- exp(pmin(validation$relative_logposterior_smc, 0) / as.numeric(relevance_temperature))
  validation$refinement_score <- validation$abs_loo_residual * validation$relevance_weight
  validation <- validation[is.finite(validation$refinement_score) & validation$refinement_score > 0, , drop = FALSE]
  if (!nrow(validation)) {
    return(structure(data.frame(), class = c("theta_correction_design", "data.frame")))
  }

  center <- validation[validation$role == "center", , drop = FALSE]
  center_id <- if (nrow(center)) center$point_id[1L] else validation$point_id[which.max(validation$relative_logposterior_smc)]
  best_id <- validation$point_id[which.max(validation$relative_logposterior_smc)]
  validation <- validation[order(validation$refinement_score, decreasing = TRUE), , drop = FALSE]

  existing_z <- z
  selected_z <- NULL
  rows <- list()
  point_id <- as.integer(next_point_id %||% (max(training$point_id) + 1L))
  n_points <- as.integer(max(0L, n_points))
  if (n_points == 0L) {
    return(structure(data.frame(), class = c("theta_correction_design", "data.frame")))
  }

  add_edge <- function(seed_row, target_id, suffix) {
    seed_id <- as.character(seed_row$point_id)
    target_id <- as.character(target_id)
    if (!seed_id %in% rownames(z) || !target_id %in% rownames(z) || identical(seed_id, target_id)) {
      return(FALSE)
    }
    z_mid <- matrix((z[seed_id, ] + z[target_id, ]) / 2, nrow = 1L)
    if (!is.null(selected_z)) {
      if (min(sqrt(rowSums(sweep(selected_z, 2L, as.numeric(z_mid), "-")^2))) < min_distance) return(FALSE)
    }
    if (min(sqrt(rowSums(sweep(existing_z, 2L, as.numeric(z_mid), "-")^2))) < min_distance) return(FALSE)
    theta_mid <- .theta_sur_from_z(z_mid, surrogate$whitener, model)
    label <- paste0("refine_", seed_row$label, "_to_", suffix)
    rows[[length(rows) + 1L]] <<- data.frame(
      point_id = point_id + length(rows),
      source = "active_refinement",
      role = "surrogate_residual",
      label = label,
      matrix(theta_mid[1L, ], nrow = 1L, dimnames = list(NULL, model$hyper_names)),
      parent_point_id = seed_row$point_id,
      target_point_id = as.integer(target_id),
      parent_loo_residual = seed_row$loo_residual,
      parent_relative_logposterior_smc = seed_row$relative_logposterior_smc,
      parent_refinement_score = seed_row$refinement_score,
      check.names = FALSE
    )
    selected_z <<- if (is.null(selected_z)) z_mid else rbind(selected_z, z_mid)
    TRUE
  }

  for (i in seq_len(nrow(validation))) {
    seed <- validation[i, , drop = FALSE]
    add_edge(seed, center_id, "center")
    if (length(rows) >= n_points) break
    add_edge(seed, best_id, "best")
    if (length(rows) >= n_points) break
  }

  out <- if (length(rows)) do.call(rbind, rows) else data.frame()
  if (nrow(out) > n_points) out <- out[seq_len(n_points), , drop = FALSE]
  rownames(out) <- NULL
  structure(out, class = c("theta_correction_design", "data.frame"))
}
