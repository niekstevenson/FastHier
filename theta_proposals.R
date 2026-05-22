#!/usr/bin/env Rscript
# ============================================================================
# Theta proposal objects
# - Continuous, normalized q0(theta) proposals fitted from scored candidates
# - Defensive Student-t mixtures with an optional hyperprior floor
# ============================================================================

if (!exists("%||%", mode = "function") ||
    !exists("logsumexp", mode = "function") ||
    !exists(".rowLogSumExp", mode = "function") ||
    !exists("ESS", mode = "function") ||
    !exists("weighted_cov", mode = "function") ||
    !exists("regularize_cov", mode = "function") ||
    !exists("dmvt_chol_log", mode = "function")) {
  source("smc_core.R")
}
if (!exists("normalize_population_model", mode = "function") ||
    !exists("population_model_log_hyperprior", mode = "function") ||
    !exists("population_model_sample_hyper", mode = "function")) {
  source("population_models.R")
}

.theta_proposal_align <- function(theta, population_model) {
  model <- normalize_population_model(population_model)
  .as_hyper_matrix(theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
}

.theta_proposal_normalize_log_weights <- function(log_weight) {
  log_weight <- as.numeric(log_weight)
  ok <- is.finite(log_weight)
  if (!any(ok)) {
    stop("No finite theta candidate scores.")
  }
  w <- rep(0, length(log_weight))
  w[ok] <- exp(log_weight[ok] - logsumexp(log_weight[ok]))
  w / sum(w)
}

.theta_proposal_rmvt_chol <- function(n, mu, L, df) {
  n <- as.integer(n)
  d <- length(mu)
  z <- matrix(stats::rnorm(n * d), nrow = n, ncol = d)
  s <- sqrt(stats::rchisq(n, df = df) / df)
  out <- sweep((z %*% L) / s, 2L, mu, "+")
  colnames(out) <- names(mu)
  out
}

theta_proposal_weighted_quantile <- function(x, w, probs) {
  x <- as.numeric(x)
  w <- pmax(as.numeric(w), 0)
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) {
    return(stats::quantile(x, probs = probs, names = FALSE, type = 8))
  }
  x <- x[ok]
  w <- w[ok]
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  unique_x <- unique(x)
  if (length(unique_x) < length(x)) {
    w <- as.numeric(rowsum(w, group = match(x, unique_x), reorder = FALSE))
    x <- unique_x
  }
  w <- w / sum(w)
  stats::approx(cumsum(w), x, xout = probs, method = "constant", f = 1, rule = 2, ties = "ordered")$y
}

.theta_proposal_global_geometry <- function(theta, w) {
  d <- ncol(theta)
  mu <- colSums(theta * w)
  S <- tryCatch(weighted_cov(theta, w), error = function(e) NULL)
  if (is.null(S) || any(!is.finite(S))) {
    S <- stats::cov(theta)
  }
  if (is.null(S) || any(!is.finite(S))) {
    S <- diag(d)
  }
  S <- regularize_cov(S, min_eig = 1e-8, cond_cap = 1e8)
  L <- tryCatch(chol(S), error = function(e) chol(diag(d)))
  Z <- t(backsolve(L, t(sweep(theta, 2L, mu, "-")), transpose = TRUE))
  Z[!is.finite(Z)] <- 0
  list(mean = mu, cov = S, chol = L, z = Z)
}

.theta_proposal_assign_clusters <- function(Z, centers) {
  Z <- as.matrix(Z)
  centers <- as.matrix(centers)
  dist <- matrix(0, nrow = nrow(Z), ncol = nrow(centers))
  for (k in seq_len(nrow(centers))) {
    diff <- sweep(Z, 2L, centers[k, ], "-")
    dist[, k] <- rowSums(diff * diff)
  }
  max.col(-dist, ties.method = "first")
}

.theta_proposal_cluster <- function(theta,
                                    w,
                                    max_components,
                                    cluster_sample_size,
                                    min_component_weight,
                                    seed = NULL) {
  n <- nrow(theta)
  d <- ncol(theta)
  ess <- ESS(w)
  K <- as.integer(min(max_components, n, max(1L, floor(ess / max(4L * d, 20L)))))
  if (K <= 1L) {
    return(rep.int(1L, n))
  }

  geometry <- .theta_proposal_global_geometry(theta, w)
  if (!is.null(seed)) set.seed(as.integer(seed))
  sample_n <- as.integer(max(K, min(cluster_sample_size, n)))
  idx <- sample.int(n, sample_n, replace = TRUE, prob = w)
  km <- tryCatch(
    stats::kmeans(geometry$z[idx, , drop = FALSE], centers = K, nstart = 8L, iter.max = 80L),
    error = function(e) NULL
  )
  if (is.null(km)) {
    return(rep.int(1L, n))
  }

  cluster <- .theta_proposal_assign_clusters(geometry$z, km$centers)
  mass <- as.numeric(tapply(w, cluster, sum))
  keep <- which(mass >= min_component_weight)
  if (!length(keep)) {
    return(rep.int(1L, n))
  }
  if (length(keep) < K) {
    centers <- km$centers[keep, , drop = FALSE]
    cluster <- .theta_proposal_assign_clusters(geometry$z, centers)
  }
  cluster
}

.theta_proposal_component_from_cluster <- function(theta,
                                                  w,
                                                  global_cov,
                                                  df,
                                                  scale,
                                                  role,
                                                  min_eig,
                                                  cond_cap) {
  d <- ncol(theta)
  sw <- sum(w)
  w <- w / sw
  mu <- colSums(theta * w)
  S <- if (nrow(theta) > d + 1L && ESS(w) > d + 1L) {
    tryCatch(weighted_cov(theta, w), error = function(e) global_cov)
  } else {
    global_cov
  }
  if (is.null(S) || any(!is.finite(S))) S <- global_cov
  S <- regularize_cov(scale * S, min_eig = min_eig, cond_cap = cond_cap)
  L <- tryCatch(chol(S), error = function(e) chol(regularize_cov(S, min_eig = min_eig, cond_cap = cond_cap)))
  names(mu) <- colnames(theta)
  list(
    family = "student_t",
    role = role,
    mean = mu,
    cov = S,
    chol = L,
    df = as.numeric(df),
    base_mass = as.numeric(sw)
  )
}

fit_theta_q0_proposal <- function(theta,
                                  log_weight,
                                  population_model,
                                  max_components = 6L,
                                  core_weight = 0.85,
                                  tail_weight = 0.10,
                                  prior_weight = 0.05,
                                  df = 7,
                                  tail_df = 3,
                                  core_scale = 1.25,
                                  tail_scale = 9,
                                  cluster_sample_size = 2000L,
                                  min_component_weight = 0.02,
                                  min_eig = 1e-8,
                                  cond_cap = 1e8,
  seed = NULL,
  label = "theta_q0") {
  model <- normalize_population_model(population_model)
  candidate_source <- attr(theta, "source")
  theta <- .theta_proposal_align(theta, model)
  log_weight <- as.numeric(log_weight)
  if (length(log_weight) != nrow(theta)) {
    stop("log_weight must have one value per theta candidate.")
  }

  ok <- apply(theta, 1L, function(x) all(is.finite(x))) & is.finite(log_weight)
  if (!any(ok)) {
    stop("No finite theta candidates with finite scores.")
  }
  theta <- theta[ok, , drop = FALSE]
  source <- candidate_source %||% rep.int("candidate", length(ok))
  source <- source[ok]
  w <- .theta_proposal_normalize_log_weights(log_weight[ok])

  global <- .theta_proposal_global_geometry(theta, w)
  cluster <- .theta_proposal_cluster(
    theta = theta,
    w = w,
    max_components = as.integer(max_components),
    cluster_sample_size = as.integer(cluster_sample_size),
    min_component_weight = as.numeric(min_component_weight),
    seed = seed
  )

  cluster_ids <- sort(unique(cluster))
  base_components <- lapply(cluster_ids, function(k) {
    idx <- which(cluster == k)
    .theta_proposal_component_from_cluster(
      theta = theta[idx, , drop = FALSE],
      w = w[idx],
      global_cov = global$cov,
      df = df,
      scale = core_scale,
      role = "core",
      min_eig = min_eig,
      cond_cap = cond_cap
    )
  })
  base_mass <- vapply(base_components, `[[`, numeric(1), "base_mass")
  base_mass <- base_mass / sum(base_mass)

  components <- list()
  weights <- numeric(0)
  if (core_weight > 0) {
    components <- c(components, base_components)
    weights <- c(weights, as.numeric(core_weight) * base_mass)
  }
  if (tail_weight > 0) {
    tail_components <- lapply(base_components, function(component) {
      S <- regularize_cov(as.numeric(tail_scale) * component$cov, min_eig = min_eig, cond_cap = cond_cap)
      component$role <- "tail"
      component$cov <- S
      component$chol <- chol(S)
      component$df <- as.numeric(tail_df)
      component
    })
    components <- c(components, tail_components)
    weights <- c(weights, as.numeric(tail_weight) * base_mass)
  }

  prior_weight <- max(as.numeric(prior_weight), 0)
  total <- sum(weights) + prior_weight
  if (!is.finite(total) || total <= 0) {
    stop("Theta proposal mixture weights must have positive mass.")
  }
  weights <- weights / total
  prior_weight <- prior_weight / total
  components <- Map(function(component, weight) {
    component$weight <- as.numeric(weight)
    component
  }, components, weights)

  structure(
    list(
      name = as.character(label),
      population_model = model,
      components = components,
      prior_weight = as.numeric(prior_weight),
      hyper_names = model$hyper_names,
      hyper_dim = model$hyper_dim,
      diagnostics = list(
        n_candidates = nrow(theta),
        candidate_ess_frac = ESS(w) / length(w),
        max_candidate_weight = max(w),
        n_components = length(components),
        n_base_components = length(base_components),
        core_weight = as.numeric(core_weight) / total,
        tail_weight = as.numeric(tail_weight) / total,
        prior_weight = as.numeric(prior_weight),
        source_mass = sort(tapply(w, source, sum), decreasing = TRUE)
      )
    ),
    class = "theta_q0_proposal"
  )
}

normalize_theta_proposal <- function(proposal, population_model = NULL) {
  if (!inherits(proposal, "theta_q0_proposal")) {
    stop("proposal must inherit from 'theta_q0_proposal'.")
  }
  model <- if (is.null(population_model)) {
    normalize_population_model(proposal$population_model)
  } else {
    normalize_population_model(population_model)
  }
  proposal$population_model <- model
  proposal$hyper_names <- model$hyper_names
  proposal$hyper_dim <- model$hyper_dim

  weights <- vapply(proposal$components, `[[`, numeric(1), "weight")
  total <- sum(weights) + as.numeric(proposal$prior_weight %||% 0)
  if (!is.finite(total) || total <= 0) {
    stop("Theta proposal has no positive mixture mass.")
  }
  proposal$components <- lapply(proposal$components, function(component) {
    component$weight <- as.numeric(component$weight) / total
    component$mean <- as.numeric(component$mean)
    names(component$mean) <- model$hyper_names
    component$cov <- regularize_cov(component$cov, min_eig = 1e-8, cond_cap = 1e8)
    dimnames(component$cov) <- list(model$hyper_names, model$hyper_names)
    component$chol <- chol(component$cov)
    component$df <- as.numeric(component$df)
    component
  })
  proposal$prior_weight <- as.numeric(proposal$prior_weight %||% 0) / total
  proposal
}

combine_theta_q0_proposals <- function(proposals,
                                       weights = NULL,
                                       population_model = NULL,
                                       label = "combined_theta_q0") {
  if (!is.list(proposals) || !length(proposals)) {
    stop("proposals must be a non-empty list.")
  }
  base_model <- population_model %||% proposals[[1L]]$population_model
  model <- normalize_population_model(base_model)
  proposals <- lapply(proposals, normalize_theta_proposal, population_model = model)

  weights <- pmax(as.numeric(weights %||% rep(1, length(proposals))), 0)
  if (length(weights) != length(proposals)) {
    stop("weights must match proposals.")
  }
  sw <- sum(weights)
  if (!is.finite(sw) || sw <= 0) {
    stop("weights must have positive finite mass.")
  }
  weights <- weights / sw

  components <- list()
  prior_weight <- 0
  for (k in seq_along(proposals)) {
    proposal <- proposals[[k]]
    if (length(proposal$components)) {
      scaled <- lapply(proposal$components, function(component) {
        component$weight <- as.numeric(component$weight) * weights[k]
        component$source_proposal <- proposal$name %||% paste0("proposal_", k)
        component
      })
      components <- c(components, scaled)
    }
    prior_weight <- prior_weight + weights[k] * as.numeric(proposal$prior_weight %||% 0)
  }

  structure(
    list(
      name = as.character(label),
      population_model = model,
      components = components,
      prior_weight = as.numeric(prior_weight),
      hyper_names = model$hyper_names,
      hyper_dim = model$hyper_dim,
      diagnostics = list(
        source_proposals = vapply(proposals, function(x) x$name %||% "theta_q0", character(1)),
        source_weights = weights,
        n_components = length(components),
        prior_weight = as.numeric(prior_weight)
      )
    ),
    class = "theta_q0_proposal"
  )
}

theta_proposal_log_density <- function(proposal, theta) {
  proposal <- normalize_theta_proposal(proposal)
  theta <- .theta_proposal_align(theta, proposal$population_model)
  G <- length(proposal$components) + as.integer(proposal$prior_weight > 0)
  logs <- matrix(-Inf, nrow = nrow(theta), ncol = G)
  col <- 0L
  if (length(proposal$components)) {
    for (component in proposal$components) {
      col <- col + 1L
      logs[, col] <- log(component$weight) +
        dmvt_chol_log(theta, component$mean, component$chol, df = component$df)
    }
  }
  if (proposal$prior_weight > 0) {
    col <- col + 1L
    logs[, col] <- log(proposal$prior_weight) +
      population_model_log_hyperprior(proposal$population_model, theta)
  }
  as.numeric(.rowLogSumExp(logs[, seq_len(col), drop = FALSE]))
}

theta_proposal_sample <- function(proposal, n, seed = NULL) {
  proposal <- normalize_theta_proposal(proposal)
  if (!is.null(seed)) set.seed(as.integer(seed))
  n <- as.integer(n)
  if (n < 1L) {
    return(matrix(numeric(0), nrow = 0L, ncol = proposal$hyper_dim, dimnames = list(NULL, proposal$hyper_names)))
  }

  component_weights <- vapply(proposal$components, `[[`, numeric(1), "weight")
  weights <- c(component_weights, prior = proposal$prior_weight)
  weights <- weights / sum(weights)
  idx <- sample.int(length(weights), n, replace = TRUE, prob = weights)
  out <- matrix(NA_real_, nrow = n, ncol = proposal$hyper_dim, dimnames = list(NULL, proposal$hyper_names))

  if (length(proposal$components)) {
    for (k in seq_along(proposal$components)) {
      take <- which(idx == k)
      if (!length(take)) next
      component <- proposal$components[[k]]
      out[take, ] <- .theta_proposal_rmvt_chol(length(take), component$mean, component$chol, component$df)
    }
  }

  prior_idx <- length(weights)
  take_prior <- which(idx == prior_idx)
  if (length(take_prior)) {
    out[take_prior, ] <- population_model_sample_hyper(proposal$population_model, length(take_prior))
  }
  colnames(out) <- proposal$hyper_names
  out
}

summarize_theta_proposal_reference <- function(proposal,
                                               reference_theta,
                                               n_draws = 5000L,
                                               seed = NULL,
                                               label = proposal$name) {
  proposal <- normalize_theta_proposal(proposal)
  reference_theta <- .theta_proposal_align(reference_theta, proposal$population_model)
  draws <- theta_proposal_sample(proposal, n = n_draws, seed = seed)
  logq_ref <- theta_proposal_log_density(proposal, reference_theta)
  rows <- lapply(seq_len(ncol(reference_theta)), function(j) {
    q_prop <- stats::quantile(draws[, j], probs = c(0.01, 0.05, 0.95, 0.99), names = FALSE, type = 8)
    q_ref <- stats::quantile(reference_theta[, j], probs = c(0.01, 0.05, 0.95, 0.99), names = FALSE, type = 8)
    data.frame(
      method = as.character(label),
      parameter = colnames(reference_theta)[j],
      reference_inside_q01_q99 = mean(reference_theta[, j] >= q_prop[1L] & reference_theta[, j] <= q_prop[4L]),
      reference_inside_q05_q95 = mean(reference_theta[, j] >= q_prop[2L] & reference_theta[, j] <= q_prop[3L]),
      lower_q01_miss = max(q_prop[1L] - q_ref[1L], 0),
      upper_q99_miss = max(q_ref[4L] - q_prop[4L], 0),
      width99_ratio = (q_prop[4L] - q_prop[1L]) / max(q_ref[4L] - q_ref[1L], .Machine$double.eps),
      q_wasserstein = mean(abs(
        stats::quantile(draws[, j], probs = seq(0.01, 0.99, length.out = 99L), names = FALSE, type = 8) -
          stats::quantile(reference_theta[, j], probs = seq(0.01, 0.99, length.out = 99L), names = FALSE, type = 8)
      )),
      reference_logq_min = min(logq_ref),
      reference_logq_q01 = as.numeric(stats::quantile(logq_ref, 0.01, names = FALSE, type = 8)),
      reference_logq_median = median(logq_ref),
      check.names = FALSE
    )
  })
  do.call(rbind, rows)
}
