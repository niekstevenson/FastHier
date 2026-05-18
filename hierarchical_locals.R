#!/usr/bin/env Rscript
# ============================================================================
# Hierarchical local-reference workflow
# - Broad pilot locals
# - Pilot population fit
# - One population-informed full local pass with planned anchor components
# - Outer SMC with stratified audit and targeted factor repair
# ============================================================================

if (!exists("%||%", mode = "function") ||
    !exists("weighted_cov", mode = "function") ||
    !exists("regularize_cov", mode = "function")) {
  source("smc_core.R")
}
if (!exists("normalize_reference_prior", mode = "function") ||
    !exists("make_broad_reference_prior", mode = "function") ||
    !exists("make_reference_prior_gaussian", mode = "function") ||
    !exists("make_reference_prior_gaussian_mixture", mode = "function") ||
    !exists("combine_reference_priors", mode = "function") ||
    !exists("reference_prior_logpdf", mode = "function")) {
  source("reference_priors.R")
}
if (!exists("enhanced_smc_elite", mode = "function")) {
  source("SMC_super_fast.R")
}
if (!exists("normalize_population_model", mode = "function") ||
    !exists("population_model_prepare_theta", mode = "function") ||
    !exists("population_model_reference_components_from_theta", mode = "function")) {
  source("population_models.R")
}
if (!exists("build_population_factor_set", mode = "function") ||
    !exists("outer_population_smc", mode = "function")) {
  source("outer_population_smc.R")
}

suppressPackageStartupMessages({
  library(parallel)
})

normalize_particle_weights <- function(w) {
  w <- pmax(as.numeric(w), 0)
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) {
    rep(1 / length(w), length(w))
  } else {
    w / sw
  }
}

weighted_particle_moments <- function(Theta, w) {
  Theta <- as.matrix(Theta)
  w <- normalize_particle_weights(w)
  mu <- colSums(Theta * w)
  Sigma <- regularize_cov(weighted_cov(Theta, w), min_eig = 1e-8, cond_cap = 1e8)
  names(mu) <- colnames(Theta)
  dimnames(Sigma) <- list(colnames(Theta), colnames(Theta))
  list(mu = mu, Sigma = Sigma)
}

summarize_local_reference_fit <- function(local_fit) {
  if (is.null(local_fit$Theta) || is.null(local_fit$w)) {
    stop("local_fit must contain Theta and w.")
  }
  moments <- weighted_particle_moments(local_fit$Theta, local_fit$w)
  list(
    local_id = local_fit$local_id %||% NA_integer_,
    mu = moments$mu,
    Sigma = moments$Sigma,
    log_evidence = as.numeric(local_fit$log_evidence %||% NA_real_),
    mcse_log_evidence = as.numeric(local_fit$mcse_logZ %||% NA_real_)
  )
}

aggregate_pilot_local_moments <- function(pilot_fits) {
  if (!length(pilot_fits)) stop("pilot_fits must be non-empty.")
  summaries <- lapply(pilot_fits, summarize_local_reference_fit)
  mu_list <- lapply(summaries, `[[`, "mu")
  Sigma_list <- lapply(summaries, `[[`, "Sigma")
  mu_bar <- Reduce(`+`, mu_list) / length(mu_list)
  Sigma_bar <- Reduce(
    `+`,
    Map(
      function(mu_i, Sigma_i) {
        dm <- as.numeric(mu_i - mu_bar)
        Sigma_i + tcrossprod(dm)
      },
      mu_list,
      Sigma_list
    )
  ) / length(Sigma_list)
  Sigma_bar <- regularize_cov(Sigma_bar, min_eig = 1e-8, cond_cap = 1e8)
  names(mu_bar) <- names(mu_list[[1L]])
  dimnames(Sigma_bar) <- list(names(mu_bar), names(mu_bar))
  list(mu = mu_bar, Sigma = Sigma_bar, local_summaries = summaries)
}

default_local_feature_vector <- function(data_i) {
  if (is.data.frame(data_i)) {
    n_obs <- nrow(data_i)
    numeric_part <- data_i[, vapply(data_i, is.numeric, logical(1)), drop = FALSE]
    vals <- unlist(numeric_part, use.names = FALSE)
  } else if (is.matrix(data_i)) {
    n_obs <- nrow(data_i)
    vals <- as.numeric(data_i[is.finite(data_i)])
  } else {
    flat <- unlist(data_i, recursive = TRUE, use.names = FALSE)
    n_obs <- length(flat)
    vals <- suppressWarnings(as.numeric(flat))
    vals <- vals[is.finite(vals)]
  }

  vals <- vals[is.finite(vals)]
  if (!length(vals)) vals <- 0
  value_sd <- stats::sd(vals)
  if (!is.finite(value_sd)) value_sd <- 0

  c(
    n_obs = as.numeric(max(1L, n_obs)),
    value_mean = mean(vals),
    value_sd = value_sd,
    value_q90 = as.numeric(stats::quantile(vals, probs = 0.90, names = FALSE, na.rm = TRUE))
  )
}

compute_local_feature_matrix <- function(data_list, feature_fn = default_local_feature_vector) {
  features <- lapply(data_list, feature_fn)
  lens <- vapply(features, length, integer(1))
  if (length(unique(lens)) != 1L) {
    stop("feature_fn must return feature vectors with a constant length across locals.")
  }
  names_ref <- names(features[[1L]]) %||% paste0("feature_", seq_len(lens[1L]))
  out <- do.call(rbind, lapply(features, function(x) {
    x <- as.numeric(x)
    names(x) <- names_ref
    x
  }))
  colnames(out) <- names_ref
  out[!is.finite(out)] <- 0
  out
}

make_pilot_strata <- function(data_list = NULL,
                              feature_matrix = NULL,
                              n_strata = NULL,
                              feature_fn = default_local_feature_vector,
                              seed = 123L) {
  feature_matrix <- if (is.null(feature_matrix)) {
    compute_local_feature_matrix(data_list, feature_fn = feature_fn)
  } else {
    as.matrix(feature_matrix)
  }
  L <- nrow(feature_matrix)
  if (L <= 1L) return(rep.int(1L, L))

  n_strata <- as.integer(n_strata %||% min(max(2L, floor(sqrt(L))), L))
  n_strata <- max(1L, min(n_strata, L))
  X <- scale(feature_matrix)
  X[!is.finite(X)] <- 0

  if (n_strata == 1L || nrow(unique(X)) == 1L) {
    return(rep.int(1L, L))
  }

  set.seed(as.integer(seed))
  km <- stats::kmeans(
    X,
    centers = n_strata,
    nstart = min(10L, n_strata),
    iter.max = 50L
  )
  as.integer(km$cluster)
}

.allocate_stratified_counts <- function(strata, subset_size) {
  strata <- as.integer(strata)
  tab <- table(strata)
  counts <- as.integer(tab)
  H <- length(counts)
  subset_size <- min(as.integer(subset_size), sum(counts))
  if (subset_size <= 0L) stop("subset_size must be positive.")

  raw <- subset_size * counts / sum(counts)
  alloc <- floor(raw)
  min_alloc <- if (subset_size >= H) as.integer(counts > 0) else rep.int(0L, H)
  alloc <- pmax(alloc, min_alloc)
  alloc <- pmin(alloc, counts)

  while (sum(alloc) > subset_size) {
    reducible <- which(alloc > min_alloc)
    if (!length(reducible)) break
    idx <- reducible[which.min(raw[reducible] - alloc[reducible])]
    alloc[idx] <- alloc[idx] - 1L
  }
  while (sum(alloc) < subset_size) {
    expandable <- which(alloc < counts)
    if (!length(expandable)) break
    idx <- expandable[which.max(raw[expandable] - alloc[expandable])]
    alloc[idx] <- alloc[idx] + 1L
  }

  stats::setNames(alloc, names(tab))
}

select_stratified_pilot_subset <- function(data_list,
                                           pilot_size = 20L,
                                           strata = NULL,
                                           feature_fn = default_local_feature_vector,
                                           n_strata = NULL,
                                           seed = 123L) {
  L <- length(data_list)
  pilot_size <- min(as.integer(pilot_size), L)
  if (pilot_size <= 0L) stop("pilot_size must be positive.")

  feature_matrix <- compute_local_feature_matrix(data_list, feature_fn = feature_fn)
  strata <- strata %||% make_pilot_strata(
    feature_matrix = feature_matrix,
    n_strata = n_strata,
    seed = seed
  )
  strata <- as.integer(strata)
  alloc <- .allocate_stratified_counts(strata, subset_size = pilot_size)

  set.seed(as.integer(seed))
  groups <- split(seq_len(L), strata)
  selected <- unlist(
    Map(
      function(idx, n_take) {
        if (n_take <= 0L) return(integer(0))
        sample(idx, size = n_take, replace = FALSE)
      },
      groups,
      alloc[names(groups)]
    ),
    use.names = FALSE
  )
  selected <- sort(as.integer(selected))

  list(
    indices = selected,
    strata = strata,
    feature_matrix = feature_matrix,
    allocation = alloc
  )
}

run_reference_local_smc <- function(data_list,
                                    loglik_fn,
                                    reference_prior,
                                    indices = seq_along(data_list),
                                    M = 2000L,
                                    n_jobs = 1L,
                                    local_n_cores = 1L,
                                    base_seed = 123L,
                                    verbose = TRUE,
                                    ...) {
  reference_prior <- normalize_reference_prior(reference_prior = reference_prior)
  extra_args <- list(...)
  blocked <- intersect(names(extra_args), c("data", "loglik_fn", "reference_prior", "M", "n_cores", "seed", "verbose"))
  if (length(blocked)) {
    stop("Pass ", paste(blocked, collapse = ", "), " via the dedicated run_reference_local_smc arguments.")
  }
  indices <- sort(unique(as.integer(indices)))
  fits <- parallel::mclapply(
    indices,
    function(i) {
      fit_args <- modifyList(
        list(
          data = data_list[[i]],
          loglik_fn = loglik_fn,
          reference_prior = reference_prior,
          M = as.integer(M),
          n_cores = as.integer(local_n_cores),
          seed = as.integer(base_seed + i - 1L),
          verbose = verbose
        ),
        extra_args
      )
      fit <- do.call(enhanced_smc_elite, fit_args)
      fit$local_id <- as.integer(i)
      fit$reference_prior_label <- reference_prior$label
      fit
    },
    mc.cores = as.integer(max(1L, n_jobs))
  )
  names(fits) <- as.character(indices)
  fits
}

fit_population_model_from_local_objects <- function(local_objects,
                                                    population_model,
                                                    outer_control = list(),
                                                    n_cores = 1L,
                                                    seed = 123L,
                                                    verbose = TRUE) {
  population_model <- normalize_population_model(population_model)
  factor_set <- build_population_factor_set(local_objects, population_model)
  fit_args <- modifyList(
    list(
      factor_set = factor_set,
      N = 2000L,
      n_mcmc_moves = 3L,
      max_rounds = 80L,
      n_cores = as.integer(n_cores),
      seed = as.integer(seed),
      verbose = verbose
    ),
    outer_control
  )
  do.call(outer_population_smc, fit_args)
}

fit_pilot_population_model <- function(pilot_fits,
                                       reference_prior,
                                       population_model,
                                       outer_control = list(),
                                       n_cores = 1L,
                                       seed = 123L,
                                       verbose = TRUE) {
  pilot_objects <- build_local_reference_objects(pilot_fits, reference_prior = reference_prior)
  fit_population_model_from_local_objects(
    local_objects = pilot_objects,
    population_model = population_model,
    outer_control = outer_control,
    n_cores = n_cores,
    seed = seed,
    verbose = verbose
  )
}

.sample_weighted_indices <- function(w, size, seed = NULL) {
  w <- normalize_particle_weights(w)
  n <- length(w)
  size <- as.integer(max(1L, min(size, n)))
  if (!is.null(seed)) set.seed(as.integer(seed))
  sample.int(n, size = size, replace = FALSE, prob = w)
}

build_refined_reference_prior_from_population_fit <- function(population_fit,
                                                              population_model,
                                                              inflation = 1.5,
                                                              defensive_weight = 0.10,
                                                              broad_reference = NULL,
                                                              defensive_scale = 4,
                                                              support_size = 8L,
                                                              support_seed = 123L) {
  model <- normalize_population_model(population_model)
  theta <- .as_hyper_matrix(population_fit$theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  w <- normalize_particle_weights(population_fit$w)
  comps <- population_model_reference_components_from_theta(model, theta)
  idx <- .sample_weighted_indices(w, size = support_size, seed = support_seed)

  comp_means <- comps$component_means[idx]
  comp_covs <- lapply(comps$component_covs[idx], function(S) as.numeric(inflation) * S)
  comp_weights <- normalize_particle_weights(w[idx])
  core_prior <- make_reference_prior_gaussian_mixture(
    component_means = comp_means,
    component_covs = comp_covs,
    weights = comp_weights,
    param_names = model$alpha_names,
    label = "pilot_population_mixture"
  )

  defensive_prior <- if (is.null(broad_reference)) {
    inflate_reference_prior(core_prior, scale = defensive_scale, label = "pilot_population_defensive_component")
  } else {
    normalize_reference_prior(reference_prior = broad_reference)
  }

  if (!is.finite(defensive_weight) || defensive_weight <= 0) {
    return(core_prior)
  }
  combine_reference_priors(
    priors = list(core_prior, defensive_prior),
    weights = c(1 - defensive_weight, defensive_weight),
    label = "pilot_population_defensive_mixture"
  )
}

build_local_reference_component <- function(local_fit,
                                           reference_prior,
                                           label = NULL) {
  if (is.null(local_fit$Theta) || is.null(local_fit$w)) {
    stop("local_fit must contain Theta and w.")
  }
  particles <- as.matrix(local_fit$Theta)
  reference_prior <- normalize_reference_prior(reference_prior = reference_prior)
  list(
    label = as.character(label %||% reference_prior$label %||% ""),
    particles = particles,
    weights = normalize_particle_weights(local_fit$w),
    proposal_type = "posterior_reference",
    log_likelihood = as.numeric(local_fit$loglik %||% rep(NA_real_, nrow(particles))),
    log_evidence = as.numeric(local_fit$log_evidence %||% NA_real_),
    mcse_log_evidence = as.numeric(local_fit$mcse_logZ %||% NA_real_),
    reference_prior = reference_prior,
    log_reference_density = reference_prior_logpdf(reference_prior, particles),
    diagnostics = list(
      rounds = as.integer(local_fit$meta$rounds %||% NA_integer_),
      final_lambda = as.numeric(local_fit$final_lambda %||% NA_real_)
    )
  )
}

build_local_reference_object <- function(local_fit,
                                         reference_prior,
                                         label = NULL) {
  component <- build_local_reference_component(
    local_fit = local_fit,
    reference_prior = reference_prior,
    label = label
  )
  structure(
    list(
      local_id = as.integer(local_fit$local_id %||% NA_integer_),
      particles = component$particles,
      weights = component$weights,
      proposal_type = component$proposal_type,
      log_likelihood = component$log_likelihood,
      log_evidence = component$log_evidence,
      mcse_log_evidence = component$mcse_log_evidence,
      reference_prior = component$reference_prior,
      log_reference_density = component$log_reference_density,
      diagnostics = component$diagnostics,
      components = list(component),
      mixture_weights = 1
    ),
    class = "reference_local_object"
  )
}

build_local_reference_objects <- function(local_fits, reference_prior) {
  objs <- lapply(local_fits, build_local_reference_object, reference_prior = reference_prior)
  names(objs) <- names(local_fits)
  objs
}

make_population_reference_prior_from_theta <- function(population_model,
                                                       theta,
                                                       label = "population_reference") {
  model <- normalize_population_model(population_model)
  theta <- .as_hyper_matrix(theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  if (nrow(theta) != 1L) {
    stop("make_population_reference_prior_from_theta expects one theta row.")
  }
  comps <- population_model_reference_components_from_theta(model, theta)
  make_reference_prior_gaussian(
    mu = comps$component_means[[1L]],
    Sigma = comps$component_covs[[1L]],
    param_names = model$alpha_names,
    label = label
  )
}

build_qmc_population_anchor_local_objects <- function(data_list,
                                                      loglik_fn,
                                                      population_model,
                                                      theta,
                                                      local_ids,
                                                      qmc_size = 8192L,
                                                      qmc_randomizations = 2L,
                                                      seed = 123L,
                                                      n_jobs = 1L,
                                                      local_n_cores = 1L,
                                                      label = "qmc_population_anchor",
                                                      diagnostics = list()) {
  model <- normalize_population_model(population_model)
  theta <- .as_hyper_matrix(theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  if (nrow(theta) != 1L) stop("QMC anchor components require one theta row.")
  theta_prepared <- population_model_prepare_theta(model, theta)
  if (!identical(theta_prepared$family, "gaussian") ||
      !identical(theta_prepared$quadratic_kind, "diag")) {
    stop("QMC anchor components require a diagonal Gaussian population model.")
  }

  local_ids <- sort(unique(as.integer(local_ids)))
  local_ids <- local_ids[local_ids >= 1L & local_ids <= length(data_list)]
  if (!length(local_ids)) stop("No valid locals selected for QMC anchor components.")

  qmc_size <- as.integer(max(16L, qmc_size))
  qmc_randomizations <- as.integer(max(1L, qmc_randomizations))
  reference_prior <- make_population_reference_prior_from_theta(
    population_model = model,
    theta = theta,
    label = label
  )
  sigma <- sqrt(1 / theta_prepared$quadratic_coef[1L, ])
  mu <- theta_prepared$mean[1L, ]

  alpha_reps <- lapply(seq_len(qmc_randomizations), function(r) {
    set.seed(as.integer(seed + r - 1L))
    U <- qrng::sobol(n = qmc_size, d = model$alpha_dim, randomize = TRUE)
    z <- qnorm(pmin(pmax(U, 1e-12), 1 - 1e-12))
    alpha <- sweep(z, 2L, sigma, "*")
    alpha <- sweep(alpha, 2L, mu, "+")
    colnames(alpha) <- model$alpha_names
    alpha
  })
  alpha_all <- do.call(rbind, alpha_reps)
  log_reference_density <- reference_prior_logpdf(reference_prior, alpha_all)

  objects <- parallel::mclapply(
    local_ids,
    function(local_id) {
      ll_reps <- lapply(alpha_reps, function(alpha) {
        ll <- ll_parallel(alpha, data_list[[local_id]], loglik_fn, n_cores = local_n_cores)
        ll[!is.finite(ll)] <- -Inf
        ll
      })
      logm_rep <- vapply(ll_reps, function(ll) logsumexp(ll) - log(length(ll)), numeric(1))
      ll_all <- unlist(ll_reps, use.names = FALSE)
      log_evidence <- logsumexp(ll_all) - log(length(ll_all))
      mcse <- if (length(logm_rep) > 1L) stats::sd(logm_rep) / sqrt(length(logm_rep)) else NA_real_
      component_diagnostics <- modifyList(
        list(
          method = "qmc_population_anchor",
          qmc_size = qmc_size,
          qmc_randomizations = qmc_randomizations,
          theta = stats::setNames(as.numeric(theta[1L, ]), model$hyper_names)
        ),
        diagnostics
      )
      component <- list(
        label = as.character(label),
        particles = alpha_all,
        weights = rep(1 / nrow(alpha_all), nrow(alpha_all)),
        proposal_type = "prior_reference",
        log_likelihood = as.numeric(ll_all),
        log_evidence = as.numeric(log_evidence),
        mcse_log_evidence = as.numeric(mcse),
        reference_prior = reference_prior,
        log_reference_density = log_reference_density,
        diagnostics = component_diagnostics
      )
      structure(
        list(
          local_id = as.integer(local_id),
          components = list(component),
          mixture_weights = 1,
          diagnostics = list(mode = "qmc_population_anchor")
        ),
        class = "reference_local_object"
      )
    },
    mc.cores = as.integer(max(1L, min(n_jobs, length(local_ids))))
  )
  names(objects) <- as.character(local_ids)
  list(
    local_objects = objects,
    reference_prior = reference_prior,
    theta = theta,
    local_ids = local_ids
  )
}

select_population_anchor_support <- function(population_fit,
                                             population_model,
                                             n_anchors = 4L,
                                             seed = 123L) {
  model <- normalize_population_model(population_model)
  theta <- .as_hyper_matrix(population_fit$theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  w <- normalize_particle_weights(population_fit$w)
  n_anchors <- as.integer(max(1L, min(n_anchors, nrow(theta))))
  idx <- .certification_theta_design(theta = theta, w = w, size = n_anchors, seed = seed)
  theta[idx, , drop = FALSE]
}

build_planned_anchor_local_objects <- function(data_list,
                                               loglik_fn,
                                               population_model,
                                               anchor_theta,
                                               local_ids = seq_along(data_list),
                                               qmc_size = 2048L,
                                               qmc_randomizations = 2L,
                                               seed = 123L,
                                               n_jobs = 1L,
                                               local_n_cores = 1L,
                                               label = "planned_anchor") {
  local_ids <- sort(unique(as.integer(local_ids)))
  local_ids <- local_ids[local_ids >= 1L & local_ids <= length(data_list)]
  if (!length(local_ids)) stop("No valid locals selected for planned anchors.")
  anchor_theta <- as.matrix(anchor_theta)
  if (!nrow(anchor_theta)) stop("anchor_theta must contain at least one row.")

  local_objects <- NULL
  anchors <- vector("list", nrow(anchor_theta))
  for (a in seq_len(nrow(anchor_theta))) {
    anchor <- build_qmc_population_anchor_local_objects(
      data_list = data_list,
      loglik_fn = loglik_fn,
      population_model = population_model,
      theta = anchor_theta[a, , drop = FALSE],
      local_ids = local_ids,
      qmc_size = qmc_size,
      qmc_randomizations = qmc_randomizations,
      seed = seed + 1000L * a,
      n_jobs = n_jobs,
      local_n_cores = local_n_cores,
      label = sprintf("%s_%d", label, a),
      diagnostics = list(
        role = "planned_anchor",
        anchor_rank = a
      )
    )
    local_objects <- if (is.null(local_objects)) {
      anchor$local_objects
    } else {
      add_local_reference_components_dmis(local_objects, anchor$local_objects)
    }
    anchors[[a]] <- anchor
  }

  list(
    local_objects = local_objects,
    anchors = anchors,
    theta = anchor_theta,
    local_ids = local_ids
  )
}

.default_reference_component_weights <- function(components) {
  proposal_types <- vapply(components, function(component) {
    as.character(component$proposal_type %||% "posterior_reference")
  }, character(1))
  posterior <- proposal_types == "posterior_reference"
  prior <- proposal_types == "prior_reference"
  if (any(posterior) && any(prior)) {
    # Posterior-reference SMC components carry local logZ error; calibrated QMC anchors should dominate the mixture density.
    posterior_total <- 0.01
    out <- numeric(length(components))
    out[posterior] <- posterior_total / sum(posterior)
    out[prior] <- (1 - posterior_total) / sum(prior)
    return(out)
  }
  rep(1 / length(components), length(components))
}

rebalance_local_reference_object_components <- function(local_object,
                                                        component_weights = NULL) {
  if (!is.list(local_object) || is.null(local_object$components)) {
    local_object <- validate_reference_local_object(local_object)
  } else {
    local_object$components <- lapply(local_object$components, validate_reference_local_component)
  }
  n_components <- length(local_object$components)
  weights <- normalize_particle_weights(component_weights %||% .default_reference_component_weights(local_object$components))
  if (length(weights) != n_components) {
    stop("component_weights must match the number of local components.")
  }
  local_object$mixture_weights <- weights
  local_object
}

.reference_local_component_count <- function(local_object) {
  if (.is_compressed_population_local_object(local_object)) {
    return(as.integer(local_object$diagnostics$n_source_components %||% 1L))
  }
  length(validate_reference_local_object(local_object)$components)
}

.reference_local_particle_count <- function(local_object) {
  if (.is_compressed_population_local_object(local_object)) {
    return(as.integer(nrow(local_object$factor_particles)))
  }
  components <- validate_reference_local_object(local_object)$components
  as.integer(sum(vapply(components, function(component) nrow(component$particles), integer(1))))
}

.reference_local_budget_summary <- function(local_objects) {
  data.frame(
    local = vapply(local_objects, function(x) as.integer(x$local_id %||% NA_integer_), integer(1)),
    components = vapply(local_objects, .reference_local_component_count, integer(1)),
    particles = vapply(local_objects, .reference_local_particle_count, integer(1)),
    check.names = FALSE
  )
}

.stratified_resample_n <- function(w, n, deterministic = TRUE) {
  n <- as.integer(n)
  if (n <= 0L) return(integer(0))
  w <- normalize_reference_local_weights(w)
  cw <- c(0, cumsum(w))
  u0 <- if (isTRUE(deterministic)) 0.5 / n else stats::runif(1L) / n
  u <- u0 + (0:(n - 1L)) / n
  as.integer(findInterval(u, cw, rightmost.closed = TRUE))
}

.local_particle_order <- function(particles) {
  particles <- as.matrix(particles)
  if (nrow(particles) <= 1L) return(seq_len(nrow(particles)))
  Z <- scale(particles)
  Z[!is.finite(Z)] <- 0
  order(rowSums(Z))
}

.reference_component_theta <- function(component, population_model) {
  diagnostics <- component$diagnostics %||% list()
  theta <- diagnostics$theta %||% diagnostics$anchor_theta
  if (is.null(theta)) return(NULL)
  theta <- as.numeric(theta)
  model <- normalize_population_model(population_model)
  if (length(theta) != model$hyper_dim || any(!is.finite(theta))) return(NULL)
  names(theta) <- model$hyper_names
  theta
}

.local_object_anchor_theta_design <- function(local_objects,
                                              population_model,
                                              extra_theta = NULL,
                                              max_rows = 32L) {
  model <- normalize_population_model(population_model)
  rows <- list()
  for (obj in local_objects) {
    if (.is_compressed_population_local_object(obj) || is.null(obj$components)) next
    components <- validate_reference_local_object(obj)$components
    for (component in components) {
      theta <- .reference_component_theta(component, model)
      if (!is.null(theta)) rows[[length(rows) + 1L]] <- theta
    }
  }
  if (!is.null(extra_theta)) {
    extra_theta <- .as_hyper_matrix(extra_theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
    for (i in seq_len(nrow(extra_theta))) {
      rows[[length(rows) + 1L]] <- extra_theta[i, ]
    }
  }
  if (!length(rows)) stop("No theta design is available for local factor compression.")
  theta <- unique(as.data.frame(do.call(rbind, rows), check.names = FALSE))
  theta <- as.matrix(theta)
  colnames(theta) <- model$hyper_names
  max_rows <- as.integer(max(1L, max_rows))
  if (nrow(theta) > max_rows) {
    theta <- theta[seq_len(max_rows), , drop = FALSE]
  }
  theta
}

.compress_population_local_factor <- function(factor,
                                              theta_design,
                                              max_particles,
                                              defensive_weight = 0.05) {
  stopifnot(inherits(factor, "population_local_factor"))
  n <- nrow(factor$particles)
  max_particles <- as.integer(max_particles)
  if (!is.finite(max_particles) || max_particles <= 0L || n <= max_particles) {
    return(structure(list(
      local_id = factor$local_id,
      factor_particles = factor$particles,
      factor_log_base = factor$log_base,
      factor_log_constant = factor$log_constant,
      diagnostics = list(
        mode = "compressed_population_factor",
        n_source_components = as.integer(length(unique(factor$component_id))),
        compressed_from = as.integer(n),
        compressed_to = as.integer(n),
        compression_draws = as.integer(n)
      )
    ), class = "compressed_population_local_object"))
  }

  theta_design <- .as_hyper_matrix(
    theta_design,
    hyper_names = factor$population_model$hyper_names,
    hyper_dim = factor$population_model$hyper_dim
  )
  ord <- .local_particle_order(factor$particles)
  alpha <- factor$particles[ord, , drop = FALSE]
  log_base <- as.numeric(factor$log_base[ord])

  theta_prepared <- population_model_prepare_theta(factor$population_model, theta_design)
  log_terms <- .population_log_alpha_given_theta_many(
    model = factor$population_model,
    alpha = alpha,
    theta_prepared = theta_prepared
  )
  log_terms <- sweep(log_terms, 2L, log_base, "+")
  log_norm <- apply(log_terms, 1L, logsumexp)
  log_influence <- sweep(log_terms, 1L, log_norm, "-")
  log_q_influence <- .rowLogSumExp(t(log_influence)) - log(nrow(log_influence))
  log_q_base <- log_base - logsumexp(log_base)

  defensive_weight <- pmin(pmax(as.numeric(defensive_weight), 0), 1)
  log_q <- rlogsumexp2(
    log1p(-defensive_weight) + log_q_influence,
    log(defensive_weight) + log_q_base
  )
  q <- exp(log_q - logsumexp(log_q))
  q[!is.finite(q)] <- 0
  q <- pmax(q, .Machine$double.xmin)
  q <- q / sum(q)

  idx <- .stratified_resample_n(q, max_particles, deterministic = TRUE)
  counts <- tabulate(idx, nbins = n)
  keep <- which(counts > 0L)
  compressed_log_base <- log_base[keep] + log(counts[keep]) - log(max_particles) - log(q[keep])

  structure(list(
    local_id = factor$local_id,
    factor_particles = alpha[keep, , drop = FALSE],
    factor_log_base = as.numeric(compressed_log_base),
    factor_log_constant = factor$log_constant,
    diagnostics = list(
      mode = "compressed_population_factor",
      n_source_components = as.integer(length(unique(factor$component_id))),
      compressed_from = as.integer(n),
      compressed_to = as.integer(length(keep)),
      compression_draws = as.integer(max_particles),
      compression_design_rows = as.integer(nrow(theta_design)),
      compression_defensive_weight = as.numeric(defensive_weight)
    )
  ), class = "compressed_population_local_object")
}

compress_reference_local_objects_to_factors <- function(local_objects,
                                                        population_model,
                                                        theta_design,
                                                        data_list = NULL,
                                                        loglik_fn = NULL,
                                                        local_ids = NULL,
                                                        existing = NULL,
                                                        max_particles = 5000L,
                                                        defensive_weight = 0.05,
                                                        local_n_cores = 1L) {
  out <- existing %||% vector("list", length(local_objects))
  names(out) <- names(local_objects)
  object_ids <- vapply(local_objects, function(x) as.integer(x$local_id %||% NA_integer_), integer(1))
  idx <- if (is.null(local_ids)) {
    seq_along(local_objects)
  } else {
    match(sort(unique(as.integer(local_ids))), object_ids)
  }
  idx <- idx[!is.na(idx)]
  local_names <- names(local_objects)
  for (i in idx) {
    factor <- build_population_local_factor(
      local_objects[[i]],
      population_model = population_model,
      data_i = .data_for_local_object(local_objects[[i]], i, data_list, local_names),
      loglik_fn = loglik_fn,
      local_n_cores = local_n_cores
    )
    out[[i]] <- .compress_population_local_factor(
      factor = factor,
      theta_design = theta_design,
      max_particles = max_particles,
      defensive_weight = defensive_weight
    )
  }
  out
}

budget_reference_local_object <- function(local_object,
                                          max_components = Inf) {
  local_object <- validate_reference_local_object(local_object)
  components <- local_object$components
  if (!length(components)) stop("local object must contain at least one component.")

  priority <- vapply(seq_along(components), function(j) {
    component <- components[[j]]
    if (identical(component$proposal_type, "posterior_reference")) return(Inf)
    diagnostics <- component$diagnostics %||% list()
    if (identical(as.character(diagnostics$role %||% ""), "planned_anchor")) return(1e6 - j)
    as.numeric(diagnostics$anchor_score %||% 0)
  }, numeric(1))

  if (!is.null(max_components) && is.finite(max_components) && length(components) > max_components) {
    keep <- order(priority, decreasing = TRUE)[seq_len(as.integer(max_components))]
    keep <- sort(keep)
    components <- components[keep]
  }

  local_object$components <- components
  local_object <- rebalance_local_reference_object_components(local_object)

  local_object$diagnostics <- modifyList(
    local_object$diagnostics %||% list(),
    list(mode = "component_budget", n_components = length(components))
  )
  rebalance_local_reference_object_components(local_object)
}

budget_reference_local_objects <- function(local_objects,
                                           max_components = Inf,
                                           local_ids = NULL) {
  out <- local_objects
  object_ids <- vapply(out, function(x) as.integer(x$local_id %||% NA_integer_), integer(1))
  idx <- if (is.null(local_ids)) {
    seq_along(out)
  } else {
    match(sort(unique(as.integer(local_ids))), object_ids)
  }
  idx <- idx[!is.na(idx)]
  for (i in idx) {
    out[[i]] <- budget_reference_local_object(
      out[[i]],
      max_components = max_components
    )
  }
  out
}

add_local_reference_components_dmis <- function(local_objects,
                                                added_objects,
                                                component_weights = NULL) {
  if (!length(added_objects)) return(local_objects)
  out <- local_objects
  object_names <- names(out)
  if (is.null(object_names)) object_names <- as.character(seq_along(out))

  for (nm in names(added_objects) %||% as.character(seq_along(added_objects))) {
    added <- added_objects[[nm]]
    added <- validate_reference_local_object(added)
    local_id <- as.integer(added$local_id %||% suppressWarnings(as.integer(nm)))
    target_name <- if (nm %in% object_names) {
      nm
    } else {
      hit <- which(vapply(out, function(x) as.integer(x$local_id %||% NA_integer_) == local_id, logical(1)))
      if (!length(hit)) stop("No existing local object for added local ", local_id, ".")
      object_names[hit[1L]]
    }

    current <- validate_reference_local_object(out[[target_name]])
    current$components <- c(current$components, added$components)
    current$diagnostics <- modifyList(
      current$diagnostics %||% list(),
      list(mode = "component_merge", n_components = length(current$components))
    )
    out[[target_name]] <- rebalance_local_reference_object_components(
      current,
      component_weights = component_weights
    )
  }
  out
}

.weighted_theta_moments <- function(theta, w) {
  theta <- as.matrix(theta)
  w <- normalize_particle_weights(w)
  center <- colSums(theta * w)
  S <- tryCatch(weighted_cov(theta, w), error = function(e) stats::cov(theta))
  if (is.null(S) || any(!is.finite(S))) {
    S <- diag(ncol(theta))
  }
  S <- regularize_cov(S, min_eig = 1e-8, cond_cap = 1e8)
  list(center = center, cov = S)
}

select_population_certification_design <- function(population_fit,
                                                   factor_set = NULL,
                                                   population_model = NULL,
                                                   size = 8L,
                                                   seed = NULL,
                                                   stress_pool_size = 0L,
                                                   stress_scale = 2.5,
                                                   stress_weight = 0.35,
                                                   stress_log_drop = 30,
                                                   n_cores = 1L) {
  if (is.null(population_fit$theta) || is.null(population_fit$w)) {
    stop("population_fit must contain theta and w.")
  }
  model <- normalize_population_model(population_model %||% population_fit$population_model %||% factor_set$population_model)
  theta <- .as_hyper_matrix(population_fit$theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  w <- normalize_particle_weights(population_fit$w)

  theta_all <- theta
  source <- rep.int("particle", nrow(theta))
  source_row <- seq_len(nrow(theta))
  weight_all <- w * (1 - if (stress_pool_size > 0L && !is.null(factor_set)) as.numeric(stress_weight) else 0)

  if (stress_pool_size > 0L && !is.null(factor_set)) {
    if (!is.null(seed)) set.seed(as.integer(seed) + 7919L)
    moments <- .weighted_theta_moments(theta, w)
    pool <- as.integer(max(size, stress_pool_size))
    Z <- matrix(stats::rnorm(pool * ncol(theta)), nrow = pool, ncol = ncol(theta))
    L <- tryCatch(chol(moments$cov), error = function(e) diag(ncol(theta)))
    stress_theta <- sweep(Z %*% L * as.numeric(stress_scale), 2L, moments$center, "+")
    eig <- eigen(moments$cov, symmetric = TRUE)
    pc_levels <- c(1.5, 2.25, 3.0)
    pc_theta <- do.call(rbind, lapply(seq_len(ncol(theta)), function(j) {
      step <- sqrt(max(eig$values[j], 0)) * eig$vectors[, j]
      do.call(rbind, lapply(pc_levels, function(level) {
        rbind(moments$center + level * step, moments$center - level * step)
      }))
    }))
    stress_theta <- rbind(stress_theta, pc_theta)
    colnames(stress_theta) <- model$hyper_names
    stress_lp <- population_factor_set_logposterior(
      factor_set,
      stress_theta,
      include_constant = FALSE,
      n_cores = n_cores
    )
    keep <- is.finite(stress_lp)
    if (any(keep)) {
      keep[keep] <- stress_lp[keep] >= max(stress_lp[keep]) - as.numeric(stress_log_drop)
    }
    stress_theta <- stress_theta[keep, , drop = FALSE]
    stress_lp <- stress_lp[keep]
    if (nrow(stress_theta)) {
      stress_w <- exp(pmin(stress_lp - max(stress_lp), 0) / 2)
      stress_w <- normalize_particle_weights(stress_w) * as.numeric(stress_weight)
      theta_all <- rbind(theta_all, stress_theta)
      source <- c(source, rep.int("stress", nrow(stress_theta)))
      source_row <- c(source_row, seq_len(nrow(stress_theta)))
      weight_all <- c(weight_all, stress_w)
    }
  }

  stress_rows <- which(source == "stress")
  particle_rows <- which(source == "particle")
  n_stress <- if (length(stress_rows)) {
    min(length(stress_rows), ceiling(size * as.numeric(stress_weight)))
  } else {
    0L
  }
  n_particle <- max(0L, size - n_stress)
  idx <- integer(0)
  if (n_particle > 0L && length(particle_rows)) {
    idx <- c(
      idx,
      particle_rows[.certification_theta_design(
        theta = theta_all[particle_rows, , drop = FALSE],
        w = normalize_particle_weights(weight_all[particle_rows]),
        size = n_particle,
        seed = seed
      )]
    )
  }
  if (n_stress > 0L) {
    idx <- c(
      idx,
      stress_rows[.certification_theta_design(
        theta = theta_all[stress_rows, , drop = FALSE],
        w = rep(1 / length(stress_rows), length(stress_rows)),
        size = n_stress,
        seed = if (is.null(seed)) NULL else seed + 104729L
      )]
    )
  }
  idx <- unique(idx)
  if (length(idx) < size) {
    fill_from <- setdiff(seq_len(nrow(theta_all)), idx)
    if (length(fill_from)) {
      fill_idx <- .certification_theta_design(
        theta = theta_all[fill_from, , drop = FALSE],
        w = normalize_particle_weights(weight_all[fill_from]),
        size = min(length(fill_from), size - length(idx)),
        seed = if (is.null(seed)) NULL else seed + 1299709L
      )
      idx <- c(idx, fill_from[fill_idx])
    }
  }
  idx <- idx[seq_len(min(length(idx), size))]

  list(
    theta = theta_all[idx, , drop = FALSE],
    source_rows = source_row[idx],
    source = source[idx],
    source_weights = normalize_particle_weights(weight_all[idx])
  )
}

.certification_theta_design <- function(theta,
                                        w,
                                        size,
                                        seed = NULL) {
  theta <- as.matrix(theta)
  n <- nrow(theta)
  size <- as.integer(max(1L, min(size, n)))
  w <- normalize_particle_weights(w)
  if (!is.null(seed)) set.seed(as.integer(seed))
  if (size >= n) return(seq_len(n))

  Z <- tryCatch(.outer_whiten_theta(theta, w = w), error = function(e) scale(theta))
  Z[!is.finite(Z)] <- 0
  selected <- integer(0)

  if (n > 1L) {
    D2 <- as.matrix(stats::dist(Z))^2
    selected <- which.min(as.numeric(D2 %*% w))
  } else {
    selected <- 1L
  }

  while (length(selected) < size) {
    dist_to_selected <- vapply(seq_len(n), function(i) {
      min(rowSums((Z[selected, , drop = FALSE] -
                     matrix(Z[i, ], nrow = length(selected), ncol = ncol(Z), byrow = TRUE))^2))
    }, numeric(1))
    dist_to_selected[selected] <- 0
    score <- w * dist_to_selected
    if (!any(is.finite(score) & score > 0)) break
    selected <- c(selected, which.max(score))
  }

  if (length(selected) < size) {
    available <- setdiff(seq_len(n), selected)
    fill <- sample(available, size = size - length(selected), prob = w[available])
    selected <- c(selected, fill)
  }
  unique(as.integer(selected))[seq_len(size)]
}

select_population_certification_points <- function(population_fit,
                                                   size = 8L,
                                                   seed = NULL) {
  if (is.null(population_fit$theta) || is.null(population_fit$w)) {
    stop("population_fit must contain theta and w.")
  }
  idx <- .certification_theta_design(
    theta = population_fit$theta,
    w = population_fit$w,
    size = size,
    seed = seed
  )
  theta <- population_fit$theta[idx, , drop = FALSE]
  w <- normalize_particle_weights(population_fit$w)
  list(
    theta = theta,
    source_rows = idx,
    source_weights = normalize_particle_weights(w[idx])
  )
}

run_local_marginal_qmc_certification <- function(data_list,
                                                 loglik_fn,
                                                 factor_set,
                                                 population_model,
                                                 theta,
                                                 qmc_size = 8192L,
                                                 qmc_randomizations = 2L,
                                                 local_subset = NULL,
                                                 n_jobs = 1L,
                                                 local_n_cores = 1L,
                                                 seed = 123L) {
  model <- normalize_population_model(population_model)
  theta <- .as_hyper_matrix(theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  theta_prepared <- population_model_prepare_theta(model, theta)
  if (!identical(theta_prepared$family, "gaussian") ||
      !identical(theta_prepared$quadratic_kind, "diag")) {
    stop("QMC certification requires a diagonal Gaussian population model.")
  }
  local_subset <- sort(unique(as.integer(local_subset %||% seq_along(data_list))))
  local_subset <- local_subset[local_subset >= 1L & local_subset <= length(data_list)]
  if (!length(local_subset)) stop("local_subset selected no valid locals.")

  base_by_local <- population_factor_set_loglik_by_local(
    factor_set,
    theta = theta,
    include_constant = TRUE,
    n_cores = n_jobs
  )

  qmc_size <- as.integer(max(16L, qmc_size))
  qmc_randomizations <- as.integer(max(1L, qmc_randomizations))
  z_list <- lapply(seq_len(qmc_randomizations), function(r) {
    set.seed(as.integer(seed + r - 1L))
    U <- qrng::sobol(n = qmc_size, d = model$alpha_dim, randomize = TRUE)
    qnorm(pmin(pmax(U, 1e-12), 1 - 1e-12))
  })

  local_results <- parallel::mclapply(
    local_subset,
    function(local_id) {
      rows <- vector("list", nrow(theta))
      for (k in seq_len(nrow(theta))) {
        sigma <- sqrt(1 / theta_prepared$quadratic_coef[k, ])
        logm_rep <- vapply(seq_len(qmc_randomizations), function(r) {
          alpha <- sweep(z_list[[r]], 2L, sigma, "*")
          alpha <- sweep(alpha, 2L, theta_prepared$mean[k, ], "+")
          colnames(alpha) <- model$alpha_names
          ll <- ll_parallel(alpha, data_list[[local_id]], loglik_fn, n_cores = local_n_cores)
          ll[!is.finite(ll)] <- -Inf
          logsumexp(ll) - log(length(ll))
        }, numeric(1))
        certified <- logsumexp(logm_rep) - log(length(logm_rep))
        mcse <- if (length(logm_rep) > 1L) stats::sd(logm_rep) / sqrt(length(logm_rep)) else NA_real_
        rows[[k]] <- data.frame(
          theta_row = k,
          local = as.integer(local_id),
          base_logm = as.numeric(base_by_local[k, local_id]),
          certified_logm = as.numeric(certified),
          delta = as.numeric(certified - base_by_local[k, local_id]),
          certified_mcse = as.numeric(mcse),
          rounds = NA_integer_,
          final_lambda = 1,
          check.names = FALSE
        )
      }
      do.call(rbind, rows)
    },
    mc.cores = as.integer(max(1L, n_jobs))
  )
  summary <- do.call(rbind, local_results)
  summary <- summary[order(summary$theta_row, summary$local), , drop = FALSE]
  total <- do.call(rbind, lapply(seq_len(nrow(theta)), function(k) {
    rows <- summary[summary$theta_row == k, , drop = FALSE]
    data.frame(
      theta_row = k,
      base_loglik = sum(rows$base_logm),
      certified_loglik = sum(rows$certified_logm),
      delta = sum(rows$delta),
      certified_mcse = sqrt(sum(rows$certified_mcse^2, na.rm = TRUE)),
      n_locals = nrow(rows),
      check.names = FALSE
    )
  }))
  list(theta = theta, local = summary, total = total)
}

.append_certification_bank <- function(bank, certification, weights = NULL) {
  weights <- normalize_particle_weights(weights %||% rep(1, nrow(certification$theta)))
  add <- list(
    theta = certification$theta,
    delta = certification$total$delta,
    weights = weights,
    total = certification$total,
    local = certification$local
  )
  if (is.null(bank)) return(add)
  list(
    theta = rbind(bank$theta, add$theta),
    delta = c(bank$delta, add$delta),
    weights = normalize_particle_weights(c(bank$weights, add$weights)),
    total = rbind(bank$total, add$total),
    local = rbind(bank$local, add$local)
  )
}

fit_certified_population_model <- function(data_list,
                                           loglik_fn,
                                           local_objects,
                                           population_model,
                                           outer_control = list(),
                                           certification_control = list(),
                                           n_cores = 1L,
                                           seed = 123L,
                                           verbose = TRUE) {
  model <- normalize_population_model(population_model)
  control <- modifyList(
    list(
      max_repairs = 1L,
      audit_theta_size = 24L,
      audit_local_size = NULL,
      qmc_size = 8192L,
      qmc_randomizations = 2L,
      max_components_per_local = 10L,
      max_particles_per_local = 5000L,
      factor_compression_defensive_weight = 0.05,
      factor_compression_design_size = 32L,
      repair_anchors = NULL,
      repair_particles = NULL,
      repair_qmc_randomizations = NULL,
      repair_theta_delta_coverage = 0.90,
      repair_min_anchors = 1L,
      repair_local_delta_coverage = 0.90,
      repair_local_delta_min = 0.10,
      repair_local_ess_threshold = 50,
      repair_min_locals_per_anchor = 1L,
      repair_max_locals_per_anchor = Inf,
      local_n_cores = 1L,
      design_stress_pool_size = 384L,
      design_stress_scale = 2.5,
      design_stress_weight = 0.35,
      design_stress_log_drop = 30,
      validation_rmse_tol = 0.75,
      validation_median_abs_tol = 0.50,
      require_validation = FALSE
    ),
    certification_control
  )

  proposal_objects_current <- budget_reference_local_objects(
    local_objects,
    max_components = control$max_components_per_local
  )
  compression_theta <- .local_object_anchor_theta_design(
    proposal_objects_current,
    population_model = model,
    max_rows = control$factor_compression_design_size
  )
  local_objects_current <- compress_reference_local_objects_to_factors(
    proposal_objects_current,
    population_model = model,
    theta_design = compression_theta,
    data_list = data_list,
    loglik_fn = loglik_fn,
    max_particles = control$max_particles_per_local,
    defensive_weight = control$factor_compression_defensive_weight,
    local_n_cores = control$local_n_cores
  )
  current_factor_set <- build_population_factor_set(
    local_objects_current,
    model
  )
  outer_args <- modifyList(
    list(
      factor_set = current_factor_set,
      N = 2000L,
      n_mcmc_moves = 3L,
      max_rounds = 80L,
      n_cores = as.integer(n_cores),
      seed = as.integer(seed),
      verbose = verbose
    ),
    outer_control
  )
  current_fit <- do.call(outer_population_smc, outer_args)
  base_fit <- current_fit
  base_factor_set <- current_factor_set
  history <- list()
  certification_bank <- NULL
  max_repairs <- as.integer(max(0L, control$max_repairs))
  repair_qmc_randomizations <- as.integer(control$repair_qmc_randomizations %||% control$qmc_randomizations)
  repair_particle_ceiling <- as.integer(control$repair_particles %||% min(as.integer(control$qmc_size), 2048L))

  choose_audit_locals <- function(iter) {
    if (!is.null(control$local_subset)) {
      out <- sort(unique(as.integer(control$local_subset)))
    } else {
      audit_size <- control$audit_local_size %||% if (length(data_list) <= 64L) length(data_list) else 64L
      audit_size <- as.integer(max(1L, min(audit_size, length(data_list))))
      out <- if (audit_size >= length(data_list)) {
        seq_along(data_list)
      } else {
        select_stratified_pilot_subset(
          data_list = data_list,
          pilot_size = audit_size,
          seed = seed + 70000L + iter
        )$indices
      }
    }
    out[out >= 1L & out <= length(data_list)]
  }

  audit_current_fit <- function(iter, factor_set, fit) {
    audit_locals <- choose_audit_locals(iter)
    design <- select_population_certification_design(
      population_fit = fit,
      factor_set = factor_set,
      population_model = model,
      size = control$audit_theta_size,
      seed = seed + 10000L * (iter + 1L),
      stress_pool_size = control$design_stress_pool_size,
      stress_scale = control$design_stress_scale,
      stress_weight = control$design_stress_weight,
      stress_log_drop = control$design_stress_log_drop,
      n_cores = n_cores
    )
    audit <- run_local_marginal_qmc_certification(
      data_list = data_list,
      loglik_fn = loglik_fn,
      factor_set = factor_set,
      population_model = model,
      theta = design$theta,
      qmc_size = control$qmc_size,
      qmc_randomizations = control$qmc_randomizations,
      local_subset = audit_locals,
      n_jobs = n_cores,
      local_n_cores = control$local_n_cores,
      seed = seed + 20000L * (iter + 1L)
    )
    scale <- length(data_list) / length(audit_locals)
    residual <- audit$total$delta * scale
    local_residual <- audit$local$delta
    summary <- data.frame(
      iteration = iter,
      audit_locals = length(audit_locals),
      rmse = sqrt(mean(residual^2)),
      median_abs = median(abs(residual)),
      max_abs = max(abs(residual)),
      mean = mean(residual),
      local_rmse = sqrt(mean(local_residual^2)),
      certified_mcse = sqrt(mean((audit$total$certified_mcse * scale)^2)),
      check.names = FALSE
    )
    list(
      design = design,
      audit = audit,
      audit_locals = audit_locals,
      scaled_residual = residual,
      validation_summary = summary
    )
  }

  local_objects_for_ids <- function(local_objects_in, local_ids) {
    local_ids <- sort(unique(as.integer(local_ids)))
    object_ids <- vapply(local_objects_in, function(x) as.integer(x$local_id %||% NA_integer_), integer(1))
    idx <- match(local_ids, object_ids)
    idx <- idx[!is.na(idx)]
    local_objects_in[idx]
  }

  repair_anchor_budget <- function(audited, local_objects_in) {
    theta_score <- abs(audited$audit$total$delta)
    if (all(!is.finite(theta_score)) || !any(theta_score > 0)) {
      theta_score <- rep(1, length(theta_score))
    }
    ord <- order(theta_score, decreasing = TRUE)
    target <- as.numeric(control$repair_theta_delta_coverage) * sum(theta_score[ord])
    residual_need <- if (target > 0) {
      which(cumsum(theta_score[ord]) >= target)[1L]
    } else {
      as.integer(control$repair_min_anchors)
    }
    residual_need <- max(as.integer(control$repair_min_anchors), residual_need %||% 1L)

    component_capacity <- length(theta_score)
    if (!is.null(control$max_components_per_local) && is.finite(control$max_components_per_local)) {
      audited_objects <- local_objects_for_ids(local_objects_in, audited$audit_locals)
      if (!length(audited_objects)) return(0L)
      reserved_count <- vapply(audited_objects, function(obj) {
        components <- validate_reference_local_object(obj)$components
        as.integer(sum(vapply(components, function(component) {
          diagnostics <- component$diagnostics %||% list()
          identical(component$proposal_type, "posterior_reference") ||
            identical(as.character(diagnostics$role %||% ""), "planned_anchor")
        }, logical(1))))
      }, integer(1))
      component_capacity <- max(0L, min(as.integer(control$max_components_per_local) - reserved_count))
    }
    requested <- control$repair_anchors
    if (!is.null(requested)) residual_need <- min(residual_need, as.integer(requested))
    as.integer(max(0L, min(length(theta_score), residual_need, component_capacity)))
  }

  repair_particle_budget <- function(local_objects_in, changed_locals, n_anchors) {
    if (!length(changed_locals) || n_anchors <= 0L) return(0L)
    qmc_size <- repair_particle_ceiling
    if (!is.finite(qmc_size) || qmc_size < 16L) return(0L)
    as.integer(qmc_size)
  }

  select_repair_plan <- function(certification, factor_set, max_anchors) {
    theta <- certification$theta
    total_delta <- certification$total$delta
    score <- abs(total_delta)
    if (all(!is.finite(score)) || !any(score > 0)) {
      score <- rep(1, length(total_delta))
    }
    Z <- tryCatch(.outer_whiten_theta(theta, w = rep(1 / nrow(theta), nrow(theta))), error = function(e) scale(theta))
    Z[!is.finite(Z)] <- 0
    max_anchors <- as.integer(max(0L, min(max_anchors, length(score))))
    if (max_anchors <= 0L) {
      return(list(
        anchors = data.frame(),
        locals = data.frame()
      ))
    }

    first <- which.max(score)
    anchor_order <- first
    while (length(anchor_order) < min(length(score), as.integer(max_anchors))) {
      remaining <- setdiff(seq_along(score), anchor_order)
      dist_to_selected <- vapply(remaining, function(i) {
        min(rowSums((Z[anchor_order, , drop = FALSE] -
                       matrix(Z[i, ], nrow = length(anchor_order), ncol = ncol(Z), byrow = TRUE))^2))
      }, numeric(1))
      anchor_order <- c(anchor_order, remaining[which.max(score[remaining] * sqrt(pmax(dist_to_selected, 1e-8)))])
    }

    anchor_rows <- vector("list", length(anchor_order))
    local_rows <- list()
    for (a in seq_along(anchor_order)) {
      theta_row <- anchor_order[a]
      theta_one <- theta[theta_row, , drop = FALSE]
      local_delta <- certification$local$delta[certification$local$theta_row == theta_row]
      local_id <- certification$local$local[certification$local$theta_row == theta_row]
      local_ess <- population_factor_set_local_ess(factor_set, theta_one)
      local_ess <- local_ess[local_id]

      local_score <- abs(local_delta)
      ord <- order(local_score, decreasing = TRUE)
      total_abs <- sum(local_score)
      max_locals <- if (is.finite(control$repair_max_locals_per_anchor)) {
        as.integer(control$repair_max_locals_per_anchor)
      } else {
        length(ord)
      }
      cover_n <- if (total_abs > 0) {
        which(cumsum(local_score[ord]) >= as.numeric(control$repair_local_delta_coverage) * total_abs)[1L]
      } else {
        as.integer(control$repair_min_locals_per_anchor)
      }
      cover_n <- max(as.integer(control$repair_min_locals_per_anchor), cover_n %||% 0L)
      cover_n <- min(length(ord), cover_n, max_locals)
      selected <- ord[seq_len(max(1L, cover_n))]
      low_ess <- which(local_ess < as.numeric(control$repair_local_ess_threshold))
      high_delta <- which(local_score >= as.numeric(control$repair_local_delta_min))
      selected <- sort(unique(c(selected, low_ess, high_delta)))
      if (is.finite(control$repair_max_locals_per_anchor)) {
        selected <- selected[order(local_score[selected], decreasing = TRUE)]
        selected <- selected[seq_len(min(length(selected), as.integer(control$repair_max_locals_per_anchor)))]
      }
      selected_local <- local_id[selected]
      reason <- ifelse(local_ess[selected] < as.numeric(control$repair_local_ess_threshold), "low_ess", "delta")

      anchor_rows[[a]] <- data.frame(
        anchor = a,
        theta_row = theta_row,
        total_delta = total_delta[theta_row],
        score = score[theta_row],
        n_locals = length(selected_local),
        median_selected_ess = median(local_ess[match(selected_local, local_id)]),
        check.names = FALSE
      )
      local_rows[[a]] <- data.frame(
        anchor = a,
        theta_row = theta_row,
        local = selected_local,
        base_logm = certification$local$base_logm[certification$local$theta_row == theta_row][match(selected_local, local_id)],
        certified_logm = certification$local$certified_logm[certification$local$theta_row == theta_row][match(selected_local, local_id)],
        certified_mcse = certification$local$certified_mcse[certification$local$theta_row == theta_row][match(selected_local, local_id)],
        local_delta = local_delta[match(selected_local, local_id)],
        local_ess = local_ess[match(selected_local, local_id)],
        reason = reason,
        check.names = FALSE
      )
    }

    list(
      anchors = do.call(rbind, anchor_rows),
      locals = do.call(rbind, local_rows)
    )
  }

  repair_from_plan <- function(local_objects_in, theta, plan, iter, qmc_size) {
    local_objects_out <- local_objects_in
    components_by_anchor <- vector("list", nrow(plan$anchors))
    for (a in seq_len(nrow(plan$anchors))) {
      theta_row <- plan$anchors$theta_row[a]
      theta_anchor <- theta[theta_row, , drop = FALSE]
      local_ids <- sort(unique(plan$locals$local[plan$locals$anchor == plan$anchors$anchor[a]]))
      anchor_components <- build_qmc_population_anchor_local_objects(
        data_list = data_list,
        loglik_fn = loglik_fn,
        population_model = model,
        theta = theta_anchor,
        local_ids = local_ids,
        qmc_size = qmc_size,
        qmc_randomizations = repair_qmc_randomizations,
        seed = seed + 60000L * iter + 1000L * a,
        n_jobs = n_cores,
        local_n_cores = control$local_n_cores,
        label = sprintf("audit_repair%d_anchor%d", iter, a),
        diagnostics = list(
          role = "audit_repair",
          audit_iteration = iter,
          theta_row = theta_row,
          anchor_score = plan$anchors$score[a],
          total_delta = plan$anchors$total_delta[a]
        )
      )
      local_objects_out <- add_local_reference_components_dmis(
        local_objects_out,
        anchor_components$local_objects
      )
      components_by_anchor[[a]] <- list(
        theta_row = theta_row,
        theta = theta_anchor,
        local_ids = local_ids,
        reference_prior = anchor_components$reference_prior,
        local_objects = anchor_components$local_objects
      )
    }
    list(local_objects = local_objects_out, components = components_by_anchor)
  }

  selected_iteration <- NA_integer_
  validated <- FALSE

  for (iter in seq.int(0L, max_repairs)) {
    audited <- audit_current_fit(iter, current_factor_set, current_fit)
    audit <- audited$audit
    validation_summary <- audited$validation_summary
    if (isTRUE(verbose)) {
      cat(sprintf(
        "Population audit %d: rmse=%.3f median_abs=%.3f max_abs=%.3f across %d locals.\n",
        iter,
        validation_summary$rmse,
        validation_summary$median_abs,
        validation_summary$max_abs,
        validation_summary$audit_locals
      ))
    }

    history[[iter + 1L]] <- list(
      iteration = iter,
      design = audited$design,
      audit = audit,
      audit_locals = audited$audit_locals,
      local_objects = local_objects_current,
      factor_set = current_factor_set,
      fit = current_fit,
      validation_residual = audited$scaled_residual,
      validation_summary = validation_summary
    )
    certification_bank <- .append_certification_bank(
      certification_bank,
      audit,
      weights = audited$design$source_weights
    )
    passed <- validation_summary$rmse <= control$validation_rmse_tol &&
      validation_summary$median_abs <= control$validation_median_abs_tol
    if (passed) {
      selected_iteration <- iter + 1L
      validated <- TRUE
      break
    }
    if (iter >= max_repairs) break

    max_repair_anchors <- repair_anchor_budget(
      audited,
      local_objects_in = proposal_objects_current
    )
    if (max_repair_anchors <= 0L) {
      if (isTRUE(verbose)) {
        cat(sprintf(
          "Population audit %d: no repair capacity left under the local factor budget.\n",
          iter
        ))
      }
      break
    }

    repair_plan <- select_repair_plan(
      audit,
      factor_set = current_factor_set,
      max_anchors = max_repair_anchors
    )
    changed_locals <- sort(unique(repair_plan$locals$local))
    repair_particles <- repair_particle_budget(
      local_objects_current,
      changed_locals = changed_locals,
      n_anchors = nrow(repair_plan$anchors)
    )
    if (!nrow(repair_plan$anchors) || !length(changed_locals) || repair_particles <= 0L) {
      if (isTRUE(verbose)) {
        cat(sprintf(
          "Population audit %d: repair skipped because the local factor particle budget is exhausted.\n",
          iter
        ))
      }
      break
    }
    if (isTRUE(verbose)) {
      cat(sprintf(
        "Population audit %d: repairing %d anchors and %d local factor shards (%d QMC particles x %d randomizations each).\n",
        iter,
        nrow(repair_plan$anchors),
        length(changed_locals),
        repair_particles,
        repair_qmc_randomizations
      ))
    }
    repair <- repair_from_plan(
      local_objects_in = proposal_objects_current,
      theta = audit$theta,
      plan = repair_plan,
      iter = iter + 1L,
      qmc_size = repair_particles
    )
    proposal_objects_current <- budget_reference_local_objects(
      repair$local_objects,
      max_components = control$max_components_per_local,
      local_ids = changed_locals
    )
    compression_theta <- .local_object_anchor_theta_design(
      proposal_objects_current,
      population_model = model,
      extra_theta = audit$theta[repair_plan$anchors$theta_row, , drop = FALSE],
      max_rows = control$factor_compression_design_size
    )
    local_objects_current <- compress_reference_local_objects_to_factors(
      proposal_objects_current,
      population_model = model,
      theta_design = compression_theta,
      data_list = data_list,
      loglik_fn = loglik_fn,
      local_ids = changed_locals,
      existing = local_objects_current,
      max_particles = control$max_particles_per_local,
      defensive_weight = control$factor_compression_defensive_weight,
      local_n_cores = control$local_n_cores
    )
    old_factor_set <- current_factor_set
    current_factor_set <- update_population_factor_set_locals(
      factor_set = current_factor_set,
      local_objects = local_objects_current,
      local_ids = changed_locals,
      data_list = data_list,
      loglik_fn = loglik_fn,
      local_n_cores = control$local_n_cores
    )
    current_fit <- update_outer_population_fit(
      fit = current_fit,
      old_factor_set = old_factor_set,
      new_factor_set = current_factor_set,
      n_mcmc_moves = outer_control$n_mcmc_moves %||% 3L,
      n_cores = n_cores,
      seed = seed + 30000L * (iter + 1L),
      verbose = verbose
    )
    history[[iter + 1L]]$repair_plan <- repair_plan
    history[[iter + 1L]]$repair <- repair
  }

  validation_scores <- vapply(history, function(h) {
    h$validation_summary$rmse + h$validation_summary$median_abs
  }, numeric(1))
  if (is.na(selected_iteration)) {
    selected_iteration <- which.min(validation_scores)
  }
  selected_history <- history[[selected_iteration]]
  status <- if (isTRUE(validated)) {
    "validated"
  } else if (isTRUE(control$require_validation)) {
    current_fit <- selected_history$fit
    current_factor_set <- selected_history$factor_set
    local_objects_current <- selected_history$local_objects
    "failed_validation_selected_best"
  } else {
    current_fit <- selected_history$fit
    current_factor_set <- selected_history$factor_set
    local_objects_current <- selected_history$local_objects
    "best_unvalidated"
  }

  list(
    mode = "planned_anchor_audit",
    base_factor_set = base_factor_set,
    factor_set = current_factor_set,
    local_objects = local_objects_current,
    base_fit = base_fit,
    fit = current_fit,
    certification_bank = certification_bank,
    history = history,
    selected_iteration = selected_iteration,
    validated = validated,
    status = status,
    control = control
  )
}

compact_certified_population_result <- function(result) {
  compact_repair <- function(repair) {
    if (is.null(repair)) return(NULL)
    list(
      components = lapply(repair$components %||% list(), function(component) {
        list(
          theta_row = component$theta_row,
          theta = component$theta,
          local_ids = component$local_ids
        )
      })
    )
  }

  compact_history <- lapply(result$history %||% list(), function(h) {
    list(
      iteration = h$iteration,
      design = h$design,
      audit = h$audit,
      audit_locals = h$audit_locals,
      validation_residual = h$validation_residual,
      validation_summary = h$validation_summary,
      repair_plan = h$repair_plan,
      repair = compact_repair(h$repair)
    )
  })

  list(
    mode = result$mode,
    status = result$status,
    validated = result$validated,
    selected_iteration = result$selected_iteration,
    control = result$control,
    local_objects = result$local_objects,
    fit = result$fit,
    base_fit = result$base_fit,
    certification_bank = result$certification_bank,
    history = compact_history,
    factor_budget = .reference_local_budget_summary(result$local_objects),
    base_factor_budget = .reference_local_budget_summary(result$history[[1L]]$local_objects)
  )
}

prepare_reference_local_stage <- function(data_list,
                                          loglik_fn,
                                          base_mu,
                                          base_Sigma,
                                          population_model,
                                          pilot_size = 10L,
                                          broad_scale = 1,
                                          broad_defensive = FALSE,
                                          pilot_particles = 800L,
                                          full_particles = 2000L,
                                          inflation = 1.5,
                                          defensive_weight = 0.10,
                                          defensive_scale = 4,
                                          n_jobs = 1L,
                                          pilot_local_n_cores = 1L,
                                          full_local_n_cores = 1L,
                                          base_seed = 123L,
                                          feature_fn = default_local_feature_vector,
                                          n_strata = NULL,
                                          pilot_reference_support_size = 8L,
                                          max_components_per_local = 10L,
                                          planned_anchor_count = NULL,
                                          planned_anchor_qmc_size = 2048L,
                                          planned_anchor_qmc_randomizations = 2L,
                                          pilot_outer_control = list(N = 1000L, n_mcmc_moves = 2L, max_rounds = 50L),
                                          verbose = TRUE,
                                          pilot_smc_control = list(max_rounds = 40L),
                                          full_smc_control = list()) {
  population_model <- normalize_population_model(population_model)

  broad_reference <- make_broad_reference_prior(
    mu = base_mu,
    Sigma = base_Sigma,
    scale = broad_scale,
    defensive = broad_defensive,
    defensive_scale = broad_scale * defensive_scale,
    defensive_weight = defensive_weight
  )

  pilot <- select_stratified_pilot_subset(
    data_list = data_list,
    pilot_size = pilot_size,
    feature_fn = feature_fn,
    n_strata = n_strata,
    seed = base_seed
  )

  pilot_args <- modifyList(
    list(
      data_list = data_list,
      loglik_fn = loglik_fn,
      reference_prior = broad_reference,
      indices = pilot$indices,
      M = pilot_particles,
      n_jobs = n_jobs,
      local_n_cores = pilot_local_n_cores,
      base_seed = base_seed,
      verbose = verbose
    ),
    pilot_smc_control
  )
  pilot_fits <- do.call(run_reference_local_smc, pilot_args)

  pilot_population_fit <- fit_pilot_population_model(
    pilot_fits = pilot_fits,
    reference_prior = broad_reference,
    population_model = population_model,
    outer_control = pilot_outer_control,
    n_cores = n_jobs,
    seed = base_seed + 50000L,
    verbose = verbose
  )

  refined_reference <- build_refined_reference_prior_from_population_fit(
    population_fit = pilot_population_fit,
    population_model = population_model,
    inflation = inflation,
    defensive_weight = defensive_weight,
    broad_reference = broad_reference,
    defensive_scale = defensive_scale,
    support_size = pilot_reference_support_size,
    support_seed = base_seed + 50001L
  )

  if (is.null(planned_anchor_count)) {
    planned_anchor_count <- if (!is.null(max_components_per_local) && is.finite(max_components_per_local)) {
      min(4L, max(1L, floor(0.45 * (as.integer(max_components_per_local) - 1L))))
    } else {
      4L
    }
  }
  planned_anchor_count <- as.integer(max(1L, planned_anchor_count))
  planned_anchor_qmc_size <- as.integer(planned_anchor_qmc_size)
  planned_anchor_qmc_size <- as.integer(max(16L, planned_anchor_qmc_size))

  planned_anchor_theta <- select_population_anchor_support(
    population_fit = pilot_population_fit,
    population_model = population_model,
    n_anchors = planned_anchor_count,
    seed = base_seed + 50002L
  )
  planned_anchor_library <- build_planned_anchor_local_objects(
    data_list = data_list,
    loglik_fn = loglik_fn,
    population_model = population_model,
    anchor_theta = planned_anchor_theta,
    local_ids = seq_along(data_list),
    qmc_size = planned_anchor_qmc_size,
    qmc_randomizations = planned_anchor_qmc_randomizations,
    seed = base_seed + 60000L,
    n_jobs = n_jobs,
    local_n_cores = full_local_n_cores,
    label = "planned_population_anchor"
  )

  full_args <- modifyList(
    list(
      data_list = data_list,
      loglik_fn = loglik_fn,
      reference_prior = refined_reference,
      indices = seq_along(data_list),
      M = full_particles,
      n_jobs = n_jobs,
      local_n_cores = full_local_n_cores,
      base_seed = base_seed + 100000L,
      verbose = verbose
    ),
    full_smc_control
  )
  local_fits <- do.call(run_reference_local_smc, full_args)
  local_objects <- add_local_reference_components_dmis(
    build_local_reference_objects(local_fits, reference_prior = refined_reference),
    planned_anchor_library$local_objects
  )
  local_objects <- budget_reference_local_objects(
    local_objects,
    max_components = max_components_per_local
  )

  list(
    broad_reference = broad_reference,
    pilot = list(
      selection = pilot,
      fits = pilot_fits,
      summary = aggregate_pilot_local_moments(pilot_fits),
      population_fit = pilot_population_fit,
      refinement_source = "pilot_population_fit"
    ),
    refined_reference = refined_reference,
    planned_anchor_library = planned_anchor_library,
    local_fits = local_fits,
    local_objects = local_objects,
    outer_fit = NULL
  )
}
