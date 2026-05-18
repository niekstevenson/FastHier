#!/usr/bin/env Rscript
# ============================================================================
# Hierarchical local-reference workflow
# - Broad pilot reference
# - Population-informed refined defensive mixture
# - Full local pass, optionally followed by exactly one DMIS refresh
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

.local_reference_object_components <- function(local_object) {
  if (!is.null(local_object$components)) {
    return(local_object$components)
  }
  list(
    list(
      label = "",
      particles = as.matrix(local_object$particles),
      weights = normalize_particle_weights(local_object$weights),
      proposal_type = as.character(local_object$proposal_type %||% "posterior_reference"),
      log_likelihood = as.numeric(local_object$log_likelihood %||% rep(NA_real_, nrow(local_object$particles))),
      log_evidence = as.numeric(local_object$log_evidence %||% NA_real_),
      mcse_log_evidence = as.numeric(local_object$mcse_log_evidence %||% NA_real_),
      reference_prior = normalize_reference_prior(reference_prior = local_object$reference_prior),
      log_reference_density = as.numeric(
        local_object$log_reference_density %||%
          reference_prior_logpdf(local_object$reference_prior, local_object$particles)
      ),
      diagnostics = local_object$diagnostics %||% list()
    )
  )
}

merge_local_reference_objects_dmis <- function(primary_objects,
                                               secondary_objects,
                                               mixture_weights = c(0.5, 0.5),
                                               labels = c("initial", "refreshed")) {
  if (!length(primary_objects) || !length(secondary_objects)) {
    stop("Both primary_objects and secondary_objects must be non-empty.")
  }
  if (!identical(names(primary_objects), names(secondary_objects))) {
    stop("primary_objects and secondary_objects must have identical names for DMIS merging.")
  }
  mixture_weights <- normalize_particle_weights(mixture_weights)
  if (length(mixture_weights) != 2L) {
    stop("mixture_weights must have length 2 for DMIS merging.")
  }
  labels <- rep_len(as.character(labels), 2L)

  merged <- lapply(names(primary_objects), function(nm) {
    primary_components <- .local_reference_object_components(primary_objects[[nm]])
    secondary_components <- .local_reference_object_components(secondary_objects[[nm]])
    components <- c(primary_components, secondary_components)

    group_labels <- c(
      rep.int(labels[1L], length(primary_components)),
      rep.int(labels[2L], length(secondary_components))
    )
    if (length(components) > 1L) {
      for (idx in seq_along(components)) {
        components[[idx]]$label <- if (sum(group_labels == group_labels[idx]) > 1L) {
          paste0(group_labels[idx], "_", idx)
        } else {
          group_labels[idx]
        }
      }
    } else {
      components[[1L]]$label <- labels[1L]
    }

    component_weights <- c(
      rep(mixture_weights[1L] / length(primary_components), length(primary_components)),
      rep(mixture_weights[2L] / length(secondary_components), length(secondary_components))
    )

    structure(
      list(
        local_id = as.integer(primary_objects[[nm]]$local_id %||% secondary_objects[[nm]]$local_id %||% NA_integer_),
        components = components,
        mixture_weights = component_weights,
        diagnostics = list(mode = "dmis_merge")
      ),
      class = "reference_local_object"
    )
  })
  names(merged) <- names(primary_objects)
  merged
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
                                                      label = "qmc_population_anchor") {
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
        diagnostics = list(
          method = "qmc_population_anchor",
          qmc_size = qmc_size,
          qmc_randomizations = qmc_randomizations
        )
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

rebalance_local_reference_object_components <- function(local_object,
                                                        component_weights = NULL) {
  if (!is.list(local_object) || is.null(local_object$components)) {
    local_object <- validate_reference_local_object(local_object)
  } else {
    local_object$components <- lapply(local_object$components, validate_reference_local_component)
  }
  n_components <- length(local_object$components)
  weights <- normalize_particle_weights(component_weights %||% rep(1, n_components))
  if (length(weights) != n_components) {
    stop("component_weights must match the number of local components.")
  }
  local_object$mixture_weights <- weights
  local_object
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
    component_types <- vapply(current$components, function(component) {
      as.character(component$proposal_type %||% "posterior_reference")
    }, character(1))
    if (any(component_types == "prior_reference")) {
      current$components <- current$components[component_types == "prior_reference"]
    }
    current$diagnostics <- modifyList(
      current$diagnostics %||% list(),
      list(mode = "adaptive_dmis", n_components = length(current$components))
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

run_local_marginal_certification <- function(data_list,
                                             loglik_fn,
                                             factor_set,
                                             population_model,
                                             theta,
                                             particles = 4000L,
                                             local_subset = NULL,
                                             n_jobs = 1L,
                                             local_n_cores = 1L,
                                             base_seed = 123L,
                                             smc_control = list(),
                                             verbose = TRUE) {
  model <- normalize_population_model(population_model)
  theta <- .as_hyper_matrix(theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  local_subset <- sort(unique(as.integer(local_subset %||% seq_along(data_list))))
  local_subset <- local_subset[local_subset >= 1L & local_subset <= length(data_list)]
  if (!length(local_subset)) stop("local_subset selected no valid locals.")

  base_by_local <- population_factor_set_loglik_by_local(
    factor_set,
    theta = theta,
    include_constant = TRUE,
    n_cores = n_jobs
  )

  tasks <- expand.grid(
    theta_row = seq_len(nrow(theta)),
    local = local_subset,
    KEEP.OUT.ATTRS = FALSE
  )
  tasks <- tasks[order(tasks$theta_row, tasks$local), , drop = FALSE]
  particles <- as.integer(particles)
  if (particles <= 1L) stop("particles must be greater than one.")

  task_results <- parallel::mclapply(
    seq_len(nrow(tasks)),
    function(task_id) {
      task <- tasks[task_id, , drop = FALSE]
      theta_row <- theta[task$theta_row, , drop = FALSE]
      reference_prior <- make_population_reference_prior_from_theta(
        population_model = model,
        theta = theta_row,
        label = sprintf("cert_theta%d_local%d", task$theta_row, task$local)
      )
      fit_args <- modifyList(
        list(
          data = data_list[[task$local]],
          loglik_fn = loglik_fn,
          reference_prior = reference_prior,
          M = particles,
          n_cores = as.integer(local_n_cores),
          seed = as.integer(base_seed + task_id - 1L),
          verbose = verbose
        ),
        smc_control
      )
      fit <- do.call(enhanced_smc_elite, fit_args)
      data.frame(
        theta_row = as.integer(task$theta_row),
        local = as.integer(task$local),
        base_logm = as.numeric(base_by_local[task$theta_row, task$local]),
        certified_logm = as.numeric(fit$log_evidence),
        delta = as.numeric(fit$log_evidence - base_by_local[task$theta_row, task$local]),
        certified_mcse = as.numeric(fit$mcse_logZ %||% NA_real_),
        rounds = as.integer(fit$meta$rounds %||% NA_integer_),
        final_lambda = as.numeric(fit$final_lambda %||% NA_real_),
        check.names = FALSE
      )
    },
    mc.cores = as.integer(max(1L, n_jobs))
  )

  summary <- do.call(rbind, task_results)
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

  list(
    theta = theta,
    local = summary,
    total = total
  )
}

.gauss_hermite_rule_standard_normal <- function(order) {
  order <- as.integer(order)
  if (order < 2L) stop("Gauss-Hermite order must be at least two.")
  J <- matrix(0, order, order)
  off <- sqrt(seq_len(order - 1L) / 2)
  J[cbind(seq_len(order - 1L), 2:order)] <- off
  J[cbind(2:order, seq_len(order - 1L))] <- off
  eig <- eigen(J, symmetric = TRUE)
  ord <- order(eig$values)
  nodes <- sqrt(2) * eig$values[ord]
  weights <- (eig$vectors[1L, ord]^2)
  weights <- weights / sum(weights)
  list(nodes = as.numeric(nodes), weights = as.numeric(weights))
}

.standard_normal_tensor_grid <- function(dim, order) {
  rule <- .gauss_hermite_rule_standard_normal(order)
  grid_index <- expand.grid(rep(list(seq_along(rule$nodes)), dim))
  z <- as.matrix(data.frame(lapply(grid_index, function(idx) rule$nodes[idx])))
  logw <- rowSums(as.matrix(data.frame(lapply(grid_index, function(idx) log(rule$weights[idx])))))
  colnames(z) <- paste0("z", seq_len(dim))
  list(z = z, logw = as.numeric(logw))
}

run_local_marginal_quadrature_certification <- function(data_list,
                                                        loglik_fn,
                                                        factor_set,
                                                        population_model,
                                                        theta,
                                                        order = 7L,
                                                        local_subset = NULL,
                                                        n_jobs = 1L,
                                                        local_n_cores = 1L) {
  model <- normalize_population_model(population_model)
  theta <- .as_hyper_matrix(theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  theta_prepared <- population_model_prepare_theta(model, theta)
  if (!identical(theta_prepared$family, "gaussian") ||
      !identical(theta_prepared$quadratic_kind, "diag")) {
    stop("Gauss-Hermite certification requires a diagonal Gaussian population model.")
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

  grid <- .standard_normal_tensor_grid(model$alpha_dim, order = order)
  local_results <- parallel::mclapply(
    local_subset,
    function(local_id) {
      rows <- vector("list", nrow(theta))
      for (k in seq_len(nrow(theta))) {
        sigma <- sqrt(1 / theta_prepared$quadratic_coef[k, ])
        alpha <- sweep(grid$z, 2L, sigma, "*")
        alpha <- sweep(alpha, 2L, theta_prepared$mean[k, ], "+")
        colnames(alpha) <- model$alpha_names
        ll <- ll_parallel(alpha, data_list[[local_id]], loglik_fn, n_cores = local_n_cores)
        ll[!is.finite(ll)] <- -Inf
        certified <- logsumexp(grid$logw + ll)
        rows[[k]] <- data.frame(
          theta_row = k,
          local = as.integer(local_id),
          base_logm = as.numeric(base_by_local[k, local_id]),
          certified_logm = as.numeric(certified),
          delta = as.numeric(certified - base_by_local[k, local_id]),
          certified_mcse = 0,
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
      certified_mcse = 0,
      n_locals = nrow(rows),
      check.names = FALSE
    )
  }))
  list(theta = theta, local = summary, total = total)
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
      max_iterations = 3L,
      calibration_size = 32L,
      validation_size = 24L,
      certification_estimator = "qmc",
      particles = 4000L,
      enrichment_particles = NULL,
      enrichment_qmc_randomizations = NULL,
      enrichment_anchors = 4L,
      enrichment_local_delta_coverage = 0.90,
      enrichment_local_delta_min = 0.10,
      enrichment_local_ess_threshold = 50,
      enrichment_min_locals_per_anchor = 1L,
      enrichment_max_locals_per_anchor = Inf,
      enrich_all_locals = FALSE,
      quadrature_order = 7L,
      qmc_size = 8192L,
      qmc_randomizations = 2L,
      local_subset = NULL,
      local_n_cores = 1L,
      smc_control = list(max_rounds = 80L, n_mcmc_moves = 3L, G_mix = 12L),
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

  local_objects_current <- local_objects
  base_factor_set <- build_population_factor_set(
    local_objects_current,
    model,
    data_list = data_list,
    loglik_fn = loglik_fn,
    local_n_cores = control$local_n_cores
  )
  outer_args <- modifyList(
    list(
      factor_set = base_factor_set,
      N = 2000L,
      n_mcmc_moves = 3L,
      max_rounds = 80L,
      n_cores = as.integer(n_cores),
      seed = as.integer(seed),
      verbose = verbose
    ),
    outer_control
  )
  base_fit <- do.call(outer_population_smc, outer_args)
  current_fit <- base_fit
  current_factor_set <- base_factor_set
  history <- list()
  certification_bank <- NULL
  max_iterations <- as.integer(max(1L, control$max_iterations))
  enrichment_particles <- as.integer(control$enrichment_particles %||% min(as.integer(control$qmc_size), 2048L))
  enrichment_qmc_randomizations <- as.integer(control$enrichment_qmc_randomizations %||% control$qmc_randomizations)

  certify_points <- function(theta, factor_set, seed_offset) {
    estimator <- match.arg(
      as.character(control$certification_estimator),
      choices = c("qmc", "gh_quadrature", "smc")
    )
    if (identical(estimator, "qmc")) {
      run_local_marginal_qmc_certification(
        data_list = data_list,
        loglik_fn = loglik_fn,
        factor_set = factor_set,
        population_model = model,
        theta = theta,
        qmc_size = control$qmc_size,
        qmc_randomizations = control$qmc_randomizations,
        local_subset = control$local_subset,
        n_jobs = n_cores,
        local_n_cores = control$local_n_cores,
        seed = seed + seed_offset
      )
    } else if (identical(estimator, "gh_quadrature")) {
      run_local_marginal_quadrature_certification(
        data_list = data_list,
        loglik_fn = loglik_fn,
        factor_set = factor_set,
        population_model = model,
        theta = theta,
        order = control$quadrature_order,
        local_subset = control$local_subset,
        n_jobs = n_cores,
        local_n_cores = control$local_n_cores
      )
    } else {
      run_local_marginal_certification(
        data_list = data_list,
        loglik_fn = loglik_fn,
        factor_set = factor_set,
        population_model = model,
        theta = theta,
        particles = control$particles,
        local_subset = control$local_subset,
        n_jobs = n_cores,
        local_n_cores = control$local_n_cores,
        base_seed = seed + seed_offset,
        smc_control = control$smc_control,
        verbose = verbose
      )
    }
  }

  select_enrichment_plan <- function(certification, factor_set, max_anchors) {
    theta <- certification$theta
    total_delta <- certification$total$delta
    score <- abs(total_delta)
    if (all(!is.finite(score)) || !any(score > 0)) {
      score <- rep(1, length(total_delta))
    }
    Z <- tryCatch(.outer_whiten_theta(theta, w = rep(1 / nrow(theta), nrow(theta))), error = function(e) scale(theta))
    Z[!is.finite(Z)] <- 0
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

      if (isTRUE(control$enrich_all_locals)) {
        selected_local <- local_id
        reason <- rep("all", length(selected_local))
      } else {
        local_score <- abs(local_delta)
        ord <- order(local_score, decreasing = TRUE)
        total_abs <- sum(local_score)
        max_locals <- if (is.finite(control$enrichment_max_locals_per_anchor)) {
          as.integer(control$enrichment_max_locals_per_anchor)
        } else {
          length(ord)
        }
        cover_n <- if (total_abs > 0) {
          which(cumsum(local_score[ord]) >= as.numeric(control$enrichment_local_delta_coverage) * total_abs)[1L]
        } else {
          as.integer(control$enrichment_min_locals_per_anchor)
        }
        cover_n <- max(as.integer(control$enrichment_min_locals_per_anchor), cover_n %||% 0L)
        cover_n <- min(length(ord), cover_n, max_locals)
        selected <- ord[seq_len(max(1L, cover_n))]
        low_ess <- which(local_ess < as.numeric(control$enrichment_local_ess_threshold))
        high_delta <- which(local_score >= as.numeric(control$enrichment_local_delta_min))
        selected <- sort(unique(c(selected, low_ess, high_delta)))
        if (is.finite(control$enrichment_max_locals_per_anchor)) {
          selected <- selected[order(local_score[selected], decreasing = TRUE)]
          selected <- selected[seq_len(min(length(selected), as.integer(control$enrichment_max_locals_per_anchor)))]
        }
        selected_local <- local_id[selected]
        reason <- ifelse(local_ess[selected] < as.numeric(control$enrichment_local_ess_threshold), "low_ess", "delta")
      }

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

  enrich_from_plan <- function(local_objects_in, theta, plan, iter) {
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
        qmc_size = enrichment_particles,
        qmc_randomizations = enrichment_qmc_randomizations,
        seed = seed + 60000L * iter + 1000L * a,
        n_jobs = n_cores,
        local_n_cores = control$local_n_cores,
        label = sprintf("adaptive_iter%d_anchor%d", iter, a)
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

  for (iter in seq_len(max_iterations)) {
    design <- select_population_certification_design(
      population_fit = current_fit,
      factor_set = current_factor_set,
      population_model = model,
      size = control$calibration_size,
      seed = seed + 10000L * iter,
      stress_pool_size = control$design_stress_pool_size,
      stress_scale = control$design_stress_scale,
      stress_weight = control$design_stress_weight,
      stress_log_drop = control$design_stress_log_drop,
      n_cores = n_cores
    )
    if (isTRUE(verbose)) {
      cat(sprintf(
        "Adaptive DMIS iteration %d: certifying %d theta points across %d locals.\n",
        iter,
        nrow(design$theta),
        length(control$local_subset %||% data_list)
      ))
    }
    calibration <- certify_points(
      design$theta,
      factor_set = current_factor_set,
      seed_offset = 20000L * iter
    )
    certification_bank <- .append_certification_bank(
      certification_bank,
      calibration,
      weights = design$source_weights
    )
    enrichment_plan <- select_enrichment_plan(
      calibration,
      factor_set = current_factor_set,
      max_anchors = control$enrichment_anchors
    )
    if (isTRUE(verbose)) {
      cat(sprintf(
        "Adaptive DMIS iteration %d: enriching %d anchors and %d local components.\n",
        iter,
        nrow(enrichment_plan$anchors),
        nrow(enrichment_plan$locals)
      ))
    }
    enrichment <- enrich_from_plan(
      local_objects_in = local_objects_current,
      theta = calibration$theta,
      plan = enrichment_plan,
      iter = iter
    )
    local_objects_current <- enrichment$local_objects
    current_factor_set <- build_population_factor_set(
      local_objects_current,
      model,
      data_list = data_list,
      loglik_fn = loglik_fn,
      local_n_cores = control$local_n_cores
    )
    fit_args <- modifyList(
      list(
        factor_set = current_factor_set,
        N = 2000L,
        n_mcmc_moves = 3L,
        max_rounds = 80L,
        n_cores = as.integer(n_cores),
        seed = as.integer(seed + 30000L * iter),
        verbose = verbose
      ),
      outer_control
    )
    current_fit <- do.call(outer_population_smc, fit_args)

    validation_design <- select_population_certification_design(
      population_fit = current_fit,
      factor_set = current_factor_set,
      population_model = model,
      size = control$validation_size,
      seed = seed + 40000L * iter,
      stress_pool_size = max(control$validation_size * 8L, control$design_stress_pool_size %/% 2L),
      stress_scale = control$design_stress_scale,
      stress_weight = control$design_stress_weight,
      stress_log_drop = control$design_stress_log_drop,
      n_cores = n_cores
    )
    validation <- certify_points(
      validation_design$theta,
      factor_set = current_factor_set,
      seed_offset = 50000L * iter
    )
    residual <- validation$total$delta
    validation_summary <- data.frame(
      iteration = iter,
      rmse = sqrt(mean(residual^2)),
      median_abs = median(abs(residual)),
      max_abs = max(abs(residual)),
      mean = mean(residual),
      certified_mcse = sqrt(mean(validation$total$certified_mcse^2)),
      check.names = FALSE
    )
    if (isTRUE(verbose)) {
      cat(sprintf(
        "Adaptive DMIS validation %d: rmse=%.3f median_abs=%.3f max_abs=%.3f\n",
        iter,
        validation_summary$rmse,
        validation_summary$median_abs,
        validation_summary$max_abs
      ))
    }

    history[[iter]] <- list(
      iteration = iter,
      design = design,
      calibration = calibration,
      enrichment_plan = enrichment_plan,
      enrichment = enrichment,
      local_objects = local_objects_current,
      factor_set = current_factor_set,
      fit = current_fit,
      validation_design = validation_design,
      validation = validation,
      validation_residual = residual,
      validation_summary = validation_summary
    )

    certification_bank <- .append_certification_bank(
      certification_bank,
      validation,
      weights = validation_design$source_weights
    )
    passed <- validation_summary$rmse <= control$validation_rmse_tol &&
      validation_summary$median_abs <= control$validation_median_abs_tol
    if (passed) {
      selected_iteration <- iter
      validated <- TRUE
      break
    }
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
    mode = "adaptive_dmis",
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

.resolve_population_model <- function(population_model, pilot_population_model) {
  if (!is.null(population_model) && !is.null(pilot_population_model) && !identical(population_model, pilot_population_model)) {
    stop("Pass only population_model. pilot_population_model is kept as a compatibility alias.")
  }
  model <- population_model %||% pilot_population_model
  if (is.null(model)) {
    stop("prepare_reference_local_stage now requires population_model; pooled local-moment refinement was removed.")
  }
  normalize_population_model(model)
}

prepare_reference_local_stage <- function(data_list,
                                          loglik_fn,
                                          base_mu,
                                          base_Sigma,
                                          population_model = NULL,
                                          pilot_population_model = NULL,
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
                                          pilot_outer_control = list(N = 1000L, n_mcmc_moves = 2L, max_rounds = 50L),
                                          refresh_once = FALSE,
                                          refresh_outer_control = list(),
                                          verbose = TRUE,
                                          pilot_smc_control = list(max_rounds = 40L),
                                          full_smc_control = list()) {
  population_model <- .resolve_population_model(population_model, pilot_population_model)
  refresh_once <- isTRUE(refresh_once)

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

  run_full_local_pass <- function(reference_prior, seed_offset) {
    full_args <- modifyList(
      list(
        data_list = data_list,
        loglik_fn = loglik_fn,
        reference_prior = reference_prior,
        indices = seq_along(data_list),
        M = full_particles,
        n_jobs = n_jobs,
        local_n_cores = full_local_n_cores,
        base_seed = base_seed + as.integer(seed_offset),
        verbose = verbose
      ),
      full_smc_control
    )
    local_fits <- do.call(run_reference_local_smc, full_args)
    list(
      local_fits = local_fits,
      local_objects = build_local_reference_objects(local_fits, reference_prior = reference_prior)
    )
  }

  full_pass <- run_full_local_pass(
    reference_prior = refined_reference,
    seed_offset = 100000L
  )
  local_fits <- full_pass$local_fits
  local_objects <- full_pass$local_objects
  outer_fit <- NULL
  pre_refresh <- NULL
  refresh_merge_method <- NULL

  if (refresh_once) {
    initial_outer_fit <- fit_population_model_from_local_objects(
      local_objects = local_objects,
      population_model = population_model,
      outer_control = refresh_outer_control,
      n_cores = n_jobs,
      seed = base_seed + 150000L,
      verbose = verbose
    )

    refreshed_reference <- build_refined_reference_prior_from_population_fit(
      population_fit = initial_outer_fit,
      population_model = population_model,
      inflation = inflation,
      defensive_weight = defensive_weight,
      broad_reference = broad_reference,
      defensive_scale = defensive_scale,
      support_size = pilot_reference_support_size,
      support_seed = base_seed + 150001L
    )

    pre_refresh <- list(
      refined_reference = refined_reference,
      local_fits = local_fits,
      local_objects = local_objects,
      outer_fit = initial_outer_fit
    )

    refreshed_pass <- run_full_local_pass(
      reference_prior = refreshed_reference,
      seed_offset = 200000L
    )
    refined_reference <- refreshed_reference
    local_fits <- refreshed_pass$local_fits
    local_objects <- merge_local_reference_objects_dmis(
      primary_objects = pre_refresh$local_objects,
      secondary_objects = refreshed_pass$local_objects,
      mixture_weights = c(0.5, 0.5),
      labels = c("initial", "refreshed")
    )
    refresh_merge_method <- "dmis"

    outer_fit <- fit_population_model_from_local_objects(
      local_objects = local_objects,
      population_model = population_model,
      outer_control = refresh_outer_control,
      n_cores = n_jobs,
      seed = base_seed + 250000L,
      verbose = verbose
    )
  }

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
    local_fits = local_fits,
    local_objects = local_objects,
    outer_fit = outer_fit,
    pre_refresh = pre_refresh,
    refresh_merge_method = refresh_merge_method
  )
}
