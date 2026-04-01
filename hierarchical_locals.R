#!/usr/bin/env Rscript
# ============================================================================
# Hierarchical local-reference workflow
# - Broad reference prior construction
# - Stratified pilot-subset selection
# - Refined reference-prior construction from pilot local fits
# - Full local refits under the refined reference prior
# ============================================================================

if (!exists("%||%", mode = "function") ||
    !exists("weighted_cov", mode = "function") ||
    !exists("regularize_cov", mode = "function")) {
  source("smc_core.R")
}
if (!exists("normalize_reference_prior", mode = "function") ||
    !exists("make_broad_reference_prior", mode = "function")) {
  source("reference_priors.R")
}
if (!exists("enhanced_smc_elite", mode = "function")) {
  source("SMC_super_fast.R")
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
  list(
    mu = mu_bar,
    Sigma = Sigma_bar,
    local_summaries = summaries
  )
}

.ensure_population_refinement_helpers <- function() {
  if (!exists("build_population_factor_set", mode = "function") ||
      !exists("outer_population_smc", mode = "function") ||
      !exists("population_model_reference_components_from_theta", mode = "function")) {
    source("outer_population_smc.R")
  }
}

.sample_weighted_indices <- function(w, size, seed = NULL) {
  w <- normalize_particle_weights(w)
  n <- length(w)
  size <- as.integer(max(1L, min(size, n)))
  if (!is.null(seed)) set.seed(as.integer(seed))
  sample.int(n, size = size, replace = FALSE, prob = w)
}

fit_pilot_population_model <- function(pilot_fits,
                                       reference_prior,
                                       population_model,
                                       outer_control = list(),
                                       n_cores = 1L,
                                       seed = 123L,
                                       verbose = TRUE) {
  .ensure_population_refinement_helpers()
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

fit_population_model_from_local_objects <- function(local_objects,
                                                    population_model,
                                                    outer_control = list(),
                                                    n_cores = 1L,
                                                    seed = 123L,
                                                    verbose = TRUE) {
  .ensure_population_refinement_helpers()
  factor_set <- build_population_factor_set(local_objects, population_model)
  fit_args <- modifyList(
    list(
      factor_set = factor_set,
      N = 1000L,
      n_mcmc_moves = 2L,
      max_rounds = 50L,
      n_cores = as.integer(n_cores),
      seed = as.integer(seed),
      verbose = verbose
    ),
    outer_control
  )
  do.call(outer_population_smc, fit_args)
}

build_refined_reference_prior_from_population_fit <- function(population_fit,
                                                              population_model,
                                                              method = c("defensive_mixture", "broadened_gaussian"),
                                                              inflation = 1.5,
                                                              defensive_weight = 0.10,
                                                              broad_reference = NULL,
                                                              defensive_scale = 4,
                                                              support_size = 8L,
                                                              support_seed = 123L) {
  .ensure_population_refinement_helpers()
  method <- match.arg(method)
  model <- normalize_population_model(population_model)
  theta <- .as_hyper_matrix(population_fit$theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  w <- normalize_particle_weights(population_fit$w)
  comps <- population_model_reference_components_from_theta(model, theta)

  if (identical(method, "broadened_gaussian")) {
    mu_bar <- Reduce(
      `+`,
      Map(function(w_i, mu_i) w_i * mu_i, as.list(w), comps$component_means)
    )
    Sigma_bar <- Reduce(
      `+`,
      Map(
        function(w_i, mu_i, Sigma_i) {
          dm <- as.numeric(mu_i - mu_bar)
          w_i * (Sigma_i + tcrossprod(dm))
        },
        as.list(w),
        comps$component_means,
        comps$component_covs
      )
    )
    Sigma_bar <- regularize_cov(Sigma_bar, min_eig = 1e-8, cond_cap = 1e8)
    core_prior <- make_reference_prior_gaussian(
      mu = mu_bar,
      Sigma = Sigma_bar,
      scale = inflation,
      param_names = model$alpha_names,
      label = "pilot_population_gaussian"
    )
  } else {
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
  }

  if (identical(method, "broadened_gaussian")) {
    return(core_prior)
  }

  defensive_prior <- if (is.null(broad_reference)) {
    inflate_reference_prior(core_prior, scale = defensive_scale, label = "pilot_population_defensive_component")
  } else {
    normalize_reference_prior(reference_prior = broad_reference)
  }

  combine_reference_priors(
    priors = list(core_prior, defensive_prior),
    weights = c(1 - defensive_weight, defensive_weight),
    label = "pilot_population_defensive_mixture"
  )
}

build_refined_reference_prior <- function(pilot_fits,
                                          method = c("defensive_mixture", "broadened_gaussian"),
                                          inflation = 1.5,
                                          defensive_weight = 0.10,
                                          broad_reference = NULL,
                                          defensive_scale = 4) {
  method <- match.arg(method)
  pilot_summary <- aggregate_pilot_local_moments(pilot_fits)

  core_prior <- make_reference_prior_gaussian(
    mu = pilot_summary$mu,
    Sigma = pilot_summary$Sigma,
    scale = inflation,
    param_names = names(pilot_summary$mu),
    label = "pilot_centered_gaussian"
  )

  if (identical(method, "broadened_gaussian")) {
    return(core_prior)
  }

  defensive_prior <- if (is.null(broad_reference)) {
    inflate_reference_prior(core_prior, scale = defensive_scale, label = "pilot_defensive_component")
  } else {
    normalize_reference_prior(reference_prior = broad_reference)
  }

  combine_reference_priors(
    priors = list(core_prior, defensive_prior),
    weights = c(1 - defensive_weight, defensive_weight),
    label = "pilot_defensive_mixture"
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

prepare_reference_local_stage <- function(data_list,
                                          loglik_fn,
                                          base_mu,
                                          base_Sigma,
                                          pilot_size = 20L,
                                          broad_scale = 4,
                                          broad_defensive = FALSE,
                                          pilot_particles = 1000L,
                                          full_particles = 4000L,
                                          refined_method = c("defensive_mixture", "broadened_gaussian"),
                                          inflation = 1.5,
                                          defensive_weight = 0.10,
                                          defensive_scale = 4,
                                          n_jobs = 1L,
                                          pilot_local_n_cores = 1L,
                                          full_local_n_cores = 1L,
                                          base_seed = 123L,
                                          feature_fn = default_local_feature_vector,
                                          n_strata = NULL,
                                          pilot_population_model = NULL,
                                          pilot_reference_support_size = 8L,
                                          pilot_outer_control = list(),
                                          refresh_once = FALSE,
                                          refresh_outer_control = list(),
                                          verbose = TRUE,
                                          pilot_smc_control = list(
                                            n_mcmc_moves = 1L,
                                            max_rounds = 50L,
                                            G_mix = 8L,
                                            hist_mix_enable = FALSE,
                                            gss_enable = FALSE,
                                            da_enable = FALSE
                                          ),
                                          full_smc_control = list()) {
  refined_method <- match.arg(refined_method)
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

  pilot_population_fit <- NULL
  refined_reference <- if (is.null(pilot_population_model)) {
    build_refined_reference_prior(
      pilot_fits = pilot_fits,
      method = refined_method,
      inflation = inflation,
      defensive_weight = defensive_weight,
      broad_reference = broad_reference,
      defensive_scale = defensive_scale
    )
  } else {
    pilot_population_fit <- fit_pilot_population_model(
      pilot_fits = pilot_fits,
      reference_prior = broad_reference,
      population_model = pilot_population_model,
      outer_control = pilot_outer_control,
      n_cores = n_jobs,
      seed = base_seed + 50000L,
      verbose = verbose
    )
    build_refined_reference_prior_from_population_fit(
      population_fit = pilot_population_fit,
      population_model = pilot_population_model,
      method = refined_method,
      inflation = inflation,
      defensive_weight = defensive_weight,
      broad_reference = broad_reference,
      defensive_scale = defensive_scale,
      support_size = pilot_reference_support_size,
      support_seed = base_seed + 50001L
    )
  }

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
    if (is.null(pilot_population_model)) {
      stop("refresh_once requires pilot_population_model so the refreshed outer fit can be constructed.")
    }

    initial_outer_fit <- fit_population_model_from_local_objects(
      local_objects = local_objects,
      population_model = pilot_population_model,
      outer_control = refresh_outer_control,
      n_cores = n_jobs,
      seed = base_seed + 150000L,
      verbose = verbose
    )

    refreshed_reference <- build_refined_reference_prior_from_population_fit(
      population_fit = initial_outer_fit,
      population_model = pilot_population_model,
      method = refined_method,
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
      population_model = pilot_population_model,
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
      refinement_source = if (is.null(pilot_population_model)) "pooled_local_moments" else "pilot_population_fit"
    ),
    refined_reference = refined_reference,
    local_fits = local_fits,
    local_objects = local_objects,
    outer_fit = outer_fit,
    pre_refresh = pre_refresh,
    refresh_merge_method = refresh_merge_method
  )
}
