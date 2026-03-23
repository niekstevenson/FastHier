#!/usr/bin/env Rscript
# ============================================================================
# Collapsed hierarchical SMC variants.
#
# This file contains the subject-sequential collapsed outer samplers only:
# - collapsed_subject_smc(): add one subject at a time
# - collapsed_bridge_smc(): add one subject through within-subject bridges
#
# It expects the generic utilities from hierarchical_smc.R to already be
# available in the current session.
#
# Exact pseudo-marginal correctness is available only under the usual local
# assumptions: each subject-level evidence estimator must be nonnegative and
# unbiased for the corresponding marginal likelihood term. The outer code now
# carries the current evidence estimates as explicit auxiliary state; optional
# storage of full local fits is for diagnostics/reuse and is not required for
# the baseline pseudo-marginal argument.
# ============================================================================

if (!exists(".estimate_local_log_evidence_at_phi", mode = "function")) {
  stop("Source hierarchical_smc.R before hierarchical_collapsed_smc.R")
}

.collapsed_validate_subject_order <- function(data_list, subject_order = NULL) {
  n_subjects <- length(data_list)
  if (n_subjects <= 0L) stop("data_list must not be empty.")
  if (is.null(subject_order)) {
    return(seq_len(n_subjects))
  }
  subject_order <- as.integer(subject_order)
  if (!setequal(subject_order, seq_len(n_subjects))) {
    stop("subject_order must be a permutation of seq_along(data_list).")
  }
  subject_order
}

.collapsed_prior_dim <- function(prior) {
  seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  seed_old <- if (seed_exists) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
  on.exit({
    if (seed_exists) {
      assign(".Random.seed", seed_old, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  probe <- as.matrix(prior$rprior(1L))
  ncol(probe)
}

.normalize_bridge_schedule <- function(bridge_schedule = NULL, n_bridge_steps = 4L) {
  if (is.null(bridge_schedule)) {
    n_bridge_steps <- as.integer(n_bridge_steps)
    if (n_bridge_steps <= 0L) stop("n_bridge_steps must be positive when bridge_schedule is NULL.")
    bridge_schedule <- seq_len(n_bridge_steps) / n_bridge_steps
  }
  bridge_schedule <- sort(unique(as.numeric(bridge_schedule)))
  bridge_schedule <- bridge_schedule[is.finite(bridge_schedule) & bridge_schedule > 0 & bridge_schedule <= 1]
  if (!length(bridge_schedule)) {
    stop("bridge_schedule must contain at least one value in (0, 1].")
  }
  if (abs(tail(bridge_schedule, 1L) - 1.0) > 1e-12) {
    bridge_schedule <- c(bridge_schedule, 1.0)
  }
  bridge_schedule
}

.collapsed_stage_plan <- function(subject_order, bridge_schedule = 1.0) {
  bridge_schedule <- as.numeric(bridge_schedule)
  n_subjects <- length(subject_order)
  n_bridge <- length(bridge_schedule)
  data.frame(
    stage_index = seq_len(n_subjects * n_bridge),
    subject_position = rep(seq_len(n_subjects), each = n_bridge),
    subject = rep(as.integer(subject_order), each = n_bridge),
    lambda = rep(bridge_schedule, times = n_subjects),
    bridge_position = rep(seq_len(n_bridge), times = n_subjects),
    bridge_count = rep(n_bridge, n_subjects * n_bridge),
    subject_complete = rep(seq_len(n_bridge) == n_bridge, times = n_subjects),
    stringsAsFactors = FALSE
  )
}

.collapsed_identity_transform <- function() {
  list(
    name = "identity",
    forward = function(phi) as.numeric(phi),
    inverse = function(varphi) as.numeric(varphi),
    log_abs_det_jacobian = function(phi = NULL, varphi = NULL) 0
  )
}

.collapsed_normalize_proposal_transform <- function(proposal_transform, d_phi = NULL) {
  if (is.null(proposal_transform)) {
    proposal_transform <- .collapsed_identity_transform()
  }
  needed <- c("forward", "inverse", "log_abs_det_jacobian")
  if (!all(needed %in% names(proposal_transform))) {
    stop("proposal_transform must provide forward, inverse, and log_abs_det_jacobian.")
  }
  probe_phi <- rep(0, d_phi %||% 1L)
  probe_varphi <- as.numeric(proposal_transform$forward(probe_phi))
  if (!is.numeric(probe_varphi)) {
    stop("proposal_transform$forward must return a numeric vector.")
  }
  list(
    name = proposal_transform$name %||% "custom_transform",
    forward = proposal_transform$forward,
    inverse = proposal_transform$inverse,
    log_abs_det_jacobian = proposal_transform$log_abs_det_jacobian
  )
}

.collapsed_transform_matrix <- function(phi, proposal_transform) {
  phi <- as.matrix(phi)
  out <- t(vapply(
    seq_len(nrow(phi)),
    function(i) as.numeric(proposal_transform$forward(phi[i, , drop = TRUE])),
    numeric(ncol(phi))
  ))
  colnames(out) <- colnames(phi)
  out
}

.collapsed_make_evidence_cache <- function(enabled = TRUE,
                                           digits = 8L,
                                           keep_local_fit = FALSE,
                                           stochastic_by_seed = TRUE) {
  cache_env <- new.env(parent = emptyenv())
  index_env <- new.env(parent = emptyenv())
  surrogate_env <- new.env(parent = emptyenv())
  stats_env <- new.env(parent = emptyenv())
  stats_env$enabled <- isTRUE(enabled)
  stats_env$digits <- as.integer(digits)
  stats_env$keep_local_fit <- isTRUE(keep_local_fit)
  stats_env$stochastic_by_seed <- isTRUE(stochastic_by_seed)
  stats_env$hits <- 0L
  stats_env$misses <- 0L
  stats_env$stores <- 0L
  stats_env$size <- 0L
  list(env = cache_env, index = index_env, surrogate = surrogate_env, stats = stats_env)
}

.collapsed_cache_tag <- function(subject, lambda, digits = 8L) {
  sprintf(
    "s=%d|l=%s",
    as.integer(subject),
    format(round(as.numeric(lambda), digits = digits), scientific = FALSE, trim = TRUE)
  )
}

.collapsed_cache_key <- function(subject,
                                 lambda,
                                 phi,
                                 seed = NULL,
                                 aux_signature = NULL,
                                 digits = 8L,
                                 stochastic_by_seed = TRUE) {
  phi_txt <- paste(format(round(as.numeric(phi), digits = digits), scientific = FALSE, trim = TRUE), collapse = ",")
  aux_signature <- aux_signature %||% .collapsed_auxiliary_signature(seed = seed)
  seed_txt <- if (isTRUE(stochastic_by_seed) && !is.null(aux_signature) && nzchar(aux_signature)) {
    sprintf("|aux=%s", aux_signature)
  } else {
    ""
  }
  sprintf(
    "%s|phi=%s%s",
    .collapsed_cache_tag(subject = subject, lambda = lambda, digits = digits),
    phi_txt,
    seed_txt
  )
}

.collapsed_cache_get <- function(cache_store, subject, lambda, phi, seed = NULL, aux_signature = NULL) {
  if (is.null(cache_store) || !isTRUE(cache_store$stats$enabled)) {
    return(NULL)
  }
  key <- .collapsed_cache_key(
    subject = subject,
    lambda = lambda,
    phi = phi,
    seed = seed,
    aux_signature = aux_signature,
    digits = cache_store$stats$digits,
    stochastic_by_seed = cache_store$stats$stochastic_by_seed
  )
  if (exists(key, envir = cache_store$env, inherits = FALSE)) {
    cache_store$stats$hits <- cache_store$stats$hits + 1L
    out <- get(key, envir = cache_store$env, inherits = FALSE)
    out$cache_key <- key
    return(out)
  }
  cache_store$stats$misses <- cache_store$stats$misses + 1L
  NULL
}

.collapsed_cache_put <- function(cache_store, subject, lambda, phi, seed = NULL, aux_signature = NULL, value) {
  if (is.null(cache_store) || !isTRUE(cache_store$stats$enabled)) {
    return(invisible(NULL))
  }
  key <- .collapsed_cache_key(
    subject = subject,
    lambda = lambda,
    phi = phi,
    seed = seed,
    aux_signature = aux_signature,
    digits = cache_store$stats$digits,
    stochastic_by_seed = cache_store$stats$stochastic_by_seed
  )
  tag <- .collapsed_cache_tag(subject = subject, lambda = lambda, digits = cache_store$stats$digits)
  is_new <- !exists(key, envir = cache_store$env, inherits = FALSE)
  if (is_new) {
    cache_store$stats$stores <- cache_store$stats$stores + 1L
    cache_store$stats$size <- cache_store$stats$size + 1L
    recs <- if (exists(tag, envir = cache_store$index, inherits = FALSE)) {
      get(tag, envir = cache_store$index, inherits = FALSE)
    } else {
      list(
        key = character(0),
        phi = list(),
        log_evidence = numeric(0),
        seed = integer(0),
        aux_signature = character(0),
        has_local_fit = logical(0),
        version = 0L
      )
    }
    recs$key <- c(recs$key, key)
    recs$phi[[length(recs$phi) + 1L]] <- as.numeric(phi)
    recs$log_evidence <- c(recs$log_evidence, as.numeric(value$log_evidence))
    recs$seed <- c(recs$seed, as.integer(seed %||% NA_integer_))
    recs$aux_signature <- c(recs$aux_signature, as.character(aux_signature %||% ""))
    recs$has_local_fit <- c(recs$has_local_fit, !is.null(value$local_fit))
    recs$version <- as.integer(recs$version + 1L)
    assign(tag, recs, envir = cache_store$index)
    if (exists(tag, envir = cache_store$surrogate, inherits = FALSE)) {
      rm(list = tag, envir = cache_store$surrogate)
    }
  }
  value$cache_key <- key
  assign(key, value, envir = cache_store$env)
  invisible(NULL)
}

.collapsed_cache_records <- function(cache_store,
                                     subject,
                                     lambda,
                                     require_local_fit = FALSE) {
  if (is.null(cache_store) || !isTRUE(cache_store$stats$enabled)) {
    return(list())
  }
  tag <- .collapsed_cache_tag(subject = subject, lambda = lambda, digits = cache_store$stats$digits)
  if (!exists(tag, envir = cache_store$index, inherits = FALSE)) {
    return(list())
  }
  recs <- get(tag, envir = cache_store$index, inherits = FALSE)
  out <- vector("list", length(recs$key))
  keep <- logical(length(out))
  for (ii in seq_along(recs$key)) {
    val <- get(recs$key[ii], envir = cache_store$env, inherits = FALSE)
    if (isTRUE(require_local_fit) && is.null(val$local_fit)) {
      next
    }
    val$cache_key <- recs$key[ii]
    val$phi <- recs$phi[[ii]]
    val$seed <- recs$seed[ii]
    val$aux_signature <- recs$aux_signature[ii] %||% ""
    out[[ii]] <- val
    keep[ii] <- TRUE
  }
  out[keep]
}

.collapsed_cache_stats <- function(cache_store) {
  if (is.null(cache_store)) {
    return(list(
      enabled = FALSE,
      keep_local_fit = FALSE,
      stochastic_by_seed = TRUE,
      hits = 0L,
      misses = 0L,
      stores = 0L,
      size = 0L
    ))
  }
  list(
    enabled = isTRUE(cache_store$stats$enabled),
    keep_local_fit = isTRUE(cache_store$stats$keep_local_fit),
    stochastic_by_seed = isTRUE(cache_store$stats$stochastic_by_seed),
    hits = as.integer(cache_store$stats$hits),
    misses = as.integer(cache_store$stats$misses),
    stores = as.integer(cache_store$stats$stores),
    size = as.integer(cache_store$stats$size)
  )
}

.collapsed_make_aux_entry <- function(phi,
                                      subject,
                                      lambda,
                                      log_evidence,
                                      seed = NA_integer_,
                                      aux_latent = NULL,
                                      auxiliary_state = NULL,
                                      local_fit = NULL,
                                      cache_key = NULL,
                                      aux_signature = NULL,
                                      checkpoint_lambdas = numeric(0),
                                      checkpoint_log_evidence = numeric(0),
                                      checkpoint_log_increments = numeric(0),
                                      source = "fresh",
                                      warm_start_key = NULL) {
  seed_plan <- auxiliary_state$seed_plan %||% NULL
  list(
    phi = as.numeric(phi),
    subject = as.integer(subject),
    lambda = as.numeric(lambda),
    log_evidence = as.numeric(log_evidence),
    seed = as.integer(seed %||% NA_integer_),
    aux_latent = if (is.null(aux_latent)) NULL else as.numeric(aux_latent),
    auxiliary_state = auxiliary_state,
    seed_plan = seed_plan,
    local_fit = local_fit,
    cache_key = cache_key,
    aux_signature = as.character(aux_signature %||% .collapsed_auxiliary_signature(
      seed = seed,
      seed_plan = seed_plan,
      auxiliary_state = auxiliary_state
    )),
    checkpoint_lambdas = as.numeric(checkpoint_lambdas %||% numeric(0)),
    checkpoint_log_evidence = as.numeric(checkpoint_log_evidence %||% numeric(0)),
    checkpoint_log_increments = as.numeric(checkpoint_log_increments %||% numeric(0)),
    source = as.character(source),
    warm_start_key = warm_start_key
  )
}

.collapsed_trim_aux_entry <- function(aux_entry, keep_local_fit = FALSE) {
  if (is.null(aux_entry)) {
    return(NULL)
  }
  aux_entry$local_fit <- if (isTRUE(keep_local_fit)) aux_entry$local_fit else NULL
  aux_entry
}

.collapsed_new_subject_aux <- function(N, n_subjects, store_local_fits = FALSE) {
  replicate(as.integer(n_subjects), replicate(as.integer(N), NULL, simplify = FALSE), simplify = FALSE)
}

.collapsed_resample_subject_aux <- function(subject_aux, idx) {
  if (is.null(subject_aux)) {
    return(NULL)
  }
  lapply(subject_aux, function(x) {
    if (is.null(x)) {
      NULL
    } else {
      x[idx]
    }
  })
}

.collapsed_set_subject_aux <- function(subject_aux, subject, aux_entries) {
  if (is.null(subject_aux)) {
    return(NULL)
  }
  subject_aux[[as.integer(subject)]] <- aux_entries
  subject_aux
}

.collapsed_subject_aux_row <- function(subject_aux, j) {
  if (is.null(subject_aux)) {
    return(NULL)
  }
  lapply(subject_aux, function(x) {
    if (is.null(x) || length(x) < j) NULL else x[[j]]
  })
}

.collapsed_rebuild_subject_aux <- function(prop_list, n_subjects) {
  if (!length(prop_list) || is.null(prop_list[[1]]$subject_aux_row)) {
    return(NULL)
  }
  lapply(seq_len(as.integer(n_subjects)), function(s) {
    lapply(prop_list, function(x) x$subject_aux_row[[s]])
  })
}

.collapsed_random_int_seeds <- function(n, seed_base) {
  set.seed(as.integer(seed_base))
  sample.int(.Machine$integer.max - 1L, size = as.integer(n), replace = TRUE)
}

.collapsed_local_seed_plan_shape <- function(local_smc_control = list()) {
  max_rounds <- as.integer(local_smc_control$max_rounds %||% 200L)
  if (!is.finite(max_rounds) || max_rounds <= 0L) {
    max_rounds <- 200L
  }
  list(
    max_rounds = max_rounds,
    header_names = c("global", "init", "warm_subsample", "warm_elite_fit"),
    round_names = c("da_calib", "pre_move", "resample", "jitter1", "jitter2", "elite_fit", "main_move")
  )
}

.collapsed_auxiliary_dim <- function(local_smc_control = list()) {
  shape <- .collapsed_local_seed_plan_shape(local_smc_control)
  length(shape$header_names) + shape$max_rounds * length(shape$round_names)
}

.collapsed_draw_auxiliary_latents <- function(n, seed_base, local_smc_control = list()) {
  set.seed(as.integer(seed_base))
  dim_aux <- .collapsed_auxiliary_dim(local_smc_control)
  matrix(stats::rnorm(as.integer(n) * dim_aux), nrow = as.integer(n), ncol = dim_aux)
}

.collapsed_seed_plan_from_latent <- function(aux_latent, local_smc_control = list()) {
  shape <- .collapsed_local_seed_plan_shape(local_smc_control)
  z <- as.numeric(aux_latent)
  need <- length(shape$header_names) + shape$max_rounds * length(shape$round_names)
  if (length(z) != need) {
    stop("aux_latent has wrong length for the local seed plan shape.")
  }
  u <- pmin(pmax(stats::pnorm(z), 1e-12), 1 - 1e-12)
  seeds <- as.integer(floor(u * (.Machine$integer.max - 2L)) + 1L)
  ptr <- 0L
  take <- function(n) {
    idx <- seq.int(ptr + 1L, ptr + as.integer(n))
    ptr <<- ptr + as.integer(n)
    seeds[idx]
  }
  header_vals <- take(length(shape$header_names))
  names(header_vals) <- shape$header_names
  round_vals <- matrix(
    take(shape$max_rounds * length(shape$round_names)),
    nrow = shape$max_rounds,
    ncol = length(shape$round_names),
    byrow = TRUE,
    dimnames = list(NULL, shape$round_names)
  )
  list(
    latent = z,
    seeds = seeds,
    global = as.integer(header_vals[["global"]]),
    init = as.integer(header_vals[["init"]]),
    warm_subsample = as.integer(header_vals[["warm_subsample"]]),
    warm_elite_fit = as.integer(header_vals[["warm_elite_fit"]]),
    round = round_vals
  )
}

.collapsed_seed_plan_signature <- function(seed_plan) {
  seeds <- as.integer(seed_plan$seeds %||% integer(0))
  if (!length(seeds)) {
    return("none")
  }
  paste(seeds, collapse = ",")
}

.collapsed_auxiliary_signature <- function(seed = NULL, seed_plan = NULL, auxiliary_state = NULL) {
  if (!is.null(auxiliary_state$seed_plan)) {
    return(.collapsed_seed_plan_signature(auxiliary_state$seed_plan))
  }
  if (!is.null(seed_plan)) {
    return(.collapsed_seed_plan_signature(seed_plan))
  }
  if (!is.null(seed) && is.finite(seed)) {
    return(sprintf("seed=%d", as.integer(seed)))
  }
  "none"
}

.collapsed_propose_auxiliary_state <- function(current_aux_entry,
                                               mode = c("independent", "sticky", "frozen", "correlated"),
                                               stickiness = 0.9,
                                               correlation = 0.99,
                                               local_smc_control = list()) {
  mode <- match.arg(mode)
  dim_aux <- .collapsed_auxiliary_dim(local_smc_control)
  current_latent <- as.numeric(current_aux_entry$aux_latent %||% numeric(0))
  if (length(current_latent) != dim_aux) {
    current_latent <- NULL
  }
  latent <- switch(
    mode,
    independent = stats::rnorm(dim_aux),
    sticky = if (!is.null(current_latent) && stats::runif(1L) < as.numeric(stickiness)) current_latent else stats::rnorm(dim_aux),
    frozen = if (!is.null(current_latent)) current_latent else stats::rnorm(dim_aux),
    correlated = {
      rho <- min(max(as.numeric(correlation), 0), 0.999999)
      if (is.null(current_latent)) {
        stats::rnorm(dim_aux)
      } else {
        rho * current_latent + sqrt(1 - rho^2) * stats::rnorm(dim_aux)
      }
    }
  )
  seed_plan <- .collapsed_seed_plan_from_latent(latent, local_smc_control = local_smc_control)
  list(
    aux_latent = latent,
    seed_plan = seed_plan,
    signature = .collapsed_seed_plan_signature(seed_plan)
  )
}

.collapsed_aux_entry_matches_phi <- function(aux_entry, phi, tol = 1e-12) {
  if (is.null(aux_entry) || is.null(aux_entry$phi)) {
    return(FALSE)
  }
  aux_phi <- as.numeric(aux_entry$phi)
  phi <- as.numeric(phi)
  length(aux_phi) == length(phi) && all(abs(aux_phi - phi) <= tol)
}

.collapsed_aux_entry_value_at_lambda <- function(aux_entry, lambda) {
  if (is.null(aux_entry)) {
    return(NULL)
  }
  lambda <- as.numeric(lambda)
  cp <- as.numeric(aux_entry$checkpoint_lambdas %||% numeric(0))
  logz_cp <- as.numeric(aux_entry$checkpoint_log_evidence %||% numeric(0))
  logz_inc <- as.numeric(aux_entry$checkpoint_log_increments %||% numeric(0))
  if (length(cp)) {
    hit <- which(abs(cp - lambda) <= 1e-12)
    if (length(hit)) {
      idx <- hit[1L]
      return(list(
        log_evidence = as.numeric(logz_cp[idx]),
        log_increment = as.numeric(logz_inc[idx] %||% NA_real_)
      ))
    }
  }
  if (isTRUE(abs(as.numeric(aux_entry$lambda %||% NA_real_) - lambda) <= 1e-12) &&
      is.finite(aux_entry$log_evidence %||% NA_real_)) {
    return(list(
      log_evidence = as.numeric(aux_entry$log_evidence),
      log_increment = NA_real_
    ))
  }
  NULL
}

.collapsed_lookup_bridge_schedule <- function(bridge_schedule_lookup, subject, default = NULL) {
  if (is.null(bridge_schedule_lookup)) {
    return(default)
  }
  key <- as.character(as.integer(subject))
  if (!is.null(names(bridge_schedule_lookup)) && key %in% names(bridge_schedule_lookup)) {
    return(bridge_schedule_lookup[[key]])
  }
  idx <- as.integer(subject)
  if (length(bridge_schedule_lookup) >= idx) {
    return(bridge_schedule_lookup[[idx]] %||% default)
  }
  default
}

.collapsed_nearest_cached_local_fit <- function(cache_store,
                                                subject,
                                                lambda,
                                                phi,
                                                proposal_transform,
                                                max_candidates = 64L) {
  recs <- .collapsed_cache_records(
    cache_store = cache_store,
    subject = subject,
    lambda = lambda,
    require_local_fit = TRUE
  )
  if (!length(recs)) {
    return(NULL)
  }
  phi_tgt <- as.numeric(proposal_transform$forward(phi))
  dists <- vapply(
    recs,
    function(rec) {
      rec_phi <- rec$phi %||% (rec$local_fit$phi_anchor %||% NULL)
      if (is.null(rec_phi)) {
        return(Inf)
      }
      diff <- as.numeric(proposal_transform$forward(rec_phi)) - phi_tgt
      sum(diff * diff)
    },
    numeric(1)
  )
  ord <- order(dists, decreasing = FALSE, na.last = TRUE)
  ord <- ord[seq_len(min(length(ord), as.integer(max_candidates)))]
  for (ii in ord) {
    if (!is.null(recs[[ii]]$local_fit)) {
      return(recs[[ii]])
    }
  }
  NULL
}

.collapsed_pick_warm_start <- function(subject_aux_entry,
                                       cache_store,
                                       subject,
                                       lambda,
                                       phi,
                                       proposal_transform,
                                       warm_start_mode = c("none", "nearest", "current_then_nearest")) {
  warm_start_mode <- match.arg(warm_start_mode)
  if (identical(warm_start_mode, "none")) {
    return(NULL)
  }
  if (identical(warm_start_mode, "current_then_nearest") &&
      !is.null(subject_aux_entry) &&
      !is.null(subject_aux_entry$local_fit)) {
    return(list(local_fit = subject_aux_entry$local_fit, cache_key = subject_aux_entry$cache_key %||% NA_character_))
  }
  rec <- .collapsed_nearest_cached_local_fit(
    cache_store = cache_store,
    subject = subject,
    lambda = lambda,
    phi = phi,
    proposal_transform = proposal_transform
  )
  if (is.null(rec)) {
    return(NULL)
  }
  list(local_fit = rec$local_fit, cache_key = rec$cache_key %||% NA_character_)
}

.collapsed_surrogate_feature_matrix <- function(phi_mat, degree = 2L) {
  phi_mat <- as.matrix(phi_mat)
  cols <- list(`(Intercept)` = rep(1, nrow(phi_mat)))
  for (jj in seq_len(ncol(phi_mat))) {
    cols[[sprintf("x%d", jj)]] <- phi_mat[, jj]
  }
  if (as.integer(degree) >= 2L) {
    for (jj in seq_len(ncol(phi_mat))) {
      cols[[sprintf("x%d_sq", jj)]] <- phi_mat[, jj]^2
    }
  }
  as.matrix(as.data.frame(cols, check.names = FALSE))
}

.collapsed_build_surrogate <- function(cache_store,
                                       subject,
                                       lambda,
                                       proposal_transform,
                                       min_points = NULL,
                                       ridge = 1e-4,
                                       degree = 2L,
                                       neighbors = NULL) {
  recs <- .collapsed_cache_records(cache_store, subject = subject, lambda = lambda, require_local_fit = FALSE)
  if (!length(recs)) {
    return(NULL)
  }
  tag <- .collapsed_cache_tag(subject = subject, lambda = lambda, digits = cache_store$stats$digits)
  rec_index <- get(tag, envir = cache_store$index, inherits = FALSE)
  cached_model <- if (exists(tag, envir = cache_store$surrogate, inherits = FALSE)) {
    get(tag, envir = cache_store$surrogate, inherits = FALSE)
  } else {
    NULL
  }
  if (!is.null(cached_model) && identical(as.integer(cached_model$version), as.integer(rec_index$version))) {
    return(cached_model$model)
  }
  phi_mat <- t(vapply(recs, function(rec) as.numeric(rec$phi), numeric(length(recs[[1L]]$phi))))
  varphi <- .collapsed_transform_matrix(phi_mat, proposal_transform)
  y <- vapply(recs, function(rec) as.numeric(rec$log_evidence), numeric(1))
  p_basic <- ncol(varphi)
  min_points <- as.integer(min_points %||% max(6L, 2L * p_basic + 4L))
  neighbors <- as.integer(neighbors %||% max(8L, min(length(y), 3L * p_basic + 6L)))
  X <- .collapsed_surrogate_feature_matrix(varphi, degree = degree)
  beta <- NULL
  if (nrow(X) >= min_points && nrow(X) >= ncol(X)) {
    XtX <- crossprod(X)
    penalty <- diag(as.numeric(ridge), ncol(X))
    penalty[1L, 1L] <- 0
    beta <- tryCatch(
      solve(XtX + penalty, crossprod(X, y)),
      error = function(e) NULL
    )
  }
  model <- list(
    subject = as.integer(subject),
    lambda = as.numeric(lambda),
    version = as.integer(rec_index$version),
    proposal_transform = proposal_transform,
    varphi = varphi,
    y = y,
    beta = beta,
    degree = as.integer(degree),
    neighbors = max(1L, min(as.integer(neighbors), nrow(varphi)))
  )
  assign(tag, list(version = as.integer(rec_index$version), model = model), envir = cache_store$surrogate)
  model
}

.collapsed_predict_surrogate <- function(model, phi) {
  if (is.null(model)) {
    return(NA_real_)
  }
  varphi <- matrix(as.numeric(model$proposal_transform$forward(phi)), nrow = 1L)
  d2 <- rowSums((model$varphi - matrix(varphi, nrow(model$varphi), ncol(model$varphi), byrow = TRUE))^2)
  ord <- order(d2, decreasing = FALSE, na.last = TRUE)
  k <- min(length(ord), as.integer(model$neighbors))
  ord <- ord[seq_len(k)]
  local_pred <- if (!length(ord)) {
    NA_real_
  } else if (all(d2[ord] <= 1e-16)) {
    mean(model$y[ord])
  } else {
    ww <- 1 / pmax(sqrt(d2[ord]), 1e-8)
    sum(ww * model$y[ord]) / sum(ww)
  }
  if (is.null(model$beta)) {
    return(as.numeric(local_pred))
  }
  X <- .collapsed_surrogate_feature_matrix(varphi, degree = model$degree)
  global_pred <- as.numeric(drop(X %*% model$beta))
  if (!is.finite(local_pred)) {
    return(global_pred)
  }
  0.5 * global_pred + 0.5 * local_pred
}

.collapsed_estimate_log_evidence <- function(phi,
                                             subject,
                                             lambda,
                                             data_i,
                                             loglik_fn,
                                             gaussian_map_fn,
                                             d_theta,
                                             local_particles,
                                             local_log_evidence_fn = NULL,
                                             theta_names = NULL,
                                             local_n_cores = 1L,
                                             seed = NULL,
                                             seed_plan = NULL,
                                             aux_latent = NULL,
                                             checkpoint_lambdas = NULL,
                                             auxiliary_state = NULL,
                                             verbose = FALSE,
                                             warm_start_fit = NULL,
                                             warm_start_key = NULL,
                                             local_smc_control = list(),
                                             cache_store = NULL) {
  aux_signature <- .collapsed_auxiliary_signature(
    seed = seed,
    seed_plan = seed_plan,
    auxiliary_state = auxiliary_state
  )
  cached <- .collapsed_cache_get(
    cache_store = cache_store,
    subject = subject,
    lambda = lambda,
    phi = phi,
    seed = seed,
    aux_signature = aux_signature
  )
  if (!is.null(cached)) {
    return(cached)
  }
  out <- .estimate_local_log_evidence_at_phi(
    phi = phi,
    data_i = data_i,
    loglik_fn = loglik_fn,
    gaussian_map_fn = gaussian_map_fn,
    d_theta = d_theta,
    local_particles = local_particles,
    local_log_evidence_fn = local_log_evidence_fn,
    lambda = lambda,
    theta_names = theta_names,
    local_n_cores = local_n_cores,
    seed = seed,
    seed_plan = seed_plan,
    checkpoint_lambdas = checkpoint_lambdas,
    auxiliary_state = auxiliary_state,
    verbose = verbose,
    warm_start_fit = warm_start_fit,
    local_smc_control = local_smc_control
  )
  aux_entry <- .collapsed_make_aux_entry(
    phi = phi,
    subject = subject,
    lambda = lambda,
    log_evidence = out$log_evidence,
    seed = seed,
    aux_latent = aux_latent,
    auxiliary_state = out$auxiliary_state %||% auxiliary_state,
    local_fit = out$local_fit,
    source = if (is.null(warm_start_fit)) "fresh" else "warm_start",
    warm_start_key = warm_start_key,
    aux_signature = aux_signature,
    checkpoint_lambdas = out$checkpoint_lambdas,
    checkpoint_log_evidence = out$checkpoint_log_evidence,
    checkpoint_log_increments = out$checkpoint_log_increments
  )
  trimmed <- .collapsed_trim_aux_entry(
    aux_entry = aux_entry,
    keep_local_fit = !is.null(cache_store) && isTRUE(cache_store$stats$keep_local_fit)
  )
  .collapsed_cache_put(
    cache_store = cache_store,
    subject = subject,
    lambda = lambda,
    phi = phi,
    seed = seed,
    aux_signature = aux_signature,
    value = trimmed
  )
  trimmed
}

.collapsed_initialize_state <- function(prior, N, base_seed, n_subjects, store_local_fits = FALSE) {
  set.seed(as.integer(base_seed))
  phi <- as.matrix(prior$rprior(as.integer(N)))
  if (nrow(phi) != as.integer(N)) stop("prior$rprior returned wrong number of particles.")
  list(
    phi = phi,
    w = rep(1 / nrow(phi), nrow(phi)),
    logprior = apply(phi, 1L, prior$lprior),
    subject_loglik = matrix(NA_real_, nrow = nrow(phi), ncol = as.integer(n_subjects)),
    subject_aux = .collapsed_new_subject_aux(N = nrow(phi), n_subjects = n_subjects, store_local_fits = store_local_fits),
    completed_subjects = integer(0),
    current_subject = NULL,
    log_evidence = 0
  )
}

.collapsed_finalize_current_subject <- function(state) {
  if (is.null(state$current_subject)) {
    return(state)
  }
  state$completed_subjects <- c(state$completed_subjects, as.integer(state$current_subject))
  state$current_subject <- NULL
  state
}

.collapsed_enter_subject <- function(state, subject) {
  subject <- as.integer(subject)
  if (!is.null(state$current_subject) && identical(as.integer(state$current_subject), subject)) {
    return(state)
  }
  state <- .collapsed_finalize_current_subject(state)
  state$current_subject <- subject
  state
}

.collapsed_resample_state <- function(state) {
  idx <- .resample_indices(state$w)
  state$phi <- state$phi[idx, , drop = FALSE]
  state$w <- rep(1 / nrow(state$phi), nrow(state$phi))
  state$logprior <- state$logprior[idx]
  state$subject_loglik <- state$subject_loglik[idx, , drop = FALSE]
  state$subject_aux <- .collapsed_resample_subject_aux(state$subject_aux, idx)
  state
}

.collapsed_should_rejuvenate <- function(rejuvenation,
                                         rejuvenate_every,
                                         rejuvenate_after_resample,
                                         resampled,
                                         stage_index) {
  if (identical(rejuvenation, "none")) {
    return(FALSE)
  }
  every_hit <- is.finite(rejuvenate_every) &&
    rejuvenate_every > 0 &&
    (stage_index %% as.integer(rejuvenate_every) == 0L)
  isTRUE(rejuvenate_after_resample && resampled) || every_hit
}

.collapsed_verbose_stage_header <- function(stage_row, n_subjects, verbose) {
  if (!isTRUE(verbose)) return(invisible(NULL))
  if (stage_row$bridge_count <= 1L) {
    cat(sprintf(
      "Subject stage %d/%d: subject %d\n",
      stage_row$subject_position,
      n_subjects,
      stage_row$subject
    ))
  } else {
    if (stage_row$bridge_position == 1L) {
      cat(sprintf(
        "Subject %d/%d: subject %d\n",
        stage_row$subject_position,
        n_subjects,
        stage_row$subject
      ))
    }
    cat(sprintf(
      "  Bridge stage %d/%d: lambda %.4f\n",
      stage_row$bridge_position,
      stage_row$bridge_count,
      stage_row$lambda
    ))
  }
}

.collapsed_verbose_resample <- function(stage_row, ess_now, verbose) {
  if (!isTRUE(verbose)) return(invisible(NULL))
  prefix <- if (stage_row$bridge_count <= 1L) "  " else "    "
  cat(sprintf("%sresampled at ESS=%.1f\n", prefix, ess_now))
}

.collapsed_verbose_rejuvenation <- function(stage_row, accept_rate, verbose) {
  if (!isTRUE(verbose)) return(invisible(NULL))
  prefix <- if (stage_row$bridge_count <= 1L) "  " else "    "
  cat(sprintf("%srejuvenation accept=%.3f\n", prefix, accept_rate))
}

.collapsed_verbose_stage_footer <- function(stage_row, ess_now, logz_inc, logz_total, verbose) {
  if (!isTRUE(verbose)) return(invisible(NULL))
  prefix <- if (stage_row$bridge_count <= 1L) "  " else "    "
  cat(sprintf("%sESS=%.1f | logZ+=%.4f -> %.4f\n", prefix, ess_now, logz_inc, logz_total))
}

.collapsed_rejuvenation_cov <- function(phi, w, diag_floor = 1e-8) {
  phi <- as.matrix(phi)
  d_phi <- ncol(phi)
  S <- tryCatch(weighted_cov(phi, w), error = function(e) NULL)
  if (is.null(S) || !all(dim(S) == c(d_phi, d_phi)) || any(!is.finite(S))) {
    S <- tryCatch(stats::cov(phi), error = function(e) NULL)
  }
  if (is.null(S) || !all(dim(S) == c(d_phi, d_phi)) || any(!is.finite(S))) {
    vars <- apply(phi, 2L, stats::var)
    vars[!is.finite(vars)] <- 1
    S <- diag(pmax(vars, diag_floor), d_phi)
  }
  S <- (S + t(S)) / 2
  if (exists(".regularize_cov_safe", mode = "function")) {
    S <- .regularize_cov_safe(S, jitter = diag_floor)
  } else {
    S <- tryCatch(as.matrix(Matrix::nearPD(S, conv.tol = 1e-7)$mat), error = function(e) S)
    S <- S + diag(diag_floor, d_phi)
  }
  diag(S) <- pmax(diag(S), diag_floor)
  S
}

.collapsed_stage_increment <- function(phi,
                                       subject,
                                       data_i,
                                       loglik_fn,
                                       gaussian_map_fn,
                                       d_theta,
                                       local_particles,
                                       local_log_evidence_fn = NULL,
                                       lambda = 1.0,
                                       theta_names = NULL,
                                       local_n_cores = 1L,
                                       outer_n_cores = 1L,
                                       seed_base = 1L,
                                       verbose = FALSE,
                                       subject_aux = NULL,
                                       proposal_transform = .collapsed_identity_transform(),
                                       warm_start_mode = c("none", "nearest", "current_then_nearest"),
                                       full_bridge_schedule = NULL,
                                       local_smc_control = list(),
                                       cache_store = NULL,
                                       store_local_fits = FALSE) {
  phi <- as.matrix(phi)
  N <- nrow(phi)
  warm_start_mode <- match.arg(warm_start_mode)
  full_bridge_schedule <- sort(unique(as.numeric(full_bridge_schedule %||% lambda)))
  full_bridge_schedule <- full_bridge_schedule[
    is.finite(full_bridge_schedule) &
      full_bridge_schedule > 0 &
      full_bridge_schedule <= 1
  ]
  if (!length(full_bridge_schedule)) {
    full_bridge_schedule <- as.numeric(lambda)
  }
  use_checkpoint_path <- length(full_bridge_schedule) > 1L
  aux_latents <- .collapsed_draw_auxiliary_latents(N, seed_base = seed_base, local_smc_control = local_smc_control)
  out <- vector("list", N)
  log_mhat <- rep(NA_real_, N)
  miss_idx <- integer(0)
  for (j in seq_len(N)) {
    current_aux <- NULL
    if (!is.null(subject_aux) &&
        length(subject_aux) >= as.integer(subject) &&
        length(subject_aux[[as.integer(subject)]]) >= j) {
      current_aux <- subject_aux[[as.integer(subject)]][[j]]
    }
    current_val <- if (.collapsed_aux_entry_matches_phi(current_aux, phi[j, , drop = TRUE])) {
      .collapsed_aux_entry_value_at_lambda(current_aux, lambda)
    } else {
      NULL
    }
    if (!is.null(current_val) && is.finite(current_val$log_evidence)) {
      out[[j]] <- current_aux
      log_mhat[j] <- as.numeric(current_val$log_evidence)
      next
    }
    aux_prop <- list(
      aux_latent = aux_latents[j, ],
      seed_plan = .collapsed_seed_plan_from_latent(aux_latents[j, ], local_smc_control = local_smc_control),
      signature = .collapsed_seed_plan_signature(
        .collapsed_seed_plan_from_latent(aux_latents[j, ], local_smc_control = local_smc_control)
      )
    )
    cached <- .collapsed_cache_get(
      cache_store = cache_store,
      subject = subject,
      lambda = if (use_checkpoint_path) max(full_bridge_schedule) else lambda,
      phi = phi[j, , drop = TRUE],
      aux_signature = aux_prop$signature
    )
    cached_val <- .collapsed_aux_entry_value_at_lambda(cached, lambda)
    if (!is.null(cached) && !is.null(cached_val) && is.finite(cached_val$log_evidence)) {
      out[[j]] <- cached
      log_mhat[j] <- as.numeric(cached_val$log_evidence)
    } else {
      miss_idx <- c(miss_idx, j)
    }
  }
  if (length(miss_idx)) {
    miss_out <- parallel::mclapply(
      miss_idx,
      function(j) {
        current_aux <- NULL
        if (!is.null(subject_aux) &&
            length(subject_aux) >= as.integer(subject) &&
            length(subject_aux[[as.integer(subject)]]) >= j) {
          current_aux <- subject_aux[[as.integer(subject)]][[j]]
        }
        aux_latent <- aux_latents[j, ]
        seed_plan <- .collapsed_seed_plan_from_latent(aux_latent, local_smc_control = local_smc_control)
        warm <- .collapsed_pick_warm_start(
          subject_aux_entry = current_aux,
          cache_store = cache_store,
          subject = subject,
          lambda = lambda,
          phi = phi[j, , drop = TRUE],
          proposal_transform = proposal_transform,
          warm_start_mode = warm_start_mode
        )
        .collapsed_estimate_log_evidence(
          phi = phi[j, , drop = TRUE],
          subject = subject,
          lambda = if (use_checkpoint_path) max(full_bridge_schedule) else lambda,
          data_i = data_i,
          loglik_fn = loglik_fn,
          gaussian_map_fn = gaussian_map_fn,
          d_theta = d_theta,
          local_particles = local_particles,
          local_log_evidence_fn = local_log_evidence_fn,
          theta_names = theta_names,
          local_n_cores = local_n_cores,
          seed_plan = seed_plan,
          aux_latent = aux_latent,
          checkpoint_lambdas = if (use_checkpoint_path) full_bridge_schedule else NULL,
          auxiliary_state = list(seed_plan = seed_plan),
          verbose = verbose,
          warm_start_fit = warm$local_fit %||% NULL,
          warm_start_key = warm$cache_key %||% NA_character_,
          local_smc_control = local_smc_control,
          cache_store = NULL
        )
      },
      mc.cores = as.integer(max(1L, min(outer_n_cores, length(miss_idx))))
    )
    for (k in seq_along(miss_idx)) {
      j <- miss_idx[k]
      out[[j]] <- miss_out[[k]]
      .collapsed_cache_put(
        cache_store = cache_store,
        subject = subject,
        lambda = if (use_checkpoint_path) max(full_bridge_schedule) else lambda,
        phi = phi[j, , drop = TRUE],
        aux_signature = out[[j]]$aux_signature %||% .collapsed_auxiliary_signature(
          seed_plan = out[[j]]$seed_plan,
          auxiliary_state = out[[j]]$auxiliary_state
        ),
        value = .collapsed_trim_aux_entry(
          aux_entry = out[[j]],
          keep_local_fit = !is.null(cache_store) && isTRUE(cache_store$stats$keep_local_fit)
        )
      )
      out_val <- .collapsed_aux_entry_value_at_lambda(out[[j]], lambda)
      log_mhat[j] <- as.numeric(out_val$log_evidence)
    }
  }
  list(
    log_mhat = log_mhat,
    aux_entries = lapply(out, function(x) .collapsed_trim_aux_entry(x, keep_local_fit = store_local_fits)),
    n_local_evals = as.integer(length(miss_idx))
  )
}

.collapsed_log_abs_det_jacobian <- function(phi, proposal_transform, varphi = NULL) {
  as.numeric(proposal_transform$log_abs_det_jacobian(
    phi = as.numeric(phi),
    varphi = as.numeric(varphi %||% proposal_transform$forward(phi))
  ))
}

.collapsed_build_surrogate_bundle <- function(cache_store,
                                              active_subjects,
                                              active_lambdas,
                                              proposal_transform,
                                              surrogate_min_points = NULL,
                                              surrogate_ridge = 1e-4,
                                              surrogate_degree = 2L,
                                              surrogate_neighbors = NULL) {
  bundle <- vector("list", length(active_subjects))
  names(bundle) <- paste(active_subjects, format(active_lambdas, digits = 8), sep = "@")
  for (ii in seq_along(active_subjects)) {
    bundle[[ii]] <- .collapsed_build_surrogate(
      cache_store = cache_store,
      subject = active_subjects[ii],
      lambda = active_lambdas[ii],
      proposal_transform = proposal_transform,
      min_points = surrogate_min_points,
      ridge = surrogate_ridge,
      degree = surrogate_degree,
      neighbors = surrogate_neighbors
    )
  }
  bundle
}

.collapsed_rejuvenate_particles <- function(phi,
                                            w,
                                            logprior,
                                            subject_loglik,
                                            subject_aux,
                                            completed_subjects,
                                            data_list,
                                            loglik_fn,
                                            prior,
                                            gaussian_map_fn,
                                            d_theta,
                                            local_particles,
                                            local_log_evidence_fn = NULL,
                                            current_subject = NULL,
                                            current_lambda = 0,
                                            rejuvenation = c("recompute", "surrogate"),
                                            proposal = c("mix", "rw", "independence"),
                                            proposal_indep_prob = 0.7,
                                            rw_scale = 0.8,
                                            indep_scale = 1.0,
                                            proposal_transform = NULL,
                                            auxiliary_proposal = c("independent", "sticky", "frozen", "correlated"),
                                            auxiliary_stickiness = 0.9,
                                            auxiliary_correlation = 0.99,
                                            warm_start_mode = c("none", "nearest", "current_then_nearest"),
                                            bridge_schedule_lookup = NULL,
                                            surrogate_min_points = NULL,
                                            surrogate_ridge = 1e-4,
                                            surrogate_degree = 2L,
                                            surrogate_neighbors = NULL,
                                            theta_names = NULL,
                                            n_moves = 1L,
                                            local_n_cores = 1L,
                                            outer_n_cores = 1L,
                                            seed_base = 1L,
                                            verbose = FALSE,
                                            local_smc_control = list(),
                                            cache_store = NULL,
                                            store_local_fits = FALSE) {
  N <- nrow(phi)
  d_phi <- ncol(phi)
  rejuvenation <- match.arg(rejuvenation)
  proposal <- match.arg(proposal)
  auxiliary_proposal <- match.arg(auxiliary_proposal)
  warm_start_mode <- match.arg(warm_start_mode)
  proposal_transform <- .collapsed_normalize_proposal_transform(proposal_transform, d_phi = d_phi)
  current_lambda <- as.numeric(current_lambda %||% 0)
  has_current <- !is.null(current_subject) && is.finite(current_lambda) && current_lambda > 0
  active_subjects <- completed_subjects
  active_lambdas <- rep(1.0, length(completed_subjects))
  if (has_current) {
    active_subjects <- c(active_subjects, as.integer(current_subject))
    active_lambdas <- c(active_lambdas, as.numeric(current_lambda))
  }
  if (N <= 0L || (length(completed_subjects) <= 0L && !has_current)) {
    return(list(
      phi = phi,
      logprior = logprior,
      subject_loglik = subject_loglik,
      subject_aux = subject_aux,
      accept_rate = 0
    ))
  }

  surrogate_bundle <- if (identical(rejuvenation, "surrogate")) {
    .collapsed_build_surrogate_bundle(
      cache_store = cache_store,
      active_subjects = active_subjects,
      active_lambdas = active_lambdas,
      proposal_transform = proposal_transform,
      surrogate_min_points = surrogate_min_points,
      surrogate_ridge = surrogate_ridge,
      surrogate_degree = surrogate_degree,
      surrogate_neighbors = surrogate_neighbors
    )
  } else {
    NULL
  }

  accepted_total <- 0L
  proposed_total <- 0L
  proposed_rw <- 0L
  proposed_indep <- 0L
  accepted_rw <- 0L
  accepted_indep <- 0L
  local_eval_total <- 0L

  for (mv in seq_len(as.integer(n_moves))) {
    varphi <- .collapsed_transform_matrix(phi, proposal_transform)
    log_jac <- vapply(
      seq_len(N),
      function(i) .collapsed_log_abs_det_jacobian(phi[i, ], proposal_transform, varphi = varphi[i, ]),
      numeric(1)
    )
    S <- .collapsed_rejuvenation_cov(varphi, w)
    center <- .weighted_mean_vec(varphi, w)
    rw_cov <- (as.numeric(rw_scale)^2 / max(1L, d_phi)) * S + diag(1e-10, d_phi)
    indep_cov <- (as.numeric(indep_scale)^2) * S + diag(1e-10, d_phi)

    prop_list <- parallel::mclapply(
      seq_len(N),
      function(j) {
        set.seed(as.integer(seed_base + 100000L * mv + j))
        kernel <- proposal
        if (identical(proposal, "mix")) {
          kernel <- if (stats::runif(1L) < as.numeric(proposal_indep_prob)) "independence" else "rw"
        }
        if (identical(kernel, "independence")) {
          varphi_prop <- as.numeric(mvtnorm::rmvnorm(1L, mean = center, sigma = indep_cov))
          log_q_ratio <- mvtnorm::dmvnorm(varphi[j, ], mean = center, sigma = indep_cov, log = TRUE) -
            mvtnorm::dmvnorm(varphi_prop, mean = center, sigma = indep_cov, log = TRUE)
        } else {
          varphi_prop <- as.numeric(varphi[j, ] + mvtnorm::rmvnorm(1L, sigma = rw_cov))
          log_q_ratio <- 0
        }
        phi_prop <- tryCatch(as.numeric(proposal_transform$inverse(varphi_prop)), error = function(e) rep(NA_real_, d_phi))
        if (length(phi_prop) != d_phi || any(!is.finite(phi_prop))) {
          return(list(
            accepted = FALSE,
            phi = phi[j, ],
            logprior = logprior[j],
            subject_loglik_row = subject_loglik[j, ],
            subject_aux_row = .collapsed_subject_aux_row(subject_aux, j),
            kernel = kernel,
            n_local = 0L
          ))
        }
        lp_prop <- prior$lprior(phi_prop)
        if (!is.finite(lp_prop)) {
          return(list(
            accepted = FALSE,
            phi = phi[j, ],
            logprior = logprior[j],
            subject_loglik_row = subject_loglik[j, ],
            subject_aux_row = .collapsed_subject_aux_row(subject_aux, j),
            kernel = kernel,
            n_local = 0L
          ))
        }
        log_jac_prop <- .collapsed_log_abs_det_jacobian(phi_prop, proposal_transform, varphi = varphi_prop)

        prop_row <- subject_loglik[j, ]
        prop_aux_row <- .collapsed_subject_aux_row(subject_aux, j)
        n_local <- 0L
        curr_total <- 0
        for (k in seq_along(active_subjects)) {
          subj <- active_subjects[k]
          lambda_subj <- as.numeric(active_lambdas[k])
          current_aux <- if (!is.null(prop_aux_row) && length(prop_aux_row) >= subj) prop_aux_row[[subj]] else NULL
          surrogate_model <- if (!is.null(surrogate_bundle)) surrogate_bundle[[k]] else NULL
          curr_component <- if (identical(rejuvenation, "surrogate") && !is.null(surrogate_model)) {
            .collapsed_predict_surrogate(surrogate_model, phi[j, ])
          } else {
            subject_loglik[j, subj]
          }
          if (!is.finite(curr_component)) {
            curr_component <- subject_loglik[j, subj]
          }
          curr_total <- curr_total + curr_component
          if (identical(rejuvenation, "surrogate") && !is.null(surrogate_model)) {
            prop_component <- .collapsed_predict_surrogate(surrogate_model, phi_prop)
            if (!is.finite(prop_component)) {
              prop_component <- subject_loglik[j, subj]
            }
            prop_row[subj] <- prop_component
            if (!is.null(prop_aux_row)) {
              prop_aux_row[[subj]] <- .collapsed_trim_aux_entry(
                .collapsed_make_aux_entry(
                  phi = phi_prop,
                  subject = subj,
                  lambda = lambda_subj,
                  log_evidence = prop_component,
                  seed = NA_integer_,
                  local_fit = NULL,
                  source = "surrogate_prediction"
                ),
                keep_local_fit = store_local_fits
              )
            }
          } else {
            aux_prop <- .collapsed_propose_auxiliary_state(
              current_aux_entry = current_aux,
              mode = auxiliary_proposal,
              stickiness = auxiliary_stickiness,
              correlation = auxiliary_correlation,
              local_smc_control = local_smc_control
            )
            checkpoint_schedule <- NULL
            eval_lambda <- lambda_subj
            if (!is.null(bridge_schedule_lookup)) {
              sched <- .collapsed_lookup_bridge_schedule(bridge_schedule_lookup, subj, default = NULL)
              sched <- sort(unique(as.numeric(sched %||% numeric(0))))
              sched <- sched[is.finite(sched) & sched > 0 & sched <= 1]
              if (length(sched) > 1L && lambda_subj < 1 - 1e-12) {
                checkpoint_schedule <- sched
                eval_lambda <- max(sched)
              }
            }
            warm <- .collapsed_pick_warm_start(
              subject_aux_entry = current_aux,
              cache_store = cache_store,
              subject = subj,
              lambda = lambda_subj,
              phi = phi_prop,
              proposal_transform = proposal_transform,
              warm_start_mode = warm_start_mode
            )
            ev <- .collapsed_estimate_log_evidence(
              phi = phi_prop,
              subject = subj,
              lambda = eval_lambda,
              data_i = data_list[[subj]],
              loglik_fn = loglik_fn,
              gaussian_map_fn = gaussian_map_fn,
              d_theta = d_theta,
              local_particles = local_particles,
              local_log_evidence_fn = local_log_evidence_fn,
              theta_names = theta_names,
              local_n_cores = local_n_cores,
              seed = aux_prop$seed_plan$init %||% NA_integer_,
              seed_plan = aux_prop$seed_plan,
              aux_latent = aux_prop$aux_latent,
              checkpoint_lambdas = checkpoint_schedule,
              auxiliary_state = list(seed_plan = aux_prop$seed_plan),
              verbose = FALSE,
              warm_start_fit = warm$local_fit %||% NULL,
              warm_start_key = warm$cache_key %||% NA_character_,
              local_smc_control = local_smc_control,
              cache_store = NULL
            )
            n_local <- n_local + 1L
            prop_val <- .collapsed_aux_entry_value_at_lambda(ev, lambda_subj) %||% list(log_evidence = ev$log_evidence)
            prop_row[subj] <- prop_val$log_evidence
            if (!is.null(prop_aux_row)) {
              prop_aux_row[[subj]] <- .collapsed_trim_aux_entry(ev, keep_local_fit = store_local_fits)
            }
          }
        }

        prop_total <- sum(prop_row[active_subjects], na.rm = TRUE)
        log_ratio <- (lp_prop + prop_total + log_jac_prop) - (logprior[j] + curr_total + log_jac[j]) + log_q_ratio
        accept <- is.finite(log_ratio) && (log(runif(1L)) < min(0, log_ratio))
        if (!accept) {
          return(list(
            accepted = FALSE,
            phi = phi[j, ],
            logprior = logprior[j],
            subject_loglik_row = subject_loglik[j, ],
            subject_aux_row = .collapsed_subject_aux_row(subject_aux, j),
            kernel = kernel,
            n_local = n_local
          ))
        }
        list(
          accepted = TRUE,
          phi = phi_prop,
          logprior = lp_prop,
          subject_loglik_row = prop_row,
          subject_aux_row = prop_aux_row,
          kernel = kernel,
          n_local = n_local
        )
      },
      mc.cores = as.integer(max(1L, min(outer_n_cores, N)))
    )

    if (!is.null(cache_store) && isTRUE(cache_store$stats$enabled)) {
      for (jj in seq_along(prop_list)) {
        row_aux <- prop_list[[jj]]$subject_aux_row %||% NULL
        if (is.null(row_aux)) next
        for (subj in active_subjects) {
          aux_entry <- row_aux[[subj]] %||% NULL
          if (is.null(aux_entry)) next
          cache_lambda <- if (length(aux_entry$checkpoint_lambdas %||% numeric(0))) {
            max(aux_entry$checkpoint_lambdas)
          } else {
            aux_entry$lambda
          }
          .collapsed_cache_put(
            cache_store = cache_store,
            subject = subj,
            lambda = cache_lambda,
            phi = aux_entry$phi,
            aux_signature = aux_entry$aux_signature %||% .collapsed_auxiliary_signature(
              seed_plan = aux_entry$seed_plan,
              auxiliary_state = aux_entry$auxiliary_state
            ),
            value = .collapsed_trim_aux_entry(
              aux_entry = aux_entry,
              keep_local_fit = !is.null(cache_store) && isTRUE(cache_store$stats$keep_local_fit)
            )
          )
        }
      }
    }

    proposed_total <- proposed_total + N
    acc <- vapply(prop_list, function(x) isTRUE(x$accepted), logical(1))
    kernels <- vapply(prop_list, `[[`, character(1), "kernel")
    local_eval_total <- local_eval_total + sum(vapply(prop_list, `[[`, integer(1), "n_local"))
    proposed_rw <- proposed_rw + sum(kernels == "rw")
    proposed_indep <- proposed_indep + sum(kernels == "independence")
    accepted_rw <- accepted_rw + sum(acc & kernels == "rw")
    accepted_indep <- accepted_indep + sum(acc & kernels == "independence")
    accepted_total <- accepted_total + sum(acc)
    phi <- do.call(rbind, lapply(prop_list, `[[`, "phi"))
    logprior <- vapply(prop_list, `[[`, numeric(1), "logprior")
    subject_loglik <- do.call(rbind, lapply(prop_list, `[[`, "subject_loglik_row"))
    subject_aux <- .collapsed_rebuild_subject_aux(prop_list, ncol(subject_loglik))
  }

  list(
    phi = phi,
    logprior = logprior,
    subject_loglik = subject_loglik,
    subject_aux = subject_aux,
    accept_rate = accepted_total / max(1L, proposed_total),
    accept_rate_rw = accepted_rw / max(1L, proposed_rw),
    accept_rate_indep = accepted_indep / max(1L, proposed_indep),
    n_local_evals = as.integer(local_eval_total)
  )
}

.run_collapsed_schedule <- function(stage_plan,
                                    data_list,
                                    loglik_fn,
                                    prior,
                                    gaussian_map_fn,
                                    d_theta,
                                    N,
                                    local_particles,
                                    resample_threshold,
                                    local_log_evidence_fn,
                                    rejuvenation,
                                    rejuvenate_every,
                                    rejuvenate_after_resample,
                                    n_rejuvenation_moves,
                                    rejuvenation_proposal,
                                    rejuvenation_indep_prob,
                                    rejuvenation_rw_scale,
                                    rejuvenation_indep_scale,
                                    proposal_transform,
                                    auxiliary_proposal,
                                    auxiliary_stickiness,
                                    auxiliary_correlation,
                                    warm_start_mode,
                                    bridge_schedule_lookup,
                                    surrogate_min_points,
                                    surrogate_ridge,
                                    surrogate_degree,
                                    surrogate_neighbors,
                                    theta_names,
                                    local_n_cores,
                                    outer_n_cores,
                                    base_seed,
                                    cache_store,
                                    store_local_fits,
                                    verbose,
                                    local_smc_control) {
  state <- .collapsed_initialize_state(
    prior = prior,
    N = N,
    base_seed = base_seed,
    n_subjects = length(data_list),
    store_local_fits = store_local_fits
  )
  ess_hist <- c(.ess_from_weights(state$w))
  stage_logZ <- numeric(nrow(stage_plan))
  resampled_hist <- logical(nrow(stage_plan))
  rejuvenation_accept_hist <- rep(NA_real_, nrow(stage_plan))
  n_local_runs <- 0L
  n_subjects <- length(unique(stage_plan$subject_position))

  for (ss in seq_len(nrow(stage_plan))) {
    stage_row <- stage_plan[ss, , drop = FALSE]
    state <- .collapsed_enter_subject(state, stage_row$subject)
    .collapsed_verbose_stage_header(stage_row = stage_row, n_subjects = n_subjects, verbose = verbose)

    inc <- .collapsed_stage_increment(
      phi = state$phi,
      subject = stage_row$subject,
      data_i = data_list[[stage_row$subject]],
      loglik_fn = loglik_fn,
      gaussian_map_fn = gaussian_map_fn,
      d_theta = d_theta,
      local_particles = local_particles,
      local_log_evidence_fn = local_log_evidence_fn,
      lambda = stage_row$lambda,
      theta_names = theta_names,
      local_n_cores = local_n_cores,
      outer_n_cores = outer_n_cores,
      seed_base = as.integer(base_seed + 1000000L * ss),
      verbose = FALSE,
      subject_aux = state$subject_aux,
      proposal_transform = proposal_transform,
      warm_start_mode = warm_start_mode,
      full_bridge_schedule = .collapsed_lookup_bridge_schedule(bridge_schedule_lookup, stage_row$subject, default = stage_row$lambda),
      local_smc_control = local_smc_control,
      cache_store = cache_store,
      store_local_fits = store_local_fits
    )
    n_local_runs <- n_local_runs + (inc$n_local_evals %||% nrow(state$phi))

    old_loglik <- state$subject_loglik[, stage_row$subject]
    old_loglik[!is.finite(old_loglik)] <- 0
    log_inc <- inc$log_mhat - old_loglik
    norm <- .normalize_logweights(log(pmax(state$w, .Machine$double.eps)) + log_inc)
    state$w <- norm$w
    state$log_evidence <- state$log_evidence + norm$log_norm
    stage_logZ[ss] <- norm$log_norm
    state$subject_loglik[, stage_row$subject] <- inc$log_mhat
    state$subject_aux <- .collapsed_set_subject_aux(state$subject_aux, stage_row$subject, inc$aux_entries)
    ess_now <- .ess_from_weights(state$w)
    ess_hist <- c(ess_hist, ess_now)

    resampled <- FALSE
    if (ess_now < resample_threshold * nrow(state$phi)) {
      state <- .collapsed_resample_state(state)
      resampled <- TRUE
      ess_now <- .ess_from_weights(state$w)
      ess_hist[length(ess_hist)] <- ess_now
      .collapsed_verbose_resample(stage_row = stage_row, ess_now = ess_now, verbose = verbose)
    }
    resampled_hist[ss] <- resampled

    if (.collapsed_should_rejuvenate(
      rejuvenation = rejuvenation,
      rejuvenate_every = rejuvenate_every,
      rejuvenate_after_resample = rejuvenate_after_resample,
      resampled = resampled,
      stage_index = ss
    )) {
      rej <- .collapsed_rejuvenate_particles(
        phi = state$phi,
        w = state$w,
        logprior = state$logprior,
        subject_loglik = state$subject_loglik,
        subject_aux = state$subject_aux,
        completed_subjects = state$completed_subjects,
        data_list = data_list,
        loglik_fn = loglik_fn,
        prior = prior,
        gaussian_map_fn = gaussian_map_fn,
        d_theta = d_theta,
        local_particles = local_particles,
        local_log_evidence_fn = local_log_evidence_fn,
        current_subject = state$current_subject,
        current_lambda = stage_row$lambda,
        rejuvenation = rejuvenation,
        proposal = rejuvenation_proposal,
        proposal_indep_prob = rejuvenation_indep_prob,
        rw_scale = rejuvenation_rw_scale,
        indep_scale = rejuvenation_indep_scale,
        proposal_transform = proposal_transform,
        auxiliary_proposal = auxiliary_proposal,
        auxiliary_stickiness = auxiliary_stickiness,
        auxiliary_correlation = auxiliary_correlation,
        warm_start_mode = warm_start_mode,
        bridge_schedule_lookup = bridge_schedule_lookup,
        surrogate_min_points = surrogate_min_points,
        surrogate_ridge = surrogate_ridge,
        surrogate_degree = surrogate_degree,
        surrogate_neighbors = surrogate_neighbors,
        theta_names = theta_names,
        n_moves = n_rejuvenation_moves,
        local_n_cores = local_n_cores,
        outer_n_cores = outer_n_cores,
        seed_base = as.integer(base_seed + 2000000L * ss),
        verbose = FALSE,
        local_smc_control = local_smc_control,
        cache_store = cache_store,
        store_local_fits = store_local_fits
      )
      state$phi <- rej$phi
      state$logprior <- rej$logprior
      state$subject_loglik <- rej$subject_loglik
      state$subject_aux <- rej$subject_aux
      rejuvenation_accept_hist[ss] <- rej$accept_rate
      n_local_runs <- n_local_runs + as.integer(rej$n_local_evals %||% 0L)
      .collapsed_verbose_rejuvenation(stage_row = stage_row, accept_rate = rej$accept_rate, verbose = verbose)
    }

    .collapsed_verbose_stage_footer(
      stage_row = stage_row,
      ess_now = ess_now,
      logz_inc = stage_logZ[ss],
      logz_total = state$log_evidence,
      verbose = verbose
    )

    if (isTRUE(stage_row$subject_complete)) {
      state <- .collapsed_finalize_current_subject(state)
    }
  }

  list(
    phi = state$phi,
    w = state$w,
    log_evidence = state$log_evidence,
    subject_loglik = state$subject_loglik,
    subject_aux = state$subject_aux,
    ess_hist = ess_hist,
    stage_logZ = stage_logZ,
    resampled = resampled_hist,
    rejuvenation_accept_hist = rejuvenation_accept_hist,
    n_local_runs = as.integer(n_local_runs),
    stage_subject = stage_plan$subject,
    stage_lambda = stage_plan$lambda
  )
}

.collapsed_eval_lambda_candidate <- function(state,
                                             subject,
                                             lambda,
                                             data_i,
                                             loglik_fn,
                                             gaussian_map_fn,
                                             d_theta,
                                             local_particles,
                                             local_log_evidence_fn,
                                             theta_names,
                                             local_n_cores,
                                             outer_n_cores,
                                             seed_base,
                                             proposal_transform,
                                             warm_start_mode,
                                             local_smc_control,
                                             cache_store,
                                             store_local_fits = FALSE) {
  inc <- .collapsed_stage_increment(
    phi = state$phi,
    subject = subject,
    data_i = data_i,
    loglik_fn = loglik_fn,
    gaussian_map_fn = gaussian_map_fn,
    d_theta = d_theta,
    local_particles = local_particles,
    local_log_evidence_fn = local_log_evidence_fn,
    lambda = lambda,
    theta_names = theta_names,
    local_n_cores = local_n_cores,
    outer_n_cores = outer_n_cores,
    seed_base = seed_base,
    verbose = FALSE,
    subject_aux = state$subject_aux,
    proposal_transform = proposal_transform,
    warm_start_mode = warm_start_mode,
    local_smc_control = local_smc_control,
    cache_store = cache_store,
    store_local_fits = store_local_fits
  )
  old_loglik <- state$subject_loglik[, subject]
  old_loglik[!is.finite(old_loglik)] <- 0
  logw <- log(pmax(state$w, .Machine$double.eps)) + (inc$log_mhat - old_loglik)
  norm <- .normalize_logweights(logw)
  list(
    inc = inc,
    ess = .ess_from_weights(norm$w),
    log_norm = norm$log_norm,
    n_local_evals = inc$n_local_evals %||% nrow(state$phi)
  )
}

.collapsed_select_next_lambda <- function(state,
                                          subject,
                                          lambda_prev,
                                          data_i,
                                          loglik_fn,
                                          gaussian_map_fn,
                                          d_theta,
                                          local_particles,
                                          local_log_evidence_fn,
                                          theta_names,
                                          local_n_cores,
                                          outer_n_cores,
                                          seed_base,
                                          proposal_transform,
                                          warm_start_mode,
                                          local_smc_control,
                                          cache_store,
                                          store_local_fits,
                                          ess_target,
                                          lambda_tol,
                                          max_bisect) {
  eval_full <- .collapsed_eval_lambda_candidate(
    state = state,
    subject = subject,
    lambda = 1.0,
    data_i = data_i,
    loglik_fn = loglik_fn,
    gaussian_map_fn = gaussian_map_fn,
    d_theta = d_theta,
    local_particles = local_particles,
    local_log_evidence_fn = local_log_evidence_fn,
    theta_names = theta_names,
    local_n_cores = local_n_cores,
    outer_n_cores = outer_n_cores,
    seed_base = seed_base,
    proposal_transform = proposal_transform,
    warm_start_mode = warm_start_mode,
    local_smc_control = local_smc_control,
    cache_store = cache_store,
    store_local_fits = store_local_fits
  )
  if (eval_full$ess >= ess_target || (1 - lambda_prev) <= lambda_tol) {
    return(list(lambda = 1.0, eval = eval_full))
  }

  left <- lambda_prev
  right <- 1.0
  best_lambda <- lambda_prev
  best_eval <- NULL

  for (iter in seq_len(as.integer(max_bisect))) {
    mid <- (left + right) / 2
    if ((right - left) <= lambda_tol) {
      break
    }
    eval_mid <- .collapsed_eval_lambda_candidate(
      state = state,
      subject = subject,
      lambda = mid,
      data_i = data_i,
      loglik_fn = loglik_fn,
      gaussian_map_fn = gaussian_map_fn,
      d_theta = d_theta,
      local_particles = local_particles,
      local_log_evidence_fn = local_log_evidence_fn,
      theta_names = theta_names,
      local_n_cores = local_n_cores,
      outer_n_cores = outer_n_cores,
      seed_base = seed_base + iter,
      proposal_transform = proposal_transform,
      warm_start_mode = warm_start_mode,
      local_smc_control = local_smc_control,
      cache_store = cache_store,
      store_local_fits = store_local_fits
    )
    if (eval_mid$ess >= ess_target) {
      best_lambda <- mid
      best_eval <- eval_mid
      left <- mid
    } else {
      right <- mid
    }
  }

  if (is.null(best_eval)) {
    best_lambda <- max(lambda_prev + lambda_tol, (lambda_prev + right) / 2)
    best_eval <- .collapsed_eval_lambda_candidate(
      state = state,
      subject = subject,
      lambda = best_lambda,
      data_i = data_i,
      loglik_fn = loglik_fn,
      gaussian_map_fn = gaussian_map_fn,
      d_theta = d_theta,
      local_particles = local_particles,
      local_log_evidence_fn = local_log_evidence_fn,
      theta_names = theta_names,
      local_n_cores = local_n_cores,
      outer_n_cores = outer_n_cores,
      seed_base = seed_base + 500L,
      proposal_transform = proposal_transform,
      warm_start_mode = warm_start_mode,
      local_smc_control = local_smc_control,
      cache_store = cache_store,
      store_local_fits = store_local_fits
    )
  }

  list(lambda = best_lambda, eval = best_eval)
}

.run_collapsed_bridge_adaptive <- function(data_list,
                                           subject_order,
                                           loglik_fn,
                                           prior,
                                           gaussian_map_fn,
                                           d_theta,
                                           N,
                                           local_particles,
                                           resample_threshold,
                                           local_log_evidence_fn,
                                           rejuvenation,
                                           rejuvenate_every,
                                           rejuvenate_after_resample,
                                           n_rejuvenation_moves,
                                           rejuvenation_proposal,
                                           rejuvenation_indep_prob,
                                           rejuvenation_rw_scale,
                                           rejuvenation_indep_scale,
                                           proposal_transform,
                                           auxiliary_proposal,
                                           auxiliary_stickiness,
                                           auxiliary_correlation,
                                           warm_start_mode,
                                           surrogate_min_points,
                                           surrogate_ridge,
                                           surrogate_degree,
                                           surrogate_neighbors,
                                           theta_names,
                                           local_n_cores,
                                           outer_n_cores,
                                           base_seed,
                                           cache_store,
                                           store_local_fits,
                                           adaptive_ess_target,
                                           adaptive_lambda_tol,
                                           adaptive_max_steps_per_subject,
                                           adaptive_max_bisect,
                                           verbose,
                                           local_smc_control) {
  state <- .collapsed_initialize_state(
    prior = prior,
    N = N,
    base_seed = base_seed,
    n_subjects = length(data_list),
    store_local_fits = store_local_fits
  )
  ess_hist <- c(.ess_from_weights(state$w))
  stage_logZ <- numeric(0)
  stage_subject <- integer(0)
  stage_lambda <- numeric(0)
  resampled_hist <- logical(0)
  rejuvenation_accept_hist <- numeric(0)
  n_local_runs <- 0L
  bridge_schedule_realized <- list()
  stage_counter <- 0L

  for (tt in seq_along(subject_order)) {
    subj <- subject_order[tt]
    state <- .collapsed_enter_subject(state, subj)
    lambda_prev <- 0
    realized_subject_schedule <- numeric(0)
    if (isTRUE(verbose)) {
      cat(sprintf("Subject %d/%d: subject %d\n", tt, length(subject_order), subj))
    }
    for (ll in seq_len(as.integer(adaptive_max_steps_per_subject))) {
      selected <- .collapsed_select_next_lambda(
        state = state,
        subject = subj,
        lambda_prev = lambda_prev,
        data_i = data_list[[subj]],
        loglik_fn = loglik_fn,
        gaussian_map_fn = gaussian_map_fn,
        d_theta = d_theta,
        local_particles = local_particles,
        local_log_evidence_fn = local_log_evidence_fn,
        theta_names = theta_names,
        local_n_cores = local_n_cores,
        outer_n_cores = outer_n_cores,
        seed_base = as.integer(base_seed + 1000000L * (stage_counter + 1L)),
        proposal_transform = proposal_transform,
        warm_start_mode = warm_start_mode,
        local_smc_control = local_smc_control,
        cache_store = cache_store,
        store_local_fits = store_local_fits,
        ess_target = adaptive_ess_target * N,
        lambda_tol = adaptive_lambda_tol,
        max_bisect = adaptive_max_bisect
      )
      lambda_now <- selected$lambda
      inc <- selected$eval$inc
      stage_counter <- stage_counter + 1L
      realized_subject_schedule <- c(realized_subject_schedule, lambda_now)
      stage_subject <- c(stage_subject, subj)
      stage_lambda <- c(stage_lambda, lambda_now)
      rejuvenation_accept_hist <- c(rejuvenation_accept_hist, NA_real_)
      if (isTRUE(verbose)) {
        cat(sprintf("  Adaptive bridge stage %d: lambda %.4f\n", ll, lambda_now))
      }
      n_local_runs <- n_local_runs + (inc$n_local_evals %||% nrow(state$phi))

      old_loglik <- state$subject_loglik[, subj]
      old_loglik[!is.finite(old_loglik)] <- 0
      norm <- .normalize_logweights(log(pmax(state$w, .Machine$double.eps)) + (inc$log_mhat - old_loglik))
      state$w <- norm$w
      state$log_evidence <- state$log_evidence + norm$log_norm
      stage_logZ <- c(stage_logZ, norm$log_norm)
      state$subject_loglik[, subj] <- inc$log_mhat
      state$subject_aux <- .collapsed_set_subject_aux(state$subject_aux, subj, inc$aux_entries)
      ess_now <- .ess_from_weights(state$w)
      ess_hist <- c(ess_hist, ess_now)

      resampled <- FALSE
      if (ess_now < resample_threshold * nrow(state$phi)) {
        state <- .collapsed_resample_state(state)
        resampled <- TRUE
        ess_now <- .ess_from_weights(state$w)
        ess_hist[length(ess_hist)] <- ess_now
        if (isTRUE(verbose)) cat(sprintf("    resampled at ESS=%.1f\n", ess_now))
      }
      resampled_hist <- c(resampled_hist, resampled)

      if (.collapsed_should_rejuvenate(
        rejuvenation = rejuvenation,
        rejuvenate_every = rejuvenate_every,
        rejuvenate_after_resample = rejuvenate_after_resample,
        resampled = resampled,
        stage_index = stage_counter
      )) {
        rej <- .collapsed_rejuvenate_particles(
          phi = state$phi,
          w = state$w,
          logprior = state$logprior,
          subject_loglik = state$subject_loglik,
          subject_aux = state$subject_aux,
          completed_subjects = state$completed_subjects,
          data_list = data_list,
          loglik_fn = loglik_fn,
          prior = prior,
          gaussian_map_fn = gaussian_map_fn,
          d_theta = d_theta,
          local_particles = local_particles,
          local_log_evidence_fn = local_log_evidence_fn,
          current_subject = subj,
          current_lambda = lambda_now,
          rejuvenation = rejuvenation,
          proposal = rejuvenation_proposal,
          proposal_indep_prob = rejuvenation_indep_prob,
          rw_scale = rejuvenation_rw_scale,
          indep_scale = rejuvenation_indep_scale,
          proposal_transform = proposal_transform,
          auxiliary_proposal = auxiliary_proposal,
          auxiliary_stickiness = auxiliary_stickiness,
          auxiliary_correlation = auxiliary_correlation,
          warm_start_mode = warm_start_mode,
          bridge_schedule_lookup = setNames(list(realized_subject_schedule), as.character(subj)),
          surrogate_min_points = surrogate_min_points,
          surrogate_ridge = surrogate_ridge,
          surrogate_degree = surrogate_degree,
          surrogate_neighbors = surrogate_neighbors,
          theta_names = theta_names,
          n_moves = n_rejuvenation_moves,
          local_n_cores = local_n_cores,
          outer_n_cores = outer_n_cores,
          seed_base = as.integer(base_seed + 2000000L * stage_counter),
          verbose = FALSE,
          local_smc_control = local_smc_control,
          cache_store = cache_store,
          store_local_fits = store_local_fits
        )
        state$phi <- rej$phi
        state$logprior <- rej$logprior
        state$subject_loglik <- rej$subject_loglik
        state$subject_aux <- rej$subject_aux
        rejuvenation_accept_hist[stage_counter] <- rej$accept_rate
        n_local_runs <- n_local_runs + as.integer(rej$n_local_evals %||% 0L)
        if (isTRUE(verbose)) cat(sprintf("    rejuvenation accept=%.3f\n", rej$accept_rate))
      }

      if (isTRUE(verbose)) {
        cat(sprintf("    ESS=%.1f | logZ+=%.4f -> %.4f\n", ess_now, stage_logZ[stage_counter], state$log_evidence))
      }

      lambda_prev <- lambda_now
      if (lambda_now >= 1 - 1e-12) {
        break
      }
    }
    bridge_schedule_realized[[tt]] <- realized_subject_schedule
    state <- .collapsed_finalize_current_subject(state)
  }

  list(
    phi = state$phi,
    w = state$w,
    log_evidence = state$log_evidence,
    subject_loglik = state$subject_loglik,
    subject_aux = state$subject_aux,
    ess_hist = ess_hist,
    stage_logZ = stage_logZ,
    stage_subject = stage_subject,
    stage_lambda = stage_lambda,
    resampled = resampled_hist,
    rejuvenation_accept_hist = rejuvenation_accept_hist,
    n_local_runs = as.integer(n_local_runs),
    bridge_schedule = bridge_schedule_realized
  )
}

.collapsed_exactness_meta <- function(assume_unbiased_local_evidence,
                                      local_log_evidence_fn,
                                      store_local_fits = FALSE,
                                      rejuvenation = "recompute",
                                      bridge_mode = c("subject", "checkpoint_path", "adaptive_repeated_partial"),
                                      auxiliary_proposal = "independent",
                                      warm_start_mode = "none",
                                      proposal_transform = .collapsed_identity_transform()) {
  bridge_mode <- match.arg(bridge_mode)
  has_nested_aux_contract <- !is.function(local_log_evidence_fn)
  pm_viable <- isTRUE(assume_unbiased_local_evidence) &&
    !identical(rejuvenation, "surrogate") &&
    has_nested_aux_contract &&
    !identical(bridge_mode, "adaptive_repeated_partial")
  list(
    mode = if (identical(rejuvenation, "surrogate")) {
      "approximate_surrogate_rejuvenation"
    } else if (pm_viable) {
      "pseudo_marginal_under_local_evidence_assumption"
    } else {
      "approximate"
    },
    requires_positive_unbiased_local_evidence = isTRUE(assume_unbiased_local_evidence),
    local_backend = if (is.function(local_log_evidence_fn)) {
      "custom_local_log_evidence_fn"
    } else {
      "nested_local_smc"
    },
    carries_current_auxiliary_estimates = TRUE,
    auxiliary_state = if (has_nested_aux_contract) {
      "per_subject_local_seed_plan_and_optional_checkpoint_path"
    } else {
      "custom_backend_defined"
    },
    bridge_mode = bridge_mode,
    stores_full_local_fits = isTRUE(store_local_fits),
    rejuvenation = rejuvenation,
    auxiliary_proposal = auxiliary_proposal,
    warm_start_mode = warm_start_mode,
    transformed_proposal_space = proposal_transform$name %||% "identity"
  )
}

collapsed_subject_smc <- function(data_list,
                                  loglik_fn,
                                  prior,
                                  gaussian_map_fn,
                                  d_theta,
                                  N = 500L,
                                  local_particles = 1000L,
                                  resample_threshold = 0.5,
                                  local_log_evidence_fn = NULL,
                                  subject_order = NULL,
                                  rejuvenation = c("none", "recompute", "surrogate"),
                                  rejuvenate_every = Inf,
                                  rejuvenate_after_resample = FALSE,
                                  n_rejuvenation_moves = 1L,
                                  rejuvenation_proposal = c("mix", "rw", "independence"),
                                  rejuvenation_indep_prob = 0.7,
                                  rejuvenation_rw_scale = 0.8,
                                  rejuvenation_indep_scale = 1.0,
                                  proposal_transform = NULL,
                                  auxiliary_proposal = c("independent", "sticky", "frozen", "correlated"),
                                  auxiliary_stickiness = 0.9,
                                  auxiliary_correlation = 0.99,
                                  warm_start_mode = c("none", "nearest", "current_then_nearest"),
                                  surrogate_min_points = NULL,
                                  surrogate_ridge = 1e-4,
                                  surrogate_degree = 2L,
                                  surrogate_neighbors = NULL,
                                  theta_names = NULL,
                                  local_n_cores = 1L,
                                  outer_n_cores = 1L,
                                  base_seed = 123L,
                                  assume_unbiased_local_evidence = FALSE,
                                  evidence_cache = TRUE,
                                  evidence_cache_digits = 8L,
                                  evidence_cache_stochastic_by_seed = NULL,
                                  store_local_fits = FALSE,
                                  verbose = TRUE,
                                  local_smc_control = list()) {
  rejuvenation <- match.arg(rejuvenation)
  rejuvenation_proposal <- match.arg(rejuvenation_proposal)
  auxiliary_proposal <- match.arg(auxiliary_proposal)
  warm_start_mode <- match.arg(warm_start_mode)
  proposal_transform <- .collapsed_normalize_proposal_transform(proposal_transform, d_phi = .collapsed_prior_dim(prior))
  cache_keep_local_fit <- isTRUE(store_local_fits) || !identical(warm_start_mode, "none")
  evidence_cache_stochastic_by_seed <- isTRUE(evidence_cache_stochastic_by_seed %||% is.null(local_log_evidence_fn))
  subject_order <- .collapsed_validate_subject_order(data_list, subject_order)
  stage_plan <- .collapsed_stage_plan(subject_order = subject_order, bridge_schedule = 1.0)
  cache_store <- .collapsed_make_evidence_cache(
    enabled = evidence_cache,
    digits = evidence_cache_digits,
    keep_local_fit = cache_keep_local_fit,
    stochastic_by_seed = evidence_cache_stochastic_by_seed
  )

  fit <- .run_collapsed_schedule(
    stage_plan = stage_plan,
    data_list = data_list,
    loglik_fn = loglik_fn,
    prior = prior,
    gaussian_map_fn = gaussian_map_fn,
    d_theta = d_theta,
    N = N,
    local_particles = local_particles,
    resample_threshold = resample_threshold,
    local_log_evidence_fn = local_log_evidence_fn,
    rejuvenation = rejuvenation,
    rejuvenate_every = rejuvenate_every,
    rejuvenate_after_resample = rejuvenate_after_resample,
    n_rejuvenation_moves = n_rejuvenation_moves,
    rejuvenation_proposal = rejuvenation_proposal,
    rejuvenation_indep_prob = rejuvenation_indep_prob,
    rejuvenation_rw_scale = rejuvenation_rw_scale,
    rejuvenation_indep_scale = rejuvenation_indep_scale,
    proposal_transform = proposal_transform,
    auxiliary_proposal = auxiliary_proposal,
    auxiliary_stickiness = auxiliary_stickiness,
    auxiliary_correlation = auxiliary_correlation,
    warm_start_mode = warm_start_mode,
    bridge_schedule_lookup = list(),
    surrogate_min_points = surrogate_min_points,
    surrogate_ridge = surrogate_ridge,
    surrogate_degree = surrogate_degree,
    surrogate_neighbors = surrogate_neighbors,
    theta_names = theta_names,
    local_n_cores = local_n_cores,
    outer_n_cores = outer_n_cores,
    base_seed = base_seed,
    cache_store = cache_store,
    store_local_fits = store_local_fits,
    verbose = verbose,
    local_smc_control = local_smc_control
  )

  list(
    phi = fit$phi,
    w = fit$w,
    log_evidence = fit$log_evidence,
    subject_loglik = fit$subject_loglik,
    subject_aux = fit$subject_aux,
    meta = list(
      method = "collapsed_subject",
      subject_order = subject_order,
      n_subjects = length(subject_order),
      ess_hist = fit$ess_hist,
      stage_logZ = fit$stage_logZ,
      resampled = fit$resampled,
      rejuvenation = rejuvenation,
      rejuvenation_proposal = rejuvenation_proposal,
      proposal_transform = proposal_transform$name,
      auxiliary_proposal = auxiliary_proposal,
      auxiliary_stickiness = auxiliary_stickiness,
      auxiliary_correlation = auxiliary_correlation,
      warm_start_mode = warm_start_mode,
      surrogate = list(
        min_points = surrogate_min_points,
        ridge = surrogate_ridge,
        degree = surrogate_degree,
        neighbors = surrogate_neighbors
      ),
      rejuvenation_accept_hist = fit$rejuvenation_accept_hist,
      n_local_runs = fit$n_local_runs,
      evidence_cache = .collapsed_cache_stats(cache_store),
      extended_space = TRUE,
      exactness = .collapsed_exactness_meta(
        assume_unbiased_local_evidence = assume_unbiased_local_evidence,
        local_log_evidence_fn = local_log_evidence_fn,
        store_local_fits = cache_keep_local_fit,
        rejuvenation = rejuvenation,
        bridge_mode = "subject",
        auxiliary_proposal = auxiliary_proposal,
        warm_start_mode = warm_start_mode,
        proposal_transform = proposal_transform
      )
    )
  )
}

recover_collapsed_subject_posteriors <- function(fit,
                                                 data_list,
                                                 loglik_fn,
                                                 gaussian_map_fn,
                                                 d_theta,
                                                 subjects = seq_along(data_list),
                                                 theta_names = NULL,
                                                 local_particles = 1000L,
                                                 max_outer_particles = 64L,
                                                 local_n_cores = 1L,
                                                 outer_n_cores = 1L,
                                                 base_seed = 9001L,
                                                 local_smc_control = list()) {
  phi <- as.matrix(fit$phi)
  w <- pmax(as.numeric(fit$w), 0)
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) {
    w <- rep(1 / nrow(phi), nrow(phi))
  } else {
    w <- w / sw
  }
  keep <- seq_len(nrow(phi))
  if (length(keep) > as.integer(max_outer_particles)) {
    set.seed(as.integer(base_seed))
    keep <- sample.int(nrow(phi), size = as.integer(max_outer_particles), replace = FALSE, prob = w)
  }
  phi_keep <- phi[keep, , drop = FALSE]
  w_keep <- w[keep]
  w_keep <- w_keep / sum(w_keep)

  out <- vector("list", length(subjects))
  names(out) <- paste0("subject_", subjects)
  for (ii in seq_along(subjects)) {
    subj <- as.integer(subjects[ii])
    local_parts <- vector("list", nrow(phi_keep))
    for (jj in seq_len(nrow(phi_keep))) {
      est <- NULL
      if (!is.null(fit$subject_aux) &&
          length(fit$subject_aux) >= subj &&
          !is.null(fit$subject_aux[[subj]]) &&
          length(fit$subject_aux[[subj]]) >= keep[jj]) {
        est_aux <- fit$subject_aux[[subj]][[keep[jj]]]
        if (!is.null(est_aux) && !is.null(est_aux$local_fit)) {
          est <- list(local_fit = est_aux$local_fit)
        }
      }
      if (is.null(est)) {
        est <- .estimate_local_log_evidence_at_phi(
          phi = phi_keep[jj, ],
          data_i = data_list[[subj]],
          loglik_fn = loglik_fn,
          gaussian_map_fn = gaussian_map_fn,
          d_theta = d_theta,
          local_particles = local_particles,
          theta_names = theta_names,
          local_n_cores = local_n_cores,
          seed = as.integer(base_seed + 10000L * ii + jj),
          verbose = FALSE,
          local_smc_control = local_smc_control
        )
      }
      if (is.null(est$local_fit) || is.null(est$local_fit$Theta) || is.null(est$local_fit$w)) {
        stop("Posterior recovery requires stored local fits or nested local SMC output with Theta and w.")
      }
      Theta <- as.matrix(est$local_fit$Theta)
      w_local <- pmax(as.numeric(est$local_fit$w), 0)
      sw_local <- sum(w_local)
      if (!is.finite(sw_local) || sw_local <= 0) {
        w_local <- rep(1 / nrow(Theta), nrow(Theta))
      } else {
        w_local <- w_local / sw_local
      }
      local_parts[[jj]] <- list(alpha = Theta, w = w_keep[jj] * w_local)
    }
    alpha_all <- do.call(rbind, lapply(local_parts, `[[`, "alpha"))
    w_all <- unlist(lapply(local_parts, `[[`, "w"))
    w_all <- w_all / sum(w_all)
    out[[ii]] <- list(
      alpha = alpha_all,
      w = w_all,
      mean = colSums(alpha_all * w_all),
      subject = subj
    )
  }
  out
}

collapsed_bridge_smc <- function(data_list,
                                 loglik_fn,
                                 prior,
                                 gaussian_map_fn,
                                 d_theta,
                                 N = 500L,
                                 local_particles = 1000L,
                                 resample_threshold = 0.5,
                                 local_log_evidence_fn = NULL,
                                 subject_order = NULL,
                                 bridge_schedule = NULL,
                                 n_bridge_steps = 4L,
                                 rejuvenation = c("none", "recompute", "surrogate"),
                                 rejuvenate_every = Inf,
                                 rejuvenate_after_resample = FALSE,
                                 n_rejuvenation_moves = 1L,
                                 rejuvenation_proposal = c("mix", "rw", "independence"),
                                 rejuvenation_indep_prob = 0.7,
                                 rejuvenation_rw_scale = 0.8,
                                 rejuvenation_indep_scale = 1.0,
                                 proposal_transform = NULL,
                                 auxiliary_proposal = c("independent", "sticky", "frozen", "correlated"),
                                 auxiliary_stickiness = 0.9,
                                 auxiliary_correlation = 0.99,
                                 warm_start_mode = c("none", "nearest", "current_then_nearest"),
                                 surrogate_min_points = NULL,
                                 surrogate_ridge = 1e-4,
                                 surrogate_degree = 2L,
                                 surrogate_neighbors = NULL,
                                 theta_names = NULL,
                                 local_n_cores = 1L,
                                 outer_n_cores = 1L,
                                 base_seed = 123L,
                                 assume_unbiased_local_evidence = FALSE,
                                  evidence_cache = TRUE,
                                 evidence_cache_digits = 8L,
                                 evidence_cache_stochastic_by_seed = NULL,
                                 store_local_fits = FALSE,
                                 adaptive_bridge = FALSE,
                                 adaptive_ess_target = 0.8,
                                 adaptive_lambda_tol = 1e-3,
                                 adaptive_max_steps_per_subject = 32L,
                                 adaptive_max_bisect = 12L,
                                 verbose = TRUE,
                                 local_smc_control = list()) {
  rejuvenation <- match.arg(rejuvenation)
  rejuvenation_proposal <- match.arg(rejuvenation_proposal)
  auxiliary_proposal <- match.arg(auxiliary_proposal)
  warm_start_mode <- match.arg(warm_start_mode)
  proposal_transform <- .collapsed_normalize_proposal_transform(proposal_transform, d_phi = .collapsed_prior_dim(prior))
  cache_keep_local_fit <- isTRUE(store_local_fits) || !identical(warm_start_mode, "none")
  evidence_cache_stochastic_by_seed <- isTRUE(evidence_cache_stochastic_by_seed %||% is.null(local_log_evidence_fn))
  subject_order <- .collapsed_validate_subject_order(data_list, subject_order)
  cache_store <- .collapsed_make_evidence_cache(
    enabled = evidence_cache,
    digits = evidence_cache_digits,
    keep_local_fit = cache_keep_local_fit,
    stochastic_by_seed = evidence_cache_stochastic_by_seed
  )
  if (!isTRUE(adaptive_bridge)) {
    bridge_schedule <- .normalize_bridge_schedule(
      bridge_schedule = bridge_schedule,
      n_bridge_steps = n_bridge_steps
    )
    stage_plan <- .collapsed_stage_plan(
      subject_order = subject_order,
      bridge_schedule = bridge_schedule
    )
  } else {
    bridge_schedule <- NULL
    stage_plan <- NULL
  }

  fit <- if (!isTRUE(adaptive_bridge)) {
    .run_collapsed_schedule(
      stage_plan = stage_plan,
      data_list = data_list,
      loglik_fn = loglik_fn,
      prior = prior,
      gaussian_map_fn = gaussian_map_fn,
      d_theta = d_theta,
      N = N,
      local_particles = local_particles,
      resample_threshold = resample_threshold,
      local_log_evidence_fn = local_log_evidence_fn,
      rejuvenation = rejuvenation,
      rejuvenate_every = rejuvenate_every,
      rejuvenate_after_resample = rejuvenate_after_resample,
      n_rejuvenation_moves = n_rejuvenation_moves,
      rejuvenation_proposal = rejuvenation_proposal,
      rejuvenation_indep_prob = rejuvenation_indep_prob,
      rejuvenation_rw_scale = rejuvenation_rw_scale,
      rejuvenation_indep_scale = rejuvenation_indep_scale,
      proposal_transform = proposal_transform,
      auxiliary_proposal = auxiliary_proposal,
      auxiliary_stickiness = auxiliary_stickiness,
      auxiliary_correlation = auxiliary_correlation,
      warm_start_mode = warm_start_mode,
      bridge_schedule_lookup = setNames(
        replicate(length(subject_order), bridge_schedule, simplify = FALSE),
        as.character(subject_order)
      ),
      surrogate_min_points = surrogate_min_points,
      surrogate_ridge = surrogate_ridge,
      surrogate_degree = surrogate_degree,
      surrogate_neighbors = surrogate_neighbors,
      theta_names = theta_names,
      local_n_cores = local_n_cores,
      outer_n_cores = outer_n_cores,
      base_seed = base_seed,
      cache_store = cache_store,
      store_local_fits = store_local_fits,
      verbose = verbose,
      local_smc_control = local_smc_control
    )
  } else {
    .run_collapsed_bridge_adaptive(
      data_list = data_list,
      subject_order = subject_order,
      loglik_fn = loglik_fn,
      prior = prior,
      gaussian_map_fn = gaussian_map_fn,
      d_theta = d_theta,
      N = N,
      local_particles = local_particles,
      resample_threshold = resample_threshold,
      local_log_evidence_fn = local_log_evidence_fn,
      rejuvenation = rejuvenation,
      rejuvenate_every = rejuvenate_every,
      rejuvenate_after_resample = rejuvenate_after_resample,
      n_rejuvenation_moves = n_rejuvenation_moves,
      rejuvenation_proposal = rejuvenation_proposal,
      rejuvenation_indep_prob = rejuvenation_indep_prob,
      rejuvenation_rw_scale = rejuvenation_rw_scale,
      rejuvenation_indep_scale = rejuvenation_indep_scale,
      proposal_transform = proposal_transform,
      auxiliary_proposal = auxiliary_proposal,
      auxiliary_stickiness = auxiliary_stickiness,
      auxiliary_correlation = auxiliary_correlation,
      warm_start_mode = warm_start_mode,
      surrogate_min_points = surrogate_min_points,
      surrogate_ridge = surrogate_ridge,
      surrogate_degree = surrogate_degree,
      surrogate_neighbors = surrogate_neighbors,
      theta_names = theta_names,
      local_n_cores = local_n_cores,
      outer_n_cores = outer_n_cores,
      base_seed = base_seed,
      cache_store = cache_store,
      store_local_fits = store_local_fits,
      adaptive_ess_target = adaptive_ess_target,
      adaptive_lambda_tol = adaptive_lambda_tol,
      adaptive_max_steps_per_subject = adaptive_max_steps_per_subject,
      adaptive_max_bisect = adaptive_max_bisect,
      verbose = verbose,
      local_smc_control = local_smc_control
    )
  }

  list(
    phi = fit$phi,
    w = fit$w,
    log_evidence = fit$log_evidence,
    subject_loglik = fit$subject_loglik,
    subject_aux = fit$subject_aux,
    meta = list(
      method = "collapsed_bridge",
      subject_order = subject_order,
      bridge_schedule = bridge_schedule %||% fit$bridge_schedule,
      adaptive_bridge = isTRUE(adaptive_bridge),
      n_subjects = length(subject_order),
      n_stages = length(fit$stage_lambda),
      ess_hist = fit$ess_hist,
      stage_logZ = fit$stage_logZ,
      stage_subject = fit$stage_subject,
      stage_lambda = fit$stage_lambda,
      resampled = fit$resampled,
      rejuvenation = rejuvenation,
      rejuvenation_proposal = rejuvenation_proposal,
      proposal_transform = proposal_transform$name,
      auxiliary_proposal = auxiliary_proposal,
      auxiliary_stickiness = auxiliary_stickiness,
      auxiliary_correlation = auxiliary_correlation,
      warm_start_mode = warm_start_mode,
      surrogate = list(
        min_points = surrogate_min_points,
        ridge = surrogate_ridge,
        degree = surrogate_degree,
        neighbors = surrogate_neighbors
      ),
      rejuvenation_accept_hist = fit$rejuvenation_accept_hist,
      n_local_runs = fit$n_local_runs,
      evidence_cache = .collapsed_cache_stats(cache_store),
      extended_space = TRUE,
      exactness = .collapsed_exactness_meta(
        assume_unbiased_local_evidence = assume_unbiased_local_evidence,
        local_log_evidence_fn = local_log_evidence_fn,
        store_local_fits = cache_keep_local_fit,
        rejuvenation = rejuvenation,
        bridge_mode = if (isTRUE(adaptive_bridge)) "adaptive_repeated_partial" else "checkpoint_path",
        auxiliary_proposal = auxiliary_proposal,
        warm_start_mode = warm_start_mode,
        proposal_transform = proposal_transform
      )
    )
  )
}
