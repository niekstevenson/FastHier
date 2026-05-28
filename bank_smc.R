#!/usr/bin/env Rscript
# ============================================================================
# Local SMC bank nodes
# - Self-contained local proposal nodes for amortized hierarchical factors
# ============================================================================

if (!exists("%||%", mode = "function") ||
    !exists(".rowLogSumExp", mode = "function") ||
    !exists("logsumexp", mode = "function") ||
    !exists("ll_parallel", mode = "function") ||
    !exists("ESS", mode = "function") ||
    !exists("weighted_cov", mode = "function") ||
    !exists("regularize_cov", mode = "function")) {
  source("smc_core.R")
}
if (!exists("make_reference_prior_gaussian", mode = "function") ||
    !exists("reference_prior_geometry", mode = "function")) {
  source("reference_priors.R")
}
if (!exists("run_tempered_smc", mode = "function")) {
  source("SMC_super_fast.R")
}
if (!exists("build_population_factor_set", mode = "function") ||
    !exists("outer_population_smc", mode = "function") ||
    !exists("population_sufficient_stats_from_alpha", mode = "function")) {
  source("outer_population_smc.R")
}
if (!exists("normalize_population_model", mode = "function") ||
    !exists("population_model_log_alpha_given_theta", mode = "function") ||
    !exists("population_model_log_alpha_given_theta_many", mode = "function")) {
  source("population_models.R")
}

.bank_normalize_weights <- function(w) {
  w <- pmax(as.numeric(w), 0)
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) {
    rep(1 / length(w), length(w))
  } else {
    w / sw
  }
}

.bank_logsumexp2 <- function(a, b) {
  m <- pmax(a, b)
  both_zero <- !is.finite(a) & !is.finite(b)
  out <- m + log(exp(a - m) + exp(b - m))
  out[both_zero] <- -Inf
  out
}

.bank_mean_or_na <- function(x) {
  x <- as.numeric(x)
  if (!length(x) || all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

.bank_assert_smc_complete <- function(fit, label, tol = 1e-10) {
  final_lambda <- as.numeric(fit$final_lambda %||% NA_real_)
  if (!is.finite(final_lambda) || final_lambda < 1 - tol) {
    stop(label, " SMC stopped before lambda_target; log_marginal_anchor is not a completed local normalizer.")
  }
  invisible(TRUE)
}

.bank_align_alpha <- function(alpha, population_model) {
  model <- normalize_population_model(population_model)
  alpha <- as.matrix(alpha)
  if (ncol(alpha) != model$alpha_dim) {
    stop("Bank node alpha dimension does not match the population model.")
  }
  if (!is.null(colnames(alpha)) && setequal(colnames(alpha), model$alpha_names)) {
    alpha <- alpha[, model$alpha_names, drop = FALSE]
  } else {
    colnames(alpha) <- model$alpha_names
  }
  alpha
}

.bank_align_theta <- function(theta, population_model) {
  model <- normalize_population_model(population_model)
  theta <- .as_hyper_matrix(theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  if (nrow(theta) != 1L) {
    stop("Bank node theta_anchor must contain exactly one theta row.")
  }
  theta
}

new_bank_smc_node <- function(local_id,
                              theta_anchor,
                              alpha,
                              weights,
                              log_likelihood,
                              log_prior_anchor,
                              log_marginal_anchor,
                              mcse_log_marginal_anchor = NA_real_,
                              diagnostics = list(),
                              bridge_provenance = list(type = "initial"),
                              alpha_names = colnames(alpha),
                              hyper_names = colnames(theta_anchor)) {
  alpha <- as.matrix(alpha)
  theta_anchor <- as.matrix(theta_anchor)
  n <- nrow(alpha)
  if (n < 1L) {
    stop("Bank node alpha must contain at least one particle.")
  }
  weights <- .bank_normalize_weights(weights)
  log_likelihood <- as.numeric(log_likelihood)
  log_prior_anchor <- as.numeric(log_prior_anchor)

  if (length(weights) != n ||
      length(log_likelihood) != n ||
      length(log_prior_anchor) != n) {
    stop("Bank node particle fields must all have length nrow(alpha).")
  }
  if (!is.finite(log_marginal_anchor)) {
    stop("Bank node log_marginal_anchor must be finite.")
  }

  colnames(alpha) <- alpha_names %||% colnames(alpha)
  colnames(theta_anchor) <- hyper_names %||% colnames(theta_anchor)

  structure(
    list(
      local_id = as.integer(local_id %||% NA_integer_),
      theta_anchor = theta_anchor,
      alpha = alpha,
      weights = weights,
      log_likelihood = log_likelihood,
      log_prior_anchor = log_prior_anchor,
      log_marginal_anchor = as.numeric(log_marginal_anchor),
      mcse_log_marginal_anchor = as.numeric(mcse_log_marginal_anchor %||% NA_real_),
      diagnostics = diagnostics %||% list(),
      bridge_provenance = bridge_provenance %||% list(type = "initial"),
      alpha_names = as.character(colnames(alpha)),
      hyper_names = as.character(colnames(theta_anchor)),
      n_particles = as.integer(n)
    ),
    class = "bank_smc_node"
  )
}

validate_bank_smc_node <- function(node, population_model = NULL) {
  if (!inherits(node, "bank_smc_node")) {
    stop("node must inherit from 'bank_smc_node'.")
  }
  required <- c(
    "local_id",
    "theta_anchor",
    "alpha",
    "weights",
    "log_likelihood",
    "log_prior_anchor",
    "log_marginal_anchor",
    "diagnostics",
    "bridge_provenance"
  )
  missing <- setdiff(required, names(node))
  if (length(missing)) {
    stop("Bank node is missing: ", paste(missing, collapse = ", "))
  }

  if (!is.null(population_model)) {
    model <- normalize_population_model(population_model)
    node$alpha <- .bank_align_alpha(node$alpha, model)
    node$theta_anchor <- .bank_align_theta(node$theta_anchor, model)
  } else {
    node$alpha <- as.matrix(node$alpha)
    node$theta_anchor <- as.matrix(node$theta_anchor)
  }

  n <- nrow(node$alpha)
  if (n < 1L) {
    stop("Bank node alpha must contain at least one particle.")
  }
  node$weights <- .bank_normalize_weights(node$weights)
  node$log_likelihood <- as.numeric(node$log_likelihood)
  node$log_prior_anchor <- as.numeric(node$log_prior_anchor)
  if (length(node$weights) != n ||
      length(node$log_likelihood) != n ||
      length(node$log_prior_anchor) != n) {
    stop("Bank node particle fields must all have length nrow(alpha).")
  }
  if (!is.finite(node$log_marginal_anchor)) {
    stop("Bank node log_marginal_anchor must be finite.")
  }

  node$alpha_names <- as.character(colnames(node$alpha))
  node$hyper_names <- as.character(colnames(node$theta_anchor))
  node$n_particles <- as.integer(n)
  node
}

build_bank_smc_node_from_local_fit <- function(local_fit,
                                               theta_anchor,
                                               population_model,
                                               local_id = NULL,
                                               diagnostics = list(),
                                               bridge_provenance = list(type = "initial")) {
  model <- normalize_population_model(population_model)
  if (is.null(local_fit$Theta) || is.null(local_fit$w)) {
    stop("local_fit must contain Theta and w.")
  }
  if (is.null(local_fit$loglik)) {
    stop("local_fit must contain loglik for bank node construction.")
  }

  alpha <- .bank_align_alpha(local_fit$Theta, model)
  theta_anchor <- .bank_align_theta(theta_anchor, model)
  log_prior_anchor <- population_model_log_alpha_given_theta(
    model,
    alpha = alpha,
    theta = theta_anchor
  )

  fit_diagnostics <- local_fit$meta %||% list()
  fit_diagnostics$source <- fit_diagnostics$source %||% "local_smc_fit"
  fit_diagnostics$final_lambda <- local_fit$final_lambda %||% NA_real_
  fit_diagnostics$rounds <- local_fit$meta$rounds %||% NA_integer_

  new_bank_smc_node(
    local_id = local_id %||% local_fit$local_id %||% NA_integer_,
    theta_anchor = theta_anchor,
    alpha = alpha,
    weights = local_fit$w,
    log_likelihood = local_fit$loglik,
    log_prior_anchor = log_prior_anchor,
    log_marginal_anchor = local_fit$log_evidence,
    mcse_log_marginal_anchor = local_fit$mcse_logZ %||% NA_real_,
    diagnostics = modifyList(fit_diagnostics, diagnostics),
    bridge_provenance = bridge_provenance,
    alpha_names = model$alpha_names,
    hyper_names = model$hyper_names
  )
}

new_bank_smc_local_bank <- function(local_id,
                                    nodes,
                                    eta = NULL,
                                    max_nodes = Inf,
                                    max_particles = Inf,
                                    support_diagnostics = list(),
                                    stack = NULL) {
  if (!is.list(nodes) || !length(nodes)) {
    stop("Local bank nodes must be a non-empty list.")
  }
  if (!all(vapply(nodes, inherits, logical(1), what = "bank_smc_node"))) {
    stop("Every local bank node must inherit from 'bank_smc_node'.")
  }

  eta <- .bank_normalize_weights(eta %||% rep(1, length(nodes)))
  if (length(eta) != length(nodes)) {
    stop("Local bank eta must match the number of nodes.")
  }

  structure(
    list(
      local_id = as.integer(local_id %||% NA_integer_),
      nodes = nodes,
      eta = eta,
      max_nodes = as.numeric(max_nodes),
      max_particles = as.numeric(max_particles),
      support_diagnostics = support_diagnostics %||% list(),
      stack = stack,
      n_nodes = as.integer(length(nodes)),
      n_particles = as.integer(sum(vapply(nodes, `[[`, integer(1), "n_particles")))
    ),
    class = "bank_smc_local_bank"
  )
}

validate_bank_smc_local_bank <- function(bank, population_model = NULL) {
  if (!inherits(bank, "bank_smc_local_bank")) {
    stop("bank must inherit from 'bank_smc_local_bank'.")
  }
  required <- c(
    "local_id",
    "nodes",
    "eta",
    "max_nodes",
    "max_particles",
    "support_diagnostics"
  )
  missing <- setdiff(required, names(bank))
  if (length(missing)) {
    stop("Local bank is missing: ", paste(missing, collapse = ", "))
  }
  if (!is.list(bank$nodes) || !length(bank$nodes)) {
    stop("Local bank nodes must be a non-empty list.")
  }

  bank$nodes <- lapply(bank$nodes, validate_bank_smc_node, population_model = population_model)
  node_ids <- vapply(bank$nodes, function(node) as.integer(node$local_id %||% NA_integer_), integer(1))
  finite_node_ids <- node_ids[is.finite(node_ids)]
  if (is.finite(bank$local_id) && length(finite_node_ids) && any(finite_node_ids != bank$local_id)) {
    stop("Local bank node local_id values do not match the bank local_id.")
  }
  if (!is.finite(bank$local_id) && length(finite_node_ids)) {
    bank$local_id <- finite_node_ids[1L]
  }

  bank$eta <- .bank_normalize_weights(bank$eta)
  if (length(bank$eta) != length(bank$nodes)) {
    stop("Local bank eta must match the number of nodes.")
  }

  bank$max_nodes <- as.numeric(bank$max_nodes)
  bank$max_particles <- as.numeric(bank$max_particles)
  bank$n_nodes <- as.integer(length(bank$nodes))
  bank$n_particles <- as.integer(sum(vapply(bank$nodes, `[[`, integer(1), "n_particles")))

  if (is.finite(bank$max_nodes) && bank$n_nodes > bank$max_nodes) {
    stop("Local bank exceeds max_nodes.")
  }
  if (is.finite(bank$max_particles) && bank$n_particles > bank$max_particles) {
    stop("Local bank exceeds max_particles.")
  }
  if (is.null(bank$support_diagnostics)) {
    bank$support_diagnostics <- list()
  }
  bank
}

bank_smc_local_bank_theta_anchors <- function(bank, population_model = NULL) {
  bank <- validate_bank_smc_local_bank(bank, population_model = population_model)
  do.call(rbind, lapply(bank$nodes, `[[`, "theta_anchor"))
}

bank_smc_local_bank_add_node <- function(bank,
                                         node,
                                         eta = NULL,
                                         population_model = NULL) {
  bank <- validate_bank_smc_local_bank(bank, population_model = population_model)
  node <- validate_bank_smc_node(node, population_model = population_model)

  if (is.finite(bank$local_id) && is.finite(node$local_id) && node$local_id != bank$local_id) {
    stop("Added node local_id does not match the bank local_id.")
  }
  if (!is.finite(node$local_id)) node$local_id <- bank$local_id

  nodes <- c(bank$nodes, list(node))
  eta <- if (is.null(eta)) {
    c(bank$eta, mean(bank$eta))
  } else {
    c(bank$eta, as.numeric(eta))
  }

  validate_bank_smc_local_bank(
    new_bank_smc_local_bank(
      local_id = bank$local_id,
      nodes = nodes,
      eta = eta,
      max_nodes = bank$max_nodes,
      max_particles = bank$max_particles,
      support_diagnostics = bank$support_diagnostics,
      stack = NULL
    ),
    population_model = population_model
  )
}

bank_smc_calibrate_local_bank_normalizers <- function(bank,
                                                      population_model,
                                                      max_iter = 200L,
                                                      tol = 1e-8,
                                                      anchor = c("mean_smc", "first_smc")) {
  model <- normalize_population_model(population_model)
  bank <- validate_bank_smc_local_bank(bank, population_model = model)
  anchor <- match.arg(anchor)
  S <- length(bank$nodes)
  if (S <= 1L) {
    bank$support_diagnostics$normalizer_calibration <- list(
      method = "reverse_logistic_mbar",
      converged = TRUE,
      iterations = 0L,
      max_abs_delta = 0,
      max_abs_shift_from_smc = 0
    )
    return(bank)
  }

  theta_anchors <- bank_smc_local_bank_theta_anchors(bank, model)
  log_eta <- ifelse(bank$eta > 0, log(bank$eta), -Inf)
  logZ_smc <- vapply(bank$nodes, `[[`, numeric(1), "log_marginal_anchor")
  if (any(!is.finite(logZ_smc))) {
    stop("Cannot calibrate a local bank with non-finite node normalizers.")
  }

  logp_parts <- vector("list", S)
  log_sample_weight_parts <- vector("list", S)
  for (s in seq_len(S)) {
    node <- bank$nodes[[s]]
    logp_parts[[s]] <- t(population_model_log_alpha_given_theta_many(
      model,
      alpha = node$alpha,
      theta = theta_anchors
    ))
    log_sample_weight_parts[[s]] <- log_eta[s] + ifelse(node$weights > 0, log(node$weights), -Inf)
  }
  logp <- do.call(rbind, logp_parts)
  log_sample_weight <- unlist(log_sample_weight_parts, use.names = FALSE)
  if (nrow(logp) != length(log_sample_weight) || ncol(logp) != S) {
    stop("Internal normalizer calibration dimensions are inconsistent.")
  }

  ref <- logZ_smc[1L]
  logZ_rel <- logZ_smc - ref
  max_delta <- Inf
  iter <- 0L
  for (iter in seq_len(as.integer(max_iter))) {
    log_den <- .rowLogSumExp(sweep(sweep(logp, 2L, log_eta, "+"), 2L, logZ_rel, "-"))
    logZ_new <- vapply(
      seq_len(S),
      function(j) logsumexp(log_sample_weight + logp[, j] - log_den),
      numeric(1)
    )
    logZ_new <- logZ_new - logZ_new[1L]
    max_delta <- max(abs(logZ_new - logZ_rel))
    logZ_rel <- logZ_new
    if (is.finite(max_delta) && max_delta <= as.numeric(tol)) break
  }

  shift <- switch(
    anchor,
    mean_smc = sum(bank$eta * (logZ_smc - logZ_rel)),
    first_smc = logZ_smc[1L] - logZ_rel[1L]
  )
  logZ_calibrated <- logZ_rel + shift
  for (s in seq_len(S)) {
    bank$nodes[[s]]$diagnostics$raw_log_marginal_anchor <- bank$nodes[[s]]$diagnostics$raw_log_marginal_anchor %||%
      bank$nodes[[s]]$log_marginal_anchor
    bank$nodes[[s]]$diagnostics$normalizer_calibration_shift <- as.numeric(logZ_calibrated[s] - logZ_smc[s])
    bank$nodes[[s]]$log_marginal_anchor <- as.numeric(logZ_calibrated[s])
  }
  bank$stack <- NULL
  bank$support_diagnostics$normalizer_calibration <- list(
    method = "reverse_logistic_mbar",
    anchor = anchor,
    converged = is.finite(max_delta) && max_delta <= as.numeric(tol),
    iterations = as.integer(iter),
    max_abs_delta = as.numeric(max_delta),
    max_abs_shift_from_smc = as.numeric(max(abs(logZ_calibrated - logZ_smc))),
    mean_abs_shift_from_smc = as.numeric(mean(abs(logZ_calibrated - logZ_smc)))
  )
  validate_bank_smc_local_bank(bank, population_model = model)
}

bank_smc_calibrate_banks_normalizers <- function(banks,
                                                 population_model,
                                                 max_iter = 200L,
                                                 tol = 1e-8,
                                                 anchor = "mean_smc",
                                                 n_jobs = 1L) {
  model <- normalize_population_model(population_model)
  out <- parallel::mclapply(
    banks,
    bank_smc_calibrate_local_bank_normalizers,
    population_model = model,
    max_iter = max_iter,
    tol = tol,
    anchor = anchor,
    mc.cores = as.integer(max(1L, n_jobs))
  )
  names(out) <- names(banks)
  out
}

bank_smc_local_bank_build_stack <- function(bank, population_model) {
  model <- normalize_population_model(population_model)
  bank <- validate_bank_smc_local_bank(bank, population_model = model)
  theta_anchors <- bank_smc_local_bank_theta_anchors(bank, model)
  log_eta <- ifelse(bank$eta > 0, log(bank$eta), -Inf)
  log_marginal_anchor <- vapply(bank$nodes, `[[`, numeric(1), "log_marginal_anchor")
  anchor_offset <- log_eta - log_marginal_anchor

  alpha_parts <- vector("list", length(bank$nodes))
  log_base_parts <- vector("list", length(bank$nodes))
  log_denominator_parts <- vector("list", length(bank$nodes))
  node_index_parts <- vector("list", length(bank$nodes))

  for (b in seq_along(bank$nodes)) {
    node <- bank$nodes[[b]]
    alpha <- node$alpha
    anchor_log_prior <- population_model_log_alpha_given_theta_many(
      model,
      alpha = alpha,
      theta = theta_anchors
    )
    log_denominator <- .rowLogSumExp(t(sweep(anchor_log_prior, 1L, anchor_offset, "+")))
    if (any(!is.finite(log_denominator))) {
      stop("Bank deterministic-mixture denominator is non-finite for local ", bank$local_id, ".")
    }

    log_node_weight <- log_eta[b] + ifelse(node$weights > 0, log(node$weights), -Inf)
    alpha_parts[[b]] <- alpha
    log_base_parts[[b]] <- as.numeric(log_node_weight - log_denominator)
    log_denominator_parts[[b]] <- as.numeric(log_denominator)
    node_index_parts[[b]] <- rep.int(b, nrow(alpha))
  }

  alpha <- do.call(rbind, alpha_parts)
  colnames(alpha) <- model$alpha_names
  structure(
    list(
      alpha = alpha,
      log_base = unlist(log_base_parts, use.names = FALSE),
      log_denominator = unlist(log_denominator_parts, use.names = FALSE),
      node_index = unlist(node_index_parts, use.names = FALSE),
      theta_anchors = theta_anchors,
      eta = bank$eta,
      local_id = bank$local_id,
      n_nodes = bank$n_nodes,
      n_particles = nrow(alpha)
    ),
    class = "bank_smc_local_bank_stack"
  )
}

.bank_smc_stack_usable <- function(stack, bank) {
  inherits(stack, "bank_smc_local_bank_stack") &&
    !is.null(stack$alpha) &&
    !is.null(stack$log_base) &&
    nrow(stack$alpha) == bank$n_particles &&
    length(stack$log_base) == bank$n_particles
}

.bank_smc_local_bank_log_marginal_from_stack <- function(stack,
                                                        population_model,
                                                        theta,
                                                        block_size = 1024L) {
  model <- normalize_population_model(population_model)
  theta <- .as_hyper_matrix(theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  block_size <- as.integer(max(1L, block_size))
  accum <- rep.int(-Inf, nrow(theta))

  for (start in seq.int(1L, nrow(stack$alpha), by = block_size)) {
    idx <- seq.int(start, min(start + block_size - 1L, nrow(stack$alpha)))
    log_prior <- population_model_log_alpha_given_theta_many(
      model,
      alpha = stack$alpha[idx, , drop = FALSE],
      theta = theta
    )
    block_lse <- .rowLogSumExp(sweep(log_prior, 2L, stack$log_base[idx], "+"))
    accum <- .bank_logsumexp2(accum, block_lse)
  }
  as.numeric(accum)
}

bank_smc_local_bank_log_marginal_many <- function(bank,
                                                  population_model,
                                                  theta,
                                                  block_size = 1024L,
                                                  use_cached_stack = TRUE) {
  bank <- validate_bank_smc_local_bank(bank, population_model = population_model)
  stack <- if (isTRUE(use_cached_stack) && .bank_smc_stack_usable(bank$stack, bank)) {
    bank$stack
  } else {
    bank_smc_local_bank_build_stack(bank, population_model)
  }
  .bank_smc_local_bank_log_marginal_from_stack(
    stack = stack,
    population_model = population_model,
    theta = theta,
    block_size = block_size
  )
}

bank_smc_local_bank_log_marginal <- function(bank,
                                             population_model,
                                             theta,
                                             block_size = 1024L,
                                             use_cached_stack = TRUE) {
  bank_smc_local_bank_log_marginal_many(
    bank = bank,
    population_model = population_model,
    theta = theta,
    block_size = block_size,
    use_cached_stack = use_cached_stack
  )[1L]
}

bank_smc_local_bank_log_terms <- function(bank,
                                          population_model,
                                          theta,
                                          use_cached_stack = TRUE) {
  model <- normalize_population_model(population_model)
  bank <- validate_bank_smc_local_bank(bank, population_model = model)
  theta <- .as_hyper_matrix(theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  stack <- if (isTRUE(use_cached_stack) && .bank_smc_stack_usable(bank$stack, bank)) {
    bank$stack
  } else {
    bank_smc_local_bank_build_stack(bank, model)
  }
  log_prior <- population_model_log_alpha_given_theta_many(
    model,
    alpha = stack$alpha,
    theta = theta
  )
  log_terms <- sweep(log_prior, 2L, stack$log_base, "+")
  list(
    theta = theta,
    log_terms = log_terms,
    log_marginal = as.numeric(.rowLogSumExp(log_terms)),
    stack = stack
  )
}

.bank_gaussian_prior_sufficient_mean <- function(population_model, theta) {
  model <- normalize_population_model(population_model)
  if (!identical(model$fast_family %||% NULL, "gaussian")) return(NULL)
  prepared <- population_model_prepare_theta(model, theta)
  if (!identical(prepared$family, "gaussian") ||
      !identical(prepared$quadratic_kind, "diag")) {
    return(NULL)
  }
  sigma2 <- 1 / pmax(prepared$quadratic_coef, 1e-12)
  out <- cbind(prepared$mean, prepared$mean * prepared$mean + sigma2)
  colnames(out) <- c(
    paste0("linear_", model$alpha_names),
    paste0("quadratic_", model$alpha_names)
  )
  out
}

bank_smc_local_bank_score_summary <- function(bank,
                                              population_model,
                                              theta,
                                              omit_node = NULL,
                                              use_cached_stack = TRUE) {
  model <- normalize_population_model(population_model)
  bank <- validate_bank_smc_local_bank(bank, population_model = model)
  theta <- .bank_align_theta(theta, model)
  if (!is.null(omit_node)) {
    omit_node <- as.integer(omit_node)
    if (length(omit_node) != 1L || omit_node < 1L || omit_node > bank$n_nodes) {
      stop("omit_node is out of range.")
    }
    keep <- setdiff(seq_along(bank$nodes), omit_node)
    if (!length(keep)) stop("Cannot score a local bank with every node omitted.")
    bank <- new_bank_smc_local_bank(
      local_id = bank$local_id,
      nodes = bank$nodes[keep],
      eta = bank$eta[keep],
      max_nodes = bank$max_nodes,
      max_particles = bank$max_particles
    )
    use_cached_stack <- FALSE
  }
  terms <- bank_smc_local_bank_log_terms(
    bank = bank,
    population_model = model,
    theta = theta,
    use_cached_stack = use_cached_stack
  )
  lw <- as.numeric(terms$log_terms[1L, ])
  log_marginal <- as.numeric(terms$log_marginal[1L])
  w <- exp(lw - log_marginal)
  out <- list(log_marginal = log_marginal)
  prior_mean <- .bank_gaussian_prior_sufficient_mean(model, theta)
  if (!is.null(prior_mean)) {
    Talpha <- .bank_gaussian_sufficient_stat_matrix(terms$stack$alpha, model)
    Et <- as.numeric(colSums(Talpha * w))
    score_eta <- Et - as.numeric(prior_mean[1L, ])
    names(score_eta) <- colnames(Talpha)
    out$score_eta <- score_eta
  }
  out
}

.bank_surface_local_stability <- function(bank,
                                          population_model,
                                          theta,
                                          use_cached_stack = TRUE) {
  model <- normalize_population_model(population_model)
  bank <- validate_bank_smc_local_bank(bank, population_model = model)
  theta <- .bank_align_theta(theta, model)
  full <- bank_smc_local_bank_score_summary(
    bank = bank,
    population_model = model,
    theta = theta,
    use_cached_stack = use_cached_stack
  )
  anchors <- bank_smc_local_bank_theta_anchors(bank, model)
  natural <- bank_smc_local_bank_natural_diagnostic(bank, model, theta)
  nearest <- as.integer((natural %||% list())$best_node %||%
    which.min(rowSums(abs(sweep(anchors, 2L, as.numeric(theta), "-")))))
  loo <- if (bank$n_nodes > 1L) {
    bank_smc_local_bank_score_summary(
      bank = bank,
      population_model = model,
      theta = theta,
      omit_node = nearest,
      use_cached_stack = FALSE
    )
  } else {
    NULL
  }
  log_sensitivity <- if (!is.null(loo)) {
    abs(as.numeric(full$log_marginal) - as.numeric(loo$log_marginal))
  } else {
    Inf
  }
  score_sensitivity <- if (!is.null(loo) &&
      !is.null(full$score_eta) &&
      !is.null(loo$score_eta)) {
    sqrt(mean((full$score_eta - loo$score_eta)^2))
  } else {
    NA_real_
  }
  data.frame(
    local_id = as.integer(bank$local_id),
    nearest_node = nearest,
    log_sensitivity = as.numeric(log_sensitivity),
    score_sensitivity = as.numeric(score_sensitivity),
    natural_ess_frac = as.numeric((natural %||% list())$best_ess_frac %||% NA_real_),
    natural_distance = as.numeric((natural %||% list())$best_distance %||% NA_real_),
    log_marginal = as.numeric(full$log_marginal),
    check.names = FALSE
  )
}

bank_smc_score_surface_theta_design_with_banks <- function(theta,
                                                           banks,
                                                           population_model,
                                                           local_count = 20L,
                                                           log_weight = 1,
                                                           score_weight = 0.25,
                                                           distance_weight = 0.05,
                                                           n_jobs = 1L) {
  model <- normalize_population_model(population_model)
  theta <- .as_hyper_matrix(theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  if (is.null(banks) || !length(banks) || !nrow(theta)) {
    return(list(score = rep(0, nrow(theta)), diagnostics = data.frame()))
  }
  n_score <- min(length(banks), as.integer(max(1L, local_count)))
  local_ids <- unique(as.integer(round(seq(1L, length(banks), length.out = n_score))))
  parts <- parallel::mclapply(
    seq_len(nrow(theta)),
    function(j) {
      local_rows <- lapply(local_ids, function(i) {
        row <- .bank_surface_local_stability(
          bank = banks[[i]],
          population_model = model,
          theta = theta[j, , drop = FALSE],
          use_cached_stack = TRUE
        )
        row$theta_id <- as.integer(j)
        row
      })
      do.call(rbind, local_rows)
    },
    mc.cores = as.integer(max(1L, n_jobs))
  )
  diagnostics <- do.call(rbind, parts)
  diagnostics$component_score <-
    as.numeric(log_weight) * pmax(diagnostics$log_sensitivity, 0) +
    as.numeric(score_weight) * pmax(diagnostics$score_sensitivity, 0) +
    as.numeric(distance_weight) * pmax(diagnostics$natural_distance, 0)
  score <- vapply(seq_len(nrow(theta)), function(j) {
    x <- diagnostics$component_score[diagnostics$theta_id == j]
    x <- x[is.finite(x)]
    if (!length(x)) return(0)
    sum(x) + max(x)
  }, numeric(1))
  list(score = score, diagnostics = diagnostics)
}

.bank_log_weight_summary <- function(log_weights, use_psis = TRUE) {
  log_weights <- as.numeric(log_weights)
  finite <- is.finite(log_weights)
  n_total <- length(log_weights)
  n_finite <- sum(finite)

  if (!n_total || !n_finite) {
    return(list(
      ess = 0,
      ess_frac = 0,
      max_weight = NA_real_,
      q99_weight = NA_real_,
      log_weight_var = Inf,
      pareto_k = NA_real_,
      n_finite = as.integer(n_finite)
    ))
  }

  lw <- log_weights[finite]
  lse <- logsumexp(lw)
  w <- exp(lw - lse)
  ess <- 1 / sum(w * w)
  lw_center <- sum(w * lw)

  pareto_k <- NA_real_
  if (isTRUE(use_psis) && requireNamespace("loo", quietly = TRUE) && length(lw) > 2L) {
    psis <- tryCatch(
      suppressWarnings(loo::psis(matrix(lw - max(lw), ncol = 1L))),
      error = function(e) NULL
    )
    if (!is.null(psis)) {
      pareto_k <- as.numeric(loo::pareto_k_values(psis)[1L])
    }
  }

  list(
    ess = as.numeric(ess),
    ess_frac = as.numeric(ess / n_total),
    max_weight = as.numeric(max(w)),
    q99_weight = as.numeric(stats::quantile(w, 0.99, names = FALSE, type = 8)),
    log_weight_var = as.numeric(sum(w * (lw - lw_center)^2)),
    pareto_k = pareto_k,
    n_finite = as.integer(n_finite)
  )
}

.bank_single_node_log_weights <- function(node, population_model, theta) {
  model <- normalize_population_model(population_model)
  node <- validate_bank_smc_node(node, population_model = model)
  theta <- .bank_align_theta(theta, model)
  log_prior_theta <- population_model_log_alpha_given_theta(
    model,
    alpha = node$alpha,
    theta = theta
  )
  ifelse(node$weights > 0, log(node$weights), -Inf) +
    log_prior_theta -
    node$log_prior_anchor
}

.bank_weighted_var <- function(x, w) {
  x <- as.numeric(x)
  w <- .bank_normalize_weights(w)
  mu <- sum(w * x)
  as.numeric(sum(w * (x - mu)^2))
}

.bank_gaussian_natural_eta <- function(population_model, theta) {
  model <- normalize_population_model(population_model)
  if (!identical(model$fast_family %||% NULL, "gaussian")) {
    return(NULL)
  }
  prepared <- population_model_prepare_theta(model, theta)
  if (!identical(prepared$family, "gaussian") ||
      !identical(prepared$quadratic_kind, "diag")) {
    return(NULL)
  }
  out <- cbind(prepared$eta, -0.5 * prepared$quadratic_coef)
  colnames(out) <- c(
    paste0("linear_", model$alpha_names),
    paste0("quadratic_", model$alpha_names)
  )
  out
}

.bank_gaussian_sufficient_stat_matrix <- function(alpha, population_model) {
  model <- normalize_population_model(population_model)
  alpha <- .bank_align_alpha(alpha, model)
  out <- cbind(alpha, alpha * alpha)
  colnames(out) <- c(
    paste0("linear_", model$alpha_names),
    paste0("quadratic_", model$alpha_names)
  )
  out
}

bank_smc_local_bank_natural_diagnostic <- function(bank,
                                                   population_model,
                                                   theta) {
  model <- normalize_population_model(population_model)
  bank <- validate_bank_smc_local_bank(bank, population_model = model)
  theta <- .bank_align_theta(theta, model)
  eta_target <- .bank_gaussian_natural_eta(model, theta)
  if (is.null(eta_target)) return(NULL)
  eta_target <- as.numeric(eta_target[1L, ])

  rows <- lapply(seq_along(bank$nodes), function(i) {
    node <- bank$nodes[[i]]
    eta_anchor <- as.numeric(.bank_gaussian_natural_eta(model, node$theta_anchor)[1L, ])
    delta_eta <- eta_target - eta_anchor
    Talpha <- .bank_gaussian_sufficient_stat_matrix(node$alpha, model)
    log_ratio_projection <- as.numeric(Talpha %*% delta_eta)
    variance <- .bank_weighted_var(log_ratio_projection, node$weights)
    data.frame(
      node = as.integer(i),
      natural_variance = variance,
      natural_distance = sqrt(pmax(variance, 0)),
      natural_ess_frac = exp(-pmax(variance, 0)),
      check.names = FALSE
    )
  })

  single <- do.call(rbind, rows)
  best <- which.min(single$natural_variance)
  list(
    single = single,
    best_node = as.integer(single$node[best]),
    best_variance = as.numeric(single$natural_variance[best]),
    best_distance = as.numeric(single$natural_distance[best]),
    best_ess_frac = as.numeric(single$natural_ess_frac[best])
  )
}

bank_smc_local_bank_single_node_diagnostics <- function(bank,
                                                        population_model,
                                                        theta,
                                                        use_psis = TRUE) {
  model <- normalize_population_model(population_model)
  bank <- validate_bank_smc_local_bank(bank, population_model = model)
  theta <- .bank_align_theta(theta, model)

  rows <- lapply(seq_along(bank$nodes), function(i) {
    node <- bank$nodes[[i]]
    log_weights <- .bank_single_node_log_weights(node, model, theta)
    summary <- .bank_log_weight_summary(log_weights, use_psis = use_psis)
    log_marginal_single <- if (summary$n_finite) {
      node$log_marginal_anchor + logsumexp(log_weights[is.finite(log_weights)])
    } else {
      -Inf
    }
    data.frame(
      node = as.integer(i),
      ess = summary$ess,
      ess_frac = summary$ess_frac,
      max_weight = summary$max_weight,
      q99_weight = summary$q99_weight,
      log_weight_var = summary$log_weight_var,
      pareto_k = summary$pareto_k,
      n_finite = summary$n_finite,
      log_marginal_single = as.numeric(log_marginal_single)
    )
  })

  do.call(rbind, rows)
}

bank_smc_local_bank_mixture_diagnostic <- function(bank,
                                                   population_model,
                                                   theta,
                                                   use_psis = TRUE,
                                                   use_cached_stack = TRUE) {
  bank <- validate_bank_smc_local_bank(bank, population_model = population_model)
  terms <- bank_smc_local_bank_log_terms(
    bank = bank,
    population_model = population_model,
    theta = theta,
    use_cached_stack = use_cached_stack
  )
  if (nrow(terms$theta) != 1L) {
    stop("Bank mixture diagnostic expects exactly one theta row.")
  }
  summary <- .bank_log_weight_summary(terms$log_terms[1L, ], use_psis = use_psis)
  node_particle_scale <- max(1L, max(vapply(bank$nodes, `[[`, integer(1), "n_particles")))
  data.frame(
    ess = summary$ess,
    ess_frac = as.numeric(summary$ess / node_particle_scale),
    total_ess_frac = summary$ess_frac,
    max_weight = summary$max_weight,
    q99_weight = summary$q99_weight,
    log_weight_var = summary$log_weight_var,
    pareto_k = summary$pareto_k,
    n_finite = summary$n_finite,
    log_marginal = as.numeric(terms$log_marginal[1L])
  )
}

bank_smc_local_bank_coverage_diagnostic <- function(bank,
                                                    population_model,
                                                    theta,
                                                    min_mixture_ess_frac = 0.2,
                                                    min_single_ess_frac = 0.2,
                                                    max_pareto_k = 0.7,
                                                    use_psis = TRUE,
                                                    use_cached_stack = TRUE) {
  model <- normalize_population_model(population_model)
  bank <- validate_bank_smc_local_bank(bank, population_model = model)
  theta <- .bank_align_theta(theta, model)

  single <- bank_smc_local_bank_single_node_diagnostics(
    bank = bank,
    population_model = model,
    theta = theta,
    use_psis = use_psis
  )
  mixture <- bank_smc_local_bank_mixture_diagnostic(
    bank = bank,
    population_model = model,
    theta = theta,
    use_psis = use_psis,
    use_cached_stack = use_cached_stack
  )
  natural <- bank_smc_local_bank_natural_diagnostic(
    bank = bank,
    population_model = model,
    theta = theta
  )

  best <- which.max(single$ess_frac)
  pareto_ok <- is.na(mixture$pareto_k[1L]) || mixture$pareto_k[1L] <= max_pareto_k
  natural_threshold <- -log(pmax(as.numeric(min_mixture_ess_frac), .Machine$double.eps))
  covered <- mixture$ess_frac[1L] >= min_mixture_ess_frac && pareto_ok
  best_single_covered <- single$ess_frac[best] >= min_single_ess_frac

  list(
    theta = theta,
    covered = isTRUE(covered),
    best_node = as.integer(single$node[best]),
    best_single_covered = isTRUE(best_single_covered),
    mixture = mixture,
    single = single,
    natural = natural,
    thresholds = list(
      min_mixture_ess_frac = min_mixture_ess_frac,
      min_single_ess_frac = min_single_ess_frac,
      max_natural_variance = natural_threshold,
      max_pareto_k = max_pareto_k
    )
  )
}

.bank_gaussian_theta_from_natural_eta <- function(population_model, eta) {
  model <- normalize_population_model(population_model)
  if (!identical(model$fast_family %||% NULL, "gaussian")) return(NULL)
  eta <- as.matrix(eta)
  d <- model$alpha_dim
  if (ncol(eta) != 2L * d) {
    stop("Gaussian natural parameter matrix has incompatible dimension.")
  }
  linear <- eta[, seq_len(d), drop = FALSE]
  quadratic <- eta[, d + seq_len(d), drop = FALSE]
  sigma2 <- -0.5 / pmin(quadratic, -1e-12)
  mu <- linear * sigma2
  theta <- cbind(mu, log(sigma2))
  colnames(theta) <- model$hyper_names
  .as_hyper_matrix(theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
}

.bank_interpolate_theta_natural <- function(population_model, theta_a, theta_b, t = 0.5) {
  model <- normalize_population_model(population_model)
  theta_a <- .bank_align_theta(theta_a, model)
  theta_b <- .bank_align_theta(theta_b, model)
  eta_a <- .bank_gaussian_natural_eta(model, theta_a)
  eta_b <- .bank_gaussian_natural_eta(model, theta_b)
  if (!is.null(eta_a) && !is.null(eta_b)) {
    eta <- (1 - as.numeric(t)) * eta_a + as.numeric(t) * eta_b
    theta <- .bank_gaussian_theta_from_natural_eta(model, eta)
  } else {
    theta <- (1 - as.numeric(t)) * theta_a + as.numeric(t) * theta_b
    colnames(theta) <- model$hyper_names
  }
  .bank_align_theta(theta, model)
}

.bank_theta_anchor_index <- function(bank, population_model, theta, tol = 1e-8) {
  model <- normalize_population_model(population_model)
  bank <- validate_bank_smc_local_bank(bank, population_model = model)
  theta <- .bank_align_theta(theta, model)
  anchors <- bank_smc_local_bank_theta_anchors(bank, model)
  distance <- rowSums(abs(sweep(anchors, 2L, as.numeric(theta), "-")))
  hit <- which(distance <= as.numeric(tol))
  if (length(hit)) as.integer(hit[1L]) else NA_integer_
}

.bank_overlap_ok <- function(summary,
                             min_ess_frac,
                             max_pareto_k) {
  ess_ok <- is.finite(summary$ess_frac) && summary$ess_frac >= as.numeric(min_ess_frac)
  k <- as.numeric(summary$pareto_k %||% NA_real_)
  k_ok <- is.na(k) || !is.finite(k) || k <= as.numeric(max_pareto_k)
  isTRUE(ess_ok && k_ok)
}

.bank_worst_pareto_k <- function(a, b) {
  k <- as.numeric(c(a, b))
  k <- k[is.finite(k)]
  if (length(k)) max(k) else NA_real_
}

bank_smc_local_bank_overlap_graph <- function(bank,
                                              population_model,
                                              min_edge_ess_frac = 0.05,
                                              max_pareto_k = 0.7,
                                              use_psis = TRUE) {
  model <- normalize_population_model(population_model)
  bank <- validate_bank_smc_local_bank(bank, population_model = model)
  S <- bank$n_nodes
  adjacency <- matrix(FALSE, nrow = S, ncol = S)
  diag(adjacency) <- TRUE
  if (S <= 1L) {
    return(list(
      adjacency = adjacency,
      edges = data.frame(),
      components = rep.int(1L, S),
      n_components = 1L,
      certified = TRUE
    ))
  }

  edge_rows <- vector("list", S * (S - 1L) / 2L)
  row_id <- 0L
  for (a in seq_len(S - 1L)) {
    for (b in seq.int(a + 1L, S)) {
      forward <- .bank_log_weight_summary(
        .bank_single_node_log_weights(bank$nodes[[a]], model, bank$nodes[[b]]$theta_anchor),
        use_psis = use_psis
      )
      reverse <- .bank_log_weight_summary(
        .bank_single_node_log_weights(bank$nodes[[b]], model, bank$nodes[[a]]$theta_anchor),
        use_psis = use_psis
      )
      pass <- .bank_overlap_ok(forward, min_edge_ess_frac, max_pareto_k) &&
        .bank_overlap_ok(reverse, min_edge_ess_frac, max_pareto_k)
      adjacency[a, b] <- adjacency[b, a] <- isTRUE(pass)
      row_id <- row_id + 1L
      edge_rows[[row_id]] <- data.frame(
        node_a = as.integer(a),
        node_b = as.integer(b),
        pass = isTRUE(pass),
        forward_ess_frac = as.numeric(forward$ess_frac),
        reverse_ess_frac = as.numeric(reverse$ess_frac),
        worst_ess_frac = as.numeric(min(forward$ess_frac, reverse$ess_frac)),
        forward_pareto_k = as.numeric(forward$pareto_k %||% NA_real_),
        reverse_pareto_k = as.numeric(reverse$pareto_k %||% NA_real_),
        worst_pareto_k = .bank_worst_pareto_k(forward$pareto_k %||% NA_real_, reverse$pareto_k %||% NA_real_),
        check.names = FALSE
      )
    }
  }
  edges <- do.call(rbind, edge_rows)

  components <- rep.int(NA_integer_, S)
  component_id <- 0L
  for (start in seq_len(S)) {
    if (!is.na(components[start])) next
    component_id <- component_id + 1L
    queue <- start
    components[start] <- component_id
    while (length(queue)) {
      current <- queue[1L]
      queue <- queue[-1L]
      next_nodes <- which(adjacency[current, ] & is.na(components))
      if (length(next_nodes)) {
        components[next_nodes] <- component_id
        queue <- c(queue, next_nodes)
      }
    }
  }

  list(
    adjacency = adjacency,
    edges = edges,
    components = components,
    n_components = as.integer(component_id),
    certified = component_id == 1L
  )
}

.bank_overlap_cut_pairs <- function(graph, root_component = 1L) {
  edges <- graph$edges
  if (is.null(edges) || !nrow(edges)) return(list())
  components <- graph$components
  a_root <- components[edges$node_a] == root_component
  b_root <- components[edges$node_b] == root_component
  candidates <- edges[xor(a_root, b_root), , drop = FALSE]
  if (!nrow(candidates)) return(list())
  candidates <- candidates[order(-candidates$worst_ess_frac), , drop = FALSE]
  lapply(seq_len(nrow(candidates)), function(i) {
    row <- candidates[i, , drop = FALSE]
    a <- as.integer(row$node_a)
    b <- as.integer(row$node_b)
    if (components[a] == root_component) {
      list(source = a, target = b, edge = row)
    } else {
      list(source = b, target = a, edge = row)
    }
  })
}

bank_smc_anchor_node <- function(local_id,
                                 population_model,
                                 theta_anchor,
                                 data_i,
                                 loglik_fn,
                                 M,
                                 target_cess = 0.9,
                                 resample_threshold = 0.5,
                                 n_mcmc_moves = 2L,
                                 rw_scale = 0.75,
                                 G_mix = 8L,
                                 da_enable = TRUE,
                                 refit_every = 2L,
                                 max_steps = 64L,
                                 deterministic_resampling = FALSE,
                                 n_cores = 1L,
                                 seed = NULL,
                                 verbose = FALSE) {
  model <- normalize_population_model(population_model)
  theta_anchor <- .bank_align_theta(theta_anchor, model)
  M <- as.integer(M)
  if (M <= 0L) stop("Anchor particle count M must be positive.")

  reference_prior <- .bank_population_reference_prior_from_theta(
    model,
    theta_anchor,
    label = "bank_anchor_endpoint"
  )
  fit <- run_tempered_smc(
    reference_prior = reference_prior,
    bridge_stat_fn = function(alpha) ll_parallel(
      alpha,
      data_i,
      loglik_fn,
      n_cores = as.integer(n_cores)
    ),
    M = M,
    resample_threshold = resample_threshold,
    n_mcmc_moves = as.integer(n_mcmc_moves),
    max_rounds = as.integer(max_steps),
    lambda_target = 1,
    cess_target = as.numeric(target_cess),
    G_mix = as.integer(G_mix),
    refit_every = as.integer(refit_every),
    rw_scale_init = rw_scale,
    post_adapt_n_mcmc_moves = as.integer(n_mcmc_moves),
    da_enable = isTRUE(da_enable),
    deterministic_resampling = deterministic_resampling,
    n_cores = n_cores,
    seed = seed,
    verbose = verbose
  )

  .bank_assert_smc_complete(fit, "Bank anchor")

  node <- build_bank_smc_node_from_local_fit(
    local_fit = fit,
    theta_anchor = theta_anchor,
    population_model = model,
    local_id = local_id,
    diagnostics = list(
      source = "bank_anchor_endpoint",
      final_ess = as.numeric(ESS(fit$w)),
      final_ess_frac = as.numeric(ESS(fit$w) / length(fit$w)),
      mean_accept_rate = .bank_mean_or_na(fit$meta$accept_rate_hist)
    ),
    bridge_provenance = list(
      type = "endpoint_smc",
      target_theta = theta_anchor,
      target_cess = as.numeric(target_cess),
      resample_threshold = as.numeric(resample_threshold),
      n_mcmc_moves = as.integer(n_mcmc_moves),
      rw_scale = as.numeric(rw_scale),
      kernel = "run_tempered_smc"
    )
  )
  node
}

bank_smc_bridge_node <- function(local_id,
                                 population_model,
                                 source_node,
                                 theta_target,
                                 data_i,
                                 loglik_fn,
                                 M = NULL,
                                 target_cess = 0.9,
                                 resample_threshold = 0.5,
                                 n_mcmc_moves = 2L,
                                 rw_scale = 0.75,
                                 G_mix = 8L,
                                 da_enable = TRUE,
                                 refit_every = 2L,
                                 max_steps = 64L,
                                 deterministic_resampling = FALSE,
                                 n_cores = 1L,
                                 seed = NULL,
                                 verbose = FALSE) {
  model <- normalize_population_model(population_model)
  source_node <- validate_bank_smc_node(source_node, population_model = model)
  theta_target <- .bank_align_theta(theta_target, model)
  theta_source <- source_node$theta_anchor
  old_log_marginal <- as.numeric(source_node$log_marginal_anchor)
  if (!is.finite(old_log_marginal)) {
    stop("Cannot bridge from a source node without finite log_marginal_anchor.")
  }

  reference_prior <- .bank_population_reference_prior_from_theta(
    model,
    theta_target,
    label = "bank_bridge_endpoint"
  )
  base_logpdf_fn <- function(alpha) {
    ll_parallel(
      alpha,
      data_i,
      loglik_fn,
      n_cores = as.integer(n_cores)
    ) +
      population_model_log_alpha_given_theta(model, alpha = alpha, theta = theta_source) -
      old_log_marginal
  }
  bridge_stat_fn <- function(alpha) {
    population_model_log_alpha_given_theta(model, alpha = alpha, theta = theta_target) -
      population_model_log_alpha_given_theta(model, alpha = alpha, theta = theta_source)
  }

  fit <- run_tempered_smc(
    reference_prior = reference_prior,
    bridge_stat_fn = bridge_stat_fn,
    base_logpdf_fn = base_logpdf_fn,
    initial_particles = source_node$alpha,
    initial_weights = source_node$weights,
    initial_log_normalizer = old_log_marginal,
    M = as.integer(M %||% source_node$n_particles),
    resample_threshold = resample_threshold,
    n_mcmc_moves = as.integer(n_mcmc_moves),
    max_rounds = as.integer(max_steps),
    lambda_target = 1,
    cess_target = as.numeric(target_cess),
    G_mix = as.integer(G_mix),
    refit_every = as.integer(refit_every),
    rw_scale_init = rw_scale,
    post_adapt_n_mcmc_moves = as.integer(n_mcmc_moves),
    da_enable = isTRUE(da_enable),
    deterministic_resampling = deterministic_resampling,
    n_cores = n_cores,
    seed = seed,
    verbose = verbose
  )
  .bank_assert_smc_complete(fit, "Bank bridge")
  fit$loglik <- ll_parallel(
    fit$Theta,
    data_i,
    loglik_fn,
    n_cores = as.integer(n_cores)
  )

  build_bank_smc_node_from_local_fit(
    local_fit = fit,
    theta_anchor = theta_target,
    population_model = model,
    local_id = local_id,
    diagnostics = list(
      source = "bank_bridge_endpoint",
      final_ess = as.numeric(ESS(fit$w)),
      final_ess_frac = as.numeric(ESS(fit$w) / length(fit$w)),
      mean_accept_rate = .bank_mean_or_na(fit$meta$accept_rate_hist)
    ),
    bridge_provenance = list(
      type = "endpoint_bridge",
      source_theta = theta_source,
      target_theta = theta_target,
      source_log_marginal = old_log_marginal,
      target_cess = as.numeric(target_cess),
      resample_threshold = as.numeric(resample_threshold),
      n_mcmc_moves = as.integer(n_mcmc_moves),
      rw_scale = as.numeric(rw_scale),
      kernel = "run_tempered_smc"
    )
  )
}

bank_smc_local_bank_add_anchor <- function(bank,
                                           population_model,
                                           theta_target,
                                           data_i,
                                           loglik_fn,
                                           M = NULL,
                                           eta = NULL,
                                           bridge_min_single_ess_frac = 0.05,
                                           bridge_max_pareto_k = 0.7,
                                           ...,
                                           use_psis = TRUE) {
  model <- normalize_population_model(population_model)
  bank <- validate_bank_smc_local_bank(bank, population_model = model)
  theta_target <- .bank_align_theta(theta_target, model)

  coverage <- bank_smc_local_bank_coverage_diagnostic(
    bank = bank,
    population_model = model,
    theta = theta_target,
    use_psis = use_psis
  )
  source_diag <- bank_smc_local_bank_single_node_diagnostics(
    bank = bank,
    population_model = model,
    theta = theta_target,
    use_psis = TRUE
  )
  source_index <- which.max(source_diag$ess_frac)
  source_pareto_k <- as.numeric(source_diag$pareto_k[source_index])
  source_overlap_ok <- source_diag$ess_frac[source_index] >= as.numeric(bridge_min_single_ess_frac) &&
    (!is.finite(source_pareto_k) || source_pareto_k <= as.numeric(bridge_max_pareto_k))

  M <- as.integer(M %||% bank$nodes[[source_index]]$n_particles)
  node <- if (isTRUE(source_overlap_ok)) {
    bank_smc_bridge_node(
      local_id = bank$local_id,
      population_model = model,
      source_node = bank$nodes[[source_index]],
      theta_target = theta_target,
      data_i = data_i,
      loglik_fn = loglik_fn,
      M = M,
      ...
    )
  } else {
    bank_smc_anchor_node(
      local_id = bank$local_id,
      population_model = model,
      theta_anchor = theta_target,
      data_i = data_i,
      loglik_fn = loglik_fn,
      M = M,
      ...
    )
  }
  updated <- bank_smc_local_bank_add_node(
    bank = bank,
    node = node,
    eta = eta,
    population_model = model
  )

  list(
    bank = updated,
    node = node,
    start_node_index = as.integer(source_index),
    start_coverage = coverage,
    source_overlap = source_diag[source_index, , drop = FALSE],
    used_bridge = isTRUE(source_overlap_ok)
  )
}

bank_smc_local_bank_certify_bridge_graph <- function(bank,
                                                     population_model,
                                                     data_i,
                                                     loglik_fn,
                                                     M = NULL,
                                                     min_edge_ess_frac = 0.05,
                                                     max_pareto_k = 0.7,
                                                     use_psis = TRUE,
                                                     max_rounds = 50L,
                                                     target_cess = 0.9,
                                                     n_mcmc_moves = 2L,
                                                     n_cores = 1L,
                                                     seed = 123L,
                                                     verbose = FALSE,
                                                     ...) {
  model <- normalize_population_model(population_model)
  bank <- validate_bank_smc_local_bank(bank, population_model = model)
  max_rounds <- as.integer(max(0L, max_rounds))
  rows <- vector("list", 0L)

  graph <- bank_smc_local_bank_overlap_graph(
    bank = bank,
    population_model = model,
    min_edge_ess_frac = min_edge_ess_frac,
    max_pareto_k = max_pareto_k,
    use_psis = use_psis
  )
  if (isTRUE(graph$certified) || max_rounds <= 0L) {
    bank$support_diagnostics$overlap_graph <- list(
      method = "bidirectional_importance_overlap_graph",
      certified = isTRUE(graph$certified),
      n_components = as.integer(graph$n_components),
      rounds = 0L,
      inserted = 0L,
      min_edge_ess_frac = as.numeric(min_edge_ess_frac),
      max_pareto_k = as.numeric(max_pareto_k),
      use_psis = isTRUE(use_psis),
      edges = graph$edges
    )
    return(list(bank = bank, graph = graph, insertions = data.frame()))
  }

  for (round in seq_len(max_rounds)) {
    if (isTRUE(graph$certified)) break
    if (bank$n_nodes >= bank$max_nodes) break
    pair_candidates <- .bank_overlap_cut_pairs(graph, root_component = graph$components[1L])
    if (!length(pair_candidates)) break
    pair <- NULL
    theta_mid <- NULL
    for (candidate in pair_candidates) {
      source_node <- bank$nodes[[candidate$source]]
      target_node <- bank$nodes[[candidate$target]]
      candidate_mid <- .bank_interpolate_theta_natural(
        population_model = model,
        theta_a = source_node$theta_anchor,
        theta_b = target_node$theta_anchor,
        t = 0.5
      )
      if (is.na(.bank_theta_anchor_index(bank, model, candidate_mid, tol = 1e-8))) {
        pair <- candidate
        theta_mid <- candidate_mid
        break
      }
    }
    if (is.null(pair)) break

    source_node <- bank$nodes[[pair$source]]
    particles <- as.integer(M %||% source_node$n_particles)
    if (is.finite(bank$max_particles) && bank$n_particles + particles > bank$max_particles) break
    node <- bank_smc_bridge_node(
      local_id = bank$local_id,
      population_model = model,
      source_node = source_node,
      theta_target = theta_mid,
      data_i = data_i,
      loglik_fn = loglik_fn,
      M = particles,
      target_cess = target_cess,
      n_mcmc_moves = n_mcmc_moves,
      n_cores = n_cores,
      seed = as.integer(seed + round),
      verbose = verbose,
      ...
    )
    bank <- bank_smc_local_bank_add_node(
      bank = bank,
      node = node,
      eta = NULL,
      population_model = model
    )
    inserted_node <- bank$n_nodes

    rows[[length(rows) + 1L]] <- data.frame(
      local_id = as.integer(bank$local_id),
      round = as.integer(round),
      source_node = as.integer(pair$source),
      target_node = as.integer(pair$target),
      inserted_node = as.integer(inserted_node),
      pre_components = as.integer(graph$n_components),
      pre_worst_ess_frac = as.numeric(pair$edge$worst_ess_frac),
      pre_forward_ess_frac = as.numeric(pair$edge$forward_ess_frac),
      pre_reverse_ess_frac = as.numeric(pair$edge$reverse_ess_frac),
      nodes_after = as.integer(bank$n_nodes),
      particles_after = as.integer(bank$n_particles),
      check.names = FALSE
    )

    graph <- bank_smc_local_bank_overlap_graph(
      bank = bank,
      population_model = model,
      min_edge_ess_frac = min_edge_ess_frac,
      max_pareto_k = max_pareto_k,
      use_psis = use_psis
    )
  }

  insertions <- if (length(rows)) do.call(rbind, rows) else data.frame()
  bank$support_diagnostics$overlap_graph <- list(
    method = "bidirectional_importance_overlap_graph",
    certified = isTRUE(graph$certified),
    n_components = as.integer(graph$n_components),
    rounds = as.integer(length(rows)),
    inserted = as.integer(nrow(insertions)),
    min_edge_ess_frac = as.numeric(min_edge_ess_frac),
    max_pareto_k = as.numeric(max_pareto_k),
    use_psis = isTRUE(use_psis),
    edges = graph$edges
  )

  list(bank = validate_bank_smc_local_bank(bank, population_model = model), graph = graph, insertions = insertions)
}

bank_smc_overlap_graph_summary <- function(banks) {
  rows <- lapply(banks, function(bank) {
    diag <- bank$support_diagnostics$overlap_graph %||% list()
    certified <- if (!is.null(diag$certified)) isTRUE(diag$certified) else bank$n_nodes <= 1L
    n_components <- if (!is.null(diag$n_components)) as.integer(diag$n_components) else {
      if (bank$n_nodes <= 1L) 1L else NA_integer_
    }
    data.frame(
      local_id = as.integer(bank$local_id),
      certified = certified,
      n_components = n_components,
      inserted = as.integer(diag$inserted %||% NA_integer_),
      min_edge_ess_frac = as.numeric(diag$min_edge_ess_frac %||% NA_real_),
      max_pareto_k = as.numeric(diag$max_pareto_k %||% NA_real_),
      check.names = FALSE
    )
  })
  do.call(rbind, rows)
}

.bank_population_reference_prior_from_theta <- function(population_model,
                                                        theta,
                                                        label = "bank_population_reference") {
  model <- normalize_population_model(population_model)
  theta <- .bank_align_theta(theta, model)
  components <- population_model_reference_components_from_theta(model, theta)
  make_reference_prior_gaussian(
    mu = components$component_means[[1L]],
    Sigma = components$component_covs[[1L]],
    param_names = model$alpha_names,
    label = label
  )
}

.bank_weighted_mean <- function(x, w) {
  x <- as.matrix(x)
  w <- .bank_normalize_weights(w)
  as.numeric(colSums(x * w))
}

.bank_weighted_quantile <- function(x, w, probs) {
  x <- as.numeric(x)
  w <- .bank_normalize_weights(w)
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  cw <- cumsum(w)
  vapply(
    as.numeric(probs),
    function(p) x[which(cw >= p)[1L]],
    numeric(1)
  )
}

.bank_unique_theta_rows <- function(theta, model, digits = 8L) {
  theta <- .as_hyper_matrix(theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  key <- apply(round(theta, digits = digits), 1L, paste, collapse = "\r")
  theta[!duplicated(key), , drop = FALSE]
}

.bank_theta_prior_stress_candidates <- function(population_model,
                                                center,
                                                tail_probs = c(0.025, 0.975)) {
  model <- normalize_population_model(population_model)
  center <- .bank_align_theta(center, model)
  spec <- model$prior_spec %||% NULL
  rows <- list()

  if (!is.null(spec$mean_prior_mean) &&
      !is.null(spec$mean_prior_var) &&
      !is.null(spec$sigma2_prior_shape) &&
      !is.null(spec$sigma2_prior_rate)) {
    d <- model$alpha_dim
    for (j in seq_len(d)) {
      vals <- stats::qnorm(
        tail_probs,
        mean = as.numeric(spec$mean_prior_mean[j]),
        sd = sqrt(as.numeric(spec$mean_prior_var[j]))
      )
      for (val in vals) {
        row <- center
        row[1L, j] <- val
        rows[[length(rows) + 1L]] <- row
      }
    }
    for (j in seq_len(d)) {
      vals <- log(1 / stats::qgamma(
        1 - tail_probs,
        shape = as.numeric(spec$sigma2_prior_shape[j]),
        rate = as.numeric(spec$sigma2_prior_rate[j])
      ))
      for (val in vals) {
        row <- center
        row[1L, d + j] <- val
        rows[[length(rows) + 1L]] <- row
      }
    }
    for (j in seq_len(d)) {
      mean_vals <- stats::qnorm(
        tail_probs,
        mean = as.numeric(spec$mean_prior_mean[j]),
        sd = sqrt(as.numeric(spec$mean_prior_var[j]))
      )
      high_sigma2 <- log(1 / stats::qgamma(
        1 - max(tail_probs),
        shape = as.numeric(spec$sigma2_prior_shape[j]),
        rate = as.numeric(spec$sigma2_prior_rate[j])
      ))
      for (val in mean_vals) {
        row <- center
        row[1L, j] <- val
        row[1L, d + j] <- high_sigma2
        rows[[length(rows) + 1L]] <- row
      }
    }
  } else {
    draws <- population_model_sample_hyper(model, 512L)
    for (j in seq_len(model$hyper_dim)) {
      vals <- stats::quantile(draws[, j], probs = tail_probs, na.rm = TRUE, names = FALSE)
      for (val in vals) {
        row <- center
        row[1L, j] <- val
        rows[[length(rows) + 1L]] <- row
      }
    }
  }

  out <- do.call(rbind, rows)
  colnames(out) <- model$hyper_names
  .bank_unique_theta_rows(out, model)
}

.bank_score_theta_design_with_banks <- function(theta,
                                                banks,
                                                population_model,
                                                target_ess_frac = 0.3,
                                                score_local_count = 20L) {
  model <- normalize_population_model(population_model)
  theta <- .as_hyper_matrix(theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  if (is.null(banks) || !length(banks) || !nrow(theta)) {
    return(rep(0, nrow(theta)))
  }
  n_score <- min(length(banks), as.integer(max(1L, score_local_count)))
  local_ids <- unique(as.integer(round(seq(1L, length(banks), length.out = n_score))))

  vapply(seq_len(nrow(theta)), function(j) {
    coverage <- vapply(local_ids, function(i) {
      diag <- bank_smc_local_bank_coverage_diagnostic(
        bank = banks[[i]],
        population_model = model,
        theta = theta[j, , drop = FALSE],
        min_mixture_ess_frac = target_ess_frac,
        use_psis = FALSE
      )
      if (!is.null(diag$natural)) {
        as.numeric(diag$natural$best_ess_frac)
      } else {
        as.numeric(diag$mixture$ess_frac)
      }
    }, numeric(1))
    sum(pmax(0, target_ess_frac - coverage)) + (target_ess_frac - stats::median(coverage))
  }, numeric(1))
}

bank_smc_sample_alpha_configs <- function(banks,
                                          population_model,
                                          n_configs = 128L,
                                          theta = NULL,
                                          seed = NULL) {
  model <- normalize_population_model(population_model)
  n_configs <- as.integer(max(1L, n_configs))
  if (!is.null(seed)) set.seed(as.integer(seed))
  if (!is.null(theta)) theta <- .bank_align_theta(theta, model)

  out <- array(
    NA_real_,
    dim = c(n_configs, length(banks), model$alpha_dim),
    dimnames = list(NULL, names(banks) %||% as.character(seq_along(banks)), model$alpha_names)
  )

  for (i in seq_along(banks)) {
    bank <- validate_bank_smc_local_bank(banks[[i]], population_model = model)
    if (is.null(theta)) {
      node_idx <- sample.int(bank$n_nodes, n_configs, replace = TRUE, prob = bank$eta)
      alpha_i <- matrix(NA_real_, nrow = n_configs, ncol = model$alpha_dim)
      for (b in unique(node_idx)) {
        rows <- which(node_idx == b)
        node <- bank$nodes[[b]]
        pick <- sample.int(node$n_particles, length(rows), replace = TRUE, prob = node$weights)
        alpha_i[rows, ] <- node$alpha[pick, , drop = FALSE]
      }
    } else {
      terms <- bank_smc_local_bank_log_terms(
        bank = bank,
        population_model = model,
        theta = theta,
        use_cached_stack = FALSE
      )
      lw <- as.numeric(terms$log_terms[1L, ])
      w <- exp(lw - logsumexp(lw))
      pick <- sample.int(nrow(terms$stack$alpha), n_configs, replace = TRUE, prob = w)
      alpha_i <- terms$stack$alpha[pick, , drop = FALSE]
    }
    out[, i, ] <- alpha_i
  }

  out
}

.bank_diag_group_log_ell_marginal <- function(ell,
                                              sum_y,
                                              sum_y2,
                                              n,
                                              m0,
                                              s0,
                                              a,
                                              b) {
  ell <- as.numeric(ell)
  inv_sigma2 <- exp(-ell)
  inv_s0 <- 1 / as.numeric(s0)
  precision <- inv_s0 + as.numeric(n) * inv_sigma2
  h <- as.numeric(m0) * inv_s0 + as.numeric(sum_y) * inv_sigma2
  cterm <- as.numeric(m0)^2 * inv_s0 + as.numeric(sum_y2) * inv_sigma2

  log_ig_jac <- as.numeric(a) * log(as.numeric(b)) -
    lgamma(as.numeric(a)) -
    as.numeric(a) * ell -
    as.numeric(b) * inv_sigma2
  log_norm <- -0.5 * (log(2 * pi * as.numeric(s0)) + as.numeric(n) * log(2 * pi) + as.numeric(n) * ell)

  log_ig_jac + log_norm - 0.5 * cterm + 0.5 * h * h / precision + 0.5 * (log(2 * pi) - log(precision))
}

.bank_inv_gamma_quantile <- function(p, shape, rate) {
  p <- pmin(pmax(as.numeric(p), .Machine$double.eps), 1 - .Machine$double.eps)
  1 / stats::qgamma(1 - p, shape = as.numeric(shape), rate = as.numeric(rate))
}

.bank_sample_diag_group_dim <- function(y,
                                        m0,
                                        s0,
                                        a,
                                        b,
                                        n_draws,
                                        ell_grid_size = 96L,
                                        tail_prob = 1e-4) {
  y <- as.numeric(y)
  n <- length(y)
  sum_y <- sum(y)
  sum_y2 <- sum(y * y)
  empirical_var <- if (n > 1L) stats::var(y) else as.numeric(s0)
  empirical_var <- pmax(empirical_var, 1e-8)

  prior_low <- log(.bank_inv_gamma_quantile(tail_prob, a, b))
  prior_high <- log(.bank_inv_gamma_quantile(1 - tail_prob, a, b))
  lower <- min(prior_low, log(empirical_var) - 4)
  upper <- max(prior_high, log(empirical_var) + 4)

  mode <- tryCatch(
    stats::optimize(
      function(x) -.bank_diag_group_log_ell_marginal(x, sum_y, sum_y2, n, m0, s0, a, b),
      interval = c(lower, upper)
    )$minimum,
    error = function(e) NA_real_
  )
  if (is.finite(mode)) {
    lower <- min(lower, mode - 6)
    upper <- max(upper, mode + 6)
  }

  ell_grid <- seq(lower, upper, length.out = as.integer(max(32L, ell_grid_size)))
  logp <- .bank_diag_group_log_ell_marginal(ell_grid, sum_y, sum_y2, n, m0, s0, a, b)
  logp[!is.finite(logp)] <- -Inf
  if (!any(is.finite(logp))) {
    ell <- rep(log(empirical_var), n_draws)
  } else {
    w <- exp(logp - logsumexp(logp))
    ell <- sample(ell_grid, as.integer(n_draws), replace = TRUE, prob = w)
  }

  inv_sigma2 <- exp(-ell)
  inv_s0 <- 1 / as.numeric(s0)
  precision <- inv_s0 + as.numeric(n) * inv_sigma2
  h <- as.numeric(m0) * inv_s0 + as.numeric(sum_y) * inv_sigma2
  mu <- stats::rnorm(as.integer(n_draws), mean = h / precision, sd = sqrt(1 / precision))

  list(mu = mu, ell = ell)
}

bank_smc_conditional_group_theta_draws <- function(alpha_configs,
                                                   population_model,
                                                   draws_per_config = 1L,
                                                   seed = NULL,
                                                   ell_grid_size = 96L) {
  model <- normalize_population_model(population_model)
  if (!identical(model$fast_family %||% NULL, "gaussian")) {
    stop("Analytic bank shape currently requires the diagonal Gaussian population model.")
  }
  spec <- model$prior_spec
  if (is.null(spec)) stop("Population model is missing prior_spec.")

  dims <- dim(alpha_configs)
  if (length(dims) != 3L || dims[3L] != model$alpha_dim) {
    stop("alpha_configs must be n_config x n_subject x alpha_dim.")
  }
  if (!is.null(seed)) set.seed(as.integer(seed))

  n_config <- dims[1L]
  draws_per_config <- as.integer(max(1L, draws_per_config))
  n_out <- n_config * draws_per_config
  theta <- matrix(NA_real_, nrow = n_out, ncol = model$hyper_dim)
  colnames(theta) <- model$hyper_names

  row_start <- 1L
  for (k in seq_len(n_config)) {
    rows <- seq.int(row_start, row_start + draws_per_config - 1L)
    alpha_k <- alpha_configs[k, , , drop = FALSE]
    alpha_k <- matrix(alpha_k, ncol = model$alpha_dim)
    colnames(alpha_k) <- model$alpha_names

    for (j in seq_len(model$alpha_dim)) {
      draw <- .bank_sample_diag_group_dim(
        y = alpha_k[, j],
        m0 = spec$mean_prior_mean[j],
        s0 = spec$mean_prior_var[j],
        a = spec$sigma2_prior_shape[j],
        b = spec$sigma2_prior_rate[j],
        n_draws = draws_per_config,
        ell_grid_size = ell_grid_size
      )
      theta[rows, j] <- draw$mu
      theta[rows, model$alpha_dim + j] <- draw$ell
    }
    row_start <- row_start + draws_per_config
  }

  theta
}

bank_smc_build_analytic_shape <- function(banks,
                                          population_model,
                                          theta_reference = NULL,
                                          n_alpha_configs = 128L,
                                          draws_per_config = 1L,
                                          ell_grid_size = 96L,
                                          seed = 123L) {
  model <- normalize_population_model(population_model)
  alpha_configs <- bank_smc_sample_alpha_configs(
    banks = banks,
    population_model = model,
    n_configs = n_alpha_configs,
    theta = theta_reference,
    seed = seed
  )
  theta <- bank_smc_conditional_group_theta_draws(
    alpha_configs = alpha_configs,
    population_model = model,
    draws_per_config = draws_per_config,
    seed = seed + 1L,
    ell_grid_size = ell_grid_size
  )
  structure(
    list(
      theta = theta,
      w = rep(1 / nrow(theta), nrow(theta)),
      population_model = model,
      meta = list(
        source = "bank_imputed_conditional_group_shape",
        n_alpha_configs = as.integer(n_alpha_configs),
        draws_per_config = as.integer(draws_per_config),
        ell_grid_size = as.integer(ell_grid_size)
      )
    ),
    class = "bank_smc_analytic_shape"
  )
}

.bank_population_broad_reference_prior_from_theta <- function(population_model,
                                                              theta,
                                                              scale = 16,
                                                              defensive = TRUE,
                                                              defensive_scale = 64,
                                                              defensive_weight = 0.20,
                                                              label = "bank_alpha_first_reference") {
  model <- normalize_population_model(population_model)
  theta <- .bank_align_theta(theta, model)
  components <- population_model_reference_components_from_theta(model, theta)
  make_broad_reference_prior(
    mu = components$component_means[[1L]],
    Sigma = components$component_covs[[1L]],
    scale = scale,
    defensive = defensive,
    defensive_scale = defensive_scale,
    defensive_weight = defensive_weight,
    param_names = model$alpha_names,
    label = label
  )
}

bank_smc_alpha_first_local_draws <- function(data_list,
                                             loglik_fn,
                                             population_model,
                                             theta_reference,
                                             M = 300L,
                                             scale = 16,
                                             defensive = TRUE,
                                             defensive_scale = 64,
                                             defensive_weight = 0.20,
                                             resample_threshold = 0.6,
                                             n_mcmc_moves = 2L,
                                             max_rounds = 80L,
                                             cess_target = 0.9,
                                             G_mix = 8L,
                                             rw_scale = 0.9,
                                             n_jobs = 1L,
                                             local_n_cores = 1L,
                                             seed = 123L,
                                             verbose = FALSE) {
  model <- normalize_population_model(population_model)
  theta_reference <- .bank_align_theta(theta_reference, model)
  reference_prior <- .bank_population_broad_reference_prior_from_theta(
    population_model = model,
    theta = theta_reference,
    scale = scale,
    defensive = defensive,
    defensive_scale = defensive_scale,
    defensive_weight = defensive_weight,
    label = sprintf("alpha_first_scale_%s", format(scale, trim = TRUE))
  )
  ids <- seq_along(data_list)
  fits <- parallel::mclapply(
    ids,
    function(i) {
      fit <- run_tempered_smc(
        reference_prior = reference_prior,
        bridge_stat_fn = function(alpha) ll_parallel(
          alpha,
          data_list[[i]],
          loglik_fn,
          n_cores = as.integer(local_n_cores)
        ),
        M = as.integer(M),
        resample_threshold = as.numeric(resample_threshold),
        n_mcmc_moves = as.integer(n_mcmc_moves),
        post_adapt_n_mcmc_moves = as.integer(n_mcmc_moves),
        max_rounds = as.integer(max_rounds),
        cess_target = as.numeric(cess_target),
        G_mix = as.integer(G_mix),
        rw_scale_init = as.numeric(rw_scale),
        n_cores = as.integer(local_n_cores),
        seed = as.integer(seed + i - 1L),
        verbose = verbose
      )
      structure(
        list(
          local_id = as.integer(i),
          alpha = .bank_align_alpha(fit$Theta, model),
          weights = .bank_normalize_weights(fit$w),
          diagnostics = list(
            final_ess_frac = as.numeric(ESS(fit$w) / length(fit$w)),
            rounds = as.integer(fit$meta$rounds %||% NA_integer_),
            mean_accept_rate = .bank_mean_or_na(fit$meta$accept_rate_hist)
          )
        ),
        class = "bank_smc_alpha_first_local_fit"
      )
    },
    mc.cores = as.integer(max(1L, n_jobs))
  )
  names(fits) <- names(data_list) %||% as.character(ids)

  structure(
    list(
      fits = fits,
      population_model = model,
      meta = list(
        source = "alpha_first_local_draws",
        scale = as.numeric(scale),
        M = as.integer(M),
        defensive = isTRUE(defensive),
        defensive_scale = as.numeric(defensive_scale),
        defensive_weight = as.numeric(defensive_weight),
        median_ess_frac = stats::median(vapply(
          fits,
          function(fit) fit$diagnostics$final_ess_frac,
          numeric(1)
        ), na.rm = TRUE)
      )
    ),
    class = "bank_smc_alpha_first_local_draws"
  )
}

bank_smc_sample_alpha_configs_from_alpha_first <- function(draws,
                                                          population_model = NULL,
                                                          n_configs = 128L,
                                                          seed = NULL) {
  if (!inherits(draws, "bank_smc_alpha_first_local_draws")) {
    stop("draws must inherit from 'bank_smc_alpha_first_local_draws'.")
  }
  model <- normalize_population_model(population_model %||% draws$population_model)
  n_configs <- as.integer(max(1L, n_configs))
  if (!is.null(seed)) set.seed(as.integer(seed))

  fits <- draws$fits
  out <- array(
    NA_real_,
    dim = c(n_configs, length(fits), model$alpha_dim),
    dimnames = list(NULL, names(fits) %||% as.character(seq_along(fits)), model$alpha_names)
  )
  for (i in seq_along(fits)) {
    alpha_i <- .bank_align_alpha(fits[[i]]$alpha, model)
    w_i <- .bank_normalize_weights(fits[[i]]$weights)
    pick <- sample.int(nrow(alpha_i), n_configs, replace = TRUE, prob = w_i)
    out[, i, ] <- alpha_i[pick, , drop = FALSE]
  }
  out
}

bank_smc_rank_coherent_alpha_configs_from_alpha_first <- function(draws,
                                                                 population_model = NULL,
                                                                 probs = c(0.025, 0.05, 0.1, 0.2, 0.8, 0.9, 0.95, 0.975)) {
  if (!inherits(draws, "bank_smc_alpha_first_local_draws")) {
    stop("draws must inherit from 'bank_smc_alpha_first_local_draws'.")
  }
  model <- normalize_population_model(population_model %||% draws$population_model)
  probs <- sort(unique(pmin(pmax(as.numeric(probs), 1e-4), 1 - 1e-4)))
  fits <- draws$fits

  quantiles <- array(
    NA_real_,
    dim = c(length(probs) + 1L, length(fits), model$alpha_dim),
    dimnames = list(c("median", paste0("q", probs)), names(fits) %||% as.character(seq_along(fits)), model$alpha_names)
  )
  for (i in seq_along(fits)) {
    alpha_i <- .bank_align_alpha(fits[[i]]$alpha, model)
    w_i <- .bank_normalize_weights(fits[[i]]$weights)
    for (j in seq_len(model$alpha_dim)) {
      quantiles[1L, i, j] <- .bank_weighted_quantile(alpha_i[, j], w_i, 0.5)
      quantiles[-1L, i, j] <- .bank_weighted_quantile(alpha_i[, j], w_i, probs)
    }
  }

  n_configs <- 1L + model$alpha_dim * length(probs)
  out <- array(
    NA_real_,
    dim = c(n_configs, length(fits), model$alpha_dim),
    dimnames = list(NULL, names(fits) %||% as.character(seq_along(fits)), model$alpha_names)
  )
  out[1L, , ] <- quantiles[1L, , ]
  ptr <- 2L
  for (j in seq_len(model$alpha_dim)) {
    for (p_idx in seq_along(probs)) {
      out[ptr, , ] <- quantiles[1L, , ]
      out[ptr, , j] <- quantiles[p_idx + 1L, , j]
      ptr <- ptr + 1L
    }
  }
  out
}

.bank_smc_bind_alpha_config_arrays <- function(parts, population_model) {
  parts <- parts[vapply(parts, function(x) length(dim(x)) == 3L && dim(x)[1L] > 0L, logical(1))]
  if (!length(parts)) stop("No alpha configuration arrays to bind.")
  model <- normalize_population_model(population_model)
  n_total <- sum(vapply(parts, function(x) dim(x)[1L], integer(1)))
  n_local <- dim(parts[[1L]])[2L]
  out <- array(
    NA_real_,
    dim = c(n_total, n_local, model$alpha_dim),
    dimnames = list(NULL, dimnames(parts[[1L]])[[2L]], model$alpha_names)
  )
  ptr <- 1L
  for (part in parts) {
    n <- dim(part)[1L]
    out[seq.int(ptr, ptr + n - 1L), , ] <- part
    ptr <- ptr + n
  }
  out
}

bank_smc_alpha_first_anchor_candidates <- function(data_list,
                                                  loglik_fn,
                                                  population_model,
                                                  theta_reference,
                                                  n_alpha_configs = 128L,
                                                  draws_per_config = 1L,
                                                  ell_grid_size = 96L,
                                                  M = 300L,
                                                  scales = c(16, 4),
                                                  coherent_probs = c(0.025, 0.05, 0.1, 0.2, 0.8, 0.9, 0.95, 0.975),
                                                  defensive = TRUE,
                                                  defensive_scale = 64,
                                                  defensive_weight = 0.20,
                                                  resample_threshold = 0.6,
                                                  n_mcmc_moves = 2L,
                                                  max_rounds = 80L,
                                                  cess_target = 0.9,
                                                  G_mix = 8L,
                                                  rw_scale = 0.9,
                                                  n_jobs = 1L,
                                                  local_n_cores = 1L,
                                                  seed = 123L,
                                                  verbose = FALSE) {
  model <- normalize_population_model(population_model)
  scales <- as.numeric(scales)
  scales <- scales[is.finite(scales) & scales > 0]
  if (!length(scales)) stop("Alpha-first scales must contain at least one positive value.")

  n_alpha_configs <- as.integer(max(1L, n_alpha_configs))
  draws_per_config <- as.integer(max(1L, draws_per_config))
  config_counts <- rep(floor(n_alpha_configs / length(scales)), length(scales))
  config_counts[seq_len(n_alpha_configs %% length(scales))] <-
    config_counts[seq_len(n_alpha_configs %% length(scales))] + 1L

  theta_parts <- vector("list", length(scales))
  meta_rows <- vector("list", length(scales))
  for (k in seq_along(scales)) {
    local_draws <- bank_smc_alpha_first_local_draws(
      data_list = data_list,
      loglik_fn = loglik_fn,
      population_model = model,
      theta_reference = theta_reference,
      M = M,
      scale = scales[k],
      defensive = defensive,
      defensive_scale = defensive_scale,
      defensive_weight = defensive_weight,
      resample_threshold = resample_threshold,
      n_mcmc_moves = n_mcmc_moves,
      max_rounds = max_rounds,
      cess_target = cess_target,
      G_mix = G_mix,
      rw_scale = rw_scale,
      n_jobs = n_jobs,
      local_n_cores = local_n_cores,
      seed = as.integer(seed + 100000L * k),
      verbose = verbose
    )
    coherent_configs <- bank_smc_rank_coherent_alpha_configs_from_alpha_first(
      draws = local_draws,
      population_model = model,
      probs = coherent_probs
    )
    if (dim(coherent_configs)[1L] > config_counts[k]) {
      keep <- unique(round(seq(1L, dim(coherent_configs)[1L], length.out = config_counts[k])))
      coherent_configs <- coherent_configs[keep, , , drop = FALSE]
    }
    n_random <- as.integer(config_counts[k] - dim(coherent_configs)[1L])
    alpha_configs <- if (n_random > 0L) {
      random_configs <- bank_smc_sample_alpha_configs_from_alpha_first(
        draws = local_draws,
        population_model = model,
        n_configs = n_random,
        seed = as.integer(seed + 200000L * k)
      )
      .bank_smc_bind_alpha_config_arrays(list(coherent_configs, random_configs), model)
    } else {
      coherent_configs
    }
    theta_parts[[k]] <- bank_smc_conditional_group_theta_draws(
      alpha_configs = alpha_configs,
      population_model = model,
      draws_per_config = draws_per_config,
      seed = as.integer(seed + 300000L * k),
      ell_grid_size = ell_grid_size
    )
    meta_rows[[k]] <- data.frame(
      scale = as.numeric(scales[k]),
      n_alpha_configs = as.integer(config_counts[k]),
      coherent_configs = as.integer(dim(coherent_configs)[1L]),
      median_local_ess_frac = as.numeric(local_draws$meta$median_ess_frac),
      check.names = FALSE
    )
  }

  theta <- do.call(rbind, theta_parts)
  structure(
    list(
      theta = theta,
      w = rep(1 / nrow(theta), nrow(theta)),
      population_model = model,
      meta = list(
        source = "alpha_first_anchor_candidates",
        scales = scales,
        n_alpha_configs = as.integer(n_alpha_configs),
        draws_per_config = as.integer(draws_per_config),
        ell_grid_size = as.integer(ell_grid_size),
        local = do.call(rbind, meta_rows)
      )
    ),
    class = "bank_smc_anchor_candidates"
  )
}

bank_smc_anchor_design_fit <- function(anchor_candidates,
                                       reference_fit,
                                       population_model = NULL,
                                       anchor_weight = 0.70) {
  if (is.null(anchor_candidates)) return(reference_fit)
  model <- normalize_population_model(population_model %||% anchor_candidates$population_model)
  anchor_theta <- .as_hyper_matrix(anchor_candidates$theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  anchor_w <- .bank_normalize_weights(anchor_candidates$w %||% rep(1, nrow(anchor_theta)))
  ref_theta <- .as_hyper_matrix(reference_fit$theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  ref_w <- .bank_normalize_weights(reference_fit$w)

  anchor_weight <- min(max(as.numeric(anchor_weight), 0), 1)
  theta <- rbind(anchor_theta, ref_theta)
  colnames(theta) <- model$hyper_names
  w <- .bank_normalize_weights(c(anchor_weight * anchor_w, (1 - anchor_weight) * ref_w))

  structure(
    list(
      theta = theta,
      w = w,
      population_model = model,
      meta = list(
        source = "alpha_first_anchor_design",
        anchor_weight = anchor_weight,
        anchor_rows = as.integer(nrow(anchor_theta)),
        reference_rows = as.integer(nrow(ref_theta)),
        anchor_source = anchor_candidates$meta$source %||% NA_character_,
        reference_source = reference_fit$meta$source %||% NA_character_
      )
    ),
    class = "bank_smc_analytic_shape"
  )
}

bank_smc_select_theta_design <- function(population_fit,
                                         population_model,
                                         max_points = 6L,
                                         support_points = 3L,
                                         profile_dims = 1L,
                                         tail_probs = c(0.1, 0.9),
                                         include_center = TRUE,
                                         banks = NULL,
                                         initial_theta = NULL,
                                         extra_theta = NULL,
                                         prior_stress = FALSE,
                                         prior_tail_probs = c(0.025, 0.975),
                                         paired_effect_profiles = TRUE,
                                         protect_profile_tails = TRUE,
                                         target_ess_frac = 0.3,
                                         score_local_count = 20L) {
  model <- normalize_population_model(population_model)
  theta <- .as_hyper_matrix(population_fit$theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  w <- .bank_normalize_weights(population_fit$w)
  center <- .bank_weighted_mean(theta, w)
  names(center) <- model$hyper_names

  max_points <- as.integer(max(1L, max_points))
  profile_dims <- as.integer(max(0L, profile_dims))
  S <- if (nrow(theta) > 1L) {
    tryCatch(weighted_cov(theta, w), error = function(e) stats::cov(theta))
  } else {
    matrix(0, nrow = model$hyper_dim, ncol = model$hyper_dim)
  }
  S <- regularize_cov(S, min_eig = 1e-8, cond_cap = 1e8)
  spread <- sqrt(pmax(diag(S), 0))
  spread[!is.finite(spread)] <- 0

  protected <- list()
  if (isTRUE(include_center)) protected[[length(protected) + 1L]] <- center

  rows <- list()

  support_points <- as.integer(max(0L, support_points))
  if (support_points > 0L && nrow(theta) > 1L) {
    eig <- eigen(S, symmetric = TRUE)
    score <- as.numeric(scale(theta, center = center, scale = FALSE) %*% eig$vectors[, 1L])
    q <- .bank_weighted_quantile(score, w, seq(0.15, 0.85, length.out = support_points))
    for (val in q) {
      idx <- which.min(abs(score - val))
      rows[[length(rows) + 1L]] <- theta[idx, ]
    }
  }

  if (profile_dims > 0L) {
    dims <- head(order(spread, decreasing = TRUE), profile_dims)
    for (j in dims) {
      vals <- .bank_weighted_quantile(theta[, j], w, tail_probs)
      if (isTRUE(protect_profile_tails) && length(vals)) {
        for (val in range(vals, na.rm = TRUE)) {
          row <- center
          row[j] <- val
          protected[[length(protected) + 1L]] <- row
        }
      }
      for (val in vals) {
        row <- center
        row[j] <- val
        rows[[length(rows) + 1L]] <- row
      }
    }
  }

  if (isTRUE(paired_effect_profiles) &&
      identical(model$fast_family %||% NULL, "gaussian") &&
      model$hyper_dim == 2L * model$alpha_dim &&
      profile_dims > 0L) {
    d <- model$alpha_dim
    effect_score <- spread[seq_len(d)] + spread[d + seq_len(d)]
    effect_dims <- head(order(effect_score, decreasing = TRUE), min(profile_dims, d))
    pair_probs <- sort(unique(as.numeric(tail_probs)))
    for (j in effect_dims) {
      mean_vals <- .bank_weighted_quantile(theta[, j], w, pair_probs)
      var_vals <- .bank_weighted_quantile(theta[, d + j], w, pair_probs)
      for (mean_val in mean_vals) {
        for (var_val in var_vals) {
          row <- center
          row[j] <- mean_val
          row[d + j] <- var_val
          rows[[length(rows) + 1L]] <- row
        }
      }
    }
  }

  if (isTRUE(prior_stress)) {
    stress_center <- initial_theta %||% matrix(center, nrow = 1L)
    stress <- .bank_theta_prior_stress_candidates(model, stress_center, tail_probs = prior_tail_probs)
    rows <- c(rows, lapply(seq_len(nrow(stress)), function(i) stress[i, , drop = FALSE]))
  }
  if (!is.null(extra_theta)) {
    extra_theta <- .as_hyper_matrix(extra_theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
    rows <- c(rows, lapply(seq_len(nrow(extra_theta)), function(i) extra_theta[i, , drop = FALSE]))
  }

  protected_mat <- if (length(protected)) {
    .bank_unique_theta_rows(do.call(rbind, protected), model)
  } else {
    matrix(numeric(0), nrow = 0L, ncol = model$hyper_dim, dimnames = list(NULL, model$hyper_names))
  }
  n_protected <- nrow(protected_mat)
  out <- if (length(rows)) {
    rbind(protected_mat, do.call(rbind, rows))
  } else {
    protected_mat
  }
  colnames(out) <- model$hyper_names
  out <- .bank_unique_theta_rows(out, model)
  if (!is.null(banks)) {
    keep_first <- seq_len(min(n_protected, nrow(out)))
    score <- .bank_score_theta_design_with_banks(
      theta = out,
      banks = banks,
      population_model = model,
      target_ess_frac = target_ess_frac,
      score_local_count = score_local_count
    )
    candidates <- setdiff(seq_len(nrow(out)), keep_first)
    ranked <- candidates[order(score[candidates], decreasing = TRUE)]
    take <- unique(c(keep_first, ranked))
    out <- out[take, , drop = FALSE]
  }
  out[seq_len(min(nrow(out), max_points)), , drop = FALSE]
}

bank_smc_select_joint_tail_theta <- function(population_fit,
                                             population_model,
                                             max_points = 4L,
                                             profile_dims = 1L,
                                             tail_probs = c(0.01, 0.05, 0.95, 0.99)) {
  model <- normalize_population_model(population_model)
  theta <- .as_hyper_matrix(population_fit$theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  w <- .bank_normalize_weights(population_fit$w)
  if (nrow(theta) < 1L || max_points <= 0L) {
    return(matrix(numeric(0), nrow = 0L, ncol = model$hyper_dim, dimnames = list(NULL, model$hyper_names)))
  }

  profile_dims <- as.integer(max(1L, profile_dims))
  tail_probs <- sort(unique(pmin(pmax(as.numeric(tail_probs), 1e-4), 1 - 1e-4)))
  S <- if (nrow(theta) > 1L) {
    tryCatch(weighted_cov(theta, w), error = function(e) stats::cov(theta))
  } else {
    matrix(0, nrow = model$hyper_dim, ncol = model$hyper_dim)
  }
  spread <- sqrt(pmax(diag(regularize_cov(S, min_eig = 1e-8, cond_cap = 1e8)), 0))
  spread[!is.finite(spread)] <- 0
  dims <- head(order(spread, decreasing = TRUE), min(profile_dims, model$hyper_dim))

  idx <- integer(0)
  for (j in dims) {
    vals <- .bank_weighted_quantile(theta[, j], w, tail_probs)
    idx <- c(idx, vapply(vals, function(v) which.min(abs(theta[, j] - v)), integer(1)))
  }
  idx <- unique(idx)
  if (length(idx) > as.integer(max_points)) idx <- idx[seq_len(as.integer(max_points))]
  out <- theta[idx, , drop = FALSE]
  colnames(out) <- model$hyper_names
  .bank_unique_theta_rows(out, model)
}

bank_smc_select_inflated_profile_theta <- function(population_fit,
                                                   population_model,
                                                   max_points = 8L,
                                                   profile_dims = 4L,
                                                   tail_probs = c(0.05, 0.95),
                                                   paired_mean_tail_probs = NULL,
                                                   paired_variance_tail_probs = NULL,
                                                   inflation_scale = 2.5,
                                                   paired_effect_profiles = TRUE,
                                                   paired_effect_grid = TRUE) {
  model <- normalize_population_model(population_model)
  theta <- .as_hyper_matrix(population_fit$theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  w <- .bank_normalize_weights(population_fit$w)
  empty <- matrix(numeric(0), nrow = 0L, ncol = model$hyper_dim, dimnames = list(NULL, model$hyper_names))
  max_points <- as.integer(max_points)
  if (nrow(theta) < 2L || max_points <= 0L) return(empty)

  center <- .bank_weighted_mean(theta, w)
  names(center) <- model$hyper_names
  S <- tryCatch(weighted_cov(theta, w), error = function(e) stats::cov(theta))
  spread <- sqrt(pmax(diag(regularize_cov(S, min_eig = 1e-8, cond_cap = 1e8)), 0))
  spread[!is.finite(spread)] <- 0
  tail_probs <- sort(unique(pmin(pmax(as.numeric(tail_probs), 1e-4), 1 - 1e-4)))
  z <- stats::qnorm(tail_probs)
  mean_tail_probs <- sort(unique(pmin(pmax(as.numeric(paired_mean_tail_probs %||% tail_probs), 1e-4), 1 - 1e-4)))
  variance_tail_probs <- sort(unique(pmin(pmax(as.numeric(paired_variance_tail_probs %||% tail_probs), 1e-4), 1 - 1e-4)))
  z_mean_grid <- stats::qnorm(mean_tail_probs)
  z_variance_grid <- stats::qnorm(variance_tail_probs)
  profile_dims <- as.integer(max(1L, profile_dims))
  rows <- list()

  add_axis_rows <- function(j) {
    if (!is.finite(spread[j]) || spread[j] <= 0) return(NULL)
    for (zz in z) {
      row <- center
      row[j] <- center[j] + as.numeric(inflation_scale) * zz * spread[j]
      rows[[length(rows) + 1L]] <<- row
    }
    NULL
  }

  if (isTRUE(paired_effect_profiles) &&
      identical(model$fast_family %||% NULL, "gaussian") &&
      model$hyper_dim == 2L * model$alpha_dim) {
    d <- model$alpha_dim
    effect_score <- spread[seq_len(d)] + spread[d + seq_len(d)]
    effects <- head(order(effect_score, decreasing = TRUE), min(profile_dims, d))
    for (j in effects) {
      if (isTRUE(paired_effect_grid) &&
          is.finite(spread[j]) && spread[j] > 0 &&
          is.finite(spread[d + j]) && spread[d + j] > 0) {
        for (z_mean in z_mean_grid) {
          for (z_var in z_variance_grid) {
            row <- center
            row[j] <- center[j] + as.numeric(inflation_scale) * z_mean * spread[j]
            row[d + j] <- center[d + j] + as.numeric(inflation_scale) * z_var * spread[d + j]
            rows[[length(rows) + 1L]] <- row
          }
        }
      } else {
        add_axis_rows(j)
        add_axis_rows(d + j)
      }
    }
  } else {
    dims <- head(order(spread, decreasing = TRUE), min(profile_dims, model$hyper_dim))
    for (j in dims) add_axis_rows(j)
  }

  if (!length(rows)) return(empty)
  out <- do.call(rbind, rows)
  colnames(out) <- model$hyper_names
  out <- .bank_unique_theta_rows(out, model)
  log_prior <- population_model_log_hyperprior(model, out)
  out <- out[is.finite(log_prior), , drop = FALSE]
  if (!nrow(out)) return(empty)
  out[seq_len(min(nrow(out), max_points)), , drop = FALSE]
}

bank_smc_select_surface_theta <- function(population_fit,
                                          population_model,
                                          max_points = 8L,
                                          profile_dims = 4L,
                                          tail_probs = c(0.05, 0.95),
                                          paired_mean_tail_probs = NULL,
                                          paired_variance_tail_probs = NULL,
                                          inflation_scale = 2.5,
                                          include_inflated_profile = TRUE,
                                          include_posterior_tail = TRUE,
                                          paired_effect_profiles = TRUE,
                                          paired_effect_grid = TRUE,
                                          candidate_multiplier = 4L,
                                          banks = NULL,
                                          score_local_count = 20L,
                                          score_log_weight = 1,
                                          score_score_weight = 0.25,
                                          score_distance_weight = 0.05,
                                          n_jobs = 1L,
                                          return_diagnostics = FALSE) {
  model <- normalize_population_model(population_model)
  max_points <- as.integer(max_points)
  empty <- matrix(numeric(0), nrow = 0L, ncol = model$hyper_dim, dimnames = list(NULL, model$hyper_names))
  if (max_points <= 0L) return(empty)
  pool_points <- max_points * as.integer(max(1L, candidate_multiplier))

  rows <- list()
  if (isTRUE(include_inflated_profile)) {
    inflated <- bank_smc_select_inflated_profile_theta(
      population_fit = population_fit,
      population_model = model,
      max_points = pool_points,
      profile_dims = profile_dims,
      tail_probs = tail_probs,
      paired_mean_tail_probs = paired_mean_tail_probs,
      paired_variance_tail_probs = paired_variance_tail_probs,
      inflation_scale = inflation_scale,
      paired_effect_profiles = paired_effect_profiles,
      paired_effect_grid = paired_effect_grid
    )
    if (nrow(inflated)) rows[[length(rows) + 1L]] <- inflated
  }
  if (isTRUE(include_posterior_tail)) {
    posterior_tail <- bank_smc_select_joint_tail_theta(
      population_fit = population_fit,
      population_model = model,
      max_points = pool_points,
      profile_dims = profile_dims,
      tail_probs = tail_probs
    )
    if (nrow(posterior_tail)) rows[[length(rows) + 1L]] <- posterior_tail
  }
  if (!length(rows)) return(empty)
  out <- .bank_unique_theta_rows(do.call(rbind, rows), model)
  score_info <- NULL
  if (!is.null(banks) && length(banks) && nrow(out) > max_points) {
    score_info <- bank_smc_score_surface_theta_design_with_banks(
      theta = out,
      banks = banks,
      population_model = model,
      local_count = score_local_count,
      log_weight = score_log_weight,
      score_weight = score_score_weight,
      distance_weight = score_distance_weight,
      n_jobs = n_jobs
    )
    ranked <- order(score_info$score, decreasing = TRUE)
    out <- out[ranked, , drop = FALSE]
    if (!is.null(score_info$diagnostics) && nrow(score_info$diagnostics)) {
      score_info$ranking <- data.frame(
        theta_id = seq_along(score_info$score),
        surface_score = score_info$score,
        rank = rank(-score_info$score, ties.method = "first"),
        check.names = FALSE
      )
    }
  }
  selected <- out[seq_len(min(nrow(out), max_points)), , drop = FALSE]
  if (isTRUE(return_diagnostics)) {
    return(list(theta = selected, score = score_info))
  }
  selected
}

bank_smc_initial_banks <- function(data_list,
                                   loglik_fn,
                                   population_model,
                                   initial_theta,
                                   M = 1000L,
                                   max_nodes = 6L,
                                   max_particles = Inf,
                                   n_jobs = 1L,
                                   local_n_cores = 1L,
                                   seed = 123L,
                                   verbose = FALSE,
                                   ...) {
  model <- normalize_population_model(population_model)
  initial_theta <- .bank_align_theta(initial_theta, model)
  reference_prior <- .bank_population_reference_prior_from_theta(
    model,
    initial_theta,
    label = "bank_initial_anchor"
  )
  extra <- list(...)
  ids <- seq_along(data_list)
  banks <- parallel::mclapply(
    ids,
    function(i) {
      fit <- do.call(
        run_tempered_smc,
        modifyList(
          list(
            reference_prior = reference_prior,
            bridge_stat_fn = function(alpha) ll_parallel(
              alpha,
              data_list[[i]],
              loglik_fn,
              n_cores = as.integer(local_n_cores)
            ),
            M = as.integer(M),
            n_cores = as.integer(local_n_cores),
            seed = as.integer(seed + i - 1L),
            verbose = verbose
          ),
          extra
        )
      )
      node <- build_bank_smc_node_from_local_fit(
        fit,
        theta_anchor = initial_theta,
        population_model = model,
        local_id = i,
        bridge_provenance = list(type = "initial")
      )
      new_bank_smc_local_bank(
        local_id = i,
        nodes = list(node),
        max_nodes = max_nodes,
        max_particles = max_particles
      )
    },
    mc.cores = as.integer(max(1L, n_jobs))
  )
  names(banks) <- names(data_list) %||% as.character(ids)
  banks
}

bank_smc_local_bank_factor_object <- function(bank, population_model) {
  model <- normalize_population_model(population_model)
  bank <- validate_bank_smc_local_bank(bank, population_model = model)
  stack <- bank_smc_local_bank_build_stack(bank, model)
  structure(
    list(
      local_id = bank$local_id,
      factor_sufficient_stats = population_sufficient_stats_from_alpha(model, stack$alpha),
      factor_log_base = stack$log_base,
      factor_log_constant = 0,
      diagnostics = list(
        mode = "bank_smc",
        representation = "sufficient_stats",
        n_nodes = bank$n_nodes,
        n_particles = bank$n_particles
      )
    ),
    class = "compressed_population_local_object"
  )
}

bank_smc_banks_to_factor_set <- function(banks,
                                         population_model,
                                         particle_block_size = 1024L) {
  model <- normalize_population_model(population_model)
  local_objects <- bank_smc_banks_to_local_objects(banks, model)
  build_population_factor_set(
    local_objects = local_objects,
    population_model = model,
    particle_block_size = particle_block_size
  )
}

bank_smc_banks_to_local_objects <- function(banks, population_model) {
  model <- normalize_population_model(population_model)
  out <- lapply(banks, bank_smc_local_bank_factor_object, population_model = model)
  names(out) <- names(banks)
  out
}

bank_smc_bank_budget_summary <- function(banks) {
  data.frame(
    local_id = vapply(banks, function(bank) as.integer(bank$local_id), integer(1)),
    nodes = vapply(banks, function(bank) as.integer(bank$n_nodes), integer(1)),
    particles = vapply(banks, function(bank) as.integer(bank$n_particles), integer(1)),
    max_nodes = vapply(banks, function(bank) as.numeric(bank$max_nodes), numeric(1)),
    max_particles = vapply(banks, function(bank) as.numeric(bank$max_particles), numeric(1)),
    check.names = FALSE
  )
}

bank_smc_refine_banks_to_design <- function(banks,
                                            theta_design,
                                            data_list,
                                            loglik_fn,
                                            population_model,
                                            target_ess_frac = 0.3,
                                            max_pareto_k = 0.7,
                                            bridge_particles = NULL,
                                            target_cess = 0.9,
                                            n_mcmc_moves = 2L,
                                            max_steps = 128L,
                                            bridge_min_single_ess_frac = 0.05,
                                            bridge_max_pareto_k = 0.7,
                                            force = FALSE,
                                            anchor_tol = 1e-8,
                                            graph_control = list(),
                                            n_jobs = 1L,
                                            local_n_cores = 1L,
                                            seed = 123L,
                                            verbose = FALSE) {
  model <- normalize_population_model(population_model)
  theta_design <- .as_hyper_matrix(theta_design, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  graph_defaults <- list(
    enabled = TRUE,
    min_edge_ess_frac = 0.05,
    max_pareto_k = max_pareto_k,
    use_psis = FALSE,
    max_rounds = 50L
  )
  graph_control <- modifyList(graph_defaults, graph_control)
  ids <- seq_along(banks)
  rows <- vector("list", 0L)
  graph_rows <- vector("list", 0L)
  refined <- parallel::mclapply(
    ids,
    function(pos) {
      bank <- validate_bank_smc_local_bank(banks[[pos]], population_model = model)
      local_rows <- vector("list", 0L)
      for (j in seq_len(nrow(theta_design))) {
        theta_j <- theta_design[j, , drop = FALSE]
        diag <- bank_smc_local_bank_coverage_diagnostic(
          bank,
          model,
          theta_j,
          min_mixture_ess_frac = target_ess_frac,
          max_pareto_k = max_pareto_k,
          use_psis = FALSE
        )
        anchors <- bank_smc_local_bank_theta_anchors(bank, model)
        has_anchor <- any(rowSums(abs(sweep(anchors, 2L, as.numeric(theta_j), "-"))) <= as.numeric(anchor_tol))
        added <- FALSE
        log_marginal_before <- as.numeric(diag$mixture$log_marginal)
        log_marginal_anchor <- NA_real_
        predicted_anchor_error <- NA_real_
        forced <- isTRUE(force) && !has_anchor
        if ((!isTRUE(diag$covered) || forced) && bank$n_nodes < bank$max_nodes) {
          if (is.finite(bank$max_particles) &&
              bank$n_particles + as.integer(bridge_particles %||% bank$nodes[[diag$best_node]]$n_particles) > bank$max_particles) {
            added <- FALSE
          } else {
            bridged <- bank_smc_local_bank_add_anchor(
              bank = bank,
              population_model = model,
              theta_target = theta_j,
              data_i = data_list[[pos]],
              loglik_fn = loglik_fn,
              M = bridge_particles,
              target_cess = target_cess,
              n_mcmc_moves = n_mcmc_moves,
              max_steps = max_steps,
              bridge_min_single_ess_frac = bridge_min_single_ess_frac,
              bridge_max_pareto_k = bridge_max_pareto_k,
              n_cores = local_n_cores,
              seed = as.integer(seed + 10000L * pos + j),
              verbose = verbose,
              use_psis = FALSE
            )
            bank <- bridged$bank
            added <- TRUE
            log_marginal_anchor <- as.numeric(bridged$node$log_marginal_anchor)
            predicted_anchor_error <- log_marginal_anchor - log_marginal_before
          }
        }
        local_rows[[length(local_rows) + 1L]] <- data.frame(
          local_id = as.integer(bank$local_id),
          theta_id = as.integer(j),
          covered_before = isTRUE(diag$covered),
          forced = isTRUE(forced),
          added = isTRUE(added),
          ess_frac = as.numeric(diag$mixture$ess_frac),
          natural_ess_frac = as.numeric(diag$natural$best_ess_frac %||% NA_real_),
          natural_variance = as.numeric(diag$natural$best_variance %||% NA_real_),
          best_node = as.integer(diag$best_node),
          log_marginal_before = log_marginal_before,
          log_marginal_anchor = log_marginal_anchor,
          predicted_anchor_error = predicted_anchor_error,
          nodes_after = as.integer(bank$n_nodes),
          particles_after = as.integer(bank$n_particles),
          check.names = FALSE
        )
      }
      graph_insertions <- data.frame()
      if (isTRUE(graph_control$enabled) && bank$n_nodes > 1L) {
        certified <- bank_smc_local_bank_certify_bridge_graph(
          bank = bank,
          population_model = model,
          data_i = data_list[[pos]],
          loglik_fn = loglik_fn,
          M = bridge_particles,
          min_edge_ess_frac = graph_control$min_edge_ess_frac,
          max_pareto_k = graph_control$max_pareto_k,
          use_psis = graph_control$use_psis,
          max_rounds = graph_control$max_rounds,
          target_cess = target_cess,
          n_mcmc_moves = n_mcmc_moves,
          max_steps = max_steps,
          n_cores = local_n_cores,
          seed = as.integer(seed + 500000L + 10000L * pos),
          verbose = verbose
        )
        bank <- certified$bank
        graph_insertions <- certified$insertions
      }
      list(bank = bank, rows = do.call(rbind, local_rows), graph_rows = graph_insertions)
    },
    mc.cores = as.integer(max(1L, n_jobs))
  )
  failed <- vapply(refined, inherits, logical(1), what = "try-error")
  if (any(failed)) {
    msg <- conditionMessage(attr(refined[[which(failed)[1L]]], "condition"))
    stop("Local bank design refinement failed for worker ", which(failed)[1L], ": ", msg)
  }
  for (i in seq_along(refined)) {
    banks[[i]] <- refined[[i]]$bank
    rows[[i]] <- refined[[i]]$rows
    graph_rows[[i]] <- refined[[i]]$graph_rows
  }
  graph_rows <- Filter(nrow, graph_rows)
  list(
    banks = banks,
    refinements = do.call(rbind, rows),
    graph_insertions = if (length(graph_rows)) do.call(rbind, graph_rows) else data.frame()
  )
}

bank_smc_audit_banks <- function(banks,
                                 theta_audit,
                                 population_model,
                                 target_ess_frac = 0.3,
                                 max_pareto_k = 0.7,
                                 local_ids = seq_along(banks),
                                 n_jobs = 1L) {
  model <- normalize_population_model(population_model)
  theta_audit <- .as_hyper_matrix(theta_audit, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  local_ids <- sort(unique(as.integer(local_ids)))
  local_ids <- local_ids[local_ids >= 1L & local_ids <= length(banks)]
  parts <- parallel::mclapply(
    local_ids,
    function(i) {
      bank <- validate_bank_smc_local_bank(banks[[i]], population_model = model)
      rows <- vector("list", nrow(theta_audit))
      for (j in seq_len(nrow(theta_audit))) {
        diag <- bank_smc_local_bank_coverage_diagnostic(
          bank,
          model,
          theta_audit[j, , drop = FALSE],
          min_mixture_ess_frac = target_ess_frac,
          max_pareto_k = max_pareto_k,
          use_psis = FALSE
        )
        rows[[j]] <- data.frame(
          local_id = as.integer(i),
          theta_id = as.integer(j),
          covered = isTRUE(diag$covered),
          ess_frac = as.numeric(diag$mixture$ess_frac),
          natural_ess_frac = as.numeric(diag$natural$best_ess_frac %||% NA_real_),
          natural_variance = as.numeric(diag$natural$best_variance %||% NA_real_),
          natural_distance = as.numeric(diag$natural$best_distance %||% NA_real_),
          max_weight = as.numeric(diag$mixture$max_weight),
          log_weight_var = as.numeric(diag$mixture$log_weight_var),
          best_node = as.integer(diag$best_node),
          nodes = as.integer(bank$n_nodes),
          particles = as.integer(bank$n_particles),
          check.names = FALSE
        )
      }
      do.call(rbind, rows)
    },
    mc.cores = as.integer(max(1L, n_jobs))
  )
  do.call(rbind, parts)
}

bank_smc_repair_from_audit <- function(banks,
                                       audit,
                                       theta_audit,
                                       data_list,
                                       loglik_fn,
                                       population_model,
                                       max_repairs = 20L,
                                       target_ess_frac = 0.3,
                                       max_pareto_k = 0.7,
                                       n_jobs = 1L,
                                       local_n_cores = 1L,
                                       seed = 123L,
                                       verbose = FALSE,
                                       graph_control = list(),
                                       ...) {
  failed <- audit[!audit$covered, , drop = FALSE]
  empty_repairs <- data.frame(
    local_id = integer(0),
    theta_id = integer(0),
    old_ess_frac = numeric(0),
    nodes_after = integer(0),
    particles_after = integer(0),
    check.names = FALSE
  )
  if (!nrow(failed)) {
    return(list(banks = banks, repairs = empty_repairs))
  }
  coverage_key <- if ("natural_ess_frac" %in% names(failed) && any(is.finite(failed$natural_ess_frac))) {
    ifelse(is.finite(failed$natural_ess_frac), failed$natural_ess_frac, Inf)
  } else {
    failed$ess_frac
  }
  failed <- failed[order(coverage_key, -failed$max_weight), , drop = FALSE]
  failed <- failed[seq_len(min(nrow(failed), as.integer(max_repairs))), , drop = FALSE]

  model <- normalize_population_model(population_model)
  graph_defaults <- list(
    enabled = TRUE,
    min_edge_ess_frac = 0.05,
    max_pareto_k = max_pareto_k,
    use_psis = FALSE,
    max_rounds = 50L
  )
  graph_control <- modifyList(graph_defaults, graph_control)
  dots <- list(...)
  failed_by_local <- split(failed, failed$local_id)
  repaired <- parallel::mclapply(
    failed_by_local,
    function(local_failed) {
      i <- as.integer(local_failed$local_id[1L])
      bank <- validate_bank_smc_local_bank(banks[[i]], population_model = model)
      local_rows <- vector("list", 0L)
      for (r in seq_len(nrow(local_failed))) {
        j <- local_failed$theta_id[r]
        diag <- bank_smc_local_bank_coverage_diagnostic(
          bank = bank,
          population_model = model,
          theta = theta_audit[j, , drop = FALSE],
          min_mixture_ess_frac = target_ess_frac,
          max_pareto_k = max_pareto_k,
          use_psis = FALSE
        )
        if (isTRUE(diag$covered)) next
        if (bank$n_nodes >= bank$max_nodes) next
        M <- as.integer(dots$M %||% bank$nodes[[diag$best_node]]$n_particles)
        if (is.finite(bank$max_particles) && bank$n_particles + M > bank$max_particles) next
        bridged <- do.call(
          bank_smc_local_bank_add_anchor,
          c(
            list(
              bank = bank,
              population_model = model,
              theta_target = theta_audit[j, , drop = FALSE],
              data_i = data_list[[i]],
              loglik_fn = loglik_fn,
              n_cores = local_n_cores,
              seed = as.integer(seed + 10000L * i + r),
              verbose = verbose,
              use_psis = FALSE
            ),
            dots
          )
        )
        bank <- bridged$bank
        local_rows[[length(local_rows) + 1L]] <- data.frame(
          local_id = as.integer(i),
          theta_id = as.integer(j),
          old_ess_frac = as.numeric(diag$mixture$ess_frac),
          nodes_after = as.integer(bridged$bank$n_nodes),
          particles_after = as.integer(bridged$bank$n_particles),
          check.names = FALSE
        )
      }
      graph_insertions <- data.frame()
      if (isTRUE(graph_control$enabled) && bank$n_nodes > 1L) {
        certified <- bank_smc_local_bank_certify_bridge_graph(
          bank = bank,
          population_model = model,
          data_i = data_list[[i]],
          loglik_fn = loglik_fn,
          M = dots$M %||% NULL,
          min_edge_ess_frac = graph_control$min_edge_ess_frac,
          max_pareto_k = graph_control$max_pareto_k,
          use_psis = graph_control$use_psis,
          max_rounds = graph_control$max_rounds,
          target_cess = dots$target_cess %||% 0.9,
          n_mcmc_moves = dots$n_mcmc_moves %||% 2L,
          max_steps = dots$max_steps %||% 128L,
          n_cores = local_n_cores,
          seed = as.integer(seed + 600000L + 10000L * i),
          verbose = verbose
        )
        bank <- certified$bank
        graph_insertions <- certified$insertions
      }
      list(
        bank = bank,
        rows = if (length(local_rows)) do.call(rbind, local_rows) else empty_repairs,
        graph_rows = graph_insertions
      )
    },
    mc.cores = as.integer(max(1L, n_jobs))
  )
  failed_workers <- vapply(repaired, inherits, logical(1), what = "try-error")
  if (any(failed_workers)) {
    msg <- conditionMessage(attr(repaired[[which(failed_workers)[1L]]], "condition"))
    stop("Local bank audit repair failed for worker ", which(failed_workers)[1L], ": ", msg)
  }
  repair_rows <- vector("list", length(repaired))
  graph_rows <- vector("list", length(repaired))
  for (k in seq_along(repaired)) {
    i <- as.integer(names(failed_by_local)[k])
    banks[[i]] <- repaired[[k]]$bank
    repair_rows[[k]] <- repaired[[k]]$rows
    graph_rows[[k]] <- repaired[[k]]$graph_rows
  }
  repair_rows <- Filter(nrow, repair_rows)
  graph_rows <- Filter(nrow, graph_rows)
  list(
    banks = banks,
    repairs = if (length(repair_rows)) do.call(rbind, repair_rows) else empty_repairs,
    graph_insertions = if (length(graph_rows)) do.call(rbind, graph_rows) else data.frame()
  )
}

fit_bank_smc_population_model <- function(data_list,
                                          loglik_fn,
                                          population_model,
                                          initial_theta,
                                          local_control = list(),
                                          anchor_control = list(),
                                          shape_control = list(),
                                          outer_control = list(),
                                          design_control = list(),
                                          calibration_control = list(),
                                          graph_control = list(),
                                          surface_control = list(),
                                          audit_control = list(),
                                          n_cores = 1L,
                                          seed = 123L,
                                          verbose = TRUE) {
  model <- normalize_population_model(population_model)
  local_defaults <- list(
    M = 600L,
    max_nodes = 5L,
    max_particles = Inf,
    target_ess_frac = 0.3,
    bridge_particles = NULL,
    target_cess = 0.9,
    n_mcmc_moves = 2L,
    max_steps = 128L,
    bridge_min_single_ess_frac = Inf,
    bridge_max_pareto_k = 0.7
  )
  anchor_defaults <- list(
    enabled = TRUE,
    M = NULL,
    n_alpha_configs = 128L,
    scales = c(16, 4),
    coherent_probs = c(0.025, 0.05, 0.1, 0.2, 0.8, 0.9, 0.95, 0.975),
    anchor_weight = 0.70,
    defensive = TRUE,
    defensive_scale = 64,
    defensive_weight = 0.20,
    draws_per_config = 1L,
    ell_grid_size = 96L,
    n_mcmc_moves = 2L,
    max_rounds = 80L,
    cess_target = 0.9,
    resample_threshold = 0.6,
    G_mix = 8L,
    rw_scale = 0.9
  )
  design_defaults <- list(
    support_points = 6L,
    profile_dims = 4L,
    tail_probs = c(0.025, 0.1, 0.9, 0.975),
    initial_tail_probs = c(0.01, 0.025, 0.1, 0.9, 0.975, 0.99),
    paired_effect_profiles = TRUE
  )
  shape_defaults <- list(n_alpha_configs = 128L, draws_per_config = 1L, ell_grid_size = 96L)
  calibration_defaults <- list(
    enabled = TRUE,
    max_iter = 3000L,
    tol = 1e-8,
    anchor = "mean_smc"
  )
  graph_defaults <- list(
    enabled = TRUE,
    min_edge_ess_frac = 0.05,
    max_pareto_k = 0.7,
    use_psis = FALSE,
    max_rounds = 50L,
    require_connected = TRUE
  )
  surface_defaults <- list(
    enabled = FALSE,
    max_points = 0L,
    refine_rounds = 1L,
    profile_dims = 4L,
    tail_probs = c(0.05, 0.95),
    paired_mean_tail_probs = NULL,
    paired_variance_tail_probs = NULL,
    inflation_scale = 2.5,
    include_inflated_profile = TRUE,
    include_posterior_tail = TRUE,
    paired_effect_profiles = TRUE,
    paired_effect_grid = TRUE,
    candidate_multiplier = 4L,
    score_local_count = 20L,
    score_log_weight = 1,
    score_score_weight = 0.25,
    score_distance_weight = 0.05
  )
  audit_defaults <- list(
    max_points = NULL,
    target_ess_frac = 0.05,
    max_pareto_k = 0.7,
    max_repairs = NULL,
    refine_rounds = 3L,
    force_refine_points = 0L,
    force_profile_dims = 1L,
    force_tail_probs = c(0.01, 0.05, 0.95, 0.99)
  )
  local_control <- modifyList(local_defaults, local_control)
  anchor_control <- modifyList(anchor_defaults, anchor_control)
  anchor_control$M <- as.integer(anchor_control$M %||% min(400L, local_control$M))
  shape_control <- modifyList(shape_defaults, shape_control)
  calibration_control <- modifyList(calibration_defaults, calibration_control)
  graph_control <- modifyList(graph_defaults, graph_control)
  surface_control <- modifyList(surface_defaults, surface_control)
  design_control <- modifyList(design_defaults, design_control)
  surface_control$max_points <- as.integer(surface_control$max_points %||% 0L)
  surface_control$refine_rounds <- as.integer(surface_control$refine_rounds %||% 0L)
  surface_control$profile_dims <- as.integer(surface_control$profile_dims %||% 1L)
  surface_control$enabled <- isTRUE(surface_control$enabled) && surface_control$max_points > 0L
  node_budget <- if (is.finite(local_control$max_nodes)) as.integer(local_control$max_nodes) else 5L
  design_control$max_points <- as.integer(
    design_control$max_points %||%
      max(2L, node_budget - 1L)
  )
  audit_control <- modifyList(audit_defaults, audit_control)
  audit_control$max_points <- as.integer(
    audit_control$max_points %||%
      max(5L, node_budget - 1L)
  )
  audit_control$max_repairs <- as.integer(
    audit_control$max_repairs %||%
      max(20L, length(data_list) * audit_control$max_points)
  )

  if (isTRUE(verbose)) cat("Bank SMC: initial local banks\n")
  banks <- bank_smc_initial_banks(
    data_list = data_list,
    loglik_fn = loglik_fn,
    population_model = model,
    initial_theta = initial_theta,
    M = local_control$M,
    max_nodes = local_control$max_nodes,
    max_particles = local_control$max_particles,
    cess_target = local_control$target_cess,
    n_mcmc_moves = local_control$n_mcmc_moves,
    max_rounds = local_control$max_steps,
    n_jobs = n_cores,
    local_n_cores = 1L,
    seed = seed,
    verbose = FALSE
  )

  anchor_candidates <- NULL
  if (isTRUE(anchor_control$enabled)) {
    if (isTRUE(verbose)) cat("Bank SMC: alpha-first anchor candidates\n")
    anchor_candidates <- bank_smc_alpha_first_anchor_candidates(
      data_list = data_list,
      loglik_fn = loglik_fn,
      population_model = model,
      theta_reference = initial_theta,
      n_alpha_configs = anchor_control$n_alpha_configs,
      draws_per_config = anchor_control$draws_per_config,
      ell_grid_size = anchor_control$ell_grid_size,
      M = anchor_control$M,
      scales = anchor_control$scales,
      coherent_probs = anchor_control$coherent_probs,
      defensive = anchor_control$defensive,
      defensive_scale = anchor_control$defensive_scale,
      defensive_weight = anchor_control$defensive_weight,
      resample_threshold = anchor_control$resample_threshold,
      n_mcmc_moves = anchor_control$n_mcmc_moves,
      max_rounds = anchor_control$max_rounds,
      cess_target = anchor_control$cess_target,
      G_mix = anchor_control$G_mix,
      rw_scale = anchor_control$rw_scale,
      n_jobs = n_cores,
      local_n_cores = 1L,
      seed = seed + 50000L,
      verbose = FALSE
    )
  }

  if (isTRUE(verbose)) cat("Bank SMC: bank-imputed population shape\n")
  shape_fit <- bank_smc_build_analytic_shape(
    banks = banks,
    population_model = model,
    theta_reference = initial_theta,
    n_alpha_configs = shape_control$n_alpha_configs,
    draws_per_config = shape_control$draws_per_config,
    ell_grid_size = shape_control$ell_grid_size,
    seed = seed + 100000L
  )

  design_fit <- bank_smc_anchor_design_fit(
    anchor_candidates = anchor_candidates,
    reference_fit = shape_fit,
    population_model = model,
    anchor_weight = anchor_control$anchor_weight
  )
  initial_design_control <- design_control
  initial_design_control$tail_probs <- initial_design_control$initial_tail_probs
  initial_design_control$initial_tail_probs <- NULL
  theta_design <- do.call(
    bank_smc_select_theta_design,
    c(
      list(
        population_fit = design_fit,
        population_model = model,
        banks = banks,
        initial_theta = initial_theta,
        target_ess_frac = local_control$target_ess_frac
      ),
      initial_design_control
    )
  )

  if (isTRUE(verbose)) cat(sprintf("Bank SMC: initial design expansion (%d theta points)\n", nrow(theta_design)))
  design_refinement <- bank_smc_refine_banks_to_design(
    banks = banks,
    theta_design = theta_design,
    data_list = data_list,
    loglik_fn = loglik_fn,
    population_model = model,
    target_ess_frac = local_control$target_ess_frac,
    bridge_particles = local_control$bridge_particles,
    target_cess = local_control$target_cess,
    n_mcmc_moves = local_control$n_mcmc_moves,
    max_steps = local_control$max_steps,
    bridge_min_single_ess_frac = local_control$bridge_min_single_ess_frac,
    bridge_max_pareto_k = local_control$bridge_max_pareto_k,
    force = TRUE,
    graph_control = graph_control,
    n_jobs = n_cores,
    local_n_cores = 1L,
    seed = seed + 200000L,
    verbose = FALSE
  )
  banks <- design_refinement$banks
  graph_insertions <- list(initial_design = design_refinement$graph_insertions)

  if (isTRUE(graph_control$enabled) && isTRUE(graph_control$require_connected)) {
    graph_summary <- bank_smc_overlap_graph_summary(banks)
    if (any(!graph_summary$certified)) {
      stop(
        "Local bank overlap graph certification failed before outer SMC for ",
        sum(!graph_summary$certified),
        " locals."
      )
    }
  }

  audits <- list()
  repairs <- list()
  surface_refinements <- list()
  surface_designs <- list()
  surface_scores <- list()
  factor_updates <- list()
  surface_rounds_done <- 0L
  for (round in seq_len(as.integer(audit_control$refine_rounds) + 1L)) {
    if (isTRUE(graph_control$enabled) && isTRUE(graph_control$require_connected)) {
      graph_summary <- bank_smc_overlap_graph_summary(banks)
      if (any(!graph_summary$certified)) {
        stop(
          "Local bank overlap graph certification failed before outer SMC for ",
          sum(!graph_summary$certified),
          " locals."
        )
      }
    }
    if (isTRUE(calibration_control$enabled)) {
      if (isTRUE(verbose)) cat("Bank SMC: calibrating local normalizers\n")
      banks <- bank_smc_calibrate_banks_normalizers(
        banks = banks,
        population_model = model,
        max_iter = calibration_control$max_iter,
        tol = calibration_control$tol,
        anchor = calibration_control$anchor,
        n_jobs = n_cores
      )
    }
    if (isTRUE(verbose)) {
      label <- if (round == 1L) "outer population fit" else "outer population refit"
      cat(sprintf("Bank SMC: %s\n", label))
    }
    factor_set <- bank_smc_banks_to_factor_set(banks, model)
    fit <- do.call(
      outer_population_smc,
      modifyList(
        list(
          factor_set = factor_set,
          N = 1200L,
          n_mcmc_moves = 3L,
          max_rounds = 80L,
          n_cores = n_cores,
          seed = seed + 300000L + 10000L * (round - 1L),
          verbose = FALSE
        ),
        outer_control
      )
    )

    if (isTRUE(surface_control$enabled) &&
        surface_rounds_done < as.integer(surface_control$refine_rounds) &&
        round <= as.integer(audit_control$refine_rounds)) {
      surface_selection <- bank_smc_select_surface_theta(
        population_fit = fit,
        population_model = model,
        max_points = surface_control$max_points,
        profile_dims = surface_control$profile_dims,
        tail_probs = surface_control$tail_probs,
        paired_mean_tail_probs = surface_control$paired_mean_tail_probs,
        paired_variance_tail_probs = surface_control$paired_variance_tail_probs,
        inflation_scale = surface_control$inflation_scale,
        include_inflated_profile = surface_control$include_inflated_profile,
        include_posterior_tail = surface_control$include_posterior_tail,
        paired_effect_profiles = surface_control$paired_effect_profiles,
        paired_effect_grid = surface_control$paired_effect_grid,
        candidate_multiplier = surface_control$candidate_multiplier,
        banks = banks,
        score_local_count = surface_control$score_local_count,
        score_log_weight = surface_control$score_log_weight,
        score_score_weight = surface_control$score_score_weight,
        score_distance_weight = surface_control$score_distance_weight,
        n_jobs = n_cores,
        return_diagnostics = TRUE
      )
      surface_theta <- surface_selection$theta
      surface_designs[[paste0("surface_", round)]] <- surface_theta
      surface_scores[[paste0("surface_", round)]] <- surface_selection$score
      if (nrow(surface_theta)) {
        surface_refinement <- bank_smc_refine_banks_to_design(
          banks = banks,
          theta_design = surface_theta,
          data_list = data_list,
          loglik_fn = loglik_fn,
          population_model = model,
          target_ess_frac = audit_control$target_ess_frac,
          bridge_particles = local_control$bridge_particles,
          target_cess = local_control$target_cess,
          n_mcmc_moves = local_control$n_mcmc_moves,
          max_steps = local_control$max_steps,
          bridge_min_single_ess_frac = local_control$bridge_min_single_ess_frac,
          bridge_max_pareto_k = local_control$bridge_max_pareto_k,
          graph_control = graph_control,
          force = TRUE,
          n_jobs = n_cores,
          local_n_cores = 1L,
          seed = seed + 325000L + round,
          verbose = FALSE
        )
        surface_rounds_done <- surface_rounds_done + 1L
        surface_added <- surface_refinement$refinements[
          as.logical(surface_refinement$refinements$added), ,
          drop = FALSE
        ]
        graph_added <- surface_refinement$graph_insertions
        changed_locals <- sort(unique(c(
          surface_added$local_id,
          graph_added$local_id
        )))
        surface_refinements[[paste0("surface_", round)]] <- surface_refinement$refinements
        if (length(changed_locals)) {
          banks <- surface_refinement$banks
          graph_insertions[[paste0("surface_", round)]] <- surface_refinement$graph_insertions
          factor_updates[[paste0("surface_", round)]] <- list(
            local_ids = changed_locals,
            action = "surface_profile_anchor"
          )
          if (isTRUE(verbose)) {
            finite_error <- surface_added$predicted_anchor_error[
              is.finite(surface_added$predicted_anchor_error)
            ]
            max_error <- if (length(finite_error)) max(abs(finite_error)) else NA_real_
            cat(sprintf(
              "Bank SMC: forced %d surface anchors | max prediction error %.3f\n",
              nrow(surface_added),
              max_error
            ))
          }
          next
        }
      }
    }

    if (as.integer(audit_control$force_refine_points) > 0L &&
        round <= as.integer(audit_control$refine_rounds)) {
      force_theta <- bank_smc_select_joint_tail_theta(
        population_fit = fit,
        population_model = model,
        max_points = audit_control$force_refine_points,
        profile_dims = audit_control$force_profile_dims,
        tail_probs = audit_control$force_tail_probs
      )
      if (nrow(force_theta)) {
        forced_refinement <- bank_smc_refine_banks_to_design(
          banks = banks,
          theta_design = force_theta,
          data_list = data_list,
          loglik_fn = loglik_fn,
          population_model = model,
          target_ess_frac = audit_control$target_ess_frac,
          bridge_particles = local_control$bridge_particles,
          target_cess = local_control$target_cess,
          n_mcmc_moves = local_control$n_mcmc_moves,
          max_steps = local_control$max_steps,
          bridge_min_single_ess_frac = local_control$bridge_min_single_ess_frac,
          bridge_max_pareto_k = local_control$bridge_max_pareto_k,
          graph_control = graph_control,
          force = TRUE,
          n_jobs = n_cores,
          local_n_cores = 1L,
          seed = seed + 350000L + round,
          verbose = FALSE
        )
        forced_added <- forced_refinement$refinements[
          as.logical(forced_refinement$refinements$added), ,
          drop = FALSE
        ]
        graph_added <- forced_refinement$graph_insertions
        changed_locals <- sort(unique(c(
          forced_added$local_id,
          graph_added$local_id
        )))
        if (length(changed_locals)) {
          banks <- forced_refinement$banks
          repairs[[paste0("force_", round)]] <- forced_refinement$refinements
          graph_insertions[[paste0("force_", round)]] <- forced_refinement$graph_insertions
          factor_updates[[paste0("force_", round)]] <- list(
            local_ids = changed_locals,
            action = "force_outer_tail_anchor"
          )
          if (isTRUE(verbose)) {
            cat(sprintf("Bank SMC: forced %d posterior-boundary anchors\n", nrow(forced_added)))
          }
          next
        }
      }
    }

    audit_design_control <- modifyList(design_control, list(max_points = audit_control$max_points))
    audit_design_control$initial_tail_probs <- NULL
    theta_audit <- do.call(
      bank_smc_select_theta_design,
      c(
        list(
          population_fit = fit,
          population_model = model,
          banks = banks,
          target_ess_frac = audit_control$target_ess_frac
        ),
        audit_design_control
      )
    )
    audit <- bank_smc_audit_banks(
      banks = banks,
      theta_audit = theta_audit,
      population_model = model,
      target_ess_frac = audit_control$target_ess_frac,
      max_pareto_k = audit_control$max_pareto_k,
      n_jobs = n_cores
    )
    audits[[round]] <- list(theta = theta_audit, audit = audit)
    if (!any(!audit$covered) || round > as.integer(audit_control$refine_rounds)) {
      break
    }
    if (isTRUE(verbose)) cat(sprintf("Bank SMC: repairing %d audit failures\n", sum(!audit$covered)))
    repair <- bank_smc_repair_from_audit(
      banks = banks,
      audit = audit,
      theta_audit = theta_audit,
      data_list = data_list,
      loglik_fn = loglik_fn,
      population_model = model,
      max_repairs = audit_control$max_repairs,
      target_ess_frac = audit_control$target_ess_frac,
      max_pareto_k = audit_control$max_pareto_k,
      M = local_control$bridge_particles,
      target_cess = local_control$target_cess,
      n_mcmc_moves = local_control$n_mcmc_moves,
      max_steps = local_control$max_steps,
      bridge_min_single_ess_frac = local_control$bridge_min_single_ess_frac,
      bridge_max_pareto_k = local_control$bridge_max_pareto_k,
      n_jobs = n_cores,
      local_n_cores = 1L,
      seed = seed + 400000L + round,
      verbose = FALSE,
      graph_control = graph_control
    )
    banks <- repair$banks
    repairs[[paste0("repair_", round)]] <- repair$repairs
    graph_insertions[[paste0("repair_", round)]] <- repair$graph_insertions
    changed <- sort(unique(c(repair$repairs$local_id, repair$graph_insertions$local_id)))
    if (!length(changed)) break
    factor_updates[[paste0("repair_", round)]] <- list(local_ids = changed, action = "rerun_outer")
  }

  list(
    banks = banks,
    factor_set = factor_set,
    shape_fit = shape_fit,
    anchor_candidates = anchor_candidates,
    fit = fit,
    theta_design = theta_design,
    design_refinement = design_refinement$refinements,
    audits = audits,
    repairs = repairs,
    surface_designs = surface_designs,
    surface_scores = surface_scores,
    surface_refinements = surface_refinements,
    graph_insertions = graph_insertions,
    factor_updates = factor_updates,
    settings = list(
      local_control = local_control,
      shape_control = shape_control,
      anchor_control = anchor_control,
      design_control = design_control,
      calibration_control = calibration_control,
      graph_control = graph_control,
      surface_control = surface_control,
      audit_control = audit_control,
      seed = seed
    )
  )
}
