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

bank_smc_node_log_prior <- function(node, population_model, theta) {
  node <- validate_bank_smc_node(node, population_model = population_model)
  population_model_log_alpha_given_theta(
    population_model,
    alpha = node$alpha,
    theta = theta
  )
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

bank_smc_local_bank_node_count <- function(bank) {
  validate_bank_smc_local_bank(bank)$n_nodes
}

bank_smc_local_bank_particle_count <- function(bank) {
  validate_bank_smc_local_bank(bank)$n_particles
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

bank_smc_local_bank_replace_node <- function(bank,
                                             index,
                                             node,
                                             eta = NULL,
                                             population_model = NULL) {
  bank <- validate_bank_smc_local_bank(bank, population_model = population_model)
  index <- as.integer(index)
  if (length(index) != 1L || index < 1L || index > length(bank$nodes)) {
    stop("Replacement index is out of range.")
  }
  node <- validate_bank_smc_node(node, population_model = population_model)
  if (is.finite(bank$local_id) && is.finite(node$local_id) && node$local_id != bank$local_id) {
    stop("Replacement node local_id does not match the bank local_id.")
  }
  if (!is.finite(node$local_id)) node$local_id <- bank$local_id

  bank$nodes[[index]] <- node
  if (!is.null(eta)) {
    bank$eta[index] <- as.numeric(eta)
  }
  bank$stack <- NULL
  validate_bank_smc_local_bank(bank, population_model = population_model)
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

bank_smc_local_bank_with_stack <- function(bank, population_model) {
  bank <- validate_bank_smc_local_bank(bank, population_model = population_model)
  stack <- bank$stack
  if (!inherits(stack, "bank_smc_local_bank_stack") ||
      is.null(stack$alpha) ||
      is.null(stack$log_base) ||
      nrow(stack$alpha) != bank$n_particles ||
      length(stack$log_base) != bank$n_particles) {
    stack <- bank_smc_local_bank_build_stack(bank, population_model)
  }
  bank$stack <- stack
  bank
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
  data.frame(
    ess = summary$ess,
    ess_frac = summary$ess_frac,
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

  best <- if (!is.null(natural)) natural$best_node else which.max(single$ess_frac)
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

  if (fit$final_lambda < 1 - 1e-8) {
    stop("Anchor SMC did not reach the local posterior within max_steps.")
  }

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

bank_smc_local_bank_add_anchor <- function(bank,
                                           population_model,
                                           theta_target,
                                           data_i,
                                           loglik_fn,
                                           M = NULL,
                                           eta = NULL,
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
  M <- as.integer(M %||% bank$nodes[[coverage$best_node]]$n_particles)
  node <- bank_smc_anchor_node(
    local_id = bank$local_id,
    population_model = model,
    theta_anchor = theta_target,
    data_i = data_i,
    loglik_fn = loglik_fn,
    M = M,
    ...
  )
  updated <- bank_smc_local_bank_add_node(
    bank = bank,
    node = node,
    eta = eta,
    population_model = model
  )

  list(
    bank = updated,
    node = node,
    start_node_index = as.integer(coverage$best_node),
    start_coverage = coverage
  )
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
      for (val in vals) {
        row <- center
        row[j] <- val
        rows[[length(rows) + 1L]] <- row
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
                                            n_jobs = 1L,
                                            local_n_cores = 1L,
                                            seed = 123L,
                                            verbose = FALSE) {
  model <- normalize_population_model(population_model)
  theta_design <- .as_hyper_matrix(theta_design, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
  ids <- seq_along(banks)
  rows <- vector("list", 0L)
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
        added <- FALSE
        if (!isTRUE(diag$covered) && bank$n_nodes < bank$max_nodes) {
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
              n_cores = local_n_cores,
              seed = as.integer(seed + 10000L * pos + j),
              verbose = verbose,
              use_psis = FALSE
            )
            bank <- bridged$bank
            added <- TRUE
          }
        }
        local_rows[[length(local_rows) + 1L]] <- data.frame(
          local_id = as.integer(bank$local_id),
          theta_id = as.integer(j),
          covered_before = isTRUE(diag$covered),
          added = isTRUE(added),
          ess_frac = as.numeric(diag$mixture$ess_frac),
          natural_ess_frac = as.numeric(diag$natural$best_ess_frac %||% NA_real_),
          natural_variance = as.numeric(diag$natural$best_variance %||% NA_real_),
          best_node = as.integer(diag$best_node),
          nodes_after = as.integer(bank$n_nodes),
          particles_after = as.integer(bank$n_particles),
          check.names = FALSE
        )
      }
      list(bank = bank, rows = do.call(rbind, local_rows))
    },
    mc.cores = as.integer(max(1L, n_jobs))
  )
  for (i in seq_along(refined)) {
    banks[[i]] <- refined[[i]]$bank
    rows[[i]] <- refined[[i]]$rows
  }
  list(
    banks = banks,
    refinements = do.call(rbind, rows)
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
  dots <- list(...)
  repair_rows <- vector("list", 0L)
  for (r in seq_len(nrow(failed))) {
    i <- failed$local_id[r]
    j <- failed$theta_id[r]
    bank <- validate_bank_smc_local_bank(banks[[i]], population_model = model)
    if (bank$n_nodes >= bank$max_nodes) next
    M <- as.integer(dots$M %||% bank$nodes[[failed$best_node[r]]]$n_particles)
    if (is.finite(bank$max_particles) && bank$n_particles + M > bank$max_particles) next
    bridged <- bank_smc_local_bank_add_anchor(
      bank = bank,
      population_model = model,
      theta_target = theta_audit[j, , drop = FALSE],
      data_i = data_list[[i]],
      loglik_fn = loglik_fn,
      ...,
      use_psis = FALSE
    )
    banks[[i]] <- bridged$bank
    repair_rows[[r]] <- data.frame(
      local_id = as.integer(i),
      theta_id = as.integer(j),
      old_ess_frac = failed$ess_frac[r],
      nodes_after = as.integer(bridged$bank$n_nodes),
      particles_after = as.integer(bridged$bank$n_particles),
      check.names = FALSE
    )
  }
  list(
    banks = banks,
    repairs = if (length(repair_rows)) do.call(rbind, repair_rows) else empty_repairs
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
    n_mcmc_moves = 2L
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
    initial_tail_probs = c(0.01, 0.025, 0.1, 0.9, 0.975, 0.99)
  )
  shape_defaults <- list(n_alpha_configs = 128L, draws_per_config = 1L, ell_grid_size = 96L)
  audit_defaults <- list(
    max_points = 5L,
    target_ess_frac = 0.05,
    max_pareto_k = 0.7,
    max_repairs = 20L,
    refine_rounds = 3L
  )
  local_control <- modifyList(local_defaults, local_control)
  anchor_control <- modifyList(anchor_defaults, anchor_control)
  anchor_control$M <- as.integer(anchor_control$M %||% min(400L, local_control$M))
  shape_control <- modifyList(shape_defaults, shape_control)
  design_control <- modifyList(design_defaults, design_control)
  node_budget <- if (is.finite(local_control$max_nodes)) as.integer(local_control$max_nodes) else 5L
  design_control$max_points <- as.integer(
    design_control$max_points %||%
      max(2L, min(5L, node_budget - 1L))
  )
  audit_control <- modifyList(audit_defaults, audit_control)

  if (isTRUE(verbose)) cat("Bank SMC: initial local banks\n")
  banks <- bank_smc_initial_banks(
    data_list = data_list,
    loglik_fn = loglik_fn,
    population_model = model,
    initial_theta = initial_theta,
    M = local_control$M,
    max_nodes = local_control$max_nodes,
    max_particles = local_control$max_particles,
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
    n_jobs = n_cores,
    local_n_cores = 1L,
    seed = seed + 200000L,
    verbose = FALSE
  )
  banks <- design_refinement$banks

  if (isTRUE(verbose)) cat("Bank SMC: outer population fit\n")
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
        seed = seed + 300000L,
        verbose = FALSE
      ),
      outer_control
    )
  )

  audits <- list()
  repairs <- list()
  factor_updates <- list()
  for (round in seq_len(as.integer(audit_control$refine_rounds) + 1L)) {
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
      M = local_control$bridge_particles,
      target_cess = local_control$target_cess,
      n_mcmc_moves = local_control$n_mcmc_moves,
      n_cores = 1L,
      seed = seed + 400000L + round,
      verbose = FALSE
    )
    banks <- repair$banks
    repairs[[round]] <- repair$repairs
    changed <- sort(unique(repair$repairs$local_id))
    if (!length(changed)) break

    old_factor_set <- factor_set
    local_objects <- bank_smc_banks_to_local_objects(banks, model)
    factor_set <- update_population_factor_set_locals(
      factor_set = factor_set,
      local_objects = local_objects,
      local_ids = changed
    )
    fit <- update_outer_population_fit(
      fit = fit,
      old_factor_set = old_factor_set,
      new_factor_set = factor_set,
      n_mcmc_moves = as.integer(outer_control$n_mcmc_moves %||% 3L),
      min_mcmc_moves = as.integer(outer_control$min_mcmc_moves %||% 1L),
      n_cores = n_cores,
      seed = seed + 500000L + round,
      verbose = FALSE
    )
    factor_updates[[round]] <- list(local_ids = changed)
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
    factor_updates = factor_updates,
    settings = list(
      local_control = local_control,
      shape_control = shape_control,
      anchor_control = anchor_control,
      design_control = design_control,
      audit_control = audit_control,
      seed = seed
    )
  )
}
