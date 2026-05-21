#!/usr/bin/env Rscript
# ============================================================================
# Theta-level local-SMC correction observations
# - Phase 5 support-engine building block
# - Evaluates exact local SMC evidence at selected theta support points
# ============================================================================

if (!exists("%||%", mode = "function") ||
    !exists("ll_parallel", mode = "function") ||
    !exists("ESS", mode = "function") ||
    !exists("weighted_cov", mode = "function") ||
    !exists("regularize_cov", mode = "function")) {
  source("smc_core.R")
}

if (!exists("run_tempered_smc", mode = "function")) {
  source("SMC_super_fast.R")
}

if (!exists("make_reference_prior_gaussian", mode = "function")) {
  source("reference_priors.R")
}

if (!exists("normalize_population_model", mode = "function") ||
    !exists("population_model_reference_components_from_theta", mode = "function")) {
  source("population_models.R")
}

if (!exists("population_factor_set_loglik_by_local", mode = "function")) {
  source("outer_population_smc.R")
}

if (!exists("theta_pathfinder_whitener", mode = "function") &&
    file.exists("theta_support_pathfinder.R")) {
  source("theta_support_pathfinder.R")
}

suppressPackageStartupMessages({
  library(parallel)
})

.theta_corr_as_matrix <- function(theta, model) {
  .as_hyper_matrix(theta, hyper_names = model$hyper_names, hyper_dim = model$hyper_dim)
}

.theta_corr_normalize_weights <- function(w, n) {
  if (is.null(w)) return(rep(1 / n, n))
  if (length(w) != n) stop("weights length does not match theta rows.")
  w <- pmax(as.numeric(w), 0)
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) rep(1 / n, n) else w / sw
}

.theta_corr_weighted_quantile <- function(x, w, probs) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) {
    return(as.numeric(stats::quantile(as.numeric(x), probs = probs, na.rm = TRUE, names = FALSE)))
  }
  x <- as.numeric(x[ok])
  w <- as.numeric(w[ok])
  ord <- order(x)
  x <- x[ord]
  w <- .theta_corr_normalize_weights(w[ord], length(w))
  cw <- c(0, cumsum(w))
  qx <- c(x[1L], x)
  keep <- !duplicated(cw)
  as.numeric(stats::approx(cw[keep], qx[keep], xout = probs, rule = 2)$y)
}

.theta_corr_prob_label <- function(prob) {
  gsub("\\s+", "", formatC(100 * as.numeric(prob), format = "fg", digits = 3))
}

.theta_corr_support_whitener <- function(theta, w, model) {
  if (exists("theta_pathfinder_whitener", mode = "function")) {
    theta_pathfinder_whitener(theta, w = w, model = model)
  } else {
    theta <- .theta_corr_as_matrix(theta, model)
    w <- .theta_corr_normalize_weights(w, nrow(theta))
    center <- colSums(theta * w)
    S <- tryCatch(weighted_cov(theta, w), error = function(e) stats::cov(theta))
    S <- regularize_cov(S, min_eig = 1e-8, cond_cap = 1e8)
    list(center = center, cov = S, chol = chol(S), names = model$hyper_names)
  }
}

.theta_corr_to_z <- function(theta, whitener, model) {
  theta <- .theta_corr_as_matrix(theta, model)
  centered <- sweep(theta[, whitener$names, drop = FALSE], 2L, whitener$center, "-")
  z <- t(backsolve(whitener$chol, t(centered), transpose = TRUE))
  colnames(z) <- whitener$names
  z
}

.theta_corr_unique_design <- function(design, theta, model, support_theta, support_w, min_distance = 1e-6) {
  if (!nrow(design)) return(design)
  whitener <- .theta_corr_support_whitener(support_theta, support_w, model)
  z <- .theta_corr_to_z(theta, whitener, model)
  keep <- logical(nrow(z))
  kept <- NULL
  for (i in seq_len(nrow(z))) {
    zi <- as.numeric(z[i, ])
    if (is.null(kept)) {
      keep[i] <- TRUE
      kept <- matrix(zi, nrow = 1L)
      next
    }
    dist <- sqrt(rowSums(sweep(kept, 2L, zi, "-")^2))
    if (min(dist) > min_distance) {
      keep[i] <- TRUE
      kept <- rbind(kept, zi)
    }
  }
  design[keep, , drop = FALSE]
}

select_theta_correction_points <- function(final_fit,
                                           population_model = final_fit$population_model %||% final_fit$factor_set$population_model,
                                           support_theta = NULL,
                                           support_w = NULL,
                                           pathfinder = NULL,
                                           tail_probs = c(0.05, 0.95)) {
  model <- normalize_population_model(population_model)
  theta <- .theta_corr_as_matrix(final_fit$theta, model)
  w <- .theta_corr_normalize_weights(final_fit$w, nrow(theta))
  support_theta <- .theta_corr_as_matrix(support_theta %||% theta, model)
  support_w <- .theta_corr_normalize_weights(support_w, nrow(support_theta))

  rows <- list()
  theta_rows <- list()
  point_id <- 0L

  add_point <- function(source, role, label, theta_row) {
    point_id <<- point_id + 1L
    theta_row <- .theta_corr_as_matrix(theta_row, model)[1L, , drop = FALSE]
    theta_rows[[point_id]] <<- theta_row
    rows[[point_id]] <<- data.frame(
      point_id = point_id,
      source = source,
      role = role,
      label = label,
      matrix(theta_row[1L, ], nrow = 1L, dimnames = list(NULL, model$hyper_names)),
      check.names = FALSE
    )
  }

  center <- colSums(theta * w)
  z_wh <- .theta_corr_support_whitener(theta, w, model)
  z <- .theta_corr_to_z(theta, z_wh, model)
  center_z <- .theta_corr_to_z(matrix(center, nrow = 1L, dimnames = list(NULL, model$hyper_names)), z_wh, model)
  center_idx <- which.min(rowSums(sweep(z, 2L, center_z[1L, ], "-")^2))
  add_point("rho_smc", "center", "center", theta[center_idx, , drop = FALSE])

  tail_probs <- sort(unique(as.numeric(tail_probs)))
  for (nm in model$hyper_names) {
    qs <- .theta_corr_weighted_quantile(theta[, nm], w, tail_probs)
    for (j in seq_along(tail_probs)) {
      idx <- which.min(abs(theta[, nm] - qs[j]))
      side <- if (tail_probs[j] < 0.5) "lower" else "upper"
      add_point(
        "rho_smc",
        paste0("marginal_", side),
        paste0(nm, "_q", .theta_corr_prob_label(tail_probs[j])),
        theta[idx, , drop = FALSE]
      )
    }
  }

  if (!is.null(pathfinder) &&
      !is.null(pathfinder$candidates) &&
      nrow(pathfinder$candidates) &&
      "nonduplicate_support" %in% names(pathfinder$candidates)) {
    candidates <- pathfinder$candidates[pathfinder$candidates$nonduplicate_support, , drop = FALSE]
    if (nrow(candidates)) {
      for (j in seq_len(nrow(candidates))) {
        label <- paste0(
          "pathfinder_",
          candidates$path_id[j],
          "_",
          candidates$point[j],
          "_p",
          .theta_corr_prob_label(candidates$coverage_prob[j])
        )
        add_point(
          "pathfinder",
          "stress",
          label,
          candidates[j, model$hyper_names, drop = FALSE]
        )
      }
    }
  }

  design <- do.call(rbind, rows)
  theta_mat <- do.call(rbind, theta_rows)
  design <- .theta_corr_unique_design(
    design = design,
    theta = theta_mat,
    model = model,
    support_theta = support_theta,
    support_w = support_w
  )
  design$point_id <- seq_len(nrow(design))
  rownames(design) <- NULL
  structure(design, class = c("theta_correction_design", "data.frame"))
}

theta_local_smc_log_marginal <- function(theta,
                                         data_i,
                                         loglik_fn,
                                         population_model,
                                         M = 600L,
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
  theta <- .theta_corr_as_matrix(theta, model)
  if (nrow(theta) != 1L) stop("theta_local_smc_log_marginal requires one theta row.")
  components <- population_model_reference_components_from_theta(model, theta)
  reference_prior <- make_reference_prior_gaussian(
    mu = components$component_means[[1L]],
    Sigma = components$component_covs[[1L]],
    param_names = model$alpha_names,
    label = "theta_correction_endpoint"
  )

  fit <- run_tempered_smc(
    reference_prior = reference_prior,
    bridge_stat_fn = function(alpha) ll_parallel(
      alpha,
      data_i,
      loglik_fn,
      n_cores = as.integer(n_cores)
    ),
    M = as.integer(M),
    resample_threshold = resample_threshold,
    n_mcmc_moves = as.integer(n_mcmc_moves),
    max_rounds = as.integer(max_steps),
    lambda_target = 1,
    cess_target = as.numeric(target_cess),
    G_mix = as.integer(G_mix),
    refit_every = as.integer(refit_every),
    rw_scale_init = as.numeric(rw_scale),
    post_adapt_n_mcmc_moves = as.integer(n_mcmc_moves),
    da_enable = isTRUE(da_enable),
    deterministic_resampling = deterministic_resampling,
    n_cores = as.integer(n_cores),
    seed = seed,
    verbose = verbose
  )

  if (fit$final_lambda < 1 - 1e-8) {
    stop("Local correction SMC did not reach lambda_target.")
  }

  list(
    log_marginal = as.numeric(fit$log_evidence),
    mcse_log_marginal = as.numeric(fit$mcse_logZ %||% NA_real_),
    final_lambda = as.numeric(fit$final_lambda),
    rounds = as.integer(fit$meta$rounds %||% NA_integer_),
    final_ess = as.numeric(ESS(fit$w)),
    final_ess_frac = as.numeric(ESS(fit$w) / length(fit$w)),
    mean_accept_rate = {
      acc <- as.numeric(fit$meta$accept_rate_hist %||% NA_real_)
      if (!length(acc) || all(is.na(acc))) NA_real_ else mean(acc, na.rm = TRUE)
    }
  )
}

.theta_correction_replicate_ids <- function(design) {
  ids <- design$point_id[design$role == "center"]
  stress <- design$point_id[design$role == "stress"]
  if (length(stress)) ids <- c(ids, stress[1L])
  sort(unique(ids))
}

fit_theta_smc_corrections <- function(design,
                                      factor_set,
                                      data_list,
                                      loglik_fn,
                                      M = 600L,
                                      n_jobs = 1L,
                                      local_n_cores = 1L,
                                      base_seed = 123L,
                                      replicate_ids = NULL,
                                      replicate_count = 2L,
                                      local_control = list(),
                                      verbose = TRUE) {
  stopifnot(inherits(factor_set, "population_factor_set"))
  model <- normalize_population_model(factor_set$population_model)
  design <- as.data.frame(design)
  theta <- .theta_corr_as_matrix(design[, model$hyper_names, drop = FALSE], model)
  n_points <- nrow(theta)
  n_local <- length(data_list)
  if (n_points < 1L) stop("design must contain at least one theta point.")

  sketch_by_subject <- population_factor_set_loglik_by_local(
    factor_set = factor_set,
    theta = theta,
    include_constant = FALSE,
    n_cores = as.integer(n_jobs)
  )
  log_prior <- population_model_log_hyperprior(model, theta)
  if (ncol(sketch_by_subject) != n_local) {
    stop("factor_set local dimension does not match data_list length.")
  }

  replicate_ids <- sort(unique(as.integer(replicate_ids %||% .theta_correction_replicate_ids(design))))
  replicate_ids <- intersect(replicate_ids, design$point_id)
  replicate_count <- as.integer(max(1L, replicate_count))

  eval_plan <- do.call(rbind, lapply(seq_len(n_points), function(j) {
    n_rep <- if (design$point_id[j] %in% replicate_ids) replicate_count else 1L
    data.frame(point_row = j, point_id = design$point_id[j], replicate_id = seq_len(n_rep))
  }))

  point_results <- vector("list", nrow(eval_plan))
  by_subject_rows <- vector("list", nrow(eval_plan))
  for (r in seq_len(nrow(eval_plan))) {
    point_row <- eval_plan$point_row[r]
    point_id <- eval_plan$point_id[r]
    replicate_id <- eval_plan$replicate_id[r]
    if (isTRUE(verbose)) {
      cat(sprintf(
        "Correction point %d/%d | %s | replicate %d\n",
        point_row,
        n_points,
        design$label[point_row],
        replicate_id
      ))
    }
    theta_row <- theta[point_row, , drop = FALSE]
    local_out <- parallel::mclapply(
      seq_len(n_local),
      function(i) {
        args <- modifyList(
          list(
            theta = theta_row,
            data_i = data_list[[i]],
            loglik_fn = loglik_fn,
            population_model = model,
            M = as.integer(M),
            n_cores = as.integer(local_n_cores),
            seed = as.integer(base_seed + 100000L * point_id + 1000L * replicate_id + i),
            verbose = FALSE
          ),
          local_control
        )
        fit <- do.call(theta_local_smc_log_marginal, args)
        fit$local_id <- i
        fit
      },
      mc.cores = as.integer(max(1L, n_jobs))
    )

    log_m_smc <- vapply(local_out, `[[`, numeric(1), "log_marginal")
    mcse <- vapply(local_out, `[[`, numeric(1), "mcse_log_marginal")
    log_m_tilde <- as.numeric(sketch_by_subject[point_row, ])
    delta <- log_m_smc - log_m_tilde
    by_subject <- data.frame(
      point_id = point_id,
      replicate_id = replicate_id,
      local_id = seq_len(n_local),
      log_m_smc = log_m_smc,
      log_m_tilde = log_m_tilde,
      delta = delta,
      mcse_log_m_smc = mcse,
      final_ess_frac = vapply(local_out, `[[`, numeric(1), "final_ess_frac"),
      rounds = vapply(local_out, `[[`, integer(1), "rounds"),
      mean_accept_rate = vapply(local_out, `[[`, numeric(1), "mean_accept_rate"),
      check.names = FALSE
    )
    by_subject_rows[[r]] <- by_subject
    point_results[[r]] <- data.frame(
      point_id = point_id,
      replicate_id = replicate_id,
      source = design$source[point_row],
      role = design$role[point_row],
      label = design$label[point_row],
      matrix(theta_row[1L, ], nrow = 1L, dimnames = list(NULL, model$hyper_names)),
      log_m_smc_total = sum(log_m_smc),
      log_m_tilde_total = sum(log_m_tilde),
      log_prior = log_prior[point_row],
      logposterior_smc = log_prior[point_row] + sum(log_m_smc),
      logposterior_tilde = log_prior[point_row] + sum(log_m_tilde),
      delta_total = sum(delta),
      mcse_total = sqrt(sum(mcse^2, na.rm = TRUE)),
      max_subject_abs_delta = max(abs(delta), na.rm = TRUE),
      max_subject_mcse = max(mcse, na.rm = TRUE),
      min_subject_log_m_smc = min(log_m_smc, na.rm = TRUE),
      degenerate_subjects = sum(by_subject$final_ess_frac >= 0.999 & by_subject$mcse_log_m_smc <= 0, na.rm = TRUE),
      min_final_ess_frac = min(by_subject$final_ess_frac, na.rm = TRUE),
      median_final_ess_frac = stats::median(by_subject$final_ess_frac, na.rm = TRUE),
      check.names = FALSE
    )
  }

  observations <- do.call(rbind, point_results)
  finite_smc <- is.finite(observations$logposterior_smc)
  finite_tilde <- is.finite(observations$logposterior_tilde)
  observations$relative_logposterior_smc <- NA_real_
  observations$relative_logposterior_tilde <- NA_real_
  if (any(finite_smc)) {
    observations$relative_logposterior_smc[finite_smc] <-
      observations$logposterior_smc[finite_smc] - max(observations$logposterior_smc[finite_smc])
  }
  if (any(finite_tilde)) {
    observations$relative_logposterior_tilde[finite_tilde] <-
      observations$logposterior_tilde[finite_tilde] - max(observations$logposterior_tilde[finite_tilde])
  }
  by_subject <- do.call(rbind, by_subject_rows)
  replicate_summary <- do.call(rbind, lapply(split(observations, observations$point_id), function(x) {
    data.frame(
      point_id = x$point_id[1L],
      n_replicates = nrow(x),
      delta_mean = mean(x$delta_total),
      delta_sd = if (nrow(x) > 1L) stats::sd(x$delta_total) else NA_real_,
      mcse_mean = mean(x$mcse_total),
      label = x$label[1L],
      role = x$role[1L],
      source = x$source[1L],
      check.names = FALSE
    )
  }))
  rownames(observations) <- NULL
  rownames(by_subject) <- NULL
  rownames(replicate_summary) <- NULL

  structure(
    list(
      design = design,
      observations = observations,
      by_subject = by_subject,
      replicate_summary = replicate_summary,
      factor_set = factor_set,
      settings = list(
        M = as.integer(M),
        n_jobs = as.integer(n_jobs),
        local_n_cores = as.integer(local_n_cores),
        base_seed = as.integer(base_seed),
        replicate_ids = replicate_ids,
        replicate_count = replicate_count,
        local_control = local_control
      )
    ),
    class = "theta_smc_correction_set"
  )
}

combine_theta_smc_correction_sets <- function(...) {
  sets <- list(...)
  sets <- sets[!vapply(sets, is.null, logical(1))]
  if (!length(sets)) stop("At least one correction set is required.")
  if (!all(vapply(sets, inherits, logical(1), "theta_smc_correction_set"))) {
    stop("All inputs must inherit from 'theta_smc_correction_set'.")
  }
  rbind_fill <- function(xs) {
    cols <- unique(unlist(lapply(xs, names), use.names = FALSE))
    xs <- lapply(xs, function(x) {
      missing <- setdiff(cols, names(x))
      for (nm in missing) x[[nm]] <- NA
      x[, cols, drop = FALSE]
    })
    do.call(rbind, xs)
  }
  factor_set <- sets[[1L]]$factor_set
  design <- rbind_fill(lapply(sets, function(x) as.data.frame(x$design)))
  observations <- rbind_fill(lapply(sets, `[[`, "observations"))
  by_subject <- rbind_fill(lapply(sets, `[[`, "by_subject"))
  replicate_summary <- rbind_fill(lapply(sets, `[[`, "replicate_summary"))
  if ("logposterior_smc" %in% names(observations)) {
    finite_smc <- is.finite(observations$logposterior_smc)
    observations$relative_logposterior_smc <- NA_real_
    if (any(finite_smc)) {
      observations$relative_logposterior_smc[finite_smc] <-
        observations$logposterior_smc[finite_smc] -
        max(observations$logposterior_smc[finite_smc])
    }
  }
  if ("logposterior_tilde" %in% names(observations)) {
    finite_tilde <- is.finite(observations$logposterior_tilde)
    observations$relative_logposterior_tilde <- NA_real_
    if (any(finite_tilde)) {
      observations$relative_logposterior_tilde[finite_tilde] <-
        observations$logposterior_tilde[finite_tilde] -
        max(observations$logposterior_tilde[finite_tilde])
    }
  }
  rownames(design) <- NULL
  rownames(observations) <- NULL
  rownames(by_subject) <- NULL
  rownames(replicate_summary) <- NULL
  structure(
    list(
      design = structure(design, class = c("theta_correction_design", "data.frame")),
      observations = observations,
      by_subject = by_subject,
      replicate_summary = replicate_summary,
      factor_set = factor_set,
      settings = list(
        combined_sets = length(sets),
        components = lapply(sets, `[[`, "settings")
      )
    ),
    class = "theta_smc_correction_set"
  )
}
