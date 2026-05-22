rm(list = ls())

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[1L]))
} else {
  normalizePath("benchmarks/run_shifted_gamma_support_diagnostics.R")
}
repo_dir <- dirname(dirname(script_path))
setwd(repo_dir)

parse_cli_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (!startsWith(arg, "--")) next
    arg <- sub("^--", "", arg)
    parts <- strsplit(arg, "=", fixed = TRUE)[[1L]]
    key <- gsub("-", "_", parts[1L])
    out[[key]] <- if (length(parts) > 1L) paste(parts[-1L], collapse = "=") else "true"
  }
  out
}

`%||%` <- function(x, y) if (is.null(x)) y else x
arg_chr <- function(args, key, default) {
  val <- args[[key]]
  if (is.null(val) || !nzchar(val)) default else as.character(val)
}
arg_int <- function(args, key, default) {
  val <- args[[key]]
  if (is.null(val) || !nzchar(val)) return(as.integer(default))
  as.integer(val)
}
arg_num <- function(args, key, default) {
  val <- args[[key]]
  if (is.null(val) || !nzchar(val)) return(as.numeric(default))
  as.numeric(val)
}
arg_lgl <- function(args, key, default = FALSE) {
  val <- args[[key]]
  if (is.null(val) || !nzchar(val)) return(isTRUE(default))
  tolower(as.character(val)) %in% c("1", "true", "t", "yes", "y")
}

if (!file.exists("smc_core.R")) {
  stop("Run this script from the FastHierarchical repository root.")
}

suppressPackageStartupMessages({
  library(parallel)
  library(qrng)
})

source("smc_core.R")
source("reference_priors.R")
source("SMC_super_fast.R")
source("population_models.R")
source("local_likelihood_sketches.R")
source("theta_proposals.R")
source("local_predictive_smc.R")

cli_args <- parse_cli_args(commandArgs(trailingOnly = TRUE))
detected_cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
if (!is.finite(detected_cores) || detected_cores < 1L) detected_cores <- 1L

stan_results_file <- arg_chr(
  cli_args,
  "stan_results_file",
  file.path("benchmarks", "samples", "shifted_gamma_hierarchy_stan_results.rds")
)
results_file <- arg_chr(
  cli_args,
  "results_file",
  file.path("benchmarks", "results", "shifted_gamma_support_diagnostics.rds")
)
summary_csv <- arg_chr(
  cli_args,
  "summary_csv",
  file.path("benchmarks", "results", "shifted_gamma_support_summary.csv")
)
q0_csv <- arg_chr(
  cli_args,
  "q0_csv",
  file.path("benchmarks", "results", "shifted_gamma_q0_summary.csv")
)
theta_scout_csv <- arg_chr(
  cli_args,
  "theta_scout_csv",
  file.path("benchmarks", "results", "shifted_gamma_theta_scout_summary.csv")
)
fitted_q0_csv <- arg_chr(
  cli_args,
  "fitted_q0_csv",
  file.path("benchmarks", "results", "shifted_gamma_fitted_q0_summary.csv")
)
fitted_q0_importance_csv <- arg_chr(
  cli_args,
  "fitted_q0_importance_csv",
  file.path("benchmarks", "results", "shifted_gamma_fitted_q0_importance.csv")
)
predictive_audit_csv <- arg_chr(
  cli_args,
  "predictive_audit_csv",
  file.path("benchmarks", "results", "shifted_gamma_predictive_audit.csv")
)
predictive_repair_csv <- arg_chr(
  cli_args,
  "predictive_repair_csv",
  file.path("benchmarks", "results", "shifted_gamma_predictive_repairs.csv")
)
plot_file <- arg_chr(
  cli_args,
  "plot_file",
  file.path("benchmarks", "results", "shifted_gamma_support_diagnostics.png")
)

seed <- arg_int(cli_args, "seed", 20260522L)
n_cores <- arg_int(cli_args, "cores", min(4L, detected_cores))
local_limit_arg <- cli_args[["local_limit"]]
local_limit <- if (is.null(local_limit_arg) || !nzchar(local_limit_arg)) Inf else as.integer(local_limit_arg)
theta_points <- arg_int(cli_args, "theta_points", 11L)
smc_particles <- arg_int(cli_args, "smc_particles", 400L)
laplace_draws <- arg_int(cli_args, "laplace_draws", 800L)
qmc_size <- arg_int(cli_args, "qmc_size", 4096L)
qmc_randomizations <- arg_int(cli_args, "qmc_randomizations", 2L)
q0_configs <- arg_int(cli_args, "q0_configs", 1000L)
theta_candidate_configs <- arg_int(cli_args, "theta_candidate_configs", q0_configs)
predictive_enabled <- arg_lgl(cli_args, "predictive_enabled", TRUE)
predictive_q0_method <- arg_chr(cli_args, "predictive_q0_method", "dmis_broad_laplace")
predictive_particles <- arg_int(cli_args, "predictive_particles", smc_particles)
predictive_components <- arg_int(cli_args, "predictive_components", 64L)
predictive_bridge <- arg_lgl(cli_args, "predictive_bridge", FALSE)
predictive_repair_enabled <- arg_lgl(cli_args, "predictive_repair_enabled", FALSE)
predictive_audit_points <- arg_int(cli_args, "predictive_audit_points", 24L)
predictive_audit_ess <- arg_num(cli_args, "predictive_audit_ess", 0.05)
predictive_max_repairs <- arg_int(cli_args, "predictive_max_repairs", 20L)
predictive_repair_particles <- arg_int(cli_args, "predictive_repair_particles", predictive_particles)
laplace_df <- arg_num(cli_args, "laplace_df", 5)
laplace_scale <- arg_num(cli_args, "laplace_scale", 6)
broad_scale <- arg_num(cli_args, "broad_scale", 16)
defensive_scale <- arg_num(cli_args, "defensive_scale", 64)
verbose <- arg_lgl(cli_args, "verbose", TRUE)

dir.create(dirname(results_file), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(summary_csv), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(q0_csv), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(theta_scout_csv), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(fitted_q0_csv), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(fitted_q0_importance_csv), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(predictive_audit_csv), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(predictive_repair_csv), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(plot_file), showWarnings = FALSE, recursive = TRUE)

if (!file.exists(stan_results_file)) {
  stop("Missing Stan benchmark results: ", stan_results_file)
}

set.seed(seed)
bundle <- readRDS(stan_results_file)
y <- bundle$data$y
data_list <- lapply(seq_len(nrow(y)), function(i) y[i, ])
if (is.finite(local_limit)) {
  data_list <- data_list[seq_len(min(length(data_list), as.integer(local_limit)))]
}

alpha_names <- c("eta_shape", "eta_scale", "eta_shift")
m0 <- stats::setNames(as.numeric(bundle$priors$m0), alpha_names)
s0 <- stats::setNames(as.numeric(bundle$priors$s0), alpha_names)
a0 <- stats::setNames(as.numeric(bundle$priors$a0), alpha_names)
b0 <- stats::setNames(as.numeric(bundle$priors$b0), alpha_names)

loglik_shifted_gamma <- function(Theta, y_i) {
  Theta <- as.matrix(Theta)
  colnames(Theta) <- alpha_names
  eps <- 1e-9
  shape <- exp(Theta[, "eta_shape"]) + eps
  scale <- exp(Theta[, "eta_scale"]) + eps
  shift <- exp(Theta[, "eta_shift"]) + eps
  min_y <- min(y_i)
  out <- rep(-1e12, nrow(Theta))
  ok <- shift < min_y
  if (!any(ok)) return(out)
  for (i in which(ok)) {
    out[i] <- sum(stats::dgamma(y_i - shift[i], shape = shape[i], scale = scale[i], log = TRUE))
  }
  out[!is.finite(out)] <- -1e12
  out
}

population_model <- make_population_model_diag_gaussian(
  alpha_names = alpha_names,
  mean_prior_mean = m0,
  mean_prior_var = s0,
  sigma2_prior_shape = a0,
  sigma2_prior_rate = b0,
  label = "shifted_gamma_hierarchy"
)

initial_theta <- matrix(c(m0, log(b0 / (a0 - 1))), nrow = 1L)
colnames(initial_theta) <- population_model$hyper_names
initial_ref <- make_reference_prior_gaussian(
  mu = m0,
  Sigma = diag(b0 / (a0 - 1), length(alpha_names), length(alpha_names)),
  param_names = alpha_names,
  label = "initial_population_reference"
)
broad_ref <- make_broad_reference_prior(
  mu = m0,
  Sigma = diag(b0 / (a0 - 1), length(alpha_names), length(alpha_names)),
  scale = broad_scale,
  defensive = TRUE,
  defensive_scale = defensive_scale,
  defensive_weight = 0.20,
  param_names = alpha_names,
  label = "broad_likelihood_scout_reference"
)

stan_theta <- cbind(bundle$draws$mu, log(bundle$draws$sigma2))
colnames(stan_theta) <- population_model$hyper_names
stan_theta <- as.matrix(stan_theta)

select_theta_stress_points <- function(theta, n_points) {
  theta <- as.matrix(theta)
  n_points <- as.integer(max(1L, min(n_points, nrow(theta))))
  center <- apply(theta, 2L, stats::median)
  S <- stats::cov(theta)
  S <- regularize_cov(S, min_eig = 1e-10, cond_cap = 1e8)
  L <- chol(S)
  Z <- t(backsolve(L, t(sweep(theta, 2L, center, "-")), transpose = TRUE))
  r2 <- rowSums(Z * Z)
  pc <- tryCatch({
    eig <- eigen(S, symmetric = TRUE)
    as.numeric(scale(theta, center = center, scale = FALSE) %*% eig$vectors[, 1L])
  }, error = function(e) theta[, 1L])

  idx <- integer(0)
  idx <- c(idx, which.min(r2))
  probs <- c(0.02, 0.05, 0.20, 0.50, 0.80, 0.95, 0.98)
  vals <- stats::quantile(pc, probs = probs, names = FALSE, type = 8)
  idx <- c(idx, vapply(vals, function(v) which.min(abs(pc - v)), integer(1)))
  tail_idx <- order(r2, decreasing = TRUE)
  idx <- unique(c(idx, tail_idx[seq_len(min(length(tail_idx), n_points))]))
  idx <- idx[seq_len(min(length(idx), n_points))]
  out <- theta[idx, , drop = FALSE]
  rownames(out) <- paste0("theta_", seq_len(nrow(out)))
  attr(out, "source_rows") <- idx
  out
}

theta_stress <- select_theta_stress_points(stan_theta, theta_points)

fit_local_smc <- function(data_i, reference_prior, M, seed_i) {
  run_tempered_smc(
    reference_prior = reference_prior,
    bridge_stat_fn = function(alpha) ll_parallel(alpha, data_i, loglik_shifted_gamma, n_cores = 1L),
    M = as.integer(M),
    resample_threshold = 0.6,
    n_mcmc_moves = 2L,
    post_adapt_n_mcmc_moves = 2L,
    max_rounds = 100L,
    cess_target = 0.9,
    G_mix = 8L,
    rw_scale_init = 0.9,
    n_cores = 1L,
    seed = as.integer(seed_i),
    verbose = FALSE
  )
}

rmvt_chol <- function(n, mu, L, df) {
  z <- matrix(stats::rnorm(n * length(mu)), nrow = n)
  s <- sqrt(stats::rchisq(n, df = df) / df)
  out <- sweep((z %*% L) / s, 2L, mu, "+")
  colnames(out) <- names(mu)
  out
}

fit_local_support_objects <- function(local_id) {
  data_i <- data_list[[local_id]]
  init_fit <- fit_local_smc(
    data_i,
    reference_prior = initial_ref,
    M = smc_particles,
    seed_i = seed + 1000L + local_id
  )
  broad_fit <- fit_local_smc(
    data_i,
    reference_prior = broad_ref,
    M = smc_particles,
    seed_i = seed + 2000L + local_id
  )
  sketch <- fit_local_likelihood_sketch(
    data_i = data_i,
    loglik_fn = loglik_shifted_gamma,
    alpha_names = alpha_names,
    local_id = local_id,
    population_model = population_model,
    theta_reference = initial_theta,
    n_starts = 12L,
    start_scale = 9,
    reference_scale = 50,
    seed = seed + 3000L + local_id
  )
  lap_mu <- as.numeric(sketch$component_means[1L, ])
  names(lap_mu) <- alpha_names
  lap_cov <- regularize_cov(laplace_scale * sketch$component_covs[[1L]], min_eig = 1e-8, cond_cap = 1e8)
  dimnames(lap_cov) <- list(alpha_names, alpha_names)
  lap_L <- chol(lap_cov)
  lap_alpha <- rmvt_chol(laplace_draws, lap_mu, lap_L, df = laplace_df)
  lap_loglik <- ll_parallel(lap_alpha, data_i, loglik_shifted_gamma, n_cores = 1L)
  lap_logq <- dmvt_chol_log(lap_alpha, lap_mu, lap_L, df = laplace_df)

  list(
    local_id = local_id,
    data_i = data_i,
    initial = list(
      alpha = as.matrix(init_fit$Theta),
      w = normalize_weights(init_fit$w),
      loglik = as.numeric(init_fit$loglik),
      reference_prior = initial_ref,
      log_ref = reference_prior_logpdf(initial_ref, init_fit$Theta),
      logZ = as.numeric(init_fit$log_evidence),
      mcse_logZ = as.numeric(init_fit$mcse_logZ)
    ),
    broad = list(
      alpha = as.matrix(broad_fit$Theta),
      w = normalize_weights(broad_fit$w),
      loglik = as.numeric(broad_fit$loglik),
      reference_prior = broad_ref,
      log_ref = reference_prior_logpdf(broad_ref, broad_fit$Theta),
      logZ = as.numeric(broad_fit$log_evidence),
      mcse_logZ = as.numeric(broad_fit$mcse_logZ)
    ),
    laplace = list(
      alpha = lap_alpha,
      loglik = lap_loglik,
      logq = lap_logq,
      mean = lap_mu,
      cov = lap_cov,
      chol = lap_L,
      df = laplace_df,
      sketch = sketch
    )
  )
}

normalize_weights <- function(w) {
  w <- pmax(as.numeric(w), 0)
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) rep(1 / length(w), length(w)) else w / sw
}

log_ess_frac <- function(log_terms) {
  log_terms <- as.numeric(log_terms)
  ok <- is.finite(log_terms)
  if (!any(ok)) return(0)
  lw <- log_terms[ok]
  lse <- logsumexp(lw)
  w <- exp(lw - lse)
  as.numeric(1 / sum(w * w) / length(log_terms))
}

single_posterior_logm <- function(obj, theta) {
  logp_theta <- population_model_log_alpha_given_theta(population_model, obj$alpha, theta)
  log_terms <- log(obj$w) + logp_theta - obj$log_ref
  list(
    logm = obj$logZ + logsumexp(log_terms),
    ess_frac = log_ess_frac(log_terms)
  )
}

laplace_logm <- function(local, theta) {
  logp_theta <- population_model_log_alpha_given_theta(population_model, local$laplace$alpha, theta)
  log_terms <- local$laplace$loglik + logp_theta - local$laplace$logq - log(nrow(local$laplace$alpha))
  list(
    logm = logsumexp(log_terms),
    ess_frac = log_ess_frac(log_terms)
  )
}

dmis_broad_laplace_logm <- function(local, theta, eta_broad = 0.5) {
  broad <- local$broad
  lap <- local$laplace
  eta_lap <- 1 - eta_broad

  alpha <- rbind(broad$alpha, lap$alpha)
  source <- c(rep.int("broad", nrow(broad$alpha)), rep.int("laplace", nrow(lap$alpha)))
  loglik <- c(broad$loglik, lap$loglik)
  logp_theta <- population_model_log_alpha_given_theta(population_model, alpha, theta)
  logq_broad <- loglik + reference_prior_logpdf(broad$reference_prior, alpha) - broad$logZ
  logq_lap <- dmvt_chol_log(alpha, lap$mean, lap$chol, df = lap$df)
  logq_mix <- rlogsumexp2(log(eta_broad) + logq_broad, log(eta_lap) + logq_lap)

  sample_logw <- numeric(nrow(alpha))
  sample_logw[source == "broad"] <- log(eta_broad) + log(broad$w)
  sample_logw[source == "laplace"] <- log(eta_lap) - log(nrow(lap$alpha))
  log_terms <- sample_logw + loglik + logp_theta - logq_mix
  list(
    logm = logsumexp(log_terms),
    ess_frac = log_ess_frac(log_terms)
  )
}

posterior_obj_from_component <- function(component) {
  alpha <- as.matrix(component$particles)
  list(
    alpha = alpha,
    w = normalize_weights(component$weights),
    loglik = as.numeric(component$log_likelihood),
    reference_prior = component$reference_prior,
    log_ref = as.numeric(component$log_reference_density %||% reference_prior_logpdf(component$reference_prior, alpha)),
    logZ = as.numeric(component$log_evidence),
    mcse_logZ = as.numeric(component$mcse_log_evidence %||% NA_real_),
    component = component
  )
}

component_from_posterior_obj <- function(obj, local_id, source) {
  structure(
    list(
      local_id = as.integer(local_id),
      particles = as.matrix(obj$alpha),
      weights = normalize_weights(obj$w),
      reference_prior = obj$reference_prior,
      log_reference_density = as.numeric(obj$log_ref),
      proposal_type = "posterior_reference",
      log_likelihood = as.numeric(obj$loglik),
      log_evidence = as.numeric(obj$logZ),
      mcse_log_evidence = as.numeric(obj$mcse_logZ %||% NA_real_),
      diagnostics = list(source = source)
    ),
    class = "reference_local_component"
  )
}

local_dmis_stack <- function(local, component_names, eta = NULL) {
  eta <- normalize_weights(eta %||% rep(1, length(component_names)))
  parts <- lapply(seq_along(component_names), function(k) {
    name <- component_names[k]
    if (identical(name, "initial")) {
      obj <- local$initial
      return(list(
        name = name,
        kind = "posterior",
        alpha = obj$alpha,
        loglik = obj$loglik,
        sample_logw = log(eta[k]) + log(obj$w),
        reference_prior = obj$reference_prior,
        logZ = obj$logZ
      ))
    }
    if (identical(name, "broad")) {
      obj <- local$broad
      return(list(
        name = name,
        kind = "posterior",
        alpha = obj$alpha,
        loglik = obj$loglik,
        sample_logw = log(eta[k]) + log(obj$w),
        reference_prior = obj$reference_prior,
        logZ = obj$logZ
      ))
    }
    if (identical(name, "predictive")) {
      obj <- local$predictive
      return(list(
        name = name,
        kind = "posterior",
        alpha = obj$alpha,
        loglik = obj$loglik,
        sample_logw = log(eta[k]) + log(obj$w),
        reference_prior = obj$reference_prior,
        logZ = obj$logZ
      ))
    }
    if (identical(name, "laplace")) {
      obj <- local$laplace
      return(list(
        name = name,
        kind = "direct",
        alpha = obj$alpha,
        loglik = obj$loglik,
        sample_logw = rep(log(eta[k]) - log(nrow(obj$alpha)), nrow(obj$alpha)),
        mean = obj$mean,
        chol = obj$chol,
        df = obj$df
      ))
    }
    stop("Unknown local DMIS component: ", name)
  })

  alpha <- do.call(rbind, lapply(parts, `[[`, "alpha"))
  colnames(alpha) <- alpha_names
  loglik <- unlist(lapply(parts, `[[`, "loglik"), use.names = FALSE)
  sample_logw <- unlist(lapply(parts, `[[`, "sample_logw"), use.names = FALSE)
  logq <- matrix(-Inf, nrow = nrow(alpha), ncol = length(parts))
  for (k in seq_along(parts)) {
    part <- parts[[k]]
    if (identical(part$kind, "posterior")) {
      logq[, k] <- loglik + reference_prior_logpdf(part$reference_prior, alpha) - part$logZ
    } else {
      logq[, k] <- dmvt_chol_log(alpha, part$mean, part$chol, df = part$df)
    }
  }
  logq_mix <- .rowLogSumExp(sweep(logq, 2L, log(eta), "+"))
  list(
    alpha = alpha,
    loglik = loglik,
    sample_logw = sample_logw,
    logq_mix = as.numeric(logq_mix),
    component_names = component_names,
    eta = eta
  )
}

dmis_local_logm <- function(local, theta, component_names, eta = NULL) {
  stack <- local_dmis_stack(local, component_names = component_names, eta = eta)
  logp_theta <- population_model_log_alpha_given_theta(population_model, stack$alpha, theta)
  log_terms <- stack$sample_logw + stack$loglik + logp_theta - stack$logq_mix
  list(
    logm = logsumexp(log_terms),
    ess_frac = log_ess_frac(log_terms)
  )
}

qmc_local_logm <- function(data_i, theta, qmc_n, R, seed_i) {
  theta_prepared <- population_model_prepare_theta(population_model, theta)
  sigma <- sqrt(1 / theta_prepared$quadratic_coef[1L, ])
  mu <- as.numeric(theta_prepared$mean[1L, ])
  names(mu) <- alpha_names
  reps <- vapply(seq_len(R), function(r) {
    set.seed(as.integer(seed_i + r))
    U <- qrng::sobol(n = qmc_n, d = length(alpha_names), randomize = TRUE)
    alpha <- sweep(qnorm(pmin(pmax(U, 1e-12), 1 - 1e-12)), 2L, sigma, "*")
    alpha <- sweep(alpha, 2L, mu, "+")
    colnames(alpha) <- alpha_names
    ll <- ll_parallel(alpha, data_i, loglik_shifted_gamma, n_cores = 1L)
    logsumexp(ll) - log(length(ll))
  }, numeric(1))
  list(
    logm = logsumexp(reps) - log(length(reps)),
    mcse = if (length(reps) > 1L) stats::sd(reps) / sqrt(length(reps)) else NA_real_
  )
}

evaluate_local_support <- function(local) {
  rows <- vector("list", nrow(theta_stress))
  for (k in seq_len(nrow(theta_stress))) {
    theta <- theta_stress[k, , drop = FALSE]
    ref <- qmc_local_logm(
      data_i = local$data_i,
      theta = theta,
      qmc_n = qmc_size,
      R = qmc_randomizations,
      seed_i = seed + 50000L + 1000L * local$local_id + k
    )
    estimates <- list(
      initial_anchor = single_posterior_logm(local$initial, theta),
      broad_scout = single_posterior_logm(local$broad, theta),
      laplace_t = laplace_logm(local, theta),
      dmis_broad_laplace = dmis_broad_laplace_logm(local, theta)
    )
    if (!is.null(local$predictive)) {
      estimates$q0_predictive <- single_posterior_logm(local$predictive, theta)
      estimates$dmis_broad_predictive_laplace <- dmis_local_logm(
        local,
        theta,
        component_names = c("broad", "predictive", "laplace"),
        eta = c(0.35, 0.45, 0.20)
      )
    }
    if (!is.null(local$predictive_repaired_factor)) {
      estimates$q0_predictive_bank_dmis <- list(
        logm = as.numeric(population_local_factor_log_marginal(local$predictive_repaired_factor, theta)),
        ess_frac = as.numeric(population_local_factor_ess(local$predictive_repaired_factor, theta) /
                                local$predictive_repaired_factor$n_particles)
      )
    }
    rows[[k]] <- do.call(rbind, lapply(names(estimates), function(method) {
      est <- estimates[[method]]
      data.frame(
        local_id = as.integer(local$local_id),
        theta_id = as.integer(k),
        method = method,
        reference_logm = as.numeric(ref$logm),
        reference_mcse = as.numeric(ref$mcse),
        estimate_logm = as.numeric(est$logm),
        error = as.numeric(est$logm - ref$logm),
        abs_error = abs(as.numeric(est$logm - ref$logm)),
        ess_frac = as.numeric(est$ess_frac),
        check.names = FALSE
      )
    }))
  }
  do.call(rbind, rows)
}

sample_local_alpha <- function(local, method, n) {
  n <- as.integer(n)
  if (identical(method, "initial_anchor")) {
    idx <- sample.int(nrow(local$initial$alpha), n, replace = TRUE, prob = local$initial$w)
    return(local$initial$alpha[idx, , drop = FALSE])
  }
  if (identical(method, "broad_scout")) {
    idx <- sample.int(nrow(local$broad$alpha), n, replace = TRUE, prob = local$broad$w)
    return(local$broad$alpha[idx, , drop = FALSE])
  }
  if (identical(method, "laplace_t")) {
    return(rmvt_chol(n, local$laplace$mean, local$laplace$chol, df = local$laplace$df))
  }
  if (identical(method, "dmis_broad_laplace")) {
    take_broad <- stats::runif(n) < 0.5
    out <- matrix(NA_real_, nrow = n, ncol = length(alpha_names), dimnames = list(NULL, alpha_names))
    if (any(take_broad)) {
      idx <- sample.int(nrow(local$broad$alpha), sum(take_broad), replace = TRUE, prob = local$broad$w)
      out[take_broad, ] <- local$broad$alpha[idx, , drop = FALSE]
    }
    if (any(!take_broad)) {
      out[!take_broad, ] <- rmvt_chol(sum(!take_broad), local$laplace$mean, local$laplace$chol, df = local$laplace$df)
    }
    return(out)
  }
  stop("Unknown method: ", method)
}

inv_gamma_quantile <- function(p, shape, rate) {
  1 / stats::qgamma(1 - p, shape = shape, rate = rate)
}

log_ell_marginal <- function(ell, y_alpha, m0_j, s0_j, a_j, b_j) {
  inv_sigma2 <- exp(-ell)
  n <- length(y_alpha)
  sum_y <- sum(y_alpha)
  sum_y2 <- sum(y_alpha * y_alpha)
  inv_s0 <- 1 / s0_j
  precision <- inv_s0 + n * inv_sigma2
  h <- m0_j * inv_s0 + sum_y * inv_sigma2
  cterm <- m0_j^2 * inv_s0 + sum_y2 * inv_sigma2
  log_ig_jac <- a_j * log(b_j) - lgamma(a_j) - a_j * ell - b_j * inv_sigma2
  log_norm <- -0.5 * (log(2 * pi * s0_j) + n * log(2 * pi) + n * ell)
  log_ig_jac + log_norm - 0.5 * cterm + 0.5 * h * h / precision +
    0.5 * (log(2 * pi) - log(precision))
}

sample_theta_given_alpha_config <- function(alpha_config, seed_i = NULL) {
  if (!is.null(seed_i)) set.seed(as.integer(seed_i))
  alpha_config <- as.matrix(alpha_config)
  theta <- numeric(population_model$hyper_dim)
  names(theta) <- population_model$hyper_names
  for (j in seq_along(alpha_names)) {
    yj <- alpha_config[, j]
    empirical_var <- if (length(yj) > 1L) stats::var(yj) else b0[j] / (a0[j] - 1)
    empirical_var <- pmax(empirical_var, 1e-8)
    lower <- min(log(inv_gamma_quantile(1e-4, a0[j], b0[j])), log(empirical_var) - 4)
    upper <- max(log(inv_gamma_quantile(1 - 1e-4, a0[j], b0[j])), log(empirical_var) + 4)
    ell_grid <- seq(lower, upper, length.out = 128L)
    logp <- log_ell_marginal(ell_grid, yj, m0[j], s0[j], a0[j], b0[j])
    logp[!is.finite(logp)] <- -Inf
    p <- exp(logp - logsumexp(logp))
    ell <- sample(ell_grid, 1L, prob = p)
    inv_sigma2 <- exp(-ell)
    precision <- 1 / s0[j] + length(yj) * inv_sigma2
    h <- m0[j] / s0[j] + sum(yj) * inv_sigma2
    mu <- stats::rnorm(1L, mean = h / precision, sd = sqrt(1 / precision))
    theta[j] <- mu
    theta[length(alpha_names) + j] <- ell
  }
  theta
}

sample_q0_theta <- function(locals, method, n_configs) {
  theta <- matrix(NA_real_, nrow = n_configs, ncol = population_model$hyper_dim)
  colnames(theta) <- population_model$hyper_names
  for (k in seq_len(n_configs)) {
    alpha_config <- do.call(rbind, lapply(locals, sample_local_alpha, method = method, n = 1L))
    colnames(alpha_config) <- alpha_names
    theta[k, ] <- sample_theta_given_alpha_config(alpha_config)
  }
  theta
}

quantile_distance <- function(x, y, probs = seq(0.01, 0.99, length.out = 99L)) {
  qx <- stats::quantile(x, probs = probs, names = FALSE, type = 8)
  qy <- stats::quantile(y, probs = probs, names = FALSE, type = 8)
  mean(abs(qx - qy))
}

summarize_q0 <- function(theta_draws, method) {
  rows <- lapply(seq_len(ncol(stan_theta)), function(j) {
    q_method <- stats::quantile(theta_draws[, j], probs = c(0.01, 0.05, 0.95, 0.99), names = FALSE, type = 8)
    q_stan <- stats::quantile(stan_theta[, j], probs = c(0.01, 0.05, 0.95, 0.99), names = FALSE, type = 8)
    data.frame(
      method = method,
      parameter = colnames(stan_theta)[j],
      stan_inside_q01_q99 = mean(stan_theta[, j] >= q_method[1L] & stan_theta[, j] <= q_method[4L]),
      stan_inside_q05_q95 = mean(stan_theta[, j] >= q_method[2L] & stan_theta[, j] <= q_method[3L]),
      lower_q01_miss = max(q_method[1L] - q_stan[1L], 0),
      upper_q99_miss = max(q_stan[4L] - q_method[4L], 0),
      width99_ratio = (q_method[4L] - q_method[1L]) / max(q_stan[4L] - q_stan[1L], .Machine$double.eps),
      q_wasserstein = quantile_distance(theta_draws[, j], stan_theta[, j]),
      check.names = FALSE
    )
  })
  do.call(rbind, rows)
}

posterior_logm_many <- function(obj, theta) {
  logp_theta <- population_model_log_alpha_given_theta_many(population_model, obj$alpha, theta)
  offset <- log(obj$w) - obj$log_ref
  as.numeric(obj$logZ + .rowLogSumExp(sweep(logp_theta, 2L, offset, "+")))
}

laplace_logm_many <- function(local, theta) {
  logp_theta <- population_model_log_alpha_given_theta_many(population_model, local$laplace$alpha, theta)
  offset <- local$laplace$loglik - local$laplace$logq - log(nrow(local$laplace$alpha))
  as.numeric(.rowLogSumExp(sweep(logp_theta, 2L, offset, "+")))
}

dmis_broad_laplace_logm_many <- function(local, theta, eta_broad = 0.5) {
  broad <- local$broad
  lap <- local$laplace
  eta_lap <- 1 - eta_broad

  alpha <- rbind(broad$alpha, lap$alpha)
  source <- c(rep.int("broad", nrow(broad$alpha)), rep.int("laplace", nrow(lap$alpha)))
  loglik <- c(broad$loglik, lap$loglik)
  logq_broad <- loglik + reference_prior_logpdf(broad$reference_prior, alpha) - broad$logZ
  logq_lap <- dmvt_chol_log(alpha, lap$mean, lap$chol, df = lap$df)
  logq_mix <- rlogsumexp2(log(eta_broad) + logq_broad, log(eta_lap) + logq_lap)

  sample_logw <- numeric(nrow(alpha))
  sample_logw[source == "broad"] <- log(eta_broad) + log(broad$w)
  sample_logw[source == "laplace"] <- log(eta_lap) - log(nrow(lap$alpha))
  offset <- sample_logw + loglik - logq_mix
  logp_theta <- population_model_log_alpha_given_theta_many(population_model, alpha, theta)
  as.numeric(.rowLogSumExp(sweep(logp_theta, 2L, offset, "+")))
}

dmis_local_logm_many <- function(local, theta, component_names, eta = NULL) {
  stack <- local_dmis_stack(local, component_names = component_names, eta = eta)
  logp_theta <- population_model_log_alpha_given_theta_many(population_model, stack$alpha, theta)
  offset <- stack$sample_logw + stack$loglik - stack$logq_mix
  as.numeric(.rowLogSumExp(sweep(logp_theta, 2L, offset, "+")))
}

local_logm_many <- function(local, theta, method) {
  if (identical(method, "initial_anchor")) return(posterior_logm_many(local$initial, theta))
  if (identical(method, "broad_scout")) return(posterior_logm_many(local$broad, theta))
  if (identical(method, "laplace_t")) return(laplace_logm_many(local, theta))
  if (identical(method, "dmis_broad_laplace")) return(dmis_broad_laplace_logm_many(local, theta))
  if (identical(method, "q0_predictive")) return(posterior_logm_many(local$predictive, theta))
  if (identical(method, "dmis_broad_predictive_laplace")) {
    return(dmis_local_logm_many(
      local,
      theta,
      component_names = c("broad", "predictive", "laplace"),
      eta = c(0.35, 0.45, 0.20)
    ))
  }
  stop("Unknown theta scoring method: ", method)
}

make_theta_candidate_pool <- function(q0_draws, n_extra, seed_i) {
  set.seed(as.integer(seed_i))
  raw <- do.call(rbind, q0_draws)
  raw_source <- rep(names(q0_draws), vapply(q0_draws, nrow, integer(1L)))
  colnames(raw) <- population_model$hyper_names

  prior <- population_model_sample_hyper(population_model, n_extra)
  recombined <- matrix(
    NA_real_,
    nrow = n_extra,
    ncol = ncol(raw),
    dimnames = list(NULL, colnames(raw))
  )
  for (j in seq_len(ncol(raw))) {
    recombined[, j] <- sample(raw[, j], n_extra, replace = TRUE)
  }

  mu_cols <- seq_len(population_model$alpha_dim)
  sigma_cols <- population_model$alpha_dim + seq_len(population_model$alpha_dim)
  support_methods <- intersect(c("broad_scout", "laplace_t", "dmis_broad_laplace"), names(q0_draws))
  exploratory <- do.call(rbind, q0_draws[support_methods])
  anchor <- q0_draws$initial_anchor %||% raw
  mu_splice <- matrix(
    NA_real_,
    nrow = n_extra,
    ncol = ncol(raw),
    dimnames = list(NULL, colnames(raw))
  )
  mu_idx <- sample.int(nrow(exploratory), n_extra, replace = TRUE)
  sigma_idx <- sample.int(nrow(anchor), n_extra, replace = TRUE)
  mu_splice[, mu_cols] <- exploratory[mu_idx, mu_cols, drop = FALSE]
  mu_splice[, sigma_cols] <- anchor[sigma_idx, sigma_cols, drop = FALSE]

  theta <- rbind(raw, prior, recombined, mu_splice)
  source <- c(
    raw_source,
    rep.int("hyperprior", nrow(prior)),
    rep.int("marginal_recombine", nrow(recombined)),
    rep.int("mu_support_variance_anchor", nrow(mu_splice))
  )
  ok <- apply(theta, 1L, function(x) all(is.finite(x)))
  theta <- theta[ok, , drop = FALSE]
  attr(theta, "source") <- source[ok]
  theta
}

normalize_log_weights <- function(logw) {
  logw <- as.numeric(logw)
  ok <- is.finite(logw)
  if (!any(ok)) return(rep(1 / length(logw), length(logw)))
  out <- rep(0, length(logw))
  out[ok] <- exp(logw[ok] - logsumexp(logw[ok]))
  out
}

weighted_quantile_no_warn <- function(x, w, probs) {
  x <- as.numeric(x)
  w <- pmax(as.numeric(w), 0)
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(stats::quantile(x, probs = probs, names = FALSE, type = 8))
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
  cw <- cumsum(w)
  stats::approx(cw, x, xout = probs, method = "constant", f = 1, rule = 2, ties = "ordered")$y
}

weighted_quantile_distance <- function(x, w, y, probs = seq(0.01, 0.99, length.out = 99L)) {
  qx <- weighted_quantile_no_warn(x, w, probs = probs)
  qy <- stats::quantile(y, probs = probs, names = FALSE, type = 8)
  mean(abs(qx - qy))
}

score_theta_candidate_pool <- function(theta_candidates, locals, method) {
  theta <- as.matrix(theta_candidates)
  log_score <- population_model_log_hyperprior(population_model, theta)
  for (local in locals) {
    log_score <- log_score + local_logm_many(local, theta, method)
  }
  w <- normalize_log_weights(log_score)
  source <- attr(theta_candidates, "source") %||% rep.int("candidate", nrow(theta))
  list(
    log_weight = log_score,
    ess_frac = as.numeric(1 / sum(w * w) / length(w)),
    max_weight = max(w),
    source_mass = tapply(w, source, sum)
  )
}

score_theta_matrix <- function(theta, locals, method) {
  theta <- as.matrix(theta)
  log_score <- population_model_log_hyperprior(population_model, theta)
  for (local in locals) {
    log_score <- log_score + local_logm_many(local, theta, method)
  }
  as.numeric(log_score)
}

summarize_weighted_theta_support <- function(theta_candidates, theta_score, method) {
  theta <- as.matrix(theta_candidates)
  w <- normalize_log_weights(theta_score$log_weight)
  rows <- lapply(seq_len(ncol(stan_theta)), function(j) {
    q_method <- weighted_quantile_no_warn(theta[, j], w, probs = c(0.01, 0.05, 0.95, 0.99))
    q_stan <- stats::quantile(stan_theta[, j], probs = c(0.01, 0.05, 0.95, 0.99), names = FALSE, type = 8)
    data.frame(
      method = paste0("theta_score_", method),
      parameter = colnames(stan_theta)[j],
      stan_inside_q01_q99 = mean(stan_theta[, j] >= q_method[1L] & stan_theta[, j] <= q_method[4L]),
      stan_inside_q05_q95 = mean(stan_theta[, j] >= q_method[2L] & stan_theta[, j] <= q_method[3L]),
      lower_q01_miss = max(q_method[1L] - q_stan[1L], 0),
      upper_q99_miss = max(q_stan[4L] - q_method[4L], 0),
      width99_ratio = (q_method[4L] - q_method[1L]) / max(q_stan[4L] - q_stan[1L], .Machine$double.eps),
      q_wasserstein = weighted_quantile_distance(theta[, j], w, stan_theta[, j]),
      theta_ess_frac = as.numeric(theta_score$ess_frac),
      theta_max_weight = as.numeric(theta_score$max_weight),
      n_candidates = nrow(theta),
      check.names = FALSE
    )
  })
  do.call(rbind, rows)
}

fit_scored_q0_proposal <- function(theta_candidates, theta_score, method, seed_i) {
  fit_theta_q0_proposal(
    theta = theta_candidates,
    log_weight = theta_score$log_weight,
    population_model = population_model,
    max_components = 6L,
    core_weight = 0.85,
    tail_weight = 0.10,
    prior_weight = 0.05,
    df = 7,
    tail_df = 3,
    core_scale = 1.25,
    tail_scale = 9,
    cluster_sample_size = 2500L,
    min_component_weight = 0.02,
    seed = seed_i,
    label = paste0("fitted_q0_", method)
  )
}

summarize_fitted_q0_importance <- function(proposal, locals, method, n_draws, seed_i) {
  theta <- theta_proposal_sample(proposal, n_draws, seed = seed_i)
  log_target <- score_theta_matrix(theta, locals = locals, method = method)
  logq <- theta_proposal_log_density(proposal, theta)
  logw <- log_target - logq
  ok <- is.finite(logw)
  if (!any(ok)) {
    return(data.frame(
      method = paste0("fitted_q0_", method),
      n = n_draws,
      finite_frac = 0,
      importance_ess_frac = 0,
      max_weight = NA_real_,
      log_evidence_estimate = -Inf,
      log_weight_q01 = NA_real_,
      log_weight_median = NA_real_,
      log_weight_q99 = NA_real_,
      check.names = FALSE
    ))
  }
  lw <- logw[ok]
  w <- exp(lw - logsumexp(lw))
  data.frame(
    method = paste0("fitted_q0_", method),
    n = n_draws,
    finite_frac = mean(ok),
    importance_ess_frac = as.numeric(1 / sum(w * w) / length(logw)),
    max_weight = max(w),
    log_evidence_estimate = logsumexp(lw) - log(length(logw)),
    log_weight_q01 = as.numeric(stats::quantile(lw, 0.01, names = FALSE, type = 8)),
    log_weight_median = median(lw),
    log_weight_q99 = as.numeric(stats::quantile(lw, 0.99, names = FALSE, type = 8)),
    check.names = FALSE
  )
}

add_predictive_component_to_local <- function(local) {
  q0 <- fitted_q0_proposals[[predictive_q0_method]]
  if (is.null(q0)) {
    stop("Unknown predictive_q0_method: ", predictive_q0_method)
  }
  predictive_ref <- make_population_predictive_reference_prior(
    population_model = population_model,
    theta_proposal = q0,
    n_components = predictive_components,
    defensive_weight = 0.05,
    defensive_scale = 16,
    seed = seed + 830000L + local$local_id,
    label = paste0("q0_predictive_", predictive_q0_method)
  )
  broad_component <- component_from_posterior_obj(
    local$broad,
    local_id = local$local_id,
    source = "broad_scout"
  )
  predictive_component <- if (isTRUE(predictive_bridge)) {
    bridge_local_reference_smc_component(
      component = broad_component,
      data_i = local$data_i,
      loglik_fn = loglik_shifted_gamma,
      new_reference_prior = predictive_ref,
      M = predictive_particles,
      resample_threshold = 0.6,
      n_mcmc_moves = 1L,
      max_rounds = 80L,
      cess_target = 0.9,
      G_mix = 8L,
      rw_scale_init = 0.7,
      n_cores = 1L,
      seed = seed + 840000L + local$local_id,
      diagnostics = list(
        q0_method = predictive_q0_method,
        predictive_components = as.integer(predictive_components)
      ),
      verbose = FALSE
    )
  } else {
    fit_local_reference_smc_component(
      local_id = local$local_id,
      data_i = local$data_i,
      loglik_fn = loglik_shifted_gamma,
      reference_prior = predictive_ref,
      M = predictive_particles,
      resample_threshold = 0.6,
      n_mcmc_moves = 2L,
      max_rounds = 100L,
      cess_target = 0.9,
      G_mix = 8L,
      rw_scale_init = 0.9,
      n_cores = 1L,
      seed = seed + 840000L + local$local_id,
      diagnostics = list(
        source = "q0_predictive_reference",
        q0_method = predictive_q0_method,
        predictive_components = as.integer(predictive_components)
      ),
      verbose = FALSE
    )
  }
  local$predictive <- posterior_obj_from_component(predictive_component)
  local$predictive_component <- predictive_component
  local
}

predictive_reference_object_from_local <- function(local) {
  broad_component <- component_from_posterior_obj(
    local$broad,
    local_id = local$local_id,
    source = "broad_scout"
  )
  components <- list(broad_component)
  if (!is.null(local$predictive_component)) {
    components <- c(components, list(local$predictive_component))
  }
  reference_local_object_from_components(
    local_id = local$local_id,
    components = components
  )
}

attach_predictive_repaired_factor <- function(local, local_object) {
  local$predictive_repaired_object <- local_object
  local$predictive_repaired_factor <- reference_local_object_factor(
    local_object = local_object,
    population_model = population_model,
    data_i = local$data_i,
    loglik_fn = loglik_shifted_gamma,
    local_n_cores = 1L
  )
  local
}

plot_support_summary <- function(summary_df, q0_df, theta_scout_df, fitted_q0_df, file) {
  methods <- unique(summary_df$method)
  med_abs <- tapply(summary_df$abs_error, summary_df$method, median, na.rm = TRUE)
  q90_abs <- tapply(summary_df$abs_error, summary_df$method, function(x) stats::quantile(x, 0.90, na.rm = TRUE))
  med_ess <- tapply(summary_df$ess_frac, summary_df$method, median, na.rm = TRUE)
  q0_cov <- tapply(q0_df$stan_inside_q01_q99, q0_df$method, mean, na.rm = TRUE)
  scout_cov <- tapply(theta_scout_df$stan_inside_q01_q99, theta_scout_df$method, mean, na.rm = TRUE)
  fitted_cov <- tapply(fitted_q0_df$reference_inside_q01_q99, fitted_q0_df$method, mean, na.rm = TRUE)
  methods <- methods[order(med_abs[methods])]

  grDevices::png(file, width = 1600, height = 1000)
  graphics::par(mfrow = c(2, 2), mar = c(8, 4, 3, 1))
  graphics::barplot(
    med_abs[methods],
    names.arg = methods,
    las = 2,
    ylab = "median |log m error|",
    main = "Local Marginal Error"
  )
  graphics::barplot(
    q90_abs[methods],
    names.arg = methods,
    las = 2,
    ylab = "q90 |log m error|",
    main = "Tail Local Marginal Error"
  )
  graphics::barplot(
    med_ess[methods],
    names.arg = methods,
    las = 2,
    ylab = "median ESS fraction",
    main = "Support ESS"
  )
  support_methods <- intersect(methods, names(med_abs))
  proposal_methods <- intersect(methods, names(q0_cov))
  graphics::barplot(
    c(
      q0_cov[proposal_methods],
      scout_cov[paste0("theta_score_", proposal_methods)],
      fitted_cov[paste0("fitted_q0_", proposal_methods)]
    ),
    names.arg = c(paste0("raw_", proposal_methods), paste0("score_", proposal_methods), paste0("q0_", proposal_methods)),
    las = 2,
    ylab = "mean Stan coverage by 1%-99%",
    ylim = c(0, 1),
    main = "Raw, Scored, and Fitted q0 Support"
  )
  grDevices::dev.off()
}

if (isTRUE(verbose)) {
  cat(sprintf("Shifted-gamma support diagnostics\n"))
  cat(sprintf("Locals: %d | theta stress points: %d\n", length(data_list), nrow(theta_stress)))
  cat(sprintf("SMC particles: %d | Laplace draws: %d | QMC: %d x %d\n",
              smc_particles, laplace_draws, qmc_size, qmc_randomizations))
  cat(sprintf("Raw q0 configs: %d | extra theta candidates: %d\n",
              q0_configs, theta_candidate_configs))
}

local_ids <- seq_along(data_list)
locals <- parallel::mclapply(
  local_ids,
  fit_local_support_objects,
  mc.cores = as.integer(max(1L, n_cores))
)

methods <- c("initial_anchor", "broad_scout", "laplace_t", "dmis_broad_laplace")
q0_draws <- lapply(methods, function(method) {
  set.seed(seed + match(method, methods) * 100000L)
  sample_q0_theta(locals, method, q0_configs)
})
names(q0_draws) <- methods
q0_df <- do.call(rbind, Map(summarize_q0, q0_draws, names(q0_draws)))

theta_candidates <- make_theta_candidate_pool(
  q0_draws = q0_draws,
  n_extra = theta_candidate_configs,
  seed_i = seed + 700000L
)
theta_scout_scores <- lapply(methods, function(method) {
  score_theta_candidate_pool(theta_candidates, locals, method)
})
names(theta_scout_scores) <- methods
theta_scout_df <- do.call(rbind, Map(
  summarize_weighted_theta_support,
  theta_score = theta_scout_scores,
  method = names(theta_scout_scores),
  MoreArgs = list(theta_candidates = theta_candidates)
))

fitted_q0_proposals <- lapply(methods, function(method) {
  fit_scored_q0_proposal(
    theta_candidates = theta_candidates,
    theta_score = theta_scout_scores[[method]],
    method = method,
    seed_i = seed + 800000L + match(method, methods)
  )
})
names(fitted_q0_proposals) <- methods
fitted_q0_df <- do.call(rbind, Map(
  summarize_theta_proposal_reference,
  proposal = fitted_q0_proposals,
  label = paste0("fitted_q0_", names(fitted_q0_proposals)),
  MoreArgs = list(
    reference_theta = stan_theta,
    n_draws = q0_configs,
    seed = seed + 810000L
  )
))
fitted_q0_importance_df <- do.call(rbind, Map(
  summarize_fitted_q0_importance,
  proposal = fitted_q0_proposals,
  method = names(fitted_q0_proposals),
  seed_i = seed + 820000L + seq_along(fitted_q0_proposals),
  MoreArgs = list(
    locals = locals,
    n_draws = q0_configs
  )
))

predictive_theta_audit <- NULL
predictive_audit_df <- data.frame()
predictive_repair_df <- data.frame()
predictive_local_objects <- NULL

if (isTRUE(predictive_enabled)) {
  if (isTRUE(verbose)) {
    cat(sprintf(
      "Adding q0-predictive local components from fitted q0 '%s' (%d alpha components, %d particles)\n",
      predictive_q0_method,
      predictive_components,
      predictive_particles
    ))
  }
  locals <- parallel::mclapply(
    locals,
    add_predictive_component_to_local,
    mc.cores = as.integer(max(1L, n_cores))
  )
  predictive_local_objects <- lapply(locals, predictive_reference_object_from_local)

  if (isTRUE(predictive_repair_enabled)) {
    if (isTRUE(verbose)) {
      cat(sprintf(
        "Auditing q0-predictive local objects over %d fitted-q0 theta draws (ESS threshold %.3f)\n",
        predictive_audit_points,
        predictive_audit_ess
      ))
    }
    predictive_theta_audit <- theta_proposal_sample(
      fitted_q0_proposals[[predictive_q0_method]],
      n = as.integer(predictive_audit_points),
      seed = seed + 850000L
    )
    predictive_audit_df <- audit_reference_local_objects(
      local_objects = predictive_local_objects,
      theta_audit = predictive_theta_audit,
      population_model = population_model,
      data_list = data_list,
      loglik_fn = loglik_shifted_gamma,
      min_ess_frac = predictive_audit_ess,
      n_jobs = n_cores,
      local_n_cores = 1L
    )
    if (isTRUE(verbose)) {
      cat(sprintf("Predictive audit failures: %d / %d\n", sum(!predictive_audit_df$covered), nrow(predictive_audit_df)))
    }
    repair <- bridge_repair_reference_local_objects(
      local_objects = predictive_local_objects,
      audit = predictive_audit_df,
      theta_audit = predictive_theta_audit,
      data_list = data_list,
      loglik_fn = loglik_shifted_gamma,
      population_model = population_model,
      max_repairs = predictive_max_repairs,
      M = predictive_repair_particles,
      min_ess_frac = predictive_audit_ess,
      reference_scale = 1,
      resample_threshold = 0.6,
      n_mcmc_moves = 1L,
      max_rounds = 80L,
      cess_target = 0.9,
      G_mix = 8L,
      rw_scale_init = 0.7,
      n_jobs = n_cores,
      local_n_cores = 1L,
      seed = seed + 860000L,
      verbose = FALSE
    )
    predictive_local_objects <- repair$local_objects
    predictive_repair_df <- repair$repairs
    if (isTRUE(verbose)) {
      cat(sprintf("Predictive bridge repairs added: %d\n", nrow(predictive_repair_df)))
    }
    locals <- Map(attach_predictive_repaired_factor, locals, predictive_local_objects)
  }
}

support_parts <- parallel::mclapply(
  locals,
  evaluate_local_support,
  mc.cores = as.integer(max(1L, n_cores))
)
support_df <- do.call(rbind, support_parts)

support_summary <- do.call(rbind, lapply(split(support_df, support_df$method), function(x) {
  data.frame(
    method = x$method[1L],
    n = nrow(x),
    median_abs_error = median(x$abs_error, na.rm = TRUE),
    q90_abs_error = as.numeric(stats::quantile(x$abs_error, 0.90, na.rm = TRUE, names = FALSE)),
    max_abs_error = max(x$abs_error, na.rm = TRUE),
    median_ess_frac = median(x$ess_frac, na.rm = TRUE),
    q10_ess_frac = as.numeric(stats::quantile(x$ess_frac, 0.10, na.rm = TRUE, names = FALSE)),
    failure_ess_lt_01 = mean(x$ess_frac < 0.01, na.rm = TRUE),
    failure_abs_error_gt_1 = mean(x$abs_error > 1, na.rm = TRUE),
    check.names = FALSE
  )
}))
support_summary <- support_summary[order(support_summary$median_abs_error), , drop = FALSE]

utils::write.csv(support_df, summary_csv, row.names = FALSE)
utils::write.csv(q0_df, q0_csv, row.names = FALSE)
utils::write.csv(theta_scout_df, theta_scout_csv, row.names = FALSE)
utils::write.csv(fitted_q0_df, fitted_q0_csv, row.names = FALSE)
utils::write.csv(fitted_q0_importance_df, fitted_q0_importance_csv, row.names = FALSE)
utils::write.csv(predictive_audit_df, predictive_audit_csv, row.names = FALSE)
utils::write.csv(predictive_repair_df, predictive_repair_csv, row.names = FALSE)
plot_support_summary(support_df, q0_df, theta_scout_df, fitted_q0_df, plot_file)

saveRDS(
  list(
    settings = list(
      stan_results_file = normalizePath(stan_results_file, winslash = "/", mustWork = TRUE),
      seed = seed,
      n_cores = n_cores,
      locals = length(data_list),
      theta_points = nrow(theta_stress),
      smc_particles = smc_particles,
      laplace_draws = laplace_draws,
      qmc_size = qmc_size,
      qmc_randomizations = qmc_randomizations,
      q0_configs = q0_configs,
      theta_candidate_configs = theta_candidate_configs,
      predictive_enabled = predictive_enabled,
      predictive_q0_method = predictive_q0_method,
      predictive_particles = predictive_particles,
      predictive_components = predictive_components,
      predictive_bridge = predictive_bridge,
      predictive_repair_enabled = predictive_repair_enabled,
      predictive_audit_points = predictive_audit_points,
      predictive_audit_ess = predictive_audit_ess,
      predictive_max_repairs = predictive_max_repairs,
      predictive_repair_particles = predictive_repair_particles,
      laplace_df = laplace_df,
      laplace_scale = laplace_scale,
      broad_scale = broad_scale,
      defensive_scale = defensive_scale
    ),
    theta_stress = theta_stress,
    theta_stress_source_rows = attr(theta_stress, "source_rows"),
    support = support_df,
    support_summary = support_summary,
    q0 = q0_df,
    q0_draws = q0_draws,
    theta_candidates = theta_candidates,
    theta_scout = theta_scout_df,
    theta_scout_scores = theta_scout_scores,
    fitted_q0 = fitted_q0_df,
    fitted_q0_importance = fitted_q0_importance_df,
    fitted_q0_proposals = fitted_q0_proposals,
    predictive_theta_audit = predictive_theta_audit,
    predictive_audit = predictive_audit_df,
    predictive_repairs = predictive_repair_df,
    predictive_local_objects = predictive_local_objects,
    locals = locals
  ),
  results_file
)

cat("Saved support rows to:", summary_csv, "\n")
cat("Saved q0 rows to:", q0_csv, "\n")
cat("Saved theta scout rows to:", theta_scout_csv, "\n")
cat("Saved fitted q0 rows to:", fitted_q0_csv, "\n")
cat("Saved fitted q0 importance rows to:", fitted_q0_importance_csv, "\n")
cat("Saved predictive audit rows to:", predictive_audit_csv, "\n")
cat("Saved predictive repair rows to:", predictive_repair_csv, "\n")
cat("Saved plot to:", plot_file, "\n")
cat("Saved RDS to:", results_file, "\n")
cat("\nSupport summary:\n")
print(support_summary, row.names = FALSE)
cat("\nTheta scout summary:\n")
print(theta_scout_df[order(theta_scout_df$method, theta_scout_df$parameter), ], row.names = FALSE)
cat("\nFitted q0 summary:\n")
print(fitted_q0_df[order(fitted_q0_df$method, fitted_q0_df$parameter), ], row.names = FALSE)
cat("\nFitted q0 importance summary:\n")
print(fitted_q0_importance_df, row.names = FALSE)
if (nrow(predictive_audit_df)) {
  cat("\nPredictive audit summary:\n")
  print(data.frame(
    n = nrow(predictive_audit_df),
    failures = sum(!predictive_audit_df$covered),
    median_ess_frac = median(predictive_audit_df$ess_frac, na.rm = TRUE),
    q10_ess_frac = as.numeric(stats::quantile(predictive_audit_df$ess_frac, 0.10, na.rm = TRUE, names = FALSE)),
    check.names = FALSE
  ), row.names = FALSE)
}
if (nrow(predictive_repair_df)) {
  cat("\nPredictive repairs:\n")
  print(predictive_repair_df, row.names = FALSE)
}
