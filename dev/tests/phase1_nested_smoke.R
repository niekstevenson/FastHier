#!/usr/bin/env Rscript

source("smc_core.R")
source("SMC_super_fast.R")
source("new_SMC_cache.R")
source("nested_population_SMC.R")
source("make_prior.R")

set.seed(20260311L)

S <- 3L
N_obs <- 20L
sigma_y <- 1.0
mu_true <- 1.2
tau2_true <- 0.4

alpha_true <- rnorm(S, mean = mu_true, sd = sqrt(tau2_true))
data_list <- lapply(seq_len(S), function(i) rnorm(N_obs, mean = alpha_true[i], sd = sigma_y))

loglik_fn <- function(Theta, y_i) {
  Theta <- as.matrix(Theta)
  alpha <- Theta[, 1L]
  vapply(alpha, function(a) sum(dnorm(y_i, mean = a, sd = sigma_y, log = TRUE)), numeric(1L))
}

gaussian_map_fn <- function(phi, d = 1L) {
  if (!is.null(d) && as.integer(d) != 1L) stop("This smoke test expects d = 1.")
  mu <- as.numeric(phi[1L])
  tau2 <- exp(as.numeric(phi[2L]))
  list(
    mu = mu,
    Sigma_inv = matrix(1 / tau2, nrow = 1L, ncol = 1L),
    logdet = log(tau2),
    const = NULL
  )
}

mu_ref <- setNames(mu_true, "alpha")
Sigma_ref <- matrix(tau2_true, nrow = 1L, ncol = 1L, dimnames = list("alpha", "alpha"))
Sigma_ref_inv <- chol2inv(chol(Sigma_ref))
Sigma_ref_logdet <- as.numeric(determinant(Sigma_ref, logarithm = TRUE)$modulus)

outer_subject <- lapply(
  seq_len(S),
  function(i) {
    out <- enhanced_smc_elite(
      data = data_list[[i]],
      loglik_fn = loglik_fn,
      mu_ref = mu_ref,
      Sigma_ref = Sigma_ref,
      M = 300L,
      max_rounds = 40L,
      seed = as.integer(100L + i),
      verbose = FALSE
    )
    out$working_prior <- list(
      mu = as.numeric(mu_ref),
      Sigma = Sigma_ref,
      Sigma_inv = Sigma_ref_inv,
      logdet = Sigma_ref_logdet
    )
    out
  }
)

local_objs <- lapply(
  seq_len(S),
  function(i) {
    build_local_exact_object(
      smc_out = outer_subject[[i]],
      data = data_list[[i]],
      subj_id = i,
      loglik_fn = loglik_fn
    )
  }
)

prior <- make_prior_phi_diag(
  m0 = mu_true,
  s0 = 1.0,
  a = 5.0,
  b = (5.0 - 1.0) * tau2_true,
  d = 1L
)

fit <- nested_population_smc(
  rprior_phi = prior$rprior,
  local_objs = local_objs,
  logprior_phi = prior$lprior,
  gaussian_map_fn = gaussian_map_fn,
  N = 48L,
  M_local = 64L,
  rho_step = NULL,
  rho_res = 0.5,
  rho_local = 0.5,
  n_population_moves = 1L,
  n_local_moves = 0L,
  n_population_refresh_moves = 0L,
  max_rounds = 10L,
  seed = 77L,
  verbose = FALSE
)

lambda_final <- tail(fit$meta$lambda_hist, 1L)
eval_counts <- fit$meta$local_loglik_evals
stopifnot(is.finite(lambda_final), lambda_final >= 1 - 1e-8)
stopifnot(all(is.finite(fit$phi)))
stopifnot(all(is.finite(fit$w)))
stopifnot(abs(sum(fit$w) - 1) < 1e-8)
stopifnot(identical(fit$meta$implementation, "shared_local_banks_baseline"))
stopifnot(eval_counts["total"] >= eval_counts["initialization"])
stopifnot(eval_counts["enrichment"] == 0L)

cat("PHASE1_NESTED_SMOKE_OK\n")
cat(sprintf("lambda_final=%.6f\n", lambda_final))
cat(sprintf("rounds=%d\n", fit$meta$rounds))
cat(sprintf("local_loglik_total=%d\n", eval_counts["total"]))
cat(sprintf("local_loglik_enrichment=%d\n", eval_counts["enrichment"]))
