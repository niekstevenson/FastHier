#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) > 0L) {
  script_path <- normalizePath(sub("^--file=", "", file_arg[1L]), winslash = "/", mustWork = TRUE)
  setwd(normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE))
}

source("smc_core.R")
source("SMC_super_fast.R")
source("new_SMC_cache.R")

results <- list()

run_test <- function(name, expr) {
  err <- NULL
  tryCatch(
    force(expr),
    error = function(e) err <<- e
  )
  if (is.null(err)) {
    cat(sprintf("[PASS] %s\n", name))
    results[[name]] <<- TRUE
  } else {
    cat(sprintf("[FAIL] %s: %s\n", name, conditionMessage(err)))
    results[[name]] <<- FALSE
  }
}

id_tmap <- list(
  fwd = function(X) as.matrix(X),
  inv = function(Z) as.matrix(Z),
  log_jac = function(X) rep(0, nrow(X))
)

toy_ll <- function(Theta, data) {
  X <- as.matrix(Theta)
  if (ncol(X) != length(data$mu) && nrow(X) == length(data$mu)) {
    X <- t(X)
  }
  X <- sweep(X, 2L, data$mu, `-`)
  vals <- -0.5 * rowSums(X^2) / (data$sigma^2)
  as.numeric(vals)
}

run_test("1) cache regression with duplicates", {
  ll_parallel_orig <- ll_parallel
  on.exit(assign("ll_parallel", ll_parallel_orig, envir = .GlobalEnv), add = TRUE)
  assign(
    "ll_parallel",
    function(Theta, data, loglik_fn, n_cores = 1L) loglik_fn(Theta, data),
    envir = .GlobalEnv
  )

  set.seed(1)
  Theta <- matrix(rnorm(40), nrow = 20, ncol = 2)
  Theta[11:20, ] <- Theta[1:10, ]
  cache <- .ll_cache_make(digits = 8L, cap = 500L, prune_every = 2L)
  ll_fn <- function(Theta, data) {
    X <- as.matrix(Theta)
    if (ncol(X) == 1L && length(Theta) == 2L) X <- t(X)
    as.numeric(rowSums(X))
  }

  out1 <- .ll_cached_eval(Theta, NULL, ll_fn, cache = cache, expect_dups = TRUE, n_cores = 1L)
  out2 <- .ll_cached_eval(Theta, NULL, ll_fn, cache = cache, expect_dups = TRUE, n_cores = 1L)

  stopifnot(length(out1) == nrow(Theta), length(out2) == nrow(Theta))
  stopifnot(all(is.finite(out1)), all(is.finite(out2)))
  stopifnot(is.environment(cache$env))
  stopifnot(is.environment(cache$state))
})

run_test("2) pCN DA vs non-DA consistency", {
  set.seed(11)
  d <- 2L
  N <- 400L
  Z0 <- matrix(rnorm(N * d), nrow = N, ncol = d)
  mu_ref <- c(0, 0)
  prior_L <- chol(diag(d))
  w <- rep(1 / N, N)
  data <- list(mu = c(0.7, -0.3), sigma = 1.2)
  loglik0 <- toy_ll(Z0, data)
  lpz0 <- as.numeric(dmvnorm_chol_log(Z0, mu_ref, prior_L))
  elite <- .default_std_normal_mix(d)

  res_no_da <- mcmc_moves_z_mix_batched(
    Z = Z0, loglik = loglik0, lpz = lpz0, Tmap = id_tmap, lambda = 1.0,
    mu_ref = mu_ref, prior_L = prior_L, w = w,
    elite_mix = elite, hist_mix = NULL,
    data = data, loglik_fn = toy_ll,
    n_moves = 30L, rw_prob = 0, pcn_prob = 1, pcn_beta = 0.35,
    indep_t_prob = 0, seed = 9001L, da_enable = FALSE, ll_cache = NULL, n_cores = 1L
  )

  res_da <- mcmc_moves_z_mix_batched(
    Z = Z0, loglik = loglik0, lpz = lpz0, Tmap = id_tmap, lambda = 1.0,
    mu_ref = mu_ref, prior_L = prior_L, w = w,
    elite_mix = elite, hist_mix = NULL,
    data = data, loglik_fn = toy_ll,
    n_moves = 30L, rw_prob = 0, pcn_prob = 1, pcn_beta = 0.35,
    indep_t_prob = 0, seed = 9001L,
    da_enable = TRUE, da_screen_mix = elite, da_alpha = 1.0, da_calib = c(0, 1),
    ll_cache = NULL, n_cores = 1L
  )

  m_no <- colMeans(res_no_da$Z)
  m_da <- colMeans(res_da$Z)
  v_no <- apply(res_no_da$Z, 2L, var)
  v_da <- apply(res_da$Z, 2L, var)

  stopifnot(max(abs(m_no - m_da)) < 0.20)
  stopifnot(max(abs(v_no - v_da)) < 0.25)
})

run_test("3) GSS refit invalidation wiring", {
  txt <- readLines("SMC_super_fast.R", warn = FALSE)
  i_refit <- grep("if \\(isTRUE\\(ref\\$refit\\)\\)", txt)
  i_refmix_null <- grep("ref_mix <- NULL", txt)
  i_lpref_null <- grep("lp_ref <- NULL", txt)
  stopifnot(length(i_refit) >= 1L, length(i_refmix_null) >= 1L, length(i_lpref_null) >= 1L)
  stopifnot(any(i_refmix_null > i_refit[1L]))
  stopifnot(any(i_lpref_null > i_refit[1L]))
})

run_test("4) randomized count allocation is unbiased in expectation", {
  meansZ <- list(c(-1, 0), c(1, 0))
  covsZ <- list(diag(c(1, 0.5)), diag(c(0.8, 1.2)))
  wZ <- c(0.35, 0.65)

  mix_g <- list(meansZ = meansZ, covsZ = covsZ, wZ = wZ)
  mix_g$cache <- prep_mix_cache(mix_g$meansZ, mix_g$covsZ, mix_g$wZ)

  ri <- list(
    mix_g = mix_g,
    mix_t = mix_g,
    t_eps = 0.25,
    t_df = 5L,
    Tmap = id_tmap
  )

  M <- 256L
  R <- 500L
  acc_counts <- numeric(4L)
  for (r in seq_len(R)) {
    smp <- .sample_from_ri_rqmc(M, ri, sobol_seed = 1000L + r, deterministic_counts = FALSE)
    acc_counts <- acc_counts + smp$strata_counts
  }

  emp <- acc_counts / sum(acc_counts)
  target <- c((1 - ri$t_eps) * wZ, ri$t_eps * wZ)
  stopifnot(max(abs(emp - target)) < 0.03)

  d1 <- .sample_from_ri_rqmc(M, ri, sobol_seed = 999L, deterministic_counts = TRUE)$strata_counts
  d2 <- .sample_from_ri_rqmc(M, ri, sobol_seed = 999L, deterministic_counts = TRUE)$strata_counts
  stopifnot(all(d1 == d2))
})

run_test("5) resampling stochasticity toggles correctly", {
  w <- c(0.05, 0.1, 0.2, 0.25, 0.4)

  set.seed(1); a1 <- stratified_resample_sorted(w, deterministic = FALSE)
  set.seed(2); a2 <- stratified_resample_sorted(w, deterministic = FALSE)
  stopifnot(!identical(a1, a2))

  set.seed(1); d1 <- stratified_resample_sorted(w, deterministic = TRUE)
  set.seed(2); d2 <- stratified_resample_sorted(w, deterministic = TRUE)
  stopifnot(identical(d1, d2))
})

run_test("6) end-to-end smoke with DA + ll_cache", {
  data <- list(mu = c(0.2, -0.1), sigma = 1.0)
  mu_ref <- c(theta1 = 0, theta2 = 0)
  Sigma_ref <- diag(c(1, 1))

  fit <- enhanced_smc_elite(
    data = data,
    loglik_fn = toy_ll,
    mu_ref = mu_ref,
    Sigma_ref = Sigma_ref,
    M = 128L,
    resample_threshold = 0.6,
    n_mcmc_moves = 2L,
    max_rounds = 80L,
    G_mix = 6L,
    hist_mix_enable = TRUE,
    gss_enable = TRUE,
    da_enable = TRUE,
    ll_cache_enable = TRUE,
    deterministic_resampling = FALSE,
    n_cores = 1L,
    seed = 123,
    verbose = FALSE
  )

  stopifnot(is.finite(fit$log_evidence))
  stopifnot(is.finite(fit$mcse_logZ))
  stopifnot(is.finite(fit$final_lambda))
  stopifnot(fit$final_lambda > 1 - 1e-10)
})

n_ok <- sum(unlist(results))
n_all <- length(results)
cat(sprintf("\nSummary: %d/%d tests passed\n", n_ok, n_all))
if (n_ok != n_all) quit(status = 1L)
