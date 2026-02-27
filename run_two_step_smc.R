#!/usr/bin/env Rscript

# Inner mode toggle:
pm_mode <- "strict"  # options: "strict", "fast"

suppressPackageStartupMessages({
  library(parallel)
})

source("smc_core.R")
source("SMC_super_fast.R")
source("new_SMC_cache.R")
source("outer_SMC.R")
source("make_prior.R")

if (!pm_mode %in% c("strict", "fast")) stop("pm_mode must be 'strict' or 'fast'")

stan_path <- file.path("samples", "stan_full_results.rds")
if (!file.exists(stan_path)) {
  stop("Missing Stan results at ", stan_path, ". Run run_stan.R first.")
}
stan_res <- readRDS(stan_path)

dir.create("samples", showWarnings = FALSE, recursive = TRUE)
y_mat <- stan_res$data$y
S <- nrow(y_mat)
data_list <- lapply(seq_len(S), function(i) as.numeric(y_mat[i, ]))

loglik_shifted_gamma <- function(Theta, y) {
  Theta <- as.matrix(Theta)
  eps <- 1e-9
  shape <- exp(Theta[, 1L]) + eps
  scale <- exp(Theta[, 2L]) + eps
  shift <- exp(Theta[, 3L]) + eps
  out <- rep(-1e12, nrow(Theta))
  min_y <- min(y)
  ok <- shift < min_y
  if (!any(ok)) return(out)
  for (i in which(ok)) out[i] <- sum(dgamma(y - shift[i], shape = shape[i], scale = scale[i], log = TRUE))
  out[!is.finite(out)] <- -1e12
  out
}

theta_names <- c("eta_shape", "eta_scale", "eta_shift")
mu_ref <- setNames(rep(0, 3L), theta_names)
Sigma_ref <- diag(1, 3L)
colnames(Sigma_ref) <- rownames(Sigma_ref) <- theta_names

# Outer (subject-level) SMC settings
M_subject <- 300L

cat(sprintf("Running subject-level outer SMC for %d subjects...\n", S))
outer_subject <- mclapply(
  seq_len(S),
  function(i) {
    enhanced_smc_elite(
      data = data_list[[i]],
      loglik_fn = loglik_shifted_gamma,
      mu_ref = mu_ref,
      Sigma_ref = Sigma_ref,
      M = M_subject,
      resample_threshold = 0.6,
      n_mcmc_moves = 2L,
      max_rounds = 150L,
      G_mix = 8L,
      da_enable = TRUE,
      gss_enable = TRUE,
      ll_cache_enable = TRUE,
      deterministic_resampling = FALSE,
      n_cores = 1L,
      seed = 1000L + i,
      verbose = FALSE
    )
  },
  mc.cores = 1L
)
saveRDS(outer_subject, file.path("samples", "outer_smc_results.rds"))
cat("Saved: samples/outer_smc_results.rds\n")

cat("Building subject caches...\n")
caches <- mclapply(
  seq_len(S),
  function(i) {
    build_subject_cache_from_smc(
      smc_out = outer_subject[[i]],
      data = data_list[[i]],
      subj_id = i,
      loglik_fn = loglik_shifted_gamma,
      M = 256L,
      K_batches = 4L,
      deterministic_counts = FALSE,
      n_cores = 1L
    )
  },
  mc.cores = 1L
)

prior <- make_prior_phi_diag(
  m0 = stan_res$priors$m0,
  s0 = stan_res$priors$s0,
  a = stan_res$priors$a0,
  b = stan_res$priors$b0,
  d = 3L
)
gaussian_map_fn <- phi_to_gaussian_params_diag_factory()

# Inner SMC settings
M_inner <- 1200L

cat(sprintf("Running inner SMC (pm_mode=%s)...\n", pm_mode))
inner <- outer_smc_phi_batch(
  caches = caches,
  rprior_phi = prior$rprior,
  logprior_phi = prior$lprior,
  gaussian_map_fn = gaussian_map_fn,
  M = M_inner,
  cess_target = 0.95,
  resample_threshold = 0.5,
  n_moves = 2L,
  rw_scale_init = 1.2,
  pm_mode = pm_mode,
  block_refresh_every = if (pm_mode == "strict") 5L else 0L,
  block_refresh_frac = 0.10,
  auto_enrich_enable = FALSE,
  diag_enable = FALSE,
  max_rounds = 200L,
  seed = 123,
  verbose = FALSE
)

saveRDS(inner, file.path("samples", "inner_smc_results.rds"))
cat("Saved: samples/inner_smc_results.rds\n")
cat(sprintf("Inner log evidence: %.6f (MCSE %.6f)\n", inner$log_evidence, inner$mcse_log_evidence))
