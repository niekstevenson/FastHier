#!/usr/bin/env Rscript

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[1L]))
} else {
  normalizePath("benchmarks/nested_hierarchical_shifted_gamma.R")
}
repo_dir <- dirname(dirname(script_path))
setwd(repo_dir)

source("smc_core.R")
source("SMC_super_fast.R")
source("new_SMC_cache.R")
source("nested_population_SMC.R")
source("make_prior.R")

dir.create("benchmarks", showWarnings = FALSE, recursive = TRUE)
dir.create("samples", showWarnings = FALSE, recursive = TRUE)

stan_results_file <- file.path("samples", "stan_full_results.rds")
if (!file.exists(stan_results_file)) {
  stop("Missing Stan results. Run `Rscript run_stan.R` first.")
}

stan_res <- readRDS(stan_results_file)

S <- as.integer(stan_res$config$S)
P <- as.integer(stan_res$config$P)
y <- stan_res$data$y
data_list <- lapply(seq_len(S), function(i) as.numeric(y[i, ]))

prior_m0 <- as.numeric(stan_res$priors$m0)
prior_s0 <- as.numeric(stan_res$priors$s0)
prior_a0 <- as.numeric(stan_res$priors$a0)
prior_b0 <- as.numeric(stan_res$priors$b0)

local_particles <- 800L
local_rounds <- 180L
population_particles <- 1000L
local_bank_particles <- 64L
population_rounds <- 50L
mc_cores <- 1L

theta_names <- c("eta_shape", "eta_scale", "eta_shift")
mu_names <- c("mu_shape", "mu_scale", "mu_shift")
sigma2_names <- c("sigma2_shape", "sigma2_scale", "sigma2_shift")

loglik_shifted_gamma <- function(Theta, y_i) {
  Theta <- as.matrix(Theta)
  eps <- 1e-9
  shape <- exp(Theta[, 1L]) + eps
  scale <- exp(Theta[, 2L]) + eps
  shift <- exp(Theta[, 3L]) + eps
  out <- rep(-1e12, nrow(Theta))
  min_y <- min(y_i)
  ok <- shift < min_y
  if (!any(ok)) return(out)
  for (j in which(ok)) {
    out[j] <- sum(dgamma(y_i - shift[j], shape = shape[j], scale = scale[j], log = TRUE))
  }
  out[!is.finite(out)] <- -1e12
  out
}

gaussian_map_fn <- phi_to_gaussian_params_diag_factory()

alpha_ref_mean <- setNames(prior_m0, theta_names)
alpha_ref_var <- prior_s0 + prior_b0 / (prior_a0 - 1.0)
Sigma_ref <- diag(alpha_ref_var, P)
colnames(Sigma_ref) <- rownames(Sigma_ref) <- theta_names
Sigma_ref_inv <- chol2inv(chol(Sigma_ref))
Sigma_ref_logdet <- as.numeric(determinant(Sigma_ref, logarithm = TRUE)$modulus)

ll_eval_counter_reset()
local_fits <- parallel::mclapply(
  seq_len(S),
  function(i) {
    out <- enhanced_smc_elite(
      data = data_list[[i]],
      loglik_fn = loglik_shifted_gamma,
      mu_ref = alpha_ref_mean,
      Sigma_ref = Sigma_ref,
      M = local_particles,
      resample_threshold = 0.6,
      n_mcmc_moves = 3L,
      max_rounds = local_rounds,
      G_mix = 8L,
      da_enable = TRUE,
      gss_enable = TRUE,
      ll_cache_enable = TRUE,
      deterministic_resampling = FALSE,
      n_cores = 1L,
      seed = as.integer(1123L + 31L * i),
      verbose = FALSE
    )
    out$working_prior <- list(
      mu = as.numeric(alpha_ref_mean),
      Sigma = Sigma_ref,
      Sigma_inv = Sigma_ref_inv,
      logdet = Sigma_ref_logdet
    )
    out
  },
  mc.cores = mc_cores
)
local_prefit_loglik <- ll_eval_counter_get()
ll_eval_counter_disable()

local_objs <- lapply(
  seq_len(S),
  function(i) {
    build_local_exact_object(
      smc_out = local_fits[[i]],
      data = data_list[[i]],
      subj_id = i,
      loglik_fn = loglik_shifted_gamma,
      base_seed = as.integer(40000L + 1009L * i)
    )
  }
)

prior <- make_prior_phi_diag(
  m0 = prior_m0,
  s0 = prior_s0,
  a = prior_a0,
  b = prior_b0,
  d = P
)

nested_fit <- nested_population_smc(
  rprior_phi = prior$rprior,
  local_objs = local_objs,
  logprior_phi = prior$lprior,
  gaussian_map_fn = gaussian_map_fn,
  N = population_particles,
  M_local = local_bank_particles,
  rho_res = 0.5,
  rho_local = 0.5,
  n_population_moves = 1L,
  n_local_moves = 0L,
  n_population_refresh_moves = 0L,
  max_rounds = population_rounds,
  max_bank_topups = 3L,
  bank_split_tol = 0.01,
  seed = 123L,
  verbose = TRUE
)

stan_mu <- as.matrix(stan_res$draws$mu)
stan_sigma2 <- as.matrix(stan_res$draws$sigma2)
colnames(stan_mu) <- mu_names
colnames(stan_sigma2) <- sigma2_names

nested_phi <- as.matrix(nested_fit$phi)
nested_w <- pmax(as.numeric(nested_fit$w), 0)
nested_w <- nested_w / sum(nested_w)
nested_mu <- nested_phi[, seq_len(P), drop = FALSE]
nested_sigma2 <- exp(nested_phi[, P + seq_len(P), drop = FALSE])
colnames(nested_mu) <- mu_names
colnames(nested_sigma2) <- sigma2_names

set.seed(123L)
draw_idx <- sample.int(nrow(nested_phi), size = 4000L, replace = TRUE, prob = nested_w)
nested_mu_draws <- nested_mu[draw_idx, , drop = FALSE]
nested_sigma2_draws <- nested_sigma2[draw_idx, , drop = FALSE]

plot_file <- file.path("benchmarks", "nested_hierarchical_shifted_gamma_posteriors.png")
png(plot_file, width = 1200, height = 800, res = 120)
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))

for (j in seq_len(P)) {
  plot(
    density(stan_mu[, j]),
    main = mu_names[j],
    xlab = mu_names[j],
    ylab = "Density",
    lwd = 2,
    col = "black"
  )
  lines(density(nested_mu_draws[, j]), lwd = 2, col = "red3")
  legend("topright", legend = c("Stan", "Nested SMC"), col = c("black", "red3"), lwd = 2, bty = "n")
}

for (j in seq_len(P)) {
  plot(
    density(stan_sigma2[, j]),
    main = sigma2_names[j],
    xlab = sigma2_names[j],
    ylab = "Density",
    lwd = 2,
    col = "black"
  )
  lines(density(nested_sigma2_draws[, j]), lwd = 2, col = "red3")
  legend("topright", legend = c("Stan", "Nested SMC"), col = c("black", "red3"), lwd = 2, bty = "n")
}

dev.off()

results <- list(
  stan_results_file = stan_results_file,
  local_fits = local_fits,
  nested_fit = nested_fit,
  plot_file = plot_file,
  settings = list(
    local_particles = local_particles,
    local_rounds = local_rounds,
    population_particles = population_particles,
    local_bank_particles = local_bank_particles,
    population_rounds = population_rounds,
    alpha_ref_mean = unname(alpha_ref_mean),
    alpha_ref_var = unname(alpha_ref_var),
    max_bank_topups = 3L,
    bank_split_tol = 0.01,
    mc_cores = mc_cores
  )
)

outfile <- file.path("samples", "nested_hierarchical_shifted_gamma_results.rds")
saveRDS(results, outfile)

eval_counts <- nested_fit$meta$local_loglik_evals
local_end_to_end <- as.integer(local_prefit_loglik + eval_counts["total"])

stan_mu_mean <- colMeans(stan_mu)
stan_sigma2_mean <- colMeans(stan_sigma2)
nested_mu_mean <- colSums(nested_mu * nested_w)
nested_sigma2_mean <- colSums(nested_sigma2 * nested_w)

cat("NESTED_HIERARCHICAL_SHIFTED_GAMMA_OK\n")
cat(sprintf("stan_results=%s\n", stan_results_file))
cat(sprintf("results=%s\n", outfile))
cat(sprintf("plot=%s\n", plot_file))
cat(sprintf("population_final_lambda=%.4f\n", tail(nested_fit$meta$lambda_hist, 1L)))
cat(sprintf("population_rounds=%d\n", nested_fit$meta$rounds))
cat(sprintf("implementation=%s\n", nested_fit$meta$implementation))
cat(sprintf("bank_topup_points=%d\n", nested_fit$meta$bank_topup_points))
cat(sprintf("bank_mean_size_final=%.1f\n", tail(nested_fit$meta$bank_size_mean_hist, 1L)))
cat(sprintf("local_loglik_prefit=%d\n", local_prefit_loglik))
cat(sprintf("local_loglik_nested_total=%d\n", eval_counts["total"]))
cat(sprintf("local_loglik_nested_init=%d\n", eval_counts["initialization"]))
cat(sprintf("local_loglik_nested_enrichment=%d\n", eval_counts["enrichment"]))
cat(sprintf("local_loglik_nested_population_refresh=%d\n", eval_counts["population_refreshes"]))
cat(sprintf("local_loglik_end_to_end=%d\n", local_end_to_end))
for (j in seq_len(P)) {
  cat(sprintf("%s_mean_stan=%.6f %s_mean_nested=%.6f\n", mu_names[j], stan_mu_mean[j], mu_names[j], nested_mu_mean[j]))
}
for (j in seq_len(P)) {
  cat(sprintf("%s_mean_stan=%.6f %s_mean_nested=%.6f\n", sigma2_names[j], stan_sigma2_mean[j], sigma2_names[j], nested_sigma2_mean[j]))
}
