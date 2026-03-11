#!/usr/bin/env Rscript

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[1L]))
} else {
  normalizePath("benchmarks/nested_normal_normal.R")
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

stan_results_file <- file.path("samples", "normal_normal_stan_results.rds")
if (!file.exists(stan_results_file)) {
  stop("Missing Stan results. Run `Rscript run_stan_normal_normal.R` first.")
}

stan_res <- readRDS(stan_results_file)

S <- stan_res$simulation$S
N_obs <- stan_res$simulation$N
sigma_y <- stan_res$simulation$sigma_y
mu_true <- stan_res$simulation$mu_true
tau2_true <- stan_res$simulation$tau2_true
y <- stan_res$data$y
data_list <- lapply(seq_len(S), function(i) as.numeric(y[i, ]))

local_particles <- 1000L
local_rounds <- 200L
population_particles <- 1000L
local_module_particles <- 64L
population_rounds <- 50L
mc_cores <- 1L

loglik_fn <- function(Theta, y_i) {
  Theta <- as.matrix(Theta)
  alpha <- Theta[, 1L]
  vapply(alpha, function(a) sum(dnorm(y_i, mean = a, sd = sigma_y, log = TRUE)), numeric(1L))
}

gaussian_map_fn <- function(phi, d = 1L) {
  if (!is.null(d) && as.integer(d) != 1L) stop("This benchmark expects d = 1.")
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

local_fits <- parallel::mclapply(
  seq_len(S),
  function(i) {
    out <- enhanced_smc_elite(
      data = data_list[[i]],
      loglik_fn = loglik_fn,
      mu_ref = mu_ref,
      Sigma_ref = Sigma_ref,
      M = local_particles,
      max_rounds = local_rounds,
      seed = as.integer(1123L + 31L * i),
      verbose = FALSE
    )
    out$working_prior <- list(
      mu = as.numeric(mu_ref),
      Sigma = Sigma_ref,
      Sigma_inv = Sigma_ref_inv,
      logdet = Sigma_ref_logdet
    )
    out
  },
  mc.cores = mc_cores
)

local_objs <- lapply(
  seq_len(S),
  function(i) {
    build_local_exact_object(
      smc_out = local_fits[[i]],
      data = data_list[[i]],
      subj_id = i,
      loglik_fn = loglik_fn,
      base_seed = as.integer(40000L + 1009L * i)
    )
  }
)

prior <- make_prior_phi_diag(
  m0 = stan_res$prior$m0,
  s0 = stan_res$prior$s0,
  a = stan_res$prior$a0,
  b = stan_res$prior$b0,
  d = 1L
)

nested_fit <- nested_population_smc(
  rprior_phi = prior$rprior,
  local_objs = local_objs,
  logprior_phi = prior$lprior,
  gaussian_map_fn = gaussian_map_fn,
  N = population_particles,
  M_local = local_module_particles,
  rho_step = 0.4,
  rho_res = 0.5,
  rho_local = 0.5,
  n_population_moves = 1L,
  n_local_moves = 1L,
  n_population_refresh_moves = 0L,
  max_rounds = population_rounds,
  seed = 123L,
  verbose = TRUE
)

stan_mu <- stan_res$stan$mu
stan_tau2 <- stan_res$stan$tau2

nested_phi <- as.matrix(nested_fit$phi)
nested_w <- pmax(as.numeric(nested_fit$w), 0)
nested_w <- nested_w / sum(nested_w)
nested_mu <- nested_phi[, 1L]
nested_tau2 <- exp(nested_phi[, 2L])

set.seed(123L)
nested_mu_draws <- sample(nested_mu, size = 4000L, replace = TRUE, prob = nested_w)
nested_tau2_draws <- sample(nested_tau2, size = 4000L, replace = TRUE, prob = nested_w)

plot_file <- file.path("benchmarks", "nested_normal_normal_posteriors.png")
png(plot_file, width = 900, height = 420)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
plot(density(stan_mu), main = "Posterior of mu", xlab = "mu", ylab = "Density", lwd = 2, col = "black")
lines(density(nested_mu_draws), lwd = 2, col = "red3")
legend("topright", legend = c("Stan", "Nested SMC"), col = c("black", "red3"), lwd = 2, bty = "n")
plot(density(stan_tau2), main = "Posterior of tau2", xlab = "tau2", ylab = "Density", lwd = 2, col = "black")
lines(density(nested_tau2_draws), lwd = 2, col = "red3")
legend("topright", legend = c("Stan", "Nested SMC"), col = c("black", "red3"), lwd = 2, bty = "n")
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
    local_module_particles = local_module_particles,
    population_rounds = population_rounds,
    mc_cores = mc_cores
  )
)

outfile <- file.path("samples", "nested_normal_normal_results.rds")
saveRDS(results, outfile)

eval_counts <- nested_fit$meta$local_loglik_evals

cat("NESTED_NORMAL_NORMAL_OK\n")
cat(sprintf("stan_results=%s\n", stan_results_file))
cat(sprintf("results=%s\n", outfile))
cat(sprintf("plot=%s\n", plot_file))
cat(sprintf("population_final_lambda=%.4f\n", tail(nested_fit$meta$lambda_hist, 1L)))
cat(sprintf("population_rounds=%d\n", nested_fit$meta$rounds))
cat(sprintf("local_loglik_total=%d\n", eval_counts["total"]))
cat(sprintf("local_loglik_init=%d\n", eval_counts["initialization"]))
cat(sprintf("local_loglik_bridge=%d\n", eval_counts["bridge_updates"]))
cat(sprintf("local_loglik_rejuvenation=%d\n", eval_counts["local_rejuvenation"]))
cat(sprintf("local_loglik_population_refresh=%d\n", eval_counts["population_refreshes"]))
cat(sprintf("mu_mean_stan=%.6f mu_mean_nested=%.6f\n", mean(stan_mu), sum(nested_mu * nested_w)))
cat(sprintf("tau2_mean_stan=%.6f tau2_mean_nested=%.6f\n", mean(stan_tau2), sum(nested_tau2 * nested_w)))
