rm(list = ls())
source("smc_core.R")
source("SMC_super_fast.R")
source("new_SMC_cache.R")
source("outer_SMC.R")
source("make_prior.R")

stan_results_file <- file.path("samples", "normal_normal_stan_results.rds")
if (!file.exists(stan_results_file)) {
  stop("Missing Stan results. Run `Rscript run_stan_normal_normal.R` first.")
}

stan_res <- readRDS(stan_results_file)

S <- stan_res$simulation$S
N <- stan_res$simulation$N
sigma_y <- stan_res$simulation$sigma_y
mu_true <- stan_res$simulation$mu_true
tau2_true <- stan_res$simulation$tau2_true
y <- stan_res$data$y
data_list <- lapply(seq_len(S), function(i) as.numeric(y[i, ]))
prior_m0 <- as.numeric(stan_res$prior$m0)
prior_s0 <- as.numeric(stan_res$prior$s0)
prior_a0 <- as.numeric(stan_res$prior$a0)
prior_b0 <- as.numeric(stan_res$prior$b0)

outer_particles <- 1000L
outer_rounds <- 200L
inner_particles <- 1000L
inner_subject_particles <- 64L
inner_rounds <- 50L

mc_cores <- 1 #max(1L, min(S, parallel::detectCores(logical = TRUE) - 1L))

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

rtheta_given_phi <- function(phi, n, aux = NULL) {
  tau2 <- exp(as.numeric(phi[2L]))
  matrix(rnorm(as.integer(n), mean = as.numeric(phi[1L]), sd = sqrt(tau2)), ncol = 1L)
}

# 1. Run the outers.
alpha_ref_mean <- setNames(prior_m0, "alpha")
alpha_ref_var <- prior_s0 + prior_b0 / (prior_a0 - 1.0)
Sigma_ref <- matrix(alpha_ref_var, nrow = 1L, ncol = 1L, dimnames = list("alpha", "alpha"))
Sigma_ref_inv <- chol2inv(chol(Sigma_ref))
Sigma_ref_logdet <- as.numeric(determinant(Sigma_ref, logarithm = TRUE)$modulus)

outer_subject <- parallel::mclapply(
  seq_len(S),
  function(i) {
    out <- enhanced_smc_elite(
      data = data_list[[i]],
      loglik_fn = loglik_fn,
      mu_ref = alpha_ref_mean,
      Sigma_ref = Sigma_ref,
      M = outer_particles,
      max_rounds = outer_rounds,
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

surrogates <- lapply(
  seq_len(S),
  function(i) {
    build_subject_surrogate_from_outer(
      smc_out = outer_subject[[i]],
      data = data_list[[i]],
      subj_id = i,
      loglik_fn = loglik_fn,
      base_seed = as.integer(40000L + 1009L * i)
    )
  }
)

# 2. Run the inner.
prior <- make_prior_phi_diag(
  m0 = prior_m0,
  s0 = prior_s0,
  a = prior_a0,
  b = prior_b0,
  d = 1L
)

inner <- outer_smc_phi_batch(
  caches = surrogates,
  rprior_phi = prior$rprior,
  logprior_phi = prior$lprior,
  gaussian_map_fn = gaussian_map_fn,
  rtheta_given_phi = rtheta_given_phi,
  data_list = data_list,
  loglik_fn = loglik_fn,
  M = inner_particles,
  pm_mode = "strict",
  inner_ll_mode = "adaptive_pm",
  adaptive_pm_control = list(
    M_default = inner_subject_particles,
    M_levels = inner_subject_particles,
    da_M = 16L
  ),
  max_rounds = inner_rounds,
  seed = 123L,
  verbose = TRUE
)

# 3. Load the Stan draws.
stan_mu <- stan_res$stan$mu
stan_tau2 <- stan_res$stan$tau2

# 4. Make a comparison plot.
smc_phi <- as.matrix(inner$phi)
smc_w <- pmax(as.numeric(inner$w), 0)
smc_w <- smc_w / sum(smc_w)
smc_mu <- smc_phi[, 1L]
smc_tau2 <- exp(smc_phi[, 2L])

set.seed(123L)
smc_mu_draws <- sample(smc_mu, size = 4000L, replace = TRUE, prob = smc_w)
smc_tau2_draws <- sample(smc_tau2, size = 4000L, replace = TRUE, prob = smc_w)

# plot_file <- file.path("samples", "normal_normal_benchmark_posteriors.png")
# png(plot_file, width = 900, height = 420)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
plot(density(stan_mu), main = "Posterior of mu", xlab = "mu", ylab = "Density", lwd = 2, col = "black")
lines(density(smc_mu_draws), lwd = 2, col = "red3")
legend("topright", legend = c("Stan", "Framework"), col = c("black", "red3"), lwd = 2, bty = "n")
plot(density(stan_tau2), main = "Posterior of tau2", xlab = "tau2", ylab = "Density", lwd = 2, col = "black")
lines(density(smc_tau2_draws), lwd = 2, col = "red3")
legend("topright", legend = c("Stan", "Framework"), col = c("black", "red3"), lwd = 2, bty = "n")
# dev.off()

results <- list(
  stan_results_file = stan_results_file,
  outer_subject = outer_subject,
  inner = inner,
  plot_file = plot_file,
  settings = list(
    outer_particles = outer_particles,
    outer_rounds = outer_rounds,
    inner_particles = inner_particles,
    inner_subject_particles = inner_subject_particles,
    inner_rounds = inner_rounds,
    mc_cores = mc_cores
  )
)

outfile <- file.path("samples", "normal_normal_benchmark_results.rds")
saveRDS(results, outfile)

cat("NORMAL_NORMAL_BENCHMARK_OK\n")
cat(sprintf("stan_results=%s\n", stan_results_file))
cat(sprintf("results=%s\n", outfile))
cat(sprintf("plot=%s\n", plot_file))
cat(sprintf("inner_final_lambda=%.4f\n", tail(inner$meta$lambda_hist, 1L)))
cat(sprintf("mu_mean_stan=%.6f mu_mean_framework=%.6f\n", mean(stan_mu), sum(smc_mu * smc_w)))
cat(sprintf("tau2_mean_stan=%.6f tau2_mean_framework=%.6f\n", mean(stan_tau2), sum(smc_tau2 * smc_w)))
