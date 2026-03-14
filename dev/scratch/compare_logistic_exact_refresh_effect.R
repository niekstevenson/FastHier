#!/usr/bin/env Rscript

source("smc_core.R")
source("SMC_super_fast.R")
source("new_SMC_cache.R")
source("nested_population_SMC.R")
source("make_prior.R")

repo_dir <- getwd()

exact_file <- file.path(repo_dir, "samples", "nested_hierarchical_logistic_results.rds")
stan_file <- file.path(repo_dir, "samples", "stan_hierarchical_logistic_results.rds")
plot_file <- file.path(repo_dir, "benchmarks", "nested_hierarchical_logistic_effect.png")
cached_file <- file.path(repo_dir, "samples", "nested_hierarchical_logistic_cached_only_results.rds")

exact_res <- readRDS(exact_file)
stan_res <- readRDS(stan_file)

theta_names <- colnames(stan_res$truth$alpha)
mu_names <- paste0("mu_", theta_names)
sigma2_names <- paste0("sigma2_", theta_names)
S <- dim(stan_res$data$X)[1L]
P <- dim(stan_res$data$X)[3L]
N_per_group <- dim(stan_res$data$X)[2L]

data_list <- lapply(
  seq_len(S),
  function(i) {
    list(
      X = stan_res$data$X[i, , , drop = FALSE][1L, , , drop = TRUE],
      y = as.numeric(stan_res$data$y[i, ])
    )
  }
)

loglik_hier_logistic <- function(Theta, data_i) {
  Theta <- as.matrix(Theta)
  X <- as.matrix(data_i$X)
  y <- as.numeric(data_i$y)
  eta <- Theta %*% t(X)
  softplus <- pmax(eta, 0) + log1p(exp(-abs(eta)))
  as.numeric(drop(eta %*% y) - rowSums(softplus))
}

local_objs <- lapply(
  seq_len(S),
  function(i) {
    build_local_exact_object(
      smc_out = exact_res$local_fits[[i]],
      data = data_list[[i]],
      subj_id = i,
      loglik_fn = loglik_hier_logistic,
      base_seed = as.integer(40000L + 1009L * i)
    )
  }
)

prior <- make_prior_phi_diag(
  m0 = stan_res$priors$m0,
  s0 = stan_res$priors$s0,
  a = stan_res$priors$a0,
  b = stan_res$priors$b0,
  d = P
)

cached_fit <- nested_population_smc(
  rprior_phi = prior$rprior,
  local_objs = local_objs,
  logprior_phi = prior$lprior,
  gaussian_map_fn = phi_to_gaussian_params_diag_factory(),
  N = as.integer(exact_res$settings$population_particles),
  M_local = as.integer(exact_res$settings$local_bank_particles),
  rho_res = 0.5,
  rho_local = 0.5,
  n_population_moves = 1L,
  n_local_moves = 0L,
  n_population_refresh_moves = 0L,
  max_rounds = as.integer(exact_res$settings$population_rounds),
  max_bank_topups = as.integer(exact_res$settings$max_bank_topups),
  bank_split_tol = as.numeric(exact_res$settings$bank_split_tol),
  enable_exact_local_refresh = FALSE,
  seed = 123L,
  verbose = FALSE
)

saveRDS(
  list(
    stan_results_file = stan_file,
    local_fits = exact_res$local_fits,
    nested_fit = cached_fit,
    settings = modifyList(exact_res$settings, list(enable_exact_local_refresh = FALSE))
  ),
  cached_file
)

.weighted_draws <- function(fit, P, n_draw) {
  phi <- as.matrix(fit$phi)
  w <- pmax(as.numeric(fit$w), 0)
  w <- w / sum(w)
  mu <- phi[, seq_len(P), drop = FALSE]
  sigma2 <- exp(phi[, P + seq_len(P), drop = FALSE])
  draw_idx <- sample.int(nrow(phi), size = n_draw, replace = TRUE, prob = w)
  list(
    mu = mu[draw_idx, , drop = FALSE],
    sigma2 = sigma2[draw_idx, , drop = FALSE],
    mean_mu = colSums(mu * w),
    mean_sigma2 = colSums(sigma2 * w)
  )
}

set.seed(123L)
stan_mu <- as.matrix(stan_res$draws$mu)
stan_sigma2 <- as.matrix(stan_res$draws$sigma2)
colnames(stan_mu) <- mu_names
colnames(stan_sigma2) <- sigma2_names

cached_draws <- .weighted_draws(cached_fit, P = P, n_draw = nrow(stan_mu))
exact_draws <- .weighted_draws(exact_res$nested_fit, P = P, n_draw = nrow(stan_mu))
colnames(cached_draws$mu) <- mu_names
colnames(cached_draws$sigma2) <- sigma2_names
colnames(exact_draws$mu) <- mu_names
colnames(exact_draws$sigma2) <- sigma2_names

png(plot_file, width = 1400, height = 900, res = 120)
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))

for (j in seq_len(P)) {
  plot(density(stan_mu[, j]), main = mu_names[j], xlab = mu_names[j], ylab = "Density", lwd = 2, col = "black")
  lines(density(cached_draws$mu[, j]), lwd = 2, col = "grey40")
  lines(density(exact_draws$mu[, j]), lwd = 2, col = "red3")
  abline(v = stan_res$truth$mu[j], col = "grey70", lty = 3)
  legend("topright",
         legend = c("Stan", "Cached-only", "Exact refresh", "Truth"),
         col = c("black", "grey40", "red3", "grey70"),
         lwd = c(2, 2, 2, 1),
         lty = c(1, 1, 1, 3),
         bty = "n")
}

for (j in seq_len(P)) {
  plot(density(stan_sigma2[, j]), main = sigma2_names[j], xlab = sigma2_names[j], ylab = "Density", lwd = 2, col = "black")
  lines(density(cached_draws$sigma2[, j]), lwd = 2, col = "grey40")
  lines(density(exact_draws$sigma2[, j]), lwd = 2, col = "red3")
  abline(v = stan_res$truth$sigma2[j], col = "grey70", lty = 3)
  legend("topright",
         legend = c("Stan", "Cached-only", "Exact refresh", "Truth"),
         col = c("black", "grey40", "red3", "grey70"),
         lwd = c(2, 2, 2, 1),
         lty = c(1, 1, 1, 3),
         bty = "n")
}

dev.off()

cat("COMPARE_LOGISTIC_EFFECT_OK\n")
cat(sprintf("cached_results=%s\n", cached_file))
cat(sprintf("plot=%s\n", plot_file))
cat(sprintf("cached_only_nested_total=%d\n", cached_fit$meta$local_loglik_evals["total"]))
cat(sprintf("exact_refresh_nested_total=%d\n", exact_res$nested_fit$meta$local_loglik_evals["total"]))
