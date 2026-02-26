rm(list = ls())
library(EMC2)
library(parallel)
set.seed(123)

load("single_DDM30.RData")
emc <- single[[1]]

# ---------- prior ----------
mu_ref <- emc[[1]]$prior$theta_mu_mean
Sigma_ref <- emc[[1]]$prior$theta_mu_var
data <- lapply(single, function(x) return(x[[1]]$data[[1]]))
loglik_fn <- emc[[1]]$model

source("SMC_super_fast.R")
source("smc_diagnostics.R")

smc <- mclapply(data, enhanced_smc_elite, loglik_fn = loglik_fn, mu_ref = mu_ref, Sigma_ref = Sigma_ref,
                seed = 123, M = 5000, mc.cores = 12)


for(i in 1:length(smc)){
  pars <- parameters(single[[i]], selection = "alpha")
  pars <- pars[,-1]
  recovery_stats <- diagnose_recovery(
    true_values = pars,
    smc_res = smc[[i]],
    save_plot = FALSE,
    density_method = "smooth"  # Use smooth density estimation
  )
}

source("phi_cache.R")
source("make_prior.R")

d <- ncol(smc[[1]]$Theta)

# dimension-wise hyperparameters (scalars broadcast or give length-d vectors)
m0 <- 0         # or rep(0, d)         # prior means for μ_j
s0 <- 1        # or rep(10, d)        # prior variances for μ_j
a  <- 2         # or rep(2, d)         # InvGamma shape for σ²_j
b  <- .3         # or rep(2, d)         # InvGamma rate  for σ²_j

# build prior and Gaussian map
pri_diag <- make_prior_phi_diag(m0, s0, a, b, d)
gauss_map_fn_diag <- phi_to_gaussian_params_diag_factory()

caches <- build_phi_aware_caches_from_base_proxy(
  smc_list = smc,
  data_list = data,
  subjects = 1:length(smc),
  rprior_phi = pri_diag$rprior,
  logprior_phi = pri_diag$lprior,
  loglik_fn = loglik_fn,
  mc.cores = 12,Npilot = 4000, Nselect = 2000
)

caches_ppmh <- make_pmmh_cache(caches$caches_phi_aware)

source("PMMH_RW_works.R")
res <- run_pmmh(caches_ppmh,
                n_iter = 20000L,
                burn    = 1000L,
                thin    = 25L,
                m_use   = 100L,
                rho_sub = 0.999)

matplot(res$mu_draws, type = "l")

mu_draws <- res$mu_draws
# ---------------------------- Run outer SMC ------------------------------
# Suppose you already built per-subject caches:
#   caches <- mcmapply(smc_list, data_list, FUN = build_subject_cache_from_smc, ...)
# with each cache having fields: Theta, Z, log_py, log_r, meta$d, etc.
source("outer_SMC.R")


set.seed(1)
out <- outer_smc_phi_batch(
  caches = caches$caches_phi_aware,                      # list of length S
  rprior_phi = pri_diag$rprior,
  logprior_phi = pri_diag$lprior,
  gaussian_map_fn = gauss_map_fn_diag,       # fast θ|φ evaluator
  data_list = data,
  loglik_fn = loglik_fn,
  # use_independence_moves = TRUE,
  # auto_enrich_M_add_total = 128,
  # use_control_variate = TRUE,
  M = 5000,
  cess_target = 0.95,
  auto_enrich_enable = FALSE,
  refresh_batches_after_resample = FALSE,
  resample_threshold = 0.6,
  n_moves = 2,
  rw_scale_init = 1.0,
  max_rounds = 200,
  seed = 42,
  verbose = TRUE,
)

source("outer_SMC_check.R")
phi_draws <- posterior_draws_phi(out, n = 10000, method = "systematic")

d <- ncol(phi_draws) / 2
mu_draws     <- phi_draws[, 1:d, drop = FALSE]
sigma2_draws <- exp(phi_draws[, d + seq_len(d), drop = FALSE])



load("full_DDM30.RData")
credint(full)

par_names <- names(sampled_pars(full))
par(mfrow = c(3,2))

pars <- parameters(full)

for(i in 1:ncol(mu_draws)){
  plot(density(mu_draws[,i]), main = par_names[i])
  lines(density(pars[,i]), col = "red")
}

pars <- parameters(full, selection = "sigma2")

for(i in 1:ncol(mu_draws)){
  plot(density(sigma2_draws[,i]), main = par_names[i])
  lines(density(pars[,i]), col = "red")
}


#
#  # 3a) Unbiased estimator at φ (Gaussian prior fast path)
# phi_to_gauss <- function(phi) {
#   S <- Sigma_of_phi(phi); K <- solve(S)
#   list(mu = mu_of_phi(phi), Sigma_inv = K,
#        logdet = as.numeric(determinant(S, TRUE)$modulus),
#        const = 0.0)
# }
# cache_i <- register_gaussian_prior_map(cache_i, phi_to_gauss)
# log_p_hat <- log_marginal_unbiased_gaussian(cache_i, phi)
#
# # 3b) DA surrogate (SNIS from base SMC)
# log_p_tilde <- log_marginal_proxy_from_base(cache_i, phi,
#                                             log_prior_theta_given_phi_mat = function(Theta, phi, aux) .log_prior_gauss_mat(Theta, phi_to_gauss(phi))
# )
#
# # Optional PSIS-k to monitor proxy stability
# k_base <- pareto_k_from_base_proxy(cache_i, phi,
#                                    log_prior_theta_given_phi_mat = function(Theta, phi, aux) .log_prior_gauss_mat(Theta, phi_to_gauss(phi))
# )
