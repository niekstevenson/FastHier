rm(list = ls())
library(EMC2)
source("SMC_super_fast.R")
source("PMMH.R")
source("build_cache.R")

set.seed(1)
load("~/Documents/2025/TwoStep/samples/single_DDM12.RData")

emc <- single[[1]]
mu_ref <- emc[[1]]$prior$theta_mu_mean
Sigma_ref <- emc[[1]]$prior$theta_mu_var
ll <- emc[[1]]$model

data <- lapply(single, function(x) return(x[[1]]$data[[1]]))


library(parallel)
caches<- mclapply(data, FUN = build_inner_cache, verbose = FALSE, M = 5000,
                    loglik_fn = ll, mu0 = mu_ref, Sigma0 = Sigma_ref,
                    mc.cores = 12)

refresh_fn <- function(i, mu, tau2, old_cache) {
  Sigma_diag <- diag(as.numeric(tau2)); colnames(Sigma_diag) <- rownames(Sigma_diag) <- colnames(old_cache$theta)
  build_inner_cache(
    data = old_cache$meta$data %||% stop("store data handle in cache$meta$data for refresh"),
    loglik_fn = old_cache$meta$loglik_fn %||% stop("store loglik_fn in cache$meta$loglik_fn"),
    mu0 = setNames(mu, colnames(old_cache$theta)),
    Sigma0 = Sigma_diag,
    M = old_cache$meta$M_in,
    target_size = nrow(old_cache$theta), # keep size
    dedup = TRUE,
    store_transport = FALSE,
    store_mixZ = FALSE,
    verbose = FALSE
  )
}

source("PMMH.R")

d <- length(mu_ref)
fit <- pmmh_outer_cpm_da(
  caches,
  m0 = rep(0, d), s0 = rep(1, d),
  a0 = rep(2, d), b0 = rep(.3, d),
  n_iter = 4000,
  K1 = 16, K2 = 128, rho = 0.995,
  refresh_fn = refresh_fn,
  verbose = TRUE
)

mu <- fit$mu
colnames(fit$mu) <- names(mu_ref)
mu <- mu[2000:nrow(mu),]
load("~/Documents/2025/TwoStep/samples/full_DDM12.RData")
plot(full, selection = "mu")
matplot(mu, type = "l")
colnames(mu) <- names(mu_ref)
t(apply(mu, 2, quantile, probs = c(0.025, .5, .975)))
credint(full)

source("smc_diagnostics.R")
pars <- parameters(single[[1]], selection = "alpha")[,-1]
scm_res <-
diagnose_recovery(pars, caches[[1]])

source("PMMH_health.R")
check_outer_consistency_once(fit$final_state$mu, fit$final_state$log_tau2, caches)

