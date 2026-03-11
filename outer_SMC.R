#!/usr/bin/env Rscript

# Legacy benchmark compatibility shim.
# The maintained Phase 1 implementation lives in nested_population_SMC.R.

source("nested_population_SMC.R")

if (!exists("plot_file", inherits = FALSE)) {
  plot_file <- NA_character_
}

outer_smc_phi_batch <- function(caches,
                                rprior_phi,
                                logprior_phi,
                                gaussian_map_fn,
                                rtheta_given_phi = NULL,
                                data_list = NULL,
                                loglik_fn = NULL,
                                M = 500L,
                                pm_mode = NULL,
                                inner_ll_mode = NULL,
                                adaptive_pm_control = list(),
                                max_rounds = 100L,
                                seed = 123L,
                                verbose = TRUE,
                                ...) {
  M_local <- as.integer(
    adaptive_pm_control$M_default %||%
      adaptive_pm_control$M_levels[1L] %||%
      64L
  )
  nested_population_smc(
    rprior_phi = rprior_phi,
    local_objs = caches,
    logprior_phi = logprior_phi,
    gaussian_map_fn = gaussian_map_fn,
    N = as.integer(M),
    M_local = M_local,
    rho_step = 0.4,
    rho_res = 0.5,
    rho_local = 0.5,
    n_population_moves = 1L,
    n_local_moves = 1L,
    n_population_refresh_moves = 0L,
    max_rounds = as.integer(max_rounds),
    seed = as.integer(seed),
    verbose = verbose,
    ...
  )
}

