#!/usr/bin/env Rscript

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[1L]))
} else {
  normalizePath("dev/scratch/profile_nested_population_benchmark.R")
}
repo_dir <- dirname(dirname(dirname(script_path)))
setwd(repo_dir)

source("smc_core.R")
source("SMC_super_fast.R")
source("new_SMC_cache.R")
source("nested_population_SMC.R")
source("make_prior.R")

dir.create("dev/scratch", showWarnings = FALSE, recursive = TRUE)

stan_results_file <- file.path("samples", "normal_normal_stan_results.rds")
if (!file.exists(stan_results_file)) {
  stop("Missing Stan results. Run `Rscript run_stan_normal_normal.R` first.")
}

stan_res <- readRDS(stan_results_file)

S <- stan_res$simulation$S
sigma_y <- stan_res$simulation$sigma_y
y <- stan_res$data$y
data_list <- lapply(seq_len(S), function(i) as.numeric(y[i, ]))

prior_m0 <- as.numeric(stan_res$prior$m0)
prior_s0 <- as.numeric(stan_res$prior$s0)
prior_a0 <- as.numeric(stan_res$prior$a0)
prior_b0 <- as.numeric(stan_res$prior$b0)

local_particles <- 1000L
local_rounds <- 200L
population_particles <- 1000L
local_bank_particles <- 64L
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

alpha_ref_mean <- setNames(prior_m0, "alpha")
alpha_ref_var <- prior_s0 + prior_b0 / (prior_a0 - 1.0)
Sigma_ref <- matrix(alpha_ref_var, nrow = 1L, ncol = 1L, dimnames = list("alpha", "alpha"))
Sigma_ref_inv <- chol2inv(chol(Sigma_ref))
Sigma_ref_logdet <- as.numeric(determinant(Sigma_ref, logarithm = TRUE)$modulus)

cat("Building benchmark local prefits...\n")
ll_eval_counter_reset()
prefit_time <- system.time({
  local_fits <- parallel::mclapply(
    seq_len(S),
    function(i) {
      out <- enhanced_smc_elite(
        data = data_list[[i]],
        loglik_fn = loglik_fn,
        mu_ref = alpha_ref_mean,
        Sigma_ref = Sigma_ref,
        M = local_particles,
        max_rounds = local_rounds,
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
})
local_prefit_loglik <- ll_eval_counter_get()
ll_eval_counter_disable()

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
  m0 = prior_m0,
  s0 = prior_s0,
  a = prior_a0,
  b = prior_b0,
  d = 1L
)

profile_file <- file.path("dev", "scratch", "nested_population_benchmark_Rprof.out")
summary_file <- file.path("dev", "scratch", "nested_population_benchmark_profile.txt")

if (file.exists(profile_file)) unlink(profile_file)
if (file.exists(summary_file)) unlink(summary_file)

cat("Profiling nested population benchmark run...\n")
Rprof(profile_file, interval = 0.01)
nested_time <- system.time({
  fit <- nested_population_smc(
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
    verbose = FALSE
  )
})
Rprof(NULL)

prof <- summaryRprof(profile_file)

lines_out <- c(
  "PROFILE_NESTED_POPULATION_BENCHMARK",
  sprintf("stan_results=%s", stan_results_file),
  sprintf("groups=%d", S),
  sprintf("local_particles=%d", local_particles),
  sprintf("population_particles=%d", population_particles),
  sprintf("local_bank_particles=%d", local_bank_particles),
  sprintf("prefit_elapsed=%.3f", unname(prefit_time["elapsed"])),
  sprintf("prefit_local_loglik=%d", as.integer(local_prefit_loglik)),
  sprintf("nested_elapsed=%.3f", unname(nested_time["elapsed"])),
  sprintf("population_rounds=%d", fit$meta$rounds),
  sprintf("bank_topup_points=%d", fit$meta$bank_topup_points),
  sprintf("nested_local_loglik_total=%d", fit$meta$local_loglik_evals["total"]),
  sprintf("bank_mean_size_final=%.1f", tail(fit$meta$bank_size_mean_hist, 1L)),
  "",
  "BY_SELF",
  capture.output(print(head(prof$by.self, 20L))),
  "",
  "BY_TOTAL",
  capture.output(print(head(prof$by.total, 20L)))
)

writeLines(lines_out, con = summary_file)
cat(paste(lines_out, collapse = "\n"))
cat("\n")
