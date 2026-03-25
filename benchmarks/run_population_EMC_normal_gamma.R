rm(list = ls())

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[1L]))
} else {
  normalizePath("benchmarks/run_population_EMC_normal_gamma.R")
}
repo_dir <- dirname(dirname(script_path))
setwd(repo_dir)

suppressPackageStartupMessages({
  library(EMC2)
})

set.seed(123L)

source("hierarchical_locals.R")
source("population_models.R")
source("outer_population_smc.R")
source("utilities.R")

data_path <- file.path("benchmarks", "samples", "full_EMC2.RData")
results_path <- file.path("benchmarks", "results", "population_emc_normal_gamma_results.rds")
posterior_plot_path <- file.path("benchmarks", "results", "population_emc_normal_gamma_posteriors.pdf")

pilot_size <- 20L
pilot_particles <- 1000L
local_particles <- 5000L
outer_particles <- 4000L
n_jobs <- 12

dir.create(file.path("benchmarks", "results"), showWarnings = FALSE, recursive = TRUE)

if (!file.exists(data_path)) {
  stop("Missing data file: ", data_path)
}

load(data_path)
if (!exists("ELP_DDM", inherits = FALSE)) {
  stop("The data file does not define an object named 'ELP_DDM'.")
}

emc <- ELP_DDM[[1L]]
data_list <- emc$data
model_factory <- emc$model
param_names <- emc$par_names
d <- length(param_names)

base_mu <- stats::setNames(rep(0, d), param_names)
base_Sigma <- diag(1, nrow = d, ncol = d)
rownames(base_Sigma) <- colnames(base_Sigma) <- param_names

population_model <- make_population_model_diag_gaussian(
  alpha_names = param_names,
  mean_prior_mean = rep(0, d),
  mean_prior_var = rep(1, d),
  sigma2_prior_shape = rep(2, d),
  sigma2_prior_rate = rep(0.3, d),
  label = "emc_normal_gamma"
)

loglik_emc2 <- function(Theta, data_i) {
  Theta <- as.matrix(Theta)
  colnames(Theta) <- param_names
  out <- as.numeric(EMC2:::calc_ll_manager(Theta, data_i, model_factory, r_cores = 1L))
  out[!is.finite(out) | out > 100] <- min(out)
  out
}

if(!file.exists("benchmarks/samples/stage_ELP.RData")){
  stage <- prepare_reference_local_stage(
    data_list = data_list,
    loglik_fn = loglik_emc2,
    base_mu = base_mu,
    base_Sigma = base_Sigma,
    broad_scale = 2,
    n_jobs = n_jobs,
    pilot_population_model = population_model
  )
  save(stage, file = "benchmarks/samples/stage_ELP.RData")
} else{
  load("benchmarks/samples/stage_ELP.RData")
}

factor_set <- build_population_factor_set(stage$local_objects, population_model)

fit <- outer_population_smc(
  factor_set = factor_set,
  N = outer_particles,
  n_mcmc_moves = 2L,
  max_rounds = 100L,
  n_cores = n_jobs,
  seed = 123L,
  verbose = TRUE
)

posterior_summary <- summarize_population_posterior_diag(
  theta = fit$theta,
  w = fit$w,
  model = population_model
)

workflow_draws <- smc_posteriors(
  fit,
  population_model = population_model,
  seed = 124L
)
mu_draws <- workflow_draws$mu
sigma2_draws <- workflow_draws$sigma2

mu_EMC_draws <- parameters(ELP_DDM, selection = "mu")
sigma_EMC_draws <- parameters(ELP_DDM, selection = "sigma2")
colnames(mu_EMC_draws) <- colnames(mu_draws)
colnames(sigma_EMC_draws) <- colnames(sigma2_draws)

# grDevices::pdf(posterior_plot_path, width = 12, height = 7)
plot_posteriors(mu_draws, mu_EMC_draws, n_cols = 4L, main_prefix = "mu")
graphics::mtext("Group-level means", outer = TRUE, line = 0.5, cex = 1.1)
plot_posteriors(sigma2_draws, sigma_EMC_draws, n_cols = 4L, main_prefix = "sigma2")
graphics::mtext("Group-level variances", outer = TRUE, line = 0.5, cex = 1.1)
# grDevices::dev.off()
# 
# results <- list(
#   fit = fit,
#   posterior_summary = posterior_summary,
#   population_model = population_model,
#   local_objects = stage$local_objects,
#   broad_reference = stage$broad_reference,
#   refined_reference = stage$refined_reference,
#   pilot = stage$pilot,
#   settings = list(
#     data_path = normalizePath(data_path, winslash = "/", mustWork = TRUE),
#     results_path = results_path,
#     posterior_plot_path = posterior_plot_path,
#     pilot_size = pilot_size,
#     pilot_particles = pilot_particles,
#     local_particles = local_particles,
#     outer_particles = outer_particles,
#     n_jobs = n_jobs
#   )
# )
# 
# saveRDS(results, results_path)
# 
# 
# 
# cat("Saved results to:", results_path, "\n")
# cat("Saved posterior plots to:", posterior_plot_path, "\n")
# cat("Locals:", length(data_list), "\n")
# cat("Pilot locals:", length(stage$pilot$selection$indices), "\n")
# cat("Outer rounds:", fit$meta$rounds, "\n")
# cat("Final beta:", sprintf("%.4f", fit$beta), "\n")
# cat("Log evidence:", sprintf("%.4f +/- %.4f", fit$log_evidence, fit$mcse_log_evidence), "\n")
# 
# cat("\nGroup-level means:\n")
# print(posterior_summary$mu, row.names = FALSE)
# 
# cat("\nGroup-level variances:\n")
# print(posterior_summary$sigma2, row.names = FALSE)
