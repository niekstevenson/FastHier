rm(list = ls())
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[1L]))
} else {
  normalizePath("benchmarks/run_local_EMC.R")
}
repo_dir <- dirname(dirname(script_path))
setwd(repo_dir)

library(EMC2)
library(parallel)
set.seed(123)

# Change this if the EMC2 example data lives somewhere else.
data_path <- file.path("benchmarks", "samples", "single_DDM30.RData")

# Main run settings.
mc.cores <- 12L
M <- 4000L
base_seed <- 123L
verbose <- TRUE

# Output files.
results_file <- file.path("benchmarks", "results", "local_emc_results.rds")
posterior_plot_file <- file.path("benchmarks", "results", "local_emc_posteriors.pdf")

dir.create(file.path("benchmarks", "samples"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path("benchmarks", "results"), showWarnings = FALSE, recursive = TRUE)

source("SMC_super_fast.R")
source("utilities.R")

if (!file.exists(data_path)) {
  stop("Missing data file: ", data_path)
}

load(data_path)
emc <- single[[1]]

# ---------- prior ----------
mu_ref <- emc[[1]]$prior$theta_mu_mean
Sigma_ref <- emc[[1]]$prior$theta_mu_var
data_list <- lapply(single, function(x) x[[1]]$data[[1]])
model_factory <- emc[[1]]$model
param_names <- names(mu_ref)

loglik_emc2 <- function(Theta, data_i) {
  Theta <- as.matrix(Theta)
  colnames(Theta) <- param_names
  as.numeric(EMC2:::calc_ll_manager(Theta, data_i, model_factory, r_cores = 1L))
}

run_one_subject <- function(data_i, subject_id) {
  run_tempered_smc(
    bridge_stat_fn = function(Theta) loglik_emc2(Theta, data_i),
    reference_prior = make_reference_prior_gaussian(
      mu = mu_ref,
      Sigma = Sigma_ref,
      param_names = param_names
    ),
    M = M,
    n_cores = 1L,
    seed = base_seed + subject_id - 1L,
    verbose = verbose
  )
}

smc <- mclapply(
  seq_along(data_list),
  function(i) run_one_subject(data_list[[i]], i),
  mc.cores = mc.cores
)

recovery_stats <- vector("list", length(smc))
grDevices::pdf(posterior_plot_file, width = 12, height = 8)
for (i in seq_along(smc)) {
  reference_post <- parameters(single[[i]], selection = "alpha")
  reference_post <- reference_post[, -1L, drop = FALSE]
  posterior_draws <- smc_posteriors(smc[[i]], seed = base_seed + 1000L + i)
  reference_post <- reference_post[, colnames(posterior_draws), drop = FALSE]

  post_q <- posterior_intervals(posterior_draws)
  ref_q <- posterior_intervals(reference_post)

  recovery_stats[[i]] <- data.frame(
    parameter = colnames(posterior_draws),
    posterior_mean = colMeans(posterior_draws),
    posterior_median = post_q$q500,
    posterior_ci_lower = post_q$q025,
    posterior_ci_upper = post_q$q975,
    reference_median = ref_q$q500,
    reference_ci_lower = ref_q$q025,
    reference_ci_upper = ref_q$q975,
    row.names = NULL,
    check.names = FALSE
  )

  plot_posteriors(
    posterior_draws,
    reference_post,
    labels = c("SMC", "Reference"),
    cols = c("steelblue", "firebrick")
  )
  graphics::mtext(sprintf("Subject %d posterior overlays", i), outer = TRUE, cex = 1.1)
}
grDevices::dev.off()

saveRDS(
  list(
    smc = smc,
    recovery_stats = recovery_stats
  ),
  results_file
)

cat("Saved results to:", results_file, "\n")
cat("Saved plots to:", posterior_plot_file, "\n")

for (i in seq_along(recovery_stats)) {
  cat("\nSubject", i, "recovery summary:\n")
  print(recovery_stats[[i]], row.names = FALSE)
}
