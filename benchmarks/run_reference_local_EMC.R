rm(list = ls())
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[1L]))
} else {
  normalizePath("benchmarks/run_reference_local_EMC.R")
}
repo_dir <- dirname(dirname(script_path))
setwd(repo_dir)

suppressPackageStartupMessages({
  library(EMC2)
})

set.seed(123L)

data_path <- file.path("benchmarks", "samples", "single_DDM30.RData")
results_file <- file.path("benchmarks", "results", "reference_local_emc_results.rds")
posterior_plot_file <- file.path("benchmarks", "results", "reference_local_emc_posteriors.pdf")

mc.cores <- 12L
pilot_size <- 20L
pilot_particles <- 1000L
local_particles <- 4000L
base_seed <- 123L
verbose <- TRUE

dir.create(file.path("benchmarks", "samples"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path("benchmarks", "results"), showWarnings = FALSE, recursive = TRUE)

source("hierarchical_locals.R")
source("utilities.R")

if (!file.exists(data_path)) {
  stop("Missing data file: ", data_path)
}

load(data_path)
if (!exists("single", inherits = FALSE)) {
  stop("The data file does not define an object named 'single'.")
}

emc <- single[[1]]
base_mu <- emc[[1]]$prior$theta_mu_mean
base_Sigma <- emc[[1]]$prior$theta_mu_var
param_names <- names(base_mu)
data_list <- lapply(single, function(x) x[[1]]$data[[1]])
model_factory <- emc[[1]]$model

loglik_emc2 <- function(Theta, data_i) {
  Theta <- as.matrix(Theta)
  colnames(Theta) <- param_names
  as.numeric(EMC2:::calc_ll_manager(Theta, data_i, model_factory, r_cores = 1L))
}

population_model <- default_population_model_diag_gaussian(
  base_mu = base_mu,
  base_Sigma = base_Sigma,
  mean_var_scale = 1.0,
  sigma2_shape = 3.0,
  sigma2_mean = diag(base_Sigma)
)

stage <- prepare_reference_local_stage(
  data_list = data_list,
  loglik_fn = loglik_emc2,
  base_mu = base_mu,
  base_Sigma = base_Sigma,
  population_model = population_model,
  broad_defensive = TRUE,
  n_jobs = mc.cores
)

saveRDS(
  list(
    stage = stage,
    settings = list(
      data_path = normalizePath(data_path, winslash = "/", mustWork = TRUE),
      results_file = results_file,
      mc.cores = mc.cores,
      pilot_size = pilot_size,
      pilot_particles = pilot_particles,
      local_particles = local_particles,
      base_seed = base_seed
    )
  ),
  results_file
)

broad_summary <- reference_prior_summary(stage$broad_reference)
refined_summary <- reference_prior_summary(stage$refined_reference)
pilot_ids <- stage$pilot$selection$indices

grDevices::pdf(posterior_plot_file, width = 12, height = 8)
for (local_id in pilot_ids) {
  broad_fit <- stage$pilot$fits[[as.character(local_id)]]
  refined_fit <- stage$local_fits[[as.character(local_id)]]
  broad_post <- smc_posteriors(broad_fit, seed = base_seed + local_id)
  refined_post <- smc_posteriors(refined_fit, seed = base_seed + 1000L + local_id)
  plot_posteriors(
    broad_post,
    refined_post,
    labels = c("Broad", "Refined"),
    cols = c("firebrick", "steelblue")
  )
  graphics::mtext(sprintf("Local %d: broad vs refined reference posterior", local_id), outer = TRUE, cex = 1.1)
}
grDevices::dev.off()

cat("Saved results to:", results_file, "\n")
cat("Saved plots to:", posterior_plot_file, "\n")
cat(sprintf("Locals: %d\n", length(data_list)))
cat(sprintf("Pilot subset size: %d\n", length(stage$pilot$selection$indices)))
cat("Pilot subset:", paste(stage$pilot$selection$indices, collapse = ", "), "\n")
cat(sprintf("Broad reference: %s (%d component%s)\n",
            broad_summary$label,
            broad_summary$n_components,
            if (broad_summary$n_components == 1L) "" else "s"))
cat(sprintf("Refined reference: %s (%d component%s)\n",
            refined_summary$label,
            refined_summary$n_components,
            if (refined_summary$n_components == 1L) "" else "s"))
