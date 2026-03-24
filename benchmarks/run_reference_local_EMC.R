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

normalize_weights <- function(w) {
  w <- pmax(as.numeric(w), 0)
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) {
    rep(1 / length(w), length(w))
  } else {
    w / sw
  }
}

weighted_density_smooth <- function(x, w, n_resample = 5000L, adjust = 1) {
  idx <- sample.int(length(x), n_resample, replace = TRUE, prob = normalize_weights(w))
  stats::density(x[idx], adjust = adjust)
}

plot_reference_comparison <- function(broad_fit,
                                      refined_fit,
                                      local_id,
                                      smooth_adjust = 1,
                                      n_resample = 5000L) {
  Theta_broad <- as.matrix(broad_fit$Theta)
  Theta_refined <- as.matrix(refined_fit$Theta)
  w_broad <- normalize_weights(broad_fit$w)
  w_refined <- normalize_weights(refined_fit$w)

  d <- ncol(Theta_refined)
  n_cols <- min(3L, d)
  n_rows <- ceiling(d / n_cols)

  par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))

  for (j in seq_len(d)) {
    dens_broad <- weighted_density_smooth(
      Theta_broad[, j],
      w_broad,
      n_resample = n_resample,
      adjust = smooth_adjust
    )
    dens_refined <- weighted_density_smooth(
      Theta_refined[, j],
      w_refined,
      n_resample = n_resample,
      adjust = smooth_adjust
    )

    xlim <- range(c(dens_broad$x, dens_refined$x))
    ylim <- c(0, max(dens_broad$y, dens_refined$y))

    plot(
      dens_broad,
      xlim = xlim,
      ylim = ylim,
      lwd = 2,
      col = "firebrick",
      xlab = colnames(Theta_refined)[j],
      ylab = "Density",
      main = colnames(Theta_refined)[j]
    )
    lines(dens_refined, lwd = 2, col = "steelblue")
    legend(
      "topright",
      legend = c("Broad", "Refined"),
      col = c("firebrick", "steelblue"),
      lwd = 2,
      bty = "n",
      cex = 0.8
    )
  }

  mtext(sprintf("Local %d: broad vs refined reference posterior", local_id), outer = TRUE, cex = 1.1)
}

stage <- prepare_reference_local_stage(
  data_list = data_list,
  loglik_fn = loglik_emc2,
  base_mu = base_mu,
  base_Sigma = base_Sigma,
  pilot_size = pilot_size,
  broad_scale = 4,
  broad_defensive = TRUE,
  pilot_particles = pilot_particles,
  full_particles = local_particles,
  refined_method = "defensive_mixture",
  inflation = 1.5,
  defensive_weight = 0.10,
  defensive_scale = 4,
  n_jobs = mc.cores,
  pilot_local_n_cores = 1L,
  full_local_n_cores = 1L,
  base_seed = base_seed,
  verbose = verbose,
  pilot_smc_control = list(
    max_rounds = 50L,
    n_mcmc_moves = 1L,
    G_mix = 8L,
    hist_mix_enable = FALSE,
    gss_enable = FALSE,
    da_enable = FALSE
  ),
  full_smc_control = list()
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

# pdf(posterior_plot_file, width = 12, height = 8)
for (local_id in pilot_ids) {
  broad_fit <- stage$pilot$fits[[as.character(local_id)]]
  refined_fit <- stage$local_fits[[as.character(local_id)]]
  plot_reference_comparison(
    broad_fit = broad_fit,
    refined_fit = refined_fit,
    local_id = local_id
  )
}
# dev.off()

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
