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

normalize_weights <- function(w) {
  w <- pmax(as.numeric(w), 0)
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) {
    rep(1 / length(w), length(w))
  } else {
    w / sw
  }
}

weighted_quantile <- function(x, w, probs = c(0.025, 0.5, 0.975)) {
  ord <- order(x)
  x_ord <- x[ord]
  w_ord <- normalize_weights(w[ord])
  cw <- cumsum(w_ord)
  as.numeric(stats::approx(cw, x_ord, xout = probs, rule = 2)$y)
}

weighted_density_smooth <- function(x, w, n_resample = 5000L, adjust = 1) {
  idx <- sample.int(length(x), n_resample, replace = TRUE, prob = normalize_weights(w))
  stats::density(x[idx], adjust = adjust)
}

summarize_recovery <- function(reference_values, smc_res) {
  Theta <- as.matrix(smc_res$Theta)
  w <- normalize_weights(smc_res$w)
  post_q <- t(vapply(seq_len(ncol(Theta)), function(j) weighted_quantile(Theta[, j], w), numeric(3)))

  if (is.matrix(reference_values) || is.data.frame(reference_values)) {
    ref_mat <- as.matrix(reference_values)
    if (!is.null(colnames(ref_mat))) {
      ref_mat <- ref_mat[, colnames(Theta), drop = FALSE]
    }
    ref_q <- t(apply(ref_mat, 2L, stats::quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))
    return(data.frame(
      parameter = colnames(Theta),
      posterior_mean = colSums(Theta * w),
      posterior_median = post_q[, 2L],
      posterior_ci_lower = post_q[, 1L],
      posterior_ci_upper = post_q[, 3L],
      reference_median = ref_q[, 2L],
      reference_ci_lower = ref_q[, 1L],
      reference_ci_upper = ref_q[, 3L],
      row.names = NULL,
      check.names = FALSE
    ))
  }

  truth <- as.numeric(reference_values[colnames(Theta)])
  data.frame(
    parameter = colnames(Theta),
    truth = truth,
    posterior_mean = colSums(Theta * w),
    posterior_median = post_q[, 2L],
    posterior_ci_lower = post_q[, 1L],
    posterior_ci_upper = post_q[, 3L],
    covered = truth >= post_q[, 1L] & truth <= post_q[, 3L],
    row.names = NULL,
    check.names = FALSE
  )
}

plot_subject_posteriors <- function(reference_values, smc_res, subject_id,
                                    smooth_adjust = 1, n_resample = 5000L) {
  Theta <- as.matrix(smc_res$Theta)
  w <- normalize_weights(smc_res$w)
  d <- ncol(Theta)
  n_cols <- min(3L, d)
  n_rows <- ceiling(d / n_cols)

  ref_is_chain <- is.matrix(reference_values) || is.data.frame(reference_values)
  ref_mat <- NULL
  if (ref_is_chain) {
    ref_mat <- as.matrix(reference_values)
    if (!is.null(colnames(ref_mat))) {
      ref_mat <- ref_mat[, colnames(Theta), drop = FALSE]
    }
  }

  par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))

  for (j in seq_len(d)) {
    posterior_density <- weighted_density_smooth(
      Theta[, j],
      w,
      n_resample = n_resample,
      adjust = smooth_adjust
    )

    xlim <- range(posterior_density$x)
    ylim <- c(0, max(posterior_density$y))

    if (ref_is_chain) {
      reference_density <- stats::density(ref_mat[, j], adjust = smooth_adjust)
      xlim <- range(c(xlim, reference_density$x))
      ylim <- c(0, max(ylim[2], reference_density$y))
    }

    plot(
      posterior_density,
      xlim = xlim,
      ylim = ylim,
      lwd = 2,
      col = "steelblue",
      xlab = colnames(Theta)[j],
      ylab = "Density",
      main = colnames(Theta)[j]
    )

    if (ref_is_chain) {
      lines(reference_density, lwd = 2, col = "firebrick")
      legend(
        "topright",
        legend = c("SMC", "Reference"),
        col = c("steelblue", "firebrick"),
        lwd = 2,
        bty = "n",
        cex = 0.8
      )
    } else {
      abline(v = as.numeric(reference_values[colnames(Theta)[j]]), col = "firebrick", lwd = 2, lty = 2)
      legend(
        "topright",
        legend = c("SMC", "Truth"),
        col = c("steelblue", "firebrick"),
        lwd = 2,
        lty = c(1, 2),
        bty = "n",
        cex = 0.8
      )
    }
  }

  mtext(sprintf("Subject %d posterior overlays", subject_id), outer = TRUE, cex = 1.1)
}

loglik_emc2 <- function(Theta, data_i) {
  Theta <- as.matrix(Theta)
  colnames(Theta) <- param_names
  as.numeric(EMC2:::calc_ll_manager(Theta, data_i, model_factory, r_cores = 1L))
}

run_one_subject <- function(data_i, subject_id) {
  enhanced_smc_elite(
    data = data_i,
    loglik_fn = loglik_emc2,
    mu_ref = mu_ref,
    Sigma_ref = Sigma_ref,
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
# pdf(posterior_plot_file, width = 12, height = 8)
for (i in seq_along(smc)) {
  reference_post <- parameters(single[[i]], selection = "alpha")
  reference_post <- reference_post[, -1L, drop = FALSE]
  recovery_stats[[i]] <- summarize_recovery(reference_post, smc[[i]])
  plot_subject_posteriors(reference_post, smc[[i]], subject_id = i)
}
# dev.off()

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
