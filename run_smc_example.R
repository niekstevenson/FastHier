#!/usr/bin/env Rscript
# Example script demonstrating the SMC algorithm on a 5-parameter nonlinear regression

rm(list = ls())

# Load required libraries
library(mvtnorm)

# Source the SMC algorithm, likelihood functions, and diagnostics
source("SMC_super_fast.R")
source("hard_ll.R")
source("smc_diagnostics.R")

# Set random seed for reproducibility
set.seed(123)

# ---------- Generate synthetic data ----------
cat(paste(rep("=", 60), collapse=""), "\n")
cat("SMC ALGORITHM DEMONSTRATION\n")
cat(paste(rep("=", 60), collapse=""), "\n")

data_setup <- gen_data(N = 200, sigma_y = 0.2, prior_scale = 2, seed = 123)

cat("Data generation complete:\n")
cat(sprintf("  N = %d observations, D = %d parameters\n", data_setup$N, data_setup$D))
cat(sprintf("  True parameters: [%.3f, %.3f, %.3f, %.3f, %.3f]\n",
            data_setup$theta_true[1], data_setup$theta_true[2], data_setup$theta_true[3],
            data_setup$theta_true[4], data_setup$theta_true[5]))
cat(sprintf("  sigma_y = %.3f\n", data_setup$sigma_y))

# Extract data and prior
data <- list(X = data_setup$X, y = data_setup$y, sigma_y = data_setup$sigma_y)
mu_ref <- data_setup$mu
Sigma_ref <- data_setup$Sigma

# Test likelihood function
cat("\nTesting likelihood function...\n")
theta_test <- matrix(data_setup$theta_true, 1, 5)
ll_test <- loglik_sinmix(theta_test, data)
cat(sprintf("  Likelihood at true parameters: %.3f\n", ll_test))

# ---------- Run SMC algorithm ----------
cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("RUNNING SMC ALGORITHM\n")
cat(paste(rep("=", 60), collapse=""), "\n")

start_time <- Sys.time()

# Use automatic adaptive SMC
smc_res <- auto_smc(
  data = data,
  loglik_fn = loglik_sinmix,
  mu_ref = mu_ref,
  Sigma_ref = Sigma_ref,
  max_rounds = 100,
  seed = 12
)

end_time <- Sys.time()
runtime <- as.numeric(difftime(end_time, start_time, units = "secs"))

# ---------- Results and Diagnostics ----------
cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("RESULTS SUMMARY\n")
cat(paste(rep("=", 60), collapse=""), "\n")

cat(sprintf("Algorithm completed in %.2f seconds\n", runtime))
cat(sprintf("Converged in %d rounds (lambda = %.4f)\n",
            smc_res$meta$rounds, smc_res$final_lambda))
cat(sprintf("Final ESS: %.1f (%.1f%%)\n",
            smc_res$meta$ess, 100 * smc_res$meta$ess / length(smc_res$w)))

# Parameter estimates (weighted mean and quantiles)
theta_est <- colSums(smc_res$Theta * smc_res$w)
theta_true <- data_setup$theta_true

# Compute weighted quantiles
weighted_quantile <- function(x, w, probs) {
  ord <- order(x)
  x_ord <- x[ord]
  w_ord <- w[ord] / sum(w)
  cw <- cumsum(w_ord)
  approx(cw, x_ord, xout = probs, rule = 2)$y
}

cat("\n", paste(rep("-", 60), collapse=""), "\n")
cat("PARAMETER RECOVERY ANALYSIS\n")
cat(paste(rep("-", 60), collapse=""), "\n")
cat(sprintf("%-8s %8s %8s %8s %10s %10s %10s\n",
            "Param", "True", "Estimate", "Error", "Q2.5%", "Q50%", "Q97.5%"))
cat(paste(rep("-", 60), collapse=""), "\n")

for (i in 1:5) {
  q <- weighted_quantile(smc_res$Theta[, i], smc_res$w, c(0.025, 0.5, 0.975))
  cat(sprintf("theta%-3d %8.3f %8.3f %+8.3f %10.3f %10.3f %10.3f\n",
              i, theta_true[i], theta_est[i], theta_est[i] - theta_true[i],
              q[1], q[2], q[3]))
}

# Coverage check (true parameter within 95% credible interval)
coverage <- numeric(5)
for (i in 1:5) {
  q <- weighted_quantile(smc_res$Theta[, i], smc_res$w, c(0.025, 0.975))
  coverage[i] <- (theta_true[i] >= q[1]) && (theta_true[i] <= q[2])
}
cat(sprintf("\n95%% Credible interval coverage: %d/5 parameters\n", sum(coverage)))

# ---------- Algorithm Diagnostics ----------
cat("\n", paste(rep("-", 60), collapse=""), "\n")
cat("ALGORITHM DIAGNOSTICS\n")
cat(paste(rep("-", 60), collapse=""), "\n")

cat(sprintf("Final random walk scale: %.3f\n", smc_res$meta$rw_scale_final))

if (length(smc_res$meta$acc_hist) > 0) {
  cat("Acceptance rates by round:\n")
  for (i in seq_along(smc_res$meta$acc_hist)) {
    cat(sprintf("  Round %2d: %.3f\n", i, smc_res$meta$acc_hist[i]))
  }
  cat(sprintf("Final acceptance rate: %.3f\n", tail(smc_res$meta$acc_hist, 1)))
}

# Lambda progression
cat("\nTempering schedule (lambda progression):\n")
for (i in seq_along(smc_res$meta$lambda_hist)) {
  cat(sprintf("  Round %2d: %.4f\n", i-1, smc_res$meta$lambda_hist[i]))
}

# Transport map quality (if available)
if (!is.null(smc_res$transport)) {
  cat("\nTransport map test (Gaussianity check):\n")
  Z_final <- smc_res$Z
  # Simple normality test on each transformed dimension
  for (j in 1:ncol(Z_final)) {
    ks_test <- ks.test(Z_final[, j], "pnorm")
    cat(sprintf("  Dimension %d: KS p-value = %.3f\n", j, ks_test$p.value))
  }
}

# ---------- Visual Diagnostics ----------
cat("\n", paste(rep("-", 60), collapse=""), "\n")
cat("VISUAL DIAGNOSTICS\n")
cat(paste(rep("-", 60), collapse=""), "\n")

# Create diagnostic plots showing posterior recovery
cat("Creating parameter recovery plots (improved density estimation)...\n")
recovery_stats <- diagnose_recovery(
  true_values = data_setup$theta_true,
  smc_res = smc_res,
  main_title = "5-Parameter Nonlinear Regression: Parameter Recovery",
  save_plot = FALSE,
  filename = "parameter_recovery_diagnostics.png",
  density_method = "smooth"  # Use smooth density estimation
)
#
# # Create convergence diagnostic plots
# cat("Creating convergence diagnostic plots...\n")
# convergence_diagnostics(
#   smc_res = smc_res,
#   save_plot = FALSE,
#   filename = "smc_convergence_diagnostics.png"
# )
#
# # Print tuning advice for better posteriors
# print_smc_tuning_advice(smc_res)
#
# cat("\n", paste(rep("=", 60), collapse=""), "\n")
# cat("Analysis complete!\n")
# cat("Diagnostic plots saved:\n")
# cat("  - parameter_recovery_diagnostics.png\n")
# cat("  - smc_convergence_diagnostics.png\n")
# cat(paste(rep("=", 60), collapse=""), "\n")
