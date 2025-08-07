#!/usr/bin/env Rscript
# Example: Using the SMC algorithm with a custom likelihood function
# This demonstrates the modular structure of the SMC implementation

rm(list = ls())
library(mvtnorm)

# Source the SMC algorithm and diagnostics (no specific likelihood assumed)
source("SMC_super_fast.R")
source("smc_diagnostics.R")

# ---------- Define a custom likelihood function ----------
# Example: Simple linear regression with unknown slope, intercept, and error variance
# Parameters: theta = [beta0, beta1, log_sigma]
# Model: y = beta0 + beta1 * x + epsilon, epsilon ~ N(0, sigma^2)

custom_loglik <- function(theta, data) {
  # theta: M x 3 matrix of parameters [beta0, beta1, log_sigma]
  # data: list with elements x (N-vector) and y (N-vector)

  theta <- as.matrix(theta)
  M <- nrow(theta)
  x <- data$x
  y <- data$y
  N <- length(y)

  beta0 <- theta[, 1]     # intercept
  beta1 <- theta[, 2]     # slope
  log_sigma <- theta[, 3] # log error standard deviation
  sigma <- exp(log_sigma)

  # Compute predicted values for each parameter vector
  loglik <- numeric(M)
  for (m in 1:M) {
    y_pred <- beta0[m] + beta1[m] * x
    resid <- y - y_pred
    # Normal likelihood: sum of log densities
    loglik[m] <- sum(dnorm(resid, mean = 0, sd = sigma[m], log = TRUE))
  }

  loglik
}

# ---------- Generate synthetic data for linear regression ----------
set.seed(456)
N <- 100
x <- runif(N, -2, 2)
beta0_true <- 2.0
beta1_true <- -1.5
sigma_true <- 0.3
y <- beta0_true + beta1_true * x + rnorm(N, 0, sigma_true)

data <- list(x = x, y = y)

# ---------- Set up prior ----------
# Prior on [beta0, beta1, log_sigma]
mu_ref <- c(beta0 = 0, beta1 = 0, log_sigma = log(1))
Sigma_ref <- diag(c(4, 4, 1))  # fairly diffuse priors
param_names <- names(mu_ref)

cat("Custom likelihood example: Linear regression\n")
cat("============================================\n")
cat(sprintf("True parameters: beta0=%.2f, beta1=%.2f, sigma=%.2f\n",
            beta0_true, beta1_true, sigma_true))
cat(sprintf("Data: N=%d observations\n", N))

# Test likelihood function
theta_test <- matrix(c(beta0_true, beta1_true, log(sigma_true)), 1, 3)
ll_test <- custom_loglik(theta_test, data)
cat(sprintf("Likelihood at true parameters: %.3f\n", ll_test))

# ---------- Run SMC ----------
cat("\nRunning SMC algorithm...\n")

# Test the new automatic adaptive SMC
smc_res <- auto_smc(
  data = data,
  loglik_fn = custom_loglik,
  mu_ref = mu_ref,
  Sigma_ref = Sigma_ref,
  seed = 12345
)

# ---------- Results ----------
cat("\nResults:\n")
cat("========\n")
cat(sprintf("Converged in %d rounds\n", smc_res$meta$rounds))
cat(sprintf("Final ESS: %.1f (%.1f%%)\n",
            smc_res$meta$ess, 100 * smc_res$meta$ess / length(smc_res$w)))

# Parameter estimates
theta_est <- colSums(smc_res$Theta * smc_res$w)
theta_true_vec <- c(beta0_true, beta1_true, log(sigma_true))

cat("\nParameter estimates:\n")
cat("                  True   Estimate     Error\n")
cat(sprintf("Beta0 (intercept) %6.3f    %6.3f   %+6.3f\n",
            theta_true_vec[1], theta_est[1], theta_est[1] - theta_true_vec[1]))
cat(sprintf("Beta1 (slope)     %6.3f    %6.3f   %+6.3f\n",
            theta_true_vec[2], theta_est[2], theta_est[2] - theta_true_vec[2]))
cat(sprintf("log(sigma)        %6.3f    %6.3f   %+6.3f\n",
            theta_true_vec[3], theta_est[3], theta_est[3] - theta_true_vec[3]))
cat(sprintf("sigma             %6.3f    %6.3f   %+6.3f\n",
            sigma_true, exp(theta_est[3]), exp(theta_est[3]) - sigma_true))

# ---------- Visual Diagnostics ----------
cat("\nCreating diagnostic plots (improved density estimation)...\n")
true_vals <- c(beta0 = beta0_true, beta1 = beta1_true, log_sigma = log(sigma_true))
recovery_stats <- diagnose_recovery(
  true_values = true_vals,
  smc_res = smc_res,
  main_title = "Linear Regression: Parameter Recovery",
  save_plot = FALSE,
  filename = "linear_regression_diagnostics.png",
  density_method = "smooth"  # Use smooth density estimation
)

cat("\nThis demonstrates how to use the SMC algorithm with any likelihood function!\n")
cat("Simply define your likelihood as a function taking (parameters, data) and pass it to enhanced_smc_elite().\n")
cat("Diagnostic plot saved to: linear_regression_diagnostics.png\n")

