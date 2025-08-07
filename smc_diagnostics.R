# SMC Diagnostics and Visualization Functions
# Contains functions for analyzing SMC results and parameter recovery

# Required for better density estimation
if (!requireNamespace("KernSmooth", quietly = TRUE)) {
  cat("Note: Install KernSmooth package for even better density estimation\n")
}

#' Alternative weighted density estimation using sampling approach
#' 
#' This function creates a smooth posterior by resampling particles according to weights
weighted_density_smooth <- function(x, w, n_resample = 5000, adjust = 1) {
  # Resample particles according to weights to create smooth density
  idx <- sample(length(x), n_resample, replace = TRUE, prob = w)
  x_resampled <- x[idx]
  
  # Use R's built-in density estimation on resampled data
  dens <- density(x_resampled, adjust = adjust)
  list(x = dens$x, y = dens$y)
}

#' Diagnose Parameter Recovery from SMC Results
#'
#' Creates diagnostic plots showing posterior distributions vs true parameter values
#'
#' @param true_values Named vector of true parameter values
#' @param smc_res SMC results object from enhanced_smc_elite()
#' @param main_title Main title for the plot (optional)
#' @param save_plot Logical, whether to save plot to file (default: FALSE)
#' @param filename Filename for saved plot (default: "smc_diagnostics.png")
#' @param width Plot width in inches (default: 12)
#' @param height Plot height in inches (default: 8)
#' @param density_method Method for density estimation: "kernel" (default) or "smooth" 
#' @param smooth_adjust Bandwidth adjustment for smooth method (default: 1)
#'
#' @return Invisibly returns summary statistics for each parameter
diagnose_recovery <- function(true_values, smc_res, 
                             main_title = "SMC Parameter Recovery Diagnostics",
                             save_plot = FALSE, filename = "smc_diagnostics.png",
                             width = 12, height = 8, density_method = "kernel",
                             smooth_adjust = 1) {
  
  # Extract particles and weights
  Theta <- smc_res$Theta
  w <- smc_res$w
  D <- ncol(Theta)
  param_names <- colnames(Theta)
  
  # If no parameter names, create generic ones
  if (is.null(param_names)) {
    param_names <- paste0("theta", 1:D)
    colnames(Theta) <- param_names
  }
  
  # Ensure true_values has matching names
  if (is.null(names(true_values))) {
    names(true_values) <- param_names
  }
  
  # Check dimensions match
  if (length(true_values) != D) {
    stop("Length of true_values must match number of parameters in smc_res")
  }
  
  # Set up plotting device
  if (save_plot) {
    png(filename, width = width, height = height, units = "in", res = 300)
  }
  
  # Set up multi-panel plot
  n_cols <- min(3, D)
  n_rows <- ceiling(D / n_cols)
  par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 3, 1), oma = c(0, 0, 3, 0))
  
  # Compute weighted quantiles function
  weighted_quantile <- function(x, w, probs) {
    ord <- order(x)
    x_ord <- x[ord]
    w_ord <- w[ord] / sum(w)
    cw <- cumsum(w_ord)
    approx(cw, x_ord, xout = probs, rule = 2)$y
  }
  
  # Storage for summary statistics
  summary_stats <- data.frame(
    parameter = param_names,
    true_value = true_values[param_names],
    posterior_mean = numeric(D),
    posterior_median = numeric(D),
    ci_lower = numeric(D),
    ci_upper = numeric(D),
    coverage = logical(D),
    stringsAsFactors = FALSE
  )
  
  # Plot each parameter
  for (i in 1:D) {
    param_name <- param_names[i]
    true_val <- true_values[param_name]
    theta_i <- Theta[, i]
    
    # Compute posterior statistics
    post_mean <- sum(theta_i * w)
    quantiles <- weighted_quantile(theta_i, w, c(0.025, 0.5, 0.975))
    post_median <- quantiles[2]
    ci_lower <- quantiles[1]
    ci_upper <- quantiles[3]
    coverage <- (true_val >= ci_lower) && (true_val <= ci_upper)
    
    # Store statistics
    summary_stats[i, "posterior_mean"] <- post_mean
    summary_stats[i, "posterior_median"] <- post_median
    summary_stats[i, "ci_lower"] <- ci_lower
    summary_stats[i, "ci_upper"] <- ci_upper
    summary_stats[i, "coverage"] <- coverage
    
    # Create weighted density estimate
    if (density_method == "smooth") {
      # Use resampling-based smooth density
      dens_result <- weighted_density_smooth(theta_i, w, n_resample = 5000, adjust = smooth_adjust)
      x_grid <- dens_result$x
      density_vals <- dens_result$y
    } else {
      # Use proper weighted kernel density estimation
      x_range <- range(theta_i)
      x_expand <- diff(x_range) * 0.15
      x_range <- c(x_range[1] - x_expand, x_range[2] + x_expand)
      
      # Improved weighted kernel density estimation
      n_grid <- 200
      x_grid <- seq(x_range[1], x_range[2], length.out = n_grid)
      
      # Adaptive bandwidth using Silverman's rule of thumb, adjusted for weights
      n_eff <- sum(w)^2 / sum(w^2)  # Effective sample size
      sigma_hat <- sqrt(sum(w * (theta_i - post_mean)^2))  # Weighted std dev
      # Silverman's bandwidth with effective sample size
      bw <- 1.06 * sigma_hat * (n_eff^(-1/5))
      bw <- max(bw, diff(range(theta_i)) / 50)  # Minimum bandwidth
      
      # Weighted kernel density
      density_vals <- numeric(n_grid)
      for (j in 1:n_grid) {
        # Gaussian kernel weights
        kernel_vals <- exp(-0.5 * ((theta_i - x_grid[j]) / bw)^2) / (bw * sqrt(2 * pi))
        density_vals[j] <- sum(w * kernel_vals)
      }
      density_vals <- density_vals / sum(w)  # Normalize by total weight
    }
    
    # Plot posterior density
    plot(x_grid, density_vals, type = "l", lwd = 2, col = "steelblue",
         xlab = param_name, ylab = "Posterior Density",
         main = sprintf("%s\n(True: %.3f, Est: %.3f)", param_name, true_val, post_mean))
    
    # Add credible interval shading
    x_ci <- x_grid[x_grid >= ci_lower & x_grid <= ci_upper]
    y_ci <- density_vals[x_grid >= ci_lower & x_grid <= ci_upper]
    if (length(x_ci) > 0) {
      polygon(c(x_ci, rev(x_ci)), c(y_ci, rep(0, length(y_ci))), 
               col = adjustcolor("lightblue", alpha = 0.3), border = NA)
    }
    
    # Add true value line
    abline(v = true_val, col = "red", lwd = 3, lty = 1)
    
    # Add posterior mean line
    abline(v = post_mean, col = "darkblue", lwd = 2, lty = 2)
    
    # Add legend for first plot
    if (i == 1) {
      legend("topright", 
             legend = c("True Value", "Posterior Mean", "95% CI"),
             col = c("red", "darkblue", "lightblue"),
             lty = c(1, 2, 1), lwd = c(3, 2, 8),
             cex = 0.8, bg = "white")
    }
    
    # Add coverage indicator
    coverage_text <- ifelse(coverage, "✓", "✗")
    coverage_color <- ifelse(coverage, "darkgreen", "red")
    text(x = par("usr")[2], y = par("usr")[4], 
         labels = coverage_text, col = coverage_color, 
         cex = 1.5, font = 2, adj = c(1.1, 1.1))
  }
  
  # Add overall title
  mtext(main_title, outer = TRUE, cex = 1.2, font = 2)
  
  # Add summary text at bottom
  n_covered <- sum(summary_stats$coverage)
  coverage_pct <- round(100 * n_covered / D)
  summary_text <- sprintf("95%% CI Coverage: %d/%d parameters (%d%%)", 
                         n_covered, D, coverage_pct)
  
  mtext(summary_text, side = 1, outer = TRUE, line = -1, 
        col = ifelse(coverage_pct >= 80, "darkgreen", "red"), font = 2)
  
  if (save_plot) {
    dev.off()
    cat(sprintf("Diagnostic plot saved to: %s\n", filename))
  }
  
  # Print summary table
  cat("\nParameter Recovery Summary:\n")
  cat("===========================\n")
  cat(sprintf("%-10s %8s %8s %8s %10s %10s %8s\n", 
              "Parameter", "True", "Post.Mean", "Error", "CI.Lower", "CI.Upper", "Coverage"))
  cat(paste(rep("-", 70), collapse = ""), "\n")
  
  for (i in 1:D) {
    error <- summary_stats$posterior_mean[i] - summary_stats$true_value[i]
    coverage_symbol <- ifelse(summary_stats$coverage[i], "✓", "✗")
    cat(sprintf("%-10s %8.3f %8.3f %+8.3f %10.3f %10.3f %8s\n",
                summary_stats$parameter[i],
                summary_stats$true_value[i],
                summary_stats$posterior_mean[i],
                error,
                summary_stats$ci_lower[i],
                summary_stats$ci_upper[i],
                coverage_symbol))
  }
  
  cat(paste(rep("-", 70), collapse = ""), "\n")
  cat(sprintf("Overall 95%% CI Coverage: %d/%d (%d%%)\n", n_covered, D, coverage_pct))
  
  # Return summary statistics invisibly
  invisible(summary_stats)
}

#' Additional SMC Convergence Diagnostics
#'
#' Plots convergence diagnostics for the SMC algorithm
#'
#' @param smc_res SMC results object from enhanced_smc_elite()
#' @param save_plot Logical, whether to save plot to file
#' @param filename Filename for saved plot
convergence_diagnostics <- function(smc_res, save_plot = FALSE, 
                                   filename = "smc_convergence.png") {
  
  if (save_plot) {
    png(filename, width = 10, height = 6, units = "in", res = 300)
  }
  
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  
  # Plot 1: Lambda progression
  lambda_hist <- smc_res$meta$lambda_hist
  rounds <- 0:(length(lambda_hist) - 1)
  
  plot(rounds, lambda_hist, type = "b", pch = 19, col = "steelblue",
       xlab = "Round", ylab = "Lambda (Tempering Parameter)",
       main = "Tempering Schedule", lwd = 2)
  grid()
  abline(h = 1, col = "red", lty = 2)
  
  # Plot 2: Acceptance rates
  if (length(smc_res$meta$acc_hist) > 0) {
    acc_rounds <- 1:length(smc_res$meta$acc_hist)
    plot(acc_rounds, smc_res$meta$acc_hist, type = "b", pch = 19, col = "darkgreen",
         xlab = "Round", ylab = "Acceptance Rate",
         main = "MCMC Acceptance Rates", lwd = 2, ylim = c(0, 1))
    grid()
    abline(h = c(0.2, 0.5), col = c("red", "blue"), lty = 2)
    legend("topright", legend = c("Target Range"), 
           col = c("blue"), lty = 2, cex = 0.8)
  }
  
  if (save_plot) {
    dev.off()
    cat(sprintf("Convergence diagnostics saved to: %s\n", filename))
  }
}

#' Get Better SMC Posteriors - Tuning Suggestions
#'
#' Print recommendations for improving SMC posterior quality
#'
#' @param smc_res SMC results object from enhanced_smc_elite()
#' @param target_ess_pct Target ESS percentage (default: 70%)
print_smc_tuning_advice <- function(smc_res, target_ess_pct = 70) {
  final_ess_pct <- 100 * smc_res$meta$ess / length(smc_res$w)
  final_acc_rate <- if (length(smc_res$meta$acc_hist) > 0) tail(smc_res$meta$acc_hist, 1) else NA
  
  cat("\n", paste(rep("=", 60), collapse=""), "\n")
  cat("SMC POSTERIOR QUALITY TUNING ADVICE\n")
  cat(paste(rep("=", 60), collapse=""), "\n")
  
  cat(sprintf("Current Performance:\n"))
  cat(sprintf("  - Final ESS: %.1f%% (target: >%d%%)\n", final_ess_pct, target_ess_pct))
  if (!is.na(final_acc_rate)) {
    cat(sprintf("  - Final acceptance rate: %.3f (target: 0.2-0.6)\n", final_acc_rate))
  }
  cat(sprintf("  - Rounds to convergence: %d\n", smc_res$meta$rounds))
  
  recommendations <- c()
  
  # ESS recommendations
  if (final_ess_pct < target_ess_pct) {
    recommendations <- c(recommendations, 
      "• INCREASE M (particles): More particles → better diversity",
      "• DECREASE resample_threshold: Resample less frequently", 
      "• INCREASE n_mcmc_moves: More MCMC steps after resampling")
  }
  
  # Acceptance rate recommendations  
  if (!is.na(final_acc_rate)) {
    if (final_acc_rate < 0.15) {
      recommendations <- c(recommendations,
        "• DECREASE rw_scale_init: Smaller random walk steps",
        "• INCREASE cov_inflation: More spread in mixture proposals")
    } else if (final_acc_rate > 0.7) {
      recommendations <- c(recommendations,
        "• INCREASE rw_scale_init: Larger random walk steps for better exploration")
    }
  }
  
  # Convergence speed recommendations
  if (smc_res$meta$rounds > 20) {
    recommendations <- c(recommendations,
      "• DECREASE target ESS: Allow more aggressive tempering",
      "• INCREASE dalpha_floor_factor: Larger minimum steps")
  }
  
  # Transport map recommendations
  lambda_hist <- smc_res$meta$lambda_hist
  if (length(lambda_hist) > 2 && (lambda_hist[2] - lambda_hist[1]) < 0.01) {
    recommendations <- c(recommendations,
      "• CHECK transport map quality: Very small first step suggests poor Gaussianization")
  }
  
  # General recommendations for better posteriors
  recommendations <- c(recommendations,
    "",
    "For Consistently Filled Posteriors:",
    "• USE density_method='smooth' in diagnose_recovery() for smoother plots",
    "• INCREASE M to 2000-5000 for complex posteriors", 
    "• SET hist_mix_enable=TRUE for better late-stage exploration",
    "• DECREASE refit_every to update transport maps more often",
    "• CONSIDER multi-start: run SMC multiple times with different seeds")
  
  if (length(recommendations) > 0) {
    cat("\nRecommendations for improvement:\n")
    for (rec in recommendations) {
      if (nzchar(rec)) cat(sprintf("  %s\n", rec))
      else cat("\n")
    }
  } else {
    cat("\nExcellent performance! No tuning needed.\n")
  }
  
  cat("\nExample improved settings:\n")
  cat("enhanced_smc_elite(\n")
  cat("  M = 3000,                    # More particles\n") 
  cat("  resample_threshold = 0.3,    # Less frequent resampling\n")
  cat("  n_mcmc_moves = 3,            # More MCMC steps\n")
  cat("  G_mix = min(15, 3*D),        # More mixture components\n")
  cat("  hist_mix_enable = TRUE,      # History-based proposals\n") 
  cat("  refit_every = 3,             # More frequent transport updates\n")
  cat("  use_qmc = TRUE               # Quasi-Monte Carlo\n")
  cat(")\n")
  
  cat(paste(rep("=", 60), collapse=""), "\n")
}

#' Create Multiple SMC Runs for Robust Posteriors
#'
#' Run SMC multiple times and combine results for better posterior coverage
#'
#' @param n_runs Number of independent SMC runs (default: 3)
#' @param ... Arguments passed to enhanced_smc_elite()
#' @return Combined results with pooled particles and weights
multi_run_smc <- function(n_runs = 3, ...) {
  cat(sprintf("Running %d independent SMC chains...\n", n_runs))
  
  results <- list()
  for (i in 1:n_runs) {
    cat(sprintf("Chain %d/%d: ", i, n_runs))
    
    # Use different seed for each run
    args <- list(...)
    args$seed <- (args$seed %||% 123) + i * 1000
    
    results[[i]] <- do.call(enhanced_smc_elite, args)
    cat(sprintf("ESS=%.1f%%\n", 100 * results[[i]]$meta$ess / length(results[[i]]$w)))
  }
  
  # Combine results
  combined_Theta <- do.call(rbind, lapply(results, `[[`, "Theta"))
  combined_w <- do.call(c, lapply(results, `[[`, "w"))
  combined_w <- combined_w / sum(combined_w)  # Renormalize
  
  # Use first result as template and update key fields
  combined_result <- results[[1]]
  combined_result$Theta <- combined_Theta
  combined_result$w <- combined_w
  combined_result$meta$ess <- sum(combined_w)^2 / sum(combined_w^2)
  combined_result$meta$multi_run <- TRUE
  combined_result$meta$n_runs <- n_runs
  
  cat(sprintf("Combined: %d total particles, ESS=%.1f%%\n", 
              nrow(combined_Theta), 100 * combined_result$meta$ess / length(combined_w)))
  
  combined_result
}

# Helper for NULL coalescing
`%||%` <- function(x, y) if (is.null(x)) y else x 