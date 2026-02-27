#!/usr/bin/env Rscript

stan_path <- file.path("samples", "stan_full_results.rds")
inner_path <- file.path("samples", "inner_smc_results.rds")

if (!file.exists(stan_path)) stop("Missing ", stan_path, ". Run run_stan.R first.")
if (!file.exists(inner_path)) stop("Missing ", inner_path, ". Run run_two_step_smc.R first.")

stan_res <- readRDS(stan_path)
inner <- readRDS(inner_path)

weighted_quantile <- function(x, w, probs = c(0.025, 0.5, 0.975)) {
  o <- order(x)
  x <- x[o]
  w <- w[o]
  w <- w / sum(w)
  cw <- cumsum(w)
  sapply(probs, function(p) x[which(cw >= p)[1L]])
}

weighted_summary <- function(X, w) {
  w <- w / sum(w)
  means <- colSums(X * w)
  q <- vapply(seq_len(ncol(X)), function(j) weighted_quantile(X[, j], w), numeric(3))
  data.frame(
    parameter = colnames(X),
    mean = means,
    q2.5 = q[1, ],
    q50 = q[2, ],
    q97.5 = q[3, ],
    row.names = NULL,
    check.names = FALSE
  )
}

stan_summary <- function(X, names) {
  data.frame(
    parameter = names,
    mean = colMeans(X),
    q2.5 = apply(X, 2L, quantile, probs = 0.025),
    q50 = apply(X, 2L, quantile, probs = 0.5),
    q97.5 = apply(X, 2L, quantile, probs = 0.975),
    row.names = NULL,
    check.names = FALSE
  )
}

phi <- as.matrix(inner$phi)
w <- as.numeric(inner$w) / sum(as.numeric(inner$w))

phi_nat <- cbind(
  mu_shape = phi[, 1L],
  mu_scale = phi[, 2L],
  mu_shift = phi[, 3L],
  sigma2_shape = exp(phi[, 4L]),
  sigma2_scale = exp(phi[, 5L]),
  sigma2_shift = exp(phi[, 6L])
)
smc_sum <- weighted_summary(phi_nat, w)

stan_mu <- as.matrix(stan_res$draws$mu)
stan_sigma2 <- as.matrix(stan_res$draws$sigma2)
stan_sum <- rbind(
  stan_summary(stan_mu, c("mu_shape", "mu_scale", "mu_shift")),
  stan_summary(stan_sigma2, c("sigma2_shape", "sigma2_scale", "sigma2_shift"))
)

cmp <- merge(stan_sum, smc_sum, by = "parameter", suffixes = c("_stan", "_smc"))
cmp$mean_diff <- cmp$mean_smc - cmp$mean_stan

cat("Posterior comparison (Stan vs Inner SMC):\n")
print(cmp[, c("parameter", "mean_stan", "mean_smc", "mean_diff")], row.names = FALSE)
cat(sprintf("\nInner log evidence: %.6f (MCSE %.6f)\n", inner$log_evidence, inner$mcse_log_evidence))

# Plot posterior densities for each inner parameter (interactive sessions only)
if (interactive()) {
  par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
  for (p in colnames(phi_nat)) {
    d_smc <- density(sample(phi_nat[, p], size = min(10000L, nrow(phi_nat)), replace = TRUE, prob = w))
    if (grepl("^mu_", p)) {
      j <- match(sub("^mu_", "", p), c("shape", "scale", "shift"))
      d_stan <- density(stan_mu[, j])
    } else {
      j <- match(sub("^sigma2_", "", p), c("shape", "scale", "shift"))
      d_stan <- density(stan_sigma2[, j])
    }
    plot(d_smc, main = p, xlab = p, ylab = "density", col = "black", lwd = 2)
    lines(d_stan, col = "firebrick", lwd = 2, lty = 2)
    legend("topright", legend = c("Inner SMC", "Stan"), col = c("black", "firebrick"),
           lwd = 2, lty = c(1, 2), bty = "n", cex = 0.9)
  }
} else {
  cat("\nNon-interactive session: skipping plots (run interactively to view posterior plots).\n")
}
