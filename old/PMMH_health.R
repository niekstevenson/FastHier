# ---- helper: per-subject posterior mean under current eta using the cache ----
subject_mean_from_cache <- function(cache, mu, tau2) {
  Theta <- as.matrix(cache$theta)
  # Δ = log p_eta(θ) - log p0(θ)
  d_all <- log_prod_normal_rows(Theta, mu, tau2) - cache$log_p0
  lw <- log(pmax(cache$w, .Machine$double.eps)) + d_all
  lw <- lw - max(lw)
  w_eta <- exp(lw); w_eta <- w_eta / sum(w_eta)
  colSums(Theta * w_eta)
}

# ---- single-state check (for any mu, log_tau2) ----
check_outer_consistency_once <- function(mu, log_tau2, caches, verbose = TRUE) {
  tau2 <- exp(log_tau2)
  subj_means <- t(vapply(caches, subject_mean_from_cache, numeric(length(mu)),
                         mu = mu, tau2 = tau2))
  colnames(subj_means) <- names(mu)
  mbar <- colMeans(subj_means)
  delta <- as.numeric(mu - mbar)
  names(delta) <- names(mu)

  # Across-subject variance of subject means -> rough SE of the average
  s2 <- apply(subj_means, 2, stats::var)
  se <- sqrt(pmax(s2, 0) / nrow(subj_means) + 1e-12)
  z <- delta / se

  out <- list(
    mu = mu,
    mbar_subject_means = mbar,
    delta = delta,
    l2_norm = sqrt(sum(delta^2)),
    max_abs = max(abs(delta)),
    z = z,
    se = se
  )
  if (isTRUE(verbose)) {
    cat("Outer consistency check (μ vs mean of subject posterior means):\n")
    print(round(rbind(mu = mu,
                      mean_subject = mbar,
                      delta = delta,
                      z = z), 3))
    cat(sprintf("L2(μ - mean_subject) = %.4f | max|Δ| = %.4f\n",
                out$l2_norm, out$max_abs))
  }
  out
}

# ---- chain-wide quick diagnostic (thin and check a few iterations) ----
check_outer_consistency_chain <- function(fit, caches, iters = NULL, thin = 100, verbose = TRUE) {
  stopifnot(is.matrix(fit$mu), is.matrix(fit$log_tau2))
  T <- nrow(fit$mu)
  if (is.null(iters)) iters <- unique(pmax(1, seq(1, T, by = thin)))
  res <- lapply(iters, function(it) {
    mu <- fit$mu[it, ]; names(mu) <- colnames(fit$mu)
    lt2 <- fit$log_tau2[it, ]; names(lt2) <- colnames(fit$log_tau2)
    ch <- check_outer_consistency_once(mu, lt2, caches, verbose = FALSE)
    c(iter = it, ch$l2_norm, ch$max_abs, setNames(ch$delta, paste0("delta_", names(ch$delta))))
  })
  df <- as.data.frame(do.call(rbind, res))
  if (isTRUE(verbose)) {
    head_cols <- c("iter", "l2_norm", "max_abs")
    print(utils::head(df[ , head_cols, drop = FALSE], 10))
  }
  df
}

# --- example usage ---
# final state:
# check_outer_consistency_once(fit$final_state$mu, fit$final_state$log_tau2, caches)
