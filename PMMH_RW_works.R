# ---------------------------------------------------------------------------
#  Modular fully‑correlated pseudo‑marginal MCMC (PMMH) sampler
#  ---------------------------------------------------------------------------
#  * Independence kernel: Normal–Inverse‑Wishart (NIW)
#  * Adaptive block Gaussian random‑walk (RW) kernel
#  * Correlated pseudo‑marginal estimator with per‑subject AR(1) drivers
#
#  The main entry‑point is `run_pmmh(cache, …)` which works with *any* person‑
#  level likelihood as long as you provide a cache of importance‑sampling
#  particles (Θ), their proposal log‑densities (logq) and exact log‑likelihoods
#  (loglik).
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(mvtnorm)
  library(coda)
})

# ---------- helper utilities ------------------------------------------------
logsumexp <- function(x) {
  m <- max(x); m + log(sum(exp(x - m)))
}

vec_from_L <- function(L) {
  d <- nrow(L)
  c(log(diag(L)), L[lower.tri(L)])
}

L_from_vec <- function(vec, d = NULL) {
  if (is.null(d)) d <- (-1 + sqrt(1 + 8 * length(vec))) / 2
  stopifnot(length(vec) == d * (d + 1) / 2)
  L <- matrix(0, d, d)
  diag(L) <- exp(vec[1:d])
  if (d > 1) L[lower.tri(L)] <- vec[(d + 1):length(vec)]
  L
}

Sigma_from_L_vec <- function(vec) {
  L <- L_from_vec(vec)
  L %*% t(L)
}

rinvwishart <- function(df, Psi) {
  Wi <- rWishart(1, df, solve(Psi))[,,1]
  solve(Wi)
}

dinvwishart <- function(S, df, Psi, log = FALSE) {
  p <- nrow(S)
  log_const <- -0.5 * df * p * log(2) -
    0.25 * p * (p - 1) * log(pi) -
    sum(lgamma(0.5 * (df + 1 - 1:p))) +
    0.5 * df * log(det(Psi))
  log_dens <- log_const -
    0.5 * (df + p + 1) * log(det(S)) -
    0.5 * sum(diag(Psi %*% solve(S)))
  if (log) log_dens else exp(log_dens)
}

generate_Z <- function(M_cache) rnorm(M_cache)
update_Z   <- function(Z_list, rho) {
  lapply(Z_list, function(z) rho * z + sqrt(1 - rho^2) * rnorm(length(z)))
}

compute_log_marginals <- function(mu, L_vec, Z_list, cache, p_inc, z_cut) {
  Sig <- Sigma_from_L_vec(L_vec)
  total <- 0
  M_cache <- length(Z_list[[1]])
  for (i in seq_along(cache)) {
    inc <- Z_list[[i]] < z_cut
    if (!any(inc)) { k <- which.min(Z_list[[i]]); inc[k] <- TRUE }
    idx <- which(inc)
    logw <- dmvnorm(cache[[i]]$Theta[idx, , drop = FALSE], mu, Sig, log = TRUE) -
      cache[[i]]$logq[idx] + cache[[i]]$loglik[idx] + log(1 / p_inc)
    total <- total + (-log(M_cache) + logsumexp(logw))
  }
  total
}

log_q_NIW <- function(mu, L_vec, mu_bar, nu_p, kappa_p, Psi_p) {
  Sig <- Sigma_from_L_vec(L_vec)
  log_dens_S  <- dinvwishart(Sig, nu_p, Psi_p, log = TRUE)
  log_dens_mu <- dmvnorm(mu, mu_bar, Sig / kappa_p, log = TRUE)
  log_dens_S + log_dens_mu
}

draw_NIW <- function(mu_bar, nu_p, kappa_p, Psi_p) {
  Sig <- rinvwishart(nu_p, Psi_p)
  mu  <- rmvnorm(1, mu_bar, Sig / kappa_p)[1, ]
  list(mu = mu, L_vec = vec_from_L(t(chol(Sig))))
}

# ---------- generic priors --------------------------------------------------
log_prior_mu <- function(mu, m0 = rep(0, length(mu)), s0 = rep(3, length(mu))) {
  sum(dnorm(mu, m0, s0, log = TRUE))
}

log_prior_L <- function(L_vec, meanlog_diag = 0, sdlog_diag = .5, sd_off = 1) {
  d <- (-1 + sqrt(1 + 8 * length(L_vec))) / 2
  diag_logs <- L_vec[1:d]
  off       <- if (d > 1) L_vec[(d + 1):length(L_vec)] else numeric(0)
  sum(dlnorm(exp(diag_logs), meanlog_diag, sdlog_diag, log = TRUE) + diag_logs) +
    sum(dnorm(off, 0, sd_off, log = TRUE))
}

log_prior_phi <- function(mu, L_vec, ...) {
  log_prior_mu(mu, ...) + log_prior_L(L_vec, ...)
}

# ---------- main PMMH sampler ----------------------------------------------
run_pmmh <- function(cache,
                     n_iter      = 6000L,
                     burn        = 1000L,
                     thin        = 5L,
                     p_rw        = 0.25,
                     rho_sub     = 0.999,
                     m_use       = 100L,
                     acc_target  = 0.40,
                     adapt_start = 100L,
                     gamma_da    = 0.05,
                     t0          = 10,
                     kappa_da    = 0.75,
                     seed        = 1) {

  set.seed(seed)
  N  <- length(cache)
  d  <- ncol(cache[[1]]$Theta)
  M_cache <- nrow(cache[[1]]$Theta)
  p_inc   <- m_use / M_cache
  z_cut   <- qnorm(p_inc)

  # empirical moments for NIW hyper‑parameters
  theta_hat <- t(sapply(cache, function(el) colMeans(el$Theta)))
  mu_bar    <- colMeans(theta_hat)
  S_sample  <- cov(theta_hat)
  nu_p      <- N + d
  kappa_p   <- N
  Psi_p     <- S_sample * (nu_p - d - 1)

  # --- initial state ------------------------------------------------------
  mu_curr    <- mu_bar
  L_vec_curr <- vec_from_L(t(chol(S_sample)))
  Z_curr     <- replicate(N, generate_Z(M_cache), simplify = FALSE)

  loglike_curr  <- compute_log_marginals(mu_curr, L_vec_curr, Z_curr, cache, p_inc, z_cut)
  logprior_curr <- log_prior_phi(mu_curr, L_vec_curr)

  # --- dual‑averaging state ----------------------------------------------
  param_len   <- d + length(L_vec_curr)
  log_sd_block <- rep(log(0.05), param_len)
  log_sd_bar   <- log_sd_block
  H_bar        <- rep(0, param_len)
  mu_log_sd    <- log_sd_block

  # --- storage ------------------------------------------------------------
  n_keep <- length(seq(burn, n_iter, by = thin))
  keep_mu     <- matrix(NA_real_, n_keep, d)
  keep_L_vec  <- matrix(NA_real_, n_keep, length(L_vec_curr))
  keep_accept <- logical(n_iter)
  keep_rwstep <- logical(n_iter)

  param_hist  <- matrix(NA_real_, n_iter, param_len)

  iter_keep <- 0L

  for (t in 1:n_iter) {

    # -------- choose kernel ---------------------------------------------
    use_rw <- runif(1) < p_rw
    keep_rwstep[t] <- use_rw

    if (use_rw) {
      # ----- build empirical correlation R ------------------------------
      if (t > adapt_start + 5) {
        R <- cor(param_hist[1:(t - 1), , drop = FALSE], use = "pairwise")
        R[is.na(R)] <- 0; diag(R) <- 1
        R <- R + 1e-6 * diag(param_len)
        cholR <- try(chol(R), silent = TRUE)
        if (inherits(cholR, "try-error")) { R <- diag(param_len); cholR <- diag(param_len) }
      } else {
        R <- diag(param_len); cholR <- diag(param_len)
      }

      s_vec  <- exp(log_sd_block)
      cholS  <- diag(s_vec) %*% cholR
      step   <- drop(rmvnorm(1, sigma = cholS %*% t(cholS)))

      mu_prop     <- mu_curr + step[1:d]
      L_vec_prop  <- L_vec_curr + step[(d + 1):param_len]
      Z_prop      <- update_Z(Z_curr, rho_sub)

      log_q_prop <- log_q_curr <- 0  # symmetric RW cancels
    } else {
      # independence NIW move
      niw <- draw_NIW(mu_bar, nu_p, kappa_p, Psi_p)
      mu_prop    <- niw$mu
      L_vec_prop <- niw$L_vec
      Z_prop     <- update_Z(Z_curr, rho_sub)

      log_q_prop <- log_q_NIW(mu_prop, L_vec_prop, mu_bar, nu_p, kappa_p, Psi_p)
      log_q_curr <- log_q_NIW(mu_curr, L_vec_curr, mu_bar, nu_p, kappa_p, Psi_p)
    }

    # -------- likelihood & prior ---------------------------------------
    loglike_prop  <- compute_log_marginals(mu_prop, L_vec_prop, Z_prop, cache, p_inc, z_cut)
    logprior_prop <- log_prior_phi(mu_prop, L_vec_prop)

    # -------- acceptance -------------------------------------------------
    logacc <- (loglike_prop + logprior_prop + log_q_curr) -
      (loglike_curr + logprior_curr + log_q_prop)
    acc <- (log(runif(1)) < logacc)
    keep_accept[t] <- acc

    if (acc) {
      mu_curr       <- mu_prop
      L_vec_curr    <- L_vec_prop
      Z_curr        <- Z_prop
      loglike_curr  <- loglike_prop
      logprior_curr <- logprior_prop
    }

    # -------- dual‑averaging adaptation (RW only) -----------------------
    if (use_rw && (t > adapt_start)) {
      m   <- t - adapt_start
      eta <- 1 / (m + t0)
      H_bar        <- (1 - eta) * H_bar + eta * (acc_target - acc)
      log_sd_block <- mu_log_sd - (sqrt(m) / gamma_da) * H_bar
      w <- m^(-kappa_da)
      log_sd_bar <- w * log_sd_block + (1 - w) * log_sd_bar
    }

    # -------- save draws -------------------------------------------------
    param_hist[t, ] <- c(mu_curr, L_vec_curr)
    if (t >= burn && ((t - burn) %% thin == 0)) {
      iter_keep <- iter_keep + 1L
      keep_mu[iter_keep, ]    <- mu_curr
      keep_L_vec[iter_keep, ] <- L_vec_curr
    }

    # -------- diagnostics every 1000 iters ------------------------------
    if (t %% 1000 == 0) {
      cat(sprintf(
        "iter %d | acc=%.3f | RW-share=%.2f | loĝ=%.1f | effSize(μ1)=%.1f\n",
        t,
        mean(keep_accept[1:t]),
        mean(keep_rwstep[1:t]),
        loglike_curr,
        if (iter_keep > 20) effectiveSize(keep_mu[1:iter_keep, 1]) else NA_real_
      ))
    }
  }

  list(mu_draws     = keep_mu,
       L_vec_draws  = keep_L_vec,
       accept       = keep_accept,
       rw_step      = keep_rwstep)
}

# ---------------------------------------------------------------------------
#  Example application: gamma(shape, rate) person‑level likelihood
# ---------------------------------------------------------------------------
# set.seed(1)
# N <- 30; n_i <- rep(120, N)
# mu_true    <- c(log(2), log(3))
# Sigma_true <- matrix(c(.2^2, 0.025, 0.025, .3^2), 2, 2)
# THETA_true <- mvtnorm::rmvnorm(N, mu_true, Sigma_true)
#
# Y <- lapply(1:N, function(i) {
#   rgamma(n_i[i], shape = exp(THETA_true[i, 1]), rate = exp(THETA_true[i, 2]))
# })
#
# build_cache_gamma <- function(Y, M_cache = 400L) {
#   N <- length(Y)
#   cache <- vector("list", N)
#   theta_hat_mom <- matrix(NA_real_, N, 2)
#   for (i in seq_len(N)) {
#     yi <- Y[[i]]; m <- mean(yi); v <- var(yi)
#     theta_hat_mom[i, ] <- th <- c(log(m^2 / v), log(m / v))
#     s <- 0.5 / sqrt(length(yi) / 50)
#     Sig_q <- diag(rep(s^2, 2))
#     Theta_i <- mvtnorm::rmvnorm(M_cache, th, Sig_q)
#     cache[[i]] <- list(
#       Theta  = Theta_i,
#       loglik = apply(Theta_i, 1, function(th) {
#         shape <- exp(th[1]); rate <- exp(th[2])
#         sum(dgamma(yi, shape, rate, log = TRUE))
#       }),
#       logq   = mvtnorm::dmvnorm(Theta_i, th, Sig_q, log = TRUE)
#     )
#   }
#   cache
# }
#
# cache <- build_cache_gamma(Y, M_cache = 400L)
#
# # ---- run the sampler -------------------------------------------------------
# res <- run_pmmh(cache,
#                 n_iter = 6000L,
#                 burn    = 1000L,
#                 thin    = 5L,
#                 m_use   = 100L,
#                 rho_sub = 0.999)
#
# # Very brief post‑processing --------------------------------------------------
# print(head(res$mu_draws))
# cat(sprintf("\nOverall acceptance: %.3f\n", mean(res$accept)))
