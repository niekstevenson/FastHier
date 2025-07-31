rm(list = ls())
# ----------------------------------------------------
#
#  * Works for ANY parameter‑dimension p  (taken from the cache)
#  * Still mixes an independence‑NIW kernel with an adaptive RW kernel
#  * Correlated‑PM machinery kept untouched, but rho_sub and the
#    Bernoulli‑thinning rate are now explicit arguments
#  * Prints progress every diag_print iterations
#
# A minimal example of use sits at the bottom of the file.
# ----------------------------------------------------

library(coda)
library(EMC2)
library(mclust)      # EM for Gaussian mixtures
library(mvtnorm)     # dmvnorm / rmvnorm
library(loo)         # PSIS diagnostics
library(matrixStats) # fast ESS helpers

#' Fit a Gaussian–mixture proposal q(θ)
#' @param chain  n × p matrix/data-frame of MCMC draws
#' @param G      number of mixture components; if NULL let Mclust pick by BIC
#' @returns      list with log_q(⋅) and r_q(⋅) samplers
fit_gmm_proposal <- function(chain, G = NULL) {
  chain <- as.matrix(chain)
  gmm   <- Mclust(chain, G = G)              # EM fit

  # helper to compute mixture density
  log_q <- function(theta) {
    theta <- as.matrix(theta)
    ll <- rep(0, nrow(theta))
    pars <- gmm$parameters
    for (g in seq_len(gmm$G)) {
      ll <- ll + pars$pro[g] *
        mvtnorm::dmvnorm(theta,
                         mean  = pars$mean[, g],
                         sigma = pars$variance$sigma[, , g])
    }
    log(ll)
  }

  # helper to draw N samples
  r_q <- function(N) {
    comp <- sample(seq_len(gmm$G), N, replace = TRUE, prob = gmm$parameters$pro)
    out  <- matrix(NA_real_, N, ncol(chain))
    for (g in seq_len(gmm$G)) {
      idx <- which(comp == g)
      if (length(idx))
        out[idx, ] <- mvtnorm::rmvnorm(length(idx),
                                       mean  = gmm$parameters$mean[, g],
                                       sigma = gmm$parameters$variance$sigma[, , g])
    }
    colnames(out) <- colnames(chain)
    out
  }

  list(model = gmm, log_q = log_q, r_q = r_q)
}


#' Importance sampling with PSIS smoothing
#' @param data      the data set your likelihood needs
#' @param ll        function(data, theta) → log-likelihood  **(vectorised over rows of theta)**
#' @param proposal  list returned by fit_gmm_proposal()
#' @param N_imp     number of importance draws
#' @returns         θ, smoothed weights, raw ESS, PSIS k-hat
importance_sample <- function(data, ll, model, proposal, N_imp = 3000) {

  theta  <- proposal$r_q(N_imp)              # θ ~ q
  log_p  <- EMC2:::calc_ll_manager(theta, data, model)                  # vectorised log-likelihood
  log_q  <- proposal$log_q(theta)
  log_w  <- log_p - log_q

  ## Pareto-smoothed importance weights
  psis   <- loo::psis(log_w)
  w      <- loo::weights.importance_sampling(psis, normalize = TRUE)

  ess    <- (sum(w)^2) / sum(w^2)           # Kish effective sample size
  list(theta = theta,
       w     = w,
       ess   = ess,
       k_hat = psis$diagnostics$pareto_k)
}



# ---------- generic helpers -------------------------------------------------
logsumexp <- function(x) { m <- max(x); m + log(sum(exp(x - m))) }

# ----- Cholesky factor utilities -------------------------------------------
vec_to_L <- function(param_vec, p) {
  diag_logs <- param_vec[seq_len(p)]
  off       <- param_vec[-seq_len(p)]
  L <- matrix(0, p, p)
  diag(L) <- exp(diag_logs)
  idx <- 1L
  for (j in 1:(p - 1)) {
    for (i in (j + 1):p) {
      L[i, j] <- off[idx]
      idx <- idx + 1L
    }
  }
  L
}

L_to_vec <- function(L) {
  p <- nrow(L)
  c(log(diag(L)), L[lower.tri(L)])
}

Sigma_from_L <- function(L) L %*% t(L)

# ---------- prior densities -------------------------------------------------
log_prior_mu <- function(mu, m0 = NULL, s0 = NULL) {
  if (is.null(m0)) m0 <- rep(0, length(mu))
  if (is.null(s0)) s0 <- rep(3, length(mu))
  sum(dnorm(mu, m0, s0, log = TRUE))
}

log_prior_L <- function(L, meanlog_diag = 0, sdlog_diag = .5, sd_off = 1) {
  p <- nrow(L)
  diag_logs <- log(diag(L))
  off <- L[lower.tri(L)]
  sum(dlnorm(exp(diag_logs), meanlog_diag, sdlog_diag, log = TRUE) + diag_logs) +
    sum(dnorm(off, 0, sd_off, log = TRUE))
}

log_prior_phi <- function(mu, L) log_prior_mu(mu) + log_prior_L(L)

# ---------- correlated‑PM helpers -------------------------------------------
log_marginals_from_Z <- function(mu, L, cache, Z_list, z_cut, p_inc) {
  N <- length(cache)
  Sig <- Sigma_from_L(L)
  M_cache <- nrow(cache[[1]]$Theta)
  total <- 0
  for (i in seq_len(N)) {
    inc <- Z_list[[i]] < z_cut
    if (!any(inc)) {           # guarantee ≥1 inclusion
      inc[which.min(Z_list[[i]])] <- TRUE
    }
    idx <- which(inc)
    logw <- dmvnorm(cache[[i]]$Theta[idx, , drop = FALSE], mu, Sig, log = TRUE) -
      cache[[i]]$logq[idx] +
      cache[[i]]$loglik[idx] + log(1 / p_inc)
    total <- total + (-log(M_cache) + logsumexp(logw))
  }
  total
}

generate_Z <- function(M_cache) rnorm(M_cache)

update_Z <- function(Z_list, rho_sub) {
  lapply(Z_list, function(z) rho_sub * z + sqrt(1 - rho_sub^2) * rnorm(length(z)))
}

# ---------- NIW proposal ----------------------------------------------------
rinvwishart <- function(df, Psi) {
  Wi <- rWishart(1, df, solve(Psi))[,,1]
  solve(Wi)
}

log_multigamma <- function(a, p) {
  sum(lgamma(a + (1 - seq_len(p)) / 2)) + (p * (p - 1) / 4) * log(pi)
}

log_dinvwishart <- function(Sigma, df, Psi) {
  p <- nrow(Sigma)
  term1 <- - (df + p + 1) / 2 * determinant(Sigma, logarithm = TRUE)$modulus
  term2 <- -0.5 * sum(diag(Psi %*% solve(Sigma)))
  term3 <- -df * p / 2 * log(2) - 0.5 * determinant(Psi, logarithm = TRUE)$modulus -
    log_multigamma(df / 2, p)
  as.numeric(term1 + term2 + term3)
}

draw_from_NIW <- function(mu_bar, S_sample, nu_p, kappa_p) {
  p <- length(mu_bar)
  Sig <- rinvwishart(nu_p, S_sample * (nu_p - p - 1))
  mu  <- rmvnorm(1, mu_bar, Sig / kappa_p)[1, ]
  list(mu = mu, L = t(chol(Sig)))
}

# ---------- RW‑proposal helpers --------------------------------------------
build_rw_step <- function(param_vec, log_sd_block, cholR) {
  s_vec <- exp(log_sd_block)
  cholS <- diag(s_vec) %*% cholR
  step  <- drop(rmvnorm(1, sigma = cholS %*% t(cholS)))
  param_vec + step
}

# ---------- main PMMH driver -----------------------------------------------
run_pmmh <- function(cache,
                     n_iter      = 6000L,
                     burn        = 1000L,
                     thin        = 5L,
                     p_rw        = 0.9,
                     rho_sub     = 0.999,
                     m_use       = 100L,
                     adapt_start = 100L,
                     acc_target  = 0.40,
                     gamma_da    = 0.05,
                     t0          = 10,
                     kappa_da    = 0.75,
                     diag_print  = 1000L,
                     seed        = NULL) {

  if (!is.null(seed)) set.seed(seed)

  # --- dimensions & constants --------------------------------------------
  N <- length(cache)
  p <- ncol(cache[[1]]$Theta)
  M_cache <- nrow(cache[[1]]$Theta)
  p_param <- p + p * (p - 1) / 2           # vector length for (diag+off‑diag)
  m_use   <- as.integer(m_use)
  p_inc   <- m_use / M_cache
  z_cut   <- qnorm(p_inc)

  # --- NIW hyper‑params ----------------------------------------------------
  theta_hat_mom <- t(vapply(cache, function(li) colMeans(li$Theta), numeric(p)))
  mu_bar   <- colMeans(theta_hat_mom)
  S_sample <- cov(theta_hat_mom)
  nu_p     <- N + p
  kappa_p  <- N

  # --- initial state -------------------------------------------------------
  Sig_init <- S_sample
  L_curr   <- t(chol(Sig_init))
  mu_curr  <- mu_bar
  Z_curr   <- replicate(N, generate_Z(M_cache), simplify = FALSE)
  loglike_curr  <- log_marginals_from_Z(mu_curr, L_curr, cache, Z_curr, z_cut, p_inc)
  logprior_curr <- log_prior_phi(mu_curr, L_curr)

  # --- dual‑averaging state -------------------------------------------------
  log_sd_block <- rep(log(0.05), p_param + p)  # (mu, Lvec)
  log_sd_bar   <- log_sd_block
  H_bar        <- rep(0, p_param + p)
  mu_log_sd    <- log_sd_block

  # --- storage -------------------------------------------------------------
  draws_needed <- floor((n_iter - burn) / thin)
  keep <- list(
    mu          = matrix(NA_real_, draws_needed, p),
    L           = array(NA_real_, c(p, p, draws_needed)),
    loglike_hat = numeric(draws_needed)
  )
  accept  <- logical(n_iter)
  rw_step <- logical(n_iter)
  param_hist <- matrix(NA_real_, n_iter, p_param + p)
  draw_idx <- 0L

  # --- main loop -----------------------------------------------------------
  for (t in seq_len(n_iter)) {

    ## 1. NIW independence proposal --------------------------------------
    prop       <- draw_from_NIW(mu_bar, S_sample, nu_p, kappa_p)
    mu_prop    <- prop$mu
    L_prop     <- prop$L
    Z_prop     <- update_Z(Z_curr, rho_sub)
    loglike_p  <- log_marginals_from_Z(mu_prop, L_prop, cache, Z_prop, z_cut, p_inc)
    logprior_p <- log_prior_phi(mu_prop, L_prop)

    # proposal densities (needed because NIW is independence kernel)
    q_curr <- log_dinvwishart(Sigma_from_L(L_curr), nu_p, S_sample * (nu_p - p - 1)) +
      dmvnorm(mu_curr, mu_bar, Sigma_from_L(L_curr) / kappa_p, log = TRUE)
    q_prop <- log_dinvwishart(Sigma_from_L(L_prop), nu_p, S_sample * (nu_p - p - 1)) +
      dmvnorm(mu_prop, mu_bar, Sigma_from_L(L_prop) / kappa_p, log = TRUE)

    logacc_ind <- (loglike_p + logprior_p + q_curr) - (loglike_curr + logprior_curr + q_prop)

    ## 2. potential RW proposal -------------------------------------------
    use_rw <- runif(1) < p_rw
    rw_step[t] <- use_rw

    if (use_rw) {
      # empirical correlation matrix
      if (t > adapt_start + 5) {
        R <- cor(param_hist[1:(t - 1), , drop = FALSE], use = "pairwise")
        R[is.na(R)] <- 0; diag(R) <- 1
        R <- R + 1e-6 * diag(ncol(R))
        cholR <- try(chol(R), silent = TRUE)
        if (inherits(cholR, "try-error")) cholR <- diag(ncol(R))
      } else {
        cholR <- diag(p_param + p)
      }

      param_vec_curr <- c(mu_curr, L_to_vec(L_curr))
      param_vec_prop <- build_rw_step(param_vec_curr, log_sd_block, cholR)

      mu_prop  <- param_vec_prop[seq_len(p)]
      L_prop   <- vec_to_L(param_vec_prop[-seq_len(p)], p)
      Z_prop   <- update_Z(Z_curr, rho_sub)
      loglike_p  <- log_marginals_from_Z(mu_prop, L_prop, cache, Z_prop, z_cut, p_inc)
      logprior_p <- log_prior_phi(mu_prop, L_prop)

      logacc <- (loglike_p + logprior_p) - (loglike_curr + logprior_curr)  # symmetric
    } else {
      logacc <- logacc_ind
    }

    ## 3. MH accept/reject --------------------------------------------------
    acc <- (log(runif(1)) < logacc)
    accept[t] <- acc

    if (acc) {
      mu_curr       <- mu_prop
      L_curr        <- L_prop
      Z_curr        <- Z_prop
      loglike_curr  <- loglike_p
      logprior_curr <- logprior_p
    }

    ## 4. dual‑averaging update (only RW steps) ----------------------------
    if (use_rw && (t > adapt_start)) {
      m   <- t - adapt_start
      eta <- 1 / (m + t0)
      H_bar <- (1 - eta) * H_bar + eta * (acc_target - acc)
      log_sd_block <- mu_log_sd - (sqrt(m) / gamma_da) * H_bar
      w <- m^(-kappa_da)
      log_sd_bar <- w * log_sd_block + (1 - w) * log_sd_bar
    }

    ## 5. bookkeeping ------------------------------------------------------
    param_hist[t, ] <- c(mu_curr, L_to_vec(L_curr))

    if (t > burn && ((t - burn) %% thin == 0)) {
      draw_idx <- draw_idx + 1L
      keep$mu[draw_idx, ]        <- mu_curr
      keep$L[,, draw_idx]        <- L_curr
      keep$loglike_hat[draw_idx] <- loglike_curr
    }

    if (t %% diag_print == 0L) {
      cat(sprintf("iter %d | acc=%.3f | RW‑share=%.2f | loglikê=%.1f\n",
                  t, mean(accept[1:t]), mean(rw_step[1:t]), loglike_curr))
    }
  }

  list(draws   = keep,
       accept  = accept,
       rw_step = rw_step,
       p_inc   = p_inc,
       rho_sub = rho_sub)
}

# ---------------------------- Example usage ---------------------------------
# The following illustrates how you might call the sampler on the toy Gamma
# example that was in the original script.  Commented‑out so sourcing this file
# does not run anything by default.

set.seed(1)
# N <- 30; n_i <- rep(120, N)
# mu_true    <- c(log(2), log(3))
# Sigma_true <- matrix(c(.2^2, 0.025, 0.025, .3^2), 2, 2)
# THETA_true <- rmvnorm(N, mu_true, Sigma_true)
# Y <- lapply(1:N, function(i) rgamma(n_i[i],
#                                     shape = exp(THETA_true[i, 1]),
#                                     rate  = exp(THETA_true[i, 2])))

library(EMC2)
load("~/Documents/2025/TwoStep/samples/single.RData")

## ------------------------------------------------------------------
## 1.  Utility that turns *one* element of the single list
##     into the structure the PMMH sampler expects
## ------------------------------------------------------------------
make_is_cache <- function(single_subject,
                          G        = 5,       # GMM components for q(θ)
                          M_cache  = 1000L) { # particles per subject

  ## 1a. Posterior draws of subject-level parameters (α’s)
  chain <- parameters(single_subject, selection = "alpha")
  chain <- chain[, -1, drop = FALSE]          # drop iteration index if present

  ## 1b. Fit Gaussian-mixture proposal  q(θ)
  prop   <- fit_gmm_proposal(chain, G = G)

  ## 1c. Draw Θ₁:ₘ  ~  q(θ)
  Theta  <- prop$r_q(M_cache)

  ## 1d. Exact log-likelihood *for each particle*
  mdl    <- single_subject[[1]]$model
  dat    <- single_subject[[1]]$data[[1]]
  loglik <- EMC2:::calc_ll_manager(Theta, dat, mdl)  # vectorised over rows of Θ

  ## 1e. Proposal log-density  log q(Θ)
  logq   <- prop$log_q(Theta)

  list(Theta  = Theta,
       loglik = loglik,
       logq   = logq)
}

## ------------------------------------------------------------------
## 2.  Build the cache list for *all* subjects in one line
## ------------------------------------------------------------------
M_cache <- 1000L
cache   <- lapply(single, make_is_cache, G = 5, M_cache = M_cache)

# ## Optional: quick diagnostics
# rel_ess <- sapply(single, \(subj) {
#   ch   <- parameters(subj, selection = "alpha")[, -1, drop = FALSE]
#   prop <- fit_gmm_proposal(ch, G = 5)
#   is   <- importance_sample(subj[[1]]$data[[1]],
#                             ll    = subj[[1]]$model$ll,
#                             model = subj[[1]]$model,
#                             proposal = prop,
#                             N_imp = 3000)
#   is$ess / length(is$w)
# })
# print(rel_ess)        # relative ESS per subject

## ------------------------------------------------------------------
## 3.  Run the fully-correlated PMMH sampler
## ------------------------------------------------------------------
res <- run_pmmh(cache,
                n_iter = 5000,
                burn   = 1000,
                thin   = 5,
                seed   = 1234)

sel <- seq_len(length(res$draws$loglike_hat))
matplot(res$draws$mu[sel, ], type = "l")
plot(res$draws$L[2,1,sel], type = "l")
