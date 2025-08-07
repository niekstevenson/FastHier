# minimal deps
library(mvtnorm)
library(loo)
library(coda)
library(ggplot2)

logsumexp <- function(x) { m <- max(x); m + log(sum(exp(x - m))) }
kish_ess  <- function(w) { s <- sum(w); (s*s)/sum(w*w) }

# 1) Build posterior weights for one subject at given phi=(mu,Sigma)
subject_posterior_weights <- function(cache, mu, Sigma) {
  stopifnot(!is.null(cache$Theta), !is.null(cache$loglik))
  den <- if (!is.null(cache$logqmix)) cache$logqmix else cache$logq
  stopifnot(!is.null(den))
  r <- cache$loglik + mvtnorm::dmvnorm(cache$Theta, mean = mu, sigma = Sigma, log = TRUE) - den
  w <- exp(r - logsumexp(r))                       # self-normalized weights (for inspection)
  list(logw_raw = r, w = w)
}

# 2) Quick diagnostics (ESS and PSIS k-hat)
subject_is_diagnostics <- function(cache, mu, Sigma) {
  ww <- subject_posterior_weights(cache, mu, Sigma)
  ess  <- kish_ess(ww$w)
  rel_ess <- ess / length(ww$w)
  ps <- loo::psis(ww$logw_raw)                     # PSIS diagnostic (for assessment only)
  k_hat <- ps$diagnostics$pareto_k
  k_thr <- min(1 - 1/log10(length(ww$w)), 0.7)     # sample-size-specific threshold
  list(ESS = ess, rel_ESS = rel_ess, k_hat = k_hat, k_thr = k_thr)
}

# 3) Get posterior draws for plotting/summaries (simple multinomial resample)
subject_posterior_draws <- function(cache, mu, Sigma, R = 4000L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  ww <- subject_posterior_weights(cache, mu, Sigma)
  idx <- sample.int(nrow(cache$Theta), R, replace = TRUE, prob = ww$w)
  draws <- cache$Theta[idx, , drop = FALSE]
  attr(draws, "weights") <- ww$w                     # keep around if needed
  draws
}

# 4) Summaries + HPD intervals
subject_posterior_summary <- function(draws, prob = 0.95) {
  # resampled draws -> ordinary summaries
  sm <- as.data.frame(cbind(mean = colMeans(draws),
                            t(apply(draws, 2, median)),
                            t(apply(draws, 2, sd))))
  colnames(sm) <- c("mean","median","sd")
  # HPD via coda
  hpd <- coda::HPDinterval(coda::as.mcmc(draws), prob = prob)
  cbind(sm, hpd)
}

# 5) Plots (marginals, and overlay true theta if provided)
plot_subject_posterior <- function(draws, theta_true = NULL, names = NULL) {
  D <- as.data.frame(draws)
  if (!is.null(names)) colnames(D) <- names
  Dlong <- reshape2::melt(D, variable.name = "param", value.name = "value")
  p <- ggplot(Dlong, aes(x = value)) +
    geom_histogram(aes(y = after_stat(density)), bins = 80) +
    facet_wrap(~ param, scales = "free", ncol = 2) +
    theme_bw() + labs(title = "Subject posterior (IS resampled)")
  if (!is.null(theta_true)) {
    tf <- data.frame(param = colnames(D), value = as.numeric(theta_true))
    p <- p + geom_vline(data = tf, aes(xintercept = value), linetype = 2)
  }
  p
}
