rm(list = ls())

suppressPackageStartupMessages({
  library(parallel)
})

set.seed(20260323L)

stan_results_file <- "old/samples/normal_normal_meanvar_hierarchy_compare_results.rds"
results_file <- "normal_normal_meanvar_hierarchy_current_results.rds"
plot_file <- "normal_normal_meanvar_hierarchy_current_posteriors.png"

detected_cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
if (!is.finite(detected_cores) || detected_cores < 1L) {
  detected_cores <- 1L
}

mc.cores <- as.integer(max(1L, min(4L, detected_cores)))
pilot_size <- 10L
pilot_particles <- 800L
full_particles <- 2000L
outer_particles <- 2000L
outer_mcmc_moves <- 3L
outer_max_rounds <- 80L
base_seed <- 20260323L
verbose <- TRUE

source("hierarchical_locals.R")
source("population_models.R")
source("outer_population_smc.R")

if (!file.exists(stan_results_file)) {
  stop("Missing Stan benchmark results: ", stan_results_file)
}

bundle <- readRDS(stan_results_file)

y <- bundle$simulation$y
stan_draws <- bundle$stan_draws
phi_prior_mean <- bundle$prior$mean
phi_prior_sd <- bundle$prior$sd

S <- nrow(y)
data_list <- lapply(seq_len(S), function(i) y[i, ])

alpha_names <- c("subject_mean", "subject_log_var")
hyper_names <- c("mu_mean", "mu_log_var", "log_tau2_mean", "log_tau2_log_var")

base_mu <- c(
  subject_mean = unname(phi_prior_mean["mu_mean"]),
  subject_log_var = unname(phi_prior_mean["mu_log_var"])
)
base_Sigma <- structure(
  diag(c(4.0, 2.0), 2L),
  dimnames = list(alpha_names, alpha_names)
)

normalize_weights <- function(w) {
  w <- pmax(as.numeric(w), 0)
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) {
    rep(1 / length(w), length(w))
  } else {
    w / sw
  }
}

resample_weighted_rows <- function(x, w, n_draws, seed) {
  set.seed(as.integer(seed))
  x <- as.matrix(x)
  idx <- sample.int(
    nrow(x),
    size = as.integer(n_draws),
    replace = TRUE,
    prob = normalize_weights(w)
  )
  x[idx, , drop = FALSE]
}

theta_draws_to_named_list <- function(theta_draws) {
  theta_draws <- as.matrix(theta_draws)
  colnames(theta_draws) <- hyper_names
  list(
    mu_mean = theta_draws[, "mu_mean"],
    mu_log_var = theta_draws[, "mu_log_var"],
    tau2_mean = exp(theta_draws[, "log_tau2_mean"]),
    tau2_log_var = exp(theta_draws[, "log_tau2_log_var"])
  )
}

plot_density_overlay <- function(stan_x, workflow_x, main, xlab) {
  d_stan <- stats::density(stan_x)
  d_workflow <- stats::density(workflow_x)
  xlim <- range(c(d_stan$x, d_workflow$x))
  ylim <- c(0, 1.05 * max(d_stan$y, d_workflow$y))

  plot(
    d_stan,
    lwd = 2,
    col = "black",
    xlim = xlim,
    ylim = ylim,
    main = main,
    xlab = xlab,
    ylab = "Density"
  )
  lines(d_workflow, lwd = 2, col = "firebrick3")
  legend(
    "topright",
    legend = c("Stan", "Current workflow"),
    col = c("black", "firebrick3"),
    lwd = 2,
    bty = "n",
    cex = 0.85
  )
}

loglik_fn <- function(Theta, y_i) {
  Theta <- as.matrix(Theta)
  colnames(Theta) <- alpha_names
  vapply(
    seq_len(nrow(Theta)),
    function(i) {
      mu_i <- Theta[i, "subject_mean"]
      sd_i <- pmax(exp(0.5 * Theta[i, "subject_log_var"]), 1e-9)
      sum(stats::dnorm(y_i, mean = mu_i, sd = sd_i, log = TRUE))
    },
    numeric(1L)
  )
}

make_population_model_normal_logvar <- function(prior_mean, prior_sd) {
  prior_mean <- as.numeric(prior_mean[hyper_names])
  prior_sd <- as.numeric(prior_sd[hyper_names])

  sample_hyper <- function(n) {
    out <- matrix(
      stats::rnorm(
        n * length(hyper_names),
        mean = rep(prior_mean, each = n),
        sd = rep(prior_sd, each = n)
      ),
      nrow = n,
      ncol = length(hyper_names),
      byrow = FALSE
    )
    colnames(out) <- hyper_names
    out
  }

  log_hyperprior <- function(theta) {
    theta <- as.matrix(theta)
    colnames(theta) <- hyper_names
    rowSums(
      vapply(
        seq_along(hyper_names),
        function(j) {
          stats::dnorm(theta[, j], mean = prior_mean[j], sd = prior_sd[j], log = TRUE)
        },
        numeric(nrow(theta))
      )
    )
  }

  log_alpha_given_theta <- function(alpha, theta_row) {
    alpha <- as.matrix(alpha)
    colnames(alpha) <- alpha_names
    theta_row <- as.matrix(theta_row)
    colnames(theta_row) <- hyper_names

    mu <- c(
      subject_mean = theta_row[1L, "mu_mean"],
      subject_log_var = theta_row[1L, "mu_log_var"]
    )
    sigma2 <- c(
      subject_mean = exp(theta_row[1L, "log_tau2_mean"]),
      subject_log_var = exp(theta_row[1L, "log_tau2_log_var"])
    )

    centered <- sweep(alpha, 2L, mu, "-")
    as.numeric(
      -0.5 * (length(alpha_names) * log(2 * pi) + sum(log(sigma2)) +
                rowSums(sweep(centered^2, 2L, sigma2, "/")))
    )
  }

  log_alpha_given_theta_many <- function(alpha, theta) {
    alpha <- as.matrix(alpha)
    colnames(alpha) <- alpha_names
    theta <- as.matrix(theta)
    colnames(theta) <- hyper_names

    mu_mean <- theta[, "mu_mean"]
    mu_log_var <- theta[, "mu_log_var"]
    tau2_mean <- exp(theta[, "log_tau2_mean"])
    tau2_log_var <- exp(theta[, "log_tau2_log_var"])

    out <- matrix(
      -0.5 * log(2 * pi * tau2_mean) - 0.5 * log(2 * pi * tau2_log_var),
      nrow = nrow(theta),
      ncol = nrow(alpha)
    )
    out <- out - 0.5 * sweep(outer(mu_mean, alpha[, "subject_mean"], "-")^2, 1L, tau2_mean, "/")
    out <- out - 0.5 * sweep(outer(mu_log_var, alpha[, "subject_log_var"], "-")^2, 1L, tau2_log_var, "/")
    out
  }

  normalize_population_model(
    list(
      name = "normal_normal_meanvar_hierarchy",
      alpha_dim = length(alpha_names),
      hyper_dim = length(hyper_names),
      alpha_names = alpha_names,
      hyper_names = hyper_names,
      sample_hyper = sample_hyper,
      log_hyperprior = log_hyperprior,
      log_alpha_given_theta = log_alpha_given_theta,
      log_alpha_given_theta_many = log_alpha_given_theta_many
    )
  )
}

population_model <- make_population_model_normal_logvar(
  prior_mean = phi_prior_mean,
  prior_sd = phi_prior_sd
)

cat(sprintf("Loaded Stan benchmark bundle: %s\n", stan_results_file))
cat(sprintf("Data: %d subjects x %d trials\n", nrow(y), ncol(y)))
cat("Running local-reference stage...\n")

stage <- prepare_reference_local_stage(
  data_list = data_list,
  loglik_fn = loglik_fn,
  base_mu = base_mu,
  base_Sigma = base_Sigma,
  pilot_size = min(pilot_size, S),
  broad_scale = 1,
  broad_defensive = FALSE,
  pilot_particles = pilot_particles,
  full_particles = full_particles,
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
    max_rounds = 40L,
    n_mcmc_moves = 1L,
    G_mix = 8L,
    hist_mix_enable = FALSE,
    gss_enable = FALSE,
    da_enable = FALSE
  ),
  full_smc_control = list(
    hist_mix_enable = FALSE,
    gss_enable = FALSE,
    da_enable = FALSE
  )
)

cat("Running outer population SMC...\n")

factor_set <- build_population_factor_set(stage$local_objects, population_model)
fit <- outer_population_smc(
  factor_set = factor_set,
  N = outer_particles,
  n_mcmc_moves = outer_mcmc_moves,
  max_rounds = outer_max_rounds,
  n_cores = mc.cores,
  seed = base_seed,
  verbose = verbose
)

workflow_draws <- theta_draws_to_named_list(
  resample_weighted_rows(
    fit$theta,
    fit$w,
    n_draws = length(stan_draws$mu_mean),
    seed = base_seed + 1L
  )
)

png(plot_file, width = 1200, height = 900)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
plot_density_overlay(stan_draws$mu_mean, workflow_draws$mu_mean, "mu_mean", "mu_mean")
plot_density_overlay(stan_draws$tau2_mean, workflow_draws$tau2_mean, "tau2_mean", "tau2_mean")
plot_density_overlay(stan_draws$mu_log_var, workflow_draws$mu_log_var, "mu_log_var", "mu_log_var")
plot_density_overlay(stan_draws$tau2_log_var, workflow_draws$tau2_log_var, "tau2_log_var", "tau2_log_var")
dev.off()

saveRDS(
  list(
    stan_source = stan_results_file,
    stage = stage,
    fit = fit,
    stan_draws = stan_draws,
    workflow_draws = workflow_draws,
    settings = list(
      mc.cores = mc.cores,
      pilot_size = min(pilot_size, S),
      pilot_particles = pilot_particles,
      full_particles = full_particles,
      outer_particles = outer_particles,
      outer_mcmc_moves = outer_mcmc_moves,
      outer_max_rounds = outer_max_rounds,
      base_seed = base_seed
    ),
    plot_file = plot_file
  ),
  results_file
)

cat("Saved results to:", results_file, "\n")
cat("Saved plot to:", plot_file, "\n")
