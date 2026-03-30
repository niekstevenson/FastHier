rm(list = ls())
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[1L]))
} else {
  normalizePath("benchmarks/run_shifted_gamma_hierarchy_compare.R")
}
repo_dir <- dirname(dirname(script_path))
setwd(repo_dir)

parse_cli_args <- function(args) {
  out <- list()
  if (!length(args)) return(out)
  for (arg in args) {
    if (!startsWith(arg, "--")) next
    arg <- sub("^--", "", arg)
    parts <- strsplit(arg, "=", fixed = TRUE)[[1L]]
    key <- gsub("-", "_", parts[1L])
    value <- if (length(parts) > 1L) paste(parts[-1L], collapse = "=") else "true"
    out[[key]] <- value
  }
  out
}

arg_chr <- function(args, key, default = NULL) {
  val <- args[[key]]
  if (is.null(val) || !nzchar(val)) default else as.character(val)
}

arg_int <- function(args, key, default) {
  val <- args[[key]]
  if (is.null(val) || !nzchar(val)) return(as.integer(default))
  as.integer(val)
}

arg_lgl <- function(args, key, default = FALSE) {
  val <- args[[key]]
  if (is.null(val) || !nzchar(val)) return(isTRUE(default))
  tolower(as.character(val)) %in% c("1", "true", "t", "yes", "y")
}

config_label_default <- function(refined_method, transport_method, gss_enable, da_enable) {
  refined_tag <- if (identical(refined_method, "defensive_mixture")) "defmix" else "broad"
  transport_tag <- if (identical(transport_method, "gaussian_copula")) "gcop" else "tri"
  mode_tag <- if (gss_enable && da_enable) {
    "gss_da"
  } else if (gss_enable) {
    "gss"
  } else if (da_enable) {
    "da"
  } else {
    "base"
  }
  paste(refined_tag, transport_tag, mode_tag, sep = "_")
}

suppressPackageStartupMessages({
  library(parallel)
})

set.seed(20260324L)

cli_args <- parse_cli_args(commandArgs(trailingOnly = TRUE))

transport_method <- arg_chr(cli_args, "transport_method", "sparse_triangular")
refined_method <- arg_chr(cli_args, "refined_method", "defensive_mixture")
gss_enable <- arg_lgl(cli_args, "gss_enable", FALSE)
da_enable <- arg_lgl(cli_args, "da_enable", FALSE)
hist_mix_enable <- arg_lgl(cli_args, "hist_mix_enable", gss_enable || da_enable)
run_label <- arg_chr(cli_args, "label", NULL)

stan_results_file <- file.path("benchmarks", "samples", "shifted_gamma_hierarchy_stan_results.rds")
results_file <- arg_chr(cli_args, "results_file", file.path("benchmarks", "results", "shifted_gamma_hierarchy_current_results.rds"))
plot_file <- arg_chr(cli_args, "plot_file", file.path("benchmarks", "results", "shifted_gamma_hierarchy_current_posteriors.png"))

detected_cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
if (!is.finite(detected_cores) || detected_cores < 1L) {
  detected_cores <- 1L
}

mc.cores <- as.integer(max(1L, min(arg_int(cli_args, "mc_cores", 4L), detected_cores)))
pilot_size <- arg_int(cli_args, "pilot_size", 10L)
pilot_particles <- arg_int(cli_args, "pilot_particles", 800L)
full_particles <- arg_int(cli_args, "full_particles", 2000L)
outer_particles <- arg_int(cli_args, "outer_particles", 2000L)
outer_mcmc_moves <- arg_int(cli_args, "outer_mcmc_moves", 3L)
outer_max_rounds <- arg_int(cli_args, "outer_max_rounds", 80L)
base_seed <- arg_int(cli_args, "base_seed", 20260324L)
verbose <- arg_lgl(cli_args, "verbose", TRUE)

if (is.null(run_label)) {
  if (length(commandArgs(trailingOnly = TRUE))) {
    run_label <- config_label_default(refined_method, transport_method, gss_enable, da_enable)
  } else {
    run_label <- "current"
  }
}

if (is.null(cli_args[["results_file"]])) {
  if (!identical(run_label, "current")) {
    results_file <- file.path("benchmarks", "results", sprintf("shifted_gamma_%s_results.rds", run_label))
  }
}
if (is.null(cli_args[["plot_file"]])) {
  if (!identical(run_label, "current")) {
    plot_file <- file.path("benchmarks", "results", sprintf("shifted_gamma_%s_posteriors.png", run_label))
  }
}

dir.create(file.path("benchmarks", "samples"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path("benchmarks", "results"), showWarnings = FALSE, recursive = TRUE)

source("smc_core.R")
source("reference_priors.R")
source("utilities.R")
source("SMC_super_fast.R")
if (identical(transport_method, "gaussian_copula")) {
  fit_copula_transform <- .fit_gaussian_copula_transform
}
source("hierarchical_locals.R")
source("population_models.R")
source("outer_population_smc.R")

if (!file.exists(stan_results_file)) {
  stop("Missing Stan benchmark results: ", stan_results_file)
}

bundle <- readRDS(stan_results_file)

y <- bundle$data$y
m0 <- as.numeric(bundle$priors$m0)
s0 <- as.numeric(bundle$priors$s0)
a0 <- as.numeric(bundle$priors$a0)
b0 <- as.numeric(bundle$priors$b0)

alpha_names <- c("eta_shape", "eta_scale", "eta_shift")

data_list <- lapply(seq_len(nrow(y)), function(i) y[i, ])

names(m0) <- alpha_names
names(s0) <- alpha_names
names(a0) <- alpha_names
names(b0) <- alpha_names

base_mu <- m0
base_var <- s0 + b0 / (a0 - 1)
base_Sigma <- diag(base_var, nrow = length(alpha_names))
dimnames(base_Sigma) <- list(alpha_names, alpha_names)

loglik_shifted_gamma <- function(Theta, y_i) {
  Theta <- as.matrix(Theta)
  colnames(Theta) <- alpha_names

  eps <- 1e-9
  shape <- exp(Theta[, "eta_shape"]) + eps
  scale <- exp(Theta[, "eta_scale"]) + eps
  shift <- exp(Theta[, "eta_shift"]) + eps
  min_y <- min(y_i)

  out <- rep(-1e12, nrow(Theta))
  ok <- shift < min_y
  if (!any(ok)) {
    return(out)
  }

  for (i in which(ok)) {
    out[i] <- sum(stats::dgamma(y_i - shift[i], shape = shape[i], scale = scale[i], log = TRUE))
  }
  out[!is.finite(out)] <- -1e12
  out
}

population_model <- make_population_model_diag_gaussian(
  alpha_names = alpha_names,
  mean_prior_mean = m0,
  mean_prior_var = s0,
  sigma2_prior_shape = a0,
  sigma2_prior_rate = b0,
  label = "shifted_gamma_hierarchy"
)

stan_draws <- data.frame(
  mu_shape = bundle$draws$mu[, 1L],
  mu_scale = bundle$draws$mu[, 2L],
  mu_shift = bundle$draws$mu[, 3L],
  sigma2_shape = bundle$draws$sigma2[, 1L],
  sigma2_scale = bundle$draws$sigma2[, 2L],
  sigma2_shift = bundle$draws$sigma2[, 3L],
  check.names = FALSE
)

cat(sprintf("Loaded Stan benchmark bundle: %s\n", stan_results_file))
cat(sprintf("Data: %d subjects x %d trials\n", nrow(y), ncol(y)))
cat(sprintf("Configuration: label=%s | refined=%s | transport=%s | hist_mix=%s | gss=%s | da=%s\n",
            run_label,
            refined_method,
            transport_method,
            hist_mix_enable,
            gss_enable,
            da_enable))
cat("Running local-reference stage...\n")

stage <- prepare_reference_local_stage(
  data_list = data_list,
  loglik_fn = loglik_shifted_gamma,
  base_mu = base_mu,
  base_Sigma = base_Sigma,
  pilot_size = min(pilot_size, length(data_list)),
  broad_scale = 1,
  pilot_particles = pilot_particles,
  full_particles = full_particles,
  refined_method = refined_method,
  n_jobs = mc.cores,
  base_seed = base_seed,
  pilot_smc_control = list(
    max_rounds = 40L,
    hist_mix_enable = hist_mix_enable,
    gss_enable = gss_enable,
    da_enable = da_enable
  ),
  full_smc_control = list(
    hist_mix_enable = hist_mix_enable,
    gss_enable = gss_enable,
    da_enable = da_enable
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

workflow_theta <- smc_posteriors(
  fit,
  n_draws = nrow(as.matrix(bundle$draws$mu)),
  seed = base_seed + 1L
)

workflow_draws <- data.frame(
  mu_shape = workflow_theta$mu_eta_shape,
  mu_scale = workflow_theta$mu_eta_scale,
  mu_shift = workflow_theta$mu_eta_shift,
  sigma2_shape = exp(workflow_theta$log_sigma2_eta_shape),
  sigma2_scale = exp(workflow_theta$log_sigma2_eta_scale),
  sigma2_shift = exp(workflow_theta$log_sigma2_eta_shift),
  check.names = FALSE
)

grDevices::png(plot_file, width = 1400, height = 900)
plot_posteriors(
  stan_draws,
  workflow_draws,
  labels = c("Stan", run_label),
  cols = c("black", "firebrick3"),
  n_cols = 3L
)
grDevices::dev.off()

saveRDS(
  list(
    stan_source = stan_results_file,
    stage = stage,
    fit = fit,
    stan_draws = stan_draws,
    workflow_draws = workflow_draws,
    settings = list(
      label = run_label,
      transport_method = transport_method,
      refined_method = refined_method,
      hist_mix_enable = hist_mix_enable,
      gss_enable = gss_enable,
      da_enable = da_enable,
      mc.cores = mc.cores,
      pilot_size = min(pilot_size, length(data_list)),
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
