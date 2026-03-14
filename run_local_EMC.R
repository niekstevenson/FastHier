rm(list = ls())
library(EMC2)
library(parallel)


env_chr <- function(name, default) {
  val <- Sys.getenv(name, unset = "")
  if (!nzchar(val)) default else val
}

env_int <- function(name, default) {
  val <- Sys.getenv(name, unset = "")
  if (!nzchar(val)) default else as.integer(val)
}

weighted_quantile <- function(x, w, probs = c(0.025, 0.5, 0.975)) {
  ord <- order(x)
  x_ord <- x[ord]
  w_ord <- pmax(as.numeric(w[ord]), 0)
  sw <- sum(w_ord)
  if (!is.finite(sw) || sw <= 0) {
    w_ord <- rep(1 / length(w_ord), length(w_ord))
  } else {
    w_ord <- w_ord / sw
  }
  cw <- cumsum(w_ord)
  as.numeric(stats::approx(cw, x_ord, xout = probs, rule = 2)$y)
}

summarize_recovery <- function(true_values, smc_res) {
  Theta <- as.matrix(smc_res$Theta)
  w <- pmax(as.numeric(smc_res$w), 0)
  sw <- sum(w)
  if (!is.finite(sw) || sw <= 0) {
    w <- rep(1 / nrow(Theta), nrow(Theta))
  } else {
    w <- w / sw
  }

  param_names <- colnames(Theta)
  if (is.null(param_names)) {
    param_names <- paste0("theta", seq_len(ncol(Theta)))
    colnames(Theta) <- param_names
  }

  if (is.matrix(true_values) || is.data.frame(true_values)) {
    ref_mat <- as.matrix(true_values)
    if (!is.null(colnames(ref_mat))) {
      ref_mat <- ref_mat[, param_names, drop = FALSE]
    }
    ref_q <- t(apply(ref_mat, 2L, stats::quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))
    out <- data.frame(
      parameter = param_names,
      posterior_mean = colSums(Theta * w),
      posterior_median = vapply(seq_len(ncol(Theta)), function(j) weighted_quantile(Theta[, j], w)[2L], numeric(1)),
      posterior_ci_lower = vapply(seq_len(ncol(Theta)), function(j) weighted_quantile(Theta[, j], w)[1L], numeric(1)),
      posterior_ci_upper = vapply(seq_len(ncol(Theta)), function(j) weighted_quantile(Theta[, j], w)[3L], numeric(1)),
      reference_median = ref_q[, 2L],
      reference_ci_lower = ref_q[, 1L],
      reference_ci_upper = ref_q[, 3L],
      row.names = NULL,
      check.names = FALSE
    )
    return(out)
  }

  truth <- as.numeric(true_values[param_names])
  post_q <- t(vapply(seq_len(ncol(Theta)), function(j) weighted_quantile(Theta[, j], w), numeric(3)))
  data.frame(
    parameter = param_names,
    truth = truth,
    posterior_mean = colSums(Theta * w),
    posterior_median = post_q[, 2L],
    posterior_ci_lower = post_q[, 1L],
    posterior_ci_upper = post_q[, 3L],
    covered = truth >= post_q[, 1L] & truth <= post_q[, 3L],
    row.names = NULL,
    check.names = FALSE
  )
}

data_path <- env_chr("FH_LOCAL_EMC_DATA", "single_DDM30.RData")
save_path <- env_chr("FH_LOCAL_EMC_SAVE", "local_emc_results.rds")
mc_cores <- env_int("FH_LOCAL_EMC_MC_CORES", 12L)
M_particles <- env_int("FH_LOCAL_EMC_M", 4000L)
seed_base <- env_int("FH_LOCAL_EMC_SEED", 123L)
verbose <- identical(toupper(env_chr("FH_LOCAL_EMC_VERBOSE", "TRUE")), "TRUE")

if (!file.exists(data_path)) {
  stop(
    "Missing data file at ", data_path, ". ",
    "Set FH_LOCAL_EMC_DATA to the path of single_DDM30.RData."
  )
}

source("SMC_super_fast.R")

load(data_path)
if (!exists("single", inherits = FALSE)) {
  stop("The data file does not define an object named 'single'.")
}

emc <- single[[1L]]
mu_ref <- emc[[1L]]$prior$theta_mu_mean
Sigma_ref <- emc[[1L]]$prior$theta_mu_var
data_list <- lapply(single, function(x) x[[1L]]$data[[1L]])
model_factory <- emc[[1L]]$model
param_names <- names(mu_ref)

loglik_emc2 <- function(Theta, data_i) {
  Theta <- as.matrix(Theta)
  if (is.null(colnames(Theta))) {
    colnames(Theta) <- param_names
  }
  as.numeric(EMC2:::calc_ll_manager(Theta, data_i, model_factory, r_cores = 1L))
}

run_one_subject <- function(data_i, subject_idx) {
  enhanced_smc_elite(
    data = data_i,
    loglik_fn = loglik_emc2,
    mu_ref = mu_ref,
    Sigma_ref = Sigma_ref,
    M = M_particles,
    n_cores = 1L,
    seed = seed_base + as.integer(subject_idx) - 1L,
    verbose = verbose
  )
}

smc <- mclapply(
  seq_along(data_list),
  function(i) run_one_subject(data_list[[i]], i),
  mc.cores = mc_cores
)

recovery_stats <- lapply(seq_along(smc), function(i) {
  truth_i <- parameters(single[[i]], selection = "alpha")
  truth_i <- truth_i[, -1L, drop = FALSE]
  summarize_recovery(truth_i, smc[[i]])
})

out <- list(
  config = list(
    data_path = normalizePath(data_path, winslash = "/", mustWork = TRUE),
    save_path = save_path,
    mc_cores = mc_cores,
    M_particles = M_particles,
    seed_base = seed_base,
    verbose = verbose
  ),
  smc = smc,
  recovery_stats = recovery_stats
)

saveRDS(out, save_path)

cat("Saved local EMC results to:", save_path, "\n")
cat("Subjects run:", length(smc), "\n")
for (i in seq_along(recovery_stats)) {
  cat("\nSubject", i, "recovery summary:\n")
  print(recovery_stats[[i]], row.names = FALSE)
}
