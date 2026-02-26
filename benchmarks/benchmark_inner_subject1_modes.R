#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(EMC2)
})

source("smc_core.R")
source("SMC_super_fast.R")

load("single_DDM30.RData")

subject_idx <- 1L
seeds <- 123:127
M <- 5000L

emc_subj <- single[[subject_idx]][[1L]]
data_i <- emc_subj$data[[1L]]
mu_ref <- emc_subj$prior$theta_mu_mean
Sigma_ref <- emc_subj$prior$theta_mu_var

loglik_fn <- function(Theta, data) {
  EMC2:::calc_ll_manager(as.matrix(Theta), dadm = data, model = emc_subj$model, r_cores = 1L)
}

mode_overrides <- list(
  strict = list(),
  fast = list(post_adapt_n_mcmc_moves = 2L)
)

weighted_mean <- function(Theta, w) {
  w <- as.numeric(w)
  w <- w / sum(w)
  colSums(Theta * w)
}

run_one <- function(mode, seed) {
  args <- c(
    list(
      data = data_i,
      loglik_fn = loglik_fn,
      mu_ref = mu_ref,
      Sigma_ref = Sigma_ref,
      M = M,
      seed = seed,
      verbose = FALSE
    ),
    mode_overrides[[mode]]
  )

  t0 <- proc.time()[["elapsed"]]
  fit <- do.call(enhanced_smc_elite, args)
  elapsed <- proc.time()[["elapsed"]] - t0

  muw <- weighted_mean(fit$Theta, fit$w)
  out <- data.frame(
    mode = mode,
    seed = seed,
    elapsed_sec = as.numeric(elapsed),
    log_evidence = as.numeric(fit$log_evidence),
    mcse_logZ = as.numeric(fit$mcse_logZ),
    rounds = as.integer(fit$meta$rounds),
    stringsAsFactors = FALSE
  )
  for (j in seq_along(muw)) out[[paste0("theta_mean_", names(muw)[j])]] <- as.numeric(muw[j])
  out
}

rows <- list()
k <- 1L
for (mode in names(mode_overrides)) {
  for (seed in seeds) {
    cat(sprintf("Running mode=%s seed=%d ...\n", mode, seed))
    rows[[k]] <- run_one(mode, seed)
    k <- k + 1L
  }
}

res <- do.call(rbind, rows)

summ <- do.call(rbind, lapply(split(res, res$mode), function(df) {
  data.frame(
    mode = df$mode[1L],
    n = nrow(df),
    mean_elapsed_sec = mean(df$elapsed_sec),
    sd_elapsed_sec = sd(df$elapsed_sec),
    mean_log_evidence = mean(df$log_evidence),
    sd_log_evidence = sd(df$log_evidence),
    spread_log_evidence = diff(range(df$log_evidence)),
    mean_rounds = mean(df$rounds),
    sd_rounds = sd(df$rounds),
    mean_mcse_logZ = mean(df$mcse_logZ),
    stringsAsFactors = FALSE
  )
}))

out_csv <- file.path("benchmarks", "ddm_subject1_strict_vs_fast_5rep.csv")
out_rds <- file.path("benchmarks", "ddm_subject1_strict_vs_fast_5rep.rds")
write.csv(res, out_csv, row.names = FALSE)
saveRDS(list(per_run = res, summary = summ), out_rds)

cat("\nSummary:\n")
print(summ, row.names = FALSE)
cat(sprintf("\nSaved per-run CSV: %s\n", out_csv))
cat(sprintf("Saved RDS: %s\n", out_rds))
