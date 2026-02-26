rm(list = ls())
library(EMC2)
# library(cmdstanr)
# load("~/Documents/2025/TwoStep/samples/single30.RData")
set.seed(123)

# Source the SMC algorithm and diagnostics (no specific likelihood assumed)


Smat <- cbind(d = c(-1, 1))

# emc <- single[[1]]

new_model <- DDM()
new_model$bound$minmax[1,'sv'] <- 1e-2
new_DDM <- function(){return(new_model)}

# des <- design(factors=list(subjects=1,C=c("HARD", "EASY")),
#                       Rlevels = c(0, 1),
#                       formula =list(v~C,a~1, t0~1, s~1, Z~1, sv~1),
#                       constants=c(s=log(1)),
#                       model = new_DDM)
# p_vector <- sampled_pars(des, doMap = F)
# p_vector[] <- c(1.5, 0.1, log(1), log(.25), qnorm(.5), log(.3))

des <- design(factors=list(subjects=1,C=c("HARD", "EASY")),
              Rlevels = c(0, 1),
              formula =list(v~C,B~1, t0~C, A ~ 1, sv ~ lR),
              constants = c(sv = 1),
              model = LBA)
p_vector <- sampled_pars(des, doMap = F)
p_vector[] <- c(1.5, 0.1, log(1), log(.25), .1, log(.3), 1.2)



dat <- make_data(p_vector, des, n_trials = 50)
# plot_density(dat, factors = "S")

emc <- make_emc(dat, des, type = "single")
emc <- fit(emc)
emc <- save(emc, file = "emc.RData")
load("emc.RData")
# file5 <- file.path("../EMC2/paper/script/Comparison_Stan/ddm5.stan")
#
# mod5 <- cmdstan_model(file5)
# standata5 <- list(N=nrow(dat), cnd = ifelse(dat$C == "HARD", 1, 2),
#                   Ncnds = length(unique(dat$C)), rt = dat$rt,
#                   resp = as.numeric(as.character(dat$R)), parallel = 0)
# fit5 <- mod5$sample(data = standata5, seed = 123, chains = 4, parallel_chains = 4,
#                     iter_warmup = 300, iter_sampling = 1000, refresh = 0)
#
# draws_df <- fit5$draws(format = "df")  # returns a data frame of draws
# save(draws_df, file = "stan.RData")
# load("stan.RData")
# apply(draws_df, 2, quantile, probs = c(.025, .5, .975))[,"v_s_tr"]
# credint(emc, use_par = "sv")


# ---------- prior ----------
mu_ref <- emc[[1]]$prior$theta_mu_mean
Sigma_ref <- emc[[1]]$prior$theta_mu_var
data <- emc[[1]]$data[[1]]
ll <- emc[[1]]$model

# ---------- Run SMC ----------
cat("\nRunning SMC algorithm...\n")
set.seed(NULL)
source("SMC_super_fast.R")
source("smc_diagnostics.R")

system.time(
  smc_res <- enhanced_smc_elite(
    data = data, loglik_fn = ll, mu_ref = mu_ref, Sigma_ref = Sigma_ref,
    seed = sample.int(150, 1), verbose = TRUE, da_enable = TRUE,
    M = 5000, ll_cache_enable = TRUE, n_cores = 1)
)


source("smc_health.R")
res <- smc_quick_check_deterministic(smc_res, mu_ref, Sigma_ref)
# library(profvis)

# res <- profvis({
  # put the *entire* optimization here so the stack is captured
  # smc_res <- enhanced_smc_elite(
  #   data = data, loglik_fn = ll, mu_ref = mu_ref, Sigma_ref = Sigma_ref,
  #   seed = sample.int(75, 1), verbose = TRUE,
  #   resample_threshold = .6,
  #   M = 5000,max_rounds = 150
  # )
# })
# res

load("emc.RData")
theta_true_vec <- credint(emc, probs = .5)[[1]]

pars <- parameters(emc, selection = "alpha")
pars <- pars[,-1]

recovery_stats <- diagnose_recovery(
  true_values = pars,
  smc_res = smc_res,
  main_title = "Linear Regression: Parameter Recovery",
  save_plot = FALSE,
  density_method = "smooth"  # Use smooth density estimation
)

credint(emc)

# for(q in 2:ncol(pars)){
#   lines(density(pars[,q]))
# }


