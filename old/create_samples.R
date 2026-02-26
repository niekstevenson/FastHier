rm(list = ls())
library(EMC2)
set.seed(123)


N_subj <- 30
N_trials <- 100
# des <- design(formula = list(m ~ 1, t0 ~ 1, s ~ 1),
#               factors = list(subjects = 1:N_subj),
#               Rlevels = 1:2,
#               model = LNR)

des <- design(factors=list(subjects=1:N_subj,C=c("HARD", "EASY")),
              Rlevels = c(0, 1),
              formula =list(v~C,a~1, t0~1, s~1, Z~1, sv~1),
              constants=c(s=log(1)),
              model = DDM)
p_vector <- sampled_pars(des, doMap = F)
p_vector[] <- c(1.5, 0.1, log(1), log(.25), qnorm(.5), log(.3))

mu <- p_vector
cov <- diag(c(.25, .1, .1, .05, .05, .05)^2)

REs <- make_random_effects(des, mu, N_subj, covariances = cov)
dat <- make_data(REs, des, n_trials = N_trials)
# I mean sure
plot_density(dat, factors = c("C"))
# pri_single <- prior(des, type = "single", pmean = mu, theta_mu_var = 5*cov)
# pri <- prior(des, type = "standard", mu_mean=mu, mu_var = cov)

run_single <- function(dat, design){
  emc <- make_emc(dat, design, type = "single")
  emc <- fit(emc, cores_for_chains = 1)
  return(emc)
}
library(parallel)

dat_single <- split(dat, dat$subjects)
single <- mclapply(dat_single, run_single, des, mc.cores = 12)

full <- make_emc(dat, des, type = "diagonal-gamma")
full <- fit(full, cores_per_chain = 4)

save(full, file = "full_DDM30.RData")
save(single, file = "single_DDM30.RData")




