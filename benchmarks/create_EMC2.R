# We're going to analyze a subset of the English Lexicon Project (ELP) data, where people completed
# a lexical decision-making task
rm(list = ls())
library(EMC2)
load("~/Documents/2026/Erasmus/data/ELP.RData")

# Which is 0 for non-words

Smat <- cbind(d = c(-1, 1))
des_ELP <- design(list(v ~ S + LogFreq,
                       a ~ 1, t0 ~ 1, Z ~ 1, sv ~ 1, s ~ S),
                  data = data, model = DDM,
                  contrasts = list(v = list(S = Smat)),
                  constants = c(s = log(1)))


ELP_DDM <- make_emc(data, design = des_ELP, type = "diagonal-gamma")

ELP_DDM <- fit(ELP_DDM, cores_per_chain = 3)
save(ELP_DDM, file = "benchmarks/samples/full_EMC2")

summary(get_prior(ELP_DDM))
