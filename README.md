# FastHier
Fast hierarchical modeling

The sequential Monte Carlo sampler automatically increases the multiple-try
Metropolis factor (`mtm_K`) as the inferred difficulty of the problem grows,
starting at 1 and interpolating toward a higher cap for challenging targets.
