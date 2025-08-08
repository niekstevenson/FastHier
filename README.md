# FastHier
Fast hierarchical modeling

The adaptive SMC sampler now skips elite mixture fitting when the difficulty
estimate (`difficulty_ema`) falls below a threshold. This avoids unnecessary
mixture updates on easy stages while re-enabling them once the difficulty
rises again. The threshold defaults to `0.2` and can be tuned via the
`mix_difficulty_threshold` argument in `enhanced_smc_elite()` or `auto_smc()`.


The sequential Monte Carlo sampler automatically increases the multiple-try
Metropolis factor (`mtm_K`) as the inferred difficulty of the problem grows,
starting at 1 and interpolating toward a higher cap for challenging targets.

