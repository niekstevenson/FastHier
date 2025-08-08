# FastHier
Fast hierarchical modeling

## Elite mixture fitting

The adaptive SMC sampler now skips elite mixture fitting when the difficulty
estimate (`difficulty_ema`) falls below a threshold. This avoids unnecessary
mixture updates on easy stages while re-enabling them once the difficulty
rises again. The threshold defaults to `0.2` and can be tuned via the
`mix_difficulty_threshold` argument in `enhanced_smc_elite()` or `auto_smc()`.
