# Adjusted Nested SMC Implementation Plan

## Goal

Keep the same statistical target as the original nested SMC note, but change the implementation design so it is computationally natural for this project.

The central adjustment is:

- the expensive object is the local likelihood `p(y_i | alpha)`
- so local likelihood evaluations should be treated as a cached resource
- changes in population state `phi` and bridge level `lambda` should mostly trigger cheap reweighting, not fresh local likelihood work

This document uses:

- `population` for the hyperparameter layer
- `local` for the per-group latent layer
- the intercept-only formulation for now

So for now:

\[
\phi = (\mu, \Sigma)
\]

with

\[
\alpha_i \mid \phi \sim \mathcal N(\mu,\Sigma),
\qquad
y_i \sim p(y_i \mid \alpha_i),
\qquad
i = 1,\dots,p.
\]

This is the same nested idea as before, just with `X_i = 1`.

---

## Ground Rules

### Keep the framework clean

- maintain one clear path
- do not accumulate modes
- do not add bells and whistles unless they earn their keep
- if something is tried and it does not help, remove it afterwards
- scaffolding does not become permanent just because it exists

### Be brutally honest about progress

Something counts as complete only if it is:

- implemented in the maintained path
- wired into the actual algorithm
- sourceable and runnable
- exercised by a benchmark or test that we actually use

Something does not count as complete if it is only:

- a helper with no maintained caller
- a scratch script
- a compatibility shim
- a partial approximation that has not been benchmarked
- a design note without code

### Optimize the metric that matters

The main efficiency metric is:

- total number of local likelihood evaluations of `log p(y_i | alpha)`

This is the primary score.

Other metrics are secondary:

- wall-clock time
- number of population bridge rounds
- number of population resamples
- number of population mutation proposals

If a method uses fewer bridge rounds but more local likelihood evaluations, that does not count as an improvement.

### Keep the repo tidy

Core code stays in the main folder.

Non-core scripts go in dedicated folders:

- `benchmarks/`
- `dev/tests/`
- `dev/scratch/`

---

## 1. Statistical Target

The full posterior is

\[
\pi_1(\phi,\alpha_{1:p})
\propto
p(\phi)\prod_{i=1}^p
p(y_i\mid \alpha_i)\,
\mathcal N(\alpha_i;\mu,\Sigma).
\]

Use the continuous bridge

\[
\pi_\lambda(\phi,\alpha_{1:p})
\propto
p(\phi)\prod_{i=1}^p
p(y_i\mid \alpha_i)\,
\mathcal N(\alpha_i;\mu,\Sigma)^\lambda\,
g_i(\alpha_i)^{1-\lambda},
\qquad \lambda \in [0,1].
\]

where `g_i` is a working local density.

Special cases:

- `lambda = 0`: local targets are `p(y_i | alpha_i) g_i(alpha_i)`
- `lambda = 1`: full hierarchical posterior

This target family does not change.

What changes is the implementation strategy for approximating it.

---

## 2. Adjustment to the Architecture

### 2.1 What the Phase 1 scaffold does now

The current Phase 1 code uses:

- population particles over `phi`
- a separate mutable `local_module` for every `(population particle, group)` pair
- optional local rejuvenation inside the nested run

This is a valid scaffold because it is easy to reason about and it proves the bridge machinery can run.

But it is probably not the computationally natural final design.

### 2.2 Why the naive nested implementation is too expensive

In the current naive design, local work is duplicated across population particles.

That means if there are:

- `N` population particles
- `p` groups
- `M` local particles

then local rejuvenation can scale like:

\[
O(N p M)
\]

in fresh local likelihood evaluations per update if done naively.

That is exactly what produced the large benchmark count in the current Phase 1 implementation.

The issue is not the target.

The issue is duplicating local candidate generation and local likelihood reevaluation across many population particles when the expensive local likelihood itself has not changed.

### 2.3 The central design correction

The correct computational principle is:

> never reevaluate a local likelihood unless a genuinely new local candidate `alpha` has been created

When `phi` or `lambda` changes, the term that changes is the Gaussian bridge factor, not the local likelihood term.

So local likelihood values should be cached and reused.

---

## 3. Adjusted Inner Representation

The maintained inner object should no longer be a separate local particle cloud for every population particle.

Instead, for each group `i`, maintain one shared adaptive local bank:

\[
\mathcal B_i =
\left\{
\alpha_{i,m},
\ell_{i,m},
\log g_i(\alpha_{i,m})
\right\}_{m=1}^{M_i}
\]

where

\[
\ell_{i,m} = \log p(y_i \mid \alpha_{i,m})
\]

is evaluated once and then cached.

Optional additional cached quantities may include:

- proposal density values
- bank provenance or refresh batch id
- transport coordinates if useful

The important point is:

- the bank is shared across all population particles
- local likelihoods are attached to bank elements, not recomputed per population particle

---

## 4. Adjusted Population Representation

Each population particle still carries:

\[
\phi^{(n)} = (\mu^{(n)}, \Sigma^{(n)}).
\]

But it does not carry its own separate full local system for every group.

Instead, the population layer queries the shared local banks to compute group contributions.

So the maintained state becomes:

- population particles `phi^(1:N)`
- population weights
- current bridge value `lambda`
- one shared local bank `B_i` per group

This is still nested in spirit, but it is bank-based rather than per-particle duplicated.

---

## 5. Local Conditional Target and Shared-Bank Approximation

For fixed population particle `phi^(n)` and bridge level `lambda`, the local target remains

\[
\gamma_{i,\lambda}(\alpha \mid \phi^{(n)})
\propto
p(y_i \mid \alpha)\,
\mathcal N(\alpha;\mu^{(n)},\Sigma^{(n)})^\lambda\,
g_i(\alpha)^{1-\lambda}.
\]

The normalizing constant remains

\[
Z_{i,\lambda}(\phi^{(n)}) =
\int
p(y_i \mid \alpha)\,
\mathcal N(\alpha;\mu^{(n)},\Sigma^{(n)})^\lambda\,
g_i(\alpha)^{1-\lambda}\,
d\alpha.
\]

The difference is how this is approximated.

Using a shared bank with proposal density `h_i`, a natural estimator is

\[
\hat Z_{i,\lambda}(\phi^{(n)})
=
\frac{1}{M_i}
\sum_{m=1}^{M_i}
\frac{
p(y_i \mid \alpha_{i,m})\,
\mathcal N(\alpha_{i,m};\mu^{(n)},\Sigma^{(n)})^\lambda\,
g_i(\alpha_{i,m})^{1-\lambda}
}{
h_i(\alpha_{i,m})
}.
\]

The important practical point is:

- `p(y_i | alpha_{i,m})` is cached
- only the Gaussian bridge term changes with `phi^(n)` and `lambda`

This turns the inner layer into shared adaptive importance banks rather than per-particle local particle systems.

---

## 6. Adjusted Design Principle for Updates

The right update order is:

1. reweight existing bank points under the new `phi` and `lambda`
2. measure bank adequacy
3. only enrich the bank if adequacy is poor

This replaces the current Phase 1 rule of always doing a fixed number of local rejuvenation moves.

In other words:

- reweight first
- refresh only if diagnostics say reweighting is failing

That is the natural SMC logic for this problem.

---

## 7. Adaptive Bridge in the Adjusted Design

Let

\[
0=\lambda_0 < \lambda_1 < \cdots < \lambda_T = 1.
\]

As before, do not fix `T` in advance.

Choose the next bridge value adaptively using population ESS or rCESS.

For a candidate `lambda > lambda_t`, compute hypothetical population increments using the current shared local banks.

The bridge schedule should still be driven by:

- population ESS
- or population rCESS

The adjustment is only that the group contributions are computed from shared banks, not from separate local modules.

---

## 8. Bank Reweighting Under a Population Particle

For each group `i`, bank point `m`, population particle `n`, and bridge level `lambda`, define the cheap log weight contribution

\[
\tilde w_{i,m}^{(n)}(\lambda)
=
\ell_{i,m}
+
\lambda \log \mathcal N(\alpha_{i,m};\mu^{(n)},\Sigma^{(n)})
+
(1-\lambda)\log g_i(\alpha_{i,m})
-
\log h_i(\alpha_{i,m}).
\]

Everything here is cheap except the cached local log likelihood `\ell_{i,m}`, which has already been paid for.

These bank weights are used to:

- estimate `Z_{i,\lambda}(\phi^(n))`
- estimate incremental ratios
- compute local ESS diagnostics

No fresh local likelihood evaluation is needed for this reweighting step.

---

## 9. Adequacy Diagnostics for Shared Banks

Each group bank needs its own adequacy diagnostic.

The basic one is local ESS after reweighting.

For group `i`, evaluate bank adequacy across the current population particles.

Practical options:

- minimum local ESS across population particles
- a low quantile of local ESS across population particles
- coefficient of variation of importance weights
- instability in estimated `log Z_i`

Recommended default:

- use a low quantile of local ESS across population particles, not the minimum

This avoids overreacting to one pathological population particle.

The bank is enriched only if this adequacy measure falls below threshold.

---

## 10. Bank Enrichment

When group `i` needs refresh, enrich `B_i` with a small number of new local candidates.

The rule should be adaptive:

- first try `0` new points
- if adequacy is poor, add a small batch such as `4` or `8`
- only escalate further if still poor

Do not default to large fixed refreshes such as `64` new local points every time.

Each new candidate `alpha` requires:

- one fresh evaluation of `log p(y_i | alpha)`
- one cached value of `log g_i(alpha)`
- bookkeeping for proposal density `h_i(alpha)` if needed

Old candidates remain in the bank unless there is a clear reason to prune them.

This makes bank growth incremental and measurable.

---

## 11. Candidate Sources for Bank Enrichment

New local candidates can come from a mixture of sources.

Recommended simple mixture:

- draws from the working density `g_i`
- draws from the saved local fit bank from the current local SMC engine
- draws from a refresh proposal centered around the current population region

This should stay simple at first.

Do not add many proposal families unless they clearly help.

The aim is for banks to remain useful over a wider range of `phi` values, so that frequent enrichment is unnecessary.

---

## 12. Optional Later Compression: Population Clustering

If one shared bank per group is still not adequate enough, the next approximation is not to return to one local system per population particle.

The next step is:

- cluster population particles in `phi` space
- maintain one bank refresh context per cluster
- let particles within a cluster share those refreshed bank updates

This is a later optimization, not Phase 1 or Phase 2.

It should only be added if the simpler shared-bank design is not good enough.

---

## 13. Adjusted Population Mutation

Population mutation still targets `phi`.

For a proposed `phi'`, acceptance uses the bank-based estimates of the group contributions at the current `lambda`.

The important correction is:

- a `phi` proposal should not automatically trigger fresh local likelihood evaluations
- it should first reuse the shared banks
- only if bank adequacy is poor under the proposed region should enrichment happen

So population mutation becomes:

1. propose `phi'`
2. evaluate approximate group contributions using shared banks
3. accept or reject
4. if accepted and some group banks are inadequate in the newly occupied region, enrich those banks

This preserves the idea that local likelihood calls are scarce and precious.

---

## 14. What the Current Code Already Gives Us

The current codebase already has useful pieces.

### 14.1 What is genuinely useful

- [SMC_super_fast.R](/Users/nstevenson/Documents/2025/FastHierarchical/SMC_super_fast.R) gives the maintained local SMC engine
- [smc_core.R](/Users/nstevenson/Documents/2025/FastHierarchical/smc_core.R) gives shared SMC numerics such as ESS and stepping utilities
- [new_SMC_cache.R](/Users/nstevenson/Documents/2025/FastHierarchical/new_SMC_cache.R) already wraps local fits and preserves working-prior information and local particle banks
- [benchmarks/local_shifted_gamma_vs_stan.R](/Users/nstevenson/Documents/2025/FastHierarchical/benchmarks/local_shifted_gamma_vs_stan.R) is a clean local benchmark
- [benchmarks/nested_normal_normal.R](/Users/nstevenson/Documents/2025/FastHierarchical/benchmarks/nested_normal_normal.R) is a dedicated nested benchmark

### 14.2 What is only a scaffold

- [nested_population_SMC.R](/Users/nstevenson/Documents/2025/FastHierarchical/nested_population_SMC.R) currently implements the naive per-particle local-module design
- that file is useful as a proof of concept for bridge stepping and bookkeeping
- but it should not be mistaken for the final adjusted architecture

### 14.3 Honest status

What counts today:

- there is a runnable population bridge scaffold
- there is a real local SMC engine
- there is benchmark infrastructure
- local likelihood evaluations are now counted inside the current nested run

What does not count today:

- end-to-end optimal local likelihood accounting
- shared-bank inner machinery
- adaptive bank enrichment
- bank adequacy diagnostics
- a final computationally natural nested algorithm

---

## 15. What Must Change in Code Terms

The main code shift is:

- from `local_module` attached to each `(population particle, group)`
- to `local_bank` attached only to each group

### 15.1 New maintained object

Introduce a first-class `local_bank` object with methods such as:

- `predict_group_increment(phi, lambda_new)`
- `estimate_logZ(phi, lambda)`
- `ess_under(phi, lambda)`
- `enrich(phi_context, lambda, n_new)`
- `clone()` only if needed for temporary work, not as core design

The key invariant:

- cached local likelihoods stay with bank points permanently

### 15.2 Population driver

The maintained population driver should:

- hold `phi`, weights, and `lambda`
- query group banks for increments
- choose `lambda` adaptively
- resample population particles
- mutate population particles
- trigger bank enrichment only when adequacy diagnostics say so

### 15.3 Benchmark accounting

Benchmark output should report:

- local prefit likelihood count from the initial local SMC runs
- nested bank initialization count
- bank enrichment count
- any population-proposal-triggered enrichment count
- grand total

The grand total is the number that matters.

---

## 16. Recommended Implementation Order

### Phase A: honest accounting

1. keep the current scaffold runnable
2. add full end-to-end local likelihood accounting, including the up-front local SMC fits
3. report this in the nested benchmark

Completion criteria:

- benchmark prints nested-only and end-to-end local likelihood counts
- counts are traceable to code paths

### Phase B: shared-bank baseline

1. create `local_bank`
2. initialize one bank per group from the current local SMC outputs
3. replace per-particle local modules with bank-based reweighting
4. disable unconditional per-particle local rejuvenation

Completion criteria:

- dedicated nested benchmark reaches `lambda = 1`
- shared-bank code is the maintained path
- total local likelihood count is substantially below the naive per-particle design

What does not count:

- leaving the old per-particle design as the maintained default
- only adding a bank helper while still running the naive algorithm

### Phase C: adaptive enrichment

1. add adequacy diagnostics for each group bank
2. enrich banks only when diagnostics fail
3. top up banks in small adaptive batches

Completion criteria:

- enrichment is diagnostic-driven, not fixed-schedule driven
- benchmark reports how many fresh local likelihood calls came from enrichment
- benchmark shows whether enrichment actually improves posterior quality enough to justify its cost

### Phase D: only if needed

Possible later steps:

- population clustering for bank refresh contexts
- more sophisticated bank proposals
- pruning or compression of stale banks

These should only be added if the simpler shared-bank design is inadequate.

---

## 17. Honest Bottom Line

The original nested SMC target is still the right target.

The naive implementation with one local system per population particle is useful as a scaffold, but it is probably not the right maintained algorithm for this problem.

The adjusted plan is:

- keep population SMC over `phi`
- keep the continuous bridge in `lambda`
- replace duplicated local systems with one shared adaptive bank per group
- cache local likelihoods permanently
- reweight first
- enrich only when diagnostics say the bank is failing
- measure success primarily by total local likelihood evaluations

This is less pure than the fully literal nested construction.

But it is much more likely to be the computationally natural algorithm for the actual problem we care about.
