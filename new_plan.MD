# New Plan: Extended-State Two-Stage Hierarchical SMC

## Purpose

Build a fast two-step hierarchical inference pipeline where:

- expensive subject-level likelihood work is done once in independent local analyses;
- the group-level stage only uses cheap Gaussian hierarchy corrections;
- memory stays bounded by storing bank indices rather than full replicated local states;
- the method scales to many subjects;
- Stan remains the calibration and debugging reference.

This plan is best viewed as a new mode for the repository: a fast approximate-first hierarchical engine with a clear path to later strictification if needed.

## Model

For groups `i = 1, ..., p`:

- group parameter: `alpha_i in R^d`
- hyperparameters: `B in R^(q x d)`, `Sigma in S_{++}^d`
- group design row/object: `X_i`
- group data: `y_i`
- subject likelihood: `l_i(alpha_i) = p(y_i | alpha_i)`

Hierarchical model:

`alpha_i | B, Sigma ~ N(X_i B, Sigma)`

Posterior:

`pi(B, Sigma, alpha_1:p) propto p(B, Sigma) prod_i l_i(alpha_i) N(alpha_i; X_i B, Sigma)`

## Core Idea

For each subject, choose a working prior `g_i(alpha_i)` and define the local target

`q_i(alpha_i) propto l_i(alpha_i) g_i(alpha_i)`

with normalizing constant

`Z_i = int l_i(alpha_i) g_i(alpha_i) d alpha_i`.

Define the correction ratio

`r_i(alpha_i; B, Sigma) = N(alpha_i; X_i B, Sigma) / g_i(alpha_i)`.

Then

`l_i(alpha_i) N(alpha_i; X_i B, Sigma) = Z_i q_i(alpha_i) r_i(alpha_i; B, Sigma)`,

so the full posterior can be rewritten, up to the constant `prod_i Z_i`, as

`pi(B, Sigma, alpha_1:p) propto p(B, Sigma) prod_i q_i(alpha_i) r_i(alpha_i; B, Sigma)`.

After Stage 1, the expensive `l_i(alpha_i)` evaluations are not used in Stage 2.

## Stage 1: Independent Local Analyses

For each subject `i`:

1. Choose a working prior `g_i(alpha_i)`.
2. Run a local SMC analysis targeting `q_i(alpha_i) propto l_i(alpha_i) g_i(alpha_i)`.
3. Store a local bank representing `q_i`.

Preferred first representation:

- resampled approximately unweighted draws `a_{i,1}, ..., a_{i,M_i} ~ q_i`

Minimum stored artifact per subject:

- local bank `a_{i,1:M_i}`
- seed/provenance for reproducibility
- local diagnostics
- optional local mean/covariance summaries for screening

For the first prototype, Stage 1 should be allowed to spend extra work on hard subjects. This is the right place to pay for weak identification or boundary-sensitive groups.

## Stage 2: Outer Bridge SMC on the Extended State

Define the tempered sequence

`pi_t(B, Sigma, alpha_1:p) propto p(B, Sigma) prod_i q_i(alpha_i) r_i(alpha_i; B, Sigma)^(beta_t)`

with

`0 = beta_0 < beta_1 < ... < beta_T = 1`.

Interpretation:

- at `beta = 0`, subjects are independent under `q_i`
- at `beta = 1`, the hierarchy is fully enforced

## Particle State

Conceptual particle:

`X^(n) = (B^(n), Sigma^(n), alpha_1^(n), ..., alpha_p^(n))`

Implementation should use bank indices:

- store `k_i^(n) in {1, ..., M_i}`
- define `alpha_i^(n) = a_{i, k_i^(n)}`

Actual particle state:

`X^(n) = (B^(n), Sigma^(n), k_1^(n), ..., k_p^(n))`

This is the default implementation because it bounds memory cleanly.

## Initialization

For each outer particle:

1. Draw `B, Sigma` from the hyperprior or a pilot-informed initialization.
2. For each subject `i`, draw `k_i` uniformly from the local bank if the bank is approximately unweighted.

At `beta = 0`, all particle weights are equal.

## Reweighting

When moving from `beta_{t-1}` to `beta_t`, update

`log w_t^(n) = log w_(t-1)^(n) + Delta beta_t sum_i log r_i(alpha_i^(n); B^(n), Sigma^(n))`

where `Delta beta_t = beta_t - beta_{t-1}`.

Choose temperatures adaptively using a CESS or ESS target.

## Resampling

If particle ESS falls below threshold, resample.

Default:

- systematic or stratified resampling

After resampling:

- reset weights to `1 / N`

## Mutation Kernel

At each temperature, apply a kernel invariant to `pi_t`.

Recommended blocked kernel:

1. Update `B | Sigma, alpha_1:p`
2. Update `Sigma | B, alpha_1:p`
3. Refresh a subset of subject latents

### Update `B`

If `p(B)` is Gaussian and the hierarchy is parameterized conjugately, use Gibbs.

Otherwise:

- use MH, MALA, or HMC on `B`

### Update `Sigma`

If `p(Sigma)` is chosen conjugately, use Gibbs.

Otherwise:

- update a Cholesky parameterization by MH or HMC

### Refresh Subject Latents

For a chosen subset `J_t` of subjects:

1. Propose `alpha_i' ~ q_i` by drawing a new bank index `k_i'`
2. Accept with probability

`min(1, [r_i(alpha_i'; B, Sigma) / r_i(alpha_i; B, Sigma)]^(beta_t))`

This is the key rejuvenation move. It uses the stored local law directly and does not call the expensive likelihood.

## Sufficient-Statistics Acceleration

For each particle maintain:

- `S_xx = sum_i X_i^T X_i`
- `S_xa = sum_i X_i^T alpha_i^T`
- `S_aa = sum_i alpha_i alpha_i^T`

Then the quadratic form

`sum_i (alpha_i - X_i B)^T Sigma^(-1) (alpha_i - X_i B)`

can be evaluated without looping over all subjects every time.

When a subject latent is refreshed:

- update `S_xa` incrementally
- update `S_aa` incrementally

This is essential for scale.

## Recommended Build Order

1. Implement the approximate version with fixed finite local banks.
2. Validate against Stan on the current gamma example.
3. Add hard-subject refinement in Stage 1.
4. Add multi-anchor local analysis if shift-tail coverage is weak.
5. Decide whether this mode remains approximate-only or is wrapped in a strict final layer.

## Intended Role in This Repository

This plan should be treated as a new algorithmic mode, not as a small patch to the current `H1` path.

Current `H1`:

- strict adaptive pseudo-marginal over hyperparameters
- can rediscover missed subject-tail regions because it re-evaluates the true likelihood in the final stage

This new plan:

- avoids inner likelihood reevaluation in Stage 2
- is likely faster and stabler than the current fixed-cache path
- is attractive as an approximate-first engine and adaptation layer

## Caveats

### Caveat A: Not Just Tuning

The fixed-bank version is not strict for the original continuous model.

If Stage 2 only uses a finite stored bank per subject, then the final target is the posterior induced by the empirical local laws, not the exact hierarchical posterior. That is not merely a tuning detail; it is a definition-level limitation of the fixed-bank construction.

Implication:

- if strict correctness for the published `phi` posterior is mandatory, this must be handled explicitly

Three ways to handle that:

1. Use this method as an approximate-first engine, then run a strict final method afterward.
2. Add a mechanism to regenerate exact or asymptotically exact draws from `q_i` on demand in Stage 2.
3. Accept and label this mode as approximate.

### Caveat B: Important, But Can Be Added Later

These are real improvements, but they do not need to block a first prototype:

1. Multi-anchor or subject-specific `g_i`

   A single shared working prior is likely too brittle for weakly identified shift subjects. A pilot-based subject-specific or anchor-specific `g_i` should improve bank coverage substantially.

2. Adaptive subject refresh fraction

   Refreshing only `1%` to `5%` of subjects per sweep may be too conservative once `Sigma` starts moving. The refresh rate should eventually depend on temperature, recent movement in `Sigma`, or acceptance diagnostics.

3. Selective refinement of hard subjects

   Hard subjects should get larger local banks, extra anchors, or local reruns. Easy subjects should remain cheap.

4. Bank coverage diagnostics

   The method needs subject-level diagnostics that detect when the local bank fails to cover the hierarchical region actually used in Stage 2.

5. Continuous local surrogates

   A fitted local surrogate can reduce disk and memory costs, but it should be introduced only after the bank-based version works.

## Exact Variant

This section describes how to regain exactness while keeping the same overall two-stage architecture.

### Principle

To regain exactness, the Stage 1 artifacts must be treated as proposal mechanisms, not as the final local law itself.

That means:

- Stage 1 still does expensive local work once and produces bank-informed proposal objects.
- Stage 2 still carries the extended state `(B, Sigma, alpha_1:p)`.
- But every exact mutation of `alpha_i` is corrected using the true local likelihood `l_i(alpha_i)`.

With that change, the final Stage 2 target is again the exact continuous hierarchical posterior, provided initialization and mutation are both corrected exactly.

### Exact Initialization

At `beta = 0`, the target is

`pi_0(B, Sigma, alpha_1:p) propto p(B, Sigma) prod_i q_i(alpha_i)`

with

`q_i(alpha_i) propto l_i(alpha_i) g_i(alpha_i)`.

There are two exact ways to initialize:

1. Draw exactly from each `q_i`

   This is exact in principle, but usually not practical unless the local sampler can be resumed on demand.

2. Draw from a continuous proposal `s_i(alpha_i)` and importance-correct

   Sample `alpha_i ~ s_i` and assign initial weight proportional to

   `prod_i [l_i(alpha_i) g_i(alpha_i) / s_i(alpha_i)]`.

   The unknown local constants `Z_i` cancel because they do not depend on the particle state.

In practice, the second route is the realistic exact version.

### Exact Mutation of Subject Latents

For any refreshed subject `i`, use a bank-informed continuous proposal `s_i(alpha_i' | current state)` and accept using the exact MH ratio under the tempered target

`pi_t(alpha_i | B, Sigma) propto l_i(alpha_i) g_i(alpha_i)^(1 - beta_t) N(alpha_i; X_i B, Sigma)^(beta_t)`.

If the proposal is independent of the current `alpha_i`, the log acceptance ratio is

`log a = min(0, log l_i(alpha_i') - log l_i(alpha_i) + (1 - beta_t)[log g_i(alpha_i') - log g_i(alpha_i)] + beta_t[log N(alpha_i'; X_i B, Sigma) - log N(alpha_i; X_i B, Sigma)] + log s_i(alpha_i) - log s_i(alpha_i'))`.

This is the key tradeoff:

- exactness is recovered;
- but exact local likelihood evaluations return during Stage 2, at least for refreshed subjects.

### Delayed-Acceptance Option

To control cost, the exact latent-refresh move should use delayed acceptance:

1. First screen the proposal using only cheap quantities, for example the hierarchy correction term and proposal density ratio.
2. Only if it passes, evaluate the exact local likelihood and do the final MH correction.

This preserves exactness while reducing wasted expensive likelihood calls on obviously bad proposals.

### How This Ties Into the Current Repository

The current code already contains most of the proposal machinery needed for the exact variant.

Stage 1 reuse:

- `run_outer_subjects()` already performs the independent local analyses.
- `build_subject_surrogate_from_outer()` already converts a subject outer fit into a continuous proposal object.
- `make_theta_proposal_adaptive()` already builds a mixture proposal with local and anchor-informed components and exposes `draw()` and `log_q()`.

Exact local likelihood reuse:

- `loglik_shifted_gamma()` is already the exact local likelihood for the current gamma example.

Hyperprior mapping reuse:

- `gaussian_map_fn` already maps the current hyperparameter vector to a Gaussian prior over subject parameters.
- In the current gamma model, `phi` is effectively the intercept-only special case of `(B, Sigma)`.
- For a more general design matrix, the same interface can be extended so that the prior mean is `X_i B` rather than a shared intercept.

Where the new code would live:

- keep the existing `run_outer_subjects()` Stage 1 path;
- keep `build_subject_surrogate_from_outer()` as the Stage 1 proposal-builder;
- add a new inner mode, for example `extended_exact`, rather than modifying `H1` in place;
- implement a new outer kernel that stores `(B, Sigma, k_1:p)` or direct `alpha_1:p`, plus cached sufficient statistics.

Recommended reuse pattern:

1. Use the existing local outer fits to create per-subject proposal objects.
2. Use those proposal objects to initialize the extended-state SMC exactly by weighting with `l_i g_i / s_i`.
3. Reuse the same proposal objects during Stage 2 latent refreshes.
4. Use the existing exact likelihood function only inside initialization and accepted-screen exact refresh proposals.

### Expected Cost

Relative to the pure fixed-bank approximate version, the exact variant adds back likelihood cost in two places:

1. Initialization

   If initialized from continuous proposals `s_i`, this costs about `N * p` local likelihood evaluations.

2. Latent refreshes during Stage 2

   If `rho_t` is the fraction of subjects refreshed at round `t` and `m_t` is the number of proposals per refreshed subject, the added cost is roughly

   `sum_t N * rho_t * p * m_t`

   exact local likelihood evaluations.

This is the dominant extra cost of regaining exactness.

### Expected Cost Relative to Existing Modes

Qualitatively:

- more expensive than the pure fixed-bank extended-state approximation;
- usually much cheaper than full strict pseudo-marginal if refresh rates remain sparse;
- likely competitive if proposals are good and delayed acceptance is used well.

For the current gamma example, a reasonable expectation is:

- the exact variant will cost a few times more than the pure bank version;
- but it may still be cheaper and much stabler than the current noisy strict PM path, because it spends exact likelihood work only on refreshed subject states rather than on Monte Carlo evidence estimation for every subject at every hyperparameter proposal.

That last point is a hypothesis, not a benchmark. The actual crossover depends on:

- number of subjects;
- refresh fraction;
- number of temperatures;
- local likelihood cost;
- proposal quality from the Stage 1 surrogates.

## Repository Execution Plan

This section is the concrete implementation plan for this repository as it exists today.

The key constraint is that the current codebase does not yet represent the general hierarchy as explicit `(B, Sigma, X_i)` objects inside the inner engine. The current maintained path uses a packed hyperparameter vector `phi` and a Gaussian map `gaussian_map_fn(phi, d)`.

So the exact execution plan should be staged:

1. First implement an exact extended-state engine for the current intercept-only diagonal-Gaussian hierarchy.
2. Then, once that works, generalize the same engine to explicit `(B, Sigma, X_i)`.

### Phase 1 Goal

Add a new exact inner mode for the current gamma example with state

`(phi, alpha_1:p)`

where:

- `phi` is the current packed hyperparameter vector already used in the repository;
- `alpha_i` is the subject latent carried explicitly in the particle state.

This is already enough to recover exactness for the current model while reusing the present code.

### Existing Code To Reuse As-Is

From `run_two_step_smc.R`:

- `run_outer_subjects()`
- `estimate_phi_anchor_from_outer()`
- `loglik_shifted_gamma()`
- current experiment wiring and scorecard flow

From `new_SMC_cache.R`:

- `build_subject_surrogate_from_outer()`
- `make_theta_proposal_adaptive()`
- Gaussian prior helpers such as `.log_prior_gauss_mat()`

From `make_prior.R`:

- `make_prior_phi_diag()`
- `phi_to_gaussian_params_diag_factory()`

From `smc_core.R` and `outer_SMC.R`:

- resampling helpers
- temperature adaptation helpers
- generic random-walk covariance adaptation patterns

### Existing Code That Should Be Extended, Not Replaced

1. `run_two_step_smc.R`

   Add a new execution mode, for example `X1` or `extended_exact`, without changing `E4` or `H1`.

2. `run_inner()`

   Either:

   - dispatch to a new `run_inner_extended_exact()`, or
   - add a third `inner_ll_mode`

   The cleaner path is a new wrapper so the current PM code stays untouched.

3. `build_subject_surrogate_from_outer()`

   Reuse it directly, but the exact mode also needs the working prior `g_i` stored explicitly alongside the surrogate.

### New Objects To Add

#### 1. Exact local proposal object

File:

- `new_SMC_cache.R`

New function:

- `build_subject_exact_local_from_outer(smc_out, data, subj_id, loglik_fn, working_prior, base_seed = NULL)`

Purpose:

- wrap the existing surrogate built from `build_subject_surrogate_from_outer()`;
- store the subject working prior `g_i`;
- store exact local likelihood handles needed during Stage 2;
- provide a single object used by initialization and latent refreshes.

Suggested contents:

- `surrogate`
- `working_prior`
- `data`
- `loglik_fn`
- `subj_id`
- `base_seed`

#### 2. Exact working-prior evaluator

File:

- `new_SMC_cache.R` or `smc_core.R`

New function:

- `log_working_prior_gaussian_mat(Theta, working_prior)`

Purpose:

- evaluate `log g_i(alpha)` exactly for the stored local working prior.

This is required because exact initialization and exact latent refresh use `l_i(alpha)` and `g_i(alpha)` explicitly.

#### 3. Initialization proposal builder

File:

- `new_SMC_cache.R`

New function:

- `make_theta_proposal_exact_init(local_obj, init_ctl = list())`

Purpose:

- build a continuous proposal close to `q_i`;
- reuse `make_theta_proposal_adaptive()` internally, but typically with local-heavy weights and little or no anchor dependence.

This is for exact initialization at `beta = 0`.

#### 4. Exact refresh proposal builder

File:

- `new_SMC_cache.R`

New function:

- `make_theta_proposal_exact_refresh(phi, local_obj, refresh_ctl = list(), callbacks = list())`

Purpose:

- build the continuous proposal used in exact latent refresh moves;
- reuse `make_theta_proposal_adaptive()` internally, but now allow current-hyperparameter anchor information.

This is for exact Stage 2 mutation.

### New Inner Engine To Add

File:

- `outer_SMC.R`

Main new function:

- `outer_smc_phi_extended_exact(local_objs, rprior_phi, logprior_phi, gaussian_map_fn, data_list, loglik_fn, M = ..., ...)`

This should be a separate engine, not a branch inside the current pseudo-marginal code path.

Its particle state for the current repository should be:

- `phi`
- subject latent matrix or index-backed latent representation for `alpha_1:p`
- cached subject log-likelihoods `log l_i(alpha_i)`
- cached sufficient statistics for the Gaussian hierarchy

### New Supporting Functions To Add

#### 1. Exact particle initialization

File:

- `outer_SMC.R`

New function:

- `init_extended_exact_particles(local_objs, phi_particles, phi_anchor, init_ctl = list())`

What it does:

1. sample `phi` from the current hyperprior;
2. for each subject and particle, propose `alpha_i` from the subject initialization proposal;
3. compute exact initial log weights using

   `log l_i(alpha_i) + log g_i(alpha_i) - log s_i(alpha_i)`

4. normalize the particle weights;
5. cache subject-level exact log-likelihood values.

Note:

- this is where exactness is recovered at `beta = 0` if exact draws from `q_i` are not available.

#### 2. Exact tempered log target contribution

File:

- `outer_SMC.R`

New function:

- `log_tempered_hier_contrib(phi_row, alpha_mat, beta, gaussian_map_fn, working_priors)`

Purpose:

- evaluate the cheap Stage 2 hierarchical correction term

  `sum_i [(1 - beta) log g_i(alpha_i) + beta log p(alpha_i | phi)]`

- without recomputing subject likelihoods.

For the current repo, `p(alpha_i | phi)` can be evaluated from `gaussian_map_fn`.

#### 3. Exact `phi` rejuvenation kernel

File:

- `outer_SMC.R`

New function:

- `rejuvenate_phi_given_alpha_exact(phi, alpha_state, logprior, w, lambda, logprior_phi, gaussian_map_fn, rw_scale = ..., n_moves = ...)`

Purpose:

- update `phi` while holding explicit subject latents fixed;
- use the exact conditional target induced by the current `alpha_1:p`;
- start with RW-MH by copying the current covariance-adaptation pattern from the maintained inner SMC code.

Important:

- do not try to introduce Gibbs on day one;
- first get the exact intercept-only RW kernel working and benchmarked.

#### 4. Exact latent-refresh block

File:

- `outer_SMC.R`

New function:

- `refresh_alpha_block_exact(phi, alpha_state, local_objs, beta, refresh_frac = ..., delayed_accept = TRUE, refresh_ctl = list())`

Purpose:

- choose a subset of subjects to refresh per particle;
- propose new `alpha_i` from the exact refresh proposal;
- use delayed acceptance;
- if the proposal survives screening, evaluate the exact subject likelihood and perform the final MH correction.

This is the core exact move that replaces the current PM likelihood estimator.

#### 5. Sufficient-statistics cache

File:

- `outer_SMC.R`

New functions:

- `compute_alpha_suff_stats(alpha_state, X_list = NULL)`
- `update_alpha_suff_stats_one(stats_obj, alpha_old, alpha_new, x_i = NULL)`

Purpose:

- support cheap `phi` updates and later generalization to `(B, Sigma, X_i)`.

For the current intercept-only model, `X_i = 1`, so the first implementation can use the simpler special case:

- `sum_i alpha_i`
- `sum_i alpha_i alpha_i^T`

Only after this works should the code be generalized to arbitrary `X_i`.

#### 6. Optional latent storage indirection

File:

- `outer_SMC.R`

New function or helper set:

- `materialize_alpha_from_indices()`
- `refresh_alpha_indices_exact()`

Purpose:

- allow the first version to store explicit `alpha` values if that is simpler;
- later swap to bank-index or mixed representations if memory becomes a bottleneck.

For the current problem size, storing explicit `alpha` values is acceptable for the first exact prototype.

### Changes To Experiment Wiring

File:

- `run_two_step_smc.R`

Concrete changes:

1. Add a new experiment mode entry, for example

   - `X1 = list(mode = "extended_exact")`

2. In `run_single_experiment()`:

   - reuse the current pilot-anchor logic;
   - reuse `run_outer_subjects()` unchanged;
   - build `local_objs` using the new `build_subject_exact_local_from_outer()`;
   - call `run_inner_extended_exact()` instead of `run_inner()` for this mode.

3. Add output summaries:

   - runtime split for exact initialization vs exact mutation;
   - average refreshed subjects per round;
   - number of exact local likelihood calls;
   - acceptance rate of delayed-acceptance screens and final exact acceptances.

### Minimal Current-Repo Deliverable

The first exact deliverable should not attempt the full generic `(B, Sigma, X_i)` hierarchy.

Definition of the first shippable exact mode:

- current gamma model only;
- current `phi` parameterization only;
- exact extended-state inner SMC on `(phi, alpha_1:p)`;
- exact initialization via proposal correction;
- exact latent refresh via delayed acceptance;
- RW update of `phi`;
- comparison against Stan using the existing `compare_results.R`.

If that works, the repository has already regained exactness for the current example.

### Phase 2 Generalization

Only after the exact current-model mode works:

1. unpack `phi` into explicit `(B, Sigma)`;
2. allow nontrivial `X_i`;
3. replace the intercept-only sufficient statistics with the full `S_xx`, `S_xa`, `S_aa` form;
4. add conjugate or semi-conjugate `B` and `Sigma` updates where available.

This second phase is a model-generalization step, not the first exactness step.

### Implementation Order

1. Modify `run_outer_subjects()` outputs or wrappers so each subject carries its working prior `g_i`.
2. Add `build_subject_exact_local_from_outer()`.
3. Add exact initialization helpers.
4. Add the new exact inner engine in `outer_SMC.R`.
5. Add the exact `phi` update kernel.
6. Add the exact latent-refresh block with delayed acceptance.
7. Wire the new mode into `run_two_step_smc.R`.
8. Add diagnostics and call counters.
9. Validate against Stan.

### Definition Of Done

The exact mode is done only when all of the following are true:

- the final target is exact for the current continuous gamma model;
- no fixed-bank approximation remains in the target definition;
- all expensive Stage 2 calls are limited to explicit exact initialization and exact latent refresh proposals;
- Stan comparison is at least as good as or better than the current `H1` central tendency on `mu_shift` and `sigma2_shift`;
- runtime and exact likelihood-call counts are reported explicitly.

## Practical Position

If the immediate goal is to get a better scalable inner engine than the current degenerate fixed-cache estimator, these caveats should not stop the prototype.

If the immediate goal is a strict final published method for `phi`, Caveat A cannot be deferred. In that case the right interpretation is:

- build this algorithm as the fast exploratory and adaptation layer
- keep a strict final layer for publication-grade results
