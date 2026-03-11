# Nested SMC Implementation Plan

## Goal

Move from the current disconnected two-stage workflow to one continuously evolving nested SMC algorithm.

In this document:

- `population` means the hyperparameter layer
- `local` means the per-group latent layer
- we use the intercept-only formulation for now

So the maintained target is not yet `X_i B`. It is:

\[
\phi = (\mu, \Sigma)
\]

with

\[
\alpha_i \mid \phi \sim \mathcal N(\mu, \Sigma),
\qquad
y_i \sim p(y_i \mid \alpha_i),
\qquad i = 1,\dots,p.
\]

This is the same nested idea as the original note, just with `X_i = 1`.

---

## Ground Rules

### Keep the framework clean

- the first implementation should be minimal and easy to reason about
- prefer one clear maintained path over multiple modes
- do not add bells and whistles unless they are needed for the core algorithm
- if something is tried and does not work or does not help, remove it afterwards
- dead options and half-maintained side paths should not accumulate

### Keep the population layer thin

The population layer should not try to copy the full complexity of the current local SMC engine.

The local layer is where the richer machinery should live. The population layer should mostly orchestrate:

- adaptive `lambda` choice
- population reweighting
- module resampling
- simple population mutation

### Be honest about progress

Something counts as complete only if it is:

- implemented in the maintained path
- wired into the actual nested algorithm
- sourceable and runnable
- exercised by a reproducible benchmark or check

Something does **not** count as complete if it is only:

- a scaffold
- a helper with no caller
- a scratch script
- an optional mode that is not used
- a partial implementation that cannot run end-to-end

### Optimize the metric that matters

The main efficiency metric is:

- total number of local-level likelihood evaluations

That means evaluations of terms of the form:

\[
\log p(y_i \mid \alpha)
\]

inside local modules, including:

- initialization
- bridge updates
- local rejuvenation
- population-proposal refreshes

This is the number to drive down over time.

Other metrics still matter, but they are secondary:

- wall-clock time
- number of bridge rounds
- number of population resamples
- number of population MH proposals

A method does not count as better just because it uses fewer bridge rounds if it uses more local likelihood evaluations.

### Keep the repository tidy

Core modules should stay in the main folder.

Non-core scripts should not pollute the repo root.

Planned layout:

- `benchmarks/` for stable benchmark scripts
- `dev/scratch/` for exploratory scripts
- `dev/tests/` for repeatable test or verification scripts

---

## 1. Target Family

The full posterior is

\[
\pi_1(\phi,\alpha_{1:p})
\propto
p(\phi)\prod_{i=1}^p p(y_i\mid \alpha_i)\,\mathcal N(\alpha_i;\mu,\Sigma).
\]

Define the bridge

\[
\pi_\lambda(\phi,\alpha_{1:p})
\propto
p(\phi)\prod_{i=1}^p
p(y_i\mid \alpha_i)\,
\mathcal N(\alpha_i;\mu,\Sigma)^\lambda\,
g_i(\alpha_i)^{1-\lambda},
\qquad \lambda \in [0,1].
\]

Here `g_i` is a working local density for group `i`.

Special cases:

- `lambda = 0`: local groups are independent under `p(y_i | alpha_i) g_i(alpha_i)`
- `lambda = 1`: full hierarchical posterior

This is the path the nested algorithm should track continuously.

---

## 2. Architecture and Nested State

The implementation should be organized around:

- a thin `population_driver`
- a mutable `local_module` object

This is the main architectural adjustment to the original plan.

### 2.1 Population driver

The `population_driver` owns:

- the population particles `phi^(1:N)`
- the population weights
- the current bridge value `lambda`
- adaptive step-size selection
- population resampling
- population mutation

It should **not** own complicated local proposal logic directly.

### 2.2 Local module

Each local module is a persistent, changeable object attached to one `(population particle, group)` pair.

It should hold at least:

- local particles `alpha_{i,1:M}`
- local weights
- current `lambda`
- current attached population state `phi`
- cached `log p(y_i | alpha)`
- cached `log g_i(alpha)`
- cached `log N(alpha; phi)`
- local transport / proposal state

The current `enhanced_smc_elite()` output is not the final local module object, but it is the main source for what state and proposal machinery the module should carry.

The local module should expose operations like:

- `predict_increment(phi, lambda_new)`
- `advance(phi, lambda_new)`
- `move_population(phi_new)`
- `clone()`

The first of these must be side-effect free. The others may mutate module state.

### 2.3 Nested state

Each population particle `n = 1, ..., N` carries:

\[
\phi^{(n)} = (\mu^{(n)}, \Sigma^{(n)})
\]

and, attached to it, one local particle system per group:

\[
\mathcal A_i^{(n)} =
\left\{
\alpha_{i,1}^{(n)},\dots,\alpha_{i,M}^{(n)}
\right\}.
\]

Each local system also carries local weights:

\[
\omega_{i,1:M}^{(n)}.
\]

So the full nested state is

\[
\left(
\phi^{(n)},
\left\{
\alpha_{i,1:M}^{(n)},\omega_{i,1:M}^{(n)}
\right\}_{i=1}^p
\right),
\qquad n=1,\dots,N.
\]

This is the natural starting point. Later, parts of this may be compressed, but not at the first implementation.

---

## 3. Local Conditional Target

For fixed population particle `phi^(n)` and bridge level `lambda`, define the local target for group `i`:

\[
\gamma_{i,\lambda}(\alpha \mid \phi^{(n)})
\propto
p(y_i \mid \alpha)\,
\mathcal N(\alpha; \mu^{(n)}, \Sigma^{(n)})^\lambda\,
g_i(\alpha)^{1-\lambda}.
\]

Its normalizing constant is

\[
Z_{i,\lambda}(\phi^{(n)}) =
\int
p(y_i \mid \alpha)\,
\mathcal N(\alpha; \mu^{(n)}, \Sigma^{(n)})^\lambda\,
g_i(\alpha)^{1-\lambda}
\,
d\alpha.
\]

Then the marginal population target is

\[
\overline\pi_\lambda(\phi)
\propto
p(\phi)\prod_{i=1}^p Z_{i,\lambda}(\phi).
\]

This is the core nesting relation:

- local systems estimate `Z_{i,lambda}(phi)`
- population SMC targets the product over groups

---

## 4. Adaptive Bridge Schedule

Let

\[
0 = \lambda_0 < \lambda_1 < \cdots < \lambda_T = 1.
\]

Do not fix `T` in advance.

At step `t`, choose `lambda_{t+1}` adaptively as the largest value with acceptable population ESS under the predicted incremental weights.

This should reuse the same basic logic as the current local CESS code:

- compute particle-wise hypothetical log increments
- apply an ESS / rCESS criterion
- solve for the largest admissible step

The difference is that the hypothetical increments now come from the attached local modules rather than directly from a simple likelihood increment.

For candidate `lambda > lambda_t`, define

\[
\tilde w^{(n)}(\lambda)
=
w_t^{(n)}
\prod_{i=1}^p
u_i^{(n)}(\lambda),
\]

where `u_i^(n)(lambda)` is the estimated local normalizing-constant increment from the current local system.

In code terms, this should come from a side-effect free local-module operation such as:

- `predict_increment(phi = phi^(n), lambda_new = lambda)`

Then compute

\[
\mathrm{ESS}(\lambda)
=
\frac{
\left(\sum_n \tilde w^{(n)}(\lambda)\right)^2
}{
\sum_n \tilde w^{(n)}(\lambda)^2
}.
\]

Choose the largest `lambda_{t+1}` such that

\[
\mathrm{ESS}(\lambda_{t+1}) \ge \rho_{\text{step}} N.
\]

This gives the natural number of bridge updates automatically.

---

## 5. Initialization at Lambda = 0

At `lambda = 0`, the local targets are

\[
\gamma_{i,0}(\alpha)
\propto
p(y_i\mid \alpha)\,g_i(\alpha).
\]

These do not depend on `phi`.

Initialization:

1. draw `phi_0^(n) ~ p(phi)` for `n = 1, ..., N`
2. for each population particle `n` and group `i`, initialize a local module targeting `gamma_{i,0}`
3. set `w_0^(n) = 1/N`

This is the correct place to exploit the current local SMC machinery and local working densities.

The intent is not to rerun the full local SMC engine wholesale at every later step. The intent is:

- use the current local engine to help build good initial local modules
- then evolve those modules continuously as `lambda` and `phi` change

---

## 6. Population Incremental Weighting

Suppose we move from `lambda_t` to `lambda_{t+1}` with `Delta lambda = lambda_{t+1} - lambda_t`.

For group `i` and population particle `n`, under current local particles targeting `gamma_{i,\lambda_t}(. | phi^(n))`, the natural incremental factor is

\[
u_{i,t+1}^{(n)}
=
\frac{1}{M}
\sum_{m=1}^M
\left[
\frac{
\mathcal N(\alpha_{i,m}^{(n)};\mu^{(n)},\Sigma^{(n)})
}{
g_i(\alpha_{i,m}^{(n)})
}
\right]^{\Delta\lambda}.
\]

Then update population weights by

\[
w_{t+1}^{(n)}
\propto
w_t^{(n)}
\prod_{i=1}^p u_{i,t+1}^{(n)}.
\]

In log form:

\[
\log w_{t+1}^{(n)}
=
\log w_t^{(n)}
+
\sum_{i=1}^p \log u_{i,t+1}^{(n)}.
\]

This is the main weighting rule for the nested implementation.

In code, the population driver should obtain these local increment estimates by asking each local module for its predicted increment under the candidate `lambda`.

---

## 7. Population Resampling

After reweighting, compute

\[
\mathrm{ESS}_{t+1}
=
\frac{1}{\sum_n (w_{t+1}^{(n)})^2}.
\]

If

\[
\mathrm{ESS}_{t+1} < \rho_{\text{res}} N,
\]

resample the population particles and duplicate their attached local systems.

The resampling object is the whole module:

\[
\left(
\phi^{(n)},
\left\{
\mathcal A_i^{(n)}
\right\}_{i=1}^p
\right).
\]

Use systematic or stratified resampling.

After resampling, reset weights to `1/N`.

This is where the `clone()` behavior of the local module matters.

---

## 8. Mutation Step

After each bridge move, apply a mutation kernel targeting the new `pi_{lambda_{t+1}}`.

This mutation has two parts.

### 8.1 Population mutation

For each population particle `n`, propose

\[
\phi' \sim K_\phi(\cdot \mid \phi^{(n)}).
\]

For the intercept-only first version:

- `phi = (mu, log sigma2)` in the diagonal case
- population proposals can be random-walk MH on `mu` and `log sigma2`
- later this can be replaced by better blocked updates if needed

The acceptance ratio is based on the estimated product of local normalizing constants:

\[
a
=
\min\left(
1,
\frac{
p(\phi') \prod_{i=1}^p Z_{i,\lambda}(\phi')
}{
p(\phi^{(n)}) \prod_{i=1}^p Z_{i,\lambda}(\phi^{(n)})
}
\cdot
\frac{K_\phi(\phi^{(n)}\mid \phi')}{K_\phi(\phi' \mid \phi^{(n)})}
\right).
\]

In practice, those `Z_i` terms are estimated by lightly refreshing the local systems under `phi'`.

This is another place where the population layer should stay thin:

- propose a new `phi'`
- ask the attached local modules for a cheap evidence update under `phi'`
- accept or reject

Do not build a second large standalone SMC engine for the population layer unless the simple version clearly fails.

### 8.2 Local mutation

For each group `i` within population particle `n`, update the local system from

\[
\gamma_{i,\lambda_t}(\cdot \mid \phi_{\text{old}}^{(n)})
\]

to

\[
\gamma_{i,\lambda_{t+1}}(\cdot \mid \phi_{\text{new}}^{(n)}).
\]

The natural local update is:

1. reweight local particles by the incremental ratio
2. resample locally if local ESS is low
3. apply one or a few local rejuvenation moves targeting the new local target

The exact incremental ratio is

\[
v_{i,m}^{(n)}
=
\frac{
\gamma_{i,\lambda_{t+1}}(\alpha_{i,m}^{(n)} \mid \phi_{\text{new}}^{(n)})
}{
\gamma_{i,\lambda_t}(\alpha_{i,m}^{(n)} \mid \phi_{\text{old}}^{(n)})
}.
\]

Under the current bridge this becomes the change in the Gaussian-vs-working-prior part, plus any change induced by `phi`.

Algorithmically, this should live inside the local module, not in the population driver.

---

## 9. Full Algorithm

### Inputs

- `N`: number of population particles
- `M`: number of local particles per group per population particle
- `g_i`: working local density for each group
- `rho_step`: ESS threshold for choosing the next `lambda`
- `rho_res`: ESS threshold for resampling population particles
- `rho_local`: ESS threshold for resampling local particles

### Initialization

For `n = 1, ..., N`:

1. draw `phi_0^(n) ~ p(phi)`
2. for each group `i`, initialize a local module targeting `p(y_i | alpha_i) g_i(alpha_i)`
3. set population weight `w_0^(n) = 1/N`

Set `lambda_0 = 0`.

### Main loop

Repeat until `lambda_t = 1`:

1. choose `lambda_{t+1}` adaptively using predicted population ESS
2. for each population particle and group:
   advance the attached local module from `lambda_t` to `lambda_{t+1}`
3. update population weights using the product of estimated local normalizing-constant increments
4. if population ESS is low, resample whole population modules
5. mutate each population particle
6. after a population move, lightly refresh local modules under the accepted population state
7. set `t <- t + 1`

### Output

- weighted population approximation to `p(phi | y)`
- attached local systems approximating `p(alpha_i | phi, y_i)`

---

## 10. What Already Exists

The current codebase already contains useful pieces for this nested plan.

It is important to distinguish carefully between what already helps and what actually counts as nested-SMC progress.

### 10.1 Local SMC engine

`SMC_super_fast.R`

- `enhanced_smc_elite()`

This is a maintained per-group SMC sampler. It already provides:

- adaptive bridge stepping by CESS
- transport-map fitting
- local resampling
- local rejuvenation
- likelihood caching

It is not the final nested local updater as-is, but it is the main source of local proposal machinery and bridge design ideas.

What this **does** count as:

- a proven single-group SMC engine
- a proven source of transport and rejuvenation ideas
- a proven source of local proposal state

What this **does not** count as:

- a completed local module for the nested algorithm
- any completed population-level nested logic

### 10.2 Shared SMC numerics

`smc_core.R`

Useful functions already available:

- `ESS()`
- `rCESS()`
- `rCESS_stat()`
- `next_lambda_via_rCESS()`
- `next_lambda_via_rCESS_stat()`
- `residual_resample()`
- `stratified_resample_sorted()`
- `ll_parallel()`
- Gaussian-mixture helper code for local proposals

The population bridge-step solver for the nested algorithm will likely be built on the same `rCESS_stat()` pattern, but with nested incremental log-weights.

### 10.3 Population prior and Gaussian map

`make_prior.R`

Useful functions already available:

- `make_prior_phi_diag()`
- `phi_to_gaussian_params_diag_factory()`

These already support the intercept-only diagonal-Gaussian population layer:

- `phi = (mu, log sigma2)`
- `alpha_i | phi ~ N(mu, diag(sigma2))`

This is the right starting point for the first implementation.

### 10.4 Local working densities and local proposal objects

`new_SMC_cache.R`

Useful functions already available:

- `build_local_surrogate()`
- `build_local_exact_object()`
- `log_working_prior_gaussian_mat()`
- `make_local_theta_proposal_init()`
- `make_local_theta_proposal_refresh()`
- `make_local_theta_proposal()`

These already give:

- local working Gaussian densities `g_i`
- local proposal draws for initialization
- local proposal draws for refresh under a new population state

These are directly relevant to the nested implementation.

What this **does** count as:

- working local-density and local-proposal support code

What this **does not** count as:

- a complete local module lifecycle
- a complete nested local update kernel

### 10.5 Explicit population/local helper skeleton

`hierarchical_state_helpers.R`

Useful functions already available:

- `compute_local_suff_stats()`
- `update_local_suff_stats_one()`
- `log_tempered_population_contrib()`
- `init_population_particles_exact()`
- `rejuvenate_population_given_local_exact()`
- `refresh_local_block_exact()`

These are not yet the final nested driver, but they already encode:

- explicit local-state storage
- intercept-only sufficient statistics
- exact initialization from local proposal objects
- a population mutation step conditional on local states
- a local refresh step conditional on a population particle

This file is the closest existing scaffold for the nested implementation.

What this **does** count as:

- useful intercept-only scaffolding
- proof that explicit population/local state can be represented

What this **does not** count as:

- a finished nested algorithm
- a finished population driver
- a finished local module abstraction

---

## 11. What Is Missing

The current code does **not** yet implement the true continuously evolving nested SMC algorithm.

The main missing algorithmic pieces are:

### 11.1 A real thin population driver

Need a top-level driver, for example:

- `nested_population_smc()`

This function should own:

- population particles
- attached local modules
- adaptive bridge stepping
- resampling of whole population modules
- mutation of both population and local states

Nothing in the current codebase yet does that end-to-end.

### 11.2 A real local module abstraction

The current helper code stores local states, but not yet a proper module object with a stable interface.

Need a real `local_module` representation containing:

- current particles
- current local weights
- current `lambda`
- current attached population state
- cached local log-likelihoods
- cached `log g_i(alpha)`
- cached `log N(alpha; mu, Sigma)`
- transport / proposal state

and methods:

- `predict_increment()`
- `advance()`
- `move_population()`
- `clone()`

### 11.3 Population-step lambda solver using nested increments

Need a solver for `lambda_{t+1}` based on predicted population ESS where the hypothetical incremental weights are built from current local systems:

- compute `u_i^(n)(lambda)` from each local system
- aggregate across groups
- solve for largest feasible `lambda`

This should deliberately reuse the style of `rCESS_stat()` rather than inventing a new step-size logic.

### 11.4 Local incremental update kernel

Need a dedicated function, for example:

- `advance_local_module()`

For one group `i` and one population particle `n`, it should:

1. compute local incremental log weights
2. estimate the local normalizing-constant ratio
3. update local weights
4. resample if local ESS is low
5. rejuvenate locally under the new target

This is a central missing component.

### 11.5 Population weight update from local estimates

Need a dedicated function, for example:

- `compute_population_increment()`

Given all local normalizing-constant increment estimates, return:

- incremental log weights
- predicted ESS
- diagnostics for the lambda selection step

### 11.6 Population resampling of full modules

Need a clean module-resampling function that duplicates:

- one population particle
- all attached local systems

The general resampling utilities already exist, but the nested state wrapper does not.

### 11.7 Consistent population MH with refreshed local evidence

`rejuvenate_population_given_local_exact()` is a useful scaffold, but the true nested version should use current local modules to estimate the required local normalizing-constant ratios under proposed `phi'`.

So the final population mutation step still needs a proper nested implementation.

### 11.8 Diagnostics and bookkeeping

Need to store:

- `lambda` history
- population ESS history
- local ESS summaries
- number of local resamples
- number of local rejuvenation moves
- number of population acceptances
- total local-level likelihood evaluations
- local-level likelihood evaluations split by source:
  initialization, bridge updates, local rejuvenation, population-proposal refreshes

This is necessary for tuning and for comparing against the old workflow.

The primary efficiency report for any run should include:

- total local-level likelihood evaluations
- local-level likelihood evaluations per completed run to `lambda = 1`
- local-level likelihood evaluations by phase of the algorithm

---

## Honest Status Review

Current status, without inflation:

### What genuinely exists

- a real single-group local SMC engine
- real local working-prior and proposal helpers
- real intercept-only population prior helpers
- a useful explicit-state scaffold in `hierarchical_state_helpers.R`
- a simple local-vs-Stan benchmark in `benchmarks/`

### What does not yet exist

- a working continuously evolving nested SMC driver
- persistent local modules attached to each population particle
- adaptive population bridge stepping driven by local-module predictions
- a module-level population resampling mechanism
- a correct end-to-end nested benchmark

### What the main score should be

For this project, the main score is not just:

- number of bridge rounds
- number of refresh events

The main score is:

- total local-level likelihood evaluations needed to reach `lambda = 1`

That is the quantity future changes should be judged against.

### What should not be overstated

- `hierarchical_state_helpers.R` is scaffold, not completion
- `nested_SMC.MD` is a plan, not progress
- a local benchmark is useful validation, but it is not nested-SMC progress by itself
- any old two-stage script does not count toward completion of the new method

---

## Repository Hygiene

The repo should stay clean while this is being built.

### Root folder policy

Only keep these in the main folder:

- core maintained modules
- minimal project entry points that are genuinely part of the maintained workflow

### Non-core script policy

All non-core scripts should live under dedicated folders:

- `benchmarks/` for stable benchmark scripts
- `dev/scratch/` for exploratory scripts
- `dev/tests/` for repeatable checks that are not core modules

That means:

- no new ad hoc scratch script in the root
- no temporary debug runner in the root
- no benchmark helper dumped in the root

If a script stops being useful, delete it rather than leaving it around.

---

## 12. Recommended Implementation Order

The safest implementation order is:

### Phase 1. Intercept-only exact skeleton

Build the smallest correct nested driver using:

- `phi = (mu, log sigma2)`
- diagonal Gaussian `Sigma`
- local working Gaussian densities
- explicit local particle systems per population particle and group

At this stage, correctness matters more than efficiency.

Needed deliverables:

- thin population driver
- real local-module object
- local-module update function
- population lambda-step solver
- population reweight/resample/mutate loop

### Phase 2. Reuse current local proposal machinery

For local initialization and refresh proposals, reuse:

- `build_local_exact_object()`
- `make_local_theta_proposal_init()`
- `make_local_theta_proposal_refresh()`

This avoids inventing new proposal code too early.

### Phase 3. Replace scaffold refresh logic with true local-module updates

The current helper refresh logic is only scaffold.

Replace it with a true local nested update that:

- tracks local weights continuously in `lambda`
- estimates local normalizing-constant increments
- uses local ESS-based resampling

### Phase 4. Improve population mutation

Once the bridge and local updates are correct:

- improve the population proposal
- possibly exploit intercept-only sufficient statistics
- reduce the number of local refreshes needed per population proposal

Only add complexity here if the simple population MH step is a real bottleneck.

### Phase 5. Generalize beyond intercept-only

Only after the intercept-only nested method is stable:

- replace `mu` by `X_i B`
- replace current diagonal helpers as needed

This should be a later extension, not part of the first implementation.

---

## Phase Completion Criteria

These are the standards for claiming that a phase is done.

### Phase 1 counts as complete only if

- the nested driver can run end-to-end on a small intercept-only example
- the run reaches `lambda = 1`
- population resampling duplicates full modules correctly
- local modules persist across bridge steps rather than being thrown away each time
- population mutation is wired into the maintained path
- the run reports total local-level likelihood evaluations

Phase 1 does **not** count as complete if there is only:

- a driver skeleton
- a lambda solver with no real caller
- a module structure with no bridge updates
- a benchmark that still relies on the old disconnected workflow

### Phase 2 counts as complete only if

- initialization and refresh proposals from the current local machinery are actually used by the nested driver
- the maintained path can build local modules from those proposal helpers

Phase 2 does **not** count as complete if the helpers merely exist but the nested driver does not rely on them.

### Phase 3 counts as complete only if

- local weights evolve continuously with `lambda`
- local ESS triggers local resampling
- local rejuvenation is applied inside the module update
- the driver no longer depends on scaffold refresh logic for the maintained path

Phase 3 does **not** count as complete if local modules are effectively rebuilt from scratch each step under a different name.

### Phase 4 counts as complete only if

- the simple population mutation baseline already works
- the improved mutation is shown to help the total local-level likelihood evaluation count, or at least wall-clock time without increasing that count badly, in a reproducible check
- the added complexity is kept in the maintained path only if it earns its keep

Phase 4 does **not** count as complete if it only adds knobs, modes, or proposal variants without evidence.

### Phase 5 counts as complete only if

- the intercept-only path still works after generalization
- the generalized path is exercised on a reproducible example
- the new complexity is still compatible with the thin population-driver design

Phase 5 does **not** count as complete if the generalization breaks the intercept-only maintained path.

---

## 13. Proposed First Driver API

The first concrete maintained objects should be:

- one `local_module`
- one `nested_population_smc()` driver

### 13.1 Local module API

The local module should support functions along the lines of:

```r
build_local_module <- function(local_obj, phi, lambda = 0, M_local, seed = 123)
predict_increment_local_module <- function(module, phi, lambda_new)
advance_local_module <- function(module, phi, lambda_new, rho_local, n_local_moves)
move_population_local_module <- function(module, phi_new, n_local_moves = 1L)
clone_local_module <- function(module)
```

The contract should be:

- `predict_increment_*` is side-effect free
- `advance_*` mutates the module state after the new `lambda` is chosen
- `move_population_*` adjusts the module after a population move
- `clone_*` makes resampling straightforward

### 13.2 Population driver API

The first population driver should look something like:

```r
nested_population_smc <- function(
    rprior_phi,
    local_objs,
    logprior_phi,
    gaussian_map_fn,
    N,
    M_local,
    rho_step = 0.8,
    rho_res = 0.5,
    rho_local = 0.5,
    n_population_moves = 1L,
    n_local_moves = 1L,
    seed = 123,
    verbose = TRUE
)
```

Expected inputs:

- `rprior_phi`: prior sampler for population particles
- `local_objs`: one local object per group, containing data, likelihood, and working density
- `logprior_phi`: population prior density
- `gaussian_map_fn`: map from population `phi` to local Gaussian parameters

Expected output:

- population particles and weights
- full nested local systems
- diagnostics and histories

---

## 14. Minimal New Functions To Add

The first implementation should add only a small number of new functions:

- `build_local_module()`
- `predict_increment_local_module()`
- `advance_local_module()`
- `move_population_local_module()`
- `clone_local_module()`
- `predict_population_log_increment()`
- `choose_next_lambda_nested()`
- `resample_population_modules()`
- `mutate_population_particle_nested()`
- `nested_population_smc()`

Everything else should be reused from the current code where possible.

---

## 15. Summary

The current codebase already has:

- a strong local SMC engine
- the working-prior machinery needed for `g_i`
- intercept-only population priors and Gaussian maps
- an exact population/local scaffold

What it does **not** yet have is the actual nested algorithm:

- persistent local modules attached to each population particle
- adaptive population bridge stepping based on nested ESS
- population reweighting from local normalizing-constant increments
- consistent module resampling
- a full continuously evolving nested SMC driver

So the immediate task is not to redesign everything. It is to build the smallest correct nested implementation around a thin `population_driver` and a rich mutable `local_module`, in the intercept-only case first, and then measure it primarily by total local-level likelihood evaluations.
