# Top-Level Support Engine Plan

## Problem Statement

The hierarchical target is

```math
p(\theta \mid y)
\propto
p(\theta)
\prod_{i=1}^S m_i(\theta),
```

where

```math
m_i(\theta)
=
\int L_i(\alpha_i)p(\alpha_i \mid \theta)d\alpha_i.
```

The current bank approach tries to discover useful population values by first
building local alpha banks. That is backwards for support discovery. The hard
object is the posterior-relevant theta region

```math
\mathcal R_\epsilon
=
\left\{
\theta:
\log p(\theta)+\sum_i \log m_i(\theta)
\ge
\max_\vartheta
\left[
\log p(\vartheta)+\sum_i \log m_i(\vartheta)
\right]
-\epsilon
\right\}.
```

The implementation should therefore build a conservative top-level support
engine first, then spend local SMC only inside that discovered region.

## Design Direction

The new architecture is:

```text
local likelihood sketch
-> analytic approximate local marginal
-> conservative theta support discovery
-> expensive local SMC correction
-> final theta SMC
-> certification at posterior functional stress points
```

The local bank is not responsible for discovering theta support. It can remain a
baseline, but the next architecture should move support discovery to the
population level.

## Recent Methods To Use

### 1. INLA-like support map

For each subject, approximate the local likelihood by a Gaussian or a small
Gaussian mixture:

```math
L_i(\alpha)
\approx
\widetilde L_i(\alpha)
=
\sum_{k=1}^{K_i}
c_{ik}N(\alpha \mid a_{ik},V_{ik}).
```

For a normal population prior,

```math
p(\alpha \mid \theta)=N(\alpha \mid \mu,\Sigma),
```

the approximate marginal is analytic:

```math
\widetilde m_i(\theta)
=
\int \widetilde L_i(\alpha)p(\alpha \mid \theta)d\alpha
=
\sum_k
c_{ik}N(a_{ik}\mid \mu,\Sigma+V_{ik}).
```

Then

```math
\widetilde \ell(\theta)
=
\log p(\theta)+\sum_i \log \widetilde m_i(\theta)
```

is cheap to evaluate. Use it for support discovery, not final correctness.

Make the support map deliberately conservative:

```math
\widetilde \ell_\rho(\theta)
=
\log p(\theta)
+
\rho\sum_i \log \widetilde m_i(\theta),
\qquad
0<\rho\le 1.
```

Small rho gives wider support. A useful SMC sequence is

```text
rho = 0.2 -> 0.4 -> 0.7 -> 1.0
```

Recent references:

- Laplace Matching for latent Gaussian models, 2021:
  https://arxiv.org/abs/2105.03109
- Fast scalable approximations for extended latent Gaussian models, 2021:
  https://arxiv.org/abs/2103.07425
- INLA extensions for double hierarchical models, 2022:
  https://link.springer.com/article/10.1007/s11222-022-10122-1

### 2. Pathfinder-style theta support proposals

Run multiple quasi-Newton optimizations on the cheap approximate target

```math
\widetilde \ell(\theta).
```

Each path gives a local Gaussian proposal

```math
q_r(\theta)=N(\theta\mid \hat\theta_r,H_r^{-1}).
```

Use the mixture

```math
q_{\mathrm{PF}}(\theta)=\sum_r\omega_r q_r(\theta)
```

as a support proposal, including inflated components

```math
N(\theta\mid \hat\theta_r,c^2H_r^{-1}), \qquad c>1.
```

This is not final inference. It is a systematic way to find ridges, skewed
regions, and tail directions in theta-space.

Reference:

- Pathfinder, Zhang, Carpenter, Gelman, Vehtari, JMLR 2022:
  https://jmlr.org/papers/v23/21-0889.html

### 3. Active correction over theta

After the cheap support map is built, estimate the correction

```math
\Delta(\theta)
=
\ell(\theta)-\widetilde\ell(\theta)
=
\sum_i
\left[
\log m_i(\theta)-\log \widetilde m_i(\theta)
\right].
```

Use local SMC only at selected theta points. Then fit a correction model over
theta:

```math
\Delta(\theta)\sim \mathcal{GP}
```

or a simpler local regression if the dimension is too high for a GP.

The acquisition rule must be posterior-focused:

```math
a(\theta)
=
\widetilde p(\theta\mid y)
\cdot
\operatorname{Var}\{\widehat\Delta(\theta)\}
\cdot
s(\theta),
```

where `s(theta)` is high for tail and posterior functional stress points.

Recent references:

- GPry, fast Bayesian inference with GP log-posterior surrogate, 2022:
  https://arxiv.org/abs/2211.02045
- Constrained GP active learning for expensive Bayesian inverse problems, 2023:
  https://arxiv.org/abs/2312.08085

### 4. Transport SMC as an optional later upgrade

Flow or transport SMC can help move theta particles through complicated support:

```math
T_t\#N(0,I)\approx \pi_t(\theta).
```

This is useful only after there is a reliable top-level approximate target. It
does not solve noisy local evidence by itself.

Relevant work:

- Annealed Flow Transport Monte Carlo, 2021:
  https://arxiv.org/abs/2102.07501
- Flow Annealed Importance Sampling Bootstrap, 2022:
  https://arxiv.org/abs/2208.01893
- Transport Score Climbing, 2022:
  https://arxiv.org/abs/2202.01841
- flowMC, 2023:
  https://joss.theoj.org/papers/10.21105/joss.05021

This is not phase one. It is a later upgrade if vanilla theta SMC still has
transport problems after the support map exists.

## Clean Implementation Plan

The implementation must start from the clean bank-SMC baseline. Do not build on
the failed residual-corrector or bridge-graph experiments.

Each phase below has a narrow deliverable and a hard progress gate. If a phase
fails, first classify the failure as:

```text
implementation bug
tuning / numerical stability issue
bad approximation family
bad architecture
```

Do not discard a method after one bad posterior plot. Do not keep a method after
it fails its own diagnostic gates.

### Phase 0: Baseline and cleanup

Goal: define the clean starting point.

Tasks:

1. Confirm the working tree contains only the clean bank-SMC baseline.
2. Run the shifted-gamma benchmark once with current defaults.
3. Record:
   - elapsed time;
   - factor size;
   - theta posterior quantile errors versus Stan;
   - known weak points.
4. Do not add new abstractions in this phase.

Progress gate:

```text
The baseline benchmark must be reproducible from one script.
The result file and plot path must be recorded.
```

Failure handling:

```text
If the clean baseline cannot be reproduced, fix that first.
Do not start the new support-engine work until the baseline is stable.
```

### Phase 1: Local likelihood sketch object

Goal: represent a cheap approximation to `L_i(alpha)` independently of theta.

Object:

```text
local_likelihood_sketch
  local_id
  alpha_names
  components:
    weight_log_c[k]
    mean_alpha[k,]
    cov_alpha[k,,]
  diagnostics:
    optimizer status
    Hessian condition number
    local mode log likelihood
    profile checks
```

Start with one Gaussian component:

```math
\widetilde L_i(\alpha)=c_iN(\alpha\mid a_i,V_i).
```

Use a mixture only after the one-component path is verified.

Implementation steps:

1. Add a function to fit one local mode for subject `i`.
2. Compute the Hessian or observed curvature at the mode.
3. Regularize the covariance only when the Hessian is numerically invalid.
4. Store the normalizing constant term `c_i` honestly. If it is approximate,
   name it as approximate.
5. Add a diagnostic that evaluates the sketch against direct local likelihood
   values along the main local directions.

Progress gate:

```text
For each subject, the sketch must place its mode near a high-likelihood alpha.
The covariance must be positive definite.
Sketch log density profiles must be directionally sensible.
```

Failure handling:

```text
If profiles are bad, inspect whether the local likelihood is skewed,
bounded, or multimodal.
If the issue is skew/multimodality, add a second component deliberately.
Do not jump directly to a flexible mixture without proving why one component
failed.
```

### Phase 2: Analytic theta marginal from sketches

Goal: compute

```math
\widetilde m_i(\theta)
=
\sum_k c_{ik}N(a_{ik}\mid \mu,\Sigma+V_{ik})
```

for all subjects.

Implementation steps:

1. Implement the analytic convolution for diagonal Gaussian population priors.
2. Keep the interface specific and honest:

   ```text
   build_sketch_factor_set(sketch_set, population_model)
   ```

3. The factor set evaluates

   ```text
   theta -> sum_i log m_tilde_i(theta)
   ```

4. Avoid storing per-theta values. Evaluation should be vectorized over theta.

Progress gate:

```text
For simple normal-normal test cases, analytic sketch factors must match the
known marginal likelihood.
For shifted gamma, selected theta values must be compared against fresh local
SMC to measure approximation error, not to tune the final posterior yet.
```

Failure handling:

```text
If the analytic formula fails on the normal-normal case, it is a bug.
If it passes normal-normal but is poor on shifted gamma, it is approximation
error. Diagnose which subjects and which theta directions cause it.
```

### Phase 3: Conservative theta support SMC

Goal: sample a conservative approximate theta posterior before local SMC
correction.

Target:

```math
\widetilde\pi_\rho(\theta)
\propto
p(\theta)
\prod_i \widetilde m_i(\theta)^\rho.
```

Implementation steps:

1. Reuse the existing outer SMC where possible.
2. Add a simple `rho` tempering wrapper around the sketch factor set.
3. Run SMC for a sequence:

   ```text
   rho = 0.2, 0.4, 0.7, 1.0
   ```

4. Save theta particles and weights at each rho.
5. Summarize support expansion:
   - marginal quantiles;
   - covariance;
   - tail ranges;
   - ESS;
   - number of SMC rounds.

Progress gate:

```text
The rho=1 approximate posterior must be cheap and stable.
The rho<1 particles must form a wider support envelope than rho=1.
For shifted gamma diagnostics only, Stan posterior tails should lie inside the
rho<1 support envelope. This is not allowed as an algorithmic dependency.
```

Failure handling:

```text
If rho<1 is not wider, the factor tempering is implemented incorrectly.
If the support is still too narrow, inflate the local sketch covariances before
adding more algorithmic machinery.
```

### Phase 4: Pathfinder-style support proposals

Goal: improve support discovery on ridges and skewed regions of theta.

Implementation steps:

1. Run multiple optimizers on

   ```math
   \widetilde\ell(\theta)
   ```

   from dispersed starts.

2. Store for each path:
   - final theta;
   - approximate inverse Hessian;
   - objective value;
   - convergence code;
   - gradient norm if available.

3. Construct Gaussian probability-contour proposal components:

   ```math
   \{\theta:(\theta-\hat\theta_r)^\top H_r(\theta-\hat\theta_r)
   =\chi^2_d(p)\}
   ```

   with several probability levels.

4. Merge Pathfinder candidates with SMC candidates from Phase 3.

Progress gate:

```text
The Pathfinder proposal must add non-duplicate tail/ridge candidates not already
covered by rho-SMC.
It must not replace rho-SMC.
```

Failure handling:

```text
If optimizers collapse to one point, check starts and parameter scaling.
If Hessians are unstable, use covariance from optimizer path history or the
rho-SMC covariance instead of pretending the Hessian is valid.
```

### Phase 5: Expensive local SMC correction points

Goal: choose a small set of theta points for accurate local SMC evaluation.

Candidate sources:

```text
rho-SMC posterior center
rho-SMC posterior tails
inflated rho-SMC draws
Pathfinder components
prior stress points
subject-profile stress points
```

For each candidate theta, estimate

```math
\widehat\Delta(\theta)
=
\sum_i
\left[
\log \widehat m_i(\theta)-\log \widetilde m_i(\theta)
\right].
```

Implementation steps:

1. Start with a small fixed evaluation design, for example 10 to 20 theta
   points.
2. At each theta, run local SMC for every subject.
3. Use replicated local SMC at a subset of points to estimate log evidence
   noise and bias.
4. Store correction observations as theta-level objects:

   ```text
   theta_correction_observation
     theta
     log_m_smc_by_subject
     log_m_tilde_by_subject
     delta_total
     delta_by_subject
     mcse_by_subject
     replicate diagnostics
   ```

Progress gate:

```text
At selected theta points, correction MCSE must be small enough that posterior
ranking errors are interpretable.
If aggregate correction uncertainty is larger than 1-2 log units in relevant
regions, increase local SMC precision before fitting any correction surrogate.
```

Failure handling:

```text
If replicated local SMC disagrees with single-run estimates, the issue is local
evidence precision. Do not fit a correction surrogate to noisy logZ values.
If only a few subjects dominate the error, diagnose those subjects separately.
```

### Phase 6: Theta-level correction surrogate

Goal: model only

```math
\Delta(\theta)
```

over theta. Do not model per-subject alpha banks.

Start simple:

```text
local linear / quadratic regression in whitened theta coordinates
```

Only add GP after the regression baseline is understood.

Implementation steps:

1. Whiten theta using the conservative support covariance from Phase 3.
2. Fit a correction model with observation noise from Phase 5.
3. Validate by leave-one-out prediction over correction points.
4. Compare:
   - constant correction;
   - linear correction;
   - quadratic correction;
   - GP correction, only if needed.

Progress gate:

```text
The selected surrogate must win by held-out predictive error, not posterior
appearance.
If constant correction wins, keep constant correction.
```

Failure handling:

```text
If all surrogates fail, do not tune blindly. Add correction points in the region
with highest posterior relevance and highest correction uncertainty.
```

### Phase 7: Corrected top-level SMC

Goal: sample

```math
\pi_{\mathrm{corr}}(\theta)
\propto
\exp\{
\widetilde\ell(\theta)+\widehat\Delta(\theta)
\}.
```

Implementation steps:

1. Build a corrected analytic factor set that evaluates the sketch likelihood
   plus the theta correction.
2. Run outer SMC on this corrected target.
3. Compare against:
   - sketch-only posterior;
   - clean bank baseline;
   - Stan benchmark where available.

Progress gate:

```text
The corrected posterior must improve over sketch-only on held-out fresh-SMC
diagnostics.
It does not need to beat the clean bank baseline on the first run.
It must not be slower and worse without a clear diagnostic reason.
```

Failure handling:

```text
If corrected SMC worsens posterior shape, inspect correction surrogate
validation before changing SMC settings.
If SMC has low ESS or bad transport, then consider transport/flow SMC later.
```

### Phase 8: Active refinement

Goal: add expensive local SMC evaluations only where they can change posterior
functionals.

Acquisition:

```math
a(\theta)
=
\widetilde p_{\mathrm{corr}}(\theta\mid y)
\cdot
\operatorname{Var}\{\widehat\Delta(\theta)\}
\cdot
s(\theta).
```

Implementation steps:

1. Define posterior functionals:
   - marginal quantiles;
   - tail probabilities;
   - evidence if needed.
2. Draw candidate theta points from the corrected posterior and inflated support
   proposals.
3. Score candidates by acquisition.
4. Add a small batch of correction points.
5. Refit the correction surrogate.
6. Rerun corrected SMC.

Progress gate:

```text
Stop when relevant posterior functionals are stable under added correction
points.
Stability means changes are below a predeclared tolerance, not merely visually
acceptable.
```

Failure handling:

```text
If refinement keeps adding points in the same region, the local sketch is too
weak there. Improve the sketch for the responsible subjects.
If refinement jumps across unrelated regions, the support proposal is too
narrow or the correction uncertainty model is wrong.
```

### Phase 9: Certification

Goal: prove the approximate posterior is trustworthy enough for the benchmark.

Certification points:

```text
posterior center
each marginal 5%, 50%, 95% theta quantile
inflated tail points
Pathfinder ridge points
points with high correction uncertainty
points with high functional sensitivity
```

At each point, run fresh local SMC with enough precision to estimate

```math
\ell(\theta)-\widehat\ell_{\mathrm{corr}}(\theta).
```

Progress gate:

```text
For posterior-relevant theta points, aggregate log posterior error must be small
enough not to change declared posterior functionals.
For irrelevant theta points, large absolute error is acceptable only if fresh
log posterior remains far below the posterior region.
```

Failure handling:

```text
If certification fails in a posterior-relevant region, add correction points
there and rerun active refinement.
If certification fails because fresh local SMC itself is unstable, improve local
SMC precision or local proposal quality before changing the top-level method.
```

## Anti-Bloat Rules

1. Do not keep a new abstraction unless it survives a progress gate.
2. Do not keep bridge-graph, MBAR, or residual-bank code in this path.
3. Do not add command-line knobs for every internal choice.
4. Expose only:
   - sketch complexity;
   - support inflation;
   - correction budget;
   - SMC particle budgets.
5. Every stored result must be reloadable and sufficient to continue the run.
6. Do not save full local SMC particles by default unless they are required for
   continuation.

## First Minimal Experiment

The first real experiment should be deliberately small:

```text
shifted gamma benchmark
one Gaussian local sketch per subject
analytic sketch factor
rho-SMC support discovery
fresh local SMC correction at 10 theta points
constant vs linear correction comparison
corrected theta SMC
certification at 5 theta points
```

Success is not defined as matching Stan immediately. Success is:

```text
The sketch posterior finds the relevant theta support.
The correction estimates are numerically stable.
The correction surrogate passes held-out validation.
The corrected posterior moves in the direction predicted by fresh SMC audits.
```

If this fails, the next action is determined by diagnostics:

```text
bad support envelope -> improve/inflate sketches
bad correction MCSE -> improve local SMC precision
bad surrogate validation -> add correction points or simplify surrogate
bad corrected SMC transport -> consider transport SMC
```

No posterior plot alone is allowed to decide the architecture.

## Phase 0 Result: Clean Baseline

Date: 2026-05-21.

Code state:

```text
git status --short
?? update.md
```

There are no code diffs relative to `HEAD`; the only worktree change is this
planning document. The benchmark therefore ran from the clean bank-SMC baseline.

Command:

```text
Rscript benchmarks/run_shifted_gamma_hierarchy_compare.R
```

Benchmark:

```text
script:  benchmarks/run_shifted_gamma_hierarchy_compare.R
label:   bank_smc_full
seed:    20260519
cores:   4
data:    20 subjects x 200 trials
local:   5 bank nodes x 600 particles
```

Artifacts:

```text
results: /Users/nstevenson/Documents/2025/FastHierarchical/benchmarks/results/shifted_gamma_bank_smc_full_results.rds
plot:    /Users/nstevenson/Documents/2025/FastHierarchical/benchmarks/results/shifted_gamma_bank_smc_full_posteriors.png
```

Timing and size:

```text
elapsed_sec:             39.88014
bank nodes per subject:  min=5, median=5, max=5
particles per subject:   min=3000, median=3000, max=3000
```

Audit summary:

```text
round  failures  min_ess_frac  median_ess_frac
1      20        0.00879096    0.1155578
2      1         0.04888068    0.1129694
```

The remaining uncovered audit case is:

```text
local_id:          15
theta_id:          2
ess_frac:          0.04888068
natural_ess_frac:  1.778591e-12
natural_distance:  5.201461
max_weight:        0.02369464
best_node:         1
nodes:             5
particles:         3000
```

Posterior quantile error versus Stan:

```text
parameter       q05            q50            q95
mu_shape       -0.001046151    0.006468681    0.003528335
mu_scale       -0.007269630   -0.008384840   -0.002074302
mu_shift        0.012263740    0.024696410   -0.003930265
sigma2_shape    0.002357050    0.002355234    0.002353435
sigma2_scale    0.000063830    0.000921614    0.000765555
sigma2_shift    0.000152727    0.000167961   -0.003452949
```

Largest absolute quantile error:

```text
0.02469641
```

`mu_shift` quantiles:

```text
Stan:      q05=-2.089517, q50=-1.678471, q95=-1.416751
Workflow:  q05=-2.077253, q50=-1.653775, q95=-1.420681
```

Phase 0 assessment:

```text
The clean bank-SMC baseline is reproducible and currently stronger than the
discarded residual-corrector experiments.

Accuracy is good on the shifted-gamma posterior plot and quantiles.

The weak point is certification, not visible posterior accuracy: one
subject/theta audit case remains uncovered after repair, with extremely poor
natural-parameter ESS. This is the baseline limitation that the new support
engine must improve or explain.
```

Phase 0 gate:

```text
passed for reproducibility and benchmark recording.
not fully certified because one audit point remains uncovered.
```

## Phase 1 Execution: Local Likelihood Sketch Object

Date: 2026-05-21

Implemented:

```text
/Users/nstevenson/Documents/2025/FastHierarchical/local_likelihood_sketches.R
```

The Phase 1 object represents a theta-independent local likelihood sketch:

```math
\widetilde L_i(\alpha)=c_iN(\alpha\mid a_i,V_i).
```

The implementation fits one local likelihood mode per subject, estimates local
curvature with Nelder-Mead plus BFGS, stores the Laplace constant
approximately, and evaluates direct likelihood profile checks along principal
covariance directions.

Important implementation correction during validation:

```text
Do not discard all observed Hessian curvature when one direction is weak.
Keep valid Hessian directions and borrow fallback precision only for invalid
directions.
```

Without that correction, nearly flat shift directions caused the sketch to use
the broad population covariance for all coordinates, creating artificial profile
errors in shape and scale.

Shifted-gamma Phase 1 diagnostic:

```text
result: /Users/nstevenson/Documents/2025/FastHierarchical/benchmarks/results/shifted_gamma_phase1_local_likelihood_sketches.rds
elapsed_sec: 0.649
result_size: 192.5 Kb
sketch_set_size: 166.5 Kb
subjects: 20
starts_per_subject: 16
start_scale: 4
```

Mode and curvature summary:

```text
hessian_valid: 14 / 20
regularized:   6 / 20
median_condition: 110.334
max_condition:    156123645
mode_loglik:      min=-388.5249, median=-303.2778, max=-216.6185
```

Profile diagnostic:

```text
profile_points:       300
valid_profile_points: 276
invalid_points:       24
median_abs_error:     0.02221915
q90_abs_error:        1.581997
max_abs_error:        127.2475
abs_error_gt_1:       30
abs_error_gt_3:       12
```

Worst valid profile errors:

```text
local_id  direction  radius  error
19        1          2       127.24749
12        1          2        84.65554
3         1          2        70.71889
5         1          2        55.42238
7         1          2        38.27456
```

Interpretation:

```text
Phase 1 implementation is in place, compact, and fast.

The one-Gaussian local likelihood sketch does not pass the shifted-gamma
profile gate globally. It is accurate near the mode for most directions, but
some subjects have bounded/skewed shift geometry where the Gaussian principal
direction walks into regions with much sharper likelihood decay than the sketch
predicts. This is a real approximation limitation, not only a tuning issue.

This validates the planned sequencing: Phase 2 can build the analytic marginal
machinery, but the one-component sketch should be treated as a support map, not
as a final calibrated likelihood representation. If Phase 2 support is too
narrow or too centered, the next surgical improvement is a deliberate second
component or bounded/skew-aware local sketch for the subjects whose profile
check fails, not a generic flexible mixture.
```

Phase 1 gate:

```text
passed for object layer and curvature handling.
not passed as a universally adequate one-Gaussian approximation for shifted
gamma local likelihoods.
```

## Phase 2 Execution: Analytic Tempered Sketch Marginal

Date: 2026-05-21

Implemented in:

```text
/Users/nstevenson/Documents/2025/FastHierarchical/local_likelihood_sketches.R
```

For a one-component sketch

```math
\widetilde L_i(\alpha)=c_iN(\alpha\mid a_i,V_i),
```

the tempered support sketch is evaluated exactly as

```math
\widetilde L_{i,\rho}(\alpha)
=
\widetilde L_i(\alpha)^\rho
=
c_{i,\rho}N(\alpha\mid a_i,V_i/\rho),
```

with the correct Gaussian power constant stored in `c_{i,\rho}`.

For a Gaussian population prior

```math
p(\alpha\mid\theta)=N(\alpha\mid\mu_\theta,\Sigma_\theta),
```

the analytic marginal is

```math
\widetilde m_{i,\rho}(\theta)
=
\int \widetilde L_{i,\rho}(\alpha)p(\alpha\mid\theta)d\alpha
=
c_{i,\rho}N(a_i\mid\mu_\theta,\Sigma_\theta+V_i/\rho).
```

The implementation evaluates this through the population-model reference
component interface, so it is not hard-wired to the shifted-gamma benchmark or
to diagonal covariance parsing in the sketch layer.

Added functions:

```text
local_likelihood_sketch_tempered_components()
local_likelihood_sketch_log_marginal()
local_likelihood_sketch_set_log_marginal_matrix()
local_likelihood_sketch_set_loglik()
local_likelihood_sketch_set_logposterior()
```

Validation:

```text
1D normal-normal analytic marginal checked against direct numerical integration
for rho = 1, 0.5, 0.2.

Saved shifted-gamma Phase 1 sketches evaluated at three theta stress points for
rho = 1, 0.5, 0.2. All subject marginals and aggregate log posteriors were
finite, and aggregate loglik matched row sums of subject marginals.
```

Important scope note:

```text
For one-component sketches, rho tempering is exact.

For future multi-component sketches, the current formula intentionally applies
component-wise power tempering. That is a conservative support-sketch
definition, not the exact power of a Gaussian mixture.
```

## Phase 3 Execution: Tempered Theta Support Engine

Date: 2026-05-21

Implemented:

```text
/Users/nstevenson/Documents/2025/FastHierarchical/benchmarks/run_shifted_gamma_sketch_support_compare.R
```

Supporting framework changes:

```text
build_local_likelihood_sketch_factor_set()
population_factor_set_loglik() dispatch for local likelihood sketch factor sets
optional cess_target in outer_population_smc()
```

The benchmark runs the analytic sketch target through the existing outer SMC:

```math
\widetilde\pi_\rho(\theta)
\propto
p(\theta)\prod_i\widetilde m_{i,\rho}(\theta).
```

For efficiency it does not run each rho independently from the hyperprior.
It runs SMC once at a broad target and then moves through the rho ladder by
factor updates:

```text
rho ladder: 0.05, 0.10, 0.20, 0.35, 0.50, 0.75, 1.00
outer particles: 1000
```

Artifacts:

```text
results: /Users/nstevenson/Documents/2025/FastHierarchical/benchmarks/results/shifted_gamma_sketch_support_results.rds
plot:    /Users/nstevenson/Documents/2025/FastHierarchical/benchmarks/results/shifted_gamma_sketch_support_rho_ladder.png
```

Timing:

```text
elapsed_sec: 17.8

rho   elapsed_sec  rounds  update_ess  accept_rate
0.05  11.495       18      NA          0.2625
0.10   0.888        0      0.4997999   0.2450
0.20   0.850        0      0.7785876   0.2625
0.35   0.850        0      0.4641690   0.2435
0.50   0.852        0      0.8434598   0.2485
0.75   0.858        0      0.4716356   0.2520
1.00   0.859        0      0.8471942   0.2290
```

Support recall against Stan:

```text
parameter       stan_q01       stan_q99       support_q005   support_q995   covers
mu_shape         1.301524       1.576896       1.215717       1.718927      yes
mu_scale        -0.623656      -0.391371      -0.764786      -0.305563      yes
mu_shift        -2.333599      -1.338451      -2.966177      -0.488483      yes
sigma2_shape     0.016634       0.068739       0.014834       0.118959      yes
sigma2_scale     0.014599       0.063044       0.011257       0.081133      yes
sigma2_shift     0.007077       0.062934      18.101719      82.422869      no
```

Critical diagnosis:

```text
Phase 3 does not pass the support-recall gate.

The failure is not outer SMC. The rho ladder sampled the analytic sketch target
efficiently.

The failure is that raw local likelihood sketches are the wrong support source
for weakly identified bounded shift effects. Several subject-level likelihood
modes are at eta_shift around -14 to -16, while the real hierarchical posterior
has mu_shift around -1.68 and sigma2_shift around 0.017. The analytic support
engine then explains those incompatible raw local modes by inflating
sigma2_shift to roughly 20-80.
```

Concrete evidence:

```text
local eta_shift sketch modes:
1=-14.35, 2=-15.70, 9=-13.37, 11=-15.39, 14=-14.29, 15=-16.13, 20=-14.88

Stan sigma2_shift quantiles:
q01=0.007077, q50=0.016672, q99=0.062934

Sketch support sigma2_shift quantiles by rho:
rho=0.05: q01=22.24, q50=37.16, q99=72.70
rho=0.20: q01=19.72, q50=37.22, q99=69.89
rho=1.00: q01=19.65, q50=37.96, q99=72.09
```

Median theta comparison under the sketch posterior:

```text
rho   logpost(Stan median theta)  logpost(sketch-supported theta)  difference
0.05  -1649.089                  -415.623                         -1233.466
0.20  -5966.153                  -1365.822                        -4600.331
1.00  -21122.960                 -6314.474                        -14808.480

Stan median shift:    mu=-1.678471, sigma2=0.016672
Sketch support shift: mu=-1.789728, sigma2=38.56552
```

Phase 3 gate:

```text
passed for implementation and efficient rho-ladder execution.
failed for support recall.
```

Consequence:

```text
The next architectural correction is not more rho tempering and not random
inflation. The top-level support map must not be built from raw local likelihood
modes alone.

The local sketch should be a weakly regularized or population-aware likelihood
site, not the unconstrained MLE/Laplace approximation to L_i(alpha). The support
engine needs local sketches of L_i(alpha) times a deliberately broad but finite
reference prior, so non-identifiable local ridges do not generate impossible
population variance support.
```

## Reference-Prior Site Replacement

Date: 2026-05-21

Implemented in:

```text
/Users/nstevenson/Documents/2025/FastHierarchical/local_likelihood_sketches.R
```

The raw likelihood-MLE sketch was replaced by a reference-prior likelihood site.
For each subject, fit

```math
q_i^0(\alpha)
\propto
L_i(\alpha)r_0(\alpha),
```

where `r0` is a broad finite Gaussian reference prior derived from the
population-model reference theta. The likelihood sketch is then represented by
dividing out the reference prior in natural-parameter form:

```math
\log \widetilde L_i(\alpha)
=
c_i + h_i^\top \alpha - \frac12\alpha^\top P_i\alpha.
```

This is the important change: flat local likelihood directions become near-zero
site precision instead of arbitrary precise pseudo-observations at the raw MLE.

The analytic marginal is now evaluated directly from the quadratic site:

```math
\widetilde m_{i,\rho}(\theta)
=
\int
\exp\{\rho(c_i+h_i^\top\alpha-\frac12\alpha^\top P_i\alpha)\}
N(\alpha\mid\mu_\theta,\Sigma_\theta)
d\alpha.
```

Diagnostic effect on shifted gamma:

```text
old raw eta_shift sketch centers included:
-14.35, -15.70, -13.37, -15.39, -14.29, -16.13, -14.88

new reference-site eta_shift centers:
min=-1.894, median=-1.495, max=-0.479

old support sigma2_shift q005-q995:
18.10 to 82.42

new support sigma2_shift q005-q995:
0.00648 to 0.12609
```

The first replacement run fixed the impossible variance but still undercovered
the left `mu_shift` tail. That was not the old boundary-MLE bug; it was a support
ladder problem. With 20 subjects, `rho=0.05` is already one subject-equivalent
likelihood. The support ladder must start substantially below `1 / n_subjects`
when the goal is high recall rather than posterior accuracy.

Final recorded support ladder:

```text
rho: 0.0005, 0.001, 0.005, 0.01, 0.02, 0.05,
     0.10, 0.20, 0.35, 0.50, 0.75, 1.00
```

Final support benchmark:

```text
results: /Users/nstevenson/Documents/2025/FastHierarchical/benchmarks/results/shifted_gamma_sketch_support_results.rds
plot:    /Users/nstevenson/Documents/2025/FastHierarchical/benchmarks/results/shifted_gamma_sketch_support_rho_ladder.png
elapsed_sec: 13.1
```

Support recall after replacement:

```text
parameter       stan_q01       stan_q99       support_q005   support_q995   covers
mu_shape         1.301524       1.576896       0.704191       2.165247      yes
mu_scale        -0.623656      -0.391371      -1.161467       0.221223      yes
mu_shift        -2.333599      -1.338451      -2.366883      -0.436539      yes
sigma2_shape     0.016634       0.068739       0.013003       0.107548      yes
sigma2_scale     0.014599       0.063044       0.009763       0.071786      yes
sigma2_shift     0.007077       0.062934       0.006480       0.126087      yes
```

Assessment:

```text
This fixes the diagnosed Phase 3 failure mode. The support engine is now a
conservative theta-support generator, not a posterior approximation.

The rho=1 sketch posterior is still not expected to match Stan. That is not the
target of this phase. The correct next use is to select local-SMC correction
anchors from the recorded support envelope.
```

## Phase 4 Execution: Pathfinder-Style Support Proposals

Date: 2026-05-21

Implemented:

```text
/Users/nstevenson/Documents/2025/FastHierarchical/theta_support_pathfinder.R
```

The implementation is deliberately narrower than full Pathfinder. It uses the
parts needed for theta-support discovery:

```text
rho-SMC support cloud
-> dispersed starts in whitened theta coordinates
-> multi-start BFGS on the cheap sketch log posterior
-> final Hessian-based Gaussian approximations
-> probability-contour stress candidates
-> empirical novelty score against rho-SMC support spacing
```

The grounding is:

```text
Pathfinder: multi-start quasi-Newton paths plus inverse-Hessian Gaussian
approximations. See Zhang, Carpenter, Gelman, Vehtari, JMLR 2022.

Gaussian contour candidates: ellipsoids defined by chi-square probability
levels under the local Gaussian approximation.

Novelty: distance in whitened theta coordinates compared to the empirical
nearest-neighbor spacing of the existing rho-SMC support cloud. This is a
support-design diagnostic, not a posterior-correctness claim.
```

Important cleanup during implementation:

```text
Discarded arbitrary candidate-producing rules:
- no fixed absolute "1 SD" novelty threshold;
- no arbitrary gradient-norm cutoff for accepting paths;
- no arbitrary log-posterior drop cutoff for accepting paths;
- no hand inflation factors such as 1, 2, 3 as the primary contour rule.
```

Current acceptance rule for a path is intentionally strict:

```text
finite log posterior
optimizer convergence code 0
positive-definite final Hessian
```

Bad optimizer endpoints are retained in `path_summary` for diagnosis but do not
create proposal components.

Candidate contours are:

```text
coverage probabilities: 0.50, 0.90, 0.99
radius: sqrt(qchisq(probability, df = dim(theta)))
```

Benchmark command:

```text
Rscript benchmarks/run_shifted_gamma_sketch_support_compare.R
```

Artifacts:

```text
results: /Users/nstevenson/Documents/2025/FastHierarchical/benchmarks/results/shifted_gamma_sketch_support_results.rds
plot:    /Users/nstevenson/Documents/2025/FastHierarchical/benchmarks/results/shifted_gamma_sketch_support_rho_ladder.png
```

Timing:

```text
elapsed_sec: 27.3
```

Support recall stayed unchanged from Phase 3:

```text
All Stan q01-q99 intervals are covered by the rho-ladder support q005-q995.
```

Pathfinder result:

```text
starts:                    24
finite optimizer paths:    24
converged paths:           23
unique endpoint clusters:  2
usable paths:              1
components:                3
candidates:                39
nonduplicate candidates:   1
support distance cutoff:   1.288083
max candidate distance:    1.495158
```

The single nonduplicate candidate is a 99% contour point on the only usable
Hessian component:

```text
candidate_id: 29
point:        pc1_minus
mu_shift:     -0.8600031
log_sigma2_shift: -5.886551
nearest_support_distance: 1.495158
relative sketch log posterior: -19.00749
```

Assessment:

```text
Phase 4 is now a conservative support-proposal diagnostic, not a way to force
extra anchors.

For this shifted-gamma sketch target, Pathfinder mostly confirms that the
rho-SMC ladder already covers the dominant approximate posterior basin. It adds
one novel 99% Gaussian-contour stress candidate, but it does not discover a new
high-posterior mode or ridge.

This is useful negative evidence. The next accuracy work should not expect
Pathfinder alone to fix local evidence calibration. Phase 5 should evaluate
expensive local-SMC corrections at rho-SMC tail/center points plus the small
Pathfinder stress set.
```

## Phase 5 Execution: Expensive Local-SMC Correction Points

Date: 2026-05-21

Implemented:

```text
/Users/nstevenson/Documents/2025/FastHierarchical/theta_smc_corrections.R
```

The Phase 5 object is a theta-level correction observation set:

```text
theta point
-> exact local SMC log marginal by subject
-> analytic sketch log marginal by subject
-> total correction Delta(theta)
-> MCSE and degeneracy diagnostics
```

This deliberately does not create bank nodes and does not store local particles.
The local evidence evaluations reuse the existing refined local SMC machinery
through `run_tempered_smc()` with endpoint population priors built from each
theta point.

Correction design:

```text
rho-SMC center
rho-SMC marginal 5% and 95% points for each theta coordinate
nonduplicate Pathfinder stress candidates
```

For the shifted-gamma run this produced:

```text
points:             14
point evaluations:  16
local SMC runs:     320
particles/run:      600
correction time:    26.569 sec
total benchmark:    51.6 sec
```

Replicated points:

```text
center
pathfinder_9_pc1_minus_p99
```

Artifacts:

```text
results: /Users/nstevenson/Documents/2025/FastHierarchical/benchmarks/results/shifted_gamma_sketch_support_results.rds
plot:    /Users/nstevenson/Documents/2025/FastHierarchical/benchmarks/results/shifted_gamma_sketch_support_rho_ladder.png
```

Main correction observations:

```text
label                         Delta(theta)        MCSE      rel logpost SMC
center rep1                   -3.317794           0.238     -3.559676
center rep2                   -2.507367           0.238     -2.749249
mu_eta_shape_q5              -31.05046            1.145    -36.11929
mu_eta_shape_q95              -4.274400           0.247     -8.841911
mu_eta_scale_q5              -14.58428            0.295    -21.48750
mu_eta_scale_q95              -1.232299           0.235     -5.084044
mu_eta_shift_q5                1.899180           0.226      0.000000
mu_eta_shift_q95             -23.43075            0.482    -26.51119
log_sigma2_eta_shape_q5       -9.477545           0.271    -16.08492
log_sigma2_eta_shape_q95      -3.716944           0.242     -6.809885
log_sigma2_eta_scale_q5        2.585814           0.225     -0.206009
log_sigma2_eta_scale_q95      -2.837206           0.244     -8.004267
log_sigma2_eta_shift_q5        0.747557           0.232     -2.859212
log_sigma2_eta_shift_q95      -4.403122           0.249     -9.190739
pathfinder stress rep1        -5.0e12             0.245     -5.0e12
pathfinder stress rep2        -4.0e12             1.029     -4.0e12
```

The Pathfinder stress point is not a normal calibration point. It is exact-local
support failure:

```text
pathfinder_9_pc1_minus_p99:
  degenerate local subjects: 5 in replicate 1, 4 in replicate 2
  sketch relative log posterior: -18.72
  exact local-SMC relative log posterior: effectively impossible
```

This is useful because it shows the cheap sketch/Pathfinder approximation can
place a stress contour outside the feasible local likelihood support. Those
points should be used as rejection/certification evidence, not as ordinary
surrogate training data.

Assessment:

```text
Phase 5 is implemented and produces reusable correction observations.

The correction surface is not close to constant. Some posterior-relevant
rho-SMC tail points change by several log units, and the exact local SMC ranking
differs strongly from the sketch ranking:

  sketch center is best;
  exact local SMC prefers the mu_shift lower-tail point in this design.

This means Phase 6 cannot honestly use a constant correction. It needs a
theta-level correction surrogate fitted and validated against these local-SMC
observations, while excluding exact-impossible stress points from smooth
calibration fits.
```

## Phase 6 Execution: Theta-Level Correction Surrogate

Date: 2026-05-21

Implemented:

```text
/Users/nstevenson/Documents/2025/FastHierarchical/theta_correction_surrogate.R
```

The Phase 6 object fits

```math
\Delta(\theta)
=
\ell(\theta)-\widetilde\ell(\theta)
```

from Phase 5 exact local-SMC correction observations.

Training construction:

```text
aggregate replicated correction observations by theta point
exclude exact-impossible / degenerate stress points from smooth calibration
whiten theta using the conservative support cloud
fit constant, linear, and quadratic ridge regressions in whitened theta
select by leave-one-point-out MCSE-aware predictive score
```

Observation noise for a correction point is estimated from local-SMC MCSE and
replicate variation:

```math
s_j
=
\max\left(
  \sqrt{\sum_r \mathrm{MCSE}_{jr}^2}/R_j,\,
  \mathrm{sd}(\Delta_{jr})/\sqrt{R_j}
\right).
```

The regression loss is weighted by this `s_j`. The selected model minimizes the
leave-one-out Gaussian negative predictive score.

Benchmark command:

```text
Rscript benchmarks/run_shifted_gamma_sketch_support_compare.R
```

Artifacts:

```text
results: /Users/nstevenson/Documents/2025/FastHierarchical/benchmarks/results/shifted_gamma_sketch_support_results.rds
plot:    /Users/nstevenson/Documents/2025/FastHierarchical/benchmarks/results/shifted_gamma_sketch_support_rho_ladder.png
```

Run summary:

```text
total elapsed_sec:    52.7
training points:      13
rejected points:       1
selected model:        quadratic ridge
selected lambda:       31.62278
LOO NLP:              121.5703
LOO RMSE:               5.954904
LOO MAE:                4.508979
LOO weighted RMSE:     15.66939
LOO max abs error:     15.17210
saved object size:     12.037 MB
```

Model comparison:

```text
constant correction:
  best LOO NLP  = 259.4102
  best LOO RMSE = 10.6960

linear correction:
  best LOO NLP  = 144.3156
  best LOO RMSE = 5.8946

quadratic ridge correction:
  best LOO NLP  = 121.5703
  best LOO RMSE = 5.9549
```

The selected model is quadratic by MCSE-aware predictive score, not by visual
posterior appearance. The best unweighted RMSE is a near-interpolating quadratic
fit, but the MCSE-aware score rejects that in favor of stronger ridge
regularization.

Rejected from smooth calibration:

```text
pathfinder_9_pc1_minus_p99
reason: exact-local support failure / degenerate subjects
```

Important validation residuals:

```text
label                       Delta     LOO pred   LOO residual
center                     -2.913     -2.789     -0.124
mu_eta_shape_q5           -31.050    -21.409     -9.642
mu_eta_scale_q5           -14.584    -19.292      4.708
mu_eta_shift_q5             1.899     -0.544      2.443
mu_eta_shift_q95          -23.431     -8.259    -15.172
log_sigma2_eta_scale_q5     2.586     -0.975      3.561
```

Assessment:

```text
Phase 6 is implemented, but the current correction design is not yet strong
enough to justify Phase 7 as a final corrected posterior.

The result is useful because it proves:

1. constant correction is decisively wrong;
2. linear correction is better but still weak;
3. quadratic structure helps, but leave-one-out errors remain several log units
   in posterior-relevant directions.

## Phase 6b: Active Correction Refinement

The correction surrogate now supports active refinement. The refinement rule is
deliberately narrow:

```text
fit correction surrogate
compute leave-one-out residuals at correction points
weight residuals by exact-local-SMC posterior relevance
add midpoint correction points between bad residual points and:
  1. the center point;
  2. the current best exact-local-SMC point
rerun exact local SMC at those points
refit the surrogate
repeat
```

This is not a new posterior approximation. It is a validation-driven way to
add expensive correction points only where the current correction surface is
provably weak.

One implementation detail matters: after combining correction batches, relative
log posterior values are recomputed globally. They cannot be inherited from each
small batch, because batch-local relevance would corrupt later refinement
rounds.

Benchmark command:

```text
Rscript benchmarks/run_shifted_gamma_sketch_support_compare.R
```

Artifacts:

```text
results: /Users/nstevenson/Documents/2025/FastHierarchical/benchmarks/results/shifted_gamma_sketch_support_results.rds
plot:    /Users/nstevenson/Documents/2025/FastHierarchical/benchmarks/results/shifted_gamma_sketch_support_rho_ladder.png
```

Run summary after two refinement rounds:

```text
total elapsed_sec:          79.5
initial correction points:  14
refinement points:           8
final correction points:    22
point evaluations:          32
local SMC runs:            640
particles per local run:   600
saved object size:          13.8 MB
```

Surrogate validation improved substantially:

```text
                           initial      refined
training points                 13           21
rejected points                  1            1
selected degree                  2            2
selected lambda           31.62278          0.1
LOO NLP                   121.5703      53.6660
LOO RMSE                    5.9549       2.9302
LOO MAE                     4.5090       2.0962
LOO weighted RMSE          15.6694      10.4749
LOO max abs error          15.1721       8.0303
```

The selected refinement points were:

```text
refine_log_sigma2_eta_scale_q5_to_center
refine_log_sigma2_eta_scale_q5_to_best
refine_mu_eta_shift_q5_to_center
refine_log_sigma2_eta_scale_q95_to_center
refine_mu_eta_scale_q95_to_center
refine_mu_eta_scale_q95_to_best
refine_log_sigma2_eta_shift_q5_to_center
refine_log_sigma2_eta_shift_q5_to_best
```

Important remaining validation misses:

```text
label                                 Delta   LOO pred   residual  rel logpost
log_sigma2_eta_shift_q5               0.748    -4.046      4.794       -3.076
log_sigma2_eta_scale_q5               2.586     4.355     -1.769       -0.423
refine_log_sigma2_eta_shift_q5_best   1.536    -0.240      1.775       -0.564
```

Assessment:

```text
Active refinement is doing the right kind of work: it reduced large smooth
surrogate errors without adding arbitrary global particles or anchors.

It is still not sufficient evidence for final corrected posterior sampling.
The correction surface has a posterior-relevant miss in the eta_shift variance
lower-tail direction. That is the next point to resolve before Phase 7.
```
```
