# Uncertainty Propagation Plan

## Problem

The framework needs the population posterior and total model evidence

```text
p(theta | y) proportional to p(theta) prod_i m_i(theta)
m_i(theta) = int L_i(alpha_i) p(alpha_i | theta) d alpha_i
```

The current chart-atlas machinery approximates each `log m_i(theta)` with a
frozen deterministic local evidence surface. EMC reference diagnostics show
that support is not the current blocker: reference theta points are certified.
The blocker is coherent numerical error in the product evidence surface:

```text
e_i(theta) = log mhat_i(theta) - log m_i(theta)
E(theta) = sum_i e_i(theta)
```

A constant `E(theta)` only shifts total evidence. A varying `E(theta)` reshapes
the posterior. The observed centered `E(theta)` range over EMC reference
profiles is about 25 log units, so the surface error is posterior-shaping.

## Required Change

The atlas must stop pretending that its local evidence surface is exact. It is
a proposal/surrogate layer. The framework needs an uncertainty-aware evidence
layer that can:

1. estimate local evidence error at posterior-relevant theta points;
2. model local error as a function of theta;
3. propagate uncertainty into outer posterior and evidence;
4. optionally produce a delayed-acceptance or pseudo-marginal correction for
   final inference.

Do not add more generic anchors as the main fix. More anchors address support
and interpolation gaps. The current failure is calibrated-but-biased local
evidence inside the certified region.

## Architecture

```text
theta design/proposal
  -> local chart atlas surrogate
  -> reference-theta calibration design
  -> replicated local SMC evidence probes
  -> local error model e_i(theta)
  -> corrected stochastic factor set
  -> outer SMC with numerical uncertainty
  -> delayed-acceptance final correction when required
```

The atlas remains useful because it cheaply proposes and evaluates plausible
theta regions. It is not final evidence until the local surface error is either
bounded or corrected.

## Phase 1: Evidence Audit Dataset

Create a first-class object:

```text
local_evidence_audit
```

It stores rows:

```text
local_id
theta_id
theta
atlas_log_m
atlas_se
atlas_status
fresh_log_m_replicates
fresh_log_m_center
fresh_log_m_sd
fresh_log_m_se
reference_source
error_fresh_minus_atlas
```

Rules:

- `atlas_status != certified` is a support/certification failure.
- `fresh_log_m_sd` is empirical replicate uncertainty and dominates reported
  path MCSE when larger.
- A single fresh SMC endpoint is screening information, not ground truth.
- Reference posterior draws such as EMC/Stan may be used only for benchmarks,
  never for production calibration design.

Implementation owner:

```text
local_charts.R
```

Expected functions:

```text
build_local_evidence_audit()
summarize_local_evidence_audit()
select_audit_failures()
extend_local_evidence_audit()
```

## Phase 2: Posterior-Relevant Calibration Design

Calibration theta points must be chosen from the workflow itself:

```text
outer posterior profiles
inflated posterior tails
high surface-uncertainty directions
high posterior leverage points
```

Benchmark-only diagnostics may additionally use reference theta points from
EMC/Stan to reveal error, but those points must not train the production
correction.

Selection score should include:

```text
posterior_weight(theta)
predicted numerical uncertainty
local sensitivity |d log m_i / d theta|
ensemble/member disagreement
fresh-probe screening error
```

The selected design must include profiles along variance dimensions. EMC
diagnostics show large tilt in:

```text
log_sigma2_v_LogFreq
log_sigma2_v
log_sigma2_sv
mu_sv
```

Production entry point:

```text
build_local_evidence_calibration_design()
```

It accepts a workflow posterior cloud or `fit`, optional frozen atlas
`factor_sets`, and produces a normalized `local_evidence_theta_design`.
Reference posterior draws may be passed only by benchmark scripts and should be
marked through `theta_source`.

## Phase 3: Replicated Local SMC Probes

For selected `(local, theta)` pairs, run independent local SMC replicates.

Default policy:

```text
screening: cheap, broad, one replicate
confirmation: 2-4 independent replicates for selected high-impact pairs
promotion/correction: use only confirmed estimates
```

The empirical uncertainty is:

```text
s_i(theta)^2 = max(
  var(log fresh_m_replicates),
  mean(path_mcse^2) / R,
  floor^2
)
```

Do not treat SMC-reported MCSE alone as reliable. The EMC diagnostics showed
replicate SDs up to about 1.5-1.8 log units on hard local/theta pairs.

Production entry points:

```text
run_local_evidence_replicates()
build_local_evidence_audit()
extend_local_evidence_audit()
```

`build_local_evidence_audit()` runs broad screening probes. `select_audit_failures()`
chooses posterior-weighted local/theta failures. `extend_local_evidence_audit()`
adds independent confirmation replicates for those pairs and recomputes the
audit table without changing the frozen atlas surface.

## Phase 4: Local Error Model

Fit the error:

```text
delta_i(theta) = log m_i(theta) - log mhat_i(theta)
```

Initial model:

```text
delta_i(theta) = b_i^T phi(theta) + epsilon_i(theta)
epsilon_i(theta) ~ N(0, s_i(theta)^2)
```

Use a small, explicit basis in whitened theta:

```text
1
linear terms in active hyperparameters
selected quadratic terms
interaction terms only when diagnostics justify them
```

Regularization is mandatory. Hard locals should shrink toward zero correction
unless the replicated evidence supports a stable trend.

A later version can replace this with a sparse GP or Bayesian MBAR-style model,
but the first implementation should be small and auditable.

Expected functions:

```text
fit_local_evidence_error_model()
predict_local_evidence_error()
```

Production entry points:

```text
fit_local_evidence_error_model()
predict_local_evidence_error()
```

By default, the fitter uses only confirmed replicated probes
(`n_replicates >= 2` and `n_ok >= 2`). Screening rows are retained in the audit
but are not promoted into the correction surface unless explicitly requested.
The fitted object carries local diagnostics, theta-level calibration residuals,
and posterior-weighted total-surface summaries.

Acceptance diagnostics:

```text
leave-one-theta-out error
leave-one-local-out contribution
posterior-weighted RMSE
centered total E(theta) range
calibration residual z-scores
```

## Phase 5: Corrected Factor Set

Create a factor set wrapper:

```text
local_evidence_corrected_factor_set
```

For each theta:

```text
log m_corrected_i(theta)
  = log mhat_i(theta) + E[delta_i(theta)]
```

and expose uncertainty:

```text
Var_delta_i(theta)
Cov_delta_i(theta, theta')
```

The deterministic corrected factor set may be used for proposal generation and
diagnostics. It is not final evidence unless uncertainty is negligible.

Expected functions:

```text
build_corrected_local_atlas_factor_set()
corrected_factor_set_loglik()
corrected_factor_set_uncertainty()
```

Production entry points:

```text
build_corrected_local_atlas_factor_set()
corrected_factor_set_loglik()
corrected_factor_set_loglik_by_local()
corrected_factor_set_uncertainty()
```

The corrected factor set is a standard `population_factor_set` and can be passed
to `population_factor_set_loglik()` and `outer_population_smc()` directly. It
adds the fitted mean correction to the frozen atlas surface and keeps numerical
uncertainty separate for phase 6 accounting.

## Phase 6: Propagate Numerical Uncertainty

Outer SMC must report two things:

```text
posterior under corrected mean surface
numerical uncertainty sensitivity
```

Use posterior draws to compute:

```text
S(theta) = sum_i Var_delta_i(theta)
```

and report:

```text
posterior mean uncertainty contribution
log evidence numerical SE
worst posterior-shaping direction
ranked local contributors
```

If uncertainty is not negligible, the result is not final evidence.

Production entry points:

```text
summarize_corrected_outer_uncertainty()
run_corrected_outer_smc_with_uncertainty()
```

`outer_population_smc()` also attaches `numerical_uncertainty` automatically
when the supplied factor set inherits from `local_evidence_corrected_factor_set`.
The phase 6 object reports evidence SE, posterior log-factor uncertainty,
parameter-mean sensitivity, and ranked local contributors. The uncertainty
calculation uses quadratic forms of the fitted local error model rather than
materializing per-local dense covariance matrices.

Acceptance target:

```text
centered total error range over audit theta < 2 log units
posterior-weighted local evidence RMSE materially below atlas-only
frozen outer reruns stable
reference benchmark posterior improves without using reference draws for training
```

## Phase 7: Delayed-Acceptance Correction

For final posterior/evidence, use the corrected atlas as a cheap first stage.

Stage 1:

```text
theta proposal accepted/rejected under corrected surrogate
```

Stage 2:

```text
only promising theta values get replicated local evidence correction
```

If an unbiased nonnegative estimator of `prod_i m_i(theta)` is available with
manageable variance, use a pseudo-marginal correction. Otherwise use delayed
acceptance as a bias-reduction and certification layer, and report the remaining
numerical error explicitly.

Do not run full fresh local SMC for every outer particle. That gives the right
target only in principle and is not scalable.

Expected functions:

```text
outer_population_delayed_acceptance()
estimate_theta_correction()
cache_local_evidence_corrections()
```

Production entry points:

```text
outer_population_delayed_acceptance()
cache_local_evidence_corrections()
estimate_theta_correction()
```

The delayed-acceptance implementation treats the corrected outer posterior as
stage 1, selects posterior theta points by weighted resampling by default, runs
replicated fresh local evidence probes only at those theta points, and estimates
the stage-2 evidence/posterior correction by importance reweighting. Subset
local correction is allowed for diagnostics, but the evidence account marks it
as partial and will not label the Bayes factor precise.

## Phase 8: Evidence Accounting

Total evidence should be decomposed as:

```text
log Z_total
  = log Z_outer_corrected_mean
  + local numerical correction
  +/- numerical uncertainty
```

Report:

```text
log evidence estimate
outer Monte Carlo SE
local evidence numerical SE
posterior-shaping error diagnostics
number of replicated probes
hard local contributors
```

If the numerical uncertainty dominates the model-comparison scale, do not report
a precise Bayes factor.

Production entry point:

```text
build_local_evidence_accounting()
```

The accounting object combines the corrected mean-surface outer evidence, the
stage-2 local numerical correction, outer SMC MCSE, phase-6 surrogate numerical
SE, correction SE, posterior-shaping diagnostics, replicated probe counts, and
hard-local contributors.

## Phase 9: Benchmark Acceptance

Use EMC/Stan reference only as an external benchmark.

A benchmark run must include:

1. atlas-only posterior;
2. corrected posterior;
3. validation-weighted ensemble posterior if multiple atlas members are used;
4. reference-theta diagnostic;
5. held-out local evidence probes not used for fitting correction;
6. posterior comparison against EMC/Stan.

Acceptance criteria:

```text
corrected posterior improves over atlas-only on max and mean standardized error
corrected posterior improves shape error
held-out local evidence RMSE improves
centered total log surface error range shrinks materially
fresh endpoint failures do not concentrate in a small hard-local set
```

The validation-weighted ensemble is allowed as a diagnostic of numerical bias,
but it is not the final architecture. The final architecture must model or
correct the bias that the ensemble currently cancels accidentally.

## Non-Goals

- Do not use reference posterior draws to train production correction.
- Do not promote single fresh SMC endpoints as hard truth.
- Do not add generic anchors as the main repair.
- Do not hide numerical uncertainty by averaging posterior draws.
- Do not report final model evidence without a local evidence uncertainty term.
- Do not preserve old q0-predictive local machinery as a parallel solution.

## Minimal Implementation Order

1. Implement `local_evidence_audit`.
2. Implement replicated probe selection and caching.
3. Fit a regularized local error model.
4. Build `local_evidence_corrected_factor_set`.
5. Run outer SMC on corrected mean surface.
6. Add held-out audit diagnostics.
7. Add delayed-acceptance correction for final posterior/evidence.

The first useful milestone is not a perfect posterior. It is a run where the
held-out centered total log surface error range is much smaller than the current
25-log-unit range and the posterior improves without reference-trained tuning.
