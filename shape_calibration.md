# Shape Calibration

This document defines the production layer for detecting and fixing
certified-but-biased local evidence surfaces.

It is deliberately separate from raw support certification.

Raw certification asks:

```text
Can the local atlas evaluate ell_i(theta) with adequate PMIS support?
```

Shape calibration asks:

```text
When the atlas says it can evaluate ell_i(theta), is the value accurate enough
to preserve the posterior and model evidence shape?
```

Those are different problems. The current EMC residuals show why this layer is
needed: raw PMIS support is now mostly present, compression is not the main
failure, and outer SMC appears stable, but `mu_sv`, `sigma2_sv`,
`mu_v_LogFreq`, and `sigma2_v_LogFreq` still show coherent posterior shape
error. That points to small but systematic errors in local evidence curvature.

## Target

For each local `i`, the object of interest is:

```text
ell_i(theta) = log m_i(theta)
             = log integral L_i(alpha_i) p(alpha_i | theta) d alpha_i
```

The framework already builds a local atlas estimator:

```text
ell_hat_i(theta)
```

Shape calibration estimates and controls:

```text
delta_i(theta) = ell_i(theta) - ell_hat_i(theta)
```

only over posterior- or proposal-relevant theta regions.

The goal is not to make every local exact everywhere. The goal is:

```text
posterior-weighted local evidence residuals must be small enough that the
product evidence surface and final posterior are stable.
```

## Non-Goals

Do not run nested SMC for every local/theta pair.

Do not hardcode named axes such as `sv` or `v_LogFreq` into the algorithm.
Axis labels are allowed for diagnostics only.

Do not use residual correction as the first repair mechanism. First try to make
the atlas itself explain the probes through better chart placement and graph
calibration.

Do not let compression participate in shape certification. Shape probes compare
against the raw local atlas. Compression remains a downstream acceleration
layer that must be certified separately against the raw surface.

## Inputs

### Factor Set

A frozen local atlas factor set after support certification and optional
compression. Shape calibration must evaluate raw PMIS internally, even if the
outer run uses compressed PMIS.

### Theta Cloud

A weighted theta cloud with provenance:

```text
theta_id
theta
theta_weight
theta_source
```

Allowed production sources:

```text
initial q0 proposal samples
theta design anchors
inflated proposal tails
pilot outer posterior samples
pilot outer posterior tails
automatically selected high-leverage profiles
```

Benchmark-only sources such as EMC reference draws may be used for diagnosis,
but never for production repair selection.

### Raw Certification Table

One row per local/theta query:

```text
local
theta_id
theta_weight
raw_log_m
raw_se
status
reason
particle_mis_ess
particle_mis_ess_frac
particle_mis_psis_k
min_covering_distance
nearest_charts
local_score_norm
local_curvature_norm
impact_score
```

This table is used to rank candidate probes. Unlike support repair, candidates
must include certified rows.

## Probe Selection

Shape probes should be expensive and sparse. They are selected by posterior
impact and diagnostic risk, not by named dimensions.

For a local/theta pair `(i, j)`, define a probe priority:

```text
P_ij =
  W_j
  * A_ij
  * U_ij
  * G_i
```

where:

```text
W_j = normalized theta weight
```

`A_ij` is an atlas-sensitivity score:

```text
A_ij =
  1
  + c1 * log1p(local_score_norm_ij)
  + c2 * log1p(local_curvature_norm_ij)
  + c3 * log1p(min_covering_distance_ij)
```

`U_ij` is a soft instability score:

```text
U_ij =
  1
  + c4 * max(0, min_ess_frac - ess_frac_ij) / min_ess_frac
  + c5 * max(0, psis_k_ij - psis_threshold)
  + c6 * raw_se_ij
```

`G_i` is local graph/compression risk:

```text
G_i =
  1
  + c7 * graph_edge_residual_i
  + c8 * graph_logZ_se_i
  + c9 * raw_compression_holdout_rmse_i
```

The selected set must include:

1. high-priority certified pairs;
2. a small number of uncertified high-impact pairs, if any remain;
3. weighted posterior tails found automatically from the theta cloud;
4. holdout pairs not used for repair.

The split is mandatory:

```text
repair probes != holdout probes
```

## Direct Probe Estimator

At each selected pair `(i, theta_j)`, run `R` independent direct local SMC
paths targeting:

```text
pi_ij(alpha) proportional to L_i(alpha) p(alpha | theta_j)
```

Each path returns:

```text
z_hat_r = log Z_r
s_r     = path standard error
```

Combine replicates on the evidence scale:

```text
z_probe = log mean_r exp(z_hat_r)
```

Do not average log normalizers directly.

Estimate probe uncertainty from both:

```text
replicate variability
path standard errors
```

A practical initial uncertainty is:

```text
s_probe^2 =
  var_r(z_hat_r)
  + mean_r(s_r^2)
```

with bootstrap on `exp(z_hat_r)` when `R >= 3`.

## Residual Table

For each probe:

```text
atlas_z     = raw atlas PMIS estimate
atlas_se    = raw atlas PMIS standard error
probe_z     = direct SMC probe estimate
probe_se    = direct SMC probe standard error
residual    = probe_z - atlas_z
residual_se = sqrt(probe_se^2 + atlas_se^2)
z_residual  = residual / residual_se
```

The table must also carry theta, local id, certification diagnostics, nearest
charts, and posterior weight.

## Residual Geometry

Shape calibration learns geometry from residuals, not from named axes.

For each local and for the total evidence residual, whiten theta:

```text
x_j = Sigma_theta^{-1/2} (theta_j - center_theta)
```

The shape error that matters for posterior distortion is not the absolute
offset of one local. Additive local offsets mostly shift total log evidence.
Posterior shape error comes from variation of the residual field over theta:

```text
delta_i(theta) = ell_i(theta) - ell_hat_i(theta)
```

So Phase 4 estimates the active directions of the residual field. The ideal
object is the residual active-subspace matrix:

```text
A = sum_i E[ grad_x delta_i(x) grad_x delta_i(x)^T ]
    + E[ grad_x Delta(x) grad_x Delta(x)^T ]
```

where:

```text
Delta(theta) = sum_i delta_i(theta)
```

and `x` is the whitened theta coordinate. We do not have analytic residual
gradients, so the production estimator fits robust ridge local linear models:

```text
delta_ij - offset_i = a_i + g_i^T x_j + noise_ij
```

with uncertainty:

```text
noise_ij variance approximately residual_se_ij^2 + tau_local^2
```

and Student-t style robust weights. The total residual field is fitted
similarly:

```text
centered_Delta_j = a_T + g_T^T x_j + noise_j
```

using the observed local coverage fraction as a weight. The estimated active
matrix is:

```text
A_hat =
  c_grad * ( sum_i rho_i g_i g_i^T + rho_T g_T g_T^T )
  + c_loc  * ( sum_ij e_ij x_j x_j^T + sum_j e_Tj x_j x_j^T )
```

where `rho` is a reliability factor from leave-one-out residual fit error and
`e` is robust residual energy. The second term is deliberately a stabilizer,
not the main signal: residual magnitude tells us where something went wrong;
gradients tell us which theta directions control the wrongness.

Eigenvectors of `A_hat` define the learned residual directions. Named
hyperparameters enter only after the eigensystem is learned, as loadings used
for human interpretation.

After directions are found, Phase 4 also fits small Taylor diagnostics in the
learned low-dimensional coordinates:

```text
u_j = V^T x_j
delta_ij - offset_i approximately
  a_i + b_i^T u_j + 0.5 u_j^T H_i u_j
```

These Taylor fits are diagnostic only. They are used to tell whether the
residual pattern is smooth enough for stencil chart repair or too noisy to
trust. They are not final local evidence evaluators.

Use it to identify:

```text
top residual directions
local contributors
coherent total-evidence bias
whether residuals are isolated or smooth
exact and direction-stencil chart candidates
```

Named hyperparameters may be reported as loadings after the direction is
learned.

## Repair Policy

Shape repair is triggered only when residual diagnostics show material error:

```text
posterior-weighted centered total residual RMSE > threshold
or max abs centered total residual > threshold
or coherent residual direction has non-negligible posterior mass
```

The first repair action is chart placement.

For isolated local/theta residuals:

```text
add a chart exactly at theta_j
```

For coherent residual directions:

```text
add a local stencil:
  theta_j
  theta_j + epsilon * v
  theta_j - epsilon * v
```

where `v` is a learned whitened residual direction transformed back into theta
space. `epsilon` is chosen from local posterior scale and capped by chart
distance.

New charts are support charts first. They become normalizer-certified only
through the strict graph/edge certification machinery.

## Holdout Gate

After shape repair, rerun raw PMIS and direct probes on holdout pairs.

The repair is accepted only if:

```text
holdout local centered RMSE decreases
holdout total centered RMSE decreases
max absolute centered total residual does not increase materially
raw certification does not regress
normalizer graph residuals do not regress
outer posterior reweight diagnostics do not worsen
```

If holdout fails, the repair is rejected or quarantined. Do not rerun outer SMC
from an unvalidated shape repair.

## Optional Residual Correction

A correction model is allowed only after chart repair has been tried and
holdout residuals show a smooth, low-dimensional, reproducible pattern.

The corrected local factor would be:

```text
ell_corrected_i(theta) =
  ell_hat_i(theta) + delta_hat_i(theta)
```

with propagated uncertainty.

This is dangerous if used too early. It can hide a bad atlas. It should be an
explicit mode with a benchmark gate, not the default path.

## Production Phases

### Phase 1: Probe Pair Selection

Implement a selector that takes a factor set, weighted theta cloud, raw
certification table, graph summary, and compression summary. It returns repair
and holdout probe pairs.

The selector must include certified pairs. This is the key difference from
support repair.

### Phase 2: Direct Probe Runner

Implement a runner that executes replicated direct local SMC at exact selected
local/theta pairs and combines log normalizers on the evidence scale.

### Phase 3: Residual Diagnostics

Build the residual table and summarize local, total, and posterior-weighted
residuals.

### Phase 4: Residual Geometry

Learn low-dimensional residual directions from weighted probe residuals.
Report axis/family loadings only as interpretation.

Implementation:

```text
learn_shape_residual_geometry()
```

This function consumes a shape-probe result, residual diagnostics, or a
residual table plus theta cloud. It returns a
`local_evidence_shape_residual_geometry` object with:

```text
active_energy_matrix
gradient_energy_matrix
location_energy_matrix
directions
direction_vectors_whitened
direction_vectors_theta
direction_loadings
gradient_diagnostics
local_contributors
theta_scores
local_taylor
total_taylor
stencil_candidates
local_stencil_candidates
```

`active_energy_matrix` is the matrix that owns direction discovery. It is built
from robust ridge residual gradients. `location_energy_matrix` is retained
separately and only enters direction discovery through `location_energy_weight`.
This prevents the algorithm from confusing “large residual at a point” with
“axis along which residual error changes.”

`local_stencil_candidates` is the Phase 5 handoff: it contains exact and
directional local/theta chart candidates, with no named-axis logic and no EMC
reference information.

### Phase 5: Shape Repair

Add exact charts or small stencils at residual-driven theta locations. Use the
strict chart insertion and graph certification path. Do not admit single-SMC
direct normalizers without certification.

Implementation:

```text
repair_shape_residual_geometry()
```

This function consumes the Phase 4 geometry, selects from
`geometry$local_stencil_candidates`, constructs an augmented theta cloud, and
calls:

```text
local_atlas_repair_certification_pairs()
```

The selector is deliberately shape-specific and narrow:

```text
exact residual candidates are selected for isolated local/theta misses
directional stencil candidates are selected only when the learned direction
has enough global energy, the local contributes to that direction, and the
local Taylor diagnostic is reliable
```

The repair executor remains the strict chart primitive. Every new chart is
inserted as a support chart first:

```text
normalizer_certified = FALSE
direct_observation_role = support
```

It can become evidence-bearing only through a certified relative-edge path to
the active local graph, or through the existing optional replicated/direct
normalizer gate. Phase 5 does not use derivatives as evidence and does not call
the support-failure selector.

After insertion, Phase 5 can rerun raw PMIS certification on the original theta
cloud for repaired locals or all locals. The default audit scope is repaired
locals; the holdout decision still belongs to Phase 6.

### Phase 6: Holdout Validation

Probe held-out pairs and accept the repair only if residual metrics improve
without certification or graph regressions.

Implementation:

```text
validate_shape_repair_holdout()
```

This function compares the repaired atlas against held-out shape probes. When
the original shape-probe result already contains holdout direct probes, it
reuses their direct SMC estimates and recomputes only the raw PMIS atlas side
after repair. That is the default because the holdout truth should be fixed
while the atlas changes. If requested, it can rerun held-out direct probes.

The gate reports:

```text
pair centered RMSE before/after
posterior-weighted local centered RMSE before/after
posterior-weighted total centered RMSE before/after
max centered total residual before/after
holdout certification before/after
graph residual before/after where available
```

Acceptance requires no holdout certification regression, no material graph
regression, and by default improvement or non-worsening in the centered
residual metrics. This is intentionally a validation gate only; it does not
run outer SMC.

### Phase 7: Outer Reweight Gate

Before rerunning outer SMC, reweight existing outer particles against the
repaired factor set. If reweight ESS is high but posterior diagnostics worsen,
do not rerun outer.

Implementation:

```text
shape_repair_outer_reweight_gate()
```

This function evaluates the repaired factor set on the existing outer particles
and computes:

```text
log_ratio(theta) = log L_repaired(theta) - log L_old(theta)
```

using `fit$loglik_dynamic` as the old value when available. It then forms the
reweighted posterior:

```text
w_new(theta) proportional to w_old(theta) exp(log_ratio(theta))
```

and reports:

```text
reweight ESS and ESS fraction
PSIS k of the log ratios
log evidence delta
self posterior shift between old and reweighted particles
optional reference-posterior worsening
```

If holdout validation is required, missing or failed holdout validation rejects
the repair before reweighting. A low-ESS or high-PSIS result means a full outer
rerun is required before trusting the repaired surface. A high-ESS result with
worse reference posterior diagnostics rejects the repair rather than rerunning
outer. A high-ESS result with small posterior shift accepts the reweighted
state and explicitly avoids an automatic outer rerun.

### Phase 8: Explicit Outer Rerun

Only after holdout and reweight gates pass, rerun outer SMC from the frozen
shape-calibrated factor set.

## Acceptance Metrics

At minimum, report:

```text
number of probe pairs
number of locals probed
repair/holdout split
probe replicate SD
local centered residual RMSE
total centered residual RMSE
max absolute centered total residual
raw certification fraction before/after
graph residual before/after
outer reweight ESS
posterior comparison before/after
log evidence before/after
```

The shape layer is successful only if held-out residuals improve and the final
posterior moves closer to the reference benchmark without degrading unrelated
parameters.

## Implemented Phases

### Phase 1: Probe Selection

```text
select_shape_probe_pairs()
```

This function should be pure and cheap. It should not run SMC or mutate the
factor set.

Inputs:

```text
factor_set
theta
theta_weights
raw_certification_table = NULL
graph_summary = NULL
compression_summary = NULL
control = list()
```

Outputs:

```text
repair_pairs
holdout_pairs
scored_table
selection_summary
```

The immediate benchmark is to run it on the current EMC checkpoint and verify
that it selects a small, interpretable set of certified high-risk pairs,
including hard locals already identified by diagnostics, without using any
hardcoded parameter names.

### Phase 2: Direct Probe Residuals

```text
run_shape_probe_pairs()
```

This function consumes the selected repair and holdout pairs, runs replicated
direct local SMC only at those exact local/theta pairs, combines replicated
normalizers on the evidence scale, and returns:

```text
residuals
replicates
probe_summary
summary
cloud
settings
```

The residual table carries the raw atlas PMIS estimate, the direct probe
estimate, probe uncertainty, atlas uncertainty, standardized residuals, local
diagnostics, theta provenance, posterior weight, and the repair/holdout split.

This phase is still measurement-only. It does not mutate the atlas, repair
charts, use EMC reference draws, or call the old posterior-region calibration
machinery.

### Phase 3: Residual Diagnostics

```text
diagnose_shape_probe_residuals()
```

This function consumes a shape-probe result or residual table and returns:

```text
probe_set_summary
local_summary
theta_summary
```

The local summary centers residuals within each local to separate additive
local offsets from posterior-shape error. The theta summary sums observed local
residuals by theta to diagnose the probed product-evidence error and records
coverage, so sparse probes are not mistaken for a complete total-evidence
audit. The probe-set summary reports posterior-weighted local centered RMSE,
posterior-weighted total centered RMSE, max centered total residual, total
coverage, and replicate uncertainty summaries.

`run_shape_probe_pairs()` now attaches these diagnostics directly as:

```text
summary
local_diagnostics
theta_diagnostics
diagnostics
```

### Phase 4: Residual Geometry

```text
learn_shape_residual_geometry()
```

This is now called by `run_shape_probe_pairs()`, so a shape-probe result carries
its geometry as:

```text
geometry
```

The implementation uses the active-subspace estimator derived above:

```text
robust local residual gradients + robust total residual gradient
```

with residual-location energy retained only as a stabilizing term. The returned
object exposes both pieces separately, reports learned direction loadings,
attributes each direction to local contributors, fits low-dimensional Taylor
diagnostics, and emits exact/directional chart candidates for the next repair
phase.

This phase is still measurement and design only. It does not alter local
atlases, correct local evidence values, rerun outer SMC, or consult the EMC
posterior benchmark.

### Phase 5: Shape Repair

```text
repair_shape_residual_geometry()
```

This function is now the only production Phase 5 entry point. It:

```text
1. validates the shape-probe geometry;
2. scores exact and directional local stencil candidates;
3. applies strict total, per-local, per-base-theta, and per-direction caps;
4. builds an augmented theta cloud for selected candidate locations;
5. converts candidates to exact `(local_pos, theta_row)` repair pairs;
6. delegates chart insertion/certification to
   local_atlas_repair_certification_pairs();
7. reruns raw PMIS certification on the original theta cloud for the requested
   audit scope;
8. returns the updated factor set, selected candidates, repair probes,
   certification before/after summaries, and graph summary.
```

The function does not re-rank with support diagnostics, does not use compressed
PMIS for certification, and does not turn residual derivatives into final
local evidence values. Its output is ready for Phase 6 holdout validation.

### Phase 6: Holdout Validation

```text
validate_shape_repair_holdout()
```

This function validates a `local_evidence_shape_repair_result` against held-out
probe pairs. It reuses held-out direct SMC estimates from the original
`run_shape_probe_pairs()` result by default and recomputes raw PMIS under the
repaired factor set. It returns:

```text
accepted
summary
pre_residuals
post_residuals
pre_diagnostics
post_diagnostics
```

The summary contains the exact before/after centered residual metrics and
certification counts used by the gate. It does not mutate the atlas and does
not rerun outer.

### Phase 7: Outer Reweight Gate

```text
shape_repair_outer_reweight_gate()
```

This function is the mandatory gate before any full outer rerun. It uses the
existing outer particles, evaluates the repaired factor set, computes the
importance reweighting diagnostics, and returns a decision:

```text
accept_reweight_no_outer_rerun
accept_reweight_outer_rerun_optional
outer_rerun_recommended_moderate_reweight_ess
rerun_outer_required_low_reweight_quality
reject_repair_holdout_failed
reject_repair_missing_holdout_validation
reject_repair_reference_worsened
```

It does not run outer SMC. It either accepts the reweighted posterior as
diagnostically adequate, marks a full outer rerun as necessary, or rejects the
repair when holdout/reference diagnostics say the calibrated surface got worse.
