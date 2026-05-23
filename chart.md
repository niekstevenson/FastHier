# Certified Local Evidence Charts

This document defines the replacement architecture for Step 3/4 of the
hierarchical evidence framework. It is intentionally strict: a local state is
not allowed to influence the outer posterior until its local evidence
contribution is certified.

The current failure mode was:

1. a local SMC run created a node with a wrong endpoint log normalizer;
2. ESS and final particle quality looked acceptable;
3. MBAR shifted the node into apparent consistency with the bank;
4. deterministic-mixture evaluation used the node as if it were valid;
5. the outer SMC correctly sampled a distorted local evidence surface.

The fix is not more blind anchors. The fix is to change the object being built:
each local must build a certified atlas for

```text
ell_i(theta) = log m_i(theta)
             = log integral L_i(alpha) p(alpha | theta) d alpha.
```

## Design Contract

For each local `i`, the output of Step 3/4 is an evaluable surface for
`ell_i(theta)` with explicit certification status.

The surface must satisfy all of the following.

1. It has at least one absolute reference chart.
2. Every active non-root chart is connected to the root by accepted relative
   evidence edges.
3. The graph of active charts has small cycle residuals.
4. Local value, score, and curvature predictions agree on posterior-relevant
   regions.
5. Uncertified regions return `uncertified`, not a fabricated local evidence
   value.
6. MBAR is never allowed to promote an uncertified node into the active surface.
7. Outer SMC consumes only frozen active certified surfaces.

If these conditions are not met, the local surface is incomplete. The framework
must refine or quarantine; it must not proceed as if the posterior is certified.

## Core Objects

### Chart

A chart is a local approximation to `ell_i(theta)` around an anchor `theta_a`.

It stores:

```text
local_id
chart_id
theta_anchor
status: candidate | active | quarantined | retired
logZ_abs
logZ_abs_se
alpha_particles
weights
score
curvature
local_geometry
source_run_ids
edge_ids
diagnostics
```

Only `active` charts may be used in final local evidence evaluation.

### Edge

An edge certifies a relative local evidence ratio between two charts:

```text
Delta_ab = ell_i(theta_b) - ell_i(theta_a).
```

It stores:

```text
local_id
from_chart
to_chart
delta
se
method: BAR | bridge | SMC-path | Taylor-crosscheck
overlap_ess
psis_k
forward_delta
reverse_delta
taylor_delta
cycle_residual
status: candidate | active | rejected
diagnostics
```

Only `active` edges may be used in the chart graph normalizer solve.

### Atlas

An atlas is the certified chart graph for one local:

```text
local_id
root_chart_id
charts
edges
normalizer_solution
coverage_region
surface_uncertainty_model
diagnostics
```

The atlas is the local object consumed by the outer population sampler.

## Mathematical Information Stored Per Chart

For a diagonal Gaussian population model,

```text
alpha_j | theta ~ Normal(mu_j, sigma2_j)
rho_j = log sigma2_j.
```

The local evidence score is available from particles by Fisher's identity:

```text
grad_theta ell_i(theta)
  = E_{alpha | y_i, theta}[grad_theta log p(alpha | theta)].
```

The useful per-coordinate score terms are:

```text
d/d mu_j log p(alpha | theta)
  = (alpha_j - mu_j) / sigma2_j

d/d rho_j log p(alpha | theta)
  = -1/2 + (alpha_j - mu_j)^2 / (2 sigma2_j).
```

The local curvature is available from Louis' identity:

```text
H ell_i(theta)
  = E[H log p(alpha | theta)]
    + Cov[grad log p(alpha | theta)].
```

This turns each local SMC run into more than one scalar logZ. It gives a local
quadratic chart of the evidence surface. The `sv` failure was a surface-shape
failure, so score and curvature must be first-class diagnostics.

## Strict Step 3/4 Algorithm

### 1. Build One Absolute Root Per Local

For each local `i`, choose a central posterior design point `theta_root`.

The root chart is the only place where the framework pays for an absolute
normalizer by default.

Requirements:

1. run a strong local SMC at `theta_root`;
2. compute `logZ_abs`, score, curvature, and particle diagnostics;
3. run at least one cheap independent confirmation if the root has high
   evidence leverage or weak path diagnostics;
4. reject the root if logZ confirmation, path CESS, or move diagnostics fail.

The root is active only when its absolute evidence is credible.

### 2. Add New Anchors As Candidate Charts

New anchors are proposed by theta design logic:

1. posterior quantiles and ridge directions;
2. inflated posterior axes;
3. score/curvature disagreement locations;
4. outer posterior support not covered by active charts;
5. high evidence-uncertainty locations.

A candidate chart may run SMC, but its raw endpoint logZ is not trusted as an
absolute fact. It is just one noisy measurement.

### 3. Certify Candidate Charts By Relative Edges

A candidate chart `b` becomes active only through accepted relative evidence
edges to active neighbor charts.

For an active chart `a` and candidate `b`, estimate:

```text
Delta_ab = ell_i(theta_b) - ell_i(theta_a).
```

Use paired bridge/BAR style estimators from both endpoint particle clouds.
Accept the edge only if:

1. forward and reverse estimates agree within tolerance;
2. PSIS `k` is below threshold;
3. overlap ESS is above threshold;
4. Taylor-predicted delta from chart `a` is not materially inconsistent;
5. Taylor-predicted delta from chart `b` is not materially inconsistent;
6. adding the edge does not create a large graph cycle residual.

If no edge passes, insert intermediate anchors between `a` and `b`. Do not
promote `b` by MBAR.

### 4. Solve Chart Normalizers By Graph Least Squares

Given active edge measurements:

```text
z_b - z_a = Delta_ab + epsilon_ab,
```

solve for chart log normalizers `z_s` using weighted least squares, fixing the
root chart as the absolute reference.

Diagnostics are mandatory:

1. edge residuals;
2. cycle residuals;
3. per-chart normalizer uncertainty;
4. leverage of each edge;
5. sensitivity of `ell_i(theta)` to removing each chart.

If graph residuals fail, the atlas is not certified.

### 5. Evaluate `ell_i(theta)` From Active Charts Only

Evaluation uses certified local information in this order:

1. choose nearby active charts whose local geometry covers `theta`;
2. compute quadratic chart predictions from value, score, and curvature;
3. combine predictions by uncertainty-weighted local interpolation;
4. use particle bridge/DMIS evaluation only when overlap diagnostics pass;
5. return uncertainty with the value.

If no active chart covers `theta`, return:

```text
status = uncertified
reason = no_active_chart_coverage
```

Do not extrapolate silently.

### 6. Outer SMC Uses Frozen Certified Atlases

The outer sampler receives a factor set made only from certified local atlases.

Rules:

1. no candidate or quarantined chart enters the outer target;
2. no local refinement happens during an outer SMC run;
3. if outer particles enter uncertified theta regions, stop, expand atlases,
   freeze again, and rerun outer SMC;
4. final evidence requires both posterior certification and absolute local
   normalizer uncertainty accounting.

## What MBAR Is Allowed To Do

MBAR is useful only after node and edge validity are established.

Allowed:

1. audit active chart graph consistency;
2. refine relative normalizers within a certified connected component;
3. estimate uncertainty when overlap is demonstrably good.

Forbidden:

1. promoting a candidate chart;
2. overriding failed edge diagnostics;
3. hiding a bad raw SMC normalizer behind a global shift;
4. entering any uncertified chart into the final deterministic mixture.

## Certification Metrics

The framework must report these per local:

```text
n_active_charts
n_candidate_charts
n_quarantined_charts
n_active_edges
root_logZ_se
max_edge_residual
max_cycle_residual
max_chart_leave_one_out_delta
max_surface_uncertainty_on_outer_particles
fraction_outer_particles_uncertified
posterior_leverage_of_uncertified_regions
estimated_local_log_evidence_uncertainty
```

The global run must report:

```text
max_local_surface_uncertainty
sum_local_log_evidence_uncertainty
outer_particles_uncertified_fraction
posterior_marginal_differences_vs_reference_if_available
frozen_outer_rerun_stability
```

The posterior is not certified if any high posterior mass region is uncertified.

## Expected Posterior Improvement

This design should improve the posterior because it blocks the observed failure:

```text
bad local endpoint logZ
  -> MBAR shift
  -> bad active DMIS denominator
  -> distorted local evidence surface
  -> biased outer posterior
```

Under the chart design, the same bad endpoint logZ becomes:

```text
bad local endpoint logZ
  -> edge/cycle/Taylor failure
  -> chart quarantined or bridged through intermediates
  -> no contribution to active local surface
```

The target improvement criteria for EMC are:

1. `sigma2_sv` standardized mean error moves toward zero without worsening
   other variance dimensions;
2. `log_sigma2_sv` posterior width and left/right ridge geometry match EMC
   better;
3. frozen outer reruns remain stable;
4. fresh endpoint SMC probes agree with chart evaluation in posterior-relevant
   regions;
5. total log evidence uncertainty is explicitly bounded.

Improvement is not measured by fewer audit failures alone. It is measured by
posterior shape and local evidence calibration.

## Implementation Plan

### Phase 1: New Data Model

Add compact chart objects:

```text
new_local_chart()
new_local_edge()
new_local_atlas()
validate_local_chart()
validate_local_edge()
validate_local_atlas()
```

Do not wrap existing bank nodes as charts unless all required fields are real.
If a field is not computed, the chart is not active.

### Phase 2: Score And Curvature

Implement diagonal Gaussian chart derivatives:

```text
local_chart_score(theta, alpha, weights, population_model)
local_chart_curvature(theta, alpha, weights, population_model)
```

Use weighted particle expectations. Keep this population-model specific at
first; do not add a fake generic interface.

### Phase 3: Root Construction

Implement:

```text
build_local_root_chart()
```

The root builder runs local SMC, computes value/score/curvature, and optionally
performs a cheap confirmation. It returns `active` or fails loudly.

### Phase 4: Candidate Construction

Implement:

```text
build_candidate_chart()
```

Candidate charts may contain raw SMC logZ, but `status` must remain
`candidate` until edge certification passes.

### Phase 5: Edge Certification

Implement:

```text
estimate_chart_edge()
certify_chart_edge()
```

Start with two-sided bridge/BAR-style ratio estimates and Taylor checks.
Required outputs are `delta`, `se`, overlap ESS, PSIS `k`, forward/reverse
disagreement, Taylor disagreement, and status.

### Phase 6: Graph Solve

Implement:

```text
solve_atlas_normalizers()
atlas_cycle_diagnostics()
```

Use weighted least squares with root fixed. Reject graph updates that create
large residuals.

### Phase 7: Atlas Evaluation

Implement:

```text
evaluate_local_atlas(theta)
evaluate_local_atlas_many(theta)
```

Return:

```text
log_marginal
se
status
nearest_charts
diagnostics
```

Uncovered theta must return `uncertified`.

### Phase 8: Outer Integration

Build a factor-set equivalent that calls atlas evaluators instead of active
bank DMIS over uncertified nodes.

Outer SMC must stop if any proposed or retained high-weight particle has an
uncertified local factor.

### Phase 9: Benchmark Gate

The first benchmark gate is EMC.

Required checks:

1. compare posterior marginals against EMC;
2. run frozen outer rerun against the same atlases;
3. run fresh endpoint SMC probes at selected EMC and workflow theta points;
4. report chart graph residuals and uncertainty;
5. verify that removing quarantined nodes does not materially change the
   active posterior.

The implementation is not accepted unless EMC posterior shape improves relative
to the current clean state, especially in `mu_sv`, `sigma2_sv`,
`sigma2_v_LogFreq`, and `sigma2_v`.

## Non-Goals

1. Do not add another global MBAR pass as the main fix.
2. Do not increase all particle counts as the main fix.
3. Do not replicate every anchor.
4. Do not allow ESS-only certification.
5. Do not hide uncertified theta regions behind broad defensive mixtures.
6. Do not optimize for matching one benchmark by special-casing EMC.

## Final Architecture Summary

Step 3/4 becomes:

```text
for each local:
  build one certified absolute root chart
  propose candidate theta anchors
  run local SMC for candidate charts
  compute score and curvature
  certify relative edges to active charts
  solve active chart graph normalizers
  quarantine failed candidates
  expose only certified atlas evaluator to outer SMC

outer SMC:
  sample theta using frozen certified atlas factors
  stop if posterior mass reaches uncertified atlas regions
  refine atlases outside the run
  rerun outer against frozen certified factors
```

This is the strict replacement for the current local bank path.
