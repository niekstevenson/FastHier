# Raw Local Evidence Certification

## Problem

The current framework can now run a full EMC hierarchy with raw local
particle-MIS and no compression. That run is much better than the compressed
run, but it still misses the EMC posterior:

```text
raw 1500 max standardized posterior mean error: 0.692
raw 1500 mean standardized posterior mean error: 0.211
```

The remaining failure is not outer SMC Monte Carlo noise. The outer SMC MCSE is
small. It is also not primarily compression, because this diagnostic is with
compression disabled.

The remaining failure is:

```text
outer SMC is consuming raw local evidence surfaces that are still uncertified
over posterior-relevant theta regions.
```

On a diagnostic cloud of 800 EMC-reference theta points and 800 workflow
posterior theta points, the raw atlas had:

```text
EMC-reference theta:
  any uncertified local evidence: 69.6%
  mean uncertified locals per theta: 1.595

workflow-posterior theta:
  any uncertified local evidence: 37.3%
  mean uncertified locals per theta: 1.054
```

So the algorithm's own posterior sits in a better-certified region than the
reference posterior. That is a certification-induced selection bias.

## Critical Revision Of Previous Suggestion

The earlier suggestion to repair specific named axes such as `mu_sv` and
`log_sigma2_sv` is not acceptable as an algorithm.

It is useful for diagnosis only. It explains the EMC miss, but it would bake
benchmark knowledge into the framework. A production method cannot know that an
axis named `sv` is important.

The general algorithm must infer hard directions from:

1. the theta distribution currently being trusted;
2. raw particle-MIS diagnostics;
3. local evidence sensitivity;
4. posterior or proposal impact;
5. coherent failure geometry across locals.

No step may require a hand-written list of hard parameters. Manual focus axes
are allowed only for diagnostic scripts.

The second weakness in the earlier suggestion is that "repair all uncertified
pairs" is too expensive and not principled. Raw PMIS diagnostics can be
conservative; a high PSIS or low ESS point is not automatically a posterior
error. Repair must be prioritized by estimated impact on the total evidence
surface.

The third weakness is that a pre-outer gate alone is not enough. If a pilot
outer run moves mass into a region not covered by the original proposal audit,
the framework needs an explicit rebuild/rerun loop. It should not silently
continue with local surfaces that were certified for the wrong theta
distribution.

## Certification Target

For each local:

```text
ell_i(theta) = log m_i(theta)
             = log int L_i(alpha_i) p(alpha_i | theta) d alpha_i
```

The target is not to certify every possible theta. The target is:

```text
certify ell_i(theta) wherever errors can materially affect
the posterior or total model evidence.
```

The certified object is the raw particle-MIS local evidence estimator. Any
compressed representation must later be certified against this raw surface.
Compression is downstream and cannot certify the local evidence surface.

## Core Principle

A local/theta query is trusted only if both are true:

1. it has enough local chart support under the raw particle-MIS estimator;
2. its uncertainty or instability would not materially reshape the total
   evidence over the theta distribution being used.

This is different from a smoke-test gate. It is a posterior-shape gate.

## Objects

### `local_evidence_certification_cloud`

A theta cloud with weights and provenance:

```text
theta_id
theta
weight
source
round
```

Allowed production sources:

```text
q0 proposal samples
theta design anchors
inflated q0/proposal tails
pilot outer posterior samples
pilot outer posterior tails
high-leverage synthetic profiles selected automatically
```

Benchmark-only sources such as EMC/Stan reference draws may be used for
diagnostics, but never for production repair selection.

### `local_evidence_certification_table`

One row per local/theta query:

```text
local
theta_id
theta_source
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
selected_for_repair
repair_status
```

This table is the source of truth for whether outer SMC may consume the factor
set.

### `certification_geometry`

An automatically learned description of where certification is weak:

```text
whitened theta center and scale
failure-weighted covariance
top failure directions
local contributions by direction
axis/family loadings for interpretation only
```

The important point: directions are vectors in theta space. Named axes are only
labels attached after the fact.

## Failure Score

For query `(i, theta_j)`, define a raw diagnostic severity:

```text
D_ij =
  a1 * I(status != certified)
  + a2 * max(0, min_ess_frac - ess_frac) / min_ess_frac
  + a3 * max(0, psis_k - psis_threshold)
  + a4 * log1p(min_covering_distance)
  + a5 * raw_se
```

The posterior/proposal impact weight is:

```text
W_j = normalized theta weight
```

The local sensitivity weight should use chart derivative information when
available:

```text
S_ij = 1 + || score_i(theta_near) ||_G + curvature_scale_i(theta_j)
```

The repair score is:

```text
R_ij = W_j * D_ij * S_ij
```

This score decides what to repair. It does not decide what is true.

## Automatic Failure Geometry

Given a certification cloud, whiten theta:

```text
z_j = Sigma_theta^{-1/2} (theta_j - center_theta)
```

Aggregate per-theta failure burden:

```text
B_j = sum_i R_ij
```

Then compute a weighted failure covariance:

```text
C_fail = sum_j B_j z_j z_j^T / sum_j B_j
```

The top eigenvectors of `C_fail` are the hard directions. They are not chosen by
parameter name.

Use these directions to create additional audit/probe profiles:

```text
theta_center +/- q * sqrt(lambda_k) * v_k
```

where `q` comes from weighted posterior/proposal quantiles, not from hard-coded
axis names.

For interpretability only, print the largest coordinate loadings:

```text
direction 1: +0.72 log_sigma2_sv -0.41 mu_sv ...
```

Those labels must not drive the algorithm.

## Certification Gate

A factor set is certified for a theta cloud only if:

```text
1. weighted fraction of theta with any uncertified local <= threshold_any
2. weighted mean uncertified locals per theta <= threshold_mean
3. weighted q90 uncertified locals per theta <= threshold_q90
4. no high-weight theta has more than threshold_max uncertified locals
5. top failure-direction burden is below threshold_direction
6. worst local posterior-impact burden is below threshold_local
```

Suggested starting thresholds for production-quality runs:

```text
threshold_any = 0.05
threshold_mean = 0.10
threshold_q90 = 0
threshold_max = 1
threshold_direction = 0.05 * total_burden_round0
threshold_local = 0.02 * total_burden_round0
```

For loose benchmark runs, thresholds may be relaxed, but the result must be
reported as uncertified.

## Repair Policy

Repair is not "add more anchors everywhere."

Each certification round:

1. evaluate raw PMIS on the certification cloud;
2. compute `R_ij`;
3. learn failure directions from `C_fail`;
4. add direction-profile theta points to the audit cloud;
5. rank local/theta pairs by `R_ij`;
6. select a capped batch that covers:
   - top local burdens;
   - top theta burdens;
   - top failure directions;
   - distinct locals, not just repeated repairs for one local;
7. add repair charts exactly at selected theta points;
8. recalibrate the local graph;
9. rerun raw certification.

Stop only when the certification gate passes or a declared budget is exhausted.
If the budget is exhausted, the factor set is returned with status
`uncertified`, and outer SMC should not be called by default.

## Outer Loop

The production workflow should be:

```text
1. build chart design scout
2. build local atlases
3. raw pre-outer certification against q0/proposal/tail cloud
4. repair until the pre-outer gate passes
5. run pilot outer SMC
6. raw post-pilot certification against pilot posterior/tails
7. repair until the post-pilot gate passes
8. rerun final outer SMC from the certified factor set
9. optionally compress and certify compression against raw PMIS
```

The current loose workflow skipped the hard version of steps 4 and 7. That is
why the posterior was allowed to settle in a better-certified but biased region.

## Compression Rule

Compression is disabled until raw certification passes.

After raw certification:

1. build function-aware sparse quadrature per local;
2. evaluate compressed vs raw PMIS on the same certification cloud plus
   held-out direction profiles;
3. compute total compression error:

```text
Delta(theta) =
  sum_i ell_i_compressed(theta) - sum_i ell_i_raw(theta)
```

4. reject compression or increase per-local `K` if:

```text
weighted sd Delta(theta) > threshold_total_sd
weighted range Delta(theta) > threshold_total_range
any top failure direction has nontrivial correlation with Delta(theta)
any local compression error exceeds local threshold
```

Compression must never reduce uncertified raw PMIS queries to certified
queries. It may only approximate already-certified raw evidence.

## Concrete Implementation Plan

### Phase 1: Certification Table

Implement:

```text
build_local_evidence_certification_cloud()
evaluate_raw_local_evidence_certification()
summarize_raw_local_evidence_certification()
```

Owner:

```text
local_charts.R
```

Required behavior:

- uses raw particle-MIS only;
- ignores compression fields;
- returns per-local/per-theta diagnostics;
- computes per-theta and per-local burden summaries;
- writes CSV/RDS artifacts from the EMC benchmark.

### Phase 2: Automatic Failure Geometry

Implement:

```text
learn_certification_failure_geometry()
add_failure_direction_profiles()
```

Required behavior:

- whiten theta by the weighted certification cloud covariance;
- compute burden-weighted covariance;
- return top eigen-directions;
- produce new theta profiles along those directions;
- attach human-readable loadings, but do not require hard-coded names.

### Phase 3: Impact-Based Repair Selection

Replace the current tiny repair selector with:

```text
select_certification_repairs()
```

Inputs:

```text
certification_table
failure_geometry
max_repairs
max_repairs_per_local
max_repairs_per_direction
```

Selection must cover:

- highest local burden;
- highest theta burden;
- top failure directions;
- high-weight theta points;
- repeated hard locals with a per-local cap.

This prevents the current failure where thousands of bad pairs are found but
only a narrow set of repairs is installed.

### Phase 4: Raw Certification Loop

Implement:

```text
certify_and_repair_raw_atlas()
```

Pseudo-flow:

```text
cloud <- build_local_evidence_certification_cloud(...)
for round in seq_len(max_rounds):
  table <- evaluate_raw_local_evidence_certification(factor_set, cloud)
  summary <- summarize_raw_local_evidence_certification(table)
  if gate_passed(summary):
    return(certified factor_set)

  geometry <- learn_certification_failure_geometry(table, cloud)
  cloud <- add_failure_direction_profiles(cloud, geometry)
  repairs <- select_certification_repairs(table, geometry)
  factor_set <- add_repair_charts_and_recalibrate(factor_set, repairs)

return uncertified factor_set with table and summary
```

Default production behavior:

```text
if certification_status != certified:
  stop before outer SMC
```

Loose benchmark behavior:

```text
allow_outer_on_uncertified = TRUE
```

but the result must report the certification failure.

### Phase 5: Pilot Outer And Rerun

Add a two-stage outer driver:

```text
run_certified_outer_workflow()
```

Flow:

```text
pre_outer certified factor set
pilot outer SMC
post_pilot certification cloud from pilot posterior and inflated tails
repair if needed
final outer SMC from rebuilt certified factor set
```

No automatic infinite rerun. The maximum number of certify/outer cycles must be
explicit.

### Phase 6: Compression Re-Entry

Only after Phase 5:

```text
compress_certified_local_evidence()
certify_compressed_against_raw()
```

If compression fails, the framework should fall back to raw PMIS for hard
locals and compress only easy locals.

The first useful version is mixed:

```text
hard locals: raw PMIS
easy certified locals: compressed quadrature
```

This is more principled than forcing a global `K`.

## First Step

The first code change should be Phase 1 only:

```text
build and save a raw local evidence certification table for a factor set
and theta cloud.
```

Do not change repair yet.
Do not change compression yet.
Do not run another full benchmark as the first step.

Reason:

The current diagnosis came from scratch scripts. The framework needs this as a
first-class artifact before repair policy can be made principled. Once the table
exists, every subsequent benchmark can answer:

```text
Did this run pass raw local evidence certification?
If not, which locals and theta regions failed?
Did repair reduce posterior-impact burden?
Did compression preserve the certified raw surface?
```

Without this table, we will keep tuning budgets blindly.

## Expected Improvement Mechanism

The posterior should improve because the current error mechanism is:

```text
uncertified local evidence is more common in the reference-like theta region
than in the workflow posterior region.
```

Certification repair changes the local atlas before outer SMC so that the
posterior is not rewarded for moving into better-certified but wrong regions.

This is not a guarantee that the posterior will match EMC. It is the next
necessary condition: outer SMC must consume a local evidence surface whose
certification burden is low and balanced over the theta region that can affect
the posterior and evidence.
