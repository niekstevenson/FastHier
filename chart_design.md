# Adaptive Chart Design

This document defines the production replacement for all-axis local chart
design.

The all-axis design is not a default algorithm. It is a diagnostic stress test.
It answers one question:

```text
If every population axis gets some chart coverage, does the atlas recover?
```

For the EMC run, the answer was no. The design created roughly 1256 active
charts, made outer SMC expensive, and still missed `sigma2_sv`. That means the
main problem is not missing a named axis. The problem is inefficient allocation:
too many charts are spent on easy directions while hard local evidence curvature
and normalizer error remain under-resolved.

The production design must be adaptive, anisotropic, local-specific, and cheap
to evaluate.

## Role Of All-Axis Runs

All-axis runs are allowed only as diagnostics.

Use them to:

1. rule out a trivial omitted-axis explanation;
2. measure worst-case chart-count scaling;
3. expose which posterior dimensions remain wrong under broad coverage;
4. provide a negative control for more targeted designs.

Do not use them as the default benchmark setting. They are too expensive and
they dilute chart budget across directions that are not equally important.

## Design Target

For each local `i`, we need an evaluable surface

```text
ell_i(theta) = log m_i(theta)
             = log integral L_i(alpha) p(alpha | theta) d alpha.
```

The design problem is not "cover every hyperparameter axis." The design problem
is:

```text
Place local charts where errors in ell_i(theta) can materially change the
outer posterior or total model evidence.
```

That means chart allocation should depend on:

1. local sensitivity of `ell_i(theta)`;
2. local curvature of `ell_i(theta)`;
3. proposal/posterior-relevant theta mass;
4. particle-MIS instability;
5. graph normalizer uncertainty;
6. expected contribution to total evidence error.

## Mathematical Principle

At a chart anchor `theta_a`, each local stores value, score, and curvature:

```text
ell_i(theta_a)
g_i(theta_a) = grad_theta ell_i(theta_a)
H_i(theta_a) = grad_theta^2 ell_i(theta_a)
```

For nearby theta,

```text
ell_i(theta) approx ell_i(theta_a)
  + g_i(theta_a)^T delta
  + 0.5 delta^T H_i(theta_a) delta

delta = theta - theta_a.
```

The chart design should use this derivative information to decide where charts
are needed. Derivatives may guide design, but they do not certify final values.
Final values still come from certified chart normalizers and particle-MIS.

## Hyperparameter Families

For the diagonal Gaussian population model, hyperparameters naturally group by
alpha dimension:

```text
family k = { mu_k, log_sigma2_k }
```

The design unit should usually be a family, not a single scalar coordinate.

This matters because the failures are often manifold-shaped. For example,
`mu_sv` and `log_sigma2_sv` interact: moving one changes which alpha particles
matter under the other. Selecting only `mu_sv` or only `log_sigma2_sv` can still
miss the local evidence shape.

## Production Algorithm

### 1. Root Scout

Build only a root/core chart for every local, or for a pilot subset if the local
count is very large.

This scout stage must be cheap:

```text
root_M << production_M
candidate charts = none or very few
outer SMC = not run
validation gate = not run
```

The output is local derivative information and local SMC reliability:

```text
local_id
root theta
score vector
curvature matrix
logZ_se
path ESS
acceptance
local likelihood cost
```

### 2. Direction Risk Score

For each local `i` and hyperparameter family `f`, estimate a design risk score.

A practical score is:

```text
R_if =
  q_weighted_taylor_variation_if
  * local_reliability_weight_i
  * posterior_design_weight_f
```

where `q_weighted_taylor_variation_if` is computed over the current theta design
cloud, initially q0:

```text
delta_f = projected theta movement in family f

taylor_if(theta) =
  g_if^T delta_f + 0.5 delta_f^T H_iff delta_f

q_weighted_taylor_variation_if =
  weighted_sd_q0(taylor_if(theta))
```

The reliability weight should increase when the local root chart is noisy or
unstable:

```text
local_reliability_weight_i =
  1
  + c1 * logZ_se_i
  + c2 * weak_path_indicator_i
  + c3 * likelihood_cost_i
```

This is not a final accuracy certificate. It is only a budget allocator.

### 3. Select Active Families

Aggregate over locals:

```text
R_f = sum_i R_if
```

Select families by one of two rules:

```text
top_k families
or
all families with R_f >= tau * max_f R_f
```

The default should be threshold-based with a maximum cap. For example:

```text
min_families = 1
max_families = 4
relative_threshold = 0.25
```

This replaces hand-set `focus_hyper_names` as the production default.

Manual focus parameters are still useful for diagnostics, but they should not
be required for an independent hierarchical run.

### 4. Allocate Local Chart Budgets

Do not give every local the same number of charts.

Allocate per-local budget:

```text
B_i = B_min + round((B_total - n_locals * B_min) * R_i / sum_i R_i)
R_i = sum_f R_if over selected families
```

Hard locals get more charts. Easy locals stay small.

This directly attacks the current waste: all-axis designs spend chart budget on
easy locals and easy axes while hard locals still underfit the variance surface.

### 5. Generate Candidate Anchors

For each selected family, generate candidate anchors from the theta cloud:

```text
family quantile profiles: 0.05, 0.5, 0.95
local high-risk profiles
farthest points in Fisher/theta metric restricted to selected families
optional tail points from defensive q0 mass
```

Then choose each local's anchors from the global candidate pool using its local
risk scores. The output is a local-specific theta design, not one global design
blindly copied to every local.

### 6. Build Production Atlases

For each local:

1. keep the scout root chart if it passed root certification;
2. add the local-specific selected anchors;
3. activate charts only through certified edges;
4. solve normalizers robustly;
5. quarantine charts that fail graph or normalizer checks.

Fresh SMC remains the default for chart creation. Bridge/transport can be added
later only as a strict acceleration path, not as a certificate.

### 7. Sparse Final Evaluation

Outer SMC should not evaluate every theta against every active chart if most
charts are irrelevant.

Final local evaluation should use a sparse responsible-chart set:

```text
nearest charts in selected family metric
plus exact certified anchors
plus any chart required by coverage/PSIS guard
```

The sparse estimate must be audited against the full particle-MIS estimate on a
small theta sample. If sparse-vs-full differences exceed tolerance, increase
responsible chart count or mark the theta/local query uncertified.

Full particle-MIS remains the audit/reference evaluator. Sparse particle-MIS is
the production evaluator only after it passes this audit.

### 8. Validation Is External

Validation is not part of fitting.

The fitting script should produce:

```text
frozen atlases
outer posterior
posterior comparison
checkpoint files
build/calibration history
```

Separate validation scripts should run:

```text
frozen outer rerun
fresh endpoint probes
nested SMC comparison on selected hard locals
posterior overlap diagnostics
```

This keeps fitting fast and makes validation explicit.

## What To Delete Or Demote

Demote:

```text
all-axis design
```

to diagnostic-only.

Do not use:

```text
axis_count = all hyperparameters
```

as a production setting.

Do not run:

```text
benchmark gate
frozen outer rerun
fresh endpoint probes
```

inside the fitting script.

Do not automatically rerun outer SMC after posterior calibration. Calibration
creates a new frozen factor set; outer rerun is a separate explicit operation.

## Implementation Plan

### Phase 1: Chart Design Plan Object

Add a `chart_design_plan` object with:

```text
theta_root
theta_cloud
family_scores
local_scores
selected_families
local_anchor_budgets
candidate_theta
local_theta_designs
diagnostics
```

This object should be serializable and saved in checkpoints.

Acceptance criterion:

```text
The benchmark output records why each family and local received its chart
budget.
```

### Phase 2: Root Scout Stage

Split atlas construction into:

```text
root scout build
design planning
production atlas expansion
```

The scout stage builds root charts only and returns score/curvature diagnostics.

Acceptance criterion:

```text
An EMC run can stop after scout and write family/local risk tables without
running outer SMC.
```

### Phase 3: Automatic Family Selection

Implement derivative-based family risk scoring over the theta cloud.

Inputs:

```text
root scout charts
theta_cloud
population model
q0/proposal weights if available
```

Output:

```text
selected hyperparameter families
```

Acceptance criterion:

```text
The hard EMC run selects the dominant difficult families without using EMC
posterior comparison or hand-set focus parameters.
```

### Phase 4: Local-Specific Anchor Budgets

Replace one global theta design copied to every local with local-specific
designs.

Acceptance criterion:

```text
Median active charts per local decreases, while hard locals still receive enough
anchors.
```

### Phase 5: Sparse Particle-MIS Evaluator

Add sparse responsible-chart evaluation with full-mixture audit.

Acceptance criterion:

```text
For sampled theta/local queries, sparse-vs-full local log m_i(theta) differences
stay below tolerance, while outer factor evaluation time decreases.
```

### Phase 6: Production Defaults

Change default EMC hierarchical benchmark settings:

```text
design_method = adaptive_scout
manual focus_hyper_names = optional diagnostic override
all-axis = diagnostic preset only
validation gate = external script only
outer rerun = external script only
```

Acceptance criterion:

```text
A full 50-local run completes with fewer active charts than all-axis and does
not worsen the posterior comparison against EMC.
```

### Phase 7: Speed Optimizations

After the design is correct, optimize evaluation:

1. use Gaussian natural-parameter matrix products for
   `log_alpha_given_theta_many`;
2. cache per-atlas particle arrays, alpha squares, chart log weights, and
   anchor denominators;
3. invalidate caches only when charts are added or normalizers change.

Acceptance criterion:

```text
Full factor-set evaluation time drops without changing log-likelihood values
beyond numerical tolerance.
```

## First Implementation Step

The first step is Phase 1 plus the smallest useful part of Phase 2:

```text
Add a root-scout design planner that builds root charts, computes per-local
score/curvature family risk over theta_cloud, writes a chart_design_plan, and
does not run outer SMC.
```

This is the right first step because it changes the design source of truth. It
does not require tuning the full hierarchical benchmark, and it avoids adding
another patch after a bad atlas is already built.

## Implemented Interfaces

The root-scout stage is exposed as:

```text
build_chart_design_plan()
```

and through:

```text
fit_chart_atlas_population_model(..., design_control = list(
  design_method = "adaptive_scout",
  scout_only = TRUE
))
```

This returns a `chart_design_plan_fit` with:

```text
fit = NULL
factor_set = NULL
atlases = root_atlases
chart_design_plan = chart_design_plan
```

It checkpoints stage:

```text
chart_design_plan
```

The EMC runner is:

```text
benchmarks/run_emc_chart_design_scout.R
```

It writes:

```text
*_results.rds
*_checkpoint_chart_design_plan.rds
*_family_scores.csv
*_local_scores.csv
*_local_family_scores.csv
*_local_anchor_budgets.csv
```

No outer SMC is run in this stage.

The production expansion stage is exposed through the same runner:

```text
benchmarks/run_emc_chart_design_scout.R --stage=atlas
```

and through:

```text
fit_chart_atlas_population_model(..., design_control = list(
  design_method = "adaptive_scout",
  scout_only = FALSE,
  stop_after_atlas_build = TRUE,
  refine_rounds = 0
))
```

This consumes `chart_design_plan$local_theta_designs`, reuses the root scout
atlases, builds only local-specific production charts, writes a frozen
`factor_set`, and checkpoints stage:

```text
post_atlas_build
```

It writes the scout tables plus:

```text
*_atlas_build_history.csv
*_design_certification.csv
*_graph_summary.csv
```

No outer SMC is run in this stage. The next integration step is to run an
explicit outer command from the frozen factor-set checkpoint.

The explicit outer handoff is:

```text
benchmarks/rerun_chart_atlas_outer_from_checkpoint.R \
  --checkpoint_file=..._checkpoint.rds
```

By default this accepts only `post_atlas_build` checkpoints and preserves the
checkpoint's evaluator strictness. The handoff is not a validation gate, but it
is still a certification boundary. It performs these fitting operations before
outer SMC:

```text
normalize the theta proposal once for the population model
run pre-outer proposal-audit repair, if requested
certify the exact initial outer theta cloud
write a post_pre_outer checkpoint
```

The exact initial-cloud certification is mandatory for strict runs. If any
local/theta query remains uncertified after the configured repair rounds, the
script stops before outer SMC and leaves the `post_pre_outer` checkpoint for
inspection. This prevents the outer sampler from starting from an evidence
surface that the local atlases themselves reject.

If repeated rounds fail on the same local/theta pair, treat that as a local
endpoint evidence failure, not an outer-SMC problem. The repair chart either
could not connect through a certified bridge or its direct SMC path was too weak
to certify. The first response should be to increase the certification particle
budget or local SMC quality for that endpoint, not to loosen the final
particle-MIS gate.

During outer MCMC, uncertified proposed theta values are treated as reject-only
points. They must not be evaluated through derivative surfaces or uncertified
particle-MIS estimates. This keeps the current particle cloud strict while
allowing proposals outside the certified region to be rejected without aborting
the run.

For exploratory runs the script can loosen evaluator controls explicitly, for
example:

```text
--stop_on_uncertified=false
--use_uncertified_estimates=true
--min_particle_mis_ess=0
--max_particle_mis_psis_k=Inf
```

The script writes:

```text
*_outer_rerun_results.rds
*_checkpoint_post_outer.rds
*_outer_rerun_posterior_comparison.csv
*_outer_rerun_posteriors.png
```

The posterior comparison and plot are written when `--compare_emc=true`.

The proposal object must be idempotent under normalization. A repeated
`normalize_theta_proposal()` call must not change seeded proposal samples,
otherwise charts added during exact-cloud repair can miss their own later outer
queries by numerical drift.

The one-command EMC fitting entry point is:

```text
benchmarks/run_emc_chart_atlas_gate.R
```

Despite the historical name, this script now uses the same adaptive-scout
framework by default:

```text
design_method = adaptive_scout
calibration_rounds = 0
strict particle-MIS final evaluation
pre-outer audit/repair
exact initial-cloud certification
```

Small `--max_locals` runs are integration checks only. Their posterior
comparison is against the full EMC reference posterior and should not be read as
a scientific accuracy metric unless all locals are included.
