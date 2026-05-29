# Panel-Calibrated Nested Atlas Workflow

## Purpose

This document describes a replacement workflow for scalable hierarchical evidence estimation.

The goal is to approximate the full hierarchical posterior and marginal likelihood while avoiding the cost of full nested SMC at every outer particle. The method builds reusable local evidence surfaces, validates those surfaces on posterior-relevant theta panels, and only then runs the final outer SMC.

This is not a post-hoc calibration layer for the existing EMC benchmark. It is a new main workflow. Some existing code primitives remain useful, but the orchestration changes.

## Problem Setting

We have many independent local data sets indexed by `i = 1, ..., n`.

For each local:

```text
data_i
local parameter alpha_i
local likelihood L_i(alpha_i)
```

The population model is:

```text
alpha_i | theta ~ p(alpha | theta)
theta ~ p(theta)
```

The local marginal likelihood contribution is:

```text
m_i(theta) = integral L_i(alpha) p(alpha | theta) d alpha
```

The full posterior over theta is:

```text
pi(theta | data) proportional to p(theta) product_i m_i(theta)
```

The full marginal likelihood is:

```text
Z = integral p(theta) product_i m_i(theta) d theta
```

The hard computational problem is evaluating `m_i(theta)` accurately for many locals and many theta values.

## Reference Method: Full Nested SMC

The brute-force nested SMC reference is:

```text
for each outer theta particle:
  for each local i:
    run a fresh local SMC targeting L_i(alpha) p(alpha | theta)
    estimate m_i(theta)
  combine all local factors in the outer SMC
```

This is accurate in principle, but too expensive because local SMC is repeated for every outer theta particle.

The proposed workflow is a reusable approximation to nested SMC:

```text
run selected expensive local SMCs
store them as reusable local atlas charts
evaluate new theta by particle-MIS over certified charts
validate the summed evidence error on theta panels
run final outer only after local surfaces are certified
```

## Central Quantity

The local atlas approximates:

```text
ell_i(theta) = log m_i(theta)
```

Let the atlas estimate be:

```text
ellhat_i(theta)
```

The local error is:

```text
e_i(theta) = ell_i(theta) - ellhat_i(theta)
```

The total log-evidence-surface error seen by the outer posterior is:

```text
Delta(theta) = sum_i e_i(theta)
```

A constant shift in `Delta(theta)` changes the marginal likelihood estimate, but mostly does not change posterior shape. Posterior shape is driven by the centered field:

```text
Delta_c(theta) = Delta(theta) - E_pi[Delta(theta)]
```

The main validation target is therefore not local ESS, not PSIS alone, and not pointwise local residuals. The main target is held-out posterior-weighted error in `Delta_c(theta)`.

## Empirical Motivation

The complete EMC theta-panel diagnostic found large posterior-shaping local evidence error:

```text
pre holdout Delta_c RMSE: 4.288
pre repair  Delta_c RMSE: 3.981
max centered Delta_c:     7.735
```

This means the local evidence surface error is large enough to distort the posterior.

The current sparse repair operator barely changed this:

```text
holdout Delta_c RMSE: 4.288 -> 4.256
repair  Delta_c RMSE: 3.981 -> 4.029
```

That diagnostic supports the new architecture:

```text
local evidence surfaces need panel-calibrated chart promotion,
not post-hoc residual correction and not sparse scattered repair.
```

## Architectural Principle

The workflow has one central rule:

```text
Only a frozen, panel-validated local factor set may feed the final outer SMC.
```

Local atlases may be adapted between outer runs. They must not be mutated during an outer SMC run.

Direct local SMC is expensive, so when it is run for a panel audit it must produce reusable chart candidates. A direct local SMC run should not be discarded after only recording a scalar log marginal likelihood.

## Core Data Objects

### HierarchicalModel

Owns the statistical model definition.

Required fields:

```text
data_list
loglik_fn(alpha_matrix, data_i)
alpha_names
population_model
local_names
```

The population model owns:

```text
p(alpha | theta)
log p(alpha | theta)
p(theta)
theta names
alpha names
theta transformations if needed
```

For now, the population model can be diagonal Gaussian. The workflow should not hard-code EMC-specific assumptions.

### ThetaDesign

Initial theta points used to seed local atlases before any final posterior is known.

Contains:

```text
theta matrix
theta source labels
theta weights if available
design diagnostics
```

The design should be conservative and defensive. It is not evidence and not final posterior inference.

### LocalChart

A reusable local SMC state anchored at one theta value.

Contains:

```text
local_id
theta_anchor
alpha_particles
particle_weights
logZ
logZ_uncertainty
SMC path diagnostics
normalizer certification status
normalizer certification diagnostics
```

Mathematical target:

```text
pi_i,a(alpha) proportional to L_i(alpha) p(alpha | theta_a)
Z_i,a = m_i(theta_a)
```

### LocalAtlas

The collection of charts for one local.

Contains:

```text
local_id
charts
normalizer graph
chart status
particle-MIS evaluator configuration
diagnostics
```

### LocalAtlasFactorSet

The collection of local atlases for all locals.

This is the object consumed by outer SMC.

Contains:

```text
population_model
list of LocalAtlas objects
evaluator control
compression state if present
validation state
```

### PilotOuterFit

An outer SMC fit run against an initial or provisional factor set.

This is a design object, not final inference.

Contains:

```text
theta particles
weights
log likelihood/factor values
outer SMC diagnostics
```

### ThetaPanel

A small set of posterior-relevant theta points used to audit local evidence surfaces.

It has at least two subsets:

```text
training theta:
  used for diagnosis and chart promotion

heldout theta:
  used only for validation
```

Contains:

```text
theta matrix
theta weights
theta role: train or heldout
theta source
selection diagnostics
```

### PanelProbeSet

Direct local SMC results for a theta panel.

For each `(local_i, theta_j)` pair, it contains:

```text
ProbeSummary:
  logZ
  logZ uncertainty
  replicate diagnostics
  SMC path diagnostics

ChartCandidate:
  alpha_particles
  particle_weights
  theta_anchor
  logZ
  logZ_uncertainty
  local_id
```

The `ChartCandidate` is mandatory for repair-grade probes. Summary-only probes are acceptable only for cheap diagnostics.

### PanelDeltaDiagnostics

Computes the local and total evidence errors:

```text
e_i(theta)
Delta(theta)
Delta_c(theta)
```

Contains:

```text
local residual table
theta total error table
local contribution table
train summary
heldout summary
certified-but-biased local counts
uncertified local counts
replicate uncertainty summary
```

### PanelRepairPlan

Specifies which chart candidates to promote and which additional local theta anchors to run.

Contains:

```text
selected local ids
selected theta rows
promotion candidates
replication requests
extra local panel requests
reason codes
expected impact
```

### PanelValidation

Evaluates whether the repaired factor set improved held-out posterior-shaping evidence error.

Contains:

```text
pre and post heldout Delta_c metrics
pre and post train Delta_c metrics
certification changes
PMIS/PSIS changes
outer reweight diagnostics if run
accept/reject decision
```

## Full Workflow

### Stage 0: Define The Model

The user defines local likelihoods and the population model.

Inputs:

```text
data_list
loglik_fn
alpha_names
population_model
```

Output:

```text
HierarchicalModel
```

Failure modes:

```text
wrong alpha/theta naming
invalid local likelihood
population model does not align with alpha names
```

Stage 0 should be boring and explicit. It should not contain SMC logic.

### Stage 1: Build Initial Theta Design

Construct a defensive theta design before local atlas construction.

Sources may include:

```text
hyperprior draws
cheap local sketches
Laplace or VI approximations if available
small pilot samples
tail inflation points
```

Output:

```text
ThetaDesign
```

Purpose:

```text
give local atlases initial support
seed roots and early chart placement
avoid starting entirely from a biased outer posterior
```

This stage should not try to solve the posterior. It only builds a useful initial design.

### Stage 2: Build Initial Local Atlases

For each local, run local SMC at selected initial theta anchors.

Each initial chart targets:

```text
L_i(alpha) p(alpha | theta_a)
```

The local SMC result must store:

```text
particles
weights
logZ
uncertainty
diagnostics
```

Bridging can be used between nearby anchors if overlap diagnostics pass. Fresh SMC is the default when overlap is uncertain.

Output:

```text
LocalAtlasFactorSet
```

Failure modes:

```text
root chart unstable
normalizer uncertainty too high
chart graph disconnected
support too narrow
```

### Stage 3: Run Pilot Outer SMC

Run outer SMC against the initial local factor set.

Output:

```text
PilotOuterFit
```

Purpose:

```text
discover posterior-relevant theta regions
provide weights for panel design
expose obvious support gaps
```

This is not the final posterior unless the later panel validation passes without needing repair.

### Stage 4: Design Theta Panels

Select posterior-relevant theta points for complete local evidence auditing.

Training panel:

```text
used to diagnose local evidence error
used to decide chart promotions
```

Held-out panel:

```text
never used to choose repair
used only to validate whether repair generalized
```

Production selection uses internal signals:

```text
posterior mass
posterior tails
low particle-MIS ESS
high PSIS k
leave-chart-out fragility
large local graph residuals
high normalizer uncertainty
outer reweight sensitivity
uncertified local/theta queries
```

Benchmark selection may additionally use reference mismatch axes to stress the framework.

Output:

```text
ThetaPanel
```

Failure modes:

```text
panel too narrow
heldout points too similar to training points
panel ignores posterior tails
panel chosen from biased posterior only
```

### Stage 5: Run Complete Panel Direct SMC

For each selected theta point, run direct local SMC for all locals.

For the first production-quality implementation, prefer all locals over a subset. A subset is only acceptable if it has explicit uncertainty accounting for the missing locals.

Output:

```text
PanelProbeSet
```

Each probe must produce both:

```text
logZ summary
reusable chart candidate
```

Replicated probes:

```text
combine normalizers on the Z scale
estimate empirical uncertainty
flag disagreement
```

Failure modes:

```text
direct SMC itself too noisy
probe result lacks particles
replicates disagree
SMC path degeneracy
```

If probes are too noisy, the workflow cannot decide whether the local atlas is wrong. The correct response is more reliable direct probes, not weaker validation.

### Stage 6: Diagnose Panel Free-Energy Error

Evaluate the current local atlas at every panel theta and compare to direct SMC.

Compute:

```text
e_i(theta) = ell_i^direct(theta) - ell_i^atlas(theta)
Delta(theta) = sum_i e_i(theta)
Delta_c(theta) = Delta(theta) - E_pi Delta(theta)
```

Primary summaries:

```text
train Delta_c RMSE
heldout Delta_c RMSE
max abs Delta_c
local contribution to Delta_c
certified-but-biased local count
uncertified local count
replicate uncertainty
```

A local value is certified-but-biased when:

```text
atlas reports certified
direct SMC disagrees materially
disagreement is larger than uncertainty
```

Output:

```text
PanelDeltaDiagnostics
```

Interpretation:

```text
large Delta_c:
  local evidence surface can explain posterior shape error

small Delta_c but wrong posterior:
  investigate outer SMC, model mismatch, theta prior, or reference comparison

large local errors but small Delta_c:
  errors mostly cancel or are posterior-constant
```

### Stage 7: Plan Repairs By Contribution

Do not select repairs by largest pointwise residual alone.

Rank locals by contribution to posterior-shaping error:

```text
C_i = Cov_theta(e_i(theta), Delta_c(theta))
```

Also classify each local's error:

```text
mostly constant bias:
  affects evidence more than posterior shape

coherent shape bias:
  affects posterior shape and needs multi-theta repair

isolated support failure:
  needs exact chart/support repair

noisy direct probe:
  needs replication before promotion

graph conflict:
  needs normalizer graph repair
```

Output:

```text
PanelRepairPlan
```

Failure modes:

```text
planner selects too few high-contribution locals
planner chases noisy probes
planner repairs constant shifts while ignoring shape
planner proposes stencils without value anchors
```

### Stage 8: Promote Direct Probes Into Charts

This is the key replacement for current sparse repair.

Promotion uses the direct probe's own particles and weights. It does not run a new small SMC unless the original probe did not save a usable chart candidate.

Certification requirements:

```text
replicate logZ agreement on Z scale
acceptable SMC path diagnostics
normalizer uncertainty below threshold
edge or graph consistency with nearby charts
leave-chart-out stability
bounded shift to existing certified values unless strongly justified
```

Promotion result:

```text
new active chart
support-only chart
quarantined chart
rejected candidate
```

Output:

```text
UpdatedLocalAtlasFactorSet
```

Failure modes:

```text
promoted chart normalizer not reliable
promoted chart shifts graph incoherently
too few charts promoted for a high-contribution local
chart promotion improves support but not value accuracy
```

### Stage 9: Add Local Multi-Theta Repairs If Needed

If one local contributes coherent shape error across theta, promote several theta anchors for that local.

One exact chart is usually insufficient for curvature in:

```text
ell_i(theta)
```

Multi-theta local repair can include:

```text
panel theta anchors already probed
nearby tail anchors
small directional stencils
additional replicated probes
```

Derivative information may guide where to add anchors, but value anchors certify the surface.

Output:

```text
UpdatedLocalAtlasFactorSet
```

### Stage 10: Validate Held-Out Delta_c

Re-evaluate the updated factor set against the held-out panel.

Acceptance criteria:

```text
heldout Delta_c RMSE decreases
max abs heldout Delta_c decreases or stays within strict tolerance
certified-but-biased local count decreases
uncertified local count does not materially increase
PMIS/PSIS diagnostics remain acceptable
normalizer graph diagnostics remain acceptable
```

Training-only improvement is not enough.

Output:

```text
PanelValidation
```

Interpretation:

```text
train improves, heldout improves:
  repair generalized

train improves, heldout fails:
  repair is too local or overfit

train fails, heldout fails, Delta_c large:
  repair operator is inadequate

Delta_c small, posterior wrong:
  local evidence is not the main failure
```

### Stage 11: Outer Update

Only after panel validation passes, update outer inference.

First run outer reweighting:

```text
old outer particles
new frozen factor set
importance ratio = new local factor product / old local factor product
```

If reweight ESS is high and PSIS is acceptable:

```text
use reweighted fit as a diagnostic or accepted approximation
```

If reweight ESS is poor but panel validation is strong:

```text
rerun outer SMC from the frozen updated factor set
```

If the new outer posterior moves into theta regions not covered by panels:

```text
return to theta panel design
```

Output:

```text
OuterFit
```

### Stage 12: Compression

Compression happens after local validation, not before.

Compression must preserve local particle-MIS estimates on:

```text
training panel theta
heldout panel theta
posterior theta samples
tail theta samples
```

Compression output is accepted only if compressed estimates match full atlas estimates within the local evidence uncertainty budget.

Output:

```text
CompressedLocalAtlasFactorSet
```

### Stage 13: Final Evidence Report

Final evidence is reported only from:

```text
frozen local factor set
validated local panel diagnostics
outer SMC run against that factor set
```

Report:

```text
outer log evidence
outer SMC diagnostics
local panel validation summary
heldout Delta_c metrics
remaining certified-but-biased local count
compression diagnostics if used
runtime and local SMC budget
```

## Main API

The workflow should expose explicit steps:

```r
model <- define_hierarchical_model(...)

theta_design <- design_initial_theta(
  model,
  control
)

factor_set <- build_initial_local_atlases(
  model,
  theta_design,
  control
)

pilot <- run_pilot_outer(
  factor_set,
  control
)

panel <- design_theta_panel(
  model,
  factor_set,
  pilot,
  control
)

panel_probes <- run_panel_local_smc(
  model,
  factor_set,
  panel,
  control
)

diagnostics <- diagnose_panel_delta(
  factor_set,
  panel_probes
)

repair_plan <- plan_panel_chart_promotions(
  diagnostics,
  panel_probes,
  factor_set,
  control
)

factor_set <- promote_panel_charts(
  factor_set,
  panel_probes,
  repair_plan,
  control
)

validation <- validate_panel_delta(
  factor_set,
  panel_probes,
  control
)

fit <- run_validated_outer(
  factor_set,
  pilot,
  validation,
  control
)
```

A convenience orchestration function can exist:

```r
fit <- run_panel_calibrated_nested_atlas(model, control)
```

but it should only call the explicit steps. It should not hide local repair or mutate local atlases during outer SMC.

## Existing Code Mapping

The new workflow reuses primitives, not the old workflow.

| Existing piece | New role |
|---|---|
| local chart object | Keep as chart primitive |
| particle-MIS evaluator | Keep as final non-anchor evaluator |
| local normalizer graph | Keep as consistency diagnostic |
| relative edge checks | Keep as chart connection diagnostics |
| PMIS ESS / PSIS | Keep as support and tail diagnostics |
| direct local SMC replicates | Extend to return chart candidates |
| sparse shape probe selector | Not main path |
| residual geometry repair | Not main path |
| residual correction | Removed from main path |
| outer SMC | Keep, but final only after panel validation |
| compression | Keep, but only after panel validation |

## Iteration Policy

Allowed iteration:

```text
pilot outer
theta panel design
complete panel probes
chart promotion
heldout panel validation
outer update
```

If outer update changes posterior support:

```text
new outer support -> new panel -> more probes -> more chart promotion -> outer update
```

Forbidden in the main workflow:

```text
mutating local atlases during outer SMC
applying residual corrections to factor set loglik
declaring success from ESS/PSIS while heldout Delta_c worsens
discarding direct SMC particles when the run was intended for repair
reporting final evidence from an unvalidated factor set
```

## Why This Is The Right Level Of Nestedness

Full nested SMC pays:

```text
outer theta particles * locals * local SMC cost
```

The panel-calibrated atlas pays:

```text
posterior theta panels * locals * local SMC cost
```

and then reuses those local SMC states through particle-MIS.

This does not pretend local SMC is cheap. It spends local SMC where it certifies the posterior-relevant local evidence surface.

## Implementation Order

### Phase 1: Create The New Workflow Module

Create:

```text
nested_atlas_workflow.R
```

This module owns the new workflow. It may call existing low-level primitives, but it should not route through the old sparse shape-calibration benchmark.

### Phase 2: Make Direct SMC Probes Reusable

Modify direct local SMC execution so repair-grade probes return:

```text
particles
weights
logZ
uncertainty
replicate diagnostics
theta anchor
local id
```

This is the first essential code change. Without it, the new workflow collapses back into the old pattern: diagnose correctly, then repair indirectly.

### Phase 3: Implement Panel Data Objects

Add constructors and validators for:

```text
ThetaPanel
PanelProbeSet
PanelDeltaDiagnostics
PanelRepairPlan
PanelValidation
```

Keep them explicit and small.

### Phase 4: Implement Panel Delta Diagnostics

Compute:

```text
e_i(theta)
Delta(theta)
Delta_c(theta)
local contribution scores
```

These diagnostics replace sparse residual summaries as the main local evidence diagnostic.

### Phase 5: Implement Contribution-Based Repair Planning

Rank local repairs by contribution to `Delta_c`.

Separate:

```text
constant local bias
posterior-shape local bias
support failure
normalizer conflict
noisy probe
```

### Phase 6: Implement Chart Promotion

Promote reusable direct probe states into local atlases.

Combine replicated normalizers on the `Z` scale and attach uncertainty.

### Phase 7: Implement Held-Out Panel Gate

Accept repair only if held-out `Delta_c` improves or stays inside strict tolerance.

This is the main protection against changing the posterior without improving the evidence surface.

### Phase 8: Implement Validated Outer Update

Run outer reweighting or outer SMC only after local panel validation passes.

### Phase 9: Integrate Compression After Validation

Compress only validated local atlases. Validate compressed estimates against full atlas estimates on panel and posterior theta points.

## First Concrete Step

The first implementation step is:

```text
make direct local SMC panel probes return reusable chart candidates
```

Reason:

```text
the panel diagnostic already proves Delta_c is large,
but repair cannot use the expensive direct SMC states.
```

Until direct probes are promotable, any new panel workflow will remain diagnostic rather than corrective.
