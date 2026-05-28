# EMC Hierarchical Framework Steps

This document describes the framework we should be converging to. It is not a list of benchmark scripts.

The target user contract is:

1. The user supplies local data, a vectorized local likelihood, and local parameter names.
2. The user supplies or accepts a group model. For now the group model is diagonal Gaussian:
   `theta = (mu_1, ..., mu_d, log_sigma2_1, ..., log_sigma2_d)`.
3. The framework builds certified local marginal likelihood surfaces
   `m_i(theta) = integral L_i(alpha) p(alpha | theta) d alpha`.
4. The framework runs outer SMC on the frozen product
   `p(theta) prod_i m_i(theta)`.
5. The framework certifies posterior-relevant local evidence shape, repairs the local surfaces if needed, and reruns or reweights the outer posterior only when justified.

The central invariant is simple: local SMC estimates local evidence; outer SMC consumes frozen local evidence. The outer sampler should not discover or patch local evidence failures mid-run.

## Required User Inputs

| Input | Meaning | Current Representation | Necessary? |
| --- | --- | --- | --- |
| `data_list` | One local dataset per subject/item. | List of local data frames or arrays. | Yes |
| `loglik_fn(Theta, data_i)` | Vectorized local likelihood over rows of local parameter matrix `Theta`. | Function returning one log likelihood per row. | Yes |
| `alpha_names` | Names of local parameters. | Character vector. | Yes |
| `population_model` | Group prior and conditional `p(alpha | theta)`. | Currently `make_population_model_diag_gaussian()`. | Yes |
| `initial_theta_design` or proposal | Initial region of theta space to cover before outer SMC. | Posterior/sketch/q0 cloud or generated design. | Yes |
| controls | Particle counts, chart limits, certification thresholds, parallelism. | Lists. | Yes |
| reference posterior or nested SMC gold | Benchmark-only comparison target. | EMC2/Stan/nested-SMC draws. | No |

## Step 1. Define The Model

What it does:

- Wraps the local likelihood and validates that it evaluates many `alpha` rows for one local dataset.
- Builds the group model object, currently diagonal Gaussian.
- Defines names and dimensions used by all later steps.

Current owner:

- `make_population_model_diag_gaussian()` in `population_models.R`.
- The user-provided `loglik_fn`.

Runtime:

- Seconds.

Absolutely necessary:

- Yes. This is the model specification.

Failure mode:

- If this layer is unclear, every later diagnostic becomes ambiguous because `alpha`, `theta`, and parameter names are not owned by one clean interface.

## Step 2. Define Initial Theta Support

What it does:

- Builds a theta design/proposal region that locals must cover before outer SMC.
- This is where q0 belongs: a design/proposal object, not local evidence.
- The design should include central mass, posterior-like guesses if available, and defensive tails.

Current owners:

- `build_local_evidence_calibration_design()`.
- `build_local_evidence_certification_cloud()`.
- Script-level theta cloud builders in the EMC benchmark should eventually be replaced by one framework API function.

Runtime:

- Seconds unless theta support is built from a pilot run.

Absolutely necessary:

- Yes. Without an initial theta support region, local atlases are forced to extrapolate.

Failure mode:

- Too narrow: posterior tails are under-covered and outer SMC gets biased.
- Too broad: local chart budget is wasted on irrelevant theta regions.

## Step 3. Build Root Local Charts

What it does:

- For every local `i`, runs SMC at a central/root theta.
- The root chart target is
  `pi_i,root(alpha) proportional to L_i(alpha) p(alpha | theta_root)`.
- The root normalizer is an absolute local evidence value and must be explicitly certified.

Current owners:

- `build_local_root_chart()`.
- `build_local_atlas()` starts here when no root atlas is supplied.

Runtime:

- Dominated by `n_locals * root_particles * local_likelihood_cost`.
- For EMC, this is one of the expensive stages.

Absolutely necessary:

- Yes. Every local atlas needs at least one trusted normalizer anchor.

Failure mode:

- A bad root normalizer corrupts the whole local surface for that local.

## Step 4. Expand Local Atlases Over Theta Design

What it does:

- Adds theta-anchored charts for each local where the theta design requires coverage.
- Each chart target is
  `pi_i,a(alpha) proportional to L_i(alpha) p(alpha | theta_a)`.
- Bridge/edge checks are allowed to certify overlap, but adding a chart is not itself evidence certification.

Current owners:

- `build_local_atlas()`.
- `build_local_atlas_from_root()`.
- `build_candidate_chart()`.
- `estimate_chart_edge()` and `certify_chart_edge()`.

Runtime:

- Usually the main local construction cost.
- Scales roughly as `n_locals * n_charts_per_local * particles_per_chart * local_likelihood_cost`.

Absolutely necessary:

- Yes. This is how `m_i(theta)` gets support away from the root theta.

Failure mode:

- Sparse charts in high-curvature theta directions produce locally certified-looking but biased evidence surfaces.

## Step 5. Solve And Certify Local Normalizers

What it does:

- Solves the local chart normalizer graph using certified roots and certified relative edges.
- Marks which charts can act as trusted evidence anchors.
- Keeps support charts separate from certified normalizer anchors.

Current owners:

- `solve_atlas_normalizers()`.
- `atlas_cycle_diagnostics()`.
- `certify_chart_edge()`.

Runtime:

- Seconds compared with SMC, unless the graph is huge.

Absolutely necessary:

- Yes. This is the difference between particles being present and local evidence being trustworthy.

Failure mode:

- The dangerous failure is promoting a single noisy chart as a normalizer-certified evidence anchor.

## Step 6. Build The Local Factor Set

What it does:

- Converts all local atlases into a population factor set.
- The factor set exposes `log m_i(theta)` and `sum_i log m_i(theta)` to the outer sampler.
- Evaluation should use particle-MIS / calibrated chart mixtures as the final non-anchor estimator.

Current owners:

- `build_local_atlas_factor_set()`.
- `population_factor_set_loglik()`.
- `population_factor_set_loglik_by_local()`.

Runtime:

- Seconds to minutes depending on number of theta evaluations.

Absolutely necessary:

- Yes. Outer SMC should only see this factor-set interface, not raw local construction details.

Failure mode:

- If the factor set silently returns uncertified or derivative-surface-only values, the outer posterior can look stable while being wrong.

## Step 7. Optional Particle-MIS Compression

What it does:

- Compresses full local particle mixtures into smaller weighted quadrature rules for faster outer evaluation.
- The normalizers remain separate; compression only replaces the empirical measure used for particle-MIS integration.
- Compression must be certified against held-out theta points.

Current owners:

- `compress_local_atlas_particle_mis()`.
- `compress_local_atlas_factor_set()`.
- `local_atlas_compression_summary()`.

Runtime:

- Seconds to minutes.
- Saves time later if outer/reweight/certification repeatedly evaluates many theta points.

Absolutely necessary:

- No. It is an efficiency layer.

Failure mode:

- Compression can preserve average evidence while damaging posterior-shape directions if holdout certification is too weak.

## Step 8. Pre-Outer Certification

What it does:

- Audits whether the initial outer proposal/theta design can be evaluated by every local factor.
- Repairs unsupported local/theta pairs before outer SMC starts.
- This is still local evidence work, not outer posterior work.

Current owners:

- `local_atlas_pre_outer_certify()`.
- `local_atlas_certify_theta_cloud()`.
- `local_atlas_repair_certification_pairs()`.

Runtime:

- Seconds to minutes if few failures.
- Can become expensive if the initial theta design is poorly matched to the local atlases.

Absolutely necessary:

- Yes for production use. The outer sampler should not start on an uncertified factor set.

Failure mode:

- Repairing too broadly wastes time.
- Repairing with weak normalizer certification creates false confidence.

## Step 9. Run Outer SMC

What it does:

- Samples the group-level posterior:
  `pi(theta) proportional to p(theta) prod_i m_i(theta)`.
- Uses only the frozen local factor set.
- Produces posterior particles, weights, log evidence, and numerical diagnostics.

Current owners:

- `outer_population_smc()`.
- `local_atlas_frozen_outer_rerun()` is a small wrapper for rerunning outer against a frozen factor set.

Runtime:

- Scales roughly as `outer_particles * n_outer_rounds * n_locals * local_factor_eval_cost`.
- With compressed factors this should be much cheaper than local SMC construction.

Absolutely necessary:

- Yes. This is the hierarchical posterior and total model evidence stage.

Failure mode:

- Running more outer particles cannot fix biased local `m_i(theta)` surfaces. It only samples the wrong target more accurately.

## Step 10. Posterior-Region Certification

What it does:

- Takes outer posterior theta particles and checks whether every local factor is certified in the posterior region.
- Reports unsupported, low-ESS, high-PSIS, sparse-neighborhood, or high-risk local/theta pairs.

Current owners:

- `evaluate_raw_local_evidence_certification()`.
- `summarize_raw_local_evidence_certification()`.
- `build_local_evidence_certification_cloud()`.

Runtime:

- Seconds to a few minutes.
- Usually cheaper than SMC repair, but can be costly if every local/theta pair is audited.

Absolutely necessary:

- Yes. This is the main protection against plausible but wrong posteriors.

Failure mode:

- ESS/PSIS support checks alone do not prove normalizer accuracy. They only say the estimator is not obviously unsupported.

## Step 11. Select Local Evidence Probes

What it does:

- Chooses a small batch of local/theta pairs where direct local SMC probes are most informative.
- The selector should target posterior-shape loss, not just the largest local residual.
- It should rank cheaply first, then spend SMC only on the shortlist.

Current owners:

- `select_shape_probe_pairs()`.
- `select_shape_active_probe_pairs()`.
- `build_shape_residual_pool()`.
- `fit_shape_residual_feature_model()`.

Runtime:

- Seconds.

Absolutely necessary:

- Conditional. It is necessary when posterior-region certification shows remaining local evidence risk.

Failure mode:

- A selector that probes every candidate before ranking is backwards and too expensive.
- A selector that only chases pointwise local residuals can miss coherent posterior-shape error.

## Step 12. Run Direct Local Evidence Probes

What it does:

- Runs independent local SMC at selected `(local, theta)` pairs.
- Estimates the residual between the atlas estimate and a direct local evidence estimate.
- This is a diagnostic observation, not automatically a repair.

Current owners:

- `run_shape_probe_pairs()`.

Runtime:

- Scales as `n_selected_pairs * replicates * particles * local_likelihood_cost`.
- This is usually the best place to spend extra budget when the posterior shape is still off.

Absolutely necessary:

- Conditional. It is the cleanest way to detect certified-but-biased local surfaces.

Failure mode:

- Single direct probes can be noisy. High-impact probes need uncertainty estimates or replication.

## Step 13. Repair Local Atlases

What it does:

- Adds charts at selected theta points or selected stencil points.
- Treats new charts as support first.
- Promotes normalizers only after certified edge/direct/graph checks.
- Re-solves local normalizer graphs.

Current owners:

- `repair_shape_selected_probe_pairs()`.
- `repair_shape_residual_geometry()`.
- `local_atlas_repair_certification_pairs()`.

Runtime:

- Seconds to minutes for small selected batches.
- Can become local-SMC dominated if many locals are repaired.

Absolutely necessary:

- Conditional. Necessary only when probes show actionable local evidence error.

Failure mode:

- The old bad pattern was: probe, add chart, trust chart. The correct pattern is: probe, add support, certify normalizer, then trust.

## Step 14. Holdout And Reweight Gate

What it does:

- Tests whether repairs improve held-out local evidence shape.
- Reweights existing outer particles under the repaired factor set.
- Accepts repair if reweight ESS/PSIS and evidence/posterior diagnostics are acceptable.
- Decides whether a fresh outer SMC rerun is needed.

Current owners:

- `validate_shape_repair_holdout()`.
- `shape_repair_outer_reweight_gate()`.

Runtime:

- Holdout: seconds to minutes.
- Reweight gate: often one of the slower non-SMC-local stages because it reevaluates the full factor set across outer particles.

Absolutely necessary:

- Yes after repair. Repairs should not be accepted just because they were added.

Failure mode:

- Posterior overlap can improve while absolute local evidence calibration worsens. The gate must track both shape and evidence.

## Step 15. Iterate Or Stop

What it does:

- If local repairs are accepted and reweighting is stable, save the repaired checkpoint.
- If reweighting is unstable but local repairs are credible, rerun outer SMC from the repaired frozen factor set.
- If holdout or reweight diagnostics worsen, reject the repair and return to probe selection.
- Stop when posterior-region certification, holdout probes, and reweight diagnostics are all stable.

Current owner:

- Currently script-level orchestration in `benchmarks/run_emc_shape_calibration_from_checkpoint.R`.
- This should become a dedicated framework-level orchestration function.

Runtime:

- Depends on whether it reruns outer SMC.

Absolutely necessary:

- Yes as workflow logic. It should not remain scattered across benchmark scripts.

Failure mode:

- Automatic reruns hide cost and make debugging hard. Reruns should happen from checkpoints and be explicit or gate-triggered.

## Benchmark-Only Validation

These are not framework steps. They test the framework.

| Benchmark | What It Tests | Required For Production? |
| --- | --- | --- |
| EMC posterior comparison | Whether the full EMC hierarchy matches EMC2 posterior draws. | No |
| Shifted-gamma hierarchy | Whether the hierarchy behaves on a cheaper synthetic model with reference posterior. | No |
| Single-local EMC nested-SMC comparison | Whether one hard local `m_i(theta)` surface matches a strong nested-SMC reference. | No, but currently the best diagnostic for local evidence bias |

## Implemented API Layer

The framework API layer lives in `hierarchical_framework.R`. Benchmark scripts should load data, define model-specific likelihoods, choose output files, and call these functions. They should not own local-evidence calibration logic.

1. `define_hierarchical_model(data_list, loglik_fn, alpha_names, population_model)`
2. `build_theta_support(model, theta, weights, fit, population_model, ...)`
3. `build_local_evidence_atlases(model, theta_support, theta_design, ...)`
4. `certify_local_evidence(factor_set, theta_support, ...)`
5. `compress_local_evidence(factor_set, theta_support, control, ...)`
6. `run_outer_hierarchy(factor_set, initial_proposal, control, ...)`
7. `certify_posterior_region(factor_set, fit, ...)`
8. `select_local_evidence_probes(model, factor_set, certification, ...)`
9. `run_local_evidence_probes(model, factor_set, selection, ...)`
10. `repair_local_evidence(model, factor_set, probes, ...)`
11. `validate_local_evidence_repair(model, repair, probes, ...)`
12. `gate_repair_and_update_outer(repair, old_fit, old_factor_set, holdout_validation, ...)`
13. `save_hierarchy_checkpoint(state)`

For the current checkpointed EMC posterior-calibration workflow, `run_hierarchy_shape_calibration()` orchestrates steps 7 through 12 and emits stage callbacks so benchmark scripts can save intermediate objects without owning the algorithm.

Remaining cleanup: `fit_chart_atlas_population_model()` still contains too much from-scratch build orchestration. It should eventually be split internally around the same API steps rather than remain a second monolithic path.
