# Framework Steps And Runtime

This is a framework audit, not a progress summary. It separates the core hierarchy from validation and exploratory repair code.

## Active Benchmark Entry Points

Only these scripts should remain active under `benchmarks/`:

| Script | Purpose | Status |
| --- | --- | --- |
| `benchmarks/run_emc_shape_calibration_from_checkpoint.R` | Current hierarchical EMC ELP benchmark from a saved post-outer atlas checkpoint. Runs raw certification, shape probe selection, direct probes, strict repair, holdout validation, reweight gate, optional frozen outer rerun, and posterior plots. | Active EMC hierarchy entrypoint. |
| `benchmarks/run_shifted_gamma_hierarchy_compare.R` | Current hierarchical shifted-gamma benchmark using the bank-SMC population framework and Stan posterior comparison. | Active synthetic hierarchy entrypoint. |
| `benchmarks/run_emc_single_local_mi_quality.R` | Single hard EMC local comparison against nested-SMC reference values for `m_i(theta)`, with cached gold references and diagnostic plots. | Active local evidence validation entrypoint. |

Everything else in `benchmarks/` was an experiment, diagnostic, old pipeline, or helper around a superseded pipeline and should live under `old/benchmarks/`.

## EMC ELP Hierarchy From Checkpoint

This is the current usable hierarchical EMC ELP workflow. It starts from a saved post-atlas/post-outer checkpoint rather than rebuilding the initial atlas from raw data.

| Step | What It Does | Rough Time | Absolutely Necessary? | Notes |
| --- | --- | ---: | --- | --- |
| 0. Load inputs | Loads post-outer checkpoint, EMC data, baseline draws if present, and reconstructs likelihood/model handles. | Seconds | Yes | If the checkpoint is missing, this script is not a from-scratch replacement. |
| 1. Theta audit cloud | Selects posterior-relevant theta rows plus tail points from the current outer particles. | Seconds | Yes | This defines where the frozen local surfaces will be challenged. |
| 2. Raw PMIS certification | Evaluates local `m_i(theta)` support/certification for every selected local/theta query using the raw atlas evaluator. | Seconds to minutes | Yes for calibration | This is the core check that prevents unsupported extrapolation. |
| 3. Shape probe selection | Ranks local/theta pairs for direct evidence probing. Active mode uses previous residual probes and current certification diagnostics. | Seconds | Conditional | Needed for adaptive calibration. Not a mathematical part of the final posterior. |
| 4. Direct shape probes | Runs direct local SMC at selected local/theta pairs to observe local evidence residuals. | Seconds to minutes | Conditional | This is expensive enough to matter, but far cheaper than probing every pair. |
| 5. Strict shape repair | Adds charts only at selected probed locations and certifies them through the local chart machinery. | Seconds to minutes | Conditional | Needed only if probes identify actionable local surface error. |
| 6. Holdout validation | Runs held-out direct probes and checks whether repair improved shape error. | Seconds to minutes | Yes if repair is accepted | This is the guard against repairing one region while damaging another. |
| 7. Outer reweight gate | Reweights existing outer particles under the repaired factor set and checks ESS/PSIS/evidence impact. | Often 1 to 2 minutes | Yes if repair is accepted | Recent runs showed this can dominate runtime. It is cheaper than a fresh outer rerun and should remain mandatory before rerun. |
| 8. Frozen outer rerun | Optionally reruns outer SMC against the repaired frozen local factor set. | Minutes or worse | No by default | This should be explicit or gate-triggered, not automatic fitting overhead. |
| 9. Posterior comparison/plot | Compares workflow posterior draws to EMC2 posterior draws and writes plots/CSV. | Seconds | Benchmark only | Not part of production inference. |

Recent checkpointed EMC calibration runs spent most wall time in the reweight gate. Direct probes and repairs were relatively small once the selector stopped probing every local/theta pair.

## Shifted-Gamma Hierarchy

This is the synthetic hierarchy benchmark. It is useful because the likelihood is cheap enough to stress the population workflow without EMC likelihood cost.

| Step | What It Does | Rough Time | Absolutely Necessary? | Notes |
| --- | --- | ---: | --- | --- |
| 0. Load Stan benchmark | Loads synthetic data, priors, and reference Stan posterior. | Seconds | Yes | The benchmark is meaningless without the reference. |
| 1. Define local likelihood and population model | Builds shifted-gamma local log likelihood and diagonal Gaussian population prior. | Seconds | Yes | This is the model specification. |
| 2. Build local banks | Runs local SMC banks for each subject, with theta anchors, bridge checks, and local evidence surfaces. | Minutes depending on `subjects * max_bank_nodes * local_particles` | Yes | This is the main local-evidence cost. |
| 3. Design and audit refinement | Adds or audits theta profiles when local bank ESS/support diagnostics fail. | Minutes | Yes for the updated benchmark | This is where the framework checks whether local `m_i(theta)` surfaces are usable over posterior-relevant theta. |
| 4. Outer SMC | Samples the group-level posterior using frozen local evidence factors. | Seconds to minutes | Yes | This is the actual hierarchical posterior. |
| 5. Stan comparison and posterior plot | Compares workflow posterior to Stan and writes plot/CSV. | Seconds | Benchmark only | Useful for regression checking, not part of inference. |

The shifted-gamma script is not a local evidence gold-standard script. Its purpose is end-to-end hierarchy behavior against a known posterior.

## Single-Local EMC Gold Comparison

This script isolates the most important mathematical object: the local marginal evidence surface

`m_i(theta) = integral L_i(alpha) p(alpha | theta) d alpha`.

| Step | What It Does | Rough Time | Absolutely Necessary? | Notes |
| --- | --- | ---: | --- | --- |
| 0. Load EMC data | Loads EMC2 ELP data and selects one local, either by index or `local_pos=auto`. | Seconds | Yes | This avoids full 50-local hierarchy cost. |
| 1. Build theta designs | Chooses atlas anchor theta points and audit theta points from the EMC posterior cloud. | Seconds | Yes | The design determines what part of `m_i(theta)` is tested. |
| 2. Build local atlas | Runs chart SMC at atlas anchors and calibrates the chart bank. | Seconds to minutes | Yes | This is the approximation being tested. |
| 3. Evaluate atlas on audit points | Evaluates `m_i(theta)` through the particle-MIS/chart evaluator. | Seconds | Yes | This is the framework estimate. |
| 4. Nested-SMC reference | Runs independent nested SMC at audit theta points, optionally cached. | Minutes to much longer | Yes for gold comparison | This is the expensive reference. Cache reuse is essential for iteration. |
| 5. Error diagnostics | Computes local log marginal errors, centered errors, ESS/PSIS/certification diagnostics, and plots. | Seconds | Yes | This is the cleanest place to detect local evidence bias before full hierarchy. |

This benchmark is the right place to test compression, chart placement, local repair, and normalizer certification. It should not depend on the full EMC posterior benchmark except for using EMC posterior draws as a theta design cloud.

## Core Versus Non-Core Machinery

Core framework pieces:

- local SMC chart/bank construction;
- explicit local normalizer handling;
- raw particle-MIS local evidence evaluation;
- strict certification before trusting local `m_i(theta)`;
- outer SMC over frozen certified local factors;
- checkpointed factor sets and fits.

Conditional calibration pieces:

- active residual probe selection;
- direct local SMC probes;
- selected exact chart repair;
- holdout validation;
- outer reweight gate.

Benchmark-only pieces:

- comparison to EMC2 posterior draws;
- comparison to Stan posterior draws;
- nested-SMC gold references for selected locals;
- plotting scripts and CSV reports.

Not core and should not drive final posterior directly:

- derivative surfaces as final evaluators;
- broad residual correction as evidence;
- automatic full outer reruns after every repair;
- probing every local/theta pair before selection;
- old query-native visual scripts and one-off certification stress tests.

## Current Weak Point

The framework still depends on calibration quality of local evidence surfaces. The cleanest test is not the full posterior plot; it is the single-local EMC gold comparison. If a local atlas cannot match nested-SMC `m_i(theta)` on posterior-relevant theta points, the full hierarchy will inherit that shape error.

The active EMC hierarchy script should therefore be treated as a consumer of certified factor sets, not as the place where every local evidence bug gets hidden.
