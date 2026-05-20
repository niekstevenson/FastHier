# Bank SMC Plan

This document describes the proposed local SMC bank architecture for amortizing
subject-level marginal likelihoods in hierarchical models.

## 1. Target

For subjects `i = 1, ..., S`, the full hierarchical posterior is

```math
p(theta, alpha_{1:S} | y)
propto
p(theta) prod_i L_i(alpha_i) p(alpha_i | theta),
```

where

```math
L_i(alpha_i) = p(y_i | alpha_i).
```

The population sampler only needs the subject marginal likelihood functions

```math
m_i(theta)
=
int L_i(alpha) p(alpha | theta) d alpha.
```

Then

```math
p(theta | y)
propto
p(theta) prod_i m_i(theta).
```

The local computation should therefore not be framed as "sample one good local
posterior." It should be framed as building a reusable approximation to the
function

```math
theta -> m_i(theta).
```

## 2. Single-Reference Identity

For a reference population value `theta*`, define

```math
pi_i^*(alpha)
=
p(alpha | y_i, theta^*)
=
L_i(alpha) p(alpha | theta^*) / m_i(theta^*).
```

Then, for any other `theta`,

```math
m_i(theta)
=
m_i(theta^*)
E_{pi_i^*}
[
  p(alpha | theta) / p(alpha | theta^*)
].
```

An SMC approximation at `theta*` gives

```math
hat m_i(theta)
=
hat m_i(theta^*)
sum_n w_n^*
  p(alpha_n | theta) / p(alpha_n | theta^*).
```

This estimator is exact in expectation. Its failure mode is variance and
support. It fails when the local posterior induced by `theta` places mass in
alpha-space that is rare under the local posterior induced by `theta*`.

The relevant instability is approximately

```math
int pi_i(alpha | theta)^2 / pi_i^*(alpha) d alpha.
```

When this quantity is large, ESS collapses and the estimator becomes dominated
by rare or missing particles.

## 3. Local SMC Bank

Replace the single reference point with a subject-specific bank:

```math
B_i = {b = 1, ..., B_i}.
```

Each bank node stores

```math
theta_b,
pi_{ib}(alpha) propto L_i(alpha) p(alpha | theta_b),
hat m_i(theta_b),
{alpha_{ibn}, w_{ibn}}.
```

For arbitrary `theta`, evaluate `m_i(theta)` using all bank nodes as a
deterministic mixture proposal:

```math
q_i(alpha)
=
sum_b eta_b pi_{ib}(alpha).
```

Since

```math
pi_{ib}(alpha)
=
L_i(alpha) p(alpha | theta_b) / m_i(theta_b),
```

the proposal density can be evaluated up to stored normalizing constants:

```math
q_i(alpha)
=
sum_b eta_b
  L_i(alpha) p(alpha | theta_b) / hat m_i(theta_b).
```

Then

```math
hat m_i(theta)
=
sum_{b,n}
eta_b w_{ibn}
[
  L_i(alpha_{ibn}) p(alpha_{ibn} | theta)
  / q_i(alpha_{ibn})
].
```

This is the amortized local factor used by the outer population SMC.

## 4. Warm-Started Bank Expansion

Given a candidate population point `theta_new`, first try to reuse the existing
bank.

For each existing bank node `b`, diagnose the single-node transfer

```math
r_b(alpha)
=
p(alpha | theta_new) / p(alpha | theta_b).
```

A useful local distance is

```math
d_i^2(theta_b, theta_new)
approx
Var_{pi_i(. | theta_b)}
[
  log p(alpha | theta_new) - log p(alpha | theta_b)
].
```

Equivalently, compute the implied ESS of the reweighting. Pick the best
existing node:

```math
b^* = argmax_b ESS_b(theta_new).
```

If the ESS is acceptable, do not add a new node.

If the ESS is poor, add a new node by bridging from the nearest existing bank
node to `theta_new`.

Use the bridge

```math
pi_{i,t}(alpha)
propto
L_i(alpha)
p(alpha | theta_{b^*})^{1 - t}
p(alpha | theta_new)^t,
t in [0, 1].
```

The incremental log weight is

```math
Delta log w
=
(t_{r+1} - t_r)
[
  log p(alpha | theta_new)
  - log p(alpha | theta_{b^*})
].
```

The likelihood term does not change during this bridge. Only the population
prior contribution changes. This makes the bridge cheaper and more stable than
starting a fresh local SMC from a broad prior.

At the end, store the endpoint as a new bank node:

```math
theta_new,
hat m_i(theta_new),
{alpha_{in}, w_{in}}.
```

Optionally store selected intermediate bridge states as additional proposal
components if they improve coverage.

## 5. Population Approximation Loop

The true posterior support of `theta` is unknown, because it depends on all
`m_i(theta)`. The bank should therefore be built iteratively.

Use a sequence of approximations

```math
h_0(theta), h_1(theta), ...
```

where `h_r(theta)` approximates the current population posterior support.

Workflow:

1. Build a cheap pilot approximation `h_0(theta)`.
2. Choose initial population design points from `h_0(theta)`.
3. Build initial local banks.
4. Run outer SMC using the bank-based factors.
5. Use the resulting population fit as `h_1(theta)`.
6. Audit bank coverage over `h_1(theta)`.
7. Add bank nodes where local coverage fails.
8. Re-run or reweight the outer population SMC.

The bank should approximate

```math
q_i(alpha)
approx
int pi_i(alpha | theta) h(theta) d theta,
```

but in practice this integral is represented by a finite adaptive mixture over
bank nodes.

## 6. Diagnostics

For each query `theta` and subject `i`, compute mixture importance weights

```math
r_{in}(theta)
=
L_i(alpha_{in}) p(alpha_{in} | theta) / q_i(alpha_{in}).
```

Useful diagnostics:

```math
ESS_i(theta)
=
(sum_n r_{in})^2 / sum_n r_{in}^2,
```

plus:

- maximum normalized weight;
- Pareto `k`;
- variance of log weights;
- local contribution to population-level error.

The global distortion relevant to the outer posterior is

```math
Delta(theta)
=
sum_i [
  log hat m_i(theta) - log m_i(theta)
].
```

Low local ESS alone should not automatically trigger expensive repair. Repair
should prioritize failures that influence posterior mass or posterior-tail
functionals.

## 7. Why This Is Not SMC2

SMC2 would maintain local particles for every population particle:

```math
cost ~ O(S * N_theta * N_alpha).
```

The bank method maintains a small reusable set of local nodes:

```math
cost ~ O(S * B_i * N_alpha) + O(S * B_i * N_theta),
```

where `B_i << N_theta`.

The goal is to keep `B_i` small by adding nodes only where local geometry is
missing.

## 8. Failure Modes

The bank approach can still fail if:

1. The local posterior changes discontinuously with `theta`.
2. The bridge path crosses a severe multimodal region.
3. Local mutation kernels do not mix.
4. Diagnostics are noisy and add unnecessary nodes.
5. The population posterior has many separated regions.
6. The bank grows without pruning or compression.

The implementation must therefore include strict bank budgets, reliable ESS or
Pareto diagnostics, and a strategy for pruning or merging redundant nodes.

## 9. Implementation Steps

### Step 1: Define Bank Node Objects

Create a local bank node representation containing:

- `local_id`;
- `theta_anchor`;
- local particles `alpha`;
- normalized local weights;
- `log_likelihood` at particles;
- `log_prior_anchor = log p(alpha | theta_anchor)`;
- `log_marginal_anchor = log hat m_i(theta_anchor)`;
- local SMC diagnostics;
- bridge provenance.

The node must be self-contained enough to evaluate its contribution to any
future `theta`.

### Step 2: Define Local Bank Objects

Create a subject-level bank object:

- `local_id`;
- list of bank nodes;
- node mixture weights `eta`;
- particle budget;
- node budget;
- cached support diagnostics;
- optional compressed stacked representation.

This object should be the new local factor input to the population layer.

### Step 3: Implement Bank Factor Evaluation

For a subject bank, evaluate

```math
hat m_i(theta)
=
sum_{b,n}
eta_b w_{ibn}
L_i(alpha_{ibn}) p(alpha_{ibn} | theta) / q_i(alpha_{ibn}).
```

Use deterministic-mixture denominator

```math
q_i(alpha)
=
sum_c eta_c
L_i(alpha) p(alpha | theta_c) / hat m_i(theta_c).
```

In log form, for a particle from node `b`,

```math
log numerator(theta)
=
log L_i(alpha)
+
log p(alpha | theta).
```

and

```math
log q_i(alpha)
=
logsumexp_c [
  log eta_c
  + log L_i(alpha)
  + log p(alpha | theta_c)
  - log hat m_i(theta_c)
].
```

The `log L_i(alpha)` term appears in every mixture component and can be
cancelled when appropriate, but store it explicitly first for clarity.

### Step 4: Implement Bank Coverage Diagnostics

For a candidate `theta`, compute:

- single-node ESS from every bank node to `theta`;
- mixture ESS using the full bank proposal;
- Pareto `k` if available;
- max normalized weight;
- log-weight variance.

Return the best existing node and whether coverage is acceptable.

### Step 5: Implement Warm-Started Bridge

Given a local bank, a start node `b*`, and `theta_new`, run a bridge:

```math
pi_t(alpha)
propto
L_i(alpha)
p(alpha | theta_{b*})^{1-t}
p(alpha | theta_new)^t.
```

Use adaptive CESS to choose the next `t`.

Mutation should reuse the existing local SMC move machinery, but with the
bridge log target:

```math
log target_t(alpha)
=
log L_i(alpha)
+
(1 - t) log p(alpha | theta_{b*})
+
t log p(alpha | theta_new).
```

At `t = 1`, store a new bank node at `theta_new`.

### Step 6: Build Initial Population Design

Use a cheap pilot population approximation to pick initial theta points:

- center point;
- weighted support points;
- weak-direction profile points;
- optional stress points.

Do not use a long list of user-facing tuning parameters. Expose only a small
budget, such as `max_bank_nodes` and `target_local_ess`.

### Step 7: Build Initial Banks

For each subject:

1. Create the first node from a broad or pilot-informed reference.
2. Visit initial design points.
3. For each design point, reuse an existing node if ESS is good.
4. Otherwise bridge from the nearest node and store the endpoint.

This replaces planned anchors as the main local geometry construction.

### Step 8: Run Outer Population SMC

Build a population factor set from the local banks and run outer SMC over
`theta`:

```math
p(theta | y)
propto
p(theta) prod_i hat m_i^{bank}(theta).
```

The outer SMC should treat the bank factors as fixed during a run.

### Step 9: Audit Posterior-Relevant Theta Points

After outer SMC, select audit points from:

- posterior weighted support;
- weak-direction profile tails;
- any regions with poor local ESS diagnostics.

For each audit point, estimate or certify

```math
Delta(theta)
=
sum_i [
  log hat m_i(theta) - log m_i(theta)
].
```

Use exact QMC only as a certification tool, not as the main representation.

### Step 10: Adaptive Bank Refinement

If audit fails:

1. Identify influential failing `theta` points.
2. Identify subjects contributing most to the error or low ESS.
3. Add bank nodes only for those subject/theta pairs.
4. Rebuild only changed local factor shards.
5. Reweight or rerun outer SMC.

Stop when posterior-tail functionals and profile audits are stable.

### Step 11: Bank Pruning and Compression

To prevent bloat:

- keep a hard `max_bank_nodes` budget;
- merge nodes with small bridge distance;
- prune nodes with negligible mixture weight and no audit influence;
- keep endpoint nodes preferred over transient bridge nodes unless an
  intermediate node materially improves coverage;
- store compact results by default, with enough bank state to resume.

### Step 12: Extension to Latent Time Processes

If `L_i(alpha)` is itself estimated by an inner state-space SMC, the same bank
architecture still applies. Each bank node stores the local SMC state and
normalizing constant estimate for

```math
L_i(alpha) p(alpha | theta_b).
```

The outer population layer still only sees the bank estimator of

```math
m_i(theta).
```

This keeps the population algorithm separate from the time-dependent local
integration machinery.

## 13. Sufficient-Statistic Bank Evaluation

For exponential-family population priors,

```math
p(alpha | theta)
=
h(alpha)
exp { eta(theta)^T T(alpha) - A(theta) },
```

the outer population layer should not repeatedly evaluate the population prior
from raw bank particles. Each bank particle should store the sufficient
statistics needed by the population family:

```math
T(alpha_n).
```

Then every new population query only needs

```math
log hat m_i(theta)
=
-A(theta)
+
logsumexp_n [
  c_{in} + eta(theta)^T T(alpha_{in})
],
```

where `c_in` contains the deterministic-mixture bank weight and all terms that
do not depend on `theta`.

For diagonal Gaussian random effects,

```math
alpha_j | mu_j, sigma_j^2 ~ Normal(mu_j, sigma_j^2),
```

the sufficient statistics are only

```math
T(alpha) = (alpha_1, ..., alpha_d, alpha_1^2, ..., alpha_d^2).
```

The evaluation becomes

```math
sum_j [
  alpha_j mu_j / sigma_j^2
  - 0.5 alpha_j^2 / sigma_j^2
  - 0.5 log(2 pi sigma_j^2)
  - 0.5 mu_j^2 / sigma_j^2
].
```

For a full-covariance Gaussian prior the statistic is `alpha` and
`alpha alpha^T`; for restricted covariance structures it should be reduced to
the matching lower-dimensional statistic.

The bank should also keep the same bank particles and random numbers fixed
during an outer SMC run. This makes

```math
theta -> prod_i hat m_i(theta)
```

a deterministic smooth approximate likelihood surface. That is preferable to
independent noisy likelihood estimates at each `theta`, because stable
likelihood ratios are what the outer sampler needs.

This is the scalable representation:

```text
local bank construction: stores particles because bridges still need alpha
outer factor evaluation: stores and evaluates only T(alpha)
```

The raw particles remain in the bank state so the computation can be resumed or
expanded. The derived population factor set should use sufficient statistics
and should be rebuildable from the bank state.

## 14. Analytic Population Shape Layer

Population anchors should not be chosen only from the outer posterior induced by
an under-built bank. That creates a circular failure: if the initial bank misses
a tail, the pilot population fit also misses that tail, and the anchor selector
never asks the local banks to represent it.

Before expensive bank expansion, build a cheap independent approximation to the
shape of `theta`.

For diagonal Gaussian random effects, draw plausible subject-level effects from
the current local banks:

```math
alpha_{1:S}^{(k)} ~ prod_i q_i(alpha_i).
```

For each imputed configuration, approximate the conditional group posterior

```math
p(theta | alpha_{1:S}^{(k)}).
```

With independent priors

```math
mu_j ~ Normal(m_{0j}, s_{0j}),
sigma_j^2 ~ InvGamma(a_{0j}, b_{0j}),
```

the conditional group posterior factorizes by dimension. For each dimension,
integrate over `mu_j` analytically to get a one-dimensional density for
`ell_j = log sigma_j^2`:

```math
p(ell_j | alpha_{1:S,j})
propto
p(sigma_j^2) sigma_j^2
int
p(mu_j)
prod_i Normal(alpha_{ij} | mu_j, sigma_j^2)
d mu_j.
```

Then sample or quadrature-select `ell_j`; conditional on `ell_j`, sample

```math
mu_j | ell_j, alpha_{1:S,j}
```

from its exact Gaussian conditional. The resulting theta draws form an analytic
shape cloud. Initial population anchor points are chosen from this cloud, then
scored against current bank coverage.

This is deliberately not the final posterior target. It is a cheap geometric
proposal for anchor placement. The final outer sampler still targets

```math
p(theta) prod_i hat m_i^{bank}(theta).
```

## 15. Natural-Parameter Coverage

Coverage should be measured in the geometry that controls importance-weight
variance, not in raw hyperparameter coordinates.

For an exponential-family population prior,

```math
log p(alpha | theta)
=
eta(theta)^T T(alpha) - A(theta) + log h(alpha),
```

the log transfer from anchor `theta_b` to candidate `theta` is, up to constants,

```math
[eta(theta) - eta(theta_b)]^T T(alpha).
```

The relevant local distance is therefore

```math
d_b^2(theta)
=
Delta eta_b(theta)^T
Var_{bank_b}(T(alpha))
Delta eta_b(theta).
```

For Gaussian random effects with diagonal covariance,

```math
T(alpha) = (alpha, alpha^2),
eta(theta) = (mu / sigma^2, -1 / (2 sigma^2)).
```

The bank implementation computes this scalar variance directly as

```math
Var_{bank_b}( T(alpha)^T Delta eta_b(theta) )
```

without materializing a full covariance matrix. Under a log-normal
approximation to the transfer weights,

```math
ESS_frac approx exp(-d_b^2).
```

So the usual ESS target becomes a natural-variance threshold:

```math
d_b^2 <= -log(target_ess_frac).
```

This natural geometry chooses bridge start nodes and scores candidate population
anchors. Empirical ESS remains a diagnostic, but for Gaussian population priors
the primary coverage decision is natural-parameter variance.
