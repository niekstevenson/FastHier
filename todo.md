# TODO

## Current Status

Implemented:

- reusable local reference fits
- reference-prior pilot and refinement workflow
- population-level factor construction
- outer population SMC
- benchmark scripts for comparison against Stan

Still missing:

- automatic local repair when recycled local factors become unreliable

## Local Repair

The current workflow assumes that each local reference fit has enough overlap with the hierarchical local posterior induced by the outer hyperparameter `theta`.

That is not guaranteed.

For local `i`, the recycled factor is

`m_i(theta) = Z_i^q * E_{pi_i^q}[ p_M(alpha | theta) / q_i(alpha) ]`

and it is only reliable when the stored local particle cloud still covers the `alpha` region that matters under `theta`.

If that overlap is poor, the local reweighting ESS collapses and the recycled factor becomes noisy or unstable.

## What Is Missing

The code currently exposes local ESS diagnostics, but it does not yet use them to trigger repair.

So the remaining gap is:

1. detect when a local factor is no longer trustworthy for a proposed or sampled `theta`
2. mark that local as needing repair
3. rebuild only that local under a better reference prior
4. swap the repaired local back into the factor set

## What Repair Means

Local repair does not mean rerunning the whole hierarchy.

It means rerunning or bridging only the problematic local(s), for example:

- rerun local SMC under a better centered reference prior
- bridge from the current local reference prior to a new one
- maintain more than one reference fit for difficult locals and choose the best one at evaluation time

## Suggested Trigger

Start simple.

For a local `i`, flag repair when one or more of these hold repeatedly over posterior `theta` draws:

- `ESS_i(theta) / N_i` falls below a threshold such as `0.05` or `0.10`
- the maximum normalized reweighting weight is too large
- `log m_i(theta)` is dominated by a very small number of particles

This should be based on repeated failure, not on a single isolated `theta`.

## Minimal Repair Plan

The clean first implementation is:

1. after an outer run, evaluate local ESS on a weighted sample of posterior `theta` particles
2. identify locals with repeated ESS collapse
3. for each bad local, build a repaired reference prior using the problematic posterior region
4. rerun only those locals
5. rebuild the factor set
6. rerun the outer population SMC

This is the main remaining algorithmic gap in the current framework.
