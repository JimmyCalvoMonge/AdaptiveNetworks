# Research proposals

Candidate research directions building on the adaptive-network behavioral epidemic
framework implemented in this repository.

## Active

**Profiles, not parameters** — making risk tolerance (`ν`) and adaptive capacity
endogenous, so that behavioral profiles become dynamic state determined by a node's
circumstances and history rather than parameters assigned by the modeller. Merges two
earlier candidates (constrained capacity, and fatigue/learning) on their shared thesis,
and adds the coupling between them: constrained nodes accumulate *exposure* where
unconstrained nodes accumulate *fatigue*, which predicts that the burden ordering
between the two groups may invert between epidemic waves.

The full proposal is maintained separately as a document. This folder holds the
directions that were considered and set aside.

## Shelved

Each of these was developed to proposal depth and ranked against the others. They are
kept here so the reasoning is not lost, not because they were judged unworkable — most
were set aside because a different direction was chosen first, and several remain
strong candidates for later.

| File | Direction | Changes |
| --- | --- | --- |
| [`imperfect-detection.md`](imperfect-detection.md) | Symptom-based rather than infection-based risk perception, with bilateral edge consent | What the agent can see |
| [`venue-level-decisions.md`](venue-level-decisions.md) | Decisions over venues/groups rather than dyadic edges; higher-order substrate | What the agent decides over |
| [`prosociality-across-health-states.md`](prosociality-across-health-states.md) | Infected and recovered nodes get their own decision problems | Who decides at all |
| [`network-scarring.md`](network-scarring.md) | Dropped ties decay rather than being restored for free | What the choice costs later |
| [`minor-candidates.md`](minor-candidates.md) | Four shorter candidates, with reasons for the lower ranking | — |

## Also considered, earlier

Three directions were developed and then set aside before the shelf above was written.
Summarised here for the record rather than kept at full depth:

- **Identifiability.** Using the cross-scale asynchrony Δ as a statistic to separate
  behavior change from susceptible depletion, which incidence data alone confounds.
  Would need a simulation-based inference layer and an emulator over the simulator.
- **Incentive design.** A planner who perturbs the utility landscape (subsidising the
  cost of contact reduction) under a budget, rather than controlling contacts directly.
- **Vaccination and measles.** Extending the action space to `{edges, vaccinate}`, so
  free-riding and risk compensation emerge rather than being assumed, applied to
  clustered vaccine hesitancy.

## Common ground

Every direction here exploits the same three properties of the implementation:

1. **The decision is modular.** `vertxmaxer` in `SIR_adaptive_net_.py` is a
   self-contained per-node dynamic program. Its utility, information set, action space
   and state space are each replaceable without touching the epidemic layer.
2. **Parameters are preferences, not fitted multipliers.** An agent facing changed
   circumstances re-optimises rather than falling silent, which is what a prescribed
   behavioral rule cannot do.
3. **Δ is a sensitive diagnostic.** The lag between population-level and
   individual-level maximal effort is not robust by construction — it can be widened,
   narrowed or abolished. Every proposal reports Δ under its extension against the
   published baseline (Δ = 9 days at `n=500, p=0.05, ν=0.05, τ=7, β=0.05, γ=0.04`).
