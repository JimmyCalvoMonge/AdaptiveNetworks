# When you cannot tell who is infectious

**Changes:** what the agent can see
**Status:** shelved — strong candidate, and the most direct robustness check on the
published results

## Motivation

The decision in `vertxmaxer` operates on `susinfs`, a list of a node's infectious
neighbours, and drops edges precisely to them. Two assumptions are buried in that:

- **Perfect detection.** A susceptible node knows exactly which of its neighbours are
  infectious.
- **Unilateral action.** Only the susceptible end of an S–I edge ever acts, and S–S
  edges are never dropped at all — even though in a real epidemic people withdraw from
  everyone, not only from the visibly sick.

The SI already concedes half of this: the framework "does not account for the
possibility that the susceptible individual at the opposite end of an edge may
independently choose to drop that same edge." The detectability half is not named, and
it is the larger of the two. Measles is infectious roughly four days before the rash
appears. Presymptomatic and asymptomatic transmission defined COVID-19. For every
pathogen of interest, the model assumes information people do not have.

Chai and Karaliopoulos published the prescribed-rule version of this question in
*Scientific Reports* (2026): with an asymptomatic infectious period, the retarding
effect of link-breaking is "far less effective," because infectious individuals are
mistakenly perceived as safe contacts. That result is individual-based mean-field with
a heuristic breaking rule. The derived-decision version should behave qualitatively
differently, and that difference is the paper.

## Question

When risk is perceived through symptoms rather than infection, does the asynchrony
survive? An agent who knows it cannot identify the source of its risk does not stop
responding — it withdraws *indiscriminately*. Does that blanket response widen Δ,
because it reaches S–S edges that targeted avoidance never touched, or collapse it,
because everyone now behaves alike?

## What to build

- Split I into detectable and undetectable with a detection probability `d`. Perceived
  local prevalence counts only detectable infectious neighbours; the transmission
  calculation still uses all of them.
- Indiscriminate reduction: when perceived risk is nonzero but sources cannot be
  identified, the optimal `e*` applies to the whole neighbourhood and dropped edges are
  sampled from all neighbours. This is the mechanism change, and it brings S–S edges
  into play for the first time.
- Bilateral consent: an edge survives only if both endpoints retain it. The graph
  operation is already symmetric; what is missing is that both endpoints get a decision
  and the drop sets union rather than one side acting alone.
- Sweep `d` from 1 (the current model) to 0 (pure prevalence-driven blanket
  withdrawal), reporting Δ, peak size and final size across the range.

## Why this framework

A prescribed link-breaking rule has nothing to say when the link cannot be identified —
the rule simply fails to fire. A utility-maximising agent re-optimises over a coarser
information set and arrives somewhere new. This is the sharpest available demonstration
that deriving behavior is more than an aesthetic preference, and it is also the honest
robustness check the current results need.

## Risk

At low `d` the model may collapse toward uniform contact reduction, which is close to
mean-field distancing and would wash out the network-specific result. That is
informative rather than fatal — it would locate the detectability threshold below which
network structure stops mattering — but it should be anticipated in the design, not
discovered in review.

## Practicalities

- **New machinery:** detectable/undetectable split; indiscriminate drop sampling;
  bilateral consent
- **External data:** none; detection rates from existing serosurvey literature
- **Compute:** low
- **Venue:** PLOS Computational Biology · J. R. Soc. Interface · Physical Review E

## References

- Chai & Karaliopoulos. Epidemic spread with asymptomatic infectious period in contact
  adaptive networks. *Scientific Reports* (2026).
  https://www.nature.com/articles/s41598-026-36212-y
