# What the epidemic leaves behind

**Changes:** what the choice costs later
**Status:** shelved — most original output; shares tie-strength machinery with the
active proposal's capacity channel

## Motivation

Every dropped edge in the current model returns for free the next day —
`rednet = net.copy()` restores the baseline at every step. That is the right first
approximation, and it is exactly what makes "temporary reduction" a different object
from rewiring. It also means the epidemic leaves no trace on the social fabric, and
empirically it does.

German panel data spanning the pandemic found around a third of respondents lost
contact with acquaintances and one in four lost contact with a friend, with a
measurable restructuring toward kin ties and a sharp decline in friendships. Older
adults' close networks churned substantially across the same period.

No adaptive-network epidemic model captures this, because they all either restore edges
completely or rewire them under a conserved degree. This framework is one parameter
away from being the first that can.

## Question

When unused ties decay, what network does the epidemic leave behind — and does the
scarred network make the *next* outbreak worse? Does the loss concentrate where the
panel data says it does, in weak non-kin ties held by high-degree nodes?

## What to build

- A per-edge reinstatement probability that declines with the number of consecutive
  days the edge has been dropped, so a briefly suppressed tie returns and a
  long-suppressed one may not. One parameter.
- Tie strength as an edge property — strong (kin, high utility, high reinstatement)
  versus weak (acquaintance, low utility, low reinstatement). This is what makes the
  decision's *choice* of edge matter, which today it does not: dropped edges are
  sampled uniformly at random from the infectious set.
- Post-epidemic network metrics against baseline — degree distribution, clustering,
  component structure, and specifically kin/non-kin composition, tested against the
  direction of the empirical restructuring.
- Sequential outbreaks on the scarred network: the second epidemic's peak and final size
  relative to the first.

## Why this framework

The simulation already computes, per node per day, exactly which edges were dropped and
why — and then discards it. Decay is a read on data that already exists. More
importantly, once ties differ in strength, the currently arbitrary implementation
detail of *which* infectious edge gets dropped becomes a modelled behavior with
consequences that outlast the epidemic.

## Risk

Two failure modes. The reinstatement probability has no direct empirical anchor, so
calibrate the *aggregate* tie-loss fraction against the panel estimates rather than
trying to justify a per-edge rate. And the second-outbreak effect is only interesting
if the scarring is large enough to matter; if it is not, report the null and keep the
restructuring result, which stands on its own.

## Practicalities

- **New machinery:** edge memory and reinstatement; tie-strength attribute; sequential
  outbreaks
- **External data:** panel estimates of pandemic tie loss, for aggregate calibration
- **Compute:** moderate
- **Venue:** Nature Human Behaviour · PNAS · Social Networks

## References

- The Making and Breaking of Social Ties During the Pandemic. *Frontiers in Sociology*
  (2022).
  https://www.frontiersin.org/journals/sociology/articles/10.3389/fsoc.2022.837968/full
- Social ties in old age: the effect of the COVID-19 pandemic.
  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12586840/
