# People skip venues, not friends

**Changes:** what the agent decides over
**Status:** shelved — best structural fit to the Manassas digital twin

## Motivation

The current decision drops individual edges to infectious neighbours. Real avoidance is
coarser and lumpier: you skip the restaurant, cancel the choir practice, stop going to
the office. You do not selectively drop three of your seven colleagues.

This matters because the digital twin is *already* built venue-first — synthetic
individuals are assigned activity sequences, mapped to locations, and edges are derived
from co-presence. The dyadic representation throws that structure away immediately
after using it.

The timing is good. Higher-order contagion is one of the most active areas in network
science, and the key recent finding is that hyperedge overlap — not the mere presence
of group interactions — drives explosive transitions and bistability (*Nature
Communications*, 2024). A 2026 preprint reports that adaptive behavior neutralises
those bistable explosive transitions, but the adaptation there is a prescribed
risk-perception rule. Nobody has asked whether *derived* group-avoidance does the same
thing, or whether utility-maximising agents avoid exactly the hyperedges whose overlap
drives explosivity.

## Question

When the action is "which groups do I attend" rather than "which edges do I drop," does
the cross-scale asynchrony survive? And does derived avoidance neutralise explosive
transitions the way prescribed avoidance does — or does it target overlapping
hyperedges preferentially and neutralise them more efficiently?

## What to build

- Reformulate the decision over venues: each node has an activity set, each venue has a
  size and a local prevalence, and the choice is a subset of venues to attend. Infection
  risk becomes venue-level, `1 − (1−β)^(infectious co-attendees)`.
- Utility over venues rather than edges — venues differ in benefit, which is where
  "I can skip the gym but not work" enters naturally, and which connects directly to the
  capacity channel of the active proposal.
- Run on synthetic hypergraphs with controlled hyperedge overlap first, then on the
  Manassas activity structure.
- Measure whether the avoided hyperedges are overlap-enriched relative to random
  avoidance of equal volume.

## Why this framework

The venue formulation is closer to the twin's native data *and* closer to how people
actually behave, and it removes an assumption the current model needs but cannot
justify — that a susceptible node knows which specific neighbours are infectious.
Venue-level prevalence is a far more defensible information set, which makes this the
natural companion to the imperfect-detection proposal.

## Risk

The action space grows combinatorially over venue subsets. Restrict to a ranked greedy
choice (attend the top `k` venues by benefit-to-risk), which stays tractable and is
arguably more behaviourally realistic than exact subset optimisation anyway.

## Practicalities

- **New machinery:** hypergraph substrate; venue-level decision; greedy action selection
- **External data:** Manassas activity–location assignments (already in hand at UVA)
- **Compute:** high — sweeps crossed with ensembles need cluster time
- **Venue:** Nature Communications · PNAS · PRX Life

## References

- Hyperedge overlap drives explosive transitions in systems with higher-order
  interactions. *Nature Communications* (2024).
  https://www.nature.com/articles/s41467-024-55506-1
- Adaptive behaviors neutralize bistable explosive transitions in higher-order
  contagion. arXiv:2601.05801 (2026). https://arxiv.org/pdf/2601.05801
