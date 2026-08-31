# When the infected and the recovered also decide

**Changes:** who decides at all
**Status:** shelved — completes the model to all three health classes

## Motivation

In the baseline, only susceptibles decide. The Methods state that infected and
recovered individuals "have no incentive to modulate their social activity," and the
SI's version of infected responsiveness is a blunt instrument — a random half of
infected nodes drop three quarters of their edges. The manuscript flags the proper
extension itself: "prosociality of infected individuals."

Two live literatures have made this urgent. PNAS has published a theory of epidemics
with altruism finding that even *extremely* weak altruism — an individual valuing their
own life at roughly a hundred thousand others' — is enough for rational self-isolation
when infected to change outcomes almost as much as strong altruism does. Separately,
npj Complexity reports the emergence of shield immunity in spatial contagions,
following Weitz and colleagues' interaction-substitution idea, in which recovered
individuals deliberately absorb risky interactions that would otherwise fall on
susceptibles.

Both are decisions under a utility. Neither has been done with a derived, node-level,
network-embedded optimisation — which is precisely what this framework already computes
for one health class and could compute for all three.

## Question

How much prosociality by the infected substitutes for how much avoidance by the
susceptible? And does recovered-node shielding change the asynchrony — given that the
recovered pool is exactly what grows through the peak, and is exactly the population
the current model leaves passively holding every edge it has?

## What to build

- An altruism weight `α` in the infected node's utility: it internalises a fraction `α`
  of the expected harm it imposes on susceptible neighbours. At `α=0` this reduces
  exactly to the current model; at `α=1` the node is a social planner over its own
  edges. Sweep it.
- A shielding motive for recovered nodes — utility from substituting into contacts that
  would otherwise be S–I.
- The substitution surface: for each `α`, the susceptible risk sensitivity `ν` that
  yields the same peak. That trade-off is the central figure, and it is a quantity no
  aggregate model can produce.
- Serology gating — shielding only works if recovered status is known, so sweep the
  fraction of recovered nodes aware of their status. This is the same axis as the
  imperfect-detection proposal's detection probability, applied to a different class.

## Why this framework

The externality is currently structural: agents are naive to the harm they cause, by
construction. Making `α` a dial converts that from an assumption into a measurement, so
the distance between decentralised behavior and the social optimum can be quantified
*inside* the model — with no planner, no budget and no departure from the decentralised
setting the paper argues for.

## Risk

Shielding requires recovered nodes to *seek* edges, and edge addition is the one thing
the framework has deliberately refused, since its whole case against rewiring models
rests on it. Frame it precisely as temporary, motivated substitution rather than
topology-altering rewiring, and defend the distinction explicitly. A reviewer who reads
it as abandoning the paper's own premise is not being unreasonable.

## Practicalities

- **New machinery:** infected and recovered MDPs; altruism weight; substitution mechanism
- **External data:** none
- **Compute:** moderate
- **Venue:** PNAS · Proc. R. Soc. B · Journal of Theoretical Biology

## References

- The theory of epidemics with altruism. *PNAS*.
  https://www.pnas.org/doi/abs/10.1073/pnas.2518893123
- Emergence of shield immunity during spatial contagions. *npj Complexity* (2025).
  https://www.nature.com/articles/s44260-025-00044-0
- Weitz et al. Intervention Serology and Interaction Substitution: Modeling the Role of
  "Shield Immunity" in Reducing COVID-19 Epidemic Spread.
  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7276032/
