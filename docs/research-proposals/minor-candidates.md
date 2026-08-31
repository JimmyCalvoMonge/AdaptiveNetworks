# Minor candidates

Four directions developed to a shallower depth and ranked below the shelved proposals,
with the reasons recorded so the judgement can be revisited.

## Risk perception as its own contagion on a second layer

Perception spreads through a different network than the pathogen does. The multiplex
awareness literature has a concrete theoretical target in the *metacritical point* — the
threshold above which awareness diffusion controls epidemic onset — and recent work
extends this to awareness cascades driven by risk perception and social interaction.
The repository's existing layer machinery in `SIR_adaptive_net_layers.py` would carry
it.

**Why lower:** the field is crowded, and the marginal contribution reduces to "we
replaced the awareness rule with an optimisation." The imperfect-detection proposal
attacks the information set from a more fundamental direction, and the group's prior
work already covers risk misperception in the mean-field setting.

- Epidemic risk perception and social interactions lead to awareness cascades on
  multiplex networks. *J. Phys. Complexity* (2025).
  https://iopscience.iop.org/article/10.1088/2632-072X/adb897

## An interpretable baseline for LLM-agent epidemic simulations

Generative-agent epidemic models are proliferating — a KDD 2026 paper on LLM
decision-making in disease spread, LLM-powered social digital twins, coupled
epidemic–economic multi-agent frameworks. They are compelling and almost entirely
uncalibratable. There is a real paper in fitting `(ν, τ)` to the decisions an LLM agent
actually makes, and asking whether that behavior is even internally consistent enough
to be described by a coherent utility. If it is, this framework is the tractable
surrogate; if it is not, that is a finding the generative-ABM field needs.

**Why lower:** it is a methods paper about someone else's method, and the group's
comparative advantage is the mechanism itself.

- An Infectious Disease Spread Simulation Based on Large Language Model Decision
  Making. *KDD* (2026). https://arxiv.org/pdf/2606.06360
- LLM Powered Social Digital Twins. arXiv:2601.06111. https://arxiv.org/html/2601.06111

## Faster exact simulation of adaptive-network outbreaks

PLOS Complex Systems published a high-acceptance-sampling method for exact simulation of
adaptive-contact outbreaks. Adopting it would materially help every sweep-heavy
direction here, the venue-level and fatigue work in particular.

**Why lower:** infrastructure rather than science. It belongs inside whichever
sweep-heavy project starts first, not as a paper of its own.

- Efficient and accurate simulation of infectious diseases on adaptive networks. *PLOS
  Complex Systems*.
  https://journals.plos.org/complexsystems/article?id=10.1371%2Fjournal.pcsy.0000049

## Multi-pathogen and strain competition

One risk perception, several pathogens: behavior adapted to the salient threat
modulates transmission of the unnoticed one. The machinery is a modest extension of the
imperfect-detection proposal's detectability split.

**Why lower:** it needs a concrete empirical anchor to avoid being a purely theoretical
exercise, and none of the current data holdings supply one.
