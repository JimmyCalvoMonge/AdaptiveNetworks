# A taxonomy in four axes

Almost every behavior-coupled network epidemic model can be located by answering four
questions independently.

## Axis 1 — the action: what does behavior actually change?

1. **Susceptibility** — masks, hygiene, distancing modelled as a reduction in per-contact
   risk. Contact structure untouched.
2. **Rewiring** — drop a risky edge, form a replacement. Degree conserved.
3. **Deletion / deactivation** — drop without replacement. Degree falls. Permanent or
   temporary.
4. **Activity modulation** — how often you participate at all, or which venues you attend.
5. **Discrete protective act** — vaccinate, isolate, adopt an app. Usually irreversible.
6. **Mobility** — travel and migration between subpopulations.

## Axis 2 — the generator: what produces the action? *(the decisive axis)*

1. **None** — static network, no behavior.
2. **Assigned rate** — rewire with probability `w`. The rate is a parameter.
3. **Prevalence function** — a chosen functional form mapping observed prevalence to response.
4. **Behavioral contagion** — fear or awareness spreads as its own epidemic process.
5. **Threshold / social reinforcement** — act when enough neighbours do.
6. **Imitation / evolutionary game** — copy the neighbour with the better payoff.
7. **Strategic equilibrium** — Nash, differential game, mean-field game.
8. **Utility maximisation / dynamic programming** — solve an individual optimal-control problem.
9. **Learned policy** — reinforcement learning.
10. **Generative agent** — an LLM decides.

Classes 2–5 are phenomenological in the strict sense: the response is described, not
derived. Classes 6–8 derive it from something, differing in what the something is
(social comparison, equilibrium, private optimisation). Classes 9–10 replace the
derivation with a fitted or pretrained policy, which relocates the prescription rather
than removing it.

### The generator axis in detail

| Class | Behavior comes from | Representative work | What it buys | What it cannot do |
| --- | --- | --- | --- | --- |
| 2 · Assigned rate | A parameter: rewire w.p. `w`, delete w.p. `κ` | Gross 2006; Shaw & Schwartz 2008; Kiss 2012; Ball & Britton 2022 | Analytic tractability; bifurcation structure; rigorous thresholds | Say why the rate has its value, or predict it under new conditions |
| 3 · Prevalence function | A chosen form, e.g. `N′ = N(1−Θ^α)` | Perra 2011; Weitz 2020; Maharaj & Kleczkowski 2012; Arthur 2021 | Rich dynamics — plateaus, shoulders, oscillations — with few parameters | Separate the functional form from the phenomenon it was chosen to produce |
| 4 · Behavioral contagion | A second spreading process | Epstein 2008; Funk 2009; Granell 2013; Moinet 2018 | Awareness has its own topology and timescale; metacritical points | Explain why an individual adopts — transmission replaces decision |
| 5 · Threshold | Act once a fraction of neighbours act | Centola & Macy; Morsky 2023 | Social reinforcement, critical mass, clustering effects | Handle heterogeneous private costs; the threshold is assigned |
| 6 · Imitation | Copy the better-performing neighbour | Bauch & Earn 2004; Fu 2011; Salathé & Bonhoeffer 2008 | Free-riding, opinion clustering, herd-immunity dilemmas emerge | Represent forward-looking choice; payoffs must be observable |
| 7 · Strategic equilibrium | Nash / differential / mean-field game | Reluga 2010; Eksin 2017; Schnyder 2026 | Welfare statements; private versus social optimum | Scale to heterogeneous individuals on an explicit graph |
| **8 · Utility maximisation** | **An individual optimal-control problem, solved per agent** | **Fenichel 2011; Morin 2013; Espinoza 2021–2025; this framework** | **Preferences are the primitive, so agents re-optimise under new constraints, incentives or information** | **Analytic tractability; the utility form is still a modelling choice** |
| 9 · Learned policy | Reinforcement learning over a reward | Cognitively-plausible RL 2025; reward-engineering platforms 2026 | No functional form assumed; complex policies discoverable | Interpretability; the reward is the new prescription |
| 10 · Generative agent | An LLM reasons about the situation | Williams 2023; GABLE 2026; KDD 2026 | Rich, demographically conditioned behavior | Calibration, identifiability, reproducibility |

## Axis 3 — the information: what does the decision see?

1. **Global prevalence** — aggregate case counts, media, official reports.
2. **Local neighbourhood** — health states of direct contacts.
3. **n-step neighbourhood** — a wider radius (exact risk computation is NP-hard).
4. **Own state / own history** — memory, fatigue, private experience.
5. **Neighbours' behavior** — not their health, but what they are doing.
6. **A separate belief layer** — awareness held and transmitted independently.

Each admits **delayed** and **noisy / partially observed** variants, which are modelled
far less often than the idealised versions.

## Axis 4 — the substrate: what does the process run on?

1. **Mean-field** — well-mixed compartments.
2. **Static network** — ER, BA, WS, configuration model.
3. **Adaptive network** — topology and state coevolve.
4. **Temporal / activity-driven** — edges exist momentarily.
5. **Multiplex** — several layers, usually contact plus information.
6. **Higher-order** — hypergraphs, simplicial complexes, group interactions.
7. **Metapopulation** — subpopulations coupled by mobility.
8. **Empirical / synthetic** — measured contacts or a digital twin.

## Where the cells are empty

Crossing the action axis with the generator axis shows the literature's actual shape.

| Action ↓ Generator → | Assigned rate | Prevalence fn | Contagion | Imitation | Equilibrium | Utility / DP | RL / LLM |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Susceptibility | dense | dense | dense | some | dense | dense (mean-field) | some |
| Rewiring | dense | some | some | sparse | sparse | sparse | sparse |
| **Deletion / deactivation** | some | sparse | sparse | — | sparse | **this framework** | — |
| Activity / venue choice | dense | some | some | sparse | sparse | — | sparse |
| Discrete act (vaccinate) | some | some | some | dense | dense | sparse | sparse |
| Mobility | dense | some | sparse | — | sparse | — | sparse |

"Dense" means a substantial sub-literature with reviews of its own; "sparse" means
isolated papers; "—" means nothing found in this survey.

The framework's cell — temporary deletion driven by a per-node dynamic program on a
coevolving graph — has close neighbours on each side but no direct occupant:

- **Same action, assigned rate:** Tunc, Shkarayev & Shaw (2013) temporary link
  deactivation; Ball, Britton, Leung & Sirl (2019) preventive dropping of edges.
- **Same generator, mean-field substrate:** Fenichel et al. (2011); Espinoza et al.
  (2021–2025).
- **Network plus utility, different action:** Eksin, Shamma & Weitz (2017) protective
  and pre-emptive effort; Qiu et al. (2022) mask adoption.
