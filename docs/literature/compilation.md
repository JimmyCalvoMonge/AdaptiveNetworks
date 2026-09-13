# Compilation

Fourteen families. Entries marked ★ are the closest prior art or define the state of a
question a new framework would have to answer. Citation details compiled from web
search — verify against publisher records before use in a manuscript.

---

## 01 · Reviews and agenda-setting

The two most recent reviews face opposite directions: one finds behavioral models
proliferating, the other finds the behavioral variables they use do not correlate with
the behaviors they claim to represent.

- Funk, Salathé & Jansen. Modelling the influence of human behaviour on the spread of
  infectious diseases: a review. *J. R. Soc. Interface* 7, 1247–1256 (2010).
  The founding review. Introduced the distinction between behavior affecting disease
  state, contact structure, and model parameters; and between local and global information.
- ★ Verelst, Willem & Beutels. Behavioural change models for infectious disease
  transmission: a systematic review (2010–2015). *J. R. Soc. Interface* 13, 20160820 (2016).
  178 papers under PRISMA. The head count matters: 76 encode behavior as information
  entering as a dynamic parameter, 37 use an economic objective function without
  imitation, 26 with imitation. The economic-objective family was roughly a third of the
  field a decade ago, and almost entirely mean-field.
- Wang, Andrews, Wu, Wang & Bauch. Coupled disease–behavior dynamics on complex networks:
  A review. *Physics of Life Reviews* 15, 1–29 (2015). Published with a full round of
  commentaries. Organises by behavior type rather than by generator.
- Pastor-Satorras, Castellano, Van Mieghem & Vespignani. Epidemic processes in complex
  networks. *Rev. Mod. Phys.* 87, 925–979 (2015). The substrate reference; behavior
  largely out of scope, which is itself informative.
- Bedson et al. A review and agenda for integrated disease models including social and
  behavioural factors. *Nature Human Behaviour* 5, 834–846 (2021).
- ★ Recent trends in socio-epidemic modelling: behaviours and their determinants.
  *Boll. Unione Mat. Ital.* / arXiv:2506.13837 (2025). Distinguishes modelling
  *behaviours* from modelling *determinants* (awareness, beliefs, trust), notes most
  papers conflate them, then tests on Italian regional COVID data and finds behavioural
  responses **poorly explained** by awareness, beliefs or trust. A direct challenge to
  the awareness-contagion family.
- Game-Theoretic Frameworks for Epidemic Spreading and Human Decision-Making: A Review.
  arXiv:2106.00214 (2021).
- Keeling & Eames. Networks and epidemic models. *J. R. Soc. Interface* 2, 295–307 (2005).
  Kiss, Miller & Simon. *Mathematics of Epidemics on Networks*. Springer (2017).

---

## 02 · Adaptive networks with assigned rewiring rates

The largest and most mathematically developed family. Behavior is a rate; the payoff is
rigorous analysis. This is the family the framework argues against — and it should be
argued against on what the rate cannot represent, not on the quality of the work.

- ★ Gross, D'Lima & Blasius. Epidemic dynamics on an adaptive network. *Phys. Rev. Lett.*
  96, 208701 (2006). The origin. Susceptibles rewire away from infected at rate `w`,
  reconnecting to a random susceptible. Produces assortative degree correlation,
  oscillations, hysteresis, first-order transitions. Low-dimensional moment-closure model
  plus full local bifurcation analysis.
- Gross & Kevrekidis. Robust oscillations in SIS epidemics on adaptive networks:
  coarse-graining by automated moment closure. *EPL* (2008).
- Shaw & Schwartz. Fluctuating epidemics on adaptive networks. *Phys. Rev. E* 77, 066101
  (2008); Enhanced vaccine control of epidemics in adaptive networks. *Phys. Rev. E* 81,
  046120 (2010).
- Marceau, Noël, Hébert-Dufresne, Allard & Dubé. Adaptive networks: coevolution of disease
  and topology. *Phys. Rev. E* 82, 036116 (2010). Improved compartmental formalism
  tracking node state *and* neighbourhood composition.
- Guo, Trajanovski, van de Bovenkamp, Wang & Van Mieghem. Epidemic threshold and
  topological structure of SIS epidemics in adaptive networks. *Phys. Rev. E* 88, 042802
  (2013). Source of the evolving-versus-adaptive distinction.
- Risau-Gusman & Zanette. Contact switching as a control strategy for epidemic outbreaks.
  *J. Theor. Biol.* 257, 52–60 (2009).
- Jolad, Liu, Schmittmann & Zia. Epidemic spreading on preferred degree adaptive networks.
  *PLOS ONE* 7, e48686 (2012). Nodes adjust toward a target degree — an early attempt to
  anchor the rule individually, but the target is still assigned.
- Kiss, Berthouze, Taylor & Simon. Modelling approaches for simple dynamic networks.
  *Proc. R. Soc. A* 468, 1332–1355 (2012). Random link activation-deletion (RLAD).
- Szabó-Solticzky, Berthouze, Kiss & Simon. Oscillating epidemics in a dynamic network
  model. *J. Math. Biol.* (2016).
- Dong, Yin, Liu, Yan & Shi. Can rewiring strategy control the epidemic spreading?
  *Physica A* 438, 169–177 (2015). Zhu, Zhi, Guo & Wang. Analysis of epidemic spreading
  process in adaptive networks. *IEEE TCAS-II* 66, 1252–1256 (2019).
- Vazquez, Serrano & San Miguel. Rescue of endemic states in interconnected networks with
  adaptive coupling. *Sci. Rep.* 6, 29342 (2016).
- Gross & Blasius. Adaptive coevolutionary networks: a review. *J. R. Soc. Interface* 5,
  259–271 (2008). Sayama et al. Modeling complex systems with adaptive networks (2013).

---

## 03 · Deletion and deactivation without rewiring

The closest prior art on the **action**. These papers already establish that
dropping-without-replacing differs from rewiring, and two already have temporary
deactivation. What none has is a decision behind the drop.

- ★ **Tunc, Shkarayev & Shaw. Epidemics in adaptive social networks with temporary link
  deactivation. *J. Stat. Phys.* 151, 355–366 (2013).** The nearest neighbour on
  mechanism. Healthy individuals temporarily deactivate links to sick ones; links
  reactivate once both are healthy. Two regimes: slow network dynamics, where
  deactivation merely reduces effective contacts, and fast dynamics, where it efficiently
  targets dangerous connections. The deactivation and reactivation rates are parameters —
  that is the entire gap between this and a derived decision.
- ★ **Ball, Britton, Leung & Sirl. A stochastic SIR network epidemic model with preventive
  dropping of edges. *J. Math. Biol.* 78, 1875–1951 (2019).** The rigorous treatment of
  dropping without rewiring on configuration-model networks: effective-degree
  formulation, law of large numbers and functional central limit theorems, both
  Molloy–Reed and Newman–Strogatz–Watts variants. The mathematical standard any new
  dropping model will be measured against.
- Ball & Britton. Epidemics on networks with preventive rewiring. *Random Struct.
  Algorithms* 61, 250–297 (2022). Britton, Juher & Saldaña. A network epidemic model with
  preventive rewiring: comparative analysis of the initial phase. *Bull. Math. Biol.* 78
  (2016). Rewire-to-susceptible / -recovered / -non-infected compared, with `R₀` by both
  branching-process and pair approximation.
- ★ Sahneh, Vajdi, Melander & Scoglio. Contact adaption during epidemics: a multilayer
  network formulation approach. *IEEE Trans. Netw. Sci. Eng.* 6, 16–30 (2019). Each agent
  has a default contact set and an alternative, switching when alerted to infection among
  defaults. Key warning: adaptation that is **not fast enough** can *lower* network
  robustness. Non-monotonic in the adaptation rate.
- Preserving system activity while controlling epidemic spreading in adaptive temporal
  networks. *Phys. Rev. Research* 6, 033159 (2024). Frames adaptation as a trade-off
  between suppression and retained social activity — the same tension this framework puts
  in a utility function, handled here as constrained optimisation over the network
  process rather than per node.
- Balancing quarantine and self-distancing measures in adaptive epidemic networks.
  *Bull. Math. Biol.* (2022).

---

## 04 · Prevalence-driven response functions

Behavior as a chosen function of observed prevalence. Analytically productive; the
function is selected, and the phenomena it produces are downstream of that selection.

- ★ Perra, Balcan, Gonçalves & Vespignani. Towards a characterization of behavior-disease
  models. *PLOS ONE* 6, e23084 (2011). Prototypical mechanisms for self-initiated
  distancing driven by local versus non-local prevalence information, as transitions into
  behavioral classes. Rich phase space with multiple peaks and tipping points. Behavior
  alters susceptibility **without** altering contact patterns — the opposite design choice
  from this framework.
- ★ Weitz, Park, Eksin & Dushoff. Awareness-driven behavior changes can shift the shape of
  epidemics away from peaks and toward plateaus, shoulders, and oscillations. *PNAS* 117,
  32764–32771 (2020). The reference result for "behavior changes the shape, not just the
  size."
- ★ Manrubia & Zanette. Individual risk-aversion responses tune epidemics to critical
  transmissibility (R = 1). *R. Soc. Open Sci.* 9, 211667 (2022). Time-varying individual
  risk aversion generates multiple waves of decreasing amplitude tuning `R` toward 1.
  Successive waves infect individuals of gradually lower risk propensity, shaping a
  well-defined risk-aversion profile across the population. The final state is
  self-organised and parameter-independent. **The strongest existing claim that
  heterogeneous risk attitudes plus selection produce the observed multi-wave structure.**
- Maharaj & Kleczkowski. Controlling epidemic spread by social distancing: do it well or
  not at all. *BMC Public Health* 12, 679 (2012). Spatial rewiring with an individual risk
  parameter and explicit economic cost: `N_new = N₀(1 − Θ^α)`. The closest phenomenological
  analogue with individual-level sensitivity.
- Arthur, Jones, Bonds, Ram & Feldman. Adaptive social contact rates induce complex
  dynamics during epidemics. *PLOS Comput. Biol.* 17, e1008639 (2021).
- Poletti, Ajelli & Merler. Risk perception and effectiveness of uncoordinated behavioral
  responses in an emerging epidemic. *Math. Biosci.* (2012); Lattice model for influenza
  spreading with spontaneous behavioral changes. *PLOS ONE* (2013).

---

## 05 · Awareness and fear as a second contagion

Behavior is a state that spreads. Strength: information has its own topology and
timescale. Weakness: adoption is transmission rather than decision — and the 2025
socio-epidemic review's empirical finding is aimed squarely here.

- Epstein, Parker, Cummings & Hammond. Coupled contagion dynamics of fear and disease.
  *PLOS ONE* 3, e3955 (2008). Fear is caught from the infected and the fearful alike; the
  fearful remove themselves from circulation.
- Funk, Gilad, Watkins & Jansen. The spread of awareness and its impact on epidemic
  outbreaks. *PNAS* 106, 6872–6877 (2009).
- ★ Granell, Gómez & Arenas. Dynamical interplay between awareness and epidemic spreading
  in multiplex networks. *Phys. Rev. Lett.* 111, 128701 (2013). UAU-SIS on a two-layer
  multiplex solved by microscopic Markov chain approximation. Identifies a **metacritical
  point** at which epidemic onset begins to depend on the awareness process.
- Granell, Gómez & Arenas. Competing spreading processes on multiplex networks. *Phys.
  Rev. E* 90, 012808 (2014). Mass media makes the metacritical point **disappear** — a
  clean demonstration that adding a global information channel changes the qualitative
  result.
- ★ Epidemic paradox induced by awareness driven network dynamics. *Phys. Rev. Research*
  7, L012061 (2025). Awareness by susceptible-only, infected-only, or all nodes on
  scale-free networks. Epidemic size scales linearly in the susceptible-aware and
  all-aware cases but **sublinearly** in the infected-aware case — fewer aware nodes can
  reduce the epidemic more, a Braess-like paradox. Evidence that **who** adapts matters
  more than how many.
- Moinet, Pastor-Satorras & Barrat. Effect of risk perception on epidemic spreading in
  temporal networks. *Phys. Rev. E* 97, 012313 (2018).
- Epidemic risk perception and social interactions lead to awareness cascades on multiplex
  networks. *J. Phys. Complexity* (2025). Contagion dynamics on adaptive multiplex
  networks with awareness-dependent rewiring. arXiv:2012.14073. Coupled epidemic dynamics
  with awareness heterogeneity in multiplex networks. *Chaos Solitons Fractals* (2024).

---

## 06 · Games, imitation and evolutionary dynamics

Owns the vaccination problem almost completely, and owns the welfare vocabulary —
private versus social optimum — that a utility framework can borrow without a planner.

- Bauch & Earn. Vaccination and the theory of games. *PNAS* 101, 13391–13394 (2004);
  Bauch, Galvani & Earn. Group interest versus self-interest in smallpox vaccination
  policy. *PNAS* 100, 10564–10567 (2003). Self-interest can preclude eradication; coverage
  is much harder to restore after a scare than it was to lose.
- ★ Reluga. Game theory of social distancing in response to an epidemic. *PLOS Comput.
  Biol.* 6, e1000793 (2010). Differential game solved by Pontryagin's maximum principle.
  Distancing is most valuable around `R₀ ≈ 2`; optimal distancing never recovers more than
  30% of the cost of infection.
- Fu, Rosenbloom, Wang & Nowak. Imitation dynamics of vaccination behaviour on social
  networks. *Proc. R. Soc. B* 278, 42–49 (2011); Zhang et al. *PLOS Comput. Biol.* 8,
  e1002469 (2012). Imitation can **exacerbate** transmission when vaccination is cheap,
  because non-vaccinators cluster socially.
- ★ Salathé & Bonhoeffer. The effect of opinion clustering on disease outbreaks.
  *J. R. Soc. Interface* 5, 1505–1508 (2008). Opinion clustering dramatically changes
  outbreak size **at constant overall coverage**.
- ★ Eksin, Shamma & Weitz. Disease dynamics in a stochastic network game: a little empathy
  goes a long way in averting outbreaks. *Sci. Rep.* 7, 44122 (2017). Both healthy and
  sick take costly action — protective and **pre-emptive** respectively. A critical level
  of empathy by the sick above which disease is eradicated rapidly; risk-averse behavior
  by the healthy alone cannot eradicate without it.
- Eksin, Paarporn & Weitz. Systematic biases in disease forecasting — the role of behavior
  change. *Epidemics* 27, 96–105 (2019).
- Morsky, Magpantay, Day & Akçay. The impact of threshold decision mechanisms of
  collective behavior on disease spread. *PNAS* 120, e2221479120 (2023).
- ★ Young, Silk, Pritchard & Fefferman. Diversity in valuing social contact and risk
  tolerance leading to the emergence of homophily in populations facing infectious
  threats. *Phys. Rev. E* 105, 044315 (2022). **The closest existing work on heterogeneous
  risk preferences shaping network structure.** Diversity along just two dimensions —
  value of contact, risk tolerance — is sufficient for self-organising homophily. Essential
  prior art for any endogenous-profile framework.
- ★ The interplay of social constraints and individual variation in risk tolerance in the
  emergence of superspreaders. *J. R. Soc. Interface* 20, 20230077 (2023). Explicitly
  separates **social constraint** from **risk tolerance**. The rare paper that already
  treats "cannot reduce contacts" as distinct from "chooses not to."
- The theory of epidemics with altruism. *PNAS* (2026). Even extremely weak altruism
  suffices for rational self-isolation to change outcomes.
- Schnyder, Molina, Miller, Yamamoto, Kobayashi & Turner. Self-organized social distancing
  when the force of infection depends on susceptible and infectious behavior. *Math.
  Biosci. Eng.* 23 (2026). Mean-field game with both parties' behavior entering the force
  of infection.
- Game-theoretic behavioral adaptation in non-Markovian epidemic spreading on networks.
  *Chaos Solitons Fractals* (2026). Hota & Sundaram. *IFAC* (2019). A game theoretical
  analysis of voluntary mask wearing over complex networks. *Dyn. Games Appl.* (2025).

---

## 07 · Utility maximisation and economic epidemiology

This framework's own family. Note how much of it is mean-field — that is the gap the
network formulation fills.

- ★ Fenichel et al. Adaptive human behavior in epidemiological models. *PNAS* 108,
  6306–6311 (2011). The origin: individuals choose contact rates to maximise expected
  utility over a planning horizon, solved as a dynamic program, coupled back into SIR.
- Fenichel. *J. Health Econ.* 32, 440–451 (2013). Morin, Fenichel & Castillo-Chavez. SIR
  dynamics with economically driven contact rates. *Nat. Resour. Model.* 26, 505–525
  (2013). Perrings et al. *EcoHealth* 11, 464–475 (2014).
- ★ Espinoza, Marathe, Swarup & Thakur. Adaptive human behavior in epidemics: the impact
  of risk misperception. *Sci. Rep.* 11, 19744 (2021); Asymptomatic individuals can
  increase the final epidemic size under adaptive human behavior. *Sci. Rep.* 11 (2021).
  The asymptomatic result — undetectable infection *increases* final size specifically
  because behavior is adaptive — is the mean-field precursor to any detectability work on
  networks.
- Espinoza, Swarup, Barrett & Marathe. Heterogeneous adaptive behavioral responses may
  increase epidemic burden. *Sci. Rep.* 12, 11276 (2022).
- Espinoza, Saad-Roy, Grenfell, Levin & Marathe. *Proc. R. Soc. B* 291, 20241772 (2024);
  The impact of risk compensation adaptive behavior on the final epidemic size. *Math.
  Biosci.* 380, 109370 (2025).
- Traulsen, Levin & Saad-Roy. Individual costs and societal benefits of interventions
  during the COVID-19 pandemic. *PNAS* 120, e2303546120 (2023).
- ★ Qiu, Espinoza, Vasconcelos, Chen, Constantino, Crabtree, Yang, Vullikanti, Chen,
  Weibull, Basu, Dixit, Levin & Marathe. Understanding the coevolution of mask wearing and
  epidemics: a network perspective. *PNAS* 119, e2123355119 (2022). Behavior-epidemic
  coevolution on a network with an economic decision layer, from an overlapping author
  group. Robust non-monotonic relation between attack rate and transmission probability,
  with an abrupt drop at a critical threshold; regimes producing multiple waves of both
  infection and mask adoption. **The most direct network-scale precedent within this
  research programme.**
- Endogenous social distancing and its underappreciated impact on the epidemic curve.
  *Sci. Rep.* 11, 3093 (2021). Endogenous social distancing and containment policies in
  social networks. *Natl. Inst. Econ. Rev.* Epidemic spreading and equilibrium social
  distancing in heterogeneous networks. arXiv:2007.04210.

---

## 08 · Temporal and activity-driven substrates

Behavioral work here is thinner than the substrate literature — an opportunity.

- Perra, Gonçalves, Pastor-Satorras & Vespignani. Activity driven modeling of time varying
  networks. *Sci. Rep.* 2, 469 (2012).
- Epidemic spreading with awareness diffusion on activity-driven networks. *Phys. Rev. E*
  98, 062322 (2018). Impact of individual behavioral changes on epidemic spreading in
  time-varying networks. arXiv:2107.14143. Epidemic criticality in temporal networks.
  *Phys. Rev. Research* 6, L022017 (2024). Local perception and the **duration** of the
  preventive effect are the decisive parameters.
- Burstiness in activity-driven networks and the epidemic threshold. arXiv:1903.11308.
  Spatiotemporal activity-driven networks. *Phys. Rev. E* (2025).

---

## 09 · Higher-order structure and group interactions

The fastest-moving substrate family. Adaptation arrived only in the last two years and is
still entirely rule-based — the clearest open cell in the grid.

- Iacopini, Petri, Barrat & Latora. Simplicial models of social contagion. *Nat. Commun.*
  10, 2485 (2019).
- ★ Hyperedge overlap drives explosive transitions in systems with higher-order
  interactions. *Nat. Commun.* 16 (2024). Higher-order interaction alone does not
  guarantee abrupt transitions — explosivity and bistability require low intra-order
  hyperedge overlap.
- Characteristic scales and adaptation in higher-order contagions. *Nat. Commun.* 16
  (2025). Adaptive behaviors neutralize bistable explosive transitions in higher-order
  contagion. arXiv:2601.05801 (2026) — risk perception conveyed through the same
  interactions along which contagion occurs, but the rule is prescribed.
- Adaptive epidemic dynamics on hypergraphs with group-level immunization and rewiring.
  arXiv:2606.15578 (2026). An adaptive simplicial SIS where node states and hyperedge
  **activity** coevolve — structurally very close to a venue-attendance formulation, with
  an assigned rule where a decision could sit. Epidemic dynamics driven by adaptive
  rewiring on higher-order networks. *Chaos Solitons Fractals* (2025). Coupled dynamics of
  vaccination behavior and epidemic spreading on multilayer higher-order networks (2026).
- Iacopini, Karsai & Barrat. The temporal dynamics of group interactions in higher-order
  social networks. *Nat. Commun.* 15, 7391 (2024). Higher-order interactions shape
  collective human behaviour. *Nat. Hum. Behav.* (2025).

---

## 10 · Empirical networks, agent-based models and digital twins

Behavior here is usually scripted at the activity level rather than decided.

- ★ Eubank, Guclu, Kumar, Marathe, Srinivasan, Toroczkai & Wang. Modelling disease
  outbreaks in realistic urban social networks. *Nature* 429, 180–184 (2004). EpiSims:
  dynamic bipartite people-location graphs from census and land-use data. The origin of
  the digital-twin lineage the Manassas network comes from — and note the substrate is
  **venue-first** by construction.
- Barrett et al. EpiSimdemics. *SC* (2008); Generation and analysis of large synthetic
  social contact networks. *WSC* (2009). Mortveit et al. Synthetic populations and
  interaction networks for the U.S. NSSAC (2020).
- ★ Chang, Pierson, Koh, Gerardin, Redbird, Grusky & Leskovec. Mobility network models of
  COVID-19 explain inequities and inform reopening. *Nature* 589, 82–87 (2021).
  Metapopulation SEIR on hourly mobility networks for 98M people across POIs in ten US
  metros. Explains higher infection rates among disadvantaged groups mechanistically:
  lower-income neighbourhoods could not reduce mobility as much, and the POIs they visited
  were more crowded. **The empirical anchor for capacity-constrained behavior.**
- Sánchez, Calvo, García, Vásquez & Barboza. An implementation of a multilayer network
  model for the COVID-19 pandemic: a Costa Rica study. *Math. Biosci. Eng.* 20, 534–551
  (2023); A multilayer network model of COVID-19: implications in public health policy in
  Costa Rica. *Epidemics* (2022). Household, social and sporadic layers with layer-specific
  engagement from fixed distributions.
- SocioPatterns collaboration (RFID proximity: workplace, school, conference). Copenhagen
  Networks Study (~1000 students; Bluetooth proximity, calls, SMS, online ties, 2012–2016).
  Aleta, Ferraz de Arruda & Moreno. Data-driven contact structures. *PLOS Comput. Biol.*
  16 (2020).
- CitySEIRCast: an agent-based city digital twin (2024). Chen et al. Prioritizing
  allocation of COVID-19 vaccines based on social contacts (2021).

---

## 11 · Fatigue, adherence decay and the long run

The newest active area. All of these treat fatigue at the **aggregate** level — a
parameter on a population adherence curve — rather than as a node-level state.

- ★ Fatigue and adherence can challenge the prevailing wisdom on the response to severe
  epidemic outbreaks. *J. R. Soc. Interface* 23, 20250287 (2026). Moderate control
  priority generates less intense actions, mitigating fatigue accumulation and keeping
  adherence high — a weaker intervention can outperform a stronger one.
- ★ Mohammed & Alsammani. Long-term coexistence of epidemics and risk awareness: impacts
  of adaptive human response and fatigue. arXiv:2607.18301 (2026). Awareness as a dynamic
  behavioral state governing accessibility for interaction, emerging endogenously from
  prevalence and decaying through fatigue. Produces transient waves converging to
  long-term coexistence; responsiveness and fatigue jointly regulate transient structure
  and endemic burden **without** altering the invasion threshold.
- Mahmud, Eshun, Espinoza & Kadelka. Adaptive human behavior and delays in information
  availability autonomously modulate epidemic waves. *PNAS Nexus* 4, pgaf145 (2025). If
  response is either too prompt or too delayed, multiple waves do not emerge; minimal
  final size occurs in the damped-oscillation regime.
- Controlling multiple COVID-19 epidemic waves: a multi-scale model linking behaviour
  change to transmission (2022). The tension between awareness and fatigue shapes COVID-19
  spread (Georgia Tech Quantitative Biosciences).

---

## 12 · Inequality, constraint and homophily

Establishes that adaptive capacity is unequally distributed — and that models treating
this as a difference in preference are misrepresenting it.

- ★ Addressing the socioeconomic divide in computational modeling for infectious diseases.
  *Nat. Commun.* 13, 2897 (2022). The field-level critique. Cites the Santiago de Chile
  analysis in which deep inequality and disparities in **achievable** mobility reduction
  significantly delayed the end of the first wave — a constraint result, not a preference
  result.
- Importance of social inequalities to contact patterns, vaccine uptake, and epidemic
  dynamics. *Nat. Commun.* 15 (2024).
- Socioeconomic inequality in compliance with precautions and health behavior changes
  during the COVID-19 outbreak (Korean Community Health Survey 2020). Relative inequality
  indices 1.20–3.05 for non-compliance by education and income — concrete magnitudes for
  calibrating a capacity floor.
- Impact of homophily in adherence to anti-epidemic measures. arXiv:2507.13848 (2025).
  Homophily impacts the success of vaccine roll-outs. *Commun. Phys.* 5 (2022). Herd
  immunity and epidemic size in networks with vaccination homophily. *Phys. Rev. E* 105,
  L052301 (2022).
- The Making and Breaking of Social Ties During the Pandemic. *Front. Sociol.* (2022).
  Social ties in old age: the effect of the COVID-19 pandemic (2025). Roughly a third lost
  contact with acquaintances, one in four with a friend; restructuring toward kin ties.

---

## 13 · Calibration, identifiability and evaluation

Read before claiming predictive value for any behavioral mechanism.

- ★ Comparative evaluation of behavioral epidemic models using COVID-19 data. *PNAS* 122,
  e2421993122 (2025). Head-to-head benchmark against real surveillance data across
  heterogeneous locations.
- ★ Parameter estimation in behavioral epidemic models with endogenous societal
  risk-response. *PLOS Comput. Biol.* 20, e1011992 (2024). Systematic biases arise even
  with complete, accurate disease data and a correctly specified model when data are
  limited to the first wave — because of the delay between evolving risk and societal
  reaction. Small amounts of public behavior data substantially improve accuracy.
- Mathematical analysis of simple behavioral epidemic models. *Math. Biosci.* (2024).
  Heterogeneous behavioral mechanisms in epidemiological models. arXiv:2606.14902 (2026) —
  Bayesian mixture partitioning a population into distinct behavioral patterns.
- Efficient and accurate simulation of infectious diseases on adaptive networks. *PLOS
  Complex Systems* (2025). High-acceptance sampling for exact simulation of
  adaptive-contact outbreaks.

---

## 14 · Learned and generative agents

Both families replace a prescribed rule with a fitted or pretrained policy, relocating
the prescription rather than removing it — but the LLM work is producing operationally
competitive results.

- Cognitively-plausible reinforcement learning in epidemiological agent-based simulations.
  *Front. Epidemiol.* (2025). Reward engineering for spatial epidemic simulations.
  arXiv:2511.18000 (2025). Reinforcement learning for policymaking in epidemic control: a
  scoping review (2026). Reward design dominates the result.
- ★ GABLE: Integrating adaptive human behavior into epidemic models with large language
  models. arXiv:2608.29535 (2026). An LLM infers behavioral responses to epidemic and
  policy conditions and emits age-structured contact matrices coupled to a mechanistic
  model. Applied to COVID-19 in France, LLM-generated contact matrices **outperformed
  mobility-derived matrices** in short-term forecasting, with the largest gains at longer
  horizons. The benchmark a mechanistic behavioral framework now has to beat or explain.
- Williams et al. Epidemic modeling with generative agents. arXiv:2307.04986 (2023). An
  infectious disease spread simulation based on LLM decision making. *KDD* (2026). LLM
  powered social digital twins. arXiv:2601.06111 (2026). An LLM-driven multi-agent
  simulation framework for coupled epidemic–economic dynamics (2026).
