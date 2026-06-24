# µODE — Limitations: what we can and cannot model

This document is a deliberately complete and honest inventory of what the muODE
modelling approach can represent and what it cannot. It exists so that results
are never over-interpreted.

muODE's core is **community dynamic FBA** (the Static Optimization Approach):
each organism is a flux-balance LP solved at pseudo-steady state every time step,
the organisms share one pool of extracellular metabolites, and that pool plus the
biomasses evolve by ODEs. On top of that core sits an optional **ecology layer**
([ecology.md](ecology.md)) of non-metabolic processes.

The right way to read the categories below:

- **Modelled** — represented by the FBA core or an ecology layer.
- **Approximation** — represented, but with a coarse-grained model the user
  should understand and can refine.
- **Out of scope (paradigm)** — *cannot* be represented faithfully in a
  metabolism-ODE framework. We deliberately do **not** fake these; the honest
  answer is to use (or couple to) a tool built for that paradigm. Pretending a
  dFBA run captures them would produce confident, wrong numbers.

---

## 1. Metabolism and ecology — Modelled

| Process | How | Where |
|---------|-----|-------|
| Nutrient competition for shared substrates | shared extracellular pool, per-step uptake bounds | [engine.md](engine.md) |
| Cross-feeding / syntrophy | emergent from the shared pool | [engine.md](engine.md) |
| Michaelis–Menten / Monod uptake | per-step bound from concentration | [kinetics.md](kinetics.md) |
| Enzyme-capacity limits (proteome) | GECKO-lite kcat caps | [kinetics.md](kinetics.md) |
| Diet / medium, open-system influx | `Diet` concentrations + influx | [assembly.md](assembly.md) |
| Chemostat washout / dilution | `dilution_rate` | [engine.md](engine.md) |
| First-order cell death | `death_rate` | [engine.md](engine.md) |
| Perturbations (antibiotic-as-bound, knockout, removal) | `Perturbation` | [perturbation.md](perturbation.md) |
| Timed biomass events (FMT, probiotic dose) | `Injection` | [injection.md](injection.md) |
| Spatial structure (2-D reaction–diffusion) | `SpatialDynamicFBA` | [spatial.md](spatial.md) |
| **pH dynamics + weak-acid (SCFA) inhibition** | `WeakAcidInhibition` | [ecology.md](ecology.md) |
| **Bile-acid transformation (BSH, 7α-dehydroxylation)** | `BileAcidTransform` | [ecology.md](ecology.md) |
| **Secondary-bile-acid growth inhibition** | `BileAcidInhibition` | [ecology.md](ecology.md) |
| **Sporulation / germination life cycle** | `SporeForming` | [ecology.md](ecology.md) |
| **Antibiotic pharmacokinetics/dynamics (time-varying)** | `Antibiotic` | [ecology.md](ecology.md) |
| **Direct antagonism (bacteriocins)** | `Bacteriocin` | [ecology.md](ecology.md) |
| **Multi-kingdom metabolism (fungi / protists)** | any GEM via the same `OrganismModel` core | [kingdoms.md](kingdoms.md) |
| **Oxygen tolerance gating + scavenging** | `OxygenSensitivity` | [kingdoms.md](kingdoms.md) |
| **Bacteriophage predation (lytic / temperate)** | `PhageInfection` (coupled infection ODE) | [kingdoms.md](kingdoms.md) |

The **bold** rows were previously listed as limitations and are now implemented.
Eukaryote *metabolism* needs no engine change (the core is paradigm-defined, not
taxon-defined); the only kingdom-specific gap is reconstruction tooling (§2/§3).
Phage *population* dynamics are now modelled; phage/host *evolution* is not (§4).

---

## 2. Modelled, but as a documented approximation

These are represented; the model is coarse and the parameters are
literature-informed defaults, not fits. Treat outputs as qualitative unless you
calibrate.

- **Euler integration.** First-order; accuracy is `O(dt)`. Use `dt ≤ 0.1 h`.
  A higher-order integrator (RK4) would let larger steps be taken safely. The
  per-step LP makes this non-trivial (multiple solves per step); not yet done.
- **pH model.** A buffered linear titration (`pH = pH0 − acid/buffer_capacity`),
  not a full charge-balance with bicarbonate/phosphate buffering. Good for the
  acidification *trend*; calibrate `ph0`/`buffer_capacity` per compartment for
  quantitative pH.
- **Weak-acid toxicity / cardinal-pH window.** Generic Hill and Rosso forms with
  per-species `Ki`/`pH_opt`. Species-specific values must be supplied to be
  quantitative.
- **Bile-acid network.** Lumped (taurocholate→cholate→deoxycholate); the full
  human bile-acid set (chenodeoxycholate, lithocholate, ursodeoxycholate,
  glyco-/tauro-conjugates) and host enterohepatic recirculation are collapsed.
- **Sporulation trigger.** Driven by a growth-rate (nutrient-stress) threshold;
  real sporulation integrates quorum, Spo0A phosphorelay and other signals.
- **Antibiotic PK.** One-compartment, first-order elimination; no absorption
  phase, tissue compartments, protein binding or gut-lumen specifics.
- **Enzyme constraints.** GECKO-*lite*: a single global enzyme budget, not a
  measured proteome allocation with per-enzyme masses.
- **Stoichiometric maintenance.** Uses each GEM's ATP-maintenance reactions as
  reconstructed; non-growth maintenance is not separately calibrated.
- **Thermodynamic directionality.** Reaction reversibility is taken from the GEM;
  there is no ΔG/concentration-dependent feasibility check, so a secretion the
  LP finds optimal may be thermodynamically marginal. (See §3.)
- **Phage infection model.** `PhageInfection` uses an exponential latent period
  (first-order lysis at `1/latent_period`), not a fixed delay, and a single host
  per layer with a **fixed host range**. Burst size is an effective virion yield
  per gDW lysed. Captures the predation/amplification dynamics; calibrate per
  phage–host pair for quantitative titers. Temperate behaviour is a coarse
  lysogeny fraction, not an explicit prophage/induction state.
- **Oxygen layer.** `OxygenSensitivity` gates growth by tolerance class and lets
  facultatives scavenge O₂ (Monod), capturing the anaerobiosis-maintenance
  mechanism; it is not a full redox/electron-acceptor balance (see §3) and, for
  real GEMs that already exchange O₂, should run with `consume=False`.
- **Eukaryote cell-scale parameters.** Fungal/protist cells are ~10–100× the dry
  mass of bacteria; muODE works in gDW/L so this is consistent, but per-cell
  counts and growth-rate priors must be set for the organism, not inherited from
  bacterial defaults.

---

## 3. Implementable extensions (not yet built)

Tractable within the paradigm; would extend the ecology/engine layer. Listed so
the roadmap is explicit.

- **Thermodynamic FBA (ΔG constraints).** Couple to group-contribution data
  (eQuilibrator/component-contribution) to forbid thermodynamically infeasible
  flux directions and make reversibility concentration-dependent. Hooks: the
  `OrganismModel` bound-setting and the metabolite pool.
- **Gas / redox coupling.** O₂ tolerance + scavenging is now modelled
  (`OxygenSensitivity`, §1). The remaining extension is the full gas/redox set —
  explicit H₂, CO₂, CH₄ balances, an electron-acceptor hierarchy and a spatial
  oxygen gradient (`uptake_factor` + the spatial engine are the hooks) — important
  for hydrogenotrophy and methanogenesis.
- **Eukaryote / virus reconstruction tooling.** *Simulating* fungi/protists needs
  no engine change, but *building* their GEMs is an external-tool gap: CarveMe is
  prokaryote-only, so a fungal route (Yeast8 template, CarveFungi, gapseq fungal
  mode) must be wired into the reconstruction phase. `reconstruction_route()`
  already flags the correct path per `Domain`; the automated workflow rule is not
  yet built. (Phages have no GEM by definition — they are `PhageInfection` layers,
  not reconstructions.)
- **Higher-order / adaptive integration** (see §2).
- **Diauxie / catabolite repression.** A substrate-preference layer to capture
  sequential substrate use; SOA-dFBA otherwise consumes substrates concurrently.
- **Host factors as dynamic layers.** Mucin secretion as a renewable nutrient,
  antimicrobial-peptide pulses, bile flow / peristalsis as time-varying
  influx/washout, IgA. Each fits the ecology-layer pattern.
- **Time-varying / multi-drug perturbations** beyond single-antibiotic PK.

---

## 4. Out of scope for this paradigm — use a different tool

These are **not** approximations we are choosing to skip; they are processes a
metabolism-ODE model cannot represent without becoming a different kind of model.
We do not fake them. Where relevant, couple muODE to the appropriate framework.

- **Gene regulation & signalling networks.** Transcriptional regulation, two-
  component systems, quorum sensing as *regulatory logic* (not just a diffusible
  molecule). Needs regulatory-FBA (rFBA/PROM), Boolean/ODE gene networks, or
  whole-cell models. muODE can model a signalling *metabolite's* concentration,
  not the regulatory program it triggers.
- **Evolution / strain dynamics.** Mutation, selection, horizontal gene transfer,
  pangenome flux over time. Needs population-genetic / eco-evolutionary models.
  muODE assumes fixed genomes for the run. **This explicitly includes phage
  host-range evolution and host phage-resistance evolution**: `PhageInfection`
  models the *ecological* predator–prey dynamics with a fixed host range, so a
  toy phage clears its host rather than co-evolving toward the coexistence /
  arms-race oscillation seen in vivo. It also does not model phage auxiliary
  metabolic genes reprogramming host flux.
- **Demographic stochasticity & extinction at small N.** The ODEs are
  deterministic and continuous; they cannot capture noise-driven extinction or
  founder effects at low cell numbers. Needs stochastic/agent-based simulation
  (e.g. BacArena, individual-based models). The spatial engine is structured but
  still deterministic.
- **Full 3-D biofilm mechanics.** EPS matrix, mechanical stress, detachment,
  channel formation. The 2-D reaction–diffusion engine captures gradients, not
  mechanics.
- **Host immune dynamics & physiology.** Adaptive immunity, inflammation
  feedback, epithelial turnover, motility as an organ-level process. Out of
  scope beyond the coarse "host factor" layers of §3.
- **Absolute clinical prediction.** Even with every layer on, parameters are
  literature defaults and GEMs are automated reconstructions. muODE produces
  **mechanistic, qualitative-to-semi-quantitative** hypotheses (who wins, which
  interaction drives an outcome, what a perturbation does), not validated
  patient-level predictions. Calibration against measured data is required before
  any quantitative claim.

---

## How to keep this honest

When adding a feature, put it in §1 only if it is a real mechanism with a
defensible model; if it is coarse, note it in §2; if it cannot be done in this
paradigm, it belongs in §4, not faked. Every ecology layer ships with tests that
assert the *mechanism*, and the FMT/CDI example (`examples/fmt_cdiff/`) is the
end-to-end check that the layers compose into the intended biology.
