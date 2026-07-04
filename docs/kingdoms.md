# µODE — Multiple kingdoms: fungi, protists and phages

A gut metagenome is not only bacteria. It contains **archaea**, **eukaryotes**
(fungi, protists — the *mycobiome*), and **viruses** (overwhelmingly phages — the
*virome*). This document explains how muODE represents the non-bacterial members
and, just as importantly, where it draws the line and why.

The governing principle: **muODE's simulation core is paradigm-defined, not
taxon-defined.** The dynamic-FBA engine ([engine.md](engine.md)) only ever drives
the `OrganismModel` protocol (`exchange_metabolites / set_uptake_bound /
optimize`). It has no concept of taxonomy. That single fact splits the question
into two completely different cases.

---

## 1. Fungi & protists (eukaryotes) — the same paradigm

A fungus is a metabolic organism: it has a stoichiometric network, a biomass
objective and exchange reactions. A fungal genome-scale model therefore satisfies
the **same** `OrganismModel` interface as a bacterial one and simulates with the
**same** engine code — compartmentalisation (mitochondria, ER) is internal to the
GEM and invisible to the shared-pool community. **No engine change is required.**

What actually differs between kingdoms is two things, and muODE handles each:

### a. Reconstruction routing (an external-tool concern)
CarveMe — muODE's default — has universes `bacteria / archaea / grampos /
gramneg` only. It **cannot** build a eukaryotic GEM. So the kingdom must route the
reconstruction:

```python
from muode import Domain, reconstruction_route
reconstruction_route(Domain.EUKARYOTE)
# (None, 'CarveMe does NOT support eukaryotes. Use a fungal/protist route ...')
```

The recommended eukaryote route is a curated template (e.g. **Yeast8** for
ascomycetes) or a eukaryote-aware reconstructor (**CarveFungi**, **AuReMe**,
**gapseq** fungal mode). The resulting SBML loads as an ordinary `CobraOrganism`.

muODE can also build these automatically: `reconstruct_mag(engine="carvefungi")`
(fungi) or `engine="eukaryote_generic"` (other eukaryotes) run MetaEuk gene calling
followed by CarveFungi or eggNOG-mapper + ModelSEEDpy, and the Snakemake workflow
routes eukaryote MAGs to them when `euk_ref_db` is set (see
[reconstruction.md](reconstruction.md)). A curated model supplied per-MAG in
`eukaryote_models` overrides automated reconstruction.

### b. The oxygen relationship (a real dynamical difference)
The reason a mycobiome member changes community behaviour is mostly **oxygen**.
The gut lumen is anaerobic, but it is kept that way *biologically*: facultative
microbes — many yeasts, and the Enterobacteriaceae — consume the O₂ diffusing in
from the mucosa, protecting the obligate-anaerobe majority. muODE captures this
with the [`OxygenSensitivity`](#oxygensensitivity) layer.

### Traits
Kingdom and oxygen relationship are carried as optional metadata
(`muode.traits`), kept *beside* the metabolic model so the core never needs to
know about them:

```python
from muode import Community, MicrobeTraits, Domain, OxygenTolerance
comm = Community(organisms, abundances, traits={
    "C_albicans": MicrobeTraits(domain=Domain.EUKARYOTE,
                                oxygen=OxygenTolerance.FACULTATIVE),
})
comm.domain_of("C_albicans")    # Domain.EUKARYOTE
```

The abundance input stays the plain two-column TSV; traits live in a side file
(see `examples/multikingdom/traits.tsv`).

---

## 2. Phages (viruses) — a different paradigm, coupled in

A phage has **no metabolism**: no biomass objective, no exchange fluxes, nothing
for Flux Balance Analysis to optimise. Representing it as an `OrganismModel`
would be a category error, and faking it as a metabolic "species" is exactly the
kind of confident-but-wrong modelling [LIMITATIONS.md](LIMITATIONS.md) exists to
prevent.

But phages matter — they lyse hosts and reshape the community. The principled
home for them is **population/infection dynamics** (the Levin–Stewart paradigm),
*coupled* to the metabolic ODEs through the one quantity they share with the FBA
world: **host biomass**. muODE implements this as the
[`PhageInfection`](#phageinfection) ecology layer, using the same `integrate`
hook that the sporulation life cycle uses to move biomass between pools.

```python
from muode import PhageInfection, EcologyModel
eco = EcologyModel([
    PhageInfection(host="K_pneumoniae", name="vB_Kpn",
                   adsorption_rate=11.0, burst_size=60.0,
                   latent_period=0.4, decay_rate=0.1, initial_titer=0.5),
])
```

Per step it advances free phage `P` and infected biomass `I`: adsorption removes
susceptible host into the infected pool, which lyses after a latent period into a
burst of new virions. A lytic phage crashes a host bloom and amplifies; a
`lysogeny_fraction > 0` diverts adsorptions into surviving carriers (a coarse
temperate model). See the layer reference in [ecology.md](ecology.md).

### What is deliberately *not* modelled
Per the honest-scope rule, the things that belong to a *different* paradigm are
not faked (LIMITATIONS §4):

- **Phage host-range evolution** and **host resistance evolution** — these are
  evolution. The host range is fixed for the run, which is why a toy phage clears
  its host completely instead of settling into a coexistence oscillation.
- **Auxiliary metabolic genes** (phage reprogramming host flux) and explicit
  prophage induction state.

---

## 3. Archaea

Archaea (e.g. the methanogen *Methanobrevibacter smithii*) are already first-class
in the paradigm: CarveMe reconstructs them with `universe=archaea`, and they
simulate unchanged. Their distinctive metabolism — hydrogenotrophic methanogenesis
— is captured to the extent the GEM and the H₂/CO₂ pool allow; the full gas/redox
coupling that would make methanogenesis quantitative is the §3 extension in
LIMITATIONS.

---

## Worked example

`examples/multikingdom/` puts all of this together — an obligate-anaerobe
bacterium, a facultative fungus and a phage-preyed pathobiont — and shows, on
toy models that run anywhere, that removing the fungus collapses the anaerobe
(oxygen) and removing the phage unleashes the pathobiont (predation), both
emergently:

```bash
PYTHONPATH=$(git rev-parse --show-toplevel) python examples/multikingdom/mechanistic_demo.py
```

---

## Layer reference

See [ecology.md](ecology.md) for the full parameter list of `OxygenSensitivity`
and `PhageInfection`, and how they compose with the other ecology layers.
