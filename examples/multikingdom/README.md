# Multi-kingdom gut community (bacteria + fungus + phage)

A real gut metagenome is not only bacteria. This example shows how muODE handles
the **other kingdoms** — a fungus (eukaryote) and a bacteriophage (virus) —
alongside bacteria, and *why* including them changes the dynamics.

The headline point: muODE's simulation core is **paradigm-defined, not
taxon-defined**. The dynamic-FBA engine only ever drives the `OrganismModel`
interface, so a fungal genome-scale model simulates with the *same* code as a
bacterial one. What differs between kingdoms is not the engine — it is (a) how
you *reconstruct* the model and (b) a couple of physiological traits the
metabolic LP does not encode. muODE addresses both explicitly.

## The community

| Member | Kingdom | Oxygen | Role |
|--------|---------|--------|------|
| `B_thetaiotaomicron` | bacteria | obligate anaerobe | keystone fermenter |
| `C_albicans` | **eukaryote** (fungus) | facultative | oxygen scavenger |
| `K_pneumoniae` | bacteria | facultative | pathobiont (phage host) |
| `vB_Kpn` | **virus** (phage) | — | lytic predator of *K. pneumoniae* |

`abundance.tsv` keeps the plain two-column `mag_id <TAB> rel_abundance` contract.
The kingdom/oxygen annotation lives in the side file `traits.tsv` — and note the
**phage is not in either table**: it has no genome-scale model, so it is declared
as a `PhageInfection` layer, not as a community member.

## What the two kingdoms add

### Fungi / protists — same paradigm, plus oxygen
Eukaryotes are metabolic organisms; their GEM drops straight into the community.
The dynamical reason their kingdom matters is **oxygen**: the gut lumen is
anaerobic, and facultative microbes (yeasts, Enterobacteriaceae) are its main O₂
scavengers. The `OxygenSensitivity` layer (a) gates growth by O₂ tolerance
(obligate anaerobe poisoned by O₂, obligate aerobe requires it, facultative
indifferent) and (b) lets facultatives *consume* O₂, drawing it down and
**protecting the obligate-anaerobe majority**.

### Phages — a different paradigm, coupled in
A phage has no metabolism, no biomass objective, no exchange fluxes — Flux
Balance Analysis has nothing to optimise. Forcing it into the FBA community would
be a category error. Instead muODE models it where it belongs: as a
**Levin–Stewart infection ODE** (`PhageInfection`) — adsorption removes host
biomass, a latent infected pool lyses into a burst of new virions — coupled to
the metabolic world only through host biomass.

## Run the mechanistic demo (any machine)

```bash
PYTHONPATH=$(git rev-parse --show-toplevel) python examples/multikingdom/mechanistic_demo.py
```

It runs the community and two ablations. Expected (toy-model) result:

```
scenario                              Bacteroides  Candida Klebsiella O2_final
full community (fungus + phage)             4.890    2.594      0.000    0.000
no fungus (O2 not scavenged)                0.032      nan      0.000    5.648
no phage (Klebsiella unchecked)             1.022    0.176     13.817    0.000
```

- Drop the **fungus** → O₂ is no longer scavenged → the obligate-anaerobe
  *Bacteroides* collapses (cross-kingdom protection lost).
- Drop the **phage** → the *Klebsiella* pathobiont blooms unchecked and crowds
  out the commensals.

## Scaling to real genomes

`download_genomes.sh` fetches representative genomes. Reconstruction is
**kingdom-specific** (muODE's `reconstruction_route()` encodes this):

- **bacteria** (`B_thetaiotaomicron`, `K_pneumoniae`) → CarveMe (`-u gramneg`).
- **eukaryote** (`C_albicans`) → CarveMe **cannot** build it. Use a fungal route:
  a curated template (Yeast8 for ascomycetes) or a eukaryote-aware reconstructor
  (CarveFungi / AuReMe / gapseq fungal mode), then load the SBML as a
  `CobraOrganism`.
- **virus** (`vB_Kpn`) → no GEM. Declare a `PhageInfection(host="K_pneumoniae", …)`
  layer with its adsorption rate, burst size and latent period.

## What is *not* modelled

Honest scope (see `docs/LIMITATIONS.md`): phage **host-range evolution** and host
**resistance evolution** are out-of-paradigm (evolution) and are *not* faked —
the host range is fixed for the run, which is why the toy phage clears its host
completely instead of settling into a coexistence oscillation. Auxiliary
metabolic genes (phage reprogramming host metabolism) and within-host strain
structure are likewise not represented.
