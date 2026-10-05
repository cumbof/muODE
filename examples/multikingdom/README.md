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
| `S_cerevisiae` | **eukaryote** (fungus) | facultative | oxygen scavenger |
| `K_pneumoniae` | bacteria | facultative | pathobiont (phage host) |
| `vB_Kpn` | **virus** (phage) | — | lytic predator of *K. pneumoniae* |

> **The fungus is *S. cerevisiae*, via Yeast8.** Yeast8 (SysBioChalmers/yeast-GEM) is the
> gold-standard curated fungal GEM, and *S. cerevisiae* is a genuine gut-resident yeast
> (*S. boulardii*, a probiotic, is a strain of it). *C. albicans* would be more gut-native
> but has no comparably curated model — and the point here is the *cross-kingdom oxygen
> mechanism*, which any respiring facultative fungus exercises. See "Scaling to real
> genomes" below.

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
scenario                              Bacteroides    yeast Klebsiella O2_final
full community (fungus + phage)             4.890    2.594      0.000    0.000
no fungus (O2 not scavenged)                0.032      nan      0.000    5.648
no phage (Klebsiella unchecked)             1.022    0.176     13.817    0.000
```

- Drop the **fungus** → O₂ is no longer scavenged → the obligate-anaerobe
  *Bacteroides* collapses (cross-kingdom protection lost).
- Drop the **phage** → the *Klebsiella* pathobiont blooms unchecked and crowds
  out the commensals.

## Scaling to real genomes

Reconstruction is **kingdom-specific** (muODE's `_euk_engine` routing in the Snakefile
encodes this). The fungus is handled *here*, locally; only the two bacteria need a
reconstruction host.

**The fungus — curated, local, no gapseq/CarveFungi.** CarveMe cannot build a eukaryote.
Rather than stand up MetaEuk + CarveFungi for one organism, use the gold-standard curated
model. The only obstacle is namespace — the community shares one pool keyed by
exchange-metabolite ids, and Yeast8 speaks its own dialect (`s_0565` = glucose) — so its
exchanges are relabelled to BiGG from the `bigg.metabolite` annotations Yeast8 already
carries. Both steps run on any machine:

```bash
bash   examples/multikingdom/fetch_yeast8.sh              # Yeast8 v9.1.0 (~12 MB, not committed)
python examples/multikingdom/harmonize_fungal_gem.py \
    --in  examples/multikingdom/data/yeast-GEM.xml \
    --out examples/multikingdom/data/eukaryote_models/S_cerevisiae_bigg.xml.gz
```

`config.yaml` already points `eukaryote_models: {S_cerevisiae: …}` at that output, so the
pipeline drops it into the community in place of an automated reconstruction. Verified:
the harmonized model grows at ~0.47/h on `mucosal_aerobic.csv` and 0 without O₂ — so the
oxygen scavenging is metabolically load-bearing, not a trait bolted on.

**The bacteria — CarveMe (needs a reconstruction host).**

```bash
bash examples/multikingdom/download_genomes.sh           # B. theta + K. pneumoniae genomes
snakemake --use-conda --cores 8 --configfile examples/multikingdom/config.yaml
```

`engine: carveme -u gramneg` builds both; the pipeline skips the fungus (curated) and the
phage (no GEM). CarveMe and Yeast8 are both BiGG, so the harmonized fungus and the
bacteria cross-feed from one pool on `mucosal_aerobic.csv` — a **defined aerobic** medium
(muODE's `western_gut` preset is anaerobic; this example needs the O₂ the fungus
scavenges).

**The phage — no GEM.** Declare a `PhageInfection(host="K_pneumoniae", …)` layer with its
adsorption rate, burst size and latent period (as in `mechanistic_demo.py`).

> **Status: runs end to end (toy and genome-scale).** All three GEMs are committed (the
> two bacteria via CarveMe, the fungus via a harmonized Yeast8 model). Over a 48 h
> genome-scale simulation (`genome_scale_result.json`) both cross-kingdom mechanisms
> emerge: the phage crashes the Klebsiella pathobiont (0.507 → 0.000), and the fungus's
> O₂ respiration keeps the lumen anoxic and protects the anaerobe (~1.25× — modest and
> apt for *B. theta*'s aerotolerance). If a strain starves, `muode.qc`'s no-growth
> diagnosis names the missing nutrient — supplement it in the methods.

## What is *not* modelled

Honest scope: phage **host-range evolution** and host
**resistance evolution** are out-of-paradigm (evolution) and are *not* faked —
the host range is fixed for the run, which is why the toy phage clears its host
completely instead of settling into a coexistence oscillation. Auxiliary
metabolic genes (phage reprogramming host metabolism) and within-host strain
structure are likewise not represented.
