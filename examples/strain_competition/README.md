# Strain competition — one species, three strains, one environment

*Which strain of the same species survives when they share a habitat and a
carbon source?* This example puts three **strains of a single species** (toy
*Escherichia coli*) into one well-mixed environment and lets the dynamic-FBA
engine decide the winner — and shows that the winner **depends on which kind of
competition dominates**, not on any single "fittest" strain.

muODE's core is paradigm-defined, not taxon-defined: the engine only ever drives
the `OrganismModel` interface and a shared extracellular pool, so N strains of
one species compete with the *same* code as N different species. Strain identity
lives entirely in how the members differ — accessory metabolism and interference
traits — which is exactly what distinguishes clones of one species.

## The strains

All three are the same species (bacteria, facultative — see `traits.tsv`). They
differ only in a few accessory traits:

| Strain (`mag_id`) | Glucose | Extra trait |
|---|---|---|
| `Ecoli__glc_specialist` | **best** (high yield + affinity, lowest R\*) | — |
| `Ecoli__colicinogenic` | weaker | secretes a **colicin** (immune to it) |
| `Ecoli__arabinose_user` | weak | can use a **private sugar** the others can't |

They all start at **equal biomass** on the same glucose pool — nobody is handed
the win.

## Three ways to decide a winner

### 1. Resource competition → competitive exclusion
On a single limiting sugar, the classic result (Gause's law; Tilman's R\*) is
that the strain able to grow down to the lowest residual glucose concentration
wins and suppresses the rest. Here that is the glucose specialist.

### 2. Interference competition → the metabolic winner is overturned
Add a **colicin** (`muode.Bacteriocin`): the producer strain secretes a
diffusible toxin that suppresses the specialist's growth, `f = 1/(1 + C/k_i)`.
The producer is a *weaker grower* but now wins — chemical warfare beats
exploitative efficiency. This is the strain-warfare case colicins/microcins
actually play out in the gut.

### 3. Niche differentiation → coexistence
Give one strain a **private accessory catabolic pathway** — a sugar the others
cannot touch. It persists on its own resource regardless of the glucose fight, so
competitive exclusion no longer applies and the community **coexists**. This is
Freter's nutrient-niche idea: strains coexist by not sharing a limiting resource.

## Run the mechanistic demo (any machine)

Dependency-light toy models (numpy/scipy only), so the *mechanisms* run anywhere
— including hosts where CarveMe cannot:

```bash
PYTHONPATH=$(git rev-parse --show-toplevel) python examples/strain_competition/mechanistic_demo.py
```

Expected (toy-model) result — final biomass, gDW/L:

```
scenario                                    specialist  colicin    niche
--------------------------------------------------------------------------
1. resource competition (glucose only)          19.717    2.050    1.014
2. + colicin (interference)                      0.563    8.297    3.413
3. + private niche (glucose+arabinose)           0.402    0.954   22.871
```

- **(1)** the glucose specialist dominates and drives the other two down;
- **(2)** the colicin flips it — the *weaker* producer takes over, the
  specialist collapses;
- **(3)** the niche strain lives on its private sugar and comes to dominate,
  while the two glucose competitors stay suppressed → coexistence.

All three are emergent outcomes of the shared-pool dynamics, not scripted.

> The losers do not fall to exactly zero because the diet feeds a steady glucose
> influx (a chemostat-like trickle) that supports a residual population. Add a
> dilution/washout term (or run pure batch) to see exclusion go all the way to
> extinction — the *ordering* is the point.

## Scaling to real genomes

`download_genomes.sh` fetches three real *E. coli* strain genomes; `config.yaml`
reconstructs one GEM per strain with **CarveMe** (all gram-negative bacteria, one
BiGG namespace → one shared pool) and simulates:

```bash
bash examples/strain_competition/download_genomes.sh
conda activate muode && pip install -e .
snakemake --snakefile workflow/Snakefile --use-conda --cores 8 \
          --configfile examples/strain_competition/config.yaml
```

**Prerequisite — strain-resolved genomes.** This only works if the inputs are
genuinely different strains. Do **not** source them from the dereplicated
catalogue produced by [`../../metagenomics/`](../../metagenomics/): dRep at 0.95
ANI **merges same-species strains into one representative**, erasing exactly the
variation this example is about. Get strain resolution from **reference isolate
genomes** (as here) or **strain-resolved metagenomics** (inStrain, strainGE)
upstream.

## What is / isn't modelled — honest scope

- **In the toy demo the strain traits are assigned** (who is efficient, who makes
  the colicin, who has the private sugar) to isolate each mechanism. **With real
  genomes those differences come out of the GEMs**, so *which strain wins is an
  output, not an input* — and it may differ from the toy labels.
- **Closely related strains can reconstruct to near-identical GEMs.** CarveMe's
  top-down gap-filling can smooth over accessory-genome differences, in which case
  the metabolic competition is degenerate (neutral coexistence). Verify the strain
  GEMs actually differ (reaction/gene content, auxotrophies) before trusting a
  metabolic winner; differentiate further with strain-specific kinetics
  (`muode.predict`) or accessory-gene knockouts if needed.
- **No evolution.** Colicin gene gain/loss, immunity acquisition, resistance
  mutations and horizontal transfer are out of paradigm and are *not* faked — the
  strain traits are fixed for the run (see `docs/LIMITATIONS.md`). That is why the
  toy interference case settles to a fixed winner rather than a coevolutionary
  arms race.
- **Contact-dependent killing (T6SS)** is spatial and belongs in
  `muode.spatial`, not the well-mixed `Bacteriocin` layer used here.
