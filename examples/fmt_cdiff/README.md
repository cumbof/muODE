# Example: fecal microbiota transplant (FMT) for recurrent *C. difficile*

Simulates a **fecal microbiota transplant** as treatment for **recurrent
*Clostridioides difficile* infection (rCDI)** — the clinical setting where FMT has
its strongest evidence base (~85–90% cure rate). The question the scenario asks is
mechanistic: *when a donor community clears the pathogen, which arm of colonization
resistance is actually doing the work* — nutrient competition, SCFA acidification, or
secondary bile acids?

Everything here runs on **real genome-scale models** reconstructed with **gapseq**.
No growth yield is a dial: each member's yield is the biomass stoichiometry that came
out of its genome, so if the pathogen is out-competed it is because the
reconstructions say so. This is a deliberate rebuild of an earlier toy version whose
yields were hand-tuned to produce clearance (that version is gone; the reasoning is
recorded in `tests/test_provenance.py`).

> **This is a workstation study, not a `muode` CLI walkthrough.** The community is
> five genome-scale LPs integrated over ~4 simulated days; a full run is ~30 min per
> arm × 5 arms. The reconstruction step needs gapseq (not available on aarch64). The
> scenario *logic* — readiness, per-member verification, the ablation arms — is
> covered by `tests/test_fmt_scenario.py` and the fast subset runs anywhere.

---

## Why gapseq, and why the whole community (pathogen included)

CarveMe — muODE's **default** engine — carves top-down from the BiGG universe, which
has **no connected butyrate pathway** for gut anaerobes: the CarveMe Roseburia and
F. prausnitzii carried no butyrate exchange at all, so the SCFA/acidification arm
could not emerge from stoichiometry. It would have had to be a dial, which is exactly
what this scenario is trying to eliminate.

gapseq builds bottom-up from MetaCyc and, on the *same* genomes, reconstructs growing
models that secrete butyrate and acetate natively (`EX_cpd00211_e0`, `EX_cpd00029_e0`).
So this scenario reconstructs the **entire** community with gapseq — the four donors
**and the pathogen** — in one namespace (ModelSEED). One namespace means the members
share a metabolite pool and cross-feed with **no id translation**; the only
BiGG↔ModelSEED boundary in the whole simulation is the diet, translated once.

muODE keeps both engines (`config engine: carveme | gapseq`, CarveMe the default).
This study chooses gapseq because its result *rests* on a pathway CarveMe cannot
supply — not as a blanket policy. See [`gems/PROVENANCE.md`](gems/PROVENANCE.md).

---

## The community (5 members)

| model | organism | role | secretes |
|---|---|---|---|
| `C_difficile_630` | *Clostridioides difficile* 630 | **the pathogen**: Stickland fermenter, spore-former, vancomycin target | — |
| `B_thetaiotaomicron_VPI5482` | *Bacteroides thetaiotaomicron* | generalist carbohydrate degrader — broad carbon competition | — |
| `R_intestinalis_L182` | *Roseburia intestinalis* L1-82 | butyrate producer — the SCFA/acidification arm | butyrate |
| `F_prausnitzii_A2165` | *Faecalibacterium prausnitzii* A2-165 | butyrate producer + published-GEM cross-check | butyrate |
| `C_scindens_ATCC35704` | *Clostridium scindens* ATCC 35704 | *bai⁺* 7α-dehydroxylase (cholate → deoxycholate) — the secondary-bile-acid effector | — |

The *bai* step is an **ecology-layer** process keyed on guild membership, not an FBA
reaction, so C. scindens carries the role without needing a bile reaction in its GEM.

---

## What muODE captures (and what it does not)

| Mechanism | In muODE | How |
|---|---|---|
| **Nutrient-niche competition** — commensals consume the pathogen's sugars / Stickland pairs / succinate | **Yes (FBA core)** | shared-pool dynamic FBA over the gapseq GEMs |
| **SCFA → pH inhibition** — a dense fermenting community acidifies the lumen | **Yes** | `WeakAcidInhibition` (ecology layer) |
| **Secondary bile acids** — *bai* 7α-dehydroxylation → deoxycholate inhibits growth and blocks germination | **Yes** | `BileAcidTransform` + `BileAcidInhibition` |
| **Sporulation / germination** — a spore reservoir survives antibiotics, germinates only when bile permits | **Yes** | `SporeForming` |
| **Antibiotic PK/PD** — kills vegetative cells, not spores | **Yes** | `Antibiotic` |
| **Gene regulation, evolution, immune dynamics** | **No** | out of paradigm — see [LIMITATIONS.md](../../docs/LIMITATIONS.md) |

The ecology parameters are the provenance-registered layer defaults (germination
Km = 15.9 mM, inhibition Ki = 0.5 mM — see `muode/provenance.py`). No parameter is
overridden per-arm; the **only** knob the experiment turns is *which layers are
present*. That is what makes the ablation a clean decomposition.

---

## Files

| file | what it is |
|---|---|
| `scenario.py` | the scenario as an importable module: `readiness()`, `verify()`, `build_community()` (refuses to assemble without the pathogen), `build_scenario(fmt, ablate=…)` |
| `run.py` | the workstation driver — runs the five ablation arms, writes CSVs + `summary.json`, prints the decomposition table |
| `derive_medium_modelseed.py` | translates the BiGG `western_gut` diet into ModelSEED once, from the models' own annotations → `gems/western_gut_modelseed.csv` |
| `gems/` | the five gapseq GEMs, `reconstruct_gapseq.sh`, the ModelSEED diet, and `PROVENANCE.md` |
| `benchmark.yaml` | the FMT-success validation spec (checked once the study has run) |

---

## Prerequisites

The GEMs are committed under `gems/`, so **you do not need to reconstruct them to run
the study** — only muODE:

```bash
conda activate muode
```

To reconstruct the GEMs from genomes yourself (workstation, gapseq 2.1.0 +
`ncbi-datasets-cli` on PATH):

```bash
bash examples/fmt_cdiff/gems/reconstruct_gapseq.sh     # ~30–90 min per genome
python examples/fmt_cdiff/derive_medium_modelseed.py   # rebuild the ModelSEED diet
```

---

## Step 1 — check the community is ready and see the yields that used to be dials

```bash
python examples/fmt_cdiff/scenario.py
```

Prints per-member readiness, then a verification line for each member: its growth
rate on the ModelSEED `western_gut` (through the diet's flux bounds — a physiological
number, not the open-medium fiction) and, for the two butyrate producers, whether the
reconstruction encodes a **connected** butyrate export pathway. This is the check that
justified moving off CarveMe: the CarveMe producers would report the pathway *absent*.

---

## Step 2 — run the mechanism-decomposition study (workstation)

```bash
python examples/fmt_cdiff/run.py --outdir results/fmt
```

Five arms, differing only in the transplant and which ecology layers are active:

| arm | fmt? | layers removed | what it isolates |
|---|---|---|---|
| `no_fmt` | no | — | recurrence baseline (spores germinate in the wiped gut) |
| `fmt_full` | yes | — | the claim: donor community clears the pathogen |
| `fmt_no_bile` | yes | bile | clearance **without** the bile arm — is bile load-bearing? |
| `fmt_no_ph` | yes | pH | clearance **without** SCFA acidification |
| `fmt_competition` | yes | bile + pH | pure nutrient competition |

`run.py` refuses to start if the pathogen GEM is not ready, writes `<arm>_biomass.csv`
/ `<arm>_metabolites.csv` / `<arm>_spores.csv` and `summary.json`, and prints:

```
MECHANISM DECOMPOSITION (pathogen final biomass per arm)
  no_fmt              …   PERSISTS
  fmt_full            …   CLEARED
  fmt_no_bile         …   …
  fmt_no_ph           …   …
  fmt_competition     …   …
```

**The reading of the table is the result:**

- If `fmt_full` clears but `fmt_competition` does **not**, the layered mechanism is
  load-bearing and muODE has decomposed colonization resistance into its arms — the
  paper.
- If `fmt_competition` clears just as well, clearance is nutrient competition and the
  bile/pH story is decorative — which is an honest and publishable finding too, not a
  failure.

Either way the conclusion comes from stoichiometry plus one binary choice per layer,
with every kinetic parameter carrying its provenance. Run the subset first to sanity
check timing: `python examples/fmt_cdiff/run.py --arms fmt_full no_fmt`.

---

## Step 3 — validate against the benchmark

```bash
muode validate --results results/fmt --expected examples/fmt_cdiff/benchmark.yaml
```

Checks that the pathogen collapsed in `fmt_full`, that butyrate is restored, and that
donor cross-feeding is present. The threshold numbers in `benchmark.yaml` are refreshed
once the workstation run has produced `summary.json`.

---

## Reproducing / auditing

- **Provenance of the GEMs** (tool version, sequence DB, accessions, sha256, roles):
  [`gems/PROVENANCE.md`](gems/PROVENANCE.md).
- **Provenance of the parameters** (which are MEASURED / DERIVED / INVENTED, and which
  the conclusion is sensitive to): `muode/provenance.py` and `tests/test_provenance.py`.
- **The honest caveats** (what dynamic FBA and this scenario cannot represent):
  [LIMITATIONS.md](../../docs/LIMITATIONS.md). Treat every output as a mechanistic
  hypothesis, not a patient-level prediction.
