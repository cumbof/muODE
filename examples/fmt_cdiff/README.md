# Example: fecal microbiota transplant (FMT) for recurrent *C. difficile*

Simulates a **fecal microbiota transplant** as treatment for **recurrent
*Clostridioides difficile* infection (rCDI)** — the clinical setting where FMT has
its strongest evidence base (~85–90% cure rate). The question the scenario asks is
mechanistic: *when a donor community clears the pathogen, which arm of colonization
resistance is actually doing the work* — nutrient competition, SCFA acidification, or
secondary bile acids?

> **Result:** on this community, **secondary bile acids** (the *C. scindens* → deoxycholate
> arm; Buffie et al. 2015). The ablation arms show nutrient competition contributes
> nothing *here*, for a structural reason worth understanding before you trust it — see
> [What the study found](#what-the-study-found-bile-not-competition). Read that section's
> limitation note: this is a food-replete, well-mixed toy that can show the bile/SCFA/spore
> arms but not competition.

Everything here runs on **real genome-scale models** reconstructed with **gapseq**.
No growth yield is a dial: each member's yield is the biomass stoichiometry that came
out of its genome, so if the pathogen is out-competed it is because the
reconstructions say so. This is a deliberate rebuild of an earlier toy version whose
yields were hand-tuned to produce clearance (that version is gone; the reasoning is
recorded in `tests/test_provenance.py`).

> **This is a workstation study, not a `muode` CLI walkthrough.** The community is
> five genome-scale LPs integrated over ~4 simulated days; a full run is ~30 min per
> arm × 5 arms. The reconstruction step needs gapseq. The scenario *logic* —
> readiness, per-member verification, the ablation arms — is covered by
> `tests/test_fmt_scenario.py` and the fast subset runs anywhere.

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
| **Nutrient-niche competition** — commensals consume the pathogen's sugars / Stickland pairs / succinate | **Machinery yes; not expressed here** | shared-pool dynamic FBA over the gapseq GEMs — but this community stays food-replete over a physiological horizon, so competition does not engage (see "Why nutrient competition does *not* show" below) |
| **SCFA → pH inhibition** — a dense fermenting community acidifies the lumen | **Yes** | `WeakAcidInhibition` (ecology layer) |
| **Secondary bile acids** — *bai* 7α-dehydroxylation → deoxycholate inhibits growth and blocks germination | **Yes** | `BileAcidTransform` + `BileAcidInhibition` |
| **Sporulation / germination** — a spore reservoir survives antibiotics, germinates only when bile permits | **Yes** | `SporeForming` |
| **Antibiotic PK/PD** — kills vegetative cells, not spores | **Yes** | `Antibiotic` |
| **Gene regulation, evolution, immune dynamics** | **No** | out of paradigm |

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
| `figures.py` | draws the two figures that carry the claim, from `summary.json` + the biomass CSVs → `figures/` |
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
python examples/fmt_cdiff/run.py --core --outdir results/fmt   # the decisive four
python examples/fmt_cdiff/run.py --outdir results/fmt          # all seven
```

Seven arms, differing only in the transplant, the drug, and which ecology layers are
active:

| arm | fmt? | layers removed | what it isolates |
|---|---|---|---|
| `untreated` | no | abx | **read this first**: does the pathogen colonize at all? |
| `abx_only` | no | — | vancomycin alone — the standard of care and the rCDI baseline |
| `fmt_only` | yes | abx | can donors displace the pathogen with no drug? |
| `fmt_full` | yes | — | drug + FMT: the treatment |
| `fmt_no_bile` | yes | bile | the treatment **without** the bile arm — is bile load-bearing? |
| `fmt_no_ph` | yes | pH | the treatment **without** SCFA acidification |
| `fmt_competition` | yes | bile + pH | pure nutrient competition |

> **Why `untreated` exists.** The first genome-scale run dosed vancomycin in *every*
> arm, so its "control" was really drug-without-FMT. The drug kills at ~4.5/h against a
> pathogen growing at 0.036/h, so it cleared the infection everywhere and all five arms
> reported CLEARED — a fact about pharmacology, read as a fact about the FMT. A study of
> colonization resistance needs an arm where nothing is done at all.

`run.py` refuses to start if the pathogen GEM is not ready, writes `<arm>_biomass.csv`
/ `<arm>_metabolites.csv` / `<arm>_spores.csv` and `summary.json`, and prints the
`MECHANISM DECOMPOSITION` table.

### The outcome measure is λ, not "cleared"

`rebound_rate` — the net specific growth rate of the pathogen over the last 30h,
`λ = dln(X)/dt` — is the result:

| λ | meaning |
|---|---|
| **< 0** | the community **excludes** the pathogen: it washes out. Colonization resistance. |
| **≈ 0** | held at a steady state |
| **> 0** | still growing at `t_end` — the arm **delays** recurrence, it does not prevent it |

`cleared` is kept as a coarse flag but is **deliberately demoted**: it is an invented
threshold (0.1) crossed with a finite horizon, so it reports *when we stopped looking*.
In the first working run every treated arm read `CLEARED` while the pathogen was growing
exponentially; at the measured rates those arms cross 0.1 at t≈194h and t≈436h. Same
run, same biology, opposite verdicts — chosen by the horizon. λ has neither dial in it.

**Clearance counts spores.** `pathogen_final` is vegetative **+** spore biomass, because
sporulation is not clearance — it is exactly how *C. difficile* survives a drug course
and comes back. A metric reading only the vegetative pool scores dormancy as a cure and
gets rCDI precisely backwards. Watch the spore *fraction*: an arm that is 78% spores has
filled the recurrence reservoir, not emptied it.

**How to read the table:**

- If `untreated` does **not** show the pathogen establishing, stop — nothing else in the
  run means anything.
- Compare `fmt_competition` against `abx_only`: equal λ means the donors contribute
  nothing through nutrient competition.
- Compare `fmt_full` against `fmt_competition`: the gap is what the bile/pH arms buy.

### What the study found: bile, not competition

Core arms (`--core`, `t_end=120h`, `D=0.025`):

| arm | pathogen total | vegetative | spores | λ (1/h) | reading |
|---|---|---|---|---|---|
| `untreated` | 0.332 | 0.332 | 0.0002 | +0.0109 | colonizes and persists — the infection is real |
| `abx_only` | 0.038 | 0.037 | 0.0002 | +0.0115 | vancomycin knocks it down, but it rebounds (λ>0) |
| `fmt_full` | 0.022 | 0.004 | 0.018 | +0.0088 | **lowest burden — and 83% of it is dormant spores** |
| `fmt_competition` | 0.039 | 0.039 | 0.0002 | +0.0119 | indistinguishable from `abx_only` |

Two comparisons carry the result:

1. **`fmt_competition` ≡ `abx_only`.** With bile and pH ablated, the transplant adds
   nothing — the donors do not suppress the pathogen by eating its food.
2. **`fmt_full` is the lowest burden, and 83% of that residual is spores.** The gap
   between `fmt_full` and `fmt_competition` is the bile/pH layers, and the mechanism is
   the bile arm holding spores from germinating — the *C. scindens* → deoxycholate story
   of **Buffie et al. 2015 (Nature)**. No parameter is tuned per arm; this comes out of
   the stoichiometry plus the ecology layers.

So the defensible, literature-consistent finding is **bile-acid-mediated colonization
resistance**, not nutrient competition.

### Why nutrient competition does *not* show — and why it is not the roster

The obvious guess (the donors are poor competitors) is wrong, and so was our first fix.

The donors are fibre degraders and the ModelSEED diet had lost its starch in
translation, so we restored it (a genuine correctness fix — see
`_CURATED_POLYMER_BIGG_TO_MODELSEED` in `muode/media.py`; BiGG `starch1200` is 1200
glucose units and ModelSEED `cpd90003` is 27, so the row is rescaled by 1200/27 to
conserve monomer flux, not renamed). With fibre in the medium *B. thetaiotaomicron*
overtakes the pathogen in monoculture (0.042 vs 0.036/h) and *C. difficile* gains
nothing, because it carries no starch exchange. **But it did not move
`fmt_competition`** — and the reason is structural, not about the roster:

**Competition needs the community to be food-limited, and at physiological parameters it
never gets there within a defensible horizon.** Two facts collide:

- The diet's uptake bounds are real dietary intake, so the community grows at μ ≈
  0.03/h. From a post-antibiotic inoculum it takes **~275 h** to reach the density where
  the medium's carbon becomes limiting (~1 gDW/L) and **~500 h** to reach a real colon's
  density (20–40 gDW/L — which the diet's own 20 gDW/L provisioning correctly targets).
- The run is 120 h — about three colonic transits, already at the edge of what a
  well-mixed model can honestly represent. Over 120 h the community stays well below
  carrying capacity, so every substrate is replete and nobody competes. (This is also
  why λ > 0 in *every* arm.)

You can have physiological density **or** a defensible horizon, not both. Forcing
competition inside 120 h would require provisioning the gut for ~100× fewer microbes
than reality — an unphysical medium chosen to manufacture the result, which we do not
do. **Bile is unaffected** because it acts as a growth-rate modifier, not a resource: it
works at any density, which is exactly why it is the mechanism this example can show.

> **Limitation, stated plainly for the reader.** This five-member community, well-mixed
> and integrated over a physiological horizon, sits in the *food-replete* regime. It
> demonstrates the **bile / SCFA / spore** arms of colonization resistance, which are
> rate-modifying and act at any density — not **nutrient-competition-based** resistance,
> which needs a food-limited community this model does not reach at physiological density
> and horizon. Read `fmt_competition` as a negative control that came out negative *for a
> defensible reason*, not as a verdict on how well the donors compete.

(Some fibre rows remain untranslated on purpose: `amylose300`, `pullulan1200`, `lmn30`
have no compound in any of these GEMs, and `xylan4`/`xylan8` do but *C. difficile*
carries them too — a real gap in the medium, not a competition asymmetry.)

---

## Step 3 — draw the figures

```bash
python examples/fmt_cdiff/figures.py --results results/fmt
```

Writes into `figures/` (created by this step; the PNGs are committed once the study
has actually run, so what is in the repo is always a picture of a real run):

| figure | what it shows |
|---|---|
| `mechanism_decomposition.png` | pathogen final biomass per arm against the clearance threshold — **the headline**: read `fmt_competition` against `fmt_full` |
| `pathogen_trajectories.png` | *C. difficile* over time per arm, so *when* the arms diverge is visible (an arm that only separates after the transplant is evidence the transplant did it) |

The figures plot the CSVs and nothing else — an arm that was not run is simply absent,
and no curve is smoothed or extrapolated.

---

## Step 4 — validate against the benchmark

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
- **The honest caveats.** Dynamic FBA and this scenario produce mechanistic
  hypotheses, not patient-level predictions; treat every output as a hypothesis to test.
