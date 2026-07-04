# Example: fecal microbiota transplant (FMT) for recurrent *C. difficile*

This example simulates a **fecal microbiota transplant** as a treatment for
**recurrent *Clostridioides difficile* infection (rCDI)** — the clinical setting
where FMT has its strongest evidence base (~85–90% cure rate, well above
vancomycin for recurrent disease).

It is the showcase for muODE's **timed biomass injection** feature
([docs](../../docs/engine.md)): two communities are assembled independently — a
dysbiotic recipient and a healthy donor — and the donor is *transplanted* into
the recipient mid-simulation. Whether the transplant works (donor commensals
engraft and displace *C. difficile*) is an **emergent** result of metabolic
competition, not something the simulation hard-codes.

> Run all commands from the repository root. The reconstruction steps require
> CarveMe + Prodigal; on aarch64 or machines without CarveMe substitute
> `--engine stub` to exercise the *plumbing* (the biology will not be meaningful).
> The injection engine itself is covered by `tests/test_inject.py` and runs
> anywhere.

---

## The clinical scenario

1. **Antibiotics → dysbiosis.** Broad-spectrum antibiotics clear the protective
   commensal community. The vacated gut is colonized by *C. difficile* and a
   bloom of facultative anaerobes / pathobionts (Enterococcus, Enterobacteriaceae).
2. **Recurrence.** Because the protective community is gone, standard-of-care
   antibiotics knock *C. difficile* back but it rebounds — the infection recurs.
   In the simulation this is the **control** arm: the recipient community alone
   keeps *C. difficile* dominant.
3. **FMT.** A healthy donor stool community is delivered as a bolus. Donor
   commensals re-establish **colonization resistance** and *C. difficile* is
   outcompeted. This is the **treatment** arm: `simulate --inject donor.json`.

---

## What muODE captures (and what it does not)

Colonization resistance against *C. difficile* has several mechanisms. The
nutrient-competition arm is in the dynamic-FBA core; the others are now
implemented as the [ecology layer](../../docs/ecology.md) and can be switched on
(the `mechanistic_demo.py` below uses all of them):

| Mechanism | In muODE? | How |
|-----------|-----------|-----|
| **Nutrient-niche competition** — commensals consume the sugars, sugar alcohols, sialic acid and succinate *C. difficile* needs | **Yes (core)** | shared-pool dynamic FBA |
| **SCFA production → pH inhibition** — a dense fermenting community acidifies the lumen; undissociated SCFA inhibits the pathogen | **Yes** | `WeakAcidInhibition` |
| **Secondary bile acids** — donor BSH + *bai* (7α-dehydroxylation) turn primary bile acids into deoxycholate, which inhibits *C. diff* and blocks spore germination | **Yes** | `BileAcidTransform` + `BileAcidInhibition` |
| **Sporulation / germination** — a spore reservoir survives antibiotics and germinates only when bile permits | **Yes** | `SporeForming` |
| **Antibiotic PK/PD** — a dosing course that kills vegetative cells but not spores | **Yes** | `Antibiotic` |
| **Direct antagonism** — bacteriocins / antimicrobial peptides | **Yes** | `Bacteriocin` |
| **Gene-regulatory logic, evolution, demographic stochasticity, immune dynamics** | **No** | out of paradigm — see [LIMITATIONS.md](../../docs/LIMITATIONS.md) §4 |

The genome-scale walkthrough below (Steps 1–6) demonstrates the
**nutrient-competition arm** with real GEMs. The
**[mechanistic demo](#mechanistic-demo-the-full-stack)** at the end turns on the
whole ecology stack (bile acids, spores, antibiotic PK, pH) to reproduce the
recurrence-vs-cure contrast — see [LIMITATIONS.md](../../docs/LIMITATIONS.md) for
exactly what is and isn't modelled, and treat all outputs as mechanistic
hypotheses, not patient-level predictions.

---

## The two communities

### Recipient — dysbiotic gut (recurrent CDI), `recipient_abundance.tsv`

| MAG ID | Organism | Role | Rel. abundance |
|--------|----------|------|---------------|
| `C_difficile_630` | *Clostridioides difficile* 630 | **The pathogen.** Generalist: ferments simple sugars, Stickland amino-acid pairs, sialic acid, succinate | 0.50 |
| `E_faecium_Aus0004` | *Enterococcus faecium* Aus0004 | Post-antibiotic pathobiont bloom | 0.22 |
| `K_pneumoniae_HS11286` | *Klebsiella pneumoniae* HS11286 | Enterobacteriaceae bloom (facultative anaerobe) | 0.16 |
| `E_coli_MG1655` | *Escherichia coli* K-12 MG1655 | Enterobacteriaceae bloom | 0.12 |

### Donor — healthy stool, `donor_abundance.tsv`

| MAG ID | Organism | Role | Rel. abundance |
|--------|----------|------|---------------|
| `B_thetaiotaomicron_VPI5482` | *Bacteroides thetaiotaomicron* | Primary carbohydrate degrader; secretes acetate/succinate | 0.20 |
| `F_prausnitzii_A2165` | *Faecalibacterium prausnitzii* | Anti-inflammatory butyrate producer | 0.15 |
| `B_fragilis_NCTC9343` | *Bacteroides fragilis* | Carbohydrate degrader; secretes propionate | 0.14 |
| `E_rectale_ATCC33656` | *Eubacterium rectale* | Butyrate producer (cross-feeds acetate) | 0.13 |
| `R_intestinalis_L182` | *Roseburia intestinalis* | Butyrate producer | 0.10 |
| `B_longum_NCC2705` | *Bifidobacterium longum* | Ferments oligosaccharides → acetate + lactate | 0.09 |
| `R_bromii_L263` | *Ruminococcus bromii* | Keystone starch degrader | 0.08 |
| `Bl_obeum_A2162` | *Blautia obeum* | Hydrogen-consuming acetogen | 0.06 |
| `A_muciniphila_BAA835` | *Akkermansia muciniphila* | Mucin specialist | 0.05 |

The donor is essentially the healthy community from the
[`gut_western`](../gut_western/) example — the same species that establish a
dense, carbon-competitive, SCFA-producing ecosystem.

---

## Prerequisites

```bash
conda activate muode                              # muODE
conda activate carveme                            # CarveMe + Prodigal (reconstruction)
conda install -c conda-forge ncbi-datasets-cli    # genome download
```

---

## Step 0 — download genomes

```bash
bash examples/fmt_cdiff/download_genomes.sh
```

Writes 4 recipient genomes to `data/mags/recipient/` and 9 donor genomes to
`data/mags/donor/` (under `examples/fmt_cdiff/`).

---

## Step 1 — reconstruct GEMs (both communities)

```bash
# recipient
muode build --mags examples/fmt_cdiff/data/mags/recipient \
            --engine carveme \
            --outdir examples/fmt_cdiff/models/recipient_draft

# donor
muode build --mags examples/fmt_cdiff/data/mags/donor \
            --engine carveme \
            --outdir examples/fmt_cdiff/models/donor_draft
```

---

## Step 2 — refine (gap-fill + QC)

```bash
muode refine --models examples/fmt_cdiff/models/recipient_draft \
             --outdir examples/fmt_cdiff/models/recipient

muode refine --models examples/fmt_cdiff/models/donor_draft \
             --outdir examples/fmt_cdiff/models/donor
```

---

## Step 3 — assemble the two communities

The recipient and donor are assembled **separately** into two manifests. The
recipient sets the dysbiotic starting state; the donor describes the transplant
inoculum.

```bash
# recipient: dysbiotic starting community (this is the gut at t=0)
muode assemble --models examples/fmt_cdiff/models/recipient \
               --abundance examples/fmt_cdiff/recipient_abundance.tsv \
               --diet western_gut \
               --total-biomass 0.02 \
               --outdir examples/fmt_cdiff/recipient_env

# donor: the transplant inoculum (a community.json describing what is delivered)
muode assemble --models examples/fmt_cdiff/models/donor \
               --abundance examples/fmt_cdiff/donor_abundance.tsv \
               --diet western_gut \
               --total-biomass 0.05 \
               --outdir examples/fmt_cdiff/donor_env
```

`--total-biomass 0.05` for the donor sets the size of the transplant bolus
(gDW/L) — a dense inoculum, larger than the 0.02 resident load, reflecting the
high biomass delivered by an FMT.

---

## Step 4a — CONTROL: no transplant (recurrence)

Simulate the recipient community on its own for 96 h. With no protective
community, *C. difficile* keeps the niche it has — the model of a recurrence.

```bash
muode simulate --community examples/fmt_cdiff/recipient_env/community.json \
               --time 96.0 --step 0.1 \
               --outdir examples/fmt_cdiff/results_control
```

---

## Step 4b — TREATMENT: FMT at t = 24 h

Same recipient, but the donor community is **transplanted in at 24 h**. After the
bolus, the dense donor community competes for the shared carbon pool.

```bash
muode simulate --community examples/fmt_cdiff/recipient_env/community.json \
               --inject examples/fmt_cdiff/donor_env/community.json \
               --inject-time 24.0 \
               --time 96.0 --step 0.1 \
               --outdir examples/fmt_cdiff/results_fmt
```

`--inject` merges the donor `community.json` into the simulation: the 9 donor
species join the run dormant (zero biomass) and are introduced as a bolus at
`--inject-time`. From 24 h on they compete with *C. difficile* and the
pathobionts for glucose/fructose/succinate.

### What to look for

**`results_fmt/biomass.png`:**
- **0–24 h:** *C. difficile* dominates (identical to the control).
- **At 24 h:** a vertical jump as the donor biomass is injected.
- **24–96 h:** donor *Bacteroides* and butyrate producers bloom; *C. difficile*
  biomass turns over and declines as its carbon niche is consumed by the donor
  community. Compare directly against `results_control/biomass.png`, where
  *C. difficile* never declines.

**`results_fmt/metabolites.png`:**
- Glucose/fructose, flat-ish in the control once *C. diff* saturates, are drawn
  down hard after the transplant by the much denser donor community.
- **Butyrate** (`but_e`) rises only in the FMT arm — restored SCFA production is
  the metabolic signature of a healthy engrafted community.

**`results_fmt/cross_feeding.png`:** donor-driven edges appear after 24 h
(*B. thetaiotaomicron* → acetate → *E. rectale* / *F. prausnitzii*).

---

## Step 5 — validate the transplant succeeded

```bash
muode validate --results examples/fmt_cdiff/results_fmt \
               --expected examples/fmt_cdiff/benchmark.yaml
```

The benchmark checks that *C. difficile* collapsed (low final relative
abundance), donor *Bacteroides* dominates, butyrate is restored, and at least
one donor cross-feeding edge is present.

---

## Step 6 — quantify the effect (control vs treatment)

```python
import pandas as pd

ctrl = pd.read_csv("examples/fmt_cdiff/results_control/biomass.csv", index_col=0)
fmt  = pd.read_csv("examples/fmt_cdiff/results_fmt/biomass.csv",     index_col=0)

def cdiff_fraction(df):
    rel = df.div(df.sum(axis=1), axis=0)
    return rel["C_difficile_630"].iloc[-1]

print(f"C. difficile final fraction — control: {cdiff_fraction(ctrl):.3f}")
print(f"C. difficile final fraction — FMT:     {cdiff_fraction(fmt):.3f}")
# Expect: control stays high (~0.5+), FMT collapses (~0.0-0.05).
```

The contrast between the two arms is the result: the *only* difference is the
timed donor injection, and it is enough to flip the community from
*C. difficile*-dominated to commensal-dominated.

---

## Mechanistic demo — the full stack

Steps 1–6 demonstrate the **nutrient-competition arm** with real GEMs (CarveMe
required). To see the *other* mechanisms of FMT — secondary bile acids, spore
survival, antibiotic pharmacokinetics and pH — there is a self-contained demo
that runs **anywhere** (dependency-light toy models, no GEMs, aarch64-friendly):

```bash
python examples/fmt_cdiff/mechanistic_demo.py --outdir results/fmt_mechanistic
```

It runs two arms that differ only by the transplant, driving the whole
[ecology layer](../../docs/ecology.md):

- **Recurrence (antibiotic only).** A vancomycin-like course clears vegetative
  *C. difficile*, but the **spore reservoir survives** (spores aren't killed).
  With the protective guild gone, the bile pool stays germinant-rich and
  secondary-bile-acid-poor, so the spores **germinate and the infection recurs**.
- **FMT (antibiotic + transplant).** The donor competes for carbon, **acidifies**
  via SCFA (pH ≈ 4.3), and restores **7α-dehydroxylation** (cholate →
  deoxycholate). Deoxycholate both **blocks germination** and **inhibits**
  vegetative *C. difficile*, so the spores stay dormant and the pathogen is
  cleared.

Typical output: *C. difficile* final burden ≈ **6.3 gDW/L (recurrence)** vs
**≈ 0 (FMT)** — purely emergent from the layered mechanisms. The genome-scale
parameters live in `muode/scenarios.py`; the run writes `biomass.csv`,
`metabolites.csv`, `spores.csv`, `environment.csv` (pH, germination signal, drug
concentration) and figures for each arm.

This is a *mechanistic illustration with literature-default parameters*, not a
calibrated clinical model — see [LIMITATIONS.md](../../docs/LIMITATIONS.md).

### Going further: the full ecological-dynamics analysis

The [`dynamics/`](dynamics/README.md) subfolder extends this two-arm demo into a
richer four-scenario study on a defined 12-guild community: antibiotic collapse +
spore relapse, FMT rescue, rational design of a 12-member therapeutic consortium
(naive core that fails vs. full consortium that cures, with the inferred
cross-feeding network), and an in-silico autopsy of two failed-FMT bottlenecks
(donor metabolic insufficiency; bacteriophage predation). Its README also
describes how to swap the guild stand-ins for real reconstructed MAGs.

```bash
python examples/fmt_cdiff/dynamics/run_all.py
```

---

## Variations to try

- **Transplant timing.** Sweep `--inject-time` (e.g. 12, 24, 48 h). Earlier
  transplants should clear *C. difficile* faster.
- **Dose.** Lower the donor `--total-biomass` (e.g. 0.01) to model an
  underdosed / partial transplant and find the engraftment threshold.
- **Donor quality.** Drop the keystone degrader *R. bromii* or the butyrate
  producers from `donor_abundance.tsv` to test which donor functions matter for
  clearing the pathogen.
- **Antibiotic pre-treatment + FMT.** Combine `--inject` with a perturbation
  (`muode perturb` supports injections too) to knock the community back before
  transplanting.

---

## Advanced: a single self-contained manifest

Instead of `--inject`, a community manifest may carry an `injections` block, so
one `community.json` fully describes the transplant. Add the donor models to the
`models` list and append:

```json
{
  "models": ["...recipient and donor .xml paths..."],
  "abundances": { "C_difficile_630": 0.5, "...recipient only...": 0.0 },
  "diet": "western_gut",
  "total_biomass": 0.02,
  "injections": [
    {
      "name": "FMT",
      "time": 24.0,
      "total_biomass": 0.05,
      "abundances": { "B_thetaiotaomicron_VPI5482": 0.20, "...donor...": 0.05 }
    }
  ]
}
```

Then `muode simulate --community that_manifest.json` runs the same experiment.
The `--inject` two-manifest form above is usually easier to drive.
