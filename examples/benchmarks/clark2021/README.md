# Benchmark: Clark et al. 2021 — 1,850 measured synthetic gut communities

The reference validation for muODE. Everything is measured, and everything is known.

> Clark RL, Connors BM, Stevenson DM, Hromada SE, Hamilton JJ, Amador-Noguez D,
> Venturelli OS. **Design of synthetic human gut microbiome assembly and butyrate
> production.** *Nat Commun* **12**, 3254 (2021).
> [doi:10.1038/s41467-021-22938-y](https://doi.org/10.1038/s41467-021-22938-y) · CC-BY 4.0

## Why this dataset

| | |
|---|---|
| **Strains** | 26 sequenced human gut isolates (type strains; genomes public) |
| **Communities** | 631 distinct, richness 1–23; 1,850 wells |
| **Medium** | **DM38** — chemically defined, composition published |
| **Growth** | 48 h, anaerobic, sealed batch |
| **Measured** | acetate, butyrate, lactate, succinate (mM, HPLC) + species relative abundance (16S) + OD + pH |

Membership is known rather than inferred, the medium is published rather than
guessed, and the incubation is finite and defined. So a prediction error is a
**model error**, not an uncertainty in the input — which is the whole point, and
what a stool cohort cannot give you.

There is also a published precedent: community-scale metabolic models were
benchmarked on this same data in [*Nat Microbiol* (2024)](https://doi.org/10.1038/s41564-024-01728-4).
muODE's numbers therefore land next to a method already in the literature.

## Running it

```bash
# 1. genomes — resolved against NCBI and matched on STRAIN, not species
python examples/benchmarks/clark2021/fetch_genomes.py --outdir data/clark2021/genomes
#    read data/clark2021/genomes/manifest.tsv before continuing

# 2. reconstruct GEMs *on DM38* — the medium the cells were actually in
snakemake --snakefile workflow/Snakefile --use-conda --cores 16 \
  --config mags_dir=data/clark2021/genomes mag_extension=fna \
           diet=dm38 gapfill_media=diet outdir=results/clark2021

# 3. score
python examples/benchmarks/clark2021/run_benchmark.py \
  --models results/clark2021/models/refined \
  --outdir results/clark2021/benchmark
```

`derive_dm38.py` regenerates `muode/data/diets/dm38.csv` from the publisher's
Supplementary Data 4. It is checked in as a derived file so a referee can diff it.

## Three tiers, reported separately

Collapsing these into one accuracy number would hide *which layer* is broken,
which is the entire diagnostic value of the dataset.

**Tier 1 — monoculture phenotype.** Grow each strain alone on DM38: which
fermentation products does it secrete? No community dynamics, no abundance
weighting, no kinetics to argue about. A *false positive* here is a reconstruction
failure nothing downstream can repair — a GEM of *B. thetaiotaomicron* that secretes
butyrate is simply wrong. A *false negative on acid secretion*, however, is **not**:
max-biomass FBA is degenerate on these GEMs and routes carbon into overflow sinks
instead of the measured acids, so it under-reports secretion for a reason that is a
property of the method, not the reconstruction. See
[A known limitation](#a-known-limitation-tier-1-acid-secretion-under-alternate-optima)
— and read Tier-1 with that asymmetry in mind.

The ground truth (derived from the measurements, not asserted from memory):

| | butyrate | succinate |
|---|---|---|
| **producers** | AC, CC, ER, RI | the *Bacteroides* + *Prevotella*: BC, BF, BO, BT, BU, BV, BY, PC, PJ |
| **non-producers** | every *Bacteroides*, *Bifidobacterium*, *Dorea*, … | the butyrate producers |

This cross-checks exactly against the authors' own designation, with one
exception — FP, below.

**Tier 2 — pairwise.** Two strains: is the interaction right?

**Tier 3 — assembly.** Up to 23 strains: composition and butyrate.

## Four traps, each with a test

**Lactate is a medium component.** DM38 supplies **28.3 mM sodium lactate**, so
lactate is a substrate as well as a product. Measured endpoint lactate averages
29.8 mM, but *net production is only 1.49 mM* — **~95% of the raw signal is just
the medium**. Scored on the endpoint, a model that predicts *nothing happens*
looks accurate. All scoring therefore defaults to **net of the DM38 baseline**.

**HB is a 26th strain.** `HB` (*Holdemanella biformis*) is distinct from `BH`
(*Blautia hydrogenotrophica*), and is not in the authors' 25-species design space.
Keying on that design space silently discards every community containing it.

***F. prausnitzii* does not grow in DM38.** FP reaches median **OD 0.02** in
monoculture, so it measures as a butyrate *non*-producer despite being a canonical
producer. A model that correctly grows FP would be scored **wrong**. So
`non_growers()` makes the exclusion explicit, and `score_phenotypes()` echoes the
excluded set back inside the score — it cannot be applied quietly.

Related: FP A2-165 has since been **reclassified** as *Faecalibacterium duncaniae*
A2-165. Searching NCBI under the name in the 2021 paper returns the species
reference (M21/2) — a different strain, with different metabolism. `fetch_genomes.py`
carries a taxon override so it resolves the strain that was actually grown.

**A zero can be the medium, not the organism — and muODE now says which.** DM38
supplies ferrous iron only (`fe2_e`, no `fe3_e`). That is chemically *correct* — Fe(III)
is not stable in a reducing anaerobic broth — but iJO1366 needs *cytoplasmic* Fe(III)
and, with no oxygen to run `FEROpp`, has no route to it. It scores **0.0 on DM38 at
any Vmax**, and that zero says nothing whatsoever about *E. coli*.

Do **not** go looking for `fe3_e` specifically. Which nutrient goes missing depends
entirely on which genomes you feed in, so it is not a fact anyone can remember to
check — it has to be computed every run. `muode.qc.diagnose_no_growth` does that, and
`refine` calls it automatically for any model that fails to grow. Two new columns in
`reconstruction_summary.tsv`:

| column | meaning |
|---|---|
| `no_growth_verdict` | `medium_gap` (the medium starved it — **not** a result) or `model_cannot_grow` (broken reconstruction; the medium is exonerated) |
| `no_growth_fixed_by_any_of` | every metabolite that **alone** restores growth |

On iJO1366/DM38 that menu contains **both** `EX_fe3_e` and `EX_no_e` — and the second
one is why it is a menu. iJO1366's only anaerobic route to cytoplasmic Fe(III) besides
importing it is `FESD2s` (`4fe4s_c + no_c → 3fe4s_c + fe3_c`), in which nitric oxide
**destroys the cell's own iron–sulfur clusters** to liberate iron. FBA will happily use
it. A diagnosis that shrinks to a single minimal answer picks between the honest fix and
the damage reaction on iteration order — ours picked NO. Enumerating the whole
equivalence class puts both on the table and lets a human tell chemistry from artefact.
The tool cannot make that judgement; it can refuse to hide it.

When a `medium_gap` appears, decide **explicitly**: supply the metabolite (a departure
from the published medium — say so in the methods) or exclude the strain the way
`non_growers()` excludes FP. Never let it be scored as a failed prediction. Pinned in
`tests/test_dm38_bounds.py` and `tests/test_no_growth_diagnosis.py`.

## A known limitation: Tier-1 acid secretion under alternate optima

Tier-1 asks a monoculture GEM which fermentation acids it secretes, and in the BiGG
namespace muODE scores this **poorly — F1 0.10 for acetate and 0.00 for butyrate,
lactate and succinate** (`tier1_secretion_report.json`, n=24). That is documented
here rather than tuned away, because the cause is a property of the *method*, not of
the reconstructions.

Flux-balance analysis maximises biomass, and for these GEMs on DM38 the optimum is
massively **degenerate**: many flux distributions give the same maximal growth, and
the measured fermentation acids are only *one* of them. The per-step max-biomass LP
in the 48 h dynamic run is free to pick any, and it routes fixed carbon into whichever
overflow sink the simplex lands on — dominantly **acetaldehyde** (a toxic
intracellular intermediate, never a real bulk product) and **branched-chain acids**
(isobutyrate, 2-methylbutyrate, from amino-acid catabolism). Net production of the
measured acids collapses toward zero.

Traced for *B. adolescentis* (BA): the cells grow normally (0.01 → 0.39 gDW/L by 8 h)
and consume glucose and maltose, but dump **+45 mM acetaldehyde** and net only
**0.2 mM acetate**. Acetate flux is positive for only **7 of 481** integration steps —
early steps land on an acetate-producing vertex, then the LP flips to an
acetaldehyde-dumping vertex of identical biomass and never returns. Butyrate is a
separate, harder problem: CarveMe's BiGG universe cannot build the pathway at all, so
the four measured producers (AC, CC, ER, RI) are **structurally unreachable**.

This is not repairable by any clean, uniform, defensible constraint. Capping the
escape exchanges (acetaldehyde / BCFA / GABA / ethanol) and per-step parsimonious FBA
were both **falsified in the 48 h dynamic run** — they lower biomass or reroute to an
unwatched sink without recovering the acids; the static pFBA signal that suggested
otherwise does not survive integration. NGAM/ATPM, amino-acid uptake caps and loopless
FBA likewise do not move the products.

**Consequence for the claims.** muODE's quantitative validation rests on **Tiers 2 and
3** — community composition, relative abundance and cross-feeding — which is the
framework's actual contribution, and where the FBA degeneracy is partly broken by
inter-species competition for shared substrates. Tier-1 acid secretion is reported
with this caveat attached to every `benchmark_report.json` (the `tier1_acid_caveat`
field) and read as a **known FBA limitation**, not a muODE or reconstruction defect. A
namespace with curated fermentation and a butyrate pathway (gapseq/ModelSEED) is the
route to a positive Tier-1 acid result; it is scaffolded (`--namespace modelseed`) but
is not part of the headline claims.

## Kinetics: the endpoint bounds Vmax, it does not identify it

```bash
python examples/benchmarks/clark2021/fit_kinetics.py \
  --models results/clark2021/models/refined --outdir results/clark2021/kinetics
```

muODE's Vmax values are invented, so the obvious move is to fit them to this data.
Mostly, **you cannot** — and the calibration says so instead of pretending.

A 48 h sealed batch culture runs to substrate exhaustion. Its endpoint is therefore
set by **yield** (biomass and product per mole of substrate — i.e. the
reconstruction's stoichiometry) and not by **rate**. Raising Vmax makes the culture
finish *sooner*, not *differently*. On E. coli core in DM38, endpoint acetate moves
22.4 → 20.4 mM while Vmax ranges over 6 → 100 mmol/gDW/h: a 16× change in the
parameter for a 9% change in the observable, most of it integration noise. Hand that
objective to an optimiser and it will return a confident minimum sitting on a
plateau.

So `muode.calibrate` classifies each strain rather than fitting it blindly:

| verdict | meaning |
|---|---|
| `identified` | the endpoint really does move with Vmax, and the optimum is interior. A number is reported. |
| `bounded_below` | Vmax must exceed some value to grow/finish; above that, nothing changes. **A bound is reported and the point estimate is withheld.** |
| `unidentifiable` | Vmax moves nothing here. |
| `no_growth` | the model never grew — that indicts the reconstruction or the medium, and must *not* be papered over by tuning Vmax. |

**This is a result, not a disappointment.** If the predictions are insensitive to
Vmax above the bound, then the Clark endpoint scores are testing the *reconstructions
and their stoichiometry*, and muODE's invented kinetics are demonstrably **not**
contaminating them. That is a stronger claim than a fitted parameter would have been,
and this script is the evidence for it.

Where kinetics genuinely do bite: **competition** — who wins a shared substrate, i.e.
community composition — and any time-resolved prediction. An endpoint measures
neither. Pinning Vmax down properly needs growth curves (Clark reports none) or the
16S fractions; fitting on the fractions would burn the very data the benchmark
scores, so it is left alone.

## What to report

From `benchmark_report.json`:

- **Tier 1:** per-metabolite accuracy / F1 across strains, and *which strains are
  misclassified* — name them.
- **Tiers 2–3:** Pearson *r*, Spearman ρ, MAE (mM) and **bias** of net production
  per metabolite, plus relative-abundance MAE.
- The **excluded** strains and the **skipped** communities. A benchmark run on half
  the panel must not be mistakable for a good one, so the runner counts both.

Report the bias, not just the correlation. A model can track butyrate across
communities (high *r*) while being systematically 3× too high — and for a
therapeutic-design claim, the offset is the part that matters.
