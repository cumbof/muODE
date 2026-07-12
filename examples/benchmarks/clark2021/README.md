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
weighting, no kinetics to argue about. **A failure here is a reconstruction
failure and nothing downstream can repair it.** A GEM of *B. thetaiotaomicron*
that secretes butyrate is simply wrong. Get this right first.

The ground truth (derived from the measurements, not asserted from memory):

| | butyrate | succinate |
|---|---|---|
| **producers** | AC, CC, ER, RI | the *Bacteroides* + *Prevotella*: BC, BF, BO, BT, BU, BV, BY, PC, PJ |
| **non-producers** | every *Bacteroides*, *Bifidobacterium*, *Dorea*, … | the butyrate producers |

This cross-checks exactly against the authors' own designation, with one
exception — FP, below.

**Tier 2 — pairwise.** Two strains: is the interaction right?

**Tier 3 — assembly.** Up to 23 strains: composition and butyrate.

## Three traps, each with a test

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
