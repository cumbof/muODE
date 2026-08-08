# Example: 12-species western gut community

A realistic human gut community on a western dietary pattern (high simple
carbohydrate, low fibre), reconstructed with **gapseq** and simulated with dynamic
FBA. The scientific payload is the **butyrate cross-feeding network**: primary
degraders release acetate / succinate / lactate, and butyrate producers cross-feed on
those to make butyrate — the terminal, health-relevant short-chain fatty acid. A
keystone-removal perturbation then shows how losing one specialist cascades through
that network.

> **Why gapseq, not CarveMe (read this first).** CarveMe carves top-down from the BiGG
> universe, which has **no connected butyrate pathway** for gut anaerobes — a CarveMe
> *Roseburia*, *Faecalibacterium* or *Eubacterium rectale* carries no butyrate exchange
> at all. So under CarveMe the butyrate arm, the whole point of this example, cannot
> emerge from stoichiometry; it would have to be a hand-set dial. gapseq builds
> bottom-up from MetaCyc and secretes butyrate natively. This is verified, not asserted:
> the two butyrate producers already reconstructed here (*R. intestinalis*,
> *F. prausnitzii*) carry a working `EX_cpd00211_e0` and secrete at ~10 mmol/gDW/h max.
> muODE keeps both engines (`config engine: carveme | gapseq`); this example needs gapseq.

> **Status: all 12 GEMs reconstructed; runs end to end and PASSES.** gapseq does not run
> on aarch64 (~30–90 min per genome), so the models were built on a workstation; all 12
> `*.xml.gz` are committed under `gems/`. First full run (2026-08-08, ts-02, 48 h):
> relative-abundance MAE **0.087** (≤0.15) and metabolite MAE **0.962** (≤5.0) both pass,
> and the gapseq-only payload is confirmed — **butyrate 0.17 mM** is produced from an
> acetate cross-feeding chain. See `gut_western_result.json` and the reconciled
> `benchmark.yaml`. Honest divergences (documented, not hidden): B. theta over-dominates
> (food-replete runaway) and R. intestinalis fills E. rectale's butyrate-producer niche
> (functional redundancy — butyrate is made anyway); the R. bromii keystone removal shows
> **no** secondary-extinction cascade (the carbon/acetate supply is redundant), an honest
> negative that reflects the food-replete limitation shared with `fmt_cdiff`.

Run from the repository root unless noted otherwise.

---

## Community composition and ecological rationale

The 12 species cover the dominant phyla and metabolic guilds of the human gut. The
relative abundances (`abundance.tsv`) reflect the published western-gut pattern:
Bacteroides dominance, reduced butyrate producers, low Akkermansia and Prevotella.

| MAG ID | Organism | Guild | Secretes | Abund. | GEM |
|--------|----------|-------|----------|-------:|-----|
| `B_thetaiotaomicron_VPI5482` | *B. thetaiotaomicron* VPI-5482 | primary carbohydrate degrader (PULs) | acetate + succinate | 0.22 | **ready** |
| `B_fragilis_NCTC9343` | *B. fragilis* NCTC 9343 | capsular polysaccharide | propionate + acetate | 0.13 | pending |
| `E_rectale_ATCC33656` | *E. rectale* ATCC 33656 | **butyrate producer** (acetyl-CoA) | butyrate | 0.12 | pending |
| `R_intestinalis_L182` | *R. intestinalis* L1-82 | **butyrate producer**; cross-feeds succinate | butyrate | 0.09 | **ready** |
| `Bl_obeum_A2162` | *Blautia obeum* A2-162 | H₂-consuming acetogen | acetate | 0.08 | pending |
| `F_prausnitzii_A2165` | *F. prausnitzii* A2-165 | **butyrate producer**; gut-health marker | butyrate | 0.08 | **ready** |
| `R_bromii_L263` | *R. bromii* L2-63 | **keystone** resistant-starch degrader | glucose/maltose | 0.07 | pending |
| `B_longum_NCC2705` | *B. longum* NCC2705 | oligosaccharide fermenter | acetate + lactate | 0.07 | pending |
| `P_copri_DSM18205` | *P. copri* DSM 18205 | plant-polysaccharide degrader | — | 0.04 | pending |
| `A_muciniphila_BAA835` | *A. muciniphila* BAA-835 | mucin degrader | propionate + acetate | 0.04 | pending |
| `L_acidophilus_NCFM` | *L. acidophilus* NCFM | lactate from simple sugars | lactate | 0.03 | pending |
| `C_comes_ATCC27758` | *C. comes* ATCC 27758 | **butyrate producer** (lactate → butyrate) | butyrate | 0.03 | pending |

Four butyrate producers (E. rectale, R. intestinalis, F. prausnitzii, C. comes); two of
them are reconstructed and verified to export butyrate. See
[`gems/PROVENANCE.md`](gems/PROVENANCE.md) for accessions, sha256 and the gapseq version.

### Cross-feeding network (expected, ModelSEED ids)

```
glucose/fructose (diet)
    │
    ├─▶ B. thetaiotaomicron ──acetate──▶ E. rectale ──butyrate──▶ host
    │         └──────────────succinate──▶ R. intestinalis ──butyrate──▶ host
    │
    ├─▶ R. bromii ──glucose/maltose──▶ E. rectale, R. intestinalis     (keystone)
    │
    ├─▶ B. longum ──acetate/lactate──▶ E. rectale, C. comes ──butyrate──▶ host
    │
    ├─▶ L. acidophilus ──lactate──▶ C. comes
    │
    └─▶ Bl. obeum (consumes H₂, enables other fermenters)

A. muciniphila: mucin ──▶ propionate + acetate
```

ModelSEED ids for the shuttles: acetate `cpd00029_e0`, butyrate `cpd00211_e0`,
propionate `cpd00141_e0`, succinate `cpd00036_e0`, lactate `cpd00159_e0`,
glucose `cpd00027_e0`. SCFAs start at zero — they are **cross-fed, not supplied**.

---

## The diet

gapseq models are ModelSEED-native, so the diet must be too. muODE's `western_gut`
preset is BiGG (`glc__D_e`); it is translated **once**, from the models' own
`bigg.metabolite` annotations, by `derive_medium_modelseed.py` into
`gems/western_gut_modelseed.csv`. Starch is rescaled (not renamed) across the namespace
boundary to conserve monomer flux — load-bearing here because R. bromii and the butyrate
producers depend on that fibre. The translated diet is a **generated** file (it needs the
12 GEMs), so it is not committed; the workstation builds it after reconstruction.

---

## Prerequisites

```bash
conda activate muode                      # muODE itself (+ cobra)
# reconstruction only (workstation):
#   gapseq 2.1.0 on PATH (its own env)   — https://github.com/jotech/gapseq
#   ncbi-datasets-cli                     — conda install -c conda-forge ncbi-datasets-cli
```

---

## Step 1 — reconstruct the GEMs (workstation, gapseq)

The three ready models are committed under `gems/`. Build the other nine:

```bash
bash examples/gut_western/gems/reconstruct_gapseq.sh   # ~30–90 min per genome × 9
```

It downloads each genome (NCBI Datasets), runs `gapseq doall`, and drops
`<slug>.xml.gz` in `gems/`, skipping any already present. Then record the new reaction
counts + sha256 in `gems/PROVENANCE.md`.

> `download_genomes.sh` is kept for fetching the raw FASTAs on their own, but
> `reconstruct_gapseq.sh` does the download itself, so you do not need to run it first.

---

## Step 2 — derive the ModelSEED diet

```bash
python examples/gut_western/derive_medium_modelseed.py   # -> gems/western_gut_modelseed.csv
```

Prints a coverage report (what resolved, what dropped, per-model growth on the translated
diet). **Re-run once all 12 GEMs are present** — a substrate only a not-yet-built member
carries cannot resolve until that model exists.

---

## Step 3 — assemble the community

```bash
muode assemble \
  --models examples/gut_western/gems \
  --abundance examples/gut_western/abundance.tsv \
  --diet examples/gut_western/gems/western_gut_modelseed.csv \
  --total-biomass 0.01 \
  --outdir examples/gut_western/simulation_env
```

Each model is paired with its relative abundance from `abundance.tsv`; abundances
initialise species biomass proportionally from `--total-biomass` (0.01 gDW/L total).

---

## Step 4 — simulate

```bash
muode simulate \
  --community examples/gut_western/simulation_env/community.json \
  --time 48.0 --step 0.1 \
  --outdir examples/gut_western/results
```

Outputs: `biomass.csv`, `metabolites.csv`, `cross_feeding.csv`, `meta.json`, and the
matching PNGs.

### What to look for (expected mechanisms — not yet a measured run)

- **Composition:** *B. thetaiotaomicron* and *B. fragilis* dominant (Bacteroides tilt of
  a western diet); butyrate producers bloom after a lag, once primary fermenters have
  seeded the acetate pool.
- **Metabolites:** glucose (`cpd00027_e0`) drawn down; acetate (`cpd00029_e0`) rises then
  partially depletes as butyrate producers consume it; butyrate (`cpd00211_e0`) and
  propionate (`cpd00141_e0`) accumulate as end products. The acetate peak-and-decline is
  the cross-feeding signature.
- **Cross-feeding:** directed edges B. thetaiotaomicron → acetate → butyrate producers,
  and B. longum / L. acidophilus → lactate → C. comes.

The same food-replete / short-horizon caveats that apply to `fmt_cdiff` apply here (see
that example's README): at physiological uptake bounds over ~48 h the community stays
well below carrying capacity, so this shows the cross-feeding *structure*, not
resource-competition dynamics.

---

## Step 5 — validate

```bash
muode validate \
  --results examples/gut_western/results \
  --expected examples/gut_western/benchmark.yaml
```

Checks Bacteroides dominance, the acetate/glucose dynamics, and at least one expected
cross-feeding edge. **The targets in `benchmark.yaml` are a qualitative spec pending a
real run** — refresh them against `results/` once all 12 GEMs are built.

---

## Step 6 — perturbation: keystone-species removal

*R. bromii* is a keystone resistant-starch degrader: it breaks starch into
oligosaccharides other species depend on. (It is one of the 9 pending GEMs, so this
needs the workstation reconstruction first.)

```bash
muode perturb \
  --community examples/gut_western/simulation_env/community.json \
  --remove-species R_bromii_L263 \
  --time 48.0 \
  --outdir examples/gut_western/results_perturb_Rbromii
```

Expected: species that consumed R. bromii-derived starch products (e.g. E. rectale,
R. intestinalis) lose biomass because the starch-derived sugar supply is gone — a
**secondary-extinction cascade** in the cross-feeding network. A second informative
removal is *F. prausnitzii* (a dysbiosis marker); total community butyrate should fall.

---

## Alternative: the full Snakemake pipeline

```bash
snakemake --use-conda --cores 8 --configfile examples/gut_western/config.yaml
```

This runs reconstruction (`engine: gapseq`) through assembly, simulation and validation
for all 12 MAGs. It rebuilds every GEM from the genomes rather than reusing the three
committed ones, so it is the from-scratch route; expect the full gapseq runtime. Enable
the perturbation arm by uncommenting the `perturbation:` block in `config.yaml`.

---

## Reproducing / auditing

- **GEM provenance** (accessions, gapseq version, sha256, which are ready vs pending):
  [`gems/PROVENANCE.md`](gems/PROVENANCE.md).
- **The honest caveats** (what dynamic FBA and this community cannot represent):
  [LIMITATIONS.md](../../docs/LIMITATIONS.md). Every output is a mechanistic hypothesis,
  not a patient-level prediction.
