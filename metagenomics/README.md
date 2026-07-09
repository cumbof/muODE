# `metagenomics/` — reads → MAGs: the upstream companion to muODE

A standalone Snakemake pipeline that turns raw paired-end shotgun metagenomes
into exactly the inputs muODE consumes. It is **not part of the `muode` package**:
muODE *starts* from genomes; this *produces* them.

```
 raw FASTQ ──▶ QC ──▶ assembly ──▶ binning ──▶ QC/derep ──▶ taxonomy ──▶ abundance
                                     │                                        ▲
              ┌──────────────────────┴──────────────┐          MetaPhlAn4 ────┤
       eukaryotes (Tiara→EukCC)          viruses (geNomad→CheckV)   (profile) │
                                                                              ▼
                                              muode_inputs/<sample>/
                                                mags/                 → muODE mags_dir
                                                domains.tsv           → muODE traits:
                                                abundance.tsv         → muODE abundance
                                                catalogue_quality.tsv
                                                reconciliation.tsv
                                                MANIFEST.md
```

It recovers **bacteria + archaea** (co-classified by GTDB-Tk), optionally
**eukaryotes** (Tiara → MetaBAT2 → EukCC), and detects **viruses/phages +
plasmids** (geNomad → CheckV).

**Samples are never pooled.** Each sample gets its own complete, self-contained
muODE input bundle, built only from genomes assembled out of its own reads. In a
case/control design the samples are different subjects, and a genome recovered
from one subject has no business appearing in another's simulated community.
Dereplication therefore runs *within* each sample. This is what lets you take one
rCDI bundle as the recipient community and one healthy bundle as the donor, and
inject the donor into the recipient (the `examples/fmt_cdiff/` experiment).

> **Why in this repo, not nf-core/mag or ATLAS?** Mature read→MAG pipelines exist
> (nf-core/mag, ATLAS, MetaWRAP, SqueezeMeta) and are excellent for the bacterial/
> archaeal core — if you already run one, keep it and just format its outputs like
> the hand-off below. This pipeline exists because it (a) adds first-class
> eukaryote + virus tracks that most of those under-serve, and (b) emits muODE's
> exact per-sample `mags/` + `domains.tsv` + `abundance.tsv` contract so the
> hand-off is zero-friction. It deliberately mirrors `../workflow/`'s conventions
> (per-tool conda envs in `envs/`, a SLURM profile in `profiles/slurm/`).

---

## Pipeline stages and tools

| Phase | Rule(s) | Tool | Notes |
|---|---|---|---|
| 0 · QC | `fastp`, `host_removal` | fastp, bowtie2 | adapter/quality trim; optional host read removal |
| 1 · Assembly | `assemble`, `map_reads` | **MEGAHIT** / metaSPAdes, seqkit, bowtie2 | per-sample; MEGAHIT default (low-mem) |
| 2 · Binning | `metabat2`, `maxbin2`, `concoct`, `das_tool` | MetaBAT2 + MaxBin2 + CONCOCT → **DAS_Tool** | 3 binners, consolidated |
| 3 · Prok QC | `checkm2`, `select_prok_mags` | CheckM2 | completeness/contamination gate (bacteria **and** archaea) |
| 4 · Dereplicate | `dereplicate` | dRep | **within-sample** non-redundant catalogue |
| 5 · Taxonomy | `gtdbtk` | GTDB-Tk | bacteria **and** archaea in one pass; per sample |
| 6 · Eukaryotes *(opt.)* | `tiara_classify`, `euk_bin`, `eukcc`, `select_euk_mags` | Tiara, MetaBAT2, EukCC | harder, lower-yield — see caveat |
| 7a · Viruses *(opt.)* | `genomad`, `checkv`, `virus_summary` | geNomad, CheckV | viral + plasmid; taxonomy + completeness |
| 7b · Profiling *(opt.)* | `metaphlan`, `metaphlan_sgb2gtdb_table` | MetaPhlAn4 | SGB-level quantitative profile, incl. unknown fraction |
| 8 · Abundance | `coverm_abundance`, `reconcile_abundance` | CoverM + (stdlib) | reconstructed ∩ profiled → `(mag_id, rel_abundance)` |
| — · Hand-off | `muode_bundle`, `muode_manifest`, `reconciliation_summary` | (stdlib) | one bundle per sample |

Per-sample assembly (not co-assembly) is the default: it scales linearly across a
large cohort and keeps strain variation resolvable. dRep runs **within** each
sample, collapsing near-identical bins produced by the three-binner ensemble —
it does **not** pool samples. Pooling would be wrong twice over: a subject's
community would gain genomes it never yielded, and because dRep keeps only one
representative per 0.95-ANI cluster, one subject's strain would silently stand in
for another's.

---

## Quick start

```bash
# 0. build a sample sheet from your run directory (e.g. an SRA project download)
python scripts/make_samplesheet.py \
    /home/cumbof/isilon/cumbof/git/muODE-results/PRJNA426573 \
    > config/samples.tsv

# 1. edit config/config.yaml: set the reference-DB paths for the tracks you want
#    (gtdbtk_db, genomad_db, checkv_db, and eukcc_db if recover_eukaryotes)

# 2. validate the plan without running anything (no tools needed)
snakemake -n --configfile config/config.yaml

# 3. run locally on a capable host (conda envs are created on first use)
snakemake --use-conda --cores 32 --configfile config/config.yaml

# ...or on SLURM (needs: pip install snakemake-executor-plugin-slurm)
snakemake --workflow-profile profiles/slurm --use-conda
```

Run `snakemake` **from this `metagenomics/` directory** (the `conda:`/`scripts:`
paths in the rules are relative to it).

### Key config toggles (`config/config.yaml`)

| key | default | effect |
|---|---|---|
| `assembler` | `megahit` | `metaspades` for higher-quality, higher-memory assembly |
| `remove_host` / `host_index` | `false` | bowtie2 host-read removal (set a human index for gut data) |
| `binners:` | all 3 on | enable/disable MetaBAT2 / MaxBin2 / CONCOCT |
| `min_completeness` / `max_contamination` | 50 / 10 | prokaryotic MAG quality gate |
| `derep_ani` | 0.95 | within-sample dRep clustering |
| `classify_gtdbtk` | `true` | GTDB-Tk taxonomy (needs `gtdbtk_db`); required by `run_profiling` |
| `recover_eukaryotes` | `false` | enable the Tiara→EukCC euk track (needs `eukcc_db`) |
| `detect_viruses` | `true` | geNomad→CheckV virus track (needs `genomad_db`, `checkv_db`) |
| `run_profiling` | `false` | MetaPhlAn4 profiling + reconciliation (needs `metaphlan_db`) |
| `abundance_source` | `reconciled` | `reconciled` (MetaPhlAn over matched MAGs) or `coverm` |

Toggles work from the CLI too, e.g. `--config recover_eukaryotes=true detect_viruses=false`.

---

## Reference databases

Fill the matching `*_db` path in `config/config.yaml` for each track you enable.
Rules whose DB is required but unset fail immediately with a clear message.

| DB | Used by | Get it |
|---|---|---|
| GTDB-Tk reference (`gtdbtk_db`) | `gtdbtk` | `gtdbtk download-db.sh`, or download the GTDB release and point `GTDBTK_DATA_PATH` at it |
| MetaPhlAn4 SGB DB (`metaphlan_db`) | `metaphlan` | `metaphlan --install --bowtie2db <dir>` (~20 GB; must contain a `*_SGB2GTDB.tsv`) |
| CheckM2 diamond DB (`checkm2_db`) | `checkm2` | `checkm2 database --download` (blank ⇒ CheckM2's default location) |
| geNomad DB (`genomad_db`) | `genomad` | `genomad download-database .` |
| CheckV DB (`checkv_db`) | `checkv` | `checkv download_database .` |
| EukCC2 DB (`eukcc_db`) | `eukcc` | `wget` the EukCC2 DB from the EukCC docs |
| Host bowtie2 index (`host_index`) | `host_removal` | `bowtie2-build` a human reference (e.g. GRCh38) |

---

## The muODE hand-off

Each sample gets a complete bundle at `results/muode_inputs/<sample>/` (see its
generated `MANIFEST.md`). Nothing is shared between bundles.

| file | muODE consumer |
|---|---|
| `mags/` | `../workflow` `mags_dir` (or `muode build`) — one `<mag_id>.fa` per genome |
| `domains.tsv` | `../workflow` `traits:` — `(mag_id, domain)`; routes bacteria/archaea → CarveMe/gapseq, eukaryote → MetaEuk route |
| `abundance.tsv` | `muode simulate` abundance / the rCDI adapter — `(mag_id, rel_abundance)` |
| `catalogue_quality.tsv` | provenance: completeness, contamination, GTDB lineage |
| `reconciliation.tsv` | reconstructed vs profiled (below) |

Because `mag_id` is derived once (`<sample>__bin.<k>` / `<sample>__euk.<k>`) and
reused verbatim as the FASTA stem, the abundance id, and the domain id, the three
files muODE cross-references are consistent by construction — the same
"same-string-id" rule the rCDI adapter enforces. `mags/` holds *exactly* the
genomes named in `abundance.tsv`.

Then, from the muODE repo root (muode env active):

```bash
S=metagenomics/results/muode_inputs/<sample>
snakemake --snakefile workflow/Snakefile --use-conda --cores 16 \
  --config mags_dir=$S/mags mag_extension=fa \
           traits=$S/domains.tsv \
           abundance=$S/abundance.tsv
```

Because bundles are independent, a two-community experiment is just two bundles —
reconstruct each, then inject one into the other:

```bash
muode simulate --community <rCDI_sample>/community.json \
               --inject <healthy_sample>/community.json \
               --inject-time 24 --time 96 --outdir results_fmt
```

### Reconstructed vs profiled genomes

Assembly answers *"which genomes can I reconstruct?"*; a quantitative profiler
answers *"which taxa are here, and how much?"*. **These sets do not coincide**, and
pretending otherwise silently misstates what a simulation represents:

- a species can be profiled from a handful of marker reads yet **never assemble**
  (binning realistically needs ≈5–10× even coverage; low-abundance and
  high-microdiversity taxa fail);
- a novel MAG can **assemble yet carry no markers**, so the profiler misses it.

With `run_profiling: true`, MetaPhlAn4 profiles each sample (SGB-level, so
*uncharacterised* taxa are covered as uSGBs, and `--unclassified_estimation`
reports the unassignable fraction). `reconcile_abundance.py` joins the two sets on
**GTDB species** — GTDB-Tk supplies it for each MAG, and the SGB→GTDB table shipped
with the MetaPhlAn DB supplies it for each SGB — then writes muODE's `abundance.tsv`
over the **intersection**, taking the profiler's relative abundance (comparable
across samples and studies, unlike within-sample assembly coverage).

Nothing is dropped silently. Every genome and every taxon lands in
`reconciliation.tsv`:

| `status` | meaning | simulated? |
|---|---|---|
| `matched` | assembled **and** profiled | ✅ yes |
| `reconstructed_only` | assembled, not profiled (novel / below markers / DB-version skew) | ❌ no |
| `profiled_only` | profiled, never assembled | ❌ no — substitute a reference genome to include it |

`results/muode_inputs/reconciliation_summary.tsv` gives one row per sample:
how many genomes matched, **what fraction of the profiled community the simulated
genomes account for** (`pct_profiled_abundance_captured`), MetaPhlAn's unknown
fraction, and which abundance source was actually used. Report that captured
fraction — it is the honest statement of how much of each ecosystem the simulation
covers. It is a statistics table only; no genomes are pooled across samples.

If a MAG has no species-level GTDB assignment it cannot be matched, and is recorded
as `reconstructed_only`. With `run_profiling: false` the pipeline falls back to
CoverM abundance over every MAG in the sample, and the summary says so
(`muode_abundance_source`), as it does for every other fallback.

**Metabolite namespace.** CarveMe emits BiGG exchange ids; gapseq emits
ModelSEED (`cpd*****`). Whichever you choose downstream, the ecology adapters
(e.g. `examples/fmt_cdiff/dynamics/real_data_adapter.py`) must use the matching
metabolite ids or members share no pool and no cross-feeding emerges. That check
lives on the muODE side; this pipeline only produces the genomes.

**Viruses/phages** are under `results/viruses/` (geNomad contigs + CheckV
quality). They are not reconstructed into GEMs (viruses have no metabolism) —
wire them into muODE's `PhageInfection` ecology layer (host, burst size, etc.).

---

## Honest scope / caveats

- **Eukaryote recovery is hard.** Euk MAGs from shotgun gut metagenomes are large,
  repeat-rich, and usually low-coverage; expect few, partial genomes unless a
  eukaryote is genuinely abundant. The track is real and standard (Tiara → MetaBAT2
  → EukCC) but off by default.
- **The simulated community is a subset of the real one.** Only `matched` genomes
  are simulated. `pct_profiled_abundance_captured` states how much of the profiled
  community that is; it is routinely well below 100%, and it belongs in any
  write-up. Species that never assembled are absent from the dynamics entirely —
  their metabolic contribution (and any cross-feeding through them) is simply not
  modelled.
- **Species matching is a taxonomy join, not sequence identity.** A MAG is tied to
  an SGB through their shared GTDB species name. If GTDB-Tk's release and the
  MetaPhlAn DB's underlying GTDB release disagree on a renamed species, the pair
  will not match and both sides show up as "only" rows. Check
  `reconstructed_only` for taxa you expected to match before concluding a genome
  is novel.
- **Strain substitution is avoided, not solved.** Not pooling samples stops one
  subject's strain from standing in for another's. Within a sample, dRep at
  `derep_ani` still collapses ≥95% ANI bins, so co-existing strains of one species
  are represented by a single genome (see `examples/strain_competition/`).
- **This is a wrapper around established tools**, chosen for the gut-metagenome
  case; it is not a novel assembler. Swap tools freely — the rules are small.
- **Validated here:** the full Snakemake DAG builds and prunes correctly across
  configurations (prokaryote-only, +virus, +eukaryote, host-removal on/off,
  profiling on/off), the `run_profiling`-without-`classify_gtdbtk` guard fails
  fast, a fabricated-upstream dry run confirms that restructuring reruns only
  `dereplicate`/`gtdbtk`/`coverm_abundance` and the new profiling+hand-off rules,
  and the Python glue scripts (`make_samplesheet`, `select_mags`,
  `collect_quality`, `write_domains`, `reconcile_abundance`)
  are exercised — `reconcile_abundance` against synthetic fixtures covering the
  matched / reconstructed-only / profiled-only split, multi-MAG species, and every
  fallback path. **Not run here:** the heavy bioinformatics tools themselves,
  MetaPhlAn included — they need the reference DBs and a capable x86_64 host (the
  muODE dev box is aarch64). In particular the exact MetaPhlAn CLI flags and the
  `*_SGB2GTDB.tsv` filename have **not** been executed against a real DB. Do a
  `snakemake -n` dry run on your cluster first, then a single-sample real run
  before launching the whole cohort.
