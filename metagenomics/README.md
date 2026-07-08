# `metagenomics/` — reads → MAGs: the upstream companion to muODE

A standalone Snakemake pipeline that turns raw paired-end shotgun metagenomes
into exactly the inputs muODE consumes. It is **not part of the `muode` package**:
muODE *starts* from genomes; this *produces* them.

```
 raw FASTQ ──▶ QC ──▶ assembly ──▶ binning ──▶ QC/derep ──▶ taxonomy ──▶ profiling
                                     │                                        │
              ┌──────────────────────┴──────────────┐                        │
       eukaryotes (Tiara→EukCC)          viruses/phages (geNomad→CheckV)      │
                                                                              ▼
                                                       muode_inputs/
                                                         mags/            → muODE mags_dir
                                                         domains.tsv      → muODE traits:
                                                         abundance/*.tsv  → muODE abundance
                                                         catalogue_quality.tsv
```

It recovers **bacteria + archaea** (co-classified by GTDB-Tk), optionally
**eukaryotes** (Tiara → MetaBAT2 → EukCC), and detects **viruses/phages +
plasmids** (geNomad → CheckV), then profiles every sample against the
dereplicated catalogue (CoverM) to emit muODE's 2-column abundance contract.

> **Why in this repo, not nf-core/mag or ATLAS?** Mature read→MAG pipelines exist
> (nf-core/mag, ATLAS, MetaWRAP, SqueezeMeta) and are excellent for the bacterial/
> archaeal core — if you already run one, keep it and just format its outputs like
> the hand-off below. This pipeline exists because it (a) adds first-class
> eukaryote + virus tracks that most of those under-serve, and (b) emits muODE's
> exact `mags/` + `domains.tsv` + per-sample `abundance/*.tsv` contract so the
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
| 4 · Dereplicate | `dereplicate` | dRep | cross-sample non-redundant catalogue |
| 5 · Taxonomy | `gtdbtk` | GTDB-Tk | bacteria **and** archaea in one pass |
| 6 · Eukaryotes *(opt.)* | `tiara_classify`, `euk_bin`, `eukcc`, `select_euk_mags` | Tiara, MetaBAT2, EukCC | harder, lower-yield — see caveat |
| 7 · Viruses *(opt.)* | `genomad`, `checkv`, `virus_summary` | geNomad, CheckV | viral + plasmid; taxonomy + completeness |
| 8 · Profiling | `coverm_abundance`, `abundance_to_muode` | CoverM | per-sample `(mag_id, rel_abundance)` |
| — · Hand-off | `muode_domains`, `muode_manifest` | (stdlib) | `mags/`, `domains.tsv`, quality, MANIFEST |

Per-sample assembly (not co-assembly) is the default: it scales linearly across a
large cohort and keeps strain variation resolvable; dRep then collapses the same
genome recovered from many samples into one catalogue entry.

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
| `derep_ani` | 0.95 | dRep clustering (0.95 ≈ species-level catalogue) |
| `classify_gtdbtk` | `true` | GTDB-Tk taxonomy (needs `gtdbtk_db`) |
| `recover_eukaryotes` | `false` | enable the Tiara→EukCC euk track (needs `eukcc_db`) |
| `detect_viruses` | `true` | geNomad→CheckV virus track (needs `genomad_db`, `checkv_db`) |

Toggles work from the CLI too, e.g. `--config recover_eukaryotes=true detect_viruses=false`.

---

## Reference databases

Fill the matching `*_db` path in `config/config.yaml` for each track you enable.
Rules whose DB is required but unset fail immediately with a clear message.

| DB | Used by | Get it |
|---|---|---|
| GTDB-Tk reference (`gtdbtk_db`) | `gtdbtk` | `gtdbtk download-db.sh`, or download the GTDB release and point `GTDBTK_DATA_PATH` at it |
| CheckM2 diamond DB (`checkm2_db`) | `checkm2` | `checkm2 database --download` (blank ⇒ CheckM2's default location) |
| geNomad DB (`genomad_db`) | `genomad` | `genomad download-database .` |
| CheckV DB (`checkv_db`) | `checkv` | `checkv download_database .` |
| EukCC2 DB (`eukcc_db`) | `eukcc` | `wget` the EukCC2 DB from the EukCC docs |
| Host bowtie2 index (`host_index`) | `host_removal` | `bowtie2-build` a human reference (e.g. GRCh38) |

---

## The muODE hand-off

Everything muODE needs lands in `results/muode_inputs/` (see its generated
`MANIFEST.md`):

| file | muODE consumer |
|---|---|
| `mags/` | `../workflow` `mags_dir` (or `muode build`) — one `<mag_id>.fa` per genome |
| `domains.tsv` | `../workflow` `traits:` — `(mag_id, domain)`; routes bacteria/archaea → CarveMe/gapseq, eukaryote → MetaEuk route |
| `abundance/<sample>.tsv` | `muode simulate` abundance / the rCDI adapter — `(mag_id, rel_abundance)` |
| `catalogue_quality.tsv` | provenance: completeness, contamination, GTDB lineage, source sample |

Because `mag_id` is derived once (`<sample>__bin.<k>` / `<sample>__euk.<k>`) and
reused verbatim as the FASTA stem, the abundance id, and the domain id, the three
files muODE cross-references are consistent by construction — the same
"same-string-id" rule the rCDI adapter enforces.

Then, from the muODE repo root (muode env active):

```bash
snakemake --snakefile workflow/Snakefile --use-conda --cores 16 \
  --config mags_dir=metagenomics/results/muode_inputs/mags mag_extension=fa \
           traits=metagenomics/results/muode_inputs/domains.tsv \
           abundance=metagenomics/results/muode_inputs/abundance/<sample>.tsv
```

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
- **This is a wrapper around established tools**, chosen for the gut-metagenome
  case; it is not a novel assembler. Swap tools freely — the rules are small.
- **Validated here:** the full Snakemake DAG builds and prunes correctly across
  configurations (prokaryote-only, +virus, +eukaryote, host-removal on/off), and
  the Python glue scripts (`make_samplesheet`, `select_mags`, `collect_quality`,
  `coverm_to_muode`, `write_domains`) are unit-exercised. **Not run here:** the
  heavy bioinformatics tools themselves — they need the reference DBs and a capable
  x86_64 host (the muODE dev box is aarch64). Do a `snakemake -n` dry run on your
  cluster first, then a single-sample real run before launching the whole cohort.
