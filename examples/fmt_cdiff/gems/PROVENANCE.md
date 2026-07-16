# FMT community GEMs — provenance

Genome-scale models for the whole rCDI/FMT community — **the pathogen included** —
reconstructed with **gapseq**, not CarveMe.  All five share ONE namespace
(ModelSEED, `cpd00211_e0`), so the members cross-feed with no id translation; the
only BiGG↔ModelSEED boundary in the simulation is the diet, translated once by
`../derive_medium_modelseed.py`.

## Why gapseq (and not CarveMe)

CarveMe carves top-down from the BiGG universe, which has **no connected butyrate
pathway** for gut anaerobes: the CarveMe Roseburia / F. prausnitzii carried no
butyrate exchange at all, so the SCFA/acidification arm of colonization resistance
could not emerge from stoichiometry — it would have had to be a hand-set dial.
gapseq builds bottom-up from MetaCyc and, on the *same* genomes, reconstructs
growing models that secrete butyrate and acetate natively
(`EX_cpd00211_e0`, `EX_cpd00029_e0`).  muODE keeps both engines
(`config engine: carveme | gapseq`), CarveMe the default; this study uses gapseq
because its result *rests* on a pathway CarveMe cannot supply.

## How they were built

| field | value |
|---|---|
| tool | **gapseq** `2.1.0` |
| sequence DB | `1.5 (Bacteria, 2026-05-29)` |
| route | `find → transport → draft → predict medium → fill` (`gapseq doall`) |
| script | `reconstruct_gapseq.sh` (workstation; the three not built earlier) |
| input | nucleotide genome FASTA (NCBI Datasets) |
| gap-fill medium | gapseq's own predicted medium (see note) |

> **Gap-fill medium.** `doall` gap-fills each model on gapseq's own predicted
> medium, not on `western_gut`.  We accept that here because the expensive `find`
> step is medium-independent and its RDS intermediates are cached
> (`gapseq_work/`), so the models can later be cheaply re-`fill`ed on the ModelSEED
> `western_gut` if we choose to.  The simulation itself always closes every
> exchange the diet does not name (`muode.media.diet_medium`), so growth in the run
> comes from the diet, not from whatever medium the model was gap-filled against.

## Members

| model | species / strain | accession | reactions | sha256 | role |
|---|---|---|---:|---|---|
| `C_difficile_630.gapseq.xml.gz` | *Clostridioides difficile* 630 | GCF_000009205.2 | 2051 | `4d76afe76abd57ee` | **THE PATHOGEN** — reconstructed with gapseq too, so it shares the ModelSEED pool (iCN900 stays a BiGG cross-check, not a member: no namespace mixing) |
| `B_thetaiotaomicron_VPI5482.gapseq.xml.gz` | *Bacteroides thetaiotaomicron* VPI-5482 | GCF_000011065.1 | 1980 | `e849fd25f32fe9be` | generalist carbohydrate degrader: broad carbon competition, reshapes the nutrient pool |
| `R_intestinalis_L182.xml.gz` | *Roseburia intestinalis* L1-82 | GCA_900537995.1 | 1825 | `006156894c8a43fb` | butyrate producer: the SCFA/acidification arm of colonization resistance |
| `C_scindens_ATCC35704.gapseq.xml.gz` | *Clostridium scindens* ATCC 35704 | GCA_004295125.1 | 1785 | `64472dee35c0c32a` | *bai⁺* 7α-dehydroxylase: cholate → deoxycholate (the secondary-bile-acid effector) |
| `F_prausnitzii_A2165.xml.gz` | *Faecalibacterium prausnitzii* A2-165 | GCA_002734145.1 | 1559 | `9632d4705aa8832a` | butyrate producer + published-GEM cross-check |

> The two files without the `.gapseq.` infix (Roseburia, F. prausnitzii) were the
> first two members reconstructed, before the naming convention settled; they are
> gapseq 2.1.0 models all the same.  `scenario.py` loads either name.

## Reproducing

```bash
# on a workstation with gapseq 2.1.0 + ncbi-datasets-cli on PATH:
bash examples/fmt_cdiff/gems/reconstruct_gapseq.sh     # ~30–90 min per genome
python examples/fmt_cdiff/derive_medium_modelseed.py   # rebuild the ModelSEED diet
```

The sha256 above will differ if a newer gapseq or sequence DB is used — that is a
different reconstruction, and the medium and any downstream result must be
re-derived against it.
