# 12-species western gut community GEMs — provenance

Genome-scale models for the western-gut community, reconstructed with **gapseq**, not
CarveMe. All members share ONE namespace (ModelSEED, `cpd00211_e0`), so they cross-feed
with no id translation; the only BiGG↔ModelSEED boundary is the diet, translated once by
`../derive_medium_modelseed.py`.

## Why gapseq (and not CarveMe)

This example's scientific payload is a **butyrate cross-feeding network** — primary
degraders release acetate/succinate/lactate, and butyrate producers cross-feed on those
to make butyrate, the terminal health-relevant SCFA. CarveMe carves top-down from the
BiGG universe, which has **no connected butyrate pathway** for gut anaerobes: a CarveMe
*Roseburia* / *Faecalibacterium* / *Eubacterium rectale* carries no butyrate exchange at
all, so the butyrate arm — the whole point of the example — cannot emerge from
stoichiometry. It would have had to be a hand-set dial.

gapseq builds bottom-up from MetaCyc and, on the *same* genomes, reconstructs growing
models that secrete butyrate and acetate natively (`EX_cpd00211_e0`, `EX_cpd00029_e0`).
This is not a claim: the two butyrate producers already reconstructed here
(*R. intestinalis*, *F. prausnitzii*) carry a working butyrate export reaction and
secrete at ~10 mmol/gDW/h maximum (checked with cobra). muODE keeps both engines
(`config engine: carveme | gapseq`); this example uses gapseq because its result *rests*
on a pathway CarveMe cannot supply.

## Status: 3 of 12 reconstructed, 9 pending (workstation)

gapseq cannot run on aarch64, and reconstruction is ~30–90 min per genome, so the GEMs
are built on the workstation with `reconstruct_gapseq.sh`. Three members are already
reconstructed — they are the **identical** gapseq models used by `examples/fmt_cdiff`
(same genome accession, same gapseq/DB version, byte-identical `.xml.gz`), reused here
rather than rebuilt. The other nine still need the workstation step; until then this
example cannot be simulated end-to-end.

## How they were built

| field | value |
|---|---|
| tool | **gapseq** `2.1.0` |
| sequence DB | `1.5 (Bacteria, 2026-05-29)` |
| route | `find → transport → draft → predict medium → fill` (`gapseq doall`) |
| script | `reconstruct_gapseq.sh` (workstation; builds the 9 missing, skips the 3 present) |
| input | nucleotide genome FASTA (NCBI Datasets) |
| gap-fill medium | gapseq's own predicted medium (see note) |

> **Gap-fill medium.** `doall` gap-fills each model on gapseq's own predicted medium, not
> on `western_gut`. The simulation itself always closes every exchange the diet does not
> name (`muode.media.diet_medium`), so growth in the run comes from the diet, not from
> whatever medium the model was gap-filled against.

## Members

READY (reconstructed; sha256 of the committed `.xml.gz`):

| model | species / strain | accession | reactions | sha256 | role |
|---|---|---|---:|---|---|
| `B_thetaiotaomicron_VPI5482.xml.gz` | *Bacteroides thetaiotaomicron* VPI-5482 | GCF_000011065.1 | 1980 | `e849fd25f32fe9be` | primary carbohydrate degrader (PULs); secretes acetate + succinate |
| `R_intestinalis_L182.xml.gz` | *Roseburia intestinalis* L1-82 | GCA_900537995.1 | 1825 | `006156894c8a43fb` | **butyrate producer** (verified export); cross-feeds on succinate |
| `F_prausnitzii_A2165.xml.gz` | *Faecalibacterium prausnitzii* A2-165 | GCA_002734145.1 | 1559 | `9632d4705aa8832a` | **butyrate producer** (verified export); gut-health marker |

> These three are byte-identical to `examples/fmt_cdiff/gems/` (same accession, same
> gapseq 2.1.0 + DB 1.5). The accessions above are the ones the committed GEMs were
> **actually** built from; `download_genomes.sh` is aligned to them.

PENDING (accession fixed; `reconstruct_gapseq.sh` builds these on the workstation):

| model | species / strain | accession | role |
|---|---|---|---|
| `B_fragilis_NCTC9343.xml.gz` | *Bacteroides fragilis* NCTC 9343 | GCF_000025985.1 | capsular polysaccharide; secretes propionate + acetate |
| `E_rectale_ATCC33656.xml.gz` | *Eubacterium rectale* ATCC 33656 | GCF_000020605.1 | **butyrate producer** (acetyl-CoA route); cross-feeds acetate |
| `R_bromii_L263.xml.gz` | *Ruminococcus bromii* L2-63 | GCF_000154245.1 | keystone resistant-starch degrader; releases glucose/maltose |
| `Bl_obeum_A2162.xml.gz` | *Blautia obeum* A2-162 | GCF_000154405.1 | H₂-consuming acetogen |
| `C_comes_ATCC27758.xml.gz` | *Coprococcus comes* ATCC 27758 | GCF_000154325.1 | **butyrate producer** (lactate → butyrate) |
| `B_longum_NCC2705.xml.gz` | *Bifidobacterium longum* NCC2705 | GCF_000007525.1 | ferments oligosaccharides; secretes acetate + lactate |
| `A_muciniphila_BAA835.xml.gz` | *Akkermansia muciniphila* ATCC BAA-835 | GCF_000020225.1 | mucin degrader; secretes propionate + acetate |
| `P_copri_DSM18205.xml.gz` | *Prevotella copri* DSM 18205 | GCF_000155875.1 | plant-polysaccharide degrader (low in western diet) |
| `L_acidophilus_NCFM.xml.gz` | *Lactobacillus acidophilus* NCFM | GCF_000011985.1 | produces lactate from simple sugars |

Filenames are uniformly `<slug>.xml.gz`: every model here is gapseq. What a file *is* is
settled by this table and, at load time, by the namespace check — never by its name.

## Reproducing

```bash
# on a workstation with gapseq 2.1.0 + ncbi-datasets-cli on PATH:
bash examples/gut_western/gems/reconstruct_gapseq.sh     # ~30–90 min per genome × 9
python examples/gut_western/derive_medium_modelseed.py   # build the ModelSEED diet
```

The sha256 above will differ if a newer gapseq or sequence DB is used — that is a
different reconstruction, and the diet and any downstream result must be re-derived
against it. When the 9 land, add their reaction counts and sha256 to the READY table.
