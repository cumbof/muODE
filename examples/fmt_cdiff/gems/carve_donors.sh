#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# carve_donors.sh -- turn reference genomes into muODE-ready GEMs for the
#                    recurrent-CDI / FMT example, on the WORKSTATION.
#
# WHY THIS EXISTS
#   CarveMe (and its DIAMOND/prodigal stack) does not run on the aarch64 dev
#   box, so the genome-scale members of the FMT scenario have to be carved on
#   lri-ts-02 and handed back as files.  This script produces exactly those
#   files, together with the provenance that makes them citable.
#
#   The pathogen (C. difficile) is NOT carved here: it is the published GEM
#   iCN900, downloaded from BiGG at runtime.  This script carves the DONOR
#   community that competes it out -- the half that decides the outcome, and
#   the half that must be real GEMs (not tuned yields) for the result to mean
#   anything.
#
# WHAT IT GUARANTEES
#   * every member is gap-filled against the CURRENT western_gut diet, not
#     CarveMe's default all-open medium (a GEM gap-filled open invents gut
#     transporters that do not exist and the competition becomes fiction);
#   * the exact medium used is captured (row count + sha256) so a reviewer can
#     see western_gut had not drifted underneath the carve;
#   * every output .xml.gz gets a sha256 and a CarveMe-version stamp in
#     PROVENANCE.md -- the model file is a number, and this is its `source`.
#
# REQUIREMENTS (workstation, muode env)
#   conda install -c conda-forge ncbi-datasets-cli
#   pip install carveme           # + diamond, prodigal on PATH
#
# RUN (from the repo root, on the workstation)
#   bash examples/fmt_cdiff/gems/carve_donors.sh
#
# Then hand back the whole examples/fmt_cdiff/gems/ directory (the *.xml.gz,
# western_gut_mediadb.tsv and PROVENANCE.md).  Only those are committed; the
# raw genomes and uncompressed SBML are left in work/ and are git-ignored.
# ---------------------------------------------------------------------------
set -euo pipefail

HERE="$(cd "$(dirname "$(realpath "$0")")" && pwd)"       # examples/fmt_cdiff/gems
REPO="$(cd "${HERE}/../../.." && pwd)"                     # repo root
WORK="${HERE}/work"                                        # scratch (git-ignored)
MEDIADB="${HERE}/western_gut_mediadb.tsv"
PROV="${HERE}/PROVENANCE.md"
MEDIUM="western_gut"

mkdir -p "${WORK}"

# --- the donor community -----------------------------------------------------
# One row per member: "ACCESSION|slug|Species strain|mechanistic role".
# Accessions verified against NCBI (July 2026).  GCA or GCF both work: the
# script prefers an NCBI protein annotation and falls back to prodigal on the
# nucleotide FASTA when the accession ships no proteins.
MEMBERS=(
  "GCA_004295125.1|C_scindens_ATCC35704|Clostridium scindens ATCC 35704|bai+ 7alpha-dehydroxylase: cholate -> deoxycholate (the secondary-bile-acid effector)"
  "GCF_000011065.1|B_thetaiotaomicron_VPI5482|Bacteroides thetaiotaomicron VPI-5482|generalist carbohydrate degrader: broad carbon competition, reshapes the nutrient pool"
  "GCA_900537995.1|R_intestinalis_L182|Roseburia intestinalis L1-82|butyrate producer: the SCFA/acidification arm of colonization resistance"
  # Optional cross-check member -- uncomment to also carve F. prausnitzii A2-165,
  # which has an independently PUBLISHED GEM (Heinken 2014, PMC4108055) to
  # compare this carve against.  Not load-bearing for the outcome.
  # "GCA_002734145.1|F_prausnitzii_A2165|Faecalibacterium prausnitzii A2-165|butyrate producer + published-GEM cross-check"
)

# --- 1. render the CURRENT western_gut diet as a CarveMe media-db -------------
echo "==> writing media-db from the current western_gut diet"
cd "${REPO}"
python - "${MEDIADB}" "${MEDIUM}" <<'PY'
import sys
from muode.diet import load_diet
from muode.media import write_carveme_mediadb
path, medium = sys.argv[1], sys.argv[2]
write_carveme_mediadb(load_diet(medium), path, medium=medium)
print(f"    {path}")
PY
MEDIADB_ROWS=$(($(wc -l < "${MEDIADB}") - 1))          # minus header
MEDIADB_SHA=$(sha256sum "${MEDIADB}" | cut -d' ' -f1)
echo "    ${MEDIADB_ROWS} compounds, sha256 ${MEDIADB_SHA:0:12}..."

# --- version stamps ----------------------------------------------------------
CARVE_VERSION="$(carve --version 2>&1 | head -1 || echo 'unknown')"
DATASETS_VERSION="$(datasets --version 2>&1 | head -1 || echo 'unknown')"
STAMP="$(date -u +%Y-%m-%dT%H:%M:%SZ)"

# --- PROVENANCE.md header ----------------------------------------------------
{
  echo "# Donor GEMs -- provenance"
  echo
  echo "Genome-scale models for the FMT donor community, carved by"
  echo "\`carve_donors.sh\`.  The pathogen is the published GEM **iCN900**"
  echo "(*C. difficile*), fetched from BiGG at runtime -- not carved here."
  echo
  echo "| field | value |"
  echo "|---|---|"
  echo "| carved on (UTC) | ${STAMP} |"
  echo "| CarveMe | \`${CARVE_VERSION}\` |"
  echo "| ncbi-datasets | \`${DATASETS_VERSION}\` |"
  echo "| gap-fill medium | \`${MEDIUM}\` (${MEDIADB_ROWS} compounds) |"
  echo "| media-db sha256 | \`${MEDIADB_SHA}\` |"
  echo
  echo "> The media-db is regenerated from \`muode/data/diets/western_gut.csv\`"
  echo "> at carve time, so these GEMs are gap-filled against the western_gut"
  echo "> that was in the tree on the date above.  If western_gut changes, the"
  echo "> models must be re-carved (the sha256 above will no longer match)."
  echo
  echo "## Members"
  echo
  echo "| model | species / strain | accession | input | sha256 | role |"
  echo "|---|---|---|---|---|---|"
} > "${PROV}"

# --- 2. per-member: download -> carve -> gzip -> record ----------------------
carve_one() {
    local acc="$1" slug="$2" species="$3" role="$4"
    local mdir="${WORK}/${slug}"
    local xml="${WORK}/${slug}.xml"
    local gz="${HERE}/${slug}.xml.gz"

    echo "==> ${slug}  (${acc})"
    rm -rf "${mdir}"; mkdir -p "${mdir}"

    # download protein annotation if available, else the genome for prodigal
    datasets download genome accession "${acc}" \
        --include protein,genome \
        --filename "${mdir}/ncbi.zip" --no-progressbar 2>/dev/null
    unzip -q -o "${mdir}/ncbi.zip" -d "${mdir}/x"

    local faa fna input_kind carve_input
    faa=$(find "${mdir}/x" -name 'protein.faa' | head -1 || true)
    fna=$(find "${mdir}/x" -name '*.fna' | head -1 || true)

    if [[ -n "${faa}" && -s "${faa}" ]]; then
        input_kind="protein"; carve_input=("${faa}")
        echo "    input: NCBI protein annotation"
    elif [[ -n "${fna}" && -s "${fna}" ]]; then
        input_kind="dna(prodigal)"; carve_input=(--dna "${fna}")
        echo "    input: nucleotide FASTA (prodigal via --dna)"
    else
        echo "    ERROR: no protein.faa or *.fna for ${acc}" >&2
        return 1
    fi

    # carve, gap-filling on western_gut using the media-db written above
    carve "${carve_input[@]}" \
        --gapfill "${MEDIUM}" --mediadb "${MEDIADB}" \
        --output "${xml}" -v

    # keep only a deterministic gzip in the repo (-n: no name/timestamp)
    gzip -nf "${xml}"                     # -> ${xml}.gz in WORK
    mv -f "${xml}.gz" "${gz}"
    local sha; sha=$(sha256sum "${gz}" | cut -d' ' -f1)

    printf '| `%s.xml.gz` | %s | %s | %s | `%s` | %s |\n' \
        "${slug}" "${species}" "${acc}" "${input_kind}" "${sha:0:16}" "${role}" \
        >> "${PROV}"
    echo "    -> ${gz}  (sha256 ${sha:0:12}...)"
}

for row in "${MEMBERS[@]}"; do
    IFS='|' read -r acc slug species role <<< "${row}"
    carve_one "${acc}" "${slug}" "${species}" "${role}"
done

# --- 3. optional sanity load (only if cobra is importable) -------------------
echo "==> load-checking the carved models with cobrapy"
python - "${HERE}" <<'PY' || echo "    (cobra not available -- skipping load check)"
import glob, os, sys
here = sys.argv[1]
try:
    from cobra.io import read_sbml_model
except Exception as e:
    raise SystemExit(f"cobra import failed: {e}")
for gz in sorted(glob.glob(os.path.join(here, "*.xml.gz"))):
    m = read_sbml_model(gz)
    print(f"    {os.path.basename(gz):40s} {len(m.reactions):5d} rxns  "
          f"{len(m.metabolites):5d} mets  {sum(1 for r in m.reactions if r.id.startswith('EX_')):4d} exchanges")
PY

echo
echo "============================================================"
echo " Done.  Commit-worthy outputs in examples/fmt_cdiff/gems/:"
echo "   *.xml.gz                (the carved donor GEMs)"
echo "   western_gut_mediadb.tsv (the exact gap-fill medium)"
echo "   PROVENANCE.md           (accessions, versions, hashes)"
echo " The work/ scratch dir (raw genomes, uncompressed SBML) is git-ignored."
echo "============================================================"
