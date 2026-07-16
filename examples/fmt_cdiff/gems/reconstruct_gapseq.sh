#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# reconstruct_gapseq.sh -- build the remaining FMT community GEMs with gapseq,
#                          on the WORKSTATION.  (Replaces the CarveMe route.)
#
# WHY GAPSEQ, NOT CARVEME
#   CarveMe carves from the BiGG universe, which has no connected butyrate
#   pathway for gut anaerobes: the CarveMe Roseburia/F. prausnitzii carried NO
#   butyrate exchange, and the SCFA/acidification arm of colonization resistance
#   could not emerge from stoichiometry.  gapseq builds bottom-up from MetaCyc;
#   on the SAME genomes it reconstructs a growing model that secretes butyrate
#   and acetate natively (EX_cpd00211_e0, EX_cpd00029_e0).  So the whole
#   community moves to gapseq -- ONE namespace (ModelSEED), so the members
#   cross-feed with no id translation.
#
# WHAT THIS BUILDS
#   The three members not yet gapseq-reconstructed.  Roseburia L1-82 and
#   F. prausnitzii A2-165 are already done (gapseq_butyrate_test.sh) and are in
#   the repo -- do NOT redo them here.
#     * B. thetaiotaomicron VPI-5482  -- generalist carbohydrate degrader
#     * C. scindens ATCC 35704        -- bai arm (called by markers.py, not the GEM)
#     * C. difficile 630              -- THE PATHOGEN, reconstructed with gapseq
#                                        too, so it shares the ModelSEED pool.
#                                        (iCN900 stays a BiGG cross-check, not a
#                                        community member -- no namespace mixing.)
#
# GAP-FILL MEDIUM
#   gapseq's `doall` gap-fills on its own predicted medium.  We accept that here
#   because the expensive `find` step is medium-independent and gapseq caches its
#   intermediates (*-draft.RDS, *-rxnWeights.RDS, ...).  Once the ModelSEED
#   western_gut medium is finalised from the 5-model union, all five are cheaply
#   re-`fill`ed on it -- the HMM search is NOT repeated.  So DO NOT delete the
#   per-member work dirs: the RDS intermediates are needed for that re-fill.
#
# REQUIREMENTS (workstation)
#   * gapseq on PATH (its own env)
#   * ncbi-datasets-cli
#
# RUNTIME  ~30-90 min PER genome (the find step).  Five genomes: most of a day.
#          Members whose .xml.gz is already present are skipped, so this is the
#          full recipe for a fresh clone and a no-op on a complete checkout.
#
# RUN
#   bash examples/fmt_cdiff/gems/reconstruct_gapseq.sh
#   # then hand back the *.xml.gz files it drops in this dir.
# ---------------------------------------------------------------------------
set -euo pipefail

HERE="$(cd "$(dirname "$(realpath "$0")")" && pwd)"      # examples/fmt_cdiff/gems
WORK="${HERE}/gapseq_work"                               # scratch + gapseq intermediates
mkdir -p "${WORK}"

# "ACCESSION|slug|role" -- the whole community, so this script IS the recipe.  The
# accessions are the ones in PROVENANCE.md; each member already present is skipped.
MEMBERS=(
  "GCF_000011065.1|B_thetaiotaomicron_VPI5482|generalist carbohydrate degrader"
  "GCA_004295125.1|C_scindens_ATCC35704|Stickland fermenter; bai arm via markers.py"
  "GCF_000009205.2|C_difficile_630|THE PATHOGEN (gapseq; iCN900 is the cross-check)"
  "GCA_900537995.1|R_intestinalis_L182|butyrate producer (the SCFA arm)"
  "GCA_002734145.1|F_prausnitzii_A2165|butyrate producer + published-GEM cross-check"
)

for row in "${MEMBERS[@]}"; do
    IFS='|' read -r acc slug role <<< "${row}"
    mdir="${WORK}/${slug}"
    fna="${mdir}/${slug}.fna"
    gz="${HERE}/${slug}.xml.gz"
    mkdir -p "${mdir}"

    if [[ -s "${gz}" ]]; then echo "[skip] ${slug} (have ${gz})"; continue; fi
    echo "==> ${slug}  (${acc})  -- ${role}"

    # 1. nucleotide FASTA (gapseq consumes DNA, not proteins)
    if [[ ! -s "${fna}" ]]; then
        datasets download genome accession "${acc}" --include genome \
            --filename "${mdir}/ncbi.zip" --no-progressbar 2>/dev/null
        unzip -q -o "${mdir}/ncbi.zip" -d "${mdir}/x"
        cp "$(find "${mdir}/x" -name '*.fna' | head -1)" "${fna}"
    fi

    # 2. reconstruct (find -> transport -> draft -> predict medium -> fill).
    #    Run inside mdir so gapseq's intermediates land there and are kept.
    ( cd "${mdir}" && gapseq doall "${fna}" 2>&1 | tee "${slug}.gapseq.log" )

    # 3. the gap-filled SBML is <slug>.xml (NOT <slug>-draft.xml)
    xml="${mdir}/${slug}.xml"
    if [[ ! -s "${xml}" ]]; then
        echo "    ERROR: gapseq produced no ${xml} -- check ${mdir}/${slug}.gapseq.log" >&2
        exit 1
    fi
    gzip -nc "${xml}" > "${gz}"
    echo "    -> ${gz}  (sha256 $(sha256sum "${gz}" | cut -d' ' -f1 | cut -c1-12)...)"
done

echo
echo "============================================================"
echo " Done.  The community is the *.xml.gz files in examples/fmt_cdiff/gems/;"
echo " compare their sha256 against PROVENANCE.md."
echo " Keep gapseq_work/ -- its RDS intermediates are needed to cheaply"
echo " re-gap-fill on the ModelSEED western_gut medium later."
echo "============================================================"
