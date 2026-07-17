#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# reconstruct_gapseq.sh -- build the 12-species western gut community GEMs with
#                          gapseq, on the WORKSTATION.  (Replaces the CarveMe route.)
#
# WHY GAPSEQ, NOT CARVEME
#   This example's payload is a BUTYRATE cross-feeding network.  CarveMe carves
#   from the BiGG universe, which has no connected butyrate pathway for gut
#   anaerobes: a CarveMe Roseburia / F. prausnitzii / E. rectale carries NO
#   butyrate exchange, so the butyrate arm cannot emerge from stoichiometry.
#   gapseq builds bottom-up from MetaCyc; on the SAME genomes it reconstructs a
#   growing model that secretes butyrate and acetate natively (EX_cpd00211_e0,
#   EX_cpd00029_e0).  So the whole community moves to gapseq -- ONE namespace
#   (ModelSEED), so the members cross-feed with no id translation.
#
# WHAT THIS BUILDS
#   The nine members not yet gapseq-reconstructed.  Three are already done and in
#   the repo -- they are the IDENTICAL gapseq models used by examples/fmt_cdiff
#   (same accession, same gapseq/DB version), so this script SKIPS them:
#     * B_thetaiotaomicron_VPI5482   (present)
#     * R_intestinalis_L182          (present -- a butyrate producer)
#     * F_prausnitzii_A2165          (present -- a butyrate producer)
#
# GAP-FILL MEDIUM
#   gapseq's `doall` gap-fills on its own predicted medium.  The simulation always
#   closes every exchange the diet does not name (muode.media.diet_medium), so the
#   run's growth comes from western_gut, not from the gap-fill medium.  Keep the
#   per-member work dirs: their RDS intermediates make a later re-`fill` on the
#   ModelSEED western_gut cheap (the HMM `find` step is not repeated).
#
# REQUIREMENTS (workstation)
#   * gapseq 2.1.0 on PATH (its own env)
#   * ncbi-datasets-cli
#
# RUNTIME  ~30-90 min PER genome (the find step).  Nine genomes: ~a day.
#          Members whose .xml.gz is already present are skipped.
#
# RUN
#   bash examples/gut_western/gems/reconstruct_gapseq.sh
#   # then update gems/PROVENANCE.md with the reaction counts + sha256,
#   # and run derive_medium_modelseed.py to build the diet.
# ---------------------------------------------------------------------------
set -euo pipefail

HERE="$(cd "$(dirname "$(realpath "$0")")" && pwd)"      # examples/gut_western/gems
WORK="${HERE}/gapseq_work"                               # scratch + gapseq intermediates
mkdir -p "${WORK}"

# "ACCESSION|slug|role" -- the whole community, so this script IS the recipe.  The
# accessions match gems/PROVENANCE.md.  The three shared with fmt_cdiff are listed so
# the recipe is complete, but their .xml.gz is already present, so they are skipped.
MEMBERS=(
  "GCF_000011065.1|B_thetaiotaomicron_VPI5482|primary carbohydrate degrader (present)"
  "GCA_900537995.1|R_intestinalis_L182|butyrate producer (present)"
  "GCA_002734145.1|F_prausnitzii_A2165|butyrate producer (present)"
  "GCF_000025985.1|B_fragilis_NCTC9343|propionate + acetate producer"
  "GCF_000020605.1|E_rectale_ATCC33656|butyrate producer (acetyl-CoA route)"
  "GCF_000154245.1|R_bromii_L263|keystone resistant-starch degrader"
  "GCF_000154405.1|Bl_obeum_A2162|H2-consuming acetogen"
  "GCF_000154325.1|C_comes_ATCC27758|butyrate producer (lactate -> butyrate)"
  "GCF_000007525.1|B_longum_NCC2705|oligosaccharide fermenter; acetate + lactate"
  "GCF_000020225.1|A_muciniphila_BAA835|mucin degrader; propionate + acetate"
  "GCF_000155875.1|P_copri_DSM18205|plant-polysaccharide degrader"
  "GCF_000011985.1|L_acidophilus_NCFM|lactate from simple sugars"
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
echo " Done.  The community is the *.xml.gz files in examples/gut_western/gems/;"
echo " record their sha256 + reaction counts in PROVENANCE.md, then build the diet:"
echo "   python examples/gut_western/derive_medium_modelseed.py"
echo " Keep gapseq_work/ -- its RDS intermediates allow a cheap re-fill later."
echo "============================================================"
