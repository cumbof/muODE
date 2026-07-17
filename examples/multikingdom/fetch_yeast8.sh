#!/usr/bin/env bash
# Fetch Yeast8 -- the consensus genome-scale model of Saccharomyces cerevisiae -- for
# the multi-kingdom example's fungus.  It is the gold-standard curated fungal GEM
# (S. cerevisiae is a facultative fungus that respires, so it can scavenge O2, which is
# the fungus's role here).  C. albicans would be more gut-realistic but has no comparably
# curated model; see README.
#
#   Lu H, Li F, Sanchez BJ, et al. A consensus S. cerevisiae metabolic model Yeast8 and
#   its ecosystem for comprehensively probing cellular metabolism.
#   Nat Commun 10, 3586 (2019).  https://github.com/SysBioChalmers/yeast-GEM
#
# The model (~12 MB) is a downloaded artifact, NOT committed -- like the bacterial
# genomes.  After fetching, harmonize its exchange ids to BiGG so it shares the pool:
#   python examples/multikingdom/harmonize_fungal_gem.py \
#       --in examples/multikingdom/data/yeast-GEM.xml \
#       --out examples/multikingdom/data/eukaryote_models/S_cerevisiae_bigg.xml.gz
#
# Run from the repository root:
#   bash examples/multikingdom/fetch_yeast8.sh
set -euo pipefail

OUTDIR="$(dirname "$(realpath "$0")")/data"
OUT="${OUTDIR}/yeast-GEM.xml"
# Pin a released tag rather than a moving branch so the model is reproducible.
VERSION="v9.1.0"
URL="https://raw.githubusercontent.com/SysBioChalmers/yeast-GEM/${VERSION}/model/yeast-GEM.xml"

mkdir -p "${OUTDIR}"
if [[ -s "${OUT}" ]]; then
    echo "[skip] ${OUT} already exists"
    exit 0
fi
echo "[fetch] Yeast8 ${VERSION} -> ${OUT}"
curl -sSL -o "${OUT}" "${URL}"
printf "[ok]    %s  (%s)\n" "${OUT}" "$(du -h "${OUT}" | cut -f1)"
