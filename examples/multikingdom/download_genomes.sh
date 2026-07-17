#!/usr/bin/env bash
# Download the multi-kingdom example genomes from NCBI.
#
# Two bacteria + one fungus.  The phage (vB_Kpn) is intentionally NOT downloaded:
# it has no genome-scale metabolic model and is declared as a PhageInfection
# layer, not reconstructed.  See README.md.
#
# Requirements:
#   conda install -c conda-forge ncbi-datasets-cli
#
# Run from the repository root:
#   bash examples/multikingdom/download_genomes.sh

set -euo pipefail

OUTDIR="$(dirname "$(realpath "$0")")/data/mags"
mkdir -p "$OUTDIR"

download_genome() {
    local accession="$1"
    local name="$2"
    local organism="$3"
    local outfile="${OUTDIR}/${name}.fna"

    if [[ -f "$outfile" ]]; then
        printf "  [skip]  %-30s  (already exists)\n" "$name"
        return
    fi
    printf "  [fetch] %-30s  %s\n" "$name" "$accession"

    local tmpdir
    tmpdir=$(mktemp -d)
    datasets download genome accession "$accession" \
        --include genome --filename "${tmpdir}/genome.zip" --no-progressbar 2>/dev/null
    unzip -q "${tmpdir}/genome.zip" -d "${tmpdir}/extracted"
    local fasta
    fasta=$(find "${tmpdir}/extracted" -name "*.fna" | head -1)
    if [[ -z "$fasta" ]]; then
        printf "  ERROR: no FASTA for %s (%s)\n" "$name" "$accession" >&2
        rm -rf "$tmpdir"; return 1
    fi
    cp "$fasta" "$outfile"; rm -rf "$tmpdir"
    printf "  [ok]    %-30s  -> %s\n" "$name" "$outfile"
}

echo "============================================================"
echo " Multi-kingdom example — the two BACTERIAL genomes"
echo " Output: ${OUTDIR}"
echo "============================================================"

# --- bacteria → CarveMe (-u gramneg) ---------------------------------------
# These are the only genomes to fetch: the fungus uses a curated model (Yeast8, via
# fetch_yeast8.sh) and the phage has no GEM (a PhageInfection layer).
download_genome "GCF_000011065.1"  "B_thetaiotaomicron"  \
    "Bacteroides thetaiotaomicron VPI-5482 (obligate anaerobe)"
download_genome "GCF_000240185.1"  "K_pneumoniae"        \
    "Klebsiella pneumoniae HS11286 (facultative pathobiont)"

cat <<'NOTE'

------------------------------------------------------------
Reconstruction is kingdom-specific (muode's _euk_engine routing):
  * B_thetaiotaomicron, K_pneumoniae  -> CarveMe  (-u gramneg), from these genomes.
  * S_cerevisiae (fungus)             -> CarveMe CANNOT build a eukaryote.  Use the
        curated Yeast8 model, harmonized to BiGG (NO genome download needed):
          bash   examples/multikingdom/fetch_yeast8.sh
          python examples/multikingdom/harmonize_fungal_gem.py \
              --in  examples/multikingdom/data/yeast-GEM.xml \
              --out examples/multikingdom/data/eukaryote_models/S_cerevisiae_bigg.xml.gz
  * vB_Kpn (phage)                    -> no genome-scale model; declare a
        PhageInfection(host="K_pneumoniae", ...) layer instead.
------------------------------------------------------------
NOTE
