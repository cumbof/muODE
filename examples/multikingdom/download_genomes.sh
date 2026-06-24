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
echo " Multi-kingdom example genomes (2 bacteria + 1 fungus)"
echo " Output: ${OUTDIR}"
echo "============================================================"

# --- bacteria → CarveMe (-u gramneg) ---------------------------------------
download_genome "GCF_000011065.1"  "B_thetaiotaomicron"  \
    "Bacteroides thetaiotaomicron VPI-5482 (obligate anaerobe)"
download_genome "GCF_000240185.1"  "K_pneumoniae"        \
    "Klebsiella pneumoniae HS11286 (facultative pathobiont)"

# --- eukaryote → NOT CarveMe; reconstruct via a fungal route (Yeast8, etc.) -
download_genome "GCF_000182965.3"  "C_albicans"          \
    "Candida albicans SC5314 (facultative fungus)"

cat <<'NOTE'

------------------------------------------------------------
Reconstruction is kingdom-specific (muode.reconstruction_route):
  * B_thetaiotaomicron, K_pneumoniae  -> CarveMe  (-u gramneg)
  * C_albicans (fungus)               -> CarveMe CANNOT build this.
        Use a fungal route: a Yeast8-derived template or a eukaryote-aware
        reconstructor (CarveFungi / AuReMe / gapseq fungal mode), then load
        the SBML as a CobraOrganism.
  * vB_Kpn (phage)                    -> no genome-scale model; declare a
        PhageInfection(host="K_pneumoniae", ...) layer instead.
------------------------------------------------------------
NOTE
