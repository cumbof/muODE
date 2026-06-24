#!/usr/bin/env bash
# Download reference genomes for the muODE fecal microbiota transplant (FMT)
# example: a dysbiotic recipient with recurrent Clostridioides difficile
# infection (rCDI) and a healthy donor stool community.
#
# Requirements:
#   conda install -c conda-forge ncbi-datasets-cli
#
# Run from the repository root:
#   bash examples/fmt_cdiff/download_genomes.sh
#
# Genomes are split into two directories so they assemble as two independent
# communities (recipient and donor), which is what the transplant simulation
# composes via `muode simulate --inject`.
#
# Accession numbers were verified against NCBI RefSeq at the time of writing.
# If a download fails, look the strain up at
#   https://www.ncbi.nlm.nih.gov/datasets/genome/
# and update the table below.

set -euo pipefail

ROOT="$(dirname "$(realpath "$0")")/data/mags"
mkdir -p "${ROOT}/recipient" "${ROOT}/donor"

download_genome() {
    local accession="$1" name="$2" subdir="$3"
    local outfile="${ROOT}/${subdir}/${name}.fna"

    if [[ -f "$outfile" ]]; then
        printf "  [skip]  %-32s (already exists)\n" "$name"
        return
    fi
    printf "  [fetch] %-32s %s\n" "$name" "$accession"

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
    cp "$fasta" "$outfile"
    rm -rf "$tmpdir"
    printf "  [ok]    %-32s -> %s\n" "$name" "$outfile"
}

echo "============================================================"
echo " RECIPIENT: dysbiotic gut, recurrent C. difficile infection"
echo "============================================================"
# The pathogen.
download_genome "GCF_000009205.2"  "C_difficile_630"        recipient
# Facultative anaerobes / pathobionts that bloom after broad-spectrum
# antibiotics clear the protective commensal community.
download_genome "GCF_000250945.1"  "E_faecium_Aus0004"      recipient
download_genome "GCF_000240185.1"  "K_pneumoniae_HS11286"   recipient
download_genome "GCF_000005845.2"  "E_coli_MG1655"          recipient

echo ""
echo "============================================================"
echo " DONOR: healthy stool (colonization-resistance restorers)"
echo "============================================================"
# Bacteroidetes -- primary fibre/carbohydrate degraders; reseed the carbon
# competition that denies C. difficile its nutrient niche.
download_genome "GCF_000011065.1"  "B_thetaiotaomicron_VPI5482"  donor
download_genome "GCF_000025985.1"  "B_fragilis_NCTC9343"          donor
# Butyrate producers -- SCFA production and acetate cross-feeding.
download_genome "GCF_000154385.1"  "F_prausnitzii_A2165"          donor
download_genome "GCF_000020605.1"  "E_rectale_ATCC33656"          donor
download_genome "GCF_000209895.1"  "R_intestinalis_L182"          donor
# Keystone starch degrader + acetogen + oligosaccharide fermenter.
download_genome "GCF_000154245.1"  "R_bromii_L263"                donor
download_genome "GCF_000154405.1"  "Bl_obeum_A2162"               donor
download_genome "GCF_000007525.1"  "B_longum_NCC2705"             donor
# Mucin specialist.
download_genome "GCF_000020225.1"  "A_muciniphila_BAA835"         donor

echo ""
echo "============================================================"
printf " Done.  recipient: %d genome(s), donor: %d genome(s)\n" \
    "$(find "${ROOT}/recipient" -name '*.fna' | wc -l)" \
    "$(find "${ROOT}/donor" -name '*.fna' | wc -l)"
echo "============================================================"
