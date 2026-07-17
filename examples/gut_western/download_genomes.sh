#!/usr/bin/env bash
# Download 12 human gut reference genomes from NCBI RefSeq for the
# muODE western gut community example.
#
# Requirements:
#   conda install -c conda-forge ncbi-datasets-cli
#   (or: pip install ncbi-datasets-pylib)
#
# Run from the repository root:
#   bash examples/gut_western/download_genomes.sh
#
# Accession numbers were verified against NCBI RefSeq at the time of writing.
# If a download fails, verify the accession at:
#   https://www.ncbi.nlm.nih.gov/datasets/genome/
# and update the table below.

set -euo pipefail

OUTDIR="$(dirname "$(realpath "$0")")/data/mags"
mkdir -p "$OUTDIR"

# ---------------------------------------------------------------------------
# Download helper: datasets download → unzip → rename to <name>.fna
# ---------------------------------------------------------------------------
download_genome() {
    local accession="$1"
    local name="$2"
    local organism="$3"
    local outfile="${OUTDIR}/${name}.fna"

    if [[ -f "$outfile" ]]; then
        printf "  [skip]  %-35s  (already exists)\n" "$name"
        return
    fi

    printf "  [fetch] %-35s  %s\n" "$name" "$accession"

    local tmpdir
    tmpdir=$(mktemp -d)

    datasets download genome accession "$accession" \
        --include genome \
        --filename "${tmpdir}/genome.zip" \
        --no-progressbar 2>/dev/null

    unzip -q "${tmpdir}/genome.zip" -d "${tmpdir}/extracted"

    local fasta
    fasta=$(find "${tmpdir}/extracted" -name "*.fna" | head -1)

    if [[ -z "$fasta" ]]; then
        printf "  ERROR: no FASTA found for %s (%s)\n" "$name" "$accession" >&2
        rm -rf "$tmpdir"
        return 1
    fi

    cp "$fasta" "$outfile"
    rm -rf "$tmpdir"
    printf "  [ok]    %-35s  -> %s\n" "$name" "$outfile"
}

echo "============================================================"
echo " Downloading 12 human gut reference genomes"
echo " Output: ${OUTDIR}"
echo "============================================================"
echo ""

# ---------------------------------------------------------------------------
# Bacteroidetes — dominant in western gut; ferment complex carbohydrates,
# produce acetate, succinate and propionate.
# ---------------------------------------------------------------------------
download_genome "GCF_000011065.1"  "B_thetaiotaomicron_VPI5482"  \
    "Bacteroides thetaiotaomicron VPI-5482"
download_genome "GCF_000025985.1"  "B_fragilis_NCTC9343"           \
    "Bacteroides fragilis NCTC 9343"

# ---------------------------------------------------------------------------
# Firmicutes — butyrate producers; cross-feed on acetate / succinate /
# lactate from primary fermenters to produce butyrate.
# ---------------------------------------------------------------------------
download_genome "GCF_000020605.1"  "E_rectale_ATCC33656"           \
    "Eubacterium rectale ATCC 33656"
# NB: the two accessions below are the ones the committed gapseq GEMs were ACTUALLY
# built from (shared byte-identical with examples/fmt_cdiff; see gems/PROVENANCE.md).
# They differ from the RefSeq assemblies an earlier draft of this script used, so the
# recipe stays consistent with the reconstruction artifact.
download_genome "GCA_900537995.1"  "R_intestinalis_L182"           \
    "Roseburia intestinalis L1-82"
download_genome "GCA_002734145.1"  "F_prausnitzii_A2165"           \
    "Faecalibacterium prausnitzii A2-165"
download_genome "GCF_000154245.1"  "R_bromii_L263"                 \
    "Ruminococcus bromii L2-63"
download_genome "GCF_000154405.1"  "Bl_obeum_A2162"                \
    "Blautia obeum A2-162"
download_genome "GCF_000154325.1"  "C_comes_ATCC27758"             \
    "Coprococcus comes ATCC 27758"

# ---------------------------------------------------------------------------
# Actinobacteriota — ferments oligosaccharides (GOS, FOS); produces acetate
# and lactate that feed butyrate producers downstream.
# ---------------------------------------------------------------------------
download_genome "GCF_000007525.1"  "B_longum_NCC2705"              \
    "Bifidobacterium longum NCC2705"

# ---------------------------------------------------------------------------
# Verrucomicrobiota — mucin layer degrader; produces propionate and acetate;
# abundance is reduced in western diet (high-fat diet).
# ---------------------------------------------------------------------------
download_genome "GCF_000020225.1"  "A_muciniphila_BAA835"          \
    "Akkermansia muciniphila ATCC BAA-835"

# ---------------------------------------------------------------------------
# Bacteroidetes (Prevotellaceae) — plant polysaccharide degrader;
# typically low in western / fibre-poor diets.
# ---------------------------------------------------------------------------
download_genome "GCF_000155875.1"  "P_copri_DSM18205"              \
    "Prevotella copri DSM 18205"

# ---------------------------------------------------------------------------
# Lactobacillales — produces lactate from simple sugars; trace abundance in
# healthy western gut, important cross-feeder for Coprococcus.
# ---------------------------------------------------------------------------
download_genome "GCF_000011985.1"  "L_acidophilus_NCFM"            \
    "Lactobacillus acidophilus NCFM"

echo ""
echo "============================================================"
printf " Done.  %d FASTA file(s) in %s\n" \
    "$(find "$OUTDIR" -name '*.fna' | wc -l)" "$OUTDIR"
echo "============================================================"
ls -lh "${OUTDIR}"/*.fna 2>/dev/null || echo "(no files found)"
