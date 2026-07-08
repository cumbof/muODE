#!/usr/bin/env bash
# Download three STRAINS OF ONE SPECIES (Escherichia coli) for the strain-
# competition example. The point is strain-level variation within a species, so
# these must be strain-resolved genomes -- reference isolates are the clean way
# to get that (a dereplicated metagenome catalogue would have merged them).
#
# The strain->role mapping below is ILLUSTRATIVE: it picks real strains whose
# biology motivates the toy roles in mechanistic_demo.py. With real GEMs, which
# strain actually wins is an OUTPUT of reconstruction, not assigned here. Swap in
# any strains you like -- keep the output filenames matching abundance.tsv /
# traits.tsv (Ecoli__glc_specialist, Ecoli__colicinogenic, Ecoli__arabinose_user).
#
# Accessions are representative RefSeq assemblies; verify/substitute on NCBI if
# one has been superseded.
#
# Requirements:
#   conda install -c conda-forge ncbi-datasets-cli
#
# Run from the repository root:
#   bash examples/strain_competition/download_genomes.sh

set -euo pipefail

OUTDIR="$(dirname "$(realpath "$0")")/data/mags"
mkdir -p "$OUTDIR"

download_genome() {
    local accession="$1"
    local name="$2"
    local note="$3"
    local outfile="${OUTDIR}/${name}.fna"

    if [[ -f "$outfile" ]]; then
        printf "  [skip]  %-28s  (already exists)\n" "$name"
        return
    fi
    printf "  [fetch] %-28s  %s  (%s)\n" "$name" "$accession" "$note"

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
    printf "  [ok]    %-28s  -> %s\n" "$name" "$outfile"
}

echo "============================================================"
echo " Strain-competition example genomes (3 strains of E. coli)"
echo " Output: ${OUTDIR}"
echo "============================================================"

# glucose specialist  <- lab generalist, workhorse central metabolism
download_genome "GCF_000005845.2"  "Ecoli__glc_specialist"  \
    "E. coli K-12 MG1655 -- generalist / glucose"
# colicin producer    <- probiotic strain known for bacteriocin/microcin output
download_genome "GCF_000714595.1"  "Ecoli__colicinogenic"   \
    "E. coli Nissle 1917 -- bacteriocin producer"
# private-niche user  <- a divergent strain with distinct accessory catabolism
download_genome "GCF_000008865.2"  "Ecoli__arabinose_user"  \
    "E. coli O157:H7 Sakai -- accessory sugar catabolism"

cat <<'NOTE'

------------------------------------------------------------
All three are gram-negative bacteria -> CarveMe (-u gramneg), one per strain,
sharing the BiGG namespace so they compete in one extracellular pool. See
config.yaml. Reconstruct + simulate with:

  snakemake --snakefile workflow/Snakefile --use-conda --cores 8 \
            --configfile examples/strain_competition/config.yaml

Reminder: keep the strains strain-resolved. Do NOT source them from the
dereplicated ../../metagenomics/ catalogue -- dRep (0.95 ANI) merges same-species
strains into one representative, which erases exactly the variation this example
is about.
------------------------------------------------------------
NOTE
