#!/usr/bin/env bash
# Genomes for the genome-scale phage-therapy scenario: a pathobiont host and a commensal
# competitor.  Both facultative gram-negatives -> CarveMe (-u gramneg), one BiGG pool.
#
#   host      Klebsiella pneumoniae HS11286  -- a real phage-therapy target pathobiont
#   commensal Escherichia coli K-12 MG1655   -- benign gut resident, competes for carbon,
#                                               NOT a host for the Klebsiella phage
#
# The phage itself has no genome-scale model: it is a PhageInfection ecology layer (see
# genome_scale.py), so nothing is downloaded for it.
#
# Run from the repository root:
#   bash examples/phage_therapy/download_genomes.sh
set -euo pipefail

OUTDIR="$(dirname "$(realpath "$0")")/data/mags"
mkdir -p "${OUTDIR}"

download_genome() {
    local accession="$1" name="$2" note="$3"
    local outfile="${OUTDIR}/${name}.fna"
    if [[ -f "${outfile}" ]]; then printf "  [skip] %-16s (exists)\n" "${name}"; return; fi
    printf "  [fetch] %-16s %s  (%s)\n" "${name}" "${accession}" "${note}"
    local tmp; tmp=$(mktemp -d)
    datasets download genome accession "${accession}" --include genome \
        --filename "${tmp}/g.zip" --no-progressbar 2>/dev/null
    unzip -q -o "${tmp}/g.zip" -d "${tmp}/x"
    cp "$(find "${tmp}/x" -name '*.fna' | head -1)" "${outfile}"
    rm -rf "${tmp}"
    printf "  [ok]    %-16s -> %s\n" "${name}" "${outfile}"
}

echo "============================================================"
echo " Phage-therapy example genomes (host + commensal)"
echo "============================================================"
download_genome "GCF_000240185.1" "K_pneumoniae" "pathobiont host of the phage"
download_genome "GCF_000005845.2" "E_coli"       "commensal competitor (phage-resistant)"

cat <<'NOTE'

------------------------------------------------------------
  * K_pneumoniae, E_coli  -> CarveMe (-u gramneg), from these genomes.
  * the phage             -> no GEM; a PhageInfection layer (genome_scale.py).
Then reconstruct + run:
  snakemake --use-conda --cores 8 --configfile examples/phage_therapy/config.yaml
  python examples/phage_therapy/genome_scale.py --models <refined GEMs dir>
------------------------------------------------------------
NOTE
