# =============================================================================
# Phase 8 -- per-sample abundance profiling (CoverM) -> muODE abundance TSVs
# =============================================================================
# Each sample's reads are mapped to the dereplicated catalogue and CoverM reports
# per-genome relative abundance. coverm_to_muode.py converts each into muODE's
# 2-column `(mag_id, rel_abundance)` contract -- the same file muode's community
# loader and the rCDI adapter read.

rule coverm_abundance:
    input:
        catalogue=f"{OUT}/catalogue/mags",
        done=f"{OUT}/catalogue/dereplicate.done",
        r1=lambda w: clean_reads(w)["r1"],
        r2=lambda w: clean_reads(w)["r2"],
    output:
        raw=f"{OUT}/abundance/{{sample}}.coverm.tsv",
    threads: config["threads"]
    conda:
        "../envs/coverm.yaml"
    shell:
        r"""
        set -euo pipefail
        if ! compgen -G "{input.catalogue}/*.fa" > /dev/null; then
          printf 'Genome\t{wildcards.sample}\nunmapped\t100\n' > {output.raw}
          exit 0
        fi
        coverm genome \
          --coupled {input.r1} {input.r2} \
          --genome-fasta-directory {input.catalogue} -x fa \
          --methods relative_abundance \
          --threads {threads} \
          --output-file {output.raw}
        """


rule abundance_to_muode:
    input:
        raw=f"{OUT}/abundance/{{sample}}.coverm.tsv",
    output:
        tsv=f"{OUT}/muode_inputs/abundance/{{sample}}.tsv",
    conda:
        "../envs/base.yaml"
    shell:
        r"""
        python scripts/coverm_to_muode.py --coverm {input.raw} --out {output.tsv}
        """
