# =============================================================================
# Phase 8 -- abundance: ONE CoverM pass over the final genome set
# =============================================================================
# muODE takes exactly one abundance table, so this pipeline produces exactly one
# measurement -- it never merges two.
#
# The temptation is to concatenate CoverM's numbers for the MAGs with MetaPhlAn's
# numbers for the fetched reference genomes. That would be wrong: CoverM's
# `relative_abundance` is a fraction of READS (with an `unmapped` remainder),
# MetaPhlAn's is a marker-normalised fraction of CELLS. Renormalising across both
# silently mixes two different quantities.
#
# So the sets are merged at the GENOME level (rules/handoff.smk builds
# MAGs + references), and CoverM is then run ONCE over that final set: one tool,
# one denominator, one table that sums to 1. MetaPhlAn's role was discovery --
# telling us which genomes were missing -- and reporting, not quantification.

rule coverm_abundance:
    input:
        genomes=f"{OUT}/muode_inputs/{{sample}}/mags",
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
        if ! compgen -G "{input.genomes}/*.fa" > /dev/null; then
          printf 'Genome\t{wildcards.sample}\nunmapped\t100\n' > {output.raw}
          exit 0
        fi
        coverm genome \
          --coupled {input.r1} {input.r2} \
          --genome-fasta-directory {input.genomes} -x fa \
          --methods relative_abundance \
          --threads {threads} \
          --output-file {output.raw}
        """


rule abundance_to_muode:
    """CoverM percentages -> muODE's `(mag_id, rel_abundance)` contract."""
    input:
        raw=f"{OUT}/abundance/{{sample}}.coverm.tsv",
    output:
        tsv=f"{OUT}/muode_inputs/{{sample}}/abundance.tsv",
    conda:
        "../envs/base.yaml"
    shell:
        r"""
        python scripts/coverm_to_muode.py --coverm {input.raw} --out {output.tsv}
        """
