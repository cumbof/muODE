# =============================================================================
# Phase 8 -- abundance: CoverM coverage + reconciliation with the profile
# =============================================================================
# Two abundance estimates per sample, over that sample's OWN catalogue:
#
#   CoverM     -- maps the sample's reads back to its own MAGs. Direct and exact
#                 for genomes that assembled, but blind to everything that didn't,
#                 and it is a within-sample coverage share (not comparable across
#                 samples with different assembly success).
#   MetaPhlAn  -- marker-based SGB relative abundance. Sees taxa that never
#                 assembled, and is comparable across samples and studies.
#
# `reconcile_abundance.py` joins them on GTDB species and writes muODE's abundance
# over the INTERSECTION (using the profiler's value), plus a reconciliation table
# naming everything excluded from either side. See scripts/reconcile_abundance.py.

rule coverm_abundance:
    input:
        catalogue=f"{OUT}/catalogue/{{sample}}/mags",
        done=f"{OUT}/catalogue/{{sample}}/dereplicate.done",
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


def _reconcile_inputs(wildcards):
    """Only depend on the profile / taxonomy that the enabled tracks produce."""
    inputs = {"coverm": f"{OUT}/abundance/{wildcards.sample}.coverm.tsv"}
    if DO_GTDBTK:
        inputs["gtdbtk"] = f"{OUT}/taxonomy/{wildcards.sample}/gtdbtk.summary.tsv"
    if DO_PROFILE:
        inputs["profile"] = f"{OUT}/profile/{wildcards.sample}.metaphlan.tsv"
        inputs["sgb2gtdb"] = f"{OUT}/profile/sgb2gtdb.tsv"
    return inputs


rule reconcile_abundance:
    """Intersect this sample's reconstructed genomes with its profiled taxa."""
    input:
        unpack(_reconcile_inputs),
    output:
        abundance=f"{OUT}/muode_inputs/{{sample}}/abundance.tsv",
        reconciliation=f"{OUT}/muode_inputs/{{sample}}/reconciliation.tsv",
        summary=f"{OUT}/muode_inputs/{{sample}}/reconciliation_summary.tsv",
    conda:
        "../envs/base.yaml"
    params:
        gtdbtk=lambda w, input: f"--gtdbtk {input.gtdbtk}" if DO_GTDBTK else "",
        profile=lambda w, input: (
            f"--profile {input.profile} --sgb2gtdb {input.sgb2gtdb}" if DO_PROFILE else ""),
        source=config.get("abundance_source", "reconciled"),
    shell:
        r"""
        python scripts/reconcile_abundance.py \
          --sample {wildcards.sample} \
          --coverm {input.coverm} {params.gtdbtk} {params.profile} \
          --source {params.source} \
          --out-abundance {output.abundance} \
          --out-reconciliation {output.reconciliation} \
          --out-summary {output.summary}
        """


rule reconciliation_summary:
    """One row per sample: how much of each community the simulated genomes cover.

    This is a statistics table only -- no genomes are pooled across samples.
    """
    input:
        rows=[f"{OUT}/muode_inputs/{s}/reconciliation_summary.tsv" for s in SAMPLE_IDS],
    output:
        summary=f"{OUT}/muode_inputs/reconciliation_summary.tsv",
    conda:
        "../envs/base.yaml"
    shell:
        r"""
        set -euo pipefail
        head -n1 {input.rows[0]} > {output.summary}
        for f in {input.rows}; do tail -n +2 "$f" >> {output.summary}; done
        """
