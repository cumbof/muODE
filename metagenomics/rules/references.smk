# =============================================================================
# Phase 7c -- reconcile the two views, and fill the gap from NCBI
# =============================================================================
# `reconcile_profile` decides which genomes exist for this sample: assembled,
# profiled, or both. It computes no abundance -- quantification happens once,
# later, when CoverM maps the reads against the final genome set.
#
# `fetch_reference_genomes` then closes the most important gap: species the
# profiler saw but assembly never recovered. Without a genome they simply cannot
# be simulated, so a public reference genome of the same GTDB species stands in.
# This EXPANDS the simulated community rather than shrinking it to the
# intersection -- and, because those reads now have somewhere correct to map, it
# also removes the inflation they were causing on assembled relatives.
#
# Everything stays per-sample: the reference genomes fetched for one sample are
# chosen from that sample's own profile and never shared with another.

def _reconcile_inputs(wildcards):
    inputs = {"catalogue": f"{OUT}/catalogue/{wildcards.sample}/mags",
              "done": f"{OUT}/catalogue/{wildcards.sample}/dereplicate.done"}
    if DO_GTDBTK:
        inputs["gtdbtk"] = f"{OUT}/taxonomy/{wildcards.sample}/gtdbtk.summary.tsv"
    if DO_PROFILE:
        inputs["profile"] = f"{OUT}/profile/{wildcards.sample}.metaphlan.tsv"
        inputs["sgb2gtdb"] = f"{OUT}/profile/sgb2gtdb.tsv"
    return inputs


rule reconcile_profile:
    """Classify this sample's genomes: matched / reconstructed_only / profiled_only."""
    input:
        unpack(_reconcile_inputs),
    output:
        status=f"{OUT}/reconcile/{{sample}}/status.tsv",
        missing=f"{OUT}/reconcile/{{sample}}/missing_species.tsv",
        stats=f"{OUT}/reconcile/{{sample}}/profile_stats.tsv",
    conda:
        "../envs/base.yaml"
    params:
        gtdbtk=lambda w, input: f"--gtdbtk {input.gtdbtk}" if DO_GTDBTK else "",
        profile=lambda w, input: (
            f"--profile {input.profile} --sgb2gtdb {input.sgb2gtdb}"
            if DO_PROFILE else ""),
    shell:
        r"""
        python scripts/reconcile_profile.py \
          --sample {wildcards.sample} \
          --catalogue {input.catalogue} {params.gtdbtk} {params.profile} \
          --out-status {output.status} \
          --out-missing {output.missing} \
          --out-stats {output.stats}
        """


rule fetch_reference_genomes:
    """Download an NCBI reference genome for each profiled-but-unassembled species.

    Accessions come from `gtdb_taxonomy.tsv` inside the GTDB-Tk reference data, so
    they live in the same taxonomy the MAGs were classified with. The fetched
    genome is a PROXY for the species, not the sample's strain -- recorded as
    `source=reference` in every downstream table.
    """
    input:
        missing=f"{OUT}/reconcile/{{sample}}/missing_species.tsv",
    output:
        genomes=directory(f"{OUT}/references/{{sample}}/genomes"),
        table=f"{OUT}/references/{{sample}}/fetched.tsv",
        classification=f"{OUT}/references/{{sample}}/classification.tsv",
    conda:
        "../envs/datasets.yaml"
    params:
        taxonomy=lambda w: (config.get("gtdb_taxonomy")
                            or f"{config.get('gtdbtk_db', '').rstrip('/')}/taxonomy/gtdb_taxonomy.tsv"),
        min_ab=config.get("reference_min_abundance", 0.1),
        max_n=config.get("max_reference_genomes", 50),
    shell:
        r"""
        set -euo pipefail
        if [ ! -f "{params.taxonomy}" ]; then
          echo "fetch_reference_genomes: no GTDB taxonomy at {params.taxonomy}" >&2
          echo "Set config 'gtdb_taxonomy' to the gtdb_taxonomy.tsv of your GTDB" >&2
          echo "release, or set 'fetch_missing_genomes: false'." >&2
          exit 1
        fi
        python scripts/fetch_reference_genomes.py \
          --sample {wildcards.sample} \
          --missing {input.missing} \
          --gtdb-taxonomy "{params.taxonomy}" \
          --outdir {output.genomes} \
          --out-table {output.table} \
          --out-classification {output.classification} \
          --min-abundance {params.min_ab} \
          --max-genomes {params.max_n}
        """
