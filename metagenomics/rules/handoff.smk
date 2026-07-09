# =============================================================================
# Final phase -- assemble ONE muODE hand-off bundle PER SAMPLE
# =============================================================================
# Every sample is a self-contained muODE input set. Nothing is shared between
# samples: a bundle's genomes were assembled from that sample's own reads, plus
# reference genomes chosen from that sample's own profile. That is what lets you
# take one rCDI bundle as the recipient community and one healthy bundle as the
# donor, and inject the donor into the recipient -- without either community
# containing genomes it never yielded.
#
#   muode_inputs/<sample>/mags/                  -> muODE `mags_dir`
#   muode_inputs/<sample>/domains.tsv            -> muODE `traits:`
#   muode_inputs/<sample>/abundance.tsv          -> muODE `abundance`
#   muode_inputs/<sample>/genome_sources.tsv        mag vs reference provenance
#   muode_inputs/<sample>/catalogue_quality.tsv     completeness/contamination/taxonomy
#   muode_inputs/<sample>/reconciliation.tsv        what was simulated, and what wasn't
#   muode_inputs/<sample>/MANIFEST.md               how to run muODE on it
#
# The genome set is built FIRST (below), because CoverM must quantify the same set
# muODE simulates -- see rules/abundance.smk for why there is no table merge.

def _genome_set_inputs(wildcards):
    inputs = {
        "catalogue": f"{OUT}/catalogue/{wildcards.sample}/mags",
        "done": f"{OUT}/catalogue/{wildcards.sample}/dereplicate.done",
        "status": f"{OUT}/reconcile/{wildcards.sample}/status.tsv",
    }
    if DO_GTDBTK:
        inputs["gtdbtk"] = f"{OUT}/taxonomy/{wildcards.sample}/gtdbtk.summary.tsv"
    if DO_FETCH:
        inputs["references"] = f"{OUT}/references/{wildcards.sample}/genomes"
        inputs["fetched"] = f"{OUT}/references/{wildcards.sample}/fetched.tsv"
    return inputs


rule build_genome_set:
    """The genomes muODE will simulate: this sample's MAGs + its fetched references."""
    input:
        unpack(_genome_set_inputs),
    output:
        mags=directory(f"{OUT}/muode_inputs/{{sample}}/mags"),
        sources=f"{OUT}/muode_inputs/{{sample}}/genome_sources.tsv",
        classification=f"{OUT}/taxonomy/{{sample}}/classification.combined.tsv",
    conda:
        "../envs/base.yaml"
    params:
        gtdbtk=lambda w, input: f"--gtdbtk {input.gtdbtk}" if DO_GTDBTK else "",
        refs=lambda w, input: (
            f"--references {input.references} --fetched {input.fetched}"
            if DO_FETCH else ""),
    shell:
        r"""
        python scripts/build_genome_set.py \
          --sample {wildcards.sample} \
          --catalogue {input.catalogue} {params.gtdbtk} {params.refs} \
          --out-mags {output.mags} \
          --out-sources {output.sources} \
          --out-classification {output.classification}
        """


rule muode_domains:
    """(mag_id, domain) for the final genome set -- MAGs and references alike."""
    input:
        mags=f"{OUT}/muode_inputs/{{sample}}/mags",
        classification=f"{OUT}/taxonomy/{{sample}}/classification.combined.tsv",
        merged_quality=f"{OUT}/catalogue/{{sample}}/all_mags.quality.tsv",
    output:
        domains=f"{OUT}/muode_inputs/{{sample}}/domains.tsv",
        quality=f"{OUT}/muode_inputs/{{sample}}/catalogue_quality.tsv",
    conda:
        "../envs/base.yaml"
    shell:
        r"""
        python scripts/write_domains.py \
          --catalogue {input.mags} \
          --quality {input.merged_quality} \
          --gtdbtk {input.classification} \
          --out-domains {output.domains} \
          --out-quality {output.quality}
        """


rule muode_reconciliation:
    """Per-sample accounting: which organisms the simulation is actually about."""
    input:
        status=f"{OUT}/reconcile/{{sample}}/status.tsv",
        stats=f"{OUT}/reconcile/{{sample}}/profile_stats.tsv",
        sources=f"{OUT}/muode_inputs/{{sample}}/genome_sources.tsv",
        abundance=f"{OUT}/muode_inputs/{{sample}}/abundance.tsv",
        coverm=f"{OUT}/abundance/{{sample}}.coverm.tsv",
        fetched=(lambda w: f"{OUT}/references/{w.sample}/fetched.tsv") if DO_FETCH else [],
    output:
        reconciliation=f"{OUT}/muode_inputs/{{sample}}/reconciliation.tsv",
        summary=f"{OUT}/muode_inputs/{{sample}}/reconciliation_summary.tsv",
    conda:
        "../envs/base.yaml"
    params:
        fetched=lambda w, input: f"--fetched {input.fetched}" if DO_FETCH else "",
    shell:
        r"""
        python scripts/write_reconciliation.py \
          --sample {wildcards.sample} \
          --status {input.status} --stats {input.stats} \
          --sources {input.sources} --abundance {input.abundance} \
          --coverm {input.coverm} {params.fetched} \
          --out-reconciliation {output.reconciliation} \
          --out-summary {output.summary}
        """


rule reconciliation_summary:
    """One row per sample. A statistics table only -- no genomes are pooled."""
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


rule muode_manifest:
    input:
        mags=f"{OUT}/muode_inputs/{{sample}}/mags",
        domains=f"{OUT}/muode_inputs/{{sample}}/domains.tsv",
        quality=f"{OUT}/muode_inputs/{{sample}}/catalogue_quality.tsv",
        abundance=f"{OUT}/muode_inputs/{{sample}}/abundance.tsv",
        sources=f"{OUT}/muode_inputs/{{sample}}/genome_sources.tsv",
        reconciliation=f"{OUT}/muode_inputs/{{sample}}/reconciliation.tsv",
    output:
        manifest=f"{OUT}/muode_inputs/{{sample}}/MANIFEST.md",
    run:
        import csv as _csv

        sample = wildcards.sample
        srcs = list(_csv.DictReader(open(input.sources), delimiter="\t"))
        n_mag = sum(1 for r in srcs if r["source"] == "mag")
        n_ref = sum(1 for r in srcs if r["source"] == "reference")

        counts, not_simulated = {}, 0
        for row in _csv.DictReader(open(input.reconciliation), delimiter="\t"):
            counts[row["status"]] = counts.get(row["status"], 0) + 1
            if row["simulated"] == "no":
                not_simulated += 1

        with open(output.manifest, "w") as fh:
            fh.write(
                f"# muODE input bundle -- sample `{sample}`\n\n"
                "Produced by the `metagenomics/` pipeline. Every genome here comes from\n"
                f"`{sample}` alone -- assembled from its reads, or a reference genome for a\n"
                "species its own profile detected. Nothing is shared with other samples.\n\n"
                "| file | muODE consumer |\n|---|---|\n"
                f"| `mags/` | workflow `mags_dir` ({n_mag} MAG + {n_ref} reference = {len(srcs)}) |\n"
                "| `domains.tsv` | workflow `traits:` (domain routing) |\n"
                "| `abundance.tsv` | workflow `abundance` -- `(mag_id, rel_abundance)` |\n"
                "| `genome_sources.tsv` | provenance: `mag` vs `reference` (+ accession) |\n"
                "| `catalogue_quality.tsv` | completeness/contamination/taxonomy |\n"
                "| `reconciliation.tsv` | what was simulated, and what wasn't |\n\n"
                "## Reconstructed vs profiled\n\n"
                f"- `matched` ({counts.get('matched', 0)}): assembled **and** profiled.\n"
                f"- `reconstructed_only` ({counts.get('reconstructed_only', 0)}): assembled, "
                "not profiled (novel / below markers). Still simulated.\n"
                f"- `profiled_only` ({counts.get('profiled_only', 0)}): profiled, never "
                "assembled. Simulated **iff** a reference genome was fetched for it.\n\n"
                f"**{not_simulated} entr(ies) are not in the simulation.** See the\n"
                "`simulated` column of `reconciliation.tsv`, and\n"
                "`../reconciliation_summary.tsv` for the fraction of the profiled\n"
                "community these genomes account for.\n\n"
                f"> {n_ref} genome(s) here are **reference proxies** for species that never\n"
                "> assembled -- a public isolate of the same species, not the strain in this\n"
                "> sample. Their accessory gene content, and hence their GEM, may differ.\n\n"
                "## Reconstruct + simulate\n\n"
                "```bash\n"
                "# from the muODE repo root, muode env active\n"
                "snakemake --snakefile workflow/Snakefile --use-conda --cores 16 \\\n"
                "  --config mags_dir=<this>/mags mag_extension=fa \\\n"
                "           traits=<this>/domains.tsv \\\n"
                "           abundance=<this>/abundance.tsv\n"
                "```\n\n"
                "## Two-community experiments (e.g. FMT)\n\n"
                "Build one bundle per sample, then treat one as recipient and one as donor:\n\n"
                "```bash\n"
                "muode simulate --community <rCDI_sample>/community.json \\\n"
                "               --inject <healthy_sample>/community.json \\\n"
                "               --inject-time 24 --time 96 --outdir results_fmt\n"
                "```\n\n"
                "Viral/phage contigs are under `../../viruses/<sample>/` (geNomad + CheckV);\n"
                "they are not GEMs -- wire them into muODE's `PhageInfection` ecology layer.\n"
            )
