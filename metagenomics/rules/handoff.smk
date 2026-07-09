# =============================================================================
# Final phase -- assemble ONE muODE hand-off bundle PER SAMPLE
# =============================================================================
# Every sample is a self-contained muODE input set. Nothing is shared between
# samples: a bundle contains only genomes assembled from that sample's own reads.
# That is what lets you take, say, one rCDI bundle as the recipient community and
# one healthy bundle as the donor, and inject the donor into the recipient --
# without either community containing genomes it never yielded.
#
#   muode_inputs/<sample>/mags/                  -> muODE `mags_dir`
#   muode_inputs/<sample>/domains.tsv            -> muODE `traits:`
#   muode_inputs/<sample>/abundance.tsv          -> muODE `abundance`
#   muode_inputs/<sample>/catalogue_quality.tsv     provenance for the above
#   muode_inputs/<sample>/reconciliation.tsv        reconstructed vs profiled
#   muode_inputs/<sample>/MANIFEST.md               how to run muODE on it
#
# `mags/` holds exactly the genomes named in `abundance.tsv` -- i.e. the simulated
# community. Genomes assembled but excluded (and taxa profiled but never
# assembled) are not lost: they are itemised in `reconciliation.tsv`.

def _gtdbtk_input(wildcards):
    return ([f"{OUT}/taxonomy/{wildcards.sample}/gtdbtk.summary.tsv"]
            if DO_GTDBTK else [])


rule muode_bundle:
    input:
        catalogue=f"{OUT}/catalogue/{{sample}}/mags",
        done=f"{OUT}/catalogue/{{sample}}/dereplicate.done",
        merged_quality=f"{OUT}/catalogue/{{sample}}/all_mags.quality.tsv",
        abundance=f"{OUT}/muode_inputs/{{sample}}/abundance.tsv",
        gtdbtk=_gtdbtk_input,
    output:
        mags=directory(f"{OUT}/muode_inputs/{{sample}}/mags"),
        domains=f"{OUT}/muode_inputs/{{sample}}/domains.tsv",
        quality=f"{OUT}/muode_inputs/{{sample}}/catalogue_quality.tsv",
    conda:
        "../envs/base.yaml"
    params:
        gtdbtk_arg=lambda w, input: f"--gtdbtk {input.gtdbtk[0]}" if DO_GTDBTK else "",
    shell:
        r"""
        set -euo pipefail
        rm -rf {output.mags}; mkdir -p {output.mags}

        # Copy exactly the genomes muODE will simulate: those named in the
        # reconciled abundance table (its first column, minus the header).
        list="$(dirname {output.mags})/.mag_list.txt"
        tail -n +2 {input.abundance} | cut -f1 > "$list"
        n=0
        while read -r mag; do
          [ -n "$mag" ] || continue
          src="{input.catalogue}/${{mag}}.fa"
          if [ ! -f "$src" ]; then
            echo "muode_bundle: {wildcards.sample}: no FASTA for '$mag'" >&2; exit 1
          fi
          cp "$src" {output.mags}/; n=$((n+1))
        done < "$list"
        rm -f "$list"
        echo "muode_bundle: {wildcards.sample}: $n genome(s)" >&2

        python scripts/write_domains.py \
          --catalogue {output.mags} \
          --quality {input.merged_quality} {params.gtdbtk_arg} \
          --out-domains {output.domains} \
          --out-quality {output.quality}
        """


rule muode_manifest:
    input:
        mags=f"{OUT}/muode_inputs/{{sample}}/mags",
        domains=f"{OUT}/muode_inputs/{{sample}}/domains.tsv",
        quality=f"{OUT}/muode_inputs/{{sample}}/catalogue_quality.tsv",
        abundance=f"{OUT}/muode_inputs/{{sample}}/abundance.tsv",
        reconciliation=f"{OUT}/muode_inputs/{{sample}}/reconciliation.tsv",
    output:
        manifest=f"{OUT}/muode_inputs/{{sample}}/MANIFEST.md",
    run:
        import csv as _csv
        from pathlib import Path as _Path

        sample = wildcards.sample
        n_mags = len(list(_Path(input.mags).glob("*.fa")))
        counts = {}
        with open(input.reconciliation) as fh:
            for row in _csv.DictReader(fh, delimiter="\t"):
                counts[row["status"]] = counts.get(row["status"], 0) + 1

        with open(output.manifest, "w") as fh:
            fh.write(
                f"# muODE input bundle -- sample `{sample}`\n\n"
                "Produced by the `metagenomics/` pipeline. Every genome here was\n"
                f"assembled from `{sample}`'s own reads; nothing is shared with other\n"
                "samples, so this bundle is a self-contained community.\n\n"
                "| file | muODE consumer |\n|---|---|\n"
                f"| `mags/` | workflow `mags_dir` ({n_mags} genome(s)) |\n"
                "| `domains.tsv` | workflow `traits:` (domain routing) |\n"
                "| `abundance.tsv` | workflow `abundance` -- `(mag_id, rel_abundance)` |\n"
                "| `catalogue_quality.tsv` | provenance (completeness/contamination/taxonomy) |\n"
                "| `reconciliation.tsv` | reconstructed vs profiled (see below) |\n\n"
                "## Reconstructed vs profiled\n\n"
                f"- `matched` ({counts.get('matched', 0)}): assembled **and** profiled -- simulated.\n"
                f"- `reconstructed_only` ({counts.get('reconstructed_only', 0)}): assembled, not "
                "profiled (novel / below markers). **Not simulated.**\n"
                f"- `profiled_only` ({counts.get('profiled_only', 0)}): profiled, never assembled. "
                "**Not simulated** -- substitute a reference genome if you need it.\n\n"
                "`../reconciliation_summary.tsv` reports what fraction of the profiled\n"
                "community the simulated genomes account for.\n\n"
                "## Reconstruct + simulate\n\n"
                "```bash\n"
                "# from the muODE repo root, muode env active\n"
                "snakemake --snakefile workflow/Snakefile --use-conda --cores 16 \\\n"
                f"  --config mags_dir=<this>/mags mag_extension=fa \\\n"
                f"           traits=<this>/domains.tsv \\\n"
                f"           abundance=<this>/abundance.tsv\n"
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
