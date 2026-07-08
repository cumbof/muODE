# =============================================================================
# Final phase -- assemble the muODE hand-off bundle
# =============================================================================
# Gathers everything muODE ingests into one directory:
#   muode_inputs/mags/             non-redundant catalogue (-> muODE mags_dir)
#   muode_inputs/domains.tsv       (mag_id, domain)        (-> muODE traits:)
#   muode_inputs/abundance/*.tsv   (mag_id, rel_abundance) (-> muODE abundance / adapters)
#   muode_inputs/catalogue_quality.tsv   completeness/contamination/taxonomy
#   muode_inputs/MANIFEST.md       how to run muODE on this bundle

def _gtdbtk_input(_):
    return [f"{OUT}/taxonomy/gtdbtk.summary.tsv"] if DO_GTDBTK else []


rule muode_domains:
    input:
        catalogue=f"{OUT}/catalogue/mags",
        done=f"{OUT}/catalogue/dereplicate.done",
        merged_quality=f"{OUT}/catalogue/all_mags.quality.tsv",
        gtdbtk=_gtdbtk_input,
    output:
        mags=directory(f"{OUT}/muode_inputs/mags"),
        domains=f"{OUT}/muode_inputs/domains.tsv",
        quality=f"{OUT}/muode_inputs/catalogue_quality.tsv",
    conda:
        "../envs/base.yaml"
    params:
        gtdbtk_arg=lambda w, input: f"--gtdbtk {input.gtdbtk[0]}" if DO_GTDBTK else "",
    shell:
        r"""
        set -euo pipefail
        mkdir -p {output.mags}
        if compgen -G "{input.catalogue}/*.fa" > /dev/null; then
          cp {input.catalogue}/*.fa {output.mags}/
        fi
        python scripts/write_domains.py \
          --catalogue {output.mags} \
          --quality {input.merged_quality} {params.gtdbtk_arg} \
          --out-domains {output.domains} \
          --out-quality {output.quality}
        """


rule muode_manifest:
    input:
        domains=f"{OUT}/muode_inputs/domains.tsv",
        quality=f"{OUT}/muode_inputs/catalogue_quality.tsv",
        abundance=[f"{OUT}/muode_inputs/abundance/{s}.tsv" for s in SAMPLE_IDS],
    output:
        manifest=f"{OUT}/muode_inputs/MANIFEST.md",
    run:
        n_samples = len(SAMPLE_IDS)
        with open(output.manifest, "w") as fh:
            fh.write(
                "# muODE input bundle\n\n"
                "Produced by the `metagenomics/` pipeline. Hand these to muODE:\n\n"
                "| file | muODE consumer |\n|---|---|\n"
                "| `mags/` | `muode build` / workflow `mags_dir` |\n"
                "| `domains.tsv` | workflow `traits:` (domain routing) |\n"
                f"| `abundance/*.tsv` | `muode simulate` abundance / rCDI adapter ({n_samples} samples) |\n"
                "| `catalogue_quality.tsv` | provenance (completeness/contamination/taxonomy) |\n\n"
                "## Reconstruct + simulate\n\n"
                "```bash\n"
                "# from the muODE repo root, muode env active\n"
                "snakemake --snakefile workflow/Snakefile --use-conda --cores 16 \\\n"
                "  --config mags_dir=<this>/mags mag_extension=fa \\\n"
                "           traits=<this>/domains.tsv \\\n"
                "           abundance=<this>/abundance/<sample>.tsv\n"
                "```\n\n"
                "Viral/phage contigs are under `../viruses/` (geNomad + CheckV); they are\n"
                "not GEMs -- wire them into muODE's `PhageInfection` ecology layer.\n"
            )
