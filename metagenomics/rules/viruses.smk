# =============================================================================
# Phase 7 -- virus / phage detection: geNomad (end-to-end) -> CheckV
# =============================================================================
# geNomad identifies viral AND plasmid contigs from an assembly and assigns viral
# taxonomy in one pass; CheckV then estimates completeness/quality and trims host
# flanks off proviruses. These are NOT reconstructed into GEMs (viruses have no
# metabolism) -- they are the substrate for muODE's PhageInfection ecology layer,
# so we surface them as a per-sample virus catalogue + summary.

rule genomad:
    """Identify viral + plasmid contigs and assign viral taxonomy."""
    input:
        contigs=f"{OUT}/assembly/{{sample}}/contigs.fa",
    output:
        virus=f"{OUT}/viruses/{{sample}}/genomad/{{sample}}_virus.fna",
        summary=f"{OUT}/viruses/{{sample}}/genomad/{{sample}}_virus_summary.tsv",
    threads: config["threads"]
    conda:
        "../envs/genomad.yaml"
    params:
        outdir=lambda w: f"{OUT}/viruses/{w.sample}/genomad",
        db=config.get("genomad_db", ""),
    shell:
        r"""
        set -euo pipefail
        if [ -z "{params.db}" ]; then
          echo "genomad: detect_viruses is true but genomad_db is empty" >&2
          exit 1
        fi
        genomad end-to-end --cleanup --threads {threads} \
          {input.contigs} {params.outdir} {params.db}
        # geNomad nests outputs under <input-stem>_summary/ ; normalise the paths.
        stem=$(basename {input.contigs} .fa)
        cp {params.outdir}/${{stem}}_summary/${{stem}}_virus.fna {output.virus}
        cp {params.outdir}/${{stem}}_summary/${{stem}}_virus_summary.tsv {output.summary}
        """


rule checkv:
    """Quality/completeness of the geNomad viral contigs; trims provirus flanks."""
    input:
        virus=f"{OUT}/viruses/{{sample}}/genomad/{{sample}}_virus.fna",
    output:
        quality=f"{OUT}/viruses/{{sample}}/checkv/quality_summary.tsv",
    threads: config["threads"]
    conda:
        "../envs/checkv.yaml"
    params:
        outdir=lambda w: f"{OUT}/viruses/{w.sample}/checkv",
        db=config.get("checkv_db", ""),
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir}
        if [ ! -s {input.virus} ]; then
          printf 'contig_id\tcontig_length\tcompleteness\tcheckv_quality\n' \
            > {output.quality}
          exit 0
        fi
        dbarg=""
        [ -n "{params.db}" ] && dbarg="-d {params.db}"
        # CheckV writes quality_summary.tsv straight into {params.outdir}, which is
        # exactly {output.quality} -- no copy needed (a cp here would be file-onto-itself).
        checkv end_to_end {input.virus} {params.outdir} -t {threads} $dbarg
        """


rule virus_summary:
    """Per-sample marker gathering the virus catalogue + CheckV quality."""
    input:
        virus=f"{OUT}/viruses/{{sample}}/genomad/{{sample}}_virus.fna",
        checkv=f"{OUT}/viruses/{{sample}}/checkv/quality_summary.tsv",
    output:
        done=f"{OUT}/viruses/{{sample}}/{{sample}}.virus_summary.done",
    shell:
        "touch {output.done}"
