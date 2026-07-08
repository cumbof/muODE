# =============================================================================
# Phase 3 -- prokaryotic MAG quality gate (CheckM2)
# =============================================================================
# CheckM2 estimates completeness/contamination for every DAS_Tool bin; passing
# bins are copied under stable ids  <sample>__bin.<k>  ready for cross-sample
# dereplication. (CheckM2 scores bacteria AND archaea, so archaea need no
# separate track here; GTDB-Tk classifies both domains downstream.)

rule checkm2:
    input:
        binsdir=f"{OUT}/bins/{{sample}}/dastool_bins",
        done=f"{OUT}/bins/{{sample}}/dastool.done",
    output:
        report=f"{OUT}/checkm2/{{sample}}/quality_report.tsv",
    threads: config["threads"]
    conda:
        "../envs/checkm2.yaml"
    params:
        outdir=lambda w: f"{OUT}/checkm2/{w.sample}",
        db=config.get("checkm2_db", ""),
        local_tmp=config.get("local_tmp", ""),
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir}
        if ! compgen -G "{input.binsdir}/*.fa" > /dev/null; then
          # no bins recovered for this sample -> emit an empty report header
          printf 'Name\tCompleteness\tContamination\n' > {output.report}
          exit 0
        fi
        # CheckM2's prediction phase uses a multiprocessing Manager whose scratch
        # dir lives in $TMPDIR. On an NFS $TMPDIR (e.g. an isilon mount) the manager
        # shutdown fails with EBUSY ("Device or resource busy") because NFS cannot
        # unlink a still-open file. Pin temp to local disk (config local_tmp, else
        # $SLURM_TMPDIR, else /tmp) and clean it up on exit.
        tmpbase="{params.local_tmp}"
        [ -n "$tmpbase" ] || tmpbase="${{SLURM_TMPDIR:-/tmp}}"
        tmpd="$(mktemp -d "$tmpbase/checkm2.{wildcards.sample}.XXXXXX")"
        export TMPDIR="$tmpd"
        trap 'rm -rf "$tmpd" 2>/dev/null || true' EXIT
        dbarg=""
        [ -n "{params.db}" ] && dbarg="--database_path {params.db}"
        checkm2 predict --threads {threads} -x fa \
          -i {input.binsdir} -o {params.outdir} --force $dbarg
        """


rule select_prok_mags:
    """Copy passing bins to a per-sample dir under stable <sample>__bin.<k> ids."""
    input:
        report=f"{OUT}/checkm2/{{sample}}/quality_report.tsv",
        binsdir=f"{OUT}/bins/{{sample}}/dastool_bins",
    output:
        magdir=directory(f"{OUT}/mags_per_sample/{{sample}}/prok"),
        quality=f"{OUT}/mags_per_sample/{{sample}}/prok.quality.tsv",
    conda:
        "../envs/base.yaml"
    params:
        mincomp=config["min_completeness"],
        maxcont=config["max_contamination"],
    shell:
        r"""
        python scripts/select_mags.py \
          --report {input.report} --bins {input.binsdir} \
          --sample {wildcards.sample} --tag bin \
          --outdir {output.magdir} --quality-out {output.quality} \
          --min-completeness {params.mincomp} --max-contamination {params.maxcont}
        """
