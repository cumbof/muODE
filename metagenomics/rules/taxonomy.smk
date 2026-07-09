# =============================================================================
# Phase 5 -- prokaryotic taxonomy (GTDB-Tk classifies BACTERIA + ARCHAEA)
# =============================================================================
# Runs per sample, against that sample's own catalogue. One tool covers both
# domains, so the bacteria/archaea split in domains.tsv is read straight from the
# `d__Bacteria` / `d__Archaea` prefix of the classification. Eukaryotic MAGs (ids
# containing `__euk.`) are excluded here.
#
# The GTDB species assigned here is also the KEY that joins a MAG to the
# quantitative profile (rules/profile.smk), so this must run for the
# reconstructed-vs-profiled reconciliation to be possible at all.

rule gtdbtk:
    input:
        catalogue=f"{OUT}/catalogue/{{sample}}/mags",
        done=f"{OUT}/catalogue/{{sample}}/dereplicate.done",
    output:
        summary=f"{OUT}/taxonomy/{{sample}}/gtdbtk.summary.tsv",
    threads: config["threads"]
    conda:
        "../envs/gtdbtk.yaml"
    params:
        db=config.get("gtdbtk_db", ""),
        genomes=lambda w: f"{OUT}/taxonomy/{w.sample}/prok_genomes",
        outdir=lambda w: f"{OUT}/taxonomy/{w.sample}/gtdbtk",
        local_tmp=config.get("local_tmp", ""),
    shell:
        r"""
        set -euo pipefail
        if [ -z "{params.db}" ]; then
          echo "gtdbtk: classify_gtdbtk is true but gtdbtk_db is empty" >&2
          exit 1
        fi
        export GTDBTK_DATA_PATH="{params.db}"

        # GTDB-Tk parallelises with Python multiprocessing, whose Manager scratch
        # dir lives in $TMPDIR. On an NFS $TMPDIR (e.g. an isilon mount) the
        # manager shutdown fails with EBUSY because NFS cannot unlink a file that
        # is still open -- the same failure CheckM2 hits (see rules/checkm.smk).
        # Pin temp to local disk: config local_tmp, else $SLURM_TMPDIR, else /tmp.
        tmpbase="{params.local_tmp}"
        [ -n "$tmpbase" ] || tmpbase="${{SLURM_TMPDIR:-/tmp}}"
        tmpd="$(mktemp -d "$tmpbase/gtdbtk.{wildcards.sample}.XXXXXX")"
        export TMPDIR="$tmpd"
        trap 'rm -rf "$tmpd" 2>/dev/null || true' EXIT

        rm -rf {params.genomes}; mkdir -p {params.genomes} {params.outdir}

        # Prokaryotic catalogue members only (exclude euk MAGs).
        shopt -s nullglob
        n=0
        for f in {input.catalogue}/*__bin.*.fa; do
          ln -sf "$(readlink -f "$f")" {params.genomes}/; n=$((n+1))
        done
        if [ "$n" -eq 0 ]; then
          printf 'user_genome\tclassification\n' > {output.summary}
          exit 0
        fi

        # GTDB-Tk >=2.7 does the ANI pre-screen with skani and no longer accepts
        # `--skip_ani_screen` / `--mash_db` (mash was dropped); passing them is an
        # argparse error. The env pins >=2.7, which is also the first release that
        # supports the R232 data.
        gtdbtk classify_wf --genome_dir {params.genomes} --out_dir {params.outdir} \
          -x fa --cpus {threads}

        # Merge the bacterial + archaeal summaries into one (keep header once).
        head -n1 {params.outdir}/gtdbtk.bac120.summary.tsv 2>/dev/null \
          > {output.summary} || printf 'user_genome\tclassification\n' > {output.summary}
        for s in {params.outdir}/gtdbtk.bac120.summary.tsv \
                 {params.outdir}/gtdbtk.ar53.summary.tsv; do
          [ -e "$s" ] && tail -n +2 "$s" >> {output.summary}
        done
        """
