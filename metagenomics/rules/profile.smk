# =============================================================================
# Phase 7 -- quantitative taxonomic profiling (MetaPhlAn4, SGB-level)
# =============================================================================
# Assembly answers "which genomes can I reconstruct?"; a quantitative profiler
# answers "which taxa are actually here, and at what abundance?". They disagree:
# a species can be profiled from a handful of marker reads yet never assemble
# (assembly + binning realistically need >=5-10x even coverage), and a novel MAG
# can assemble yet carry no markers. muODE needs abundances for the genomes it
# simulates, so both views are computed and then reconciled (rules/abundance.smk).
#
# MetaPhlAn4 is used because its SGB (species-level genome bin) database covers
# UNCHARACTERISED taxa as uSGBs -- exactly the uncultivated organisms that MAG
# assembly recovers -- rather than only named reference species.
#
# `--unclassified_estimation` makes MetaPhlAn report the fraction of the sample it
# could not assign, so the hand-off can state how much of the community the
# simulated genomes actually represent instead of quietly renormalising it away.

rule metaphlan_sgb2gtdb_table:
    """Locate the SGB->GTDB mapping shipped with the MetaPhlAn database.

    MetaPhlAn reports SGB ids; GTDB-Tk assigns GTDB lineages to MAGs. This table
    (the one `sgb_to_gtdb_profile.py` itself reads) is what puts both in a single
    namespace so they can be joined. Without it there is no reconciliation, so
    this fails loudly rather than letting the pipeline fall back to CoverM.
    """
    output:
        table=f"{OUT}/profile/sgb2gtdb.tsv",
    params:
        db=config.get("metaphlan_db", ""),
        index=config.get("metaphlan_index", ""),
        explicit=config.get("metaphlan_sgb2gtdb", ""),
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.table})"

        src="{params.explicit}"
        if [ -z "$src" ]; then
          if [ -z "{params.db}" ]; then
            echo "metaphlan: run_profiling is true but metaphlan_db is empty" >&2
            exit 1
          fi
          # Prefer the table matching the pinned index, else any SGB2GTDB table.
          if [ -n "{params.index}" ] && [ -f "{params.db}/{params.index}_SGB2GTDB.tsv" ]; then
            src="{params.db}/{params.index}_SGB2GTDB.tsv"
          else
            src="$(ls -1 {params.db}/*_SGB2GTDB.tsv 2>/dev/null | sort | tail -1 || true)"
          fi
        fi

        if [ -z "$src" ] || [ ! -f "$src" ]; then
          echo "metaphlan: could not find an SGB->GTDB table (*_SGB2GTDB.tsv) in" >&2
          echo "  metaphlan_db='{params.db}'" >&2
          echo "Without it, MetaPhlAn SGBs cannot be matched to GTDB-Tk MAG species." >&2
          echo "Point config 'metaphlan_sgb2gtdb' at the table shipped with your DB," >&2
          echo "or set 'run_profiling: false' to fall back to CoverM abundance." >&2
          exit 1
        fi
        cp "$src" {output.table}
        echo "metaphlan: SGB->GTDB table = $src" >&2
        """


rule metaphlan:
    """SGB-level relative abundance for one sample, from its cleaned reads."""
    input:
        unpack(clean_reads),
    output:
        profile=f"{OUT}/profile/{{sample}}.metaphlan.tsv",
        bt2=f"{OUT}/profile/{{sample}}.bowtie2.bz2",
    threads: config["threads"]
    conda:
        "../envs/metaphlan.yaml"
    params:
        db=config.get("metaphlan_db", ""),
        index=config.get("metaphlan_index", ""),
    shell:
        r"""
        set -euo pipefail
        if [ -z "{params.db}" ]; then
          echo "metaphlan: run_profiling is true but metaphlan_db is empty" >&2
          exit 1
        fi
        # MetaPhlAn refuses to overwrite an existing --bowtie2out; a killed job
        # leaves one behind and every rerun would then fail.
        rm -f {output.bt2}

        idx=""
        [ -n "{params.index}" ] && idx="--index {params.index}"

        # Paired files are passed comma-separated: MetaPhlAn maps reads
        # independently against its marker set (it does not use pairing).
        metaphlan {input.r1},{input.r2} \
          --input_type fastq \
          --bowtie2db {params.db} $idx \
          --bowtie2out {output.bt2} \
          --nproc {threads} \
          --unclassified_estimation \
          --sample_id {wildcards.sample} \
          -o {output.profile}
        """
