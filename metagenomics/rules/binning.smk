# =============================================================================
# Phase 2 -- binning: MetaBAT2 + MaxBin2 + CONCOCT, consolidated by DAS_Tool
# =============================================================================
# Three complementary binners are run and DAS_Tool picks the best,
# non-redundant set of bins across all of them -- consistently better and more
# complete than any single binner. Enable/disable individual binners in config
# (`binners:`); DAS_Tool consolidates whichever ran (>=2 recommended).

ENABLED_BINNERS = [b for b, on in config.get("binners", {}).items() if on]
if len(ENABLED_BINNERS) < 1:
    raise ValueError("config: at least one entry under `binners:` must be true")


rule contig_depth:
    """Per-contig coverage table shared by MetaBAT2 and MaxBin2 (jgi format)."""
    input:
        bam=f"{OUT}/mapping/{{sample}}.sorted.bam",
        bai=f"{OUT}/mapping/{{sample}}.sorted.bam.bai",
    output:
        depth=f"{OUT}/bins/{{sample}}/depth.txt",
    conda:
        "../envs/metabat2.yaml"
    shell:
        "jgi_summarize_bam_contig_depths --outputDepth {output.depth} {input.bam}"


rule metabat2:
    input:
        contigs=f"{OUT}/assembly/{{sample}}/contigs.fa",
        depth=f"{OUT}/bins/{{sample}}/depth.txt",
    output:
        c2b=f"{OUT}/bins/{{sample}}/metabat2.contig2bin.tsv",
        binsdir=directory(f"{OUT}/bins/{{sample}}/metabat2"),
    threads: config["threads"]
    conda:
        "../envs/metabat2.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {output.binsdir}
        metabat2 -i {input.contigs} -a {input.depth} \
          -o {output.binsdir}/bin -m 1500 -t {threads} || true
        {{ for f in {output.binsdir}/*.fa; do
             [ -e "$f" ] || continue
             b=$(basename "$f" .fa)
             grep '^>' "$f" | sed 's/^>//' | awk -v b="$b" '{{print $1"\t"b}}'
           done; }} > {output.c2b}
        """


rule maxbin2:
    input:
        contigs=f"{OUT}/assembly/{{sample}}/contigs.fa",
        depth=f"{OUT}/bins/{{sample}}/depth.txt",
    output:
        c2b=f"{OUT}/bins/{{sample}}/maxbin2.contig2bin.tsv",
        binsdir=directory(f"{OUT}/bins/{{sample}}/maxbin2"),
    threads: config["threads"]
    conda:
        "../envs/maxbin2.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {output.binsdir}
        # MaxBin abundance = contig<TAB>mean-coverage (col 4 of the jgi depth table).
        awk 'NR>1{{print $1"\t"$4}}' {input.depth} > {output.binsdir}/abund.txt
        run_MaxBin.pl -contig {input.contigs} -abund {output.binsdir}/abund.txt \
          -out {output.binsdir}/maxbin -thread {threads} -min_contig_length 1000 \
          || true
        {{ for f in {output.binsdir}/*.fasta; do
             [ -e "$f" ] || continue
             b=$(basename "$f" .fasta)
             grep '^>' "$f" | sed 's/^>//' | awk -v b="$b" '{{print $1"\t"b}}'
           done; }} > {output.c2b}
        """


rule concoct:
    input:
        contigs=f"{OUT}/assembly/{{sample}}/contigs.fa",
        bam=f"{OUT}/mapping/{{sample}}.sorted.bam",
        bai=f"{OUT}/mapping/{{sample}}.sorted.bam.bai",
    output:
        c2b=f"{OUT}/bins/{{sample}}/concoct.contig2bin.tsv",
        binsdir=directory(f"{OUT}/bins/{{sample}}/concoct"),
    threads: config["threads"]
    conda:
        "../envs/concoct.yaml"
    params:
        work=lambda w: f"{OUT}/bins/{w.sample}/concoct_work",
    shell:
        r"""
        set -euo pipefail
        rm -rf {params.work}; mkdir -p {params.work} {output.binsdir}
        cut_up_fasta.py {input.contigs} -c 10000 -o 0 --merge_last \
          -b {params.work}/contigs_10k.bed > {params.work}/contigs_10k.fa
        concoct_coverage_table.py {params.work}/contigs_10k.bed {input.bam} \
          > {params.work}/coverage_table.tsv
        concoct --composition_file {params.work}/contigs_10k.fa \
          --coverage_file {params.work}/coverage_table.tsv \
          -b {params.work}/ -t {threads}
        merge_cutup_clustering.py {params.work}/clustering_gt1000.csv \
          > {params.work}/clustering_merged.csv
        extract_fasta_bins.py {input.contigs} {params.work}/clustering_merged.csv \
          --output_path {output.binsdir}
        {{ for f in {output.binsdir}/*.fa; do
             [ -e "$f" ] || continue
             b=$(basename "$f" .fa)
             grep '^>' "$f" | sed 's/^>//' | awk -v b="$b" '{{print $1"\t"b}}'
           done; }} > {output.c2b}
        """


def dastool_inputs(wildcards):
    return {
        b: f"{OUT}/bins/{wildcards.sample}/{b}.contig2bin.tsv"
        for b in ENABLED_BINNERS
    }


rule das_tool:
    """Consolidate the enabled binners into one non-redundant, scored bin set."""
    input:
        contigs=f"{OUT}/assembly/{{sample}}/contigs.fa",
        c2b=lambda w: list(dastool_inputs(w).values()),
    output:
        binsdir=directory(f"{OUT}/bins/{{sample}}/dastool_bins"),
        done=f"{OUT}/bins/{{sample}}/dastool.done",
    threads: config["threads"]
    conda:
        "../envs/dastool.yaml"
    params:
        labels=",".join(ENABLED_BINNERS),
        tsvs=lambda w: ",".join(dastool_inputs(w).values()),
        prefix=lambda w: f"{OUT}/bins/{w.sample}/dastool",
    shell:
        r"""
        set -euo pipefail
        mkdir -p {output.binsdir}
        # Drop empty binner tables so DAS_Tool never sees a zero-bin input.
        keep_i=""; keep_l=""
        IFS=',' read -ra TS <<< "{params.tsvs}"
        IFS=',' read -ra LB <<< "{params.labels}"
        for k in "${{!TS[@]}}"; do
          if [ -s "${{TS[$k]}}" ]; then
            keep_i="${{keep_i:+$keep_i,}}${{TS[$k]}}"
            keep_l="${{keep_l:+$keep_l,}}${{LB[$k]}}"
          fi
        done
        if [ -z "$keep_i" ]; then
          echo "das_tool: no non-empty bins for {wildcards.sample}" >&2
          touch {output.done}; exit 0
        fi
        DAS_Tool -i "$keep_i" -l "$keep_l" -c {input.contigs} \
          -o {params.prefix} --write_bins --score_threshold 0.5 -t {threads} || true
        # DAS_Tool writes bins to {{prefix}}_DASTool_bins/ ; collect them.
        if compgen -G "{params.prefix}_DASTool_bins/*.fa" > /dev/null; then
          cp {params.prefix}_DASTool_bins/*.fa {output.binsdir}/
        fi
        touch {output.done}
        """
