# =============================================================================
# Phase 0 -- read QC: adapter/quality trimming (fastp) + optional host removal
# =============================================================================

rule fastp:
    """Adapter + quality trim paired reads; drop reads shorter than min_read_len."""
    input:
        unpack(raw_reads),
    output:
        r1=f"{OUT}/qc/trimmed/{{sample}}_1.fastq.gz",
        r2=f"{OUT}/qc/trimmed/{{sample}}_2.fastq.gz",
        json=f"{OUT}/qc/trimmed/{{sample}}.fastp.json",
        html=f"{OUT}/qc/trimmed/{{sample}}.fastp.html",
    threads: config["threads"]
    conda:
        "../envs/qc.yaml"
    params:
        minlen=config["min_read_len"],
    shell:
        r"""
        fastp \
          -i {input.r1} -I {input.r2} \
          -o {output.r1} -O {output.r2} \
          --detect_adapter_for_pe \
          --length_required {params.minlen} \
          --thread {threads} \
          --json {output.json} --html {output.html}
        """


rule host_removal:
    """Map trimmed reads to a host bowtie2 index; keep pairs with BOTH mates unmapped.

    Only used when config remove_host: true. `host_index` must be a bowtie2 index
    PREFIX (the path you'd pass to `bowtie2 -x`), e.g. a human GRCh38 index.

    Privacy-safe decontamination: `samtools view -f 12` retains a pair only when
    both mates are unmapped, so a pair is dropped if EITHER mate maps to the host.
    This is stricter than bowtie2 `--un-conc` (which can leak a pair when only one
    mate aligns) and is the right default for human clinical samples.
    """
    input:
        r1=f"{OUT}/qc/trimmed/{{sample}}_1.fastq.gz",
        r2=f"{OUT}/qc/trimmed/{{sample}}_2.fastq.gz",
    output:
        r1=f"{OUT}/qc/hostfree/{{sample}}_1.fastq.gz",
        r2=f"{OUT}/qc/hostfree/{{sample}}_2.fastq.gz",
    threads: config["threads"]
    conda:
        "../envs/qc.yaml"
    params:
        index=config.get("host_index", ""),
        prefix=lambda w: f"{OUT}/qc/hostfree/{w.sample}",
    shell:
        r"""
        set -euo pipefail
        if [ -z "{params.index}" ]; then
          echo "host_removal: config remove_host is true but host_index is empty" >&2
          exit 1
        fi
        # -f 12: read unmapped (0x4) AND mate unmapped (0x8);  -F 256: no secondary.
        # Name-sort so samtools fastq can split mates back into a clean pair.
        bowtie2 -p {threads} -x {params.index} \
          -1 {input.r1} -2 {input.r2} 2> {params.prefix}.bowtie2.log \
          | samtools view -@ {threads} -bS -f 12 -F 256 - \
          | samtools sort -@ {threads} -n -O bam - \
          | samtools fastq -@ {threads} \
              -1 {output.r1} -2 {output.r2} -0 /dev/null -s /dev/null -n -
        """
