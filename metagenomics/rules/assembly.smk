# =============================================================================
# Phase 1 -- per-sample assembly (MEGAHIT default / metaSPAdes) + read mapping
# =============================================================================
# Per-sample (not co-) assembly is the default: it scales linearly across a large
# cohort, keeps strain variation resolvable, and dRep later collapses the same
# genome recovered from several samples into one catalogue entry. MEGAHIT is the
# default assembler because it handles deep gut metagenomes in a fraction of the
# RAM metaSPAdes needs; switch with `assembler: metaspades` when memory allows.

rule assemble:
    """Assemble one sample; emit length-filtered contigs at a uniform path."""
    input:
        unpack(clean_reads),
    output:
        contigs=f"{OUT}/assembly/{{sample}}/contigs.fa",
    threads: config["threads"]
    conda:
        "../envs/assembly.yaml"
    params:
        assembler=config["assembler"],
        minlen=config["min_contig_len"],
        workdir=lambda w: f"{OUT}/assembly/{w.sample}/work",
        memgb=config.get("spades_mem_gb", 250),
    shell:
        r"""
        set -euo pipefail
        rm -rf {params.workdir}
        case "{params.assembler}" in
          megahit)
            megahit -1 {input.r1} -2 {input.r2} \
              -t {threads} --min-contig-len {params.minlen} \
              -o {params.workdir}
            raw={params.workdir}/final.contigs.fa
            ;;
          metaspades)
            spades.py --meta -1 {input.r1} -2 {input.r2} \
              -t {threads} -m {params.memgb} -o {params.workdir}
            raw={params.workdir}/contigs.fasta
            ;;
          *)
            echo "unknown assembler: {params.assembler}" >&2; exit 1 ;;
        esac
        # Uniform, length-filtered output with sample-prefixed contig names
        # (prefix keeps contig ids unique if assemblies are ever pooled).
        seqkit seq -m {params.minlen} "$raw" \
          | seqkit replace -p '^' -r '{wildcards.sample}__' \
          > {output.contigs}
        """


rule map_reads:
    """Map a sample's own reads back to its assembly -> sorted, indexed BAM.

    Provides the per-contig coverage every binner needs.
    """
    input:
        contigs=f"{OUT}/assembly/{{sample}}/contigs.fa",
        r1=lambda w: clean_reads(w)["r1"],
        r2=lambda w: clean_reads(w)["r2"],
    output:
        bam=f"{OUT}/mapping/{{sample}}.sorted.bam",
        bai=f"{OUT}/mapping/{{sample}}.sorted.bam.bai",
    threads: config["threads"]
    conda:
        "../envs/qc.yaml"
    params:
        idx=lambda w: f"{OUT}/mapping/{w.sample}.idx",
    shell:
        r"""
        set -euo pipefail
        bowtie2-build --threads {threads} {input.contigs} {params.idx} > /dev/null 2>&1
        bowtie2 -p {threads} -x {params.idx} -1 {input.r1} -2 {input.r2} \
          | samtools view -@ {threads} -bS - \
          | samtools sort -@ {threads} -o {output.bam} -
        samtools index {output.bam}
        rm -f {params.idx}*.bt2 {params.idx}*.bt2l
        """
