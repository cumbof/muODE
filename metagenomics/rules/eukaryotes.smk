# =============================================================================
# Phase 6 -- eukaryote recovery: Tiara (contig domain) -> MetaBAT2 -> EukCC
# =============================================================================
# Honest caveat: recovering eukaryotic MAGs from shotgun metagenomes is markedly
# harder and lower-yield than the prokaryotic track (large, repeat-rich, low-
# coverage genomes). This gives a real, standard pipeline for it -- Tiara flags
# eukaryotic contigs, MetaBAT2 bins them, EukCC estimates quality -- but expect
# few, partial euk MAGs unless a euk is genuinely abundant. Enable with
# `recover_eukaryotes: true` and set `eukcc_db`.

rule tiara_classify:
    """Classify contigs by domain; extract the eukaryotic contigs (length-filtered)."""
    input:
        contigs=f"{OUT}/assembly/{{sample}}/contigs.fa",
    output:
        classifications=f"{OUT}/eukaryotes/{{sample}}/tiara.txt",
        euk=f"{OUT}/eukaryotes/{{sample}}/euk_contigs.fa",
    threads: config["threads"]
    conda:
        "../envs/tiara.yaml"
    params:
        work=lambda w: f"{OUT}/eukaryotes/{w.sample}/tiara_work",
        minlen=config.get("euk_min_contig_len", 5000),
    shell:
        r"""
        set -euo pipefail
        rm -rf {params.work}; mkdir -p {params.work}
        abs_contigs=$(readlink -f {input.contigs})
        # tiara writes --to_fasta files into the CWD, so run it inside the workdir.
        ( cd {params.work} && tiara -i "$abs_contigs" \
            -o classifications.txt --to_fasta euk -t {threads} -m {params.minlen} )
        cp {params.work}/classifications.txt {output.classifications}
        # tiara writes the euk fasta as eukarya_*.fasta in its workdir.
        if compgen -G "{params.work}/eukarya_*.fasta" > /dev/null; then
          seqkit seq -m {params.minlen} {params.work}/eukarya_*.fasta > {output.euk}
        else
          : > {output.euk}
        fi
        """


rule euk_bin:
    """Bin the eukaryotic contigs with MetaBAT2 using coverage restricted to them."""
    input:
        euk=f"{OUT}/eukaryotes/{{sample}}/euk_contigs.fa",
        depth=f"{OUT}/bins/{{sample}}/depth.txt",
    output:
        binsdir=directory(f"{OUT}/eukaryotes/{{sample}}/bins"),
    threads: config["threads"]
    conda:
        "../envs/metabat2.yaml"
    params:
        work=lambda w: f"{OUT}/eukaryotes/{w.sample}",
    shell:
        r"""
        set -euo pipefail
        mkdir -p {output.binsdir}
        if [ ! -s {input.euk} ]; then
          echo "euk_bin: no eukaryotic contigs for {wildcards.sample}" >&2
          exit 0
        fi
        # Restrict the shared depth table to euk contig ids.
        grep '^>' {input.euk} | sed 's/^>//' | awk '{{print $1}}' \
          > {params.work}/euk_ids.txt
        awk 'NR==1{{print; next}} FNR==NR{{keep[$1]=1; next}} ($1 in keep)' \
          {params.work}/euk_ids.txt {input.depth} > {params.work}/euk_depth.txt
        metabat2 -i {input.euk} -a {params.work}/euk_depth.txt \
          -o {output.binsdir}/euk -m 100000 -t {threads} || true
        """


rule eukcc:
    """Estimate completeness/contamination of euk bins (EukCC lineage-aware)."""
    input:
        binsdir=f"{OUT}/eukaryotes/{{sample}}/bins",
    output:
        report=f"{OUT}/eukaryotes/{{sample}}/eukcc.csv",
    threads: config["threads"]
    conda:
        "../envs/eukcc.yaml"
    params:
        outdir=lambda w: f"{OUT}/eukaryotes/{w.sample}/eukcc_out",
        db=config.get("eukcc_db", ""),
    shell:
        r"""
        set -euo pipefail
        if ! compgen -G "{input.binsdir}/*.fa" > /dev/null; then
          printf 'bin\tcompleteness\tcontamination\n' > {output.report}
          exit 0
        fi
        if [ -z "{params.db}" ]; then
          echo "eukcc: recover_eukaryotes is true but eukcc_db is empty" >&2
          exit 1
        fi
        mkdir -p {params.outdir}
        eukcc folder --out {params.outdir} --threads {threads} \
          --db {params.db} --suffix .fa {input.binsdir}
        cp {params.outdir}/eukcc.csv {output.report}
        """


rule select_euk_mags:
    input:
        report=f"{OUT}/eukaryotes/{{sample}}/eukcc.csv",
        binsdir=f"{OUT}/eukaryotes/{{sample}}/bins",
    output:
        magdir=directory(f"{OUT}/mags_per_sample/{{sample}}/euk"),
        quality=f"{OUT}/mags_per_sample/{{sample}}/euk.quality.tsv",
    conda:
        "../envs/base.yaml"
    params:
        mincomp=config.get("eukcc_min_completeness", 50.0),
        maxcont=config.get("eukcc_max_contamination", 10.0),
    shell:
        r"""
        python scripts/select_mags.py \
          --report {input.report} --bins {input.binsdir} \
          --sample {wildcards.sample} --tag euk \
          --outdir {output.magdir} --quality-out {output.quality} \
          --min-completeness {params.mincomp} --max-contamination {params.maxcont}
        """
