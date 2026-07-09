# =============================================================================
# Phase 4 -- WITHIN-SAMPLE dereplication (dRep) -> one catalogue per sample
# =============================================================================
# Samples are NOT pooled. Each sample is a distinct biological community -- in a
# case/control design the samples are different patients -- so a genome recovered
# from patient B must never end up in patient A's simulated community.
#
# Pooling would also be lossy in a subtler way: dRep keeps ONE representative per
# 0.95-ANI cluster, so if two patients carry different strains of the same species
# the catalogue silently substitutes one patient's strain for the other's.
#
# Dereplication is therefore run *within* each sample, where it does the job it
# should: collapsing near-identical bins recovered by the binning ensemble (and
# the prokaryotic + eukaryotic tracks) into one non-redundant per-sample
# catalogue. Pre-computed completeness/contamination go to dRep via genomeInfo so
# it does not re-run CheckM (and can rank euk MAGs it otherwise couldn't score).

def sample_quality_tables(wildcards):
    tabs = [f"{OUT}/mags_per_sample/{wildcards.sample}/prok.quality.tsv"]
    if DO_EUK:
        tabs.append(f"{OUT}/mags_per_sample/{wildcards.sample}/euk.quality.tsv")
    return tabs


def sample_mag_dirs(wildcards):
    dirs = [f"{OUT}/mags_per_sample/{wildcards.sample}/prok"]
    if DO_EUK:
        dirs.append(f"{OUT}/mags_per_sample/{wildcards.sample}/euk")
    return dirs


rule dereplicate:
    input:
        quality=sample_quality_tables,
        magdirs=sample_mag_dirs,
    output:
        catalogue=directory(f"{OUT}/catalogue/{{sample}}/mags"),
        merged_quality=f"{OUT}/catalogue/{{sample}}/all_mags.quality.tsv",
        done=f"{OUT}/catalogue/{{sample}}/dereplicate.done",
    threads: config["threads"]
    conda:
        "../envs/drep.yaml"
    params:
        ani=config.get("derep_ani", 0.95),
        work=lambda w: f"{OUT}/catalogue/{w.sample}/drep",
        pool=lambda w: f"{OUT}/catalogue/{w.sample}/pool",
    shell:
        r"""
        set -euo pipefail
        rm -rf {params.pool} {params.work}
        mkdir -p {params.pool} {output.catalogue}

        # Pool this sample's MAGs only (prok + euk); ids are <sample>__tag.k.
        found=0
        for d in {input.magdirs}; do
          if compgen -G "$d/*.fa" > /dev/null; then
            cp "$d"/*.fa {params.pool}/ && found=1
          fi
        done

        # Merge quality tables + build dRep genomeInfo.
        python scripts/collect_quality.py \
          --quality {input.quality} \
          --merged-out {output.merged_quality} \
          --genomeinfo-out {params.pool}/genomeInfo.csv

        if [ "$found" -eq 0 ]; then
          echo "dereplicate: no MAGs passed QC for {wildcards.sample}" >&2
          touch {output.done}; exit 0
        fi

        n=$(ls {params.pool}/*.fa | wc -l)
        if [ "$n" -eq 1 ]; then
          # dRep needs >=2 genomes; a lone MAG is trivially non-redundant.
          cp {params.pool}/*.fa {output.catalogue}/
          touch {output.done}; exit 0
        fi

        dRep dereplicate {params.work} \
          -g {params.pool}/*.fa \
          -p {threads} -sa {params.ani} \
          --genomeInfo {params.pool}/genomeInfo.csv \
          -comp 0 -con 1000
        cp {params.work}/dereplicated_genomes/*.fa {output.catalogue}/
        touch {output.done}
        """
