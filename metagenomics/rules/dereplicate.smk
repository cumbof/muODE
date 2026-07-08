# =============================================================================
# Phase 4 -- cross-sample dereplication (dRep) -> non-redundant MAG catalogue
# =============================================================================
# The same genome is usually recovered from many samples. dRep collapses those
# into one representative per species-level cluster, yielding the non-redundant
# catalogue that becomes muODE's `mags_dir`. Prokaryotic and (optional) eukaryotic
# MAGs are pooled; pre-computed completeness/contamination are handed to dRep via
# genomeInfo so it does not re-run CheckM (and can rank euk MAGs it otherwise
# couldn't score).

def all_quality_tables(_):
    tabs = [f"{OUT}/mags_per_sample/{s}/prok.quality.tsv" for s in SAMPLE_IDS]
    if DO_EUK:
        tabs += [f"{OUT}/mags_per_sample/{s}/euk.quality.tsv" for s in SAMPLE_IDS]
    return tabs


def all_mag_dirs(_):
    dirs = [f"{OUT}/mags_per_sample/{s}/prok" for s in SAMPLE_IDS]
    if DO_EUK:
        dirs += [f"{OUT}/mags_per_sample/{s}/euk" for s in SAMPLE_IDS]
    return dirs


rule dereplicate:
    input:
        quality=all_quality_tables,
        magdirs=all_mag_dirs,
    output:
        catalogue=directory(f"{OUT}/catalogue/mags"),
        merged_quality=f"{OUT}/catalogue/all_mags.quality.tsv",
        done=f"{OUT}/catalogue/dereplicate.done",
    threads: config["threads"]
    conda:
        "../envs/drep.yaml"
    params:
        ani=config.get("derep_ani", 0.95),
        work=f"{OUT}/catalogue/drep",
        pool=f"{OUT}/catalogue/pool",
    shell:
        r"""
        set -euo pipefail
        rm -rf {params.pool} {params.work}
        mkdir -p {params.pool} {output.catalogue}

        # Pool every passing MAG (unique ids guaranteed by <sample>__tag.k naming).
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
          echo "dereplicate: no MAGs passed QC in any sample" >&2
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
