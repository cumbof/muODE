"""Pipeline configuration shared by the CLI and the Snakemake workflow.

A single :class:`PipelineConfig` (loadable from YAML) describes a whole muODE
run, from MAG quality gates through reconstruction, kinetics and the dynamic
simulation.  Keeping it in one typed place means the CLI and the workflow can
never drift out of sync.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Optional


@dataclass
class PipelineConfig:
    # --- inputs / outputs --------------------------------------------------
    mags_dir: str = "data/raw_mags"
    abundance: Optional[str] = None  # precomputed 2-column abundance TSV, or None
    outdir: str = "results"

    # --- MAG quality gate (Phase 0) ---------------------------------------
    run_checkm2: bool = False        # gate MAGs on CheckM2 (off by default; see config)
    min_completeness: float = 50.0   # CheckM2 completeness (%)
    max_contamination: float = 10.0  # CheckM2 contamination (%)

    # --- reconstruction (Phase 1) -----------------------------------------
    gene_caller: str = "pyrodigal"   # pyrodigal | prodigal
    engine: str = "carveme"          # carveme | gapseq | stub
    universe: str = "bacteria"       # CarveMe universe / template
    gapfill_media: Optional[str] = None  # e.g. "M9" or a media db entry

    # --- model QC (Phase 3b) ----------------------------------------------
    universal_model: Optional[str] = None  # SBML universal for LP gap-filling
    run_memote: bool = False         # generate full per-model memote reports

    # --- kinetics (Phase 3) -----------------------------------------------
    predict_kinetics: bool = False   # write per-MAG predicted Km + kcat
    predictor: str = "heuristic"     # heuristic (no deps) | dlkcat | km-ml (ml extra)
    enzyme_constraints: bool = False # apply GECKO-lite kcat caps at simulation time
    default_vmax: float = 10.0
    default_km: float = 0.01

    # --- environment -------------------------------------------------------
    # WHAT KIND OF PLACE IS THIS?  The `diet` is the environment's chemistry -- what
    # is there and how much.  `environment` is everything else: how fast anything
    # could plausibly grow here, and whether a no-growth diagnosis is allowed to
    # propose oxygen.  Both used to be bare constants in qc.py that silently assumed a
    # human colon for every sample muODE ever ran.  See muode/environment.py; muODE
    # ships only the environments it can justify with a citation.
    environment: str = "human_gut"   # human_gut | generic_anaerobic | generic_aerobic

    # --- simulation (Phase 4) ---------------------------------------------
    diet: str = "western_gut"        # preset name OR path to a diet CSV
    total_biomass: float = 0.01
    t_end: float = 24.0
    dt: float = 0.1
    death_rate: float = 0.0
    dilution_rate: float = 0.0

    # --- scale & validation -----------------------------------------------
    max_species: Optional[int] = None        # keep only the N most abundant MAGs
    min_abundance: Optional[float] = None     # drop MAGs below this relative abundance
    abundance_coverage: Optional[float] = None  # keep top MAGs reaching this cumulative abundance
    benchmark: Optional[str] = None           # benchmark expectation YAML for validation

    # --- solver / compute --------------------------------------------------
    solver: str = "highs"            # highs (free) | gurobi | cplex
    threads: int = 4

    # --- perturbation (Phase 5, optional) ---------------------------------
    perturbation: Optional[dict] = None

    # -- IO -----------------------------------------------------------------
    @classmethod
    def from_yaml(cls, path: str | Path) -> "PipelineConfig":
        import yaml

        data = yaml.safe_load(Path(path).read_text()) or {}
        known = {f for f in cls.__dataclass_fields__}  # type: ignore[attr-defined]
        return cls(**{k: v for k, v in data.items() if k in known})

    def to_yaml(self, path: str | Path) -> None:
        import yaml

        Path(path).write_text(yaml.safe_dump(asdict(self), sort_keys=False))

    def to_dict(self) -> dict:
        return asdict(self)
