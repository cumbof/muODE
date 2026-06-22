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
    abundance: Optional[str] = None
    outdir: str = "results"

    # --- MAG quality gate (Phase 0) ---------------------------------------
    min_completeness: float = 50.0   # CheckM2 completeness (%)
    max_contamination: float = 10.0  # CheckM2 contamination (%)

    # --- reconstruction (Phase 1) -----------------------------------------
    gene_caller: str = "pyrodigal"   # pyrodigal | prodigal
    engine: str = "carveme"          # carveme | gapseq
    universe: str = "bacteria"       # CarveMe universe / template
    gapfill_media: Optional[str] = None  # e.g. "M9" or a media db entry

    # --- kinetics (Phase 3) -----------------------------------------------
    predict_kinetics: bool = False   # opt-in DLKcat/Km refinement layer
    default_vmax: float = 10.0
    default_km: float = 0.01

    # --- simulation (Phase 4) ---------------------------------------------
    diet: str = "western_gut"        # preset name OR path to a diet CSV
    total_biomass: float = 0.01
    t_end: float = 24.0
    dt: float = 0.1
    death_rate: float = 0.0
    dilution_rate: float = 0.0

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
