"""Input/output helpers: abundance profiles, model directories, manifests.

These functions define the on-disk contracts that glue the Snakemake workflow to
the Python library.  They deliberately avoid heavy dependencies (no cobra import
at module load) so they can be used from any pipeline phase.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd


def read_abundance(path: str | Path) -> Dict[str, float]:
    """Read a MAG abundance profile (produced by any taxonomic profiler).

    Accepts a 2-column TSV/CSV of ``mag_id`` and ``abundance``.  A header is
    optional; if the second column of the first row is non-numeric it is treated
    as a header and skipped.  Abundances are returned as-is (normalisation is the
    community's job).
    """
    path = Path(path)
    sep = "," if path.suffix.lower() == ".csv" else "\t"
    # sniff a header
    first = path.read_text().splitlines()[0].split(sep)
    header = 0
    if len(first) >= 2:
        try:
            float(first[1])
            header = None  # numeric -> no header
        except ValueError:
            header = 0
    df = pd.read_csv(path, sep=sep, header=header)
    if df.shape[1] < 2:
        raise ValueError(f"{path}: expected at least 2 columns (mag_id, abundance)")
    ids = df.iloc[:, 0].astype(str)
    vals = pd.to_numeric(df.iloc[:, 1], errors="coerce").fillna(0.0)
    return {str(i): float(v) for i, v in zip(ids, vals)}


def write_abundance(abundance: Dict[str, float], path: str | Path) -> None:
    path = Path(path)
    sep = "," if path.suffix.lower() == ".csv" else "\t"
    pd.DataFrame(
        {"mag_id": list(abundance), "abundance": list(abundance.values())}
    ).to_csv(path, sep=sep, index=False)


def list_models(model_dir: str | Path, patterns: tuple[str, ...] = ("*.xml", "*.json", "*.xml.gz")) -> List[Path]:
    """List GEM files in a directory across the common SBML/JSON extensions."""
    model_dir = Path(model_dir)
    out: List[Path] = []
    for pat in patterns:
        out.extend(sorted(model_dir.glob(pat)))
    return out


def write_manifest(path: str | Path, **fields) -> None:
    """Write a small JSON manifest describing a pipeline step's output."""
    Path(path).write_text(json.dumps(fields, indent=2, default=str))


def read_manifest(path: str | Path) -> dict:
    return json.loads(Path(path).read_text())


def load_simulation_result(outdir: str | Path):
    """Reload a :class:`~muode.dfba.SimulationResult` written by ``to_csv``."""
    from muode.dfba import SimulationResult

    outdir = Path(outdir)
    biomass = pd.read_csv(outdir / "biomass.csv", index_col="time_h")
    metabolites = pd.read_csv(outdir / "metabolites.csv", index_col="time_h")
    growth = pd.read_csv(outdir / "growth_rates.csv", index_col="time_h")
    return SimulationResult(
        times=biomass.index.values,
        biomass=biomass,
        metabolites=metabolites,
        growth_rates=growth,
    )
