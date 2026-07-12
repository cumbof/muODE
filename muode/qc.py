"""Quality control -- for MAGs (upstream) and for reconstructed models.

Good simulations require good inputs.  Two gates:

* **MAG QC (Phase 0).**  :func:`checkm2` / :func:`filter_mags` drop low-quality
  or contaminated bins before any modelling effort is spent on them.

* **Model QC.**  :func:`sanity_check_model` runs cheap, decisive tests that
  catch the GEM artefacts most likely to wreck a dynamic simulation: mass /
  charge imbalance and -- crucially -- *energy-generating cycles* (an ability to
  make ATP from nothing, which would let a species "grow" on an empty medium).
  :func:`memote_report` wraps the community-standard memote test suite for a full
  report.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional

from muode.external import require, run


# ---------------------------------------------------------------------------
# MAG quality gate
# ---------------------------------------------------------------------------


def checkm2(mags_dir: str | Path, output_dir: str | Path, extension: str = "fna", threads: int = 4) -> Path:
    """Assess MAG completeness/contamination with CheckM2; return the report TSV."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    bin_ = require("checkm2", env_hint="mamba install -c bioconda checkm2")
    run([bin_, "predict", "--input", str(mags_dir), "--output-directory", str(output_dir),
         "-x", extension, "--threads", str(threads), "--force"],
        log=output_dir / "checkm2.log")
    return output_dir / "quality_report.tsv"


def filter_mags(report_tsv: str | Path, min_completeness: float = 50.0,
                max_contamination: float = 10.0) -> List[str]:
    """Return MAG ids passing the completeness/contamination thresholds."""
    import pandas as pd

    df = pd.read_csv(report_tsv, sep="\t")
    cols = {c.lower(): c for c in df.columns}
    name_c = cols.get("name") or cols.get("bin id") or df.columns[0]
    comp_c = cols.get("completeness")
    cont_c = cols.get("contamination")
    keep = df[(df[comp_c] >= min_completeness) & (df[cont_c] <= max_contamination)]
    return keep[name_c].astype(str).tolist()


# ---------------------------------------------------------------------------
# Reconstruction-at-scale reporting
# ---------------------------------------------------------------------------


def summarize_reconstruction(qc_json_paths):
    """Aggregate per-MAG refine/QC JSONs into one tidy reconstruction table.

    This is the deliverable you read after a large run: one row per MAG telling
    you whether its model grows, whether it carries an energy-generating cycle,
    how many reactions gap-filling added, and -- via the derived ``simulatable``
    column -- whether it is safe to include in the community.  Returns a
    :class:`pandas.DataFrame` (empty if no inputs).
    """
    import json

    import pandas as pd

    rows = []
    for path in qc_json_paths:
        path = Path(path)
        try:
            d = json.loads(path.read_text())
        except (OSError, json.JSONDecodeError):
            rows.append({"mag": path.stem.replace(".qc", ""), "grows_now": False,
                         "error": "missing or unreadable QC report"})
            continue
        qc = d.get("qc", {}) or {}
        rows.append({
            "mag": d.get("model", path.stem.replace(".qc", "")),
            "grew_initially": d.get("grew_initially"),
            "n_reactions_added": len(d.get("reactions_added") or []),
            "grows_now": bool(d.get("grows_now")),
            # growth on the *simulation medium* vs with every exchange open.  A
            # model that grows only on the latter will flatline in the community.
            "growth_on_diet": d.get("growth_on_diet"),
            "growth_on_complete_medium": d.get("growth_on_complete_medium"),
            "n_mass_unbalanced": qc.get("n_mass_unbalanced"),
            "energy_generating_cycle": qc.get("energy_generating_cycle"),
            "qc_passed": qc.get("passed"),
        })
    df = pd.DataFrame(rows)
    if not df.empty:
        # the diet columns are absent when refine ran without one: don't carry
        # a column of blanks into the report
        df = df.dropna(axis=1, how="all")
        df = df.sort_values("mag").reset_index(drop=True)
        egc = df.get("energy_generating_cycle")
        df["simulatable"] = df["grows_now"].fillna(False) & (egc.fillna(False) == False)  # noqa: E712
    return df


def simulatable_mags(summary_tsv: str | Path) -> List[str]:
    """Read a reconstruction summary TSV and return the ids safe to simulate."""
    import pandas as pd

    df = pd.read_csv(summary_tsv, sep="\t")
    if "simulatable" not in df.columns or df.empty:
        return []
    return df.loc[df["simulatable"] == True, "mag"].astype(str).tolist()  # noqa: E712


# ---------------------------------------------------------------------------
# Model quality
# ---------------------------------------------------------------------------


def sanity_check_model(model, medium: Optional[Dict[str, float]] = None) -> dict:
    """Fast structural sanity checks on a GEM.

    Returns a report dict; ``passed`` is False if a likely-fatal artefact (an
    energy-generating cycle) is detected.
    """
    report: dict = {}

    # mass / charge balance (ignore boundary and objective pseudo-reactions;
    # biomass/objective reactions are intentionally not mass-balanced)
    objective_ids = {
        getattr(v, "name", str(v)) for v in model.objective.variables
    } | {r.id for r in model.reactions if r.objective_coefficient != 0}
    unbalanced = []
    for rxn in model.reactions:
        if rxn.boundary or rxn.id in objective_ids or "biomass" in rxn.id.lower():
            continue
        try:
            imbalance = rxn.check_mass_balance()
        except Exception:
            imbalance = {"error": 1}
        if imbalance:
            unbalanced.append(rxn.id)
    report["n_mass_unbalanced"] = len(unbalanced)
    report["mass_unbalanced_examples"] = unbalanced[:10]

    # energy-generating cycle: with all uptake closed, can ATP be hydrolysed?
    report["energy_generating_cycle"] = _has_energy_generating_cycle(model)

    # baseline growth
    if medium is not None:
        with model:
            model.medium = {k: v for k, v in medium.items() if k in [r.id for r in model.exchanges]}
            report["growth_on_medium"] = float(model.slim_optimize() or 0.0)

    report["passed"] = not report["energy_generating_cycle"]
    return report


def _has_energy_generating_cycle(model) -> bool:
    """True if the model can make ATP with every nutrient uptake switched off."""
    try:
        with model:
            for ex in model.exchanges:
                ex.lower_bound = 0.0  # forbid all uptake
            atpm = None
            for rid in ("ATPM", "ATP_maintenance", "NGAM"):
                if rid in model.reactions:
                    atpm = model.reactions.get_by_id(rid)
                    break
            if atpm is None:
                return False
            atpm.lower_bound = 0.0
            model.objective = atpm
            return (model.slim_optimize() or 0.0) > 1e-6
    except Exception:
        return False


def memote_report(model_path: str | Path, output_html: str | Path) -> Path:
    """Generate a full memote quality report for a model."""
    output_html = Path(output_html)
    output_html.parent.mkdir(parents=True, exist_ok=True)
    memote = require("memote", env_hint="pip install memote")
    run([memote, "report", "snapshot", "--filename", str(output_html), str(model_path)],
        log=output_html.with_suffix(".log"))
    return output_html
