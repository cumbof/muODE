"""Quality control -- for MAGs (upstream) and for reconstructed models.

Good simulations require good inputs.  Two gates:

* **MAG QC (Phase 0).**  :func:`checkm2` / :func:`filter_mags` drop low-quality
  or contaminated bins before any modelling effort is spent on them.

* **Model QC.**  :func:`sanity_check_model` runs cheap, decisive tests that
  catch the GEM artefacts most likely to wreck a dynamic simulation: mass /
  charge imbalance, *energy-generating cycles* (an ability to make ATP from
  nothing, which would let a species "grow" on an empty medium), and growth at a
  rate no organism achieves (see :data:`MAX_PLAUSIBLE_GROWTH`).
  :func:`memote_report` wraps the community-standard memote test suite for a full
  report.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional

from muode.external import require, run


#: Growth rates above this (1/h) are not biology, they are a modelling artefact.
#:
#: Gut anaerobes run 0.1-0.5/h; the fastest organism ever measured (*Vibrio
#: natriegens*, rich aerobic medium) tops out near 4/h, and no gut commensal comes
#: close.  A GEM reporting more than this on a defined medium is being handed more
#: nutrient than a cell can physically consume -- typically because every exchange
#: in a large diet is opened at the same uniform Vmax, so the model eats 80+ carbon
#: sources at once.  The number is a *ceiling on the absurd*, deliberately well
#: above any real gut organism: tripping it means the medium is wrong, not that the
#: species is unusually fast.
MAX_PLAUSIBLE_GROWTH = 1.5


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
            "implausible_growth": qc.get("implausible_growth"),
            # Why a zero is a zero.  `medium_gap` means the model is dead because the
            # medium never fed it -- its 0.0 is an artefact, not a finding about the
            # organism -- and `no_growth_missing` names what would fix it.  A gap
            # buried in a per-MAG JSON is a gap nobody reads.
            "no_growth_verdict": (qc.get("no_growth") or {}).get("verdict"),
            "no_growth_fixed_by_any_of": ";".join(
                (qc.get("no_growth") or {}).get("rescued_by_any_of") or []
            ) or None,
            "no_growth_fixed_by_all_of": ";".join(
                (qc.get("no_growth") or {}).get("rescued_by_all_of") or []
            ) or None,
            "qc_passed": qc.get("passed"),
        })
    df = pd.DataFrame(rows)
    if not df.empty:
        # the diet columns are absent when refine ran without one: don't carry
        # a column of blanks into the report
        df = df.dropna(axis=1, how="all")
        df = df.sort_values("mag").reset_index(drop=True)
        egc = df.get("energy_generating_cycle")
        # NB: `implausible_growth` deliberately does NOT gate `simulatable`.  It
        # indicts the *medium*, not the model -- every model on a too-rich diet
        # trips it -- so dropping those models would delete most of the community
        # and hide the cause.  It is surfaced, loudly, and left in.
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


def sanity_check_model(
    model,
    medium: Optional[Dict[str, float]] = None,
    growth_rate: Optional[float] = None,
    diet=None,
    kinetics=None,
) -> dict:
    """Fast structural sanity checks on a GEM.

    Returns a report dict; ``passed`` is False if a likely-fatal artefact is
    detected: an energy-generating cycle, or -- when ``growth_rate`` is supplied
    (the model's growth on the simulation medium) -- a rate above
    :data:`MAX_PLAUSIBLE_GROWTH`.

    Pass ``diet`` and a model that does **not** grow, and the report also carries a
    :func:`diagnose_no_growth` verdict under ``no_growth``.  This is not optional
    politeness: a zero looks the same whether the genome lacks the pathways or the
    medium lacks a nutrient the biomass demands, and only one of those is a result.
    Which nutrient it is depends entirely on which genomes you fed in, so the check
    has to run every time rather than being remembered.
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
        if growth_rate is None:
            growth_rate = report["growth_on_medium"]

    # a rate no organism achieves.  This is not a slow-growing-species judgement
    # call: it is a rate faster than any cell that has been measured, so it can
    # only be the medium handing out more nutrient than a cell can consume.
    if growth_rate is not None:
        report["growth_rate"] = float(growth_rate)
        report["implausible_growth"] = bool(growth_rate > MAX_PLAUSIBLE_GROWTH)

    # A model that does not grow owes us a reason.  Without this, "0.0" in the
    # results table means either "this organism cannot do it" or "we forgot to feed
    # it", and the two are indistinguishable.
    if diet is not None and (growth_rate is None or growth_rate <= NO_GROWTH_TOL):
        report["no_growth"] = diagnose_no_growth(model, diet, kinetics)

    report["passed"] = not report["energy_generating_cycle"] and not report.get(
        "implausible_growth", False
    )
    return report


#: Never proposed as a rescue, even when it would restore growth.
#:
#: A diagnosis is allowed to say "your medium is missing a nutrient".  It is NOT
#: allowed to say "your anaerobic medium should be aerobic" -- that overturns the
#: environment's premise rather than completing it, and oxygen will "rescue" almost
#: any GEM, so it would drown every real finding.  If a diet omits o2_e that is a
#: statement about the environment, and the diagnosis must respect it.  A diet that
#: DOES supply oxygen never reaches this: it is not a missing exchange.
FORBIDDEN_RESCUES = frozenset({"o2_e"})

#: Growth below this is "does not grow" for diagnostic purposes.
NO_GROWTH_TOL = 1e-6


def diagnose_no_growth(
    model,
    diet,
    kinetics=None,
    forbid: frozenset = FORBIDDEN_RESCUES,
    rescue_bound: float = 10.0,
) -> dict:
    """Why does this model not grow on this diet -- biology, or a hole in the medium?

    A zero growth rate has two causes that look identical in a results table, and
    conflating them is how a medium artefact gets reported as a biological finding:

    * **the genome really lacks the pathways** -- a true negative, and the answer;
    * **the medium lacks a nutrient this model's biomass requires** -- an artefact.
      The model is dead on arrival, and its zero says nothing about the organism.

    The second is not hypothetical and it is not specific to one nutrient.  DM38
    supplies Fe(II) and no Fe(III), which is chemically correct for a reducing
    anaerobic broth, and iJO1366 cannot grow on it *at any Vmax* because its biomass
    demands ``fe3_e``.  Whatever genomes you feed muODE, some of them may require
    something the medium does not carry, and nobody will be looking for it.  So the
    framework has to look.

    Returns a report with ``verdict``:

    ``grows``
        It does grow; nothing to diagnose.
    ``medium_gap``
        It grows only if the medium supplies more.  ``rescued_by_any_of`` lists every
        single metabolite that *on its own* restores growth; if no single one does,
        ``rescued_by_all_of`` gives a minimal combination that does.  Either way the
        zero is an artefact: do not report it as biology.
    ``model_cannot_grow``
        It does not grow even with every one of its own exchanges wide open.  The
        medium is exonerated: this is a broken or incomplete reconstruction.

    **Why a menu and not one answer.**  The obvious implementation -- shrink a full
    rescue set until every member is load-bearing -- returns *a* minimal set, and
    when the model has redundant routes it silently picks one of them by iteration
    order.  On iJO1366/DM38 that returns nitric oxide: the model's only anaerobic
    route to cytoplasmic Fe(III) is ``FESD2s``, which lets NO *destroy its own
    iron-sulfur clusters* to liberate iron.  That is a true minimal rescue, a solver
    artefact, and useless to a human -- while the honest answer, ferric iron, sits
    right next to it in the same equivalence class and loses a coin-toss on
    alphabetical order.

    So we enumerate the whole class instead.  It is deterministic, it is complete,
    and it lets a person see that ``EX_fe3_e`` and ``EX_no_e`` fix the same hole and
    pick the one that is chemistry rather than damage.  A tool cannot make that
    judgement; it can refuse to hide it.
    """
    from muode.media import diet_medium

    exchanges = {r.id for r in model.exchanges}
    supplied = diet_medium(model, diet, kinetics)
    base = {k: v for k, v in supplied.items() if k in exchanges}

    with model:
        model.medium = base
        growth = float(model.slim_optimize() or 0.0)
    if growth > NO_GROWTH_TOL:
        return {"verdict": "grows", "growth": growth,
                "rescued_by_any_of": [], "rescued_by_all_of": []}

    # Everything this model could eat that the diet does not offer -- minus anything
    # we refuse to propose (see FORBIDDEN_RESCUES).
    open_now = {k for k, v in supplied.items() if v > 0}
    candidates = sorted(
        ex for ex in exchanges
        if ex not in open_now and ex[len("EX_"):] not in forbid
    )

    with model:
        model.medium = {**base, **{ex: rescue_bound for ex in candidates}}
        rescued = float(model.slim_optimize() or 0.0)
    if rescued <= NO_GROWTH_TOL:
        return {
            "verdict": "model_cannot_grow",
            "growth": growth,
            "rescued_by_any_of": [],
            "rescued_by_all_of": [],
            "n_candidates": len(candidates),
        }

    # Every metabolite that ALONE restores growth.  Order-independent and complete.
    singles = []
    for ex in candidates:
        with model:
            model.medium = {**base, ex: rescue_bound}
            if float(model.slim_optimize() or 0.0) > NO_GROWTH_TOL:
                singles.append(ex)

    report = {
        "verdict": "medium_gap",
        "growth": growth,
        "growth_if_rescued": rescued,
        "rescued_by_any_of": singles,
        "rescued_by_all_of": [],
        "n_candidates": len(candidates),
    }
    if singles:
        return report

    # No single metabolite is enough: the model needs a COMBINATION.  Now -- and only
    # now -- shrink the full set, because there is no equivalence class to enumerate.
    needed = list(candidates)
    for ex in candidates:
        trial = [e for e in needed if e != ex]
        with model:
            model.medium = {**base, **{e: rescue_bound for e in trial}}
            if float(model.slim_optimize() or 0.0) > NO_GROWTH_TOL:
                needed = trial
    report["rescued_by_all_of"] = needed
    report["combination_is_one_of_several"] = True
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
