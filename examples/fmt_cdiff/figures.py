#!/usr/bin/env python
"""Figures for the rCDI/FMT mechanism-decomposition study.

Reads what ``run.py`` wrote (per-arm ``summary.json`` + ``<arm>_biomass.csv``) plus the
leave-one-out robustness JSONs, and draws the panels that carry the study's claim:

  a. ``mechanism_decomposition`` -- pathogen biomass per arm, split into vegetative and
     spore pools. THE headline, and the split is half of it: a bar that is mostly spore is
     not a cure, it is a pathogen waiting for the drug to wash out.
  b. ``rebound_rate`` -- net growth rate of the pathogen at t_end. The zero line is the
     whole panel: left of it the community excludes the pathogen, right of it the arm only
     delays recurrence.
  c. ``pathogen_trajectories`` -- C. difficile over time, one line per arm, so *when* the
     arms diverge is visible (before or after the transplant).
  d. ``mechanism_specificity`` -- leave-one-donor-out: deoxycholate collapses only when the
     C. scindens keystone is removed, and only then does the pathogen rebound. The
     specificity control and the in-silico recovery of Buffie 2015.

``mechanism_specificity`` and ``robustness_envelope`` read the robustness/ result JSONs;
the others read results/fmt. Every number comes from those files -- nothing here fills in,
smooths, or extrapolates. ``plot_main_figure`` composes a-d into a single main-text figure.

    python examples/fmt_cdiff/figures.py --results results/fmt
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import scenario as gs  # noqa: E402

HERE = Path(__file__).resolve().parent

#: Arms in the order the argument runs: control, treatment, then the ablations that
#: take the treatment apart.  Colours are bound to the ARM, never to its rank -- an
#: arm keeps its hue whether or not the others were run (categorical slots 1-5 of the
#: validated palette, in fixed order).
ARMS = [
    ("untreated", "untreated (no drug, no FMT)", "#2a78d6"),
    ("abx_only", "vancomycin only", "#008300"),
    ("fmt_only", "FMT only (no drug)", "#e87ba4"),
    ("fmt_full", "vancomycin + FMT", "#eda100"),
    ("fmt_no_bile", "+ FMT, bile ablated", "#1baf7a"),
    ("fmt_no_ph", "+ FMT, pH ablated", "#eb6834"),
    ("fmt_competition", "+ FMT, competition only", "#4a3aa7"),
]

#: Imported rather than restated so the figures and the summary cannot drift apart.
from run import CLEARED_BELOW, rebound_rate  # noqa: E402

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_MUTED = "#52514e"

#: The two pools the pathogen can be in; colour bound to the POOL, not the arm.
VEG_COLOUR = "#2a78d6"
SPORE_COLOUR = "#eb6834"

#: Donors dropped one at a time in the leave-one-out specificity panel, and the single
#: 7-alpha-dehydroxylating keystone whose removal abolishes deoxycholate.
SPECIFICITY_DONORS = ["R_intestinalis_L182", "F_prausnitzii_A2165",
                      "B_thetaiotaomicron_VPI5482", "C_scindens_ATCC35704"]
SPECIFICITY_SHORT = {"R_intestinalis_L182": "R. intestinalis",
                     "F_prausnitzii_A2165": "F. prausnitzii",
                     "B_thetaiotaomicron_VPI5482": "B. theta",
                     "C_scindens_ATCC35704": "C. scindens"}
KEYSTONE = "C_scindens_ATCC35704"

#: The three perturbation sweeps behind the operating-envelope panel.
ENVELOPE_SWEEPS = [
    ("colonic washout (1/h)", "dil", [0.0125, 0.02, 0.025, 0.033, 0.05], "#2a78d6"),
    ("FMT dose (gDW/L)", "dose", [0.02, 0.035, 0.05, 0.08], "#008300"),
    ("FMT timing (h)", "time", [6.0, 12.0, 24.0, 36.0], "#4a3aa7"),
]


def _style(ax) -> None:
    """Recessive axes: the data is the ink, the frame is not."""
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("#d8d7d2")
    ax.tick_params(colors=INK_MUTED, labelsize=9, length=3)
    for lbl in ax.get_xticklabels() + ax.get_yticklabels():
        lbl.set_color(INK_MUTED)


def _title(ax, text: str, panel: str | None, size: float = 12) -> None:
    """Left-set title, optionally prefixed with a panel letter for the composite figure."""
    ax.set_title((f"{panel}   " if panel else "") + text,
                 color=INK, fontsize=size, loc="left", pad=12)


# ----------------------------------------------------------------------------------------
# Data loading: run.py writes one summary.json per arm under results/fmt/<arm>/, and the
# CSVs live in the same subdir. Merge them, and resolve either layout.
# ----------------------------------------------------------------------------------------
def load_merged_summary(results: Path) -> dict | None:
    """Merge the per-arm summary.json files into the single {'arms': {...}} shape the
    panels expect. Falls back to a flat results/summary.json if one exists."""
    flat = results / "summary.json"
    merged: dict = {"arms": {}}
    if flat.exists():
        merged.update(json.loads(flat.read_text()))
        merged.setdefault("arms", {})
    for key, _, _ in ARMS:
        p = results / key / "summary.json"
        if p.exists():
            d = json.loads(p.read_text())
            merged["arms"].update(d.get("arms", {}))
            for k in ("t_end", "dt", "dilution_rate"):
                merged.setdefault(k, d.get(k))
    return merged if merged["arms"] else None


def biomass_csv(results: Path, key: str) -> Path | None:
    """Resolve an arm's biomass CSV in either the flat or per-arm-subdir layout."""
    for cand in (results / f"{key}_biomass.csv", results / key / f"{key}_biomass.csv"):
        if cand.exists():
            return cand
    return None


# ----------------------------------------------------------------------------------------
# Panel draw-cores: each draws onto a supplied ax so the standalone PNGs and the composed
# main figure share exactly one code path and cannot drift apart.
# ----------------------------------------------------------------------------------------
def _draw_decomposition(ax, summary: dict, panel: str | None = None) -> bool:
    """Pathogen biomass per arm, split into vegetative and spore -- the headline.

    Stacked rather than totalled because the split IS the finding: a bar that is mostly
    spore is not a cure, it is a pathogen waiting for the drug to wash out.
    """
    arms = [(k, label, c) for k, label, c in ARMS if k in summary["arms"]]
    if not arms:
        return False
    a = [summary["arms"][k] for k, _, _ in arms]
    veg = [x.get("pathogen_final_vegetative", x["pathogen_final"]) for x in a]
    spo = [x.get("pathogen_final_spores", 0.0) for x in a]
    vals = [v + s for v, s in zip(veg, spo)]
    labels = [label for _, label, _ in arms]
    verdicts = [f"({100*s/(v+s):.0f}% spores)" if (v + s) > 0 else ""
                for v, s in zip(veg, spo)]
    y = range(len(arms))
    ax.barh(list(y), veg, height=0.62, color=VEG_COLOUR, edgecolor=SURFACE,
            linewidth=2, label="vegetative (growing)")
    ax.barh(list(y), spo, height=0.62, left=veg, color=SPORE_COLOUR, edgecolor=SURFACE,
            linewidth=2, label="spores (dormant reservoir)")
    ax.legend(frameon=False, fontsize=8, labelcolor=INK_MUTED, loc="lower right")
    ax.axvline(CLEARED_BELOW, color=INK_MUTED, linestyle="--", linewidth=1.2, zorder=0)
    ax.text(CLEARED_BELOW, 1.01, f"clearance threshold ({CLEARED_BELOW})",
            transform=ax.get_xaxis_transform(), color=INK_MUTED, fontsize=8,
            ha="left", va="bottom")
    span = max(vals) if max(vals) > 0 else 1.0
    for i, (v, verdict) in enumerate(zip(vals, verdicts)):
        ax.annotate(f"{v:.3f}  {verdict}", xy=(v, i), xytext=(6, 0),
                    textcoords="offset points", va="center", fontsize=9, color=INK)
    ax.set_yticks(list(y))
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlim(0, span * 1.42)
    ax.set_xlabel("C. difficile biomass at t_end (gDW/L), vegetative + spores",
                  color=INK_MUTED, fontsize=9)
    _title(ax, "How much pathogen is left — and is it awake or dormant?", panel)
    ax.grid(axis="x", color="#ecebe6", linewidth=1)
    ax.set_axisbelow(True)
    _style(ax)
    return True


def _draw_rebound(ax, results: Path, summary: dict, panel: str | None = None) -> bool:
    """Net growth rate at t_end, per arm -- THE outcome, free of the clearance threshold
    and the horizon. The zero line is the whole panel."""
    rows = []
    for key, label, colour in ARMS:
        lam = (summary["arms"].get(key) or {}).get("rebound_rate")
        if lam is None:                       # recompute for runs predating the metric
            csv = biomass_csv(results, key)
            if csv is None:
                continue
            df = pd.read_csv(csv, index_col=0)
            if gs.PATHOGEN not in df.columns:
                continue
            lam = rebound_rate(df[gs.PATHOGEN])
        rows.append((label, float(lam), colour))
    if not rows:
        return False
    y = range(len(rows))
    ax.barh(list(y), [r[1] for r in rows], height=0.62,
            color=[r[2] for r in rows], edgecolor=SURFACE, linewidth=2)
    ax.axvline(0, color=INK, linewidth=1.4, zorder=3)
    for i, (_, lam, _) in enumerate(rows):
        ax.annotate(f"{lam:+.4f}/h", xy=(lam, i),
                    xytext=(6 if lam >= 0 else -6, 0), textcoords="offset points",
                    va="center", ha="left" if lam >= 0 else "right",
                    fontsize=9, color=INK)
    ax.set_yticks(list(y))
    ax.set_yticklabels([r[0] for r in rows])
    ax.invert_yaxis()
    lo = min(0.0, min(r[1] for r in rows))
    hi = max(0.0, max(r[1] for r in rows))
    pad = 0.42 * max(abs(lo), abs(hi), 1e-6)
    ax.set_xlim(lo - pad, hi + pad)
    ax.set_xlabel("net growth rate of C. difficile at t_end (1/h)",
                  color=INK_MUTED, fontsize=9)
    _title(ax, "Is the pathogen excluded, or only delayed?", panel)
    ax.text(0, 1.01, "  <- excluded (washout)   |   growing ->",
            transform=ax.get_xaxis_transform(), color=INK_MUTED, fontsize=8,
            ha="center", va="bottom")
    ax.grid(axis="x", color="#ecebe6", linewidth=1)
    ax.set_axisbelow(True)
    _style(ax)
    return True


def _draw_trajectories(ax, results: Path, panel: str | None = None) -> bool:
    """Pathogen biomass over time, one line per arm -- shows *when* the arms separate."""
    drawn = 0
    for key, label, colour in ARMS:
        csv = biomass_csv(results, key)
        if csv is None:
            continue
        df = pd.read_csv(csv, index_col=0)
        if gs.PATHOGEN not in df.columns:
            continue
        ax.plot(df.index, df[gs.PATHOGEN], color=colour, linewidth=2, label=label)
        drawn += 1
    if not drawn:
        return False
    ax.axhline(CLEARED_BELOW, color=INK_MUTED, linestyle="--", linewidth=1.2, zorder=0)
    ax.text(0.995, CLEARED_BELOW, f"clearance threshold ({CLEARED_BELOW})",
            transform=ax.get_yaxis_transform(), color=INK_MUTED, fontsize=8,
            ha="right", va="bottom")
    ax.set_xlabel("time (h)", color=INK_MUTED, fontsize=9)
    ax.set_ylabel("C. difficile biomass (gDW/L)", color=INK_MUTED, fontsize=9)
    _title(ax, "C. difficile over time, by arm", panel)
    ax.grid(color="#ecebe6", linewidth=1)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=8, labelcolor=INK_MUTED, loc="upper left")
    ax.margins(x=0.02)
    _style(ax)
    return True


def _load_robust(robust_dir: Path):
    res = Path(robust_dir)

    def load(lbl: str) -> dict | None:
        p = res / f"{lbl}.json"
        return json.loads(p.read_text()) if p.exists() else None

    return load


def _draw_specificity(ax, robust_dir: Path, panel: str | None = None) -> bool:
    """Leave-one-donor-out: deoxycholate is entirely attributable to C. scindens, and only
    its removal frees the pathogen. Bars carry DCA; lambda rides the label."""
    load = _load_robust(robust_dir)
    ref = load("ref_bile_on")
    if ref is None:
        return False
    rows = [("all 4 donors", ref, False)]
    for d in SPECIFICITY_DONORS:
        r = load(f"loo_{d}")
        if r is not None:
            rows.append((f"– {SPECIFICITY_SHORT[d]}", r, d == KEYSTONE))
    if len(rows) < 2:
        return False
    dca = [r["dca_final"] for _, r, _ in rows]
    colours = [SPORE_COLOUR if key else VEG_COLOUR for _, _, key in rows]
    y = range(len(rows))
    ax.barh(list(y), dca, height=0.62, color=colours, edgecolor=SURFACE, linewidth=2)
    for i, (_, r, key) in enumerate(rows):
        lam, v = r["lam"], r["dca_final"]
        txt = f"{v:.3f} mM   (λ {lam:+.4f}/h)" + ("  — resistance lost" if key else "")
        ax.annotate(txt, xy=(v, i), xytext=(6, 0), textcoords="offset points",
                    va="center", fontsize=9, color=(SPORE_COLOUR if key else INK))
    ax.set_yticks(list(y))
    ax.set_yticklabels([lbl for lbl, _, _ in rows])
    ax.invert_yaxis()
    span = max(dca) if max(dca) > 0 else 1.0
    ax.set_xlim(0, span * 1.6)
    ax.set_xlabel("secondary bile acid (deoxycholate) at t_end (mM)",
                  color=INK_MUTED, fontsize=9)
    _title(ax, "Only the keystone matters — DCA collapses without C. scindens", panel)
    ax.grid(axis="x", color="#ecebe6", linewidth=1)
    ax.set_axisbelow(True)
    _style(ax)
    return True


# ----------------------------------------------------------------------------------------
# Standalone-PNG wrappers (kept so the figures can also be composed in LaTeX one per file).
# ----------------------------------------------------------------------------------------
def _save(fig, path: Path) -> Path:
    fig.tight_layout()
    fig.savefig(path, dpi=200, facecolor=SURFACE)
    plt.close(fig)
    return path


def plot_decomposition(summary: dict, out: Path) -> Path | None:
    n = sum(1 for k, _, _ in ARMS if k in summary["arms"])
    fig, ax = plt.subplots(figsize=(8.6, 1.4 + 0.52 * n), facecolor=SURFACE)
    if not _draw_decomposition(ax, summary):
        plt.close(fig)
        return None
    return _save(fig, out / "mechanism_decomposition.png")


def plot_rebound(results: Path, summary: dict, out: Path) -> Path | None:
    fig, ax = plt.subplots(figsize=(8.6, 1.4 + 0.52 * len(ARMS)), facecolor=SURFACE)
    if not _draw_rebound(ax, results, summary):
        plt.close(fig)
        return None
    return _save(fig, out / "rebound_rate.png")


def plot_trajectories(results: Path, out: Path) -> Path | None:
    fig, ax = plt.subplots(figsize=(8.2, 4.2), facecolor=SURFACE)
    if not _draw_trajectories(ax, results):
        plt.close(fig)
        return None
    return _save(fig, out / "pathogen_trajectories.png")


def plot_specificity(robust_dir: Path, out: Path) -> Path | None:
    load = _load_robust(robust_dir)
    n = 1 + sum(1 for d in SPECIFICITY_DONORS if load(f"loo_{d}") is not None)
    fig, ax = plt.subplots(figsize=(8.6, 1.4 + 0.52 * n), facecolor=SURFACE)
    if not _draw_specificity(ax, robust_dir):
        plt.close(fig)
        return None
    return _save(fig, out / "mechanism_specificity.png")


def plot_envelope(robust_dir: Path, out: Path) -> Path | None:
    """Operating envelope, small-multiples: one real-valued axis per perturbation so the
    Delta-lambda curves never collide. Supplementary to the specificity panel."""
    load = _load_robust(robust_dir)
    series = []
    for title, name, vals, colour in ENVELOPE_SWEEPS:
        pts = [(v, off["lam"] - on["lam"]) for v in vals
               for on in [load(f"{name}_{v:g}_on")] for off in [load(f"{name}_{v:g}_off")]
               if on and off]
        if pts:
            series.append((title, pts, colour))
    if not series:
        return None
    lo = min(d for _, pts, _ in series for _, d in pts)
    hi = max(d for _, pts, _ in series for _, d in pts)
    pad = 0.15 * max(hi - lo, 1e-3)
    ylim = (lo - pad, hi + pad)
    fig, axes = plt.subplots(len(series), 1, figsize=(6.6, 2.35 * len(series)),
                             facecolor=SURFACE)
    if len(series) == 1:
        axes = [axes]
    for ax, (title, pts, colour) in zip(axes, series):
        xs = [v for v, _ in pts]
        dl = [d for _, d in pts]
        ax.axhspan(0, ylim[1], color="#eef6ee", zorder=0)
        ax.axhline(0, color=INK, linewidth=1.3, zorder=3)
        ax.plot(xs, dl, "o-", color=colour, linewidth=2, ms=6)
        ax.set_ylim(*ylim)
        rng = max(xs) - min(xs) or 1.0
        ax.set_xlim(min(xs) - 0.06 * rng, max(xs) + 0.06 * rng)
        ax.set_xticks(xs)
        ax.set_ylabel("Δλ (1/h)", color=INK_MUTED, fontsize=9)
        ax.set_xlabel(title, color=INK_MUTED, fontsize=9)
        ax.grid(color="#ecebe6", linewidth=1)
        ax.set_axisbelow(True)
        _style(ax)
    _title(axes[0], "The bile-acid benefit holds across the physiological range, fading at extremes\n"
           "(Δλ = λ$_{off}$ − λ$_{on}$; shaded = benefit region, Δλ > 0)", None, size=11.5)
    return _save(fig, out / "robustness_envelope.png")


def plot_main_figure(results: Path, summary: dict, robust_dir: Path, out: Path) -> Path | None:
    """Compose the main-text FMT figure: (a) decomposition, (b) rebound, (c) trajectories,
    (d) leave-one-out specificity -- the added keystone panel -- as a single 2x2 image."""
    fig, axes = plt.subplots(2, 2, figsize=(16, 12), facecolor=SURFACE)
    ok = _draw_decomposition(axes[0, 0], summary, panel="a")
    ok &= _draw_rebound(axes[0, 1], results, summary, panel="b")
    _draw_trajectories(axes[1, 0], results, panel="c")
    ok &= _draw_specificity(axes[1, 1], robust_dir, panel="d")
    if not ok:
        plt.close(fig)
        return None
    fig.tight_layout(w_pad=3.0, h_pad=3.5)
    return _save(fig, out / "fmt_main_figure.png")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", default="results/fmt",
                    help="directory run.py wrote (default: results/fmt)")
    ap.add_argument("--outdir", default=str(HERE / "figures"),
                    help="where to write the PNGs (default: the example's figures/)")
    args = ap.parse_args()

    results = Path(args.results)
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    robust_res = HERE / "robustness" / "results"

    # The robustness panels read only the robustness JSONs, so they render with or without
    # a base run.
    for p in (plot_specificity(robust_res, out), plot_envelope(robust_res, out)):
        if p is not None:
            print(f"wrote {p}")

    summary = load_merged_summary(results)
    if summary is None:
        print(f"no per-arm summaries under {results} -- run examples/fmt_cdiff/run.py "
              "for the arm panels", file=sys.stderr)
        return 0

    for p in (plot_decomposition(summary, out),
              plot_rebound(results, summary, out),
              plot_trajectories(results, out),
              plot_main_figure(results, summary, robust_res, out)):
        if p is not None:
            print(f"wrote {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
