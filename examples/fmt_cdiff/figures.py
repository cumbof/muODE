#!/usr/bin/env python
"""Figures for the rCDI/FMT mechanism-decomposition study.

Reads what ``run.py`` wrote (``summary.json`` + ``<arm>_biomass.csv``) and draws the
two figures that carry the study's claim:

  1. ``mechanism_decomposition.png`` -- pathogen final biomass per arm.  THE result:
     if ``fmt_competition`` (bile and pH both ablated) clears the pathogen just as
     well as ``fmt_full``, then clearance is nutrient competition and the advertised
     bile mechanism is decoration.  The figure is built to make that readable at a
     glance rather than to flatter the hypothesis.
  2. ``pathogen_trajectories.png`` -- C. difficile biomass over time, one line per
     arm, so *when* the arms diverge is visible (before or after the transplant).

This plots measured output only: every number comes from the CSVs, and an arm that
was not run is simply absent.  Nothing here fills in, smooths, or extrapolates.

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
    ("no_fmt", "no FMT (control)", "#2a78d6"),
    ("fmt_full", "FMT, full mechanism", "#008300"),
    ("fmt_no_bile", "FMT, bile ablated", "#e87ba4"),
    ("fmt_no_ph", "FMT, pH ablated", "#eda100"),
    ("fmt_competition", "FMT, competition only", "#1baf7a"),
]

#: run.py's clearance rule, restated here so the figure's line and the summary's
#: verdict cannot drift apart.
CLEARED_BELOW = 0.1

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_MUTED = "#52514e"


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


def plot_decomposition(summary: dict, out: Path) -> Path:
    """Pathogen final biomass per arm -- the headline.

    One measure across five named conditions, so: one hue, no legend, and the
    verdict spelled out per bar.  Colour never carries the CLEARED/PERSISTS
    distinction on its own -- the text does.
    """
    arms = [(k, label, c) for k, label, c in ARMS if k in summary["arms"]]
    if not arms:
        raise SystemExit("summary.json contains no known arms")

    vals = [summary["arms"][k]["pathogen_final"] for k, _, _ in arms]
    labels = [label for _, label, _ in arms]
    verdicts = ["CLEARED" if summary["arms"][k]["cleared"] else "PERSISTS"
                for k, _, _ in arms]

    fig, ax = plt.subplots(figsize=(8.2, 3.6), facecolor=SURFACE)
    y = range(len(arms))
    ax.barh(list(y), vals, height=0.62, color=[c for _, _, c in arms],
            edgecolor=SURFACE, linewidth=2)

    ax.axvline(CLEARED_BELOW, color=INK_MUTED, linestyle="--", linewidth=1.2, zorder=0)
    # Blended transform: x in data units (pinned to the line), y in axes units (pinned
    # to the top), so the caption cannot fall off the axis as the bar count changes.
    ax.text(CLEARED_BELOW, 1.01, f"clearance threshold ({CLEARED_BELOW})",
            transform=ax.get_xaxis_transform(), color=INK_MUTED, fontsize=8,
            ha="left", va="bottom")

    span = max(vals) if max(vals) > 0 else 1.0
    for i, (v, verdict) in enumerate(zip(vals, verdicts)):
        # Start the label past the threshold line for bars shorter than it, so the
        # text never sits on top of the rule it is being compared against.
        ax.annotate(f"{v:.3f}  {verdict}", xy=(max(v, CLEARED_BELOW), i), xytext=(6, 0),
                    textcoords="offset points", va="center", fontsize=9, color=INK)

    ax.set_yticks(list(y))
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlim(0, span * 1.32)
    ax.set_xlabel("C. difficile biomass at t_end (gDW/L)", color=INK_MUTED, fontsize=9)
    ax.set_title("Which arm of colonization resistance clears the pathogen?",
                 color=INK, fontsize=12, loc="left", pad=12)
    ax.grid(axis="x", color="#ecebe6", linewidth=1)
    ax.set_axisbelow(True)
    _style(ax)

    fig.tight_layout()
    path = out / "mechanism_decomposition.png"
    fig.savefig(path, dpi=200, facecolor=SURFACE)
    plt.close(fig)
    return path


def plot_trajectories(results: Path, summary: dict, out: Path) -> Path | None:
    """Pathogen biomass over time, one line per arm.

    Shows *when* the arms separate, which the bar chart cannot: an arm that only
    diverges after the transplant is evidence the transplant did it.
    """
    fig, ax = plt.subplots(figsize=(8.2, 4.2), facecolor=SURFACE)
    drawn = 0
    for key, label, colour in ARMS:
        csv = results / f"{key}_biomass.csv"
        if not csv.exists():
            continue
        df = pd.read_csv(csv, index_col=0)
        if gs.PATHOGEN not in df.columns:
            continue
        ax.plot(df.index, df[gs.PATHOGEN], color=colour, linewidth=2, label=label)
        drawn += 1

    if not drawn:
        plt.close(fig)
        return None

    ax.axhline(CLEARED_BELOW, color=INK_MUTED, linestyle="--", linewidth=1.2, zorder=0)
    ax.text(0.995, CLEARED_BELOW, f"clearance threshold ({CLEARED_BELOW})",
            transform=ax.get_yaxis_transform(), color=INK_MUTED, fontsize=8,
            ha="right", va="bottom")
    ax.set_xlabel("time (h)", color=INK_MUTED, fontsize=9)
    ax.set_ylabel("C. difficile biomass (gDW/L)", color=INK_MUTED, fontsize=9)
    ax.set_title("C. difficile over time, by arm", color=INK, fontsize=12,
                 loc="left", pad=12)
    ax.grid(color="#ecebe6", linewidth=1)
    ax.set_axisbelow(True)
    # Identity rides the legend, not the line ends: the arms that work all converge on
    # ~0 and their end-labels collide.  The per-arm CSVs are the table view, which is
    # what the palette's sub-3:1 slots oblige us to provide.
    ax.legend(frameon=False, fontsize=8, labelcolor=INK_MUTED, loc="upper left")
    ax.margins(x=0.02)
    _style(ax)

    fig.tight_layout()
    path = out / "pathogen_trajectories.png"
    fig.savefig(path, dpi=200, facecolor=SURFACE)
    plt.close(fig)
    return path


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", default="results/fmt",
                    help="directory run.py wrote (default: results/fmt)")
    ap.add_argument("--outdir", default=str(HERE / "figures"),
                    help="where to write the PNGs (default: the example's figures/)")
    args = ap.parse_args()

    results = Path(args.results)
    summary_path = results / "summary.json"
    if not summary_path.exists():
        print(f"no {summary_path} -- run examples/fmt_cdiff/run.py first", file=sys.stderr)
        return 1

    summary = json.loads(summary_path.read_text())
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    for p in (plot_decomposition(summary, out), plot_trajectories(results, summary, out)):
        if p is not None:
            print(f"wrote {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
