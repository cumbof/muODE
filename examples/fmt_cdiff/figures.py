#!/usr/bin/env python
"""Figures for the rCDI/FMT mechanism-decomposition study.

Reads what ``run.py`` wrote (``summary.json`` + ``<arm>_biomass.csv``) and draws the
two figures that carry the study's claim:

  1. ``mechanism_decomposition.png`` -- pathogen biomass per arm, split into the
     vegetative and spore pools.  THE result, and the split is half of it: a bar that
     is mostly spore is not a cure, it is a pathogen waiting for the drug to wash out.
     Read ``untreated`` first (did it colonize at all?), then ``fmt_competition``
     against ``fmt_full``: the figure shows the bile/pH arm is load-bearing and
     competition is not (fmt_competition tracks abx_only) -- see the README.
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


#: The two pools the pathogen can be in.  Colour is bound to the POOL here (not to the
#: arm, as in the trajectory figure): the question this figure answers is "how much
#: pathogen is left, and is it awake or dormant?"
VEG_COLOUR = "#2a78d6"
SPORE_COLOUR = "#eb6834"


def plot_decomposition(summary: dict, out: Path) -> Path:
    """Pathogen biomass per arm, split into vegetative and spore -- the headline.

    Stacked rather than totalled because the split IS the finding: a bar that is mostly
    spore is not a cure, it is a pathogen waiting for the drug to wash out. Reading only
    the vegetative pool is what made an earlier run report CLEARED in every arm.
    """
    arms = [(k, label, c) for k, label, c in ARMS if k in summary["arms"]]
    if not arms:
        raise SystemExit("summary.json contains no known arms")

    a = [summary["arms"][k] for k, _, _ in arms]
    # Fall back to the total when an older summary.json has no per-pool split.
    veg = [x.get("pathogen_final_vegetative", x["pathogen_final"]) for x in a]
    spo = [x.get("pathogen_final_spores", 0.0) for x in a]
    vals = [v + s for v, s in zip(veg, spo)]
    labels = [label for _, label, _ in arms]
    # Annotate with the spore FRACTION, not the CLEARED verdict: the verdict is a
    # threshold crossed with a horizon (see run.CLEARED_BELOW), whereas "how much of
    # what is left is dormant" is the thing the reader has to know -- an arm that is
    # 78% spores has not cured anything, it has filled the recurrence reservoir.
    verdicts = [f"({100*s/(v+s):.0f}% spores)" if (v + s) > 0 else ""
                for v, s in zip(veg, spo)]

    # Height tracks the arm count so seven arms are not crammed into a five-arm figure.
    fig, ax = plt.subplots(figsize=(8.6, 1.4 + 0.52 * len(arms)), facecolor=SURFACE)
    y = range(len(arms))
    # 2px surface-coloured gap between the stacked segments.
    ax.barh(list(y), veg, height=0.62, color=VEG_COLOUR, edgecolor=SURFACE,
            linewidth=2, label="vegetative (growing)")
    ax.barh(list(y), spo, height=0.62, left=veg, color=SPORE_COLOUR, edgecolor=SURFACE,
            linewidth=2, label="spores (dormant reservoir)")
    ax.legend(frameon=False, fontsize=8, labelcolor=INK_MUTED, loc="lower right")

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
        ax.annotate(f"{v:.3f}  {verdict}", xy=(v, i), xytext=(6, 0),
                    textcoords="offset points", va="center", fontsize=9, color=INK)

    ax.set_yticks(list(y))
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlim(0, span * 1.42)
    ax.set_xlabel("C. difficile biomass at t_end (gDW/L), vegetative + spores",
                  color=INK_MUTED, fontsize=9)
    ax.set_title("How much pathogen is left — and is it awake or dormant?",
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


def plot_rebound(results: Path, summary: dict, out: Path) -> Path | None:
    """Net growth rate at t_end, per arm -- THE outcome.

    Free of the clearance threshold and of the horizon, both of which decide the
    `cleared` verdict without reference to the biology.  The zero line is the whole
    figure: left of it the community excludes the pathogen, right of it the arm only
    delays recurrence.  An arm can sit at a flattering final biomass and still be on
    the wrong side of this line.
    """
    rows = []
    for key, label, colour in ARMS:
        lam = (summary["arms"].get(key) or {}).get("rebound_rate")
        if lam is None:                       # recompute for runs predating the metric
            csv = results / f"{key}_biomass.csv"
            if not csv.exists():
                continue
            df = pd.read_csv(csv, index_col=0)
            if gs.PATHOGEN not in df.columns:
                continue
            lam = rebound_rate(df[gs.PATHOGEN])
        rows.append((label, float(lam), colour))
    if not rows:
        return None

    fig, ax = plt.subplots(figsize=(8.6, 1.4 + 0.52 * len(rows)), facecolor=SURFACE)
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
    ax.set_title("Is the pathogen excluded, or only delayed?", color=INK, fontsize=12,
                 loc="left", pad=12)
    ax.text(0, 1.01, "  <- excluded (washout)   |   growing ->",
            transform=ax.get_xaxis_transform(), color=INK_MUTED, fontsize=8,
            ha="center", va="bottom")
    ax.grid(axis="x", color="#ecebe6", linewidth=1)
    ax.set_axisbelow(True)
    _style(ax)

    fig.tight_layout()
    path = out / "rebound_rate.png"
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

    for p in (plot_rebound(results, summary, out),
              plot_decomposition(summary, out),
              plot_trajectories(results, summary, out)):
        if p is not None:
            print(f"wrote {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
