"""Visualisation of simulation results.

All matplotlib/networkx imports are local so that importing :mod:`muode` (and
running the engine on a headless cluster) never requires a plotting stack.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Optional

if TYPE_CHECKING:
    from muode.dfba import SimulationResult


def _mpl():
    import matplotlib

    matplotlib.use("Agg")  # headless-safe
    import matplotlib.pyplot as plt

    return plt


def plot_biomass(result: "SimulationResult", path: str | Path, top: Optional[int] = 20) -> None:
    """Stacked-area chart of community composition over time."""
    plt = _mpl()
    bio = result.biomass
    if top is not None and bio.shape[1] > top:
        keep = bio.iloc[-1].sort_values(ascending=False).index[:top]
        bio = bio[keep]
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.stackplot(bio.index, bio.T.values, labels=list(bio.columns))
    ax.set_xlabel("time (h)")
    ax.set_ylabel("biomass (gDW/L)")
    ax.set_title("Community biomass over time")
    ax.legend(loc="upper left", fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_metabolites(result: "SimulationResult", path: str | Path, top: int = 15) -> None:
    """Line plot of the most dynamic extracellular metabolites."""
    plt = _mpl()
    met = result.metabolites
    variability = (met.max() - met.min()).sort_values(ascending=False)
    keep = variability.index[:top]
    fig, ax = plt.subplots(figsize=(9, 5))
    for m in keep:
        ax.plot(met.index, met[m], label=m)
    ax.set_xlabel("time (h)")
    ax.set_ylabel("concentration (mmol/L)")
    ax.set_title("Extracellular metabolites over time")
    ax.legend(loc="best", fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_cross_feeding_network(result: "SimulationResult", path: str | Path) -> None:
    """Directed network of producer -> metabolite -> consumer interactions."""
    import networkx as nx

    plt = _mpl()
    cf = result.cross_feeding()
    if cf.empty:
        return
    g = nx.DiGraph()
    for _, row in cf.iterrows():
        g.add_edge(row["producer"], row["consumer"],
                   metabolite=row["metabolite"], weight=row["strength"])
    fig, ax = plt.subplots(figsize=(8, 7))
    pos = nx.spring_layout(g, seed=0)
    weights = [g[u][v]["weight"] for u, v in g.edges()]
    nx.draw_networkx_nodes(g, pos, ax=ax, node_color="#9ecae1", node_size=1400)
    nx.draw_networkx_labels(g, pos, ax=ax, font_size=8)
    nx.draw_networkx_edges(g, pos, ax=ax, width=[1 + 2 * w / (max(weights) or 1) for w in weights],
                           edge_color="#636363", arrowsize=18, connectionstyle="arc3,rad=0.1")
    nx.draw_networkx_edge_labels(
        g, pos, ax=ax,
        edge_labels={(u, v): g[u][v]["metabolite"] for u, v in g.edges()}, font_size=6,
    )
    ax.set_title("Cross-feeding network")
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def save_all(result: "SimulationResult", outdir: str | Path) -> None:
    """Write the standard figure set to ``outdir``."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    plot_biomass(result, outdir / "biomass.png")
    plot_metabolites(result, outdir / "metabolites.png")
    plot_cross_feeding_network(result, outdir / "cross_feeding.png")
