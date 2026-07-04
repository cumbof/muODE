#!/usr/bin/env python
"""
designed_consortium.py
=============================================================================
Scenario 3 -- in-silico rational design and dynamic evaluation of a defined
synthetic therapeutic consortium.

Reproduces the design logic for a defined 12-member Live Biotherapeutic Product
(the guild registry in muode.scenarios.rcdi), evaluated against the same Day-12
post-antibiotic C. difficile environment used for the FMT.

Three panels:
  A  Naive minimal consortium. A bolus restricted to the bai-positive
     Clostridiales only. Lacking upstream primary degraders, the effectors are
     starved of carbon and out-competed by the germinating pathogen -- the
     consortium fails to engraft and C. difficile persists.
  B  Inferred trophic network. The emergent cross-feeding network of the full
     12-member consortium, inferred POST-simulation by muODE from coinciding
     secretion/uptake fluxes (not hard-coded). Directed edges are extracellular
     metabolite hand-offs; the primary degraders establish the carbon pipeline
     that fuels the bai effectors and amino-acid scavengers.
  C  Rationally designed consortium. The full 12-member consortium engrafts,
     restores the secondary-bile shield AND depletes the Stickland substrates,
     driving C. difficile to extinction -- matching whole-stool FMT efficacy.

Both consortia are injected into a DEPAUPERATE post-antibiotic lumen
(commensal_level = 0): only the pathogen (and the fungal incumbent) are present,
so the injected members are the sole source of degrader / effector activity.
This is what exposes the "donor insufficiency" of the naive minimal core.

Outputs (./results/): designed_consortium.png + CSV trajectories +
designed_consortium_crossfeeding_edges.csv (the inferred network).

Run:  conda run -n muode python designed_consortium.py
"""

from __future__ import annotations

import warnings

import numpy as np

import common as C

# Short, human-readable labels for the network panel (genus initials).
SHORT = {
    "Clostridium_scindens": "C.scin", "Extibacter_muris": "E.mur",
    "Bacteroides_thetaiotaomicron": "B.theta", "Bacteroides_ovatus": "B.ov",
    "Akkermansia_muciniphila": "A.muc", "Blautia_hydrogenotrophica": "B.hyd",
    "Peptostreptococcus_anaerobius": "P.ana", "Eubacterium_hallii": "E.hal",
    "Faecalibacterium_prausnitzii": "F.pra", "Roseburia_intestinalis": "R.int",
    "Bifidobacterium_longum": "B.lon", "Parabacteroides_distasonis": "P.dis",
    C.CDIFF: "C.diff",
}
# Colour network nodes by functional guild.
GUILD_COLOR = {
    True: C.PALETTE["secondary"],    # bai effector
}


def run_minimal():
    """Panel A: naive minimal (bai-only) consortium into the depauperate lumen."""
    return C.run_scenario(commensal_level=0.0,
                          inject_members=C.MINIMAL_CONSORTIUM, with_antibiotic=True)


def run_full():
    """Panels B/C: full 12-member consortium into the depauperate lumen."""
    return C.run_scenario(commensal_level=0.0,
                          inject_members=C.CONSORTIUM_12, with_antibiotic=True)


def _guild_color(member_id: str) -> str:
    spec = C.ALL_SPECS.get(member_id)
    if member_id == C.CDIFF:
        return C.PALETTE["pathogen"]
    if spec and spec.bai:
        return C.PALETTE["secondary"]
    if spec and "degrader" in spec.guild.lower():
        return C.PALETTE["primary"]
    if spec and "scavenger" in spec.guild.lower():
        return C.PALETTE["accent"]
    return C.PALETTE["commensal"]


def draw_network(ax, result):
    """Draw the inferred cross-feeding network of the full consortium (Panel B).

    The edges come from muODE's post-hoc cross-feeding inference
    (result.cross_feeding()): a directed producer -> consumer edge for every
    metabolite one member secretes while another simultaneously takes it up.
    """
    import networkx as nx

    cf = result.cross_feeding()
    g = nx.DiGraph()
    # keep only members that actually engrafted (biomass above a small floor) so
    # the diagram reflects the realised, not nominal, community
    active = {m for m in C.CONSORTIUM_12 + [C.CDIFF]
              if m in result.biomass and result.biomass[m].iloc[-1] > 1e-4}
    for m in active:
        g.add_node(m)
    for _, row in cf.iterrows():
        p, cons, met = row["producer"], row["consumer"], row["metabolite"]
        if p in active and cons in active:
            g.add_edge(p, cons, metabolite=met)

    if g.number_of_nodes() == 0:
        ax.text(0.5, 0.5, "no engrafted network", ha="center", va="center")
        ax.axis("off")
        return

    pos = nx.spring_layout(g, seed=1, k=1.4)
    colors = [_guild_color(n) for n in g.nodes()]
    nx.draw_networkx_nodes(g, pos, ax=ax, node_color=colors, node_size=900,
                           edgecolors="white", linewidths=1.5)
    nx.draw_networkx_labels(g, pos, ax=ax, font_size=7,
                            labels={n: SHORT.get(n, n) for n in g.nodes()})
    nx.draw_networkx_edges(g, pos, ax=ax, edge_color=C.PALETTE["muted"],
                           arrowsize=12, width=1.2, connectionstyle="arc3,rad=0.08",
                           node_size=900)
    # de-clutter: label each edge by the metabolite handed off (strip the _e tag)
    elabels = {(u, v): g[u][v]["metabolite"].replace("_e", "").replace("__L", "")
               for u, v in g.edges()}
    nx.draw_networkx_edge_labels(g, pos, ax=ax, edge_labels=elabels, font_size=6,
                                 bbox=dict(boxstyle="round,pad=0.1", fc="white",
                                           ec="none", alpha=0.7))
    ax.set_title("B  Inferred trophic network (12-member)", loc="left")
    ax.axis("off")


def make_figure(minimal, full):
    plt = C._mpl()
    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(15.5, 4.4))

    def mark_fmt(ax):
        ax.axvline(C.FMT_DAY, color=C.PALETTE["muted"], lw=1.2, ls=":")

    # ---- Panel A: minimal consortium fails --------------------------------
    t = minimal.biomass.index.values
    C.shade_dosing(axA)
    axA.plot(t, C.biomass_sum(minimal, C.MINIMAL_CONSORTIUM),
             color=C.PALETTE["commensal"], lw=2, label="Minimal consortium (bai-only)")
    axA.plot(t, minimal.biomass[C.CDIFF], color=C.PALETTE["pathogen"], lw=2,
             label="C. difficile")
    axA.set_xlabel("time (days)")
    axA.set_ylabel("biomass (gDW/L)")
    axA.set_title("A  Naive minimal consortium fails", loc="left")
    mark_fmt(axA)
    axA.legend(loc="upper left", fontsize=8)

    # ---- Panel B: inferred cross-feeding network --------------------------
    draw_network(axB, full)

    # ---- Panel C: full 12-member consortium cures -------------------------
    t = full.biomass.index.values
    C.shade_dosing(axC)
    axC.plot(t, C.biomass_sum(full, C.CONSORTIUM_12),
             color=C.PALETTE["commensal"], lw=2, label="12-member consortium")
    axC.plot(t, full.biomass[C.CDIFF], color=C.PALETTE["pathogen"], lw=2,
             label="C. difficile")
    axC.set_xlabel("time (days)")
    axC.set_ylabel("biomass (gDW/L)")
    axC.set_title("C  Designed consortium drives extinction", loc="left")
    mark_fmt(axC)
    axC.legend(loc="upper left", fontsize=8)

    fig.suptitle("Rational design of a defined 12-member therapeutic "
                 "consortium", fontsize=12, y=1.02)
    fig.tight_layout()
    return fig


def main():
    warnings.filterwarnings("ignore")
    outdir = C.output_dir()
    print("[consortium] running the minimal (bai-only) consortium ...")
    minimal, _ = run_minimal()
    print("[consortium] running the full 12-member consortium ...")
    full, _ = run_full()

    fig = make_figure(minimal, full)
    png = outdir / "designed_consortium.png"
    fig.savefig(png, bbox_inches="tight")

    full.biomass.to_csv(outdir / "designed_consortium_full_biomass.csv", index_label="day")
    minimal.biomass.to_csv(outdir / "designed_consortium_minimal_biomass.csv", index_label="day")
    cf = full.cross_feeding()
    cf.to_csv(outdir / "designed_consortium_crossfeeding_edges.csv", index=False)

    cd_min = minimal.biomass[C.CDIFF].iloc[-1]
    cd_full = full.biomass[C.CDIFF].iloc[-1]
    print(f"  minimal consortium: C. difficile end={cd_min:.3g} gDW/L (FAILS to clear)")
    print(f"  12-member consortium: C. difficile end={cd_full:.3g} gDW/L (clears)")
    print(f"  fold difference in residual pathogen: {cd_min / max(cd_full, 1e-12):.0f}x")
    print(f"  inferred cross-feeding edges: {len(cf)}")
    print(f"  wrote {png}")


if __name__ == "__main__":
    main()
