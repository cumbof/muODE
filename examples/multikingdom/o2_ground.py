#!/usr/bin/env python3
"""Ground the O2-protection arm: find the mucosal-O2-influx regime where a SINGLE
facultative scavenger cannot hold the lumen anaerobic, so removing the fungus lets O2
rise DURING the anaerobe's growth window (not just after it, as at influx=0.15).

ki_o2 is fixed at the strict-anaerobe anchor (0.005 mmol/L ~ 5 uM); we sweep o2_influx
(sustained mucosal delivery) and, for each, run the full and no-fungus arms and report:
  * O2 at t=12 h (mid growth-window) for each arm -- does removing the fungus keep O2 up?
  * final Bacteroides for each arm and the protection ratio full/no_fungus.
One extra row uses consume=False (O2 uptake left entirely to the real GEMs' FBA), the
docstring-prescribed genome-scale mode, as a cross-check. Prints a table; no artifact.
"""
import warnings, logging
warnings.filterwarnings("ignore"); logging.disable(logging.CRITICAL)
import genome_scale as gs


def o2_at(res, t=12.0):
    s = res.environment["oxygen"]; idx = s.index
    tt = min(idx, key=lambda x: abs(float(x) - t))
    return float(s.loc[tt])


def bact(res):
    return float(res.biomass[gs.BACTEROIDES].iloc[-1])


def main():
    print(f"{'influx':>7s} {'consume':>7s} | {'O2@12 full':>10s} {'O2@12 noF':>10s} | "
          f"{'Bact full':>9s} {'Bact noF':>9s} {'ratio':>6s}", flush=True)
    print("-" * 72, flush=True)
    grid = [(0.15, True), (0.5, True), (1.0, True), (2.0, True), (4.0, True), (1.0, False)]
    for influx, consume in grid:
        kw = dict(ki_o2=0.005, o2_influx=influx, consume=consume)
        rf = gs.run(True, True, **kw)
        rn = gs.run(False, True, **kw)
        bf, bn = bact(rf), bact(rn)
        ratio = bf / bn if bn > 1e-6 else float("inf")
        print(f"{influx:7.2f} {str(consume):>7s} | {o2_at(rf):10.3f} {o2_at(rn):10.3f} | "
              f"{bf:9.3f} {bn:9.3f} {ratio:6.2f}", flush=True)


if __name__ == "__main__":
    main()
