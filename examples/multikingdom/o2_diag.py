#!/usr/bin/env python3
"""Diagnose WHY the genome-scale O2-protection arm is weak before grounding it.

Two hypotheses:
  (A) gating too weak  -- O2 is high during B. theta's growth window but ki_o2 is too
      large to bite;
  (B) timing artifact  -- O2 stays ~0 while B. theta grows (glucose-replete early
      phase), then rises only after growth stops, so gating never bites.

Dumps O2(t) and Bacteroides(t) for the full and no-fungus arms, and reports whether
B. theta's own GEM draws O2 down via FBA (which would keep O2 low without the fungus and
explain (B)). Prints; no artifact.
"""
import warnings, logging
warnings.filterwarnings("ignore"); logging.disable(logging.CRITICAL)
import genome_scale as gs


def traj(res, key, n=9):
    s = res.biomass[key] if key in res.biomass else res.environment[key]
    idx = s.index
    picks = [idx[int(round(i * (len(idx) - 1) / (n - 1)))] for i in range(n)]
    return [(round(float(t), 1), round(float(s.loc[t]), 3)) for t in picks]


def main():
    for label, fungus, phage in [("no_fungus", False, True), ("full", True, True)]:
        res = gs.run(fungus, phage)
        print(f"\n=== {label} ===", flush=True)
        print("  t,O2   :", traj(res, "oxygen"), flush=True)
        print("  t,Bact :", traj(res, gs.BACTEROIDES), flush=True)


if __name__ == "__main__":
    main()
