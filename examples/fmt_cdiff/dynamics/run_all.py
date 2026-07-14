#!/usr/bin/env python
"""
run_all.py -- Regenerate every figure in the muODE rCDI reproducibility package.

Runs, in order:
  * antibiotic_relapse.py    (antibiotic collapse + spore-driven relapse)
  * fmt_resolution.py        (FMT engraftment + cure)
  * designed_consortium.py   (naive vs. designed 12-member consortium)
  * failed_fmt_autopsy.py    (failed-FMT bottleneck autopsy)

Every panel is produced by executing the real muODE dynamic-FBA engine and its
ecology layer. The model (guilds, medium, ecology, scenario runner) lives in
designed.py and is re-exported through common.py, which also holds the
plotting style and the honest-scope note. All outputs land in ./results/.

This is deterministic: the muODE engine contains no stochasticity, so a rerun
reproduces byte-identical figures and CSVs.

Run:  conda run -n muode python run_all.py
Total runtime: on the order of 15 minutes (seven full 42-day dynamic-FBA runs:
one each for the relapse and FMT arms, two for the consortium, three for the
autopsy).
"""

from __future__ import annotations

import time
import warnings

import antibiotic_relapse as relapse
import fmt_resolution as fmt
import designed_consortium as consortium
import failed_fmt_autopsy as autopsy


def main():
    warnings.filterwarnings("ignore")
    steps = [
        ("Antibiotic collapse & spore-driven relapse", relapse.main),
        ("FMT engraftment & cure", fmt.main),
        ("Naive vs. designed synthetic consortium", consortium.main),
        ("Failed-FMT autopsy -- two bottlenecks", autopsy.main),
    ]
    t0 = time.time()
    for title, fn in steps:
        print("\n" + "=" * 74)
        print(title)
        print("=" * 74)
        fn()
    print("\n" + "=" * 74)
    print(f"All figures regenerated in {time.time() - t0:.0f}s. See ./results/.")
    print("=" * 74)


if __name__ == "__main__":
    main()
