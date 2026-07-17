#!/usr/bin/env python3
"""Relabel a curated fungal GEM's EXCHANGE metabolites into the BiGG namespace, so it
shares one metabolite pool with the CarveMe (BiGG) bacteria and can actually cross-feed.

Why this exists
---------------
muODE's community is a SHARED extracellular pool, and members touch it by matching
exchange-metabolite IDS.  The multikingdom bacteria (B. thetaiotaomicron, K. pneumoniae)
reconstruct with CarveMe -> BiGG (`glc__D_e`, `o2_e`, ...).  A curated fungal GEM speaks
its own dialect -- Yeast8 (the consensus S. cerevisiae model, SysBioChalmers/yeast-GEM)
names glucose `s_0565`, oxygen `s_1277` -- so out of the box the fungus and the bacteria
would sit in the same simulation and never touch the same molecule.  Its point in this
example (a facultative fungus SCAVENGING mucosal O2 to protect the obligate anaerobe)
requires them to share `o2_e`, so the namespaces must be reconciled.

The reconciliation is cheap and annotation-driven, not hand-typed: Yeast8 already carries
`bigg.metabolite` annotations on ~87% of its exchanges, so for each exchange we read that
annotation and rename the extracellular metabolite to `<bigg>_e`.  This is the same
annotation-reading move `muode.media` uses to translate a diet across namespaces -- here
applied to a model's exchange boundary rather than a diet.  Only the metabolite IDS
change; stoichiometry and bounds are untouched, so growth is identical (verified).

Exchanges with no `bigg.metabolite` annotation (Yeast8: 34, all yeast-specific aroma
esters, unusual nucleotides, tyrosol, ...) are left under their native ids: they are not
gut cross-feeding metabolites, so they simply do not participate in the shared pool,
which is correct.

    bash  examples/multikingdom/fetch_yeast8.sh          # -> yeast-GEM.xml
    python examples/multikingdom/harmonize_fungal_gem.py \
        --in yeast-GEM.xml --out S_cerevisiae_bigg.xml.gz

Point `eukaryote_models:` in config.yaml at the output; the pipeline then drops it into
the community in place of an automated reconstruction (see gems/_euk_engine in the
Snakefile).  The model is a downloaded/derived artifact (~12 MB), so it is NOT committed
-- fetch + harmonize locally, like the bacterial genomes.
"""
from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

#: Metabolites that MUST end up BiGG-named or the community cannot cross-feed on them.
#: Not exhaustive -- a tripwire for a wrong/annotation-poor input model.
_SHARED = ("glc__D_e", "o2_e", "co2_e", "ac_e", "etoh_e", "nh4_e", "pi_e", "h2o_e")


def harmonize(model):
    """Rename each exchange metabolite to ``<bigg>_e`` from its bigg.metabolite annotation.

    Returns ``(renamed, skipped, collisions)``.  Mutates ``model`` in place; only ids
    change, so the LP (and thus growth) is unchanged.
    """
    renamed = skipped = collisions = 0
    ids = {m.id for m in model.metabolites}
    for ex in list(model.exchanges):
        if len(ex.metabolites) != 1:
            continue
        met = next(iter(ex.metabolites))
        bigg = (met.annotation or {}).get("bigg.metabolite")
        if not bigg:
            skipped += 1
            continue
        bigg = bigg[0] if isinstance(bigg, list) else bigg
        new = f"{bigg}_e"
        if new == met.id:
            renamed += 1
            continue
        if new in ids:                       # two exchanges claim one bigg id: keep the
            collisions += 1                  # first, leave this one native (report it)
            continue
        ids.discard(met.id)
        ids.add(new)
        met.id = new
        renamed += 1
    model.repair()
    return renamed, skipped, collisions


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in", dest="inp", type=Path, required=True,
                    help="curated fungal GEM (SBML) with bigg.metabolite annotations")
    ap.add_argument("--out", type=Path, required=True,
                    help="BiGG-exchange-harmonized GEM (.xml or .xml.gz)")
    args = ap.parse_args()

    import cobra
    import logging
    logging.disable(logging.CRITICAL)

    model = cobra.io.read_sbml_model(str(args.inp))
    before = model.slim_optimize() or 0.0

    renamed, skipped, collisions = harmonize(model)
    print(f"renamed {renamed} exchange metabolites to <bigg>_e; "
          f"skipped {skipped} (no bigg annotation); collisions {collisions}", file=sys.stderr)

    # relabeling must not touch the LP
    after = model.slim_optimize() or 0.0
    assert abs(after - before) < 1e-6, (
        f"growth changed on relabel ({before:.4f} -> {after:.4f}); a rename hit a "
        "reaction it should not have")
    print(f"growth preserved: {before:.4f}/h -> {after:.4f}/h", file=sys.stderr)

    exch = {next(iter(ex.metabolites)).id for ex in model.exchanges if len(ex.metabolites) == 1}
    missing = [m for m in _SHARED if m not in exch]
    print("shared-pool metabolites now BiGG-named:", file=sys.stderr)
    for m in _SHARED:
        print(f"  {'OK ' if m in exch else '!! '} {m}", file=sys.stderr)
    if missing:
        print(f"WARNING: {missing} did not resolve -- the input model may lack bigg "
              "annotations on them; the fungus cannot cross-feed these.", file=sys.stderr)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    if args.out.suffix == ".gz":
        tmp = args.out.with_suffix("")
        cobra.io.write_sbml_model(model, str(tmp))
        with open(tmp, "rb") as fh, gzip.open(args.out, "wb") as gz:
            gz.writelines(fh)
        tmp.unlink()
    else:
        cobra.io.write_sbml_model(model, str(args.out))
    print(f"wrote {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
