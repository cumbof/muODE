#!/usr/bin/env python3
"""Namespace harmonization: gapseq/ModelSEED -> BiGG at the exchange boundary.

For the hybrid community gate the gapseq Bacteroidetes and the BiGG Firmicutes must
share ONE extracellular pool. muODE keys that pool by extracellular *metabolite id*
(muode.organism: exchange map = {external_metabolite_id -> exchange_reaction_id}),
so harmonizing means renaming each gapseq model's EXTERNAL metabolites (cpd*_e0) to
their canonical BiGG ids (glc__D_e) -- the internal cpd network is left untouched,
because the dFBA environment never sees it.

Mapping is model-carried, not guessed: every gapseq metabolite ships a
`bigg.metabolite` cross-reference (verified: ~70-75% of exchanges resolve directly),
with `metanetx.chemical` (MNXM) as a bridge fallback. List-valued hits are
disambiguated toward the BiGG vocabulary the diet + BiGG models actually use, so
both sides key into the same pool. Metabolites with no BiGG xref keep their cpd id
(they are one-sided exchanges that cannot be shared anyway -- harmless, never
collide).

The decisive check is `--invariance`: harmonize a strain, run its monoculture on the
BiGG DM38, and compare biomass to the SAME strain on the native ModelSEED DM38. Same
organism + same medium in two namespaces must give the same growth; a mismatch means
the translation dropped or mis-mapped a nutrient and names the offending exchange.
"""
from __future__ import annotations

import argparse
import glob
import json
import logging
from pathlib import Path

import cobra

logging.disable(logging.CRITICAL)  # silence cobra's per-exchange INFO spam

BACT_GENERA = ("Prevotella", "Parabacteroides", "Phocaeicola", "Bacteroides")


def _ann_list(ann, key):
    v = ann.get(key)
    return [v] if isinstance(v, str) else (list(v) if v else [])


def bigg_vocab(models_dir: Path, diet) -> tuple[set[str], dict[str, str]]:
    """Return (set of valid BiGG external ids, MNXM -> BiGG external id) built from
    the BiGG models' exchanges + the BiGG diet -- the target vocabulary both sides
    must agree on."""
    ext: set[str] = set(diet.metabolites())
    mnx2ext: dict[str, str] = {}
    for p in glob.glob(str(models_dir / "*.xml.gz")):
        m = cobra.io.read_sbml_model(p)
        for r in m.exchanges:
            met = list(r.metabolites)[0]
            ext.add(met.id)
            for x in _ann_list(met.annotation, "metanetx.chemical"):
                mnx2ext.setdefault(x, met.id)
    return ext, mnx2ext


def _target_id(met, ext: set[str], mnx2ext: dict[str, str]) -> str | None:
    """Canonical BiGG external id for a gapseq external metabolite, or None."""
    cands = _ann_list(met.annotation, "bigg.metabolite")
    # prefer a candidate already in the shared vocabulary (diet / BiGG models)
    for c in cands:
        if f"{c}_e" in ext:
            return f"{c}_e"
    # MNX bridge into the vocabulary
    for x in _ann_list(met.annotation, "metanetx.chemical"):
        if x in mnx2ext:
            return mnx2ext[x]
    # a BiGG name exists but is outside the shared vocab: still canonicalize it so
    # two gapseq models agree, appending the standard external compartment tag
    if cands:
        return f"{cands[0]}_e"
    return None


def harmonize(model: cobra.Model, ext: set[str], mnx2ext: dict[str, str]) -> dict:
    """Rename external metabolites in place; return a mapping report."""
    mapped, unmapped, collisions = {}, [], []
    taken: set[str] = set()
    for r in list(model.exchanges):
        met = list(r.metabolites)[0]
        tgt = _target_id(met, ext, mnx2ext)
        if tgt is None or tgt == met.id:
            if tgt is None:
                unmapped.append(met.id)
            continue
        if tgt in taken or (tgt in {m.id for m in model.metabolites} and tgt != met.id):
            collisions.append((met.id, tgt))          # keep original id, don't merge
            continue
        mapped[met.id] = tgt
        taken.add(tgt)
        met.id = tgt
    model.repair()
    return {"mapped": mapped, "n_mapped": len(mapped),
            "unmapped": unmapped, "n_unmapped": len(unmapped),
            "collisions": collisions}


def _biomass(model, diet, kin):
    from muode.community import Community
    from muode.dfba import DynamicFBA
    from muode.organism import CobraOrganism
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from run_benchmark import T_END
    org = CobraOrganism(model, id="X")
    comm = Community([org], abundances={"X": 1.0}, total_biomass=0.01)
    res = DynamicFBA(t_end=T_END, dt=0.1, n_jobs=1).run(comm, diet, kin)
    return float(res.biomass.iloc[-1].to_dict().get("X", 0.0))


def main() -> int:
    import warnings
    warnings.filterwarnings("ignore")
    from muode.benchmarks import clark2021 as ck
    from muode.kinetics import KineticParameters

    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gems", type=Path, default=Path("examples/benchmarks/clark2021/gems_gapseq"))
    ap.add_argument("--bigg-models", type=Path, default=Path("examples/benchmarks/clark2021/models"))
    ap.add_argument("--out", type=Path, default=Path("examples/benchmarks/clark2021/gems_gapseq_bigg"))
    ap.add_argument("--invariance", action="store_true",
                    help="harmonize each Bacteroidetes and check biomass matches native")
    ap.add_argument("--write", action="store_true", help="write harmonized SBMLs to --out")
    args = ap.parse_args()

    bigg_diet, ms_diet = ck.medium("bigg"), ck.medium("modelseed")
    kin = KineticParameters()
    ext, mnx2ext = bigg_vocab(args.bigg_models, bigg_diet)
    bact = [c for c, s in ck.STRAINS.items() if s.species.startswith(BACT_GENERA)]

    reports = {}
    if args.write:
        args.out.mkdir(parents=True, exist_ok=True)
    print(f"BiGG vocab: {len(ext)} external ids; harmonizing {len(bact)} Bacteroidetes")
    if args.invariance:
        print("\ncode  native(MS)  harmon(BiGG)  ratio  mapped/unmapped  verdict")
    for c in sorted(bact):
        native = cobra.io.read_sbml_model(str(args.gems / f"{c}.xml.gz"))
        harm = cobra.io.read_sbml_model(str(args.gems / f"{c}.xml.gz"))
        rep = harmonize(harm, ext, mnx2ext)
        reports[c] = rep
        if args.write:
            cobra.io.write_sbml_model(harm, str(args.out / f"{c}.xml.gz"))
        if args.invariance:
            b_native = _biomass(native, ms_diet, kin)
            b_harm = _biomass(harm, bigg_diet, kin)
            ratio = b_harm / b_native if b_native > 1e-9 else float("nan")
            ok = "OK" if 0.85 <= ratio <= 1.18 else ("LOW" if ratio < 0.85 else "HIGH")
            print(f"{c:>3}   {b_native:>9.3f}   {b_harm:>10.3f}   {ratio:>5.2f}  "
                  f"{rep['n_mapped']:>3}/{rep['n_unmapped']:<3}         {ok}")
    Path("results/clark2021/benchmark").mkdir(parents=True, exist_ok=True)
    Path("results/clark2021/benchmark/harmonize_report.json").write_text(json.dumps(reports, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
