#!/usr/bin/env python
"""
real_data_adapter.py -- Run the rCDI / FMT ecological-dynamics analysis on YOUR
                        assembled MAGs instead of the bundled toy models.

=============================================================================
WHAT THIS IS
=============================================================================
The scenario scripts in this folder (antibiotic_relapse.py, fmt_resolution.py,
designed_consortium.py, failed_fmt_autopsy.py) run on a *defined synthetic*
community of toy metabolic models (designed.py). That community is
hard-wired: each member's functional role (pathogen / bai effector / degrader /
spore-former / drug-susceptible) is baked into a GuildSpec registry.

This adapter replaces that hard-wired toy community with THREE real-data inputs:

  1. a directory of genome-scale models (GEMs) reconstructed from your MAGs
     (`muode build` + `muode refine`, i.e. CarveMe/gapseq -> gap-filled SBML),
  2. a 2-column abundance TSV per sample (mag_id, rel_abundance) from your
     quantitative profiling, and
  3. a ROLE table (roles.tsv, below) that maps each MAG id to the functional
     flags the ecology layers key on.

It then builds the *identical* ecology stack that designed.py uses
(WeakAcidInhibition + BileAcidTransform + BileAcidInhibition + SporeForming +
Antibiotic, plus an optional PhageInfection), and runs the same DynamicFBA
engine. The layers themselves are organism-agnostic -- they operate on sets of
member ids -- so the only thing you supply is "which MAG is in which set."

=============================================================================
THE THREE ID SPACES MUST AGREE  (this is 90% of the wiring)
=============================================================================
The SAME string id must appear in all three places for a member:
  * the GEM filename stem       ->  <model_dir>/<mag_id>.xml   (Community.from_models
                                    sets organism.id = file stem)
  * the abundance TSV           ->  first column
  * the role table              ->  first column
If they disagree, that member is silently dropped from either the community or a
role set. A consistency check below fails loudly instead.

=============================================================================
METABOLITE NAMESPACE  (read before you trust any cross-feeding)
=============================================================================
Two different kinds of metabolite id appear here:

  * CARBON / AMINO-ACID / SCFA ids (glucose, acetate, butyrate, proline, ...)
    are consumed and secreted by the GEMs THROUGH FBA. Cross-feeding emerges
    only if members share the same exchange-metabolite ids. CarveMe emits BiGG
    ids (e.g. glc__D_e, ac_e, but_e, pro__L_e, gly_e); gapseq/ModelSEED emit
    cpd***** ids. **You MUST set the CONFIG ids below to whatever namespace your
    GEMs actually use** -- otherwise every member is on its own island and no
    cross-feeding (hence no colonisation resistance) can occur.

  * BILE-ACID ids (tca/ca/dca/lca) are NOT GEM reactions. The primary->secondary
    transformation is done by the BileAcidTransform ecology layer on a tracked
    pool keyed to guild membership (see docs/ecology.md). So these ids are
    ecology-internal: leave them as-is unless your GEMs genuinely exchange bile
    acids. The medium just needs to SEED the primary-bile pool.

Run:  conda run -n muode python real_data_adapter.py \
          --models  ./gems \
          --abundance recipient_abundance.tsv \
          --roles   roles.tsv \
          --donor-abundance donor_abundance.tsv \
          --arm fmt
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set

from muode.community import Community
from muode.dfba import DynamicFBA, SimulationResult
from muode.diet import Diet
from muode.ecology import EcologyModel
from muode.inject import Injection
from muode.io_utils import read_abundance
from muode.kinetics import KineticParameters

from muode.ph import WeakAcidInhibition
from muode.bile import BileAcidTransform, BileAcidInhibition
from muode.lifecycle import SporeForming
from muode.antibiotic import Antibiotic
from muode.phage import PhageInfection

# Reuse the *timeline* and the bile-acid pool ids from the published model so the
# real-data run and the toy example share one source of truth for those.
import sys as _sys
from pathlib import Path as _Path

_sys.path.insert(0, str(_Path(__file__).resolve().parents[1]))
from designed import (  # noqa: E402
    T_END, DT, FMT_DAY, DILUTION,
    TCA, CA, DCA, LCA, SECONDARY_BILE,
)


# ===========================================================================
# CONFIG -- metabolite ids in YOUR GEMs' namespace. EDIT THESE.
# ===========================================================================
# Defaults are BiGG (CarveMe). For gapseq/ModelSEED replace with cpd***** ids.
# Any id you list here that a member does not exchange is simply ignored for that
# member, so it is safe to include the full menu.
MUCIN = "mucin_e"          # a host glycan pool, if your degraders exchange one
GLC = "glc__D_e"           # D-glucose  (shared carbon currency after degradation)
PRO = "pro__L_e"           # L-proline  (Stickland substrate -- the pathogen's fuel)
GLY = "gly_e"              # glycine    (Stickland substrate)
AC = "ac_e"                # acetate    (SCFA, acidifying)
BUT = "but_e"              # butyrate   (SCFA, colonocyte fuel)
LAC = "lac__L_e"           # L-lactate  (SCFA)
SUCC = "succ_e"            # succinate

ACIDS = (AC, BUT, LAC)     # weak acids that acidify the lumen (WeakAcidInhibition)


# ===========================================================================
# 1. THE ROLE TABLE  (maps MAG id -> functional flags the ecology keys on)
# ===========================================================================
# roles.tsv is a tab-separated file with a header. First column is mag_id; the
# rest are 0/1 flags. You populate it from functional annotation of the genomes:
#
#   mag_id            pathogen  bai  bsh  spore_former  abx_susceptible  phage_host
#   Cdiff_bin.3       1         0    0    1             1                0
#   Cscindens_bin.7   0         1    0    1             1                0
#   Emuris_bin.4      0         1    1    1             1                0
#   Btheta_bin.1      0         0    0    0             1                0
#   ...
#
# WHERE EACH FLAG COMES FROM (annotation evidence, not the abundance table):
#   pathogen         : taxonomy == C. difficile (the germinating spore-former)
#   bai              : carries the bai operon / 7alpha-dehydroxylase
#                      (baiCD/baiE/baiH etc.) -> makes the secondary-bile shield
#   bsh              : carries a bile-salt hydrolase (deconjugates taurocholate)
#   spore_former     : sporulation machinery present (spo0A, Spo V, ...)
#   abx_susceptible  : expected to be killed by the drug course. For vancomycin,
#                      clinical microbiome studies show a broad community crash,
#                      so most residents = 1; drug-tolerant taxa (many fungi) = 0.
#   phage_host       : susceptible host of a transferred phage (autopsy scenario;
#                      infer from CRISPR spacers / prophage / host-range assays).

@dataclass
class Roles:
    pathogen: Optional[str] = None            # the single pathogen id (or None)
    bai: Set[str] = field(default_factory=set)
    bsh: Set[str] = field(default_factory=set)
    spore_former: Set[str] = field(default_factory=set)
    abx_susceptible: Set[str] = field(default_factory=set)
    phage_host: Set[str] = field(default_factory=set)
    all_ids: Set[str] = field(default_factory=set)


_TRUE = {"1", "true", "yes", "y", "t"}


def read_roles(path: str | Path) -> Roles:
    """Parse roles.tsv into a :class:`Roles` bundle of id sets."""
    roles = Roles()
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        cols = {c.lower(): c for c in (reader.fieldnames or [])}
        if "mag_id" not in cols:
            raise ValueError(f"{path}: first column must be 'mag_id' "
                             f"(got {reader.fieldnames})")

        def flag(row, name):
            col = cols.get(name)
            return bool(col) and str(row[col]).strip().lower() in _TRUE

        for row in reader:
            mid = str(row[cols["mag_id"]]).strip()
            if not mid or mid.startswith("#"):
                continue                        # skip blank / commented rows
            roles.all_ids.add(mid)
            if flag(row, "pathogen"):
                if roles.pathogen and roles.pathogen != mid:
                    raise ValueError(
                        f"{path}: more than one pathogen flagged "
                        f"({roles.pathogen}, {mid}); this model expects one.")
                roles.pathogen = mid
            for name, target in (("bai", roles.bai), ("bsh", roles.bsh),
                                 ("spore_former", roles.spore_former),
                                 ("abx_susceptible", roles.abx_susceptible),
                                 ("phage_host", roles.phage_host)):
                if flag(row, name):
                    target.add(mid)
    if roles.pathogen is None:
        raise ValueError(f"{path}: no member flagged pathogen=1; "
                         "the rCDI ecology needs exactly one pathogen.")
    return roles


# ===========================================================================
# 2. MEDIUM + KINETICS  (over YOUR GEMs' exchange namespace)
# ===========================================================================
# This mirrors designed.py.gut_medium(), but the ids are the CONFIG ids
# above -- so edit CONFIG, not this, to match your GEMs. Concentrations are
# mmol/L; influx is mmol/L/day (open-system renewal). Bile acids (TCA/CA) are
# seeded here even though no GEM consumes them: they feed the ecology bile pool.

def gut_medium() -> Diet:
    return Diet(
        concentrations={
            MUCIN: 8.0, GLC: 1.0,
            PRO: 4.0, GLY: 4.0,
            TCA: 3.0, CA: 3.0,          # primary bile (ecology-tracked germinants)
            DCA: 0.0, LCA: 0.0,         # secondary bile -- made by the bai guild
            AC: 0.0, BUT: 0.0, LAC: 0.0, SUCC: 0.0,
        },
        influx={
            MUCIN: 1.0, PRO: 0.4, GLY: 0.4,
            TCA: 0.4, CA: 0.4,
        },
        name="defined_gut_medium",
    )


def cdi_kinetics() -> KineticParameters:
    """Michaelis-Menten uptake constants (per day) for the shared-pool substrates."""
    vmax, km = 14.0, 0.5
    mets = [MUCIN, GLC, PRO, GLY, AC, BUT, LAC, SUCC]
    return KineticParameters(
        default_vmax=vmax, default_km=km,
        metabolite_defaults={m: (vmax, km) for m in mets},
    )


# ===========================================================================
# 3. THE ECOLOGY STACK  (built from the ROLE table, not a GuildSpec registry)
# ===========================================================================
# Parameter values are copied from designed.py.build_ecology so the
# real-data run behaves like the published example. The ONLY change is that the
# bai / bsh / spore-former / susceptible sets come from your roles table,
# intersected with the members actually present in the community.

def build_ecology(community: Community, roles: Roles, *,
                  with_antibiotic: bool = True,
                  phage_hosts: Optional[Sequence[str]] = None) -> EcologyModel:
    present = set(community.organism_ids)
    pathogen = roles.pathogen
    bai = roles.bai & present
    bsh = roles.bsh & present
    spore_formers = roles.spore_former & present
    susceptible = roles.abx_susceptible & present

    layers: List = [
        WeakAcidInhibition(
            buffer_capacity=60.0, default_ki=20.0, ki={pathogen: 20.0},
            acids=ACIDS,
        ),
        BileAcidTransform(
            bsh_producers=bsh, bai_producers=bai,
            vmax_bsh=25.0, vmax_bai=20.0, tca=TCA, ca=CA, dca=DCA,
        ),
        BileAcidInhibition(targets={pathogen}, ki=1.0, secondary=SECONDARY_BILE),
        SporeForming(
            species=spore_formers,
            initial_spores={pathogen: 0.06},
            germinant=TCA, km_germinant=0.15,
            inhibitor=SECONDARY_BILE, ki_inhibitor=0.03,
            k_germination=0.12, k_sporulation=1.5, mu_stress=0.3,
            k_spore_decay=0.06,
        ),
    ]
    if with_antibiotic and susceptible:
        layers.append(Antibiotic(
            name="vancomycin", susceptible=set(susceptible),
            dose_times=tuple(float(d) for d in range(10)),   # 10-day course
            dose=4.0, half_life=0.4, emax=8.0, ec50=0.5,
        ))
    # Optional bacteriophage cocktail (the failed-FMT autopsy): one layer per
    # susceptible host, so a transferred phage can knock out a key guild.
    for host in (phage_hosts or []):
        if host in present:
            layers.append(PhageInfection(
                host=host, name=f"phage_{host}",
                adsorption_rate=11.0, burst_size=60.0,
                latent_period=0.4, decay_rate=0.1, initial_titer=0.5,
            ))
    return EcologyModel(layers)


# ===========================================================================
# 4. SCENARIO RUNNER
# ===========================================================================

def run_scenario(
    model_dir: str | Path,
    abundance_path: str | Path,
    roles: Roles,
    *,
    with_antibiotic: bool = True,
    inject_models: Optional[str | Path] = None,   # donor abundance TSV -> FMT bolus
    inoculum_biomass: float = 0.05,
    total_biomass: float = 0.65,
    phage_hosts: Optional[Sequence[str]] = None,
    t_end: float = T_END,
    dt: float = DT,
) -> SimulationResult:
    """Assemble a real-MAG community + its ecology and run one dynamic-FBA scenario.

    Parameters
    ----------
    model_dir       : directory of <mag_id>.xml GEMs (from `muode build/refine`).
    abundance_path  : recipient abundance TSV (mag_id, rel_abundance).
    roles           : parsed roles.tsv (see :func:`read_roles`).
    inject_models   : donor abundance TSV; if given, its members are delivered as
                      a timed FMT bolus at FMT_DAY (they must have GEMs too, so
                      keep donor + recipient GEMs in the same model_dir).
    phage_hosts     : ids to attack with a phage cocktail (autopsy scenario).
    """
    recipient_ab = read_abundance(abundance_path)
    community = Community.from_models(
        model_dir, abundance=recipient_ab, total_biomass=total_biomass,
    )
    _check_ids(community, roles, recipient_ab)

    ecology = build_ecology(community, roles,
                            with_antibiotic=with_antibiotic, phage_hosts=phage_hosts)

    injections = None
    if inject_models is not None:
        donor_ab = read_abundance(inject_models)
        injections = [Injection.from_abundances(
            FMT_DAY, donor_ab, total_biomass=inoculum_biomass, name="fmt_bolus",
        )]

    return DynamicFBA(t_end=t_end, dt=dt, dilution_rate=DILUTION).run(
        community, gut_medium(), cdi_kinetics(),
        injections=injections, ecology=ecology,
    )


def _check_ids(community: Community, roles: Roles, abundance: Dict[str, float]) -> None:
    """Fail loudly if the three id spaces disagree (the #1 real-data pitfall)."""
    gem_ids = set(community.organism_ids)
    ab_ids = set(abundance)
    problems = []
    missing_gem = ab_ids - gem_ids
    if missing_gem:
        problems.append(f"  abundance ids with no GEM in model_dir: {sorted(missing_gem)}")
    role_orphans = roles.all_ids - gem_ids
    if role_orphans:
        problems.append(f"  role-table ids with no GEM: {sorted(role_orphans)}")
    if roles.pathogen not in gem_ids:
        problems.append(f"  pathogen '{roles.pathogen}' has no GEM in model_dir")
    unroled = gem_ids - roles.all_ids
    if unroled:
        # not fatal -- an unroled member just gets no ecology flags -- but warn
        print(f"[adapter] note: {len(unroled)} GEM(s) absent from roles.tsv "
              f"(no ecology roles): {sorted(unroled)}")
    if problems:
        raise ValueError("id-space mismatch between GEMs / abundance / roles:\n"
                         + "\n".join(problems))


# ===========================================================================
# 5. CLI
# ===========================================================================

def main(argv: Optional[Sequence[str]] = None) -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--models", required=True,
                   help="directory of <mag_id>.xml GEMs (muode build/refine output)")
    p.add_argument("--abundance", required=True,
                   help="recipient abundance TSV (mag_id, rel_abundance)")
    p.add_argument("--roles", required=True, help="roles.tsv (functional flags)")
    p.add_argument("--donor-abundance", default=None,
                   help="donor abundance TSV -> delivered as an FMT bolus at FMT_DAY")
    p.add_argument("--arm", choices=["relapse", "fmt", "phage"], default="fmt",
                   help="relapse = antibiotic, no transplant; fmt = antibiotic + "
                        "donor bolus; phage = fmt + a phage cocktail on bai hosts")
    p.add_argument("--outdir", default="./results_realdata")
    args = p.parse_args(argv)

    roles = read_roles(args.roles)
    phage_hosts = sorted(roles.bai) if args.arm == "phage" else None
    donor = args.donor_abundance if args.arm in ("fmt", "phage") else None

    result = run_scenario(
        args.models, args.abundance, roles,
        with_antibiotic=True, inject_models=donor, phage_hosts=phage_hosts,
    )

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    result.biomass.to_csv(outdir / "biomass.csv", index_label="day")
    result.metabolites.to_csv(outdir / "metabolites.csv", index_label="day")
    cd = result.biomass[roles.pathogen].iloc[-1]
    print(f"[adapter] arm={args.arm}  pathogen '{roles.pathogen}' "
          f"end biomass = {cd:.3g} gDW/L  ({'CLEARS' if cd < 1e-3 else 'PERSISTS'})")
    print(f"[adapter] wrote {outdir}/biomass.csv, metabolites.csv")
    print("[adapter] plot with the same helpers the toy scenarios use "
          "(common.py: PALETTE, shade_dosing, biomass_sum).")


if __name__ == "__main__":
    main()
