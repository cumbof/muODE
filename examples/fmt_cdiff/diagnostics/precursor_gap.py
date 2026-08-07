import cobra, logging, warnings
warnings.filterwarnings("ignore"); logging.disable(logging.CRITICAL)
from pathlib import Path
from cobra.util import linear_reaction_coefficients
GEMS=Path("examples/fmt_cdiff/gems")
for code in ["C_scindens_ATCC35704"]:
    m=cobra.io.read_sbml_model(str(GEMS/f"{code}.xml.gz"))
    for r in m.exchanges: r.lower_bound=-1000.0   # everything open
    bio=list(linear_reaction_coefficients(m))[0]
    print(f"{code}: biomass rxn = {bio.id}, max growth (all open) = {m.slim_optimize():.4f}")
    # for each biomass PRECURSOR (reactant), can the cell make it? add a demand sink, maximize
    reactants=[met for met,coef in bio.metabolites.items() if coef<0]
    print(f"  biomass has {len(reactants)} precursors; testing which are BLOCKED...")
    blocked=[]
    for met in reactants:
        with m:
            try:
                dm=m.add_boundary(met, type="demand")
                m.objective=dm
                v=m.slim_optimize() or 0.0
            except Exception:
                v=-1
        if v is not None and v < 1e-6:
            blocked.append((met.id, met.name[:40]))
    print(f"  BLOCKED precursors ({len(blocked)} of {len(reactants)}):")
    for mid,nm in blocked[:25]: print(f"    {mid:16} {nm}")
