"""muODE command-line interface.

Mirrors the five pipeline phases from the README plus a self-contained ``demo``.
The CLI is a thin layer over the library: every command builds the relevant
objects and calls into :mod:`muode`.  The heavyweight phases (``build``,
``refine``) shell out to external tools and are normally driven by the Snakemake
workflow for scale; ``assemble`` / ``simulate`` / ``perturb`` / ``demo`` run the
native engine and work out of the box.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import List, Optional

import typer
from rich.console import Console
from rich.table import Table

from muode import __version__

app = typer.Typer(add_completion=False, help="µODE: ODE-based microbial community dynamics from MAGs.")
console = Console()


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def _load_diet(spec: str):
    from muode.diet import Diet, available_presets, load_preset

    if spec.endswith(".csv") and Path(spec).exists():
        return Diet.from_csv(spec)
    try:
        return load_preset(spec)
    except KeyError as exc:
        raise typer.BadParameter(f"{exc}; or pass a path to a diet CSV") from exc


def _kinetics(default_vmax: float, default_km: float):
    from muode.kinetics import KineticParameters

    return KineticParameters(default_vmax=default_vmax, default_km=default_km)


def _save_outputs(result, outdir: Path, make_figures: bool = True) -> None:
    result.to_csv(outdir)
    Path(outdir, "meta.json").write_text(json.dumps(result.meta, indent=2, default=str))
    if make_figures:
        try:
            from muode import viz

            viz.save_all(result, outdir)
        except ImportError:
            console.print("[yellow]matplotlib not installed; skipping figures "
                          "(`pip install muode[viz]`).[/yellow]")


def _summary(result) -> None:
    table = Table(title="Final community state")
    table.add_column("species")
    table.add_column("initial (gDW/L)", justify="right")
    table.add_column("final (gDW/L)", justify="right")
    table.add_column("fold", justify="right")
    init, final = result.biomass.iloc[0], result.biomass.iloc[-1]
    for sp in result.biomass.columns:
        fold = final[sp] / init[sp] if init[sp] > 0 else float("nan")
        table.add_row(sp, f"{init[sp]:.4g}", f"{final[sp]:.4g}", f"{fold:.2f}x")
    console.print(table)
    extinct = result.extinct()
    if extinct:
        console.print(f"[red]Extinct species:[/red] {', '.join(extinct)}")
    cf = result.cross_feeding()
    if not cf.empty:
        console.print(f"[green]{len(cf)} cross-feeding interaction(s) detected.[/green]")


# ---------------------------------------------------------------------------
# commands
# ---------------------------------------------------------------------------


@app.command()
def version() -> None:
    """Print the muODE version."""
    console.print(f"muODE {__version__}")


@app.command()
def demo(
    outdir: Path = typer.Option("results/demo", help="Output directory."),
    t_end: float = typer.Option(24.0, help="Simulated time (h)."),
    dt: float = typer.Option(0.05, help="Integration step (h)."),
    remove: Optional[str] = typer.Option(None, help="Remove a species to see secondary extinction, e.g. 'A_glucose'."),
) -> None:
    """Run the bundled two-species cross-feeding demo (no external data needed)."""
    from muode.dfba import DynamicFBA
    from muode.examples import build_toy_community, toy_diet, toy_kinetics
    from muode.perturb import Perturbation

    outdir.mkdir(parents=True, exist_ok=True)
    pert = Perturbation.remove_species([remove]) if remove else None
    result = DynamicFBA(t_end=t_end, dt=dt).run(build_toy_community(), toy_diet(), toy_kinetics(), pert)
    _save_outputs(result, outdir)
    _summary(result)
    console.print(f"[bold]Wrote results to {outdir}[/bold]")


@app.command()
def build(
    mags: Path = typer.Option(..., help="Directory of MAG FASTA files."),
    outdir: Path = typer.Option("models/draft_gems", help="Output directory for GEMs."),
    engine: str = typer.Option("carveme", help="carveme | gapseq."),
    universe: str = typer.Option("bacteria", help="CarveMe universe / template."),
    gapfill_media: Optional[str] = typer.Option(None, help="Medium to gap-fill on during build."),
    solver: Optional[str] = typer.Option(None, help="LP solver (gurobi/cplex)."),
    pattern: str = typer.Option("*.fna", help="Glob for MAG files."),
) -> None:
    """Phase 1 -- reconstruct draft GEMs from MAGs (one model per genome)."""
    from muode.reconstruct import reconstruct_mag

    outdir.mkdir(parents=True, exist_ok=True)
    genomes = sorted(Path(mags).glob(pattern))
    if not genomes:
        raise typer.BadParameter(f"no files matching {pattern} in {mags}")
    for g in genomes:
        out = outdir / f"{g.stem}.xml"
        console.print(f"reconstructing [cyan]{g.stem}[/cyan] -> {out}")
        reconstruct_mag(g, out, engine=engine, universe=universe,
                        gapfill_media=gapfill_media, solver=solver)
    console.print(f"[bold]Built {len(genomes)} model(s).[/bold] "
                  "For >100 MAGs use the Snakemake workflow for parallelism.")


@app.command()
def refine(
    models: Path = typer.Option(..., help="Directory of draft GEMs."),
    outdir: Path = typer.Option("models/kinetic_gems", help="Output directory."),
    universal: Optional[Path] = typer.Option(None, help="Universal model (SBML) for LP gap-filling."),
    predict_kinetics: bool = typer.Option(False, help="Predict Km + kcat and write {stem}.kinetics.json."),
    predictor: str = typer.Option("heuristic", help="heuristic (default, no deps) | dlkcat | km-ml (need the 'ml' extra)."),
) -> None:
    """Phase 2/3 -- gap-fill models and (optionally) predict kinetic parameters."""
    import cobra

    from muode.gapfill import ensure_biomass

    outdir.mkdir(parents=True, exist_ok=True)
    uni = cobra.io.read_sbml_model(str(universal)) if universal else None

    pred = None
    if predict_kinetics:
        from muode.predict import get_predictor

        try:
            pred = get_predictor(predictor)
        except Exception as exc:  # missing 'ml' extra, unknown name, ...
            raise typer.BadParameter(str(exc)) from exc

    for path in sorted(Path(models).glob("*.xml")):
        model = cobra.io.read_sbml_model(str(path))
        info = ensure_biomass(model, universal=uni)
        cobra.io.write_sbml_model(model, str(outdir / path.name))
        flag = "ok" if info["grows_now"] else "STILL DEAD"
        msg = f"{path.stem}: +{len(info['reactions_added'])} rxns [{flag}]"
        if pred is not None:
            from muode.predict import refine_kinetics

            kin = refine_kinetics(model, path.stem, predictor=pred)
            kin.to_json(outdir / f"{path.stem}.kinetics.json")
            msg += f"  ·  kinetics[{pred.name}]: {len(kin.kcat)} kcat, {len(kin.overrides)} Km"
        console.print(msg)


@app.command()
def assemble(
    models: Path = typer.Option(..., help="Directory of refined GEMs."),
    abundance: Optional[Path] = typer.Option(None, help="MAG abundance TSV (2-column or MetaSBT profile)."),
    metasbt: bool = typer.Option(False, help="Parse --abundance as a MetaSBT profile (auto-detect columns)."),
    diet: str = typer.Option("western_gut", help="Diet preset name or CSV path."),
    outdir: Path = typer.Option("simulation_env", help="Output directory."),
    total_biomass: float = typer.Option(0.01, help="Total community biomass (gDW/L)."),
    max_species: Optional[int] = typer.Option(None, help="Keep only the N most abundant MAGs."),
    min_abundance: Optional[float] = typer.Option(None, help="Drop MAGs below this relative abundance."),
    coverage: Optional[float] = typer.Option(None, help="Keep the fewest top MAGs reaching this cumulative abundance (e.g. 0.95)."),
) -> None:
    """Phase 4a -- write a community manifest (models + abundances + diet)."""
    from muode.io_utils import list_models, read_abundance

    outdir.mkdir(parents=True, exist_ok=True)
    files = [str(p) for p in list_models(models)]
    if not files:
        raise typer.BadParameter(f"no models found in {models}")

    if abundance and metasbt:
        from muode.metasbt import read_metasbt_profile

        abund = read_metasbt_profile(abundance).abundances
    elif abundance:
        abund = read_abundance(abundance)
    else:
        abund = {Path(f).stem: 1.0 for f in files}

    # abundance-aware subsampling for very large communities
    if any(v is not None for v in (max_species, min_abundance, coverage)):
        from muode.subsample import select_by_abundance

        stem_to_file = {Path(f).stem: f for f in files}
        present = {s: abund.get(s, 0.0) for s in stem_to_file}
        keep = select_by_abundance(present, top_n=max_species,
                                   min_abundance=min_abundance, coverage=coverage)
        dropped = len(files) - len(keep)
        files = [stem_to_file[s] for s in keep if s in stem_to_file]
        if dropped:
            console.print(f"[yellow]Subsampled to {len(files)} MAG(s) by abundance "
                          f"(dropped {dropped}).[/yellow]")

    manifest = {
        "models": files,
        "abundances": {Path(f).stem: abund.get(Path(f).stem, 0.0) for f in files},
        "diet": diet,
        "total_biomass": total_biomass,
    }
    out = outdir / "community.json"
    out.write_text(json.dumps(manifest, indent=2))
    console.print(f"[bold]Wrote community manifest with {len(files)} model(s) to {out}[/bold]")


def _community_from_manifest(path: Path):
    from muode.community import Community
    from muode.organism import CobraOrganism

    manifest = json.loads(Path(path).read_text())
    organisms = [CobraOrganism.from_file(p, id=Path(p).stem) for p in manifest["models"]]
    comm = Community(organisms, manifest.get("abundances", {}),
                     total_biomass=manifest.get("total_biomass", 0.01))
    return comm, manifest.get("diet", "western_gut")


def _load_kinetics(spec: Path):
    """Load and merge predicted kinetics from a file or a directory of them."""
    from muode.kinetics import KineticParameters

    spec = Path(spec)
    files = sorted(spec.glob("*.kinetics.json")) if spec.is_dir() else ([spec] if spec.exists() else [])
    merged = KineticParameters()
    for f in files:
        merged.merge(KineticParameters.from_json(f))
    return merged, len(files)


@app.command()
def simulate(
    community: Optional[Path] = typer.Option(None, help="community.json manifest from `assemble`."),
    models: Optional[Path] = typer.Option(None, help="Alternatively, a directory of GEMs."),
    abundance: Optional[Path] = typer.Option(None, help="Abundance TSV (with --models)."),
    diet: str = typer.Option("western_gut", help="Diet preset name or CSV path."),
    time: float = typer.Option(24.0, "--time", help="Simulated time (h)."),
    step: float = typer.Option(0.1, "--step", help="Integration step (h)."),
    death_rate: float = typer.Option(0.0, help="First-order biomass death (1/h)."),
    dilution_rate: float = typer.Option(0.0, help="Chemostat dilution D (1/h)."),
    total_biomass: float = typer.Option(0.01, help="Total community biomass (gDW/L)."),
    default_vmax: float = typer.Option(10.0, help="Default Vmax (mmol/gDW/h)."),
    default_km: float = typer.Option(0.01, help="Default Km (mmol/L)."),
    kinetics: Optional[Path] = typer.Option(None, help="Predicted kinetics: a *.kinetics.json file or a directory of them (from `refine --predict-kinetics`)."),
    enzyme_constraints: bool = typer.Option(False, help="Apply GECKO-lite kcat caps to intracellular reactions (needs --kinetics with kcat)."),
    outdir: Path = typer.Option("results", help="Output directory."),
) -> None:
    """Phase 4 -- run the dynamic community simulation."""
    from muode.community import Community
    from muode.dfba import DynamicFBA
    from muode.io_utils import read_abundance

    if community:
        comm, diet_spec = _community_from_manifest(community)
        diet = diet_spec
    elif models:
        abund = read_abundance(abundance) if abundance else None
        comm = Community.from_models(models, abundance=abund, total_biomass=total_biomass)
    else:
        raise typer.BadParameter("provide either --community or --models")

    kin = _kinetics(default_vmax, default_km)
    if kinetics:
        loaded, n = _load_kinetics(kinetics)
        kin.merge(loaded)
        console.print(f"[green]Loaded refined kinetics from {n} file(s) "
                      f"({len(kin.kcat)} kcat, {len(kin.overrides)} Km).[/green]")
    if enzyme_constraints:
        if not kin.kcat:
            console.print("[yellow]--enzyme-constraints set but no kcat available "
                          "(pass --kinetics); skipping.[/yellow]")
        for o in comm.organisms:
            if hasattr(o, "apply_enzyme_constraints"):
                rep = o.apply_enzyme_constraints(kin, organism_id=o.id)
                console.print(f"  {o.id}: enzyme-constrained {rep['n_constrained']} reaction(s)")

    engine = DynamicFBA(t_end=time, dt=step, death_rate=death_rate, dilution_rate=dilution_rate)
    result = engine.run(comm, _load_diet(diet), kin)
    outdir.mkdir(parents=True, exist_ok=True)
    _save_outputs(result, outdir)
    _summary(result)
    console.print(f"[bold]Wrote results to {outdir}[/bold]")


@app.command()
def perturb(
    community: Path = typer.Option(..., help="community.json manifest from `assemble`."),
    target_pathway: Optional[str] = typer.Option(None, help="Pathway/subsystem to attenuate."),
    efficacy: float = typer.Option(0.95, help="Fraction of capacity removed (1.0 = knockout)."),
    remove_species: Optional[str] = typer.Option(None, help="Comma-separated species ids to remove."),
    knockout: Optional[str] = typer.Option(None, help="Comma-separated reaction ids to knock out."),
    diet: Optional[str] = typer.Option(None, help="Override diet preset/CSV."),
    time: float = typer.Option(24.0, "--time"),
    step: float = typer.Option(0.1, "--step"),
    default_vmax: float = typer.Option(10.0),
    default_km: float = typer.Option(0.01),
    outdir: Path = typer.Option("results_perturbation", help="Output directory."),
) -> None:
    """Phase 5 -- run a simulation under a perturbation (antibiotic / knockout)."""
    from muode.dfba import DynamicFBA
    from muode.perturb import Perturbation, PerturbationTarget

    comm, diet_spec = _community_from_manifest(community)
    diet = diet or diet_spec

    targets: List[PerturbationTarget] = []
    if target_pathway:
        targets.append(PerturbationTarget(subsystem=target_pathway, efficacy=efficacy))
    if knockout:
        targets.append(PerturbationTarget(reactions=knockout.split(","), efficacy=1.0))
    if remove_species:
        targets.append(PerturbationTarget(organisms=remove_species.split(","), remove_organism=True))
    if not targets:
        raise typer.BadParameter("specify --target-pathway, --knockout and/or --remove-species")
    pert = Perturbation(name="cli_perturbation", targets=targets)

    engine = DynamicFBA(t_end=time, dt=step)
    result = engine.run(comm, _load_diet(diet), _kinetics(default_vmax, default_km), pert)
    outdir.mkdir(parents=True, exist_ok=True)
    _save_outputs(result, outdir)
    console.print(f"[bold]Perturbation:[/bold] {pert.describe()}")
    _summary(result)
    console.print(f"[bold]Wrote results to {outdir}[/bold]")


@app.command()
def qc(
    model: Path = typer.Option(..., help="SBML model to sanity-check."),
) -> None:
    """Run fast structural QC on a reconstructed model."""
    import cobra

    from muode.qc import sanity_check_model

    report = sanity_check_model(cobra.io.read_sbml_model(str(model)))
    console.print_json(json.dumps(report, default=str))
    if not report["passed"]:
        raise typer.Exit(code=1)


@app.command()
def validate(
    results: Path = typer.Option(..., help="Results directory written by `simulate`/`demo`."),
    expected: Path = typer.Option(..., help="Benchmark expectation YAML (see examples/benchmarks/)."),
    report: Optional[Path] = typer.Option(None, help="Write the JSON validation report here."),
) -> None:
    """Validate a simulation against a known/expected community (Phase 6)."""
    from muode.validate import BenchmarkExpectation, validate_outputs

    expectation = BenchmarkExpectation.from_file(expected)
    rep = validate_outputs(results, expectation)
    console.print_json(json.dumps(rep.to_dict(), default=str))
    if report:
        rep.to_json(report)
    status = "[green]PASSED[/green]" if rep.passed else "[red]FAILED[/red]"
    console.print(f"[bold]Benchmark '{rep.name}': {status}[/bold]")
    if not rep.passed:
        raise typer.Exit(code=1)


if __name__ == "__main__":
    app()
