"""Bundle a simulation results directory into one self-contained HTML report.

`muode simulate` / `demo` / `perturb` write a directory of CSVs and PNGs. This
module stitches them into a single, portable ``report.html`` — summary tables
(final composition, fold-changes, extinctions, cross-feeding edges, run
metadata) with the figures embedded as base64 — so a whole run can be reviewed
or shared as one file, with no server and no extra dependencies.

It is deliberately pure stdlib (``csv``/``json``/``base64``/``html``): a static
*report*, not a live dashboard. An interactive, server-backed dashboard (with
zoomable time sliders over the same data) remains future work; the CSVs and the
``.npz`` spatial fields are exactly the data such a dashboard would render.
"""

from __future__ import annotations

import base64
import csv
import html
import json
from pathlib import Path
from typing import List, Optional, Tuple

# figures emitted by muode.viz, in display order
_FIGURES = ["biomass.png", "metabolites.png", "cross_feeding.png",
            "spatial_total_biomass.png", "spatial_final.png"]


def _read_csv(path: Path) -> Tuple[List[str], List[List[str]]]:
    """Return ``(header, rows)`` for a CSV; ``([], [])`` if it is missing/empty."""
    if not path.exists():
        return [], []
    with path.open(newline="") as fh:
        rows = list(csv.reader(fh))
    if not rows:
        return [], []
    return rows[0], rows[1:]


def _fmt(x: float) -> str:
    return f"{x:.4g}"


def _table(headers: List[str], rows: List[List[str]]) -> str:
    head = "".join(f"<th>{html.escape(str(h))}</th>" for h in headers)
    body = "".join(
        "<tr>" + "".join(f"<td>{html.escape(str(c))}</td>" for c in r) + "</tr>"
        for r in rows
    )
    return f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"


def _biomass_summary(results: Path) -> str:
    header, rows = _read_csv(results / "biomass.csv")
    if not rows:
        return "<p><em>No biomass.csv found.</em></p>"
    species = header[1:]
    first = [float(v) for v in rows[0][1:]]
    last = [float(v) for v in rows[-1][1:]]
    total_final = sum(last) or 1.0
    table_rows = []
    extinct = []
    for sp, x0, x1 in zip(species, first, last):
        fold = (x1 / x0) if x0 > 0 else float("nan")
        rel = x1 / total_final
        gone = x1 < 1e-6 or (x0 > 0 and x1 < 1e-3 * x0)
        if gone:
            extinct.append(sp)
        table_rows.append([
            sp, _fmt(x0), _fmt(x1), f"{fold:.2f}x" if fold == fold else "—",
            f"{100 * rel:.1f}%", "yes" if gone else "",
        ])
    out = ["<h2>Community composition</h2>"]
    out.append(_table(
        ["species", "initial (gDW/L)", "final (gDW/L)", "fold", "final share", "extinct"],
        table_rows,
    ))
    if extinct:
        out.append(f"<p class='warn'>Extinct: {', '.join(html.escape(s) for s in extinct)}</p>")
    return "".join(out)


def _cross_feeding_summary(results: Path) -> str:
    header, rows = _read_csv(results / "cross_feeding.csv")
    if not rows:
        return ""
    return "<h2>Cross-feeding interactions</h2>" + _table(header, rows)


def _metabolite_summary(results: Path, top: int = 12) -> str:
    header, rows = _read_csv(results / "metabolites.csv")
    if not rows:
        return ""
    mets = header[1:]
    first = [float(v) for v in rows[0][1:]]
    last = [float(v) for v in rows[-1][1:]]
    # rank by absolute change over the run
    ranked = sorted(zip(mets, first, last), key=lambda t: abs(t[2] - t[1]), reverse=True)
    table_rows = [[m, _fmt(a), _fmt(b), _fmt(b - a)] for m, a, b in ranked[:top] if abs(b - a) > 1e-9]
    if not table_rows:
        return ""
    return ("<h2>Most dynamic metabolites</h2>"
            + _table(["metabolite", "initial (mmol/L)", "final (mmol/L)", "Δ"], table_rows))


def _meta_summary(results: Path) -> str:
    path = results / "meta.json"
    if not path.exists():
        return ""
    meta = json.loads(path.read_text())
    rows = [[k, json.dumps(v) if isinstance(v, (dict, list)) else str(v)] for k, v in meta.items()]
    return "<h2>Run parameters</h2>" + _table(["key", "value"], rows)


def _figures(results: Path) -> str:
    blocks = []
    for name in _FIGURES:
        p = results / name
        if not p.exists():
            continue
        b64 = base64.b64encode(p.read_bytes()).decode("ascii")
        title = html.escape(name.replace(".png", "").replace("_", " "))
        blocks.append(
            f"<figure><figcaption>{title}</figcaption>"
            f"<img alt='{title}' src='data:image/png;base64,{b64}'></figure>"
        )
    if not blocks:
        return ""
    return "<h2>Figures</h2><div class='figs'>" + "".join(blocks) + "</div>"


_CSS = """
body{font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;
 max-width:1000px;margin:2rem auto;padding:0 1rem;color:#222;line-height:1.5}
h1{border-bottom:2px solid #9ecae1;padding-bottom:.3rem}
h2{margin-top:2rem;color:#08519c}
table{border-collapse:collapse;margin:.5rem 0;font-size:.9rem}
th,td{border:1px solid #ddd;padding:.3rem .6rem;text-align:right}
th{background:#f0f5fb;text-align:center}td:first-child,th:first-child{text-align:left}
.warn{color:#a50f15;font-weight:600}
.figs{display:flex;flex-wrap:wrap;gap:1rem}
figure{margin:0;max-width:480px}figure img{max-width:100%;border:1px solid #eee}
figcaption{font-weight:600;color:#08519c;margin-bottom:.3rem}
footer{margin-top:3rem;color:#888;font-size:.8rem;border-top:1px solid #eee;padding-top:.5rem}
"""


def build_report(results_dir: str | Path, output_html: Optional[str | Path] = None) -> Path:
    """Render a results directory into one self-contained ``report.html``.

    Parameters
    ----------
    results_dir:
        A directory written by ``muode simulate`` / ``demo`` / ``perturb``.
    output_html:
        Destination file; defaults to ``<results_dir>/report.html``.

    Returns the path written.
    """
    results = Path(results_dir)
    if not results.exists():
        raise FileNotFoundError(results)
    out = Path(output_html) if output_html else results / "report.html"
    out.parent.mkdir(parents=True, exist_ok=True)

    title = f"muODE report — {results.name}"
    sections = [
        _biomass_summary(results),
        _metabolite_summary(results),
        _cross_feeding_summary(results),
        _figures(results),
        _meta_summary(results),
    ]
    body = "\n".join(s for s in sections if s)
    doc = (
        "<!doctype html><html lang='en'><head><meta charset='utf-8'>"
        f"<title>{html.escape(title)}</title><style>{_CSS}</style></head><body>"
        f"<h1>{html.escape(title)}</h1>"
        f"{body}"
        "<footer>Generated by muODE (muode.report). A static report; values are "
        "model outputs — mechanistic hypotheses, not validated predictions.</footer>"
        "</body></html>"
    )
    out.write_text(doc, encoding="utf-8")
    return out
