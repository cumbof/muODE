"""Tests for the static HTML report builder (`muode.report`).

Pure stdlib: reads CSV/PNG/JSON and writes one HTML file, so it runs anywhere.
"""

import base64

from muode.report import build_report


def _write_run(d):
    (d / "biomass.csv").write_text(
        "time_h,A_glucose,B_acetate\n"
        "0.0,0.01,0.01\n"
        "12.0,0.40,0.00000001\n"
    )
    (d / "metabolites.csv").write_text(
        "time_h,glc_e,ac_e\n"
        "0.0,20.0,0.0\n"
        "12.0,0.5,8.0\n"
    )
    (d / "cross_feeding.csv").write_text(
        "producer,metabolite,consumer,strength\n"
        "A_glucose,ac_e,B_acetate,1.5\n"
    )
    (d / "meta.json").write_text('{"t_end": 12.0, "dt": 0.1, "diet": "toy"}')
    # a 1x1 PNG so the embed path is exercised
    png = base64.b64decode(
        "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mNk+M9QDwADhgGAWjR9awAAAABJRU5ErkJggg=="
    )
    (d / "biomass.png").write_bytes(png)


def test_build_report_renders_sections(tmp_path):
    _write_run(tmp_path)
    out = build_report(tmp_path)
    assert out.exists() and out.name == "report.html"
    doc = out.read_text()

    # composition + fold + extinction
    assert "A_glucose" in doc and "B_acetate" in doc
    assert "40.00x" in doc                  # 0.40 / 0.01
    assert "Extinct: B_acetate" in doc      # collapsed below 1e-3 of start

    # cross-feeding and metabolites
    assert "Cross-feeding interactions" in doc
    assert "Most dynamic metabolites" in doc and "ac_e" in doc

    # run metadata and the embedded figure
    assert "Run parameters" in doc and "t_end" in doc
    assert "data:image/png;base64," in doc


def test_build_report_custom_output_path_and_sparse_dir(tmp_path):
    # only a biomass.csv present: report still builds, optional sections skipped
    (tmp_path / "biomass.csv").write_text("time_h,X\n0.0,0.01\n5.0,0.05\n")
    out = tmp_path / "sub" / "custom.html"
    result = build_report(tmp_path, out)
    assert result == out and out.exists()
    doc = out.read_text()
    assert "Community composition" in doc
    assert "Cross-feeding interactions" not in doc
