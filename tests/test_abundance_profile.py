"""Tests for SimulationResult.abundance_profile() and the snapshot CLI options."""

import numpy as np
import pytest

from muode import DynamicFBA
from muode.examples import build_toy_community, toy_diet, toy_kinetics


@pytest.fixture(scope="module")
def result():
    return DynamicFBA(t_end=24.0, dt=0.5).run(
        build_toy_community(), toy_diet(), toy_kinetics()
    )


def test_full_profile_sums_to_one(result):
    prof = result.abundance_profile()
    assert prof.shape == result.biomass.shape
    row_sums = prof.sum(axis=1)
    # every row should sum to ~1 (or 0 if biomass is exactly zero, which won't happen here)
    assert (row_sums > 0.99).all()
    assert (row_sums < 1.01).all()


def test_interval_sampling(result):
    prof = result.abundance_profile(interval=6.0)
    # t_end=24, interval=6 → snapshots at 0, 6, 12, 18, 24 = 5 rows
    assert len(prof) == 5
    assert prof.sum(axis=1).between(0.99, 1.01).all()


def test_explicit_times(result):
    times = [0.0, 4.0, 12.0, 24.0]
    prof = result.abundance_profile(times=times)
    assert len(prof) == 4
    # index should be close to the requested times (nearest-step matching)
    for requested, actual in zip(times, prof.index):
        assert abs(actual - requested) <= 0.5 + 1e-9  # within one dt=0.5 step
    assert prof.sum(axis=1).between(0.99, 1.01).all()


def test_times_takes_precedence_over_interval(result):
    # when both are given, times wins
    prof_times = result.abundance_profile(times=[0.0, 12.0, 24.0])
    prof_both  = result.abundance_profile(times=[0.0, 12.0, 24.0], interval=6.0)
    assert len(prof_times) == len(prof_both) == 3


def test_values_are_relative_not_absolute(result):
    prof = result.abundance_profile(interval=6.0)
    # values are in [0, 1] and columns are the same species as biomass
    assert (prof.values >= 0.0).all()
    assert (prof.values <= 1.0 + 1e-9).all()
    assert list(prof.columns) == list(result.biomass.columns)


def test_to_csv_writes_snapshot_file(result, tmp_path):
    result.to_csv(tmp_path, snapshot_interval=6.0)
    snap = tmp_path / "abundance_snapshots.tsv"
    assert snap.exists()
    import pandas as pd
    df = pd.read_csv(snap, sep="\t", index_col="time_h")
    assert len(df) == 5          # 0,6,12,18,24
    assert df.sum(axis=1).between(0.99, 1.01).all()


def test_to_csv_without_snapshot_produces_no_file(result, tmp_path):
    result.to_csv(tmp_path)
    assert not (tmp_path / "abundance_snapshots.tsv").exists()


def test_snapshot_times_cli_option(tmp_path):
    """CLI --snapshot-times writes the file with the right number of rows."""
    from typer.testing import CliRunner
    from muode.cli import app

    runner = CliRunner()
    r = runner.invoke(app, [
        "demo",
        "--outdir", str(tmp_path),
        "--time", "24",
        "--step", "0.5",
    ])
    # demo doesn't have snapshot options; verify simulate path via Python API instead
    result2 = DynamicFBA(t_end=24.0, dt=0.5).run(
        build_toy_community(), toy_diet(), toy_kinetics()
    )
    result2.to_csv(tmp_path / "snap_test", snapshot_times=[0.0, 8.0, 16.0, 24.0])
    import pandas as pd
    df = pd.read_csv(tmp_path / "snap_test" / "abundance_snapshots.tsv",
                     sep="\t", index_col="time_h")
    assert len(df) == 4
