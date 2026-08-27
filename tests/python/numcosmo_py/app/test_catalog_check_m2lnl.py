#!/usr/bin/env python
#
# test_catalog_check_m2lnl.py
#
# Sat Aug 1 2026
# Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
#
# numcosmo is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# numcosmo is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program. If not, see <http://www.gnu.org/licenses/>.

"""CLI-level tests for the `catalog check-m2lnl` command.

Exercises the full command through the actual Typer CLI, on a real (small)
MCMC catalog: recomputing -2ln(L) row by row and comparing against the
value the catalog already has stored.
"""

from pathlib import Path
from typing import Tuple

import pytest
from typer.testing import CliRunner

pytest.importorskip("astropy")
pytest.importorskip("getdist")
# flake8: noqa: E402
# pylint: disable=wrong-import-position

from numcosmo_py import Ncm
from numcosmo_py.app import app

pytestmark = pytest.mark.app
runner = CliRunner()

Ncm.cfg_init()


@pytest.fixture(name="mcmc_catalog")
def fixture_mcmc_catalog(tmp_path: Path) -> Tuple[Path, Path]:
    """Run a small MCMC (2-parameter MVN) and return (experiment, catalog) paths."""
    rng = Ncm.RNG.seeded_new(None, 1234)
    model_mvnd = Ncm.ModelMVND.new(2)
    mset = Ncm.MSet.new_array([model_mvnd])
    mset.param_set_all_ftype(Ncm.ParamType.FREE)

    data_mvnd = Ncm.DataGaussCovMVND.new_full(2, 1.0, 2.0, 30.0, 0.0, 0.0, rng)
    likelihood = Ncm.Likelihood.new(dset=Ncm.Dataset.new_array([data_mvnd]))

    experiment = tmp_path / "simple_experiment.yaml"
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)

    experiment_dict = Ncm.ObjDictStr.new()
    experiment_dict.add("model-set", mset)
    experiment_dict.add("likelihood", likelihood)
    ser.dict_str_to_yaml_file(experiment_dict, experiment.as_posix())

    output = experiment.with_suffix(".out.yaml")
    result = runner.invoke(
        app,
        ["run", "mcmc", "apes", experiment.as_posix(), "--output", output.as_posix()],
    )
    if result.exit_code != 0:
        raise result.exception

    catalog = output.absolute().with_suffix(".mcmc.fits")
    assert catalog.exists()

    return experiment, catalog


def test_check_m2lnl_matches_when_same_experiment(mcmc_catalog):
    """Recomputing -2ln(L) with the SAME experiment that produced the
    catalog must reproduce the stored value for every row: 0 mismatches,
    no "Worst mismatches" table."""
    experiment, catalog = mcmc_catalog

    result = runner.invoke(
        app,
        ["catalog", "check-m2lnl", experiment.as_posix(), catalog.as_posix()],
    )

    assert result.exit_code == 0, result.output
    assert "mismatched rows" in result.output
    assert "0/" in result.output
    assert "Worst mismatches" not in result.output


def test_check_m2lnl_stride_and_max_rows(mcmc_catalog):
    """--stride/--max-rows limit how many rows are actually recomputed."""
    experiment, catalog = mcmc_catalog

    result = runner.invoke(
        app,
        [
            "catalog",
            "check-m2lnl",
            experiment.as_posix(),
            catalog.as_posix(),
            "--stride",
            "5",
            "--max-rows",
            "3",
        ],
    )

    assert result.exit_code == 0, result.output
    assert "Checked 3 row(s)" in result.output


def test_check_m2lnl_reports_mismatches(mcmc_catalog, monkeypatch: pytest.MonkeyPatch):
    """A recomputed -2ln(L) that disagrees with the stored value is flagged
    as a mismatch and listed in the "Worst mismatches" table."""
    experiment, catalog = mcmc_catalog

    real_m2lnl_val = Ncm.Likelihood.m2lnL_val

    def fake_m2lnl_val(self, mset):
        return real_m2lnl_val(self, mset) + 10.0

    monkeypatch.setattr(Ncm.Likelihood, "m2lnL_val", fake_m2lnl_val)

    result = runner.invoke(
        app,
        [
            "catalog",
            "check-m2lnl",
            experiment.as_posix(),
            catalog.as_posix(),
            "--max-rows",
            "5",
        ],
    )

    assert result.exit_code == 0, result.output
    assert "mismatched rows" in result.output
    assert "5/5" in result.output
    assert "Worst mismatches" in result.output


def test_check_m2lnl_max_report_limits_table_rows(
    mcmc_catalog, monkeypatch: pytest.MonkeyPatch
):
    """--max-report caps how many mismatching rows are printed individually,
    even when more rows mismatch."""
    experiment, catalog = mcmc_catalog

    real_m2lnl_val = Ncm.Likelihood.m2lnL_val

    def fake_m2lnl_val(self, mset):
        return real_m2lnl_val(self, mset) + 10.0

    monkeypatch.setattr(Ncm.Likelihood, "m2lnL_val", fake_m2lnl_val)

    result = runner.invoke(
        app,
        [
            "catalog",
            "check-m2lnl",
            experiment.as_posix(),
            catalog.as_posix(),
            "--max-rows",
            "5",
            "--max-report",
            "2",
        ],
    )

    assert result.exit_code == 0, result.output
    assert "5/5" in result.output
    table_rows = [
        line
        for line in result.output.splitlines()
        if "│" in line and line.split("│")[1].strip().isdigit()
    ]
    assert len(table_rows) == 2


def test_check_m2lnl_worst_loop_stops_at_first_matching_row(
    mcmc_catalog, monkeypatch: pytest.MonkeyPatch
):
    """The "Worst mismatches" table is built from ALL rows sorted by
    |delta|, not just the mismatching ones: once a non-mismatching row is
    reached (impossible to see further mismatches after it, since the list
    is sorted descending) the loop stops early, before --max-report would
    have. Only the first 2 of 5 rows are made to mismatch here, with
    --max-report larger than that, so this exercises that early stop
    specifically (as opposed to the --max-report cap, covered separately
    above)."""
    experiment, catalog = mcmc_catalog

    real_m2lnl_val = Ncm.Likelihood.m2lnL_val
    call_count = {"n": 0}

    def fake_m2lnl_val(self, mset):
        call_count["n"] += 1
        value = real_m2lnl_val(self, mset)
        return value + 10.0 if call_count["n"] <= 2 else value

    monkeypatch.setattr(Ncm.Likelihood, "m2lnL_val", fake_m2lnl_val)

    result = runner.invoke(
        app,
        [
            "catalog",
            "check-m2lnl",
            experiment.as_posix(),
            catalog.as_posix(),
            "--max-rows",
            "5",
            "--max-report",
            "10",
        ],
    )

    assert result.exit_code == 0, result.output
    assert "2/5" in result.output
    table_rows = [
        line
        for line in result.output.splitlines()
        if "│" in line and line.split("│")[1].strip().isdigit()
    ]
    assert len(table_rows) == 2


def test_check_m2lnl_incompatible_fparams_raises(mcmc_catalog, tmp_path: Path):
    """check-m2lnl refuses to run when the experiment's number of free
    parameters does not match the catalog's."""
    _, catalog = mcmc_catalog

    rng = Ncm.RNG.seeded_new(None, 1234)
    model_mvnd = Ncm.ModelMVND.new(3)
    mset = Ncm.MSet.new_array([model_mvnd])
    mset.param_set_all_ftype(Ncm.ParamType.FREE)

    data_mvnd = Ncm.DataGaussCovMVND.new_full(3, 1.0, 2.0, 30.0, 0.0, 0.0, rng)
    likelihood = Ncm.Likelihood.new(dset=Ncm.Dataset.new_array([data_mvnd]))

    other_experiment = tmp_path / "other_experiment.yaml"
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    experiment_dict = Ncm.ObjDictStr.new()
    experiment_dict.add("model-set", mset)
    experiment_dict.add("likelihood", likelihood)
    ser.dict_str_to_yaml_file(experiment_dict, other_experiment.as_posix())

    result = runner.invoke(
        app,
        ["catalog", "check-m2lnl", other_experiment.as_posix(), catalog.as_posix()],
    )

    assert result.exit_code != 0
    assert "are not compatible" in result.output


def test_check_m2lnl_missing_m2lnl_column_raises(
    mcmc_catalog, monkeypatch: pytest.MonkeyPatch
):
    """check-m2lnl requires the catalog to have an M2LNL column."""
    experiment, catalog = mcmc_catalog

    monkeypatch.setattr(Ncm.MSetCatalog, "col_by_name", lambda self, name: (False, 0))

    result = runner.invoke(
        app,
        ["catalog", "check-m2lnl", experiment.as_posix(), catalog.as_posix()],
    )

    assert result.exit_code != 0
    assert "has no" in result.output


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
