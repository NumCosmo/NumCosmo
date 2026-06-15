#!/usr/bin/env python
#
# test_inspect_app.py
#
# Thu May 14 2026
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

"""Tests for inspect CLI commands."""

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
from numcosmo_py.app import generate as gen
from numcosmo_py import Nc

pytestmark = pytest.mark.app
runner = CliRunner()

Ncm.cfg_init()


@pytest.fixture(name="simple_experiment")
def fixture_simple_experiment(tmp_path: Path) -> Tuple[Path, Ncm.ObjDictStr]:
    """Create a generic simple experiment with one Gaussian data object."""
    rng = Ncm.RNG.seeded_new(None, 1234)
    model_mvnd = Ncm.ModelMVND.new(5)
    mset = Ncm.MSet.new_array([model_mvnd])
    mset.param_set_all_ftype(Ncm.ParamType.FREE)

    data_mvnd = Ncm.DataGaussCovMVND.new_full(5, 1.0, 2.0, 30.0, 0.0, 0.0, rng)
    likelihood = Ncm.Likelihood.new(dset=Ncm.Dataset.new_array([data_mvnd]))

    exp_file = tmp_path / "simple_experiment.yaml"
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)

    experiment = Ncm.ObjDictStr.new()
    experiment.add("model-set", mset)
    experiment.add("likelihood", likelihood)

    ser.dict_str_to_yaml_file(experiment, exp_file.as_posix())

    return exp_file, experiment


@pytest.fixture(name="jpas_experiment")
def fixture_jpas_experiment(tmp_path: Path) -> Path:
    """Create a small JPAS cluster-count experiment for inspect tests."""
    exp_file = tmp_path / "jpas_small.yaml"

    _ = gen.GenerateJpasForecast(
        experiment=exp_file.absolute(),
        znknots=3,
        lnMobsnknots=2,
        fitting_sky_cut=gen.JpasSSCType.NO_SSC,
        resample_sky_cut=gen.JpasSSCType.NO_SSC,
    )

    return exp_file


def test_inspect_summary_lists_dataset_terms(
    simple_experiment: Tuple[Path, Ncm.ObjDictStr],
):
    """inspect summary should list likelihood terms from dataset."""
    exp_file, _ = simple_experiment

    result = runner.invoke(app, ["inspect", "summary", exp_file.as_posix()])

    assert result.exit_code == 0, result.output
    assert "Experiment summary" in result.output
    assert "Dataset likelihood terms" in result.output
    assert "DataGaussCovMVND" in result.output
    assert "Models present" in result.output


def test_inspect_summary_show_model_params(
    simple_experiment: Tuple[Path, Ncm.ObjDictStr],
):
    """inspect summary --show-model-params should print parameter table details."""
    exp_file, _ = simple_experiment

    result = runner.invoke(
        app,
        ["inspect", "summary", exp_file.as_posix(), "--show-model-params"],
    )

    assert result.exit_code == 0, result.output
    assert "Parameters:" in result.output
    assert "Fit type" in result.output
    assert "FREE" in result.output
    assert "ParamType." not in result.output


def test_inspect_cluster_ncounts_creates_plots(jpas_experiment: Path, tmp_path: Path):
    """inspect cluster-ncounts should generate plot files when output prefix is set."""
    out_prefix = tmp_path / "inspect_plots"

    result = runner.invoke(
        app,
        [
            "inspect",
            "cluster-ncounts",
            jpas_experiment.as_posix(),
            "--no-show-plot",
            "--output-prefix",
            out_prefix.as_posix(),
            "--no-show-sij",
        ],
    )

    assert result.exit_code == 0, result.output
    assert (tmp_path / "inspect_plots_data_vector.png").exists()
    assert (tmp_path / "inspect_plots_covariance.png").exists()


def test_inspect_cluster_ncounts_no_sij_message(jpas_experiment: Path):
    """inspect cluster-ncounts should report missing S_ij when none is configured."""
    result = runner.invoke(
        app,
        [
            "inspect",
            "cluster-ncounts",
            jpas_experiment.as_posix(),
            "--no-show-plot",
        ],
    )

    assert result.exit_code == 0, result.output
    assert "No S_ij matrices available to plot." in result.output


def test_inspect_cluster_ncounts_log_data_and_plot_options(
    jpas_experiment: Path, tmp_path: Path
):
    """inspect cluster-ncounts should accept plot style options and save figures."""
    out_prefix = tmp_path / "inspect_opts"

    result = runner.invoke(
        app,
        [
            "inspect",
            "cluster-ncounts",
            jpas_experiment.as_posix(),
            "--no-show-plot",
            "--output-prefix",
            out_prefix.as_posix(),
            "--log-data",
            "--cmap",
            "plasma",
            "--dpi",
            "120",
            "--no-show-sij",
        ],
    )

    assert result.exit_code == 0, result.output
    assert (tmp_path / "inspect_opts_data_vector.png").exists()
    assert (tmp_path / "inspect_opts_covariance.png").exists()


def test_inspect_cluster_ncounts_reshape_fallback(
    jpas_experiment: Path, monkeypatch: pytest.MonkeyPatch
):
    """inspect cluster-ncounts.

    inspect cluster-ncounts should fallback to flattened vector when bin shape
    mismatches.
    """

    def fake_get_z_obs(_self):
        # Force a shape mismatch with the original data-vector length.
        return Ncm.Vector.new_array([0.1, 0.3, 0.5, 0.7])

    monkeypatch.setattr(Nc.DataClusterNCountsGauss, "get_z_obs", fake_get_z_obs)

    result = runner.invoke(
        app,
        [
            "inspect",
            "cluster-ncounts",
            jpas_experiment.as_posix(),
            "--no-show-plot",
            "--no-show-sij",
        ],
    )

    assert result.exit_code == 0, result.output
    assert (
        "Could not infer (z, lnM) 2D grid from bins; plotting flattened vector."
        in result.output
    )


def test_inspect_cluster_ncounts_sij_plot_created(
    jpas_experiment: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    """inspect cluster-ncounts should save S_ij panel when S_ij matrices are present."""

    sij = Ncm.Matrix.new_array([1.0, 0.2, 0.2, 0.5], 2)

    monkeypatch.setattr(Nc.DataClusterNCountsGauss, "get_s_matrix", lambda _self: sij)
    monkeypatch.setattr(
        Nc.DataClusterNCountsGauss,
        "get_resample_s_matrix",
        lambda _self: sij,
    )

    out_prefix = tmp_path / "inspect_sij"
    result = runner.invoke(
        app,
        [
            "inspect",
            "cluster-ncounts",
            jpas_experiment.as_posix(),
            "--no-show-plot",
            "--output-prefix",
            out_prefix.as_posix(),
            "--show-sij",
        ],
    )

    assert result.exit_code == 0, result.output
    assert (tmp_path / "inspect_sij_sij.png").exists()


def test_inspect_cluster_ncounts_requires_cluster_data(
    simple_experiment: Tuple[Path, Ncm.ObjDictStr],
):
    """inspect cluster-ncounts should fail for experiments without cluster counts."""
    exp_file, _ = simple_experiment

    result = runner.invoke(
        app,
        ["inspect", "cluster-ncounts", exp_file.as_posix(), "--no-show-plot"],
    )

    assert result.exit_code != 0
    assert result.exception is not None
    assert "No Nc.DataClusterNCountsGauss data object found in dataset." in str(
        result.exception
    )
