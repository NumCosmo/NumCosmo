#!/usr/bin/env python
#
# test_cluster_richness_count_app.py
#
# Fri Jul 10 00:00:00 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_cluster_richness_count_app.py
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
# with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Test NumCosmo app's cluster mass-richness-count generator."""

from pathlib import Path

import numpy as np
import pytest
from typer.testing import CliRunner

from numcosmo_py import Ncm
from numcosmo_py.app import app

pytestmark = pytest.mark.app

runner = CliRunner()

Ncm.cfg_init()


@pytest.fixture(name="dc2_like_file")
def fixture_dc2_like_file(tmp_path: Path) -> Path:
    """Create a small FITS file mimicking a DC2-style true-count catalog."""
    pytest.importorskip("astropy")
    # flake8: noqa: E402
    # pylint: disable=import-outside-toplevel
    from astropy.table import Table

    rng = np.random.default_rng(7)
    n = 200
    halo_mass = 10 ** rng.uniform(13.0, 15.0, n)
    redshift = rng.uniform(0.1, 0.9, n)
    lnM = np.log(halo_mass)
    mu = 4.0 + 1.0 * (lnM - np.log(3.0e14)) + 0.2 * np.log1p(redshift)
    sigma = np.maximum(
        0.5 + 0.03 * (lnM - np.log(3.0e14)) + 0.15 * np.log1p(redshift), 0.05
    )
    lam = np.exp(rng.normal(mu, sigma))
    richness = rng.poisson(lam)

    data_file = tmp_path / "dc2_like.fits"
    Table({"halo_mass": halo_mass, "redshift": redshift, "richness": richness}).write(
        data_file
    )
    return data_file


@pytest.fixture(name="experiment_file")
def fixture_experiment_file(tmp_path: Path) -> Path:
    """Fixture for the experiment file."""
    return tmp_path / "experiment.yaml"


def test_generate_cluster_richness_count_default(experiment_file: Path):
    """Test the default generation of the cluster-richness-count experiment."""
    result = runner.invoke(
        app,
        [
            "generate",
            "cluster-richness-count",
            experiment_file.as_posix(),
            "--seed",
            "42",
            "--n-clusters",
            "300",
        ],
    )

    assert result.exit_code == 0, result.output
    assert experiment_file.exists()


def test_generate_cluster_richness_count_summary(experiment_file: Path):
    """Test the --summary output of the cluster-richness-count generator."""
    result = runner.invoke(
        app,
        [
            "generate",
            "cluster-richness-count",
            experiment_file.as_posix(),
            "--seed",
            "42",
            "--n-clusters",
            "300",
            "--summary",
        ],
    )

    assert result.exit_code == 0, result.output
    assert "Cluster Mass-Richness Count Mock" in result.output


def test_generate_cluster_richness_count_bad_suffix(tmp_path: Path):
    """Test that a non-.yaml suffix is rejected."""
    bad_file = tmp_path / "experiment.txt"
    result = runner.invoke(
        app, ["generate", "cluster-richness-count", bad_file.as_posix()]
    )

    assert result.exit_code != 0
    assert not bad_file.exists()


def test_generate_cluster_richness_count_bad_parameter(experiment_file: Path):
    """Test that an invalid fit parameter is rejected."""
    result = runner.invoke(
        app,
        [
            "generate",
            "cluster-richness-count",
            experiment_file.as_posix(),
            "--parameter-list",
            "NcClusterMass:bogus",
        ],
    )

    assert result.exit_code != 0


def test_run_test_cluster_richness_count(experiment_file: Path):
    """Test that the generated experiment can be loaded and evaluated once."""
    gen_result = runner.invoke(
        app,
        [
            "generate",
            "cluster-richness-count",
            experiment_file.as_posix(),
            "--seed",
            "42",
            "--n-clusters",
            "300",
        ],
    )
    assert gen_result.exit_code == 0, gen_result.output

    result = runner.invoke(app, ["run", "test", experiment_file.as_posix()])
    assert result.exit_code == 0, result.output
    assert "NcDataClusterMassRichCount" in result.output


def test_run_fit_cluster_richness_count(experiment_file: Path):
    """Test that the generated experiment can be fit through the CLI."""
    gen_result = runner.invoke(
        app,
        [
            "generate",
            "cluster-richness-count",
            experiment_file.as_posix(),
            "--seed",
            "42",
            "--n-clusters",
            "1000",
        ],
    )
    assert gen_result.exit_code == 0, gen_result.output

    result = runner.invoke(app, ["run", "fit", experiment_file.as_posix()])
    assert result.exit_code == 0, result.output


def test_generate_cluster_richness_count_from_file(
    dc2_like_file: Path, experiment_file: Path
):
    """Test loading real (DC2-like) data via --data-file."""
    result = runner.invoke(
        app,
        [
            "generate",
            "cluster-richness-count",
            experiment_file.as_posix(),
            "--data-file",
            dc2_like_file.as_posix(),
            "--mass-column",
            "halo_mass",
            "--redshift-column",
            "redshift",
            "--n-gal-column",
            "richness",
            "--summary",
        ],
    )

    assert result.exit_code == 0, result.output
    assert experiment_file.exists()
    assert "Cluster Mass-Richness Count Data" in result.output
    assert "Clusters loaded" in result.output


def test_generate_cluster_richness_count_from_file_bad_column(
    dc2_like_file: Path, experiment_file: Path
):
    """Test that a missing data column is reported clearly."""
    result = runner.invoke(
        app,
        [
            "generate",
            "cluster-richness-count",
            experiment_file.as_posix(),
            "--data-file",
            dc2_like_file.as_posix(),
            "--n-gal-column",
            "bogus_column",
        ],
    )

    assert result.exit_code != 0
    assert not experiment_file.exists()


def test_run_test_cluster_richness_count_from_file(
    dc2_like_file: Path, experiment_file: Path
):
    """Test that data loaded from a file can be evaluated through run test."""
    gen_result = runner.invoke(
        app,
        [
            "generate",
            "cluster-richness-count",
            experiment_file.as_posix(),
            "--data-file",
            dc2_like_file.as_posix(),
            "--mass-column",
            "halo_mass",
            "--redshift-column",
            "redshift",
            "--n-gal-column",
            "richness",
        ],
    )
    assert gen_result.exit_code == 0, gen_result.output

    result = runner.invoke(app, ["run", "test", experiment_file.as_posix()])
    assert result.exit_code == 0, result.output
    assert "NcDataClusterMassRichCount" in result.output
