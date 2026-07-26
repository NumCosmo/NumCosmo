#!/usr/bin/env python
#
# test_catalog_derived.py
#
# Sun Jul 26 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_catalog_derived.py
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

"""CLI-level tests for the `catalog derived-error` command and the
`catalog plot-corner --derived-*` options.

These exercise the full `--variable`/`--expr` pipeline (parameter
resolution, safe-expression evaluation, and posterior-statistic reporting)
through the actual Typer CLI, on a real (small) MCMC catalog.
"""

from typing import Tuple
from pathlib import Path

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
def fixture_mcmc_catalog(tmp_path) -> Tuple[Path, Path]:
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


def test_derived_error_median(mcmc_catalog):
    """catalog derived-error reports a median and 1/2/3-sigma bounds."""
    experiment, catalog = mcmc_catalog
    result = runner.invoke(
        app,
        [
            "catalog",
            "derived-error",
            experiment.as_posix(),
            catalog.as_posix(),
            "-x",
            "x=mu_0",
            "-x",
            "y=mu_1",
            "--expr",
            "x + y",
        ],
    )
    if result.exit_code != 0:
        raise result.exception
    assert "Derived quantity" in result.stdout
    assert "median" in result.stdout


def test_derived_error_pow10_transform(mcmc_catalog):
    """The motivating use case: reporting 10**log10(M)-style transforms."""
    experiment, catalog = mcmc_catalog
    result = runner.invoke(
        app,
        [
            "catalog",
            "derived-error",
            experiment.as_posix(),
            catalog.as_posix(),
            "-x",
            "x=mu_0",
            "--expr",
            "10**x",
            "--symbol",
            "10^{x}",
        ],
    )
    if result.exit_code != 0:
        raise result.exception
    assert "10^{x}" in result.stdout


def test_derived_error_all_stats(mcmc_catalog):
    """--stat may be repeated to report median, mode, and best-fit together."""
    experiment, catalog = mcmc_catalog
    result = runner.invoke(
        app,
        [
            "catalog",
            "derived-error",
            experiment.as_posix(),
            catalog.as_posix(),
            "-x",
            "x=mu_0",
            "--expr",
            "x",
            "--stat",
            "median",
            "--stat",
            "mode",
            "--stat",
            "bestfit",
        ],
    )
    if result.exit_code != 0:
        raise result.exception
    assert "median" in result.stdout
    assert "mode" in result.stdout
    assert "bestfit" in result.stdout


def test_derived_error_requires_variable(mcmc_catalog):
    """--expr without any --variable binding fails with a clear error."""
    experiment, catalog = mcmc_catalog
    result = runner.invoke(
        app,
        [
            "catalog",
            "derived-error",
            experiment.as_posix(),
            catalog.as_posix(),
            "--expr",
            "1",
        ],
    )
    assert result.exit_code != 0


def test_derived_error_unsafe_expr_rejected(mcmc_catalog):
    """An expression outside the safe-eval whitelist is rejected, not run."""
    experiment, catalog = mcmc_catalog
    result = runner.invoke(
        app,
        [
            "catalog",
            "derived-error",
            experiment.as_posix(),
            catalog.as_posix(),
            "-x",
            "x=mu_0",
            "--expr",
            "__import__('os').system('echo hi')",
        ],
    )
    assert result.exit_code != 0


def test_plot_corner_with_derived_dimension(mcmc_catalog, tmp_path):
    """plot-corner accepts an extra --derived-* dimension in the triangle plot."""
    experiment, catalog = mcmc_catalog
    output = tmp_path / "corner_out"
    result = runner.invoke(
        app,
        [
            "catalog",
            "plot-corner",
            experiment.as_posix(),
            catalog.as_posix(),
            "--no-show",
            "--output",
            output.as_posix(),
            "--derived-variable",
            "x=mu_0",
            "--derived-variable",
            "y=mu_1",
            "--derived-expr",
            "x + y",
            "--derived-symbol",
            "x+y",
            "--derived-name",
            "sum_xy",
        ],
    )
    if result.exit_code != 0:
        raise result.exception
    assert output.with_suffix(".corner.pdf").exists()


def test_plot_corner_derived_expr_requires_variable(mcmc_catalog, tmp_path):
    """--derived-expr without --derived-variable fails with a clear error."""
    experiment, catalog = mcmc_catalog
    output = tmp_path / "corner_out"
    result = runner.invoke(
        app,
        [
            "catalog",
            "plot-corner",
            experiment.as_posix(),
            catalog.as_posix(),
            "--no-show",
            "--output",
            output.as_posix(),
            "--derived-expr",
            "x",
        ],
    )
    assert result.exit_code != 0
