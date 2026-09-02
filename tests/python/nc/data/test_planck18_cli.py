#
# test_planck18_cli.py
#
# Mon September 01 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_planck18_cli.py
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

"""Generating and evaluating a Planck 2018 experiment through the CLI.

The subject is #NcDataPlanckLKL -- construction, properties, filename resolution, the
Boltzmann hookup and the likelihood evaluation -- which nothing else reaches; the CLI is
only the driver. Split out of test_numcosmo_app.py, which is marked ``app`` as a whole so
the default runs skip its parallel MCMC tests, leaving that surface uninstrumented. These
need no Planck data files and take a few seconds.

It lives here rather than beside the other CLI tests because conftest gates on
``"app" in item.keywords``, and pytest puts a test's directory names in its keywords: any
file under tests/python/numcosmo_py/app/ is gated whether or not it carries the marker.
"""

from pathlib import Path

import pytest
from typer.testing import CliRunner

pytest.importorskip("astropy")
pytest.importorskip("getdist")
# flake8: noqa: E402
# pylint: disable=wrong-import-position

from numcosmo_py.app import app
from numcosmo_py.app.generate import Planck18Types

runner = CliRunner()


@pytest.fixture(name="planck18_type", params=[e.value for e in Planck18Types])
def fixture_planck18_type(request) -> str:
    """Return planck18 type."""
    return request.param


def test_generate_planck(tmp_path: Path, planck18_type):
    """Generate a Planck 2018 experiment file."""
    tmp_file = tmp_path / "planck_generated1.yaml"

    result = runner.invoke(
        app,
        ["generate", "planck18", tmp_file.as_posix(), "--data-type", planck18_type],
    )

    if result.exit_code != 0:
        raise result.exception


def test_generate_planck_test(tmp_path: Path, planck18_type):
    """Generate a Planck 2018 experiment file and evaluate it."""
    tmp_file = tmp_path / "planck_generated2.yaml"

    result = runner.invoke(
        app,
        ["generate", "planck18", tmp_file.as_posix(), "--data-type", planck18_type],
    )

    if result.exit_code != 0:
        raise result.exception

    result = runner.invoke(app, ["run", "test", tmp_file.as_posix()])

    if result.exit_code != 0:
        raise result.exception
