#!/usr/bin/env python
#
# test_xcor_view_app.py
#
# Fri Aug 08 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_xcor_view_app.py
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

"""Tests for the xcor kernel view CLI command."""

import pytest
from typer.testing import CliRunner

import matplotlib

matplotlib.use("Agg")

# flake8: noqa: E402
# pylint: disable=wrong-import-position

from numcosmo_py import Ncm
from numcosmo_py.app import app

pytestmark = pytest.mark.app
runner = CliRunner()

Ncm.cfg_init()

KERNEL = "weak-lensing survey=LSST-Y1 bin_idx=0"


def _view(*extra: str):
    """Run `xcor kernel view` on a single multipole, without displaying."""
    return runner.invoke(
        app,
        [
            "xcor",
            "kernel",
            "view",
            "--kernel",
            KERNEL,
            "--ell",
            "2",
            "--n-ell",
            "1",
            "--no-show-plot",
            *extra,
        ],
    )


@pytest.mark.parametrize("l_limber", ["-1", "0", "5"])
def test_view_kernel_l_limber_tiers(l_limber: str) -> None:
    """Each evaluation tier is reachable from the command line."""
    result = _view("--l-limber", l_limber)

    assert result.exit_code == 0, result.output
    assert "Kernel evaluation complete" in result.output


def test_view_kernel_integrator_precision_options() -> None:
    """The Levin precision knobs reach the integrator.

    Left unset each of these keeps the library default, so the test asserts
    the values it passed are the ones reported back rather than merely that
    the command succeeded.
    """
    result = _view(
        "--l-limber",
        "-1",
        "--integrator-reltol",
        "1e-2",
        "--integrator-cheb-reltol",
        "1e-2",
        "--integrator-max-order",
        "64",
    )

    assert result.exit_code == 0, result.output
    assert "reltol=1.0e-02" in result.output
    assert "cheb_reltol=1.0e-02" in result.output
    assert "max_order=64" in result.output


def test_view_kernel_defaults_leave_library_precision_alone() -> None:
    """Without the options the integrator keeps its own defaults."""
    result = _view("--l-limber", "0")

    assert result.exit_code == 0, result.output
    assert "reltol=1.0e-13" in result.output
    assert "cheb_reltol=1.0e-08" in result.output
