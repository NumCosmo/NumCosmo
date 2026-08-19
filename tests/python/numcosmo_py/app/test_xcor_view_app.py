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


def test_view_kernel_cls_single() -> None:
    """--cls computes the auto-spectrum of a single kernel."""
    result = _view("--cls")

    assert result.exit_code == 0, result.output
    assert "Computing C_ell (non-Limber)" in result.output
    assert "1 kernel(s), 1 spectra" in result.output


def test_view_kernel_cls_pairs() -> None:
    """Every auto- and cross-pair of the requested kernels is computed."""
    result = runner.invoke(
        app,
        [
            "xcor",
            "kernel",
            "view",
            "--kernel",
            KERNEL,
            "--kernel",
            "number-counts survey=LSST-Y1 bin_idx=0 bias=1.5",
            "--ell",
            "2",
            "--n-ell",
            "2",
            "--no-show-plot",
            "--cls",
        ],
    )

    assert result.exit_code == 0, result.output
    assert "2 kernel(s), 3 spectra" in result.output


def test_view_kernel_cls_compare_limber() -> None:
    """--compare-limber adds a second, Limber, C_ell computation.

    The kernels are left in Limber mode by the k-space evaluation step, so this
    also covers that the Limber mode is set explicitly before each solve rather
    than inherited from whatever ran last.
    """
    result = _view("--cls", "--compare-limber")

    assert result.exit_code == 0, result.output
    assert "Computing C_ell (non-Limber)" in result.output
    assert "Computing C_ell (Limber)" in result.output


def test_view_kernel_cls_fixed_method() -> None:
    """The fixed-knot quadrature is reachable for the C_ell computation."""
    result = _view("--cls", "--cls-method", "fixed")

    assert result.exit_code == 0, result.output
    assert "method=fixed" in result.output


def test_integrator_tolerance_must_be_set_at_construction() -> None:
    """Guards the construction pattern the view command depends on.

    An ODE operator keeps the tolerance in force when it was built, so setting a
    tolerance on an already-constructed integrator leaves the result at the
    quality it was constructed with. The view command therefore builds its
    integrator with ``new_full``. If that regressed to ``new()`` followed by
    ``set_reltol()``, ``--integrator-reltol`` would silently do nothing -- which
    the reported-value assertions above cannot catch, since the getter agrees
    either way.
    """
    args = (5.0, 1.0, 0.1, 10.0, 50.0, 2)

    def build(tol: float) -> float:
        sbi = Ncm.SBesselIntegratorLevin.new_full(
            2, 2, 1e-4, 1e6, 21, 1200, tol, 2, tol
        )
        return sbi.integrate_gaussian_ell(*args)

    loose, tight = build(1e-4), build(1e-12)

    # Construction-time tolerance reaches the computation, and matters a lot here.
    assert abs(loose / tight - 1.0) > 1.0

    # Setting a tighter tolerance afterwards does not buy the tighter answer:
    # the result stays with the tolerance the operators were built at.
    relaxed = Ncm.SBesselIntegratorLevin.new_full(
        2, 2, 1e-4, 1e6, 21, 1200, 1e-4, 2, 1e-4
    )
    relaxed.set_reltol(1e-12)
    relaxed.set_cheb_reltol(1e-12)
    after = relaxed.integrate_gaussian_ell(*args)

    assert relaxed.get_reltol() == 1e-12  # the getter reports the new value ...
    assert abs(after / loose - 1.0) < 1e-6  # ... but the answer is the loose one
    assert abs(after / tight - 1.0) > 1.0
