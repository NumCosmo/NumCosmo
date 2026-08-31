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

import numpy as np

from numcosmo_py import Ncm
from numcosmo_py.app import app
from numcosmo_py.app.xcor import view

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


@pytest.mark.parametrize("option", ["--integrator-reltol", "--integrator-cheb-reltol"])
@pytest.mark.parametrize("value", ["0", "-1e-8", "2.0"])
def test_view_kernel_rejects_out_of_range_tolerance(option: str, value: str) -> None:
    """Neither tolerance is usable outside (0, 1]; the library would abort."""
    result = _view("--l-limber", "-1", option, value)

    assert result.exit_code != 0
    assert "Kernel evaluation complete" not in result.output


def test_view_kernel_defaults_leave_library_precision_alone() -> None:
    """Without the options the integrator keeps its own defaults."""
    result = _view("--l-limber", "0")

    assert result.exit_code == 0, result.output
    assert "reltol=1.0e-13" in result.output
    assert "cheb_reltol=1.0e-08" in result.output


@pytest.mark.parametrize("closure_type", ["spline", "chebyshev"])
def test_view_kernel_plots_either_closure(closure_type: str) -> None:
    """Both representations are reachable from the command line.

    The viewer holds no NcXcor when it builds a closure for the plot, so the
    representation has to come from its own option rather than from the
    computation -- and the C_ell path has to pass the same one on to NcXcor,
    or a plot and the spectrum beside it would be showing different fits.
    """
    result = _view("--l-limber", "-1", "--closure-type", closure_type, "--cls")

    assert result.exit_code == 0, result.output
    assert "Kernel evaluation complete" in result.output
    assert f"closure={closure_type}" in result.output


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


def test_view_kernel_compare_closure() -> None:
    """--compare-closure draws the other representation beside the chosen one.

    The two closures fit the same sampled window over the same domain, so the
    comparison is of the fits and not of the physics -- both C_ell runs stay on
    the same Limber tier, and only the representation differs between them.
    """
    result = _view(
        "--l-limber", "-1", "--closure-type", "spline", "--compare-closure", "--cls"
    )

    assert result.exit_code == 0, result.output
    assert "Chebyshev kernel k range" in result.output
    assert "closure=spline" in result.output
    assert "closure=chebyshev" in result.output
    # Both spectra are non-Limber: the tier is not what is being compared.
    assert "Computing C_ell (Limber)" not in result.output


def test_view_kernel_rejects_both_comparisons() -> None:
    """There is one alternative curve, so only one thing can occupy it."""
    result = _view("--compare-limber", "--compare-closure")

    assert result.exit_code != 0
    assert "Kernel evaluation complete" not in result.output


def test_view_kernel_cls_fixed_method() -> None:
    """The fixed-knot quadrature is reachable for the C_ell computation."""
    result = _view("--cls", "--cls-method", "fixed")

    assert result.exit_code == 0, result.output
    assert "method=fixed" in result.output


def test_integrator_tolerance_setting() -> None:
    """Test the setting of integrator tolerances."""
    args = (5.0, 1.0, 0.1, 10.0, 50.0, 2)

    def build(tol: float) -> float:
        sbi = Ncm.SBesselIntegratorLevin.new_full(
            2, 2, 1e-4, 1e6, 21, 1200, tol, 2, tol
        )
        return sbi.integrate_gaussian_ell(*args)

    loose, tight = build(1e-4), build(1e-12)
    assert abs(loose / tight - 1.0) > 1e-3

    sbi = Ncm.SBesselIntegratorLevin.new_full(2, 2, 1e-4, 1e6, 21, 1200, 1e-3, 2, 1e-3)
    sbi.set_reltol(1e-4)
    sbi.set_cheb_reltol(1e-4)
    assert sbi.get_reltol() == 1e-4
    assert sbi.get_cheb_reltol() == 1e-4
    loose_after = sbi.integrate_gaussian_ell(*args)

    sbi.set_reltol(1e-12)
    sbi.set_cheb_reltol(1e-12)
    assert sbi.get_reltol() == 1e-12
    assert sbi.get_cheb_reltol() == 1e-12
    tight_after = sbi.integrate_gaussian_ell(*args)

    assert abs(loose_after / loose - 1.0) < 1e-11
    assert abs(tight_after / tight - 1.0) < 1e-11


def test_view_kernel_integrator_reltol_reaches_the_computation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Test ``--integrator-reltol`` option."""
    captured: list[np.ndarray] = []
    evaluate = view.KernelEvaluation.evaluate

    def spy(self: view.KernelEvaluation, k_Mpc: np.ndarray) -> np.ndarray:
        values = evaluate(self, k_Mpc)
        captured.append(np.asarray(values, dtype=float).copy())
        return values

    monkeypatch.setattr(view.KernelEvaluation, "evaluate", spy)

    # 1e-8 rather than something tighter: the Levin RHS Chebyshev fit uses the
    # relative criterion alone, and a request below the integrator's accuracy
    # floor cannot converge, which is fatal at max-order.
    for reltol in ("1e-4", "1e-8"):
        result = _view(
            "--integrator-reltol", reltol, "--integrator-cheb-reltol", reltol
        )
        assert result.exit_code == 0, result.output

    assert len(captured) == 2
    loose, tight = captured[0], captured[1]
    assert loose.shape == tight.shape

    # Built at the requested tolerances the two runs differ by ~5e-5 here. Were
    # the tolerance applied after construction instead, both runs would compute
    # at the library default and agree to ~1e-9.
    assert np.abs(loose / tight - 1.0).max() > 1e-6
