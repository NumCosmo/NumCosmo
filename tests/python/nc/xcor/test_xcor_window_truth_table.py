#!/usr/bin/env python
#
# test_xcor_window_truth_table.py
#
# Wed Aug 27 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_xcor_window_truth_table.py
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

"""Check the radial integral against certified Arb values.

The radial integral

.. math:: I_\\ell(k) = \\int W(\\chi)\\, j_\\ell(k\\chi)\\, \\mathrm{d}\\chi

is where three quarters of an ``NcXcorSolver`` run is spent, and it has no
independent check inside NumCosmo -- every route to it goes through the same
Levin machinery. ``data/truth_tables/xcor/xcor_window_ilk.json.gz`` holds it
computed in Arb instead, certified to a relative ball radius of 2e-26, for
every analytic window on a grid that follows the peak across multipoles. This
test is the guard to run before and after touching the sbessel integrators.

Regenerating the table needs FLINT; running this test does not. See
``tests/tools/make_xcor_window_truth_table.py``.
"""

import gzip
import json
import pathlib

import numpy as np
import pytest
from numpy.testing import assert_allclose

from numcosmo_py import Nc, Ncm

Ncm.cfg_init()

# Pinned to one worker under --dist loadgroup: this file is one of the
# xcor lane's memory peaks, and an xdist worker is its own session, so
# without this its cost is paid once per worker rather than once.
pytestmark = [pytest.mark.xcor, pytest.mark.xdist_group("window_truth")]

TRUTH_TABLE = "truth_tables/xcor/xcor_window_ilk.json.gz"

# Measured, not chosen: over all 252 entries the worst relative deviation where
# |I| is within three orders of its own peak is 6.9e-10, and the worst
# deviation measured against the peak is 5.7e-11. Both tolerances sit about a
# factor of 15 above that. The integrator's nominal reltol of 1e-13 is *not*
# the right anchor -- it does not reach it, and cannot: the documented
# conditioning floor of these integrands is around 2e-10.
RTOL = 1.0e-8

# Applied as a fraction of each multipole's own peak. Far below the peak a
# relative criterion is meaningless -- the grid reaches |I|/peak = 1e-70 -- so
# this is what those entries are actually held to.
ATOL_FRAC = 1.0e-11

# The production rule, from nc_xcor_kernel.c: the caller knows the scale the
# result feeds into, and without an absolute floor the Chebyshev fit refuses to
# converge on the deep-tail entries and aborts on max-order.
INTEG_ABSTOL_FRAC = 1.0e-16


def _make_kernel(shape, dist, ps, sbi, ctor):
    """Build one analytic kernel from the constructor arguments in the table."""
    match shape:
        case "gauss":
            return Nc.XcorKernelAnalyticGauss.new_full(dist, ps, *ctor, sbi)
        case "tophat":
            return Nc.XcorKernelAnalyticTophat.new_full(dist, ps, *ctor, sbi)
        case "tophat_smooth":
            return Nc.XcorKernelAnalyticTophatSmooth.new_full(dist, ps, *ctor, sbi)
        case "student_t":
            return Nc.XcorKernelAnalyticStudentT.new_full(dist, ps, *ctor, sbi)
        case "power_exp":
            return Nc.XcorKernelAnalyticPowerExp.new_full(dist, ps, *ctor, sbi)
        case "lensing":
            return Nc.XcorKernelAnalyticLensing.new_full(dist, ps, *ctor, sbi)
        case "multi":
            mu, sigma, weight, n_sigma = ctor

            return Nc.XcorKernelAnalyticMulti.new_full(
                dist,
                ps,
                Ncm.Vector.new_array(mu),
                Ncm.Vector.new_array(sigma),
                Ncm.Vector.new_array(weight),
                n_sigma,
                sbi,
            )
        case _:
            raise ValueError(f"truth table names an unknown shape {shape!r}")


@pytest.fixture(name="truth_table", scope="module")
def fixture_truth_table() -> dict:
    """Load the certified table once for the module."""
    path = pathlib.Path(Ncm.cfg_get_data_filename(TRUTH_TABLE, True))

    with gzip.open(path, "rt") as f:
        return json.load(f)


@pytest.fixture(name="cosmo_bits", scope="module")
def fixture_cosmo_bits() -> tuple:
    """A prepared distance and a power spectrum, both constructor filler.

    The analytic windows are closed forms in chi alone: neither object reaches
    the value this test compares. They are built here rather than through
    ``Cosmology.default()`` because that imports the FFTW wisdom file, which
    costs more than every assertion in this file put together.
    """
    cosmo = Nc.HICosmoDEXcdm.new()
    dist = Nc.Distance.new(5.0)
    dist.prepare(cosmo)
    ps = Ncm.PowspecAnalytic.new(
        Ncm.PowspecAnalyticShape.BBKS, Ncm.PowspecAnalyticGrowth.LCDM
    )

    return cosmo, dist, ps


def test_table_is_certified_far_below_the_tolerance(truth_table: dict) -> None:
    """The reference must be a reference: its own uncertainty cannot matter.

    Guards a regeneration that silently lowered the precision target -- the
    comparison would still pass while no longer checking anything.
    """
    worst = max(
        radius / abs(float(value))
        for shape in truth_table["shapes"].values()
        for values, radii in zip(shape["table"], shape["radius"])
        for value, radius in zip(values, radii)
        if float(value) != 0.0
    )

    assert worst < 1.0e-20


@pytest.mark.parametrize(
    "shape",
    ["gauss", "tophat", "tophat_smooth", "student_t", "power_exp", "lensing", "multi"],
)
def test_radial_integral_matches_arb(
    shape: str, truth_table: dict, cosmo_bits: tuple
) -> None:
    """Compare I_ell(k) against Arb for one window, over every ell and k."""
    cosmo, dist, ps = cosmo_bits
    entry = truth_table["shapes"][shape]

    kernel = _make_kernel(
        shape, dist, ps, Ncm.SBesselIntegratorLevin.new(0, 8), entry["ctor"]
    )
    kernel.set_l_limber(-1)
    kernel.prepare(cosmo)

    supports = [kernel.get_comp_support(comp) for comp in range(kernel.get_n_comps())]

    for index, ell in enumerate(truth_table["ells"]):
        expected = np.array([float(value) for value in entry["table"][index]])
        peak = np.abs(expected).max()

        integrator = Ncm.SBesselIntegratorLevin.new(ell, ell)
        integrator.set_abstol(INTEG_ABSTOL_FRAC * peak)

        got = np.array(
            [
                sum(
                    integrator.integrate_ell(
                        # The window is a hard zero outside its support and the
                        # k*chi round trip puts nodes a hair beyond the edge, a
                        # cliff the Chebyshev fit cannot resolve. The library's
                        # own path clamps for this reason; so does this one.
                        lambda _ud, chi, _k, c=comp, lo=low, hi=high: (
                            kernel.eval_W_comp(c, min(max(chi, lo), hi))
                        ),
                        low,
                        high,
                        k,
                        ell,
                        None,
                    )
                    for comp, (low, high) in enumerate(supports)
                )
                for k in entry["kvals"][index]
            ]
        )

        assert_allclose(
            got,
            expected,
            rtol=RTOL,
            atol=ATOL_FRAC * peak,
            err_msg=f"{shape} at ell = {ell}",
        )


# Worst deviation measured per shape, closure at reltol = scaled-abstol = 1e-6,
# with roughly a factor of three of headroom. Re-measure these after changing
# NC_XCOR_KERNEL_CHEB_PANEL_K_CAP: a smaller cap makes more, lower-order panels,
# which converge to the requested tolerance by a different route.
#
# Read the spread against the spline's, not on its own -- against Arb the spline
# reaches 2.5e-4, 6.2e-5, 7.4e-4, 8.6e-3, 3.1e-5, 1.2e-4, 1.5e-5 on these same
# shapes.
# Where a number here is large it is the *sampling* that binds, not the fit:
# per multipole the closure sits on the sampling floor wherever the convergence
# criterion lets it reach, which a worst-over-ell figure like this cannot show.
CLOSURE_TOL = {
    "gauss": 9.0e-4,
    "tophat": 7.0e-7,
    "student_t": 2.0e-8,
    "power_exp": 2.0e-6,
    "lensing": 7.0e-6,
    "multi": 4.0e-4,
    "tophat_smooth": 4.0e-8,
}


@pytest.mark.parametrize("shape", sorted(CLOSURE_TOL))
def test_chebyshev_closure_matches_arb(
    shape: str, truth_table: dict, cosmo_bits: tuple
) -> None:
    """The closure itself against Arb, not just the radial integral under it.

    A closure holds C * sqrt(P(k)) * I_ell(k) for a constant C, so dividing the
    certified I_ell and sqrt(P) out has to leave a constant. How far that ratio
    wanders across k is the closure's own fitting error, measured against values
    with proven radii rather than against another of NumCosmo's paths.

    """
    cosmo, dist, ps = cosmo_bits
    entry = truth_table["shapes"][shape]
    RH = Nc.HICosmo.RH_Mpc(cosmo)

    worst = 0.0
    compared = 0

    for index, ell in enumerate(truth_table["ells"]):
        kernel = _make_kernel(
            shape, dist, ps, Ncm.SBesselIntegratorLevin.new(0, 8), entry["ctor"]
        )
        kernel.set_l_limber(-1)
        kernel.set_property("reltol", 1.0e-6)
        kernel.set_property("scaled-abstol", 1.0e-6)
        kernel.prepare(cosmo)

        integrand = kernel.get_eval_vectorized_full(
            cosmo,
            ell,
            ell,
            Ncm.SBesselIntegratorLevin.new(ell, ell),
            Nc.XcorKernelClosure.CHEBYSHEV,
        )
        k_min, k_max = integrand.get_range()

        expected = np.array([float(value) for value in entry["table"][index]])
        peak = np.abs(expected).max()

        ratios = [
            integrand.eval_array(k * RH)[0] / (value * np.sqrt(ps.eval(cosmo, 0.0, k)))
            for k, value in zip(entry["kvals"][index], expected)
            # Outside the fitted domain the closure extrapolates, and where
            # I_ell is a thousandth of its peak the ratio is dominated by the
            # certified value's own smallness rather than by the fit.
            if k_min < k * RH < k_max and abs(value) > 1.0e-3 * peak
        ]

        if len(ratios) < 3:
            continue

        ratios = np.array(ratios)
        worst = max(worst, np.abs(ratios / np.median(ratios) - 1.0).max())
        compared += len(ratios)

    assert compared > 0, f"{shape} shared no k with the closure's fitted range"
    assert worst < CLOSURE_TOL[shape], (
        f"{shape}: closure wanders by {worst:.3e} against Arb over "
        f"{compared} points, above the measured {CLOSURE_TOL[shape]:.1e}"
    )
