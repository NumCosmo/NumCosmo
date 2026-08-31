#!/usr/bin/env python
#
# test_kernel_table.py
#
# Sat August 30 2026
# Copyright  2026  Sandro Dias Pinto Vitenti <vitenti@uel.br>
#
# numcosmo is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# numcosmo is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Tests for NcXcorKernelTable, the tabulated radial window.

The anchor is a round trip: sample an analytic window onto a table, rebuild it
as a table kernel, and require the C_ell back. The analytic kernel is an
independently known answer, so this measures the reconstruction rather than
comparing the class to itself.
"""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from numcosmo_py import Nc, Ncm

Ncm.cfg_init()

pytestmark = [pytest.mark.xcor]

CHI_MEAN, CHI_SIGMA, N_SIGMA = 1500.0, 300.0, 4.0
CHI_LO, CHI_HI = CHI_MEAN - N_SIGMA * CHI_SIGMA, CHI_MEAN + N_SIGMA * CHI_SIGMA
ELLS = (2, 8, 32, 100)


@pytest.fixture(name="bits", scope="module")
def fixture_bits():
    """A cosmology, a distance and a closed-form power spectrum."""
    cosmo = Nc.HICosmoDEXcdm.new()
    dist = Nc.Distance.new(5.0)
    dist.prepare(cosmo)
    ps = Ncm.PowspecAnalytic.new(
        Ncm.PowspecAnalyticShape.BBKS, Ncm.PowspecAnalyticGrowth.LCDM
    )
    return cosmo, dist, ps


def _analytic(bits):
    cosmo, dist, ps = bits
    sbi = Ncm.SBesselIntegratorLevin.new(0, 128)
    k = Nc.XcorKernelAnalyticGauss.new_full(dist, ps, CHI_MEAN, CHI_SIGMA, N_SIGMA, sbi)
    k.set_l_limber(-1)
    k.prepare(cosmo)
    return k


def _samples(kernel, n=2000, lo=CHI_LO, hi=CHI_HI):
    chi = np.linspace(lo, hi, n)
    return chi, np.array([kernel.eval_W(c) for c in chi])


def _table(bits, chi, W, kind=None, order=8, normalize=True, with_sbi=True):
    _, dist, ps = bits
    kind = Nc.XcorKernelTableKind.DENSITY if kind is None else kind
    sbi = Ncm.SBesselIntegratorLevin.new(0, 128) if with_sbi else None
    return Nc.XcorKernelTable.new_full(
        dist,
        ps,
        Ncm.Vector.new_array(chi.tolist()),
        Ncm.Vector.new_array(W.tolist()),
        kind,
        order,
        normalize,
        sbi,
    )


def _cls(bits, kernel, ells=ELLS):
    cosmo, dist, ps = bits
    xc = Nc.Xcor.new(dist, ps, Nc.XcorMethod.KERNEL_EXACT)
    xc.prepare(cosmo)
    out = []
    for ell in ells:
        v = Ncm.Vector.new(1)
        xc.compute(kernel, kernel, cosmo, ell, ell, v)
        out.append(v.get(0))
    return np.array(out)


@pytest.mark.parametrize("order,tol", [(4, 1e-10), (6, 1e-11), (8, 1e-12)])
def test_round_trip_reproduces_the_analytic_kernel(bits, order, tol):
    """A sampled analytic window, rebuilt as a table, gives back its own C_ell."""
    cosmo = bits[0]
    ana = _analytic(bits)
    chi, W = _samples(ana)

    tab = _table(bits, chi, W, order=order)
    tab.set_l_limber(-1)
    tab.prepare(cosmo)

    assert_allclose(_cls(bits, tab), _cls(bits, ana), rtol=tol)


def test_higher_order_is_not_worse(bits):
    """Degree 7 beats a cubic on the same samples, which is why it is the default."""
    cosmo = bits[0]
    ana = _analytic(bits)
    chi, W = _samples(ana, n=200)
    ref = _cls(bits, ana)

    errs = []
    for order in (4, 8):
        tab = _table(bits, chi, W, order=order)
        tab.set_l_limber(-1)
        tab.prepare(cosmo)
        errs.append(np.max(np.abs(_cls(bits, tab) / ref - 1.0)))

    assert errs[1] < errs[0]


def test_window_matches_the_sampled_one(bits):
    """eval_W reproduces the samples it was built from."""
    ana = _analytic(bits)
    chi, W = _samples(ana)
    tab = _table(bits, chi, W)

    got = np.array([tab.eval_W(c) for c in chi])
    assert_allclose(got, W, rtol=1e-10, atol=1e-14)


def test_normalized_window_integrates_to_one(bits):
    """The window is rescaled to unit integral over its support."""
    ana = _analytic(bits)
    chi, W = _samples(ana)
    tab = _table(bits, chi, 3.7 * W)

    lo, hi = tab.get_support()
    grid = np.linspace(lo, hi, 20001)
    integ = np.trapezoid([tab.eval_W(c) for c in grid], grid)
    assert_allclose(integ, 1.0, rtol=1e-8)
    assert_allclose(tab.get_norm(), 3.7, rtol=1e-6)


def test_unnormalized_keeps_the_supplied_scale(bits):
    """With normalize off the samples are used as given."""
    ana = _analytic(bits)
    chi, W = _samples(ana)
    tab = _table(bits, chi, 3.7 * W, normalize=False)

    assert tab.get_norm() == 1.0
    i = len(chi) // 2
    assert_allclose(tab.eval_W(chi[i]), 3.7 * W[i], rtol=1e-10)


def test_zero_runs_are_trimmed_from_the_support(bits):
    """Padding the table with zeros does not widen the integration domain."""
    ana = _analytic(bits)
    chi, W = _samples(ana)

    pad_lo = np.linspace(50.0, CHI_LO, 200, endpoint=False)
    pad_hi = np.linspace(CHI_HI, 6000.0, 200)[1:]
    chi_p = np.concatenate([pad_lo, chi, pad_hi])
    W_p = np.concatenate([np.zeros_like(pad_lo), W, np.zeros_like(pad_hi)])

    tab = _table(bits, chi_p, W_p)
    lo, hi = tab.get_support()

    # One padding sample is kept on each side, no more.
    assert lo == pytest.approx(pad_lo[-1], rel=1e-12)
    assert hi == pytest.approx(pad_hi[0], rel=1e-12)


def test_knots_are_the_trimmed_abscissae(bits):
    """The reconstruction reports the breakpoints a quadrature should align to."""
    ana = _analytic(bits)
    chi, W = _samples(ana)
    tab = _table(bits, chi, W)

    knots = tab.peek_knots()
    assert knots.len() == len(chi)
    assert knots.get(0) == pytest.approx(chi[0])
    assert knots.get(knots.len() - 1) == pytest.approx(chi[-1])


def test_properties_round_trip(bits):
    """The accessors and the GObject properties agree; serialization reads both."""
    ana = _analytic(bits)
    chi, W = _samples(ana)
    tab = _table(bits, chi, W, kind=Nc.XcorKernelTableKind.SHEAR, order=6)

    assert tab.get_kind() == Nc.XcorKernelTableKind.SHEAR
    assert tab.get_order() == 6
    assert tab.get_normalize() is True
    assert tab.props.kind == Nc.XcorKernelTableKind.SHEAR
    assert tab.props.order == 6
    assert tab.props.normalize is True


def test_shear_suppresses_the_spectrum(bits):
    """Shear carries 1/(k chi)^2, so its C_ell is far below the density one."""
    cosmo = bits[0]
    ana = _analytic(bits)
    chi, W = _samples(ana)

    dens = _table(bits, chi, W)
    dens.set_l_limber(-1)
    dens.prepare(cosmo)

    shear = _table(bits, chi, W, kind=Nc.XcorKernelTableKind.SHEAR)
    shear.set_l_limber(-1)
    shear.prepare(cosmo)

    cl_d, cl_s = _cls(bits, dens), _cls(bits, shear)
    assert np.all(cl_s > 0.0)
    assert np.all(cl_s < cl_d)


def test_shear_prefactor(bits):
    """sqrt((l+2)(l+1)l(l-1)), and zero where that would be negative."""
    cosmo = bits[0]
    ana = _analytic(bits)
    chi, W = _samples(ana)

    shear = _table(bits, chi, W, kind=Nc.XcorKernelTableKind.SHEAR)
    dens = _table(bits, chi, W)

    for ell in (0, 1):
        assert shear.eval_prefactor(cosmo, ell) == 0.0

    for ell in (2, 8, 32, 100):
        expect = np.sqrt((ell + 2) * (ell + 1) * ell * (ell - 1))
        assert_allclose(shear.eval_prefactor(cosmo, ell), expect, rtol=1e-14)
        assert dens.eval_prefactor(cosmo, ell) == 1.0


def test_shear_kernel_factor(bits):
    """Shear carries 1/(k chi)^2 in the kernel; density carries nothing."""
    cosmo = bits[0]
    ana = _analytic(bits)
    chi, W = _samples(ana)

    shear = _table(bits, chi, W, kind=Nc.XcorKernelTableKind.SHEAR)
    dens = _table(bits, chi, W)

    for c, k in ((1500.0, 0.1), (900.0, 0.02)):
        assert_allclose(
            shear.eval_kernel_factor(cosmo, c, k), 1.0 / (k * c) ** 2, rtol=1e-14
        )
        assert dens.eval_kernel_factor(cosmo, c, k) == 1.0


# Construct-time invariants (mismatched lengths, too few samples for the requested
# degree) are reported with g_error, as every sibling kernel reports its own. That
# is fatal by design and so cannot be exercised from here.
