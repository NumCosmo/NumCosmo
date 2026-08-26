#!/usr/bin/env python
#
# test_powspec_analytic.py
#
# Tue Aug 26 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_powspec_analytic.py
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

"""Tests for NcmPowspecAnalytic.

Every other #NcmPowspec is spline-backed, so its value has no independently
known answer. This one is defined by a formula, and the point of the tests is
that the formula is what the library evaluates: each shape is checked against
the closed form written out here from the same definition, and the derivatives
against finite differences of the library's own eval().

A certified check against Arb ball arithmetic lives in
`tests/tools/compare_powspec_analytic_arb.py`; it needs FLINT and is not part
of this suite.
"""

import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.special import hyp2f1

from numcosmo_py import Ncm

# Required for this file to run at all: conftest.py skips anything whose keywords
# contain "powspec" -- which the directory name alone supplies -- unless
# --run-powspec is given, while the py-powspec lane selects with `-m powspec`,
# which matches markers only.
pytestmark = pytest.mark.powspec

Ncm.cfg_init()

SHAPES = [
    Ncm.PowspecAnalyticShape.POWER_LAW,
    Ncm.PowspecAnalyticShape.BBKS,
    Ncm.PowspecAnalyticShape.RATIONAL,
]
GROWTHS = [
    Ncm.PowspecAnalyticGrowth.NONE,
    Ncm.PowspecAnalyticGrowth.LCDM,
    Ncm.PowspecAnalyticGrowth.RATIONAL,
]

K_GRID = np.logspace(-5.0, 1.0, 41)
Z_GRID = np.array([0.0, 0.1, 0.5, 1.0, 3.0, 10.0, 20.0])

BBKS_C = (2.34, 3.89, 16.1, 5.46, 6.71)


def _transfer(shape, k, k_eq):
    """T(k), written from the definition rather than read from the library."""
    q = k / k_eq

    if shape == Ncm.PowspecAnalyticShape.POWER_LAW:
        return np.ones_like(q)

    if shape == Ncm.PowspecAnalyticShape.BBKS:
        c1, c2, c3, c4, c5 = BBKS_C
        big_m = 1.0 + c2 * q + (c3 * q) ** 2 + (c4 * q) ** 3 + (c5 * q) ** 4
        return np.log1p(c1 * q) / (c1 * q) * big_m**-0.25

    return 1.0 / (1.0 + q**2)


def _growth(growth, z, omega_m, a_t):
    """D(z), normalized to D(0) = 1."""
    a = 1.0 / (1.0 + z)

    if growth == Ncm.PowspecAnalyticGrowth.NONE:
        return np.ones_like(np.asarray(a, dtype=float))

    if growth == Ncm.PowspecAnalyticGrowth.LCDM:

        def raw(scale):
            return scale * hyp2f1(
                1.0 / 3.0, 1.0, 11.0 / 6.0, -(1.0 - omega_m) / omega_m * scale**3
            )

        return raw(a) / raw(1.0)

    def raw_rational(scale):
        return scale * (1.0 + (scale / a_t) ** 3) ** (-1.0 / 3.0)

    return raw_rational(a) / raw_rational(1.0)


def _reference(psa, z, k):
    """A k^n_s T(k)^2 D(z)^2, from the definitions above."""
    return (
        psa.get_amplitude()
        * k ** psa.get_n_s()
        * _transfer(psa.get_shape(), k, psa.get_k_eq()) ** 2
        * _growth(psa.get_growth(), z, psa.get_Omega_m(), 1.0) ** 2
    )


@pytest.fixture(name="psa", params=SHAPES, ids=lambda s: s.value_nick)
def fixture_psa(request):
    """One prepared object per transfer-function shape, LCDM growth."""
    obj = Ncm.PowspecAnalytic.new(request.param, Ncm.PowspecAnalyticGrowth.LCDM)
    obj.prepare(None)

    return obj


@pytest.mark.parametrize("shape", SHAPES, ids=lambda s: s.value_nick)
@pytest.mark.parametrize("growth", GROWTHS, ids=lambda g: g.value_nick)
def test_eval_matches_closed_form(shape, growth):
    """eval() reproduces the formula the object documents."""
    psa = Ncm.PowspecAnalytic.new(shape, growth)
    psa.prepare(None)

    for z in Z_GRID:
        got = np.array([psa.eval(None, z, k) for k in K_GRID])
        assert_allclose(got, _reference(psa, z, K_GRID), rtol=1.0e-13)


def test_properties_round_trip():
    """new_full() stores what it is given, and the getters report it."""
    psa = Ncm.PowspecAnalytic.new_full(
        Ncm.PowspecAnalyticShape.BBKS,
        Ncm.PowspecAnalyticGrowth.LCDM,
        1.5e7,
        0.97,
        0.09,
        0.31,
    )

    assert psa.get_shape() == Ncm.PowspecAnalyticShape.BBKS
    assert psa.get_growth() == Ncm.PowspecAnalyticGrowth.LCDM
    assert psa.get_amplitude() == 1.5e7
    assert psa.get_n_s() == 0.97
    assert psa.get_k_eq() == 0.09
    assert psa.get_Omega_m() == 0.31


def test_eval_vec_matches_eval(psa):
    """The vectorized path agrees with the scalar one exactly."""
    kv = Ncm.Vector.new_array(list(K_GRID))
    pv = Ncm.Vector.new(len(K_GRID))
    psa.eval_vec(None, 0.3, kv, pv)

    got = np.array([pv.get(i) for i in range(len(K_GRID))])
    expected = np.array([psa.eval(None, 0.3, k) for k in K_GRID])

    assert_allclose(got, expected, rtol=0.0, atol=0.0)


def test_deriv_k_matches_finite_differences(psa):
    """dP/dk agrees with a central difference of eval()."""
    for k in (1.0e-4, 1.0e-2, 0.1, 1.0, 5.0):
        h = k * 1.0e-5
        fd = (psa.eval(None, 0.4, k + h) - psa.eval(None, 0.4, k - h)) / (2.0 * h)

        assert_allclose(psa.deriv_k(None, 0.4, k), fd, rtol=1.0e-7)


@pytest.mark.parametrize("growth", GROWTHS, ids=lambda g: g.value_nick)
def test_deriv_z_matches_finite_differences(growth):
    """dP/dz agrees with a central difference of eval(), every growth shape."""
    psa = Ncm.PowspecAnalytic.new(Ncm.PowspecAnalyticShape.BBKS, growth)
    psa.prepare(None)

    for z in (0.5, 2.0, 10.0):
        h = 1.0e-5
        fd = (psa.eval(None, z + h, 0.1) - psa.eval(None, z - h, 0.1)) / (2.0 * h)

        assert_allclose(psa.deriv_z(None, z, 0.1), fd, rtol=1.0e-7, atol=1.0e-30)


@pytest.mark.parametrize("growth", GROWTHS, ids=lambda g: g.value_nick)
def test_growth_normalized_at_zero(growth):
    """D(0) = 1 by construction, for every growth shape."""
    psa = Ncm.PowspecAnalytic.new(Ncm.PowspecAnalyticShape.BBKS, growth)
    psa.prepare(None)

    assert psa.eval_growth(0.0) == pytest.approx(1.0, rel=1.0e-15)


def test_growth_decreases_with_redshift(psa):
    """LCDM growth is monotonically decreasing in z."""
    values = np.array([psa.eval_growth(z) for z in Z_GRID])

    assert np.all(np.diff(values) < 0.0)


def test_bbks_asymptotic_slopes():
    """The BBKS shape has the two limits it exists to supply.

    Below the turnover P goes as k^n_s; far above it the transfer function
    contributes k^-2 ln k, so the logarithmic slope approaches n_s - 4 from
    above and has not reached it by k = 10. A pure power-law tail would sit
    at n_s - 4 exactly, which is the RATIONAL shape's contrast.
    """
    psa = Ncm.PowspecAnalytic.new(
        Ncm.PowspecAnalyticShape.BBKS, Ncm.PowspecAnalyticGrowth.NONE
    )
    psa.prepare(None)

    k = np.logspace(-8.0, 1.0, 200)
    slope = np.gradient(
        np.log(np.array([psa.eval(None, 0.0, ki) for ki in k])), np.log(k)
    )

    assert slope[2] == pytest.approx(psa.get_n_s(), abs=1.0e-4)
    assert psa.get_n_s() - 4.0 < slope[-3] < psa.get_n_s() - 2.0


def test_rational_tail_is_a_pure_power_law():
    """The RATIONAL shape reaches n_s - 4 exactly, unlike BBKS."""
    psa = Ncm.PowspecAnalytic.new(
        Ncm.PowspecAnalyticShape.RATIONAL, Ncm.PowspecAnalyticGrowth.NONE
    )
    psa.prepare(None)

    k = np.logspace(2.0, 4.0, 50)
    slope = np.gradient(
        np.log(np.array([psa.eval(None, 0.0, ki) for ki in k])), np.log(k)
    )

    assert slope[len(slope) // 2] == pytest.approx(psa.get_n_s() - 4.0, abs=1.0e-3)


def test_bao_factor_off_by_default(psa):
    """B(k) is identically one until bao-amplitude is set."""
    for k in K_GRID:
        assert psa.eval_bao(k) == 1.0


def test_bao_factor_oscillates_and_damps():
    """With the factor on, P is modulated and the modulation dies out."""
    psa = Ncm.PowspecAnalytic.new(
        Ncm.PowspecAnalyticShape.BBKS, Ncm.PowspecAnalyticGrowth.NONE
    )
    psa.set_property("bao-amplitude", 0.1)
    psa.set_property("bao-rd", 147.0)
    psa.set_property("bao-sigma", 10.0)
    psa.prepare(None)

    k = np.linspace(0.01, 0.2, 400)
    bao = np.array([psa.eval_bao(ki) for ki in k])

    assert np.any(bao > 1.0)
    assert np.any(bao < 1.0)
    # Damped away well before k = 1, where exp(-(k sigma)^2) is ~1e-44.
    assert psa.eval_bao(1.0) == pytest.approx(1.0, abs=1.0e-15)


def test_transfer_is_one_at_zero_k():
    """The removable 0/0 in the BBKS log term evaluates to its limit."""
    psa = Ncm.PowspecAnalytic.new(
        Ncm.PowspecAnalyticShape.BBKS, Ncm.PowspecAnalyticGrowth.NONE
    )
    psa.prepare(None)

    assert psa.eval_transfer(0.0) == 1.0
    assert psa.eval_transfer(1.0e-300) == pytest.approx(1.0, rel=1.0e-15)


def test_serialization_round_trip(psa):
    """The object survives a serialize/deserialize cycle unchanged."""
    ser = Ncm.Serialize.new(0)
    clone = ser.from_string(ser.to_string(psa, True))

    for z in (0.0, 1.0):
        for k in (1.0e-3, 0.1, 5.0):
            assert clone.eval(None, z, k) == psa.eval(None, z, k)


@pytest.mark.parametrize("prop", ["k-eq", "Omega-m", "a-t"])
def test_scale_parameters_are_bounded_away_from_zero(prop):
    """Zero is rejected by the property's own range, before any evaluation.

    These three parameters divide, so zero would be a division by zero rather
    than a merely odd shape. The GParamSpec minimum is what enforces it, which
    is why the object's constructed() carries no redundant check of its own.
    """
    psa = Ncm.PowspecAnalytic.new(
        Ncm.PowspecAnalyticShape.BBKS, Ncm.PowspecAnalyticGrowth.LCDM
    )

    with pytest.warns(Warning, match="out of range"):
        psa.set_property(prop, 0.0)

    assert psa.get_property(prop) > 0.0
