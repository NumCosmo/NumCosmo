#!/usr/bin/env python
#
# test_hicosmo_de_wpoly.py
#
# Mon Sep 28 10:00:00 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_hicosmo_de_wpoly.py
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

"""Unit tests for NcHICosmoDEWPoly."""

import numpy as np
from numpy.polynomial import legendre as L
from numpy.testing import assert_allclose
import pytest

from numcosmo_py import Ncm, Nc

Ncm.cfg_init()

Z_F = 2.33


def set_coeffs(cosmo: Nc.HICosmoDEWPoly, coeffs: np.ndarray) -> None:
    """Set the Legendre coefficients of w + 1."""
    for i, c in enumerate(coeffs):
        cosmo.orig_vparam_set(Nc.HICosmoDEWPolyVParams.A, i, float(c))


def w_of_z(coeffs: np.ndarray, z: np.ndarray) -> np.ndarray:
    """w(z) of the model, frozen beyond Z_F."""
    zc = np.minimum(z, Z_F)
    return -1.0 + L.legval(2.0 * zc / Z_F - 1.0, coeffs)


def e2omega_de_ref(coeffs: np.ndarray, omega_x: float, z: float) -> float:
    """Reference dark energy density by a 128-node Gauss-Legendre rule."""
    zc = min(z, Z_F)
    x, wts = L.leggauss(128)
    t = 0.5 * zc * (x + 1.0)
    val = 0.5 * zc * np.sum(wts * (1.0 + w_of_z(coeffs, t)) / (1.0 + t))
    if z > Z_F:
        val += (1.0 + w_of_z(coeffs, np.array([Z_F]))[0]) * (
            np.log1p(z) - np.log1p(Z_F)
        )
    return omega_x * np.exp(3.0 * val)


@pytest.fixture(name="rng")
def fixture_rng() -> np.random.Generator:
    """Seeded generator."""
    return np.random.default_rng(1234)


@pytest.mark.parametrize("ncoeffs", [1, 2, 8, 51])
def test_construction(ncoeffs: int) -> None:
    """Defaults are LambdaCDM."""
    cosmo = Nc.HICosmoDEWPoly.new(ncoeffs=ncoeffs, z_f=Z_F)

    assert cosmo.get_ncoeffs() == ncoeffs
    assert cosmo.get_z_f() == Z_F
    assert cosmo.vparam_len(Nc.HICosmoDEWPolyVParams.A) == ncoeffs
    for z in (0.0, 0.5, Z_F, 5.0):
        assert cosmo.w_de(z) == -1.0
        assert cosmo.dw_dz(z) == 0.0


def test_lcdm_matches_xcdm() -> None:
    """All coefficients zero reproduces XCDM with w = -1."""
    cosmo = Nc.HICosmoDEWPoly.new(ncoeffs=8, z_f=Z_F)
    xcdm = Nc.HICosmoDEXcdm()
    xcdm.param_set_by_name("w", -1.0)

    for z in np.linspace(0.0, 5.0, 23):
        assert_allclose(cosmo.E2Omega_de(z), xcdm.E2Omega_de(z), rtol=1.0e-14)
        assert_allclose(cosmo.dE2Omega_de_dz(z), xcdm.dE2Omega_de_dz(z), rtol=1.0e-13)
        assert_allclose(
            cosmo.d2E2Omega_de_dz2(z), xcdm.d2E2Omega_de_dz2(z), rtol=1.0e-13
        )


def test_constant_w_matches_xcdm() -> None:
    """Only a_0 nonzero reproduces XCDM with w = -1 + a_0."""
    cosmo = Nc.HICosmoDEWPoly.new(ncoeffs=8, z_f=Z_F)
    set_coeffs(cosmo, np.array([0.3]))
    xcdm = Nc.HICosmoDEXcdm()
    xcdm.param_set_by_name("w", -0.7)

    for z in np.linspace(0.0, 5.0, 23):
        assert cosmo.w_de(z) == pytest.approx(-0.7)
        assert_allclose(cosmo.E2Omega_de(z), xcdm.E2Omega_de(z), rtol=1.0e-13)
        assert_allclose(cosmo.dE2Omega_de_dz(z), xcdm.dE2Omega_de_dz(z), rtol=1.0e-13)
        assert_allclose(
            cosmo.d2E2Omega_de_dz2(z), xcdm.d2E2Omega_de_dz2(z), rtol=1.0e-13
        )


@pytest.mark.parametrize("ncoeffs", [8, 51])
def test_w_and_derivative(rng: np.random.Generator, ncoeffs: int) -> None:
    """w(z) and dw/dz match the Legendre series, frozen beyond z_f."""
    coeffs = rng.standard_normal(ncoeffs) * 0.1
    cosmo = Nc.HICosmoDEWPoly.new(ncoeffs=ncoeffs, z_f=Z_F)
    set_coeffs(cosmo, coeffs)
    dcoeffs = L.legder(coeffs) * 2.0 / Z_F

    z = np.linspace(0.0, Z_F, 101)
    w_ref = w_of_z(coeffs, z)
    dw_ref = L.legval(2.0 * z / Z_F - 1.0, dcoeffs)
    assert_allclose([cosmo.w_de(zi) for zi in z], w_ref, rtol=1.0e-13, atol=1.0e-15)
    # dw/dz is zero from z_f on, so the polynomial derivative is compared below z_f
    assert_allclose(
        [cosmo.dw_dz(zi) for zi in z[:-1]], dw_ref[:-1], rtol=1.0e-12, atol=1.0e-14
    )

    for zi in (Z_F, 3.0, 10.0):
        assert cosmo.w_de(zi) == pytest.approx(w_ref[-1], rel=1.0e-14)
        assert cosmo.dw_dz(zi) == 0.0


@pytest.mark.parametrize("ncoeffs", [8, 51])
def test_e2omega_de(rng: np.random.Generator, ncoeffs: int) -> None:
    """E2Omega_de and its derivatives match adaptive quadrature at both degrees."""
    coeffs = rng.standard_normal(ncoeffs) * 0.1
    cosmo = Nc.HICosmoDEWPoly.new(ncoeffs=ncoeffs, z_f=Z_F)
    set_coeffs(cosmo, coeffs)
    omega_x = cosmo.param_get_by_name("Omegax")

    for z in np.concatenate((np.linspace(0.0, Z_F, 13), [3.0, 5.0])):
        e2 = e2omega_de_ref(coeffs, omega_x, z)
        w = w_of_z(coeffs, np.array([z]))[0]
        dw = cosmo.dw_dz(z)
        assert_allclose(cosmo.E2Omega_de(z), e2, rtol=1.0e-12)
        assert_allclose(
            cosmo.dE2Omega_de_dz(z), 3.0 * (1.0 + w) * e2 / (1.0 + z), rtol=1.0e-12
        )
        d2_ref = (
            3.0 * (2.0 + w * (5.0 + 3.0 * w) + dw * (1.0 + z)) * e2 / (1.0 + z) ** 2
        )
        assert_allclose(cosmo.d2E2Omega_de_dz2(z), d2_ref, rtol=1.0e-11)


def test_parameter_change_is_seen() -> None:
    """Changing a coefficient after an evaluation updates the frozen tail."""
    cosmo = Nc.HICosmoDEWPoly.new(ncoeffs=3, z_f=Z_F)
    assert cosmo.w_de(5.0) == -1.0
    set_coeffs(cosmo, np.array([0.0, 0.2]))
    assert cosmo.w_de(5.0) == pytest.approx(-0.8)
    assert cosmo.w_de(0.0) == pytest.approx(-1.2)


def test_serialization() -> None:
    """Serialization round trip keeps z_f and the coefficients."""
    cosmo = Nc.HICosmoDEWPoly.new(ncoeffs=5, z_f=1.7)
    set_coeffs(cosmo, np.array([0.1, -0.2, 0.05, 0.0, 0.3]))
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    cosmo2 = ser.dup_obj(cosmo)
    assert isinstance(cosmo2, Nc.HICosmoDEWPoly)
    assert cosmo2.get_z_f() == 1.7
    assert cosmo2.get_ncoeffs() == 5
    for z in np.linspace(0.0, 3.0, 7):
        assert cosmo2.w_de(z) == cosmo.w_de(z)
        assert cosmo2.E2Omega_de(z) == cosmo.E2Omega_de(z)


def coeffs_from_roots(
    roots_z: list[float], pairs: list[tuple[float, float]], sign: float = 1.0
) -> np.ndarray:
    """Legendre coefficients of w + 1 whose dw/dz has the given real roots.

    ``pairs`` are (b, s) complex-conjugate root pairs of dw/dz in z, which add
    no turning point.
    """
    P = np.polynomial.polynomial
    x_roots = [2.0 * z / Z_F - 1.0 for z in roots_z]
    p = np.array([sign])
    for xr in x_roots:
        p = P.polymul(p, np.array([-xr, 1.0]))
    for b, s in pairs:
        xb, xs = 2.0 * b / Z_F - 1.0, 2.0 * s / Z_F
        p = P.polymul(p, np.array([xb * xb + xs * xs, -2.0 * xb, 1.0]))
    dw_leg = L.poly2leg(p)
    return L.legint(dw_leg) * 0.5 * Z_F


@pytest.mark.parametrize(
    "roots_z,pairs",
    [
        ([], []),
        ([0.5], []),
        ([0.3, 1.2], []),
        ([0.2, 0.9, 2.0], []),
        ([0.5], [(1.5, 0.2)]),
        ([0.3, 1.2], [(0.8, 0.05), (2.0, 0.4)]),
        ([0.1, 0.15, 2.2], [(1.0, 0.01)]),
    ],
)
def test_turning_points_constructed(
    roots_z: list[float], pairs: list[tuple[float, float]]
) -> None:
    """Prescribed real roots of dw/dz are recovered; complex pairs add none."""
    coeffs = coeffs_from_roots(roots_z, pairs)
    cosmo = Nc.HICosmoDEWPoly.new(ncoeffs=len(coeffs), z_f=Z_F)
    set_coeffs(cosmo, coeffs)

    tp = cosmo.turning_points()
    assert len(tp) == len(roots_z)
    assert_allclose(tp, sorted(roots_z), rtol=1.0e-10, atol=1.0e-12)
    for z in tp:
        assert cosmo.dw_dz(z) == pytest.approx(0.0, abs=1.0e-10)


def test_turning_points_outside_range() -> None:
    """Roots of dw/dz outside (0, z_f) are not turning points."""
    coeffs = coeffs_from_roots([-0.5, 0.7, 3.0], [])
    cosmo = Nc.HICosmoDEWPoly.new(ncoeffs=len(coeffs), z_f=Z_F)
    set_coeffs(cosmo, coeffs)
    tp = cosmo.turning_points()
    assert len(tp) == 1
    assert tp[0] == pytest.approx(0.7)


def test_turning_points_degenerate() -> None:
    """LambdaCDM, constant and linear w have no turning points."""
    cosmo = Nc.HICosmoDEWPoly.new(ncoeffs=8, z_f=Z_F)
    assert len(cosmo.turning_points()) == 0
    set_coeffs(cosmo, np.array([0.3]))
    assert len(cosmo.turning_points()) == 0
    set_coeffs(cosmo, np.array([0.3, -0.2]))
    assert len(cosmo.turning_points()) == 0
    # w quadratic: dw/dz linear, one root
    set_coeffs(cosmo, np.array([0.0, 0.0, 0.1]))
    tp = cosmo.turning_points()
    assert len(tp) == 1
    assert tp[0] == pytest.approx(0.5 * Z_F)


@pytest.mark.parametrize("ncoeffs", [8, 51])
def test_turning_points_random(rng: np.random.Generator, ncoeffs: int) -> None:
    """Random series agree with numpy's Legendre root finder."""
    for _ in range(20):
        coeffs = rng.standard_normal(ncoeffs)
        cosmo = Nc.HICosmoDEWPoly.new(ncoeffs=ncoeffs, z_f=Z_F)
        set_coeffs(cosmo, coeffs)
        r = L.legroots(L.legder(coeffs))
        real = r[np.abs(r.imag) < 1.0e-8].real
        ref = np.sort(0.5 * Z_F * (real[(real > -1.0) & (real < 1.0)] + 1.0))
        tp = np.array(cosmo.turning_points())
        assert len(tp) == len(ref)
        assert_allclose(tp, ref, rtol=1.0e-8, atol=1.0e-10)
