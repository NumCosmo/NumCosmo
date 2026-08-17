#!/usr/bin/env python
#
# test_multiplicity_func_castro.py
#
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

"""Tests for NcMultiplicityFuncCastro."""

import math

import numpy as np
import pytest
from scipy.integrate import quad

from numcosmo_py import Ncm, Nc

Ncm.cfg_init()

MODELS = [
    Nc.MultiplicityFuncCastroModel.C23,
    Nc.MultiplicityFuncCastroModel.C25,
]
HALO_FINDERS = [
    Nc.MultiplicityFuncCastroHaloFinder.ROCKSTAR,
    Nc.MultiplicityFuncCastroHaloFinder.AHF,
    Nc.MultiplicityFuncCastroHaloFinder.SUBFIND,
    Nc.MultiplicityFuncCastroHaloFinder.VELOCIRAPTOR,
]


@pytest.fixture(name="cosmo")
def fixture_cosmo() -> Nc.HICosmo:
    """Return a dark-energy cosmology, as required by the C25 calibration."""
    cosmo = Nc.HICosmoDEXcdm.new()
    cosmo.omega_x2omega_k()
    cosmo.param_set_by_name("H0", 70.0)
    cosmo.param_set_by_name("Omegac", 0.25)
    cosmo.param_set_by_name("Omegab", 0.05)
    cosmo.param_set_by_name("Omegak", 0.0)
    cosmo.param_set_by_name("w", -1.0)
    return cosmo


def _power_law_psf(n_index: float = -1.5) -> Ncm.PowspecFilter:
    """Return a tophat filter over a power law, with a known dlnsigma/dlnR."""
    ka = np.geomspace(1.0e-6, 1.0e3, 3000)
    za = np.linspace(0.0, 5.0, 20, dtype=np.float64)
    za_v, ka_v = np.meshgrid(za, ka)
    lnPk = np.log(1.0e-9 * ka_v**n_index) + 0.0 * za_v

    Pk2d = Ncm.Spline2dBicubic(
        spline=Ncm.SplineCubicNotaknot.new(),
        x_vector=Ncm.Vector.new_array(za.tolist()),
        y_vector=Ncm.Vector.new_array(np.log(ka).tolist()),
        z_matrix=Ncm.Matrix.new_array(lnPk.flatten().tolist(), len(za)),
    )
    ps = Ncm.PowspecSpline2d.new(Pk2d)
    ps.prepare()

    psf = Ncm.PowspecFilter.new(ps, Ncm.PowspecFilterType.TOPHAT)
    psf.set_best_lnr0()
    return psf


def _make_castro(
    model: Nc.MultiplicityFuncCastroModel,
    cosmo: Nc.HICosmo,
    halo_finder: Nc.MultiplicityFuncCastroHaloFinder = (
        Nc.MultiplicityFuncCastroHaloFinder.ROCKSTAR
    ),
    n_index: float = -1.5,
) -> Nc.MultiplicityFuncCastro:
    """Build a Castro multiplicity function wired to a power-law filter."""
    psf = _power_law_psf(n_index)
    psf.prepare(cosmo)

    mc = Nc.MultiplicityFuncCastro.new_full(model, halo_finder)
    mc.set_psf(psf)
    return mc


@pytest.mark.parametrize("model", MODELS)
def test_construction(model: Nc.MultiplicityFuncCastroModel, cosmo: Nc.HICosmo) -> None:
    """The object is constructible and reports what it was asked for."""
    mc = _make_castro(model, cosmo)

    assert isinstance(mc, Nc.MultiplicityFuncCastro)
    assert isinstance(mc, Nc.MultiplicityFunc)
    assert mc.get_model() == model
    # Castro is calibrated on virial masses.
    assert mc.get_mdef() == Nc.MultiplicityFuncMassDef.VIRIAL


def test_z_ta() -> None:
    """Turn-around redshift follows the Einstein-de Sitter relation."""
    mc = Nc.MultiplicityFuncCastro.new()

    for z in (0.0, 0.5, 1.0, 2.0):
        assert mc.z_ta(z) == pytest.approx((1.0 + z) * 2.0 ** (2.0 / 3.0) - 1.0)


def test_delta_c_kitayama_suto(cosmo: Nc.HICosmo) -> None:
    """delta_c follows the Kitayama & Suto fit in Omega_m(z)."""
    mc = Nc.MultiplicityFuncCastro.new()
    eds = 3.0 / 20.0 * (12.0 * math.pi) ** (2.0 / 3.0)

    for z in (0.0, 0.5, 1.0, 3.0):
        om_z = cosmo.E2Omega_m(z) / cosmo.E2(z)
        expected = eds * (1.0 + 0.012299 * math.log10(om_z))
        assert mc.delta_c(cosmo, z) == pytest.approx(expected, rel=1.0e-14)

    # Omega_m(z) < 1 at every epoch here, so the correction always lowers it.
    assert mc.delta_c(cosmo, 0.0) < eds


@pytest.mark.parametrize("model", MODELS)
@pytest.mark.parametrize("z", [0.0, 0.5, 1.0])
def test_normalization(
    model: Nc.MultiplicityFuncCastroModel, z: float, cosmo: Nc.HICosmo
) -> None:
    """The amplitude must normalise the multiplicity function exactly.

    Both A(p, q) and A(p, q, r) are closed-form normalisation constants, so
    int f(nu) dnu/nu = 1 identically. That makes this an absolute check on the
    gamma-function algebra, independent of any reference implementation.
    """
    mc = _make_castro(model, cosmo)
    lnR = math.log(8.0)
    delta_c = mc.delta_c(cosmo, z)

    def integrand(ln_nu: float) -> float:
        # eval takes sigma; invert nu = delta_c / sigma to sweep nu directly.
        sigma = delta_c / math.exp(ln_nu)
        return mc.eval(cosmo, sigma, lnR, z)

    integral, _ = quad(integrand, -40.0, 15.0, limit=500)

    assert integral == pytest.approx(1.0, rel=1.0e-6)


@pytest.mark.parametrize("halo_finder", HALO_FINDERS)
def test_c23_halo_finders_differ(
    halo_finder: Nc.MultiplicityFuncCastroHaloFinder, cosmo: Nc.HICosmo
) -> None:
    """Each halo finder has its own calibration, and all stay normalised."""
    mc = _make_castro(Nc.MultiplicityFuncCastroModel.C23, cosmo, halo_finder)
    assert mc.get_halo_finder() == halo_finder

    lnR = math.log(8.0)
    delta_c = mc.delta_c(cosmo, 0.0)

    def integrand(ln_nu: float) -> float:
        return mc.eval(cosmo, delta_c / math.exp(ln_nu), lnR, 0.0)

    integral, _ = quad(integrand, -40.0, 15.0, limit=500)
    assert integral == pytest.approx(1.0, rel=1.0e-6)


def test_c25_ignores_halo_finder(cosmo: Nc.HICosmo) -> None:
    """C25 has a single parameter set, so the halo finder must not matter."""
    lnR = math.log(8.0)
    sigma = 0.8

    values = []
    for halo_finder in HALO_FINDERS:
        mc = _make_castro(Nc.MultiplicityFuncCastroModel.C25, cosmo, halo_finder)
        values.append(mc.eval(cosmo, sigma, lnR, 0.0))

    for value in values[1:]:
        assert value == pytest.approx(values[0], rel=1.0e-14)


def test_c25_differs_from_c23(cosmo: Nc.HICosmo) -> None:
    """The two calibrations are distinct fits, not one nested in the other.

    Checked at w = -1, where the dark-energy correction vanishes: they still
    differ, because C25 carries an extra shape parameter.
    """
    lnR = math.log(8.0)
    c23 = _make_castro(Nc.MultiplicityFuncCastroModel.C23, cosmo)
    c25 = _make_castro(Nc.MultiplicityFuncCastroModel.C25, cosmo)

    differs = False
    for sigma in (0.4, 0.8, 1.2, 2.0):
        v23 = c23.eval(cosmo, sigma, lnR, 0.0)
        v25 = c25.eval(cosmo, sigma, lnR, 0.0)
        assert v23 > 0.0
        assert v25 > 0.0
        if abs(v25 / v23 - 1.0) > 1.0e-3:
            differs = True

    assert differs


def test_dark_energy_term_moves_c25(cosmo: Nc.HICosmo) -> None:
    """The C25 dark-energy correction must respond to w."""
    lnR = math.log(8.0)
    sigma = 0.8

    mc = _make_castro(Nc.MultiplicityFuncCastroModel.C25, cosmo)
    ref = mc.eval(cosmo, sigma, lnR, 0.0)

    cosmo.param_set_by_name("w", -0.8)
    assert mc.eval(cosmo, sigma, lnR, 0.0) != pytest.approx(ref, rel=1.0e-6)


def test_serialization() -> None:
    """The object round-trips through the serializer."""
    mc = Nc.MultiplicityFuncCastro.new_full(
        Nc.MultiplicityFuncCastroModel.C25,
        Nc.MultiplicityFuncCastroHaloFinder.SUBFIND,
    )
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    mc_dup = ser.dup_obj(mc)

    assert isinstance(mc_dup, Nc.MultiplicityFuncCastro)
    assert mc_dup.get_model() == Nc.MultiplicityFuncCastroModel.C25
    assert mc_dup.get_halo_finder() == Nc.MultiplicityFuncCastroHaloFinder.SUBFIND
    assert mc_dup.get_mdef() == mc.get_mdef()


@pytest.mark.parametrize("model", MODELS)
def test_lnf_matches_eval_full(model, cosmo: Nc.HICosmo) -> None:
    """The log form agrees with log(f) wherever f is representable."""
    mc = _make_castro(model, cosmo)

    for sigma in (2.0, 1.0, 0.5, 0.2, 0.1):
        f = mc.eval_full(cosmo, sigma, -0.5, 0.0)
        assert f > 0.0
        assert mc.eval_lnf(cosmo, sigma, -0.5, 0.0) == pytest.approx(
            math.log(f), rel=1.0e-13
        )


@pytest.mark.parametrize("model", MODELS)
def test_lnf_finite_where_f_underflows(model, cosmo: Nc.HICosmo) -> None:
    """The log form stays finite in the rare-peak regime where f underflows.

    Regression: taking log() of the underflowed f gave -inf, and the central
    difference of two -inf values is nan, which propagated into the bias.
    """
    mc = _make_castro(model, cosmo)

    # nu reaches ~1e8 at the smallest sigma here.
    for sigma in (1.0e-2, 1.0e-3, 1.0e-5, 1.0e-8):
        lnf = mc.eval_lnf(cosmo, sigma, -0.5, 0.0)
        assert math.isfinite(lnf), f"ln f not finite at sigma={sigma}"
        assert lnf < 0.0
        # f itself is expected to underflow; that is not an error.
        assert mc.eval_full(cosmo, sigma, -0.5, 0.0) == 0.0


@pytest.mark.parametrize("model", MODELS)
def test_lnf_decreases_with_peak_height(model, cosmo: Nc.HICosmo) -> None:
    """The log form must fall monotonically as sigma shrinks, with no jumps."""
    mc = _make_castro(model, cosmo)
    sigmas = [1.0, 0.5, 0.1, 1.0e-2, 1.0e-3, 1.0e-4]
    lnfs = [mc.eval_lnf(cosmo, sig, -0.5, 0.0) for sig in sigmas]

    for prev, cur in zip(lnfs, lnfs[1:]):
        assert math.isfinite(cur)
        assert cur < prev


@pytest.mark.parametrize("model", MODELS)
@pytest.mark.parametrize("slope", [-1.2, -1.0, -0.8, -0.5, -0.2, 0.0])
def test_normalization_over_wide_slope_range(model, slope: float, cosmo) -> None:
    """The gamma-function amplitude stays well conditioned across slopes.

    The 2025 amplitude carries (2p - q) and (q + r) in a denominator, so a
    parameter combination approaching either zero would wreck it. Normalisation
    is the sharpest probe available: it holds identically, so any conditioning
    problem shows up here immediately.
    """
    mc = _make_castro(model, cosmo)
    delta_c = mc.delta_c(cosmo, 0.0)

    def integrand(ln_nu: float) -> float:
        return mc.eval_full(cosmo, delta_c / math.exp(ln_nu), slope, 0.0)

    # f ~ nu^q at small nu, and q falls to ~0.28 at the shallow end, so the
    # lower tail needs far more room there than the default -40 provides.
    integral, _ = quad(integrand, -200.0, 15.0, limit=800)
    assert integral == pytest.approx(1.0, rel=1.0e-6)
