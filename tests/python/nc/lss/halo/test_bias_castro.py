#!/usr/bin/env python
#
# test_bias_castro.py
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

"""Tests for NcHaloBiasCastro.

How the truth table was produced
--------------------------------
``data/truth_tables/halo/nc_halo_bias_castro_pbs.bin`` holds
``[slope, z, sigma, b_pbs]``, with ``b_pbs`` obtained from CCToolkit's
``multiplicity_function_castro23`` by **complex-step differentiation**:
``dlnf/dlnnu = Im[f(nu(1 + ih))]/h / f(nu)`` with ``h = 1e-20``.

CCToolkit's own ``pbs_bias`` takes this derivative with ``np.gradient`` over a
mass grid, which is second order and converges only to ~6e-5 even at 128 points
per decade (and is off by ~1% at the 4-8 points per decade the API invites).
The complex step has no subtractive cancellation, so it is exact to round-off
and is the better reference. Our own value agrees with it to 5.7e-11.

Run this module as a script to rewrite the table.

A power law is used throughout because it makes ``dlnsigma/dlnR`` constant, so
``ds/dlnR = 0`` and the total derivative along the mass direction coincides with
the partial one -- which is what lets the table be a pure function of the slope.
``test_eval_total_equals_partial_for_power_law`` checks that the implementation
respects that.
"""

import math

import numpy as np
import pytest

from numcosmo_py import Ncm, Nc

Ncm.cfg_init()

_PBS_TABLE = "truth_tables/halo/nc_halo_bias_castro_pbs.bin"

# Central differences on the closed-form model, against an exact reference.
RTOL_PBS = 1.0e-9


def _load_table(name: str) -> np.ndarray:
    """Load a truth table as an (nrows, ncols) array."""
    path = Ncm.cfg_get_data_filename(name, True)
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    matrix = ser.from_binfile(path)
    assert isinstance(matrix, Ncm.Matrix)
    return np.array(matrix.dup_array()).reshape(matrix.nrows(), matrix.ncols())


PBS_GOLDEN = _load_table(_PBS_TABLE)


def _make_cosmo(w: float = -1.0) -> Nc.HICosmo:
    """Return the cosmology the truth table was generated with."""
    cosmo = Nc.HICosmoDEXcdm.new()
    cosmo.omega_x2omega_k()
    cosmo.param_set_by_name("H0", 70.0)
    cosmo.param_set_by_name("Omegac", 0.25)
    cosmo.param_set_by_name("Omegab", 0.05)
    cosmo.param_set_by_name("Omegak", 0.0)
    cosmo.param_set_by_name("w", w)
    return cosmo


def _build_psf(slope: float, amplitude: float) -> Ncm.PowspecFilter:
    """Return a tophat filter over a power law of the requested slope."""
    n_index = -2.0 * slope - 3.0
    ka = np.geomspace(1.0e-10, 1.0e3, 6000)
    za = np.linspace(0.0, 5.0, 20, dtype=np.float64)
    za_v, ka_v = np.meshgrid(za, ka)
    lnPk = np.log(amplitude * ka_v**n_index) + 0.0 * za_v

    Pk2d = Ncm.Spline2dBicubic(
        spline=Ncm.SplineCubicNotaknot.new(),
        x_vector=Ncm.Vector.new_array(za.tolist()),
        y_vector=Ncm.Vector.new_array(np.log(ka).tolist()),
        z_matrix=Ncm.Matrix.new_array(lnPk.flatten().tolist(), len(za)),
    )
    ps = Ncm.PowspecSpline2d.new(Pk2d)
    ps.prepare()

    psf = Ncm.PowspecFilter.new(ps, Ncm.PowspecFilterType.TOPHAT)
    psf.set_reltol(1.0e-9)
    psf.set_best_lnr0()
    return psf


def _make_psf(slope: float, cosmo: Nc.HICosmo) -> Ncm.PowspecFilter:
    """Return a power-law filter normalised to a realistic sigma_8.

    The amplitude matters here, unlike in the tests that pass ``sigma``
    directly: ``eval`` derives sigma from the filter, and an arbitrary
    amplitude puts it at ~1e-6, where nu ~ 1e5 and exp(-a nu^2 / 2) underflows.
    Since sigma scales as the square root of the amplitude, one rescaling pass
    lands it exactly.
    """
    target_sigma8 = 0.8
    lnR8 = math.log(8.0 / cosmo.h())

    psf = _build_psf(slope, 1.0e-9)
    psf.prepare(cosmo)
    sigma8 = psf.eval_sigma_lnr(0.0, lnR8)

    return _build_psf(slope, 1.0e-9 * (target_sigma8 / sigma8) ** 2)


def _make_bias(
    slope: float, cosmo: Nc.HICosmo
) -> tuple[Nc.HaloBiasCastro, Nc.HaloMassFunction]:
    """Build a Castro bias on a mass function using a Castro multiplicity."""
    psf = _make_psf(slope, cosmo)
    mulf = Nc.MultiplicityFuncCastro.new_full(
        Nc.MultiplicityFuncCastroModel.C23,
        Nc.MultiplicityFuncCastroHaloFinder.ROCKSTAR,
    )
    dist = Nc.Distance.new(5.0)
    mfp = Nc.HaloMassFunction.new(dist, psf, mulf)
    bias = Nc.HaloBiasCastro.new(mfp)
    psf.prepare(cosmo)
    return bias, mfp


def test_construction() -> None:
    """The object builds and carries the calibrated coefficients."""
    cosmo = _make_cosmo()
    bias, mfp = _make_bias(-0.5, cosmo)

    assert isinstance(bias, Nc.HaloBiasCastro)
    assert isinstance(bias, Nc.HaloBias)
    assert bias.peek_mass_function() == mfp

    assert bias.get_A0() == pytest.approx(1.150)
    assert bias.get_a1() == pytest.approx(0.0929)
    assert bias.get_b1() == pytest.approx(0.256)
    assert bias.get_b2() == pytest.approx(0.173)
    assert bias.get_c1() == pytest.approx(-0.0372)


def test_requires_second_derivative() -> None:
    """Construction raises the filter's derivative order to two."""
    cosmo = _make_cosmo()
    psf = _make_psf(-0.5, cosmo)
    assert psf.get_nderivs() == 1

    mulf = Nc.MultiplicityFuncCastro.new()
    dist = Nc.Distance.new(5.0)
    mfp = Nc.HaloMassFunction.new(dist, psf, mulf)
    Nc.HaloBiasCastro.new(mfp)

    assert psf.get_nderivs() == 2


@pytest.mark.parametrize("row", PBS_GOLDEN)
def test_pbs_matches_exact_derivative(row: np.ndarray) -> None:
    """b_PBS reproduces the complex-step reference derivative."""
    slope, z, sigma, expected = row
    cosmo = _make_cosmo()
    bias, _ = _make_bias(slope, cosmo)

    assert bias.pbs(cosmo, sigma, slope, z) == pytest.approx(expected, rel=RTOL_PBS)


def test_S8_uses_h_inverse_radius() -> None:
    """S8 is sigma(8 h^-1 Mpc) rescaled, so it must carry the h factor."""
    cosmo = _make_cosmo()
    bias, mfp = _make_bias(-0.5, cosmo)

    psf = mfp.peek_psf()
    h = cosmo.h()
    sigma8 = psf.eval_sigma_lnr(0.0, math.log(8.0 / h))
    expected = sigma8 * math.sqrt(cosmo.Omega_m0() / 0.3)

    assert bias.S8(cosmo) == pytest.approx(expected, rel=1.0e-14)
    # Using 8 Mpc instead of 8 Mpc/h would be a different number entirely.
    wrong = psf.eval_sigma_lnr(0.0, math.log(8.0)) * math.sqrt(cosmo.Omega_m0() / 0.3)
    assert bias.S8(cosmo) != pytest.approx(wrong, rel=1.0e-3)


@pytest.mark.parametrize("slope", [-0.8, -0.5, -0.2])
@pytest.mark.parametrize("z", [0.0, 1.0])
def test_correction_factor(slope: float, z: float) -> None:
    """The correction reproduces its defining product."""
    cosmo = _make_cosmo()
    bias, _ = _make_bias(slope, cosmo)

    om_z = cosmo.E2Omega_m(z) / cosmo.E2(z)
    expected = (
        bias.get_A0()
        * (1.0 + bias.get_a1() * om_z)
        * (1.0 + bias.get_b1() * slope + bias.get_b2() * slope**2)
        * (1.0 + bias.get_c1() * bias.S8(cosmo))
    )

    assert bias.correction(cosmo, slope, z) == pytest.approx(expected, rel=1.0e-14)


@pytest.mark.parametrize("slope", [-0.8, -0.5])
def test_eval_total_equals_partial_for_power_law(slope: float) -> None:
    """For a scale-free spectrum the slope does not run, killing the extra term.

    eval() adds ``dlnf/ds * ds/dlnnu`` to the partial derivative. A power law has
    ``ds/dlnR = 0`` exactly, so eval() must collapse onto pbs() times the
    correction. This is what checks that the extra term is wired with the right
    sign and vanishes when it should.
    """
    cosmo = _make_cosmo()
    bias, mfp = _make_bias(slope, cosmo)
    z = 0.5

    for lnM in (math.log(1.0e13), math.log(1.0e14), math.log(1.0e15)):
        lnR = mfp.lnM_to_lnR(cosmo, lnM)
        sigma = mfp.sigma_lnR(cosmo, lnR, z)
        s_actual = 0.5 * mfp.peek_psf().eval_dnlnvar_dlnrn(z, lnR, 1)

        expected = bias.pbs(cosmo, sigma, s_actual, z) * bias.correction(
            cosmo, s_actual, z
        )
        assert bias.eval(cosmo, sigma, lnM, z) == pytest.approx(expected, rel=1.0e-6)


def test_integrand_uses_eval() -> None:
    """The mean-bias integrand is the mass function times the bias."""
    cosmo = _make_cosmo()
    bias, mfp = _make_bias(-0.5, cosmo)
    mfp.set_eval_limits(cosmo, math.log(1.0e13), math.log(1.0e15), 0.0, 1.0)
    mfp.prepare(cosmo)

    lnM, z = math.log(1.0e14), 0.4
    sigma = mfp.sigma_lnM(cosmo, lnM, z)
    expected = mfp.d2n_dzdlnM(cosmo, lnM, z) * bias.eval(cosmo, sigma, lnM, z)

    assert bias.integrand(cosmo, lnM, z) == pytest.approx(expected, rel=1.0e-14)


def test_serialization() -> None:
    """The object round-trips through the serializer."""
    cosmo = _make_cosmo()
    bias, _ = _make_bias(-0.5, cosmo)

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    dup = ser.dup_obj(bias)

    assert isinstance(dup, Nc.HaloBiasCastro)
    assert dup.get_A0() == pytest.approx(bias.get_A0())
    assert dup.get_c1() == pytest.approx(bias.get_c1())


@pytest.mark.parametrize("sigma", [1.0e-2, 1.0e-3, 1.0e-5, 1.0e-8])
def test_pbs_finite_in_rare_peak_regime(sigma: float) -> None:
    """b_PBS stays finite where f itself underflows to zero.

    Regression: the derivative used to be taken as log(f(sigma_+)) -
    log(f(sigma_-)). Once f underflowed both logs became -inf and their
    difference nan, which then propagated through eval() into the likelihood.
    Differencing ln f directly keeps it finite.
    """
    cosmo = _make_cosmo()
    bias, _ = _make_bias(-0.5, cosmo)

    b_pbs = bias.pbs(cosmo, sigma, -0.5, 0.0)

    assert math.isfinite(b_pbs), f"b_PBS not finite at sigma={sigma}"
    assert b_pbs > 1.0


def test_pbs_follows_rare_peak_asymptote() -> None:
    """Deep in the rare-peak regime b_PBS must grow as nu^2.

    The peak-background split gives b -> 1 + a nu^2 / delta_c there, so halving
    sigma has to quadruple b - 1. This checks the limit is not merely finite but
    the physically right one.
    """
    cosmo = _make_cosmo()
    bias, _ = _make_bias(-0.5, cosmo)

    ratios = []
    for sigma in (1.0e-4, 1.0e-5, 1.0e-6):
        b_1 = bias.pbs(cosmo, sigma, -0.5, 0.0) - 1.0
        b_2 = bias.pbs(cosmo, 0.5 * sigma, -0.5, 0.0) - 1.0
        ratios.append(b_2 / b_1)

    for ratio in ratios:
        assert ratio == pytest.approx(4.0, rel=1.0e-6)


def test_eval_finite_for_extreme_masses() -> None:
    """eval() stays finite across a mass range far wider than any survey."""
    cosmo = _make_cosmo()
    bias, mfp = _make_bias(-0.5, cosmo)

    for lnM in (math.log(1.0e10), math.log(1.0e14), math.log(1.0e17)):
        lnR = mfp.lnM_to_lnR(cosmo, lnM)
        sigma = mfp.sigma_lnR(cosmo, lnR, 0.5)
        value = bias.eval(cosmo, sigma, lnM, 0.5)

        assert math.isfinite(value), f"bias not finite at lnM={lnM}"
        assert value > 0.0


def _regenerate() -> None:  # pragma: no cover - developer tool
    """Rewrite the truth table from CCToolkit via complex-step differentiation."""
    # pylint: disable=import-outside-toplevel
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "cc_hmf", "../CCToolkit/cctoolkit/hmf.py"
    )
    hmf = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(hmf)

    step = 1.0e-20
    cosmo = _make_cosmo()
    rows = []
    for slope in sorted(set(PBS_GOLDEN[:, 0])):
        for z in sorted(set(PBS_GOLDEN[:, 1])):
            om = cosmo.E2Omega_m(z) / cosmo.E2(z)
            d_c = (
                3.0
                / 20.0
                * (12.0 * math.pi) ** (2.0 / 3.0)
                * (1.0 + 0.012299 * math.log10(om))
            )
            for sigma in sorted(set(PBS_GOLDEN[:, 2])):
                nu = d_c / sigma
                f_0 = hmf.multiplicity_function_castro23(
                    nu, slope, om, hmf.best_fit_values_ROCKSTAR
                )
                f_i = hmf.multiplicity_function_castro23(
                    nu * (1.0 + 1j * step), slope, om, hmf.best_fit_values_ROCKSTAR
                )
                b_pbs = 1.0 - ((f_i.imag / step) / f_0.real) / d_c
                rows.append([slope, z, sigma, float(b_pbs)])

    arr = np.array(rows, dtype=np.float64)
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    matrix = Ncm.Matrix.new_array(arr.flatten().tolist(), arr.shape[1])
    ser.to_binfile(matrix, f"data/{_PBS_TABLE}")
    print(f"wrote data/{_PBS_TABLE}: {arr.shape[0]} x {arr.shape[1]}")


if __name__ == "__main__":  # pragma: no cover
    _regenerate()
