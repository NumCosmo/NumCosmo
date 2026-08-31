#!/usr/bin/env python
#
# test_kernel_analytic.py
#
# Fri Aug 21 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_kernel_analytic.py
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

"""Tests for the analytic xcor kernels.

Every other kernel is spline-backed, so its C_ell has no independently known
value. These are defined by a formula, and the point of the tests is that the
formula is what the library actually integrates: the window is checked against
its closed form, and the C_ell against a quadrature written here from the same
definition.
"""

import math
import typing
from decimal import Decimal, localcontext

import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.integrate import quad
from scipy.special import (  # pylint: disable=no-name-in-module
    erf,
    erfc,
    spherical_jn,
)
from scipy.stats import t as student_t

from numcosmo_py import Nc, Ncm
from numcosmo_py.cosmology import Cosmology

pytest_plugins = [
    "python.fixtures_xcor",
]

pytestmark = pytest.mark.xcor

CHI_MEAN, CHI_SIGMA, N_SIGMA = 1500.0, 300.0, 4.0
CHI_LOWER, CHI_UPPER = 500.0, 2500.0
ST_MEAN, ST_SCALE, ST_NU, ST_NSCALE = 1500.0, 200.0, 2.0, 6.0
MULTI_MEAN, MULTI_SIGMA, MULTI_WEIGHT = [1000.0, 1600.0], [300.0, 300.0], [1.0, 0.6]
MULTI_NSIGMA = 4.0
MULTI_DISJOINT_MEAN, MULTI_DISJOINT_SIGMA = [600.0, 2600.0], [100.0, 150.0]
PE_SCALE, PE_ALPHA, PE_BETA, PE_LOWER, PE_UPPER = 1200.0, 2.0, 1.5, 50.0, 4000.0
TS_LOWER, TS_UPPER, TS_SIGMA, TS_NSIGMA = 1000.0, 2000.0, 150.0, 6.0
LENS_LOWER, LENS_SRC_LOWER, LENS_SRC_UPPER = 50.0, 2000.0, 3000.0


def _integrator() -> Ncm.SBesselIntegrator:
    return Ncm.SBesselIntegratorLevin.new(0, 8)


# Every kernel here runs at the library's default tolerances. scaled-abstol
# floors W_ell(k) against its own peak, but the C_ell integrand is k^2 W_1 W_2,
# so the floor enters squared: the 1e-4 default is already 1e-8 on what is
# integrated, which is about what the outer integral can carry. Overriding it to
# 1e-8 -- as some xcor tests do -- asks for 1e-16 there, costs up to 500x the
# knots on a compactly-supported window, and measurably buys nothing: against
# the quadrature below the agreement saturates at a few times 1e-6 either way.


def _gauss(cosmology: Cosmology, l_limber: int = -1) -> Nc.XcorKernelAnalyticGauss:
    """Gaussian window, prepared and in the requested Limber tier."""
    kernel = Nc.XcorKernelAnalyticGauss.new_full(
        cosmology.dist, cosmology.ps_ml, CHI_MEAN, CHI_SIGMA, N_SIGMA, _integrator()
    )
    kernel.set_l_limber(l_limber)
    kernel.prepare(cosmology.cosmo)

    return kernel


def _tophat(cosmology: Cosmology, l_limber: int = -1) -> Nc.XcorKernelAnalyticTophat:
    """Top-hat window, prepared and in the requested Limber tier."""
    kernel = Nc.XcorKernelAnalyticTophat.new_full(
        cosmology.dist, cosmology.ps_ml, CHI_LOWER, CHI_UPPER, _integrator()
    )
    kernel.set_l_limber(l_limber)
    kernel.prepare(cosmology.cosmo)

    return kernel


def _student_t(
    cosmology: Cosmology, l_limber: int = -1
) -> Nc.XcorKernelAnalyticStudentT:
    """Power-law-tailed window, prepared and in the requested Limber tier."""
    kernel = Nc.XcorKernelAnalyticStudentT.new_full(
        cosmology.dist,
        cosmology.ps_ml,
        ST_MEAN,
        ST_SCALE,
        ST_NU,
        ST_NSCALE,
        _integrator(),
    )
    kernel.set_l_limber(l_limber)
    kernel.prepare(cosmology.cosmo)

    return kernel


def _multi(
    cosmology: Cosmology,
    means: list[float],
    sigmas: list[float],
    l_limber: int = -1,
) -> Nc.XcorKernelAnalyticMulti:
    """Multimodal window, prepared and in the requested Limber tier."""
    kernel = Nc.XcorKernelAnalyticMulti.new_full(
        cosmology.dist,
        cosmology.ps_ml,
        Ncm.Vector.new_array(means),
        Ncm.Vector.new_array(sigmas),
        Ncm.Vector.new_array(MULTI_WEIGHT),
        MULTI_NSIGMA,
        _integrator(),
    )
    kernel.set_l_limber(l_limber)
    kernel.prepare(cosmology.cosmo)

    return kernel


def _power_exp(
    cosmology: Cosmology, l_limber: int = -1
) -> Nc.XcorKernelAnalyticPowerExp:
    """Power-law rise, stretched-exponential fall: the dn/dz family."""
    kernel = Nc.XcorKernelAnalyticPowerExp.new_full(
        cosmology.dist,
        cosmology.ps_ml,
        PE_SCALE,
        PE_ALPHA,
        PE_BETA,
        PE_LOWER,
        PE_UPPER,
        _integrator(),
    )
    kernel.set_l_limber(l_limber)
    kernel.prepare(cosmology.cosmo)

    return kernel


def _tophat_smooth(
    cosmology: Cosmology, chi_sigma: float = TS_SIGMA, l_limber: int = -1
) -> Nc.XcorKernelAnalyticTophatSmooth:
    """Top-hat convolved with a Gaussian: a real tomographic bin."""
    kernel = Nc.XcorKernelAnalyticTophatSmooth.new_full(
        cosmology.dist,
        cosmology.ps_ml,
        TS_LOWER,
        TS_UPPER,
        chi_sigma,
        TS_NSIGMA,
        _integrator(),
    )
    kernel.set_l_limber(l_limber)
    kernel.prepare(cosmology.cosmo)

    return kernel


def _lensing(cosmology: Cosmology, l_limber: int = -1) -> Nc.XcorKernelAnalyticLensing:
    """Lensing efficiency of a top-hat source bin."""
    kernel = Nc.XcorKernelAnalyticLensing.new_full(
        cosmology.dist,
        cosmology.ps_ml,
        LENS_LOWER,
        LENS_SRC_LOWER,
        LENS_SRC_UPPER,
        _integrator(),
    )
    kernel.set_l_limber(l_limber)
    kernel.prepare(cosmology.cosmo)

    return kernel


def _lensing_shape(chi: float) -> float:
    """The unnormalized efficiency, written here from the geometry.

    Evaluated in ``decimal`` rather than ``math``: on the upper branch the true
    value goes as ``(b - chi)^2 / 2chi`` as chi approaches b, so the difference
    ``(b - chi) - chi log (b/chi)`` cancels away every significant digit near
    the edge -- in float it is wrong by a factor of ninety, with the wrong sign,
    one ulp inside. The library evaluates that branch through a stable
    rearrangement; a reference written the naive way would be the inaccurate
    side of the comparison. ``Decimal.ln`` is correctly rounded and in the
    standard library, so this stays an independent check rather than a copy of
    the library's own rearrangement.
    """
    a, b = LENS_SRC_LOWER, LENS_SRC_UPPER

    if chi >= b:
        return 0.0

    with localcontext() as ctx:
        ctx.prec = 50
        chi_d = Decimal(repr(chi))
        b_d = Decimal(repr(b))
        m_d = Decimal(repr(a)) if chi <= a else chi_d

        return float(
            chi_d * ((b_d - m_d) - chi_d * (b_d / m_d).ln()) / (b_d - Decimal(repr(a)))
        )


def test_gauss_support_is_the_truncation() -> None:
    """The support is the mean plus/minus n sigma, clipped at the observer."""
    kernel = Nc.XcorKernelAnalyticGauss(
        dist=Nc.Distance.new(3.0),
        powspec=Nc.PowspecMLTransfer.new(Nc.TransferFuncEH.new()),
        chi_mean=CHI_MEAN,
        chi_sigma=CHI_SIGMA,
        n_sigma=N_SIGMA,
    )

    chi_min, chi_max = kernel.get_support()

    assert chi_min == CHI_MEAN - N_SIGMA * CHI_SIGMA
    assert chi_max == CHI_MEAN + N_SIGMA * CHI_SIGMA


def test_gauss_support_clips_at_the_observer() -> None:
    """A window wider than its distance to the observer starts at chi = 0."""
    kernel = Nc.XcorKernelAnalyticGauss(
        dist=Nc.Distance.new(3.0),
        powspec=Nc.PowspecMLTransfer.new(Nc.TransferFuncEH.new()),
        chi_mean=300.0,
        chi_sigma=300.0,
        n_sigma=4.0,
    )

    chi_min, _ = kernel.get_support()

    assert chi_min == 0.0


def test_gauss_window_matches_closed_form(cosmology: Cosmology) -> None:
    """eval_W is the truncated Gaussian, renormalized over what survives."""
    kernel = _gauss(cosmology)
    chi_min, chi_max = kernel.get_support()
    chi = np.linspace(chi_min, chi_max, 257)

    s2 = np.sqrt(2.0) * CHI_SIGMA
    norm = (
        0.5
        * np.sqrt(2.0 * np.pi)
        * CHI_SIGMA
        * (math.erf((chi_max - CHI_MEAN) / s2) - math.erf((chi_min - CHI_MEAN) / s2))
    )
    expected = np.exp(-0.5 * ((chi - CHI_MEAN) / CHI_SIGMA) ** 2) / norm

    got = np.array([kernel.eval_W(c) for c in chi])

    assert_allclose(got, expected, rtol=1.0e-14)


@pytest.mark.parametrize(
    "shape",
    ["gauss", "tophat", "student_t", "power_exp", "tophat_smooth", "lensing", "multi"],
)
def test_window_vanishes_outside_support(cosmology: Cosmology, shape: str) -> None:
    """The truncation is part of the definition, not a tolerance.

    Every shape guards its support with its own expression of it -- a pair of
    limits, a mean plus a width, a smoothed interval -- so the check runs per
    shape rather than on a representative one. Inside, the window is alive
    right up to the lower edge. The upper edge has one exception: the lensing
    efficiency of a source bin is geometrically zero at the far edge of that
    bin, since nothing behind it is being lensed.
    """
    kernel = _any_shape(cosmology, shape)
    chi_min, chi_max = kernel.get_support()

    assert kernel.eval_W(chi_min - 1.0e-3) == 0.0
    assert kernel.eval_W(chi_max + 1.0e-3) == 0.0

    assert kernel.eval_W(chi_min) > 0.0
    assert kernel.eval_W(0.5 * (chi_min + chi_max)) > 0.0

    if shape == "lensing":
        assert kernel.eval_W(chi_max) == 0.0
    else:
        assert kernel.eval_W(chi_max) > 0.0


@pytest.mark.parametrize("shape", ["gauss", "tophat"])
def test_window_is_normalized(cosmology: Cosmology, shape: str) -> None:
    """int W dchi = 1 over the support, which is what fixes the C_ell scale."""
    kernel = _gauss(cosmology) if shape == "gauss" else _tophat(cosmology)
    chi_min, chi_max = kernel.get_support()

    total, _ = quad(kernel.eval_W, chi_min, chi_max, limit=400)

    assert_allclose(total, 1.0, rtol=1.0e-12)


def test_tophat_window_is_flat(cosmology: Cosmology) -> None:
    """Inside the bin the window is one constant; outside it is zero."""
    kernel = _tophat(cosmology)
    height = 1.0 / (CHI_UPPER - CHI_LOWER)

    assert kernel.get_support() == (CHI_LOWER, CHI_UPPER)
    assert_allclose(
        [kernel.eval_W(c) for c in (CHI_LOWER, 1500.0, CHI_UPPER)], height, rtol=1.0e-15
    )
    assert kernel.eval_W(CHI_LOWER - 1.0) == 0.0
    assert kernel.eval_W(CHI_UPPER + 1.0) == 0.0


def test_z_range_is_the_support_in_redshift(cosmology: Cosmology) -> None:
    """get_z_range converts the support, so mapping it back returns it."""
    kernel = _gauss(cosmology)
    cosmo = cosmology.cosmo
    chi_min, chi_max = kernel.get_support()
    zmin, zmax, zmid = kernel.get_z_range()

    rh = cosmo.RH_Mpc()
    got = [cosmology.dist.comoving(cosmo, z) * rh for z in (zmin, zmax, zmid)]

    assert_allclose(got, [chi_min, chi_max, 0.5 * (chi_min + chi_max)], rtol=1.0e-9)


@pytest.mark.parametrize(
    "shape",
    ["gauss", "tophat", "student_t", "power_exp", "tophat_smooth", "lensing", "multi"],
)
def test_serialization_roundtrip(cosmology: Cosmology, shape: str) -> None:
    """The kernel is a registered, serializable object: a spec survives a file.

    Every shape, not a representative pair: a benchmark specification is stored
    by serializing the kernel, so a shape whose properties did not round-trip
    would silently produce a spec that does not describe the run.
    """
    kernel = _any_shape(cosmology, shape)
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)

    dup = ser.dup_obj(kernel)

    assert type(dup) is type(kernel)
    assert dup.get_support() == kernel.get_support()
    assert_allclose(
        [dup.eval_W(c) for c in np.linspace(*kernel.get_support(), 32)],
        [kernel.eval_W(c) for c in np.linspace(*kernel.get_support(), 32)],
        rtol=1.0e-15,
    )


def _ctor_args(cosmology: Cosmology, shape: str) -> tuple[typing.Any, tuple]:
    """The class of a shape and the leading arguments both its constructors take."""
    head: tuple = (cosmology.dist, cosmology.ps_ml)
    cls, tail = {
        "gauss": (Nc.XcorKernelAnalyticGauss, (CHI_MEAN, CHI_SIGMA, N_SIGMA)),
        "tophat": (Nc.XcorKernelAnalyticTophat, (CHI_LOWER, CHI_UPPER)),
        "student_t": (
            Nc.XcorKernelAnalyticStudentT,
            (ST_MEAN, ST_SCALE, ST_NU, ST_NSCALE),
        ),
        "power_exp": (
            Nc.XcorKernelAnalyticPowerExp,
            (PE_SCALE, PE_ALPHA, PE_BETA, PE_LOWER, PE_UPPER),
        ),
        "tophat_smooth": (
            Nc.XcorKernelAnalyticTophatSmooth,
            (TS_LOWER, TS_UPPER, TS_SIGMA, TS_NSIGMA),
        ),
        "lensing": (
            Nc.XcorKernelAnalyticLensing,
            (LENS_LOWER, LENS_SRC_LOWER, LENS_SRC_UPPER),
        ),
        "multi": (
            Nc.XcorKernelAnalyticMulti,
            (
                Ncm.Vector.new_array(MULTI_MEAN),
                Ncm.Vector.new_array(MULTI_SIGMA),
                Ncm.Vector.new_array(MULTI_WEIGHT),
                MULTI_NSIGMA,
            ),
        ),
    }[shape]

    return cls, head + tail


@pytest.mark.parametrize(
    "shape",
    ["gauss", "tophat", "student_t", "power_exp", "tophat_smooth", "lensing", "multi"],
)
def test_the_two_constructors_differ_only_by_the_integrator(
    cosmology: Cosmology, shape: str
) -> None:
    """`new` is the Limber-only constructor, `new_full` the one that carries an sbi.

    The builders above all use `new_full`, since a kernel only accepts the
    non-Limber tiers once it holds an integrator. `new` must otherwise describe
    the same window, so the two are compared shape by shape.
    """
    cls, args = _ctor_args(cosmology, shape)

    plain = cls.new(*args)
    full = cls.new_full(*args, _integrator())

    for kernel in (plain, full):
        kernel.prepare(cosmology.cosmo)

    assert plain.peek_integrator() is None
    assert full.peek_integrator() is not None

    assert plain.get_support() == full.get_support()
    chis = np.linspace(*plain.get_support(), 32)
    assert_allclose(
        [plain.eval_W(chi) for chi in chis],
        [full.eval_W(chi) for chi in chis],
        rtol=1.0e-15,
    )


GX, GW = np.polynomial.legendre.leggauss(12)


def _gauss_with_kdep(
    cosmology: Cosmology, kdep: typing.Optional[Nc.XcorKernelRadialKDep]
) -> Nc.XcorKernelAnalyticGauss:
    """The baseline Gaussian, optionally carrying a scale dependence.

    The property is nullable -- no factor at all is a distinct state from a
    factor that evaluates to one -- but the generated stubs cannot express
    that, so the kwarg is omitted rather than passed as None.
    """
    kwargs: dict[str, typing.Any] = {}

    if kdep is not None:
        kwargs["scale_dependence"] = kdep

    kernel = Nc.XcorKernelAnalyticGauss(
        dist=cosmology.dist,
        powspec=cosmology.ps_ml,
        chi_mean=CHI_MEAN,
        chi_sigma=CHI_SIGMA,
        n_sigma=N_SIGMA,
        integrator=_integrator(),
        **kwargs,
    )
    kernel.set_l_limber(-1)
    kernel.prepare(cosmology.cosmo)

    return kernel


def _any_shape(cosmology: Cosmology, shape: str) -> Nc.XcorKernelRadial:
    """Any shape by name, including the ones the C_ell cases do not cover."""
    if shape == "multi":
        return _multi(cosmology, MULTI_MEAN, MULTI_SIGMA)

    if shape == "tophat_smooth":
        return _tophat_smooth(cosmology)

    return _cl_kernel(cosmology, shape)


def _cl_kernel(cosmology: Cosmology, shape: str) -> Nc.XcorKernelRadial:
    """The kernel a C_ell case names."""
    builder = {
        "gauss": _gauss,
        "tophat": _tophat,
        "student_t": _student_t,
        "power_exp": _power_exp,
        "lensing": _lensing,
    }[shape]

    return builder(cosmology)


def _window(shape: str) -> typing.Callable[[np.ndarray], np.ndarray]:
    """The window's closed form, written here in numpy.

    The reference must not reach the window through the library, or the check
    degenerates into comparing a routine with itself.
    """
    if shape == "tophat":
        return lambda chi: np.where(
            (chi >= CHI_LOWER) & (chi <= CHI_UPPER), 1.0 / (CHI_UPPER - CHI_LOWER), 0.0
        )

    if shape == "power_exp":
        norm = quad(
            lambda c: (c / PE_SCALE) ** PE_ALPHA * np.exp(-((c / PE_SCALE) ** PE_BETA)),
            PE_LOWER,
            PE_UPPER,
            limit=600,
        )[0]

        def power_exp_window(chi: np.ndarray) -> np.ndarray:
            x = chi / PE_SCALE

            return np.where(
                (chi >= PE_LOWER) & (chi <= PE_UPPER),
                x**PE_ALPHA * np.exp(-(x**PE_BETA)) / norm,
                0.0,
            )

        return power_exp_window

    if shape == "lensing":
        norm = (
            quad(_lensing_shape, LENS_LOWER, LENS_SRC_LOWER, limit=400)[0]
            + quad(_lensing_shape, LENS_SRC_LOWER, LENS_SRC_UPPER, limit=400)[0]
        )

        def lensing_window(chi: np.ndarray) -> np.ndarray:
            flat = np.ravel(chi)
            out = np.array([_lensing_shape(float(c)) for c in flat]) / norm

            return out.reshape(np.shape(chi))

        return lensing_window

    if shape == "student_t":
        lo = max(0.0, ST_MEAN - ST_NSCALE * ST_SCALE)
        hi = ST_MEAN + ST_NSCALE * ST_SCALE
        kept = student_t.cdf((hi - ST_MEAN) / ST_SCALE, ST_NU) - student_t.cdf(
            (lo - ST_MEAN) / ST_SCALE, ST_NU
        )

        return lambda chi: np.where(
            (chi >= lo) & (chi <= hi),
            student_t.pdf((chi - ST_MEAN) / ST_SCALE, ST_NU) / ST_SCALE / kept,
            0.0,
        )

    chi_min = max(0.0, CHI_MEAN - N_SIGMA * CHI_SIGMA)
    chi_max = CHI_MEAN + N_SIGMA * CHI_SIGMA
    s2 = np.sqrt(2.0) * CHI_SIGMA
    norm = (
        0.5
        * np.sqrt(2.0 * np.pi)
        * CHI_SIGMA
        * (math.erf((chi_max - CHI_MEAN) / s2) - math.erf((chi_min - CHI_MEAN) / s2))
    )

    return lambda chi: np.where(
        (chi >= chi_min) & (chi <= chi_max),
        np.exp(-0.5 * ((chi - CHI_MEAN) / CHI_SIGMA) ** 2) / norm,
        0.0,
    )


def _cl_reference(
    kernel: Nc.XcorKernelRadial,
    cosmology: Cosmology,
    shape: str,
    ell: int,
    *,
    n_k: int,
    per_lobe: int,
) -> float:
    """C_ell for an analytic window, from the definition, in Mpc throughout.

    C_ell = (2/pi) int dk k^2 P(k, z=0) [int dchi W(chi) j_ell(k chi)]^2. The
    radial integral is split into panels a fraction of a half-oscillation wide
    and each panel taken by Gauss-Legendre; the outer one is a log-k trapezoid
    over the range the kernel itself reports. Only that range and P(k) come
    from the library -- the window and the quadrature are written here.
    """
    cosmo, powspec = cosmology.cosmo, cosmology.ps_ml
    chi_min, chi_max = kernel.get_support()
    k_min, k_max = kernel.get_eval(cosmo, ell, Nc.XcorKernelClosure.SPLINE).get_range()
    rh = cosmo.RH_Mpc()
    window = _window(shape)

    def radial(k: float) -> float:
        n_panel = max(16, per_lobe * int((chi_max - chi_min) * k / np.pi) + 1)
        edges = np.linspace(chi_min, chi_max, n_panel + 1)
        lo, hi = edges[:-1, None], edges[1:, None]
        chi = 0.5 * (hi - lo) * GX[None, :] + 0.5 * (hi + lo)
        f = window(chi) * spherical_jn(ell, k * chi)

        return float(np.sum(0.5 * (hi - lo) * f * GW[None, :]))

    lnk = np.linspace(np.log(k_min / rh), np.log(k_max / rh), n_k)
    k = np.exp(lnk)
    pk = np.array([powspec.eval(cosmo, 0.0, ki) for ki in k])
    radial_sq = np.array([radial(ki) ** 2 for ki in k])

    return 2.0 / np.pi * np.trapezoid(k**3 * pk * radial_sq, lnk)


# (shape, ell, rtol), at the library's default tolerances.
#
# The tolerances differ per shape because the achieved accuracy does, and by two
# orders of magnitude. scaled-abstol floors W(k) against its own *peak*, so how
# much of the contributing k-range it discards depends on how broad the window
# is. Measured at ell = 8 against the quadrature below:
#
#     shape      support (Mpc)   floor 1e-4   1e-5      1e-6
#     gauss      300-2700        -4.3e-06     -2.1e-07  +5.2e-08
#     lensing    50-3000         -5.6e-05     -1.7e-06  -1.4e-06
#     power_exp  50-4000         -6.1e-04     -2.4e-05  +1.4e-06
#
# Nothing is wrong with the solver: a broad window carries weight over more of
# the k-range, so the same peak-relative floor throws away more of it. The
# numbers below are what the default delivers, not a target -- tightening the
# floor moves them, at the cost documented in NcXcorKernel:scaled-abstol.
CL_CASES = [
    ("gauss", 2, 1.0e-5),
    ("gauss", 8, 1.0e-5),
    ("gauss", 32, 1.0e-5),
    ("student_t", 2, 1.0e-4),
    ("student_t", 8, 1.0e-4),
    ("lensing", 8, 1.0e-4),
    ("power_exp", 8, 1.0e-3),
    ("tophat", 2, 1.0e-3),
]


@pytest.mark.parametrize("shape, ell, rtol", CL_CASES)
def test_cl_matches_quadrature_from_the_definition(
    cosmology: Cosmology, shape: str, ell: int, rtol: float
) -> None:
    """The non-Limber C_ell is the integral the window's formula defines.

    This is what an analytic kernel buys: every other kernel is spline-backed,
    so agreement can only ever be between two library routines sharing the same
    reconstruction. Here the right-hand side is written from the closed form,
    which also pins the whole unit convention -- Mpc-normalized window, the
    R_H factor that carries it into the internal radial measure, and the 2/pi.
    """
    kernel = _cl_kernel(cosmology, shape)

    solver = Nc.XcorSolver.new()
    kernel_id = solver.register_kernel(kernel)
    solver.request_cl(kernel_id, kernel_id, ell, ell)
    solver.plan_blocks(64)
    solver.solve(
        Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.KERNEL_EXACT),
        cosmology.cosmo,
    )
    got = solver.get_result(0).get(0)

    expected = _cl_reference(kernel, cosmology, shape, ell, n_k=512, per_lobe=2)

    assert_allclose(got, expected, rtol=rtol)


def test_limber_and_non_limber_agree_at_high_ell(cosmology: Cosmology) -> None:
    """The two paths share the window, so they must converge as ell grows.

    Limber and non-Limber reach the window through different code -- one
    evaluates it at chi = (ell+1/2)/k, the other integrates it against j_ell --
    so this catches a normalization that is right in only one of them. The
    Gaussian carries the check: the agreement is a property of the shared
    window, and the top-hat's non-Limber solve at this ell is far more
    expensive without testing anything further.
    """
    lmin, lmax = 200, 204
    out = {}

    for tier in (0, -1):
        kernel = _gauss(cosmology, tier)
        solver = Nc.XcorSolver.new()
        kernel_id = solver.register_kernel(kernel)
        solver.request_cl(kernel_id, kernel_id, lmin, lmax)
        solver.plan_blocks(64)
        solver.solve(
            Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.KERNEL_EXACT),
            cosmology.cosmo,
        )
        out[tier] = np.array(solver.get_result(0).dup_array())

    assert_allclose(out[0], out[-1], rtol=1.0e-3)


# --- power-law tails --------------------------------------------------------


def test_student_t_matches_scipy(cosmology: Cosmology) -> None:
    """eval_W is the truncated Student-t, renormalized over what is kept."""
    kernel = _student_t(cosmology)
    chi_min, chi_max = kernel.get_support()
    chi = np.linspace(chi_min, chi_max, 257)

    got = np.array([kernel.eval_W(c) for c in chi])

    assert_allclose(got, _window("student_t")(chi), rtol=1.0e-13)


def test_student_t_reports_the_mass_it_discards(cosmology: Cosmology) -> None:
    """The tail cut is not negligible here, so the class states how big it is.

    A Gaussian at four sigma leaves out 6e-5 of itself and nobody has to care.
    A nu = 2 profile at six scales leaves out 2.7%, and a benchmark that does
    not say so is not reproducible.
    """
    kernel = _student_t(cosmology)
    chi_min, chi_max = kernel.get_support()

    kept = student_t.cdf((chi_max - ST_MEAN) / ST_SCALE, ST_NU) - student_t.cdf(
        (chi_min - ST_MEAN) / ST_SCALE, ST_NU
    )

    assert_allclose(kernel.get_tail_mass(), 1.0 - kept, rtol=1.0e-12)
    assert kernel.get_tail_mass() > 1.0e-2


def test_student_t_tail_is_not_exponential(cosmology: Cosmology) -> None:
    """The tail still carries weight where a Gaussian's has vanished.

    This is the property the shape exists for. Comparing the drop from three
    widths to six: a Gaussian loses six orders of magnitude, the nu = 2 profile
    less than one. Stated as a ratio rather than as an asymptotic slope, which
    the finite support never reaches -- at six scales the "+1" in
    (1 + t^2/nu) still contributes, so the local slope is -2.67, not -3.
    """
    kernel = _student_t(cosmology)
    near, far = ST_MEAN + 3.0 * ST_SCALE, ST_MEAN + 6.0 * ST_SCALE

    drop = kernel.eval_W(far) / kernel.eval_W(near)
    gaussian_drop = np.exp(-0.5 * 6.0**2) / np.exp(-0.5 * 3.0**2)

    assert drop > 0.1
    assert drop / gaussian_drop > 1.0e4


# --- multimodal -------------------------------------------------------------


@pytest.mark.parametrize(
    "means, sigmas, n_comps",
    [
        (MULTI_MEAN, MULTI_SIGMA, 1),
        (MULTI_DISJOINT_MEAN, MULTI_DISJOINT_SIGMA, 2),
    ],
    ids=["overlapping", "disjoint"],
)
def test_multi_groups_by_support_not_by_peak(
    cosmology: Cosmology, means: list[float], sigmas: list[float], n_comps: int
) -> None:
    """Bumps that run together share a component; separated ones do not.

    The split exists so that every boundary of the window is a boundary of an
    integration domain. Cutting a continuous stretch at each bump's own n-sigma
    would do the opposite: manufacture steps inside one.
    """
    kernel = _multi(cosmology, means, sigmas)

    assert kernel.get_n_comps() == n_comps

    intervals = [kernel.get_comp_support(i) for i in range(n_comps)]

    for lo, hi in intervals:
        assert hi > lo

    for (_, hi), (lo, _) in zip(intervals[:-1], intervals[1:]):
        assert lo > hi  # a real gap, not an arbitrary cut


@pytest.mark.parametrize(
    "means, sigmas",
    [(MULTI_MEAN, MULTI_SIGMA), (MULTI_DISJOINT_MEAN, MULTI_DISJOINT_SIGMA)],
    ids=["overlapping", "disjoint"],
)
def test_multi_is_normalized_across_its_components(
    cosmology: Cosmology, means: list[float], sigmas: list[float]
) -> None:
    """The components together carry unit integral, gaps and all."""
    kernel = _multi(cosmology, means, sigmas)

    total = 0.0

    for i in range(kernel.get_n_comps()):
        lo, hi = kernel.get_comp_support(i)
        total += quad(lambda c, i=i: kernel.eval_W_comp(i, c), lo, hi, limit=400)[0]

    assert_allclose(total, 1.0, rtol=1.0e-12)


def test_multi_vanishes_in_the_gap(cosmology: Cosmology) -> None:
    """Between disjoint components the window is exactly zero, not merely small."""
    kernel = _multi(cosmology, MULTI_DISJOINT_MEAN, MULTI_DISJOINT_SIGMA)
    _, first_hi = kernel.get_comp_support(0)
    second_lo, _ = kernel.get_comp_support(1)

    assert kernel.eval_W(0.5 * (first_hi + second_lo)) == 0.0


# --- scale dependence -------------------------------------------------------


def test_kdep_is_scale_free_below_the_transition() -> None:
    """Well below k_t the factor is one, so the kernel is unchanged there."""
    alpha, k_t, chi_ref = 0.3, 0.05, 1500.0
    kdep = Nc.XcorKernelRadialKDepGrowth.new(alpha, k_t, chi_ref)

    chi = np.array([600.0, 1500.0, 2400.0])
    k = np.array([1.0e-4, 1.0e-2, 0.05, 1.0, 10.0])
    cc, kk = np.meshgrid(chi, k, indexing="ij")

    sigma = (kk / k_t) ** 2 / (1.0 + (kk / k_t) ** 2)
    expected = np.exp(-alpha * sigma * (chi_ref - cc) / chi_ref)
    got = np.array([[kdep.eval(c, kv) for kv in k] for c in chi])

    assert_allclose(got, expected, rtol=1.0e-14)

    # Scale-free well below the transition, saturating well above it, and
    # exactly one at the reference distance whatever k is.
    assert_allclose(kdep.eval(600.0, 1.0e-6), 1.0, rtol=1.0e-9)
    assert_allclose(kdep.eval(1500.0, 1.0), 1.0, rtol=1.0e-15)
    assert_allclose(kdep.eval(600.0, 1.0e4), np.exp(-alpha * 0.6), rtol=1.0e-7)
    assert kdep.eval(600.0, 1.0) < 1.0  # suppressed nearer the observer


def test_zero_amplitude_kdep_changes_nothing(cosmology: Cosmology) -> None:
    """alpha = 0 must reproduce the kernel without any factor, exactly.

    The null test for a capability slot: carrying a scale dependence that does
    nothing has to be indistinguishable from not carrying one, or every later
    comparison against the scale-free case is confounded.
    """
    lmin, lmax = 2, 8
    out = {}

    kdeps: list[tuple[str, typing.Optional[Nc.XcorKernelRadialKDep]]] = [
        ("none", None),
        ("zero", Nc.XcorKernelRadialKDepGrowth.new(0.0, 0.05, 1500.0)),
    ]

    for label, kdep in kdeps:
        kernel = _gauss_with_kdep(cosmology, kdep)

        solver = Nc.XcorSolver.new()
        kernel_id = solver.register_kernel(kernel)
        solver.request_cl(kernel_id, kernel_id, lmin, lmax)
        solver.plan_blocks(64)
        solver.solve(
            Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.KERNEL_EXACT),
            cosmology.cosmo,
        )
        out[label] = np.array(solver.get_result(0).dup_array())

    assert np.array_equal(out["none"], out["zero"])


def test_kdep_makes_the_integrand_non_separable(cosmology: Cosmology) -> None:
    """A non-zero factor changes C_ell, and not by an overall constant.

    Without it the radial integral has the same shape at every k, so any
    rescaling would be k-independent and the ratio across ell would be flat.
    A k-dependent factor has to break that.
    """
    lmin, lmax = 2, 8
    out = {}

    kdeps: list[tuple[str, typing.Optional[Nc.XcorKernelRadialKDep]]] = [
        ("off", None),
        ("on", Nc.XcorKernelRadialKDepGrowth.new(0.3, 0.05, 3000.0)),
    ]

    for label, kdep in kdeps:
        kernel = _gauss_with_kdep(cosmology, kdep)

        solver = Nc.XcorSolver.new()
        kernel_id = solver.register_kernel(kernel)
        solver.request_cl(kernel_id, kernel_id, lmin, lmax)
        solver.plan_blocks(64)
        solver.solve(
            Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.KERNEL_EXACT),
            cosmology.cosmo,
        )
        out[label] = np.array(solver.get_result(0).dup_array())

    ratio = out["on"] / out["off"]

    assert np.all(ratio != 1.0)
    assert ratio.max() / ratio.min() - 1.0 > 1.0e-3  # not an overall rescaling


# --- the remaining physically motivated shapes ------------------------------


def test_power_exp_matches_its_closed_form(cosmology: Cosmology) -> None:
    """The dn/dz family: power-law rise, stretched-exponential fall.

    The normalization is an incomplete gamma rather than a quadrature, so this
    also checks that substitution against numerical integration of the same
    unnormalized shape.
    """
    kernel = _power_exp(cosmology)
    chi = np.linspace(*kernel.get_support(), 401)

    assert_allclose(
        [kernel.eval_W(c) for c in chi], _window("power_exp")(chi), rtol=1.0e-11
    )
    assert_allclose(
        quad(kernel.eval_W, *kernel.get_support(), limit=600)[0], 1.0, rtol=1.0e-11
    )


def test_power_exp_is_skewed(cosmology: Cosmology) -> None:
    """It is not a Gaussian in disguise: the peak sits left of the mean.

    Skewness is the point of the shape -- a real redshift distribution rises as
    a power law and falls as a stretched exponential, and a symmetric window
    does not exercise the same thing.
    """
    kernel = _power_exp(cosmology)
    chi = np.linspace(*kernel.get_support(), 4001)
    w = np.array([kernel.eval_W(c) for c in chi])

    peak = chi[np.argmax(w)]
    mean = np.trapezoid(chi * w, chi)

    assert peak < mean


def _tophat_smooth_unnormalized(chi: np.ndarray) -> np.ndarray:
    """erf(tu) - erf(tl), written so neither tail cancels.

    Both error functions tend to the same limit in either tail, so the naive
    difference is accurate only to the machine epsilon of *one*, not of their
    difference -- 1e-8 relative at 6 sigma. That is above the tolerance the
    Levin expansion asks of the window, so the reference has to be written the
    way the library writes it.
    """
    s2 = np.sqrt(2.0) * TS_SIGMA
    tu = (TS_UPPER - chi) / s2
    tl = (TS_LOWER - chi) / s2

    return np.where(
        (tu >= 0.0) & (tl >= 0.0),
        erfc(tl) - erfc(tu),
        np.where((tu <= 0.0) & (tl <= 0.0), erfc(-tu) - erfc(-tl), erf(tu) - erf(tl)),
    )


def test_tophat_smooth_matches_its_closed_form(cosmology: Cosmology) -> None:
    """A difference of error functions, normalized in closed form."""
    kernel = _tophat_smooth(cosmology)
    chi_min, chi_max = kernel.get_support()
    chi = np.linspace(chi_min, chi_max, 401)

    unnormalized = _tophat_smooth_unnormalized(chi)
    norm = quad(
        lambda c: _tophat_smooth_unnormalized(np.array(c)),
        chi_min,
        chi_max,
        limit=600,
    )[0]

    assert_allclose([kernel.eval_W(c) for c in chi], unnormalized / norm, rtol=1.0e-12)
    assert_allclose(
        quad(kernel.eval_W, chi_min, chi_max, limit=600)[0], 1.0, rtol=1.0e-12
    )


def test_tophat_smooth_becomes_the_tophat(cosmology: Cosmology) -> None:
    """As sigma shrinks it reproduces the sharp bin, height and all.

    The two shapes bracket the interesting regime, so the smoothed one has to
    reduce to the sharp one rather than merely resemble it.
    """
    sharp = _tophat_smooth(cosmology, chi_sigma=1.0e-3)

    assert_allclose(sharp.eval_W(1500.0), 1.0 / (TS_UPPER - TS_LOWER), rtol=1.0e-9)


def test_lensing_matches_its_closed_form(cosmology: Cosmology) -> None:
    """Efficiency of a top-hat source bin, against the geometry written here."""
    kernel = _lensing(cosmology)
    chi_min, chi_max = kernel.get_support()
    chi = np.linspace(chi_min, chi_max, 4001)[:-1]  # both vanish at the far edge

    assert_allclose(
        [kernel.eval_W(c) for c in chi], _window("lensing")(chi), rtol=1.0e-11
    )

    # Split at the kink: the piecewise form changes where the lens enters the
    # source bin, and a single quadrature would be judging its own resolution.
    total = (
        quad(kernel.eval_W, chi_min, LENS_SRC_LOWER, limit=400)[0]
        + quad(kernel.eval_W, LENS_SRC_LOWER, chi_max, limit=400)[0]
    )

    assert_allclose(total, 1.0, rtol=1.0e-12)


def test_lensing_is_broad_and_peaks_in_front_of_its_sources(
    cosmology: Cosmology,
) -> None:
    """An integral of the source window, so it is smooth and displaced.

    This is what makes the shape worth having: it peaks around half the source
    distance, is continuous where the source distribution has hard edges, and
    vanishes at the far edge -- none of which a bin-like window does.
    """
    kernel = _lensing(cosmology)
    chi = np.linspace(*kernel.get_support(), 4001)
    w = np.array([kernel.eval_W(c) for c in chi])

    assert kernel.eval_W(LENS_SRC_UPPER) == 0.0
    assert chi[np.argmax(w)] < LENS_SRC_LOWER  # peaks in front of the sources

    # Continuous across the kink where the lens enters the source bin.
    eps = 1.0e-4
    assert_allclose(
        kernel.eval_W(LENS_SRC_LOWER - eps),
        kernel.eval_W(LENS_SRC_LOWER + eps),
        rtol=1.0e-6,
    )


# --- the spec a kernel reports back ----------------------------------------


def test_every_shape_reports_the_spec_it_was_built_from(cosmology: Cosmology) -> None:
    """Each accessor returns exactly what construction was given.

    Not bookkeeping: these accessors are how a benchmark specification is read
    back off a kernel, so a value that does not survive construction would make
    a recorded spec disagree with the object that produced the numbers.
    """
    gauss = _gauss(cosmology)
    assert (gauss.get_chi_mean(), gauss.get_chi_sigma(), gauss.get_n_sigma()) == (
        CHI_MEAN,
        CHI_SIGMA,
        N_SIGMA,
    )

    tophat = _tophat(cosmology)
    assert (tophat.get_chi_lower(), tophat.get_chi_upper()) == (CHI_LOWER, CHI_UPPER)

    student = _student_t(cosmology)
    assert (
        student.get_chi_mean(),
        student.get_chi_scale(),
        student.get_nu(),
        student.get_n_scale(),
    ) == (ST_MEAN, ST_SCALE, ST_NU, ST_NSCALE)

    power = _power_exp(cosmology)
    assert (power.get_chi_scale(), power.get_alpha(), power.get_beta()) == (
        PE_SCALE,
        PE_ALPHA,
        PE_BETA,
    )

    smooth = _tophat_smooth(cosmology)
    assert (
        smooth.get_chi_lower(),
        smooth.get_chi_upper(),
        smooth.get_chi_sigma(),
    ) == (TS_LOWER, TS_UPPER, TS_SIGMA)

    lensing = _lensing(cosmology)
    assert (
        lensing.get_chi_source_lower(),
        lensing.get_chi_source_upper(),
    ) == (LENS_SRC_LOWER, LENS_SRC_UPPER)


def test_multi_reports_its_bumps(cosmology: Cosmology) -> None:
    """The bump vectors survive construction, and n_bumps is not n_comps.

    The distinction matters: overlapping bumps share a component, so a spec
    that read the component count back as the bump count would be wrong.
    """
    kernel = _multi(cosmology, MULTI_MEAN, MULTI_SIGMA)

    assert kernel.get_n_bumps() == len(MULTI_MEAN)
    assert kernel.get_n_sigma() == MULTI_NSIGMA
    assert kernel.get_n_bumps() != kernel.get_n_comps()  # these two overlap

    for peeked, expected in (
        (kernel.peek_chi_mean(), MULTI_MEAN),
        (kernel.peek_chi_sigma(), MULTI_SIGMA),
        (kernel.peek_weight(), MULTI_WEIGHT),
    ):
        assert_allclose([peeked.get(i) for i in range(peeked.len())], expected)


def test_kdep_reports_its_parameters() -> None:
    """The scale dependence is part of the spec too.

    Through the accessors and through the GObject properties: serialization
    reads the latter, so the two have to agree.
    """
    kdep = Nc.XcorKernelRadialKDepGrowth.new(0.3, 0.05, 1500.0)

    assert (
        kdep.get_amplitude(),
        kdep.get_k_transition(),
        kdep.get_chi_ref(),
    ) == (0.3, 0.05, 1500.0)

    assert (
        kdep.props.amplitude,
        kdep.props.k_transition,
        kdep.props.chi_ref,
    ) == (0.3, 0.05, 1500.0)


def test_a_kernel_reports_the_scale_dependence_it_carries(cosmology: Cosmology) -> None:
    """peek_kdep returns the attached factor, or None when there is none.

    The distinction is part of the spec: no factor at all and a factor that
    evaluates to one are different states, and a reader of a stored kernel has
    to be able to tell them apart.
    """
    kdep = Nc.XcorKernelRadialKDepGrowth.new(0.3, 0.05, 1500.0)

    assert _gauss_with_kdep(cosmology, None).peek_kdep() is None
    assert _gauss_with_kdep(cosmology, kdep).peek_kdep() is not None


@pytest.mark.parametrize("shape", ["gauss", "tophat", "multi_disjoint"])
def test_limber_z_is_the_window_in_hubble_radius_units(
    cosmology: Cosmology, shape: str
) -> None:
    """The Limber-in-z path evaluates the same window, per unit Hubble radius.

    That path integrates over z rather than chi, so it carries the R_H factor
    of the change of variable and no ell-dependent prefactor. It also has one
    shared domain for the whole kernel, unlike the non-Limber path: the gap of
    a disjoint multimodal window is zero there only if each component refuses
    to contribute outside its own interval, so that shape is included.

    It carries one factor beyond the window. This class puts the growth in W and
    pairs it with P(k, 0), while the Limber integrand multiplies the two kernels
    by P(k, z); each kernel therefore returns sqrt(P(k,0)/P(k,z)) so the product
    restores P(k, 0). The factor is 1 for a power spectrum without growth.
    """
    if shape == "multi_disjoint":
        kernel = _multi(cosmology, MULTI_DISJOINT_MEAN, MULTI_DISJOINT_SIGMA)
    else:
        kernel = _any_shape(cosmology, shape)

    cosmo = cosmology.cosmo
    rh = cosmo.RH_Mpc()
    zmin, zmax, _ = kernel.get_z_range()

    ps = kernel.peek_powspec()

    for z in np.linspace(zmin, zmax, 16):
        chi = cosmology.dist.comoving(cosmo, z) * rh
        got = kernel.eval_limber_z_full(cosmo, z, cosmology.dist, 8)
        k = (8 + 0.5) / chi
        growth = np.sqrt(ps.eval(cosmo, 0.0, k) / ps.eval(cosmo, z, k))

        assert_allclose(got, rh * kernel.eval_W(chi) * growth, rtol=1.0e-12)

    assert kernel.eval_limber_z_prefactor(cosmo, 8) == 1.0


def test_an_analytic_kernel_is_a_window_not_a_survey(cosmology: Cosmology) -> None:
    """It has one observable, no observable parameters, and no noise.

    A kernel defined by a formula has no survey behind it, so the noise it adds
    to a spectrum is nothing at all.
    """
    kernel = _gauss(cosmology)

    assert kernel.obs_len() == 1
    assert kernel.obs_params_len() == 0

    cl = [1.0, 2.0, 3.0, 4.0]
    cl_noisy = Ncm.Vector.new_array(cl)

    kernel.add_noise(Ncm.Vector.new_array(cl), cl_noisy, 0)

    assert_allclose(cl_noisy.dup_array(), cl)


def test_limber_agrees_with_non_limber_at_high_ell():
    """The two branches must describe the same kernel.

    NcXcorKernelRadial carries the growth in W and pairs it with P(k, 0), while
    the Limber integrand multiplies the two kernels by P(k, z). Getting that
    wrong is invisible to self-convergence and to any comparison at fixed
    method: it is an ell-independent factor, so it reads as a normalization.
    Limber is an excellent approximation at these ell, so requiring the branches
    to agree is what pins it.
    """
    cosmo = Nc.HICosmoDEXcdm.new()
    dist = Nc.Distance.new(5.0)
    dist.prepare(cosmo)

    for growth in (Ncm.PowspecAnalyticGrowth.LCDM, Ncm.PowspecAnalyticGrowth.NONE):
        ps = Ncm.PowspecAnalytic.new(Ncm.PowspecAnalyticShape.BBKS, growth)

        for ell in (100, 200):
            # One integrator per multipole: new_full(lmin, lmax) allocates over
            # the whole contiguous range.
            sbi = Ncm.SBesselIntegratorLevin.new(ell, ell)
            k_nl = Nc.XcorKernelAnalyticGauss.new_full(
                dist, ps, 1500.0, 300.0, 4.0, sbi
            )
            k_nl.set_l_limber(-1)
            k_nl.prepare(cosmo)

            k_li = Nc.XcorKernelAnalyticGauss.new(dist, ps, 1500.0, 300.0, 4.0)
            k_li.set_l_limber(0)
            k_li.prepare(cosmo)

            def _cl(kernel, method):
                xc = Nc.Xcor.new(dist, ps, method)
                xc.prepare(cosmo)
                v = Ncm.Vector.new(1)
                xc.compute(kernel, kernel, cosmo, ell, ell, v)
                return v.get(0)

            non_limber = _cl(k_nl, Nc.XcorMethod.KERNEL_EXACT)
            limber = _cl(k_li, Nc.XcorMethod.LIMBER_Z_CUBATURE)

            assert limber == pytest.approx(non_limber, rel=2.0e-3)
