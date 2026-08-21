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

import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.integrate import quad
from scipy.special import spherical_jn
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
    kernel = Nc.XcorKernelAnalyticGauss(
        dist=cosmology.dist,
        powspec=cosmology.ps_ml,
        chi_mean=CHI_MEAN,
        chi_sigma=CHI_SIGMA,
        n_sigma=N_SIGMA,
        integrator=_integrator(),
    )
    kernel.set_l_limber(l_limber)
    kernel.prepare(cosmology.cosmo)

    return kernel


def _tophat(cosmology: Cosmology, l_limber: int = -1) -> Nc.XcorKernelAnalyticTophat:
    """Top-hat window, prepared and in the requested Limber tier."""
    kernel = Nc.XcorKernelAnalyticTophat(
        dist=cosmology.dist,
        powspec=cosmology.ps_ml,
        chi_lower=CHI_LOWER,
        chi_upper=CHI_UPPER,
        integrator=_integrator(),
    )
    kernel.set_l_limber(l_limber)
    kernel.prepare(cosmology.cosmo)

    return kernel


def _student_t(
    cosmology: Cosmology, l_limber: int = -1
) -> Nc.XcorKernelAnalyticStudentT:
    """Power-law-tailed window, prepared and in the requested Limber tier."""
    kernel = Nc.XcorKernelAnalyticStudentT(
        dist=cosmology.dist,
        powspec=cosmology.ps_ml,
        chi_mean=ST_MEAN,
        chi_scale=ST_SCALE,
        nu=ST_NU,
        n_scale=ST_NSCALE,
        integrator=_integrator(),
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
    kernel = Nc.XcorKernelAnalyticMulti(
        dist=cosmology.dist,
        powspec=cosmology.ps_ml,
        chi_mean=Ncm.Vector.new_array(means),
        chi_sigma=Ncm.Vector.new_array(sigmas),
        weight=Ncm.Vector.new_array(MULTI_WEIGHT),
        n_sigma=MULTI_NSIGMA,
        integrator=_integrator(),
    )
    kernel.set_l_limber(l_limber)
    kernel.prepare(cosmology.cosmo)

    return kernel


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


def test_gauss_window_vanishes_outside_support(cosmology: Cosmology) -> None:
    """The truncation is part of the definition, not a tolerance."""
    kernel = _gauss(cosmology)
    chi_min, chi_max = kernel.get_support()

    assert kernel.eval_W(chi_min - 1.0e-3) == 0.0
    assert kernel.eval_W(chi_max + 1.0e-3) == 0.0
    assert kernel.eval_W(chi_min) > 0.0
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


@pytest.mark.parametrize("shape", ["gauss", "tophat"])
def test_serialization_roundtrip(cosmology: Cosmology, shape: str) -> None:
    """The kernel is a registered, serializable object: a spec survives a file."""
    kernel = _gauss(cosmology) if shape == "gauss" else _tophat(cosmology)
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)

    dup = ser.dup_obj(kernel)

    assert type(dup) is type(kernel)
    assert dup.get_support() == kernel.get_support()
    assert_allclose(
        [dup.eval_W(c) for c in np.linspace(*kernel.get_support(), 32)],
        [kernel.eval_W(c) for c in np.linspace(*kernel.get_support(), 32)],
        rtol=1.0e-15,
    )


GX, GW = np.polynomial.legendre.leggauss(12)


def _gauss_with_kdep(
    cosmology: Cosmology, kdep: typing.Optional[Nc.XcorKernelAnalyticKDep]
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


def _cl_kernel(cosmology: Cosmology, shape: str) -> Nc.XcorKernelAnalytic:
    """The kernel a C_ell case names."""
    builder = {"gauss": _gauss, "tophat": _tophat, "student_t": _student_t}[shape]

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
    kernel: Nc.XcorKernelAnalytic,
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
    k_min, k_max = kernel.get_eval(cosmo, ell).get_range()
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


# (shape, ell, rtol). The Gaussian reference is converged -- it does not move
# between (n_k, per_lobe) = (512, 2) and (1024, 3) -- so the budget there is the
# library's own tolerance. The top-hat's discontinuous window makes I_ell(k)
# decay only like 1/k, so the reference's outer trapezoid is the loose end, and
# its non-Limber solve costs ~20x the Gaussian's at the same ell; it is checked
# at ell = 2 rather than swept.
CL_CASES = [
    ("gauss", 2, 1.0e-5),
    ("gauss", 8, 1.0e-5),
    ("gauss", 32, 1.0e-5),
    ("student_t", 2, 1.0e-4),
    ("student_t", 8, 1.0e-4),
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
        Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.KERNEL_FIXED),
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
            Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.KERNEL_FIXED),
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
    kdep = Nc.XcorKernelAnalyticKDepGrowth.new(alpha, k_t, chi_ref)

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

    kdeps: list[tuple[str, typing.Optional[Nc.XcorKernelAnalyticKDep]]] = [
        ("none", None),
        ("zero", Nc.XcorKernelAnalyticKDepGrowth.new(0.0, 0.05, 1500.0)),
    ]

    for label, kdep in kdeps:
        kernel = _gauss_with_kdep(cosmology, kdep)

        solver = Nc.XcorSolver.new()
        kernel_id = solver.register_kernel(kernel)
        solver.request_cl(kernel_id, kernel_id, lmin, lmax)
        solver.plan_blocks(64)
        solver.solve(
            Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.KERNEL_FIXED),
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

    kdeps: list[tuple[str, typing.Optional[Nc.XcorKernelAnalyticKDep]]] = [
        ("off", None),
        ("on", Nc.XcorKernelAnalyticKDepGrowth.new(0.3, 0.05, 3000.0)),
    ]

    for label, kdep in kdeps:
        kernel = _gauss_with_kdep(cosmology, kdep)

        solver = Nc.XcorSolver.new()
        kernel_id = solver.register_kernel(kernel)
        solver.request_cl(kernel_id, kernel_id, lmin, lmax)
        solver.plan_blocks(64)
        solver.solve(
            Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.KERNEL_FIXED),
            cosmology.cosmo,
        )
        out[label] = np.array(solver.get_result(0).dup_array())

    ratio = out["on"] / out["off"]

    assert np.all(ratio != 1.0)
    assert ratio.max() / ratio.min() - 1.0 > 1.0e-3  # not an overall rescaling
