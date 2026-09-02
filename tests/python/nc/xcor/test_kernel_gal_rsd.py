#!/usr/bin/env python
#
# test_kernel_gal_rsd.py
#
# Tue Sep 01 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_kernel_gal_rsd.py
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

"""Redshift-space distortions in NcXcorKernelGal against CCL.

The dorsd property adds the linear Kaiser term as a component with kernel
-f(z) dn/dz weighted by j_l''(k chi), following CCL's NumberCountsTracer
convention (transfer -f, der_bessel = 2). The exact tier integrates it with
the derivative-weighted Levin solve; the kernel-space Limber tier peak-
approximates the j_l and j_{l+1} pieces of the recurrence separately.
"""

import gc

import numpy as np
import pytest
from numpy.testing import assert_allclose

pytest.importorskip("pyccl")

# pylint: disable=wrong-import-position
import pyccl

from numcosmo_py import Nc, Ncm
from numcosmo_py.ccl.nc_ccl import create_nc_obj, CCLParams
from numcosmo_py.ccl import two_point

pytestmark = [pytest.mark.ccl, pytest.mark.xcor, pytest.mark.xdist_group("ccl_rsd")]

Ncm.cfg_init()

Z_LEN = 1000
MU, SIGMA = 0.5, 0.1
BIAS = 1.5


@pytest.fixture(name="ccl_cosmo", scope="module")
def fixture_ccl_cosmo():
    """CCL cosmology matching the NumCosmo one, linear power spectrum."""
    CCLParams.set_high_prec_params()
    return pyccl.Cosmology(
        Omega_c=0.25,
        Omega_b=0.05,
        h=0.7,
        n_s=0.96,
        sigma8=0.8,
        transfer_function="eisenstein_hu",
        matter_power_spectrum="linear",
    )


@pytest.fixture(name="cosmology", scope="module")
def fixture_cosmology(ccl_cosmo):
    """NumCosmo cosmology matching the CCL one."""
    return create_nc_obj(ccl_cosmo)


@pytest.fixture(name="dndz_arrays", scope="module")
def fixture_dndz_arrays():
    """Gaussian redshift bin shared by both codes."""
    z_a = np.linspace(0.0, 1.5, Z_LEN)
    nz_a = np.exp(-0.5 * ((z_a - MU) / SIGMA) ** 2)
    return z_a, nz_a


def _nc_gal_kernel(cosmology, dndz_arrays, *, dorsd: bool, lmax: int):
    z_a, nz_a = dndz_arrays
    dndz = Ncm.Spline.new_array(
        Ncm.SplineCubicNotaknot.new(), z_a.tolist(), nz_a.tolist(), True
    )
    kernel = Nc.XcorKernelGal(
        dist=cosmology.dist,
        powspec=cosmology.ps_ml,
        bparam_length=1,
        dndz=dndz,
        domagbias=False,
        dorsd=dorsd,
        integrator=Ncm.SBesselIntegratorLevin.new(2, lmax),
    )
    kernel.orig_vparam_set(Nc.XcorKernelGalVParams.BIAS, 0, BIAS)
    return kernel


def _ccl_gal_tracer(ccl_cosmo, dndz_arrays, *, has_rsd: bool):
    z_a, nz_a = dndz_arrays
    return pyccl.NumberCountsTracer(
        ccl_cosmo,
        has_rsd=has_rsd,
        dndz=(z_a, nz_a),
        bias=(z_a, np.ones_like(z_a) * BIAS),
        n_samples=Z_LEN,
    )


def _nc_cl(cosmology, kernel, ells, l_limber, kernel_2=None):
    kernel_2 = kernel if kernel_2 is None else kernel_2
    for k in (kernel, kernel_2):
        k.set_l_limber(l_limber)
        k.prepare(cosmology.cosmo)
    lmin, lmax = int(ells[0]), int(ells[-1])
    res = Ncm.Vector.new(lmax - lmin + 1)
    xcor = Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.KERNEL_EXACT)
    xcor.prepare(cosmology.cosmo)
    xcor.compute(kernel, kernel_2, cosmology.cosmo, lmin, lmax, res)
    return np.array(res.dup_array())[np.asarray(ells) - lmin]


def test_rsd_auto_exact_matches_ccl_limber_at_high_ell(
    ccl_cosmo, cosmology, dndz_arrays
) -> None:
    """Exact-tier C_ell with RSD against CCL's Limber method at high ell.

    CCL's Limber is the trustworthy reference here, and only at high ell:
    pyccl's FKEM non-Limber RSD measured 38% low at ell = 2 and 11% low at
    ell = 10 on this bin, while this exact tier and the CCL-kernel bridge
    (numcosmo_py.ccl.two_point.angular_cl, which integrates CCL's own kernels
    and growth with the derivative-weighted Levin solve) agree there to 1e-6.
    At ell ~ 100-300 the Limber approximation itself is good to a few 1e-4 on
    this bin, which sets the tolerance.
    """
    ells = np.array([100, 300])

    kernel = _nc_gal_kernel(cosmology, dndz_arrays, dorsd=True, lmax=int(ells[-1]))
    got = _nc_cl(cosmology, kernel, ells, l_limber=-1)

    tracer = _ccl_gal_tracer(ccl_cosmo, dndz_arrays, has_rsd=True)
    psp = ccl_cosmo.get_linear_power()
    ref = pyccl.angular_cl(
        ccl_cosmo, tracer, tracer, ells, l_limber=0, p_of_k_a_lin=psp, p_of_k_a=psp
    )

    assert_allclose(got, ref, rtol=2.0e-3)


def test_rsd_auto_exact_matches_ccl_kernel_bridge_at_low_ell(
    ccl_cosmo, cosmology, dndz_arrays
) -> None:
    """Exact-tier C_ell with RSD against the CCL-kernel bridge where RSD matters.

    numcosmo_py.ccl.two_point.angular_cl integrates CCL's own tracer kernels,
    transfer (-f for the RSD piece, der_bessel = 2) and growth with the
    derivative-weighted Levin solve, so the two sides share the solver and
    nothing else: NumCosmo's dn/dz, bias, distances and growth on one side,
    CCL's tables on the other. The bridge's support cut at 1e-4 of the kernel
    peak keeps it at a few seconds; measured deviation 3.5e-6 at ell = 2
    (2.1e-6 with the cut at 1e-5), so the tolerance below is the cut, not
    the solve.
    """
    ells = np.array([2, 5, 10])

    kernel = _nc_gal_kernel(cosmology, dndz_arrays, dorsd=True, lmax=int(ells[-1]))
    got = _nc_cl(cosmology, kernel, ells, l_limber=-1)

    tracer = _ccl_gal_tracer(ccl_cosmo, dndz_arrays, has_rsd=True)
    ref = two_point.angular_cl(cosmology, tracer, tracer, ells, support_tol=1.0e-4)

    assert_allclose(got, ref, rtol=2.0e-5)


def test_rsd_cross_spectrum_matches_ccl_kernel_bridge(
    ccl_cosmo, cosmology, dndz_arrays
) -> None:
    """Cross spectrum of an RSD bin with a second, non-overlapping bin.

    The cross of a density-plus-RSD kernel with a plain density kernel is
    negative here (the bins do not overlap), so the sign and the mixed
    j_l x j_l'' term are both checked; a second case switches RSD on in
    both bins. Measured deviations vs the bridge are at or below 1.3e-5.
    """
    ells = np.array([2, 5, 10])
    z2 = np.linspace(0.0, 2.0, Z_LEN)
    nz2 = np.exp(-0.5 * ((z2 - 0.9) / SIGMA) ** 2)

    kernel_1 = _nc_gal_kernel(cosmology, dndz_arrays, dorsd=True, lmax=int(ells[-1]))
    tracer_1 = _ccl_gal_tracer(ccl_cosmo, dndz_arrays, has_rsd=True)

    for has_rsd in (False, True):
        kernel_2 = _nc_gal_kernel(
            cosmology, (z2, nz2), dorsd=has_rsd, lmax=int(ells[-1])
        )
        tracer_2 = _ccl_gal_tracer(ccl_cosmo, (z2, nz2), has_rsd=has_rsd)

        got = _nc_cl(cosmology, kernel_1, ells, l_limber=-1, kernel_2=kernel_2)
        ref = two_point.angular_cl(
            cosmology, tracer_1, tracer_2, ells, support_tol=1.0e-4
        )

        assert np.all(ref < 0.0)
        assert_allclose(got, ref, rtol=5.0e-5)


def test_rsd_auto_limber_matches_ccl(ccl_cosmo, cosmology, dndz_arrays) -> None:
    """Kernel-space Limber C_ell with RSD against CCL's Limber method.

    Both sides use the same two-point peak approximation of the j_l''
    recurrence, so they must track each other more closely than either tracks
    the exact result.
    """
    ells = np.array([50, 100, 200, 400])

    kernel = _nc_gal_kernel(cosmology, dndz_arrays, dorsd=True, lmax=int(ells[-1]))
    got = _nc_cl(cosmology, kernel, ells, l_limber=0)

    tracer = _ccl_gal_tracer(ccl_cosmo, dndz_arrays, has_rsd=True)
    psp = ccl_cosmo.get_linear_power()
    ref = pyccl.angular_cl(
        ccl_cosmo, tracer, tracer, ells, l_limber=0, p_of_k_a_lin=psp, p_of_k_a=psp
    )

    assert_allclose(got, ref, rtol=5.0e-4)


def test_rsd_raises_the_low_ell_auto_spectrum(cosmology, dndz_arrays) -> None:
    """The Kaiser term boosts the low-ell galaxy auto spectrum."""
    ells = np.array([2, 5, 10])

    with_rsd = _nc_cl(
        cosmology,
        _nc_gal_kernel(cosmology, dndz_arrays, dorsd=True, lmax=10),
        ells,
        l_limber=-1,
    )
    without_rsd = _nc_cl(
        cosmology,
        _nc_gal_kernel(cosmology, dndz_arrays, dorsd=False, lmax=10),
        ells,
        l_limber=-1,
    )

    assert np.all(with_rsd > without_rsd)
    # at ell = 2 the boost is large, not a rounding-level difference
    assert with_rsd[0] > 1.2 * without_rsd[0]


def test_rsd_growth_rate_matches_ccl(ccl_cosmo, cosmology) -> None:
    """NumCosmo's growth rate agrees with CCL's, which the C_ell tests fold in."""
    gf = Nc.GrowthFunc.new()
    gf.prepare(cosmology.cosmo)

    z_a = np.linspace(0.0, 1.5, 16)
    f_nc = []
    for z in z_a:
        d, dd_dz = gf.eval_both(cosmology.cosmo, float(z))
        f_nc.append(-(1.0 + z) * dd_dz / d)

    f_ccl = ccl_cosmo.growth_rate(1.0 / (1.0 + z_a))

    assert_allclose(f_nc, f_ccl, rtol=1.0e-6)


def test_rsd_component_list_and_properties(cosmology, dndz_arrays) -> None:
    """dorsd adds one component of Bessel-derivative order 2 beside the density one."""
    kernel = _nc_gal_kernel(cosmology, dndz_arrays, dorsd=True, lmax=10)
    assert kernel.get_property("dorsd") is True

    comps = kernel.get_component_list()
    orders = sorted(c.get_bessel_deriv() for c in comps)
    assert orders == [0, 2]

    rsd = next(c for c in comps if c.get_bessel_deriv() == 2)
    assert rsd.get_property("bessel-deriv") == 2

    kernel_plain = _nc_gal_kernel(cosmology, dndz_arrays, dorsd=False, lmax=10)
    assert kernel_plain.get_property("dorsd") is False
    assert [c.get_bessel_deriv() for c in kernel_plain.get_component_list()] == [0]

    # dropping the kernel disposes the RSD component and its data
    del comps, rsd, kernel
    gc.collect()


def test_rsd_kernel_survives_serialization(cosmology, dndz_arrays) -> None:
    """A serialized copy keeps dorsd, rebuilds the RSD component and gives the same C_ell.

    NcXcorSolver hands each OpenMP thread NcmSerialize-duplicated objects, so
    a construct-only property that did not survive the round trip would
    silently drop the Kaiser term there.
    """
    ells = np.array([2, 5, 10])
    kernel = _nc_gal_kernel(cosmology, dndz_arrays, dorsd=True, lmax=int(ells[-1]))

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    copy = ser.dup_obj(kernel)

    assert copy.get_property("dorsd") is True
    assert sorted(c.get_bessel_deriv() for c in copy.get_component_list()) == [0, 2]

    original = _nc_cl(cosmology, kernel, ells, l_limber=-1)
    duplicate = _nc_cl(cosmology, copy, ells, l_limber=-1)

    assert_allclose(duplicate, original, rtol=1.0e-12)


def test_rsd_through_solver_matches_compute(cosmology, dndz_arrays) -> None:
    """NcXcorSolver's block path reproduces nc_xcor_compute() with the RSD component.

    One 8-ell block, exact tier; the solver runs its OpenMP team over blocks
    with per-block integrator clones, so this covers the deriv state living on
    the integrator rather than on the shared kernel.
    """
    lmin, lmax = 2, 9
    kernel = _nc_gal_kernel(cosmology, dndz_arrays, dorsd=True, lmax=lmax)
    kernel.set_l_limber(-1)
    kernel.prepare(cosmology.cosmo)

    xcor = Nc.Xcor.new(cosmology.dist, cosmology.ps_ml, Nc.XcorMethod.KERNEL_EXACT)
    xcor.prepare(cosmology.cosmo)

    expected = Ncm.Vector.new(lmax - lmin + 1)
    xcor.compute(kernel, kernel, cosmology.cosmo, lmin, lmax, expected)

    solver = Nc.XcorSolver.new()
    kid = solver.register_kernel(kernel)
    solver.request_cl(kid, kid, lmin, lmax)
    solver.plan_blocks(8)
    solver.solve(xcor, cosmology.cosmo)

    got = np.array(solver.get_result(0).dup_array())
    assert_allclose(got, np.array(expected.dup_array()), rtol=1.0e-10)


def test_deriv_one_limber_matches_exact_at_high_ell(cosmology, dndz_arrays) -> None:
    """A j_l'-weighted component: the two-point Limber peak formula vs the exact solve.

    No physical kernel uses order 1 yet, so this reconfigures the density
    component; the point is that both tiers implement the same integral.
    """
    ell = 200

    kernel = _nc_gal_kernel(cosmology, dndz_arrays, dorsd=False, lmax=ell)
    comp = kernel.get_component_list()[0]
    comp.set_bessel_deriv(1)
    assert comp.get_bessel_deriv() == 1

    exact = _nc_cl(cosmology, kernel, np.array([ell]), l_limber=-1)
    limber = _nc_cl(cosmology, kernel, np.array([ell]), l_limber=0)

    assert exact[0] > 0.0
    # measured: -8.4e-3 at ell = 200 (the Limber error of this suppressed integral)
    assert_allclose(limber, exact, rtol=2.0e-2)
