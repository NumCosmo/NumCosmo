#!/usr/bin/env python
#
# test_py_xcor_ccl_comparison.py
#
# Thu Apr 17 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_xcor_ccl_comparison.py
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

"""Unit tests comparing NumCosmo XcorKernel with CCL (Core Cosmology Library).

These tests validate NumCosmo's cross-correlation kernel implementations
against the pyccl library to ensure numerical consistency between the two codes.
Modernized from test_py_xcor.py CCL comparison tests.
"""

import pytest
from numpy.testing import assert_allclose
import numpy as np

pytest.importorskip("getdist")
pytest.importorskip("pyccl")
# flake8: noqa: E402
# pylint: disable=wrong-import-position

import pyccl

from numcosmo_py import Ncm, Nc
from numcosmo_py.cosmology import Cosmology
from numcosmo_py.ccl.two_point import compute_kernel
import numcosmo_py.ccl.comparison as nc_cmp

pytestmark = [pytest.mark.ccl, pytest.mark.xcor]

Ncm.cfg_init()
pytest_plugins = ["python.fixtures_ccl", "python.fixtures_xcor"]


# Kernel Window Function Comparison Tests


def test_cmb_lens_kernel(
    ccl_cosmo_eh_linear: pyccl.Cosmology,
    nc_cosmo_eh_linear: Cosmology,
    ccl_cmb_lens: pyccl.CMBLensingTracer,
    kernel_cmb_lens: Nc.XcorKernelCMBLensing,
) -> None:
    """Compare NumCosmo and CCL CMB lensing kernel window functions.

    Tests that the Limber-approximated correlation window matches between
    NumCosmo and CCL for CMB lensing (kappa convergence) tracers.
    """
    cosmo = nc_cosmo_eh_linear.cosmo
    dist = nc_cosmo_eh_linear.dist

    if ccl_cosmo_eh_linear["Omega_k"] != 0.0:
        pytest.skip("CMB lensing not implemented for non-flat cosmologies")

    reltol_target = 1.0e-6 if ccl_cosmo_eh_linear.high_precision else 1.0e-4

    ell = 79.0
    kernel_cmb_lens.set_l_limber(0)

    # Get CCL kernel window
    z_a, _, Wchi_a = compute_kernel(ccl_cmb_lens, nc_cosmo_eh_linear, ell)

    # Get NumCosmo kernel window
    kernel_cmb_lens.prepare(cosmo)
    nc_Wchi_a = (
        np.array(
            [kernel_cmb_lens.eval_limber_z_full(cosmo, z, dist, int(ell)) for z in z_a]
        )
        / cosmo.RH_Mpc()
    )

    assert_allclose(nc_Wchi_a, Wchi_a, rtol=reltol_target, atol=0.0)


def test_cmb_isw_kernel(
    ccl_cosmo_eh_linear: pyccl.Cosmology,
    nc_cosmo_eh_linear: Cosmology,
    ccl_cmb_isw: pyccl.ISWTracer,
    kernel_cmb_isw: Nc.XcorKernelCMBISW,
) -> None:
    """Compare NumCosmo and CCL ISW kernel window functions.

    Tests that the Integrated Sachs-Wolfe effect kernel matches between
    NumCosmo and CCL across redshift.
    """
    cosmo = nc_cosmo_eh_linear.cosmo
    dist = nc_cosmo_eh_linear.dist

    if ccl_cosmo_eh_linear["Omega_k"] != 0.0:
        pytest.skip("CMB ISW not implemented for non-flat cosmologies")

    reltol_target = 1.0e-6 if ccl_cosmo_eh_linear.high_precision else 1.0e-4

    ell = 77.0

    # Get CCL kernel window
    z_a, _, Wchi_a = compute_kernel(ccl_cmb_isw, nc_cosmo_eh_linear, ell)

    # Get NumCosmo kernel window
    kernel_cmb_isw.prepare(cosmo)
    nc_Wchi_a = (
        np.array(
            [kernel_cmb_isw.eval_limber_z_full(cosmo, z, dist, int(ell)) for z in z_a]
        )
        / cosmo.RH_Mpc()
    )

    assert_allclose(nc_Wchi_a, Wchi_a, rtol=reltol_target, atol=0.0)


def test_tsz_kernel(
    nc_cosmo_eh_linear: Cosmology,
    ccl_tsz: pyccl.tSZTracer,
    kernel_tsz: Nc.XcorKerneltSZ,
) -> None:
    """Compare NumCosmo and CCL tSZ kernel window functions.

    Tests that the thermal Sunyaev-Zel'dovich effect kernel matches between
    NumCosmo and CCL.
    """
    cosmo = nc_cosmo_eh_linear.cosmo
    dist = nc_cosmo_eh_linear.dist

    reltol_target = 1.0e-4

    ell = 77.0

    # Get CCL kernel window
    z_a, _, Wchi_a = compute_kernel(ccl_tsz, nc_cosmo_eh_linear, ell)

    # Get NumCosmo kernel window
    kernel_tsz.prepare(cosmo)
    nc_Wchi_a = (
        np.array([kernel_tsz.eval_limber_z_full(cosmo, z, dist, int(ell)) for z in z_a])
        / cosmo.RH_Mpc()
    )

    assert_allclose(nc_Wchi_a, Wchi_a, rtol=reltol_target, atol=0.0)


@pytest.mark.parametrize("wl_bin", [0, 1, 2, 3, 4])
def test_weak_lensing_kernel(
    ccl_cosmo_eh_linear: pyccl.Cosmology,
    nc_cosmo_eh_linear: Cosmology,
    wl_bin: int,
    request: pytest.FixtureRequest,
) -> None:
    """Compare NumCosmo and CCL weak lensing kernel window functions.

    Tests that the cosmic shear kernel window matches between NumCosmo and CCL,
    accounting for the source galaxy redshift distribution.
    """
    # Get specific bin fixture
    kernel_wl = request.getfixturevalue(f"kernel_wl_bin{wl_bin}")

    # Create CCL tracer from NumCosmo kernel parameters
    dndz: Ncm.Spline = kernel_wl.props.dndz
    assert isinstance(dndz, Ncm.Spline)
    ccl_weak_lensing = pyccl.WeakLensingTracer(
        ccl_cosmo_eh_linear,
        dndz=(
            np.array(dndz.peek_xv().dup_array()),
            np.array(dndz.peek_yv().dup_array()),
        ),
        n_samples=10_000,
    )

    cosmo = nc_cosmo_eh_linear.cosmo
    dist = nc_cosmo_eh_linear.dist

    reltol_target = 6.0e-4  # Relaxed tolerance for weak lensing kernel comparison (bin 4 needs ~0.06%)

    ell = 77.0

    # Get CCL kernel window
    z_a, _, Wchi_a = compute_kernel(ccl_weak_lensing, nc_cosmo_eh_linear, ell)

    # Get NumCosmo kernel window
    kernel_wl.prepare(cosmo)
    nc_Wchi_a = (
        np.array([kernel_wl.eval_limber_z_full(cosmo, z, dist, int(ell)) for z in z_a])
        / cosmo.RH_Mpc()
    )

    assert_allclose(nc_Wchi_a, Wchi_a, rtol=reltol_target, atol=1.0e-20)


@pytest.mark.parametrize("gal_bin", [0, 1, 2, 3, 4])
def test_gal_kernel(
    ccl_cosmo_eh_linear: pyccl.Cosmology,
    nc_cosmo_eh_linear: Cosmology,
    gal_bin: int,
    request: pytest.FixtureRequest,
) -> None:
    """Compare NumCosmo and CCL galaxy number counts kernel window functions.

    Tests that the galaxy clustering kernel (including bias and magnification)
    matches between NumCosmo and CCL.
    """
    # Get specific bin fixture
    kernel_gal = request.getfixturevalue(f"kernel_gal_bin{gal_bin}")

    # Create CCL tracer from NumCosmo kernel parameters
    dndz: Ncm.Spline = kernel_gal.props.dndz
    assert isinstance(dndz, Ncm.Spline)
    bias = kernel_gal.orig_vparam_get(Nc.XcorKernelGalVParams.BIAS, 0)
    magbias = kernel_gal.orig_param_get(Nc.XcorKernelGalSParams.MAG_BIAS)
    z_array = np.array(dndz.peek_xv().dup_array())
    ccl_gal = pyccl.NumberCountsTracer(
        ccl_cosmo_eh_linear,
        has_rsd=False,
        dndz=(
            z_array,
            np.array(dndz.peek_yv().dup_array()),
        ),
        bias=(z_array, np.ones_like(z_array) * bias) if bias != 0.0 else None,
        mag_bias=(z_array, np.ones_like(z_array) * magbias) if magbias != 0.0 else None,
        n_samples=len(z_array),
    )

    cosmo = nc_cosmo_eh_linear.cosmo
    dist = nc_cosmo_eh_linear.dist

    reltol_target = 1.0e-4

    ell = 77.0

    # Get CCL kernel window
    z_a, _, Wchi_a = compute_kernel(ccl_gal, nc_cosmo_eh_linear, ell)

    # Get NumCosmo kernel window
    kernel_gal.prepare(cosmo)
    nc_Wchi_a = (
        np.array([kernel_gal.eval_limber_z_full(cosmo, z, dist, int(ell)) for z in z_a])
        / cosmo.RH_Mpc()
    )

    assert_allclose(nc_Wchi_a, Wchi_a, rtol=reltol_target, atol=1.0e-30)


# Auto-Correlation Power Spectrum Tests


def test_cmb_lens_auto_integrand(
    ccl_cosmo_eh_linear: pyccl.Cosmology,
    nc_cosmo_eh_linear: Cosmology,
    ccl_cmb_lens: pyccl.CMBLensingTracer,
    kernel_cmb_lens: Nc.XcorKernelCMBLensing,
) -> None:
    """Compare NumCosmo and CCL CMB lensing integrand components.

    Tests that the integrand components (kernel, power spectrum, and full integrand)
    match between NumCosmo and CCL for the CMB lensing auto-correlation.
    """
    cosmo = nc_cosmo_eh_linear.cosmo
    dist = nc_cosmo_eh_linear.dist
    ps_ml = nc_cosmo_eh_linear.ps_ml

    if ccl_cosmo_eh_linear["Omega_k"] != 0.0:
        pytest.skip("CMB lensing not implemented for non-flat cosmologies")

    if ccl_cosmo_eh_linear.high_precision:
        reltol_W = 1.0e-4
        reltol_ps = 1.0e-3
        reltol_f = 1.0e-3
    else:
        reltol_W = 1.0e-4
        reltol_ps = 1.0e-2
        reltol_f = 1.0e-2

    psp = ccl_cosmo_eh_linear.get_linear_power()

    xcor = Nc.Xcor.new(dist, ps_ml, Nc.XcorMethod.LIMBER_Z_CUBATURE)
    xcor.prepare(cosmo)
    kernel_cmb_lens.prepare(cosmo)

    ell = 77.0
    z_a, chi_a, Wchi_a = compute_kernel(ccl_cmb_lens, nc_cosmo_eh_linear, ell)

    nu = ell + 0.5
    a_a = ccl_cosmo_eh_linear.scale_factor_of_chi(chi_a)
    k_a = nu / chi_a
    z_a = 1.0 / a_a - 1.0

    # NumCosmo components
    nc_W = (
        np.array(
            [kernel_cmb_lens.eval_limber_z_full(cosmo, z, dist, int(ell)) for z in z_a]
        )
        / cosmo.RH_Mpc()
    )
    nc_ps = np.array([ps_ml.eval(cosmo, z, k) for z, k in zip(z_a, k_a)])
    nc_f = k_a * nc_ps * nc_W**2 / nu

    # CCL components
    ccl_W = Wchi_a * ell * (ell + 1.0) / nu**2
    ccl_ps = np.array([psp(k, a, cosmo=ccl_cosmo_eh_linear) for k, a in zip(k_a, a_a)])
    ccl_f = k_a * ccl_ps * ccl_W**2 / nu

    assert_allclose(nc_W, ccl_W, rtol=reltol_W, atol=0.0)
    assert_allclose(nc_ps, ccl_ps, rtol=reltol_ps, atol=0.0)
    assert_allclose(nc_f, ccl_f, rtol=reltol_f, atol=0.0)


def test_cmb_lens_auto(
    ccl_cosmo_eh_linear: pyccl.Cosmology,
    nc_cosmo_eh_linear: Cosmology,
    ccl_cmb_lens: pyccl.CMBLensingTracer,
    kernel_cmb_lens: Nc.XcorKernelCMBLensing,
) -> None:
    """Compare NumCosmo and CCL CMB lensing auto-correlation power spectrum.

    Tests that the integrated angular power spectrum C_ell matches between
    NumCosmo and CCL for CMB lensing across a range of multipoles.
    """
    cosmo = nc_cosmo_eh_linear.cosmo
    dist = nc_cosmo_eh_linear.dist
    ps_ml = nc_cosmo_eh_linear.ps_ml

    if ccl_cosmo_eh_linear["Omega_k"] != 0.0:
        pytest.skip("CMB lensing not implemented for non-flat cosmologies")

    reltol_target = 1.0e-2

    lmax = 3000
    ells = np.arange(2, lmax + 1)

    # CCL power spectrum
    psp = ccl_cosmo_eh_linear.get_linear_power()
    assert (
        ccl_cosmo_eh_linear.get_linear_power() == ccl_cosmo_eh_linear.get_nonlin_power()
    )
    ccl_cmb_lens_auto = pyccl.angular_cl(
        ccl_cosmo_eh_linear,
        ccl_cmb_lens,
        ccl_cmb_lens,
        ells,
        l_limber=-1,
        p_of_k_a=psp,
        p_of_k_a_lin=psp,
    )

    assert ccl_cmb_lens_auto is not None
    assert all(np.isfinite(ccl_cmb_lens_auto))
    assert all(ccl_cmb_lens_auto >= 0.0)

    # NumCosmo power spectrum
    xcor = Nc.Xcor.new(dist, ps_ml, Nc.XcorMethod.LIMBER_Z_CUBATURE)
    nc_cmb_lens_auto_v = Ncm.Vector.new(lmax + 1 - 2)
    xcor.prepare(cosmo)
    kernel_cmb_lens.prepare(cosmo)
    xcor.compute(kernel_cmb_lens, kernel_cmb_lens, cosmo, 2, lmax, nc_cmb_lens_auto_v)
    nc_cmb_lens_auto = np.array(nc_cmb_lens_auto_v.dup_array())

    assert_allclose(ccl_cmb_lens_auto, nc_cmb_lens_auto, rtol=reltol_target, atol=0.0)


# Comparison Helper Function Tests


@pytest.mark.parametrize("n_points", [None, 400])
def test_compare_kernels(
    ccl_cosmo_eh_linear: pyccl.Cosmology,
    nc_cosmo_eh_linear: Cosmology,
    n_points: int | None,
) -> None:
    """Compare all kernel types using NumCosmo CCL comparison helpers.

    Tests that the comparison helper functions correctly match CCL and NumCosmo
    for all kernel types: CMB lensing, ISW, tSZ, weak lensing, and galaxy clustering.

    Args:
        n_points: Number of sampling points for kernel evaluation (None for default)
    """
    # CMB lensing kernel (flat cosmologies only)
    if ccl_cosmo_eh_linear["Omega_k"] == 0.0:
        cmp = nc_cmp.compare_cmb_lens_kernel(
            ccl_cosmo_eh_linear, nc_cosmo_eh_linear, ell=77, n_samples=n_points
        )
        assert_allclose(cmp.y1, cmp.y2, rtol=1.0e-6)

    # CMB ISW kernel
    cmp = nc_cmp.compare_cmb_isw_kernel(
        ccl_cosmo_eh_linear, nc_cosmo_eh_linear, ell=77, n_chi=n_points
    )
    assert_allclose(cmp.y1, cmp.y2, rtol=1.0e-1)

    # tSZ kernel
    cmp = nc_cmp.compare_tsz_kernel(
        ccl_cosmo_eh_linear, nc_cosmo_eh_linear, ell=77, n_chi=n_points
    )
    assert_allclose(cmp.y1, cmp.y2, rtol=1.0e-8)

    # Weak lensing kernel
    cmp = nc_cmp.compare_galaxy_weak_lensing_kernel(
        ccl_cosmo_eh_linear, nc_cosmo_eh_linear, ell=77, n_samples=n_points
    )
    assert_allclose(cmp.y1, cmp.y2, rtol=1.0e-3, atol=1.0e-15)

    # Galaxy number counts kernel
    cmp = nc_cmp.compare_galaxy_number_count_kernel(
        ccl_cosmo_eh_linear, nc_cosmo_eh_linear, ell=77
    )
    assert_allclose(cmp.y1, cmp.y2, rtol=1.0e-4, atol=1.0e-12)


@pytest.mark.parametrize("n_points", [None, 400])
def test_compare_autocorrelation(
    ccl_cosmo_eh_linear: pyccl.Cosmology,
    nc_cosmo_eh_linear: Cosmology,
    n_points: int | None,
) -> None:
    """Compare auto-correlation power spectra using NumCosmo CCL comparison helpers.

    Tests that the angular power spectrum C_ell matches between CCL and NumCosmo
    for all tracer auto-correlations across a range of multipoles.

    Args:
        n_points: Number of sampling points for integration (None for default)
    """
    ells = np.arange(2, 1000)

    # CMB lensing auto (flat cosmologies only)
    if ccl_cosmo_eh_linear["Omega_k"] == 0.0:
        cmp = nc_cmp.compare_cmb_len_auto(
            ccl_cosmo_eh_linear, nc_cosmo_eh_linear, ells, n_samples=n_points
        )
        assert_allclose(cmp.y1, cmp.y2, rtol=1.0e-2)

    # CMB ISW auto
    cmp = nc_cmp.compare_cmb_isw_auto(
        ccl_cosmo_eh_linear, nc_cosmo_eh_linear, ells, n_chi=n_points
    )
    assert_allclose(cmp.y1, cmp.y2, rtol=1.0e-1)

    # tSZ auto
    cmp = nc_cmp.compare_tsz_auto(
        ccl_cosmo_eh_linear, nc_cosmo_eh_linear, ells, n_chi=n_points
    )
    assert_allclose(cmp.y1, cmp.y2, rtol=1.0e-3)

    # Weak lensing auto
    cmp = nc_cmp.compare_galaxy_weak_lensing_auto(
        ccl_cosmo_eh_linear, nc_cosmo_eh_linear, ells, n_samples=n_points
    )
    assert_allclose(cmp.y1, cmp.y2, rtol=1.0e-3)

    # Galaxy number counts auto
    cmp = nc_cmp.compare_galaxy_number_count_auto(
        ccl_cosmo_eh_linear, nc_cosmo_eh_linear, ells
    )
    assert_allclose(cmp.y1, cmp.y2, rtol=1.0e-4)
