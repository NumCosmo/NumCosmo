#!/usr/bin/env python
#
# test_py_xcor.py
#
# Thu Aug 01 11:45:10 2024
# Copyright  2024  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_xcor.py
# Copyright (C) 2024 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

"""Unit tests for NumCosmo powwer-spectra."""

import pytest

import numpy as np
from numpy.testing import assert_allclose

import pyccl

import numcosmo_py.cosmology as ncpy
from numcosmo_py import Ncm, Nc

from .fixtures_ccl import (  # pylint: disable=unused-import # noqa: F401
    fixture_k_a,
    fixture_z_a,
    fixture_ccl_cosmo_eh_linear,
    fixture_ccl_cosmo_eh_halofit,
    fixture_nc_cosmo_eh_linear,
    fixture_nc_cosmo_eh_halofit,
)
from .fixtures_xcor import (  # pylint: disable=unused-import # noqa: F401
    fixture_ccl_cmb_lens,
    fixture_nc_cmb_lens,
    fixture_ccl_cmb_isw,
    fixture_nc_cmb_isw,
    fixture_ccl_gal,
    fixture_nc_gal,
    fixture_ccl_weak_lensing,
    fixture_nc_weak_lensing,
)

Ncm.cfg_init()


# Testing the NumCosmo tracers observables


def test_cmb_lens_obs(nc_cmb_lens: Nc.XcorLimberKernelCMBLensing) -> None:
    """Check that CMB lensing tracer has the correct number of observables."""
    assert nc_cmb_lens is not None
    assert isinstance(nc_cmb_lens, Nc.XcorLimberKernelCMBLensing)
    assert nc_cmb_lens.obs_len() == 1
    assert nc_cmb_lens.obs_params_len() == 0


def test_cmb_isw_obs(nc_cmb_isw: Nc.XcorLimberKernelCMBISW) -> None:
    """Check that CMB ISW tracer has the correct number of observables."""
    assert nc_cmb_isw is not None
    assert isinstance(nc_cmb_isw, Nc.XcorLimberKernelCMBISW)
    assert nc_cmb_isw.obs_len() == 1
    assert nc_cmb_isw.obs_params_len() == 0


def test_gal_obs(nc_gal: Nc.XcorLimberKernelGal) -> None:
    """Check that galaxy tracer has the correct number of observables."""
    assert nc_gal is not None
    assert isinstance(nc_gal, Nc.XcorLimberKernelGal)
    assert nc_gal.obs_len() == 2
    assert nc_gal.obs_params_len() == 1


def test_weak_lensing_obs(nc_weak_lensing: Nc.XcorLimberKernelWeakLensing) -> None:
    """Check that weak lensing tracer has the correct number of observables."""
    assert nc_weak_lensing is not None
    assert isinstance(nc_weak_lensing, Nc.XcorLimberKernelWeakLensing)
    assert nc_weak_lensing.obs_len() == 2
    assert nc_weak_lensing.obs_params_len() == 1


# Testing the NumCosmo tracers serialization


def test_cmb_lens_serialization(
    nc_cosmo_eh_linear: ncpy.Cosmology, nc_cmb_lens: Nc.XcorLimberKernelCMBLensing
) -> None:
    """Check that CMB lensing tracer can be serialized."""
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    nc_cmb_lens_dup = ser.dup_obj(nc_cmb_lens)
    assert nc_cmb_lens_dup is not None
    assert nc_cmb_lens_dup is not nc_cmb_lens
    assert isinstance(nc_cmb_lens_dup, Nc.XcorLimberKernelCMBLensing)

    cosmo = nc_cosmo_eh_linear.cosmo
    dist = nc_cosmo_eh_linear.dist
    nc_cmb_lens.prepare(cosmo)
    nc_cmb_lens_dup.prepare(cosmo)

    assert_allclose(
        nc_cmb_lens.eval_full(cosmo, 0.0, dist, 2),
        nc_cmb_lens_dup.eval_full(cosmo, 0.0, dist, 2),
    )


def test_cmb_isw_serialization(
    nc_cosmo_eh_linear: ncpy.Cosmology, nc_cmb_isw: Nc.XcorLimberKernelCMBISW
) -> None:
    """Check that CMB ISW tracer can be serialized."""
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    nc_cmb_isw_dup = ser.dup_obj(nc_cmb_isw)
    assert nc_cmb_isw_dup is not None
    assert nc_cmb_isw_dup is not nc_cmb_isw
    assert isinstance(nc_cmb_isw_dup, Nc.XcorLimberKernelCMBISW)

    cosmo = nc_cosmo_eh_linear.cosmo
    dist = nc_cosmo_eh_linear.dist
    nc_cmb_isw.prepare(cosmo)
    nc_cmb_isw_dup.prepare(cosmo)

    assert_allclose(
        nc_cmb_isw.eval_full(cosmo, 0.0, dist, 2),
        nc_cmb_isw_dup.eval_full(cosmo, 0.0, dist, 2),
    )


def test_gal_serialization(
    nc_cosmo_eh_linear: ncpy.Cosmology, nc_gal: Nc.XcorLimberKernelGal
) -> None:
    """Check that galaxy tracer can be serialized."""
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    nc_gal_dup = ser.dup_obj(nc_gal)
    assert nc_gal_dup is not None
    assert nc_gal_dup is not nc_gal
    assert isinstance(nc_gal_dup, Nc.XcorLimberKernelGal)

    cosmo = nc_cosmo_eh_linear.cosmo
    dist = nc_cosmo_eh_linear.dist

    nc_gal.prepare(cosmo)
    nc_gal_dup.prepare(cosmo)

    assert_allclose(
        nc_gal.eval_full(cosmo, 0.0, dist, 2),
        nc_gal_dup.eval_full(cosmo, 0.0, dist, 2),
    )


def test_weak_lensing_serialization(
    nc_cosmo_eh_linear: ncpy.Cosmology, nc_weak_lensing: Nc.XcorLimberKernelWeakLensing
) -> None:
    """Check that weak lensing tracer can be serialized."""
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    nc_wl_dup = ser.dup_obj(nc_weak_lensing)
    assert nc_wl_dup is not None
    assert nc_wl_dup is not nc_weak_lensing
    assert isinstance(nc_wl_dup, Nc.XcorLimberKernelWeakLensing)

    cosmo = nc_cosmo_eh_linear.cosmo
    dist = nc_cosmo_eh_linear.dist

    nc_weak_lensing.prepare(cosmo)
    nc_wl_dup.prepare(cosmo)

    assert_allclose(
        nc_weak_lensing.eval_full(cosmo, 0.0, dist, 2),
        nc_wl_dup.eval_full(cosmo, 0.0, dist, 2),
    )


# Testing the NumCosmo tracers noise


def test_cmb_lens_noise(nc_cmb_lens: Nc.XcorLimberKernelCMBLensing) -> None:
    """Check that CMB lensing tracer has the correct noise."""
    vp1 = Ncm.Vector.new(10)
    vp2 = Ncm.Vector.new(10)
    vp1.set_all(1.0)
    nc_cmb_lens.add_noise(vp1, vp2, 5)
    assert_allclose(vp2.dup_array(), np.arange(6, 16), atol=0.0)


def test_cmb_isw_noise(nc_cmb_isw: Nc.XcorLimberKernelCMBISW) -> None:
    """Check that CMB ISW tracer has the correct noise."""
    vp1 = Ncm.Vector.new(10)
    vp2 = Ncm.Vector.new(10)
    vp1.set_all(1.0)
    nc_cmb_isw.add_noise(vp1, vp2, 5)
    assert_allclose(vp2.dup_array(), np.arange(6, 16), atol=0.0)


def test_gal_noise(nc_gal: Nc.XcorLimberKernelGal) -> None:
    """Check that galaxy tracer has the correct noise."""
    vp1 = Ncm.Vector.new(10)
    vp2 = Ncm.Vector.new(10)
    vp1.set_all(1.0)
    nc_gal.add_noise(vp1, vp2, 5)
    assert_allclose(vp2.dup_array(), np.ones(10) * 2.234, atol=0.0)


def test_weak_lensing_noise(nc_weak_lensing: Nc.XcorLimberKernelWeakLensing) -> None:
    """Check that weak lensing tracer has the correct noise."""
    vp1 = Ncm.Vector.new(10)
    vp2 = Ncm.Vector.new(10)
    vp1.set_all(1.0)
    nc_weak_lensing.add_noise(vp1, vp2, 5)
    assert_allclose(vp2.dup_array(), np.ones(10) * (7.0) ** 2 / 3.0 + 1.0, atol=0.0)


# Testing the NumCosmo tracers kernels


def test_cmb_lens_kernel(
    ccl_cosmo_eh_linear: pyccl.Cosmology,
    nc_cosmo_eh_linear: ncpy.Cosmology,
    ccl_cmb_lens: pyccl.CMBLensingTracer,
    nc_cmb_lens: Nc.XcorLimberKernelCMBLensing,
) -> None:
    """Compare NumCosmo and CCL correlation windows."""
    cosmo = nc_cosmo_eh_linear.cosmo
    dist = nc_cosmo_eh_linear.dist
    if ccl_cosmo_eh_linear["Omega_k"] != 0.0:
        pytest.skip("CMB lensing not implemented for non-flat cosmologies")
    if ccl_cosmo_eh_linear.high_precision:
        reltol_target: float = 1.0e-7
    else:
        reltol_target = 1.0e-4

    ell = 79.0

    Wchi_list, chi_list = ccl_cmb_lens.get_kernel()
    assert chi_list is not None
    assert Wchi_list is not None
    assert len(chi_list) == 1  # Single tracer
    assert len(Wchi_list) == 1  # Single tracer

    chi_a = np.array(chi_list[0])[1:-1]
    Wchi_a = np.array(Wchi_list[0])[1:-1] * ell * (ell + 1.0) / (ell + 0.5) ** 2
    RH_Mpc = cosmo.RH_Mpc()
    nc_cmb_lens.prepare(cosmo)

    z_array = [dist.inv_comoving(cosmo, chi / RH_Mpc) for chi in chi_a]

    nc_Wchi_a = np.array(
        [
            nc_cmb_lens.eval_full(cosmo, z, dist, int(ell)) * cosmo.E(z) / RH_Mpc
            for z in z_array
        ]
    )
    assert_allclose(nc_Wchi_a, Wchi_a, rtol=reltol_target, atol=0.0)


def test_cmb_lens_auto_integrand(
    ccl_cosmo_eh_linear: pyccl.Cosmology,
    nc_cosmo_eh_linear: ncpy.Cosmology,
    ccl_cmb_lens: pyccl.CMBLensingTracer,
    nc_cmb_lens: Nc.XcorLimberKernelCMBLensing,
) -> None:
    """Compare NumCosmo and CCL transfer functions."""
    cosmo = nc_cosmo_eh_linear.cosmo
    dist = nc_cosmo_eh_linear.dist
    ps_ml = nc_cosmo_eh_linear.ps_ml
    if ccl_cosmo_eh_linear["Omega_k"] != 0.0:
        pytest.skip("CMB lensing not implemented for non-flat cosmologies")
    if ccl_cosmo_eh_linear.high_precision:
        reltol_W: float = 1.0e-6
        reltol_ps = 1.0e-4
        reltol_f = 1.0e-4
    else:
        reltol_W = 1.0e-4
        reltol_ps = 1.0e-2
        reltol_f = 1.0e-2

    psp = ccl_cosmo_eh_linear.get_linear_power()

    xcor = Nc.Xcor.new(dist, ps_ml, Nc.XcorLimberMethod.GSL)
    xcor.prepare(cosmo)
    nc_cmb_lens.prepare(cosmo)

    Wchi_list, chi_list = ccl_cmb_lens.get_kernel()
    assert chi_list is not None
    assert Wchi_list is not None
    assert len(chi_list) == 1  # Single tracer
    assert len(Wchi_list) == 1  # Single tracer

    chi_a = np.array(chi_list[0])[1:-1]
    Wchi_a = np.array(Wchi_list[0])[1:-1]
    RH_Mpc = cosmo.RH_Mpc()

    ell = 77.0
    nu = ell + 0.5
    a_a = ccl_cosmo_eh_linear.scale_factor_of_chi(chi_a)
    k_a = nu / chi_a
    z_a = 1.0 / a_a - 1.0

    nc_W = np.array(
        [
            nc_cmb_lens.eval_full(cosmo, z, dist, int(ell)) * cosmo.E(z) / RH_Mpc
            for z in z_a
        ]
    )
    nc_ps = np.array([ps_ml.eval(cosmo, z, k) for z, k in zip(z_a, k_a)])
    nc_f = k_a * nc_ps * nc_W**2 / nu

    ccl_W = Wchi_a * ell * (ell + 1.0) / nu**2
    ccl_ps = np.array([psp(k, a, cosmo=ccl_cosmo_eh_linear) for k, a in zip(k_a, a_a)])
    ccl_f = k_a * ccl_ps * ccl_W**2 / nu

    assert_allclose(nc_W, ccl_W, rtol=reltol_W, atol=0.0)
    assert_allclose(nc_ps, ccl_ps, rtol=reltol_ps, atol=0.0)
    assert_allclose(nc_f, ccl_f, rtol=reltol_f, atol=0.0)


def test_cmb_lens_auto(
    ccl_cosmo_eh_linear: pyccl.Cosmology,
    nc_cosmo_eh_linear: ncpy.Cosmology,
    ccl_cmb_lens: pyccl.CMBLensingTracer,
    nc_cmb_lens: Nc.XcorLimberKernelCMBLensing,
) -> None:
    """Compare NumCosmo and CCL transfer functions."""
    cosmo = nc_cosmo_eh_linear.cosmo
    dist = nc_cosmo_eh_linear.dist
    ps_ml = nc_cosmo_eh_linear.ps_ml
    if ccl_cosmo_eh_linear["Omega_k"] != 0.0:
        pytest.skip("CMB lensing not implemented for non-flat cosmologies")
    if ccl_cosmo_eh_linear.high_precision:
        reltol_target: float = 1.0e-2
    else:
        reltol_target = 1.0e-2

    lmax = 3000
    ells = np.arange(2, lmax + 1)

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

    xcor = Nc.Xcor.new(dist, ps_ml, Nc.XcorLimberMethod.GSL)
    nc_cmb_lens_auto_v = Ncm.Vector.new(lmax + 1 - 2)
    xcor.prepare(cosmo)
    nc_cmb_lens.prepare(cosmo)
    xcor.limber(nc_cmb_lens, nc_cmb_lens, cosmo, 2, lmax, nc_cmb_lens_auto_v)
    nc_cmb_lens_auto = np.array(nc_cmb_lens_auto_v.dup_array())

    assert_allclose(ccl_cmb_lens_auto, nc_cmb_lens_auto, rtol=reltol_target, atol=0.0)


def test_cmb_isw_kernel(
    ccl_cosmo_eh_linear: pyccl.Cosmology,
    nc_cosmo_eh_linear: ncpy.Cosmology,
    ccl_cmb_isw: pyccl.ISWTracer,
    nc_cmb_isw: Nc.XcorLimberKernelCMBISW,
) -> None:
    """Compare NumCosmo and CCL correlation windows."""
    cosmo = nc_cosmo_eh_linear.cosmo
    dist = nc_cosmo_eh_linear.dist
    if ccl_cosmo_eh_linear["Omega_k"] != 0.0:
        pytest.skip("CMB lensing not implemented for non-flat cosmologies")
    if ccl_cosmo_eh_linear.high_precision:
        reltol_target: float = 1.0e-6
    else:
        reltol_target = 1.0e-4

    RH_Mpc = cosmo.RH_Mpc()
    ell = 77.0

    Wchi_list, chi_list = ccl_cmb_isw.get_kernel()
    assert chi_list is not None
    assert Wchi_list is not None
    assert len(chi_list) == 1  # Single tracer
    assert len(Wchi_list) == 1  # Single tracer

    chi_a = np.array(chi_list[0])[1:-1]
    Wchi_a = np.array(Wchi_list[0])[1:-1] / (ell + 0.5) ** 2
    nc_cmb_isw.prepare(cosmo)

    z_array = [dist.inv_comoving(cosmo, chi / RH_Mpc) for chi in chi_a]

    nc_Wchi_a = np.array(
        [
            nc_cmb_isw.eval_full(cosmo, z, dist, int(ell)) * cosmo.E(z) / RH_Mpc
            for z in z_array
        ]
    )
    assert_allclose(nc_Wchi_a, Wchi_a, rtol=reltol_target, atol=0.0)


def test_weak_lensing_kernel(
    ccl_cosmo_eh_linear: pyccl.Cosmology,
    nc_cosmo_eh_linear: ncpy.Cosmology,
    ccl_weak_lensing: pyccl.WeakLensingTracer,
    nc_weak_lensing: Nc.XcorLimberKernelWeakLensing,
) -> None:
    """Compare NumCosmo and CCL correlation windows."""
    cosmo = nc_cosmo_eh_linear.cosmo
    dist = nc_cosmo_eh_linear.dist
    if ccl_cosmo_eh_linear.high_precision:
        reltol_target: float = 1.0e-4
    else:
        reltol_target = 1.0e-4

    RH_Mpc = cosmo.RH_Mpc()
    ell = 77.0

    Wchi_list, chi_list = ccl_weak_lensing.get_kernel()
    assert chi_list is not None
    assert Wchi_list is not None
    assert len(chi_list) == 1  # Single tracer
    assert len(Wchi_list) == 1  # Single tracer

    chi_a = np.array(chi_list[0])[1:-1]
    Wchi_a = (
        np.array(Wchi_list[0])[1:-1]
        * np.sqrt((ell + 2.0) * (ell + 1.0) * ell * (ell - 1.0))
        / (ell + 0.5) ** 2
    )
    nc_weak_lensing.prepare(cosmo)

    z_array = [dist.inv_comoving(cosmo, chi / RH_Mpc) for chi in chi_a]

    nc_Wchi_a = np.array(
        [
            nc_weak_lensing.eval_full(cosmo, z, dist, int(ell)) * cosmo.E(z) / RH_Mpc
            for z in z_array
        ]
    )
    assert_allclose(nc_Wchi_a, Wchi_a, rtol=reltol_target, atol=1.0e-30)