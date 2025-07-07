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

import itertools as it
import pytest
from pytest_lazy_fixtures import lf

import numpy as np
from numpy.testing import assert_allclose

import pyccl

import numcosmo_py.cosmology as ncpy
from numcosmo_py.ccl.two_point import compute_kernel
from numcosmo_py import Ncm, Nc
import numcosmo_py.ccl.comparison as nc_cmp

from .fixtures_ccl import (  # pylint: disable=unused-import # noqa: F401
    fixture_k_a,
    fixture_z_a,
    fixture_ccl_cosmo_eh_linear,
    fixture_ccl_cosmo_eh_halofit,
    fixture_nc_cosmo_eh_linear,
    fixture_nc_cosmo_eh_halofit,
    fixture_nc_cosmo_default,
)
from .fixtures_xcor import (  # pylint: disable=unused-import # noqa: F401
    fixture_ccl_cmb_lens,
    fixture_nc_cmb_lens,
    fixture_ccl_cmb_isw,
    fixture_nc_cmb_isw,
    fixture_ccl_tsz,
    fixture_nc_tsz,
    fixture_ccl_gal,
    fixture_nc_gal,
    fixture_ccl_weak_lensing,
    fixture_nc_weak_lensing,
)

Ncm.cfg_init()

# Testing the NumCosmo tracers observables


def test_cmb_lens_obs(nc_cmb_lens: Nc.XcorKernelCMBLensing) -> None:
    """Check that CMB lensing tracer has the correct number of observables."""
    assert nc_cmb_lens is not None
    assert isinstance(nc_cmb_lens, Nc.XcorKernelCMBLensing)
    assert nc_cmb_lens.obs_len() == 1
    assert nc_cmb_lens.obs_params_len() == 0


def test_cmb_isw_obs(nc_cmb_isw: Nc.XcorKernelCMBISW) -> None:
    """Check that CMB ISW tracer has the correct number of observables."""
    assert nc_cmb_isw is not None
    assert isinstance(nc_cmb_isw, Nc.XcorKernelCMBISW)
    assert nc_cmb_isw.obs_len() == 1
    assert nc_cmb_isw.obs_params_len() == 0


def test_tsz_obs(nc_tsz: Nc.XcorKerneltSZ) -> None:
    """Check that tSZ tracer has the correct number of observables."""
    assert nc_tsz is not None
    assert isinstance(nc_tsz, Nc.XcorKernelSZ)
    assert nc_tsz.obs_len() == 1
    assert nc_tsz.obs_params_len() == 0


def test_gal_obs(nc_gal: Nc.XcorKernelGal) -> None:
    """Check that galaxy tracer has the correct number of observables."""
    assert nc_gal is not None
    assert isinstance(nc_gal, Nc.XcorKernelGal)
    assert nc_gal.obs_len() == 2
    assert nc_gal.obs_params_len() == 1

    nc_gal.set_fast_update(True)
    assert nc_gal.get_fast_update()

    nc_gal.set_fast_update(False)
    assert not nc_gal.get_fast_update()

    if nc_gal.vparam_len(Nc.XcorKernelGalVParams.BIAS) == 1:
        nc_gal.orig_vparam_set(Nc.XcorKernelGalVParams.BIAS, 0, 3.21)
        nc_gal.set_bias_old(1.2345, 0.9876)

        bias0, bias_old0, noise_bias_old0 = nc_gal.get_bias()

        assert_allclose(bias0, 3.21, atol=0.0)
        assert_allclose(bias_old0, 1.2345, atol=0.0)
        assert_allclose(noise_bias_old0, 0.9876, atol=0.0)


def test_gal_obs_extrapolation(
    nc_cosmo_eh_linear: ncpy.Cosmology, nc_gal: Nc.XcorKernelGal
) -> None:
    """Check that galaxy tracer has the correct number of observables."""
    nc_gal.prepare(nc_cosmo_eh_linear.cosmo)
    z_a = np.array(nc_gal.props.dndz.peek_xv().dup_array())
    nc_gal.set_z_range(z_a[0], z_a[-1] * 2.0, 0.5 * (z_a[-1] + z_a[0]))
    assert np.isfinite(
        nc_gal.eval_full(nc_cosmo_eh_linear.cosmo, z_a[-1], nc_cosmo_eh_linear.dist, 77)
    )

    extra_z_a = np.linspace(z_a[-1], 2.0 * z_a[-1], 10)
    k_a = np.array(
        [
            nc_gal.eval_full(nc_cosmo_eh_linear.cosmo, z, nc_cosmo_eh_linear.dist, 77)
            for z in extra_z_a
        ]
    )
    assert np.isfinite(k_a).all()


def test_weak_lensing_obs(nc_weak_lensing: Nc.XcorKernelWeakLensing) -> None:
    """Check that weak lensing tracer has the correct number of observables."""
    assert nc_weak_lensing is not None
    assert isinstance(nc_weak_lensing, Nc.XcorKernelWeakLensing)
    assert nc_weak_lensing.obs_len() == 2
    assert nc_weak_lensing.obs_params_len() == 1


# Testing the NumCosmo tracers serialization


def test_cmb_lens_serialization(
    nc_cosmo_eh_linear: ncpy.Cosmology, nc_cmb_lens: Nc.XcorKernelCMBLensing
) -> None:
    """Check that CMB lensing tracer can be serialized."""
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    nc_cmb_lens_dup = ser.dup_obj(nc_cmb_lens)
    assert nc_cmb_lens_dup is not None
    assert nc_cmb_lens_dup is not nc_cmb_lens
    assert isinstance(nc_cmb_lens_dup, Nc.XcorKernelCMBLensing)

    cosmo = nc_cosmo_eh_linear.cosmo
    dist = nc_cosmo_eh_linear.dist
    nc_cmb_lens.prepare(cosmo)
    nc_cmb_lens_dup.prepare(cosmo)

    assert_allclose(
        nc_cmb_lens.eval_full(cosmo, 0.0, dist, 2),
        nc_cmb_lens_dup.eval_full(cosmo, 0.0, dist, 2),
    )

    # Prepare for a second time
    nc_cmb_lens.prepare(cosmo)
    assert_allclose(
        nc_cmb_lens.eval_full(cosmo, 0.0, dist, 2),
        nc_cmb_lens_dup.eval_full(cosmo, 0.0, dist, 2),
    )


def test_cmb_isw_serialization(
    nc_cosmo_eh_linear: ncpy.Cosmology, nc_cmb_isw: Nc.XcorKernelCMBISW
) -> None:
    """Check that CMB ISW tracer can be serialized."""
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    nc_cmb_isw_dup = ser.dup_obj(nc_cmb_isw)
    assert nc_cmb_isw_dup is not None
    assert nc_cmb_isw_dup is not nc_cmb_isw
    assert isinstance(nc_cmb_isw_dup, Nc.XcorKernelCMBISW)

    cosmo = nc_cosmo_eh_linear.cosmo
    dist = nc_cosmo_eh_linear.dist
    nc_cmb_isw.prepare(cosmo)
    nc_cmb_isw_dup.prepare(cosmo)

    assert_allclose(
        nc_cmb_isw.eval_full(cosmo, 0.0, dist, 2),
        nc_cmb_isw_dup.eval_full(cosmo, 0.0, dist, 2),
    )

    # Prepare for a second time
    nc_cmb_isw.prepare(cosmo)
    assert_allclose(
        nc_cmb_isw.eval_full(cosmo, 0.0, dist, 2),
        nc_cmb_isw_dup.eval_full(cosmo, 0.0, dist, 2),
    )


def test_tsz_serialization(
    nc_cosmo_eh_linear: ncpy.Cosmology, nc_tsz: Nc.XcorKerneltSZ
) -> None:
    """Check that tSZ tracer can be serialized."""
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    nc_tsz_dup = ser.dup_obj(nc_tsz)
    assert nc_tsz_dup is not None
    assert nc_tsz_dup is not nc_tsz
    assert isinstance(nc_tsz_dup, Nc.XcorKerneltSZ)

    cosmo = nc_cosmo_eh_linear.cosmo
    dist = nc_cosmo_eh_linear.dist
    nc_tsz.prepare(cosmo)
    nc_tsz_dup.prepare(cosmo)

    assert_allclose(
        nc_tsz.eval_full(cosmo, 0.0, dist, 2),
        nc_tsz_dup.eval_full(cosmo, 0.0, dist, 2),
    )

    # Prepare for a second time
    nc_tsz.prepare(cosmo)
    assert_allclose(
        nc_tsz.eval_full(cosmo, 0.0, dist, 2),
        nc_tsz_dup.eval_full(cosmo, 0.0, dist, 2),
    )


def test_gal_serialization(
    nc_cosmo_eh_linear: ncpy.Cosmology, nc_gal: Nc.XcorKernelGal
) -> None:
    """Check that galaxy tracer can be serialized."""
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    nc_gal_dup = ser.dup_obj(nc_gal)
    assert nc_gal_dup is not None
    assert nc_gal_dup is not nc_gal
    assert isinstance(nc_gal_dup, Nc.XcorKernelGal)

    cosmo = nc_cosmo_eh_linear.cosmo
    dist = nc_cosmo_eh_linear.dist

    nc_gal.prepare(cosmo)
    nc_gal_dup.prepare(cosmo)

    assert_allclose(
        nc_gal.eval_full(cosmo, 0.0, dist, 2),
        nc_gal_dup.eval_full(cosmo, 0.0, dist, 2),
    )

    # Prepare for a second time
    nc_gal.prepare(cosmo)
    assert_allclose(
        nc_gal.eval_full(cosmo, 0.0, dist, 2),
        nc_gal_dup.eval_full(cosmo, 0.0, dist, 2),
    )


def test_weak_lensing_serialization(
    nc_cosmo_eh_linear: ncpy.Cosmology, nc_weak_lensing: Nc.XcorKernelWeakLensing
) -> None:
    """Check that weak lensing tracer can be serialized."""
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    nc_wl_dup = ser.dup_obj(nc_weak_lensing)
    assert nc_wl_dup is not None
    assert nc_wl_dup is not nc_weak_lensing
    assert isinstance(nc_wl_dup, Nc.XcorKernelWeakLensing)

    cosmo = nc_cosmo_eh_linear.cosmo
    dist = nc_cosmo_eh_linear.dist

    nc_weak_lensing.prepare(cosmo)
    nc_wl_dup.prepare(cosmo)

    assert_allclose(
        nc_weak_lensing.eval_full(cosmo, 0.0, dist, 2),
        nc_wl_dup.eval_full(cosmo, 0.0, dist, 2),
    )
    # Prepare for a second time
    nc_weak_lensing.prepare(cosmo)
    assert_allclose(
        nc_weak_lensing.eval_full(cosmo, 0.0, dist, 2),
        nc_wl_dup.eval_full(cosmo, 0.0, dist, 2),
    )


# Testing the NumCosmo tracers noise


def test_cmb_lens_noise(nc_cmb_lens: Nc.XcorKernelCMBLensing) -> None:
    """Check that CMB lensing tracer has the correct noise."""
    vp1 = Ncm.Vector.new(10)
    vp2 = Ncm.Vector.new(10)
    vp1.set_all(1.0)
    nc_cmb_lens.add_noise(vp1, vp2, 5)
    assert_allclose(vp2.dup_array(), np.arange(6, 16), atol=0.0)


def test_cmb_isw_noise(nc_cmb_isw: Nc.XcorKernelCMBISW) -> None:
    """Check that CMB ISW tracer has the correct noise."""
    vp1 = Ncm.Vector.new(10)
    vp2 = Ncm.Vector.new(10)
    vp1.set_all(1.0)
    nc_cmb_isw.add_noise(vp1, vp2, 5)
    assert_allclose(vp2.dup_array(), np.arange(6, 16), atol=0.0)


def test_tsz_noise(nc_tsz: Nc.XcorKerneltSZ) -> None:
    """Check that tSZ tracer has the correct noise."""
    vp1 = Ncm.Vector.new(10)
    vp2 = Ncm.Vector.new(10)
    vp1.set_all(1.0)
    nc_tsz.add_noise(vp1, vp2, 5)
    assert_allclose(vp2.dup_array(), np.ones(10) * 1.0, atol=0.0)


def test_gal_noise(nc_gal: Nc.XcorKernelGal) -> None:
    """Check that galaxy tracer has the correct noise."""
    vp1 = Ncm.Vector.new(10)
    vp2 = Ncm.Vector.new(10)
    vp1.set_all(1.0)
    nc_gal.add_noise(vp1, vp2, 5)
    assert_allclose(vp2.dup_array(), np.ones(10) * 2.234, atol=0.0)


def test_weak_lensing_noise(nc_weak_lensing: Nc.XcorKernelWeakLensing) -> None:
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
    nc_cmb_lens: Nc.XcorKernelCMBLensing,
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

    ell = 79.0

    z_a, _, H_Mpc_a, Wchi_a = compute_kernel(ccl_cmb_lens, nc_cosmo_eh_linear, ell)
    nc_cmb_lens.prepare(cosmo)

    nc_Wchi_a = (
        np.array([nc_cmb_lens.eval_full(cosmo, z, dist, int(ell)) for z in z_a])
        * H_Mpc_a
    )
    assert_allclose(nc_Wchi_a, Wchi_a, rtol=reltol_target, atol=0.0)


def test_xcor_set_get(nc_cosmo_default: ncpy.Cosmology) -> None:
    """Test get_reltol and set_reltol."""
    nc_xcor = Nc.Xcor.new(
        nc_cosmo_default.dist, nc_cosmo_default.ps_ml, Nc.XcorMethod.GSL
    )

    nc_xcor.set_reltol(1.0e-4)
    assert nc_xcor.get_reltol() == 1.0e-4
    assert nc_xcor.props.reltol == 1.0e-4

    nc_xcor.set_reltol(1.0e-5)
    assert nc_xcor.get_reltol() == 1.0e-5
    assert nc_xcor.props.reltol == 1.0e-5

    nc_xcor.props.reltol = 1.0e-6
    assert nc_xcor.get_reltol() == 1.0e-6

    nc_xcor.props.meth = Nc.XcorMethod.GSL
    assert nc_xcor.props.meth == Nc.XcorMethod.GSL

    nc_xcor.props.meth = Nc.XcorMethod.CUBATURE
    assert nc_xcor.props.meth == Nc.XcorMethod.CUBATURE

    assert nc_xcor.props.power_spec is nc_cosmo_default.ps_ml


def test_cmb_lens_auto_integrand(
    ccl_cosmo_eh_linear: pyccl.Cosmology,
    nc_cosmo_eh_linear: ncpy.Cosmology,
    ccl_cmb_lens: pyccl.CMBLensingTracer,
    nc_cmb_lens: Nc.XcorKernelCMBLensing,
) -> None:
    """Compare NumCosmo and CCL transfer functions."""
    cosmo = nc_cosmo_eh_linear.cosmo
    dist = nc_cosmo_eh_linear.dist
    ps_ml = nc_cosmo_eh_linear.ps_ml
    if ccl_cosmo_eh_linear["Omega_k"] != 0.0:
        pytest.skip("CMB lensing not implemented for non-flat cosmologies")
    if ccl_cosmo_eh_linear.high_precision:
        reltol_W: float = 1.0e-4
        reltol_ps = 1.0e-3
        reltol_f = 1.0e-3
    else:
        reltol_W = 1.0e-4
        reltol_ps = 1.0e-2
        reltol_f = 1.0e-2

    psp = ccl_cosmo_eh_linear.get_linear_power()

    xcor = Nc.Xcor.new(dist, ps_ml, Nc.XcorMethod.CUBATURE)
    xcor.prepare(cosmo)
    nc_cmb_lens.prepare(cosmo)

    ell = 77.0
    z_a, chi_a, H_Mpc_a, Wchi_a = compute_kernel(ccl_cmb_lens, nc_cosmo_eh_linear, ell)

    nu = ell + 0.5
    a_a = ccl_cosmo_eh_linear.scale_factor_of_chi(chi_a)
    k_a = nu / chi_a
    z_a = 1.0 / a_a - 1.0

    nc_W = (
        np.array([nc_cmb_lens.eval_full(cosmo, z, dist, int(ell)) for z in z_a])
        * H_Mpc_a
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
    nc_cmb_lens: Nc.XcorKernelCMBLensing,
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

    xcor = Nc.Xcor.new(dist, ps_ml, Nc.XcorMethod.CUBATURE)
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
    nc_cmb_isw: Nc.XcorKernelCMBISW,
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

    ell = 77.0
    z_a, _, H_Mpc_a, Wchi_a = compute_kernel(ccl_cmb_isw, nc_cosmo_eh_linear, ell)
    nc_cmb_isw.prepare(cosmo)

    nc_Wchi_a = (
        np.array([nc_cmb_isw.eval_full(cosmo, z, dist, int(ell)) for z in z_a])
        * H_Mpc_a
    )
    assert_allclose(nc_Wchi_a, Wchi_a, rtol=reltol_target, atol=0.0)


def test_tsz_kernel(
    ccl_cosmo_eh_linear: pyccl.Cosmology,
    nc_cosmo_eh_linear: ncpy.Cosmology,
    ccl_tsz: pyccl.tSZTracer,
    nc_tsz: Nc.XcorKerneltSZ,
) -> None:
    """Compare NumCosmo and CCL correlation windows."""
    cosmo = nc_cosmo_eh_linear.cosmo
    dist = nc_cosmo_eh_linear.dist
    if ccl_cosmo_eh_linear.high_precision:
        reltol_target: float = 1.0e-4
    else:
        reltol_target = 1.0e-4

    ell = 77.0
    z_a, _, H_Mpc_a, Wchi_a = compute_kernel(ccl_tsz, nc_cosmo_eh_linear, ell)
    nc_tsz.prepare(cosmo)

    nc_Wchi_a = (
        np.array([nc_tsz.eval_full(cosmo, z, dist, int(ell)) for z in z_a]) * H_Mpc_a
    )
    assert_allclose(nc_Wchi_a, Wchi_a, rtol=reltol_target, atol=0.0)


def test_weak_lensing_kernel(
    ccl_cosmo_eh_linear: pyccl.Cosmology,
    nc_cosmo_eh_linear: ncpy.Cosmology,
    ccl_weak_lensing: pyccl.WeakLensingTracer,
    nc_weak_lensing: Nc.XcorKernelWeakLensing,
) -> None:
    """Compare NumCosmo and CCL correlation windows."""
    cosmo = nc_cosmo_eh_linear.cosmo
    dist = nc_cosmo_eh_linear.dist
    if ccl_cosmo_eh_linear.high_precision:
        reltol_target: float = 1.0e-4
    else:
        reltol_target = 1.0e-4

    ell = 77.0
    z_a, _, H_Mpc_a, Wchi_a = compute_kernel(ccl_weak_lensing, nc_cosmo_eh_linear, ell)
    nc_weak_lensing.prepare(cosmo)

    nc_Wchi_a = (
        np.array([nc_weak_lensing.eval_full(cosmo, z, dist, int(ell)) for z in z_a])
        * H_Mpc_a
    )
    assert_allclose(nc_Wchi_a, Wchi_a, rtol=reltol_target, atol=1.0e-20)


def test_gal_kernel(
    ccl_cosmo_eh_linear: pyccl.Cosmology,
    nc_cosmo_eh_linear: ncpy.Cosmology,
    ccl_gal: pyccl.NumberCountsTracer,
    nc_gal: Nc.XcorKernelGal,
) -> None:
    """Compare NumCosmo and CCL correlation windows."""
    cosmo = nc_cosmo_eh_linear.cosmo
    dist = nc_cosmo_eh_linear.dist
    if ccl_cosmo_eh_linear.high_precision:
        reltol_target: float = 1.0e-4
    else:
        reltol_target = 1.0e-4

    ell = 77.0
    z_a, _, H_Mpc_a, Wchi_a = compute_kernel(ccl_gal, nc_cosmo_eh_linear, ell)
    nc_gal.prepare(cosmo)

    nc_Wchi_a = (
        np.array([nc_gal.eval_full(cosmo, z, dist, int(ell)) for z in z_a]) * H_Mpc_a
    )
    assert_allclose(nc_Wchi_a, Wchi_a, rtol=reltol_target, atol=1.0e-30)


@pytest.mark.parametrize("n_points", [None, 400])
def test_compare_kernels(
    ccl_cosmo_eh_linear: pyccl.Cosmology,
    nc_cosmo_eh_linear: ncpy.Cosmology,
    n_points: int | None,
) -> None:
    """Compare CMB lensing kernel from CCL and NumCosmo."""
    # TODO: Test curvature models when CCL supports them
    if ccl_cosmo_eh_linear["Omega_k"] == 0.0:
        cmp = nc_cmp.compare_cmb_lens_kernel(
            ccl_cosmo_eh_linear, nc_cosmo_eh_linear, ell=77, n_samples=n_points
        )
        assert_allclose(cmp.y1, cmp.y2, rtol=1.0e-6)

    cmp = nc_cmp.compare_cmb_isw_kernel(
        ccl_cosmo_eh_linear, nc_cosmo_eh_linear, ell=77, n_chi=n_points
    )
    assert_allclose(cmp.y1, cmp.y2, rtol=1.0e-1)

    cmp = nc_cmp.compare_tsz_kernel(
        ccl_cosmo_eh_linear, nc_cosmo_eh_linear, ell=77, n_chi=n_points
    )
    assert_allclose(cmp.y1, cmp.y2, rtol=1.0e-8)

    cmp = nc_cmp.compare_galaxy_weak_lensing_kernel(
        ccl_cosmo_eh_linear, nc_cosmo_eh_linear, ell=77, n_samples=n_points
    )
    assert_allclose(cmp.y1, cmp.y2, rtol=1.0e-3, atol=1.0e-15)

    cmp = nc_cmp.compare_galaxy_number_count_kernel(
        ccl_cosmo_eh_linear, nc_cosmo_eh_linear, ell=77
    )
    assert_allclose(cmp.y1, cmp.y2, rtol=1.0e-4, atol=1.0e-12)


@pytest.mark.parametrize("n_points", [None, 400])
def test_compare_autocorrelation(
    ccl_cosmo_eh_linear: pyccl.Cosmology,
    nc_cosmo_eh_linear: ncpy.Cosmology,
    n_points: int | None,
) -> None:
    """Compare CCL and NumCosmo auto-correlation."""
    ells = np.arange(2, 1000)

    # TODO: Test curvature models when CCL supports them
    if ccl_cosmo_eh_linear["Omega_k"] == 0.0:
        cmp = nc_cmp.compare_cmb_len_auto(
            ccl_cosmo_eh_linear, nc_cosmo_eh_linear, ells, n_samples=n_points
        )
        assert_allclose(cmp.y1, cmp.y2, rtol=1.0e-2)

    cmp = nc_cmp.compare_cmb_isw_auto(
        ccl_cosmo_eh_linear, nc_cosmo_eh_linear, ells, n_chi=n_points
    )
    assert_allclose(cmp.y1, cmp.y2, rtol=1.0e-1)

    cmp = nc_cmp.compare_tsz_auto(
        ccl_cosmo_eh_linear, nc_cosmo_eh_linear, ells, n_chi=n_points
    )
    assert_allclose(cmp.y1, cmp.y2, rtol=1.0e-3)

    cmp = nc_cmp.compare_galaxy_weak_lensing_auto(
        ccl_cosmo_eh_linear, nc_cosmo_eh_linear, ells, n_samples=n_points
    )
    assert_allclose(cmp.y1, cmp.y2, rtol=1.0e-3)

    cmp = nc_cmp.compare_galaxy_number_count_auto(
        ccl_cosmo_eh_linear, nc_cosmo_eh_linear, ells
    )
    assert_allclose(cmp.y1, cmp.y2, rtol=1.0e-4)


@pytest.mark.parametrize(
    "k1, k2",
    it.combinations_with_replacement(
        [
            lf("nc_cmb_lens"),
            lf("nc_cmb_isw"),
            lf("nc_gal"),
            lf("nc_tsz"),
            lf("nc_weak_lensing"),
        ],
        r=2,
    ),
)
@pytest.mark.parametrize(
    "ccl_cosmo_eh_linear",
    [pytest.param((False, 0), id="high_prec_false_index_0")],
    indirect=True,
)
def test_xcor_methods(
    nc_cosmo_eh_linear: ncpy.Cosmology, k1: Nc.XcorKernel, k2: Nc.XcorKernel
) -> None:
    """Compare NumCosmo Xcor integration methods."""
    xcor_gsl = Nc.Xcor.new(
        nc_cosmo_eh_linear.dist, nc_cosmo_eh_linear.ps_ml, Nc.XcorMethod.GSL
    )

    xcor_cub = Nc.Xcor.new(
        nc_cosmo_eh_linear.dist, nc_cosmo_eh_linear.ps_ml, Nc.XcorMethod.CUBATURE
    )

    xcor_gsl.prepare(nc_cosmo_eh_linear.cosmo)
    xcor_cub.prepare(nc_cosmo_eh_linear.cosmo)

    k1.prepare(nc_cosmo_eh_linear.cosmo)
    k2.prepare(nc_cosmo_eh_linear.cosmo)

    lmin = 2
    lmax = 1000
    vp_gsl = Ncm.Vector.new(lmax - lmin + 1)
    vp_cub = Ncm.Vector.new(lmax - lmin + 1)

    xcor_gsl.set_reltol(1.0e-7)
    xcor_gsl.limber(k1, k2, nc_cosmo_eh_linear.cosmo, lmin, lmax, vp_gsl)
    xcor_cub.limber(k1, k2, nc_cosmo_eh_linear.cosmo, lmin, lmax, vp_cub)

    vp_gsl_a = np.array(vp_gsl.dup_array())
    vp_cub_a = np.array(vp_cub.dup_array())

    assert_allclose(vp_gsl_a, vp_cub_a, rtol=1.0e-5)
