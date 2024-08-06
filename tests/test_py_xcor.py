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

from numcosmo_py import Ncm, Nc
from numcosmo_py.ccl.nc_ccl import create_nc_obj

from .ccl_fixtures import (  # pylint: disable=unused-import # noqa: F401
    fixture_k_a,
    fixture_z_a,
    fixture_ccl_cosmo_eh_linear,
    fixture_ccl_cosmo_eh_halofit,
)

Ncm.cfg_init()


def test_cmb_lens_kernel(ccl_cosmo_eh_linear: pyccl.Cosmology) -> None:
    """Compare NumCosmo and CCL correlation windows."""
    cosmo, dist, _, _, _ = create_nc_obj(ccl_cosmo_eh_linear, dist_z_max=2000.0)
    if ccl_cosmo_eh_linear["Omega_k"] != 0.0:
        pytest.skip("CMB lensing not implemented for non-flat cosmologies")
    if ccl_cosmo_eh_linear.high_precision:
        reltol_target: float = 1.0e-7
    else:
        reltol_target = 1.0e-4

    dist.compute_inv_comoving(True)
    dist.prepare(cosmo)
    recomb = Nc.RecombSeager()

    lmax = 3000
    z_lss = dist.decoupling_redshift(cosmo)
    RH_Mpc = cosmo.RH_Mpc()
    ell = 77.0

    ccl_cmb_lens = pyccl.CMBLensingTracer(
        ccl_cosmo_eh_linear, z_source=z_lss, n_samples=10000
    )
    assert ccl_cmb_lens is not None

    Wchi_list, chi_list = ccl_cmb_lens.get_kernel()
    assert chi_list is not None
    assert Wchi_list is not None
    assert len(chi_list) == 1  # Single tracer
    assert len(Wchi_list) == 1  # Single tracer

    chi_a = np.array(chi_list[0])[1:-1]
    Wchi_a = np.array(Wchi_list[0])[1:-1] * ell * (ell + 1.0) / (ell + 0.5) ** 2

    nc_cmb_lens = Nc.XcorLimberKernelCMBLensing.new(
        dist, recomb, Ncm.Vector.new_array(np.zeros(lmax + 1).tolist())
    )
    nc_cmb_lens.prepare(cosmo)

    z_array = [dist.inv_comoving(cosmo, chi / RH_Mpc) for chi in chi_a]

    nc_Wchi_a = np.array(
        [
            nc_cmb_lens.eval_full(cosmo, z, dist, int(ell)) * cosmo.E(z) / RH_Mpc
            for z in z_array
        ]
    )
    assert_allclose(nc_Wchi_a, Wchi_a, rtol=reltol_target, atol=0.0)


def test_cmb_lens_serialization(ccl_cosmo_eh_linear: pyccl.Cosmology) -> None:
    """Check that CMB lensing tracer can be serialized."""
    cosmo, dist, _, _, _ = create_nc_obj(ccl_cosmo_eh_linear, dist_z_max=2000.0)

    lmax = 3000
    nc_cmb_lens = Nc.XcorLimberKernelCMBLensing.new(
        dist, Nc.RecombSeager(), Ncm.Vector.new_array(np.zeros(lmax + 1).tolist())
    )
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    nc_cmb_lens_dup = ser.dup_obj(nc_cmb_lens)
    assert nc_cmb_lens_dup is not None
    assert nc_cmb_lens_dup is not nc_cmb_lens
    assert isinstance(nc_cmb_lens_dup, Nc.XcorLimberKernelCMBLensing)

    nc_cmb_lens.prepare(cosmo)
    nc_cmb_lens_dup.prepare(cosmo)

    assert_allclose(
        nc_cmb_lens.eval_full(cosmo, 0.0, dist, 2),
        nc_cmb_lens_dup.eval_full(cosmo, 0.0, dist, 2),
    )


def test_cmb_auto_integrand(ccl_cosmo_eh_linear: pyccl.Cosmology) -> None:
    """Compare NumCosmo and CCL transfer functions."""
    cosmo, dist, ps_lin, _, _ = create_nc_obj(ccl_cosmo_eh_linear)
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
    recomb = Nc.RecombSeager()

    lmax = 3000
    z_lss = dist.decoupling_redshift(cosmo)
    RH_Mpc = cosmo.RH_Mpc()

    ccl_cmb_lens = pyccl.CMBLensingTracer(
        ccl_cosmo_eh_linear, z_source=z_lss, n_samples=10000
    )
    assert ccl_cmb_lens is not None
    psp = ccl_cosmo_eh_linear.get_linear_power()

    Nl_lensing = Ncm.Vector.new_array(np.zeros(lmax + 1).tolist())
    nc_cmb_lens = Nc.XcorLimberKernelCMBLensing.new(dist, recomb, Nl_lensing)
    xcor = Nc.Xcor.new(dist, ps_lin, Nc.XcorLimberMethod.GSL)
    xcor.prepare(cosmo)
    nc_cmb_lens.prepare(cosmo)

    Wchi_list, chi_list = ccl_cmb_lens.get_kernel()
    assert chi_list is not None
    assert Wchi_list is not None
    assert len(chi_list) == 1  # Single tracer
    assert len(Wchi_list) == 1  # Single tracer

    chi_a = np.array(chi_list[0])[1:-1]
    Wchi_a = np.array(Wchi_list[0])[1:-1]

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
    nc_ps = np.array([ps_lin.eval(cosmo, z, k) for z, k in zip(z_a, k_a)])
    nc_f = k_a * nc_ps * nc_W**2 / nu

    ccl_W = Wchi_a * ell * (ell + 1.0) / nu**2
    ccl_ps = np.array([psp(k, a, cosmo=ccl_cosmo_eh_linear) for k, a in zip(k_a, a_a)])
    ccl_f = k_a * ccl_ps * ccl_W**2 / nu

    assert_allclose(nc_W, ccl_W, rtol=reltol_W, atol=0.0)
    assert_allclose(nc_ps, ccl_ps, rtol=reltol_ps, atol=0.0)
    assert_allclose(nc_f, ccl_f, rtol=reltol_f, atol=0.0)


def test_cmb_auto(ccl_cosmo_eh_linear: pyccl.Cosmology) -> None:
    """Compare NumCosmo and CCL transfer functions."""
    cosmo, dist, ps_lin, _, _ = create_nc_obj(ccl_cosmo_eh_linear)
    if ccl_cosmo_eh_linear["Omega_k"] != 0.0:
        pytest.skip("CMB lensing not implemented for non-flat cosmologies")
    if ccl_cosmo_eh_linear.high_precision:
        reltol_target: float = 1.0e-2
    else:
        reltol_target = 1.0e-2

    lmax = 3000
    ells = np.arange(2, lmax + 1)
    z_lss = dist.decoupling_redshift(cosmo)

    ccl_cmb_lens = pyccl.CMBLensingTracer(
        ccl_cosmo_eh_linear, z_source=z_lss, n_samples=100000
    )
    assert ccl_cmb_lens is not None

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

    nc_cmb_lens = Nc.XcorLimberKernelCMBLensing.new(
        dist, Nc.RecombSeager(), Ncm.Vector.new_array(np.zeros(lmax + 1).tolist())
    )
    mset = Ncm.MSet.empty_new()
    mset.set(cosmo)
    mset.push(nc_cmb_lens)
    xcor = Nc.Xcor.new(dist, ps_lin, Nc.XcorLimberMethod.GSL)
    nc_cmb_lens_auto_v = Ncm.Vector.new(lmax + 1 - 2)
    xcor.prepare(cosmo)
    nc_cmb_lens.prepare(cosmo)
    xcor.limber(nc_cmb_lens, nc_cmb_lens, cosmo, 2, lmax, nc_cmb_lens_auto_v)
    nc_cmb_lens_auto = np.array(nc_cmb_lens_auto_v.dup_array())

    assert_allclose(ccl_cmb_lens_auto, nc_cmb_lens_auto, rtol=reltol_target, atol=0.0)


def test_cmb_isw_kernel(ccl_cosmo_eh_linear: pyccl.Cosmology) -> None:
    """Compare NumCosmo and CCL correlation windows."""
    cosmo, dist, ps_lin, _, _ = create_nc_obj(ccl_cosmo_eh_linear, dist_z_max=2000.0)
    if ccl_cosmo_eh_linear["Omega_k"] != 0.0:
        pytest.skip("CMB lensing not implemented for non-flat cosmologies")
    if ccl_cosmo_eh_linear.high_precision:
        reltol_target: float = 1.0e-7
    else:
        reltol_target = 1.0e-4

    dist.compute_inv_comoving(True)
    dist.prepare(cosmo)
    recomb = Nc.RecombSeager()

    lmax = 3000
    RH_Mpc = cosmo.RH_Mpc()
    ell = 77.0

    ccl_cmb_lens = pyccl.ISWTracer(ccl_cosmo_eh_linear, n_chi=10000)
    assert ccl_cmb_lens is not None

    Wchi_list, chi_list = ccl_cmb_lens.get_kernel()
    assert chi_list is not None
    assert Wchi_list is not None
    assert len(chi_list) == 1  # Single tracer
    assert len(Wchi_list) == 1  # Single tracer

    chi_a = np.array(chi_list[0])[1:-1]
    Wchi_a = np.array(Wchi_list[0])[1:-1] / (ell + 0.5) ** 2

    nc_cmb_isw = Nc.XcorLimberKernelCMBISW.new(
        dist, ps_lin, recomb, Ncm.Vector.new_array(np.zeros(lmax + 1).tolist())
    )
    nc_cmb_isw.prepare(cosmo)

    z_array = [dist.inv_comoving(cosmo, chi / RH_Mpc) for chi in chi_a]

    nc_Wchi_a = np.array(
        [
            nc_cmb_isw.eval_full(cosmo, z, dist, int(ell)) * cosmo.E(z) / RH_Mpc
            for z in z_array
        ]
    )
    assert_allclose(nc_Wchi_a, Wchi_a, rtol=reltol_target, atol=0.0)


def test_cmb_isw_serialization(ccl_cosmo_eh_linear: pyccl.Cosmology) -> None:
    """Check that CMB ISW tracer can be serialized."""
    cosmo, dist, ps_lin, _, _ = create_nc_obj(ccl_cosmo_eh_linear, dist_z_max=2000.0)

    lmax = 3000
    nc_cmb_isw = Nc.XcorLimberKernelCMBISW.new(
        dist,
        ps_lin,
        Nc.RecombSeager(),
        Ncm.Vector.new_array(np.zeros(lmax + 1).tolist()),
    )
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    nc_cmb_isw_dup = ser.dup_obj(nc_cmb_isw)
    assert nc_cmb_isw_dup is not None
    assert nc_cmb_isw_dup is not nc_cmb_isw
    assert isinstance(nc_cmb_isw_dup, Nc.XcorLimberKernelCMBISW)

    nc_cmb_isw.prepare(cosmo)
    nc_cmb_isw_dup.prepare(cosmo)

    assert_allclose(
        nc_cmb_isw.eval_full(cosmo, 0.0, dist, 2),
        nc_cmb_isw_dup.eval_full(cosmo, 0.0, dist, 2),
    )
