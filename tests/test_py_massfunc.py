#!/usr/bin/env python
#
# test_py_massfunc.py
#
# Fri Jan 03 12:41:47 2025
# Copyright  2025  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_massfunc.py
# Copyright (C) 2025 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

"""Unit tests for NumCosmo multiplicity and mass functions."""

import pytest
from pytest_lazy_fixtures import lf

import numpy as np
from numpy.testing import assert_allclose

import pyccl

import numcosmo_py.cosmology as ncpy
from numcosmo_py import Ncm, Nc
from numcosmo_py.ccl.nc_ccl import create_nc_obj, CCLParams


Ncm.cfg_init()


@pytest.fixture(name="cosmologies", scope="module")
def fixture_cosmologies() -> tuple[pyccl.Cosmology, ncpy.Cosmology]:
    """Fixture for cosmologies."""
    _ = CCLParams()

    cosmo_ccl = pyccl.Cosmology(
        Omega_c=0.262,
        Omega_b=0.049,
        Neff=3.046,
        h=0.6766,
        sigma8=0.8277,
        n_s=0.96,
        Omega_k=0.0,
        w0=-1.0,
        wa=0.0,
        T_CMB=2.7255,
        m_nu=[0.00, 0.0, 0.0],
        transfer_function="eisenstein_hu",
        matter_power_spectrum="linear",
    )

    cosmology = create_nc_obj(cosmo_ccl)

    return cosmo_ccl, cosmology


@pytest.fixture(name="hmd_200m")
def fixture_hmd_200m() -> pyccl.halos.MassDef:
    """Fixture for Delta=200m mass definition."""
    return pyccl.halos.MassDef200m


@pytest.fixture(name="hmd_200c")
def fixture_hmd_200c() -> pyccl.halos.MassDef:
    """Fixture for Delta=200c mass definition."""
    return pyccl.halos.MassDef200c


@pytest.fixture(name="hmd_500m")
def fixture_hmd_500m() -> pyccl.halos.MassDef:
    """Fixture for Delta=500m mass definition."""
    return pyccl.halos.MassDef(500, "matter")


@pytest.fixture(name="hmd_vir")
def fixture_hmd_vir() -> pyccl.halos.MassDef:
    """Fixture for virial overdensity mass definition."""
    return pyccl.halos.MassDefVir


@pytest.fixture(name="hmd_fof")
def fixture_hmd_fof() -> pyccl.halos.MassDef:
    """Fixture for FoF mass definition."""
    return pyccl.halos.MassDefFof


@pytest.fixture(
    name="massfunc_rho_type",
    params=[
        ("matter", Nc.MultiplicityFuncMassDef.MEAN),
        ("critical", Nc.MultiplicityFuncMassDef.CRITICAL),
    ],
    ids=["matter", "critical"],
)
def fixture_massfunc_rho_type(request) -> tuple[str, Nc.MultiplicityFuncMassDef]:
    """Fixture for mass functions and power spectra."""
    return request.param


@pytest.fixture(name="massfunc_Delta", params=[200, 300, 400, 500, 600, 800])
def fixture_massfunc_Delta(request) -> int:
    """Fixture for mass functions and power spectra."""
    return request.param


@pytest.fixture(name="massfunc_ps")
def fixture_massfunc_ps(
    cosmologies: tuple[pyccl.Cosmology, ncpy.Cosmology], hmd_fof: pyccl.halos.MassDef
) -> tuple[pyccl.halos.MassFunc, Nc.HaloMassFunction]:
    """Fixture for mass functions and power spectra."""
    _, cosmo_nc = cosmologies
    ccl_hmf_PS = pyccl.halos.MassFuncPress74(mass_def=hmd_fof)
    hmf_PS = Nc.MultiplicityFuncPS.new()
    mf_PS = Nc.HaloMassFunction.new(cosmo_nc.dist, cosmo_nc.psf, hmf_PS)

    return ccl_hmf_PS, mf_PS


@pytest.fixture(name="massfunc_st")
def fixture_massfunc_st(
    cosmologies: tuple[pyccl.Cosmology, ncpy.Cosmology], hmd_fof: pyccl.halos.MassDef
) -> tuple[pyccl.halos.MassFunc, Nc.HaloMassFunction]:
    """Fixture for mass functions and power spectra."""
    _, cosmo_nc = cosmologies
    ccl_hmf_ST = pyccl.halos.MassFuncSheth99(mass_def=hmd_fof)
    hmf_ST = Nc.MultiplicityFuncST.new()
    mf_ST = Nc.HaloMassFunction.new(cosmo_nc.dist, cosmo_nc.psf, hmf_ST)

    return ccl_hmf_ST, mf_ST


@pytest.fixture(name="massfunc_jenkins")
def fixture_massfunc_jenkins(
    cosmologies: tuple[pyccl.Cosmology, ncpy.Cosmology]
) -> tuple[pyccl.halos.MassFunc, Nc.HaloMassFunction]:
    """Fixture for mass functions and power spectra."""
    _, cosmo_nc = cosmologies
    ccl_hmf_jenkins = pyccl.halos.MassFuncJenkins01()
    hmf_jenkins = Nc.MultiplicityFuncJenkins.new()
    mf_Jenkins = Nc.HaloMassFunction.new(cosmo_nc.dist, cosmo_nc.psf, hmf_jenkins)
    return ccl_hmf_jenkins, mf_Jenkins


@pytest.fixture(name="massfunc_tinker08")
def fixture_massfunc_tinker08(
    cosmologies: tuple[pyccl.Cosmology, ncpy.Cosmology],
    massfunc_rho_type: tuple[str, Nc.MultiplicityFuncMassDef],
    massfunc_Delta: int,
) -> tuple[pyccl.halos.MassFunc, Nc.HaloMassFunction]:
    """Fixture for mass functions and power spectra."""
    _, cosmo_nc = cosmologies
    rho_type, nc_rho_type = massfunc_rho_type
    hmd = pyccl.halos.MassDef(massfunc_Delta, rho_type)
    ccl_hmf_tinker = pyccl.halos.MassFuncTinker08(mass_def=hmd)
    hmf_tinker = Nc.MultiplicityFuncTinker.new()
    hmf_tinker.set_mdef(nc_rho_type)
    hmf_tinker.set_Delta(massfunc_Delta)
    hmf_tinker.set_linear_interp(True)
    mf_tinker = Nc.HaloMassFunction.new(cosmo_nc.dist, cosmo_nc.psf, hmf_tinker)

    return ccl_hmf_tinker, mf_tinker


@pytest.fixture(name="massfunc_tinker10")
def fixture_massfunc_tinker10(
    cosmologies: tuple[pyccl.Cosmology, ncpy.Cosmology],
) -> tuple[pyccl.halos.MassFunc, Nc.HaloMassFunction]:
    """Fixture for mass functions and power spectra."""
    _, cosmo_nc = cosmologies
    ccl_hmf_tinker = pyccl.halos.MassFuncTinker10()
    hmf_tinker = Nc.MultiplicityFuncTinkerMeanNormalized.new()
    mf_tinker = Nc.HaloMassFunction.new(cosmo_nc.dist, cosmo_nc.psf, hmf_tinker)

    return ccl_hmf_tinker, mf_tinker


@pytest.fixture(name="mass_and_z_array")
def fixture_mass_and_z_array() -> tuple[np.ndarray, np.ndarray]:
    """Fixture for mass and redshift arrays."""
    # Array of masses
    m_arr = np.geomspace(1e10, 1e15, 128)

    # Array of redshifts
    z_arr = np.linspace(0.0, 1.0, 16)

    return m_arr, z_arr


@pytest.mark.parametrize(
    "massfuncs",
    [
        lf("massfunc_ps"),
        lf("massfunc_st"),
        lf("massfunc_jenkins"),
        lf("massfunc_tinker08"),
        lf("massfunc_tinker10"),
    ],
)
def test_compare_mf(
    cosmologies: tuple[pyccl.Cosmology, ncpy.Cosmology],
    massfuncs: tuple[pyccl.halos.MassFunc, Nc.HaloMassFunction],
    mass_and_z_array: tuple[np.ndarray, np.ndarray],
) -> None:
    """Test comparison of mass functions."""
    cosmo_ccl, cosmo_nc = cosmologies
    ccl_hmf, nc_hmf = massfuncs
    m_arr, z_arr = mass_and_z_array

    nc_hmf.set_eval_limits(cosmo_nc.cosmo, np.log(1.0e11), np.log(1.0e16), 0.0, 2.0)
    nc_hmf.prepare(cosmo_nc.cosmo)

    for z in z_arr:
        a = 1.0 / (1.0 + z)
        ccl_mf = ccl_hmf(cosmo_ccl, m_arr, a)
        ccl_nm = m_arr * ccl_mf
        nc_mf = np.array(
            [
                nc_hmf.dn_dlnM(cosmo_nc.cosmo, logm, z) * np.log(10.0)
                for logm in np.log(m_arr)
            ]
        )
        nc_nm = m_arr * nc_mf

        assert_allclose(ccl_nm, nc_nm, rtol=1.0e-2)


@pytest.mark.parametrize(
    "massfuncs",
    [
        lf("massfunc_ps"),
        lf("massfunc_st"),
        lf("massfunc_jenkins"),
        lf("massfunc_tinker08"),
        lf("massfunc_tinker10"),
    ],
)
def test_compare_mf_z0(
    cosmologies: tuple[pyccl.Cosmology, ncpy.Cosmology],
    massfuncs: tuple[pyccl.halos.MassFunc, Nc.HaloMassFunction],
    mass_and_z_array: tuple[np.ndarray, np.ndarray],
) -> None:
    """Test comparison of mass functions."""
    cosmo_ccl, cosmo_nc = cosmologies
    ccl_hmf, nc_hmf = massfuncs
    m_arr, _ = mass_and_z_array

    nc_hmf.set_eval_limits(cosmo_nc.cosmo, np.log(1.0e11), np.log(1.0e16), 0.0, 2.0)
    nc_hmf.prepare(cosmo_nc.cosmo)

    z = 0.0
    a = 1.0 / (1.0 + z)
    ccl_mf = ccl_hmf(cosmo_ccl, m_arr, a)
    ccl_nm = m_arr * ccl_mf
    nc_mf = np.array(
        [
            nc_hmf.dn_dlnM(cosmo_nc.cosmo, logm, z) * np.log(10.0)
            for logm in np.log(m_arr)
        ]
    )
    nc_nm = m_arr * nc_mf

    assert_allclose(ccl_nm, nc_nm, rtol=1.0e-4)


@pytest.mark.parametrize(
    "massfuncs",
    [
        lf("massfunc_ps"),
        lf("massfunc_st"),
        lf("massfunc_jenkins"),
        lf("massfunc_tinker08"),
        lf("massfunc_tinker10"),
    ],
)
def test_compare_mf_dup(
    cosmologies: tuple[pyccl.Cosmology, ncpy.Cosmology],
    massfuncs: tuple[pyccl.halos.MassFunc, Nc.HaloMassFunction],
    mass_and_z_array: tuple[np.ndarray, np.ndarray],
) -> None:
    """Test comparison of mass functions."""
    _, cosmo_nc = cosmologies
    _, nc_hmf_orig = massfuncs
    m_arr, z_arr = mass_and_z_array

    nc_hmf_orig.set_eval_limits(
        cosmo_nc.cosmo, np.log(1.0e11), np.log(1.0e16), 0.0, 2.0
    )
    nc_hmf_orig.prepare(cosmo_nc.cosmo)

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    nc_hmf = ser.dup_obj(nc_hmf_orig)
    assert isinstance(nc_hmf, Nc.HaloMassFunction)
    nc_hmf.prepare(cosmo_nc.cosmo)

    assert_allclose(nc_hmf.props.area, nc_hmf_orig.props.area, rtol=1.0e-14, atol=0.0)
    assert_allclose(nc_hmf.props.lnMi, nc_hmf_orig.props.lnMi, rtol=1.0e-14, atol=0.0)
    assert_allclose(nc_hmf.props.lnMf, nc_hmf_orig.props.lnMf, rtol=1.0e-14, atol=0.0)
    assert_allclose(nc_hmf.props.zi, nc_hmf_orig.props.zi, rtol=1.0e-14, atol=0.0)
    assert_allclose(nc_hmf.props.zf, nc_hmf_orig.props.zf, rtol=1.0e-14, atol=0.0)
    assert_allclose(nc_hmf.props.prec, nc_hmf_orig.props.prec, rtol=1.0e-14, atol=0.0)
    assert_allclose(nc_hmf.props.mf_lb, nc_hmf_orig.props.mf_lb, rtol=1.0e-14, atol=0.0)

    mf = nc_hmf.peek_multiplicity_function()
    mf_orig = nc_hmf_orig.peek_multiplicity_function()

    assert_allclose(mf.props.Delta, mf_orig.props.Delta, rtol=1.0e-14, atol=0.0)
    assert mf.props.mass_def == mf_orig.props.mass_def
    assert isinstance(mf, type(mf_orig))
    if isinstance(mf, Nc.MultiplicityFuncTinker):
        assert isinstance(mf_orig, Nc.MultiplicityFuncTinker)
        assert mf.props.linear_interp == mf_orig.props.linear_interp

    for z in z_arr:
        nc_mf_orig = np.array(
            [
                nc_hmf_orig.dn_dlnM(cosmo_nc.cosmo, logm, z) * np.log(10.0)
                for logm in np.log(m_arr)
            ]
        )
        nc_nm_orig = m_arr * nc_mf_orig
        nc_mf = np.array(
            [
                nc_hmf.dn_dlnM(cosmo_nc.cosmo, logm, z) * np.log(10.0)
                for logm in np.log(m_arr)
            ]
        )
        nc_nm = m_arr * nc_mf

        assert_allclose(nc_nm_orig, nc_nm, rtol=1.0e-14)
