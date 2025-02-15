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


TINKER08_DELTAS_CRITICAL = [200, 300, 400, 500, 600, 800]
TINKER08_DELTAS_MATTER = [200, 300, 400, 500, 600, 800, 1200, 1600, 2400, 3200]


@pytest.fixture(
    name="tinker08_mdef",
    params=[
        ("matter", Nc.MultiplicityFuncMassDef.MEAN, Delta)
        for Delta in TINKER08_DELTAS_MATTER
    ]
    + [
        ("critical", Nc.MultiplicityFuncMassDef.CRITICAL, Delta)
        for Delta in TINKER08_DELTAS_CRITICAL
    ],
    ids=[f"matter_{Delta}" for Delta in TINKER08_DELTAS_MATTER]
    + [f"critical_{Delta}" for Delta in TINKER08_DELTAS_CRITICAL],
)
def fixture_tinker08_mdef(request) -> tuple[str, Nc.MultiplicityFuncMassDef, int]:
    """Fixture for mass functions and power spectra."""
    return request.param


TINKER10_DELTAS = [200, 300, 400, 600, 800, 1200, 1600, 2400, 3200]


@pytest.fixture(
    name="tinker10_mdef",
    params=[
        ("matter", Nc.MultiplicityFuncMassDef.MEAN, Delta) for Delta in TINKER10_DELTAS
    ],
    ids=[f"mean_{Delta}" for Delta in TINKER10_DELTAS],
)
def fixture_tinker10_mdef(request) -> tuple[str, Nc.MultiplicityFuncMassDef, int]:
    """Fixture for mass functions and power spectra."""
    return request.param


BOCQUET_HYDRO = [
    (True, Nc.MultiplicityFuncBocquetSim.HYDRO),
    (False, Nc.MultiplicityFuncBocquetSim.DM),
]


@pytest.fixture(
    name="bocquet_mdef",
    params=[
        ("matter", Nc.MultiplicityFuncMassDef.MEAN) + hydro_opt + (200,)
        for hydro_opt in BOCQUET_HYDRO
    ]
    + [
        ("critical", Nc.MultiplicityFuncMassDef.CRITICAL) + hydro_opt + (200,)
        for hydro_opt in BOCQUET_HYDRO
    ]
    + [
        ("critical", Nc.MultiplicityFuncMassDef.CRITICAL) + hydro_opt + (500,)
        for hydro_opt in BOCQUET_HYDRO
    ],
    ids=[f"matter_hydro_{hydro_opt[0]}_200" for hydro_opt in BOCQUET_HYDRO]
    + [f"critical_hydro_{hydro_opt[0]}_200" for hydro_opt in BOCQUET_HYDRO]
    + [f"critical_hydro_{hydro_opt[0]}_500" for hydro_opt in BOCQUET_HYDRO],
)
def fixture_bocquet_mdef(
    request,
) -> tuple[str, Nc.MultiplicityFuncMassDef, bool, Nc.MultiplicityFuncBocquetSim, int]:
    """Fixture for mass functions and power spectra."""
    return request.param


WATSON_OPTS = [
    ("matter", Nc.MultiplicityFuncMassDef.MEAN, 178),
    ("matter", Nc.MultiplicityFuncMassDef.MEAN, 200),
    ("matter", Nc.MultiplicityFuncMassDef.MEAN, 500),
    ("matter", Nc.MultiplicityFuncMassDef.FOF, "fof"),
]


@pytest.fixture(
    name="watson_mdef",
    params=WATSON_OPTS,
    ids=[f"{mdef}_{Delta}" for mdef, _, Delta in WATSON_OPTS],
)
def fixture_watson_mdef(request) -> tuple[str, Nc.MultiplicityFuncMassDef, int | str]:
    """Fixture for mass functions and power spectra."""
    return request.param


DESPALI_DELTAS = [200, 300, 400, 500, 600]
DESPALI_OPT = (
    [
        ("matter", Nc.MultiplicityFuncMassDef.MEAN, True, Delta)
        for Delta in DESPALI_DELTAS
    ]
    + [
        ("matter", Nc.MultiplicityFuncMassDef.MEAN, False, Delta)
        for Delta in DESPALI_DELTAS
    ]
    + [
        ("critical", Nc.MultiplicityFuncMassDef.CRITICAL, True, Delta)
        for Delta in DESPALI_DELTAS
    ]
    + [
        ("critical", Nc.MultiplicityFuncMassDef.CRITICAL, False, Delta)
        for Delta in DESPALI_DELTAS
    ]
    + [
        ("critical", Nc.MultiplicityFuncMassDef.VIRIAL, True, "vir"),
        ("critical", Nc.MultiplicityFuncMassDef.VIRIAL, False, "vir"),
    ]
)


@pytest.fixture(
    name="despali_mdef",
    params=DESPALI_OPT,
    ids=[f"{mdef}_Ellip_{ellip}_{Delta}" for mdef, _, ellip, Delta in DESPALI_OPT],
)
def fixture_despali_mdef(
    request,
) -> tuple[str, Nc.MultiplicityFuncMassDef, bool, int | str]:
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
    mf_PS = Nc.HaloMassFunction.new(cosmo_nc.dist, cosmo_nc.psf_tophat, hmf_PS)

    return ccl_hmf_PS, mf_PS


@pytest.fixture(name="massfunc_st")
def fixture_massfunc_st(
    cosmologies: tuple[pyccl.Cosmology, ncpy.Cosmology], hmd_fof: pyccl.halos.MassDef
) -> tuple[pyccl.halos.MassFunc, Nc.HaloMassFunction]:
    """Fixture for mass functions and power spectra."""
    _, cosmo_nc = cosmologies
    ccl_hmf_ST = pyccl.halos.MassFuncSheth99(mass_def=hmd_fof)
    hmf_ST = Nc.MultiplicityFuncST.new()
    mf_ST = Nc.HaloMassFunction.new(cosmo_nc.dist, cosmo_nc.psf_tophat, hmf_ST)

    return ccl_hmf_ST, mf_ST


@pytest.fixture(name="massfunc_jenkins")
def fixture_massfunc_jenkins(
    cosmologies: tuple[pyccl.Cosmology, ncpy.Cosmology]
) -> tuple[pyccl.halos.MassFunc, Nc.HaloMassFunction]:
    """Fixture for mass functions and power spectra."""
    _, cosmo_nc = cosmologies
    ccl_hmf_jenkins = pyccl.halos.MassFuncJenkins01()
    hmf_jenkins = Nc.MultiplicityFuncJenkins.new()
    mf_Jenkins = Nc.HaloMassFunction.new(
        cosmo_nc.dist, cosmo_nc.psf_tophat, hmf_jenkins
    )
    return ccl_hmf_jenkins, mf_Jenkins


@pytest.fixture(name="massfunc_tinker08")
def fixture_massfunc_tinker08(
    cosmologies: tuple[pyccl.Cosmology, ncpy.Cosmology],
    tinker08_mdef: tuple[str, Nc.MultiplicityFuncMassDef, int],
) -> tuple[pyccl.halos.MassFunc, Nc.HaloMassFunction]:
    """Fixture for mass functions and power spectra."""
    _, cosmo_nc = cosmologies
    rho_type, nc_rho_type, Delta = tinker08_mdef
    hmd = pyccl.halos.MassDef(Delta, rho_type)
    ccl_hmf_tinker = pyccl.halos.MassFuncTinker08(mass_def=hmd)
    hmf_tinker = Nc.MultiplicityFuncTinker.new()
    hmf_tinker.set_mdef(nc_rho_type)
    hmf_tinker.set_Delta(Delta)
    hmf_tinker.set_linear_interp(True)
    mf_tinker = Nc.HaloMassFunction.new(cosmo_nc.dist, cosmo_nc.psf_tophat, hmf_tinker)

    return ccl_hmf_tinker, mf_tinker


@pytest.fixture(name="massfunc_tinker10")
def fixture_massfunc_tinker10(
    cosmologies: tuple[pyccl.Cosmology, ncpy.Cosmology],
    tinker10_mdef: tuple[str, Nc.MultiplicityFuncMassDef, int],
) -> tuple[pyccl.halos.MassFunc, Nc.HaloMassFunction]:
    """Fixture for mass functions and power spectra."""
    _, cosmo_nc = cosmologies
    rho_type, nc_rho_type, Delta = tinker10_mdef
    hmd = pyccl.halos.MassDef(Delta, rho_type)
    ccl_hmf_tinker = pyccl.halos.MassFuncTinker10(mass_def=hmd)
    hmf_tinker = Nc.MultiplicityFuncTinkerMeanNormalized.new()
    hmf_tinker.set_mdef(nc_rho_type)
    hmf_tinker.set_Delta(Delta)
    mf_tinker = Nc.HaloMassFunction.new(cosmo_nc.dist, cosmo_nc.psf_tophat, hmf_tinker)

    return ccl_hmf_tinker, mf_tinker


@pytest.fixture(name="massfunc_bocquet")
def fixture_massfunc_bocquet(
    cosmologies: tuple[pyccl.Cosmology, ncpy.Cosmology],
    bocquet_mdef: tuple[
        str, Nc.MultiplicityFuncMassDef, bool, Nc.MultiplicityFuncBocquetSim, int
    ],
) -> tuple[pyccl.halos.MassFunc, Nc.HaloMassFunction]:
    """Fixture for mass functions and power spectra."""
    _, cosmo_nc = cosmologies
    rho_type, nc_rho_type, ccl_hydro, nc_hydro, Delta = bocquet_mdef
    hmd = pyccl.halos.MassDef(Delta, rho_type)
    ccl_hmf_bocquet = pyccl.halos.MassFuncBocquet16(hydro=ccl_hydro, mass_def=hmd)
    hmf_bocquet = Nc.MultiplicityFuncBocquet.new_full(nc_rho_type, nc_hydro, Delta)
    mf_bocquet = Nc.HaloMassFunction.new(
        cosmo_nc.dist, cosmo_nc.psf_tophat, hmf_bocquet
    )

    return ccl_hmf_bocquet, mf_bocquet


@pytest.fixture(name="massfunc_watson")
def fixture_massfunc_watson(
    cosmologies: tuple[pyccl.Cosmology, ncpy.Cosmology],
    watson_mdef: tuple[str, Nc.MultiplicityFuncMassDef, int | str],
) -> tuple[pyccl.halos.MassFunc, Nc.HaloMassFunction]:
    """Fixture for mass functions and power spectra."""
    _, cosmo_nc = cosmologies
    rho_type, nc_rho_type, Delta = watson_mdef
    hmd = pyccl.halos.MassDef(Delta, rho_type)
    ccl_hmf_watson = pyccl.halos.MassFuncWatson13(mass_def=hmd)
    hmf_watson = Nc.MultiplicityFuncWatson.new()
    hmf_watson.set_mdef(nc_rho_type)
    if nc_rho_type != Nc.MultiplicityFuncMassDef.FOF:
        assert isinstance(Delta, int)
        hmf_watson.set_Delta(Delta)
    mf_watson = Nc.HaloMassFunction.new(cosmo_nc.dist, cosmo_nc.psf_tophat, hmf_watson)

    return ccl_hmf_watson, mf_watson


@pytest.fixture(name="massfunc_despali")
def fixture_massfunc_despali(
    cosmologies: tuple[pyccl.Cosmology, ncpy.Cosmology],
    despali_mdef: tuple[str, Nc.MultiplicityFuncMassDef, bool, int | str],
) -> tuple[pyccl.halos.MassFunc, Nc.HaloMassFunction]:
    """Fixture for mass functions and power spectra."""
    _, cosmo_nc = cosmologies
    rho_type, nc_rho_type, ellip, Delta = despali_mdef
    hmd = pyccl.halos.MassDef(Delta, rho_type)
    ccl_hmf_despali = pyccl.halos.MassFuncDespali16(ellipsoidal=ellip, mass_def=hmd)
    hmf_despali = Nc.MultiplicityFuncDespali.new()
    hmf_despali.set_mdef(nc_rho_type)
    hmf_despali.set_eo(ellip)
    if nc_rho_type != Nc.MultiplicityFuncMassDef.VIRIAL:
        assert isinstance(Delta, int)
        hmf_despali.set_Delta(Delta)
    mf_despali = Nc.HaloMassFunction.new(
        cosmo_nc.dist, cosmo_nc.psf_tophat, hmf_despali
    )

    return ccl_hmf_despali, mf_despali


@pytest.fixture(name="mass_and_z_array")
def fixture_mass_and_z_array() -> tuple[np.ndarray, np.ndarray]:
    """Fixture for mass and redshift arrays."""
    # Array of masses
    m_arr = np.geomspace(1e10, 1e15, 128)

    # Array of redshifts
    z_arr = np.linspace(0.0, 1.0, 16)

    return m_arr, z_arr


def parametrized_massfunc_tests():
    """Parametrized mass function tests."""
    return pytest.mark.parametrize(
        "massfuncs",
        [
            lf("massfunc_ps"),
            lf("massfunc_st"),
            lf("massfunc_jenkins"),
            lf("massfunc_tinker08"),
            lf("massfunc_tinker10"),
            lf("massfunc_bocquet"),
            lf("massfunc_watson"),
            lf("massfunc_despali"),
        ],
    )


@parametrized_massfunc_tests()
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


@parametrized_massfunc_tests()
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


@parametrized_massfunc_tests()
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
