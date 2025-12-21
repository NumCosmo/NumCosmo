#!/usr/bin/env python
#
# test_py_concentration_mass_relations.py
#
# Sat Jan 18 15:15:38 2025
# Copyright  2025  Mariana Penna-Lima
# <pennalima@unb.br>
#
# test_py_concentration_mass_relations.py
# Copyright (C) 2025 Mariana Penna-Lima <pennalima@unb.br>
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

"""Unit tests for NumCosmo concentration-mass relations."""

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


# Halo mass definitions - CCL


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


# Halo mass definitions - NumCosmo


DUFFY08_OPTS = [
    ("matter", Nc.HaloMassSummaryMassDef.MEAN, 200),
    ("critical", Nc.HaloMassSummaryMassDef.CRITICAL, 200),
]

BHATTACHARYA13_OPTS = [
    ("matter", Nc.HaloMassSummaryMassDef.MEAN, 200),
    ("critical", Nc.HaloMassSummaryMassDef.CRITICAL, 200),
]


@pytest.fixture(
    name="duffy08_mdef",
    params=DUFFY08_OPTS,
    ids=[f"{mdef}_{Delta}" for mdef, _, Delta in DUFFY08_OPTS],
)
def fixture_duffy08_mdef(request) -> tuple[str, Nc.HaloMassSummaryMassDef, int]:
    """Fixture for mass definition of concentration-mass relations."""
    return request.param


@pytest.fixture(
    name="bhattacharya13_mdef",
    params=BHATTACHARYA13_OPTS,
    ids=[f"{mdef}_{Delta}" for mdef, _, Delta in BHATTACHARYA13_OPTS],
)
def fixture_bhattacharya13_mdef(request) -> tuple[str, Nc.HaloMassSummaryMassDef, int]:
    """Fixture for mass definition of Bhattacharya13 concentration-mass relation."""
    return request.param


@pytest.fixture(name="cmr_duffy08")
def fixture_cmr_duffy08(
    duffy08_mdef: tuple[str, Nc.HaloMassSummaryMassDef, int],
) -> tuple[pyccl.halos.ConcentrationDuffy08, Nc.HaloMassSummary]:
    """Fixture for Duffy08 concentration-mass relation."""
    rho_type, nc_rho_type, Delta = duffy08_mdef
    hmd = pyccl.halos.MassDef(Delta, rho_type)
    ccl_cmr_Duffy08 = pyccl.halos.ConcentrationDuffy08(mass_def=hmd)
    cmr_Duffy08 = Nc.HaloCMDuffy08.new(nc_rho_type, Delta)

    return ccl_cmr_Duffy08, cmr_Duffy08


@pytest.fixture(name="cmr_duffy08_vir")
def fixture_cmr_duffy08_vir() -> (
    tuple[pyccl.halos.ConcentrationDuffy08, Nc.HaloMassSummary]
):
    """Fixture for Duffy08 concentration-mass relation for Virial mass definition."""
    hmd = pyccl.halos.MassDef("vir", "critical")
    ccl_cmr_Duffy08_vir = pyccl.halos.ConcentrationDuffy08(mass_def=hmd)
    cmr_Duffy08_vir = Nc.HaloCMDuffy08.new(Nc.HaloMassSummaryMassDef.VIRIAL, 200)

    return ccl_cmr_Duffy08_vir, cmr_Duffy08_vir


@pytest.fixture(name="cmr_klypin11")
def fixture_cmr_klypin11() -> (
    tuple[pyccl.halos.ConcentrationKlypin11, Nc.HaloMassSummary]
):
    """Fixture for Klypin11 concentration-mass relation."""
    hmd = pyccl.halos.MassDef("vir", "critical")  # pyccl.halos.MassDef(Delta, rho_type)
    ccl_cmr_Klypin11 = pyccl.halos.ConcentrationKlypin11(mass_def=hmd)
    cmr_Klypin11 = Nc.HaloCMKlypin11.new(Nc.HaloMassSummaryMassDef.VIRIAL, 200)

    return ccl_cmr_Klypin11, cmr_Klypin11


@pytest.fixture(name="cmr_bhattacharya13")
def fixture_cmr_bhattacharya13(
    cosmologies: tuple[pyccl.Cosmology, ncpy.Cosmology],
    bhattacharya13_mdef: tuple[str, Nc.HaloMassSummaryMassDef, int],
) -> tuple[pyccl.halos.ConcentrationBhattacharya13, Nc.HaloMassSummary]:
    """Fixture for Bhattacharya13 concentration-mass relation."""
    rho_type, nc_rho_type, Delta = bhattacharya13_mdef
    hmd = pyccl.halos.MassDef(Delta, rho_type)
    ccl_cmr_Bhattacharya13 = pyccl.halos.ConcentrationBhattacharya13(mass_def=hmd)

    # Get cosmology objects
    _, cosmology = cosmologies
    ps_ml = cosmology.ps_ml

    # Create mass function with matching mass definition
    mulf = Nc.MultiplicityFuncTinker.new()
    if nc_rho_type == Nc.HaloMassSummaryMassDef.MEAN:
        mulf.set_mdef(Nc.MultiplicityFuncMassDef.MEAN)
    else:
        mulf.set_mdef(Nc.MultiplicityFuncMassDef.CRITICAL)
    mulf.set_Delta(Delta)
    mfp = Nc.HaloMassFunction.new(cosmology.dist, cosmology.psf_tophat, mulf)

    # Get growth function from power spectrum
    assert isinstance(ps_ml, Nc.PowspecMLTransfer)
    gf = ps_ml.peek_gf()

    cmr_Bhattacharya13 = Nc.HaloCMBhattacharya13.new(nc_rho_type, Delta, mfp, gf)

    return ccl_cmr_Bhattacharya13, cmr_Bhattacharya13


@pytest.fixture(name="mass_and_z_array")
def fixture_mass_and_z_array() -> tuple[np.ndarray, np.ndarray]:
    """Fixture for mass and redshift arrays."""
    # Array of masses
    m_arr = np.geomspace(1e10, 1e15, 128)

    # Array of redshifts
    z_arr = np.linspace(0.0, 1.0, 16)

    return m_arr, z_arr


def parametrized_concentration_mass_tests():
    """Parametrized concentration-mass relation tests."""
    return pytest.mark.parametrize(
        "concentration_mass_relations",
        [
            lf("cmr_duffy08"),
            lf("cmr_duffy08_vir"),
            lf("cmr_klypin11"),
            lf("cmr_bhattacharya13"),
        ],
    )


def parametrized_all_cmr():
    """Parametrized all concentration-mass relation fixtures."""
    return pytest.mark.parametrize(
        "cmr_fixture",
        [
            lf("cmr_duffy08"),
            lf("cmr_duffy08_vir"),
            lf("cmr_klypin11"),
            lf("cmr_bhattacharya13"),
        ],
    )


@parametrized_all_cmr()
def test_cmr_basic_construction(
    cmr_fixture: tuple[pyccl.halos.Concentration, Nc.HaloMassSummary],
) -> None:
    """Test basic construction and type checking."""
    _, nc_hms = cmr_fixture

    # Check that object is created and is instance of HaloMassSummary
    assert nc_hms is not None
    assert isinstance(nc_hms, Nc.HaloMassSummary)
    assert isinstance(nc_hms, Ncm.Model)


@parametrized_all_cmr()
def test_cmr_properties(
    cmr_fixture: tuple[pyccl.halos.Concentration, Nc.HaloMassSummary],
) -> None:
    """Test property access and modification."""
    _, nc_hms = cmr_fixture

    # Test log10MDelta parameter access and modification
    original_log10m = nc_hms["log10MDelta"]
    assert isinstance(original_log10m, float)

    # Set a new value
    test_log10m = 14.5
    nc_hms["log10MDelta"] = test_log10m
    assert np.isclose(nc_hms["log10MDelta"], test_log10m)

    # Restore original value
    nc_hms["log10MDelta"] = original_log10m


@parametrized_all_cmr()
def test_cmr_mass_calculation(
    cmr_fixture: tuple[pyccl.halos.Concentration, Nc.HaloMassSummary],
) -> None:
    """Test mass calculation."""
    _, nc_hms = cmr_fixture

    # Set a test mass
    log10m_test = 14.0
    nc_hms["log10MDelta"] = log10m_test

    # Calculate mass
    mass = nc_hms.mass()

    # Check that mass is consistent with set log10M
    expected_mass = 10.0**log10m_test
    assert np.isclose(mass, expected_mass, rtol=1e-10)


@parametrized_all_cmr()
def test_cmr_concentration_evaluation(
    cosmologies: tuple[pyccl.Cosmology, ncpy.Cosmology],
    cmr_fixture: tuple[pyccl.halos.Concentration, Nc.HaloMassSummary],
) -> None:
    """Test concentration evaluation at different redshifts."""
    _, cosmo_nc = cosmologies
    _, nc_hms = cmr_fixture

    # Set a test mass
    nc_hms["log10MDelta"] = 14.0

    # Test at multiple redshifts
    z_test = [0.0, 0.5, 1.0, 2.0]

    for z in z_test:
        c = nc_hms.concentration(cosmo_nc.cosmo, z)

        # Check that concentration is positive and reasonable
        assert c > 0.0
        assert c < 100.0  # Sanity check: concentration should be reasonable
        assert isinstance(c, float)


@parametrized_all_cmr()
def test_cmr_concentration_mass_dependence(
    cosmologies: tuple[pyccl.Cosmology, ncpy.Cosmology],
    cmr_fixture: tuple[pyccl.halos.Concentration, Nc.HaloMassSummary],
) -> None:
    """Test that concentration changes with mass."""
    _, cosmo_nc = cosmologies
    _, nc_hms = cmr_fixture

    z = 0.0
    masses = [13.0, 14.0, 15.0]
    concentrations = []

    for log10m in masses:
        nc_hms["log10MDelta"] = log10m
        c = nc_hms.concentration(cosmo_nc.cosmo, z)
        concentrations.append(c)
        assert c > 0.0

    # Check that concentrations are different for different masses
    assert not np.allclose(concentrations[0], concentrations[1])
    assert not np.allclose(concentrations[1], concentrations[2])


@parametrized_all_cmr()
def test_cmr_concentration_redshift_dependence(
    cosmologies: tuple[pyccl.Cosmology, ncpy.Cosmology],
    cmr_fixture: tuple[pyccl.halos.Concentration, Nc.HaloMassSummary],
) -> None:
    """Test concentration evaluation at multiple redshifts."""
    _, cosmo_nc = cosmologies
    _, nc_hms = cmr_fixture

    nc_hms["log10MDelta"] = 14.0
    redshifts = [0.0, 0.5, 1.0, 2.0]
    concentrations = []

    for z in redshifts:
        c = nc_hms.concentration(cosmo_nc.cosmo, z)
        concentrations.append(c)
        assert c > 0.0
        assert c < 100.0
        assert np.isfinite(c)


def test_cmr_param_cdelta_property() -> None:
    """Test cDelta property for NcHaloCMParam."""
    cmr = Nc.HaloCMParam.new(Nc.HaloMassSummaryMassDef.CRITICAL, 200.0)

    # Test cDelta parameter access and modification
    original_cdelta = cmr["cDelta"]
    assert isinstance(original_cdelta, float)
    assert original_cdelta > 0.0

    # Set a new value
    test_cdelta = 5.5
    cmr["cDelta"] = test_cdelta
    assert np.isclose(cmr["cDelta"], test_cdelta)

    # Restore original value
    cmr["cDelta"] = original_cdelta


def test_cmr_bhattacharya13_mass_function_property(
    cosmologies: tuple[pyccl.Cosmology, ncpy.Cosmology],
) -> None:
    """Test mass-function property for NcHaloCMBhattacharya13."""
    _, cosmology = cosmologies
    ps_ml = cosmology.ps_ml

    # Create initial mass function
    mulf = Nc.MultiplicityFuncTinker.new()
    mulf.set_mdef(Nc.MultiplicityFuncMassDef.CRITICAL)
    mulf.set_Delta(200.0)
    mfp = Nc.HaloMassFunction.new(cosmology.dist, cosmology.psf_tophat, mulf)

    # Get growth function
    assert isinstance(ps_ml, Nc.PowspecMLTransfer)
    gf = ps_ml.peek_gf()

    cmr = Nc.HaloCMBhattacharya13.new(
        Nc.HaloMassSummaryMassDef.CRITICAL, 200.0, mfp, gf
    )

    # Test getting mass-function property
    retrieved_mfp = cmr.props.mass_function
    assert retrieved_mfp is not None
    assert isinstance(retrieved_mfp, Nc.HaloMassFunction)

    # Test getting growth-function property
    retrieved_gf = cmr.props.growth_function
    assert retrieved_gf is not None
    assert isinstance(retrieved_gf, Nc.GrowthFunc)


def test_cmr_diemer15_mass_function_property(
    cosmologies: tuple[pyccl.Cosmology, ncpy.Cosmology],
) -> None:
    """Test mass-function property for NcHaloCMDiemer15."""
    _, cosmology = cosmologies

    # Create mass function
    mulf = Nc.MultiplicityFuncTinker.new()
    mulf.set_mdef(Nc.MultiplicityFuncMassDef.CRITICAL)
    mulf.set_Delta(200.0)
    mfp = Nc.HaloMassFunction.new(cosmology.dist, cosmology.psf_tophat, mulf)

    cmr = Nc.HaloCMDiemer15.new(Nc.HaloMassSummaryMassDef.CRITICAL, 200.0, mfp)

    # Test getting mass-function property
    retrieved_mfp = cmr.props.mass_function
    assert retrieved_mfp is not None
    assert isinstance(retrieved_mfp, Nc.HaloMassFunction)


def test_cmr_prada12_mass_function_property(
    cosmologies: tuple[pyccl.Cosmology, ncpy.Cosmology],
) -> None:
    """Test mass-function property for NcHaloCMPrada12."""
    _, cosmology = cosmologies

    # Create mass function
    mulf = Nc.MultiplicityFuncTinker.new()
    mulf.set_mdef(Nc.MultiplicityFuncMassDef.CRITICAL)
    mulf.set_Delta(200.0)
    mfp = Nc.HaloMassFunction.new(cosmology.dist, cosmology.psf_tophat, mulf)

    cmr = Nc.HaloCMPrada12.new(Nc.HaloMassSummaryMassDef.CRITICAL, 200.0, mfp)

    # Test getting mass-function property
    retrieved_mfp = cmr.props.mass_function
    assert retrieved_mfp is not None
    assert isinstance(retrieved_mfp, Nc.HaloMassFunction)


@parametrized_all_cmr()
def test_cmr_serialize_duplicate(
    cosmologies: tuple[pyccl.Cosmology, ncpy.Cosmology],
    cmr_fixture: tuple[pyccl.halos.Concentration, Nc.HaloMassSummary],
) -> None:
    """Test serialization and duplication of concentration-mass relations."""
    _, cosmo_nc = cosmologies
    _, nc_hms = cmr_fixture

    # Set test values
    nc_hms["log10MDelta"] = 13.5

    # Serialize and duplicate
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    nc_hms_dup = ser.dup_obj(nc_hms)

    # Verify it's a different object
    assert nc_hms_dup is not nc_hms

    # Verify it's the correct type
    assert isinstance(nc_hms_dup, Nc.HaloMassSummary)
    assert type(nc_hms_dup) is type(nc_hms)

    # Verify parameter values are copied
    assert np.isclose(nc_hms_dup["log10MDelta"], nc_hms["log10MDelta"])

    # Verify mass and concentration calculations match
    assert np.isclose(nc_hms_dup.mass(), nc_hms.mass())

    z = 0.5
    c1 = nc_hms.concentration(cosmo_nc.cosmo, z)
    c2 = nc_hms_dup.concentration(cosmo_nc.cosmo, z)
    assert np.isclose(c1, c2)


@parametrized_concentration_mass_tests()
def test_compare_cmr(
    cosmologies: tuple[pyccl.Cosmology, ncpy.Cosmology],
    concentration_mass_relations: tuple[pyccl.halos.Concentration, Nc.HaloMassSummary],
    mass_and_z_array: tuple[np.ndarray, np.ndarray],
) -> None:
    """Test comparison of concentration-mass relations."""
    cosmo_ccl, cosmo_nc = cosmologies
    ccl_hmc, nc_hms = concentration_mass_relations
    m_arr, z_arr = mass_and_z_array

    # Bhattacharya13 requires looser tolerance due to implementation differences
    rtol = 3.0e-3 if isinstance(nc_hms, Nc.HaloCMBhattacharya13) else 1.0e-4

    for z in z_arr:
        a = 1.0 / (1.0 + z)
        ccl_cmr = ccl_hmc(cosmo_ccl, m_arr, a)
        nc_cmr_list = []
        for log10m in np.log10(m_arr):
            nc_hms["log10MDelta"] = log10m
            nc_cmr_list.append(nc_hms.concentration(cosmo_nc.cosmo, z))

        nc_cmr = np.array(nc_cmr_list)

        assert_allclose(ccl_cmr, nc_cmr, rtol=rtol)
