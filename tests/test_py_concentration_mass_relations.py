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


@pytest.fixture(
    name="duffy08_mdef",
    params=DUFFY08_OPTS,
    ids=[f"{mdef}_{Delta}" for mdef, _, Delta in DUFFY08_OPTS],
)
def fixture_duffy08_mdef(request) -> tuple[str, Nc.HaloMassSummaryMassDef, int]:
    """Fixture for mass definition of concentration-mass relations."""
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


@pytest.fixture(name="mass_and_z_array")
def fixture_mass_and_z_array() -> tuple[np.ndarray, np.ndarray]:
    """Fixture for mass and redshift arrays."""
    # Array of masses
    m_arr = np.geomspace(1e10, 1e15, 128)

    # Array of redshifts
    z_arr = np.linspace(0.0, 1.0, 16)

    return m_arr, z_arr


def parametrized_concmass_tests():
    """Parametrized concentration-mass relation tests."""
    return pytest.mark.parametrize(
        "concmassrels",
        [
            lf("cmr_duffy08"),
            lf("cmr_duffy08_vir"),
            lf("cmr_klypin11"),
        ],
    )


@parametrized_concmass_tests()
def test_compare_cmr(
    cosmologies: tuple[pyccl.Cosmology, ncpy.Cosmology],
    concmassrels: tuple[pyccl.halos.Concentration, Nc.HaloMassSummary],
    mass_and_z_array: tuple[np.ndarray, np.ndarray],
) -> None:
    """Test comparison of concentration-mass relations."""
    cosmo_ccl, cosmo_nc = cosmologies
    ccl_hmc, nc_hms = concmassrels
    m_arr, z_arr = mass_and_z_array

    # nc_hmf.set_eval_limits(cosmo_nc.cosmo, np.log(1.0e11), np.log(1.0e16), 0.0, 2.0)
    # nc_hmf.prepare(cosmo_nc.cosmo)

    for z in z_arr:
        a = 1.0 / (1.0 + z)
        ccl_cmr = ccl_hmc(cosmo_ccl, m_arr, a)
        nc_cmr_list = []
        for log10m in np.log10(m_arr):
            nc_hms["log10MDelta"] = log10m
            nc_cmr_list.append(nc_hms.concentration(cosmo_nc.cosmo, z))

        nc_cmr = np.array(nc_cmr_list)

        assert_allclose(ccl_cmr, nc_cmr, rtol=1.0e-4)
