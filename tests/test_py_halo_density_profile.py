#!/usr/bin/env python
#
# test_py_halo_density_profile.py
#
# Sun Nov 10 19:52:17 2024
# Copyright  2024  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_halo_density_profile.py
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

"""Tests on NcHaloDensityProfile model."""

from itertools import product
import pytest
from numpy.testing import assert_allclose
import numpy as np

from numcosmo_py import Ncm, Nc
from numcosmo_py.helper import npa_to_seq

Ncm.cfg_init()


@pytest.fixture(
    name="mass_def",
    params=[
        Nc.HaloMassSummaryMassDef.CRITICAL,
        Nc.HaloMassSummaryMassDef.MEAN,
        Nc.HaloMassSummaryMassDef.VIRIAL,
    ],
    ids=["critical", "mean", "virial"],
)
def fixture_mass_def(request) -> Nc.HaloMassSummaryMassDef:
    """Fixture for HaloMassSummaryMassDef."""
    return request.param


@pytest.fixture(
    name="Delta",
    params=[200.0, 500.0],
    ids=["200", "500"],
)
def fixture_Delta(request) -> float:
    """Fixture for Delta."""
    return request.param


@pytest.fixture(name="halo_mass_summary", params=[Nc.HaloCMParam])
def fixture_halo_mass_summary(
    request, mass_def: Nc.HaloMassSummaryMassDef, Delta: float
) -> Nc.HaloMassSummary:
    """Fixture for HaloMassSummary."""
    return request.param(mass_def=mass_def, Delta=Delta)


@pytest.fixture(name="cosmo", params=[Nc.HICosmoDEXcdm])
def fixture_cosmo(request) -> Nc.HICosmo:
    """Fixture for HICosmo."""
    return request.param()


@pytest.fixture(
    name="halo_density_profile",
    params=[
        Nc.HaloDensityProfileNFW,
        Nc.HaloDensityProfileEinasto,
        Nc.HaloDensityProfileHernquist,
    ],
)
def fixture_halo_density_profile(
    request, halo_mass_summary: Nc.HaloMassSummary
) -> Nc.HaloDensityProfile:
    """Fixture for HaloDensityProfile."""
    dp = request.param
    assert issubclass(dp, Nc.HaloDensityProfile)
    return dp.new(halo_mass_summary)


def test_halo_density_profile_basic(
    halo_density_profile: Nc.HaloDensityProfile,
    cosmo: Nc.HICosmo,
    halo_mass_summary: Nc.HaloMassSummary,
):
    """Test HaloDensityProfile basic properties."""
    R_array = np.geomspace(1e-2, 1e2, 100, dtype=np.float64)
    z_array = np.linspace(0.0, 1.0, 100, dtype=np.float64)
    assert isinstance(halo_density_profile, Nc.HaloDensityProfile)
    assert halo_density_profile.peek_mass_summary() == halo_mass_summary
    for R, z in product(R_array, z_array):
        assert halo_density_profile.eval_2d_density(cosmo, R, z) > 0.0
        assert halo_density_profile.eval_cyl_mass(cosmo, R, z) > 0.0
        assert halo_density_profile.eval_density(cosmo, R, z) > 0.0
        assert halo_density_profile.eval_spher_mass_delta(cosmo, z) > 0.0
        assert halo_density_profile.eval_dl_2d_density(R) > 0.0
        assert halo_density_profile.eval_dl_cyl_mass(R) > 0.0
        assert halo_density_profile.eval_dl_density(R) > 0.0


def test_halo_density_profile_2d_density_vectorized(
    halo_density_profile: Nc.HaloDensityProfile,
    cosmo: Nc.HICosmo,
):
    """Test HaloDensityProfile 2D density vectorized."""
    R_array = np.geomspace(1e-2, 1e2, 100, dtype=np.float64)
    z_array = np.linspace(0.0, 1.0, 100, dtype=np.float64)
    twod_density_array = [
        halo_density_profile.eval_2d_density_array(
            cosmo, npa_to_seq(R_array), 1.0, 1.0, z
        )
        for z in z_array
    ]
    twod_density_ind = []
    for z, R in product(z_array, R_array):
        twod_density_ind.append(halo_density_profile.eval_2d_density(cosmo, R, z))
    twod_density_ind_array = np.array(twod_density_ind).reshape(100, 100)
    assert_allclose(twod_density_array, twod_density_ind_array)


def test_halo_density_profile_cyl_mass_vectorized(
    halo_density_profile: Nc.HaloDensityProfile,
    cosmo: Nc.HICosmo,
):
    """Test HaloDensityProfile cylindrical mass vectorized."""
    R_array = np.geomspace(1e-2, 1e2, 100, dtype=np.float64)
    z_array = np.linspace(0.0, 1.0, 100, dtype=np.float64)
    cyl_mass_array = [
        halo_density_profile.eval_cyl_mass_array(
            cosmo, npa_to_seq(R_array), 1.0, 1.0, z
        )
        for z in z_array
    ]
    cyl_mass_ind = []
    for z, R in product(z_array, R_array):
        cyl_mass_ind.append(halo_density_profile.eval_cyl_mass(cosmo, R, z))
    cyl_mass_ind_array = np.array(cyl_mass_ind).reshape(100, 100)
    assert_allclose(cyl_mass_array, cyl_mass_ind_array)


def test_halo_density_profile_density_vectorized(
    halo_density_profile: Nc.HaloDensityProfile,
    cosmo: Nc.HICosmo,
):
    """Test HaloDensityProfile density vectorized."""
    R_array = np.geomspace(1e-2, 1e2, 100, dtype=np.float64)
    z_array = np.linspace(0.0, 1.0, 100, dtype=np.float64)
    density_array = [
        halo_density_profile.eval_density_array(cosmo, npa_to_seq(R_array), 1.0, 1.0, z)
        for z in z_array
    ]
    density_ind = []
    for z, R in product(z_array, R_array):
        density_ind.append(halo_density_profile.eval_density(cosmo, R, z))
    density_ind_array = np.array(density_ind).reshape(100, 100)
    assert_allclose(density_array, density_ind_array)


def test_halo_density_profile_get_phys_limts(
    halo_density_profile: Nc.HaloDensityProfile,
    cosmo: Nc.HICosmo,
):
    """Test HaloDensityProfile get physical limits."""
    z_array = np.linspace(0.0, 1.0, 100)
    for z in z_array:
        R_min, R_max = halo_density_profile.get_phys_limts(cosmo, z)
        assert R_min > 0.0
        assert R_max > R_min


def test_halo_density_profile_rho_s(
    halo_density_profile: Nc.HaloDensityProfile,
    cosmo: Nc.HICosmo,
    halo_mass_summary: Nc.HaloMassSummary,
):
    """Test HaloDensityProfile rho_s."""
    z_array = np.linspace(0.0, 1.0, 100)
    for z in z_array:
        rho_s = halo_density_profile.rho_s(cosmo, z)
        r_s, rho_s0 = halo_density_profile.r_s_rho_s(cosmo, z)
        mass = halo_density_profile.eval_spher_mass_delta(cosmo, z)
        cDelta = halo_mass_summary.concentration(cosmo, z)
        dl_mass = halo_density_profile.eval_dl_spher_mass(cDelta)

        assert rho_s > 0.0
        assert r_s > 0.0
        assert rho_s0 > 0.0
        assert_allclose(rho_s, rho_s0)
        assert_allclose(rho_s, mass / (4.0 * np.pi * dl_mass * r_s**3.0))


def test_halo_density_profile_get_numint_splines(
    halo_density_profile: Nc.HaloDensityProfile,
):
    """Test HaloDensityProfile get numint splines."""
    x_array = np.geomspace(1.0e-2, 1.0e1, 100)

    twod_density, cyl_mass, _ = halo_density_profile.get_numint_splines()

    assert isinstance(twod_density, Ncm.Spline)
    assert isinstance(cyl_mass, Ncm.Spline)

    twod_density_array = np.exp([twod_density.eval(np.log(x)) for x in x_array])
    cyl_mass_array = np.exp([cyl_mass.eval(np.log(x)) for x in x_array])

    twod_density_array_cmp = np.array(
        [halo_density_profile.eval_numint_dl_2d_density(x) for x in x_array]
    )

    cyl_mass_array_cmp = np.array(
        [halo_density_profile.eval_numint_dl_cyl_mass(x) for x in x_array]
    )

    assert_allclose(twod_density_array, twod_density_array_cmp)
    assert_allclose(cyl_mass_array, cyl_mass_array_cmp)
