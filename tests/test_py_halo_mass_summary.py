#!/usr/bin/env python
#
# test_py_halo_mass_summary.py
#
# Sun Nov 10 19:52:17 2024
# Copyright  2024  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_halo_mass_summary.py
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

"""Tests on NcHaloMassSummary model."""

import pytest
from numpy.testing import assert_allclose

from numcosmo_py import Ncm, Nc

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


@pytest.fixture(name="halo_mass_summary", params=[Nc.HaloMCParam])
def fixture_halo_mass_summary(
    request, mass_def: Nc.HaloMassSummaryMassDef, Delta: float
) -> Nc.HaloMassSummary:
    """Fixture for HaloMassSummary."""
    return request.param(mass_def=mass_def, Delta=Delta)


@pytest.fixture(name="halo_mass_summary", params=[Nc.HaloCMKlypin11])
def fixture_halo_mass_summary(
    request, mass_def: Nc.HaloMassSummaryMassDef, Delta: float
) -> Nc.HaloMassSummary:
    """Fixture for HaloMassSummary."""
    return request.param(mass_def=mass_def, Delta=Delta)


@pytest.fixture(name="cosmo", params=[Nc.HICosmoDEXcdm])
def fixture_cosmo(request) -> Nc.HICosmo:
    """Fixture for HICosmo."""
    return request.param()


def test_halo_mass_summary_basic(
    halo_mass_summary: Nc.HaloMassSummary, cosmo: Nc.HICosmo
):
    """Test HaloMassSummary basic properties."""
    assert isinstance(halo_mass_summary, Nc.HaloMassSummary)
    assert halo_mass_summary.concentration(cosmo) > 0.0
    assert halo_mass_summary.mass() > 0.0
    assert halo_mass_summary.Delta(cosmo, 0.0) > 0.0
    assert halo_mass_summary.Delta(cosmo, 1.0) > 0.0
    assert halo_mass_summary.rho_bg(cosmo, 0.0) > 0.0
    assert halo_mass_summary.rho_bg(cosmo, 1.0) > 0.0
    assert halo_mass_summary.Delta_rho_bg(cosmo, 0.0) > 0.0
    assert halo_mass_summary.Delta_rho_bg(cosmo, 1.0) > 0.0


def test_halo_mass_summary_mass(halo_mass_summary: Nc.HaloMassSummary):
    """Test HaloMassSummary mass."""
    match halo_mass_summary:
        case Nc.HaloMCParam():
            log10MDelta = halo_mass_summary["log10MDelta"]
            mass = 10.0**log10MDelta
        case Nc.HaloCMKlypin11():
            log10MDelta = halo_mass_summary["log10MDelta"]
            mass = 10.0**log10MDelta
        case _:
            raise ValueError("Invalid HaloMassSummary type")

    assert_allclose(halo_mass_summary.mass(), mass, rtol=1e-5)


def test_halo_mass_summary_concentration(
    halo_mass_summary: Nc.HaloMassSummary, cosmo: Nc.HICosmo
):
    """Test HaloMassSummary concentration."""
    match halo_mass_summary:
        case Nc.HaloMCParam():
            cDelta = halo_mass_summary["cDelta"]
        case Nc.HaloCMKlypin11():
            cDelta = halo_mass_summary.concentration(cosmo)
        case _:
            raise ValueError("Invalid HaloMassSummary type")

    assert_allclose(halo_mass_summary.concentration(cosmo), cDelta, rtol=1e-5)
