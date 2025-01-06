#!/usr/bin/env python
#
# test_cosmology.py
#
# Mon Jan 15 10:19:40 2024
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_serialize.py
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
#

"""Tests for the cosmology module."""

import pytest

from numcosmo_py import Nc, Ncm
from numcosmo_py.cosmology import Cosmology, create_cosmo, HIPrimModel

Ncm.cfg_init()


@pytest.fixture(name="cosmology")
def fixture_cosmo() -> Cosmology:
    """Fixture for the Cosmology class."""
    cosmology = Cosmology.default()
    return cosmology


@pytest.fixture(name="cosmology_minimal")
def fixture_cosmo_minimal() -> Cosmology:
    """Fixture for the Cosmology class with minimal objects."""
    cosmology = Cosmology.default_minimal()
    return cosmology


@pytest.fixture(
    name="prim_model_cls",
    params=[
        (HIPrimModel.ATAN, Nc.HIPrimAtan),
        (HIPrimModel.BPL, Nc.HIPrimBPL),
        (HIPrimModel.EXPC, Nc.HIPrimExpc),
        (HIPrimModel.POWER_LAW, Nc.HIPrimPowerLaw),
        (HIPrimModel.SBPL, Nc.HIPrimSBPL),
        (HIPrimModel.TWO_FLUIDS, Nc.HIPrimTwoFluids),
    ],
    ids=[
        HIPrimModel.ATAN,
        HIPrimModel.BPL,
        HIPrimModel.EXPC,
        HIPrimModel.POWER_LAW,
        HIPrimModel.SBPL,
        HIPrimModel.TWO_FLUIDS,
    ],
)
def fixture_prim_model_cls(request) -> tuple[HIPrimModel, type[Nc.HIPrim]]:
    """Fixture for the primordial model parameter."""
    return request.param


def test_cosmology_default(cosmology: Cosmology):
    """Test the default cosmology creation."""
    assert isinstance(cosmology, Cosmology)
    assert cosmology.cosmo["Omegak"] == 0.0


def test_cosmology_properties(cosmology: Cosmology):
    """Test the properties of the Cosmology class."""
    assert isinstance(cosmology.ps_ml, Nc.PowspecML)
    assert isinstance(cosmology.ps_mnl, Nc.PowspecMNL)
    assert isinstance(cosmology.psf_tophat, Ncm.PowspecFilter)
    assert isinstance(cosmology.mset, Ncm.MSet)


def test_create_cosmo():
    """Test the create_cosmo function."""
    cosmo = create_cosmo()
    assert isinstance(cosmo, Nc.HICosmo)
    assert cosmo["H0"] == 70.0
    assert cosmo["omegab"] == 0.022
    assert cosmo["omegac"] == 0.12


def test_create_cosmo_massive_nu():
    """Test the create_cosmo function with massive neutrinos."""
    cosmo = create_cosmo(massive_nu=True)
    assert isinstance(cosmo, Nc.HICosmo)
    assert cosmo["ENnu"] == 2.0328
    assert cosmo["massnu_0"] == 0.06


def test_create_cosmo_invalid_prim():
    """Test the create_cosmo function with invalid primordial model."""
    with pytest.raises(ValueError, match="Invalid primordial model"):
        _ = create_cosmo(prim_model="invalid")


def test_create_cosmo_prim_model_power_law(prim_model_cls):
    """Test the create_cosmo function with power-law primordial model."""
    prim_model, prim_cls = prim_model_cls
    cosmo_power_law = create_cosmo(prim_model=prim_model)

    assert isinstance(cosmo_power_law, Nc.HICosmo)
    assert isinstance(
        cosmo_power_law.peek_submodel_by_mid(Nc.HIPrimPowerLaw.id()), prim_cls
    )


def test_cosmology_missing_ps_ml(cosmology_minimal: Cosmology):
    """Test the Cosmology class with missing ps_ml."""
    with pytest.raises(AttributeError, match="Linear matter power spectrum not set."):
        _ = cosmology_minimal.ps_ml


def test_cosmology_missing_ps_mnl(cosmology_minimal: Cosmology):
    """Test the Cosmology class with missing ps_mnl."""
    with pytest.raises(
        AttributeError, match="Non-linear matter power spectrum not set."
    ):
        _ = cosmology_minimal.ps_mnl


def test_cosmology_missing_psf(cosmology_minimal: Cosmology):
    """Test the Cosmology class with missing psf."""
    with pytest.raises(AttributeError, match="Power spectrum filter not set."):
        _ = cosmology_minimal.psf_tophat
