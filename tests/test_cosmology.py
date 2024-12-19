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


@pytest.fixture(name="cosmo")
def fixture_cosmo() -> Cosmology:
    """Fixture for the Cosmology class."""
    cosmo = Cosmology.default()
    return cosmo


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


def test_cosmology_default(cosmo: Cosmology):
    """Test the default cosmology creation."""
    assert isinstance(cosmo, Cosmology)
    assert cosmo.cosmo["Omegak"] == 0.0


def test_cosmology_properties(cosmo: Cosmology):
    """Test the properties of the Cosmology class."""
    assert isinstance(cosmo.ps_ml, Nc.PowspecML)
    assert isinstance(cosmo.ps_mnl, Nc.PowspecMNL)
    assert isinstance(cosmo.psf, Ncm.PowspecFilter)
    assert isinstance(cosmo.mset, Ncm.MSet)


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


def test_create_cosmo_prim_model_power_law(prim_model_cls):
    """Test the create_cosmo function with power-law primordial model."""
    prim_model, prim_cls = prim_model_cls
    cosmo_power_law = create_cosmo(prim_model=prim_model)

    assert isinstance(cosmo_power_law, Nc.HICosmo)
    assert isinstance(
        cosmo_power_law.peek_submodel_by_mid(Nc.HIPrimPowerLaw.id()), prim_cls
    )
