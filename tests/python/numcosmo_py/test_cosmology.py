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


def test_cosmology_missing_dist():
    """Test the Cosmology class with missing dist."""
    cosmo = Nc.HICosmoDEXcdm()
    cosmo.omega_x2omega_k()
    cosmo["Omegak"] = 0.0
    cosmo.add_submodel(Nc.HIPrimPowerLaw.new())
    cosmo.add_submodel(Nc.HIReionCamb.new())
    dist = Nc.Distance.new(10.0)

    # Create cosmology and then manually set _dist to None to test the exception
    cosmology = Cosmology(cosmo=cosmo, dist=dist)
    # pylint: disable-next=protected-access
    cosmology._dist = None

    with pytest.raises(AttributeError, match="Distance not set."):
        _ = cosmology.dist


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
    with pytest.raises(AttributeError, match="Top-hat power spectrum filter not set."):
        _ = cosmology_minimal.psf_tophat


def test_cosmology_lazy_dist_preparation():
    """Test that dist property prepares on access."""
    # Create a cosmology without calling prepare()
    cosmo = Nc.HICosmoDEXcdm()
    cosmo.omega_x2omega_k()
    cosmo["Omegak"] = 0.0
    cosmo.add_submodel(Nc.HIPrimPowerLaw.new())
    cosmo.add_submodel(Nc.HIReionCamb.new())
    dist = Nc.Distance.new(10.0)

    # Create Cosmology object - but override prepare() to prevent automatic prep
    cosmology = Cosmology.__new__(Cosmology)
    cosmology.cosmo = cosmo
    # pylint: disable=protected-access
    cosmology._dist = dist
    cosmology._ps_ml = None
    cosmology._ps_mnl = None
    cosmology._psf_tophat = None
    cosmology._dist.compute_inv_comoving(True)
    cosmology.recomb = Nc.RecombSeager()
    cosmology._mset = Ncm.MSet.new_array([cosmo])
    # pylint: enable=protected-access

    # Access dist property should trigger preparation
    dist_obj = cosmology.dist
    assert dist_obj is not None
    assert isinstance(dist_obj, Nc.Distance)

    # Should be able to use it for calculations
    z = 1.0
    distance = dist_obj.comoving(cosmo, z)
    assert distance > 0.0


def test_cosmology_lazy_ps_ml_preparation():
    """Test that ps_ml property prepares on access."""
    cosmo = Nc.HICosmoDEXcdm()
    cosmo.omega_x2omega_k()
    cosmo["Omegak"] = 0.0
    cosmo.add_submodel(Nc.HIPrimPowerLaw.new())
    cosmo.add_submodel(Nc.HIReionCamb.new())
    dist = Nc.Distance.new(10.0)
    ps_ml = Nc.PowspecMLTransfer.new(Nc.TransferFuncEH())

    cosmology = Cosmology(cosmo=cosmo, dist=dist, ps_ml=ps_ml)

    # Access ps_ml property - should prepare if needed
    ps_ml_obj = cosmology.ps_ml
    assert ps_ml_obj is not None
    assert isinstance(ps_ml_obj, Nc.PowspecML)

    # Should be able to evaluate it
    k = 0.1
    z = 0.0
    power = ps_ml_obj.eval(cosmo, z, k)
    assert power > 0.0


def test_cosmology_lazy_ps_mnl_preparation():
    """Test that ps_mnl property prepares on access."""
    cosmo = Nc.HICosmoDEXcdm()
    cosmo.omega_x2omega_k()
    cosmo["Omegak"] = 0.0
    cosmo.add_submodel(Nc.HIPrimPowerLaw.new())
    cosmo.add_submodel(Nc.HIReionCamb.new())
    dist = Nc.Distance.new(10.0)
    ps_ml = Nc.PowspecMLTransfer.new(Nc.TransferFuncEH())
    ps_mnl = Nc.PowspecMNLHaloFit.new(ps_ml, 5.0, 1.0e-7)

    cosmology = Cosmology(cosmo=cosmo, dist=dist, ps_ml=ps_ml, ps_mnl=ps_mnl)

    # Access ps_mnl property - should prepare if needed
    ps_mnl_obj = cosmology.ps_mnl
    assert ps_mnl_obj is not None
    assert isinstance(ps_mnl_obj, Nc.PowspecMNL)

    # Should be able to evaluate it
    k = 0.1
    z = 0.0
    power = ps_mnl_obj.eval(cosmo, z, k)
    assert power > 0.0


def test_cosmology_lazy_psf_tophat_preparation():
    """Test that psf_tophat property prepares on access."""
    cosmo = Nc.HICosmoDEXcdm()
    cosmo.omega_x2omega_k()
    cosmo["Omegak"] = 0.0
    cosmo.add_submodel(Nc.HIPrimPowerLaw.new())
    cosmo.add_submodel(Nc.HIReionCamb.new())
    dist = Nc.Distance.new(10.0)
    ps_ml = Nc.PowspecMLTransfer.new(Nc.TransferFuncEH())
    psf = Ncm.PowspecFilter.new(ps_ml, Ncm.PowspecFilterType.TOPHAT)

    cosmology = Cosmology(cosmo=cosmo, dist=dist, ps_ml=ps_ml, psf_tophat=psf)

    # Access psf_tophat property - should prepare if needed
    psf_obj = cosmology.psf_tophat
    assert psf_obj is not None
    assert isinstance(psf_obj, Ncm.PowspecFilter)

    # Should be able to use it
    r = 8.0
    z = 0.0
    sigma = psf_obj.eval_sigma(z, r)
    assert sigma > 0.0


def test_cosmology_multiple_access_same_property(cosmology: Cosmology):
    """Test that accessing the same property multiple times works correctly."""
    # Access dist multiple times
    dist1 = cosmology.dist
    dist2 = cosmology.dist
    assert dist1 is dist2  # Should be the same object

    # Access ps_ml multiple times
    ps_ml1 = cosmology.ps_ml
    ps_ml2 = cosmology.ps_ml
    assert ps_ml1 is ps_ml2

    # Access ps_mnl multiple times
    ps_mnl1 = cosmology.ps_mnl
    ps_mnl2 = cosmology.ps_mnl
    assert ps_mnl1 is ps_mnl2

    # Access psf_tophat multiple times
    psf1 = cosmology.psf_tophat
    psf2 = cosmology.psf_tophat
    assert psf1 is psf2


def test_cosmology_all_properties_prepared_after_access(cosmology: Cosmology):
    """Test that all properties work correctly after being accessed."""
    # Access all properties
    dist = cosmology.dist
    ps_ml = cosmology.ps_ml
    ps_mnl = cosmology.ps_mnl
    psf = cosmology.psf_tophat

    # All should be usable for calculations
    z = 1.0
    k = 0.1
    r = 8.0

    # Test dist
    distance = dist.comoving(cosmology.cosmo, z)
    assert distance > 0.0

    # Test ps_ml
    power_ml = ps_ml.eval(cosmology.cosmo, 0.0, k)
    assert power_ml > 0.0

    # Test ps_mnl
    power_mnl = ps_mnl.eval(cosmology.cosmo, 0.0, k)
    assert power_mnl > 0.0

    # Test psf
    sigma = psf.eval_sigma(0.0, r)
    assert sigma > 0.0


def test_cosmology_prepare_method_works(cosmology_minimal: Cosmology):
    """Test that the prepare() method still works correctly."""
    # Call prepare explicitly
    cosmology_minimal.prepare()

    # After prepare, dist should be accessible
    dist = cosmology_minimal.dist
    assert dist is not None
    assert isinstance(dist, Nc.Distance)

    # Should be able to use it
    z = 1.0
    distance = dist.comoving(cosmology_minimal.cosmo, z)
    assert distance > 0.0
