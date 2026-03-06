"""Test fixtures for numcosmo_py."""

import pytest
import numpy as np
from numcosmo_py import Nc, Ncm

Ncm.cfg_init()


@pytest.fixture(name="prim")
def fixture_prim() -> Nc.HIPrim:
    """Create primordial power spectrum model."""
    prim = Nc.HIPrimPowerLaw.new()
    prim.props.n_SA = 0.967
    return prim


@pytest.fixture(name="reion")
def fixture_reion() -> Nc.HIReion:
    """Create reionization model."""
    return Nc.HIReionCamb.new()


@pytest.fixture(name="cosmo")
def fixture_cosmo(prim: Nc.HIPrim, reion: Nc.HIReion) -> Nc.HICosmo:
    """Create cosmological model with dark energy.

    Configures a flat LCDM cosmology with w = -1.
    """
    cosmo = Nc.HICosmoDEXcdm()
    cosmo.omega_x2omega_k()
    cosmo["Omegak"] = 0.0
    cosmo["H0"] = 71
    cosmo["Omegab"] = 0.0406
    cosmo["Omegac"] = 0.22
    cosmo["w"] = -1.0

    cosmo.add_submodel(prim)
    cosmo.add_submodel(reion)

    return cosmo


@pytest.fixture(name="dist")
def fixture_distribution() -> Nc.Distance:
    """Fixture for Distance."""
    return Nc.Distance.new(2000.0)


@pytest.fixture(name="psf")
def fixture_psf(cosmo: Nc.HICosmo, prim: Nc.HIPrim) -> Ncm.PowspecFilter:
    """Create power spectrum filter.

    Configures a tophat filter and normalizes to sigma8 = 0.8.
    """
    tf = Nc.TransferFuncEH()
    psml = Nc.PowspecMLTransfer.new(tf)
    psml.require_kmin(1.0e-6)
    psml.require_kmax(1.0e3)

    psf = Ncm.PowspecFilter.new(psml, Ncm.PowspecFilterType.TOPHAT)
    psf.set_best_lnr0()
    psf.prepare(cosmo)

    # Normalize to sigma8 = 0.8
    old_amplitude = np.exp(prim["ln10e10ASA"])
    prim["ln10e10ASA"] = np.log((0.8 / cosmo.sigma8(psf)) ** 2 * old_amplitude)

    return psf
